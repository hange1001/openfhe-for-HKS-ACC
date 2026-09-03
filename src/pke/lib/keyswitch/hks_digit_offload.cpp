#include "keyswitch/hks_digit_offload.h"
#include "keyswitch/hks_strategy.h"
#include "scheme/ckksrns/ckksrns-cryptoparameters.h"
#include "FpgaManager.h"
#include <algorithm>
#include <chrono>
#include <limits>
#include <mutex>
#include <stdexcept>

namespace lbcrypto {
namespace {
constexpr size_t kN = 4096, kQ = 3, kP = 2, kTotal = 5, kMeta = 18;
struct State {
    std::mutex mutex;
    HKSDigitBackend backend = HKSDigitBackend::Off;
    HKSTopFunction model = nullptr;
    std::vector<uint64_t> moduli, roots;
    HKSDigitStats stats;
};
State& GetState() { static State state; return state; }

#if NATIVEINT == 64 && !defined(WITH_REDUCED_NOISE)
bool Transfer(State& s, uint8_t opcode, const std::vector<uint64_t>& in1,
              const std::vector<uint64_t>& in2, std::vector<uint64_t>& out, int alpha, int start) {
    if (s.backend == HKSDigitBackend::CModel) {
        s.model(in1.data(), in2.data(), out.data(), opcode, alpha, start);
        return true;
    }
    return FpgaManager::GetInstance().HksDigitTransfer(opcode, in1, in2, out, alpha, start);
}

bool Init(State& s, const std::vector<uint64_t>& mods, const std::vector<uint64_t>& roots) {
    if (s.moduli == mods && s.roots == roots) return true;
    s.moduli.clear();
    s.roots.clear();
    std::vector<uint64_t> fwd(kQ * 3 + MAX_LIMBS * CG_TF_SIZE, 0);
    std::vector<uint64_t> inv(kP * 3 + MAX_LIMBS * CG_TF_SIZE, 0), dummy(1, 0);
    for (size_t i = 0; i < kTotal; ++i) {
        const uint64_t q = mods[i];
        const unsigned shift = 64 - __builtin_clzll(q) + 62;
        const uint64_t mu = (static_cast<unsigned __int128>(1) << shift) / q;
        auto& buf = i < kQ ? fwd : inv;
        const size_t idx = i < kQ ? i : i - kQ, width = i < kQ ? kQ : kP;
        buf[idx] = q; buf[width + idx] = shift; buf[2 * width + idx] = mu;
        // The Host's OpenFHE-compatible NEGACYCLIC table, not the cyclic PoC table.
        MathUtils::BuildCgTwiddle_Unified(fwd.data() + 3 * kQ + i * CG_TF_SIZE,
                                          kN, q, roots[i], true);
        MathUtils::BuildCgTwiddle_Unified(inv.data() + 3 * kP + i * CG_TF_SIZE,
                                          kN, q, MathUtils::ModInverse(roots[i], q), false);
    }
    if (!Transfer(s, OP_INIT, fwd, inv, dummy, 0, 0)) return false;
    s.moduli = mods; s.roots = roots;
    ++s.stats.initializations;
    return true;
}
#endif
}  // namespace

void ConfigureHKSDigitBackend(HKSDigitBackend backend, HKSTopFunction cmodel) {
    if (backend == HKSDigitBackend::CModel && !cmodel)
        throw std::invalid_argument("CModel requires the HLS Top function");
    auto& s = GetState();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.backend = backend; s.model = cmodel;
    s.moduli.clear(); s.roots.clear();
    // Once explicitly configured, Off means CPU, not the legacy per-op FPGA
    // hooks (whose modulus-index cache may no longer match the fused context).
    FpgaManager::GetInstance().SetHksDigitOnly(true);
}
HKSDigitStats GetHKSDigitStats() {
    auto& s = GetState(); std::lock_guard<std::mutex> lock(s.mutex); return s.stats;
}
void ResetHKSDigitStats() {
    auto& s = GetState(); std::lock_guard<std::mutex> lock(s.mutex); s.stats = {};
}

bool TryHKSDigitOffload(const DCRTPoly& c,
                       const std::shared_ptr<CryptoParametersBase<DCRTPoly>>& base,
                       std::vector<DCRTPoly>& output) {
    auto& s = GetState();
    std::lock_guard<std::mutex> lock(s.mutex);
    if (s.backend == HKSDigitBackend::Off) return false;
    auto fallback = [&](const char* reason) {
        ++s.stats.fallbacks; s.stats.last_fallback = reason; return false;
    };
#if NATIVEINT != 64 || defined(WITH_REDUCED_NOISE)
    return fallback("requires NATIVEINT=64 and standard unsigned ApproxSwitchCRTBasis");
#else
    const auto cp = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(base);
    if (!cp) return fallback("only CKKS HYBRID is enabled");
    if (GetHKSStrategy() != HKSStrategy::DC) return fallback("only DC scheduling is enabled");
    if (c.GetFormat() != Format::EVALUATION || c.GetRingDimension() != kN ||
        c.GetNumOfElements() != kQ || !cp->GetParamsP() || cp->GetParamsP()->GetParams().size() != kP)
        return fallback("requires EVAL, N=4096, current Q=3, P=2");
    const size_t alpha = cp->GetNumPerPartQ();
    const size_t digits = alpha ? (kQ + alpha - 1) / alpha : 0;
    if (alpha < 1 || alpha > kQ || digits > cp->GetNumberOfQPartitions())
        return fallback("unsupported digit partition");
    auto qp = c.GetExtendedCRTBasis(cp->GetParamsP());
    std::vector<uint64_t> mods, roots;
    for (const auto& p : qp->GetParams()) {
        const uint64_t q = p->GetModulus().ConvertToInt(), root = p->GetRootOfUnity().ConvertToInt();
        if (q <= 2 || q >= (uint64_t(1) << 62) || q % (2 * kN) != 1 ||
            root == 0 || root >= q || MathUtils::Power(root, kN, q) != q - 1)
            return fallback("unsupported modulus or root");
        mods.push_back(q); roots.push_back(root);
    }
    // Check all metadata before touching the device, including complement ordering.
    for (size_t part = 0; part < digits; ++part) {
        const size_t start = alpha * part, count = std::min(alpha, kQ - start);
        const auto& src = cp->GetParamsPartQ(part)->GetParams();
        const auto& dst = cp->GetParamsComplPartQ(kQ - 1, part)->GetParams();
        const auto& inv = cp->GetPartQlHatInvModq(part, count - 1);
        const auto& weights = cp->GetPartQlHatModp(kQ - 1, part);
        if (src.size() < count || inv.size() < count || weights.size() < count || dst.size() != kTotal - count)
            return fallback("CRT metadata dimensions mismatch");
        for (size_t i = 0; i < count; ++i) {
            if (src[i]->GetModulus().ConvertToInt() != mods[start + i] ||
                src[i]->GetRootOfUnity().ConvertToInt() != roots[start + i] ||
                inv[i].ConvertToInt() >= mods[start + i] || weights[i].size() != dst.size())
                return fallback("CRT source metadata mismatch");
            for (size_t j = 0; j < dst.size(); ++j) {
                const size_t idx = j < start ? j : j + count;
                if (dst[j]->GetModulus().ConvertToInt() != mods[idx] ||
                    dst[j]->GetRootOfUnity().ConvertToInt() != roots[idx] ||
                    weights[i][j].ConvertToInt() >= mods[idx])
                    return fallback("CRT complement metadata mismatch");
            }
        }
    }
    try {
        auto t0 = std::chrono::steady_clock::now();
        if (!Init(s, mods, roots)) return fallback("XRT backend unavailable");
        std::vector<DCRTPoly> result;
        result.reserve(digits);
        for (size_t part = 0; part < digits; ++part) {
            const size_t start = alpha * part, count = std::min(alpha, kQ - start);
            std::vector<uint64_t> in(count * kN), meta(kMeta, 0);
            std::vector<uint64_t> out(kTotal * kN, std::numeric_limits<uint64_t>::max());
            const auto& inv = cp->GetPartQlHatInvModq(part, count - 1);
            const auto& weights = cp->GetPartQlHatModp(kQ - 1, part);
            for (size_t i = 0; i < count; ++i) {
                meta[15 + i] = inv[i].ConvertToInt();
                for (size_t j = 0; j < kTotal - count; ++j)
                    meta[i * kTotal + j] = weights[i][j].ConvertToInt();
                for (size_t r = 0; r < kN; ++r)
                    in[i * kN + r] = c.GetElementAtIndex(start + i)[r].ConvertToInt();
            }
            if (!Transfer(s, OP_HKS_DIGIT, in, meta, out, count, start))
                throw std::runtime_error("HKS_DIGIT transport became unavailable");
            DCRTPoly digit(qp, Format::EVALUATION, true);
            for (size_t j = 0; j < kTotal; ++j) {
                auto tower = digit.GetElementAtIndex(j);
                for (size_t r = 0; r < kN; ++r) {
                    const auto value = out[j * kN + r];
                    if (value >= mods[j])
                        throw std::runtime_error("Invalid HKS_DIGIT result (check opcode-8 xclbin and ABI)");
                    tower[r] = value;
                }
                digit.SetElementAtIndex(j, std::move(tower));
            }
            result.push_back(std::move(digit));
            ++s.stats.digits;
            s.stats.input_bytes += in.size() * 8;
            s.stats.metadata_bytes += meta.size() * 8;
            s.stats.output_bytes += out.size() * 8;
        }
        output = std::move(result);
        auto& stats = GetHKSStats();
        stats.num_digits = digits; stats.size_ql = kQ; stats.size_p = kP;
        stats.alpha = alpha; stats.ring_dim = kN;
        stats.fused_digits += digits;
        stats.time_fused_us += std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - t0).count();
        return true;
    } catch (...) {
        s.moduli.clear(); s.roots.clear();
        ++s.stats.failures;
        throw;
    }
#endif
}
}  // namespace lbcrypto
