#include "openfhe.h"
#include "keyswitch/hks_digit_offload.h"
#include "keyswitch/hks_strategy.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>

// Keep HLS macros/types out of the OpenFHE translation unit.
extern "C" void Top(const uint64_t*, const uint64_t*, uint64_t*, uint8_t, int, int);
using namespace lbcrypto;
namespace {
size_t checks = 0, compared = 0;
void Require(bool ok, const std::string& what) {
    ++checks;
    if (!ok) throw std::runtime_error(what);
}
void EqualPoly(const DCRTPoly& a, const DCRTPoly& b, const std::string& label) {
    Require(a.GetFormat() == b.GetFormat(), label + ": format");
    Require(a.GetNumOfElements() == b.GetNumOfElements(), label + ": tower count");
    for (size_t t = 0; t < a.GetNumOfElements(); ++t) {
        const auto& x = a.GetElementAtIndex(t);
        const auto& y = b.GetElementAtIndex(t);
        Require(x.GetModulus() == y.GetModulus() && x.GetRootOfUnity() == y.GetRootOfUnity(), label + ": params");
        for (size_t i = 0; i < a.GetRingDimension(); ++i) {
            ++compared;
            if (x[i] != y[i]) throw std::runtime_error(label + ": mismatch tower=" +
                std::to_string(t) + " index=" + std::to_string(i) + " CPU=" + x[i].ToString() + " fused=" + y[i].ToString());
        }
    }
}
void EqualDigits(const std::vector<DCRTPoly>& a, const std::vector<DCRTPoly>& b, const std::string& label) {
    Require(a.size() == b.size(), label + ": digits");
    for (size_t i = 0; i < a.size(); ++i) EqualPoly(a[i], b[i], label + " digit=" + std::to_string(i));
}
void EqualCipher(const Ciphertext<DCRTPoly>& a, const Ciphertext<DCRTPoly>& b, const std::string& label) {
    Require(a->GetElements().size() == b->GetElements().size(), label + ": elements");
    for (size_t i = 0; i < a->GetElements().size(); ++i) EqualPoly(a->GetElements()[i], b->GetElements()[i], label);
}
CryptoContext<DCRTPoly> Context(uint32_t n = 4096, uint32_t scale = 50, uint32_t partitions = 2) {
    CCParams<CryptoContextCKKSRNS> p;
    p.SetSecurityLevel(HEStd_NotSet);  // small TEST parameters, not a security claim
    p.SetRingDim(n); p.SetMultiplicativeDepth(1); p.SetScalingModSize(scale);
    p.SetFirstModSize(60); p.SetBatchSize(8); p.SetKeySwitchTechnique(HYBRID);
    p.SetScalingTechnique(FLEXIBLEAUTOEXT); p.SetNumLargeDigits(partitions);
    auto cc = GenCryptoContext(p);
    cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE);
    return cc;
}
DCRTPoly Input(const CryptoContext<DCRTPoly>& cc, int pattern) {
    DCRTPoly c(cc->GetElementParams(), Format::COEFFICIENT, true);
    std::mt19937_64 rng(0x4f50454e464845ULL + pattern);
    for (size_t t = 0; t < c.GetNumOfElements(); ++t) {
        auto tower = c.GetElementAtIndex(t);
        const auto q = tower.GetModulus().ConvertToInt();
        for (size_t r = 0; r < c.GetRingDimension(); ++r)
            tower[r] = pattern == 0 ? 0 : pattern == 1 ? q - 1 : pattern == 2 ? (r == 1) : rng() % q;
        c.SetElementAtIndex(t, std::move(tower));
    }
    c.SetFormat(Format::EVALUATION);
    return c;
}
auto Precompute(const CryptoContext<DCRTPoly>& cc, const DCRTPoly& c) {
    return cc->GetScheme()->EvalKeySwitchPrecomputeCore(c, cc->GetCryptoParameters());
}
void CheckFallback(const CryptoContext<DCRTPoly>& cc, const DCRTPoly& c, HKSStrategy strategy,
                   const std::string& label) {
    SetHKSStrategy(strategy);
    ConfigureHKSDigitBackend(HKSDigitBackend::Off);
    const auto ref = Precompute(cc, c);
    ConfigureHKSDigitBackend(HKSDigitBackend::CModel, Top);
    ResetHKSDigitStats();
    const auto actual = Precompute(cc, c);
    const auto s = GetHKSDigitStats();
    Require(s.digits == 0 && s.initializations == 0 && s.fallbacks == 1, label + ": no kernel launch");
    EqualDigits(*ref, *actual, label);
    std::cout << "[PASS] fallback " << label << ": " << s.last_fallback << '\n';
    SetHKSStrategy(HKSStrategy::DC);
    ConfigureHKSDigitBackend(HKSDigitBackend::Off);
}
void OldBitstream(const uint64_t* a, const uint64_t* b, uint64_t* c, uint8_t op, int n, int start) {
    if (op == 0) Top(a, b, c, op, n, start);
}
void FailSecondDigit(const uint64_t* a, const uint64_t* b, uint64_t* c, uint8_t op, int n, int start) {
    if (op == 8 && start == 2) throw std::runtime_error("injected transport failure");
    Top(a, b, c, op, n, start);
}
std::shared_ptr<DCRTPoly::Params> transformBasis;
void ProbeTransforms(const uint64_t* a, const uint64_t* b, uint64_t* c, uint8_t op, int n, int start) {
    Top(a, b, c, op, n, start);
    if (op != 0) return;
    std::mt19937_64 rng(91);
    for (size_t t = 0; t < 5; ++t) {
        NativePoly ref(transformBasis->GetParams()[t], Format::COEFFICIENT, true);
        std::vector<uint64_t> coeff(4096), eval(4096), actual(4096), inverse(4096);
        const auto q = ref.GetModulus().ConvertToInt();
        for (size_t i = 0; i < 4096; ++i) { coeff[i] = rng() % q; ref[i] = coeff[i]; }
        ref.SetFormat(Format::EVALUATION);
        for (size_t i = 0; i < 4096; ++i) eval[i] = ref[i].ConvertToInt();
        Top(coeff.data(), coeff.data(), actual.data(), 4, 1, t);
        Top(eval.data(), eval.data(), inverse.data(), 5, 1, t);
        Require(eval == actual, "OpenFHE negacyclic forward ordering, tower=" + std::to_string(t));
        Require(coeff == inverse, "OpenFHE inverse normalization, tower=" + std::to_string(t));
        compared += 2 * 4096;
    }
    std::cout << "[PASS] NTT and INTT vs OpenFHE independently, all five moduli\n";
}
void CheckDecrypt(const CryptoContext<DCRTPoly>& cc, const PrivateKey<DCRTPoly>& key,
                  const Ciphertext<DCRTPoly>& c, const std::vector<double>& expected, const std::string& label) {
    Plaintext p;
    cc->Decrypt(key, c, &p); p->SetLength(expected.size());
    const auto values = p->GetRealPackedValue();
    if (values.size() < expected.size()) throw std::runtime_error(label + ": decoded vector too short");
    double error = 0;
    for (size_t i = 0; i < expected.size(); ++i) {
        if (!std::isfinite(values[i])) throw std::runtime_error(label + ": non-finite decoded value");
        error = std::max(error, std::abs(values[i] - expected[i]));
    }
    Require(std::isfinite(error) && error < 1e-6, label + ": decryption error");
    std::cout << "[PASS] " << label << " max_abs_error=" << error << '\n';
}
}  // namespace

int main() {
    try {
        ConfigureHKSDigitBackend(HKSDigitBackend::Off);
        SetHKSStrategy(HKSStrategy::DC);
        const auto cc = Context();
        const auto cp = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
        Require(cc->GetElementParams()->GetParams().size() == 3 && cp->GetParamsP()->GetParams().size() == 2 &&
                cp->GetNumPerPartQ() == 2, "primary test must use Q=3 P=2 alpha=2");
        std::cout << "OpenFHE + actual HLS Top C-model; N=4096 Q=3 P=2 digits=[2,1]\n";
        for (int pattern = 0; pattern < 7; ++pattern) {
            ConfigureHKSDigitBackend(HKSDigitBackend::Off);
            const auto c = Input(cc, pattern), original = c;
            const auto ref = Precompute(cc, c);
            if (pattern == 3) transformBasis = c.GetExtendedCRTBasis(cp->GetParamsP());
            ConfigureHKSDigitBackend(HKSDigitBackend::CModel, pattern == 3 ? ProbeTransforms : Top);
            ResetHKSDigitStats(); ResetHKSStats();
            const auto result = Precompute(cc, c);
            EqualDigits(*ref, *result, "ModUp pattern=" + std::to_string(pattern));
            EqualPoly(c, original, "input unchanged");
            auto s = GetHKSDigitStats();
            Require(s.digits == 2 && s.fallbacks == 0 && s.initializations == 1, "fused hook must execute");
            Require(s.input_bytes == 3 * 4096 * 8 && s.metadata_bytes == 2 * 18 * 8 &&
                    s.output_bytes == 10 * 4096 * 8, "payload accounting");
            Require(GetHKSStats().fused_digits == 2 && GetHKSStats().intt_poly == 0 &&
                    GetHKSStats().bconv == 0, "fused and CPU counters separated");
            if (pattern == 3) {
                EqualDigits(*ref, *Precompute(cc, c), "cached repeat");
                s = GetHKSDigitStats();
                Require(s.initializations == 1 && s.digits == 4, "context reused without INIT");
            }
            std::cout << "[PASS] bit-exact ModUp pattern=" << pattern << " (both digits, all five towers)\n";
        }

        // A/B/A contexts exercise automatic cache invalidation, not just Configure().
        ConfigureHKSDigitBackend(HKSDigitBackend::Off);
        const auto cc2 = Context(4096, 49);
        const auto a = Input(cc, 7), b = Input(cc2, 8);
        const auto ar = Precompute(cc, a), br = Precompute(cc2, b);
        ConfigureHKSDigitBackend(HKSDigitBackend::CModel, Top); ResetHKSDigitStats();
        EqualDigits(*ar, *Precompute(cc, a), "context A");
        EqualDigits(*br, *Precompute(cc2, b), "context B");
        EqualDigits(*ar, *Precompute(cc, a), "context A again");
        Require(GetHKSDigitStats().initializations == 3 && GetHKSDigitStats().digits == 6, "context changes reinitialize");
        std::cout << "[PASS] context A/B/A reinitialization\n";

        // Full cryptographic path. CPU and C-model consume the SAME ciphertext/key.
        ConfigureHKSDigitBackend(HKSDigitBackend::Off);
        auto keys = cc->KeyGen();
        cc->EvalRotateKeyGen(keys.secretKey, {1, -1}); cc->EvalMultKeyGen(keys.secretKey);
        const std::vector<double> x{0.25, 0.5, 0.75, 1, 2, 3, 4, 5};
        auto ct = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(x));
        for (int rotation : {1, -1}) {
            ConfigureHKSDigitBackend(HKSDigitBackend::Off);
            const auto ref = cc->EvalRotate(ct, rotation);
            ConfigureHKSDigitBackend(HKSDigitBackend::CModel, Top); ResetHKSDigitStats();
            const auto result = cc->EvalRotate(ct, rotation);
            EqualCipher(ref, result, "EvalRotate");
            Require(GetHKSDigitStats().digits == 2, "EvalRotate must use fused digits");
            std::vector<double> expected(8);
            for (int i = 0; i < 8; ++i) expected[i] = x[(i + rotation + 8) % 8];
            CheckDecrypt(cc, keys.secretKey, result, expected, "EvalRotate " + std::to_string(rotation));
        }
        ConfigureHKSDigitBackend(HKSDigitBackend::Off);
        const auto mulRef = cc->EvalMult(ct, ct);
        ConfigureHKSDigitBackend(HKSDigitBackend::CModel, Top); ResetHKSDigitStats();
        const auto mulResult = cc->EvalMult(ct, ct);
        EqualCipher(mulRef, mulResult, "EvalMult");
        std::vector<double> squared;
        for (auto v : x) squared.push_back(v * v);
        CheckDecrypt(cc, keys.secretKey, mulResult, squared, "EvalMult");
        Require(GetHKSDigitStats().digits == 0 && GetHKSDigitStats().fallbacks == 1,
                "FLEXIBLEAUTOEXT EvalMult reduces Q to 2 and must fall back");
        std::cout << "  EvalMult fused_digits=" << GetHKSDigitStats().digits
                  << " fallbacks=" << GetHKSDigitStats().fallbacks << '\n';

        // Level drop, other schedule/shape, unavailable transport, invalid output.
        ConfigureHKSDigitBackend(HKSDigitBackend::Off);
        auto reduced = a; reduced.DropLastElement();
        CheckFallback(cc, reduced, HKSStrategy::DC, "Q=2 reduced level");
        CheckFallback(cc, a, HKSStrategy::MP, "MP schedule");
        CheckFallback(cc, a, HKSStrategy::OC, "OC schedule");
        const auto large = Context(8192);
        CheckFallback(large, Input(large, 3), HKSStrategy::DC, "N=8192");
        const auto three = Context(4096, 50, 3);
        CheckFallback(three, Input(three, 3), HKSStrategy::DC, "P=1");
        ConfigureHKSDigitBackend(HKSDigitBackend::Off);
        auto coeff = a; coeff.SetFormat(Format::COEFFICIENT);
        ConfigureHKSDigitBackend(HKSDigitBackend::CModel, Top); ResetHKSDigitStats();
        std::vector<DCRTPoly> unchanged{a};
        Require(!TryHKSDigitOffload(coeff, cc->GetCryptoParameters(), unchanged), "coefficient input rejected");
        Require(GetHKSDigitStats().initializations == 0 && unchanged.size() == 1, "reject before INIT");
        EqualPoly(unchanged[0], a, "unsupported input leaves output unchanged");
        bool missingModel = false;
        try { ConfigureHKSDigitBackend(HKSDigitBackend::CModel); }
        catch (const std::invalid_argument&) { missingModel = true; }
        Require(missingModel, "missing C-model callback rejected");
        ConfigureHKSDigitBackend(HKSDigitBackend::XRT); ResetHKSDigitStats();
        EqualDigits(*ar, *Precompute(cc, a), "XRT unavailable CPU fallback");
        Require(GetHKSDigitStats().digits == 0 && GetHKSDigitStats().fallbacks == 1, "unavailable XRT");

        for (auto callback : {OldBitstream, FailSecondDigit}) {
            ConfigureHKSDigitBackend(HKSDigitBackend::CModel, callback); ResetHKSDigitStats();
            std::vector<DCRTPoly> output{a};
            bool threw = false;
            try { TryHKSDigitOffload(a, cc->GetCryptoParameters(), output); }
            catch (const std::runtime_error&) { threw = true; }
            Require(threw && GetHKSDigitStats().failures == 1, "bad transport must fail closed");
            Require(output.size() == 1, "no partial digits published");
            EqualPoly(output[0], a, "failed output unchanged");
        }
        ConfigureHKSDigitBackend(HKSDigitBackend::CModel, Top);
        EqualDigits(*ar, *Precompute(cc, a), "recovery after failed transport");
        ConfigureHKSDigitBackend(HKSDigitBackend::Off);
        std::cout << "PASS: " << checks << " checks, " << compared << " residues compared exactly.\n"
                  << "This is HLS C-model integration, NOT RTL/XRT emulation or hardware performance.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
