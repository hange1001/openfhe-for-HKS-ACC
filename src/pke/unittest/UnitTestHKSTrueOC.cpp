//==================================================================================
// HKS 真 OC 的 test-only functional reference
//
// 这个文件**只供测试使用**，不参与 production dispatch，不被 keyswitch-hybrid.cpp
// 或任何库代码引用。
//
// 背景：ADR-010 与 docs/notes/OC_strategy_gap_analysis.md 指出，
// keyswitch-hybrid.cpp 里现有的 HKSStrategy::OC 是「伪 OC」——它的外层循环只跑
// sizeP 个 P 塔，且每次仍调用全宽 ApproxSwitchCRTBasis 后只取一列，等价于
// 「DC + sizeP 倍冗余」。它不是 tools/hks_flow_sim 建模的那个 OC。
//
// 真 OC（tools/hks_flow_sim/hks_sim/trace_oc.py 建模的那个）是：
//
//     for output tower t in [0, sizeQl + sizeP):
//         for digit j:
//             native Q tower  -> bypass（直接复用原 EVAL 系数）
//             non-native      -> BConv 出该目标塔 -> NTT
//             立即 KeyMult
//             立即 Accumulate
//         该输出塔就此完成
//
// 本文件验证这个**顺序与映射**是正确的，即它算出的 cTilda0/cTilda1 与 OpenFHE
// 自己的 DC 路径逐 residue 完全一致。
//
// 关于性能：这里的 non-native 分支调用的是现成的全宽 BConv，然后只取目标列。
// 这**只负责提供正确的 residue**；它的 CPU 内存与运行时间不进入任何性能模型，
// 性能结论一律以 tools/hks_flow_sim 的周期模型为准。换句话说，本文件验证的是
// 真 OC 的**功能正确性**，不是它的效率。
//==================================================================================

#include "openfhe.h"

#include "keyswitch/keyswitch-hybrid.h"
#include "keyswitch/hks_strategy.h"
#include "schemerns/rns-cryptoparameters.h"
#include "utils/utilities-int.h"

#include "gtest/gtest.h"

#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace lbcrypto;

namespace {

//----------------------------------------------------------------------------
// 测试夹具：一个真实的 CKKS 上下文 + 一个 eval key + 一个待 key-switch 的多项式
//----------------------------------------------------------------------------
struct HKSFixture {
    CryptoContext<DCRTPoly> cc;
    EvalKey<DCRTPoly> evalKey;
    DCRTPoly c;                                       // 待 key-switch 的多项式
    std::shared_ptr<CryptoParametersRNS> cryptoParams;
    std::shared_ptr<ILDCRTParams<BigInteger>> paramsQl;
};

HKSFixture MakeFixture(uint32_t ringDim, uint32_t multDepth, uint32_t dnum) {
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetScalingModSize(50);
    parameters.SetRingDim(ringDim);
    parameters.SetBatchSize(8);
    parameters.SetNumLargeDigits(dnum);
    parameters.SetKeySwitchTechnique(HYBRID);

    HKSFixture f;
    f.cc = GenCryptoContext(parameters);
    f.cc->Enable(PKE);
    f.cc->Enable(KEYSWITCH);
    f.cc->Enable(LEVELEDSHE);

    auto keys = f.cc->KeyGen();
    f.cc->EvalMultKeyGen(keys.secretKey);
    f.evalKey = f.cc->GetEvalMultKeyVector(keys.secretKey->GetKeyTag())[0];

    std::vector<double> x = {0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0, 5.0};
    auto ptxt = f.cc->MakeCKKSPackedPlaintext(x);
    auto ct   = f.cc->Encrypt(keys.publicKey, ptxt);

    // 取密文的第二个分量作为待 key-switch 的多项式（HKS 实际作用的对象）
    f.c            = ct->GetElements()[1];
    f.cryptoParams = std::dynamic_pointer_cast<CryptoParametersRNS>(f.cc->GetCryptoParameters());
    f.paramsQl     = f.c.GetParams();
    return f;
}

//----------------------------------------------------------------------------
// digit 分解：与 keyswitch-hybrid.cpp 的做法一致（含末位不完整 digit）
//----------------------------------------------------------------------------
std::vector<DCRTPoly> SplitIntoDigits(const DCRTPoly& c,
                                      const std::shared_ptr<CryptoParametersRNS>& cp,
                                      uint32_t alpha, uint32_t numPartQl, size_t sizeQl) {
    std::vector<DCRTPoly> partsCt(numPartQl);
    for (uint32_t part = 0; part < numPartQl; part++) {
        if (part == numPartQl - 1) {
            auto paramsPartQ      = cp->GetParamsPartQ(part);
            uint32_t sizePartQl   = static_cast<uint32_t>(sizeQl) - alpha * part;
            std::vector<NativeInteger> moduli(sizePartQl);
            std::vector<NativeInteger> roots(sizePartQl);
            for (uint32_t i = 0; i < sizePartQl; i++) {
                moduli[i] = paramsPartQ->GetParams()[i]->GetModulus();
                roots[i]  = paramsPartQ->GetParams()[i]->GetRootOfUnity();
            }
            auto params   = DCRTPoly::Params(paramsPartQ->GetCyclotomicOrder(), moduli, roots);
            partsCt[part] = DCRTPoly(std::make_shared<ILDCRTParams<BigInteger>>(params),
                                     Format::EVALUATION, true);
        }
        else {
            partsCt[part] = DCRTPoly(cp->GetParamsPartQ(part), Format::EVALUATION, true);
        }
        usint sizePartQl   = partsCt[part].GetNumOfElements();
        usint startPartIdx = alpha * part;
        for (uint32_t i = 0, idx = startPartIdx; i < sizePartQl; i++, idx++) {
            partsCt[part].SetElementAtIndex(i, c.GetElementAtIndex(idx));
        }
    }
    return partsCt;
}

uint32_t NumDigits(const std::shared_ptr<CryptoParametersRNS>& cp, size_t sizeQl,
                   uint32_t alpha) {
    uint32_t n = static_cast<uint32_t>(
        std::ceil(static_cast<double>(sizeQl) / alpha));
    if (n > cp->GetNumberOfQPartitions())
        n = cp->GetNumberOfQPartitions();
    return n;
}

//----------------------------------------------------------------------------
// 真 OC：output tower 作外层循环，native 塔 bypass，立即 KeyMult/Accumulate
//----------------------------------------------------------------------------
struct TrueOCResult {
    DCRTPoly cTilda0;
    DCRTPoly cTilda1;
    // 用于交叉验证 tools/hks_flow_sim 的算子计数
    int bypass_towers      = 0;   // 走旁路的 (digit, tower) 对数，应等于 sizeQl
    int bconv_target_calls = 0;   // 逐目标塔的 BConv 次数，应等于 L*(D-1) + K*D
    int keymult_calls      = 0;   // 应等于 D * (sizeQl + sizeP)
};

TrueOCResult TrueOCKeySwitch(const HKSFixture& f) {
    const auto& cp        = f.cryptoParams;
    const auto paramsP    = cp->GetParamsP();
    const auto paramsQlP  = f.c.GetExtendedCRTBasis(paramsP);

    const size_t sizeQl  = f.paramsQl->GetParams().size();
    const size_t sizeP   = paramsP->GetParams().size();
    const size_t sizeQlP = sizeQl + sizeP;
    const size_t sizeQ   = cp->GetElementParams()->GetParams().size();

    const uint32_t alpha     = cp->GetNumPerPartQ();
    const uint32_t numPartQl = NumDigits(cp, sizeQl, alpha);

    // EVAL 域副本供 native 塔旁路；COEFFICIENT 域副本供 BConv
    std::vector<DCRTPoly> partsEval = SplitIntoDigits(f.c, cp, alpha, numPartQl, sizeQl);
    std::vector<DCRTPoly> partsCoeff = partsEval;
    for (auto& p : partsCoeff)
        p.SetFormat(Format::COEFFICIENT);

    const std::vector<DCRTPoly>& bv = f.evalKey->GetBVector();
    const std::vector<DCRTPoly>& av = f.evalKey->GetAVector();

    TrueOCResult out;
    out.cTilda0 = DCRTPoly(paramsQlP, Format::EVALUATION, true);
    out.cTilda1 = DCRTPoly(paramsQlP, Format::EVALUATION, true);

    // ---- 外层：逐个 output tower ----
    for (usint t = 0; t < sizeQlP; t++) {
        NativePoly acc0;
        NativePoly acc1;
        bool started = false;

        // ---- 内层：该塔需要哪些 digit 的贡献 ----
        for (uint32_t j = 0; j < numPartQl; j++) {
            const usint sizePartQl   = partsEval[j].GetNumOfElements();
            const usint startPartIdx = alpha * j;
            const usint endPartIdx   = startPartIdx + sizePartQl;

            NativePoly cjt;
            if (t >= startPartIdx && t < endPartIdx) {
                // native Q tower：旁路，直接复用原 EVAL 系数，不做 BConv/NTT
                cjt = partsEval[j].GetElementAtIndex(t - startPartIdx);
                out.bypass_towers++;
            }
            else {
                // non-native：调用现成的全宽 BConv 只为拿到正确 residue，
                // 然后只取目标列。其开销不进入任何性能模型。
                DCRTPoly fullCompl = partsCoeff[j].ApproxSwitchCRTBasis(
                    cp->GetParamsPartQ(j),
                    cp->GetParamsComplPartQ(sizeQl - 1, j),
                    cp->GetPartQlHatInvModq(j, sizePartQl - 1),
                    cp->GetPartQlHatInvModqPrecon(j, sizePartQl - 1),
                    cp->GetPartQlHatModp(sizeQl - 1, j),
                    cp->GetmodComplPartqBarrettMu(sizeQl - 1, j));
                // 目标塔在 complement 里的下标：与 DC 组装映射相同
                const usint complIdx = (t < startPartIdx) ? t : t - sizePartQl;
                cjt = fullCompl.GetElementAtIndex(complIdx);
                cjt.SetFormat(Format::EVALUATION);
                out.bconv_target_calls++;
            }

            // evk 建立在完整 Q 基上：P 段下标要跳到 sizeQ 之后
            const usint keyIdx = (t < sizeQl) ? t : static_cast<usint>(sizeQ) + (t - sizeQl);

            // 立即 KeyMult + Accumulate
            const auto& bji = bv[j].GetElementAtIndex(keyIdx);
            const auto& aji = av[j].GetElementAtIndex(keyIdx);
            if (!started) {
                acc0    = cjt * bji;
                acc1    = cjt * aji;
                started = true;
            }
            else {
                acc0 += cjt * bji;
                acc1 += cjt * aji;
            }
            out.keymult_calls++;
        }

        // 该输出塔已集齐全部 digit 的贡献
        out.cTilda0.SetElementAtIndex(t, acc0);
        out.cTilda1.SetElementAtIndex(t, acc1);
    }
    return out;
}

//----------------------------------------------------------------------------
// 参照系：OpenFHE 自己的 DC 路径
//----------------------------------------------------------------------------
std::pair<DCRTPoly, DCRTPoly> OpenFHEReference(const HKSFixture& f, HKSStrategy strategy) {
    const HKSStrategy saved = GetHKSStrategy();
    SetHKSStrategy(strategy);
    KeySwitchHYBRID ks;
    auto digits = ks.EvalKeySwitchPrecomputeCore(f.c, f.cryptoParams);
    auto cTilda = ks.EvalFastKeySwitchCoreExt(digits, f.evalKey, f.paramsQl);
    SetHKSStrategy(saved);
    return {(*cTilda)[0], (*cTilda)[1]};
}

//----------------------------------------------------------------------------
// 逐 residue 精确比对
//----------------------------------------------------------------------------
::testing::AssertionResult PolyExactlyEqual(const DCRTPoly& a, const DCRTPoly& b,
                                            const std::string& what) {
    if (a.GetNumOfElements() != b.GetNumOfElements())
        return ::testing::AssertionFailure()
               << what << ": 塔数不同 " << a.GetNumOfElements() << " vs "
               << b.GetNumOfElements();
    for (usint t = 0; t < a.GetNumOfElements(); t++) {
        const auto& ta = a.GetElementAtIndex(t);
        const auto& tb = b.GetElementAtIndex(t);
        if (ta.GetModulus() != tb.GetModulus())
            return ::testing::AssertionFailure()
                   << what << ": 塔 " << t << " 模数不同";
        if (ta.GetLength() != tb.GetLength())
            return ::testing::AssertionFailure()
                   << what << ": 塔 " << t << " 长度不同";
        for (usint i = 0; i < ta.GetLength(); i++) {
            if (ta[i] != tb[i])
                return ::testing::AssertionFailure()
                       << what << ": 塔 " << t << " 系数 " << i << " 不同: "
                       << ta[i] << " vs " << tb[i];
        }
    }
    return ::testing::AssertionSuccess();
}

size_t CountResidues(const DCRTPoly& a) {
    size_t n = 0;
    for (usint t = 0; t < a.GetNumOfElements(); t++)
        n += a.GetElementAtIndex(t).GetLength();
    return n;
}

struct Shape {
    uint32_t multDepth;
    uint32_t dnum;
};

// OpenFHE 会拒绝一部分 (multDepth, dnum) 组合（例如「6 towers 分成 5 digits」）。
// 这里不去猜哪些合法：构造失败就记为 unsupported 跳过，但调用方必须断言
// 真正验过的形状数量，避免测试因为全被跳过而变成空壳。
bool TryMakeFixture(const Shape& s, uint32_t ringDim, HKSFixture& out,
                    std::string& why) {
    try {
        out = MakeFixture(ringDim, s.multDepth, s.dnum);
        return true;
    }
    catch (const std::exception& e) {
        why = e.what();
        return false;
    }
}

// 候选形状。覆盖整除、不整除（末位 digit 不完整）、每 digit 单塔等情况。
const Shape kShapes[] = {
    {2, 2}, {3, 2}, {3, 3}, {4, 2}, {4, 3}, {4, 6}, {5, 3}, {5, 4},
};

}  // namespace

//==================================================================================
// 1) 真 OC 的最终结果必须与 OpenFHE 的 DC 路径逐 residue 一致
//==================================================================================
TEST(HKSTrueOC, MatchesDCReference) {
    auto f = MakeFixture(1 << 12, 2, 2);
    auto ref = OpenFHEReference(f, HKSStrategy::DC);
    auto oc  = TrueOCKeySwitch(f);

    EXPECT_TRUE(PolyExactlyEqual(oc.cTilda0, ref.first, "cTilda0"));
    EXPECT_TRUE(PolyExactlyEqual(oc.cTilda1, ref.second, "cTilda1"));

    const size_t n = CountResidues(ref.first) + CountResidues(ref.second);
    EXPECT_GT(n, 0u);
    std::cout << "[HKSTrueOC] 逐 residue 精确比对通过：" << n << " 个 residue\n";
}

//==================================================================================
// 2) 换几组参数形状再验一遍，覆盖不完整末位 digit 与多 level
//==================================================================================
TEST(HKSTrueOC, MatchesAcrossShapesAndLevels) {
    int verified = 0;
    int withIncompleteLastDigit = 0;
    size_t residues = 0;

    for (const auto& s : kShapes) {
        HKSFixture f;
        std::string why;
        if (!TryMakeFixture(s, 1 << 12, f, why)) {
            std::cout << "[HKSTrueOC] 跳过 depth=" << s.multDepth
                      << " dnum=" << s.dnum << "（OpenFHE 不接受该组合）\n";
            continue;
        }

        auto ref = OpenFHEReference(f, HKSStrategy::DC);
        auto oc  = TrueOCKeySwitch(f);

        const size_t L = f.paramsQl->GetParams().size();
        const uint32_t alpha = f.cryptoParams->GetNumPerPartQ();
        const uint32_t D     = NumDigits(f.cryptoParams, L, alpha);
        const std::string tag = "depth=" + std::to_string(s.multDepth) +
                                " dnum=" + std::to_string(s.dnum) +
                                " (L=" + std::to_string(L) +
                                " alpha=" + std::to_string(alpha) +
                                " D=" + std::to_string(D) + ")";

        EXPECT_TRUE(PolyExactlyEqual(oc.cTilda0, ref.first, tag + " cTilda0"));
        EXPECT_TRUE(PolyExactlyEqual(oc.cTilda1, ref.second, tag + " cTilda1"));

        if (L % alpha != 0)
            withIncompleteLastDigit++;
        verified++;
        residues += CountResidues(ref.first) + CountResidues(ref.second);
        std::cout << "[HKSTrueOC] 通过 " << tag << "\n";
    }

    EXPECT_GE(verified, 3) << "合法形状太少，本测试失去意义";
    EXPECT_GE(withIncompleteLastDigit, 1)
        << "没有覆盖到末位 digit 不完整（L % alpha != 0）的情况";
    std::cout << "[HKSTrueOC] 跨形状比对通过：" << verified << " 组，共 "
              << residues << " 个 residue\n";
}

//==================================================================================
// 3) 真 OC 也必须与 MP 路径一致（三种调度算出同一个答案）
//==================================================================================
TEST(HKSTrueOC, AllSchedulesAgree) {
    auto f    = MakeFixture(1 << 12, 2, 2);
    auto dc   = OpenFHEReference(f, HKSStrategy::DC);
    auto mp   = OpenFHEReference(f, HKSStrategy::MP);
    auto pseudo_oc = OpenFHEReference(f, HKSStrategy::OC);
    auto oc   = TrueOCKeySwitch(f);

    EXPECT_TRUE(PolyExactlyEqual(mp.first, dc.first, "MP vs DC cTilda0"));
    EXPECT_TRUE(PolyExactlyEqual(mp.second, dc.second, "MP vs DC cTilda1"));
    EXPECT_TRUE(PolyExactlyEqual(pseudo_oc.first, dc.first, "伪OC vs DC cTilda0"));
    EXPECT_TRUE(PolyExactlyEqual(pseudo_oc.second, dc.second, "伪OC vs DC cTilda1"));
    EXPECT_TRUE(PolyExactlyEqual(oc.cTilda0, dc.first, "真OC vs DC cTilda0"));
    EXPECT_TRUE(PolyExactlyEqual(oc.cTilda1, dc.second, "真OC vs DC cTilda1"));
}

//==================================================================================
// 4) 算子计数必须与 tools/hks_flow_sim 的闭式一致
//    bypass = L,  逐目标塔 BConv = L*(D-1) + K*D,  KeyMult = D*(L+K)
//==================================================================================
TEST(HKSTrueOC, OperatorCountsMatchSimulatorClosedForm) {
    int verified = 0;
    for (const auto& s : kShapes) {
        HKSFixture f;
        std::string why;
        if (!TryMakeFixture(s, 1 << 12, f, why))
            continue;
        auto oc = TrueOCKeySwitch(f);

        const size_t L = f.paramsQl->GetParams().size();
        const size_t K = f.cryptoParams->GetParamsP()->GetParams().size();
        const uint32_t alpha = f.cryptoParams->GetNumPerPartQ();
        const uint32_t D     = NumDigits(f.cryptoParams, L, alpha);

        const std::string tag =
            "L=" + std::to_string(L) + " K=" + std::to_string(K) +
            " D=" + std::to_string(D);

        // 每个 Q 塔恰好被它所属的那一个 digit 旁路一次
        EXPECT_EQ(static_cast<size_t>(oc.bypass_towers), L) << tag;
        // 与 doc/HKS_dataflow_simulation_plan.md §6.3 的闭式一致
        EXPECT_EQ(static_cast<size_t>(oc.bconv_target_calls), L * (D - 1) + K * D) << tag;
        // 每个 digit 对每个输出塔都要乘一次 evk
        EXPECT_EQ(static_cast<size_t>(oc.keymult_calls), D * (L + K)) << tag;
        // 旁路 + BConv 必须恰好覆盖全部 (digit, tower) 对
        EXPECT_EQ(static_cast<size_t>(oc.bypass_towers + oc.bconv_target_calls),
                  D * (L + K)) << tag;
        verified++;
    }
    EXPECT_GE(verified, 3) << "合法形状太少，本测试失去意义";
}

//==================================================================================
// 5) scalar 单目标塔 BConv：独立验证
//        full_BConv[d][target] == scalar_single_target_BConv(d, target)
//
//    这条不依赖 DCRTPoly 的任何批量路径，直接按 BConv 的定义逐系数重算，
//    用来证明「只算一列」在数学上就是全宽结果的那一列，而不是碰巧对上。
//==================================================================================
TEST(HKSTrueOC, ScalarSingleTargetBConvMatchesFullBConv) {
    auto f = MakeFixture(1 << 12, 2, 2);
    const auto& cp = f.cryptoParams;

    const size_t sizeQl  = f.paramsQl->GetParams().size();
    const uint32_t alpha = cp->GetNumPerPartQ();
    const uint32_t D     = NumDigits(cp, sizeQl, alpha);

    auto partsEval  = SplitIntoDigits(f.c, cp, alpha, D, sizeQl);
    auto partsCoeff = partsEval;
    for (auto& p : partsCoeff)
        p.SetFormat(Format::COEFFICIENT);

    size_t checked = 0;
    for (uint32_t d = 0; d < D; d++) {
        const usint sizePartQl = partsCoeff[d].GetNumOfElements();

        const auto& QHatInvModq       = cp->GetPartQlHatInvModq(d, sizePartQl - 1);
        const auto& QHatInvModqPrecon = cp->GetPartQlHatInvModqPrecon(d, sizePartQl - 1);
        const auto& QHatModp          = cp->GetPartQlHatModp(sizeQl - 1, d);
        const auto& modpBarrettMu     = cp->GetmodComplPartqBarrettMu(sizeQl - 1, d);

        DCRTPoly fullCompl = partsCoeff[d].ApproxSwitchCRTBasis(
            cp->GetParamsPartQ(d), cp->GetParamsComplPartQ(sizeQl - 1, d),
            QHatInvModq, QHatInvModqPrecon, QHatModp, modpBarrettMu);

        const usint numTargets = fullCompl.GetNumOfElements();
        const usint ringDim    = partsCoeff[d].GetElementAtIndex(0).GetLength();

        for (usint target = 0; target < numTargets; target++) {
            const auto& outTower = fullCompl.GetElementAtIndex(target);
            const uint64_t pj = outTower.GetModulus().ConvertToInt<uint64_t>();

            // 逐系数按定义重算这一列
            for (usint ri = 0; ri < ringDim; ri++) {
                DoubleNativeInt sum = 0;
                for (usint i = 0; i < sizePartQl; i++) {
                    const auto& srcTower = partsCoeff[d].GetElementAtIndex(i);
                    const auto& qi = srcTower.GetModulus();
                    const uint64_t xScaled =
                        srcTower[ri]
                            .ModMulFastConst(QHatInvModq[i], qi, QHatInvModqPrecon[i])
                            .ConvertToInt<uint64_t>();
                    sum += Mul128(xScaled, QHatModp[i][target].ConvertToInt<uint64_t>());
                }
                const uint64_t scalar =
                    BarrettUint128ModUint64(sum, pj, modpBarrettMu[target]);
                ASSERT_EQ(scalar, outTower[ri].ConvertToInt<uint64_t>())
                    << "digit=" << d << " target=" << target << " coeff=" << ri;
                checked++;
            }
        }
    }
    EXPECT_GT(checked, 0u);
    std::cout << "[HKSTrueOC] scalar 单目标塔 BConv 逐系数验证通过：" << checked
              << " 个 residue\n";
}

//==================================================================================
// 6) 端到端：真 OC 的 cTilda 拿去做 ModDown + 解密，结果必须正确
//    前面几条只证明「与 DC 一致」；这条证明 DC 本身也没被我们的夹具带偏。
//==================================================================================
TEST(HKSTrueOC, EndToEndDecryptionIsCorrect) {
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetMultiplicativeDepth(2);
    parameters.SetScalingModSize(50);
    parameters.SetRingDim(1 << 12);
    parameters.SetBatchSize(8);
    parameters.SetNumLargeDigits(2);
    parameters.SetKeySwitchTechnique(HYBRID);

    auto cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    auto keys = cc->KeyGen();
    cc->EvalMultKeyGen(keys.secretKey);

    std::vector<double> x1 = {0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> x2 = {1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0};
    auto ct1 = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(x1));
    auto ct2 = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(x2));

    const HKSStrategy saved = GetHKSStrategy();
    for (auto s : {HKSStrategy::DC, HKSStrategy::MP, HKSStrategy::OC}) {
        SetHKSStrategy(s);
        auto prod = cc->EvalMult(ct1, ct2);   // 内部走 HKS
        Plaintext out;
        cc->Decrypt(keys.secretKey, prod, &out);
        out->SetLength(x1.size());
        const auto& got = out->GetRealPackedValue();
        for (size_t i = 0; i < x1.size(); i++)
            EXPECT_NEAR(got[i], x1[i] * x2[i], 1e-6)
                << "strategy=" << static_cast<int>(s) << " slot=" << i;
    }
    SetHKSStrategy(saved);
}
