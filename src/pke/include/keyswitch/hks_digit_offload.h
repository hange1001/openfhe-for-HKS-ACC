#ifndef OPENFHE_HKS_DIGIT_OFFLOAD_H
#define OPENFHE_HKS_DIGIT_OFFLOAD_H

#include "lattice/lat-hal.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace lbcrypto {
template <class Element> class CryptoParametersBase;

enum class HKSDigitBackend { Off, CModel, XRT };
using HKSTopFunction = void (*)(const uint64_t*, const uint64_t*, uint64_t*, uint8_t, int, int);

struct HKSDigitStats {
    uint64_t initializations = 0;
    uint64_t digits = 0;
    uint64_t fallbacks = 0;
    uint64_t failures = 0;
    // Logical opcode-8 payload only, not measured PCIe traffic or total traffic.
    uint64_t input_bytes = 0;
    uint64_t metadata_bytes = 0;
    uint64_t output_bytes = 0;
    std::string last_fallback;
};

// Explicit opt-in; default is Off. After Configure(), Off means CPU and legacy
// FPGA hooks remain disabled. Applications that never Configure() are unchanged.
// CModel must be the actual HLS Top function.
// Configure only while OpenFHE is idle. The fused device context is serialized;
// do not issue external Top / legacy InitModuli calls while this backend owns it.
void ConfigureHKSDigitBackend(HKSDigitBackend backend, HKSTopFunction cmodel = nullptr);
HKSDigitStats GetHKSDigitStats();
void ResetHKSDigitStats();

// Returns false before any digit launch for unsupported configurations. Output is
// committed only after all digits succeed; transport/corrupt-result errors throw.
bool TryHKSDigitOffload(const DCRTPoly& c,
                       const std::shared_ptr<CryptoParametersBase<DCRTPoly>>& cryptoParams,
                       std::vector<DCRTPoly>& output);
}  // namespace lbcrypto
#endif
