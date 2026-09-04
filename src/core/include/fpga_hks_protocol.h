#ifndef OPENFHE_FPGA_HKS_PROTOCOL_H
#define OPENFHE_FPGA_HKS_PROTOCOL_H
#include <cstdint>

// Shared host/HLS wire constants. This header intentionally has no ap_int or
// OpenFHE types and no arithmetic precomputation. All sizes count uint64 words.
namespace fpga_hks_abi {
constexpr unsigned kVersion = 2;
constexpr unsigned kLogN = 12, kN = 1u << kLogN;
constexpr unsigned kQ = 3, kP = 2, kTotal = kQ + kP;
constexpr uint8_t kQueryOpcode = 9;  // 7 remains retired; 0..8 retain their IDs.
constexpr uint64_t kMagic = 0xC84B535F41424932ULL; // outside any supported residue range
constexpr uint64_t kCapabilities = (uint64_t(0x484B5332) << 32) |
    (uint64_t(kLogN) << 24) | (uint64_t(kQ) << 16) | (uint64_t(kP) << 8) | kTotal;
constexpr unsigned kHeaderWords = 2;
constexpr unsigned kWeightOffset = kHeaderWords, kWeightWords = kQ * kTotal;
constexpr unsigned kPreconOffset = kWeightOffset + kWeightWords;
constexpr unsigned kTailOffset = kPreconOffset + kWeightWords;
constexpr unsigned kHksInvOffset = kTailOffset, kHksWords = kHksInvOffset + kQ;
constexpr unsigned kBconvModOffset = kTailOffset, kBconvWords = kBconvModOffset + kTotal;
constexpr uint64_t kHksDescriptor = (uint64_t(kVersion) << 32) | kHksWords;
constexpr uint64_t kBconvDescriptor = (uint64_t(kVersion) << 32) | kBconvWords;
} // namespace fpga_hks_abi
#endif
