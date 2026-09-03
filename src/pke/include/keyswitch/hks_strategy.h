#ifndef LBCRYPTO_CRYPTO_KEYSWITCH_HKS_STRATEGY_H
#define LBCRYPTO_CRYPTO_KEYSWITCH_HKS_STRATEGY_H

#include <cstdint>
#include <ostream>
#include <vector>
#include <chrono>

namespace lbcrypto {

enum class HKSStrategy {
    DC,  // Digit-Centric: per-digit INTT→BConv→NTT (default, matches current code)
    MP,  // Max-Parallel:  all-INTT → all-BConv → all-NTT (global barriers between phases)
    OC,  // Output-Centric: per-output-tower BConv with sizeP=1 (minimal peak SRAM)
};

inline HKSStrategy& GetHKSStrategy() {
    static HKSStrategy s = HKSStrategy::DC;
    return s;
}

inline void SetHKSStrategy(HKSStrategy s) {
    GetHKSStrategy() = s;
}

// ---------------------------------------------------------------------------
// Per-KeySwitch operation statistics
// Reset via ResetHKSStats() before a call; read via GetHKSStats() after.
// ---------------------------------------------------------------------------
struct HKSStats {
    // --- Operation counts (EvalKeySwitchPrecomputeCore) ---
    int intt_poly   = 0;  // DCRTPoly-level INTT calls  (each covers alpha limbs)
    int ntt_poly    = 0;  // DCRTPoly-level NTT  calls
    int ntt_limb    = 0;  // single-limb NTT calls (OC: one per P-tower per digit)
    int bconv       = 0;  // ApproxSwitchCRTBasis (BConv) calls
    int fused_digits = 0;  // OP_HKS_DIGIT calls; not separate CPU INTT/BConv/NTT
    int64_t time_fused_us = 0;  // Host wall time incl. packing and context init if needed

    // --- Operation counts (EvalFastKeySwitchCoreExt) ---
    int modmul_limb = 0;  // limb-level multiply-accumulate iterations

    // --- Parameters captured at precompute time ---
    int num_digits  = 0;  // numPartQl
    int size_ql     = 0;  // current ciphertext Q limbs
    int size_p      = 0;  // auxiliary P limbs
    int alpha       = 0;  // limbs per digit
    int ring_dim    = 0;  // N

    // --- Peak SRAM: max P-tower complement ring-elements held simultaneously ---
    // One ring-element = ring_dim * 8 bytes.
    // DC:  sizeP            (one full complement per digit, one digit at a time)
    // MP:  numPartQl*sizeP  (all complements held simultaneously)
    // OC:  1                (one P-tower at a time)
    int peak_p_towers = 0;

    // --- Sub-operation timing (microseconds, accumulated over all digits/calls) ---
    // Measured on the CPU software path; FPGA path timing is in TransferStats.
    int64_t time_intt_us   = 0;
    int64_t time_bconv_us  = 0;
    int64_t time_ntt_us    = 0;
    int64_t time_modmul_us = 0;
};

inline HKSStats& GetHKSStats() {
    static HKSStats s;
    return s;
}

inline void ResetHKSStats() {
    GetHKSStats() = HKSStats{};
}

// ---------------------------------------------------------------------------
// Real-time intermediate-buffer memory tracker (software probe)
// ---------------------------------------------------------------------------
struct MemoryEvent {
    int step;
    int64_t delta;
    int64_t watermark;
    const char* op;
    int digit;
    int tower;
};

class MemoryTracker {
    int64_t current_ = 0;
    int64_t peak_    = 0;
    int     step_    = 0;
    std::vector<MemoryEvent> log_;

public:
    void alloc(int64_t bytes, const char* op, int digit, int tower = -1) {
        current_ += bytes;
        if (current_ > peak_) peak_ = current_;
        log_.push_back({step_++, bytes, current_, op, digit, tower});
    }

    void free(int64_t bytes, const char* op, int digit, int tower = -1) {
        current_ -= bytes;
        if (current_ < 0) current_ = 0;
        log_.push_back({step_++, -bytes, current_, op, digit, tower});
    }

    int64_t peak()    const { return peak_; }
    int64_t current() const { return current_; }
    const std::vector<MemoryEvent>& events() const { return log_; }

    void reset() {
        current_ = 0;
        peak_    = 0;
        step_    = 0;
        log_.clear();
    }

    void dump_csv(std::ostream& os) const {
        os << "step,op,digit,tower,delta_bytes,watermark_bytes\n";
        for (auto& e : log_) {
            os << e.step << "," << e.op << "," << e.digit << ","
               << e.tower << "," << e.delta << "," << e.watermark << "\n";
        }
    }
};

inline MemoryTracker& GetMemoryTracker() {
    static MemoryTracker t;
    return t;
}

inline void ResetMemoryTracker() {
    GetMemoryTracker().reset();
}

// ---------------------------------------------------------------------------
// PCIe transfer / kernel timing (补4-3: T_transfer)
// Populated by FpgaManager::Execute() and BConvOffload() when FPGA is active.
// ---------------------------------------------------------------------------
struct TransferStats {
    int64_t h2d_us    = 0;  // Host→Device sync time (accumulated)
    int64_t kernel_us = 0;  // Pure kernel execution time (accumulated)
    int64_t d2h_us    = 0;  // Device→Host sync time (accumulated)
    int     calls     = 0;  // Number of Execute/BConvOffload calls
};

inline TransferStats& GetTransferStats() {
    static TransferStats s;
    return s;
}

inline void ResetTransferStats() {
    GetTransferStats() = TransferStats{};
}

}  // namespace lbcrypto


#endif
