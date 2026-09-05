#ifndef CG_NTT_H
#define CG_NTT_H

//============================================================================
// File   : cg_ntt.h
// Brief  : CG-NTT (Constant Geometry NTT) — 恒定几何 NTT 头文件
//
// 设计要点：
//   - 每一层（stage）的蝶形单元永远读取间距 N/2 的两个数
//   - 每一层的写回按"完美洗牌（Perfect Shuffle）"规则：写到 2i 和 2i+1
//   - 因此地址生成极其简单：读 global_i / global_i+N/2，写 2*global_i / 2*global_i+1
//   - 旋转因子表由 Host 端预计算（模拟数据流动记录每层真实使用的 TF）
//   - FPGA 端只需顺序读取 tf_rom[stage][global_i]，无需任何运行时索引计算
//
// 函数声明：
//   CG_NTT_Kernel      - 单 limb CG-NTT/INTT 核心（1D 布局）
//   Compute_CG_NTT     - 多 limb 包装器（与 Compute_NTT 接口风格一致）
//   flatten_2d_to_1d   - [SQRT][SQRT] → [RING_DIM] 辅助转换
//   reshape_1d_to_2d   - [RING_DIM] → [SQRT][SQRT] 辅助转换
//   cg_ntt_reorder     - CG-NTT 输出重排（bit-reversal + unshuffle）
//============================================================================

#include "define.h"
#include "arithmetic.h"
#include "ntt_kernel.h"

// =========================================================
// 辅助常量（不修改 define.h，在本文件局部定义）
// =========================================================
static const int CG_HALF_N  = RING_DIM / 2;   // 2048：蝶形跨度
static const int CG_PE_NUM  = PE_PARALLEL;    // 并行 PE 数（复用 define.h 的 PE_PARALLEL）

// 乒乓缓冲 buf_A / buf_B 的 bank 数 = 2 × PE 数。
// 定这个倍数的是**写端口**，不是读端口：
//   写：NTT 每个 PE 写 [2i, 2i+1] 两个连续地址，CG_PE_NUM 个 PE 一拍共写
//       2*CG_PE_NUM 个连续地址；切成同样多的 bank 才能做到 1 写/bank。
//   读：u 与 v 相距 CG_HALF_N，而 CG_HALF_N % CG_BUF_PARTITION == 0，
//       两组读必然落在同一批 bank 上（2 访问/bank），靠 ram_2p 的双端口吸收，
//       所以读侧对 bank 数没有额外要求。
static const int CG_BUF_PARTITION = 2 * PE_PARALLEL;

// 512-bit 总线打包参数
static const int PACK_RATIO     = 512 / 64;                        // 8
static const int PACKED_RING_DIM = RING_DIM / PACK_RATIO;          // 512
static const int PACKED_TW_SIZE  = (STAGE * CG_HALF_N) / PACK_RATIO; // 3072

// =========================================================
// 核心：单 limb CG-NTT / INTT
// =========================================================
// IS_NTT = true  → 正向 NTT（Cooley-Tukey）
// IS_NTT = false → 逆向 INTT（Gentleman-Sande）
// Compatibility API for standalone tests; Top uses the runtime engine below.
template <bool IS_NTT>
void CG_NTT_Kernel(
    const uint64_t in_data[RING_DIM],
    uint64_t out_data[RING_DIM],
    const uint64_t modulus,
    const uint64_t S,
    const uint64_t M_barrett,
    const uint64_t cg_twiddle[STAGE][CG_HALF_N]
);

// One physical engine: direction is a runtime control, not a template argument.
// Both twiddle memories remain resident; only the selected direction is read.
void CG_Transform_Kernel(
    const uint64_t in_data[RING_DIM], uint64_t out_data[RING_DIM],
    uint64_t modulus, uint64_t S, uint64_t M_barrett,
    const uint64_t ntt_twiddle[STAGE][CG_HALF_N],
    const uint64_t intt_twiddle[STAGE][CG_HALF_N], bool is_ntt
);

// Transform already-loaded ping-pong banks; no intermediate flat arrays.
// The caller owns both banks and loads/stores at the existing PE width.
void CG_Transform_Banks(
    uint64_t buf_A[RING_DIM], uint64_t buf_B[RING_DIM],
    uint64_t modulus, uint64_t S, uint64_t M_barrett,
    const uint64_t ntt_twiddle[STAGE][CG_HALF_N],
    const uint64_t intt_twiddle[STAGE][CG_HALF_N], bool is_ntt
);

// Top's shared work memory is the A bank. tower is a physical slot, independent
// of the global modulus/twiddle index selected by the caller. Even STAGE leaves
// the final result in work[tower]; scratch is reused by every direction/tower.
void CG_Transform_Work(
    uint64_t work[MAX_LIMBS][SQRT][SQRT], int tower,
    uint64_t scratch[RING_DIM],
    uint64_t modulus, uint64_t S, uint64_t M_barrett,
    const uint64_t ntt_twiddle[STAGE][CG_HALF_N],
    const uint64_t intt_twiddle[STAGE][CG_HALF_N], bool is_ntt,
    bool scale_only, uint64_t scale_factor
);

// =========================================================
// 多 limb 包装器
// =========================================================
// 与现有 Compute_NTT 接口风格一致，方便在 top.cpp 中集成
extern "C" {
    void Compute_CG_NTT(
        ap_uint<512> *in_data,
        const ap_uint<512> *cg_ntt_twiddle,
        const ap_uint<512> *cg_intt_twiddle,
        const uint64_t modulus[MAX_LIMBS],
        const uint64_t S[MAX_LIMBS],
        const uint64_t M_barrett[MAX_LIMBS],
        bool is_ntt,
        int num_active_limbs,
        int mod_idx_offset
    );
}

// =========================================================
// 2D ↔ 1D 布局转换辅助（用于与现有 [SQRT][SQRT] 接口互转）
// =========================================================
extern "C" {
    // [SQRT][SQRT] → [RING_DIM]（行优先平铺）
    void flatten_2d_to_1d(
        const uint64_t src[SQRT][SQRT],
        uint64_t dst[RING_DIM]
    );
}

extern "C" {
    // [RING_DIM] → [SQRT][SQRT]（行优先还原）
    void reshape_1d_to_2d(
        const uint64_t src[RING_DIM],
        uint64_t dst[SQRT][SQRT]
    );
}

// =========================================================
// CG-NTT 输出重排
// =========================================================
// CG-NTT 完成后，数据处于经过 STAGE 次 perfect shuffle 后的乱序状态。
// 本函数通过预计算的逆排列表，将数据还原为标准 NTT 输出顺序。
// 在 FPGA 上，此操作可在 DDR 写回时顺手完成（零额外延迟）。
extern "C" {
    void cg_ntt_reorder(uint64_t data[RING_DIM]);
}

#endif // CG_NTT_H
