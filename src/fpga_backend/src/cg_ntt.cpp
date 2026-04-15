#include "../include/cg_ntt.h"

//============================================================================
// File   : cg_ntt.cpp
// Brief  : CG-NTT（恒定几何 NTT）核心实现
//
// 算法概述（CG-NTT 的三大硬件红利）：
//   ① 每层 PE 永远读取相距 N/2 的两个数：Read(i) 和 Read(i + N/2)
//   ② 蝶形运算完成后，按"完美洗牌"规则写回：Write(2i) 和 Write(2i+1)
//   ③ 由于上层已经洗过牌，下一层 PE 依然只需读固定位置的数
//   → MUX/交叉开关彻底消灭，变为固定物理硬连线
//
// 存储架构：
//   - buf_A[RING_DIM] / buf_B[RING_DIM]：1D 乒乓缓冲
//   - cyclic factor=16：保证 8 PE × (读u + 读v) 无 Bank Conflict
//   - ram_2p：读端 global_i 和 global_i+N/2 落同一 bank，需双端口
//
// 旋转因子：
//   - cg_twiddle[STAGE][CG_HALF_N]：由 Host 端预计算
//   - FPGA 顺序消费：tf = cg_twiddle[stage][global_i]，无运行时索引计算
//============================================================================

// ============================================================
// 2D ↔ 1D 布局转换辅助
// ============================================================

void flatten_2d_to_1d(
    const uint64_t src[SQRT][SQRT],
    uint64_t dst[RING_DIM]
) {
    #pragma HLS INLINE
    FLATTEN_ROW:
    for (int i = 0; i < SQRT; i++) {
        #pragma HLS PIPELINE II=1
        FLATTEN_COL:
        for (int j = 0; j < SQRT; j++) {
            #pragma HLS UNROLL factor=PE_PARALLEL
            dst[i * SQRT + j] = src[i][j];
        }
    }
}

void reshape_1d_to_2d(
    const uint64_t src[RING_DIM],
    uint64_t dst[SQRT][SQRT]
) {
    #pragma HLS INLINE
    RESHAPE_ROW:
    for (int i = 0; i < SQRT; i++) {
        #pragma HLS PIPELINE II=1
        RESHAPE_COL:
        for (int j = 0; j < SQRT; j++) {
            #pragma HLS UNROLL factor=PE_PARALLEL
            dst[i][j] = src[i * SQRT + j];
        }
    }
}

// ============================================================
// CG-NTT 输出重排：还原为标准 NTT 输出顺序
//
// CG-NTT 经过 STAGE 次 perfect shuffle 后，数据位于排列 perm[] 处。
// 逆映射：result[perm[i]] = data[i]  → 即 result[物理位置] = data[乱序位置]
//
// perfect shuffle 定义：
//   new[2*i]     = old[i]           (前半段 → 偶数位)
//   new[2*i + 1] = old[i + N/2]     (后半段 → 奇数位)
// ============================================================

void cg_ntt_reorder(uint64_t data[RING_DIM]) {
    #pragma HLS INLINE off

    // CG-NTT 输出为 bit-reversed 顺序（相对于标准 NTT 输出）。
    // 只需执行 bit-reversal 排列即可还原为标准顺序。
    uint64_t temp[RING_DIM];

    REORDER_LOOP:
    for (int i = 0; i < RING_DIM; i++) {
        #pragma HLS PIPELINE II=1
        // bit-reverse index i (STAGE bits = log2(RING_DIM) bits)
        int rev = 0;
        int x = i;
        for (int b = 0; b < STAGE; b++) {
            rev = (rev << 1) | (x & 1);
            x >>= 1;
        }
        temp[rev] = data[i];
    }

    // 写回
    WRITEBACK_REORDER:
    for (int i = 0; i < RING_DIM; i++) {
        #pragma HLS PIPELINE II=1
        data[i] = temp[i];
    }
}

// ============================================================
// CG_NTT_Kernel：单 limb CG-NTT / INTT 核心
// ============================================================

void CG_NTT_Kernel(
    uint64_t in_data[RING_DIM],
    const uint64_t modulus,
    const uint64_t K_HALF,
    const uint64_t M_barrett,
    const uint64_t cg_twiddle[STAGE][CG_HALF_N],
    bool is_ntt
) {
    // ============================================================
    // 乒乓缓冲：消除 RAW 依赖
    // NTT 方向：读 [i, i+N/2]，写 [2i, 2i+1]（perfect shuffle）
    // INTT 方向：读 [2i, 2i+1]，写 [i, i+N/2]（perfect unshuffle）
    // ============================================================
    uint64_t buf_A[RING_DIM];
    uint64_t buf_B[RING_DIM];

    #pragma HLS ARRAY_PARTITION variable=buf_A cyclic factor=16 dim=1
    #pragma HLS ARRAY_PARTITION variable=buf_B cyclic factor=16 dim=1
    #pragma HLS BIND_STORAGE variable=buf_A type=ram_2p impl=bram
    #pragma HLS BIND_STORAGE variable=buf_B type=ram_2p impl=bram

    // 旋转因子：stage 维 complete（8 副本），位置维 cyclic 按 PE_PARALLEL
    // 保证 8 个 PE 同时读 cg_twiddle[s][i*8+0..7] 无冲突
    #pragma HLS ARRAY_PARTITION variable=cg_twiddle cyclic factor=CG_PE_NUM dim=2

    // ============================================================
    // 初始化：in_data → buf_A
    // ============================================================
    INIT_LOOP:
    for (int i = 0; i < RING_DIM / CG_PE_NUM; i++) {
        #pragma HLS PIPELINE II=1
        CG_INIT_PE:
        for (int p = 0; p < CG_PE_NUM; p++) {
            #pragma HLS UNROLL
            buf_A[i * CG_PE_NUM + p] = in_data[i * CG_PE_NUM + p];
        }
    }

    // ============================================================
    // 主循环：STAGE 层 × (CG_HALF_N / CG_PE_NUM) 次迭代
    // 乒乓协议：偶数 stage 读 A 写 B，奇数 stage 读 B 写 A
    //
    // NTT（正向）：stage 0→11，每层 shuffle 写
    //   读 buf[i] 和 buf[i+N/2]，写 buf[2i] 和 buf[2i+1]
    //
    // INTT（逆向）：stage 11→0，每层 unshuffle 写
    //   读 buf[2i] 和 buf[2i+1]，写 buf[i] 和 buf[i+N/2]
    // ============================================================
    STAGE_LOOP:
    for (int stage = 0; stage < STAGE; stage++) {

        // NTT 正序（0,1,...,11），INTT 逆序（11,10,...,0）
        int actual_stage = is_ntt ? stage : (STAGE - 1 - stage);

        BUTTERFLY_LOOP:
        for (int i = 0; i < CG_HALF_N / CG_PE_NUM; i++) {
            #pragma HLS PIPELINE II=1
            #pragma HLS DEPENDENCE variable=buf_A inter false
            #pragma HLS DEPENDENCE variable=buf_B inter false

            PE_UNROLL:
            for (int p = 0; p < CG_PE_NUM; p++) {
                #pragma HLS UNROLL

                int global_i = i * CG_PE_NUM + p;  // 0 ~ CG_HALF_N-1

                uint64_t u, v;
                if (is_ntt) {
                    // NTT 读：固定跨度 [global_i, global_i + N/2]
                    if ((stage & 1) == 0) {
                        u = buf_A[global_i];
                        v = buf_A[global_i + CG_HALF_N];
                    } else {
                        u = buf_B[global_i];
                        v = buf_B[global_i + CG_HALF_N];
                    }
                } else {
                    // INTT 读：连续对 [2*global_i, 2*global_i + 1]
                    if ((stage & 1) == 0) {
                        u = buf_A[2 * global_i];
                        v = buf_A[2 * global_i + 1];
                    } else {
                        u = buf_B[2 * global_i];
                        v = buf_B[2 * global_i + 1];
                    }
                }

                // ② 顺序读取旋转因子（无运行时索引计算！）
                uint64_t tf = cg_twiddle[actual_stage][global_i];

                // ③ 蝶形运算（复用现有 Configurable_PE）
                uint64_t out_u, out_v;
                Configurable_PE(u, v, tf, out_u, out_v, modulus, K_HALF, M_barrett, is_ntt);

                if (is_ntt) {
                    // NTT 写：完美洗牌 [2*global_i, 2*global_i + 1]
                    if ((stage & 1) == 0) {
                        buf_B[2 * global_i]     = out_u;
                        buf_B[2 * global_i + 1] = out_v;
                    } else {
                        buf_A[2 * global_i]     = out_u;
                        buf_A[2 * global_i + 1] = out_v;
                    }
                } else {
                    // INTT 写：完美逆洗牌 [global_i, global_i + N/2]
                    if ((stage & 1) == 0) {
                        buf_B[global_i]              = out_u;
                        buf_B[global_i + CG_HALF_N]  = out_v;
                    } else {
                        buf_A[global_i]              = out_u;
                        buf_A[global_i + CG_HALF_N]  = out_v;
                    }
                }
            }
        }
    }

    // ============================================================
    // 回写到 in_data
    // STAGE=12（偶数）→ 最后写的 stage=11（奇）→ 写入 buf_A → 结果在 buf_A
    // 通用：偶数 STAGE 结果在 buf_A，奇数 STAGE 结果在 buf_B
    // ============================================================
    WRITEBACK_LOOP:
    for (int i = 0; i < RING_DIM / CG_PE_NUM; i++) {
        #pragma HLS PIPELINE II=1
        CG_WB_PE:
        for (int p = 0; p < CG_PE_NUM; p++) {
            #pragma HLS UNROLL
            if ((STAGE & 1) == 0) {
                in_data[i * CG_PE_NUM + p] = buf_A[i * CG_PE_NUM + p];
            } else {
                in_data[i * CG_PE_NUM + p] = buf_B[i * CG_PE_NUM + p];
            }
        }
    }
}

// ============================================================
// Compute_CG_NTT：多 limb 包装器
// 与现有 Compute_NTT 接口风格一致
// ============================================================

void Compute_CG_NTT(
    uint64_t in_data[MAX_LIMBS][RING_DIM],
    const uint64_t cg_ntt_twiddle[MAX_LIMBS][STAGE][CG_HALF_N],
    const uint64_t cg_intt_twiddle[MAX_LIMBS][STAGE][CG_HALF_N],
    const uint64_t modulus[MAX_LIMBS],
    const uint64_t K_HALF[MAX_LIMBS],
    const uint64_t M_barrett[MAX_LIMBS],
    bool is_ntt,
    int num_active_limbs,
    int mod_idx_offset
) {
    #pragma HLS ARRAY_PARTITION variable=in_data cyclic factor=CG_PE_NUM dim=2

    // 旋转因子表：dim=1 (limb) 不拆，dim=2 (stage) complete，dim=3 cyclic
    #pragma HLS ARRAY_PARTITION variable=cg_ntt_twiddle  complete dim=2
    #pragma HLS ARRAY_PARTITION variable=cg_intt_twiddle complete dim=2

    LIMB_LOOP:
    for (int l = mod_idx_offset; l < mod_idx_offset + num_active_limbs; l++) {
        #pragma HLS LOOP_TRIPCOUNT min=1 max=8 avg=3
        if (is_ntt) {
            CG_NTT_Kernel(
                in_data[l],
                modulus[l],
                K_HALF[l],
                M_barrett[l],
                cg_ntt_twiddle[l],
                true
            );
        } else {
            CG_NTT_Kernel(
                in_data[l],
                modulus[l],
                K_HALF[l],
                M_barrett[l],
                cg_intt_twiddle[l],
                false
            );
        }
    }
}
