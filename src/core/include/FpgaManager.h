#ifndef _FPGA_MANAGER_H_
#define _FPGA_MANAGER_H_

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstring>   // memset, memcpy
#include <cstdlib>   // std::getenv
#include <cmath>     // std::log2, std::ceil
#include <algorithm> // std::max, std::find
#include <numeric>   // std::gcd (C++17)
#include <chrono>    // high_resolution_clock

// =============================================================
// 1. XRT Configuration (不变)
// =============================================================
#ifdef OPENFHE_FPGA_ENABLE

    #ifndef XRT_VERSION_CODE
    #define XRT_VERSION_CODE ((2 << 16) | (0 << 8) | 0)
    #endif

    #ifndef XRT_MAJOR
    #define XRT_MAJOR(v) (((v) >> 16) & 0xffff)
    #endif

    #ifndef XRT_MINOR
    #define XRT_MINOR(v) (((v) >> 8) & 0xff)
    #endif

    #include <xrt/xrt_device.h>
    #include <xrt/xrt_kernel.h>
    #include <xrt/xrt_bo.h>

#endif

// =============================================================
// 2. Definitions 
// =============================================================
#define OP_INIT   0
#define OP_ADD    1
#define OP_SUB    2
#define OP_MULT   3
#define OP_NTT    4
#define OP_INTT   5
#define OP_BCONV  6
#define OP_AUTO   7
#define OP_HKS_DIGIT 8

// ---- 必须与 FPGA 端 define.h 保持一致 ----
// MAX_LIMBS = LIMB_Q (3) + MAX_OUT_COLS (5) = 8
// OP_INIT 时 FPGA 会按 MAX_LIMBS 个 limb 读取 CG-NTT 旋转因子表，
// 若 Host 端分配不足会发生 AXI 越界。
#define MAX_LIMBS 8
#define FPGA_LIMB_Q 3
#define FPGA_LIMB_P 2
#define FPGA_RING_DIM  4096
#define STAGE_NUM 12
#define CG_HALF_N  (FPGA_RING_DIM / 2)          // 2048
#define CG_TF_SIZE (STAGE_NUM * CG_HALF_N)      // 12 * 2048 = 24576

// =============================================================
// 2.5 PCIe Transfer Statistics (补4-3: T_transfer)
// Accumulated by Execute() and BConvOffload() when FPGA is active.
// =============================================================
struct FpgaOpStats {
    int64_t h2d_us    = 0;
    int64_t kernel_us = 0;
    int64_t d2h_us    = 0;
    // 实验 D（task.yaml step 1.2）：整个 Execute()/BConvOffload() 的墙钟。
    // wall - (h2d + kernel + d2h) = 之前完全没被计入的那一段，即
    // xrt::bo 构造/析构 x3 + bo.write() x2 + bo_out.read()。
    int64_t wall_us   = 0;
    int     calls     = 0;
};

struct FpgaTransferStats {
    FpgaOpStats by_opcode[9];  // indexed by OP_INIT..OP_HKS_DIGIT
    // Aggregate totals (sum of all opcodes)
    int64_t h2d_us    = 0;
    int64_t kernel_us = 0;
    int64_t d2h_us    = 0;
    int64_t wall_us   = 0;
    int     calls     = 0;
};

inline FpgaTransferStats& GetFpgaTransferStats() {
    static FpgaTransferStats s;
    return s;
}

inline void ResetFpgaTransferStats() {
    GetFpgaTransferStats() = FpgaTransferStats{};
}

// 实验 D 的被测量：未计时开销。预测值 30-80 us/call（task.yaml:216）。
// 其中 memcpy 分量已由 testbench/host_overhead_bench.cpp 离线测得
// （NTT 路径 1.8 us、BConv 路径 14.0 us @ i7-12700H），
// 差额即 xrt::bo 分配开销 = 「buffer 池化」能消除的收益上界。
inline int64_t FpgaUntimedUs(const FpgaOpStats& s) {
    int64_t accounted = s.h2d_us + s.kernel_us + s.d2h_us;
    return (s.wall_us > accounted) ? (s.wall_us - accounted) : 0;
}

inline void PrintFpgaTransferStats(std::ostream& os = std::cout) {
    static const char* kOpName[9] = {"INIT", "ADD", "SUB", "MULT", "NTT", "INTT", "BCONV", "AUTO", "HKS_DIGIT"};
    const auto& ts = GetFpgaTransferStats();
    os << "\n===== FPGA transfer stats (实验 D 口径) =====\n"
       << "  未计时开销 = wall - (h2d + kernel + d2h)"
          "  <- xrt::bo 分配/析构 + bo.write/read\n"
       << "  op      calls      h2d     kernel        d2h       wall    未计时   未计时占比\n";
    for (int op = 0; op < 9; ++op) {
        const auto& s = ts.by_opcode[op];
        if (s.calls == 0)
            continue;
        int64_t untimed = FpgaUntimedUs(s);
        double  pct     = (s.wall_us > 0) ? (100.0 * untimed / s.wall_us) : 0.0;
        os << "  " << kOpName[op] << "\t" << s.calls << "\t" << s.h2d_us << "us\t" << s.kernel_us
           << "us\t" << s.d2h_us << "us\t" << s.wall_us << "us\t" << untimed << "us\t" << pct
           << "%\n";
    }
    int64_t total_untimed = (ts.wall_us > ts.h2d_us + ts.kernel_us + ts.d2h_us)
                                ? ts.wall_us - (ts.h2d_us + ts.kernel_us + ts.d2h_us)
                                : 0;
    os << "  ----\n  TOTAL\t" << ts.calls << "\t" << ts.h2d_us << "us\t" << ts.kernel_us << "us\t"
       << ts.d2h_us << "us\t" << ts.wall_us << "us\t" << total_untimed << "us\n";
    if (ts.calls > 0) {
        os << "  每次 call 平均: wall " << (ts.wall_us / ts.calls) << "us，其中未计时 "
           << (total_untimed / ts.calls) << "us\n";
    }
    os << "=============================================\n";
}

inline std::string GetXclbinPath() {
    const char* env_path = std::getenv("XCLBIN_PATH");
    if (env_path) return std::string(env_path);

    const char* mode = std::getenv("XCL_EMULATION_MODE");
    std::string base = "/home/timhan/FHE/openfhe-for-HKS-ACC/src/fpga_backend/build/";
    std::string target = (mode && std::string(mode) == "hw_emu") ? "hw_emu"
                       : (mode && std::string(mode) == "hw")     ? "hw"
                       :                                           "sw_emu";
    return base + target + "/fhe_kernels_" + target + ".xclbin";
}

// =============================================================
// 3. Math Helpers 
// =============================================================
class MathUtils {
public:
    static uint64_t Power(unsigned __int128 base, unsigned __int128 exp, uint64_t mod) {
        unsigned __int128 res = 1;
        base %= mod;
        while (exp > 0) {
            if (exp % 2 == 1) res = (res * base) % mod;
            base = (base * base) % mod;
            exp /= 2;
        }
        return (uint64_t)res;
    }

    static uint64_t ModInverse(uint64_t n, uint64_t mod) {
        return Power(n, mod - 2, mod);
    }

    static uint64_t GCD(uint64_t a, uint64_t b) {
        while (b) {
            a %= b;
            std::swap(a, b);
        }
        return a;
    }

    static std::vector<uint64_t> GetPrimeFactors(uint64_t n) {
        std::vector<uint64_t> factors;
        uint64_t temp = n;
        for (uint64_t i = 2; i * i <= temp; ++i) {
            if (temp % i == 0) {
                factors.push_back(i);
                while (temp % i == 0) {
                    temp /= i;
                }
            }
        }
        if (temp > 1) {
            factors.push_back(temp);
        }
        return factors;
    }

    static bool IsPrimitiveRoot(uint64_t a, uint64_t p) {
        if (GCD(a, p) != 1) return false;
        uint64_t phi = p - 1;  
        std::vector<uint64_t> factors = GetPrimeFactors(phi);
        if (factors.empty()) return false; 

        for (uint64_t q_i : factors) {
            // 使用 128 位进行中间计算
            if (Power(a, phi / q_i, p) == 1) {
                return false;
            }
        }

        return true;
    }


    static uint64_t FindSmallestPrimitiveRoot(uint64_t p) {
        if (p <= 4) return p - 1;

     
        for (uint64_t g = 2; g < p; ++g) {
            if (IsPrimitiveRoot(g, p)) {
                return g; 
            }
        }
        return 0; 
    }


    // [修改] 使用 FindSmallestPrimitiveRoot 替换朴素搜索
    static uint64_t Find2NthRootOfUnity(uint64_t modulus, uint64_t n) {
        uint64_t required_order = 2 * n; 
        if ((modulus - 1) % required_order != 0) {
            throw std::runtime_error("Modulus does not support 2N-th root of unity");
        }
        
        // 1. 找到最小原根 g (SymPy 的 'pr')
        // WARNING: 这个调用可能非常慢，因为它包含了大数的质因数分解。
        uint64_t primitive_root_g = FindSmallestPrimitiveRoot(modulus);
        
        if (primitive_root_g == 0) {
            throw std::runtime_error("Failed to find smallest primitive root for the modulus.");
        }
        
        unsigned __int128 exponent = (modulus - 1) / required_order;
        uint64_t psi = Power(primitive_root_g, exponent, modulus);

        // 3. 验证性质: psi^N == -1 (即 p-1)
        if (Power(psi, n, modulus) == (modulus - 1)) {
            return psi;
        } else {
             throw std::runtime_error("Calculated Psi failed the N-th power test. Primitive root found but exponentiation error.");
        }
    }


    // ---------------------------------------------------------
    // CG-NTT 旋转因子预计算（Host 端）
    //
    // 输出布局：out[s * CG_HALF_N + t]，s ∈ [0, STAGE_NUM)，t ∈ [0, CG_HALF_N)
    // 算法：模拟 perfect-shuffle 路由网络，记录每层实际使用的 ψ 幂次。
    //   - 初始排列 perm[i] = bit_reverse(i, STAGE_NUM)
    //   - 每经过一层，新排列：np[2i]=perm[i]，np[2i+1]=perm[i+N/2]
    //   - 第 s 层、位置 i 的指数 = (perm[i] mod 2^s) * N / 2^s
    //
    // 与 FPGA 端 cg_ntt.cpp 中 cg_twiddle[s][t] 的顺序完全匹配。
    // root_2N 是 2N-th primitive root of unity（ψ）。
    // ---------------------------------------------------------
    static int BitReverse(int x, int bits) {
        int r = 0;
        for (int b = 0; b < bits; b++) { r = (r << 1) | (x & 1); x >>= 1; }
        return r;
    }

 // ---------------------------------------------------------
    // 统一的 CG-NTT 旋转因子生成器：通过模拟硬件洗牌路径来定位 TF
    // ---------------------------------------------------------
    static void BuildCgTwiddle_Unified(uint64_t* out, int n, uint64_t mod, uint64_t root, bool is_ntt) {
        std::vector<int> logical_idx(n);
        for(int i=0; i<n; i++) logical_idx[i] = i; // 初始逻辑索引

        for(int s=0; s<STAGE_NUM; s++) {
            // OpenFHE Negacyclic 参数：NTT 阶段 m 递增，INTT 阶段 m 递减
            int m = is_ntt ? (1 << s) : (n >> (s + 1));
            int t = is_ntt ? (n >> (s + 1)) : (1 << s);

            for(int i=0; i<CG_HALF_N; i++) {
                // 根据硬件读模式确定逻辑索引
                // NTT (DIT) 硬件读 i 和 i + N/2
                // INTT (DIF) 硬件读 2i 和 2i + 1
                int idx = is_ntt ? logical_idx[i] : logical_idx[2*i];
                int group = idx / (2 * t);
                
                uint64_t power = BitReverse(m + group, STAGE_NUM);
                uint64_t tf = Power(root, power, mod);


                int target_row = is_ntt ? s : (STAGE_NUM - 1 - s);
                out[target_row * CG_HALF_N + i] = tf;
            }

            // --- 核心：模拟硬件在 Stage 结束后的物理洗牌行为 ---
            std::vector<int> next(n);
            if(is_ntt) {
                // NTT 硬件执行 Perfect Shuffle 写回：2i <- i, 2i+1 <- i+N/2
                for(int i=0; i<CG_HALF_N; i++) {
                    next[2*i] = logical_idx[i];
                    next[2*i+1] = logical_idx[i + CG_HALF_N];
                }
            } else {
                // INTT 硬件执行 Perfect Unshuffle 写回：i <- 2i, i+N/2 <- 2i+1
                for(int i=0; i<CG_HALF_N; i++) {
                    next[i] = logical_idx[2*i];
                    next[i + CG_HALF_N] = logical_idx[2*i+1];
                }
            }
            logical_idx = next;
        }
    }
};

// =============================================================
// 4. FPGA Manager Class 
// =============================================================
class FpgaManager {
public:
    static FpgaManager& GetInstance() {
        static FpgaManager instance;
        return instance;
    }

    bool IsReady() const { return m_is_ready && !m_hks_digit_only; }

    // Fused mode owns the device context. All legacy fine-grained hooks stay on
    // CPU, including unsupported shapes, KeyMult and ModDown. Configure while idle.
    void SetHksDigitOnly(bool enabled) { m_hks_digit_only = enabled; }

    // Sized transport for INIT / HKS_DIGIT. Unlike Execute(), input, metadata and
    // output have different sizes. Errors propagate; never return partial output.
    bool HksDigitTransfer(uint8_t opcode, const std::vector<uint64_t>& in1,
                          const std::vector<uint64_t>& in2, std::vector<uint64_t>& out,
                          int alpha, int start) {
        const bool init = opcode == OP_INIT &&
            in1.size() == FPGA_LIMB_Q * 3 + MAX_LIMBS * CG_TF_SIZE &&
            in2.size() == FPGA_LIMB_P * 3 + MAX_LIMBS * CG_TF_SIZE && out.size() == 1;
        const bool digit = opcode == OP_HKS_DIGIT && alpha >= 1 && alpha <= 3 &&
            start >= 0 && start + alpha <= 3 &&
            in1.size() == static_cast<size_t>(alpha * FPGA_RING_DIM) &&
            in2.size() == 18 && out.size() == 5 * FPGA_RING_DIM;
        if (!init && !digit)
            throw std::invalid_argument("Invalid HKS_DIGIT transfer shape");
#ifdef OPENFHE_FPGA_ENABLE
        if (!m_is_ready) return false;
        using Clock = std::chrono::steady_clock;
        auto elapsed = [](Clock::time_point a, Clock::time_point b) {
            return std::chrono::duration_cast<std::chrono::microseconds>(b - a).count();
        };
        auto t0 = Clock::now();
        int64_t h2d, kernel, d2h;
        {
            auto b1 = xrt::bo(m_device, in1.size() * 8, m_kernel_top.group_id(0));
            auto b2 = xrt::bo(m_device, in2.size() * 8, m_kernel_top.group_id(1));
            auto bo = xrt::bo(m_device, out.size() * 8, m_kernel_top.group_id(2));
            b1.write(in1.data());
            b2.write(in2.data());
            bo.write(out.data());  // sentinel detects an old xclbin without opcode 8
            auto t1 = Clock::now();
            b1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
            b2.sync(XCL_BO_SYNC_BO_TO_DEVICE);
            bo.sync(XCL_BO_SYNC_BO_TO_DEVICE);
            auto t2 = Clock::now();
            auto run = m_kernel_top(b1, b2, bo, opcode, alpha, start);
            run.wait();
            auto t3 = Clock::now();
            bo.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
            auto t4 = Clock::now();
            bo.read(out.data());
            h2d = elapsed(t1, t2);
            kernel = elapsed(t2, t3);
            d2h = elapsed(t3, t4);
        }
        auto wall = elapsed(t0, Clock::now());
        auto& s = GetFpgaTransferStats();
        s.h2d_us += h2d; s.kernel_us += kernel; s.d2h_us += d2h;
        s.wall_us += wall; ++s.calls;
        auto& op = s.by_opcode[opcode];
        op.h2d_us += h2d; op.kernel_us += kernel; op.d2h_us += d2h;
        op.wall_us += wall; ++op.calls;
        return true;
#else
        return false;
#endif
    }

// ----------------------------------------------------------------------
// InitModuli: 接收 CPU 的 Roots -> Index -> Permute -> Pack
// ----------------------------------------------------------------------
    void InitModuli(const std::vector<uint64_t>& q_mods, const std::vector<uint64_t>& p_mods,
                    const std::vector<uint64_t>& q_roots, const std::vector<uint64_t>& p_roots) {
    #ifdef OPENFHE_FPGA_ENABLE
        if (!m_is_ready) return;
        
        size_t n_q = q_mods.size();
        size_t n_p = p_mods.size();
        size_t total_limbs = n_q + n_p;
        const int N = FPGA_RING_DIM;

        std::cout << "[Host] InitModuli: Q=" << n_q << ", P=" << n_p << std::endl;

        m_stored_moduli.clear();
        m_stored_moduli.insert(m_stored_moduli.end(), q_mods.begin(), q_mods.end());
        m_stored_moduli.insert(m_stored_moduli.end(), p_mods.begin(), p_mods.end());

        // 🔥 将从 CPU 传来的单位根合并
        std::vector<uint64_t> combined_roots;
        combined_roots.insert(combined_roots.end(), q_roots.begin(), q_roots.end());
        combined_roots.insert(combined_roots.end(), p_roots.begin(), p_roots.end());

        std::vector<uint64_t> S_vals(total_limbs), M_vals(total_limbs);
        for(size_t i=0; i<total_limbs; i++) {
            uint64_t p = m_stored_moduli[i];
            int pbits = 64 - __builtin_clzll(p);
            int S = pbits + 62;  // 全精度总移位量，不再除以 2
            unsigned __int128 power = (unsigned __int128)1 << S;
            uint64_t m = (uint64_t)(power / p);
            S_vals[i] = S;
            M_vals[i] = m;
            std::cout << "  [Barrett] idx=" << i << ": mod=" << p << ", S=" << S << ", m=" << m << std::endl;
        }

        // ---- CG-NTT 旋转因子预计算 ----
        // 每 limb 生成 [STAGE_NUM][CG_HALF_N] = 24576 个 TF，与 FPGA cg_twiddle[s][t] 顺序完全对齐。
        // 按 MAX_LIMBS 尺寸分配（FPGA OP_INIT 按 MAX_LIMBS × CG_TF_SIZE 搬运，不足槽位以 0 填充）。
        std::vector<uint64_t> all_ntt_twiddles(MAX_LIMBS * CG_TF_SIZE, 0);
        std::vector<uint64_t> all_intt_twiddles(MAX_LIMBS * CG_TF_SIZE, 0);

        for(size_t limb = 0; limb < total_limbs && limb < (size_t)MAX_LIMBS; limb++) {
            uint64_t mod = m_stored_moduli[limb];
            uint64_t psi = combined_roots[limb];
            uint64_t psi_inv = MathUtils::ModInverse(psi, mod);

            // 正向
            MathUtils::BuildCgTwiddle_Unified(
                all_ntt_twiddles.data() + limb * CG_TF_SIZE, N, mod, psi, true);
            
            // 逆向
            MathUtils::BuildCgTwiddle_Unified(
                all_intt_twiddles.data() + limb * CG_TF_SIZE, N, mod, psi_inv, false);
        }
        

        const int PARAMS_PER_LIMB = 3;
        // 重要：header 偏移必须与 FPGA 端 top.cpp 中 NTT_TF_BASE = LIMB_Q*3 / INTT_TF_BASE = LIMB_P*3
        // 完全一致，并按 MAX_LIMBS × CG_TF_SIZE 分配尾部 TF 区域（FPGA OP_INIT 按该尺寸搬运）。
        size_t buf1_size = FPGA_LIMB_Q * PARAMS_PER_LIMB + MAX_LIMBS * CG_TF_SIZE;
        std::vector<uint64_t> buf1_Q(buf1_size, 0);
        for(size_t i=0; i<n_q && i<(size_t)FPGA_LIMB_Q; i++) {
            buf1_Q[i]                      = m_stored_moduli[i];
            buf1_Q[FPGA_LIMB_Q + i]        = S_vals[i];
            buf1_Q[FPGA_LIMB_Q * 2 + i]    = M_vals[i];
        }
        memcpy(buf1_Q.data() + FPGA_LIMB_Q * PARAMS_PER_LIMB,
               all_ntt_twiddles.data(),
               MAX_LIMBS * CG_TF_SIZE * sizeof(uint64_t));

        size_t buf2_size = FPGA_LIMB_P * PARAMS_PER_LIMB + MAX_LIMBS * CG_TF_SIZE;
        std::vector<uint64_t> buf2_P(buf2_size, 0);
        for(size_t i=0; i<n_p && i<(size_t)FPGA_LIMB_P; i++) {
            size_t global_idx = n_q + i;
            buf2_P[i]                      = m_stored_moduli[global_idx];
            buf2_P[FPGA_LIMB_P + i]        = S_vals[global_idx];
            buf2_P[FPGA_LIMB_P * 2 + i]    = M_vals[global_idx];
        }
        memcpy(buf2_P.data() + FPGA_LIMB_P * PARAMS_PER_LIMB,
               all_intt_twiddles.data(),
               MAX_LIMBS * CG_TF_SIZE * sizeof(uint64_t));

        // 探针：如果这里打印出来了，说明准备工作完毕，即将呼叫硬件
        std::cout << "[DEBUG] Ready to allocate XRT buffers..." << std::endl;

        auto bo_1 = xrt::bo(m_device, buf1_Q.size() * sizeof(uint64_t), m_kernel_top.group_id(0));
        auto bo_2 = xrt::bo(m_device, buf2_P.size() * sizeof(uint64_t), m_kernel_top.group_id(1));
        auto bo_out = xrt::bo(m_device, 8, m_kernel_top.group_id(2)); 

        bo_1.write(buf1_Q.data());
        bo_2.write(buf2_P.data());
        
        bo_1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
        bo_2.sync(XCL_BO_SYNC_BO_TO_DEVICE);

        std::cout << "[Host] Launching Init Kernel..." << std::endl;
        auto run = m_kernel_top(bo_1, bo_2, bo_out, OP_INIT, 0, 0); 
        run.wait();
        
        std::cout << "[Host] FPGA Parameter Init Complete." << std::endl;
    #endif
    }
    
    void Execute(
        uint8_t opcode,
        const uint64_t* in1,
        const uint64_t* in2,
        uint64_t* out,
        int num_limbs,
        int mod_idx
    ) {
    #ifdef OPENFHE_FPGA_ENABLE
        if (!m_is_ready) return;
        try {
            using Clock = std::chrono::high_resolution_clock;
            using us    = std::chrono::microseconds;

            // 实验 D：整段墙钟从这里起，覆盖 xrt::bo 构造 + write + sync + kernel + read。
            // 与 h2d/kernel/d2h 三段之和的差 = 之前完全没记账的部分。
            auto t_call_start = Clock::now();
            int64_t h2d = 0, kern = 0, d2h = 0;

            // 内层作用域：让三个 xrt::bo 在 t_call_end 之前析构，墙钟才把 bo 释放
            // （驱动 unmap）也算进去。否则测出来又是一个下界，重蹈 182.8 us 的覆辙。
            {
            size_t size_bytes     = (size_t)num_limbs * FPGA_RING_DIM * sizeof(uint64_t);
            size_t out_size_bytes = size_bytes;
            auto bo_in1 = xrt::bo(m_device, size_bytes,     m_kernel_top.group_id(0));
            auto bo_out = xrt::bo(m_device, out_size_bytes, m_kernel_top.group_id(2));

            bo_in1.write(in1);

            auto bo_in2 = xrt::bo(m_device, size_bytes, m_kernel_top.group_id(1));
            if (in2) {
                bo_in2.write(in2);
            } else {
                bo_in2.write(in1);  // NTT/INTT: kernel ignores in2, but buffer must be valid
            }

            // --- Host → Device ---
            auto t_h2d_start = Clock::now();
            bo_in1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
            bo_in2.sync(XCL_BO_SYNC_BO_TO_DEVICE);
            auto t_h2d_end = Clock::now();

            // --- Kernel Exec ---
            auto t_kern_start = Clock::now();
            auto run = m_kernel_top(bo_in1, bo_in2, bo_out, opcode, num_limbs, mod_idx);
            run.wait();
            auto t_kern_end = Clock::now();

            // --- Device → Host ---
            auto t_d2h_start = Clock::now();
            bo_out.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
            auto t_d2h_end = Clock::now();

            bo_out.read(out);

            h2d  = std::chrono::duration_cast<us>(t_h2d_end  - t_h2d_start).count();
            kern = std::chrono::duration_cast<us>(t_kern_end - t_kern_start).count();
            d2h  = std::chrono::duration_cast<us>(t_d2h_end  - t_d2h_start).count();
            }  // <- bo_in1 / bo_in2 / bo_out 在此析构

            auto t_call_end  = Clock::now();
            int64_t wall     = std::chrono::duration_cast<us>(t_call_end - t_call_start).count();
            int64_t untimed  = wall - h2d - kern - d2h;

            // Accumulate into global transfer stats (total + per-opcode)
            auto& ts = GetFpgaTransferStats();
            ts.h2d_us    += h2d;
            ts.kernel_us += kern;
            ts.d2h_us    += d2h;
            ts.wall_us   += wall;
            ts.calls++;
            if (opcode < 9) {
                ts.by_opcode[opcode].h2d_us    += h2d;
                ts.by_opcode[opcode].kernel_us += kern;
                ts.by_opcode[opcode].d2h_us    += d2h;
                ts.by_opcode[opcode].wall_us   += wall;
                ts.by_opcode[opcode].calls++;
            }

            std::cout << "    [Trace] Opcode=" << (int)opcode
                      << "  H2D=" << h2d << "μs"
                      << "  Kernel=" << kern / 1000.0 << "ms"
                      << "  D2H=" << d2h << "μs"
                      << "  Wall=" << wall << "μs"
                      << "  未计时=" << untimed << "μs\n";
        } catch (const std::exception& e) {
            std::cerr << "[FPGA Exec Error] " << e.what() << std::endl;
        }
    #endif
    }

    int GetModIndex(uint64_t modulus) {
        auto it = std::find(m_stored_moduli.begin(), m_stored_moduli.end(), modulus);
        if (it != m_stored_moduli.end()) {
            return (int)std::distance(m_stored_moduli.begin(), it);
        } else {
            std::cerr << "[FPGA Warning] Modulus " << modulus << " not found! Using 0." << std::endl;
            return 0;
        }
    }

    void NttForwardOffload(
        const uint64_t* in, 
        uint64_t* out, 
        uint64_t modulus, 
        size_t n
    ) {
        std::cout << "=== [FPGA] Execute NTT ===" << std::endl;
        int mod_idx = GetModIndex(modulus);
        if (n != FPGA_RING_DIM) {
            std::cerr << "[FPGA Warning] NTT size " << n << " != FPGA_RING_DIM " << FPGA_RING_DIM << std::endl;
            return;
        }
        Execute(OP_NTT, in, nullptr, out, 1, mod_idx);

        if (std::getenv("OPENFHE_NTT_DUMP")) {
            const size_t dumpLen = std::min(n, (size_t)16);
            std::cerr << "[NTT_DUMP] NTT forward first " << dumpLen << " in/out (mod=" << modulus << "):" << std::endl;
            std::cerr << "  in :"; for (size_t i = 0; i < dumpLen; ++i) std::cerr << " " << in[i]; std::cerr << std::endl;
            std::cerr << "  out:"; for (size_t i = 0; i < dumpLen; ++i) std::cerr << " " << out[i]; std::cerr << std::endl;
        }
    }

    void NttInverseOffload(
        const uint64_t* in, 
        uint64_t* out, 
        uint64_t modulus, 
        size_t n
    ) {
        std::cout << "=== [FPGA] Execute INTT ===" << std::endl;
        int mod_idx = GetModIndex(modulus);
        Execute(OP_INTT, in, nullptr, out, 1, mod_idx);

        if (std::getenv("OPENFHE_NTT_DUMP")) {
            const size_t dumpLen = std::min(n, (size_t)16);
            std::cerr << "[NTT_DUMP] INTT first " << dumpLen << " in/out (mod=" << modulus << "):" << std::endl;
            std::cerr << "  in :"; for (size_t i = 0; i < dumpLen; ++i) std::cerr << " " << in[i]; std::cerr << std::endl;
            std::cerr << "  out:"; for (size_t i = 0; i < dumpLen; ++i) std::cerr << " " << out[i]; std::cerr << std::endl;
        }
    }

    void AutoOffload(
        const uint64_t* in,
        uint64_t* out,
        uint32_t k,
        uint32_t kinv,
        uint64_t modulus,
        size_t n
    ) {
        std::cout << "=== [FPGA] Execute Auto (k=" << k << ", kinv=" << kinv << ") ===" << std::endl;
        int mod_idx = GetModIndex(modulus);
        size_t meta_size = FPGA_RING_DIM;  // Execute expects in2 same size as poly
        std::vector<uint64_t> meta(meta_size, 0);
        meta[0] = (uint64_t)k;
        meta[1] = (uint64_t)kinv;
        Execute(OP_AUTO, in, meta.data(), out, 1, mod_idx);
    }

    // ============================================================
    // 单limb操作（保留向后兼容）
    // ============================================================
    void ModMultOffload(
        const uint64_t* a, 
        const uint64_t* b, 
        uint64_t* result, 
        uint64_t modulus, 
        size_t n
    ) {
        int mod_idx = GetModIndex(modulus);
        Execute(OP_MULT, a, b, result, 1, mod_idx);
    }

    void ModAddOffload(
        const uint64_t* a, 
        const uint64_t* b, 
        uint64_t* result, 
        uint64_t modulus, 
        size_t n
    ) {
        int mod_idx = GetModIndex(modulus);
        Execute(OP_ADD, a, b, result, 1, mod_idx);
    }
    
    void ModSubOffload(
        const uint64_t* a, 
        const uint64_t* b, 
        uint64_t* result, 
        uint64_t modulus, 
        size_t n
    ) {
        int mod_idx = GetModIndex(modulus);
        Execute(OP_SUB, a, b, result, 1, mod_idx);
    }

    // ============================================================
    // 批量limb操作（一次kernel调用处理所有limb）
    // 要求：所有limb的模数索引是连续的（从mod_idx_start开始）
    // ============================================================
    void ModOpBatchOffload(
        int opcode,             // OP_ADD, OP_SUB, or OP_MUL (注意是OP_MUL不是OP_MULT)
        const uint64_t* a,      // [numLimbs × ringDim]
        const uint64_t* b,      // [numLimbs × ringDim]
        uint64_t* result,       // [numLimbs × ringDim]
        int mod_idx_start,      // 起始模数索引
        size_t ringDim,
        size_t numLimbs
    ) {
    #ifdef OPENFHE_FPGA_ENABLE
        if (!m_is_ready) return;
        
        std::cout << "=== [FPGA] Batch Op=" << opcode << " numLimbs=" << numLimbs 
                  << " modStart=" << mod_idx_start << " ===" << std::endl;

        try {
            size_t total_size = numLimbs * ringDim * sizeof(uint64_t);
            
            auto bo_a = xrt::bo(m_device, total_size, m_kernel_top.group_id(0));
            auto bo_b = xrt::bo(m_device, total_size, m_kernel_top.group_id(1));
            auto bo_out = xrt::bo(m_device, total_size, m_kernel_top.group_id(2));

            // 一次性传输所有数据
            bo_a.write(a);
            bo_a.sync(XCL_BO_SYNC_BO_TO_DEVICE);
            
            bo_b.write(b);
            bo_b.sync(XCL_BO_SYNC_BO_TO_DEVICE);

            // 一次kernel调用处理所有limb
            // num_active_limbs = numLimbs, mod_index = mod_idx_start
            auto run = m_kernel_top(bo_a, bo_b, bo_out, opcode, (int)numLimbs, mod_idx_start);
            run.wait();

            // 一次性读取所有结果
            bo_out.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
            bo_out.read(result);
        } catch (const std::exception& e) {
            std::cerr << "[FPGA Batch Op Error] " << e.what() << std::endl;
        }
    #endif
    }

    // BConv with dynamic output moduli
    // Kernel: LIMB_Q rows × MAX_OUT_COLS columns (3×5)
    static const int KERNEL_LIMB_Q = 3;
    static const int KERNEL_MAX_OUT_COLS = 5;  // LIMB_Q + LIMB_P
    
    void BConvOffload(
        const uint64_t* x,              // 输入: [KERNEL_LIMB_Q × RING_DIM]
        const uint64_t* w,              // 权重: [KERNEL_LIMB_Q × KERNEL_MAX_OUT_COLS]
        const uint64_t* out_mod,        // 输出模数: [sizeP]
        uint64_t* result,               // 输出: [sizeP × RING_DIM]
        size_t ringDim,
        int sizeP                       // 实际输出列数
    ) {
    #ifdef OPENFHE_FPGA_ENABLE
        if (!m_is_ready) return;

        std::cout << "=== [FPGA] Execute BConv === sizeP=" << sizeP << std::endl;

        try {
            using Clock = std::chrono::high_resolution_clock;
            using us    = std::chrono::microseconds;

            // 实验 D：整段墙钟，口径与 Execute() 一致（含 bo 构造/析构 + meta 组装）
            auto t_call_start = Clock::now();
            int64_t h2d = 0, kern = 0, d2h = 0;
            {
            // Buffer sizes
            size_t in_size = KERNEL_LIMB_Q * ringDim * sizeof(uint64_t);
            size_t out_size = sizeP * ringDim * sizeof(uint64_t);

            // -------------------------------------------------------
            // FPGA Kernel (top.cpp OP_BCONV) 从 mem_in2 的布局：
            //   [0           .. LIMB_Q*MAX_OUT_COLS-1]  : weights  (15 个)
            //   [LIMB_Q*MAX_OUT_COLS .. +MAX_OUT_COLS-1]: out_mod  (5 个)
            //   [+MAX_OUT_COLS       .. +MAX_OUT_COLS-1]: S        (5 个)
            //   [+MAX_OUT_COLS       .. +MAX_OUT_COLS-1]: m_barrett(5 个)
            // 总共: 15 + 5 + 5 + 5 = 30 个 uint64_t
            // -------------------------------------------------------
            size_t weights_count = KERNEL_LIMB_Q * KERNEL_MAX_OUT_COLS;  // 15
            size_t total_meta = weights_count + 3 * KERNEL_MAX_OUT_COLS; // 15 + 15 = 30
            size_t meta_size = total_meta * sizeof(uint64_t);

            std::vector<uint64_t> meta_buffer(total_meta, 0);

            // (1) Copy weights [0..14]
            std::memcpy(meta_buffer.data(), w, weights_count * sizeof(uint64_t));

            // (2) Copy output moduli + 计算 Barrett 参数 S 和 m_barrett
            size_t mod_offset   = weights_count;                       // 15
            size_t S_offset = mod_offset + KERNEL_MAX_OUT_COLS;    // 20
            size_t m_offset     = S_offset + KERNEL_MAX_OUT_COLS;  // 25

            for (int i = 0; i < KERNEL_MAX_OUT_COLS; i++) {
                if (i < sizeP) {
                    uint64_t p = out_mod[i];
                    meta_buffer[mod_offset + i] = p;

                    // Barrett 参数：全精度总移位量 S = bitwidth(p) + 62
                    int pbits = 64 - __builtin_clzll(p);
                    int S = pbits + 62;
                    unsigned __int128 power = (unsigned __int128)1 << S;
                    uint64_t m = (uint64_t)(power / p);

                    meta_buffer[S_offset + i] = (uint64_t)S;   // S = 总移位量
                    meta_buffer[m_offset + i]     = m;              // m_barrett
                } else {
                    meta_buffer[mod_offset + i]   = 0;
                    meta_buffer[S_offset + i] = 0;
                    meta_buffer[m_offset + i]     = 0;
                }
            }

            auto bo_in   = xrt::bo(m_device, in_size,   m_kernel_top.group_id(0));
            auto bo_meta = xrt::bo(m_device, meta_size, m_kernel_top.group_id(1));
            auto bo_out  = xrt::bo(m_device, out_size,  m_kernel_top.group_id(2));

            bo_in.write(x);
            bo_meta.write(meta_buffer.data());

            // --- Host → Device ---
            auto t_h2d_start = Clock::now();
            bo_in.sync(XCL_BO_SYNC_BO_TO_DEVICE);
            bo_meta.sync(XCL_BO_SYNC_BO_TO_DEVICE);
            auto t_h2d_end = Clock::now();

            // --- Kernel Exec ---
            auto t_kern_start = Clock::now();
            auto run = m_kernel_top(bo_in, bo_meta, bo_out, OP_BCONV, sizeP, 0);
            run.wait();
            auto t_kern_end = Clock::now();

            // --- Device → Host ---
            auto t_d2h_start = Clock::now();
            bo_out.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
            auto t_d2h_end = Clock::now();

            bo_out.read(result);

            h2d  = std::chrono::duration_cast<us>(t_h2d_end  - t_h2d_start).count();
            kern = std::chrono::duration_cast<us>(t_kern_end - t_kern_start).count();
            d2h  = std::chrono::duration_cast<us>(t_d2h_end  - t_d2h_start).count();
            }  // <- bo_in / bo_meta / bo_out 在此析构

            auto t_call_end = Clock::now();
            int64_t wall    = std::chrono::duration_cast<us>(t_call_end - t_call_start).count();
            int64_t untimed = wall - h2d - kern - d2h;

            // Accumulate into global transfer stats
            // （原来 BConvOffload 只累加了总计、漏了 by_opcode，一并补上）
            auto& ts = GetFpgaTransferStats();
            ts.h2d_us    += h2d;
            ts.kernel_us += kern;
            ts.d2h_us    += d2h;
            ts.wall_us   += wall;
            ts.calls++;
            ts.by_opcode[OP_BCONV].h2d_us    += h2d;
            ts.by_opcode[OP_BCONV].kernel_us += kern;
            ts.by_opcode[OP_BCONV].d2h_us    += d2h;
            ts.by_opcode[OP_BCONV].wall_us   += wall;
            ts.by_opcode[OP_BCONV].calls++;

            std::cout << "    [Trace] BConv(sizeP=" << sizeP << ")"
                      << "  H2D=" << h2d << "μs"
                      << "  Kernel=" << kern / 1000.0 << "ms"
                      << "  D2H=" << d2h << "μs"
                      << "  Wall=" << wall << "μs"
                      << "  未计时=" << untimed << "μs\n";
        } catch (const std::exception& e) {
            std::cerr << "[FPGA BConv Error] " << e.what() << std::endl;
        }
    #endif
    }


private:
#ifdef OPENFHE_FPGA_ENABLE
    xrt::device m_device;
    xrt::kernel m_kernel_top;
#endif
    bool m_is_ready = false;                   // set to true only when FPGA is successfully initialized
    bool m_hks_digit_only = false;
    std::vector<uint64_t> m_stored_moduli;
    std::vector<uint64_t> m_stored_roots; // <--- 新增

    FpgaManager() {
#ifdef OPENFHE_FPGA_ENABLE
        try {
            m_device = xrt::device(0); 
            std::cout << "[FPGA] Device connected: " << m_device.get_info<xrt::info::device::name>() << std::endl;
            auto uuid = m_device.load_xclbin(GetXclbinPath());
            m_kernel_top = xrt::kernel(m_device, uuid, "Top");
            m_is_ready = true;
            std::cout << "[FPGA] Unified Top Kernel Loaded.\n";
        } catch (const std::exception& e) { 
            std::cerr << "[FPGA Setup Error] " << e.what() << std::endl;
            m_is_ready = false; 
        }
#endif
    }
    FpgaManager(const FpgaManager&) = delete;
    void operator=(const FpgaManager&) = delete;
};

#endif // _FPGA_MANAGER_H_
