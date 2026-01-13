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

#define MAX_LIMBS 5 
#define FPGA_RING_DIM  4096
#define STAGE_NUM 12 

inline std::string GetXclbinPath() {
    const char* mode = std::getenv("XCL_EMULATION_MODE");
    // 请根据你的实际路径修改 Base Path
    std::string base = "/home/timhan/FHE/openfhe/src/fpga_backend/";
    if (mode && std::string(mode) == "sw_emu") 
        return base + "fhe_kernels_sw_emu.xclbin";
    else if (mode && std::string(mode) == "hw_emu") 
        return base + "fhe_kernels_hw_emu.xclbin";
    else if (mode && std::string(mode) == "hw") 
        return base + "fhe_kernels_hw.xclbin";
    else return base + "fhe_kernels_sw_emu.xclbin";
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


    static std::vector<int> GenerateTwiddleIndices(int n) {
        std::vector<int> index = {0};
        for (int i = 0; i < STAGE_NUM; ++i) {
            std::vector<int> index_temp;
            int offset = n / (1 << (i + 1));
            for (int j : index) {
                index_temp.push_back(j + offset);
            }
            index.insert(index.end(), index_temp.begin(), index_temp.end());
        }
        index.erase(index.begin());
        return index;
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

    bool IsReady() const { return m_is_ready; }

 // ----------------------------------------------------------------------
    // InitModuli: Compute -> Index -> Permute -> Pack
    // ----------------------------------------------------------------------
    void InitModuli(const std::vector<uint64_t>& q_mods, const std::vector<uint64_t>& p_mods) {
    #ifdef OPENFHE_FPGA_ENABLE
        if (!m_is_ready) return;
        
        // 【核心修复】硬件被锁定为 4 (由日志 MODULUS[3]=16 证实)
        // 我们必须主动填充空洞，把 P 模数顶到 Index 4 的位置
        const size_t HARDWARE_LIMB_Q = 4; 

        size_t n_q_real = q_mods.size();
        size_t n_p = p_mods.size();
        
        // 自动计算对齐：如果真实 Q < 4，则对齐到 4
        // 这样就把 HARDWARE_LIMB_Q 用起来了，不会报 unused variable
        size_t n_q_aligned = std::max(n_q_real, HARDWARE_LIMB_Q);
        
        // 总 Limb 数 (对齐后)
        size_t total_limbs_aligned = n_q_aligned + n_p; 
        
        const int N = FPGA_RING_DIM;

        std::cout << "[Host] InitModuli (Auto-Padding): Real Q=" << n_q_real 
                  << " -> Aligned to " << n_q_aligned 
                  << ". P starts at Index " << n_q_aligned << std::endl;

        // 1. 重构 stored_moduli (用于 GetModIndex 查找)
        // 目标布局: [Q0, Q1, Q2, PAD, P0, P1] -> P0 索引变为 4
        m_stored_moduli.clear();
        
        // 1.1 填入真实 Q
        m_stored_moduli.insert(m_stored_moduli.end(), q_mods.begin(), q_mods.end());
        
        // 1.2 填充空洞 (Padding)
        // 填入 1 是为了防止某些数学库在预计算时除以 0
        if (n_q_real < n_q_aligned) {
            size_t padding = n_q_aligned - n_q_real;
            for(size_t k=0; k<padding; k++) {
                m_stored_moduli.push_back(1); 
            }
        }

        // 1.3 填入 P (此时 P 的起始索引就被顶到了 4，与硬件对齐！)
        m_stored_moduli.insert(m_stored_moduli.end(), p_mods.begin(), p_mods.end());

        // 2. Compute Constants K and M (Barrett Reduction)
        std::vector<uint64_t> K_vals(total_limbs_aligned), M_vals(total_limbs_aligned);
        for(size_t i=0; i<total_limbs_aligned; i++) {
            uint64_t p = m_stored_moduli[i];
            if (p <= 1) { 
                K_vals[i] = 0; M_vals[i] = 0; 
                continue; 
            }
            int k = (int)std::ceil(std::log2((double)p));
            unsigned __int128 power = (unsigned __int128)1 << (2 * k);
            uint64_t m = (uint64_t)(power / p);
            K_vals[i] = k;
            M_vals[i] = m;
        }

        // 3. Compute Twiddles (Base on Aligned Total Limbs)
        std::vector<uint64_t> all_ntt_twiddles(total_limbs_aligned * N);
        std::vector<uint64_t> all_intt_twiddles(total_limbs_aligned * N);

        for(size_t limb = 0; limb < total_limbs_aligned; limb++) {
            uint64_t mod = m_stored_moduli[limb];
            if (mod <= 1) continue; 

            uint64_t psi = MathUtils::Find2NthRootOfUnity(mod, N); 
            uint64_t psi_inv = MathUtils::ModInverse(psi, mod);

            uint64_t w = 1;
            uint64_t w_inv = 1;

            for(int i = 0; i < N; i++) {
                all_ntt_twiddles[limb * N + i] = w;
                all_intt_twiddles[limb * N + i] = w_inv;
                w = ((unsigned __int128)w * psi) % mod;
                w_inv = ((unsigned __int128)w_inv * psi_inv) % mod;
            }
        }

        // 4. Permute Twiddles
        std::vector<int> perm_index = MathUtils::GenerateTwiddleIndices(N);
        std::vector<uint64_t> permuted_ntt(total_limbs_aligned * N);
        std::vector<uint64_t> permuted_intt(total_limbs_aligned * N);

        for (size_t limb = 0; limb < total_limbs_aligned; limb++) {
            size_t base_offset = limb * N;
            for (size_t i = 0; i < perm_index.size(); ++i) {
                permuted_ntt[base_offset + i]  = all_ntt_twiddles[base_offset + perm_index[i]];
                permuted_intt[base_offset + i] = all_intt_twiddles[base_offset + perm_index[i]];
            }
        }

        // ============================================================
        // 5. Pack Buffers (发送给 FPGA)
        // ============================================================
        const int PARAMS_PER_LIMB = 3; // MOD, K, M
        
        // --- Buffer 1 (Q Params + NTT Twiddles) ---
        // 关键：必须按照 n_q_aligned (4) 来打包头部
        size_t buf1_size = n_q_aligned * PARAMS_PER_LIMB + total_limbs_aligned * N;
        std::vector<uint64_t> buf1_Q(buf1_size);

        for(size_t i=0; i<n_q_aligned; i++) {
            buf1_Q[i] = m_stored_moduli[i];                   // MOD
            buf1_Q[n_q_aligned + i] = K_vals[i];              // K
            buf1_Q[n_q_aligned*2 + i] = M_vals[i];            // M
        }
        // Copy NTT Twiddles (offset 会自动对齐)
        memcpy(
            buf1_Q.data() + n_q_aligned * PARAMS_PER_LIMB, 
            permuted_ntt.data(), 
            total_limbs_aligned * N * sizeof(uint64_t)
        );

        // --- Buffer 2 (P Params + INTT Twiddles) ---
        size_t buf2_size = n_p * PARAMS_PER_LIMB + total_limbs_aligned * N;
        std::vector<uint64_t> buf2_P(buf2_size);

        for(size_t i=0; i<n_p; i++) {
            // P 在 stored_moduli 中的索引从 4 开始
            size_t global_idx = n_q_aligned + i; 
            buf2_P[i] = m_stored_moduli[global_idx];        // MOD
            buf2_P[n_p + i] = K_vals[global_idx];           // K
            buf2_P[n_p*2 + i] = M_vals[global_idx];         // M
        }
        // Copy INTT Twiddles
        memcpy(buf2_P.data() + n_p * PARAMS_PER_LIMB, 
               permuted_intt.data(), 
               total_limbs_aligned * N * sizeof(uint64_t)
            );

        // 7. Send to FPGA (OP_INIT)
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
    
    // --- 其他函数 (Execute, Wrappers) 保持不变 ---
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
            size_t size_bytes = (size_t)num_limbs * FPGA_RING_DIM * sizeof(uint64_t);
            
            auto bo_in1 = xrt::bo(m_device, size_bytes, m_kernel_top.group_id(0));
            auto bo_out = xrt::bo(m_device, size_bytes, m_kernel_top.group_id(2));

            bo_in1.write(in1);
            bo_in1.sync(XCL_BO_SYNC_BO_TO_DEVICE);

            xrt::bo bo_in2;
            if (in2 && in2 != in1) {
                bo_in2 = xrt::bo(m_device, size_bytes, m_kernel_top.group_id(1));
                bo_in2.write(in2);
                bo_in2.sync(XCL_BO_SYNC_BO_TO_DEVICE);
            } else { 
                bo_in2 = bo_in1; 
            }

            auto run = m_kernel_top(bo_in1, bo_in2, bo_out, opcode, num_limbs, mod_idx);
            run.wait();
            

            bo_out.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
            bo_out.read(out);
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

    void NttForwardOffload(const uint64_t* in, uint64_t* out, uint64_t modulus, size_t n) {
        std::cout << "=== [FPGA] Execute NTT ===" << std::endl;
        int mod_idx = GetModIndex(modulus);
        Execute(OP_NTT, in, nullptr, out, 1, mod_idx); 
    }

    void NttInverseOffload(const uint64_t* in, uint64_t* out, uint64_t modulus, size_t n) {
        std::cout << "=== [FPGA] Execute INTT ===" << std::endl;
        int mod_idx = GetModIndex(modulus);
        Execute(OP_INTT, in, nullptr, out, 1, mod_idx);
    }

    void ModMultOffload(const uint64_t* a, const uint64_t* b, uint64_t* result, uint64_t modulus, size_t n) {
        std::cout << "=== [FPGA] Execute Mult ===" << std::endl;
        int mod_idx = GetModIndex(modulus);
        Execute(OP_MULT, a, b, result, 1, mod_idx);
    }

    void ModAddOffload(const uint64_t* a, const uint64_t* b, uint64_t* result, uint64_t modulus, size_t n) {
        std::cout << "=== [FPGA] Execute Add ===" << std::endl;
        int mod_idx = GetModIndex(modulus);
        Execute(OP_ADD, a, b, result, 1, mod_idx);
    }
    
    void ModSubOffload(const uint64_t* a, const uint64_t* b, uint64_t* result, uint64_t modulus, size_t n) {
        std::cout << "=== [FPGA] Execute Sub ===" << std::endl;
        int mod_idx = GetModIndex(modulus);
        Execute(OP_SUB, a, b, result, 1, mod_idx);
    }

    //TODO:testing
    // void BConvOffload(...) n是什么，要加吗
    void BConvOffload(const uint64_t* x,const uint64_t* w, uint64_t* result, uint64_t modulus,size_t n, size_t num_q_limbs) {
        std::cout << "=== [FPGA] Execute BConv ===" << std::endl;
        std::cout << " num_limbs=" << num_q_limbs << std::endl;

        int mod_idx = GetModIndex(modulus);
        Execute(OP_BCONV, x, w, result, num_q_limbs, mod_idx);
    }
    //看top.cpp的bconv_systolic函数最后一个参数是num_limbs，这里暂时写1

private:
#ifdef OPENFHE_FPGA_ENABLE
    xrt::device m_device;
    xrt::kernel m_kernel_top;
#endif
    bool m_is_ready = false;
    std::vector<uint64_t> m_stored_moduli;

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