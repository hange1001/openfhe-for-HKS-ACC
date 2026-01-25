#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <cstring>
#include "./include/bconv.h"

// ------------------------------------------------------------------
// 软件参考模型 (Golden Reference) for Compute_BConv
// ------------------------------------------------------------------
// Compute_BConv 执行的操作 (匹配 bconv_systolic):
// 对于每个 ring element 位置 (row, col) 和每个输出 limb p:
//   out[mod_idx_offset + p][row][col] = (sum_{q=0}^{LIMB_Q-1} in_x[q][row][col] * in_w[q][p]) % MODULUS[mod_idx_offset + p]
//
// 注意: 权重是 2D 数组 in_w[LIMB_Q][LIMB_P]
void bconv_ref_3d(
    uint64_t in_x[LIMB_Q + LIMB_P][SQRT][SQRT],
    const uint64_t in_w[LIMB_Q][LIMB_P],
    const uint64_t MODULUS[LIMB_Q + LIMB_P],
    int mod_idx_offset
) {
   
    // 对每个输出 limb p
    for (int p = 0; p < LIMB_P; ++p) {
        uint64_t mod = MODULUS[LIMB_Q + p];
        
        // 对每个 ring element 位置
        for (int row = 0; row < SQRT; ++row) {
            for (int col = 0; col < SQRT; ++col) {
                unsigned __int128 acc = 0;
                
                // 累加所有输入 limb q 的贡献
                for (int q = 0; q < LIMB_Q; ++q) {
                    uint64_t x_val = in_x[q][row][col];
                    uint64_t w_val = in_w[q][p];
                    acc += (unsigned __int128)x_val * (unsigned __int128)w_val;
                }
                
                // 存储到 out[mod_idx_offset + p]
                in_x[LIMB_Q + p][row][col] = (uint64_t)(acc);
            }
        }
    }
}

// ------------------------------------------------------------------
// 辅助函数：打印 3D 数组的一个 limb
// ------------------------------------------------------------------
void print_limb(const char* name, const uint64_t data[SQRT][SQRT], int max_rows = 4) {
    std::cout << name << " [" << SQRT << "x" << SQRT << "] (First " << max_rows << " rows):\n";
    for (int r = 0; r < max_rows && r < SQRT; ++r) {
        for (int c = 0; c < 8 && c < SQRT; ++c) {
            std::cout << std::setw(8) << data[r][c] << " ";
        }
        if (SQRT > 8) std::cout << "...";
        std::cout << "\n";
    }
    std::cout << std::endl;
}

int main() {
    std::cout << "=============================================" << std::endl;
    std::cout << "   Testbench for Compute_BConv (3D Arrays)   " << std::endl;
    std::cout << "   Using FIXED dimensions:" << std::endl;
    std::cout << "   LIMB_Q=" << LIMB_Q << ", LIMB_P=" << LIMB_P << std::endl;
    std::cout << "   SQRT=" << SQRT << " (RING_DIM=" << RING_DIM << ")" << std::endl;
    std::cout << "=============================================" << std::endl;

    // 1. 分配数组
    // in_x: [LIMB_Q + LIMB_P][SQRT][SQRT] - 只使用前 LIMB_Q 个 limb
    // in_w: [LIMB_Q][LIMB_P] - 2D 权重数组
    // out:  [LIMB_Q + LIMB_P][SQRT][SQRT] - 输出
    static uint64_t in_x[LIMB_Q + LIMB_P][SQRT][SQRT];
    static uint64_t in_hw[LIMB_Q + LIMB_P][SQRT][SQRT];
    static uint64_t in_w[LIMB_Q][LIMB_P];
    static uint64_t MODULUS[LIMB_Q + LIMB_P];

    // 2. 设置模数
    std::cout << "Setting up modulus array..." << std::endl;
    for (int i = 0; i < LIMB_Q + LIMB_P; ++i) {
        MODULUS[i] = 1000000007ULL + i * 1000;  // 不同的模数
    }
    std::cout << "MODULUS[0] = " << MODULUS[0] << std::endl;

    // 3. 初始化输入数据 in_x
    std::cout << "Initializing input X [" << (LIMB_Q + LIMB_P) << "][" << SQRT << "][" << SQRT << "]..." << std::endl;
    srand(42);
    for (int l = 0; l < LIMB_Q; ++l) {  // 只初始化前 LIMB_Q 个 limb
        for (int r = 0; r < SQRT; ++r) {
            for (int c = 0; c < SQRT; ++c) {
                in_x[l][r][c] = rand() % 100 + 1;  // 1-100
                in_hw[l][r][c] = in_x[l][r][c];
            }
        }
    }
    // 其余 limb 清零
    for (int l = LIMB_Q; l < LIMB_Q + LIMB_P; ++l) {
        memset(in_x[l], 0, sizeof(in_x[l]));
        memset(in_hw[l], 0, sizeof(in_hw[l]));
    }

    // 4. 初始化权重 in_w
    std::cout << "Initializing weights W [" << LIMB_Q << "][" << LIMB_P << "]..." << std::endl;
    for (int q = 0; q < LIMB_Q; ++q) {
        for (int p = 0; p < LIMB_P; ++p) {
            in_w[q][p] = ((q * LIMB_P + p) % 50) + 1;  // 1-50
        }
    }



    // 6. 运行软件参考模型
    std::cout << "\nRunning Reference Model..." << std::endl;
    int num_active_limbs = LIMB_Q;
    int mod_idx_offset = 0;
    bconv_ref_3d(in_x, in_w, MODULUS, mod_idx_offset);
    std::cout << "Reference model completed." << std::endl;

    // 7. 运行硬件 Kernel (Compute_BConv)
    std::cout << "\nRunning Hardware Kernel (Compute_BConv)..." << std::endl;
    Compute_BConv(
        in_hw,           // input
        in_w,           // weights (2D)
        MODULUS,        // modulus array
        num_active_limbs,
        mod_idx_offset
    );
    std::cout << "Hardware kernel completed." << std::endl;

    // 8. 结果比对
    std::cout << "\nVerifying Results..." << std::endl;
    int errors = 0;
    int total_checked = 0;
    
    for (int l = 0; l < LIMB_Q + LIMB_P; ++l) {
        for (int r = 0; r < SQRT; ++r) {
            for (int c = 0; c < SQRT; ++c) {
                if (in_hw[l][r][c] != in_x[l][r][c]) {
                    errors++;
                    std::cout << "Error at in_hw[" << l << "][" << r << "][" << c << "] = " << in_hw[l][r][c] << " != " << in_x[l][r][c] << std::endl;
                }
            }
        }
    }
    std::cout << "Total errors: " << errors << std::endl;
    std::cout << "Total checked: " << total_checked << std::endl;
    std::cout << "Error rate: " << (double)errors / total_checked << std::endl;

    return errors;
}
