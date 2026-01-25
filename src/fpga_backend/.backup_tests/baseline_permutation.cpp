#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <bitset>
#include <sstream> // 必须包含这个头文件

using namespace std;

// ================= 配置区域 =================
#define RING_DIM 65536      // 环维度 N
#define PARALLELISM 32      // 硬件并行度 (Ports/Lanes), 增加到32更符合实际大N场景
#define FPGA_FREQ_MHZ 300   // FPGA 频率

// ================= 1. Naive (软件参考) 算法 =================
// 这是标准的数学实现，用于生成 "Golden Reference"
void automorphism_naive(const vector<int>& input, vector<int>& output, int r, int mod, int n) {
    // 计算 k = 5^r. 注意：实际应用中要处理大数溢出，这里演示 r=1
    long long k = 1;
    for(int i=0; i<r; i++) k *= 5; 

    for (int i = 0; i < n; i++) {
        long long raw_idx = (long long)i * k;
        int dest = raw_idx % (2 * n);
        
        if (dest < n) {
            output[dest] = input[i];
        } else {
            // CKKS Automorphism 的变号性质: X^(N+j) = -X^j
            int val = input[i];
            output[dest - n] = (val == 0) ? 0 : (mod - val);
        }
    }
}

// ================= 2. 硬件网络仿真模型 =================

struct Packet {
    int value;          
    int original_idx;   
    int target_idx;     
    bool needs_negate;  
    
    string to_string() const {
        stringstream ss;
        // 如果数据量大，简化打印格式
        ss << "[" << value << "]"; 
        return ss.str();
    }
};

class MultiStagePermutationSim {
public:
    int N;
    int stages;
    vector<Packet> wires; 

    MultiStagePermutationSim(int n) : N(n), stages(log2(n)) {
        wires.resize(N);
    }

    // AGU: 计算路由信息
    void init_agu_routing(const vector<int>& input, int r) {
        long long k = 1;
        for(int i=0; i<r; i++) k *= 5;

        for (int i = 0; i < N; i++) {
            wires[i].value = input[i];
            wires[i].original_idx = i;
            
            long long raw_idx = (long long)i * k;
            int dest = raw_idx % (2 * N);
            
            if (dest < N) {
                wires[i].target_idx = dest;
                wires[i].needs_negate = false;
            } else {
                wires[i].target_idx = dest - N;
                wires[i].needs_negate = true;
            }
        }
    }

    // Perfect Shuffle 连线
    void perfect_shuffle() {
        vector<Packet> next_wires(N);
        int half = N / 2;
        for (int i = 0; i < half; i++) {
            next_wires[2 * i] = wires[i];         
            next_wires[2 * i + 1] = wires[half + i]; 
        }
        wires = next_wires;
    }

    // Switch Stage
    void switch_stage(int stage_bit) {
        for (int i = 0; i < N; i += 2) {
            Packet& p1 = wires[i];
            Packet& p2 = wires[i+1];

            // 检查目标地址的对应位
            // Omega Network / Butterfly 路由逻辑
            bool p1_wants_down = (p1.target_idx >> (stages - 1 - stage_bit)) & 1;
            bool p2_wants_up   = !((p2.target_idx >> (stages - 1 - stage_bit)) & 1);

            if (p1_wants_down || p2_wants_up) {
                swap(p1, p2); 
            }
        }
    }

    // 运行流水线
    void run_pipeline(int mod, bool verbose) {
        if (verbose) {
            cout << ">>> Initial State: "; print_wires();
        }

        for (int s = 0; s < stages; s++) {
            perfect_shuffle();
            switch_stage(s);
            if (verbose) {
                cout << "=== Stage " << s + 1 << ": "; print_wires();
            }
        }

        // 最后一级变号处理
        for (int i = 0; i < N; i++) {
            if (wires[i].needs_negate) {
                int val = wires[i].value;
                wires[i].value = (val == 0) ? 0 : (mod - val);
            }
        }

        if (verbose) {
            cout << ">>> Final Result:  "; print_wires();
        }
    }

    // 辅助函数：提取结果以便对比
    void get_results(vector<int>& output) {
        output.resize(N);
        for(int i=0; i<N; i++) {
            output[i] = wires[i].value;
        }
    }

    void print_wires() {
        // 防止 N=65536 时刷屏
        if (N > 32) {
            cout << "(Too strict to print " << N << " wires...)" << endl;
            return;
        }
        for (int i = 0; i < N; i++) cout << wires[i].to_string() << " ";
        cout << endl;
    }

    // 统计报告
    void report_stats(long long naive_cycles) {
        int pipeline_latency = stages; 
        int throughput_cycles = N / PARALLELISM;
        int hw_total_cycles = pipeline_latency + throughput_cycles;

        double hw_time_us = (double)hw_total_cycles / FPGA_FREQ_MHZ;
        // 假设 CPU 也是 300MHz (为了公平对比周期效率)，或者你可以设为 3000MHz
        double naive_time_us = (double)naive_cycles / FPGA_FREQ_MHZ; 

        cout << "\n==============================================" << endl;
        cout << "           Performance Comparison             " << endl;
        cout << "==============================================" << endl;
        cout << "Parameters:" << endl;
        cout << "  N = " << N << ", Ports = " << PARALLELISM << ", Freq = " << FPGA_FREQ_MHZ << " MHz" << endl;
        cout << "----------------------------------------------" << endl;
        cout << "[1] Naive (Sequential Software)" << endl;
        cout << "  Complexity      : O(N)" << endl;
        cout << "  Total Cycles    : " << naive_cycles << " (Est. 1 cycle/op ideal)" << endl;
        cout << "  Est. Time       : " << naive_time_us << " us" << endl;
        cout << "----------------------------------------------" << endl;
        cout << "[2] Hardware (Permutation Network)" << endl;
        cout << "  Complexity      : O(N / Ports)" << endl;
        cout << "  Pipeline Depth  : " << pipeline_latency << " cycles" << endl;
        cout << "  Throughput Time : " << throughput_cycles << " cycles" << endl;
        cout << "  Total Cycles    : " << hw_total_cycles << endl;
        cout << "  Est. Time       : " << hw_time_us << " us" << endl;
        cout << "----------------------------------------------" << endl;
        cout << ">>> Speedup       : " << fixed << setprecision(2) << naive_time_us / hw_time_us << "x" << endl;
        cout << "==============================================" << endl;
    }
};

int main() {
    int N = RING_DIM; 
    int r = 1;          
    int mod = 998244353;

    // 1. 初始化数据
    vector<int> input_data(N);
    for(int i=0; i<N; i++) input_data[i] = i;

    vector<int> naive_result(N);
    vector<int> hw_result(N);

    // 2. 运行 Naive 算法 (Golden Reference)
    cout << "Running Naive implementation..." << endl;
    automorphism_naive(input_data, naive_result, r, mod, N);

    // 3. 运行 硬件仿真
    cout << "Running Hardware simulation..." << endl;
    MultiStagePermutationSim sim(N);
    sim.init_agu_routing(input_data, r);
    
    // 如果 N 太大，关闭 verbose 打印
    bool verbose = (N <= 32); 
    sim.run_pipeline(mod, verbose);
    sim.get_results(hw_result);

    // 4. 验证结果 (Verification)
    cout << "\nVerifying results...";
    bool pass = true;
    for(int i=0; i<N; i++) {
        if (naive_result[i] != hw_result[i]) {
            pass = false;
            cout << "\n[FAIL] Mismatch at index " << i 
                 << ": Expected " << naive_result[i] 
                 << ", Got " << hw_result[i] << endl;
            // 只打印前几个错误
            if (i > 5) break; 
        }
    }

    if (pass) {
        cout << " [PASS]" << endl;
    } else {
        return -1;
    }

    // 5. 输出对比报告
    // 假设 Naive 版本的周期数大约是 N (非常理想化的串行处理)
    sim.report_stats(N);

    return 0;
}