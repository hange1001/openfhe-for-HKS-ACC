// =============================================================================
// host_overhead_bench.cpp — 实验 D 的【离线可测分量】
//
// 独立程序，不参与 csim / 不进 OpenFHE 构建。编译运行都不需要板卡和 XRT：
//   g++ -O2 -std=c++17 -o host_overhead_bench host_overhead_bench.cpp && ./host_overhead_bench
//
// -----------------------------------------------------------------------------
// 背景（task.yaml step 1.2 实验 D）
// -----------------------------------------------------------------------------
// FpgaManager::Execute() 现在只给三段计时：h2d / kernel / d2h = 40.6 / 123.0 / 19.2 us。
// 没被计进去的是：
//   (1) xrt::bo 构造 x3 + 析构 x3   —— 驱动 ioctl + pin 内存 + mmap
//   (2) bo.write(in1) / bo.write(in2) / bo_out.read(out) —— 纯 memcpy
//   (3) BConv hook 里的 gather/scatter（dcrtpoly-impl.h:1181-1223）——
//       每元素一次 ModMulFastConst，这一段连 h2d 都没算进去
//
// (1) 必须上板才能测。(2)(3) 是纯 CPU 侧，本程序把它们量出来。
//
// 它能回答的问题：预测的「未计时开销 30-80 us/call」里，memcpy 和 gather
// 各占多少？如果 (2)(3) 加起来就逼近 30 us，那 buffer 池化（只消除 (1)）
// 的收益上界比预期小；如果 (2)(3) 只有几 us，那 (1) 是大头，池化值得做。
// =============================================================================

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <random>
#include <vector>

// --- 与 src/fpga_backend/include/define.h 保持一致 -------------------------
static constexpr int RING_DIM = 1 << 12;  // 4096
static constexpr int LIMB_Q   = 3;
static constexpr int LIMB_P   = 2;

using Clock = std::chrono::steady_clock;

// 取中位数而非均值：memcpy 受调度抖动影响大，均值会被离群点拉偏
static double median_us(std::vector<double>& v) {
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

struct Result {
    const char* name;
    size_t      bytes;
    double      us;
};

// 防止 -O2 把整个循环优化掉
static volatile uint64_t g_sink = 0;

// -----------------------------------------------------------------------------
// (2) bo.write() / bo.read() 的等价物：一次 memcpy
// -----------------------------------------------------------------------------
static double bench_memcpy(size_t bytes, int iters) {
    std::vector<uint8_t> src(bytes), dst(bytes);
    std::mt19937_64 rng(12345);
    for (auto& b : src) b = (uint8_t)rng();

    std::memcpy(dst.data(), src.data(), bytes);  // warm up / 触发缺页
    std::vector<double> samples;
    samples.reserve(iters);
    for (int i = 0; i < iters; ++i) {
        auto t0 = Clock::now();
        std::memcpy(dst.data(), src.data(), bytes);
        auto t1 = Clock::now();
        g_sink += dst[i % bytes];
        samples.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    return median_us(samples);
}

// -----------------------------------------------------------------------------
// (3a) BConv hook 的输入准备：dcrtpoly-impl.h:1181-1188
//      每元素一次 ModMulFastConst（Shoup 形式的模乘）+ 一次跨 limb 的 gather 写
// -----------------------------------------------------------------------------
static inline uint64_t mul_mod_shoup(uint64_t x, uint64_t w, uint64_t w_precon, uint64_t q) {
    // OpenFHE ModMulFastConst 的等价形式：用预计算的 w_precon 省掉一次除法
    uint64_t hi = (uint64_t)(((__uint128_t)x * w_precon) >> 64);
    uint64_t r  = x * w - hi * q;
    return (r >= q) ? r - q : r;
}

static double bench_bconv_gather(int size_q, int iters) {
    std::vector<std::vector<uint64_t>> towers(size_q, std::vector<uint64_t>(RING_DIM));
    std::vector<uint64_t> flat((size_t)size_q * RING_DIM);
    std::vector<uint64_t> q(size_q), w(size_q), wp(size_q);

    std::mt19937_64 rng(999);
    for (int i = 0; i < size_q; ++i) {
        q[i]  = (1ULL << 49) - 2 * (i + 1) * 1024 + 1;  // 量级贴近真实 NTT-friendly 素数
        w[i]  = rng() % q[i];
        wp[i] = (uint64_t)(((__uint128_t)w[i] << 64) / q[i]);
        for (int r = 0; r < RING_DIM; ++r) towers[i][r] = rng() % q[i];
    }

    std::vector<double> samples;
    samples.reserve(iters);
    for (int it = 0; it < iters; ++it) {
        auto t0 = Clock::now();
        for (int i = 0; i < size_q; ++i)
            for (int r = 0; r < RING_DIM; ++r)
                flat[(size_t)i * RING_DIM + r] = mul_mod_shoup(towers[i][r], w[i], wp[i], q[i]);
        auto t1 = Clock::now();
        g_sink += flat[it % flat.size()];
        samples.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    return median_us(samples);
}

// -----------------------------------------------------------------------------
// (3b) BConv hook 的结果回写：dcrtpoly-impl.h:1218-1223
//      纯 scatter，无模运算，但跨 tower 写会打散 cache line
// -----------------------------------------------------------------------------
static double bench_bconv_scatter(int size_p, int iters) {
    std::vector<uint64_t> flat((size_t)size_p * RING_DIM);
    std::vector<std::vector<uint64_t>> towers(size_p, std::vector<uint64_t>(RING_DIM));
    std::mt19937_64 rng(4242);
    for (auto& v : flat) v = rng();

    std::vector<double> samples;
    samples.reserve(iters);
    for (int it = 0; it < iters; ++it) {
        auto t0 = Clock::now();
        for (int j = 0; j < size_p; ++j)
            for (int r = 0; r < RING_DIM; ++r) towers[j][r] = flat[(size_t)j * RING_DIM + r];
        auto t1 = Clock::now();
        g_sink += towers[0][it % RING_DIM];
        samples.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    return median_us(samples);
}

int main() {
    const int ITERS = 2000;
    const size_t LIMB_BYTES = (size_t)RING_DIM * sizeof(uint64_t);  // 32 KiB

    std::printf("host_overhead_bench — 实验 D 的离线分量\n");
    std::printf("RING_DIM=%d  LIMB_Q=%d  LIMB_P=%d  一个 limb = %zu B\n",
                RING_DIM, LIMB_Q, LIMB_P, LIMB_BYTES);
    std::printf("每项取 %d 次的中位数\n\n", ITERS);

    // --- memcpy 扫描：对应 bo.write()/bo.read() ---
    std::printf("[1] memcpy（= bo.write / bo.read 的等价代价）\n");
    std::printf("  %-8s %10s %12s %12s\n", "limbs", "bytes", "us", "GB/s");
    std::vector<Result> mc;
    for (int limbs : {1, 2, 3, 5, 8}) {
        size_t bytes = limbs * LIMB_BYTES;
        double us    = bench_memcpy(bytes, ITERS);
        mc.push_back({"memcpy", bytes, us});
        std::printf("  %-8d %10zu %12.3f %12.1f\n", limbs, bytes, us, bytes / us / 1000.0);
    }

    // --- gather / scatter：对应 BConv hook 的 host 侧准备与回写 ---
    std::printf("\n[2] BConv hook 的 host 侧 gather/scatter（dcrtpoly-impl.h:1181-1223）\n");
    double g_us = bench_bconv_gather(LIMB_Q, ITERS / 4);
    double s_us = bench_bconv_scatter(LIMB_P, ITERS / 4);
    std::printf("  gather  %d limbs + ModMulFastConst : %8.3f us\n", LIMB_Q, g_us);
    std::printf("  scatter %d limbs                   : %8.3f us\n", LIMB_P, s_us);

    // --- 汇总成「每次 offload 调用」的口径 ---
    std::printf("\n[3] 折算到单次 offload 调用（对账 task.yaml 的 182.8 us）\n");
    double ntt_memcpy = 3 * mc[0].us;  // in1 + in2(复制 in1) + out，各 1 limb
    std::printf("  OP_NTT/INTT (1 limb):  3 x memcpy(32K)        = %7.3f us\n", ntt_memcpy);

    double bconv_memcpy = bench_memcpy(LIMB_Q * LIMB_BYTES, ITERS)      // bo_in.write
                          + bench_memcpy(240, ITERS)                    // meta，30 x u64
                          + bench_memcpy(LIMB_P * LIMB_BYTES, ITERS);   // bo_out.read
    std::printf("  OP_BCONV (Q=%d -> P=%d):  memcpy               = %7.3f us\n",
                LIMB_Q, LIMB_P, bconv_memcpy);
    std::printf("                          + gather + scatter    = %7.3f us\n", g_us + s_us);
    std::printf("                          = 合计                = %7.3f us\n",
                bconv_memcpy + g_us + s_us);

    std::printf("\n注：xrt::bo 构造/析构（驱动 ioctl + pin 内存 + mmap）不在此列，必须上板测。\n");
    std::printf("    本程序量的是「未计时开销」里【不需要板卡】的那一半。\n");
    if (g_sink == 0x1234567890ABCDEFULL) std::printf(" ");  // 防优化
    return 0;
}
