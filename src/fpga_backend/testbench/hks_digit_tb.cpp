#include "../include/top.h"
#include "../include/cg_ntt.h"
#include "../include/hks_digit.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {
const uint64_t CANARY = 0xBAD0BAD0BAD0BAD0ULL;
const int IN1_DEPTH = 3 * LIMB_Q + MAX_LIMBS * STAGE * CG_HALF_N;
const int IN2_DEPTH = 3 * LIMB_P + MAX_LIMBS * STAGE * CG_HALF_N;
int failures = 0;
int cases = 0;
uint64_t primes[MAX_OUT_COLS], roots[MAX_OUT_COLS];

void expect(bool ok, const std::string& label) {
    if (!ok) {
        ++failures;
        std::cerr << "[HKS FAIL] " << label << '\n';
    }
}

uint64_t mul(uint64_t a, uint64_t b, uint64_t p) {
    return (unsigned __int128)a * b % p;
}
uint64_t power(uint64_t a, uint64_t e, uint64_t p) {
    uint64_t r = 1;
    while (e) {
        if (e & 1) r = mul(r, a, p);
        a = mul(a, a, p);
        e >>= 1;
    }
    return r;
}
void barrett(uint64_t p, uint64_t& s, uint64_t& m) {
    unsigned bits = 0;
    for (uint64_t t = p; t; t >>= 1) ++bits;
    s = bits + 62;
    m = ((unsigned __int128)1 << s) / p;
}
int reverse_bits(int x) {
    int r = 0;
    for (int b = 0; b < STAGE; ++b) {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    return r;
}

// Independent iterative radix-2 cyclic NTT (no CG kernel/MultMod calls).
// The current Top twiddle builder yields bit-reversed EVAL output.
std::vector<uint64_t> eval_reference(std::vector<uint64_t> a, int l) {
    const uint64_t p = primes[l];
    for (int i = 0; i < RING_DIM; ++i) {
        const int j = reverse_bits(i);
        if (i < j) std::swap(a[i], a[j]);
    }
    const uint64_t omega = mul(roots[l], roots[l], p);
    for (int len = 2; len <= RING_DIM; len <<= 1) {
        const uint64_t step = power(omega, RING_DIM / len, p);
        for (int base = 0; base < RING_DIM; base += len) {
            uint64_t w = 1;
            for (int k = 0; k < len / 2; ++k) {
                uint64_t u = a[base + k], v = mul(a[base + k + len / 2], w, p);
                a[base + k] = (u + v) % p;
                a[base + k + len / 2] = (u + p - v) % p;
                w = mul(w, step, p);
            }
        }
    }
    std::vector<uint64_t> native(RING_DIM);
    for (int i = 0; i < RING_DIM; ++i) native[i] = a[reverse_bits(i)];
    return native;
}

// Full-depth buffers also satisfy HLS C/RTL co-simulation pointer contracts.
std::vector<uint64_t> call(uint8_t op, int a, int start,
                          const std::vector<uint64_t>& in,
                          const std::vector<uint64_t>& meta, int nout) {
    std::vector<uint64_t> b1(IN1_DEPTH, CANARY), b2(IN2_DEPTH, CANARY);
    std::vector<uint64_t> out(MAX_LIMBS * RING_DIM, CANARY);
    std::copy(in.begin(), in.end(), b1.begin());
    std::copy(meta.begin(), meta.end(), b2.begin());
    Top(b1.data(), b2.data(), out.data(), op, a, start);
    expect(std::all_of(out.begin() + nout, out.end(),
                      [](uint64_t v) { return v == CANARY; }),
           "output guard / no-write contract, opcode=" + std::to_string(op));
    expect(std::equal(in.begin(), in.end(), b1.begin()), "input must be read-only");
    expect(std::equal(meta.begin(), meta.end(), b2.begin()), "metadata must be read-only");
    out.resize(nout);
    return out;
}

void initialize(bool invalid_modulus = false) {
    std::vector<uint64_t> in1(IN1_DEPTH, 0), in2(IN2_DEPTH, 0);
    for (int l = 0; l < MAX_OUT_COLS; ++l) {
        uint64_t s, m;
        barrett(primes[l], s, m);
        if (l < LIMB_Q) {
            in1[l] = primes[l]; in1[LIMB_Q + l] = s; in1[2 * LIMB_Q + l] = m;
        } else {
            int j = l - LIMB_Q;
            in2[j] = primes[l]; in2[LIMB_P + j] = s; in2[2 * LIMB_P + j] = m;
        }
        std::vector<int> perm(RING_DIM);
        for (int i = 0; i < RING_DIM; ++i) perm[i] = reverse_bits(i);
        uint64_t invroot = power(roots[l], primes[l] - 2, primes[l]);
        for (int stage = 0; stage < STAGE; ++stage) {
            for (int i = 0; i < CG_HALF_N; ++i) {
                const uint64_t exp = (perm[i] % (1 << stage)) * (RING_DIM >> stage);
                int offset = l * STAGE * CG_HALF_N + stage * CG_HALF_N + i;
                in1[3 * LIMB_Q + offset] = power(roots[l], exp, primes[l]);
                in2[3 * LIMB_P + offset] = power(invroot, exp, primes[l]);
            }
            std::vector<int> next(RING_DIM);
            for (int i = 0; i < CG_HALF_N; ++i) {
                next[2 * i] = perm[i]; next[2 * i + 1] = perm[i + CG_HALF_N];
            }
            perm.swap(next);
        }
    }
    if (invalid_modulus) in1[0] = 0;
    call(OP_INIT, 0, 0, in1, in2, 0);
}

bool same(const std::vector<uint64_t>& a, const std::vector<uint64_t>& b,
          const std::string& label) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) {
            std::cerr << label << " limb=" << i / RING_DIM << " coeff=" << i % RING_DIM
                      << " got=" << a[i] << " want=" << b[i] << '\n';
            return false;
        }
    }
    return true;
}

void run_case(int alpha, int start, unsigned seed, int pattern) {
    ++cases;
    const int before = failures;
    const int count = MAX_OUT_COLS - alpha;
    const std::string label = "alpha=" + std::to_string(alpha) + " start=" +
        std::to_string(start) + " seed=" + std::to_string(seed) + " pattern=" + std::to_string(pattern);
    std::mt19937_64 rng(seed);
    std::vector<std::vector<uint64_t>> coeff(alpha, std::vector<uint64_t>(RING_DIM));
    std::vector<uint64_t> input(alpha * RING_DIM), meta(HKS_META_WORDS, CANARY);
    std::vector<uint64_t> expected(MAX_OUT_COLS * RING_DIM, 0);
    std::vector<uint64_t> padded(LIMB_Q * RING_DIM, 0);
    std::vector<uint64_t> bmeta(HKS_WEIGHT_WORDS + 3 * MAX_OUT_COLS, 0);

    for (int q = 0; q < alpha; ++q) {
        const int l = start + q;
        const uint64_t qi = primes[l];
        for (int n = 0; n < RING_DIM; ++n) {
            coeff[q][n] = pattern == 1 ? 0 : pattern == 2 ? qi - 1 : rng() % qi;
        }
        std::vector<uint64_t> eval = eval_reference(coeff[q], l);
        std::copy(eval.begin(), eval.end(), input.begin() + q * RING_DIM);
        std::copy(eval.begin(), eval.end(), expected.begin() + l * RING_DIM);

        uint64_t qhat_mod_qi = 1;
        for (int k = 0; k < alpha; ++k)
            if (k != q) qhat_mod_qi = mul(qhat_mod_qi, primes[start + k] % qi, qi);
        const uint64_t inv = power(qhat_mod_qi, qi - 2, qi);
        meta[HKS_INV_OFFSET + q] = inv;
        expect(mul(qhat_mod_qi, inv, qi) == 1, "CRT inverse sanity");

        std::vector<uint64_t> recovered = call(OP_INTT, 1, l, eval, {}, RING_DIM);
        expect(same(recovered, coeff[q], label + " INTT"), "native INTT vs independent input");
        for (int n = 0; n < RING_DIM; ++n) padded[q * RING_DIM + n] = mul(recovered[n], inv, qi);

        for (int p = 0; p < count; ++p) {
            const int target = p < start ? p : p + alpha;
            uint64_t w = 1;
            for (int k = 0; k < alpha; ++k)
                if (k != q) w = mul(w, primes[start + k] % primes[target], primes[target]);
            meta[q * MAX_OUT_COLS + p] = w;
            bmeta[q * MAX_OUT_COLS + p] = w;
        }
    }
    for (int p = 0; p < count; ++p) {
        const int l = p < start ? p : p + alpha;
        uint64_t s, m;
        barrett(primes[l], s, m);
        bmeta[HKS_WEIGHT_WORDS + p] = primes[l];
        bmeta[HKS_WEIGHT_WORDS + MAX_OUT_COLS + p] = s;
        bmeta[HKS_WEIGHT_WORDS + 2 * MAX_OUT_COLS + p] = m;
        std::vector<uint64_t> converted(RING_DIM, 0);
        for (int n = 0; n < RING_DIM; ++n) {
            for (int q = 0; q < alpha; ++q) {
                uint64_t x = mul(coeff[q][n], meta[HKS_INV_OFFSET + q], primes[start + q]);
                converted[n] = (converted[n] + mul(x, meta[q * MAX_OUT_COLS + p], primes[l])) % primes[l];
            }
        }
        std::vector<uint64_t> e = eval_reference(converted, l);
        std::copy(e.begin(), e.end(), expected.begin() + l * RING_DIM);
    }

    std::vector<uint64_t> bc = call(OP_BCONV, count, 0, padded, bmeta, count * RING_DIM);
    std::vector<uint64_t> separated(MAX_OUT_COLS * RING_DIM);
    for (int q = 0; q < alpha; ++q)
        std::copy(input.begin() + q * RING_DIM, input.begin() + (q + 1) * RING_DIM,
                  separated.begin() + (start + q) * RING_DIM);
    for (int p = 0; p < count; ++p) {
        const int l = p < start ? p : p + alpha;
        std::vector<uint64_t> tower(bc.begin() + p * RING_DIM, bc.begin() + (p + 1) * RING_DIM);
        std::vector<uint64_t> e = call(OP_NTT, 1, l, tower, {}, RING_DIM);
        std::copy(e.begin(), e.end(), separated.begin() + l * RING_DIM);
    }
    // Inactive metadata is deliberately poisoned. Prior calls leave scratch dirty.
    std::vector<uint64_t> fused = call(OP_HKS_DIGIT, alpha, start, input, meta, MAX_OUT_COLS * RING_DIM);
    expect(same(fused, expected, label + " golden"), label + " scalar CRT/NTT reference");
    expect(same(fused, separated, label + " split"), label + " separate opcode reference");
    std::cout << "[HKS " << (before == failures ? "PASS" : "FAIL") << "] " << label << '\n';
}

bool is_prime(uint64_t n) {
    if (n < 2) return false;
    for (uint64_t p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL}) {
        if (n % p == 0) return n == p;
    }
    uint64_t d = n - 1;
    int s = 0;
    while (!(d & 1)) { d >>= 1; ++s; }
    for (uint64_t a : {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL}) {
        if (a % n == 0) continue;
        uint64_t x = power(a % n, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int r = 1; r < s; ++r) {
            x = mul(x, x, n);
            if (x == n - 1) { composite = false; break; }
        }
        if (composite) return false;
    }
    return true;
}
bool set_roots() {
    for (int l = 0; l < MAX_OUT_COLS; ++l) {
        expect(is_prime(primes[l]), "prime sanity");
        bool found = false;
        for (uint64_t g = 2; g < 1000; ++g) {
            uint64_t r = power(g, (primes[l] - 1) / (2 * RING_DIM), primes[l]);
            if (power(r, RING_DIM, primes[l]) == primes[l] - 1) {
                roots[l] = r; found = true; break;
            }
        }
        expect(found, "primitive 2N-th root sanity");
        if (!found) return false;
    }
    return true;
}
} // namespace

int hks_digit_uninitialized_test() {
    const int before = failures;
    call(OP_HKS_DIGIT, 2, 0, {}, {}, 0);
    return failures - before;
}

int run_hks_digit_tests() {
    const int before = failures;
    const uint64_t small[] = {998244353ULL, 1004535809ULL, 469762049ULL, 167772161ULL, 1224736769ULL};
    std::copy(small, small + MAX_OUT_COLS, primes);
    if (!set_roots()) return failures - before;
    initialize();
    for (int a = 1; a <= LIMB_Q; ++a)
        for (int s = 0; s <= LIMB_Q - a; ++s) run_case(a, s, 0x7101U + a + s, 0);
    for (unsigned seed : {0x8102U, 0x9103U, 0xA104U}) {
        run_case(2, 0, seed, 0);
        run_case(1, 2, seed, 0);
    }
    for (int pattern : {1, 2}) {
        run_case(2, 0, 1, pattern);
        run_case(1, 2, 2, pattern);
    }
    const int invalid[][2] = {{0, 0}, {-1, 0}, {4, 0}, {2, -1}, {2, 2},
                              {1, 3}, {1, 2147483647}, {2147483647, 0}};
    for (const auto& d : invalid) call(OP_HKS_DIGIT, d[0], d[1], {}, {}, 0);
    initialize(true);
    call(OP_HKS_DIGIT, 1, 2, {}, {}, 0);

    // Distinct 60-bit NTT primes exercise full-size residue arithmetic as well.
    uint64_t candidate = (uint64_t(1) << 60) - 2 * RING_DIM + 1;
    for (int l = 0; l < MAX_OUT_COLS; ++l) {
        while (!is_prime(candidate)) candidate -= 2 * RING_DIM;
        primes[l] = candidate;
        candidate -= 2 * RING_DIM;
    }
    if (!set_roots()) return failures - before;
    initialize();
    for (unsigned seed : {0xB105U, 0xC106U}) {
        run_case(2, 0, seed, 0);
        run_case(1, 2, seed, 0);
    }
    run_case(2, 0, 1, 2);
    run_case(1, 2, 2, 2);
    std::cout << "[HKS SUMMARY] " << cases << " valid cases, " << failures - before
              << " failures; invalid descriptors and canaries checked\n";
    return failures - before;
}
