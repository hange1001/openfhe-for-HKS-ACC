#include "../include/shoup.h"
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

int main() {
    std::mt19937_64 rng(0x53484F555032ULL);
    size_t checked = 0, above_target = 0, corrected_cases = 0;
    auto check = [&](uint64_t x, uint64_t w, uint64_t p) {
        const uint64_t pre = uint64_t((static_cast<unsigned __int128>(w) << 64) / p);
        const unsigned __int128 product = static_cast<unsigned __int128>(x) * w;
        const uint64_t expected = uint64_t(product % p);
        const uint64_t qhat = uint64_t((static_cast<unsigned __int128>(x) * pre) >> 64);
        const unsigned __int128 residual = product - static_cast<unsigned __int128>(qhat) * p;
        if (residual >= static_cast<unsigned __int128>(2) * p) return false;
        uint64_t actual = ~expected;
        ShoupMul(x, w, pre, p, actual);
        ++checked;
        above_target += x >= p;
        corrected_cases += residual >= p;
        if (actual != expected || actual >= p) {
            std::cerr << "FAIL x=" << x << " w=" << w << " p=" << p
                      << " got=" << actual << " expected=" << expected << '\n';
            return false;
        }
        return true;
    };
    const uint64_t maximum = std::numeric_limits<uint64_t>::max();
    const uint64_t moduli[] = {2, 3, 17, 167772161, 998244353,
        (uint64_t(1) << 60) - 16383, (uint64_t(1) << 62) - 1, (uint64_t(1) << 63) - 1};
    for (uint64_t p : moduli) {
        for (uint64_t w : {uint64_t(0), uint64_t(1), p/2, p-1})
            for (uint64_t x : {uint64_t(0), uint64_t(1), p-1, p, p+1, maximum-1, maximum})
                if (!check(x, w, p)) return 1;
        for (int i = 0; i < 2048; ++i)
            if (!check(rng(), rng() % p, p)) return 1;
    }
    // Continuously change the modulus and reciprocal to catch stale parameter state.
    for (int i = 0; i < 4096; ++i) {
        const uint64_t p = (rng() >> 1) | 3;
        if (!check(rng(), rng() % p, p)) return 1;
    }
    uint64_t disabled = 123;
    const uint64_t zero = 0;
    ShoupMul(maximum, maximum, maximum, zero, disabled);
    if (disabled != 0 || !above_target || !corrected_cases) return 1;
    std::cout << "PASS: " << checked << " exact residues, " << above_target
              << " inputs >= target modulus, " << corrected_cases
              << " one-subtraction cases; disabled lane returns zero\n";
    return 0;
}
