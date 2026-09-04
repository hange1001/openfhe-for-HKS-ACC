#ifndef FPGA_SHOUP_H
#define FPGA_SHOUP_H
#include <cstdint>

// Constant-multiplier Shoup reduction in the ordinary residue domain.
// Preconditions: 0 <= x < 2^64, 0 <= w < p < 2^63,
// w_pre = floor(w * 2^64 / p). x need NOT be smaller than p.
// p==0 denotes a disabled BConv lane and returns zero.
// Returns the canonical residue [0,p); no Montgomery conversion is involved.
extern "C" void ShoupMul(const uint64_t& x, const uint64_t& w,
                         const uint64_t& w_pre, const uint64_t& p,
                         uint64_t& result);
#endif
