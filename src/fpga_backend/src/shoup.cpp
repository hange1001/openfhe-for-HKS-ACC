#include "../include/shoup.h"
#include <ap_int.h>

void ShoupMul(const uint64_t& x, const uint64_t& w,
              const uint64_t& w_pre, const uint64_t& p, uint64_t& result) {
    #pragma HLS INLINE off
    #pragma HLS PIPELINE II=1
    // One high product and two low products: do not widen the latter to full
    // 128-bit products, and do not introduce a runtime shift or division.
    ap_uint<128> estimate_product = ap_uint<128>(x) * ap_uint<128>(w_pre);
    #pragma HLS BIND_OP variable=estimate_product op=mul impl=dsp latency=4
    const ap_uint<64> qhat = estimate_product.range(127, 64);
    ap_uint<64> product_low = ap_uint<64>(x) * ap_uint<64>(w);
    #pragma HLS BIND_OP variable=product_low op=mul impl=dsp latency=4
    ap_uint<64> multiple_low = qhat * ap_uint<64>(p);
    #pragma HLS BIND_OP variable=multiple_low op=mul impl=dsp latency=4
    ap_uint<64> remainder = product_low - multiple_low;
    #pragma HLS BIND_OP variable=remainder op=sub impl=fabric latency=2
    ap_uint<64> corrected = remainder - ap_uint<64>(p);
    #pragma HLS BIND_OP variable=corrected op=sub impl=fabric latency=2
    // Exact full residual is in [0,2p), which fits a word. Consequently low-word
    // subtraction gives the exact residual, and ONE conditional subtraction
    // suffices even when the input source residue x exceeds the target p.
    result = p == 0 ? uint64_t(0) : uint64_t(remainder >= p ? corrected : remainder);
}
