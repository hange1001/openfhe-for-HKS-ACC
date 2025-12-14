#ifndef LOAD_H
#define LOAD_H
#include <cstddef>
#include <cstdint>
#include "define.h"


void Load(
    const uint64_t *mem_in,
    uint64_t buffer[MAX_LIMBS][SQRT][SQRT],
    const int num_limbs,
    const int limb_offset // 偏移量
);

void Store(
    uint64_t buffer[MAX_LIMBS][SQRT][SQRT],
    uint64_t *mem_out,
    const int num_limbs,
    const int limb_offset // 偏移量
);

#endif // LOAD_H