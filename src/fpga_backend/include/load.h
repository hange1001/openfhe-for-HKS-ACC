#ifndef LOAD_H
#define LOAD_H
#include <cstddef>
#include <cstdint>
#include "define.h"

extern "C" {

void Load(
    const uint64_t *mem_in,
    uint64_t buffer[MAX_LIMBS][SQRT][SQRT],
    const int num_active_limbs,
    const int mod_index
);
}

extern "C" {
void Store(
    uint64_t buffer[MAX_LIMBS][SQRT][SQRT],
    uint64_t *mem_out,
    const int num_active_limbs,
    const int mod_index
);

}

#endif // LOAD_H