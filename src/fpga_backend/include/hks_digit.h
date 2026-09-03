#ifndef HKS_DIGIT_H
#define HKS_DIGIT_H

#include "define.h"

// OP_HKS_DIGIT uses Top's existing scalar arguments: alpha, digit_start.
// in1: compact native-EVAL digit; out: full Q/P basis in global order.
// in2: digit-local weights [LIMB_Q][MAX_OUT_COLS], then QHatInv[LIMB_Q].
// Active complement column p maps to p < digit_start ? p : p + alpha.
// See ../HKS_DIGIT.md for preconditions and the board-free validation scope.
static const int HKS_WEIGHT_WORDS = LIMB_Q * MAX_OUT_COLS;
static const int HKS_INV_OFFSET = HKS_WEIGHT_WORDS;
static const int HKS_META_WORDS = HKS_INV_OFFSET + LIMB_Q;

#endif
