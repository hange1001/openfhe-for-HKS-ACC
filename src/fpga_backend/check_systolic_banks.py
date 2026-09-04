#!/usr/bin/env python3
"""P2 bank access checker: verify the direct work-area systolic addressing.

Models the proposed P2 bconv_systolic_core scheduling for every cycle t in
[0, TOTAL_CYCLES):

  Feed_X  (read ): in_x[q][(t-q)>>6][(t-q)&63]          for q < active_q
  Collect (write): in_x[LIMB_Q+p][(t-3-p)>>6][(t-3-p)&63] for p < sizeP

Contract checked per cycle (based on the top-level storage view):
  - in_x is ram_2p BRAM with one physical instance per limb (dim0),
    dim3 partitioned cyclic=BANKS (BConv wrapper view, LOAD_PAR=8).
  - Reads only touch rows 0..active_q-1, writes only rows LIMB_Q..LIMB_Q+sizeP-1,
    so reads and writes never share a physical instance.
  - Within one cycle, all read banks are distinct and all write banks are
    distinct -> at most one R and one W per (row, bank) per cycle.
  - Addresses stay within [0, RING_DIM) whenever an access is guarded valid.

Usage: python3 check_systolic_banks.py [--banks 8] [--active_q 2] [--sizeP 3]
Runs the full active_q x sizeP grid when no specific values are given.
"""
import argparse
import sys

LIMB_Q = 3
MAX_OUT_COLS = 5
RING_DIM = 4096
LOG_SQRT = 6
SQRT = 64
TOTAL_CYCLES = LIMB_Q + RING_DIM + MAX_OUT_COLS - 1  # 4103


def check(active_q, sizeP, banks):
    failures = []
    if not (1 <= active_q <= LIMB_Q):
        failures.append(f"active_q={active_q} out of [1,{LIMB_Q}]")
        return failures
    if not (1 <= sizeP <= MAX_OUT_COLS):
        failures.append(f"sizeP={sizeP} out of [1,{MAX_OUT_COLS}]")
        return failures

    for t in range(TOTAL_CYCLES):
        reads = {}   # (row, bank) -> count
        writes = {}  # (row, bank) -> count
        for q in range(LIMB_Q):
            idx = t - q
            valid = q < active_q and 0 <= idx < RING_DIM
            if valid:
                row, col = idx >> LOG_SQRT, idx & (SQRT - 1)
                assert 0 <= row < SQRT and 0 <= col < SQRT, "read address out of range"
                reads[(q, col % banks)] = reads.get((q, col % banks), 0) + 1
        for p in range(MAX_OUT_COLS):
            idx = t - (LIMB_Q + p)
            if p < sizeP and 0 <= idx < RING_DIM:
                row, col = idx >> LOG_SQRT, idx & (SQRT - 1)
                assert 0 <= row < SQRT and 0 <= col < SQRT, "write address out of range"
                writes[(LIMB_Q + p, col % banks)] = writes.get((LIMB_Q + p, col % banks), 0) + 1
        for key, n in reads.items():
            if n > 1:
                failures.append(f"t={t} read bank conflict {key} x{n}")
        for key, n in writes.items():
            if n > 1:
                failures.append(f"t={t} write bank conflict {key} x{n}")
        # read/write on different limbs -> different physical instances; verify.
        if set(r[0] for r in reads) & set(w[0] for w in writes):
            failures.append(f"t={t} read and write share a limb instance")
        if failures:
            break  # one offending cycle is enough detail
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--banks", type=int, default=8, help="dim3 cyclic partition factor")
    parser.add_argument("--active_q", type=int, help="specific active_q to check")
    parser.add_argument("--sizeP", type=int, help="specific sizeP to check")
    args = parser.parse_args()

    grid = []
    if args.active_q is not None and args.sizeP is not None:
        grid = [(args.active_q, args.sizeP)]
    else:
        grid = [(q, p) for q in range(1, LIMB_Q + 1)
                for p in range(1, MAX_OUT_COLS + 1)]

    all_fail = []
    for active_q, sizeP in grid:
        fails = check(active_q, sizeP, args.banks)
        status = "PASS" if not fails else "FAIL"
        print(f"[{status}] active_q={active_q} sizeP={sizeP} banks={args.banks}")
        for f in fails[:3]:
            print("      " + f)
        all_fail.extend(fails)

    if all_fail:
        print(f"\nRESULT: FAIL ({len(all_fail)} conflict(s))")
        sys.exit(1)
    print(f"\nRESULT: PASS ({len(grid)} combos conflict-free)")


if __name__ == "__main__":
    main()
