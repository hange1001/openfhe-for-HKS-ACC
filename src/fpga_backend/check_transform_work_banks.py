#!/usr/bin/env python3
"""Independent access-count proof for direct work/scratch transform phases.

This checks the algorithmic port budget, NOT HLS scheduling or routed hardware.
Combine with check_shared_transform.py and the actual RAM binding reports.
"""
import argparse
from collections import Counter

def check(lanes):
    n, stages, slots, banks = 4096, 12, 8, 2 * lanes
    assert n % (2 * lanes) == 0 and 64 % banks == 0
    for slot in range(slots):
        for forward in (False, True):
            for stage in range(stages):
                reads, writes = set(), set()
                for base in range(0, n // 2, lanes):
                    ports = Counter()
                    for k in range(base, base + lanes):
                        r = (k, k + n // 2) if forward else (2 * k, 2 * k + 1)
                        w = (2 * k, 2 * k + 1) if forward else (k, k + n // 2)
                        for address in r:
                            assert 0 <= address < n and address not in reads
                            reads.add(address)
                            ports[(stage % 2, address % banks)] += 1
                        for address in w:
                            assert 0 <= address < n and address not in writes
                            writes.add(address)
                            ports[(1 - stage % 2, address % banks)] += 1
                    assert max(ports.values()) <= 2
                assert len(reads) == len(writes) == n
    # In-place SCALE has one read and one write per bank; addresses never repeat.
    for base in range(0, n, lanes):
        banks_used = [i % banks for i in range(base, base + lanes)]
        assert len(banks_used) == len(set(banks_used))
    print(f"PASS: {slots} slots x 2 directions x {stages} stages; {banks} T2P banks, <=2 accesses/bank/cycle")
    print("PASS: all 4096 coefficients covered once per stage; SCALE one read + one write/bank")
    print("Scope: algorithmic port budget only; actual HLS II, enables and RAM binding require separate audit")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lanes', type=int, default=4)
    check(parser.parse_args().lanes)
