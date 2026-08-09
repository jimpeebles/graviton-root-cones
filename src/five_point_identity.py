#!/usr/bin/env python3
"""Exact verification of an explicit non-decay five-graviton chamber formula."""

from __future__ import annotations

from itertools import combinations, product

from amplitude_search import add, amplitude, bracket
from chamber_enumeration import chamber_sign


TARGET_ROOT_CHAMBER = chamber_sign((-3, 2, 2, -1))
TARGET_PAIR_SIGNS = (-1, -1, 1, 1, -1, -1)  # 12,13,14,23,24,34
TARGET_TRIPLE_SIGNS = {
    (0, 1, 2): (-1, 1, 1),
    (0, 1, 3): (1, 1, -1),
    (0, 2, 3): (-1, 1, -1),
    (1, 2, 3): (1, -1, 1),
}


def sign(value: int) -> int:
    if value == 0:
        return 0
    return 1 if value > 0 else -1


def in_target_chamber(q: tuple[tuple[int, int], ...]) -> bool:
    if len(q) != 4:
        raise ValueError("Expected four positive-helicity spinors")
    total = add(*q)
    x = tuple(bracket(vector, total) for vector in q)
    if chamber_sign(x) != TARGET_ROOT_CHAMBER:
        return False

    pairs = tuple(bracket(q[i], q[j]) for i, j in combinations(range(4), 2))
    if tuple(sign(value) for value in pairs) != TARGET_PAIR_SIGNS:
        return False

    for labels, expected in TARGET_TRIPLE_SIGNS.items():
        triple_total = add(*(q[i] for i in labels))
        projected = tuple(bracket(q[i], triple_total) for i in labels)
        if tuple(sign(value) for value in projected) != expected:
            return False
    return True


def compact_formula(q: tuple[tuple[int, int], ...]) -> int:
    q5 = tuple(-value for value in add(*q))
    momenta = q + (q5,)

    def w(i: int, j: int) -> int:
        return abs(bracket(momenta[i - 1], momenta[j - 1]))

    return (
        w(1, 3) * w(1, 4) * (w(1, 2) + w(2, 3))
        - w(1, 4) * w(2, 4) * (w(2, 3) + w(3, 4))
        - w(2, 5) * w(1, 4) * w(3, 4)
        - w(3, 5) * w(1, 2) * w(2, 4)
    )


def exhaustive_check(max_flow: int = 12, max_gap: int = 8) -> int:
    checked = 0
    for a, b, c in product(range(1, max_flow + 1), repeat=3):
        d = b + c - a
        if d <= 0:
            continue
        x = (-a, b, c, -d)
        if chamber_sign(x) != TARGET_ROOT_CHAMBER:
            continue
        for gaps in product(range(1, max_gap + 1), repeat=3):
            slopes = (0, gaps[0], gaps[0] + gaps[1], sum(gaps))
            q = tuple((x_i * slope, -x_i) for x_i, slope in zip(x, slopes))
            if add(*q)[0] <= 0 or not in_target_chamber(q):
                continue
            exact = amplitude(q)
            compact = compact_formula(q)
            if exact != compact:
                raise AssertionError(
                    f"Five-point identity failed: x={x}, gaps={gaps}, "
                    f"exact={exact}, compact={compact}"
                )
            checked += 1
    return checked


if __name__ == "__main__":
    count = exhaustive_check()
    print(f"five-point compact identity: passed {count} exact chamber samples")
