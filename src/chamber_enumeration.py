#!/usr/bin/env python3
"""Exhaustive low-rank chamber classification for the graviton vertex.

Use ``--json PATH`` to write machine-readable representatives, feasible-tree
families, and arborescence roots for every realizable chamber.
"""

from __future__ import annotations

import argparse
from itertools import combinations, product
import json
from math import floor
from pathlib import Path

from amplitude_search import add, bracket, spanning_trees, vertex


Edge = tuple[int, int]
Vec = tuple[int, int]


KNOWN_RESONANCE_CHAMBERS = {4: 32, 5: 370}


def chamber_subsets(m: int) -> tuple[tuple[int, ...], ...]:
    return tuple(
        subset
        for size in range(1, m)
        for subset in combinations(range(m - 1), size)
    )


def chamber_sign(x: tuple[int, ...]) -> tuple[bool, ...] | None:
    subsets = chamber_subsets(len(x))
    values = tuple(sum(x[i] for i in subset) for subset in subsets)
    if any(value == 0 for value in values):
        return None
    return tuple(value > 0 for value in values)


def enumerate_resonance_chambers(m: int) -> dict[tuple[bool, ...], tuple[int, ...]]:
    """Grid enumeration with a completeness check against known chamber counts."""
    target = KNOWN_RESONANCE_CHAMBERS[m]
    for bound in range(1, 10):
        chambers: dict[tuple[bool, ...], tuple[int, ...]] = {}
        for coordinates in product(range(-bound, bound + 1), repeat=m - 1):
            x = coordinates + (-sum(coordinates),)
            sign = chamber_sign(x)
            if sign is not None:
                chambers.setdefault(sign, x)
        if len(chambers) == target:
            return chambers
    raise RuntimeError(f"Failed to enumerate all {target} chambers for m={m}")


def realizable_spinors(x: tuple[int, ...]) -> tuple[Vec, ...] | None:
    """Realize a chamber with increasing slopes and Q on the positive x-axis.

    For q_i=(x_i r_i,-x_i), Q=(S,0), where
    S=sum_i x_i r_i=-sum_k prefix_k(x)*(r_{k+1}-r_k).
    Positive gaps can give S>0 iff at least one proper prefix sum is negative.
    """
    prefixes = [sum(x[:k]) for k in range(1, len(x))]
    negative = [k for k, value in enumerate(prefixes) if value < 0]
    if not negative:
        return None

    gaps = [1] * (len(x) - 1)
    s_value = -sum(prefixes)
    if s_value <= 0:
        index = negative[0]
        gain = -prefixes[index]
        gaps[index] += floor((-s_value) / gain) + 1

    slopes = [0]
    for gap in gaps:
        slopes.append(slopes[-1] + gap)
    q = tuple((x_i * slope, -x_i) for x_i, slope in zip(x, slopes))
    q_total = add(*q)
    if q_total[1] != 0 or q_total[0] <= 0:
        raise AssertionError((x, prefixes, gaps, q_total))
    projected = tuple(bracket(vector, q_total) for vector in q)
    if chamber_sign(projected) != chamber_sign(x):
        raise AssertionError("Spinor realization changed resonance chamber")
    return q


def oriented_edge(q: tuple[Vec, ...], edge: Edge, advanced: bool = False) -> Edge:
    a, b = edge
    if bracket(q[a], q[b]) < 0:
        direction = (a, b)
    else:
        direction = (b, a)
    return direction[::-1] if advanced else direction


def feasible_tree_family(q: tuple[Vec, ...], advanced: bool = False) -> frozenset[tuple[Edge, ...]]:
    m = len(q)
    total = add(*q)
    netflow = tuple(bracket(vector, total) for vector in q)
    answer = []
    for tree in spanning_trees(m):
        adjacency = [set() for _ in range(m)]
        for a, b in tree:
            adjacency[a].add(b)
            adjacency[b].add(a)
        parent = {0: None}
        order = [0]
        for current in order:
            for neighbor in adjacency[current]:
                if neighbor not in parent:
                    parent[neighbor] = current
                    order.append(neighbor)
        subtree = list(netflow)
        feasible = True
        for child in reversed(order[1:]):
            parent_node = parent[child]
            assert parent_node is not None
            tail, head = oriented_edge(q, tuple(sorted((child, parent_node))), advanced)
            points_out = tail == child
            if (subtree[child] > 0) != points_out:
                feasible = False
                break
            subtree[parent_node] += subtree[child]
        if feasible:
            answer.append(tree)
    return frozenset(answer)


def is_graphic_family(family: frozenset[tuple[Edge, ...]], m: int) -> bool:
    if not family:
        return True
    edge_union = frozenset(edge for tree in family for edge in tree)
    graph_trees = frozenset(
        tree for tree in spanning_trees(m) if all(edge in edge_union for edge in tree)
    )
    return graph_trees == family


def points_to_root(tree: tuple[Edge, ...], q: tuple[Vec, ...], root: int, advanced: bool) -> bool:
    adjacency = [set() for _ in q]
    for a, b in tree:
        adjacency[a].add(b)
        adjacency[b].add(a)
    parent = {root: None}
    order = [root]
    for current in order:
        for neighbor in adjacency[current]:
            if neighbor not in parent:
                parent[neighbor] = current
                order.append(neighbor)
    return all(
        oriented_edge(q, tuple(sorted((node, parent[node]))), advanced) == (node, parent[node])
        for node in range(len(q))
        if node != root
    )


def arborescence_roots(
    family: frozenset[tuple[Edge, ...]], q: tuple[Vec, ...], advanced: bool = False
) -> tuple[int, ...]:
    """Roots for which the family is every in-arborescence of its edge union."""
    if not family:
        return ()
    edge_union = frozenset(edge for tree in family for edge in tree)
    roots = []
    for root in range(len(q)):
        generated = frozenset(
            tree
            for tree in spanning_trees(len(q))
            if all(edge in edge_union for edge in tree)
            and points_to_root(tree, q, root, advanced)
        )
        if generated == family:
            roots.append(root)
    return tuple(roots)


def bareiss_determinant(matrix: list[list[int]]) -> int:
    """Fraction-free exact determinant."""
    n = len(matrix)
    if n == 0:
        return 1
    values = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(n - 1):
        if values[pivot_index][pivot_index] == 0:
            swap = next(
                (row for row in range(pivot_index + 1, n) if values[row][pivot_index]),
                None,
            )
            if swap is None:
                return 0
            values[pivot_index], values[swap] = values[swap], values[pivot_index]
            sign *= -1
        pivot = values[pivot_index][pivot_index]
        for row in range(pivot_index + 1, n):
            for column in range(pivot_index + 1, n):
                values[row][column] = (
                    values[row][column] * pivot
                    - values[row][pivot_index] * values[pivot_index][column]
                ) // previous
        previous = pivot
    return sign * values[-1][-1]


def arborescence_determinant(
    family: frozenset[tuple[Edge, ...]], q: tuple[Vec, ...], root: int, advanced: bool = False
) -> int:
    """Directed matrix-tree cofactor for the union graph of a family."""
    edge_union = frozenset(edge for tree in family for edge in tree)
    m = len(q)
    laplacian = [[0 for _ in range(m)] for _ in range(m)]
    for edge in edge_union:
        tail, head = oriented_edge(q, edge, advanced)
        weight = abs(bracket(q[tail], q[head]))
        laplacian[tail][tail] += weight
        laplacian[tail][head] -= weight
    minor = [
        [laplacian[row][column] for column in range(m) if column != root]
        for row in range(m)
        if row != root
    ]
    return bareiss_determinant(minor)


def classify(m: int) -> dict[str, object]:
    chambers = enumerate_resonance_chambers(m)
    realizable = []
    family_histogram: dict[frozenset[tuple[Edge, ...]], int] = {}
    rows = []
    determinant_checks = 0
    for sign, x in chambers.items():
        q = realizable_spinors(x)
        if q is None:
            continue
        family = feasible_tree_family(q)
        family_histogram[family] = family_histogram.get(family, 0) + 1
        graphic = is_graphic_family(family, m)
        roots = arborescence_roots(family, q)
        if roots:
            expected = vertex(q)
            for root in roots:
                determinant = arborescence_determinant(family, q, root)
                if determinant != expected:
                    raise AssertionError(
                        f"Matrix-tree mismatch m={m}, x={x}, root={root}: "
                        f"det={determinant}, V={expected}"
                    )
                determinant_checks += 1
        rows.append((sign, x, q, len(family), graphic, roots))
        realizable.append(sign)

    nonempty = [row for row in rows if row[3] > 0]
    nontrivial = [row for row in nonempty if row[3] > 1]
    result = {
        "m": m,
        "all_chambers": len(chambers),
        "realizable_chambers": len(realizable),
        "distinct_families": len(family_histogram),
        "empty_vertices": sum(row[3] == 0 for row in rows),
        "graphic_nonempty": sum(row[4] for row in nonempty),
        "graphic_nontrivial": sum(row[4] for row in nontrivial),
        "arborescent_nonempty": sum(bool(row[5]) for row in nonempty),
        "arborescent_nontrivial": sum(bool(row[5]) for row in nontrivial),
        "max_family_size": max((row[3] for row in rows), default=0),
        "family_sizes": sorted({row[3] for row in rows}),
        "determinant_checks": determinant_checks,
        "examples_graphic_nontrivial": [row for row in rows if row[3] > 1 and row[4]][:5],
        "examples_arborescent_nontrivial": [row for row in rows if row[3] > 1 and row[5]][:5],
        "_certificates": [
            {
                "resonance_signs": [int(value) for value in sign],
                "netflow_representative": list(x),
                "spinors": [list(vector) for vector in q],
                "feasible_trees": [
                    [[a + 1, b + 1] for a, b in tree]
                    for tree in sorted(feasible_tree_family(q))
                ],
                "arborescence_roots": [root + 1 for root in roots],
            }
            for sign, x, q, _, _, roots in rows
        ],
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=Path, help="write complete chamber certificates")
    arguments = parser.parse_args()
    certificates = {
        "schema": "graviton-root-cone-chambers-v1",
        "python_minimum": "3.11",
        "bracket_convention": "[a,b] = a_0 b_1 - a_1 b_0",
        "valencies": {},
    }
    for m in (4, 5):
        result = classify(m)
        certificates["valencies"][str(m)] = result.pop("_certificates")
        print(f"m={m}")
        for key, value in result.items():
            if not key.startswith("examples_"):
                print(f"  {key}: {value}")
        for key in ("examples_graphic_nontrivial", "examples_arborescent_nontrivial"):
            print(f"  {key}:")
            for _, x, q, size, graphic, roots in result[key]:
                print(f"    x={x}, q={q}, trees={size}, graphic={graphic}, roots={roots}")
    if arguments.json:
        arguments.json.parent.mkdir(parents=True, exist_ok=True)
        arguments.json.write_text(json.dumps(certificates, indent=2) + "\n", encoding="utf-8")
        print(f"wrote chamber certificates: {arguments.json}")


if __name__ == "__main__":
    main()
