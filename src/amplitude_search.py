#!/usr/bin/env python3
"""Exact low-point checks for half-collinear single-minus graviton amplitudes.

This implements the vertices and recursion of Eqs. (28)--(33) of
arXiv:2603.04330 using integer two-spinors.  It anchors the implementation
against the paper's three-, four-, and decay-chamber formulas before testing
the root-cone reformulation and two stronger conjectures that fail at five
points.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import combinations_with_replacement, permutations, product
from fractions import Fraction
import random
import sys


Vec = tuple[int, int]
Edge = tuple[int, int]


if sys.version_info < (3, 11):
    raise RuntimeError("These verification programs require Python 3.11 or newer")


def add(*vectors: Vec) -> Vec:
    return sum(v[0] for v in vectors), sum(v[1] for v in vectors)


def bracket(a: Vec, b: Vec) -> int:
    return a[0] * b[1] - a[1] * b[0]


def set_partitions(items: tuple[int, ...], min_blocks: int = 1):
    """Yield canonical unordered set partitions as tuples of tuples."""
    if not items:
        return

    blocks: list[list[int]] = []

    def visit(index: int):
        if index == len(items):
            if len(blocks) >= min_blocks:
                yield tuple(tuple(block) for block in blocks)
            return

        item = items[index]
        for block in blocks:
            block.append(item)
            yield from visit(index + 1)
            block.pop()

        blocks.append([item])
        yield from visit(index + 1)
        blocks.pop()

    yield from visit(0)


@lru_cache(maxsize=None)
def spanning_trees(n: int) -> tuple[tuple[Edge, ...], ...]:
    """All labeled spanning trees on range(n), via Prüfer sequences."""
    if n == 1:
        return ((),)
    if n == 2:
        return (((0, 1),),)

    trees = []
    for sequence in product(range(n), repeat=n - 2):
        degree = [1] * n
        for vertex in sequence:
            degree[vertex] += 1
        edges = []
        for vertex in sequence:
            leaf = next(i for i, d in enumerate(degree) if d == 1)
            edges.append((min(leaf, vertex), max(leaf, vertex)))
            degree[leaf] -= 1
            degree[vertex] -= 1
        last = [i for i, d in enumerate(degree) if d == 1]
        edges.append((min(last), max(last)))
        trees.append(tuple(sorted(edges)))
    return tuple(trees)


def component_after_cut(n: int, edges: tuple[Edge, ...], cut: Edge, start: int) -> set[int]:
    adjacency = [set() for _ in range(n)]
    for edge in edges:
        if edge == cut:
            continue
        a, b = edge
        adjacency[a].add(b)
        adjacency[b].add(a)
    seen = {start}
    stack = [start]
    while stack:
        current = stack.pop()
        for neighbor in adjacency[current] - seen:
            seen.add(neighbor)
            stack.append(neighbor)
    return seen


def vertex(vectors: tuple[Vec, ...], advanced: bool = False) -> int:
    """The retarded V or advanced V-bar vertex away from chamber walls."""
    n = len(vectors)
    if n == 1:
        return 1
    total = 0
    for tree in spanning_trees(n):
        term = 1
        for edge in tree:
            u, v = edge
            side_a = component_after_cut(n, tree, edge, u)
            side_b = set(range(n)) - side_a
            momentum_a = add(*(vectors[i] for i in side_a))
            momentum_b = add(*(vectors[i] for i in side_b))
            edge_bracket = bracket(vectors[u], vectors[v])
            cut_bracket = bracket(momentum_a, momentum_b)
            if edge_bracket == 0 or cut_bracket == 0:
                raise ValueError("Sample lies on a chamber wall")
            theta_argument_sign = (1 if advanced else -1) * cut_bracket * edge_bracket
            if theta_argument_sign <= 0:
                term = 0
                break
            term *= abs(edge_bracket)
        total += term
    return total


def flow_cone_vertex(vectors: tuple[Vec, ...], advanced: bool = False) -> int:
    """Equivalent weighted root-cone enumerator for V or V-bar.

    Set x_i=[q_i,Q], Q=sum(q). Each pair is oriented i->j when [ij]<0
    for the retarded vertex (and oppositely for the advanced vertex). A tree
    contributes precisely when x is the divergence of a strictly positive flow
    on that oriented tree.
    """
    n = len(vectors)
    if n == 1:
        return 1
    total_momentum = add(*vectors)
    netflow = tuple(bracket(vector, total_momentum) for vector in vectors)
    if sum(netflow) != 0:
        raise AssertionError("Projected netflow must sum to zero")

    answer = 0
    for tree in spanning_trees(n):
        adjacency = [set() for _ in range(n)]
        for a, b in tree:
            adjacency[a].add(b)
            adjacency[b].add(a)

        parent = {0: None}
        order = [0]
        for current in order:
            for neighbor in adjacency[current]:
                if neighbor in parent:
                    continue
                parent[neighbor] = current
                order.append(neighbor)

        subtree_sum = list(netflow)
        feasible = True
        for child in reversed(order[1:]):
            parent_node = parent[child]
            assert parent_node is not None
            a, b = sorted((child, parent_node))
            pair_bracket = bracket(vectors[a], vectors[b])
            if pair_bracket == 0 or subtree_sum[child] == 0:
                raise ValueError("Sample lies on a chamber wall")

            # Retarded orientation is a->b iff [ab]<0; advanced reverses it.
            a_to_b = pair_bracket < 0
            if advanced:
                a_to_b = not a_to_b
            tail, head = (a, b) if a_to_b else (b, a)

            # Positive divergence is outward. The subtree divergence is positive
            # exactly when its unique boundary edge points child-side -> parent-side.
            edge_points_out = tail == child or (tail != parent_node and head == parent_node)
            expected_positive = edge_points_out
            if (subtree_sum[child] > 0) != expected_positive:
                feasible = False
                break
            subtree_sum[parent_node] += subtree_sum[child]

        if feasible:
            weight = 1
            for a, b in tree:
                weight *= abs(bracket(vectors[a], vectors[b]))
            answer += weight
    return answer


def verify_flow_cone_equivalence(samples: int = 100, seed: int = 314159) -> None:
    rng = random.Random(seed)
    for n in range(2, 7):
        checked = 0
        while checked < samples:
            vectors = tuple(
                (rng.randint(-9, 9) or 1, rng.randint(-9, 9) or -1)
                for _ in range(n)
            )
            try:
                direct = vertex(vectors)
                direct_bar = vertex(vectors, advanced=True)
                flow = flow_cone_vertex(vectors)
                flow_bar = flow_cone_vertex(vectors, advanced=True)
            except ValueError:
                continue
            if (direct, direct_bar) != (flow, flow_bar):
                raise AssertionError(
                    f"flow-cone mismatch n={n}: direct={(direct, direct_bar)}, "
                    f"flow={(flow, flow_bar)}, vectors={vectors}"
                )
            checked += 1
        print(f"flow-cone equivalence n={n}: passed {checked}/{samples}")


def verify_permutation_invariance(samples: int = 20, seed: int = 271828) -> None:
    rng = random.Random(seed)
    for n in range(3, 6):
        for _ in range(samples):
            momenta = random_sample(n, rng)
            expected = amplitude(momenta[:-1])
            tested = list(permutations(range(n)))
            rng.shuffle(tested)
            for permutation in tested[: min(30, len(tested))]:
                permuted = tuple(momenta[i] for i in permutation)
                if amplitude(permuted[:-1]) != expected:
                    raise AssertionError(
                        f"permutation mismatch n={n}, permutation={permutation}, "
                        f"momenta={momenta}"
                    )
        print(f"permutation invariance n={n}: passed {samples} samples")


def verify_published_low_points(samples: int = 50, seed: int = 161803) -> None:
    """Directly check the paper's displayed three- and four-point formulas."""
    rng = random.Random(seed)
    for _ in range(samples):
        momenta = random_sample(3, rng)
        expected = abs(bracket(momenta[0], momenta[1]))
        if amplitude(momenta[:-1]) != expected:
            raise AssertionError((momenta, amplitude(momenta[:-1]), expected))

    for _ in range(samples):
        momenta = random_sample(4, rng)
        pairings = ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2))
        expected = Fraction(
            sum(
                abs(bracket(momenta[i], momenta[j]))
                * abs(bracket(momenta[k], momenta[l]))
                for i, j, k, l in pairings
            ),
            2,
        )
        if amplitude(momenta[:-1]) != expected:
            raise AssertionError((momenta, amplitude(momenta[:-1]), expected))
    print(f"published M3 and M4 formulas: passed {2 * samples} exact samples")


def decay_sample(n: int, rng: random.Random) -> tuple[tuple[Vec, ...], tuple[Vec, ...]]:
    """Generate the ordered decay chamber used in arXiv:2603.04330.

    Positive-helicity legs have omega>0 and ordered projective coordinates
    z_1<...<z_{n-2}<z_n<z_{n-1}; the final, negative-helicity leg is fixed by
    momentum conservation and has omega_n<0.
    """
    while True:
        energies = [rng.randint(1, 6) for _ in range(n - 2)] + [rng.randint(20, 40)]
        ordered = sorted(rng.sample(range(-30, 30), n - 1))
        slopes = ordered[: n - 2] + [ordered[-1]]
        positive = tuple(
            (energy, energy * slope) for energy, slope in zip(energies, slopes)
        )
        minus = tuple(-component for component in add(*positive))
        minus_slope = Fraction(minus[1], minus[0])
        if not (slopes[n - 3] < minus_slope < slopes[-1]):
            continue
        all_momenta = positive + (minus,)
        if any(
            bracket(a, b) == 0
            for i, a in enumerate(all_momenta)
            for b in all_momenta[i + 1 :]
        ):
            continue
        try:
            amplitude(positive)
        except ValueError:
            continue
        return positive, all_momenta


def published_decay_product(all_momenta: tuple[Vec, ...]) -> int:
    """The product in Eqs. (26)/(39), transcribed in our bracket convention."""
    n = len(all_momenta)
    answer = 1
    for i in range(n - 2):
        answer *= sum(bracket(all_momenta[j], all_momenta[i]) for j in range(i + 1, n - 1))
    return answer


def verify_decay_chamber(samples: int = 50, seed: int = 141421) -> None:
    """Anchor the recursion against the paper's closed decay formula at n=4,5.

    Our determinant convention for square brackets is globally opposite to the
    spinor-component convention used in the paper.  Since M_n has bracket
    degree n-2, direct transcription differs by (-1)^(n-2)=(-1)^n.
    """
    rng = random.Random(seed)
    for n in (4, 5):
        checked = 0
        while checked < samples:
            positive, all_momenta = decay_sample(n, rng)
            expected = ((-1) ** n) * published_decay_product(all_momenta)
            actual = amplitude(positive)
            if actual != expected:
                raise AssertionError(
                    f"decay formula mismatch n={n}: actual={actual}, expected={expected}, "
                    f"momenta={all_momenta}"
                )
            checked += 1
        print(f"published decay formula n={n}: passed {checked}/{samples} exact samples")


def amplitude(positive: tuple[Vec, ...]) -> int:
    """Stripped n-point amplitude; the minus leg is minus sum(positive)."""
    @lru_cache(maxsize=None)
    def pre(labels: tuple[int, ...]) -> int:
        if len(labels) == 1:
            return 1
        if len(labels) == 2:
            return 0
        answer = 0
        for partition in set_partitions(labels, min_blocks=3):
            block_vectors = tuple(add(*(positive[i] for i in block)) for block in partition)
            factor = 1
            for block in partition:
                factor *= pre(tuple(block))
            answer -= vertex(block_vectors) * factor
        return answer

    labels = tuple(range(len(positive)))
    answer = 0
    for partition in set_partitions(labels, min_blocks=2):
        block_vectors = tuple(add(*(positive[i] for i in block)) for block in partition)
        kernel = vertex(block_vectors) - vertex(block_vectors, advanced=True)
        factor = 1
        for block in partition:
            factor *= pre(tuple(block))
        answer -= kernel * factor
    return answer


@lru_cache(maxsize=None)
def cubic_trees(n_leaves: int) -> tuple[tuple[Edge, ...], ...]:
    """All unrooted binary trees with labeled leaves 0,...,n_leaves-1.

    Internal labels are n_leaves, n_leaves+1, ... . Starting with the unique
    three-leaf tree, each new leaf is attached by subdividing every existing edge.
    """
    if n_leaves < 3:
        raise ValueError("Need at least three leaves")
    center = n_leaves
    trees: tuple[tuple[Edge, ...], ...] = (((0, center), (1, center), (2, center)),)
    for leaf in range(3, n_leaves):
        next_trees = []
        new_internal = n_leaves + leaf - 2
        for tree in trees:
            for edge in tree:
                a, b = edge
                remaining = [candidate for candidate in tree if candidate != edge]
                remaining.extend(((min(a, new_internal), max(a, new_internal)),
                                  (min(b, new_internal), max(b, new_internal)),
                                  (leaf, new_internal)))
                next_trees.append(tuple(sorted(remaining)))
        trees = tuple(next_trees)
    return trees


def component_leaves_after_vertex(
    tree: tuple[Edge, ...], removed: int, start: int, n_leaves: int
) -> set[int]:
    adjacency: dict[int, set[int]] = {}
    for a, b in tree:
        if a == removed or b == removed:
            continue
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)
    seen = {start}
    stack = [start]
    while stack:
        current = stack.pop()
        for neighbor in adjacency.get(current, set()) - seen:
            seen.add(neighbor)
            stack.append(neighbor)
    return {node for node in seen if node < n_leaves}


def cubic_tree_ansatz(all_momenta: tuple[Vec, ...]) -> Fraction:
    """2^(3-n) times the sum of products of local cubic-vertex brackets."""
    n = len(all_momenta)
    total = 0
    for tree in cubic_trees(n):
        adjacency: dict[int, set[int]] = {}
        for a, b in tree:
            adjacency.setdefault(a, set()).add(b)
            adjacency.setdefault(b, set()).add(a)
        term = 1
        for internal, neighbors in adjacency.items():
            if internal < n:
                continue
            branches = [
                component_leaves_after_vertex(tree, internal, neighbor, n)
                for neighbor in sorted(neighbors)
            ]
            if len(branches) != 3:
                raise AssertionError("Cubic tree has a non-trivalent internal vertex")
            p0 = add(*(all_momenta[i] for i in branches[0]))
            p1 = add(*(all_momenta[i] for i in branches[1]))
            term *= abs(bracket(p0, p1))
        total += term
    return Fraction(total, 2 ** (n - 3))


def random_sample(n: int, rng: random.Random) -> tuple[Vec, ...]:
    while True:
        positive = tuple(
            (rng.randint(-8, 8) or 1, rng.randint(-8, 8) or -1)
            for _ in range(n - 1)
        )
        minus = tuple(-x for x in add(*positive))
        all_momenta = positive + (minus,)
        if any(bracket(a, b) == 0 for i, a in enumerate(all_momenta) for b in all_momenta[i + 1:]):
            continue
        try:
            amplitude(positive)
        except ValueError:
            continue
        return all_momenta


def run_tests(max_n: int = 6, samples: int = 20, seed: int = 260304330) -> None:
    rng = random.Random(seed)
    for n in range(3, max_n + 1):
        print(f"n={n}: cubic trees={len(cubic_trees(n))}")
        for sample_number in range(1, samples + 1):
            momenta = random_sample(n, rng)
            actual = amplitude(momenta[:-1])
            candidate = cubic_tree_ansatz(momenta)
            if actual != candidate:
                print(f"  FAIL sample {sample_number}")
                print(f"  momenta={momenta}")
                print(f"  recursion={actual}")
                print(f"  cubic={candidate}")
                return
        print(f"  passed {samples}/{samples}")


def canonical_multigraph(edges: tuple[Edge, ...], n_vertices: int) -> tuple[int, ...]:
    """Canonical adjacency-multiplicity key under all vertex relabelings."""
    best = None
    for permutation in permutations(range(n_vertices)):
        counts = {}
        for a, b in edges:
            image = tuple(sorted((permutation[a], permutation[b])))
            counts[image] = counts.get(image, 0) + 1
        key = tuple(
            counts.get((a, b), 0)
            for a in range(n_vertices)
            for b in range(a + 1, n_vertices)
        )
        if best is None or key < best:
            best = key
    assert best is not None
    return best


@lru_cache(maxsize=None)
def symmetric_edge_monomial_orbits(n_vertices: int, degree: int):
    edges = tuple((a, b) for a in range(n_vertices) for b in range(a + 1, n_vertices))
    orbits: dict[tuple[int, ...], list[tuple[Edge, ...]]] = {}
    for monomial in combinations_with_replacement(edges, degree):
        key = canonical_multigraph(monomial, n_vertices)
        orbits.setdefault(key, []).append(monomial)
    return tuple(tuple(monomials) for monomials in orbits.values())


def symmetric_absolute_features(momenta: tuple[Vec, ...], degree: int) -> tuple[int, ...]:
    edge_value = {
        (a, b): abs(bracket(momenta[a], momenta[b]))
        for a in range(len(momenta))
        for b in range(a + 1, len(momenta))
    }
    return tuple(
        sum(
            product_value
            for monomial in orbit
            for product_value in [
                __import__("math").prod(edge_value[edge] for edge in monomial)
            ]
        )
        for orbit in symmetric_edge_monomial_orbits(len(momenta), degree)
    )


def fit_five_point_absolute_polynomial(samples: int = 40, seed: int = 5) -> None:
    """Test whether M5 is a symmetric cubic in |[ij]|."""
    rng = random.Random(seed)
    rows = []
    targets = []
    for _ in range(samples):
        momenta = random_sample(5, rng)
        rows.append(symmetric_absolute_features(momenta, degree=3))
        targets.append(amplitude(momenta[:-1]))

    def rref(values):
        matrix = [[Fraction(value) for value in row] for row in values]
        pivot_columns = []
        pivot_row = 0
        for column in range(len(matrix[0])):
            pivot = next((row for row in range(pivot_row, len(matrix)) if matrix[row][column]), None)
            if pivot is None:
                continue
            matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
            divisor = matrix[pivot_row][column]
            matrix[pivot_row] = [value / divisor for value in matrix[pivot_row]]
            for row in range(len(matrix)):
                if row == pivot_row or not matrix[row][column]:
                    continue
                factor = matrix[row][column]
                matrix[row] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(matrix[row], matrix[pivot_row])
                ]
            pivot_columns.append(column)
            pivot_row += 1
            if pivot_row == len(matrix):
                break
        return matrix, pivot_columns

    _, pivots = rref(rows)
    augmented, augmented_pivots = rref([list(row) + [target] for row, target in zip(rows, targets)])
    rank = len(pivots)
    augmented_rank = len(augmented_pivots)
    basis_size = len(rows[0])
    print(f"M5 symmetric-|bracket| basis size={basis_size}, rank={rank}")
    print(f"augmented rank={augmented_rank}")
    if rank == augmented_rank:
        if rank == basis_size:
            solution = [Fraction(0) for _ in range(basis_size)]
            for row_index, column in enumerate(augmented_pivots):
                if column < basis_size:
                    solution[column] = augmented[row_index][-1]
            print(f"unique solution={solution}")
        else:
            print("A fit exists but the sampled basis is rank-deficient.")
    else:
        print("No symmetric cubic polynomial in pairwise absolute brackets fits M5.")


if __name__ == "__main__":
    verify_published_low_points()
    verify_decay_chamber()
    run_tests()
    fit_five_point_absolute_polynomial()
    verify_flow_cone_equivalence()
    verify_permutation_invariance()
