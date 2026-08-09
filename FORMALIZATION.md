# Lean formalization scope and statement-audit guide

This document maps the manuscript's informal claims to the declarations in
`GravitonRootCones/RootCone.lean`. It is intended to make an independent
statement audit concrete. It is not itself such an audit.

## Reproducible environment

- Lean: `v4.32.2`, pinned in `lean-toolchain`
- mathlib: `v4.32.2`, pinned in `lake-manifest.json`
- Build: `lake build --wfail`
- Compiled-environment recheck in CI: Lean Action's `leanchecker`

Every principal theorem is guarded by mathlib's `assert_no_sorry` command and
followed by `#print axioms`. The expected output may contain only `propext`,
`Classical.choice`, and `Quot.sound`. It must not contain `sorryAx` or any
project-defined axiom.

## Manuscript-to-Lean map

| Manuscript content | Lean declaration | What is checked |
|---|---|---|
| Antisymmetric two-spinor bracket | `bracket`, `bracket_self` | The coordinate definition and `[a,a]=0` over any commutative ring. |
| Linearity used in the cut calculation | `bracket_add_left`, `bracket_add_right`, `bracket_finset_sum_left` | Additivity and distribution across a finite momentum sum. |
| Equation `cutisflow`: `[q_A,Q-q_A] = sum_(i in A) [q_i,Q]` | `cut_bracket_eq_sum_projected` | The exact algebraic identity, with no kinematic assumptions. |
| `sum_i x_i=[Q,Q]=0` | `sum_projected_eq_zero` | The projected vector lies in the zero-sum hyperplane. |
| Root-vector convention `e_tail-e_head` | `divergence` | Divergence is outgoing flow minus incoming flow, matching the manuscript. |
| Cancellation of internal edges across a cut | `sum_divergence_eq_sum_cutContribution` | Total divergence on any finite vertex set equals its directed cut balance. |
| Deleting an edge of a tree leaves one crossing edge | `DirectedTreeCuts` | The exact endpoint-membership and noncrossing properties used in the proof are explicit hypotheses. |
| A tree divergence fixes each edge flow uniquely | `sum_divergence_treeSide_eq_flow`, `flow_unique_of_same_divergence` | The flow on an edge is the divergence sum on its tail-side component, hence uniqueness. |
| Positive flow iff all tree-cut sums are positive | `prescribed_cut_positive_iff_positive_flow` | The central strict-inequality equivalence for a prescribed divergence. |
| Open root cone iff the tree-cut tests pass | `InOpenRootCone`, `inOpenRootCone_iff_treeCuts` | Root-cone membership is equivalent to the strict cut inequalities once a realizing tree flow is supplied. |
| Retarded factor `-cutBracket/pairBracket` | `RetardedCutPass`, `retardedCutPass_iff_positive_flow_of_orientation`, `all_retardedCutPass_iff_positive_flow` | Both signs of `[uv]` give exactly positive edge flow, individually and for all tree edges, away from walls. |

## Deliberate abstraction boundary

`DirectedTreeCuts` records the only graph-theoretic fact used by the proof: for
each deleted edge, its chosen tail-side component contains the tail, excludes
the head, and no other edge crosses that cut. The current artifact does not yet
construct this structure from mathlib's undirected `SimpleGraph.IsTree`. That
bridge is mathematically routine but important for eliminating a statement-
correspondence gap.

Likewise, the formalization does not yet encode the sum over all weighted
spanning trees, the advanced vertex, the resonance-chamber classification, the
realizability count, the determinant special cases, or the five-point formula.

## Independent statement-audit checklist

An auditor should check, without relying on the proof scripts, that:

1. `bracket` matches the manuscript's two-spinor convention.
2. `divergence` matches the root generator `e_i-e_j` for an edge `i -> j`.
3. Choosing `DirectedTreeCuts.side e` as the component containing `tail e`
   translates the manuscript's possibly opposite choice of `A_e` with the
   correct sign.
4. `DirectedTreeCuts.other_not_crossing` follows from deleting an edge of a
   spanning tree, and the formal hypotheses do not admit a relevant loophole.
5. `InOpenRootCone` is extensionally the manuscript's `C^circ(T)`.
6. The manuscript's Heaviside test is exactly the cut inequality appearing in
   `prescribed_cut_positive_iff_positive_flow` after applying the orientation
   convention `[ij] < 0`.
7. The theorem scope claimed in `VIBEMATHED-SUBMISSION.md` does not extend to the
   computational claims that remain outside Lean.

Until an independent reviewer completes and records this checklist, the honest
status is **Lean-checked, statement unaudited**.
