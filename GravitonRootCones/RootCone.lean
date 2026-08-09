/-
Copyright (c) 2026 James Kehoe. All rights reserved.
Released under MIT license as described in the file LICENSE.
-/
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Tactic.Ring
import Mathlib.Util.AssertNoSorry

set_option linter.style.header false

/-!
# Root-cone formalization for half-collinear graviton vertices

This module will formalize the central equivalence used in Theorem 3.1 of the
paper: the edge-cut sign conditions on an oriented tree are equivalent to the
existence of a strictly positive flow having the prescribed divergence.

The first milestone below records the algebraic cut identity independently of
the tree argument.  It is intentionally small so that the public project has a
kernel-checked foundation before the graph layer is added.
-/

namespace GravitonRootCones

open scoped BigOperators

/-- The antisymmetric bracket on two-component vectors. -/
def bracket {R : Type*} [CommRing R] (a b : Fin 2 → R) : R :=
  a 0 * b 1 - a 1 * b 0

theorem bracket_self {R : Type*} [CommRing R] (a : Fin 2 → R) :
    bracket a a = 0 := by
  simp only [bracket]
  ring

theorem bracket_add_right {R : Type*} [CommRing R] (a b c : Fin 2 → R) :
    bracket a (b + c) = bracket a b + bracket a c := by
  simp [bracket]
  ring

theorem bracket_add_left {R : Type*} [CommRing R] (a b c : Fin 2 → R) :
    bracket (a + b) c = bracket a c + bracket b c := by
  simp [bracket]
  ring

theorem bracket_sub_self_right {R : Type*} [CommRing R] (a Q : Fin 2 → R) :
    bracket a (Q - a) = bracket a Q := by
  simp [bracket]
  ring

theorem bracket_finset_sum_left {I R : Type*} [CommRing R]
    (q : I → Fin 2 → R) (Q : Fin 2 → R) (s : Finset I) :
    bracket (∑ i ∈ s, q i) Q = ∑ i ∈ s, bracket (q i) Q := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [bracket]
  | @insert a s ha ih =>
      simp only [Finset.sum_insert ha]
      rw [bracket_add_left, ih]

/-- Equation (3.5) of the manuscript: a cut bracket is the sum of the
projected vectors on either side of the cut. -/
theorem cut_bracket_eq_sum_projected {I R : Type*} [CommRing R]
    (q : I → Fin 2 → R) (Q : Fin 2 → R) (side : Finset I) :
    bracket (∑ i ∈ side, q i) (Q - ∑ i ∈ side, q i) =
      ∑ i ∈ side, bracket (q i) Q := by
  rw [bracket_sub_self_right, bracket_finset_sum_left]

/-- The projected vector `x i = bracket (q i) Q` lies in the zero-sum
hyperplane when `Q` is the total momentum. -/
theorem sum_projected_eq_zero {I R : Type*} [Fintype I] [CommRing R]
    (q : I → Fin 2 → R) :
    ∑ i, bracket (q i) (∑ j, q j) = 0 := by
  classical
  rw [← bracket_finset_sum_left q (∑ j, q j) Finset.univ]
  exact bracket_self _

section DirectedCuts

variable {V E R : Type*}
variable [Fintype E]
variable [DecidableEq V]
variable [CommRing R]

/-- Divergence in the manuscript's convention: outgoing flow minus incoming flow. -/
def divergence (tail head : E → V) (flow : E → R) (v : V) : R :=
  ∑ e, ((if v = tail e then flow e else 0) - (if v = head e then flow e else 0))

/-- The contribution of one directed edge to the net inflow of a vertex set. -/
def cutContribution (tail head : E → V) (flow : E → R)
    (side : Finset V) (e : E) : R :=
  (if tail e ∈ side then flow e else 0) - (if head e ∈ side then flow e else 0)

/-- Summing divergence over a set of vertices leaves only its directed cut. -/
theorem sum_divergence_eq_sum_cutContribution
    (tail head : E → V) (flow : E → R) (side : Finset V) :
    ∑ v ∈ side, divergence tail head flow v =
      ∑ e, cutContribution tail head flow side e := by
  simp only [divergence, cutContribution]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro e _
  simp only [Finset.sum_sub_distrib]
  congr 1 <;> simp

/-- Data supplied by the two components obtained after deleting each edge of a tree.
For the chosen component, the deleted edge leaves it and every other edge has both
endpoints on the same side. -/
structure DirectedTreeCuts (tail head : E → V) where
  side : E → Finset V
  tail_mem : ∀ e, tail e ∈ side e
  head_not_mem : ∀ e, head e ∉ side e
  other_not_crossing : ∀ e e', e' ≠ e → (head e' ∈ side e ↔ tail e' ∈ side e)

/-- For an edge-deletion component, its total divergence is exactly the flow on
the unique directed edge leaving that component. -/
theorem sum_divergence_treeSide_eq_flow
    (tail head : E → V) (flow : E → R) (cuts : DirectedTreeCuts tail head) (e : E) :
    ∑ v ∈ cuts.side e, divergence tail head flow v = flow e := by
  rw [sum_divergence_eq_sum_cutContribution]
  rw [Fintype.sum_eq_single e]
  · simp [cutContribution, cuts.tail_mem, cuts.head_not_mem]
  · intro e' he'
    have hsame := cuts.other_not_crossing e e' he'
    by_cases hhead : head e' ∈ cuts.side e
    · have htail : tail e' ∈ cuts.side e := hsame.mp hhead
      simp [cutContribution, hhead, htail]
    · have htail : tail e' ∉ cuts.side e := fun h ↦ hhead (hsame.mpr h)
      simp [cutContribution, hhead, htail]

/-- A flow with a prescribed divergence is unique on a directed tree. -/
theorem flow_unique_of_same_divergence
    (tail head : E → V) (cuts : DirectedTreeCuts tail head) (flow₁ flow₂ : E → R)
    (hdiv : ∀ v, divergence tail head flow₁ v = divergence tail head flow₂ v) :
    flow₁ = flow₂ := by
  funext e
  rw [← sum_divergence_treeSide_eq_flow tail head flow₁ cuts e]
  rw [← sum_divergence_treeSide_eq_flow tail head flow₂ cuts e]
  apply Finset.sum_congr rfl
  intro v _
  exact hdiv v

/-- The retarded cut inequalities for a directed tree are equivalent to strict
positivity of every edge flow.  The theorem only needs an ordering relation:
all graph-theoretic and algebraic content is in `sum_divergence_treeSide_eq_flow`. -/
theorem tree_cut_positive_iff_positive_flow [LT R]
    (tail head : E → V) (flow : E → R) (cuts : DirectedTreeCuts tail head) :
    (∀ e, 0 < ∑ v ∈ cuts.side e, divergence tail head flow v) ↔
      ∀ e, 0 < flow e := by
  constructor
  · intro h e
    rw [← sum_divergence_treeSide_eq_flow tail head flow cuts e]
    exact h e
  · intro h e
    rw [sum_divergence_treeSide_eq_flow tail head flow cuts e]
    exact h e

/-- Paper-facing form: if `x` is the prescribed divergence of a tree flow,
then the strict tree-cut inequalities for `x` hold exactly when the flow is
strictly positive. -/
theorem prescribed_cut_positive_iff_positive_flow [LT R]
    (tail head : E → V) (x : V → R) (flow : E → R)
    (cuts : DirectedTreeCuts tail head)
    (hdiv : ∀ v, divergence tail head flow v = x v) :
    (∀ e, 0 < ∑ v ∈ cuts.side e, x v) ↔ ∀ e, 0 < flow e := by
  have hsum (e : E) : ∑ v ∈ cuts.side e, x v = flow e := by
    rw [← sum_divergence_treeSide_eq_flow tail head flow cuts e]
    apply Finset.sum_congr rfl
    intro v _
    exact (hdiv v).symm
  constructor
  · intro h e
    rw [← hsum e]
    exact h e
  · intro h e
    rw [hsum e]
    exact h e

/-- Membership in the open root cone generated by the directed edges. -/
def InOpenRootCone [LT R] (tail head : E → V) (x : V → R) : Prop :=
  ∃ flow : E → R, (∀ e, 0 < flow e) ∧ ∀ v, divergence tail head flow v = x v

/-- Once a prescribed divergence is realized on a tree, root-cone membership
is exactly the family of strict inequalities on its edge-deletion cuts. -/
theorem inOpenRootCone_iff_treeCuts [LT R]
    (tail head : E → V) (x : V → R) (flow : E → R)
    (cuts : DirectedTreeCuts tail head)
    (hdiv : ∀ v, divergence tail head flow v = x v) :
    InOpenRootCone tail head x ↔ ∀ e, 0 < ∑ v ∈ cuts.side e, x v := by
  constructor
  · rintro ⟨other, hpos, hother⟩
    have heq : other = flow := flow_unique_of_same_divergence tail head cuts other flow fun v ↦
      (hother v).trans (hdiv v).symm
    rw [heq] at hpos
    exact (prescribed_cut_positive_iff_positive_flow tail head x flow cuts hdiv).mpr hpos
  · intro hcuts
    refine ⟨flow, ?_, hdiv⟩
    exact (prescribed_cut_positive_iff_positive_flow tail head x flow cuts hdiv).mp hcuts

end DirectedCuts

assert_no_sorry bracket_self
assert_no_sorry bracket_add_right
assert_no_sorry bracket_add_left
assert_no_sorry bracket_sub_self_right
assert_no_sorry bracket_finset_sum_left
assert_no_sorry cut_bracket_eq_sum_projected
assert_no_sorry sum_projected_eq_zero
assert_no_sorry sum_divergence_eq_sum_cutContribution
assert_no_sorry sum_divergence_treeSide_eq_flow
assert_no_sorry flow_unique_of_same_divergence
assert_no_sorry tree_cut_positive_iff_positive_flow
assert_no_sorry prescribed_cut_positive_iff_positive_flow
assert_no_sorry inOpenRootCone_iff_treeCuts

#print axioms bracket_self
#print axioms bracket_add_right
#print axioms bracket_add_left
#print axioms bracket_sub_self_right
#print axioms bracket_finset_sum_left
#print axioms cut_bracket_eq_sum_projected
#print axioms sum_projected_eq_zero
#print axioms sum_divergence_eq_sum_cutContribution
#print axioms sum_divergence_treeSide_eq_flow
#print axioms flow_unique_of_same_divergence
#print axioms tree_cut_positive_iff_positive_flow
#print axioms prescribed_cut_positive_iff_positive_flow
#print axioms inOpenRootCone_iff_treeCuts

end GravitonRootCones
