# Referee report

**Manuscript:** *Half-Collinear Graviton Vertices as Weighted Root-Cone Enumerators* (working draft, Aug 9, 2026)
**Reviewer stance:** skeptical mathematical-physics referee, per the enclosed review request. I did not assume correctness because tests pass; every principal claim was re-derived or re-computed with code written independently of the bundle.

## Summary of the submission

The draft takes the retarded/advanced multipoint vertices of Guevara–Lupsasca–Skinner–Strominger–Weil (arXiv:2603.04330), whose Cayley-tree weights carry global step functions, and shows each vertex equals a weighted enumerator of directed spanning trees whose root cone (in the sense of Gutekunst–Mészáros–Petersen, arXiv:1903.06595) contains the projected netflow vector x_i = [q_i, Q]. Consequences: vertex walls lie on the resonance arrangement; a criterion for which resonance chambers are realizable by slope-ordered real spinors, with the closed count (m−1)R_{m−1}/m; an exhaustive m = 4, 5 classification identifying chambers where the directed matrix-tree theorem applies; and an explicit compact five-graviton formula in a chamber outside the published decay region.

## 1. Fatal issues

None found. Specifically:

1.1. **Theorem 3.1 (root-cone form) is correct.** I re-derived the proof by hand. The key identity [A_e, B_e] = Σ_{i∈A_e} x_i is elementary and correct; the case analysis on sign([uv]) matches the orientation convention i→j ⇔ [ij] < 0; the advanced statement follows by global sign reversal. I also re-implemented both sides from scratch and confirmed equality on 300 exact random samples (60 per valency, m = 2…6), independently of the bundled code.

1.2. **The vertex and recursion conventions match the published paper exactly.** I downloaded arXiv:2603.04330 and read the equations directly (this matters: two automated summaries of the paper misreported the constraint in eq. (33) as A > 2; the PDF says A ≥ 2). The draft's eq. (1) matches the published eqs. (30)–(31) verbatim, including the Θ-argument, and the draft's five-point recursion is exactly the published eqs. (28)–(29)+(33) specialized to n = 5 using M̄_{i} = 1, M̄_{ij} = 0, M̄_{ijk} = −V_{ijk}. I confirmed empirically that the A ≥ 2 convention (used by the draft and its code) — and not A ≥ 3 — is the one that (a) reproduces the published decay-chamber closed form eq. (26)/(39) at n = 5, and (b) is fully S_n permutation invariant, which the published paper asserts below eq. (33).

1.3. **The five-graviton corollary is correct.** Working only from the sign data of Table 2 (no kinematics), I independently resolved every step function in the 16 four-block Cayley trees and the 3 trees of each triple subvertex. I recover exactly the claimed families F(V) = {{14,23,24},{14,24,34}}, F(V̄) = {{12,13,14},{13,14,23}}, V_123 = V_234 = 0, V_124 = w12·w24, V_134 = w14·w34, hence eq. (newM5) holds on the entire open chamber, confirming the draft's finite analytic argument. Numerically, my own implementation of the published recursion agrees with eq. (newM5) on 40/40 fresh chamber points.

1.4. **The counts are correct.** My own grid enumeration recovers 32 and 370 resonance chambers (matching R_3, R_4, i.e. OEIS A034997), 24 and 296 chambers with at least one negative prefix sum, and 8 = 32/4 and 74 = 370/5 chambers in the aligned positive root cone — consistent with Theorem 4.1 and eq. (physicalcount).

1.5. **Not circular, not already known in the same form (as far as I could determine).** The cited GMP paper supplies the root-cone/resonance machinery but contains no physics; the graviton paper contains no resonance-arrangement language and explicitly states both that "we cannot yet apply the matrix-tree theorem to this formula because each factor in the product depends on the global structure of the tree" and that "simplification of the general solution (B11), if possible, is left to future work." INSPIRE lists 6 papers citing arXiv:2603.04330 to date; none addresses chamber combinatorics or simplification of the half-collinear recursion. Targeted searches for prior connections (half-collinear × resonance arrangement / flow cones / Kostant / matrix-tree) returned nothing. Caveat: the source paper is ~5 months old and this cannot rule out unposted concurrent work.

## 2. Major issues

2.1. **Fit with the venue's inclusion criterion is genuinely arguable and should be addressed head-on.** VibeMathed catalogues *previously open problems now resolved*. The anchor here is soft: "simplification, if possible, is left to future work" is an open-ended request, not a precisely stated question, and the draft answers it *partially* (chamber-by-chamber determinants at low rank plus one new compact formula — not a general closed form). The strongest crisply resolvable statement the draft actually settles is the obstruction sentence: the published claim that the matrix-tree theorem "cannot yet" be applied because the step tests are global. The draft proves the global tests are exactly root-cone membership and exhibits (with proof) chambers where MTT does apply, including outside the decay region. My recommendation: restate the contribution as resolving that specific published obstruction in a precise, checkable form ("the global cut conditions are equivalent to a root-cone membership; consequently MTT applies in exactly these chambers, classified exhaustively for m ≤ 5"), and let the venue judge. Claiming a fully "answered open question" without this reframing invites rejection.

2.2. **"Physical" is doing unexamined work.** Theorem 4.1's notion — realizable by generic real two-spinors with slope-ordered labels and normalized Q — is kinematic realizability, not physicality in the scattering sense: the published paper's decay region additionally fixes energy signs (ω_n < 0, ω_a > 0), i.e. in/out assignments. A chamber counted "physical" here may or may not be reachable in a given physical channel. Rename (e.g. "realizable chambers") or add a remark separating the two notions, and state explicitly that the count 296 is of *labeled* chambers under the slope-order convention.

2.3. **The abstract overstates Corollary 3.2.** "The kinematic walls are precisely the hyperplanes of the resonance arrangement" claims equality; the text proves one direction (walls ⊆ resonance hyperplanes, plus constancy of tree support on chambers). That every resonance hyperplane is genuinely a wall (tree support changes across it, for some tournament) is plausible but not proven. Either prove the converse or weaken "precisely" to "contained among."

2.4. **All bundled numerical evidence shares one code path.** chamber_enumeration.py and graviton-five-point-identity.py both import amplitude_search.py, so every bundled check inherits any convention error in the single `vertex` implementation, and the bundle never anchors against the published paper's independent closed forms. The missing anchor is the decay-chamber product formula (eqs. (26)/(39) of arXiv:2603.04330). I performed this test myself: sampling the decay region (ω_a > 0, z̃_1 < … < z̃_{n−2} < z̃_n < z̃_{n−1}), the recursion reproduces the published product exactly at n = 4 and up to a global sign at n = 5 — the flip is the expected [ij] ↔ −[ij] convention ambiguity, invisible at even bracket degree and harmless to the draft (which uses only |[ij]| and sign ratios). Add this test to the bundle, and note the bracket-sign convention explicitly.

2.5. **Missing literature thread: generalized retarded functions and the Steinmann relations.** The resonance arrangement's chambers famously count generalized retarded functions in axiomatic QFT (A034997; Epstein–Glaser–Ruelle; more recently Norledge–Ocneanu, arXiv:1901.03243, on the adjoint braid arrangement via Steinmann relations — cited in GMP's own bibliography). The draft's vertices are literally *retarded* and *advanced* kernels, so the appearance of the resonance arrangement is unlikely to be a coincidence and has classical antecedents in spirit. This must be discussed: it simultaneously strengthens the result's plausibility and tempers the novelty of "retarded objects live on resonance chambers" as an abstract idea. The concrete new content — the flow-cone identification for *this* vertex, the realizability count, the arborescent classification, and eq. (newM5) — survives, but an expert referee will raise this immediately.

## 3. Minor issues

3.1. Three of five bibliography entries (GuevaraGluon, BerendsGiele, Hodges) are never cited in the text. Cite or cut.

3.2. Eq. (1): "Choosing u ∈ A_e" — add a half-sentence noting the Θ-argument is invariant under (u,A_e) ↔ (v,B_e), so the choice is immaterial (the published paper fixes u to the smaller label; the invariance is why the conventions agree).

3.3. The claim "reproduction of the published three- and four-point amplitudes" is only implicit in the code (via the cubic ansatz passing at n = 3, 4). Add direct assertions M_123 = |[12]| and M_1234 = ½(|12||34| + |13||24| + |14||23|).

3.4. Eq. (M5rec) is stated as "the published five-point recursion" — give the one-line derivation from eq. (33) with the M̄ base cases so the reader need not re-derive which partitions survive.

3.5. Table 1 caption: state that "distinct feasible-tree families" and all subsequent rows are counted over *realizable* chambers only, and that "22,102 exact chamber samples" is a bounded-grid enumeration, not random sampling.

3.6. The observation that at m = 5 every nonempty arborescent chamber is nontrivial (no singleton families occur at all — family sizes jump 0 → 2) is curious and worth a sentence; I confirmed it independently.

3.7. The proof of Theorem 4.1 silently uses that a generic chamber has all x_i ≠ 0 and all prefix sums nonzero (both are resonance coordinates); one clause making this explicit would help.

3.8. Reproducibility: pin a Python version, add the decay-region test (2.4), and consider emitting machine-readable chamber certificates as the draft itself proposes.

## 4. Independent checks performed

All checks below used fresh code written for this review (no imports from the bundle), with exact integer/rational arithmetic; plus direct reading of the two source papers.

- Reconstructed the published vertex (eqs. 30–31) and the full recursion (eqs. 28–29, 33) from the arXiv PDF; confirmed the draft's conventions match, and resolved the A ≥ 2 vs A > 2 ambiguity from the PDF text itself.
- Verified Theorem 3.1's proof line by line (orientation, open cones, cut identity, advanced relation) and numerically on 300 exact samples, m = 2…6.
- Verified the summation-by-parts identity and both directions of Theorem 4.1.
- Confirmed GMP Corollary 3 (R_n = (n+1)·R_n⁺, via the cyclic decomposition of Lemma 6) in the published EJC version, and that the aligned positive root cone equals {all prefix sums ≥ 0} (path-flow argument), so eq. (physicalcount) is fully supported.
- Re-enumerated resonance chambers (32, 370), realizable chambers (24, 296), and aligned-cone chambers (8, 74).
- Re-derived the entire Table 1 classification independently: distinct families 9/55, V = 0 chambers 6/38, arborescent nonempty 10/57, nontrivial 6/57, max family 5/28 — all match; all 67 matrix-tree cofactor identities re-verified with my own Bareiss determinant against my own vertex.
- Verified the directed matrix-tree convention (outgoing Laplacian, root-deleted minor = in-arborescences toward the root) on the draft's Laplacian definition.
- Symbolically re-derived the five-point chamber families from Table 2's signs alone, and numerically confirmed eq. (newM5) against my own recursion on 40 fresh chamber points. Attempted counterexample search within the listed inequalities found none (consistent with the finite symbolic argument, which is the actual proof).
- Anchored the recursion against the published decay-chamber closed form (26)/(39) at n = 4, 5 — a test absent from the bundle.
- Ran all three bundled programs to completion: outputs match every number quoted in the draft and README (including the two documented *failures* of the stronger conjectures, which reproduce as described: cubic ansatz fails at n = 5; symmetric-cubic fit has rank 7 vs augmented rank 8).
- Literature searches: INSPIRE citation list of arXiv:2603.04330 (6 citing papers, none on this topic); targeted web searches for prior resonance/flow-cone connections to half-collinear amplitudes (none found); identified the uncited Steinmann/generalized-retarded-function thread (2.5).

## 5. Verdict

**Plausible preprint after specified revisions** — with a caveat on venue fit.

On mathematical validity I found no errors: every theorem, count, and formula I tested independently is correct, and the one convention question I could manufacture (A ≥ 2 vs A ≥ 3) resolves in the draft's favor from the primary source. The work is a modest but genuine partial answer to a stated obstruction in a very recent paper, executed carefully and honestly (the draft's own framing of its limits, and its inclusion of failed conjectures, is exemplary). The required revisions are 2.1–2.5 above: reframe the resolved question precisely, fix the "physical"/"precisely" wording, add the decay-chamber anchor test, and engage the Steinmann literature. Suitability for VibeMathed specifically turns on whether a *partial* simplification against an open-ended "left to future work" counts under their criterion; the reframing in 2.1 gives it the best honest shot. Independent human expert review (amplitudes + arrangements) remains necessary before any claim of novelty is made publicly, exactly as the draft itself says.
