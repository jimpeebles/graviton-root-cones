# Response to the adversarial referee report

We thank the reviewer for independently reconstructing the published recursion, checking every main result, and identifying several important framing and reproducibility issues. The manuscript and release have been revised as follows.

## Major recommendations

### 2.1 — VibeMathed framing and the prior open problem

The abstract, introduction, and discussion now frame the contribution as a precise **partial resolution of the matrix-tree obstruction** stated in arXiv:2603.04330. We no longer imply a general closed form for unrestricted amplitudes. The new statement is that the global cut conditions are root-cone membership conditions for every chamber, and the ordinary directed matrix-tree theorem applies whenever the resulting feasible-tree family is arborescent; those chambers are classified exhaustively through five blocks.

### 2.2 — “Physical” versus kinematically realizable

“Physical chamber” has been replaced throughout by “realizable chamber.” Theorem 4.1 is now the “Realizable chamber criterion.” A new remark distinguishes generic real-spinor realizability from the energy-sign and in/out restrictions of a physical scattering channel. The text and table caption state that the count 296 concerns labeled chambers under a fixed slope ordering.

### 2.3 — Overstatement of the wall result

The abstract now says that kinematic walls are **contained among** resonance hyperplanes. No converse is claimed.

### 2.4 — Shared implementation path and missing published anchor

`amplitude_search.py` now directly asserts the displayed three- and four-point amplitudes and checks the published decay-region product at four and five points. The global bracket-sign convention and its degree-dependent sign are documented in the code and manuscript. The release also includes the referee's two independently written programs, which do not import the bundled implementation.

### 2.5 — Generalized retarded functions and Steinmann relations

The introduction now discusses the classical connection between resonance chambers and generalized retarded functions and cites Evans, Epstein, Liu–Norledge–Ocneanu, and Norledge–Ocneanu. The novelty claim is correspondingly narrowed to the flow-cone identification for this graviton vertex, the realizability result, the arborescent classification, and the five-point chamber formula.

## Minor recommendations

- **3.1:** Berends–Giele and Hodges are now cited in the introduction; the unused gluon reference was removed.
- **3.2:** The manuscript now explains why the cut-side choice in Eq. (1) is immaterial and how it matches the source convention.
- **3.3:** Direct low-point assertions were added to the verification program and described in the manuscript.
- **3.4:** A one-paragraph derivation of the five-point recursion from Eq. (33) and the preamplitude base cases was added.
- **3.5:** The Table 1 caption identifies the population over which each row is counted; the 22,102 checks are explicitly described as a bounded-grid enumeration.
- **3.6:** The absence of singleton feasible families at five blocks is noted after Table 1.
- **3.7:** Theorem 4.1 now states the relevant genericity conditions explicitly.
- **3.8:** The release is pinned to Python 3.12.3, includes the decay anchor, emits machine-readable chamber certificates, and contains a continuous-integration workflow.

## Additional release changes

- The provisional “VibeMathed Research” author line was replaced with a named human author and independent-researcher designation.
- A detailed human/AI contribution statement was added.
- The adversarial report and independent programs are included for transparency, while the manuscript continues to state that unaffiliated human expert review remains outstanding.
