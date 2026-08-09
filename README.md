# Half-collinear graviton vertices as weighted root-cone enumerators

[![Verification](https://github.com/jimpeebles/graviton-root-cones/actions/workflows/verify.yml/badge.svg)](https://github.com/jimpeebles/graviton-root-cones/actions/workflows/verify.yml)

This repository accompanies a working mathematical-physics paper on the general-region half-collinear single-minus graviton recursion of Guevara, Lupsasca, Skinner, Strominger, and Weil ([arXiv:2603.04330](https://arxiv.org/abs/2603.04330)).

## Main result

For block spinors `q_i`, let `Q = sum_i q_i`, set `x_i = [q_i,Q]`, and orient `i -> j` when `[ij] < 0`. The retarded multipoint vertex is exactly the weighted sum of directed spanning trees whose open root cone contains `x`. The advanced vertex uses `-x`.

This converts the source paper's global cut tests into root-cone membership conditions and places every possible kinematic wall among the hyperplanes of the resonance arrangement. In chambers where the feasible trees form all in-arborescences of a directed graph, the ordinary directed Matrix-Tree Theorem gives a determinant.

The release also contains:

- a real-spinor realizability theorem and the count `(m-1) R_(m-1) / m`;
- exhaustive classifications of all 32 four-block and 370 five-block resonance chambers;
- 10 four-block and 57 five-block determinant checks;
- a compact five-graviton formula in a chamber outside the published decay region;
- 22,102 exact bounded-grid checks of that formula;
- an adversarial referee report and two independent implementations written without importing the primary code.

## Status

The mathematics passed an adversarial independent LLM review, including fresh derivations and implementations. No fatal error was found. It has **not yet received review from an unaffiliated human expert** in scattering amplitudes or resonance arrangements. Treat this as an openly checkable working preprint, not an established result.

The claim is intentionally narrow: this is a partial resolution of the matrix-tree obstruction stated in arXiv:2603.04330, not a general closed form for unrestricted half-collinear amplitudes.

## Repository layout

- `paper/` — manuscript PDF and LaTeX source.
- `src/amplitude_search.py` — published recursion, direct low-point and decay anchors, flow-cone checks, and failed stronger ansatz tests.
- `src/chamber_enumeration.py` — exact chamber classification and directed-Laplacian determinants.
- `src/five_point_identity.py` — exhaustive check of the non-decay five-point formula.
- `data/chamber-certificates.json` — machine-readable representatives, feasible trees, and roots for every realizable chamber.
- `review/` — adversarial report, response, and independent verification programs.
- `VIBEMATHED-SUBMISSION.md` — proposed tracker entry and honest scope caveat.

## Reproduce

Python 3.12.3 is pinned; the programs use only the standard library.

```bash
python src/amplitude_search.py
python src/chamber_enumeration.py --json /tmp/chamber-certificates.json
cmp data/chamber-certificates.json /tmp/chamber-certificates.json
python src/five_point_identity.py
python review/independent_check.py
PYTHONPATH=review python review/independent_table_check.py
```

The GitHub Actions workflow runs this complete verification matrix on every push.

## Independent review wanted

The highest-value review tasks are:

1. audit the flow-cone proof and real-spinor normalization;
2. search for prior appearances of this specific graviton/root-cone identification;
3. verify the five-point chamber formula directly from the published recursion;
4. assess whether the result qualifies for VibeMathed as a partial answer to the source paper's stated matrix-tree obstruction.

Issues and pull requests containing counterexamples or corrected references are especially welcome.

## AI contribution disclosure

James Kehoe selected and directed the research problem and commissioned adversarial review. OpenAI Codex generated the central reformulation, proofs, verification code, and initial manuscript under that direction. A separate LLM review independently re-derived the results using fresh code. No AI system is listed as an author; the named human author is responsible for the released claims.

## License and citation

Code is MIT-licensed. The paper and prose are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Citation metadata is in `CITATION.cff`.
