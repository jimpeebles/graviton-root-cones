# Proposed VibeMathed submission

This is a draft record, not a claim that the curator has accepted the result as in scope.

## Result

**Name:** Half-collinear graviton vertices as weighted root-cone enumerators

**Short name:** Graviton root cones

**Result:** Proved

**Status:** Partial result

**Field:** Mathematical physics

**Field detail:** Scattering amplitudes; resonance arrangements; directed spanning trees

**Method:** Argument

**Solve date:** 2026-08-09

**AI contribution:** AI-discovered

**Model:** OpenAI Codex (GPT-5 family)

**Vendor:** OpenAI

**Publication:** Announced

**Verification:** Unreviewed

## Plain-language statement

The general half-collinear single-minus graviton recursion of Guevara–Lupsasca–Skinner–Strominger–Weil contains multipoint vertices whose tree weights depend on global cut tests, obstructing a direct matrix-tree formula outside a special decay region. The new result proves that these cut tests are exactly positive-flow conditions: each retarded vertex is a weighted enumerator of directed spanning-tree root cones containing a kinematic netflow vector. Consequently, the directed Matrix-Tree Theorem applies whenever the feasible trees form a complete arborescence family. These chambers are classified exhaustively through five blocks, and one gives a compact five-graviton formula outside the decay region.

## Result qualifier

Partial resolution of the matrix-tree obstruction and general-kinematics simplification problem stated in arXiv:2603.04330; it does not give a single determinant or closed form in every chamber.

## What the AI did

Under human direction, OpenAI Codex proposed the root-cone interpretation, developed the proof, wrote exact enumeration and verification programs, found the non-decay five-point chamber identity, and drafted the manuscript. A separate LLM referee independently reconstructed the published recursion and re-derived every principal claim using fresh code.

## Verification note

The theorem, chamber counts, determinant identities, and five-point formula were independently reproduced by an adversarial LLM using separately written exact-arithmetic programs. No unaffiliated human domain expert has yet endorsed the result, so the correct verification label is Unreviewed.

## Prior-question anchor

arXiv:2603.04330 says that simplification of the general solution, if possible, is left to future work, and separately explains that the ordinary matrix-tree theorem cannot yet be applied because each factor depends on the global tree structure. This work resolves that obstruction in root-cone language and obtains determinants on a classified family of chambers. Because the original prompt is open-ended and the resolution is partial, curator confirmation of scope should precede submission.

## Primary source

**Source name:** GitHub release v0.2.0 — reviewed working preprint

**Source URL:** https://github.com/jimpeebles/graviton-root-cones/releases/tag/v0.2.0

The release contains the manuscript PDF and a frozen source archive. If the manuscript is later archived on Zenodo or arXiv, update the publication field to Preprint and replace or supplement this source URL.

## Additional links

- Repository: https://github.com/jimpeebles/graviton-root-cones
- Passing verification workflow: https://github.com/jimpeebles/graviton-root-cones/actions/runs/31339397929
- Public review response: https://github.com/jimpeebles/graviton-root-cones/blob/main/review/response-to-referee.md

## Recommended curator inquiry

> We have an AI-discovered working preprint addressing the matrix-tree obstruction explicitly stated in arXiv:2603.04330. It proves that the global cut conditions are directed root-cone membership conditions, classifies the chambers where the Matrix-Tree Theorem consequently applies through five blocks, and derives a non-decay five-graviton formula. The result is public and exactly reproducible but has only adversarial machine review, not human expert endorsement. Would this qualify as an unreviewed partial result under VibeMathed's inclusion criterion?
