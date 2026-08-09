import Lake

open Lake DSL

package GravitonRootCones where
  version := v!"0.1.0"

require "leanprover-community" / "mathlib" @ git "v4.32.2"

@[default_target]
lean_lib GravitonRootCones where
  leanOptions := #[
    ⟨`pp.unicode.fun, true⟩,
    ⟨`relaxedAutoImplicit, false⟩,
    ⟨`weak.linter.mathlibStandardSet, true⟩,
    ⟨`maxSynthPendingDepth, .ofNat 3⟩
  ]
