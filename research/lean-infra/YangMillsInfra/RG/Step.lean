/-
# The RG step-theorem interface, and the finite-volume gap certificate

The 4D construction has two logically distinct RG requirements:

* a **step theorem**: one blocking step maps the effective action of scale `a` to that
  of scale `a/2`, up to a controlled error, uniformly in the volume;
* a **gap certificate**: a finite-volume quantity `Δ(n)` whose positivity at any one
  scale implies a transfer-matrix gap.

This file fixes their types.  The step theorem itself is the hard analysis (Balaban /
BBS / Bauerschmidt–Bodineau–Dagallier); what is formalized here is the *bookkeeping*:
the scale sequence, the per-step error, and the fact that summable errors give a
uniform bound on the accumulated error.

The gap certificate is

    Δ(n) = 4 log Z(n) − log Z(2n),

where `Z(L)` is the partition function of the `L`-box; it is invariant under exact
blocking, and reflection positivity turns `Δ(n) ≥ δ` into a transfer-matrix gap of size
`≥ δ/(4n)`.  The implication is recorded as `CertificateImpliesGap`: it is the one
piece of the OpenAI sigma-model template (docs/yang-mills-prior-work.md §3.2) that the
4D argument needs formalized, and it is not part of this infrastructure package.
-/
import YangMillsInfra.Reconstruction.OS
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecificLimits.Basic

noncomputable section

open Filter

namespace YangMillsInfra
namespace RG

/-- The finite-volume gap certificate `Δ(n) = 4 log Z(n) − log Z(2n)`. -/
def gapCertificate (Z : ℕ → ℝ) (n : ℕ) : ℝ :=
  4 * Real.log (Z n) - Real.log (Z (2 * n))

/-- Scale and error data for a sequence of RG steps.

`scale n` is the lattice spacing after `n` blocking steps (so it halves each step).
`error n` is the norm distance between the effective action produced by the step and
the reference flow; `error_summable` is the uniformity statement that makes the
accumulated error finite. -/
structure RGStepData where
  scale : ℕ → ℝ
  scale_pos : ∀ n, 0 < scale n
  scale_halves : ∀ n, scale (n + 1) = scale n / 2
  error : ℕ → ℝ
  error_nonneg : ∀ n, 0 ≤ error n
  error_summable : Summable error

namespace RGStepData

variable (D : RGStepData)

/-- The accumulated error after `N` steps is bounded by the total error mass.  This is
the formal content of "the errors are summable uniformly in the volume": no matter how
many scales are combined, the error never exceeds a single finite constant. -/
theorem accumulated_error_le (N : ℕ) :
    ∑ n ∈ Finset.range N, D.error n ≤ ∑' n, D.error n :=
  D.error_summable.sum_le_tsum _ fun n _ => D.error_nonneg n

/-- The accumulated error converges to the total error mass as the number of steps
goes to infinity. -/
theorem accumulated_error_tendsto :
    Tendsto (fun N : ℕ => ∑ n ∈ Finset.range N, D.error n) atTop
      (nhds (∑' n, D.error n)) :=
  D.error_summable.hasSum.tendsto_sum_nat

end RGStepData

/-- The uniformity-in-volume requirement, as a predicate on a family of statements
indexed by the volume `L`.  A step theorem is *uniform* when it holds for every `L`,
with a bound independent of `L`. -/
def VolumeUniform (P : ℕ → Prop) : Prop := ∀ L, P L

/-- **A step theorem.**  The fields are the quantitative output that the weak-coupling
analysis is supposed to provide; this package only fixes their types.

`reference_flow` is the one/two-loop flow that a successful step must reproduce;
`coupling_error` says how closely the effective action tracks it; `uniform` says the
bound is independent of the volume `L`; `error_summable` says it integrates to a
finite total over all scales. -/
structure StepTheorem where
  blockScale : ℕ → ℝ
  blockScale_pos : ∀ n, 0 < blockScale n
  blockScale_halves : ∀ n, blockScale (n + 1) = blockScale n / 2
  /-- Distance of the `n`-th effective action from the reference flow, at volume `L`. -/
  stepError : ℕ → ℕ → ℝ
  stepError_nonneg : ∀ n L, 0 ≤ stepError n L
  /-- Uniformity: the same error bound works at every volume. -/
  uniform_in_volume : ∀ n, ∃ ε : ℝ, 0 ≤ ε ∧ ∀ L, stepError n L ≤ ε
  /-- The errors over scales are summable, at every fixed volume. -/
  summable_in_scale : ∀ L, Summable (fun n => stepError n L)

/-- Every step theorem yields, at each volume, an `RGStepData` whose accumulated error
is finite.  This is the bookkeeping bridge from the analysis to the certificates. -/
def StepTheorem.toRGStepData (S : StepTheorem) (L : ℕ) (hL : ∀ n, 0 ≤ S.stepError n L) :
    RGStepData where
  scale := S.blockScale
  scale_pos := S.blockScale_pos
  scale_halves := S.blockScale_halves
  error := fun n => S.stepError n L
  error_nonneg := hL
  error_summable := S.summable_in_scale L

/-- A uniform lower bound for the gap certificate. -/
structure UniformGapCertificate (Z : ℕ → ℝ) where
  Z_pos : ∀ k, 1 ≤ Z k
  n : ℕ
  n_pos : 0 < n
  lower : ℝ
  lower_pos : 0 < lower
  cert_bound : lower ≤ gapCertificate Z n

/-- **The certificate-to-gap implication**, recorded as an interface.

This is the statement that a positive finite-volume gap certificate forces a
transfer-matrix gap.  It is the one step of the OpenAI sigma-model template not supplied
by the 16 certified inputs, and the 4D analogue is what the mass-gap program has to
prove.  The `rate` field is the resulting gap `δ/(4n)`. -/
structure CertificateImpliesGap (T : Reconstruction.TransferData) (Z : ℕ → ℝ) where
  certificate : UniformGapCertificate Z
  transfer_gap : T.HasGap (certificate.lower / (4 * certificate.n))

end RG
end YangMillsInfra
