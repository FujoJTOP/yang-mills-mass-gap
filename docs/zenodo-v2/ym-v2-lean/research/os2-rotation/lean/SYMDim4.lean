/-
  SYM-DIM4: no dimension-four anisotropy for the hypercubic group.

  The hyperoctahedral group `W₄ = (ℤ/2)⁴ ⋊ S₄` acts on symmetric 4×4 matrices by
  `A ↦ w A wᵀ`.  This file proves that the space of `W₄`-invariant symmetric
  matrices is one-dimensional, spanned by the identity.  Equivalently, the
  coefficient tensor `c_{μν}` of a hypercubic-invariant dimension-four operator
  `Σ_{μν} c_{μν} tr F_{μν}²` is forced to be `c δ_{μν}`, i.e. the operator is
  `O(4)`-invariant.  Hence lattice artifacts start at dimension six, which is the
  rigorous content of "no fine-tuning".

  Scope: this is the finite-dimensional representation theory only.  The step
  from "coefficient tensor is a multiple of `δ`" to "the dimension-four effective
  action has no anisotropy" is physics (the classification of dimension-four
  operators), and is deliberately NOT formalised here.

  Not part of any lake build.  Run `research/os2-rotation/check-symdim4.sh`.
-/
import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Real.Basic
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Logic.Equiv.Basic
import Mathlib.Tactic

open Matrix

namespace YM.OS2

/-! ## The two generating families of `W₄` -/

/-- Sign of a Boolean: `true ↦ -1`, `false ↦ 1`. -/
def sgn : Bool → ℝ := fun b => if b then -1 else 1

@[simp] lemma sgn_mul_sgn (b : Bool) : sgn b * sgn b = 1 := by
  cases b <;> norm_num [sgn]

@[simp] lemma sgn_not (b : Bool) : sgn (!b) = -sgn b := by
  cases b <;> norm_num [sgn]

/-- Sign flip `A ↦ D_s A D_s` with `D_s = diag (sgn ∘ s)`. -/
def signFlip (s : Fin 4 → Bool) (A : Matrix (Fin 4) (Fin 4) ℝ) : Matrix (Fin 4) (Fin 4) ℝ :=
  fun i j => sgn (s i) * A i j * sgn (s j)

/-- Coordinate permutation `A ↦ P_σ A P_σᵀ`, written with `σ` on the left. -/
def permute (σ : Equiv.Perm (Fin 4)) (A : Matrix (Fin 4) (Fin 4) ℝ) : Matrix (Fin 4) (Fin 4) ℝ :=
  fun i j => A (σ i) (σ j)

/-- Invariance under all sign flips. -/
def IsSignInvariant (A : Matrix (Fin 4) (Fin 4) ℝ) : Prop :=
  ∀ s : Fin 4 → Bool, signFlip s A = A

/-- Invariance under all coordinate permutations. -/
def IsPermInvariant (A : Matrix (Fin 4) (Fin 4) ℝ) : Prop :=
  ∀ σ : Equiv.Perm (Fin 4), permute σ A = A

/-- Invariance under the full hyperoctahedral group. -/
def IsInvariant (A : Matrix (Fin 4) (Fin 4) ℝ) : Prop :=
  IsSignInvariant A ∧ IsPermInvariant A

/-! ## Invariant matrices are multiples of the identity -/

/-- Flipping the sign of a single coordinate. -/
def flipAt (i : Fin 4) (s : Fin 4 → Bool) : Fin 4 → Bool :=
  Function.update s i (!s i)

lemma flipAt_involutive (i : Fin 4) : Function.Involutive (flipAt i) := by
  intro s
  funext k
  by_cases hk : k = i
  · subst hk
    simp [flipAt]
  · simp [flipAt, hk]

/-- A sum whose summand changes sign under an involution vanishes. -/
lemma sum_eq_zero_of_neg (i : Fin 4) (f : (Fin 4 → Bool) → ℝ)
    (hf : ∀ s, f (flipAt i s) = - f s) : (∑ s, f s) = 0 := by
  have h := Equiv.sum_comp (Function.Involutive.toPerm (flipAt i) (flipAt_involutive i)) f
  have h' : (∑ s, f s) = ∑ s, - f s := by
    rw [← h]
    exact Finset.sum_congr rfl (fun s _ => hf s)
  rw [Finset.sum_neg_distrib] at h'
  linarith

/-- Off-diagonal entries vanish: a single sign flip changes their sign. -/
lemma offdiag_eq_zero (A : Matrix (Fin 4) (Fin 4) ℝ) (h : IsInvariant A)
    {i j : Fin 4} (hij : i ≠ j) : A i j = 0 := by
  have h1 : signFlip (Function.update (fun _ => false) i true) A i j = A i j := by
    rw [h.1]
  have hi : (Function.update (fun _ => false) i true) i = true := by simp
  have hj : (Function.update (fun _ => false) i true) j = false := by
    simp [Function.update, Ne.symm hij]
  simp only [signFlip, hi, hj, sgn] at h1
  norm_num at h1
  linarith

/-- All diagonal entries are equal: a transposition moves one to the other. -/
lemma diag_eq (A : Matrix (Fin 4) (Fin 4) ℝ) (h : IsInvariant A) (i j : Fin 4) :
    A i i = A j j := by
  by_cases hij : i = j
  · subst hij; rfl
  · have h1 : permute (Equiv.swap i j) A i i = A i i := by rw [h.2]
    simp only [permute] at h1
    rw [Equiv.swap_apply_left i j] at h1
    exact h1.symm

/-- **Main lemma.** A `W₄`-invariant matrix is `(tr A / 4) • 1`. -/
theorem invariant_eq_trace_div_four (A : Matrix (Fin 4) (Fin 4) ℝ) (h : IsInvariant A) :
    A = (A.trace / 4) • (1 : Matrix (Fin 4) (Fin 4) ℝ) := by
  ext i j
  by_cases hij : i = j
  · subst hij
    have hconst : ∀ k : Fin 4, A k k = A i i := fun k => diag_eq A h k i
    have htr : A.trace = 4 * A i i := by
      simp only [Matrix.trace, Matrix.diag]
      rw [Finset.sum_congr rfl (fun k (_ : k ∈ Finset.univ) => hconst k)]
      simp [Finset.sum_const, Finset.card_univ]
    rw [htr]
    simp
  · have h0 : A i j = 0 := offdiag_eq_zero A h hij
    simp [hij, h0]

/-- **SYM-DIM4.** The `W₄`-invariant symmetric matrices are one-dimensional,
spanned by the identity: no dimension-four anisotropy. -/
theorem invariant_sym2_eq_smul_one (A : Matrix (Fin 4) (Fin 4) ℝ)
    (_hA : Aᵀ = A) (h : IsInvariant A) :
    ∃ c : ℝ, A = c • (1 : Matrix (Fin 4) (Fin 4) ℝ) :=
  ⟨A.trace / 4, invariant_eq_trace_div_four A h⟩

/-! ## Averaging over the two families -/

lemma signFlip_smul (s : Fin 4 → Bool) (c : ℝ) (A : Matrix (Fin 4) (Fin 4) ℝ) :
    signFlip s (c • A) = c • signFlip s A := by
  ext i j
  simp only [signFlip, Matrix.smul_apply]
  ring

lemma signFlip_sum {ι : Type*} [Fintype ι] (s : Fin 4 → Bool)
    (F : ι → Matrix (Fin 4) (Fin 4) ℝ) :
    signFlip s (∑ i, F i) = ∑ i, signFlip s (F i) := by
  ext i j
  simp only [signFlip, Matrix.sum_apply]
  rw [Finset.mul_sum, Finset.sum_mul]

lemma permute_smul (τ : Equiv.Perm (Fin 4)) (c : ℝ) (A : Matrix (Fin 4) (Fin 4) ℝ) :
    permute τ (c • A) = c • permute τ A := by
  ext i j
  simp only [permute, Matrix.smul_apply]

lemma permute_sum {ι : Type*} [Fintype ι] (τ : Equiv.Perm (Fin 4))
    (F : ι → Matrix (Fin 4) (Fin 4) ℝ) :
    permute τ (∑ i, F i) = ∑ i, permute τ (F i) := by
  ext i j
  simp only [permute, Matrix.sum_apply]

lemma permute_comp (τ σ : Equiv.Perm (Fin 4)) (A : Matrix (Fin 4) (Fin 4) ℝ) :
    permute τ (permute σ A) = permute (σ * τ) A := by
  ext i j
  rfl

lemma permute_diagonal (σ : Equiv.Perm (Fin 4)) (d : Fin 4 → ℝ) :
    permute σ (diagonal d) = diagonal (fun i => d (σ i)) := by
  ext i j
  simp only [permute, Matrix.diagonal, Matrix.of_apply]
  by_cases h : i = j
  · subst h; simp
  · have h' : σ i ≠ σ j := fun hh => h (σ.injective hh)
    simp [h, h']

/-- Average over the sign group `(ℤ/2)⁴`. -/
noncomputable def signAvg (A : Matrix (Fin 4) (Fin 4) ℝ) : Matrix (Fin 4) (Fin 4) ℝ :=
  ((Fintype.card (Fin 4 → Bool) : ℝ)⁻¹) • ∑ s : Fin 4 → Bool, signFlip s A

/-- Average over the coordinate permutations. -/
noncomputable def permAvg (A : Matrix (Fin 4) (Fin 4) ℝ) : Matrix (Fin 4) (Fin 4) ℝ :=
  ((Fintype.card (Equiv.Perm (Fin 4)) : ℝ)⁻¹) •
    ∑ σ : Equiv.Perm (Fin 4), permute σ A

lemma signFlip_diagonal (s : Fin 4 → Bool) (d : Fin 4 → ℝ) :
    signFlip s (diagonal d) = diagonal d := by
  ext i j
  by_cases hij : i = j
  · subst hij
    simp only [signFlip, Matrix.diagonal, Matrix.of_apply]
    have h : sgn (s i) * sgn (s i) = 1 := sgn_mul_sgn _
    calc sgn (s i) * d i * sgn (s i) = d i * (sgn (s i) * sgn (s i)) := by ring
      _ = d i := by rw [h, mul_one]
  · simp [signFlip, Matrix.diagonal, hij]

/-- The sign average is the diagonal. -/
lemma signAvg_eq_diagonal (A : Matrix (Fin 4) (Fin 4) ℝ) :
    signAvg A = diagonal (fun i => A i i) := by
  have hsum : ∀ i j : Fin 4,
      (∑ s : Fin 4 → Bool, signFlip s A i j)
        = A i j * (∑ s : Fin 4 → Bool, sgn (s i) * sgn (s j)) := by
    intro i j
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl (fun s _ => by simp only [signFlip]; ring)
  have hprod : ∀ i j : Fin 4,
      (∑ s : Fin 4 → Bool, sgn (s i) * sgn (s j))
        = if i = j then (Fintype.card (Fin 4 → Bool) : ℝ) else 0 := by
    intro i j
    by_cases hij : i = j
    · subst hij
      have h1 : ∀ s : Fin 4 → Bool, sgn (s i) * sgn (s i) = 1 := fun s => sgn_mul_sgn _
      simp [h1, Finset.sum_const, Finset.card_univ]
    · rw [if_neg hij]
      refine sum_eq_zero_of_neg i _ ?_
      intro s
      have h1 : flipAt i s i = !(s i) := by simp [flipAt]
      have h2 : flipAt i s j = s j := by
        simp [flipAt, Function.update, Ne.symm hij]
      simp only [h1, h2]
      rw [sgn_not]
      ring
  ext i j
  simp only [signAvg, Matrix.smul_apply, Matrix.sum_apply, hsum i j, Matrix.diagonal,
    Matrix.of_apply]
  rw [hprod i j]
  by_cases hij : i = j
  · subst hij
    have hc : (Fintype.card (Fin 4 → Bool) : ℝ) = 16 := by
      have h : Fintype.card (Fin 4 → Bool) = 16 := by decide
      rw [h]; norm_num
    simp only [hc, ite_true, smul_eq_mul]
    ring
  · rw [if_neg hij]; simp [hij]

lemma permute_trace (σ : Equiv.Perm (Fin 4)) (A : Matrix (Fin 4) (Fin 4) ℝ) :
    (permute σ A).trace = A.trace := by
  simp only [Matrix.trace, Matrix.diag, permute]
  exact Equiv.sum_comp σ (fun k => A k k)

lemma permAvg_trace (A : Matrix (Fin 4) (Fin 4) ℝ) :
    (permAvg A).trace = A.trace := by
  unfold permAvg
  rw [Matrix.trace_smul, Matrix.trace_sum]
  have h : ∀ σ : Equiv.Perm (Fin 4), (permute σ A).trace = A.trace := fun σ => permute_trace σ A
  simp only [h]
  rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
  have hc : (Fintype.card (Equiv.Perm (Fin 4)) : ℝ) = 24 := by
    have h24 : Fintype.card (Equiv.Perm (Fin 4)) = 24 := by decide
    rw [h24]; norm_num
  rw [hc, smul_eq_mul]
  ring

lemma permAvg_invariant (A : Matrix (Fin 4) (Fin 4) ℝ) : IsPermInvariant (permAvg A) := by
  intro τ
  simp only [permAvg]
  rw [permute_smul, permute_sum]
  have h : (∑ σ : Equiv.Perm (Fin 4), permute τ (permute σ A))
      = ∑ σ : Equiv.Perm (Fin 4), permute σ A := by
    rw [Finset.sum_congr rfl (fun σ (_ : σ ∈ Finset.univ) => permute_comp τ σ A)]
    exact Equiv.sum_comp (Equiv.mulRight τ) (fun σ => permute σ A)
  rw [h]

lemma signFlip_permAvg_diagonal (s : Fin 4 → Bool) (d : Fin 4 → ℝ) :
    signFlip s (permAvg (diagonal d)) = permAvg (diagonal d) := by
  simp only [permAvg]
  rw [signFlip_smul, signFlip_sum]
  have h : (∑ σ : Equiv.Perm (Fin 4), signFlip s (permute σ (diagonal d)))
      = ∑ σ : Equiv.Perm (Fin 4), permute σ (diagonal d) := by
    refine Finset.sum_congr rfl (fun σ _ => ?_)
    rw [permute_diagonal, signFlip_diagonal]
  rw [h]

/-- **Averaging identity.** Averaging a matrix over the two generating families
(equivalently over the hyperoctahedral group) projects onto `(tr A / 4) • 1`. -/
theorem average_eq_trace_div_four (A : Matrix (Fin 4) (Fin 4) ℝ) :
    permAvg (signAvg A) = (A.trace / 4) • (1 : Matrix (Fin 4) (Fin 4) ℝ) := by
  rw [signAvg_eq_diagonal]
  have hinv : IsInvariant (permAvg (diagonal (fun i => A i i))) :=
    ⟨fun s => signFlip_permAvg_diagonal s _, permAvg_invariant _⟩
  have ht : (permAvg (diagonal (fun i => A i i))).trace = A.trace := by
    rw [permAvg_trace, Matrix.trace_diagonal]
    simp [Matrix.trace, Matrix.diag]
  rw [invariant_eq_trace_div_four _ hinv, ht]

/-! ## Axiom footprint (frozen) -/

/-- info: 'YM.OS2.invariant_eq_trace_div_four' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YM.OS2.invariant_eq_trace_div_four

/-- info: 'YM.OS2.invariant_sym2_eq_smul_one' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YM.OS2.invariant_sym2_eq_smul_one

/-- info: 'YM.OS2.average_eq_trace_div_four' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YM.OS2.average_eq_trace_div_four

end YM.OS2
