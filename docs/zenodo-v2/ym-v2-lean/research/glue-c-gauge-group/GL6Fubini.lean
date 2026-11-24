/-
GL6: finite orbit Fubini / gauge-fixing independence, machine-checked.

This is the finite-lattice layer only.  No 4D nonabelian Gribov problem is
touched: we work with a finite type `X`, a finite group `Γ`, a `Γ`-invariant
weight `μ`, a `Γ`-invariant observable `F`, and a gauge slice `A` meeting every
orbit.  The theorem is the finite Faddeev--Popov identity

    ∑_x F x μ x  =  ∑_{a ∈ A} F a μ a · |orbit a| / |A ∩ orbit a| .

The factor `1/|A ∩ orbit a|` is the Gribov-copy bookkeeping: dropping it is
wrong whenever the multiplicity is non-constant (`gl6_check.py` exhibits
m ∈ {1,2,3,4,6,12} for S₃ on a 2×1 torus).

A canonical representative is chosen as the minimum of each orbit (hence the
hypothesis `[LinearOrder X]`; any finite type admits one).  This avoids
constructing the orbit-space measure, which is the one step Mathlib does not
yet package for compact group quotients.

The concrete `example` at the end is an exact rational enumeration for S₃
acting on `Fin 3` (non-free: stabilizers of order 2), with slices of
multiplicity 1 and 2.

Run: bash check-glue-c.sh  (or `lean GL6Fubini.lean` with the pinned toolchain)
-/

import Mathlib

set_option linter.unusedSectionVars false

noncomputable section

namespace GlueC

open Finset

variable {X Γ : Type*}

section Orbit

variable [Fintype X] [DecidableEq X] [Fintype Γ] [Group Γ] [MulAction Γ X]

/-- The orbit of `x`, as a finset. -/
def orbitF (x : X) : Finset X := Finset.univ.image (fun g : Γ => g • x)

omit [Fintype X] in
@[simp] lemma mem_orbitF {x y : X} :
    y ∈ orbitF (Γ := Γ) x ↔ ∃ g : Γ, g • x = y := by
  simp [orbitF]

lemma orbitF_nonempty (x : X) : (orbitF (Γ := Γ) x).Nonempty :=
  ⟨x, (mem_orbitF).mpr ⟨1, one_smul Γ x⟩⟩

lemma mem_orbitF_comm {x y : X} :
    x ∈ orbitF (Γ := Γ) y ↔ y ∈ orbitF (Γ := Γ) x := by
  rw [mem_orbitF, mem_orbitF]
  constructor
  · rintro ⟨g, hg⟩
    exact ⟨g⁻¹, by rw [← hg, inv_smul_smul]⟩
  · rintro ⟨g, hg⟩
    exact ⟨g⁻¹, by rw [← hg, inv_smul_smul]⟩

omit [Fintype X] in
lemma orbitF_eq_of_mem {x y : X} (hy : y ∈ orbitF (Γ := Γ) x) :
    orbitF (Γ := Γ) y = orbitF (Γ := Γ) x := by
  obtain ⟨h, rfl⟩ := mem_orbitF.mp hy
  ext z
  simp only [mem_orbitF]
  constructor
  · rintro ⟨g, hg⟩
    exact ⟨g * h, by rw [mul_smul, hg]⟩
  · rintro ⟨g, hg⟩
    exact ⟨g * h⁻¹, by rw [mul_smul, inv_smul_smul, hg]⟩

variable [LinearOrder X]

/-- Canonical representative of an orbit: its minimum. -/
def orbitMin (x : X) : X := (orbitF (Γ := Γ) x).min' (orbitF_nonempty x)

lemma orbitMin_mem (x : X) : orbitMin (Γ := Γ) x ∈ orbitF (Γ := Γ) x :=
  Finset.min'_mem _ _

lemma orbitMin_eq_of_mem {x y : X} (hy : y ∈ orbitF (Γ := Γ) x) :
    orbitMin (Γ := Γ) y = orbitMin (Γ := Γ) x := by
  have hEq : orbitF (Γ := Γ) y = orbitF (Γ := Γ) x := orbitF_eq_of_mem hy
  unfold orbitMin
  apply le_antisymm
  · apply Finset.min'_le
    show (orbitF (Γ := Γ) x).min' _ ∈ orbitF (Γ := Γ) y
    rw [hEq]
    exact Finset.min'_mem _ _
  · apply Finset.min'_le
    show (orbitF (Γ := Γ) y).min' _ ∈ orbitF (Γ := Γ) x
    rw [← hEq]
    exact Finset.min'_mem _ _

lemma orbitMin_idem (x : X) :
    orbitMin (Γ := Γ) (orbitMin (Γ := Γ) x) = orbitMin (Γ := Γ) x :=
  orbitMin_eq_of_mem (orbitMin_mem x)

lemma orbitMin_fixed {y : X}
    (hy : y ∈ Finset.univ.image (orbitMin (Γ := Γ))) :
    orbitMin (Γ := Γ) y = y := by
  obtain ⟨z, -, hz⟩ := Finset.mem_image.mp hy
  rw [← hz, orbitMin_idem]

lemma orbitMin_eq_iff {x y : X}
    (hy : y ∈ Finset.univ.image (orbitMin (Γ := Γ))) :
    orbitMin (Γ := Γ) x = y ↔ x ∈ orbitF (Γ := Γ) y := by
  constructor
  · intro h
    have hx : y ∈ orbitF (Γ := Γ) x := h ▸ orbitMin_mem x
    exact mem_orbitF_comm.mpr hx
  · intro hx
    have h1 : orbitMin (Γ := Γ) x = orbitMin (Γ := Γ) y :=
      orbitMin_eq_of_mem (x := y) (y := x) hx
    rw [h1, orbitMin_fixed hy]

/-- **GL6 finite orbit Fubini.**  For a finite `Γ`-invariant weight `μ`, a
`Γ`-invariant `F`, and a slice `A` meeting every orbit,

`∑_x F x μ x = ∑_{a∈A} F a μ a · |orbit a| / |A ∩ orbit a|`.

The multiplicity `|A ∩ orbit a|` is the Faddeev--Popov orbit volume; dropping it
is the standard Gribov-copy error. -/
theorem orbit_sum_identity
    (F μ : X → ℝ)
    (hF : ∀ (g : Γ) (x : X), F (g • x) = F x)
    (hμ : ∀ (g : Γ) (x : X), μ (g • x) = μ x)
    (A : Finset X)
    (hA : ∀ x : X, ∃ a ∈ A, ∃ g : Γ, g • x = a) :
    ∑ x, F x * μ x
      = ∑ a ∈ A, F a * μ a * ((orbitF (Γ := Γ) a).card : ℝ) /
          ((A.filter (fun b => b ∈ orbitF (Γ := Γ) a)).card : ℝ) := by
  have hfinv : ∀ (g : Γ) (x : X), F (g • x) * μ (g • x) = F x * μ x := by
    intro g x; rw [hF, hμ]
  -- Left side: group by the canonical representative.
  have hL : (∑ x, F x * μ x)
      = ∑ y : X, ∑ x ∈ Finset.univ.filter (fun x => orbitMin (Γ := Γ) x = y),
          F x * μ x := by
    rw [← Finset.sum_fiberwise_of_maps_to
      (s := (Finset.univ : Finset X)) (t := (Finset.univ : Finset X))
      (g := fun x => orbitMin (Γ := Γ) x) (fun x _ => Finset.mem_univ _)]
  -- Inner sum on the left, for an actual representative.
  have hinnerL : ∀ y : X, y ∈ Finset.univ.image (orbitMin (Γ := Γ)) →
      (∑ x ∈ Finset.univ.filter (fun x => orbitMin (Γ := Γ) x = y), F x * μ x)
        = ((orbitF (Γ := Γ) y).card : ℝ) * (F y * μ y) := by
    intro y hy
    have hfilter : Finset.univ.filter (fun x => orbitMin (Γ := Γ) x = y)
        = orbitF (Γ := Γ) y := by
      ext x
      simp only [Finset.mem_filter, Finset.mem_univ, true_and]
      exact orbitMin_eq_iff hy
    rw [hfilter]
    have hconst : ∀ x ∈ orbitF (Γ := Γ) y, F x * μ x = F y * μ y := by
      intro x hx
      obtain ⟨g, hg⟩ := mem_orbitF.mp hx
      rw [← hg, hfinv]
    rw [Finset.sum_congr rfl hconst, Finset.sum_const, nsmul_eq_mul]
  -- Empty fibers off the image of the representative map.
  have hempty : ∀ y : X, y ∉ Finset.univ.image (orbitMin (Γ := Γ)) →
      Finset.univ.filter (fun x => orbitMin (Γ := Γ) x = y) = ∅ ∧
      A.filter (fun a => orbitMin (Γ := Γ) a = y) = ∅ := by
    intro y hy
    constructor
    · rw [Finset.eq_empty_iff_forall_notMem]
      intro x hx
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hx
      exact hy (hx ▸ Finset.mem_image_of_mem (orbitMin (Γ := Γ)) (Finset.mem_univ x))
    · rw [Finset.eq_empty_iff_forall_notMem]
      intro a ha
      simp only [Finset.mem_filter] at ha
      exact hy (ha.2 ▸ Finset.mem_image_of_mem (orbitMin (Γ := Γ)) (Finset.mem_univ a))
  -- Right side: group `A` by the same representative.
  have hR : (∑ a ∈ A, F a * μ a * ((orbitF (Γ := Γ) a).card : ℝ) /
        ((A.filter (fun b => b ∈ orbitF (Γ := Γ) a)).card : ℝ))
      = ∑ y : X, ∑ a ∈ A.filter (fun a => orbitMin (Γ := Γ) a = y),
          F a * μ a * ((orbitF (Γ := Γ) a).card : ℝ) /
            ((A.filter (fun b => b ∈ orbitF (Γ := Γ) a)).card : ℝ) := by
    rw [← Finset.sum_fiberwise_of_maps_to
      (s := A) (t := (Finset.univ : Finset X))
      (g := fun a => orbitMin (Γ := Γ) a) (fun a _ => Finset.mem_univ _)]
  -- Inner sum on the right: every term is constant, and the fiber has exactly
  -- `|A ∩ orbit y|` terms, so the multiplicity cancels.
  have hinnerR : ∀ y : X, y ∈ Finset.univ.image (orbitMin (Γ := Γ)) →
      (∑ a ∈ A.filter (fun a => orbitMin (Γ := Γ) a = y),
          F a * μ a * ((orbitF (Γ := Γ) a).card : ℝ) /
            ((A.filter (fun b => b ∈ orbitF (Γ := Γ) a)).card : ℝ))
        = ((orbitF (Γ := Γ) y).card : ℝ) * (F y * μ y) := by
    intro y hy
    have hterm : ∀ a ∈ A.filter (fun a => orbitMin (Γ := Γ) a = y),
        F a * μ a * ((orbitF (Γ := Γ) a).card : ℝ) /
            ((A.filter (fun b => b ∈ orbitF (Γ := Γ) a)).card : ℝ)
          = F y * μ y * ((orbitF (Γ := Γ) y).card : ℝ) /
              ((A.filter (fun b => orbitMin (Γ := Γ) b = y)).card : ℝ) := by
      intro a ha
      simp only [Finset.mem_filter] at ha
      have hay : a ∈ orbitF (Γ := Γ) y := (orbitMin_eq_iff hy).mp ha.2
      have horb : orbitF (Γ := Γ) a = orbitF (Γ := Γ) y := orbitF_eq_of_mem hay
      have hfa : F a * μ a = F y * μ y := by
        obtain ⟨g, hg⟩ := mem_orbitF.mp hay
        rw [← hg, hfinv]
      have hden : A.filter (fun b => b ∈ orbitF (Γ := Γ) a)
          = A.filter (fun b => orbitMin (Γ := Γ) b = y) := by
        ext b
        simp only [Finset.mem_filter]
        rw [horb]
        constructor
        · rintro ⟨hb, hbo⟩
          exact ⟨hb, (orbitMin_eq_iff hy).mpr hbo⟩
        · rintro ⟨hb, hbo⟩
          exact ⟨hb, (orbitMin_eq_iff hy).mp hbo⟩
      rw [hfa, hden, horb]
    rw [Finset.sum_congr rfl hterm, Finset.sum_const, nsmul_eq_mul]
    have hmpos : 0 < (A.filter (fun b => orbitMin (Γ := Γ) b = y)).card := by
      obtain ⟨z, -, hz⟩ := Finset.mem_image.mp hy
      obtain ⟨a, haA, g, hg⟩ := hA z
      have hama : orbitMin (Γ := Γ) a = y := by
        have h1 : orbitMin (Γ := Γ) a = orbitMin (Γ := Γ) z :=
          orbitMin_eq_of_mem (x := z) (y := a) (mem_orbitF.mpr ⟨g, hg⟩)
        rw [h1, hz]
      exact Finset.card_pos.mpr ⟨a, Finset.mem_filter.mpr ⟨haA, hama⟩⟩
    have hmne : (((A.filter (fun b => orbitMin (Γ := Γ) b = y)).card : ℕ) : ℝ) ≠ 0 := by
      exact_mod_cast ne_of_gt hmpos
    field_simp [hmne]
  rw [hL, hR]
  refine Finset.sum_congr rfl ?_
  intro y _
  by_cases hy : y ∈ Finset.univ.image (orbitMin (Γ := Γ))
  · rw [hinnerL y hy, hinnerR y hy]
  · obtain ⟨h1, h2⟩ := hempty y hy
    rw [h1, h2]
    simp

/-- **GL6 gauge-fixing independence.**  Two slices meeting every orbit give
the same invariant expectation.  This is the "for each fixed gauge choice"
statement; neither slice needs to be a fundamental domain. -/
theorem gauge_fixing_independence
    (F μ : X → ℝ)
    (hF : ∀ (g : Γ) (x : X), F (g • x) = F x)
    (hμ : ∀ (g : Γ) (x : X), μ (g • x) = μ x)
    (A₁ A₂ : Finset X)
    (hA₁ : ∀ x : X, ∃ a ∈ A₁, ∃ g : Γ, g • x = a)
    (hA₂ : ∀ x : X, ∃ a ∈ A₂, ∃ g : Γ, g • x = a) :
    (∑ a ∈ A₁, F a * μ a * ((orbitF (Γ := Γ) a).card : ℝ) /
        ((A₁.filter (fun b => b ∈ orbitF (Γ := Γ) a)).card : ℝ))
      = (∑ a ∈ A₂, F a * μ a * ((orbitF (Γ := Γ) a).card : ℝ) /
          ((A₂.filter (fun b => b ∈ orbitF (Γ := Γ) a)).card : ℝ)) := by
  rw [← orbit_sum_identity (Γ := Γ) F μ hF hμ A₁ hA₁,
      ← orbit_sum_identity (Γ := Γ) F μ hF hμ A₂ hA₂]

end Orbit

/-! ## Exact S₃ enumeration

A concrete non-free instance: `S₃ = Equiv.Perm (Fin 3)` acts on `Fin 3`
(stabilizers of order 2).  With `F = μ = 1`, the two slices `{0}` and `{0,1}`
have multiplicities `m = 1` and `m = 2`; both weighted sums equal `3`.  This is
an exact rational identity, checked by `native_decide`. -/

example :
    (∑ a ∈ ({0} : Finset (Fin 3)),
        (1 : ℚ) * (1 : ℚ) * ((orbitF (Γ := Equiv.Perm (Fin 3)) a).card : ℚ) /
          ((({0} : Finset (Fin 3)).filter
            (fun b => b ∈ orbitF (Γ := Equiv.Perm (Fin 3)) a)).card : ℚ))
      = (∑ a ∈ ({0, 1} : Finset (Fin 3)),
          (1 : ℚ) * (1 : ℚ) * ((orbitF (Γ := Equiv.Perm (Fin 3)) a).card : ℚ) /
            ((({0, 1} : Finset (Fin 3)).filter
              (fun b => b ∈ orbitF (Γ := Equiv.Perm (Fin 3)) a)).card : ℚ)) := by
  native_decide

end GlueC

/-- info: 'GlueC.orbit_sum_identity' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms GlueC.orbit_sum_identity

/-- info: 'GlueC.gauge_fixing_independence' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms GlueC.gauge_fixing_independence
