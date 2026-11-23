/-
# Lattice gauge theory: combinatorics, gauge action, Wilson action

This is the *finite, exact* layer of "measure theory on gauge orbits".  On a finite
lattice every gauge-theoretic object is finite-dimensional, so nothing here is
approximate: the gauge group is `Λ.V → G`, the configuration space is `Λ.L → G`, and
gauge invariance of the Wilson action is an algebraic identity.

The point of formalizing this layer first is that it is where the *non-vacuous*
statement of gauge invariance lives.  The continuum story (principal bundles,
connections) has no mathlib support at all (see REPORT.md, gap G-1), while the lattice
story is complete and is exactly what the 4D construction uses.

Main results:
* `gaugeAct` is a `MulAction` of `Λ.V → G` on `Λ.L → G`;
* `holonomy_gaugeAct`: a plaquette holonomy conjugates by the value of the gauge
  transformation at the base point (the telescoping product);
* `wilsonAction_invariant`: the Wilson action with any class function is gauge invariant.

Reference: Seiler, *Gauge Theories as a Problem of Constructive Quantum Field Theory
and Statistical Mechanics*, LNP 159 (1982), Ch. II; Osterwalder–Seiler, Ann. Phys. 110
(1978) 440.
-/
import Mathlib.GroupTheory.GroupAction.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Group

noncomputable section

namespace YangMillsInfra
namespace LatticeGauge

universe u

/-- Combinatorial data of a finite lattice.

`V` are vertices, `L` directed links, `Plaq` plaquettes.  Each link has a reverse link
(`rev`), and each plaquette is an ordered 4-cycle `p1 → p2 → p3 → p4 → p1` of links,
with the four boundary conditions recorded as fields.

This is enough to define gauge transformations, holonomies and the Wilson action.  It
is *not* a theory of lattices (no embedding into `ℤ^d`); that is deliberate, since the
mass-gap statements only need the finite-volume theory at each `Λ`, plus a compatible
family of `Λ`'s for the infinite-volume limit. -/
structure LatticeData where
  V : Type u
  L : Type u
  Plaq : Type u
  instV : Fintype V
  instDecV : DecidableEq V
  instL : Fintype L
  instDecL : DecidableEq L
  instP : Fintype Plaq
  src : L → V
  tgt : L → V
  rev : L → L
  rev_src : ∀ l, src (rev l) = tgt l
  rev_tgt : ∀ l, tgt (rev l) = src l
  rev_rev : ∀ l, rev (rev l) = l
  p1 : Plaq → L
  p2 : Plaq → L
  p3 : Plaq → L
  p4 : Plaq → L
  h12 : ∀ q, src (p2 q) = tgt (p1 q)
  h23 : ∀ q, src (p3 q) = tgt (p2 q)
  h34 : ∀ q, src (p4 q) = tgt (p3 q)
  h41 : ∀ q, src (p1 q) = tgt (p4 q)

attribute [instance] LatticeData.instV LatticeData.instDecV LatticeData.instL
  LatticeData.instDecL LatticeData.instP

namespace LatticeData

variable (Λ : LatticeData.{u})

/-- Gauge field configurations: one group element per directed link. -/
abbrev Connection (G : Type u) : Type u := Λ.L → G

/-- The gauge group of a finite lattice: one group element per vertex. -/
abbrev GaugeGroup (G : Type u) : Type u := Λ.V → G

end LatticeData

variable {Λ : LatticeData.{u}} {G : Type u} [Group G]

/-- Gauge transformation of a link field: `U_l ↦ g(src l) * U_l * g(tgt l)⁻¹`. -/
def gaugeAct (g : Λ.V → G) (U : Λ.L → G) : Λ.L → G :=
  fun l => g (Λ.src l) * U l * (g (Λ.tgt l))⁻¹

instance : SMul (Λ.V → G) (Λ.L → G) := ⟨gaugeAct⟩

@[simp] lemma gaugeAct_def (g : Λ.V → G) (U : Λ.L → G) (l : Λ.L) :
    (g • U) l = g (Λ.src l) * U l * (g (Λ.tgt l))⁻¹ := rfl

instance : MulAction (Λ.V → G) (Λ.L → G) where
  one_smul U := by
    funext l
    simp [gaugeAct_def]
  mul_smul g h U := by
    funext l
    simp only [gaugeAct_def, Pi.mul_apply, mul_inv_rev]
    group

/-- The holonomy around a plaquette: the ordered product of its four link variables. -/
def holonomy (U : Λ.L → G) (q : Λ.Plaq) : G :=
  U (Λ.p1 q) * U (Λ.p2 q) * U (Λ.p3 q) * U (Λ.p4 q)

/-- **Telescoping**: under a gauge transformation a plaquette holonomy conjugates by the
value of the gauge transformation at the base point of the plaquette. -/
theorem holonomy_gaugeAct (g : Λ.V → G) (U : Λ.L → G) (q : Λ.Plaq) :
    holonomy (g • U) q = g (Λ.src (Λ.p1 q)) * holonomy U q * (g (Λ.src (Λ.p1 q)))⁻¹ := by
  have h1 : g (Λ.src (Λ.p2 q)) = g (Λ.tgt (Λ.p1 q)) := by rw [Λ.h12]
  have h2 : g (Λ.src (Λ.p3 q)) = g (Λ.tgt (Λ.p2 q)) := by rw [Λ.h23]
  have h3 : g (Λ.src (Λ.p4 q)) = g (Λ.tgt (Λ.p3 q)) := by rw [Λ.h34]
  have h4 : g (Λ.tgt (Λ.p4 q)) = g (Λ.src (Λ.p1 q)) := by rw [Λ.h41]
  simp only [holonomy, gaugeAct_def, h1, h2, h3, h4]
  group

/-- A character of the gauge group with the properties needed for the Wilson action:
normalized at `1` and invariant under conjugation.  For `G = SU(N)` this is
`χ(g) = Re tr(g) / N`. -/
structure CharacterData (G : Type u) [Group G] where
  χ : G → ℝ
  χ_one : χ 1 = 1
  χ_conj : ∀ g h : G, χ (h * g * h⁻¹) = χ g

variable (χ : CharacterData G)

/-- The Wilson action `S_β(U) = β ∑_p (1 - χ(U_p))`. -/
def wilsonAction (β : ℝ) (U : Λ.L → G) : ℝ :=
  β * ∑ q : Λ.Plaq, (1 - χ.χ (holonomy U q))

/-- **Gauge invariance of the Wilson action.**  This is the algebraic content behind
"the Wilson measure is a measure on the space of gauge orbits": the weight is constant
on orbits, so it pushes forward. -/
theorem wilsonAction_invariant (β : ℝ) (g : Λ.V → G) (U : Λ.L → G) :
    wilsonAction χ β (g • U) = wilsonAction χ β U := by
  unfold wilsonAction
  congr 1
  refine Finset.sum_congr rfl fun q _ => ?_
  rw [holonomy_gaugeAct, χ.χ_conj]

end LatticeGauge
end YangMillsInfra
