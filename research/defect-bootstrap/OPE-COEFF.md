# Bulk–defect OPE data for pure-gauge `Z_N`

> CPB-attack probe (round 4). Verdict: both quantities are **dynamical**, not symmetry-fixed —
> which is exactly why the defect-4pt → bulk-gap link is indirect.

## One-point function

Defect conformal symmetry fixes `⟨tr F²(x)⟩_def = a_R/ρ⁴` (`ρ` = transverse distance to the
line; `Δ_{trF²}=4`). `a_R` is a **dynamical** defect-EFT datum, not symmetry-fixed. Free/leading
order: the straight line in rep `R` produces the non-abelian Coulomb field
`F_{0i} ~ C₂(R)^{1/2} g/(4πρ²)`, so `a_R ∝ C₂(R) g²/(16π²)` (free `U(1)`: `q²/16π²`). At the
interacting marginal point `a_R` is a function of the coupling — **unknown**.

## "`tr F² ↔ D`" is not a direct OPE

`tr F²` is spin-0 **singlet**; the displacement `D^a` is spin-1 **adjoint**. Spin and gauge
mismatch ⇒ **no** leading-order `tr F² ↔ D` term. The physical coupling is the **defect–defect**
OPE `D^a D^b ⊃ δ^{ab} c · tr F²` (the bulk channel of `DEFECT4PT-ZN.md` / `BULK-CHANNEL.md`).
Its coefficient `c` is dynamical; free theory gives it by Wick (two transverse-field-strength
inserts fusing to the singlet), `O(g²)`, **not** symmetry-fixed.

## Honest boundary

Both quantities are dynamical CFT data, not fixed by symmetry — precisely why the defect-4pt →
bulk-gap link is indirect. The interacting values were not computed; the free-theory leading
coefficients rely on standard `C₂(R)` scaling and a convention factor not pinned exactly.
