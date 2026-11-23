# TEMPLE-WEINSTEIN — rigorous eigenvalue LOWER bounds on the `0⁺⁺` glueball: the `E₀`-lever (PF vacuum) is real and new; the `E₁`-lever is the wall, and the naive `E₀`-shift formula is *false*

> **2026-11-23 · NEW SPECIES — Temple/Weinstein/Kato lower bounds.** 判词 (verdict). This is the **cleanest sign-CORRECT species**: a genuine *lower* bound on a spectral energy from a trial — the one direction the moment/GEVP routes (GAP-G2-edge) cannot produce. Two levers sharpen it: **(a) the `E₀`-lever** — PF gives the vacuum `φ₀` *exactly* (`PROOF-GAUSS-INVARIANT-PF`), `E₀=−(1/a)ln λ₀` is the free energy, so the vacuum is *closed*; **(b) the `E₁`-lever** — a glueball trial `ψ₁`. **Verdict: the proposed formula `E₁ ≥ E_trial − Var/(E_trial − E₀^lower)` is FALSE** (verified counterexample, §1); the correct Temple bound needs `ρ≤E₁` (the *next* level); and an airtight lemma (§3) shows mean+variance can **never** lower-bound an inf-support without a spectral edge. With PF closing `E₀`, `m₁=E₁−E₀` reduces *exactly* to an `E₂` lower bound `= I2 = H*`. **#1,2,4 PASS; #3,5 FAIL = the entire Clay content. No closure, no fabrication.**

## §1 The correct inequalities — and the false one
For `R=⟨ψ,Hψ⟩`, `Var=‖(H−R)ψ‖²`, spectrum `E₀≤E₁≤…`:
- **Temple (1928):** if `ρ ≤ E₁` and `E₀ ≤ R < ρ`, then `E₀ ≥ R − Var/(ρ−R)`. *Verified*: 0 violations / 3474 random measures. **The input is a lower bound on the SECOND level, not `E₀`.**
- **Weinstein:** `E₀ ≥ R − √Var` holds **only** under `R ≤ ½(E₀+E₁)` — 49/50 over-midpoint trials it fails. Same `E₁` input.
- **The proposed `E₀`-shift** `E₁ ≥ R − Var/(R−E₀^lower)` is **not a theorem**. Counterexample `{E₀,E₁,E₂}={0,1,2}`, `ψ=√.98 φ₁+√.02 φ₂` `(⊥φ₀)`: `R=1.02`, `Var=.0196`, RHS `=1.0008 > E₁=1`. The formula **over-claims**, worsening with admixture (`c₁=.2` ⇒ RHS `=1.71 ≫ 1`). No `E₀` substitution can replace the `ρ`-slot: Temple needs `R<ρ=E₀`, contradicting `R≥E₀`.

## §2 What the `E₀`-lever really buys
`m₁=E₁−E₀`, and an upper bound on `E₀` is free (Ritz; tighter, `E₀` is *exact* via the Perron eigenvalue of `T=e^{−aH}`). So `m₁≥λ ⟺ E₁ ≥ E₀+λ` — the whole problem is an **`E₁` lower bound**. PF removes the vacuum exactly (orthogonality `ψ₁⊥φ₀` rigorous), but on `(vac)^⊥` the bottom is `E₁` and the **second is `E₂`**: the regress restarts one rung up (`PROOF-EDGE-E2` §i). Temple-on-`(vac)^⊥` needs `ρ≤E₂`; Weinstein needs `R_⊥ ≤ ½(E₁+E₂)`. **Both = `I2`** (isolation `spec(H)|_{0⁺⁺}∩(m₁,(1+δ)m₁)=∅`). The vacuum lever is real but only relabels; it cannot supply the `E₁` direction.

## §3 The impossibility lemma (airtight)
*Mean+variance do not lower-bound an inf-support.* Fix `R=0`, `Var≈1`; put mass `ε=1/N²` at `−N`, `1−ε` at `s=−εt/(1−ε)`. Then `R=0`, `Var→1`, but `inf supp=−N→−∞`. **Verified.** So no trial-only `(R,Var)` bound lower-bounds `E₁`; every valid one (Temple/Weinstein) smuggles a **hard spectral edge** `E₂`. This is the energy-space face of GAP-G2-edge's "moments cannot pin the edge": a `ψ₁` sharp enough that `√Var_⊥<R_⊥−E₀−λ` *is* a near-eigenvector — the glueball itself. **Circular.**

## §4 Free-field acid test + uniformity
Free massless: `0⁺⁺` bottom is the two-photon cut, `E₁=0`, `m₁=0`; any valid bound gives `≤0` — acid test **passes**, but *vacuously* (as F2): absent `E₂` control it gives `≤0` for YM too. **Uniform `(V,a)` FAILS**: as `V→∞,a→0` the `(vac)^⊥` spectrum becomes a continuum, `Var_⊥→O(width²)`, so `R_⊥−√Var_⊥` drops `≤0` (F2-analog `c_V→1`); a uniform `ρ≤E₂` *is* uniform `I2` — unproven.

## Five-point scorecard
1. **Nonperturbative — PASS.** Exact `H`/transfer, PF vacuum, all `β`; no `1/N`.
2. **Sign-controlling `m₁²≥λ` — PASS (direction).** A genuine *lower* bound — the right sign the moment/GEVP routes lack. But it bounds the **gap** (bottom), blind to the atom `f₁²` (F2 §2; EDGE-E2).
3. **Non-circular — FAIL.** PF closes `E₀` only; `E₁` lower ⟺ `E₂` lower `= I2` (Temple/Weinstein), or ⟺ a known glueball vector (§3).
4. **Non-SUSY + fails free field — PASS (vacuous).** `E₁=0` ⇒ `≤0`; non-SUSY. But `≤0` for YM too absent `I2`.
5. **Uniform `≥λ` in `(V,a)` — FAIL.** `Var_⊥` grows against the continuum; uniform `ρ≤E₂` `= I2`.

## Exact failing estimate
```
Temple:   ρ≤E₁ , R<ρ  ⇒  E₀ ≥ R − Var/(ρ−R)           [valid; needs E₁]
on (vac)⊥: E₁ ≥ R_⊥ − Var_⊥/(ρ−R_⊥),  ρ≤E₂            ⇒  needs E₂ = I2
naive:    E₁ ≥ R − Var/(R−E₀^lower)  is FALSE          [{0,1,2} ψ=√.98φ₁+√.02φ₂: RHS=1.0008>1]
impossibility: (R,Var) fixed, inf supp → −∞             ⇒ no trial-only E₁ lower bound
exact missing:  E₂ ≥ (E₁+δ)² (I2)  OR  Var_⊥(ψ₁) ≤ (R_⊥−E₀−λ)² (a certified glueball trial)
```

**Sharpest line.** Temple/Weinstein is the species that finally *names the wall in one line*: a rigorous lower bound on `m₁` is **equivalent** to either a lower bound on the next level `E₂` (Temple, `=I2=H*`) or a variance bound on a glueball trial sharp enough that `√Var_⊥<R_⊥−E₀−λ` — which is a known glueball, the thing being bounded (`PROOF-F1SQ`-class circularity). The `E₀`-lever is the one genuinely new asset: PF supplies the vacuum *exactly* and the free energy closes `E₀`, so unlike F2's `G″` this species does not pay on the vacuum side. But it cannot manufacture the `E₁` side, because mean+variance provably cannot pin an inf-support without a spectral edge (§3), and the edge is `I2`. The proposed combined formula is not merely weak — it is **false**. **#1,2,4 pass; #3,#5 are `I2`/`H*`; the `E₀`-lever relabels the wall, it does not cross it.**

*Refs.* Temple, *Proc. London Math. Soc.* (1928); Weinstein, *Phys. Rev.* **40** (1932) 797; Kato, *J. Phys. Soc. Japan* **4** (1949) 334; Courant–Fischer. Internal: PROOF-EDGE-E2 §i (Temple regress one rung up); F2-KREIN-RUTMAN-QUANT (bottom-not-atom, `c_V→1`, `G″`); SUB2-I2-STRONG (`Cu` vs `Be^{−δt}`, the `O(u)` non-gap); GAP-G2-edge (moments cannot pin the edge; residue⟸level-gap); PROOF-GAUSS-INVARIANT-PF (Perron = vacuum, exact `φ₀`); PROOF-F1SQ (gap ⇏ residue; gapped free `f₁²=0`). Numerics verified: Temple 0/3474 violations; naive over-claims 5/5; impossibility lemma exact; Weinstein valid iff `R≤½(E₀+E₁)`. No fabricated bound.
