# PROOF-TRANSPORT-RG — a residue-carrying block-RG step: does the off-diagonal residue escape the IR-H1 trap?

> 2026-11-23 · Attempted construction of one rigorous block-RG step mapping `β=2.5 → β_eff<1/48` that *carries the residue* `f₁²`, and of the transport inequality `f₁²(2.5) ≥ c·f₁²(β_eff)`. **Verdict: the residue-RG escapes the *specific* IR-H1 trap — the residue is mass-free and off-diagonal, so Savvidy negative modes act as mass/vertex shifts P3's weight `r₀` cannot see (**b** is a genuine mechanism) — but it does not cross the W2 wall. One factor-`b` step leaves `β_eff≈2.5`; the iterated flow to `β_eff<1/48` runs through the crossover, where the blocked measure carries an `O(1)` Savvidy deformation, so P1/P2's small-activity expansion and P3's constructive bound both stop (**a**: the atom persists only as the open `I2`; **c**: `c` unproven). The transport reduces to a *non-cancellation / operator-stability* bound on a signed-measure amplitude, which `ν≥0` gives only as an *upper* bound (B1). No fabricated proof.**

## 1. The RG map
Fine `U` on `T⁴_L`, coarse `U′` on `T⁴_{L/b}`. Fixed local covariant low-pass `K(U′,ξ)≥0`, `∫K=1`, reconstruction `U=Θ(U′,ξ)`. Blocked measure and operator
`ν_b(U′)=∫dξ K(U′,ξ)e^{−S_W[Θ(U′,ξ)]}/Z`,  `O^{(b)}=E_{ν_b}[O\,|\,U′]` (`L²`-projection of the plaquette), `β_eff(b,β)` by matching the renormalized `2βN_p` (Migdal–Kadanoff at leading order). Residue `f₁²(β)=|⟨Ω_b|O^{(b)}|G_b⟩|²`, `|G_b⟩` the lightest `0⁺⁺` state.

## 2. Exact invariance ⇒ transport *is* a two-measure comparison
Blocking is a change of variables on *one* theory, so by the tower property `C(t)=⟨O(t)O(0)⟩=⟨O^{(b)}(t)O^{(b)}(0)⟩_{ν_b}`, the KL measure `ν_O` is identical, and `f₁²` is literally the same number at every `b`: `f₁²(2.5)=f₁²_{ν_b}(β_eff)`. So the inequality cannot be won by invariance. It must compare **two measures at the same nominal coupling** — blocked `ν_b` versus *bare Wilson* at `β_eff` (where P3 applies). That comparison is the whole content.

## 3. What the residue is blind to — the real escape from IR-H1
The gap's content is the `k→0` transverse stiffness `G″(0)` (PROOF-IR-H1). The residue is **mass-free**: P3/P1P2 give `r₀=c₁(β)u^p` with no `m₁`; decorations are "a mass shift plus a positive `1+O(u)` factor, never a cancellation" (P1P2 §4). So a Savvidy mass-shift `δm_G` never enters the residue weight. Formally `f₁²=|⟨Ω|O|G⟩|²` is a *transition amplitude between distinct states*, so the vacuum's own curvature—the diagonal object failing at `k→0`—does not appear. This is the rigorous sense of **(b): negative modes resummed into `m_G` (mass) and `Z_O` (vertex), not disabled**, and `r₀` is blind to them.

## 4. The obstruction — `β_eff` is not `<1/48`, and `ΔS` is `O(1)`
`β=2.5` is `g²=2N/β≈1.6–2.4`; the radius is `β<1/48 ⇔ g²>96N` (≈200–300). **One** factor-2 step gives `β_eff≈β` (slow UV running; PROOF-H2-GLOBAL route (i)); reaching `g²~10²` spans ~`50` e-foldings ⇒ **tens** of steps, *through the crossover*. There `ν_b=Wilson_{β_eff}+ΔS` with `ΔS` the accumulated Savvidy remainder of `O(1)` activity — **not** `O(u(β_eff))` (W2: "blocking re-imports the Savvidy negativity"; `β_eff→0` does not shrink `ΔS`, which encodes the blocked-out `β=2.5` physics). P1/P2's KP condition `u·K<1` sees `ν_b`'s activity, i.e. `O(1)`.

## 5. The exact missing inequality
For `ν_b` at `β_eff`, uniform in `V` and over the iterated flow:
```
(M1)  |⟨Ω_b|O^{(b)}|G_b⟩| ≥ z₀·|A₁(β_eff)| = z₀·√(r₀(β_eff)),   z₀ = z₀(b,N) > 0 ,
(M2)  m_G(ν_b) < m_excited(ν_b) − δ ,   δ > 0 .
```
`(M1)` is operator-stability/non-cancellation (`O`'s leading-tube amplitude survives `ΔS` with fixed fraction `z₀`; then `c=z₀²`); `(M2)` is the two-mass isolation `I2` (PROOF-EDGE-E2). First-order in `ΔS`, `δ⟨Ω|O|G⟩=Σ_n V_{Gn}⟨Ω|O|n⟩/(E_G−E_n)`, an uncontrolled signed sum that can *in principle* cancel the `O(√r₀)` amplitude. Neither is provable: `(M1)` is a *lower* bound on a matrix element, and every `ν≥0`/RP total is Schur-concave ⇒ *upper* (EDGE-f1²; B1 moment lemma), while the signed `ν_b` admits no constructive lower bound; `(M2)` is open `I2`. With `(M1)∧(M2)`: `f₁²(2.5)≥z₀²r₀(β_eff)>0`.

## 6. Escape IR-H1, fall into W2⊗B1
The residue-RG **does escape IR-H1 in its specific form**: it never needs `G″(0)>0`. But it does **not** escape W2 (upstream): the controlled flow to `β_eff<1/48` is the same wall the gap-RG hit. And the residue adds a *second* wall: `(M1)` is a lower bound positivity forbids and the signed measure can't supply — B1's sign obstruction, not the soft mode. The attack converts one hard input (`H1`) into two (`W2`-flow + `(M1)/(M2)`); off-diagonality helps only against the *diagonal* failure, while the **magnitude** stays usc and unprotected.

**Sharpest line.** The residue-RG's insight is exact: `r₀` is mass-free, so Savvidy modes are genuinely resummed into `m_G`/`Z_O` and never flip the atom's square. But "resummed" needs small activity, and the `β=2.5→1/48` blocking is `O(1)`-deformed exactly where the residue's *amplitude*—not its sign—becomes a signed cancellation problem. Off-diagonality beats the soft-mode curvature; it does not beat the crossover.

*Refs.* PROOF-P3; PROOF-P1P2; PROOF-F1SQ; RESIDUE-CONDITIONAL (`H*`); PROOF-EDGE-E2 (`I1/I2`); PROOF-IR-H1; PROOF-H2-GLOBAL; REPAIR-2SCALE; INEQ-R7-rg-invariant-decay; B1-nonperturbative-residue; EDGE-f1sq-attack; W2-positive-scheme; W2-transmuted-gap; BLOCKING-DESIGN; Shen–Zhu–Zhu, CMP 400 (2023) 805.
