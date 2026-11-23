# LIEB–ROBINSON — information propagation is local *and* gap-blind: the LR commutator bound is satisfied by the gapless free field, so its rate `μ` is not a mass and cannot source `m_eff²(k→0)≥λ`; the only `LR⟹clustering` theorem re-inserts the gap. **No-go; 1 PASS, 2/3/4/5 FAIL.**

> **2026-11-23 · NEW SPECIES: Lieb–Robinson / information propagation.** Forge an LR commutator bound for the gauge transfer operator `T=e^{−aH_L}` (Kogut–Susskind), read a light cone with speed `v` and off-cone decay rate `μ`, and try to convert it into exponential clustering of the `0⁺⁺` correlator at rate `λ`, giving the target `m_eff²(k→0)=G″=E[G′]−Var(G)≥λ ⟺ Var(a)≤1/λ` (H1). **判词 (verdict): NOT ACHIEVED, a clean no-go.** LR is a genuine nonperturbative, non-SUSY dynamical species (#1 PASS), but it is **universal for finite-range local dynamics**: it holds with finite `μ>0` for the *gapless* free massless field. So `m_eff:=μ` is a **spurious gap** on the free field (#4 FAILS — the species is satisfied precisely where it must vanish); the only rigorous `LR⟹equal-time clustering` theorem (Nachtergaele–Sims) *assumes* the gap it is asked to produce (#3 FAIL); and LR upper-bounds a commutator of *time-evolved* observables, never a static curvature (#2,#5 FAIL). It re-encodes H1 as "propagation."

## 1. The forge — LR for the gauge transfer operator
`H_L=Σ_ℓ g²E_ℓ² + (1/g²)Σ_p tr(1−U_p)` on the Gauss-law-constrained slice; `T=e^{−aH_L}`. Electric and magnetic terms are nearest-neighbour / plaquette, i.e. **finite-range** with bounded interaction norm `‖Φ‖_μ<∞`, so the Lieb–Robinson theorem (1972; Nachtergaele–Sims) applies at **every β, every g**, gapped or not:
```
‖[A_X(t), B_Y]‖ ≤ 2‖A‖‖B‖ e^{−μ(dist(X,Y) − v|t|)},   v = c_μ‖Φ‖_μ.
```
A genuine light cone: signal speed `v`, off-cone decay rate `μ`. The fork: is `μ` the mass? (Infinite local dim for compact `G` is handled by a Peter–Weyl/electric cutoff that the bound descends in — UV bookkeeping, not the IR issue below.)

## 2. Acid test — `μ ≠ m_eff`; the free field is LR
The free massless lattice field is a finite-range local system; harmonic LR bounds (Nachtergaele–Raz–Schlein–Sims) give it a **finite** `μ_free>0`, `v_free` — yet its gap is exactly zero and its static `0⁺⁺` correlator is **power-law** `C(r)∼r^{−2}` (d=4), so `Var(a₀)=Σ_r C(r)=∞` and `G″=0`. Both read-offs die: (i) `m_eff:=μ` returns `μ_free>0` on a **gapless** field — a spurious gap, the exact negation of the acid test; (ii) the bound constrains the **commutator of time-evolved** operators `[A_x(t),B_y]` (off-cone-exponentially-small for the free field), *not* the **static connected correlator** `⟨O_xO_y⟩_c` — the power-law object `Var(a₀)` that defines the target. LR is blind to `C(r)`.

## 3. The `LR⟹clustering` bridge *is* the gap
The one theorem turning LR into equal-time clustering is Nachtergaele–Sims / Hastings: for a ground state with spectral gap `γ>0`, `|⟨O_xO_y⟩_c| ≤ C e^{−γ r/(2v)}`, i.e. `ξ ≤ v/γ` — it is **conditioned on `γ>0`**. There is no converse: LR (or even exponential clustering) does **not** imply a gap — and 定理 B (mass-spectrum REPORT) gives clustering over a *continuous* bottom with `f₁²=0`, so clustering would not even reach the atom. A proof of `|C(r)|≤e^{−λr}` presupposes `λ`; then `Var(a₀)≤C_dλ^{−d} ⇒ G″≥λ` merely restates it. The real-time→Euclidean-time continuation *is* the spectral-gap input. Also invisible: the Savvidy tachyon `λ_min=−g|B|` is a **stiffness (Hessian-sign)** defect, and LR bounds *speed*, not stiffness — the same blindness as the MERA causal cone (HOLOGRAPHIC-TN).

```
LR dynamics :  ‖[A_x(t),B_y]‖ ≤ 2‖A‖‖B‖ e^{−μ(r−v|t|)}     (∀ finite-range systems, free field too: μ_free>0)
target      :  |⟨O_xO_y⟩_c| ≤ e^{−λr}  ⇒  Var(a₀)=Σ_rC(r) ≤ C_dλ^{−d} ⇒ G″=E[G′]−Var(G)≥λ
bridge      :  LR⇒clustering REQUIRES γ>0 (Nachtergaele–Sims, ξ≤v/γ) ⇒ λ=γ = the gap      [circular]
free field  :  μ_free>0 but m=0, C(r)~r^{−2}, Var(a₀)=∞ ⇒ read-off μ=m is SPURIOUS (acid test fails)
```

## 4. Five-point scorecard
1. **Nonperturbative — PASS.** Exact LR bound for all β,g; no loop/`1/N` expansion; genuinely a new (dynamical) species.
2. **Sign-controlling `m_eff²(k→0)≥λ` — FAIL.** LR gives an *upper* bound on a commutator; `μ` is not `m_eff`, and no lower bound on the static curvature `G″` follows.
3. **Non-circular — FAIL.** `LR⟹clustering` assumes `γ>0`; `Σ_rC(r)≤1/λ` presupposes `|C|≤e^{−λr}` = the gap. Reduces to H1.
4. **Non-SUSY + fails free field — FAIL (decisive).** Free massless has finite LR `μ` but `m=0`: the species is *satisfied* by the gapless field. Read dynamically it emits a spurious gap; read statically it is silent. It does not cleanly vanish on the free field.
5. **Uniform `(L,a)` — FAIL.** `v,μ` are uniform bounds but on the wrong (commutator/dynamical) object; `Var(a₀)` at `k→0` is uncontrolled; no `λ`.

## Exact missing step
```
∃λ>0: |⟨O_x O_y⟩_c| ≤ Ce^{−λ|x−y|}  (equal-time, k→0)  ⟺  Var(a₀) ≤ 1/λ  is H1.
LR supplies only real-time commutator decay; the imaginary-time continuation needs γ>0.
Missing object: an equal-time clustering bound, i.e. the spectral gap itself. = H1.
```

**Sharpest line.** Information propagation has the right *shape* for a mass — a light cone with a speed and a decay length — but Lieb–Robinson's cone is a property of **locality**, shared by every finite-range lattice theory, the gapless free field included. The LR rate `μ` bounds signal *speed*, not *stiffness*: finite for a massless field, so it cannot lower-bound `G″=S−D`. Converting it to static clustering requires already holding the gap (Nachtergaele–Sims), and seeing the Savvidy tachyon requires a Hessian sign LR lacks. The species passes #1 and dies on #2/#3/#4/#5: H1 in propagation clothing, one analytic continuation short.

*Refs.* Lieb–Robinson, *CMP* **28** (1972) 251; Nachtergaele–Sims, "Lieb–Robinson bounds and the exponential clustering theorem," *CMP* **265** (2006) 119; Hastings, "Spectral gap and exponential decay of correlations," *CMP* **265** (2006) 781; Nachtergaele–Raz–Schlein–Sims, "Lieb–Robinson bounds for harmonic and anharmonic lattice systems," *CMP* **286** (2009) 1073; Bratteli–Robinson, *Operator Algebras and Quantum Statistical Mechanics 2*; Kogut–Susskind, *PRD* **12** (1975) 395; Cubitt–Pérez-García–Wolf, "Undecidability of the spectral gap," *Nature* **528** (2015) 207. Internal: PROOF-CONFINEMENT-RANDOM-SURFACE (area law ⇒ string tension but ⇏ `0⁺⁺` gap; small-λ only); G1-CONTINUUM-OS (OS gives the object `(H,T_t,ν)`, neither gap nor atom; 定理 B continuous bottom); C8-SPECIES-NOGO (only `E[G′]−Var(G)=G″` survives; torsion kills association classes); ABC-NEW-INEQ §1 (Ward–Stein `G″·Var(a)=1` iff Gaussian); PROOF-JACOBIAN-SAVES (◆: `G″≥λ⟺Var(a)≤1/λ`; `−g|B|`); BHH-TRANSLATION (same `E[G′]−Var(G)=G″`=H1 endpoint via a different device); HOLOGRAPHIC-TN (causal cone bounds speed not stiffness; same blindness); G1-EXISTENCE-THEOREM; CONDITIONAL-LATTICE-GAP.
