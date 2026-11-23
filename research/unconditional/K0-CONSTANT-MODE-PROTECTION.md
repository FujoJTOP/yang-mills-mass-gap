# K0-CONSTANT-MODE-PROTECTION — is the uniform (k=0) mode shielded from the tachyon? Three steelmen, one exact step: Var(∇S|0)

> **2026-11-23 · 判词 (verdict).** The constant (k=0) mode is **not protected** — it is the *seat* of the homogeneous chromomagnetic tachyon, and the measured `G″(k=0)≈+0.4` is carried by the geometric/entropic term `W_geo+X`, **not** by any symmetry. The three candidates collapse onto already-banked partial results and die at **one exact step**: an upper bound on the constant-mode-conditioned force variance `Var(Pᵀ∇S|0)` (the harmonic-mean stiffness `≥λ`) = **H1 verbatim**. **(a) compactness** = the β=0 Prékopa theorem `Φ″⪰N/12` (entropic, O(1), orthogonal to the dynamical gap; the `+0.4`'s positivity *is* `W_geo`, dominating `Δ_dyn` is H1). **(b) Ward** = banked PF positivity `≥0` — the tachyon is transverse-magnetic and *survives the gauge quotient* (W4: 7.1% negative), so no Ward identity nullifies it. **(c) 1-D QM** = the no-invariant-subspace wall (k0-CONSTANT-MODE-TRANSFER) + the wrong-direction reduction (k0-EXACT-ONEMODE). Best score **1–2/5**; not genuine. The constant mode is H1's sharpest face, not its escape.

## The sign, honestly stated
The `E[H]−Var` proxy at the softest `k→0` is `≈+0.043−0.208 = −0.165 < 0` (VAR-FINAL; W6); the *full* `G″=Δ_dyn+W_geo+X` is `≈+0.4–0.51` (COMPUTE-GPP; k0-CONSTRAINT-EFFPOT). The gap-shaped positivity sits **entirely in `W_geo+X`** overcoming the negative Savvidy core, so any "protection" must prove `W_geo+X > −Δ_dyn` at `k=0` — exactly (◆)=H1. None of (a),(b),(c) supplies a *mechanism* for that domination; each gives a bound strictly weaker than needed.

## (a) Compactness ⇒ Haar-convexity ⇒ bounded support
*Steelman.* `G` compact ⇒ the constant field is bounded ⇒ the per-link Haar factor `e^{−Φ}` is log-concave with `Φ″⪰N/12` (SU2 1/6, SU3 1/4). Prékopa pushes this through the linear marginal: `V_eff″(0)≥N/12` for **every** mode including `k→0`, uniform in `L` (k0-CONSTRAINT-EFFPOT §2; E1-FRG §3). Bounded support ⇒ `χ=Var(a₀)<∞`, versus the free field's `χ=∞` (D3 §3) — a genuine acid-test discriminator.
*Attack.* (i) `N/12` is the β=0 **entropic** theorem, O(1), living in `W_geo`, not the dynamical curvature. At β=2.5 the cactus action `V₀(a)=Σ_□(1/2N)tr[a_μ,a_ν]²` is **quartic and exactly flat** on the `4(N−1)`-dim Cartan valley (`V₀″(0)=0`) and **Savvidy-concave** off it; `S_W` is nonconvex so Prékopa log-concavity of the marginal **fails** (D3 §4). (ii) Bounded support ⇒ finite `χ` does **not** imply `G″(0)>0`: a symmetric marginal on a bounded interval can have a central *dip*, `−(log p)″(0)<0` — curvature sign ≠ variance magnitude (Jensen `E[1/ψ″]≥1/G″(0)`). (iii) The only bound it supplies is `Var(a₀)≤(√V c)²/4=O(V)` — no uniform λ, useless in `L`. Verdict: passes the acid *form* at β=0, fails #2/#5.

## (b) Ward identity / gauge-invariance of the uniform mode
*Steelman.* The uniform mode is gauge-covariantly special: constant `A₀` has `F=0` but a nontrivial Wilson line, and Ward/Elitzur gives `E[a]=0` plus a conserved current, suggesting a protected flat direction (the holonomy).
*Attack.* The negative Savvidy direction is **transverse-magnetic and survives the gauge quotient** (CLAUDE.md front: W4 measured 7.1% negative *after* projection) — a Ward identity cannot nullify a *physical* negative direction. PF/reflection-positivity certifies only `f₁²≥0` and `Ĥ≥0` — "**positivity, never separation**" (D3 §5): the free field is also RP with `f₁²=0`. The holonomy reading is already a no-go: `V″(Ω)=c/L²` is a one-loop **free-gluon electric** mass, the **wrong operator** (electric vs magnetic), in the center/1-form channel **provably blind** to the neutral `0⁺⁺` that carries the gap (k0-HOLONOMY-EFFPOT; PROOF-1FORM-ANOMALY). Verdict: banked ≥0 floor; #2/#3/#5 fail.

## (c) Constant mode = 1-D QM, bounded spectrum ⇒ discrete gap
*Steelman.* Integrate out all non-constant modes ⇒ a 1-D Schrödinger operator for `a₀(τ)` on a compact target ⇒ discrete spectrum with a gap ⇒ `G″(k=0)>0`.
*Attack.* Three walls, all banked (k0-CONSTANT-MODE-TRANSFER): **(A)** `[T,P₀]≠0` — the plaquette is a **momentum convolution**; at cubic order it maps `k=0` into every `k`, so `P₀TP₀` is no sub-operator. **(B)** The mode is `4(N²−1)`-dimensional, the tree potential quartic with a **flat Cartan valley** — no harmonic gap. **(C)** Any curvature from a genuine elimination *is* `G″(0)`, one-loop Savvidy-indefinite. Compact-space discreteness is the gap of the **reduced** theory, and `gap_red≥gap_true` (k0-EXACT-ONEMODE) — the **wrong direction**; deficit `D=Var(∇S|0)=H1`. Verdict: #2/#3/#5 fail; #4 reverses.

## Five-point scorecard (shared)
1. **Nonperturbative — PASS only (a)@β=0** (Prékopa+`Φ″⪰N/12`); (b),(c) tree/one-loop FAIL.
2. **Sign control `G″(k=0)>0` — FAIL all.** `+N/12` entropic @β=0; @2.5 bare `k=0` `−0.3`, `E[H]−Var≈−0.165`, `+0.4` is `W`-carried and unproven to dominate.
3. **Non-circular — FAIL at closing:** residual step is `Var(Pᵀ∇S|0)≤…` = H1.
4. **Non-SUSY + fails free — PASS (a)@β=0** (free `Φ=0⇒0`, `χ=∞`); (b) FAIL; (c) FAIL (reduced gap≥true).
5. **Uniform `≥λ` — FAIL all:** compactness `O(V)` bound; flat valley; reduced gap *upper*-bounds.

## Exact stuck step
```
Need (candidate-independent):
  W_geo + X > −Δ_dyn  at k=0
  ⇔ Var(Pᵀ∇S|0) − E[H|0] < W_geo
  ⇔ E_p[1/G″⁺] ≤ 1/λ  (harmonic-mean stiffness ≥λ)
  ⇔ Var(a₀) ≤ 1/λ          = (◆) = H1.
```
Every protection mechanism — compactness, Ward, 1-D QM — terminates at the **same** upper bound on the constant-mode-conditioned force variance. The constant mode is **not** the direction where a symmetry saves us; it is the direction where the Savvidy variance is *largest* because the tachyon is *homogeneous* (k=0).

**Sharpest line.** The constant mode is the tachyon's **home**, not its sanctuary: its off-diagonal part carries the `−0.3` bare Hessian, its symmetry-protected Cartan/holonomy part is the center modulus **provably blind** to the neutral flux sector where the gap lives, and its measured `+0.4` is **entropy** (`W_geo`, Prékopa), not symmetry. Compactness re-derives the banked β=0 `N/12`, the Ward identity re-derives banked PF positivity, the 1-D QM re-derives the wrong-direction reduction — and all three stop at `Var(G|0)≤…`: H1 verbatim. The missing machine is still the one that makes `W_geo` dominate the Savvidy variance at the constant mode without assuming it.

*Refs.* k0-CONSTRAINT-EFFPOT (`G″=Δ_dyn+W_geo+X`, β=0 `N/12`); k0-EXACT-ONEMODE (`gap_red≥gap_true`, `D=H1`); k0-CONSTANT-MODE-TRANSFER (no invariant subspace; flat valley; curvature=target); k0-HOLONOMY-EFFPOT (electric `1/L²`, center-blind); D3-SUSCEPTIBILITY-VARIATIONAL (BL `Var≤E[1/ψ″]`, PF positivity≠separation); E1-FRG (convexity `≥0`, `N/12` invariant); FOURIER-UNIFORMITY (`S=1/Var`); PROOF-COARSE-POINCARE (H1, `H_ff` indefinite, min −0.3); VAR-FINAL / W6-BE-at-2.5 (`−0.165`, 7.1% negative after quotient); PROOF-GAUSS-INVARIANT-PF (RP⟹atom≥0); PROOF-1FORM-ANOMALY (center channel). External: Prékopa 1973; Brascamp–Lieb 1976; Gross–Pisarski–Yaffe, RMP **53** (1981) 43; Savvidy, PLB **71** (1977) 133; Shen–Zhu–Zhu, CMP **400** (2023) 805.
