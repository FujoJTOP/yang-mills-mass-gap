# PROOF-JACOBIAN-SAVES — the Haar term's positivity is provable; its dominance over Savvidy is the gap

> **Verdict.** The geometric term's **sign** is now a theorem. With `Φ=−log J=Σ_{i<j}g(λ_i−λ_j)`, `g(d)=−2log(sin(d/2)/(d/2))`, the pair function is convex (`g″≥1/6` on `|d|<2π`); so `Φ` is convex, `e^{−Φ}dX` log-concave, and by Prékopa every linear marginal `a=⟨v,log U⟩` has convex free energy: `W_geo=E[Φ″]−Var(Φ′)≥0` for **every mode, every k including k→0**, uniformly in L. This closes the sign at β=0 (`G″=W_geo`) and at weak UV coupling (`G″≥δ_W·μ(typ)−Ce^{−c/g²}>0`, A1). It does **not** close at β=2.5: the exact split is `G″=Δ_dyn+W_geo+X`, `Δ_dyn=E[H]−Var(G)`, `X=−2Cov(G,Φ′)`. The measured `W_full=W_geo+X≈+0.22–0.31` beats `−Δ_dyn≈+0.02`, but certifying `W>Var−E[H]` needs an upper bound on `Var(G|0)` — the induced-mass/H1 object. The −0.165 "DEAD-lean" is a proxy that *dropped exactly +W*; the exact marginal on its mode is **+0.51**.

## (1) W, explicitly
Flat chart `X_ℓ=log U_ℓ∈su(N)`, `U_ℓ=e^{iX_ℓ}`, Haar `dU=J dX`, `Φ=−log J`. With flat Hessian `H_vv=d²S/ds²`, `G_v=dS/ds` along `X→X+sV`, the identity is exact:
`G″(0)=E[H_vv+Φ″|a=0]−Var(G_v+Φ′|a=0)=Δ_dyn+W_geo+X`,
`Δ_dyn=E[H]−Var(G)`, `W_geo=E[Φ″]−Var(Φ′)`, `X=−2Cov(G,Φ′)`.
`W_geo` is the curvature of the **pure-Haar** marginal; `X` is the cross-fluctuation of the action force with the Jacobian. W7's `W` is `W_full=W_geo+X`.

## (2) Bounds; is `W>Var−E[H]` provable?
**Lemma B (quantitative, pointwise).** SU(2) has `Φ(θ)=−2log(sin(θ/2)/(θ/2))`, `Φ″(θ)=½csc²(θ/2)−2/θ²≥1/6`, `Φ′/θ≥1/6` on `(0,2π)`; radial+transverse decomposition gives `Φ″[v]≥(1/6)|v|²` per link (exact for SU2), so at eigenvalue order `E[Φ″|0]≥(N/12)‖v‖²` (SU2 `1/6`, SU3 `1/4`; measured `0.34`,`0.46`). Crucial: `Φ″⪰0` is pointwise, so `E[Φ″]≥0` under *any* measure, any β.
**Prékopa (sign).** `Φ` convex ⇒ `e^{−Φ}dX` log-concave ⇒ pushed `a` log-concave ⇒ `V_eff` convex ⇒ `W_geo≥0`, all `k` including `k→0`. At β=0 this is the whole curvature (SU2: `V_eff″=0.248` vs `E[Φ″]−Var(Φ′)=0.246`, exact closure).
**The failing piece — `Var` from above.** The only rigorous handle on `Var(G|0)` is Lemma A (`≤β²N_□/λ_f`), but `λ_f>0` at β=2.5 is unproven and `λ_f→0` at `k→0`, so it is vacuous. Thus `W>Var−E[H]` is **not** provable at β=2.5. The exact missing inequality:
**(◆)** `Var_μ(⟨v,∇S⟩|a=0) − E_μ[H_vv|a=0] < W_full`, uniform in (L,U′). Since `G″≈1/Var(a)` to 2–5% (rerun 0.281 vs 1/Var 0.272), **(◆) ⟺ Var(a)≤1/δ ⟺ m_eff²(k→0)≥δ** — H1 on the coarse mode. Proved: `W_geo≥0` (strongest at k→0); missing: control of the dynamical force variance.

## (3) Reconciliation — three objects, not two values
The −0.165 and the +0.28/+0.43 are **different objects and regimes**, not conflicting measurements of one `G″`.
- **−0.165 (VAR-FINAL; SU(3) β=2.5, L=4):** the `b=2` block-uniform **softest blocked mode**, computed `λ_E−Var_total=+0.043−0.208`, with `λ_E` a 3-config `Ĥ` softest eigenvalue, `Var_total` the *total* force variance (`R²=0.004`), and **no +W**. A proxy, not `G″`.
- **+0.51 (W7 Part B; same mode, same β):** exact histogram, full identity: `G″=+0.51±0.08` (bracket `−0.03`, `W_full≈+0.54`, trimmed +0.56). The exact marginal on the proxy's own mode is **positive**; the sign is the +W the proxy dropped (and conditioning, `E[H|0]=0.22>λ_E=0.04`).
- **+0.28/+0.43 (W7-RERUN; SU2 β=2.5, SU3 β=5.7):** finite-`k` (`π/2,π/3`) plane waves at each group's **crossover**: `E[H]−Var≈0` (equipartition, chart-robust: ray `2.46/2.44`, flat `1.52/1.54`), `E[Φ″]=+0.34/+0.46` volume- and mode-independent, `W=+0.22–0.31` carries `G″`.

**Which regime.** SU3 β=2.5 is deep strong coupling (Haar-dominated) — so the *largest* `G″` (+0.51) sits on the **softest, most k→0-like** mode: the geometric term dominates most where the measure is least action-like. SU2 β=2.5 and SU3 β=5.7 are the crossovers. No run reaches true `k→0` (L=4,6 give `k=π/2,π/3`); the coarse-most available mode is the `b=2` block-uniform one (exact +0.51/+0.27), and its `1/Var` is volume-independent (0.272 at L=4 *and* L=6) — evidence, not proof, of a stable limit.

**Sharpest line.** The Haar-Jacobian curvature is **provably non-negative on every mode including k→0** (Φ convex ⇒ Prékopa ⇒ `W_geo≥0`); positivity is a theorem at β=0 and weak UV. At β=2.5 the measured `W_full` beats the Savvidy deficit by 0.2–0.5, but proving `W>Var−E[H]` is exactly a bound on `Var(G|0)` — H1 in curvature language. The −0.165 subtracted the full force variance and added no Jacobian; its own mode, measured exactly, gives +0.51. Jacobian saves the sign numerically and at small β; at the physical coupling the rescue is one variance bound short.

*Refs.* W7-Gpp0 (Part B +0.51); W7-RERUN (identity, β=0 closure 0.248/0.246); VAR-FINAL (−0.165 proxy); PROVE-Gpp-positive (Lemma A/B, (★)); A1-twosided-UV (UV positivity); PROOF-COARSE-POINCARE (H1 at k→0); W6-BE-at-2.5; MARGINAL (softest `PᵀHP` −0.3); Prékopa 1971; Brascamp–Lieb 1976; Savvidy 1977; Shen–Zhu–Zhu, CMP 400 (2023), arXiv:2204.12737.
