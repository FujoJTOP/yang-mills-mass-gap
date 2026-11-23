# OTHER-GPP-KZERO — proving `G″(0)>0` at `k→0` directly: +0.51 is the right sign, the Haar term is provably `≥0` on every mode, but the `k→0` lower bound is exactly one variance short

> **2026-11-23 · Verdict (判词).** The measured +0.51 is **indicative, not proof**. It is the *exact* marginal curvature on the softest accessible mode (SU(3) β=2.5, L=4, `b=2` block-uniform) and it **reverses** VAR-FINAL's −0.165: that proxy subtracted the *full* force variance and added **no** `+W`, whereas the exact identity `G″=Δ_dyn+W_geo+X` gives **+0.51±0.08** on the proxy's own mode (`E[H]−Var≈−0.03`, `W_full≈+0.54`). The carrier is the one object provably non-negative on **every** mode including `k→0` — the Haar/Jacobian term `W_geo≥0` (Φ convex ⇒ Prékopa) — volume- and mode-stable (`1/Var=0.272` at L=4 *and* L=6). But *proving* the `k→0` sign is **not reachable**: the certified half is only `W_geo≥0` (value-free), and the sign needs `Var(G|0)−E[H|0] < W_full` uniformly in `(L,U′)` — *verbatim* H1. Both attacks fail: **(a)** `W_geo`'s "full Prékopa strength beyond `N/12`" is not a per-config object; **(b)** the exact Savvidy negativity is physical and pointwise unbounded, suppressed only by `e^{−S}`. **+0.51 is a sign probe, not H1.**

## The competition, exactly
Flat chart `X=log U∈su(N)`, `a=⟨v,X⟩`, `V=S_W+Φ`, `Φ=−log J`. The identity is exact:
```
G″(0)=Δ_dyn+W_geo+X,  Δ_dyn=E[H]−Var(G),  W_geo=E[Φ″]−Var(Φ′),  X=−2Cov(G,Φ′).
```
`W_geo≥0` for **every** mode, every k, uniform in L (Prékopa; pointwise `Φ″⪰N/12`, Lemma B). At β=0 this is the whole curvature (SU2 0.248 vs 0.246, exact closure). At β=2.5 `Δ_dyn≈0` (equipartition, `E[H]≈Var(G)≈2`) and the *measured* `W_geo≈+0.32` carries the sign. Since `G″≈1/Var(a)` to 2–5%, `G″>0 ⟺ Var(a)≤1/δ ⟺` induced mass.

## Why +0.51 is indicative
- **Same mode, same β:** −0.165 (proxy: full variance, no W) vs +0.51 (exact: conditioned, +W). The gap between them *is* the geometric term; `E[H|0]=0.22 > λ_E=0.04`.
- **Stability:** `W_geo=+0.318/+0.319`, `1/Var=0.272/0.272` across L=4→6; `E[Φ″]=0.34/0.46` volume/mode-independent. SU3 β=5.7 at k=2π/L: `+0.425` (L4)/`+0.329` (L6), W-dominated (66%/94%), `W_geo` identical on both volumes.
- **Mechanism:** `corr(a,Φ′)=+0.85` — the Jacobian tracks the amplitude, so the geometric curvature *is* the transverse inverse propagator. Positive, W-carried, volume-stable: consistent with a positive `k→0` limit.
- **But:** every run stops at k=π/2 or π/3 (L=4,6); the L=4→6 drop (0.43→0.33) mixes k-decrease with volume change. A sign probe, not a limit.

## Can it be proven? Exact failure
```
(◆)  Var_μ(⟨v,∇S⟩|a=0) − E_μ[H_vv|a=0] < W_full,   uniform in (L,U′)
     ⟺ Var_μ(a)≤1/δ ⟺ m_eff²(k→0)≥δ  =  H1.
```
- **(a) Push `W_geo`'s Prékopa beyond `N/12` — fails.** `W_geo`'s gauge-invariant content is the *constant floor* `N/12`; the surplus convexity (`g″→∞`) lives at the log cut (`d→2π`, links near `U≈z∈Z_N`, so `F≈0`), **anti-correlated** with `∝|F|`. No gauge-invariant norm makes it dominate `H_ff`: spectral ratio `ρ≈3–16≫1` at β=2.5. Prékopa is a measure-level tool; a per-config Loewner bound is a category error. Scaling is the wall: `W_geo~N/12` (β-free) vs `−1.12β` (β-proportional) — factor ≳10 short.
- **(b) Exact Savvidy sign — fails.** `H=−D²+2ad F`, `λ_min=−g|B|`, `N_neg=(gB)²V/8π²`; count and depth both grow with `|B|`; only `Tr(H_-)≤8C_N∫|F|³` (Lieb–Thirring), unbounded. Flux quantizes `N_neg` (never caps); topology anti-correlates (SD and ASD both BPS ⇒ PSD; the `Q=0` equal mixture is *maximally* unstable). Only `e^{−S}` suppresses — exactly the dynamical `W>Var−E[H]`.
- **Rigorous interval arithmetic** keeping only `W_geo≥0` + statistical CIs for bracket,X gives **−0.13 (L4)/−0.27 (L6)**: 0 *not* excluded; the sign lives in `W_geo`'s magnitude, which the theorem bounds below by 0 only.

## Five-point verdict (honest)
1. **Nonperturbative — PASS.** Exact constrained identity + Prékopa; no expansion; every β.
2. **Sign-controlling — PARTIAL; FAIL at 2.5.** `W_geo≥0` on every mode incl. `k→0` is a genuine value-free lower bound on the geometric half; the total sign needs `W_full>−Δ_dyn`, uncontrolled.
3. **Non-circular — PASS for `W_geo`; FAIL at closing.** Certifying `G″≥δ` at 2.5 ⟺ bounding `Var(G|0)` ⟺ H1 — a rename.
4. **Non-SUSY + fails free field — PASS (cleanest).** Massless free: zero mode flat ⇒ `G″=0`; YM β=0 SU2 `G″=0.248`; free 0. The Haar term discriminates.
5. **Uniform lower bound `G″≥λ` — FAIL.** Only at β=0 / weak UV (`V_eff″≥N/12`, β<1/48 via SZZ). At 2.5: only `≥0`, uncontrolled.

## Sharpest line
`W_geo≥0` is a theorem at *every* k including `k→0`, and it is the measured carrier of +0.51 — so the *direction* is right, and unlike the "DEAD-lean" the exact sign on the softest mode is positive. But the only certified statement is `W_geo≥0`; the measured +0.32 that beats the Savvidy deficit has no rigorous lower bound, and closing `Var(G|0)−E[H]<W_full` is H1. **The `k→0` sign is *motivated by* and *consistent with* +0.51, but the `k→0` lower bound is not provable from the present estimates — the exact failure is (◆)=H1; +0.51 is a sign probe, not a theorem.**

*Refs.* PROOF-JACOBIAN-SAVES (identity, `W_geo≥0`, (◆)); PROOF-SAVVIDY-EXACT-SIGN (`−D²+2adF`, `λ_min=−g|B|`, LT); PROOF-SAVVIDY-WGEO-DOMINATION (`ρ≳10`, category mismatch); k0-CONSTRAINT-EFFPOT (`V_eff″=G″`, β=0 `N/12`); INVENT-VARIANCE-INEQ (`G″=S−D`); COMPUTE-GPP-SU3 (+0.425/+0.329); W7-Gpp0 (+0.51); W7-RERUN (β=0 closure); PROOF-VALIDATED-GPP-SIGN (interval −0.13/−0.27); PROOF-COARSE-POINCARE (H1, β<1/48); VAR-FINAL (−0.165); Prékopa 1971; Brascamp–Lieb 1976; Savvidy, PLB **71** (1977) 133; Shen–Zhu–Zhu, CMP **400** (2023), arXiv:2204.12737.
