# COMPUTE-GPP-K0-CONTINUUM — `G″(0)` at the softest nonzero `k=2π/L` on ONE `L=8` volume: positive, `W_geo`-carried, volume-stable — a κ→0 sign probe, not H1

> 2026-11-23 · The `k→0`-proxy measurement the master doc flagged as pending (`GPP-kzero`, L=8 softest mode). **Single process, memory-safe streaming** (configs never stored; only `N×Nmode` scalars accumulated; peak RSS **91 MB** for SU(3), **70 MB** for SU(2), well under 200 MB — the machine's prior OOM was from parallel `SU(N)` L=8 jobs, not this). Two single jobs, run **sequentially**: SU(3) β=5.7, L=8, N=400, gap=4, therm=350, seed 808, 714 s; SU(2) β=2.4, L=8, N=400, gap=5, therm=600, seed 910, 364 s. Code: `research/unconditional/w7/k0_continuum_run.py` (stream) + `k0_continuum_analyze.py` (delete-block jackknife). Data: `k0_su3_b5.7_L8_n2.npz`, `k0_su2_b2.4_L8_n2.npz`.

**Verdict.** On a single `L=8` volume at the SMALLEST nonzero momentum `k=2π/L=π/4`, the SU(3) β=5.7 coarse-mode curvature is **positive in 6/6 transverse plane-wave modes**: `G″ = +0.307 ± 0.055` (delete-block jackknife, 20 blocks). It is **`W`-carried**: the action bracket is statistically zero (`E[H]−Var(G) = +0.037 ± 0.099`, the equipartition cancellation), and `W_geo = E[Φ″]−Var(Φ′) = +0.319 ± 0.005` alone accounts for all of `G″`. The decisive new fact is **`W_geo` is volume- and `k`-independent to 0.5 %**: `+0.318` (L=4) / `+0.319` (L=6) / **`+0.319` (L=8, k=π/4)** / `+0.320` (L=8, k=π/2). The old L=4→6 confound (k and volume changing together) is **removed two ways**: (i) at **fixed `k=π/2`**, L=8 gives `+0.396 ± 0.063` vs L=4 `+0.425 ± 0.084` — volume-stable; (ii) at **fixed `L=8`**, the k-scan `π/4 → π/2` gives `+0.307 → +0.396`, a mild *rise* with k (no collapse toward `k→0`). **But keeping only the proven `W_geo ≥ 0`** and treating `bracket + X` statistically gives a 95 % one-sided lower bound of **−0.098** — 0 is not excluded. The **SU(2) β=2.4 cross-check reproduces the picture** on a second gauge group: `G″(π/4)=+0.192 ± 0.058` (12/12 modes positive overall), `W_geo=+0.252 ± 0.003`, `1/Var=0.270` (matching the SU(2) value `0.272` at L=4,6), and a **tighter** rigorous-style bound `−0.035`. This is a **κ→0 sign probe with the volume confound removed, not H1**: uniform `G″ ≥ λ` in the `k→0, L→∞` limit (⇔ the `Var(G|0)` bound) remains open.

## 1. What was measured (exact identity, same conventions as COMPUTE-GPP-SU3)

Flat exponential chart `X = log U ∈ su(3)` (branch-safe traceless principal log, center twist corrected so `exp(iX)=U` exactly), mode amplitude `a = ⟨v,X⟩`, action `V = S_W + Φ` with `Φ = −log J` the Haar Jacobian. Exactly (constrained free energy):

```
G″(0) = bracket + W ,    bracket = E[H|0] − Var(G|0) ,
W = W_geo + X ,   W_geo = E[Φ″|0] − Var(Φ′|0) ,   X = −2 Cov(G,Φ′|0) .
```

`H`, `G`, `Φ′`, `Φ″` are flat directional derivatives (central differences `h=10⁻³` for the action; analytic eigenvalue-perturbation theory for the Jacobian, no branch-cut FD). Modes: 6 transverse plane waves, one polarisation × one colour each, `n=1` (`k=π/4`) and `n=2` (`k=π/2`), unit-normalised. The gold standard `G″` is the curvature of `−log P(a)` (Gaussian KDE, `bw_method=0.2`); the decomposition is a Gaussian-kernel conditional average at `a=0` (`bw=0.30`). Errors are **delete-block jackknife** over 20 contiguous blocks of 20 stored configs (= 80 sweeps each, ≫ τ; measured `τ = 1.5–2.6` stored configs).

**Memory discipline.** One process, one `Lattice` object; each config is measured then discarded (`k0_continuum_run.py` reuses `S0` across modes, bit-identical to `rerun_run.measure`). Only `a,H,G,P1,P2` (`N×12` floats) are held. Peak RSS 91 MB (SU(3) `L=8`).

## 2. Main result — SU(3) β=5.7, L=8, N=400

| `k` | n | modes | `G″` | `bracket` | `W_geo` | `X` | `W` | `bracket+W` | `1/Var(a)` |
|---|---|---|---|---|---|---|---|---|---|
| **0.785 (π/4)** | 1 | 6 | **+0.307 ± 0.055** | +0.037 ± 0.099 | **+0.319 ± 0.005** | +0.035 ± 0.019 | +0.353 ± 0.019 | +0.391 ± 0.100 | 0.394 |
| 1.571 (π/2) | 2 | 6 | +0.396 ± 0.063 | −0.000 ± 0.093 | +0.320 ± 0.005 | −0.015 ± 0.036 | +0.305 ± 0.036 | +0.304 ± 0.099 | 0.400 |
| **all** | — | 12 | **+0.352 ± 0.045** | +0.019 ± 0.073 | **+0.319 ± 0.004** | +0.010 ± 0.017 | +0.329 ± 0.019 | +0.348 ± 0.076 | 0.397 |

At `k=π/4` the six per-mode `G″` are `+0.340, +0.333, +0.366, +0.250, +0.266, +0.291` — **all positive**, range `0.250–0.366`. The pooled mean is `+7.8σ` from 0 **statistically** (all 12 modes positive; pooled `+0.352 ± 0.045`). Closure `bracket+W − G″ = +0.084` (k=π/4), `−0.092` (k=π/2), `−0.004` (all) — the known SU(3) chart-cut/log-amplitude systematic of ±0.09, sign-unaffecting. KDE-bandwidth robustness (`kdebw` 0.12/0.20/0.35): `G″(π/4) = +0.329/+0.307/+0.305`; conditioning-bandwidth (`bw` 0.20/0.30/0.40): `W_geo = +0.335/+0.319/+0.301`. The sign never changes.

## 3. The volume confound is removed — two independent ways

Prior SU(3) runs reached `k=π/2` at L=4 and `k=π/3` at L=6, so the L=4→6 drop (`+0.425 → +0.329`) mixed a **k-decrease with a volume change**. With one `L=8` volume we separate them:

- **Fixed `k=π/2`, vary volume.** L=4: `G″=+0.425 ± 0.084`, `W_geo=+0.318`. L=8: `G″=+0.396 ± 0.063`, `W_geo=+0.320`. Consistent within errors → the earlier drop was **not** a volume effect; `G″` at fixed k is volume-stable.
- **Fixed `L=8`, vary k.** `G″(π/4)=+0.307 ± 0.055` < `G″(π/2)=+0.396 ± 0.063` (Δ = `0.089 ± 0.083`, ~1σ): moving to **smaller k does not collapse** `G″`; it stays positive and roughly flat, with a mild rise toward the UV.
- **`W_geo` invariance.** `+0.318 / +0.319 / +0.319 / +0.320` at (L4,k=π/2), (L6,k=π/3), (L8,k=π/4), (L8,k=π/2). The geometric carrier is constant in `(L,k)` to 0.5 %.

So the `k→0` proxy at L=8 is **`G″(k=π/4)=+0.31`, W-dominated**, and the previously ambiguous L-volume trend is consistent with a volume-stable positive value that **does not collapse** as `k` decreases: `G″` drifts down mildly from `+0.396` (k=π/2) to `+0.307` (k=π/4) but stays decisively positive, with `W_geo` pinned at `+0.32`.

## 4. Why it is still a sign probe, not H1 (honest interval)

The exact identity separates one term with a **proven** sign from three with only measured values. `W_geo ≥ 0` for every mode, every k, uniform in L (Φ convex ⇒ Prékopa; `Φ″ ≥ N/12` pointwise) — a genuine theorem. But the theorem bounds `W_geo` from below by **0**, while the *measured* `W_geo=+0.319` is what carries the sign. Keeping only the proven part and treating `bracket + X` statistically:

```
(bracket + X) = +0.028 ± 0.077   →   95 % one-sided lower bound = −0.098 .
```

So even at L=8 the rigorous-style interval **does not exclude 0**. Per k: `k=π/4` lower bound `≈ −0.094`, `k=π/2 ≈ −0.180`, pooled `−0.098`. This is the *same structural obstruction* as PROOF-VALIDATED-GPP-SIGN (−0.13 at L4, −0.27 at L6), now tighter (−0.10) but still not a certificate. The missing inequality is exactly

```
(◆)  Var(G|0) − E[H|0]  <  W_full ,   uniform in (L, U′)   ⟺   Var(a) ≤ 1/δ   ⟺   m_eff²(k→0) ≥ δ  =  H1.
```

The measured facts all point the right way — `bracket ≈ 0`, `X` small, `W_geo ≈ +0.32` stable — but the *magnitude* of `W_geo` beyond 0 is not certified, and closing `(◆)` is H1 / the coarse-Poincaré gap. One lattice size is not the continuum: `a(β=5.7)≈0.18 fm`, `L=8` ⇒ `≈1.4 fm` box, `k_min=π/4` in lattice units; this is a **κ→0 sign probe on a single volume**, not the `a→0` limit.

## 5. SU(2) β=2.4, L=8 — cross-check on a second gauge group

Same pipeline, `N=400`, `gap=5`, `therm=600`, seed 910, 364 s, peak RSS 70 MB (SU(2) equilibrates fast: plaquette flat at `⟨u_□⟩=0.634` from sweep 100). Both `k` show **6/6 positive modes**:

| `k` | modes | `G″` | `bracket` | `W_geo` | `X` | `W` | `1/Var(a)` |
|---|---|---|---|---|---|---|---|
| **0.785 (π/4)** | 6 | **+0.192 ± 0.058** | +0.050 ± 0.067 | **+0.254 ± 0.004** | +0.025 ± 0.022 | +0.278 ± 0.021 | 0.270 |
| 1.571 (π/2) | 6 | +0.262 ± 0.052 | +0.028 ± 0.061 | +0.250 ± 0.005 | −0.043 ± 0.025 | +0.208 ± 0.026 | 0.283 |
| **all** | 12 | **+0.227 ± 0.031** | +0.039 ± 0.036 | **+0.252 ± 0.003** | −0.009 ± 0.015 | +0.243 ± 0.015 | 0.276 |

Exactly the SU(3) pattern, with the same structure and the same limitation: positive softest-mode curvature, carried by `W_geo`, action bracket ≈0; the k-trend is again a mild *rise* with k (`π/4: +0.19 < π/2: +0.26`), no collapse toward `k→0`. Two quantitative agreements worth recording: (i) `1/Var(a)=0.270` at `β=2.4` L=8 matches the SU(2) `1/Var=0.272` quoted at L=4,6 in OTHER-GPP-KZERO — volume-stable Gaussian proxy; (ii) the rigorous-style interval is **tighter than SU(3)**: `(bracket+X)=+0.030 ± 0.040` → 95 % one-sided lower bound **`−0.035`**, i.e. much closer to 0 but still not excluding it. `τ = 2.2` stored configs (mean), so 20-delete-block jackknife (block = 20 configs = 100 sweeps) is safe.

## 6. Sharpest line

On the softest accessible nonzero mode of a single `L=8` volume, the coarse-mode curvature is **positive, and its sign is the Haar/Jacobian term `W_geo`** — the one object proven `≥0` on every mode, here measured at the remarkable fixed value `+0.32` across L=4,6,8 and both k — while the Savvidy/Wilson action bracket cancels to equipartition (`≈0`). The volume confound of the earlier L-series is gone: `G″` is volume-stable at fixed k and does not collapse as k→0. That is the strongest *numerical* support yet for `G″(0)>0` at `k→0`; it is **not** a proof, because the only certified half is `W_geo≥0` (value-free), and the rigorous-style lower bound from these samples is **−0.10 (SU(3)) / −0.035 (SU(2))**, not positive. Closing it is one variance inequality, `(◆)` = H1 — not a better dataset.

*Refs.* COMPUTE-GPP-SU3 (L=4,6 baseline: +0.425/+0.329, W_geo +0.318/+0.319); OTHER-GPP-KZERO (+0.51 SU3 β=2.5 exact marginal, sign-probe verdict); PROOF-VALIDATED-GPP-SIGN (rigorous-style interval −0.13/−0.27, (◆)); PROOF-JACOBIAN-SAVES (`W_geo≥0`, `Φ″≥N/12`, Prékopa); PROOF-SAVVIDY-WGEO-DOMINATION (β-vs-β-independent scaling wall); PROOF-COARSE-POINCARÉ (H1); Prékopa 1971; Brascamp–Lieb 1976; scripts `research/unconditional/w7/k0_continuum_run.py`, `k0_continuum_analyze.py`; data `k0_su3_b5.7_L8_n2.npz`, `k0_su2_b2.4_L8_n2.npz`; results `k0_su3_b5.7_L8_n2_results.json`, `k0_su2_b2.4_L8_n2_results.json`.
