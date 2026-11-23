# REPORT-NA3-T — Dobrushin / SDPI coefficient of the SU(3) link heat bath

> 2026-11-23 · `research/lattice-cpb/dobi.py`, `dobi_run.py`, `dobi_analyze.py`.
> **Observables only** — the heat-bath dynamics (`su3.py`, `su3_numba.py`) is
> untouched.  No cluster expansion, no defect weight is used anywhere.
>
> **Every grid point is now labelled `cert` or `est`.**  See §0.3 — a kill is
> strict only when a recorded witness pair gives a certified lower bound ≥ 1;
> a `λ̂ < 1` is evidence, never a survival certificate.

## 0. Definitions and estimators implemented

`P_e(·|U_{≠e}) ∝ exp((β/3) Re tr(U_e S_e)) dU`, `S_e` = staple sum;
`I(e,e') = sup_{U,U' differ at e'} ‖P_e(·|U) − P_e(·|U')‖_TV`;
`η = max_e Σ_{e'} I(e,e')`.

### 0.1 Structure (exact)

The staple at `e` depends on exactly **18** other links; `D.neighbor_links`
enumerates them and `D.staple_at` reproduces `su3.staple_sum` to 10⁻¹².

### 0.2 Central TV estimate — **`est` only**

Both conditionals are matrix-Fisher distributions.  The normalisation ratio is
estimated with the **Bennett acceptance-ratio (BAR)** method from samples of
each conditional (CM subgroup heat bath), then `TV = 1 − ∫min(p,q)`.  This is a
*point estimate*: it is used to report η̂ and to guide the search for the worst
link, **not** to certify anything.

### 0.3 Certified lower bound — **`cert` capable**

Witness function `f(U) = Re tr(U ΔS) / ‖ΔS‖_* ∈ [−1, 1]`, with `ΔS = S_q − S_p`.
For **any** two distributions, `TV(p,q) ≥ ½ |E_p f − E_q f|`.  This is exact
and needs no coupling argument.  We estimate each expectation by 300 CM
samples and report

```
TV_lb(e,e') = max(0, ½|E_p f − E_q f| − 4·se) ,
```

a 4σ lower confidence bound.  Then for the recorded link `e` and neighbour
change,

```
η  ≥  Σ_{e'} I(e,e')  ≥  Σ_{e'} TV_lb(e,e')  =:  row_lb .
```

If `row_lb ≥ 1` the configuration is a **witness**: `η ≥ 1` is *certified*, and
the pair `(U, {g_e'})` is saved to `runs/witness_L4_b{β}.npz` (configuration,
link coordinate `x,μ`, the 18 neighbour coordinates and replacement matrices,
and the per-pair `(lb, value, se)`).  The label is `cert`; otherwise `est`.

**Asymmetry, stated plainly.**  `cert` ⟹ Dobrushin fails (no uniform
contraction).  `est` with `η̂ < 1` would mean only *no counterexample found in
this sample* — it would **not** certify `η < 1` and would not be a survival
result for NA-3.  Only the kill direction is reachable by this construction.

## 1. Results (L=4, 30 links per β)

| L | β | η̂ (row sum) | η max | **witness lb** (max link) | cert/est | mean TV/pair | mean ‖ΔS‖_* |
|---|---|---|---|---|---|---|---|
| 4 | 2.2 | 5.74(0.23) | 8.50 | **1.18** | cert | 0.319 | 2.41 |
| 4 | 3.0 | 8.38(0.36) | 12.33 | **1.75** | cert | 0.466 | 2.43 |
| 4 | 6.0 | 17.19(0.10) | 17.98 | **1.82** | cert | 0.955 | 2.43 |

`τ_P(L=4)` = 0.52 / 0.50 / 1.61.  `m_eff` at L=4 is **not resolvable** — the
connected plaquette correlator dies within one lattice spacing at these
couplings (main report §4.3), so no `m_eff` is quoted.

Witness files: `runs/witness_L4_b2.2.npz`, `…b3.0.npz`, `…b6.0.npz`.
`results_eta.csv` holds the same rows with the `cert` column and the witness
filename in `notes`.

## 2. Reading in the agreed order (R1–R5)

* **R1 — β=6.0 is the `η→1` control and never fires D3.**  `η̂ = 17.19`, at the
  trivial cap (18 neighbours × TV≈0.95).  This is the known failure of the
  crude single-link sufficient condition, **not** an architecture verdict.
  Certified (`row_lb = 1.82`) but reported, not counted.
* **R2 — β=2.2 is the BL-6-redundant point.**  `η̂ = 5.74`, certified
  (`row_lb = 1.18`).  Reported, not counted.
* **R3 — β=3.0 is the only decisive point — and it is now a certified kill.**
  `η̂ = 8.38 ± 0.36`, and the witness pair gives
  **`η ≥ 1.75 ≥ 1`** (`runs/witness_L4_b3.0.npz`).  The strict statement is
  therefore not "no contraction was seen" but **"η ≥ 1 at β=3.0, hence the
  single-link Dobrushin bound cannot contract"**.  The certified bound
  (~1.75) is well below the point estimate (~8.4) because the witness
  function is weak; **the gap only widens the margin, it does not affect the
  direction.**
* **R4 — D4 is not assessable.**  `1−η < 0` at every β (certified), but the
  comparison object `m_eff` cannot be measured at L=4.  Per the agreed rule
  this is "the comparison carries no information", not "NA-3 is dead **for
  the wrong reason**".
* **R5 — full table reported above**; no binary D4 verdict is issued here.

## 3. What the numbers say (decision input only)

`η̂` **increases monotonically with β** (5.74 → 8.38 → 17.19) and approaches
the trivial cap 18 at β=6.  **Every measured coupling is certified to have
`η ≥ 1`.**  There is no coupling anywhere in this window where the single-link
Dobrushin coefficient could be below 1.

Consistent with the card's §6 ("most likely: D3"), the non-cluster local route
does not produce contraction at `q≈1`, and this is now a *measured fact with a
witness*, not a point-estimate extrapolation.  What is **not** established: the
`b=2` block kernel (different kernel, would need its own sampler), and the
`m_eff` comparison needed for D4.

## 4. Boundary (verbatim from the card §6)

> **通过只是决策输入，不建立 uniform 收缩，不进 Lean 链，不替代 `CPB-W/S`。**
> 判死则触发 `NEW-ARCHITECTURE.md` §8 的结构性结论（触发须是**实测事实**：D1/D2/D3 之一）。

A `cert` row is a 实测事实 and may serve as the D3 trigger input.  It does
**not** enter the Lean chain and does not replace `CPB-W/S`.

## 5. Limitations (all explicit)

1. **Statistical, not interval-arithmetic.**  `TV_lb` subtracts 4σ of
   Monte-Carlo error under a normal approximation.  It is a rigorous lower
   bound in the statistical sense, not a machine-checked interval.  Raising
   `--nsamp` tightens it; it never weakens the direction.
2. **Typical, not supremum.**  Configurations are sampled from the equilibrium
   measure; η is a supremum over configurations.  This only makes the reported
   `η̂` an *under*-estimate — harmless for a kill.
3. **One witness link is enough.**  The certificate uses a single (link,
   neighbour-change) pair; that is all a lower bound on a supremum needs.
4. **The `b=2` block kernel is not implemented** — only the link kernel here.
5. **`m_eff` at L=4 is unavailable**, so D4 cannot be tested.
6. SU(2) was not attempted (non-blocking per the card).

## 6. Files / reproduction

```
python dobi_run.py --L 4 --beta 2.2 --ncfg 15 --nlink 2 --nsamp 300 --nsweep 12 \
    --therm 600 --sep 5 --out runs/eta_L4_b2.2.npz
python dobi_run.py --L 4 --beta 3.0 --ncfg 15 --nlink 2 --nsamp 300 --nsweep 12 \
    --therm 600 --sep 5 --out runs/eta_L4_b3.0.npz
python dobi_run.py --L 4 --beta 6.0 --ncfg 15 --nlink 2 --nsamp 300 --nsweep 12 \
    --therm 600 --sep 5 --out runs/eta_L4_b6.0.npz
python dobi_analyze.py        # -> results_eta.csv, na3t_table.txt
```

Each run also writes `runs/witness_L4_b{β}.npz` (the best-certified pair).
