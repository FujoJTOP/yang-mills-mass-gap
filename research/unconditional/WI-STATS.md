# WI-STATS — statistical feasibility of the weight-invariance (WI) atom-vs-cut signature

> **VERDICT (no QED).** The WI signature (operator-reweighting leaves the analytic-functional
> mass `𝒜` invariant iff there is an isolated atom) is **statistically feasible in principle but
> decisive only if the intrinsic operator-spread `Δ_cut` exceeds the noise floor** — which is itself
> unknown. Toy calibration: zero-spread (atom) vs `30%·m₁` spread (cut) separate at `3σ` once the
> per-estimate mass error `δm ≲ Δ_cut/2`. Translated to lattice: **SU(2)/SU(3) at β=2.5, L=4–6, needs
> `~10⁴–10⁶` configs (feasible on a cluster); L=8 SU(3) is at the exponential SNR wall (likely
> infeasible).** The WI test is a **[C]/[V-num] numerical diagnostic**, not a rigorous certificate: a
> "consistent with atom" result never rules out a low-density cut (the same thin-atom blindness as
> `PNT-OPE-RHO`). **Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric.
Witness: `w7/wi_stats_witness.py` (`py -3.12 -I`, runs clean).

---

## §1 · The statistic and the estimator

The WI test uses `K` interpolating operators `O_1…O_K` in the `0⁺⁺` channel (same physical spectrum,
different `tr F²` overlaps / smearing). For each, measure the correlator `C_k(τ)`, `τ=0…T−1`, at real
Euclidean time, and extract `𝒜_k` (the analytic-functional mass via a real-`τ` fit to
`C(τ)=w₀+F(τ)`). **Discriminant: `S = std_k(𝒜_k)`.**

- **Atom** (`C_k = w₀ + f_k²e^{−m₁τ}`): `𝒜_k = m₁` for *all* operators (the mass rides on `Im τ`
  = `π/m₁`; the overlap `f_k²` only moves `Re τ`). ⇒ intrinsic spread `0`.
- **Cut** (`C_k = w₀ + ∫ρ_k(E)e^{−Eτ}dE`): the nearest-zero functional depends on the spectral
  weight distribution, which differs across operators ⇒ `𝒜_k` moves.

**The test:** reject "atom" when `S` is significantly `> 0`. Power `∝ N`.

---

## §2 · Toy calibration (`w7/wi_stats_witness.py`)

`T=6, m₁=0.8, w₀=1`; atom `f_k² ∈ {0.3,1,5}`, cut width `Δ_k ∈ {0.5,1,2}` above `m_c=1.0`;
per-point noise `σ₀`; 80 trials, `K=3`.

| quantity | value |
|---|---|
| atom **intrinsic** spread | **0.0000** |
| cut **intrinsic** spread | **0.2439** (= `30%·m₁`) |

| `σ₀` | atom `S` | cut `S` | separation `z` | `N` (`s_base=1`) |
|---|---|---|---|---|
| 0.100 | 0.297 | 0.276 | −0.1 | 10² |
| 0.050 | 0.195 | 0.238 | 0.4 | 4·10² |
| 0.020 | 0.070 | 0.240 | **3.7** | 2.5·10³ |
| 0.010 | 0.039 | 0.245 | 8.7 | 10⁴ |
| 0.005 | 0.021 | 0.244 | 19.4 | 4·10⁴ |
| 0.002 | 0.008 | 0.244 | 49.0 | 2.5·10⁵ |

The discriminator is **strong in principle** (atom spread exactly zero, cut spread `30%·m₁`). It
becomes decisive (`z≥3`) once the per-estimate mass error `δm ≲ Δ_cut/2` (here `δm≈0.1` at `z=3.7`,
`Δ_cut=0.24`). **Scaling law:** `N_eff ≳ (2 s_base/Δ_cut)²`, i.e. `z ≈ (Δ_cut/δm)/2`.

---

## §3 · Translation to SU(2)/SU(3) at β=2.5

The decisive quantity is the **per-estimate 0⁺⁺ mass error** `δm`, which must satisfy `δm ≲ Δ_cut/2`.
Real glueball technology reaches `δm/m₁ ~ 1–2%` at `N ~ 10⁴–10⁵` configs (`L=4–6`), degrading
exponentially with `L` (signal `~f²e^{−m(T−1)}`).

| target | `δm` needed for `Δ_cut`=10–30% | `N` (configs) | verdict |
|---|---|---|---|
| SU(2), L=4–6 | `0.03–0.1 m₁` (3–10%) | `10⁴–10⁵` | **feasible** (cluster, days) |
| SU(3), L=6 | `0.03–0.1 m₁` | `10⁵–10⁶` | **feasible but heavy** |
| SU(3), L=8 | ≤ `0.03 m₁` | `10⁶–10⁷` + SNR wall | **likely infeasible** |

`K` operators reuse the same configurations (no extra `N`, only more measurements).

---

## §4 · The two unknowns that decide feasibility

1. **`Δ_cut` (the intrinsic operator-spread) is UNKNOWN for the real 0⁺⁺.** The toy's `30%` comes
   from a toy continuum placed above `m₁`; the real lightest `0⁺⁺` sits *below* the two-particle
   threshold, and its continuum coupling may be **weak** (`Δ_cut` small). If `Δ_cut ≪ 10%·m₁`, the
   test needs `N ≳ 10⁶⁺` and hits the **exponential SNR wall** — the same thin-atom blindness that
   `PNT-OPE-RHO` proved for all finite-measurable data. **So feasibility is conditional on `Δ_cut`,
   which is `H1`-adjacent.**
2. **The `𝒜` extraction is an ill-posed continuation.** The toy used a `2`-parameter real-`τ` fit
   (it works only because the toy is exactly 3-parameter); the real `𝒜` needs the full continuum
   from `T=4–8` points — strictly harder. A real study would use a `T`-large effective-mass plateau
   + `2`-state fit; the systematic (model) error can exceed the statistical one at `L=8`.

---

## §5 · Honest verdict

- **Feasible as a numerical experiment** for SU(2)/SU(3) at `L≤6` with `~10⁵–10⁶` configs — a real
  cluster campaign, **not** a laptop job, and **not** a one-shot (needs `K≈3–5` operators, a
  controlled continuum fit, and a plateau study).
- **Decisive only if `Δ_cut` is not small.** If the 0⁺⁺ is a thin atom above a low-density continuum,
  the intrinsic spread falls below the noise floor and the test is **inconclusive at any feasible `N`** —
  recovering, in the WI coordinate, the `PNT-OPE-RHO` blindness.
- **Never a rigorous certificate:** a "zero spread" outcome is *consistent with* an atom but cannot
  *prove* `f₁²>0` or exclude a low-density cut; the ill-posed continuation is a [C]/[O] boundary.

**Net:** WI-STATS says the weight-invariance signature is a **genuine, measurable handle** at
`~10⁵–10⁶` configs for `L≤6`, worth a real lattice run — but its **decisiveness is gated by the
unknown intrinsic spread `Δ_cut`**, so it is a [C]/[V-num] diagnostic, not a verification of `H1`.
**Clay 0%; no QED.**

---

*Refs:* `CONSTRUCT-IR-ANALYTIC-NUM.md` (the signature), `PNT-IR-ANALYTIC.md` (the functional),
`PNT-OPE-RHO.md` (thin-atom blindness), `SLOT-BOOTSTRAP.md` (`Δ_min(K)=0` floor-blindness);
witness `w7/wi_stats_witness.py`.
