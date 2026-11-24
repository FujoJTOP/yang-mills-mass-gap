# BLAST-RESUMMATION — aggressive resummation of the strong-coupling series does NOT certify `χ ≤ 1/λ` at β=2.5: plain truncation is useless, Borel–Padé recovers only the (non-physical) PV, and the IR-renormalon transseries ambiguity (σ free) leaves an O(20–45%) uncontrolled error

> **VERDICT (no QED, Clay 0%).** Brute-force resummation of the strong-coupling series toward β=2.5 was
> attempted with Borel–Padé + optimal-truncation analysis at the physical point. **It fails as a certificate.**
> Three reasons, all real: **(1)** the strong-coupling series is a strong asymptotic series whose Borel
> singularity sits at the physical (crossover/mass-gap) scale `t_c ~ O(1)`; at `u(2.5)=0.507` the optimal
> truncation is **`n*=0`** terms for `t_c ≲ 1`, so **plain truncation gives an O(1) error** (no useful estimate).
> **(2)** Borel–Padé *does* reproduce the principal value PV to `<0.1%` (even for `t_c<u`), **but the PV is not
> the physical value**: the IR-renormalon structure makes the physical answer `PV + σ·e^{−t_c/u}` with **σ a free
> transseries coefficient** — no saddle (ROOF-RESURGENCE) — so the physical `χ` has an intrinsic
> **O(e^{−t_c/u}) ≈ 20–45% ambiguity** that no prescription (Padé, conformal, median) removes. **(3)** no
> rigorous error bound exists (asymptotic series). **So resummation is a numerical diagnostic, not a certificate.**
> The exact missing ingredient = the physical Borel-singularity location + `σ` = the condensate = `H1`.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理/启发 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric.

---

## §1 · Setup (the physical point)

SU(2) strong-coupling parameter `u(β)=I₂(β)/I₁(β)`:

| β | 2.0 | 2.3 | **2.5** | 2.7 |
|---|---|---|---|---|
| u | 0.433 | 0.479 | **0.507** | 0.533 |

The **rigorous** radius (KP/SZZ) is `β* ≈ 1/48 = 0.0208` (`u ≈ 0.02`); β=2.5 is **~2 orders in β / ~25× in u**
outside it. So the strong-coupling series is evaluated deep in its asymptotic regime.

## §2 · Plain truncation is useless at β=2.5 `[V-num]`

For an asymptotic series with a Borel singularity at `t_c` (coefficients `~Γ(n+1+2b)/Γ(1+2b) t_c^{−n}`),
the optimal-truncation order is `n* ≈ t_c/u`. At `u=0.507`:

| t_c | 0.3 | 0.5 | 0.7 | 1.0 | 1.5 | 2.0 |
|---|---|---|---|---|---|---|
| n* | 0 | 0 | 0 | 0 | 1 | 2 |
| minimal-term error | 1.00 | 1.00 | 1.00 | 1.00 | 0.68 | 0.39 |

For the physical crossover/mass-gap range `t_c ≲ 1`, the term ratios exceed 1 from the first term ⇒ **optimal
truncation is `n*=0`, minimal term O(1)** — the truncated series carries **no** usable information. (The
`e^{−t_c/u}` column is only the large-`t_c/u` asymptotic floor, a *lower* bound on the achievable error; the
real error at the physical `t_c` is worse.)

## §3 · Borel–Padé recovers the PV, but the PV is not physical `[V-num]/[C]`

Toy (Borel transform `(1−t/t_c)^{−1/2}`, `α=1/2`, summable): Borel–Padé `[8/8]` reproduces the principal-value
Borel sum to **0.00%** (`t_c=0.7>u`) and **0.02%** (`t_c=0.4<u`) — i.e. Borel–Padé **beats plain truncation** and
is numerically accurate. **But:** the real parts of the upper/lower lateral sums coincide (both = the PV); the
IR renormalon makes the **physical** value `PV + σ·e^{−t_c/u}` with **σ free** (ROOF-RESURGENCE: `u=2` IR
renormalon has no saddle, the condensate is a free transseries parameter). So the physical `χ` differs from the
PV by `σ·e^{−t_c/u}`, an **intrinsic O(20–45%) ambiguity** (σ undetermined) that **no** resummation prescription
removes, and there is **no rigorous error bound** (asymptotic series).

## §4 · Honest bottom line

- **Not a certificate.** Plain truncation: `[NG]` (useless). Borel–Padé: numerically reproduces the PV but
  `[NG]` as a certificate — the PV is not the physical `χ`, and the transseries σ-ambiguity (O(20–45%)) has no
  rigorous bound.
- **Why:** the series is asymptotic with an IR-renormalon singularity at the physical scale; the nonperturbative
  σ is exactly the gluon-condensate freedom = the `H1` content. Resummation **re-expresses** the wall, does not
  cross it.
- **Verdict:** `χ ≤ 1/λ` at β=2.5 **cannot be established by aggressive resummation** — it is a diagnostic
  `[V-num]/[C]`, never a certificate. **Clay = 0%.**
- **Exact missing step `[O]`:** the physical Borel-singularity location + the transseries coefficient `σ` (the
  condensate) — i.e. `H1` in the resummation coordinate.

*Witness:* `research/unconditional/w7/blast_resummation.py` (runs clean, `py -3.12 -I`, numpy/scipy).
