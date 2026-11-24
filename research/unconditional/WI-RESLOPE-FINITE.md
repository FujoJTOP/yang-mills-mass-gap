# WI-RESLOPE-FINITE — WI LINE 2 continued: the Re-slope cannot be moved to finite real-τ data without collapsing to P3 (READ-not-FORCE); the rung/slope proper stays A3, and even its continuum-robustness requires the gap

> **VERDICT (no QED, Clay 0%).** The Re-slope functional
> (`dRe τ/d ln c₁ = 1/E₁` along `|Im τ|`-rungs, `WI-MODIFIED-FUNCTIONAL`) is robust to a continuum
> **precisely because it lives on the complex-τ zeros** of `C_O(τ)=w₀+Σcᵢe^{−Eᵢτ}` — i.e. it is A3
> (analytic continuation, read-not-force). Pushing it to **finite real-τ data** produces a functional that
> *exists* but **collapses to P3** (the effective-mass plateau + approach law): the plateau value is
> operator-independent only as `t→∞` (P1), the operator spread decays as `e^{−gap·t}` (P2), and atom-vs-cut is
> separated by the **approach law** (exponential vs power-law, P3). All of these **READ** the gap and need
> `t→∞` **uniformly = isolation = `H1`**. The "slope `1/E₁` / rung" structure does **not** survive as an
> independent finite-data feature — it is the complex-τ dual, and it is the tautology
> `Re τ₀/(ln c₁−ln w₀)=1/E₁`, which needs `E₁` (the plateau) + `ln c₁` (the amplitude), each requiring
> `t→∞`. **So `WI-RESLOPE-FINITE` = a real read/diagnostic, NOT a forcing tool; and even the complex rung's
> robustness needs the gap.** Agent-level `[T]/[V-num]/[NG]`; Clay 0%.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric.
Acid test: must fail on free massless (`f₁²=0`, `G″=0`) *and* on the gapped-atomless `:φ²:` (`f₁²=0`).

---

## §1 · What the Re-slope is (complex-τ), and that it is A3

For an isolated atom (`w₀+ c₁e^{−E₁τ}`) the zeros are `τ_n = [ln(c₁/w₀) − i(2n+1)π]/E₁`:
- **rung** `|Im τ_n| = (2n+1)π/E₁` is **independent of `c₁`** (operator-reweighting-invariant);
- **`Re τ_n = ln(c₁/w₀)/E₁`**, so `dRe τ_n/d ln c₁ = 1/E₁` — the Re-slope.

Witness (A): `c₁∈{0.4,1,2.5,5}` → `Re τ₀ = ln(c₁/w₀)` exactly, `Im τ₀=−π` exactly
(`π/|Im|=1.000`), slope `= 1.00000 = 1/E₁`. **This is a property of the complex zeros** — you cannot read
`Re τ_n`/`Im τ_n` from real-τ values, so obtaining it requires analytic continuation = **A3**.

Crucially, witness (E): under a **gapped** continuum (`atom+gap`, edge 2 > E₁) the rung persists
(`π/|Im| = 0.944→1.000`), so the Re-slope is robust *to a gapped continuum* — i.e. its robustness already
**requires the gap** (a threshold continuum smears the rung, per `WI-MODIFIED-FUNCTIONAL`). So even the
complex Re-slope needs **isolation**.

---

## §2 · What finite real-τ data carries (and why it is P3, not the Re-slope)

At **real `t`** the recoverable objects are the `WI-PLATEAU-RIGOR` quantities:

| real-τ object | content | status |
|---|---|---|
| effective mass `m(t)=−d/dt log(C−w₀)` | `m(t)→E₁` for every operator (**P1**); finite-`t` value operator-dependent, spread decays `O(e^{−gap·t})` (**P2**) | READ (needs `t→∞`) |
| **approach law** of `m(t)` | atom/atom+gap: **exponential** `E₁+Ce^{−gap·t}`; pure/threshold cut: **power-law** `E₁+c/t` (**P3**) | READ (needs `t→∞`) |
| log-amplitude `a(t)=ln(C−w₀)+E₁t` | `→ ln c₁` only as `t→∞` | READ (needs `t→∞`) |

Witness:
- (B) `atom+gap`, `c₁∈{0.4,1,2.5,5}`: `m(t)` spread runs `0.0429→0.0115→0.0010→0.0000` (decays with `t`).
- (C) `atom 1.0000 / atom+gap →1.0000 exp / atom+thr slower / cut 1.249→1.032` (**power-law**, P3).
- (D) `a(t)→ln c₁` (exactly for a pure atom; for `atom+gap` contaminated at small `t`, clean only as `t→∞`).

---

## §3 · Why no independent finite-data Re-slope exists

The complex Re-slope content is `1/E₁ = Re τ₀/(ln c₁ − ln w₀)` (witness E: the ratio is `1.000`, a
**tautology**). To evaluate it from real data you need **both** `E₁` (the plateau) **and** `ln c₁` (the
amplitude) — and each requires `t→∞` (the plateau) / the continuum-decayed amplitude. So the finite-data
"Re-slope" is nothing but a rearrangement of `(E₁, ln c₁)`, both obtained from the `t→∞` asymptotics.
It carries **no robust finite-`t` content beyond P3**, and P3 itself certifies the atom only in the
`t→∞` (isolation) limit.

Hence: **a finite-real-τ functional exists (the plateau+approach = P3), but it does not reproduce the
Re-slope as an independent, robust, finitely-computable object** — that object is the complex-τ
dual (A3).

---

## §4 · Verdict (READ-not-FORCE, doubly blocked)

- **The rung/slope Re-slope proper is A3** (complex-τ analytic continuation), and its continuum-robustness
  **requires the isolation gap** (witness E) — so it is not even robust to a gapless continuum.
- **The finite-real-τ version collapses to P3** (plateau value + exponential-vs-power-law approach), which
  **READS** the gap and needs `t→∞` **uniformly = isolation = `H1`**.
- **It does not force.** There is no finite-real-data functional that lower-bounds `f₁²>0` (consistent with
  `ONLY-EXIT`: any forcing functional must be weak-* discontinuous = not finitely computable = A3/`H1`).
- **Acid**: the P3 diagnostic is acid-clean (free massless → no plateau at `E₁>0`; `:φ²:` → power-law
  approach, correctly *not* flagged as an atom), but it is a **read**, not a certificate.

**Bottom line.** Continuing WI LINE 2 does the same as LINE 1: it sharpens the *diagnostic* (the Re-slope is
a clean geometric picture of the isolation, realisable at finite `τ` only as P3) but **does not escape
read-not-force**. The Re-slope is the **A3 coordinate of `H1`**, finitely computable only in the `t→∞`
(isolation) limit. **Clay 0%; no QED.**

*Files:* `research/unconditional/WI-RESLOPE-FINITE.md`, `research/unconditional/w7/wi_reslope_finite.py`
(runs clean, `py -3.12 -I`). *Refs:* `WI-MODIFIED-FUNCTIONAL.md`, `WI-PLATEAU-RIGOR.md`, `WI-GENERAL.md`,
`WI-EMBEDDED.md`, `ONLY-EXIT-SYNTHESIS.md`.
