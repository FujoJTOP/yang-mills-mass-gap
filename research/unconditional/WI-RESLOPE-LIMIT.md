# WI-RESLOPE-LIMIT — the Re-slope functional is `A3` by necessity, NOT a forcing tool: any Re-slope-based inequality inherits weak-* non-computability

> **VERDICT (no QED, Clay 0%).** WI-MODIFIED-FUNCTIONAL found a robust-looking diagnostic: the **Re-slope** of
> the complex-`τ` zeros of `Z(τ)=w₀+Σfᵢ²e^{−Eᵢτ}` along `|Im τ|`-rungs (`dRe τ/d ln c₁ = 1/E₁`). This file
> pins its **exact limit**. It is **`A3` by necessity**, and **no Re-slope variant is a forcing lower bound**:
> the Re-slope reads the **simple pole at `s=E₁`** (equivalently the residue `f₁² = −Res Π` = `A3`); evaluating
> it needs the **complex-`τ` zero ladder** = analytic continuation (`A3`) or the full `Z` (`J→∞` = `H1`); a
> **finite real-`τ` sample is matched by an atomless bump to `O(h)`** (so any functional of it is blind); and
> the rung-invariance is moreover **fragile** (a weak continuum already breaks it). **A Re-slope-based
> inequality cannot force `f₁²>0` from finite data — it is `A3`/infinite, and at best a *presence/isolation*
> diagnostic.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · What the Re-slope is (recall)

For `Z(τ)=w₀+Σfᵢ²e^{−Eᵢτ}`: a simple atom at `E₁` gives zeros `τ_k = (ln(f₁²/w₀) − i(2k+1)π)/E₁`, i.e. on the
vertical line `Re τ = ln(f₁²/w₀)/E₁`, spaced `2π/E₁`. Under operator-reweighting `f₁² → λ f₁²` the real part
slides at `dRe τ/d ln λ = 1/E₁`, independent of the rung `k`. A continuum has no such clean ladder.

---

## §2 · The four exact facts (witness `w7/wi_reslope_limit.py`, clean)

1. **`[T]` Pure atom:** `Re(τ_k)=ln(f²/w₀)/E₁` exactly, `dRe/dln λ = 1/E₁` for **every** rung (verified Newton: `Re=−0.693147` at `k=0,1,2`, `Im` spacing `2π/E₁`).
2. **`[T]` Finite real-`τ` sample is blind:** over 12 real-`τ∈[0,5]` points, `max|Z_atom−Z_bump| = 8.9×10⁻³, 9.1×10⁻⁴, 9.2×10⁻⁵, 9.2×10⁻⁶` for bump width `h=10⁻¹…10⁻⁴` — **`O(h) → 0`**. Any functional of the finite real-`τ` data is blind (the atomless bump of width `h` matches it).
3. **`[V-num]` The rung-invariance is FRAGILE:** `Re(τ_k)` spread across `k=0,1,2` is `0.0` (pure atom), **`0.21` (atom + a weak continuum — broken)**, `1.9×10⁹` (pure cut — none). So the clean ladder requires an **essentially pure single level**; even a weak continuum shifts it. (This is the `τ`-coordinate form of `WI-ADMISSIBLE`'s "single-exponential dominance, not atom-vs-continuum".)
4. **`[T]` Recovering the ladder needs complex `τ`:** the zeros live at non-real `τ`; obtaining them from real-`τ` data is the ill-posed analytic continuation (= `A3`), or requires the full `Z(τ)` (= `J→∞`).

---

## §3 · The decisive classification: `A3` by necessity, no forcing

- **`[T]` The Re-slope reads the simple pole.** "Zeros on a clean vertical line with `dRe/dln λ = 1/E₁`" is exactly "the lightest singularity of the Laplace transform is a simple pole at `s=E₁`", whose residue is `f₁² = −Res Π` = **`A3` (read-not-force)**.
- **`[T]` Weak-* (ONLY-EXIT).** The atom `f₁²` is the **weak-* discontinuity** of the boundary density; the Re-slope reads it (weak-* discontinuous). Hence it is **not finitely computable**: a finite sample is matched by the atomless bump (fact 2) whose zero ladder differs (fact 3), so the discriminating content lives **entirely in the analytic continuation / full `Z`**.
- **`[NG]` No Re-slope-based forcing inequality.** Suppose `f₁² ≥ g(sample)` for a Re-slope-based `g`. To distinguish atom from bump, `g` must evaluate the zero ladder ⇒ analytic continuation (`A3`) or full `Z` (`J→∞`). Evaluated on finite real-`τ` data, `g` gives the **same** value for atom and bump (`O(h)`), so it **cannot** conclude `f₁²>0` while excluding `f₁²=0`. A **forcing** lower bound therefore does not exist.
- **`[C]` Best case is a presence/isolation diagnostic.** With complex-`τ`/full-`Z` (i.e. `A3`/infinite) the Re-slope can decide *presence* (clean ladder iff a dominant simple pole), but (i) that needs `A3`/`J→∞`, (ii) it is fragile to any continuum (fact 3), and (iii) at finite `L` all levels are poles so it is trivial — the real content is the `L→∞` survival of the isolated pole = isolation = `H1`.

---

## §4 · Honest bottom line

- **WI LINE 2 (Re-slope): classified `A3` by necessity, not a forcing tool.** Its clean signature is a *pure-single-level* property; the discriminating content is the complex-`τ` zero ladder = analytic continuation (`A3`) or the full `Z` (`J→∞`); a finite real-`τ` sample is blind (`O(h)`).
- **No Re-slope variant forces `f₁²` from finite data** — consistent with ONLY-EXIT (any forcing functional is weak-* discontinuous ⇒ `A3`/infinite = `H1`) and with WI-ADMISSIBLE (the raw complex-zero detects single-exponential dominance, not atom-vs-continuum).
- **What it is good for:** a *presence/isolation* diagnostic at the `A3`/infinite cost, fragile to continuum — i.e. a read, not a certificate.
- **Clay = 0%.** The Re-slope is one more `A3`-coordinate of `H1`, not an escape.

*Refs:* `WI-MODIFIED-FUNCTIONAL.md`, `WI-ADMISSIBLE.md`, `WI-EMBEDDED.md`, `ONLY-EXIT-SYNTHESIS.md`;
witness `w7/wi_reslope_limit.py`.
