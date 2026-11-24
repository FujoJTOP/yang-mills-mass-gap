# TOOL-ZEROCOUNT — argument-principle / Rouché zero-counting on `Z(τ)=w₀+Σfᵢ²e^{−sᵢτ}`: a READ of `s₁` and a zero-free region, but it CANNOT force the isolation gap `I2` (positivity admits `s₂=s₁`), and the zero COUNT is not a faithful level counter

> **VERDICT (no QED).** The direct non-integral attack on the isolation gap `I2`
> (`s₂−s₁≥Δ>0`) via the argument principle **does not close the blind spot.** Two rigorous
> pieces survive — the **exact single-atom zero ladder** (reads `s₁`) and a **Rouché zero-free
> right strip** (zeros live in a bounded strip) — but **(a)** the zero **count** is **not a faithful
> level counter** (a close second level adds *no* zero; an equal-weight continuum band gives a
> *constant* winding), and **(b)** **positivity admits `s₂=s₁`** (a positive measure with no atom at
> the support floor), so **no positivity constraint — hence no zero-counting built on it — can
> FORCE `I2`**. Zero-counting is a **READ** (of `s₁` and of where the zeros cannot be), not a
> **FORCE**. `I2 = H1` again. **Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric. `Z(τ)=w₀+Σᵢ fᵢ² e^{−sᵢτ}`, `w₀=|⟨O⟩|²>0`, `fᵢ²≥0` (Laplace transform of a positive KL measure `ν`, `sᵢ=mᵢ²`); `Z` is **entire**.

---

## §1 · What the argument principle can build `[T]`

**Rouché zero-free strip.** For `Re τ>0`, `|Σfᵢ²e^{−sᵢτ}| ≤ Σfᵢ²e^{−sᵢRe τ} < Σfᵢ²`, and if `w₀>Σfᵢ²` then `Z≠0` on `Re τ>0`; more generally `Z≠0` for `Re τ ≥ T₀` for some finite `T₀` (the `w₀` term dominates). Likewise `Re τ→−∞` the largest-`s` term dominates (`Z→f²_max e^{−s_maxτ}`, positive). **So the zeros of `Z` live in a bounded vertical strip** — a genuine `[T]` structural fact (useful: it bounds where zero-based extractions must look).

**Exact single-atom ladder.** `Z=w₀+f₁²e^{−s₁τ}=0 ⇒ e^{−s₁τ}=−w₀/f₁² ⇒ τ_n=(ln(f₁²/w₀) − i(2n+1)π)/s₁`. Zeros lie on **one vertical line**, `Im`-spacing `2π/s₁`; the winding in a box around `±π/s₁` is exactly 2. **This READS `s₁` (the gap) and `f₁²/w₀` (the weight).** Witness (a): spacing `6.2832=2π`, winding 2.

---

## §2 · Why the zero **count** does NOT give the isolation bound `[NG]`

The hoped-for bridge — "a bounded zero count around the lightest zero ⇒ `s₂−s₁≥Δ`" — **fails**, because the zero count is **not a faithful/monotone counter of spectral levels**:
- **A close second level adds NO zero.** Witness (d): same box `Re[−4,0.5]×Im[−14,14]` gives winding **4** for the single atom, **4** for two atoms with gap `0.05`, and **6** only for gap `0.60`. So near-degenerate levels are **invisible to the zero count**.
- **An equal-weight continuum band gives a CONSTANT winding.** Witness (d): `N` equal-weight levels on `[1,1+0.03(N−1)]` give winding `2` for `N=1,2,4,8,16` — the band collapses to one effective pole; the count neither grows nor bounds the level number.

So "count the zeros to bound the levels" is **not** a usable route: the data you'd need (the level multiplicities up to a near-degeneracy) is exactly **not encoded in the count**.

---

## §3 · The hard obstruction: positivity admits `s₂=s₁` `[NG]`

The cleanest and decisive obstruction is **logical**, needing no zeros:

> `dν = ρ(s) ds` on `[E₁,∞)`, `ρ(s)≥0`, is a **positive** measure, but the support floor `E₁` carries **no atom** — the first continuum level sits **at** `E₁`, so the isolation gap `s₂−s₁ = 0`. Hence **`ν≥0` does NOT imply `I2`.**

Since every positivity constraint used by a zero-counting argument is satisfied by this measure, **no positivity constraint can exclude it**, so **no positivity-based zero-counting can FORCE `I2`.** (This is the same atom-blind / lsc skeleton as `PNT-OPE-RHO`: positivity is compatible with no atom.)

---

## §4 · Even as a *reader*, the ladder is not a clean gap-reader `[V-num]`

Adding a second level does move the first rung, but in the crude Newton test the shift is **non-monotone in the gap** (witness (c): `|shift|=0.53, 0.69, 0.83, 0.62` for gaps `0.30,0.60,1.20,2.40` — Newton jumps rungs). So the zero-ladder robustly **reads `s₁`**, but does **not** robustly read `s₂−s₁`. (The clean gap-reader is the **approach-law** — `WI-PLATEAU-RIGOR` P3 — not the zero ladder.)

---

## §5 · Honest bottom line

- **The blind spot is NOT filled.** Zero-counting adds: an **exact read of `s₁`** and a **Rouché zero-free strip** (`[T]`), but it **cannot force `I2`** — positivity admits `s₂=s₁` (§3) — and the **zero count is not a level counter** (§2), nor is the ladder a clean gap-reader (§4).
- **What survives (real `[T]`):** the bounded-strip fact (zeros of `Z` lie in a finite strip; Rouché) and the exact single-atom ladder `2π/s₁`. These are **read** tools (like `A3`), not force tools.
- **Verdict:** the argument principle is a **READ-not-FORCE** coordinate of `I2=H1`; it does **not** fill the "non-integral force" blind spot. The invariant `I2=H1` stands.

*Witness:* `w7/tool_zerocount.py` (numpy only; `py -3.12 -I`, exit 0).
*Refs:* `WI-PLATEAU-RIGOR.md` (approach-law is the clean gap-reader), `PNT-OPE-RHO.md` (positivity
skeleton), `WI-ESCAPE.md` (A3 read-not-force).
