# STAR-NO-BZ — the "pure YM has no IR fixed point (no Banks–Zaks)" link: **[O]**, and it is *dynamical*, not a consequence of any exact property

> **VERDICT.** The (★) chain is `gapless ⟹ (★) conformal ⟹ β(g*)=0, g*>0 ⟹ [THIS LINK: no such g*] ⟹ contradiction`.
> This note pins the **second load-bearing link**: "pure YM has no IR fixed point." It is **`[T]` only
> perturbatively** (2-loop, indeed 3-/4-loop, `β<0` for all `g>0`), and **`[O]` nonperturbatively** — there is no
> rigorous proof that `β(g)` never crosses zero at strong coupling. Crucially, **no exact property forbids an
> IR fixed point**: the trace anomaly is *consistent* with `β(g*)=0` (it vanishes there), the a-theorem allows
> `a_IR<a_UV`, and unitarity allows a conformal IR. And the **conformal window** proves that "AF-UV + conformal-IR"
> is a **real phenomenon** for gauge theories — so "no BZ for pure YM" is a genuine **dynamical** statement (pure
> YM is below the window), not a kinematic theorem. **Net: (★) does not close the gap, because its no-BZ link is
> itself `[O]` and is essentially "pure YM is not in the conformal window" ≈ confining ≈ H1.** Clay 0%.

---

## §0 · Labels. [T] 真定理 · [C] 条件定理 · [O] 精确缺失步 · [NG] no-go · [V-num] numeric.

## §1 · What is actually proven

**Perturbative `[T]`.** For pure `SU(N)`: `β(g) = −β₀g³ − β₁g⁵ − …` with `β₀ = 11N/3`, `β₁ = 34N²/3`
(`>0`), so **`β(g) < 0` for all `g>0`** at 2 loops (and the 3-/4-loop coefficients keep it negative for
small/perturbative `g`). ⇒ **no perturbative zero.** But this is only valid for `g` small; the crossover is at
`α_s~O(1)`, where the perturbative series is unreliable.

**Why no exact property forbids an IR fixed point `[T]`/`[NG]`.**
- **Trace anomaly** `⟨θ^μ_μ⟩ = (β/2g)⟨tr G²⟩`: at `β(g*)=0` this is **`0`** — exactly what a CFT wants. The
  anomaly is *satisfied*, not contradicted.
- **a-theorem**: an IR CFT has `a_IR < a_UV` — *allowed* (the theorem is an inequality, not an exclusion).
- **Unitarity / KL positivity**: a conformal IR is as unitary as a gapped one.
⇒ **No exact identity excludes an IR fixed point.** The exclusion is dynamical.

**The conformal window is the decisive counterexample to any "AF ⇒ gapped" claim `[T]`.** `SU(N)` with `N_f`
fundamentals: `β₀ = 11N/3 − 2N_f/3`, `β₁ = 34N²/3 − (10N/3 + 2C_F)N_f`. For `N_f` near `11N/2` (`β₀>0` small,
`β₁<0`) the 2-loop β crosses zero — a **Banks–Zaks IR CFT**, i.e. an **asymptotically-free theory that is
gapless and conformal in the IR**. So "AF-UV + conformal-IR" is **real**, and the implication (★) is *consistent*
with it. Pure YM = the `N_f=0` corner, far **below** the window.

## §2 · The open step, pinned

> **Pure YM (`N_f=0`) is not in the conformal window — `β(g)` has no zero `g*>0` nonperturbatively.** This is
> `[O]`: believed (lattice shows confinement for all couplings; the window's lower edge is numerically known to
> be far above `N_f=0`), but **no rigorous nonperturbative proof** that the strong-coupling β never reaches zero.
> It is **not** a consequence of the trace anomaly / a-theorem / unitarity (those *allow* an IR CFT); it is a
> **dynamical** fact about pure YM's position relative to the window.

## §3 · Honest verdict on (★)

- **(★) ∧ no-BZ ⟹ pure YM gapped**, but the **no-BZ link is `[O]`**, and it says essentially **"pure YM is not in
  the conformal window" = "pure YM is confining" = `H1`.** So **(★) does not break the wall; it re-derives it**:
  the sharpest contradiction route reduces, at this link, to the open statement that pure YM is confining.
- The **only** way no-BZ could be `[T]` would be a nonperturbative sign theorem for `β` (none exists) or an exact
  obstruction to a gapless 1-form-symmetric 4D phase (also none exists) — both `[O]/[NG]`.

*Refs:* `IMPOSSIBLE-CONTRADICTION.md`, `IMPOSSIBLE-SYNTHESIS.md`, `SD-*`/`ROOF-*`. Witness `w7/star_no_bz.py`.
