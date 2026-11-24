# BLAST-EXCLUDE-BUMP — trying to exclude the atomless "flat bump" with a YM-specific property: every candidate reduces to `H1`. The blast FAILS; ONLY-EXIT is confirmed.

> **VERDICT (no QED).** The ONLY-EXIT adversary is a flat bump `ρ_h = f₁²/h` on `[m₁², m₁²+h]` (same mass,
> floor, width `h→0`) with `f₁²=0`, matching all finite moments of the atom measure. To blast `H1` one must
> exclude it by a YM-specific property. **I tried the four candidates; every one either fails to exclude it or
> reduces to `H1`.** Witness: the bump matches `M₀, χ, A` to `O(h)` (`|Δχ|,|ΔA| → 0` as `h→0`) while `f₁²=0`
> vs the atom's `0.5` — so no finite integral property separates them. **The blast confirms the wall.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric.

## §1 · The candidates, and why each fails

1. **Reflection positivity (RP) on the full field algebra `[NG]`.** RP is stronger than 0⁺⁺-measure positivity —
   but the 0⁺⁺ KL measure of **any** positive `ν` is realizable by an OS theory: a generalized-free / `:φ²:`
   -type construction gives a **gapped but atomless** 0⁺⁺ measure (`:φ²:` has `f₁²=0`). So RP alone does **not**
   exclude the bump — an atomless OS theory exists (the free-field acid adversary). To exclude it you must
   exclude the atomless *realization with YM's exact `χ, A`* — a non-perturbative existence question = `H1`.

2. **OPE / operator content `[NG]`.** Every OPE coefficient of `tr F² × tr F²` is a **moment** `M₋ₖ`; the bump
   matches them to `O(h)` (witness). Operator-specific structure is **linear in ν** (`FALSIFY-OPE`) ⇒ no
   nonlinear constraint.

3. **Exact Ward identity beyond `κχ=1−η` `[NG]`.** Ward is the tautology `η:=1−κχ`; the anomaly `A` is a moment.
   A flat bump of the same `χ, A` satisfies both.

4. **A nonlinear functional nonzero on the bump, zero on YM (or vice versa) `[NG]`.** Any functional that is a
   continuous function of the (matched) finite moments is blind; a **non-continuous** one distinguishing the
   bump is the singular-kernel/weak-* route, which detects *presence* only and needs the pole location + `ε→0`
   (= A3) — `TOOL-HARDY`, `FALSIFY-TOPOLOGICAL`.

## §2 · The one genuine cross-channel bite — and it is `H1`

The trace anomaly `θ = (β/2g) tr F²` ties the 0⁺⁺ channel to `tr T_μν`; the 2⁺⁺ (traceless `T`) is a
**separate** measure, and the dilatation Ward relates the 0⁺⁺ and 2⁺⁺ spectral densities. So a 0⁺⁺ bump
**inconsistent with the 2⁺⁺ spectrum** is excludable — **but that needs the 2⁺⁺ mass/gap**, i.e. the gap in
another channel. It does **not** give a 0⁺⁺-only exclusion; it relocates the demand to `H1` (2⁺⁺ isolation).

## §3 · Verdict

`[NG]` — no property checkable from finite `tr F²` data excludes the bump. Excluding it ≡ proving the floor
is a **discrete atom** consistent with the full operator algebra = **isolation `I2` = `H1`**. The blast fails;
**ONLY-EXIT is confirmed**: the single exit is `H1`. Clay 0%; no QED.

*Files:* `w7/blast_exclude_bump.py` (runs clean, `py -3.12 -I`).
*Refs:* `ONLY-EXIT-SYNTHESIS.md`, `FALSIFY-OPE.md`, `FALSIFY-TOPOLOGICAL.md`, `TOOL-HARDY.md`, `PROOF-P3-UNCONDITIONAL.md`.
