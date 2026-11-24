# CERT-CHI-CONDENSATE — the rank-one Savvidy form floors the DEFECT `D` renormalon-free, but NOT `S_bare`: `S_bare` is the `u=2`-ambiguous condensate = `H1`

> **VERDICT (no QED, Clay 0%).** The Savvidy-cell rank-one split `G″ = S − D`,
> `D = D₊ + D_Sav` (`D₊≤0` by Brascamp–Lieb, `D_Sav = g²Var(φ²) ≤ 2g²/λ₄`), gives a **rigorous,
> renormalon-free, `c`-independent bound on the DEFECT `D`** — a real theorem. **But it does NOT give a
> renormalon-free uniform floor on `S_bare`.** The floor `G″≥λ` reduces to `S_bare ≥ λ + 2g²/λ₄`, and
> `S_bare = ⟨tr F²⟩` carries the `u=2` IR renormalon: `a_n ~ Γ(n+1+2b)`, Borel pole at `u=2`, ambiguity
> `~Λ⁴` = the condensate's own magnitude (`ROOF-CONDENSATE-FLOOR`). The rank-one form is a **one-sided
> constant bound** (`D ≤ const`), not a **subtraction** that cancels `S`'s `u=2` pole. **The constant is still
> the ambiguous condensate = `H1`.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric.

---

## §1 · The defect IS floored (real, renormalon-free)

`[T]` From `V = −½cφ² + ¼λ₄φ⁴` (`c=gB`), the saddle gives `Var(φ²) = 2/λ₄ + 2/(gB)² + O((gB)⁻⁴) → 2/λ₄`,
**independent of the negative curvature `c`** (witness: `c = 10…160 → 2.0221, 2.0051, 2.0013, 2.0003, 2.0001`).
Hence `D_Sav = g²Var(φ²) ≤ 2g²/λ₄` — bounded, sign-positive, `c`-independent. The negative-mode contribution to
`D` is the rank-one term, and its `k→0` overlap vanishes by scale separation in the deep IR (`D→0`). **This half
is genuinely renormalon-free.**

## §2 · The curvature floor is NOT (it is the condensate)

`[T]` `G″ = S − D ≥ S − 2g²/λ₄`, so `G″ ≥ λ ⟺ S_bare ≥ λ + 2g²/λ₄`. **The floor lives entirely on `S_bare`.**

`[NG]/[T]` `S_bare = ⟨tr F²⟩` (the bare coarse curvature) is the condensate. Its perturbative series is the
`u=2` IR renormalon: `a_n ~ Γ(n+1+2b)` (`b=2`), so `a_n/a_{n−1} ~ n` grows without bound — Borel pole at
`u=2`, lateral-sum ambiguity `~Λ⁴` = the condensate's own size (witness Part 3). The rank-one form bounds `D`
(`D ≤ const`); **a one-sided inequality with a constant does not subtract off `S`'s `u=2` pole**, so the floor
on `S_bare` is not made renormalon-free. Exact same ambiguity as `ROOF-CONDENSATE-FLOOR`.

## §3 · The cancellation that would be needed, and why the rank-one form cannot supply it

`[O]` The `u=2` renormalon *does* cancel in the **exact** difference `G″ = S − D` (the physical `1/χ` is a
renormalon-free scale). But the rank-one form supplies only `D ≤ const` — **a bound, not `D` exactly**. Exact
`D` = the conditioned fiber measure (Brascamp–Lieb defect as a functional of the true measure) = **`H1`**.
So the Savvidy structure **isolates and bounds the defect** (real, renormalon-free) but **leaves the curvature
floor as the ambiguous condensate**.

---

## §4 · Honest bottom line

- **Renormalon-free floor on `S_bare`: NO.** The rank-one form bounds `D` (true, valuable), but `S_bare` stays
  `u=2`-ambiguous (`~Λ⁴`).
- **The constant is the condensate = `H1`** — same wall as `ROOF-CONDENSATE-FLOOR`, reached here via
  `G″ = S − D ≥ S − 2g²/λ₄`.
- **What is genuinely gained:** a rigorous, `c`-independent, sign-positive bound on the soft defect `D`
  (`≤ 2g²/λ₄`), which certifies "the negative mode does not blow up `D`" — the defect half of the Savvidy
  sharpening, *not* the curvature half.
- **Clay = 0%.** No QED.

*Witness:* `w7/cert_chi_condensate.py` (runs `py -3.12 -I`, numpy+scipy).
*Refs:* `BLAST-SAVVIDY-CLOSE.md`, `SAVVIDY-DTERM.md`, `SAVVIDY-VAR-PHI2.md`, `SAVVIDY-NEGMODE-EXACT.md`,
`ROOF-CONDENSATE-FLOOR.md`, `ONLY-EXIT-SYNTHESIS.md`.
