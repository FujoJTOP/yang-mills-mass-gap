# Blocked-SLC: gauge removal + Var term (refinement)

> CPB-attack probe. Verdict: **still-inconclusive, leaning dead.**

## Task 1 (gauge removal) — invalid as directed; the −0.256 is discarded

`PᵀHP` is the **pointwise** restriction of the fine Hessian. The coarse gauge image is
**not** a kernel direction of it: the thermal fine field is **not block-uniform**, so
coarse gauge transformations are not symmetries of the pointwise restriction. The
implementation confirmed this (gauge image not annihilated, relative residual ≈0.3;
the constant-ω combination failed to cancel), so the resulting "physical min −0.256" is a
**bug artifact, discarded**.

The ensemble-average `Ĥ` min **+0.043 (0/512 neg)** is therefore a **genuine, non-gauge
eigenvalue of the pointwise restriction** — but it is **not** the marginal Hessian.

## Task 2 (Var term) — uncomputed; a proxy

Within-block force variance / 16 ≈ **0.19** (mean `⟨|f̄|²⟩ ≈ 0.21`). Against the softest
E-part eigenvalue **+0.043**:

```
G'' ≲ 0.043 − 0.19 < 0     in the softest directions  →  leaning dead.
```

Caveat: this proxy is the **unconditional** within-block variance — an **upper bound** on
the conditional `Var(∂S|A')` — so the true `G''` could be **less negative**.

## Verdict

**Still-inconclusive, leaning dead.** The stiff directions have `G'' ≫ 0`; the softest
likely `G'' < 0`. The "+0.043 positive" is real for `PᵀHP` but is **not** the marginal —
the **Var term, not gauge**, is what decides, and it is still uncomputed exactly.

## Honest boundary

Gauge-image implementation was buggy; `Var` is an upper-bound proxy; L=4 (`L'=2`) is small.
Not a proof either way; the exact conditional `Var` (constrained sampling) is the only
clean decider.
