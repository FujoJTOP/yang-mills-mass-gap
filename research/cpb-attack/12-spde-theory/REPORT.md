# Method 5b — 4D critical gauge SPDE: the minimal theory

> CPB-attack probe. Verdict: **NEEDS-X** (structural, not a missing constant).

## 1. The exact renormalisation obstruction

Regularity structures / paracontrolled calculus rest on **subcriticality**: for a
nonlinearity of degree `k` in a field of regularity `α`, one needs all tree weights
positive, so the number of relevant trees is **finite** (Hairer's "subcritical" /
GIP degree condition). In `d=3`, `[g]=1/2>0` (super-renormalisable): finitely many
trees, uniform BPHZ bound. In `d=4`, `[g]=0` (marginal): the tree weights are **zero**,
the set of trees of nonnegative weight is **infinite**, and the BPHZ bound is **not
uniform** in the cutoff. Neither framework has a marginal-case theorem. 4D Φ⁴ is the
only solved marginal case and it is **trivial** (Aizenman–Duminil-Copin 2021) — no
template for a non-trivial marginal gauge theory.

## 2. Minimal "critical gauge regularity structure"

- **Scaling:** build the **running coupling** into the scaling — a `t`-dependent norm
  where the marginal coupling is absorbed as a **logarithmic** correction to tree
  weights (`weight 0` → `weight −log`). A **log-critical** structure, not a fixed-point one.
- **Counterterms:** log-divergent (not power-divergent); the RG is a **flow equation**
  (Polchinski/Wetterich-type) inside the regularity-structure setting, not a finite BPHZ
  subtraction.
- **Gauge invariance:** open — the gauge-covariant frame (parallel transport) used in 3D
  must survive the logarithmic counterterms; compatible iff the log-counterterms are
  gauge-covariant (expected via the Ward identity), but **not established**.

## 3. Literature (2024–2026)

No critical/marginal gauge SPDE theorem exists. Solved marginal cases are trivial (Φ⁴).
Critical-SPDE work is for **scalar** marginal (2D Liouville-type). The 3D YM construction
(CCHS) does not extend; the gap to 4D is structural.

## 4. Constructibility and the gap step

Constructibility: **needs-X** — a new log-critical gauge regularity structure with a
flow-equation renormalisation, not available. The gap step (Zegarlinski fictitious-time
LSI ⇒ Euclidean decay) requires supplying the LSI, which is the **global summability**
object blocked by **OB-1** — not free.

## Verdict

**NEEDS-X.** First test: write the log-critical tree-weight accounting for 4D YM and
check whether the log-counterterms close gauge-covariantly.

## Honest boundary

Reach assessment; says nothing about whether such a theory exists.
