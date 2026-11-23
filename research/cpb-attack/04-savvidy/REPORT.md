# Method 4 — effective SLC / Savvidy negative modes

> CPB-attack probe. Verdict: **dead (cheap restricted-SLC)**; **needs-X (blocked-scale SLC)**.

## Data (L=6, per-link Hessian from Gibbs sampling — ensemble average = measure-weighted)

| β | frac(min<0) | E[neg] | E[pos] | E[neg]/E[pos] | per-config mean |
|---|---|---|---|---|---|
| 2.5 | 0.493 | 0.131 | 0.139 | **0.94** | +0.008 |
| 2.6 | 0.465 | 0.124 | 0.158 | 0.79 | +0.034 |
| 3.0 | 0.362 | 0.101 | 0.249 | 0.41 | +0.148 |

**Reading:** at `β≈2.5` the negative curvature is ≈ **94%** of the positive curvature,
**49% of links have a negative eigenvalue**, and the per-link mean curvature ≈ 0. The
premise "negative modes are a small-measure defect" is **false** — it is half the links,
half the curvature.

## Literature

The only preservation theorem (Saumard–Wellner Thm 3.2.2 / Brascamp–Lieb) requires
**global** SLC. No theorem for "SLC outside a small set" was found — it would need to be
built.

## Verdict

- **dead** (cheap restricted-SLC): the defect occupies half the links / half the
  curvature, so there is no "small set" to exclude;
- **needs-X** (the one survivor): **blocked-scale** SLC — the *effective* action at the
  crossover might be SLC even though the microscopic one is not. Untested, unknown.

## Honest boundary

`link_hessian_min` is the **per-link** Hessian, not the **collective / blocked** Hessian;
the collective Savvidy mode may differ. This is a decision input, not a proof.

Sources: [Saumard–Wellner arXiv:1404.5886](https://arxiv.org/html/1404.5886v1);
[review §3 preservation theorems](https://sl2x.aimath.org/development/collectedworks/htmlpaper/1404_5886/section3.html).
