# Blocked-SLC: exact conditional Var (final proxy)

> CPB-attack probe. Verdict: **DEAD-lean, with a clearer margin.**

## Method

Exact constrained sampling of `A⊥ | A'` on SU(3) is infeasible without a constrained MCMC,
so we computed the ensemble law-of-total-variance proxy
`E_{A'}[Var(g|A')] = Var(g) − Var(E[g|A'])`, with `g = Pᵀ∇S` (coarse gradient),
`A' = Pᵀ log-link` (coarse field), and the conditional-mean variance estimated by **linear
regression** on `A'`, in the **softest eigen-direction** `v` of `Ĥ = E[PᵀHP]`. Script
`cond_var.py`; L=4, β=2.5, 40 configs (3 for `Ĥ`).

## Result (along the softest `Ĥ` direction)

- `λ_E = −0.0005` (3-config `Ĥ`; the stable 20-config reference is **+0.043**).
- `Var(g·v) = 0.209`; linear-regression residual `Var = 0.208`; `R² = 0.004`.
- **`G''_soft ≈ +0.043 − 0.208 = −0.165`.**

## Caveats

(i) The residual `Var` is an **upper bound** (linear model; `R²≈0`, so it ≈ the total force
variance). (ii) `A'` uses the absolute principal log, `g` the left-variation gradient — mixed
coordinates (but the total force variance 0.209 is coordinate-robust). (iii) No exact
constrained sampling.

## Verdict

**DEAD-lean, with a clearer margin.** The force variance along the softest direction (0.209)
is ~**5×** the softest E-part (0.043). To flip the sign, a very nonlinear conditional mean
would have to explain **>80%** of the force variance, which the linear regression
(`R²=0.004`) gives **no support** for. Leaning dead — **not a clean proof**.

## Honest boundary

The exact conditional `Var` via constrained sampling remains the only clean decider. Scripts
in this directory.
