# Radial Casimir solve — numerical attempt and the obstruction

> Fork probe. Verdict: the **conformal-cross-ratio** block is verified numerically; the **H²-variable**
> ("single ODE in `u`") reading is **not** a regular ODE — it is irregular at the OPE point. So a
> clean spectral solve of the bulk–defect radial block *in `u`* is not achieved.

## What was solved

- **(A)** 1D hypergeometric block `G(ξ)=ξ^Δ ₂F₁(Δ,Δ;2Δ;ξ)`: Chebyshev collocation (N=160) of
  `ξ(1−ξ)g''+[2Δ−(2Δ+1)ξ]g'−Δ²g=0`, `g(0)=1`, `g'(0)=Δ/2`; verified vs `scipy.hyp2f1`.
- **(B)** H²-Laplacian eigen-equation in the K₃-invariant variable `z=cosh d` (`u=(z−1)/(z+1)`):
  `(z²−1)f'' + z f' ± Δ(Δ−1) f = 0`, `solve_ivp` shooting from `z=1.01`.

## Numerical block (verified; reproduces `NUMERICS.md` digits)

`G(ξ)=ξ^Δ ₂F₁(Δ,Δ;2Δ;ξ)` at `ξ = 0.1, 0.3, 0.5, 0.7`:

| `Δ` | 0.1 | 0.3 | 0.5 | 0.7 |
|---|---|---|---|---|
| 2 | 0.011099 | 0.126948 | 0.476649 | 1.415697 |
| 3 | 0.001169 | 0.045170 | 0.327400 | 1.659044 |

Product block `𝔽_j(ξ,θ)=G(ξ)·P_j(cos θ)`. At `cos θ=0.5`: `P_0=1, P_1=0.5, P_2=−0.125`
(matches `PDE-SOLUTION.md`), e.g. `(Δ,j)=(2,0):0.12695`, `(2,1):0.06347`, `(2,2):−0.01587` at `ξ=0.3`.

## Accuracy

`Δ=2`: `max|spectral−hyp2f1| = 4.6×10⁻⁹`. `Δ=1,3`: the Chebyshev collocation is numerically
unstable (absolute error dominated by large values near `ξ→1`); the tabulated values were verified
against `scipy.hyp2f1` instead (exact to machine precision, matching `NUMERICS.md`).

## Obstruction (the honest finding)

**(B) is irregular at `z=1`.** Substituting `z=1+w`: `2w f''+f' ± Δ(Δ−1)f = 0`; for `f=w^α`,
`α(2α−1)w^{α−1}` vs `λ w^α` → no power-law balance for `λ≠0`. So the block is **not** a regular
solution of the naive H²-Laplacian eigen-equation. The `K3-RESOLUTION.md` "single ODE in `u`" is
**not** this ODE. The block is regular hypergeometric in the 1D **conformal** cross-ratio `ξ`
(object (A)); but the map `ξ↔u` (H² distance) and the correct **conformal-Casimir** operator in `u`
are **not** derived.

## Honest boundary

The product `𝔽=G(ξ)P_j(cos θ)` is verified numerically. The radial-variable (`u`) ODE in its
correct form is **open**; the obstruction is the mismatch between the naive H² Laplacian and the
conformal Casimir. Deriving the latter in `u` is the missing step.
