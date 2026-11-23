"""Conditional-variance decomposition of the coarse force.

Formula (ENTRY.md eq. (1)):
    G'' = E[d2_xx F | x] - Var(d_x F | x).

Here we measure how much of the *coarse* force f = block-average(dS_W) is
linearly explained by the *coarse* field xi = block-average(algebra coord),
at the SAME coarse link:

    f^a = sum_b M^{ab} xi^b + noise ,      R^2 per component a .

R^2 (same-link) is gauge-invariant: under a common adjoint rotation
f -> R f, xi -> R xi, the regression R^2 is unchanged.  (Neighbour predictors
would break this, so we keep the same-link regression.)

Residual fraction 1 - R^2 is the proxy for the variance term.  Compare b=1
(no blocking) vs b=2 (blocked): if the blocked force is *more* deterministic
given the field (R^2 up), the blocked field is more Gaussian-like.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lattice-cpb'))
import su3 as S  # noqa: E402


def field_coord(U):
    """xi^a = 0.5 * Im tr(T^a U)  (algebra coordinate, tr(T^a T^b)=2 delta)."""
    out = np.empty(U.shape[:5] + (8,))
    for mu in range(4):
        out[mu] = 0.5 * np.einsum('aij,...ji->...a', S.T_GEN, U[mu]).imag
    return out


def block_avg(X, b):
    """X: (D, L,L,L,L, C) -> block-average the four spatial axes."""
    D, L, C = X.shape[0], X.shape[1], X.shape[-1]
    n = L // b
    Y = X.reshape(D, n, b, n, b, n, b, n, b, C)
    return Y.mean(axis=(2, 4, 6, 8))


def r2_per_component(F, X):
    """F: (N,8) forces, X: (N,8) fields.  Returns mean R^2 over the 8 comps."""
    r2s = []
    Xd = np.concatenate([np.ones((X.shape[0], 1)), X], axis=1)
    for a in range(F.shape[1]):
        y = F[:, a]
        coef, *_ = np.linalg.lstsq(Xd, y, rcond=None)
        pred = Xd @ coef
        ss_res = np.sum((y - pred) ** 2)
        ss_tot = np.sum((y - y.mean()) ** 2)
        r2s.append(1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0)
    return float(np.mean(r2s)), r2s


def block_var_decomp(fmag, b):
    """Within-block vs between-block variance of the force magnitude |f|^2.
    fmag: (D, L,L,L,L).  Gauge-invariant (|f|^2 = sum_a f_a^2)."""
    D, L = fmag.shape[0], fmag.shape[1]
    n = L // b
    Y = fmag.reshape(D, n, b, n, b, n, b, n, b)
    bm = Y.mean(axis=(2, 4, 6, 8))                      # (D,n,n,n,n)
    bm_b = bm.reshape(D, n, 1, n, 1, n, 1, n, 1)
    within = ((Y - bm_b) ** 2).mean(axis=(2, 4, 6, 8))  # (D,n,n,n,n)
    return float(within.mean()), float(bm.var())


def run(beta, L, sweeps, therm, seed, b):
    rng = np.random.default_rng(seed)
    lat = S.Lattice(L, beta, rng, cold=True)
    for _ in range(therm):
        lat.sweep()

    Fs, Xs, fmag = [], [], []
    for _ in range(sweeps):
        lat.sweep()
        U = lat.U
        f = S.action_gradient(U, beta)      # (4,L,L,L,L,8)
        xi = field_coord(U)                 # (4,L,L,L,L,8)
        fb = block_avg(f, b) if b > 1 else f
        xb = block_avg(xi, b) if b > 1 else xi
        # flatten over (direction, site): (4*L'^4, 8)
        Fs.append(fb.reshape(-1, 8))
        Xs.append(xb.reshape(-1, 8))
        fmag.append(np.sum(f ** 2, axis=-1))   # (4,L,L,L,L)
    F = np.concatenate(Fs, axis=0)
    X = np.concatenate(Xs, axis=0)
    r2, per = r2_per_component(F, X)
    wv, bv = block_var_decomp(np.stack(fmag).reshape(-1, L, L, L, L), b)
    return r2, np.array(per), F.shape[0], wv, bv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=6)
    ap.add_argument('--sweeps', type=int, default=150)
    ap.add_argument('--therm', type=int, default=400)
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--betas', type=str, default='2.6,3.0,3.3')
    args = ap.parse_args()

    print(f"L={args.L} sweeps={args.sweeps} therm={args.therm}")
    print("beta   b   R2_same_link  within_var  between_var  within/(w+b)")
    for beta in [float(x) for x in args.betas.split(',')]:
        for b in (1, 2):
            r2, per, n, wv, bv = run(beta, args.L, args.sweeps, args.therm, args.seed, b)
            ratio = wv / (wv + bv) if (wv + bv) > 0 else 0.0
            print(f"{beta:5.2f}  {b}   {r2:+.4f}       {wv:.4f}     {bv:.4f}      {ratio:.4f}")
        print(flush=True)


if __name__ == '__main__':
    main()
