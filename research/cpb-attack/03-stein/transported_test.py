"""Validate the instrument: regress the FINE force on the 18 staple links'
algebra coordinates, PARALLEL-TRANSPORTED to the base site x.

If this gives high R^2 (unlike the raw-coordinate regression ~0), then the
instrument is right and the coarse-level question needs gauge-covariant
(transported) link blocking for b=2.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'lattice-cpb'))
import su3 as S  # noqa: E402


def xi(U):
    return 0.5 * np.einsum('aij,...ji->...a', S.T_GEN, U).imag


def ad(g, x):
    """Adjoint action: algebra coords of g (x_a T^a) g^{-1}."""
    M = np.einsum('...a,aij->...ij', x, S.T_GEN)
    M = g @ M @ S.adj(g)
    return 0.5 * np.einsum('aij,...ji->...a', S.T_GEN, M).imag


def feats(U, mu):
    """List of transported 8-vectors for the 18 staple links of direction mu."""
    xs = []
    for nu in range(4):
        if nu == mu:
            continue
        # forward plaquette at x: U_mu(x), U_nu(x+mu), U_mu(x+nu), U_nu(x)
        xs.append(xi(U[mu]))
        xs.append(ad(S.adj(U[mu]), xi(np.roll(U[nu], -1, axis=mu))))
        xs.append(ad(S.adj(U[nu]), xi(np.roll(U[mu], -1, axis=nu))))
        xs.append(xi(U[nu]))
        # backward plaquette at x-nu (relative to x)
        Un_m = np.roll(U[nu], 1, axis=nu)          # U_nu(x-nu)
        Um_m = np.roll(U[mu], 1, axis=nu)          # U_mu(x-nu)
        Un_mm = np.roll(np.roll(U[nu], 1, axis=nu), -1, axis=mu)  # U_nu(x-nu+mu)
        xs.append(ad(Un_m, xi(Un_m)))
        xs.append(ad(Un_m, xi(Um_m)))
        xs.append(ad(Un_m @ S.adj(Um_m), xi(Un_mm)))
    return xs


def r2(Y, X):
    Xd = np.concatenate([np.ones((X.shape[0], 1)), X], axis=1)
    rs = []
    for a in range(Y.shape[1]):
        y = Y[:, a]
        coef, *_ = np.linalg.lstsq(Xd, y, rcond=None)
        pred = Xd @ coef
        ss_res = np.sum((y - pred) ** 2)
        ss_tot = np.sum((y - y.mean()) ** 2)
        rs.append(1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0)
    return float(np.mean(rs))


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=6)
    ap.add_argument('--sweeps', type=int, default=8)
    ap.add_argument('--therm', type=int, default=250)
    ap.add_argument('--betas', type=str, default='6.0,3.0')
    args = ap.parse_args()
    for beta in [float(x) for x in args.betas.split(',')]:
        rng = np.random.default_rng(7)
        lat = S.Lattice(args.L, beta, rng, cold=True)
        for _ in range(args.therm):
            lat.sweep()
        Ys, Xs = [[] for _ in range(4)], [[] for _ in range(4)]
        for _ in range(args.sweeps):
            lat.sweep()
            U = lat.U
            f = S.action_gradient(U, beta)
            for mu in range(4):
                ff = feats(U, mu)
                X = np.concatenate([a.reshape(-1, 8) for a in ff], axis=1)
                Ys[mu].append(f[mu].reshape(-1, 8))
                Xs[mu].append(X)
        rs = []
        for mu in range(4):
            Y = np.concatenate(Ys[mu], 0)
            X = np.concatenate(Xs[mu], 0)
            rs.append(r2(Y, X))
        print(f"beta={beta}: R2(transported staple, b=1) = {np.mean(rs):+.4f}  per-mu={[round(r,3) for r in rs]}", flush=True)


if __name__ == '__main__':
    main()
