"""Coarse-staple neighbour regression (gauge-invariant via Coulomb gauge fixing).

Question: is the coarse force f' = block-average(dS_W) determined by the coarse
field A' (block-average algebra coord) at the *same link* (R^2 ~ 0, known) or by
the *plaquette-neighbour* coarse links?  Distinguishes 'local but off-site'
(local SLC salvageable) from 'genuinely nonlocal' (operator target).

Gauge invariance: we first Coulomb-gauge-fix each configuration, so all algebra
coordinates live in a common frame; R^2 is then a function of the gauge orbit.
(Residual Gribov-copy caveat noted in the report.)
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'lattice-cpb'))
import su3 as S  # noqa: E402


def coulomb_fix(U, n_iter=40):
    """Iterative Coulomb gauge fixing (max sum_mu Re tr U_mu)."""
    U = np.array(U, dtype=np.complex128, copy=True)
    L = U.shape[1]
    for _ in range(n_iter):
        M = np.zeros((L, L, L, L, 3, 3), dtype=np.complex128)
        for mu in range(4):
            M += U[mu] + S.adj(np.roll(U[mu], 1, axis=mu))
        g = S.project_su3(M)
        det = np.linalg.det(g)
        phase = np.exp(1j * np.angle(det) / 3.0)   # |det|=1 for the polar part
        g = g / phase[..., None, None]             # force det = 1
        for mu in range(4):
            U[mu] = g @ U[mu] @ S.adj(np.roll(g, -1, axis=mu))
    return U


def field_coord(U):
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


def design(fb, xb):
    """Build (Y, X) pooled over directions and sites.
    Predictors = the algebra coords of the 18 staple links of link (mu,x):
    for nu != mu:  U_nu(x), U_nu(x+mu), U_nu(x-nu), U_nu(x-nu+mu),
                   U_mu(x+nu), U_mu(x-nu)  (all transported; gauge fixed)."""
    Ys, Xs = [], []
    for mu in range(4):
        y = fb[mu].reshape(-1, 8)
        feats = [xb[mu].reshape(-1, 8)]                       # same link
        for nu in range(4):
            if nu == mu:
                continue
            feats.append(xb[nu].reshape(-1, 8))                           # U_nu(x)
            feats.append(np.roll(xb[nu], -1, axis=mu).reshape(-1, 8))     # U_nu(x+mu)
            feats.append(np.roll(xb[nu], 1, axis=nu).reshape(-1, 8))      # U_nu(x-nu)
            feats.append(np.roll(np.roll(xb[nu], 1, axis=nu), -1, axis=mu)
                         .reshape(-1, 8))                                 # U_nu(x-nu+mu)
            feats.append(np.roll(xb[mu], -1, axis=nu).reshape(-1, 8))     # U_mu(x+nu)
            feats.append(np.roll(xb[mu], 1, axis=nu).reshape(-1, 8))      # U_mu(x-nu)
        Xf = np.concatenate(feats, axis=1)                                # (N, 152)
        Ys.append(y)
        Xs.append(Xf)
    return np.concatenate(Ys, 0), np.concatenate(Xs, 0)


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


def run(beta, L, sweeps, therm, seed, b):
    rng = np.random.default_rng(seed)
    lat = S.Lattice(L, beta, rng, cold=True)
    for _ in range(therm):
        lat.sweep()
    Ys, Xs, Ysame, Xsame = [], [], [], []
    for _ in range(sweeps):
        lat.sweep()
        U = coulomb_fix(lat.U)
        f = S.action_gradient(U, beta)
        xi = field_coord(U)
        fb, xb = block_avg(f, b), block_avg(xi, b)
        y, x = design(fb, xb)
        Ys.append(y); Xs.append(x)
        # same-link only
        Ysame.append(np.concatenate([fb[mu].reshape(-1, 8) for mu in range(4)], 0))
        Xsame.append(np.concatenate([xb[mu].reshape(-1, 8) for mu in range(4)], 0))
    Y = np.concatenate(Ys, 0); X = np.concatenate(Xs, 0)
    r2_nb = r2(Y, X)
    r2_same = r2(np.concatenate(Ysame, 0), np.concatenate(Xsame, 0))
    return r2_nb, r2_same, Y.shape[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=6)
    ap.add_argument('--sweeps', type=int, default=25)
    ap.add_argument('--therm', type=int, default=300)
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--betas', type=str, default='2.6,3.0')
    ap.add_argument('--b', type=int, default=2)
    args = ap.parse_args()
    print(f"L={args.L} sweeps={args.sweeps} therm={args.therm} b={args.b}")
    print("beta   R2_same_link  R2_same+neighbours  n_samples")
    for beta in [float(x) for x in args.betas.split(',')]:
        r2nb, r2s, n = run(beta, args.L, args.sweeps, args.therm, args.seed, args.b)
        print(f"{beta:5.2f}   {r2s:+.4f}       {r2nb:+.4f}          {n}", flush=True)


if __name__ == '__main__':
    main()
