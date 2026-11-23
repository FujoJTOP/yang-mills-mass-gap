"""Measure the local third connected cumulant of the CPB-W perturbation density.

V = S - S_G  with S_G the quadratic (Gaussian) part of the Wilson action:
    V_x = beta * sum_{p in x} [ (1 - Re tr U_p / 3) - tr(A_p^2)/6 ],  A_p = -i log U_p
(the second term is the exact quadratic reference; sum over the 6 plaquettes
based at x, so that sum_x V_x = S_W - S_G).

Per configuration we store the two local densities S_x and Q_x with
    V_x(alpha) = S_x - alpha Q_x        (alpha = 1 is the canonical V)
and the blocked (b=2,3) unweighted block means are formed in analyze.py.

Outputs an npz with the per-configuration fields.
"""
import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', 'lattice-cpb')))
import su3_numba as N  # noqa: E402


def adj(A):
    return np.conj(np.swapaxes(A, -1, -2))


def eig_angles_batch(F):
    """Eigenvalue phases of a batch of 3x3 unitary matrices, principal branch.

    F shape (M,3,3) -> angles (M,3) in (-pi,pi].  Exact to ~1e-15.
    """
    x = np.trace(F, axis1=-2, axis2=-1)
    y = np.trace(F @ F, axis1=-2, axis2=-1)
    s2 = 0.5 * (x * x - y)
    p = s2 - x * x / 3.0
    q = -2 * x ** 3 / 27.0 + x * s2 / 3.0 - 1.0
    disc = q * q + 4.0 * (p / 3.0) ** 3
    sq = np.sqrt(disc.astype(complex))
    t1 = 0.5 * (-q + sq)
    t2 = 0.5 * (-q - sq)
    u = np.power(t1.astype(complex), 1.0 / 3.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        v = np.where(np.abs(u) > 1e-14, (-p / 3.0) / u,
                     np.power(t2.astype(complex), 1.0 / 3.0))
    om = np.exp(2j * np.pi / 3.0)
    ws = np.stack([u + v, om * u + om * om * v, om * om * u + om * v], axis=-1)
    lam = ws + (x / 3.0)[..., None]
    return np.angle(lam)


def local_densities(U, beta):
    """Return S_x, Q_x arrays of shape (L,L,L,L) (site-based action densities)."""
    L = U.shape[1]
    ps, qs = [], []
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Up = (U[mu] @ np.roll(U[nu], -1, axis=mu)
                  @ adj(np.roll(U[mu], -1, axis=nu)) @ adj(U[nu]))
            F = Up.reshape(-1, 3, 3)
            ang = eig_angles_batch(F)
            qs.append(np.sum(ang * ang, axis=-1) / 6.0)
            ps.append(1.0 - np.real(np.trace(F, axis1=-2, axis2=-1)) / 3.0)
    ps = np.array(ps); qs = np.array(qs)          # (6,V)
    Sx = beta * np.sum(ps, axis=0).reshape(L, L, L, L)
    Qx = beta * np.sum(qs, axis=0).reshape(L, L, L, L)
    return Sx, Qx, ps, qs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=6)
    ap.add_argument('--beta', type=float, required=True)
    ap.add_argument('--sweeps', type=int, default=8000)
    ap.add_argument('--therm', type=int, default=1500)
    ap.add_argument('--meas', type=int, default=2)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', type=str, required=True)
    args = ap.parse_args()

    L = args.L
    beta = args.beta
    kappa = beta / 3.0
    U = np.zeros((4, L, L, L, L, 3, 3), dtype=np.complex128)
    for mu in range(4):
        U[mu][...] = np.eye(3)
    np.random.seed(args.seed)

    t0 = time.time()
    for k in range(args.therm):
        N.sweep_par(U, L, kappa)
    print(f'therm {args.therm} in {time.time()-t0:.1f}s  P={N.plaquette_mean(U,L):.5f}',
          flush=True)

    nmeas = args.sweeps // args.meas
    Sx = np.empty((nmeas, L, L, L, L), dtype=np.float32)
    Qx = np.empty((nmeas, L, L, L, L), dtype=np.float32)
    plaq = np.empty(nmeas)
    # per-plaquette accumulators: V_p(alpha)=beta(p-alpha q), 6*L^4 values/config
    AL = (0.0, 0.25, 0.5, 0.6, 0.75, 1.0, 1.25, 1.5)
    pmom = np.empty((nmeas, len(AL), 4))      # raw sums of V_p^k, k=1..4
    pmax = np.empty(nmeas)                    # max |V_p| at alpha=1
    f3p = np.empty(nmeas)                     # per-config mean cube of field p_p
    fmom = np.empty((nmeas, 4))               # raw moments of field p_p
    t0 = time.time()
    for i in range(nmeas):
        for _ in range(args.meas):
            N.sweep_par(U, L, kappa)
        s, q, ps, qs = local_densities(U, beta)     # ps,qs shape (6,V)
        Sx[i] = s
        Qx[i] = q
        plaq[i] = N.plaquette_mean(U, L)
        pflat = ps.reshape(-1)
        qflat = qs.reshape(-1)
        for k, a in enumerate(AL):
            v = beta * (pflat - a * qflat)
            pmom[i, k] = [v.sum(), (v ** 2).sum(), (v ** 3).sum(), (v ** 4).sum()]
        pmax[i] = np.abs(beta * (pflat - qflat)).max()
        f3p[i] = np.mean(pflat ** 3)
        fmom[i] = [pflat.sum(), (pflat ** 2).sum(), (pflat ** 3).sum(),
                   (pflat ** 4).sum()]
        if (i + 1) % max(1, nmeas // 8) == 0:
            print(f'  [{i+1}/{nmeas}] {time.time()-t0:.0f}s  P={plaq[i]:.5f}', flush=True)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    np.savez(args.out, L=L, beta=beta, seed=args.seed, sweeps=args.sweeps,
             therm=args.therm, meas=args.meas, Sx=Sx, Qx=Qx, plaq=plaq,
             pmom=pmom, pmax=pmax, f3p=f3p, fmom=fmom, alphas=np.array(AL))
    print('saved', args.out, f'({time.time()-t0:.0f}s)')


if __name__ == '__main__':
    main()
