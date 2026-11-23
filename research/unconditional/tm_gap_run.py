"""
DIRECT transfer-matrix spectral gap measurement, SU(2) pure gauge, 4D.

Measures the Euclidean-time connected correlator of the lightest 0++ operator
(APE-smeared magnetic / tr F^2 density summed over each time slice).  The
exponential decay rate of C(t) = <O(t) O(0)>_c is the transfer-matrix gap
Delta = m_{0++} (lightest glueball above the vacuum), in lattice units.

Conventions (w4_blocked.py / su3.py):
  Wilson action  S = beta * sum_p (1 - 1/2 Re tr U_p),  beta = 4/g^2.
  SU(2) exact heat bath, kappa = beta/2.
  Operator on a time slice t:
     O(t) = sum_{x,y,z} sum_{1<=i<j<=3} (1/2) Re tr U_ij(x)
  Spatial links APE-smeared (alpha, n steps) before measurement; temporal
  links untouched so the transfer matrix is the genuine Euclidean one.

Output: npz with O[level] arrays shape (nconf, L) and the plaquette series.

usage:
  python tm_gap_run.py --L 6 --beta 2.5 --ntherm 800 --gapevery 4 --nconf 8000 \
      --out tm_gap_data/L6_b2.5.npz
"""
import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
import w4_blocked as W          # SU(2) sampler: sample2, staple_sum2, random_su2
from su3 import su2_heatbath, project_su3, adj  # generic, work for 2x2


def ape_spatial2(U, alpha=0.5, n=1):
    """APE-smear only the three spatial link directions (mu=1,2,3); 2x2 SU(2)."""
    U = U.copy()
    for _ in range(n):
        old = U.copy()
        for i in (1, 2, 3):
            S = np.zeros_like(U[i])
            for j in (1, 2, 3):
                if j == i:
                    continue
                # forward staple: U_j(x) U_i(x+j) U_j^dag(x+i)
                S += old[j] @ np.roll(old[i], -1, axis=j) @ adj(np.roll(old[j], -1, axis=i))
                # backward staple: U_j^dag(x-j) U_i(x-j) U_j(x-j+i)
                S += (adj(np.roll(old[j], 1, axis=j)) @ np.roll(old[i], 1, axis=j)
                      @ np.roll(np.roll(old[j], 1, axis=j), -1, axis=i))
            U[i] = project_su3((1.0 - alpha) * old[i] + alpha / 6.0 * S)
    return U


def magnetic_slice2(U):
    """O(t) = sum over spatial volume of (1/2)Re tr of the three spatial planes.
    Returns array shape (L,) indexed by time slice."""
    s = np.zeros(U.shape[1:5], dtype=np.float64)   # (t,x,y,z)
    for i in (1, 2, 3):
        for j in range(i + 1, 4):
            Up = (U[i] @ np.roll(U[j], -1, axis=i)
                  @ adj(np.roll(U[i], -1, axis=j)) @ adj(U[j]))
            s += 0.5 * np.real(np.trace(Up, axis1=-2, axis2=-1))
    return s.sum(axis=(1, 2, 3))                    # (t,)


def mean_plaq2(U):
    tot = 0.0
    n = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Up = (U[mu] @ np.roll(U[nu], -1, axis=mu)
                  @ adj(np.roll(U[mu], -1, axis=nu)) @ adj(U[nu]))
            tot += 0.5 * np.real(np.trace(Up, axis1=-2, axis2=-1)).mean()
            n += 1
    return tot / n


def run(L, beta, ntherm, gapevery, nconf, seed, alpha, levels, out):
    rng = np.random.default_rng(seed)
    ns = L ** 4
    grids = np.meshgrid(*[np.arange(L)] * 4, indexing='ij')
    par = sum(grids) % 2
    idxp = [np.nonzero(par.reshape(-1) == c)[0] for c in (0, 1)]
    kappa = beta / 2.0

    def sweep(U):
        Uf = U.reshape(4, ns, 2, 2)
        for c in (0, 1):
            ids = idxp[c]
            for mu in range(4):
                S = W.staple_sum2(U, mu).reshape(ns, 2, 2)[ids]
                A = Uf[mu, ids]
                M = A @ S
                v = su2_heatbath(kappa, M, 0, 1, rng)
                Uf[mu, ids] = v @ A

    U = W.random_su2(rng, (4, L, L, L, L))
    t0 = time.time()
    for _ in range(ntherm):
        sweep(U)
    print(f"[L={L} beta={beta}] thermalized {ntherm} sweeps in {time.time()-t0:.1f}s  "
          f"P={mean_plaq2(U):.5f}", flush=True)

    levels = sorted(levels)
    O = {lv: np.empty((nconf, L)) for lv in levels}
    plaq = np.empty(nconf)
    t0 = time.time()
    for k in range(nconf):
        for _ in range(gapevery):
            sweep(U)
        Us = U
        prev = 0
        for lv in levels:
            if lv > prev:
                Us = ape_spatial2(Us, alpha, lv - prev)
                prev = lv
            O[lv][k] = magnetic_slice2(Us)
        plaq[k] = mean_plaq2(U)
        if (k + 1) % max(1, nconf // 10) == 0:
            el = time.time() - t0
            print(f"  [{k+1}/{nconf}] P={plaq[k]:.5f}  {el:.0f}s "
                  f"(eta {el/(k+1)*(nconf-k-1):.0f}s)", flush=True)

    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    np.savez(out, L=L, beta=beta, alpha=alpha, levels=np.array(levels),
             gapevery=gapevery, nconf=nconf, plaq=plaq,
             **{f'O{lv}': O[lv] for lv in levels})
    print(f"saved {out}  total {time.time()-t0:.1f}s", flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, required=True)
    ap.add_argument('--beta', type=float, required=True)
    ap.add_argument('--ntherm', type=int, default=800)
    ap.add_argument('--gapevery', type=int, default=4)
    ap.add_argument('--nconf', type=int, default=8000)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--alpha', type=float, default=0.5)
    ap.add_argument('--levels', type=str, default='0,4,8')
    ap.add_argument('--out', type=str, required=True)
    args = ap.parse_args()
    lv = [int(x) for x in args.levels.split(',') if x.strip()]
    run(args.L, args.beta, args.ntherm, args.gapevery, args.nconf,
        args.seed, args.alpha, lv, args.out)
