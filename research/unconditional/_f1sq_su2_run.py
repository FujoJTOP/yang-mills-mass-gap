"""SU(2) connected 0++ (tr F^2) correlator for the Clay residue f1^2.

O_s(t) = sum_x (1/3) sum_{i<j spatial} (1/2) Re tr U_ij(x,t),  APE-smeared
spatial links at level s in LEVELS.  Stores per-config time-slice sums
S{n}[k,t] (shape Ncfg x L) so the connected correlator and spectral ratios
can be rebuilt offline.  Reuses lattice-cpb/w4_blocked.py (SU(2) heat bath,
su2_heatbath) and su3.project_su3 (polar projection works for 2x2 too).
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
import w4_blocked as W   # noqa: E402
import su3 as S          # noqa: E402

adj = W.adj
LEVELS = (0, 2, 4, 8, 12)


def mag_density_su2(U):
    """s(x) = (1/3) sum_{1<=i<j<=3} (1/2) Re tr U_ij(x)  (spatial plaquettes)."""
    L = U.shape[1]
    s = np.zeros((L, L, L, L))
    for i in (1, 2, 3):
        for j in range(i + 1, 4):
            Up = (U[i] @ np.roll(U[j], -1, axis=i)
                  @ adj(np.roll(U[i], -1, axis=j)) @ adj(U[j]))
            s += 0.5 * np.real(np.trace(Up, axis1=-2, axis2=-1))
    return s / 3.0


def ape_spatial2(U, alpha=0.5, n=1):
    """APE-smear the three spatial links (mu=1,2,3); temporal fixed."""
    U = U.copy()
    for _ in range(n):
        old = U.copy()
        for i in (1, 2, 3):
            Ss = np.zeros_like(U[i])
            for j in (1, 2, 3):
                if j == i:
                    continue
                Ss += old[j] @ np.roll(old[i], -1, axis=j) @ adj(np.roll(old[j], -1, axis=i))
                Ss += (adj(np.roll(old[j], 1, axis=j)) @ np.roll(old[i], 1, axis=j)
                       @ np.roll(np.roll(old[j], 1, axis=j), -1, axis=i))
            U[i] = S.project_su3((1.0 - alpha) * old[i] + alpha / 6.0 * Ss)
    return U


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=6)
    ap.add_argument('--beta', type=float, required=True)
    ap.add_argument('--N', type=int, default=2000)
    ap.add_argument('--therm', type=int, default=800)
    ap.add_argument('--gap', type=int, default=2)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--alpha', type=float, default=0.5)
    ap.add_argument('--levels', type=str, default=None,
                    help='comma list, default 0,2,4,8,12')
    ap.add_argument('--out', type=str, required=True)
    args = ap.parse_args()
    global LEVELS
    if args.levels is not None:
        LEVELS = tuple(int(x) for x in args.levels.split(','))

    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    cfgs = W.sample2(args.L, args.beta, rng, args.therm, args.N, args.gap)
    t1 = time.time()
    print(f'[L={args.L} b={args.beta}] sampled {len(cfgs)} cfgs, therm {args.therm} '
          f'gap {args.gap} in {t1 - t0:.0f}s', flush=True)

    out = {n: [] for n in LEVELS}
    plaq = []
    t1 = time.time()
    for k, U in enumerate(cfgs):
        plaq.append(float(np.mean(mag_density_su2(U))))
        prev = 0
        Us = U
        for n in LEVELS:
            if n > prev:
                Us = ape_spatial2(Us, args.alpha, n - prev)
                prev = n
            out[n].append(mag_density_su2(Us).sum(axis=(1, 2, 3)))  # zero-momentum per t
        if (k + 1) % 250 == 0:
            print(f'  [{k + 1}/{args.N}] {time.time() - t1:.0f}s', flush=True)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    np.savez(args.out, L=args.L, beta=args.beta, alpha=args.alpha,
             levels=np.array(LEVELS), plaq=np.array(plaq), N=args.N, gap=args.gap,
             therm=args.therm, seed=args.seed,
             **{f'S{n}': np.array(v) for n, v in out.items()})
    print('saved', args.out, f'total {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
