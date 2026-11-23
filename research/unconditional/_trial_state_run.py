"""Trial-state energy R(psi) for the lightest glueball, 0++ channel.

Operator: zero-momentum sum of the spatial magnetic plaquette density
    O_s(t) = sum_x Re tr F^2_{space}(x,t)/3     (a tr F^2 interpolator)
built from APE-smeared spatial links at levels s in LEVELS (incremental
smearing).  Each level is a different trial operator; the level set is the
variational basis for a GEVP.

Stores per-sweep time-slice sums S{n}[k, t] (shape nsweeps x L) so the
connected zero-momentum correlator and all its ratios can be rebuilt.

Reuses lattice-cpb/su3.py (CM heat bath, APE smearing, magnetic_density).
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
import su3 as S

LEVELS = (0, 2, 4, 8, 12)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=6)
    ap.add_argument('--beta', type=float, required=True)
    ap.add_argument('--sweeps', type=int, default=4000)
    ap.add_argument('--therm', type=int, default=1500)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--alpha', type=float, default=0.5)
    ap.add_argument('--out', type=str, required=True)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    lat = S.Lattice(args.L, args.beta, rng, cold=True)
    for _ in range(args.therm):
        lat.sweep()
    t1 = time.time()
    print(f'therm {args.therm} sweeps in {t1 - t0:.1f}s '
          f'({(t1 - t0) / max(args.therm, 1):.3f}s/sweep)', flush=True)

    out = {n: [] for n in LEVELS}
    plaq = []
    for k in range(args.sweeps):
        lat.sweep()
        Us = lat.U
        prev = 0
        for n in LEVELS:
            if n > prev:
                Us = S.ape_smear_spatial(Us, args.alpha, n - prev)
                prev = n
            d = S.magnetic_density(Us)              # (L,L,L,L) per site
            out[n].append(d.sum(axis=(1, 2, 3)))    # zero-momentum, per time
        plaq.append(S.mean_plaquette(lat.U))
        if (k + 1) % max(1, args.sweeps // 10) == 0:
            dt = time.time() - t1
            print(f'  [{k + 1}/{args.sweeps}] P={plaq[-1]:.5f} '
                  f'({dt:.0f}s)', flush=True)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    np.savez(args.out, L=args.L, beta=args.beta, alpha=args.alpha,
             levels=np.array(LEVELS), plaq=np.array(plaq),
             **{f'S{n}': np.array(v) for n, v in out.items()})
    print('saved', args.out, f'total {time.time() - t1:.0f}s', flush=True)


if __name__ == '__main__':
    main()
