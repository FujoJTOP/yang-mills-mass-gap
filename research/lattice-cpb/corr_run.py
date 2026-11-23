"""Dedicated correlation-length run: connected magnetic-plaquette correlator
with APE-smeared operators (levels n=0 and n=8).

usage: python corr_run.py --L 6 --beta 1.7 --sweeps 4000 --therm 1500 --seed 12 --out corr/L6_b1.7.npz
"""
import argparse
import os

import numpy as np

import su3 as S

LEVELS = (0, 8)


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
    lat = S.Lattice(args.L, args.beta, rng, cold=True)
    for _ in range(args.therm):
        lat.sweep()

    out = {k: [] for k in LEVELS}
    plaq = []
    for k in range(args.sweeps):
        lat.sweep()
        Us = lat.U
        prev = 0
        for n in LEVELS:
            if n > prev:
                Us = S.ape_smear_spatial(Us, args.alpha, n - prev)
                prev = n
            out[n].append(S.magnetic_density(Us).sum(axis=(1, 2, 3)))
        plaq.append(S.mean_plaquette(lat.U))
        if (k + 1) % max(1, args.sweeps // 10) == 0:
            print(f"  [{k+1}/{args.sweeps}] P={plaq[-1]:.5f}", flush=True)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    np.savez(args.out, L=args.L, beta=args.beta, alpha=args.alpha,
             levels=np.array(LEVELS), plaq=np.array(plaq),
             **{f'S{n}': np.array(v) for n, v in out.items()})
    print("saved", args.out)


if __name__ == '__main__':
    main()
