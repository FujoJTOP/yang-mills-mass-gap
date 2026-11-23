"""Pilot: APE-smeared connected magnetic-plaquette correlator on the Numba kernel.

Same operator as lattice-cpb/corr_run.py (levels n=0 and n=8, alpha=0.5),
but driven by su3_numba (sweep_par / ape_smear_num_par / mag_slice).

usage: python pilot_run.py --L 10 --beta 6.2 --sweeps 30000 --therm 2000 --seed 3 --out pilot/L10_b6.2.npz
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lattice-cpb'))
import su3_numba as N


def cold(L):
    U = np.zeros((4, L, L, L, L, 3, 3), dtype=np.complex128)
    for mu in range(4):
        for i in range(3):
            U[mu][..., i, i] = 1.0
    return U


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=10)
    ap.add_argument('--beta', type=float, required=True)
    ap.add_argument('--sweeps', type=int, default=30000)
    ap.add_argument('--therm', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--alpha', type=float, default=0.5)
    ap.add_argument('--nsmear', type=int, default=8,
                    help='largest smearing level (kept for back-compat)')
    ap.add_argument('--levels', type=str, default='',
                    help='comma-separated APE levels, e.g. 0,8,16,32 (overrides --nsmear)')
    ap.add_argument('--measure-every', type=int, default=1)
    ap.add_argument('--save-every', type=int, default=2000)
    ap.add_argument('--out', type=str, required=True)
    args = ap.parse_args()

    if args.levels:
        levels = tuple(sorted(int(x) for x in args.levels.split(',')))
    else:
        levels = (0, args.nsmear)

    L, kappa = args.L, args.beta / 3.0
    rng = np.random.default_rng(args.seed)
    np.random.seed(args.seed)

    U = cold(L)
    t0 = time.time()
    for _ in range(args.therm):
        N.sweep_par(U, L, kappa)
    print(f"thermalised in {time.time()-t0:.0f}s", flush=True)

    data = {n: [] for n in levels}
    plaq = []
    t0 = time.time()
    for k in range(args.sweeps):
        N.sweep_par(U, L, kappa)
        if (k + 1) % args.measure_every != 0:
            continue
        Us, prev = U, 0
        for n in levels:
            if n > prev:
                Us = N.ape_smear_num_par(Us, L, args.alpha, n - prev)
                prev = n
            data[n].append(N.mag_slice(Us, L))
        plaq.append(N.plaquette_mean(U, L))
        if (k + 1) % max(1, args.sweeps // 20) == 0:
            dt = time.time() - t0
            print(f"  [{k+1}/{args.sweeps}] P={plaq[-1]:.5f} "
                  f"({dt:.0f}s, {dt/(k+1):.3f}s/step)", flush=True)
        if (k + 1) % args.save_every == 0:
            save(args, plaq, data, levels)

    save(args, plaq, data, levels)
    print("saved", args.out, flush=True)


def save(args, plaq, data, levels):
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    np.savez(args.out, L=args.L, beta=args.beta, alpha=args.alpha,
             levels=np.array(levels), plaq=np.array(plaq),
             **{f'S{n}': np.array(v) for n, v in data.items()})


if __name__ == '__main__':
    main()
