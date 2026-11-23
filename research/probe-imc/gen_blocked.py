"""Fine-lattice SU(3) MC -> b=2 link blocking -> blocked-observable ensemble.

usage:
  python gen_blocked.py --L 6 --beta 2.6 --sweeps 3000 --therm 800 --meas 1 \
                        --seed 1 --out data/true_L6_b2.6.npz
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              '..', 'lattice-cpb'))

import su3 as S
import su3_numba as N
import blk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=6)
    ap.add_argument('--beta', type=float, required=True)
    ap.add_argument('--sweeps', type=int, default=3000)
    ap.add_argument('--therm', type=int, default=800)
    ap.add_argument('--meas', type=int, default=1)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--wstaple', type=float, default=1.0)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    L = args.L
    np.random.seed(args.seed)
    U = np.zeros((4, L, L, L, L, 3, 3), dtype=np.complex128)
    for mu in range(4):
        U[mu][...] = np.eye(3)
    kappa = args.beta / 3.0
    t0 = time.time()
    for k in range(args.therm):
        N.sweep_par(U, L, kappa)

    n = L // 2
    orbits = blk.orbits_for(n)
    recs = []
    t0 = time.time()
    for k in range(args.sweeps):
        N.sweep_par(U, L, kappa)
        if (k + 1) % args.meas:
            continue
        V = blk.block_links(U, args.wstaple)
        rec = blk.measure_all(V, orbits)
        rec['micro_P'] = S.mean_plaquette(U)
        recs.append(rec)
        if (len(recs) % max(1, args.sweeps // 10)) == 0:
            print(f"  [{len(recs)}] micro_P={rec['micro_P']:.5f} "
                  f"P'={rec['op_P']:.5f} {time.time()-t0:.0f}s", flush=True)

    out = {'L': L, 'beta': args.beta, 'seed': args.seed, 'nmeas': len(recs),
           'sweeps': args.sweeps, 'therm': args.therm, 'meas': args.meas,
           'wstaple': args.wstaple,
           'orbits': np.array(orbits)}
    for key in recs[0]:
        out[key] = np.array([r[key] for r in recs])
    np.savez(args.out, **out)
    print('saved', args.out, 'nmeas', len(recs), 'time', round(time.time() - t0))


if __name__ == '__main__':
    main()
