"""Pilot run at large L: parallel Numba heat bath + APE-smeared magnetic
correlator.  usage: python pilot_numba.py --L 24 --beta 6.30 --sweeps 3000 ..."""
import argparse
import os
import time

import numpy as np

import su3_numba as N


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=24)
    ap.add_argument('--beta', type=float, default=6.30)
    ap.add_argument('--sweeps', type=int, default=3000)
    ap.add_argument('--therm', type=int, default=400)
    ap.add_argument('--meas', type=int, default=10)
    ap.add_argument('--nlev', type=int, default=8)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', type=str, required=True)
    args = ap.parse_args()

    L = args.L
    U = np.zeros((4, L, L, L, L, 3, 3), dtype=np.complex128)
    for mu in range(4):
        U[mu][...] = np.eye(3)
    kappa = args.beta / 3.0
    np.random.seed(args.seed)

    t0 = time.time()
    for k in range(args.therm):
        N.sweep_par(U, L, kappa)
        if (k + 1) % max(1, args.therm // 4) == 0:
            print(f"  therm {k+1}/{args.therm}  P={N.plaquette_mean(U,L):.5f}", flush=True)

    S0, Sn, plaq, taus = [], [], [], []
    for k in range(args.sweeps):
        N.sweep_par(U, L, kappa)
        plaq.append(N.plaquette_mean(U, L))
        if (k + 1) % args.meas == 0:
            s0 = N.mag_slice(U, L)
            if args.nlev > 0:
                V = N.ape_smear_num_par(U, L, 0.5, args.nlev)
                sn = N.mag_slice(V, L)
            else:
                sn = s0
            S0.append(s0)
            Sn.append(sn)
        if (k + 1) % max(1, args.sweeps // 8) == 0:
            print(f"  [{k+1}/{args.sweeps}] {time.time()-t0:.0f}s  P={plaq[-1]:.5f}", flush=True)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    np.savez(args.out, L=L, beta=args.beta, seed=args.seed, sweeps=args.sweeps,
             therm=args.therm, meas=args.meas, nlev=args.nlev,
             plaq=np.array(plaq), S0=np.array(S0), Sn=np.array(Sn))
    print("saved", args.out, f"({time.time()-t0:.0f}s)")


if __name__ == '__main__':
    main()
