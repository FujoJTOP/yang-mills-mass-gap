"""NA-3-T run: Dobrushin/SDPI coefficients of the SU(3) link heat bath.

usage: python dobi_run.py --L 4 --beta 3.0 --ncfg 40 --nlink 3 --nsamp 80 --out eta_b3.0.npz
"""
import argparse
import os

import numpy as np

import su3 as S
import dobi as D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=4)
    ap.add_argument('--beta', type=float, default=3.0)
    ap.add_argument('--ncfg', type=int, default=40)
    ap.add_argument('--nlink', type=int, default=3)
    ap.add_argument('--nsamp', type=int, default=80)
    ap.add_argument('--nsweep', type=int, default=20)
    ap.add_argument('--therm', type=int, default=600)
    ap.add_argument('--sep', type=int, default=5)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', type=str, required=True)
    args = ap.parse_args()

    L = args.L
    kappa = args.beta / 3.0
    rng = np.random.default_rng(args.seed)
    lat = S.Lattice(L, args.beta, rng, cold=True)
    for _ in range(args.therm):
        lat.sweep()

    e = np.eye(4, dtype=int)
    rows, ub_rows, tv_pairs, lb_rows = [], [], [], []
    best_lb = (-1.0, None, None, None, None)
    for ic in range(args.ncfg):
        for _ in range(args.sep):
            lat.sweep()
        U = lat.U
        for _ in range(args.nlink):
            x = tuple(rng.integers(0, L, 4))
            mu = int(rng.integers(0, 4))
            Sp = D.staple_at(U, x, mu)
            row_tv = 0.0
            row_ub = 0.0
            row_lb = 0.0
            wit = []
            for (xp, mup) in D.neighbor_links(x, mu, L):
                g = S.random_su3(rng)
                U2 = U.copy()
                U2[mup][xp] = g @ U[mup][xp]
                Sq = D.staple_at(U2, x, mu)
                dS = Sq - Sp
                ub = (kappa / 2.0) * D.nuclear_norm(dS)
                tv = D.tv_bar(Sp, Sq, kappa, nsamp=args.nsamp,
                              nsweep=args.nsweep, rng=rng)
                lb, val, se = D.tv_witness_lb(Sp, Sq, kappa, nsamp=args.nsamp,
                                              nsweep=args.nsweep, rng=rng)
                row_lb += lb
                wit.append((xp, mup, g, lb, val, se))
                row_tv += tv
                row_ub += min(ub, 1.0)
                tv_pairs.append((np.linalg.norm(dS), tv))
            rows.append(row_tv)
            ub_rows.append(row_ub)
            lb_rows.append(row_lb)
            if row_lb > best_lb[0]:
                best_lb = (row_lb, U.copy(), x, mu, wit)
        if (ic + 1) % max(1, args.ncfg // 5) == 0:
            print(f"  cfg {ic+1}/{args.ncfg}: row-sum eta mean {np.mean(rows):.3f} "
                  f"max {np.max(rows):.3f}", flush=True)

    rows = np.array(rows); ub_rows = np.array(ub_rows); lb_rows = np.array(lb_rows)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    wfit = {'L': L, 'beta': args.beta, 'row_lb': best_lb[0],
            'x': np.array(best_lb[2]), 'mu': best_lb[3], 'U': best_lb[1]}
    if best_lb[4] is not None:
        wfit['link_coords'] = np.array([list(w[0]) + [w[1]] for w in best_lb[4]])
        wfit['link_lb'] = np.array([[w[3], w[4], w[5]] for w in best_lb[4]])
        wfit['gs'] = np.array([w[2] for w in best_lb[4]])
    np.savez(os.path.join(os.path.dirname(os.path.abspath(args.out)),
                          f'witness_L{L}_b{args.beta:.1f}.npz'), **wfit)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    np.savez(args.out, L=L, beta=args.beta, kappa=kappa,
             eta_rows=rows, ub_rows=ub_rows, lb_rows=lb_rows,
             tv_pairs=np.array(tv_pairs),
             ncfg=args.ncfg, nlink=args.nlink, nsamp=args.nsamp)
    print(f"beta={args.beta}: eta_row mean {rows.mean():.3f} max {rows.max():.3f} "
          f"| rigorous ub mean {ub_rows.mean():.2f} max {ub_rows.max():.2f}")
    print("saved", args.out)


if __name__ == '__main__':
    main()
