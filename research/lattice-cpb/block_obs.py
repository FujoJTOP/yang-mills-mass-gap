"""Block-scale measurements for the Bianchi/z(q) questions.

Per configuration:
  * character coefficients of the microscopic plaquette distribution
        c3 = <Re tr U_p>,  c6 = <((tr U)^2+tr U^2)/2>,
        c8 = <|tr U|^2 - 1>,  c3bar = <Re((tr U)^2-tr U^2)/2>  (= c3 for SU(3))
        c10 = <Re[(tr U)^3 + 3 tr U tr U^2 + 2 tr U^3]/6>  (dim-10)
    -> single-plaquette signed/absolute proxy
        r1 = |1 + sum_R c_R d_R| / (1 + sum_R |c_R| d_R)
  * defect clusters of the blocked plaquette field P_b(y) (b=1,2):
    mask = P_b < theta, 6-connectivity (non-periodic), size histogram,
    largest cluster, number of clusters, bounding-box extents.

This does NOT compute kappa_c: the polymer weights w(C) of the blocked
effective action are not accessible from a direct MC.  r1 is a single-plaquette
proxy only.
"""
import argparse
import os

import numpy as np
from scipy import ndimage

import su3 as S

HERE = os.path.dirname(os.path.abspath(__file__))


def block_mean(P, b):
    L = P.shape[0]
    n = L // b
    return P.reshape(n, b, n, b, n, b, n, b).mean(axis=(1, 3, 5, 7))


def struct6():
    st = np.zeros((3, 3, 3, 3), bool)
    st[1, 1, 1, 1] = True
    for d in range(4):
        for s in (-1, 1):
            idx = [1, 1, 1, 1]
            idx[d] = 1 + s
            st[tuple(idx)] = True
    return st


ST6 = struct6()


def cluster_stats(mask, b):
    lab, n = ndimage.label(mask, structure=ST6)
    if n == 0:
        return np.zeros(101), 0, 0
    sizes = np.bincount(lab.ravel())[1:]
    hist = np.bincount(np.clip(sizes, 0, 100), minlength=101)[:101]
    return hist, int(sizes.max()), int(n)


def chars(U):
    tr = S.plaquette_traces(U)  # Re tr/3
    out = {}
    vals3, vals6, vals8, vals3b, vals10 = [], [], [], [], []
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Up = (U[mu] @ np.roll(U[nu], -1, axis=mu)
                  @ S.adj(np.roll(U[mu], -1, axis=nu)) @ S.adj(U[nu]))
            t = np.trace(Up, axis1=-2, axis2=-1)
            t2 = np.trace(Up @ Up, axis1=-2, axis2=-1)
            vals3.append(np.real(t).mean())
            vals6.append(np.real(0.5 * (t * t + t2)).mean())
            vals8.append(np.real(np.abs(t) ** 2 - 1.0).mean())
            vals3b.append(np.real(0.5 * (t * t - t2)).mean())
            t3 = np.trace(Up @ Up @ Up, axis1=-2, axis2=-1)
            vals10.append(np.real((t * t * t + 3.0 * t * t2 + 2.0 * t3) / 6.0).mean())
    out['c3'] = float(np.mean(vals3))
    out['c6'] = float(np.mean(vals6))
    out['c8'] = float(np.mean(vals8))
    out['c3bar'] = float(np.mean(vals3b))
    out['c10'] = float(np.mean(vals10))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=6)
    ap.add_argument('--beta', type=float, required=True)
    ap.add_argument('--sweeps', type=int, default=3000)
    ap.add_argument('--therm', type=int, default=1500)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', type=str, required=True)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    lat = S.Lattice(args.L, args.beta, rng, cold=True)
    for _ in range(args.therm):
        lat.sweep()

    recs = []
    for k in range(args.sweeps):
        lat.sweep()
        U = lat.U
        tr = S.plaquette_traces(U)
        P = np.mean(np.stack([a for a in tr.values()]), axis=0)
        rec = {'plaq': S.mean_plaquette(U), **chars(U)}
        for b in (1, 2):
            Pb = block_mean(P, b)
            mp = float(Pb.mean())
            sp = float(Pb.std())
            rec[f'mean_{b}'] = mp
            rec[f'var_{b}'] = sp * sp
            rec[f'skew_{b}'] = float(np.mean(((Pb - mp) / max(sp, 1e-300)) ** 3))
            for tag, th in (('0', 0.0), ('hh', 0.5 * mp), ('h75', 0.75 * mp)):
                mask = Pb < th
                hist, big, ncl = cluster_stats(mask, b)
                rec[f'hist_{b}_{tag}'] = hist
                rec[f'big_{b}_{tag}'] = big
                rec[f'ncl_{b}_{tag}'] = ncl
                rec[f'frac_{b}_{tag}'] = float(np.mean(mask))
                rec[f'k_{b}_{tag}'] = int(mask.sum())
        recs.append(rec)
        if (k + 1) % max(1, args.sweeps // 5) == 0:
            print(f"  [{k+1}/{args.sweeps}] P={rec['plaq']:.5f}", flush=True)

    out = {'L': args.L, 'beta': args.beta, 'seed': args.seed,
           'sweeps': args.sweeps, 'therm': args.therm}
    for key in recs[0]:
        out[key] = np.array([r[key] for r in recs])
    np.savez(args.out, **out)
    print('saved', args.out)


if __name__ == '__main__':
    main()
