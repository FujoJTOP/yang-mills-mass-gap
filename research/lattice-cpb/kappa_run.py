"""Compute kappa_eff (projected Hessian) and per-link curvature statistics on
saved configurations.  usage: python kappa_run.py [L ...] [--nsnap N] [--m M]"""
import argparse
import glob
import os

import numpy as np

import su3 as S
import kappa as K

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nsnap', type=int, default=8)
    ap.add_argument('--m', type=int, default=60)
    ap.add_argument('--Lam', type=float, default=1e3)
    ap.add_argument('--Ls', type=int, nargs='*', default=[4, 6])
    ap.add_argument('--tag', type=str, default='all')
    ap.add_argument('--betas', type=float, nargs='*', default=None)
    args = ap.parse_args()

    recs = []
    for f in sorted(glob.glob(os.path.join(HERE, 'runs', '*.npz'))):
        d = dict(np.load(f, allow_pickle=True))
        if 'snapshots' not in d:
            continue
        L = int(d['L']); beta = float(d['beta'])
        if L not in args.Ls:
            continue
        if args.betas is not None and not any(abs(beta - b) < 1e-6 for b in args.betas):
            continue
        snaps = d['snapshots']
        idx = np.unique(np.linspace(0, len(snaps) - 1, min(args.nsnap, len(snaps))).astype(int))
        print(f"== {os.path.basename(f)}  L={L} beta={beta}  {len(idx)} snapshots", flush=True)
        for j in idx:
            U = snaps[j]
            hl = S.link_hessian_min(U, beta)
            res = K.lowest_eigenvalue(U, beta, m=args.m, seed=int(j), Lam=args.Lam)
            rec = dict(L=L, beta=beta, snap=int(j), kappa=res['kappa'],
                       lanczos_res=res['lanczos_res'], cg_res=res['cg_res_max'],
                       ritz=res['ritz'], hl_min=float(hl.min()),
                       hl_neg_frac=float(np.mean(hl < 0.0)),
                       hl_mean=float(hl.mean()))
            recs.append(rec)
            print(f"   snap {j:3d}: kappa={res['kappa']:+.5f} res={res['lanczos_res']:.1e} "
                  f"cg={res['cg_res_max']:.1e}  linkH min={hl.min():+.4f} "
                  f"neg={rec['hl_neg_frac']:.3f}", flush=True)

    out = {}
    for key in ('L', 'beta', 'snap', 'kappa', 'lanczos_res', 'cg_res',
                'hl_min', 'hl_neg_frac', 'hl_mean'):
        out[key] = np.array([r[key] for r in recs])
    out['ritz'] = np.array([r['ritz'] for r in recs])
    fn = os.path.join(HERE, f'kappa_results_{args.tag}.npz')
    np.savez(fn, **out)
    print("saved", fn)


if __name__ == '__main__':
    main()
