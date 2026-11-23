"""Driver: fit + validate for every dataset / operator set.

usage:
  python run_all.py --datasets data/true_L6_b2.6.npz data/true_L6_b3.0.npz \
      --ops P,R,S --iters 8 --sweeps 1500 --outdir fits
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

import imc
import validate


def tag(path, ops):
    b = os.path.basename(path).replace('true_', '').replace('.npz', '')
    return b + '_' + ''.join(ops)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--datasets', nargs='+', required=True)
    ap.add_argument('--ops', default='P,R,S')
    ap.add_argument('--iters', type=int, default=8)
    ap.add_argument('--sweeps', type=int, default=1500)
    ap.add_argument('--meas', type=int, default=2)
    ap.add_argument('--val-sweeps', type=int, default=3000)
    ap.add_argument('--outdir', default='fits')
    ap.add_argument('--no-validate', action='store_true')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    ops = args.ops.split(',')
    for path in args.datasets:
        d = dict(np.load(path))
        t = tag(path, ops)
        print(f"=== {t}", flush=True)
        fit, _ = imc.fit(d, ops, iters=args.iters, sweeps=args.sweeps,
                         meas=args.meas, seed=1)
        json.dump(fit, open(os.path.join(args.outdir, t + '.json'), 'w'), indent=1)
        print('  lam =', np.round(fit['lam'], 4), ' cond_true =',
              round(fit['cond_true'], 2), flush=True)
        if args.no_validate:
            continue
        pred = validate.sample_predictions(d, np.array(fit['lam']), fit['ops'],
                                           sweeps=args.val_sweeps)
        rows = validate.compare(d, fit['lam'], fit['ops'], pred)
        json.dump(dict(lam=fit['lam'], kinds=fit['ops'], rows=rows, pred=pred),
                  open(os.path.join(args.outdir, 'val_' + t + '.json'), 'w'),
                  indent=1)
        for r in rows:
            print(f"   {r['obs']:8s} true={r['true']:+.5f}({r['err_true']:.5f}) "
                  f"model={r['model']:+.5f}({r['err_model']:.5f}) "
                  f"rel={100*r['rel']:+.1f}% z={r['z']:+.1f}", flush=True)


if __name__ == '__main__':
    main()
