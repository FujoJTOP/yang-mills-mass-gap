"""Out-of-sample validation of a fitted blocked action.

The fit matches the per-loop densities of the loop families in `ops`.
Everything else measured on the blocked ensemble is a prediction:
  c6, c8, c10   (blocked single-plaquette character coefficients)
  c3 = 3 * op_P (identical to the fit operator, by definition)
  R4            (4th cumulant ratio of the site-averaged blocked plaquette)
  CP_<r>        (connected plaquette-plaquette correlators)
  CL_r          (link-link correlator, straight 2-link line)
"""
import json
import os
import sys

import numpy as np

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              '..', 'lattice-cpb'))

import model as M
import imc
import blk


def sample_predictions(d, lam, kinds, sweeps=4000, therm=400, meas=2, seed=7):
    L = int(d['L']) // 2          # model lives on the blocked lattice
    recs, a, _ = imc.model_run(lam, kinds, L, sweeps, therm, meas, seed,
                               verbose=True, full=True)
    keys = [k for k in recs[0] if k != 'micro_P' and k != 'hist']
    n = len(recs)
    out = {}
    for k in keys:
        v = np.array([r[k] for r in recs])
        b = max(1, int(np.ceil(2 * imc.tau_int(v))) if v.ndim == 1 else 1)
        m, e = imc.jackknife(lambda idx: float(v[idx].mean()), n, b)
        out[k] = [float(m), float(e)]
    out['_acc'] = a
    out['_nkinds'] = kinds
    return out


def compare(d, lam, kinds, pred):
    """Table of true vs predicted for the held-out observables."""
    leave = {'op_' + k for k in kinds} | {'micro_P'}
    from imc import DATA_KEY
    leave |= {DATA_KEY[k] for k in kinds}
    rows = []
    for k, val in pred.items():
        if k.startswith('_') or k not in d:
            continue
        m, e = val
        vt = np.asarray(d[k])
        b = max(1, int(np.ceil(2 * imc.tau_int(vt))) if vt.ndim == 1 else 1)
        mt, et = imc.jackknife(lambda idx: float(vt[idx].mean()), len(vt), b)
        diff = m - mt
        err = np.sqrt(e ** 2 + et ** 2)
        rows.append(dict(obs=k, true=float(mt), err_true=float(et),
                         model=float(m), err_model=float(e),
                         diff=float(diff), z=float(diff / err) if err > 0 else np.inf,
                         rel=float(diff / abs(mt)) if mt != 0 else np.inf,
                         fit=(k in leave)))
    return rows


def main():
    ap = __import__('argparse').ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--fit', required=True, help='json from imc.py')
    ap.add_argument('--sweeps', type=int, default=4000)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    d = dict(np.load(args.data))
    f = json.load(open(args.fit))
    lam = np.array(f['lam'])
    kinds = f['ops']
    pred = sample_predictions(d, lam, kinds, sweeps=args.sweeps)
    rows = compare(d, lam, kinds, pred)
    json.dump(dict(lam=lam.tolist(), kinds=kinds, rows=rows, pred=pred),
              open(args.out, 'w'), indent=1)
    for r in rows:
        print(f"{r['obs']:8s} true={r['true']:+.5f}({r['err_true']:.5f}) "
              f"model={r['model']:+.5f}({r['err_model']:.5f}) "
              f"rel={100*r['rel']:+.1f}%  z={r['z']:+.1f}")


if __name__ == '__main__':
    main()
