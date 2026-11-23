"""beta=6.0 control (weak-coupling end of the same blocking) + wstaple variant.

usage: python run_control.py [beta6|wstaple]
"""
import json
import os
import sys
import time

import numpy as np

import imc
import validate

LOG = open('runlog.txt', 'a', buffering=1)


def log(m):
    LOG.write(m + '\n')


def fitval(path, ops, iters, sweeps, val_sweeps, tag, width=1.0):
    d = dict(np.load(path))
    t0 = time.time()
    log(f"=== fit {tag} ops={ops} width={width}")
    fit, _ = imc.fit(d, ops, iters=iters, sweeps=sweeps, meas=2, seed=1,
                     verbose=True, width=width)
    json.dump(fit, open(f'fits/{tag}.json', 'w'), indent=1)
    for h in fit['hist']:
        log(f"   it{h['it']} lam={np.round(h['lam'],4)} r={np.round(h['residual'],3)}"
            f" cond={h['cond_mod']:.2f}")
    log(f"   lam={np.round(fit['lam'],4)} cond_true={fit['cond_true']:.2f} "
        f"({time.time()-t0:.0f}s)")
    if val_sweeps > 0:
        pred = validate.sample_predictions(d, np.array(fit['lam']), fit['ops'],
                                           sweeps=val_sweeps)
        rows = validate.compare(d, fit['lam'], fit['ops'], pred)
        json.dump(dict(lam=fit['lam'], kinds=fit['ops'], rows=rows, pred=pred),
                  open(f'fits/val_{tag}.json', 'w'), indent=1)
        for r in rows:
            log(f"   {r['obs']:8s} true={r['true']:+.6f}({r['err_true']:.6f}) "
                f"model={r['model']:+.6f}({r['err_model']:.6f}) "
                f"rel={100*r['rel']:+.1f}% z={r['z']:+.1f}")
    return fit


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'beta6'
    log(f"\n##### run_control {time.ctime()}")
    if which == 'beta6':
        fitval('data/true_L8_b6.0.npz', ['P', 'R', 'S'], 6, 700, 2500,
               'L8_b6.0_PRS', width=0.2)
    elif which == 'wstaple':
        d = dict(np.load('data/true_L6_b2.6_w05.npz'))
        for k in ('op_P', 'op_R', 'op_S', 'op_R3', 'c6', 'c8', 'c10', 'R4'):
            v = d[k]
            b = max(1, int(np.ceil(2 * imc.tau_int(v))))
            m, e = imc.jackknife(lambda i: float(v[i].mean()), len(v), b)
            log(f"   w05 {k:8s} {m:+.6f} +- {e:.6f}")
        fitval('data/true_L6_b2.6_w05.npz', ['P', 'R', 'S'], 6, 1200, 0,
               'L6_b2.6_w05_PRS')
    log('##### control done')
    LOG.close()
    print('done')
