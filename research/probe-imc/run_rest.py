"""Resume driver: for every pending case, reuse an existing fit JSON if present,
otherwise fit; then validate.  Robust: a failure in one case does not stop the rest.

usage: python run_rest.py [group]
  group = all | loop | char | l6 | control
"""
import json
import os
import sys
import time
import traceback

import numpy as np

import imc
import validate

LOG = open('runlog.txt', 'a', buffering=1)


def log(m):
    LOG.write(m + '\n')


def case(path, ops, tag, iters=7, sweeps=900, val_sweeps=3000, width=1.0):
    try:
        d = dict(np.load(path))
        fj = f'fits/{tag}.json'
        if os.path.exists(fj):
            fit = json.load(open(fj))
            log(f"=== reuse {tag} lam={np.round(fit['lam'],4)}")
        else:
            t0 = time.time()
            log(f"=== fit {tag} ops={ops}")
            fit, _ = imc.fit(d, ops, iters=iters, sweeps=sweeps, meas=2, seed=1,
                             verbose=True, width=width)
            json.dump(fit, open(fj, 'w'), indent=1)
            for h in fit['hist']:
                log(f"   it{h['it']} lam={np.round(h['lam'],4)} "
                    f"r={np.round(h['residual'],3)} cond={h['cond_mod']:.2f}")
            log(f"   lam={np.round(fit['lam'],4)} cond_true={fit['cond_true']:.2f} "
                f"({time.time()-t0:.0f}s)")
        if val_sweeps <= 0:
            return
        t0 = time.time()
        pred = validate.sample_predictions(d, np.array(fit['lam']), fit['ops'],
                                           sweeps=val_sweeps)
        rows = validate.compare(d, fit['lam'], fit['ops'], pred)
        json.dump(dict(lam=fit['lam'], kinds=fit['ops'], rows=rows, pred=pred),
                  open(f'fits/val_{tag}.json', 'w'), indent=1)
        for r in rows:
            log(f"   {r['obs']:9s} true={r['true']:+.6f}({r['err_true']:.6f}) "
                f"model={r['model']:+.6f}({r['err_model']:.6f}) "
                f"rel={100*r['rel']:+7.1f}% z={r['z']:+6.1f}"
                + ("  FIT" if r['fit'] else ""))
        log(f"   validate {time.time()-t0:.0f}s")
    except Exception:
        log(f"!!! {tag} FAILED\n" + traceback.format_exc())


if __name__ == '__main__':
    g = sys.argv[1] if len(sys.argv) > 1 else 'all'
    log(f"\n##### run_rest {g} {time.ctime()}")
    P0 = dict(iters=7, sweeps=1500, val_sweeps=3000)
    if g in ('all', 'core'):
        # priority order: loop ansatz (primary), then character basis
        case('data/true_L8_b2.6.npz', ['P', 'R', 'S'], 'L8_b2.6_PRS', **P0)
        case('data/true_L8_b3.0.npz', ['P', 'R', 'S'], 'L8_b3.0_PRS', **P0)
        case('data/true_L8_b2.6.npz', ['P', 'C6', 'C8'], 'L8_b2.6_PC6C8',
             iters=6, sweeps=2500, val_sweeps=3000)
        case('data/true_L8_b3.0.npz', ['P', 'C6', 'C8'], 'L8_b3.0_PC6C8',
             iters=6, sweeps=2500, val_sweeps=3000)
    if g in ('all', 'more'):
        case('data/true_L8_b2.6.npz', ['P', 'R'], 'L8_b2.6_PR',
             iters=6, sweeps=1500, val_sweeps=3000)
        case('data/true_L8_b2.6.npz', ['P', 'R', 'C6', 'C8'], 'L8_b2.6_PRC6C8', **P0)
        case('data/true_L8_b2.6.npz', ['P', 'C6', 'C8', 'C10'], 'L8_b2.6_PC6C8C10',
             iters=6, sweeps=1500, val_sweeps=0)
    if g in ('all', 'l6'):
        case('data/true_L6_b2.6.npz', ['P', 'R', 'S'], 'L6_b2.6_PRS',
             iters=7, sweeps=1500, val_sweeps=4000)
        case('data/true_L6_b3.0.npz', ['P', 'R', 'S'], 'L6_b3.0_PRS',
             iters=7, sweeps=1500, val_sweeps=4000)
    if g in ('all', 'control'):
        case('data/true_L8_b6.0.npz', ['P', 'R', 'S'], 'L8_b6.0_PRS',
             iters=6, sweeps=700, val_sweeps=2000, width=0.2)
    log('##### rest done')
    LOG.close()
    print('done')
