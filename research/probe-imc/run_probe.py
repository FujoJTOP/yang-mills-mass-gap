"""Full probe driver: fits and validations, sequential, file-logged."""
import json
import os
import sys
import time

import numpy as np

import imc
import validate

LOG = open('runlog.txt', 'a', buffering=1)


def log(msg):
    LOG.write(msg + '\n')


def do(dpath, ops, iters, sweeps, meas, val_sweeps, tag):
    d = dict(np.load(dpath))
    t0 = time.time()
    log(f"=== fit {tag} ops={ops} iters={iters} sweeps={sweeps}")
    fit, _ = imc.fit(d, ops, iters=iters, sweeps=sweeps, meas=meas, seed=1,
                     verbose=True)
    json.dump(fit, open(f'fits/{tag}.json', 'w'), indent=1)
    for h in fit['hist']:
        log(f"   it{h['it']} lam={np.round(h['lam'],4)} "
            f"r={np.round(h['residual'],3)} cond={h['cond_mod']:.1f}")
    log(f"   lam={np.round(fit['lam'],4)} cond_true={fit['cond_true']:.2f} "
        f"eig={np.round(fit['eig_true'],4)} ({time.time()-t0:.0f}s)")
    if val_sweeps <= 0:
        return fit
    t0 = time.time()
    pred = validate.sample_predictions(d, np.array(fit['lam']), fit['ops'],
                                       sweeps=val_sweeps)
    rows = validate.compare(d, fit['lam'], fit['ops'], pred)
    json.dump(dict(lam=fit['lam'], kinds=fit['ops'], rows=rows, pred=pred),
              open(f'fits/val_{tag}.json', 'w'), indent=1)
    for r in rows:
        log(f"   {r['obs']:8s} true={r['true']:+.6f}({r['err_true']:.6f}) "
            f"model={r['model']:+.6f}({r['err_model']:.6f}) "
            f"rel={100*r['rel']:+.1f}% z={r['z']:+.1f}")
    log(f"   validate ({time.time()-t0:.0f}s)")
    return fit


if __name__ == '__main__':
    os.makedirs('fits', exist_ok=True)
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    log(f"\n##### run_probe {time.ctime()}")
    if which in ('all', 'L8'):
        for b in ('2.6', '3.0'):
            f = f'data/true_L8_b{b}.npz'
            do(f, ['P', 'R', 'S'], 7, 900, 2, 3000, f'L8_b{b}_PRS')
            do(f, ['P', 'R'], 6, 900, 2, 3000, f'L8_b{b}_PR')
    if which in ('all', 'L6'):
        for b in ('2.6', '3.0'):
            f = f'data/true_L6_b{b}.npz'
            do(f, ['P', 'R', 'S'], 7, 1200, 2, 4000, f'L6_b{b}_PRS')
            do(f, ['P', 'R'], 6, 1200, 2, 4000, f'L6_b{b}_PR')
    log('##### done')
    LOG.close()
    print('done')
