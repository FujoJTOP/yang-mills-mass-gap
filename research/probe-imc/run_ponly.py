"""One-coupling (Wilson-only) fit + validation at both beta.  Tests how much of
the blocked theory a single plaquette coupling already reproduces."""
import json
import time

import numpy as np

import imc
import validate

LOG = open('runlog.txt', 'a', buffering=1)


def log(m):
    LOG.write(m + '\n')


def case(path, tag):
    d = dict(np.load(path))
    t0 = time.time()
    log(f"\n=== fit {tag} ops=[P]")
    fit, _ = imc.fit(d, ['P'], iters=6, sweeps=1500, meas=2, seed=1, verbose=True)
    json.dump(fit, open(f'fits/{tag}.json', 'w'), indent=1)
    for h in fit['hist']:
        log(f"   it{h['it']} lam={np.round(h['lam'],4)} r={np.round(h['residual'],4)}"
            f" cond={h['cond_mod']:.2f}")
    log(f"   lam={np.round(fit['lam'],4)} ({time.time()-t0:.0f}s)")
    pred = validate.sample_predictions(d, np.array(fit['lam']), fit['ops'],
                                       sweeps=3000)
    rows = validate.compare(d, fit['lam'], fit['ops'], pred)
    json.dump(dict(lam=fit['lam'], kinds=fit['ops'], rows=rows, pred=pred),
              open(f'fits/val_{tag}.json', 'w'), indent=1)
    for r in rows:
        log(f"   {r['obs']:9s} true={r['true']:+.6f}({r['err_true']:.6f}) "
            f"model={r['model']:+.6f}({r['err_model']:.6f}) "
            f"rel={100*r['rel']:+7.1f}% z={r['z']:+6.1f}"
            + ("  FIT" if r['fit'] else ""))


if __name__ == '__main__':
    log(f"\n##### run_ponly {time.ctime()}")
    for b in ('2.6', '3.0'):
        case(f'data/true_L8_b{b}.npz', f'L8_b{b}_P')
    log('##### ponly done')
    LOG.close()
    print('done')
