"""Distribution-level check: blocked-plaquette tail, true ensemble vs model.

f(th) = fraction of blocked plaquettes with (1/3)Re tr < th; c0(th) = -ln f.
Thresholds: 0, mu/2, 0.75 mu (mu = per-configuration mean), plus fixed 0.01..0.05.
Moments: <P>, Var, skew, <P^3>.

usage: python tail_check.py --data data/true_L8_b2.6_tail.npz --fit fits/L8_b2.6_PRS.json
"""
import argparse
import json

import numpy as np

import imc


KEYS = [('f(<0)', 'n_0'), ('f(<mu/2)', 'n_hh'), ('f(<0.75mu)', 'n_h75'),
        ('f(<0.01)', 'n_01'), ('f(<0.03)', 'n_03')]


def stats(d, recs=None):
    """recs: list of dicts (model); else use arrays in d."""
    if recs is None:
        g = lambda k: np.array(d[k], float)
    else:
        g = lambda k: np.array([r[k] for r in recs], float)
    mu = g('op_P').mean()
    p2 = g('P2').mean()
    p3 = g('P3').mean()
    var = p2 - mu * mu
    out = dict(mean=mu, p2=p2, var=var, p3=p3,
               skew=(p3 - 3 * mu * p2 + 2 * mu ** 3) / var ** 1.5)
    for name, key in KEYS:
        out[name] = g(key).mean()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--fit', required=True)
    ap.add_argument('--sweeps', type=int, default=1500)
    ap.add_argument('--out', default='')
    args = ap.parse_args()

    d = dict(np.load(args.data))
    f = json.load(open(args.fit))
    lam = np.array(f['lam'])
    L = int(d['L']) // 2
    st = stats(d)
    recs, a, _ = imc.model_run(lam, f['ops'], L, args.sweeps, 400, 2, 11,
                               verbose=True, full=True)
    sm = stats(None, recs)
    print(f"{'':12s} {'true':>12s} {'model':>12s} {'rel':>9s}")
    for k in ('mean', 'p2', 'var', 'p3', 'skew'):
        print(f"{k:12s} {st[k]:12.5f} {sm[k]:12.5f} "
              f"{100*(sm[k]-st[k])/abs(st[k]) if st[k] else np.nan:+8.1f}%")
    print(f"{'threshold':12s} {'f_true':>10s} {'f_model':>10s} "
          f"{'c0_true':>8s} {'c0_model':>9s}")
    for name, _ in KEYS:
        ft, fm = st[name], sm[name]
        ct = -np.log(ft) if ft > 0 else np.inf
        cm = -np.log(fm) if fm > 0 else np.inf
        print(f"{name:12s} {ft:10.5f} {fm:10.5f} {ct:8.3f} {cm:9.3f}")
    if args.out:
        json.dump(dict(true=st, model=sm, lam=lam.tolist(), acc=a),
                  open(args.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
