"""Final tables for the report (fits + validation)."""
import glob
import json
import os

import numpy as np


def fit_table():
    print("FITS")
    print(f"{'tag':18s} {'lambda':38s} {'scatter(last3)':22s} {'cond':>6s} {'nit':>4s} {'acc':>5s}")
    for f in sorted(glob.glob('fits/L*.json')):
        d = json.load(open(f))
        hist = d['hist']
        lams = np.array([h['lam'] for h in hist])
        lam = lams[-1]
        sc = lams[-3:].std(axis=0) if len(lams) >= 3 else np.zeros_like(lam)
        print(f"{os.path.basename(f)[:-5]:18s} "
              f"{np.array2string(lam, precision=4, floatmode='fixed'):38s} "
              f"{np.array2string(sc, precision=4, floatmode='fixed'):22s} "
              f"{d['cond_true']:6.2f} {len(hist):4d} {hist[-1]['acc']:5.2f}")


def val_table():
    print("\nVALIDATION  (* = held out)")
    for f in sorted(glob.glob('fits/val_*.json')):
        d = json.load(open(f))
        print(f"\n{os.path.basename(f)[4:-5]}: lam={np.round(d['lam'],4)}")
        for r in d['rows']:
            mark = '' if r['fit'] else ' *'
            print(f"   {r['obs']:9s} true={r['true']:+10.6f}({r['err_true']:.6f}) "
                  f"mod={r['model']:+10.6f}({r['err_model']:.6f}) "
                  f"rel={100*r['rel']:+7.1f}% z={r['z']:+6.1f}{mark}")


if __name__ == '__main__':
    fit_table()
    val_table()
