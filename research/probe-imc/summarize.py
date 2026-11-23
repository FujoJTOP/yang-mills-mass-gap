"""Assemble report tables from fits/*.json."""
import glob
import json
import os

import numpy as np


def main():
    print("== fits")
    for f in sorted(glob.glob('fits/L*.json')):
        d = json.load(open(f))
        hist = d['hist']
        lams = np.array([h['lam'] for h in hist])
        lam = lams[-1]
        scat = lams[len(lams) // 2:].std(axis=0) if len(lams) > 2 else np.zeros_like(lam)
        print(f"{os.path.basename(f):22s} lam={np.round(lam,4)} "
              f"scatter={np.round(scat,4)} cond_true={d['cond_true']:.2f} "
              f"nit={len(hist)} acc={hist[-1]['acc']:.2f}")
    print("\n== validation")
    for f in sorted(glob.glob('fits/val_*.json')):
        d = json.load(open(f))
        print(os.path.basename(f))
        for r in d['rows']:
            flag = '' if r['fit'] else ' *'
            print(f"   {r['obs']:9s} true={r['true']:+.6f}({r['err_true']:.6f}) "
                  f"mod={r['model']:+.6f}({r['err_model']:.6f}) "
                  f"rel={100*r['rel']:+7.1f}% z={r['z']:+6.1f}{flag}")


if __name__ == '__main__':
    main()
