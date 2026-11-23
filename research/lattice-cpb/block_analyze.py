"""Analyse block_obs.py output.

Reports, per run:
  * character coefficients c3, c6, c8, c10 with signs and errors
    -> explicit answer to "does any c_R < 0 exist?"
  * moments (mean, variance, skewness) of the blocked plaquette field P_b
  * defect counts k at thresholds P<0, P<P>/2, P<0.75P> with Poisson errors,
    the per-block fraction f, and c0 = -ln(f) with error (or a 95% lower bound
    when k=0).  This is a defect-fraction proxy, NOT w_q and NOT c0 itself.
"""
import glob
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def jack(fn, n, b):
    b = max(b, int(np.ceil(n / 400.0)))
    nb = n // b
    bs = [np.arange(i * b, (i + 1) * b) for i in range(nb)]
    full = fn(np.arange(n))
    vals = np.array([fn(np.concatenate([bs[j] for j in range(nb) if j != i]))
                     for i in range(nb)])
    err = np.sqrt((nb - 1) / nb * np.sum((vals - vals.mean(axis=0)) ** 2, axis=0))
    return full, err


def tau_int(x):
    x = np.asarray(x, float); n = len(x); x = x - x.mean()
    v = np.dot(x, x) / n
    if v <= 0:
        return 0.5
    s = 0.0
    for t in range(1, min(n // 10, 500)):
        s += np.dot(x[:-t], x[t:]) / (n - t) / v
        if t >= 6 * (0.5 + s):
            break
    return max(0.5 + s, 0.5)


def main():
    files = sorted(glob.glob(os.path.join(HERE, 'runs', 'block_*.npz')) +
                   glob.glob(os.path.join(HERE, 'runs', 'b8_*.npz')))
    out = []
    for f in files:
        d = dict(np.load(f, allow_pickle=True))
        L = int(d['L']); beta = float(d['beta']); n = len(d['plaq'])
        b = max(1, int(np.ceil(2 * tau_int(d['plaq']))))
        rows = [f"=== L={L} beta={beta:.3f}  N={n} bin={b}  "
                f"[{os.path.basename(f)}]"]

        # ---- character coefficients and their signs
        neg = []
        for k in ('c3', 'c6', 'c8', 'c10'):
            if k not in d:
                continue
            val, err = jack(lambda idx, k=k: float(np.mean(d[k][idx])), n, b)
            rows.append(f"   {k:4s} = {val:+.5f} +- {err:.5f}")
            if val < 0:
                neg.append(k)
        if neg:
            rows.append(f"   *** SIGN STRUCTURE PRESENT: {neg} < 0 ***")
        else:
            rows.append("   ALL c_R > 0  -> no c1-level sign cancellation; "
                        "consistent with BL-6 (kappa_c~1), not a measurement of it")

        # ---- moments of the blocked field
        V = {1: L ** 4, 2: (L // 2) ** 4}
        for bb in (1, 2):
            if f'mean_{bb}' not in d:
                continue
            mu, emu = jack(lambda idx, bb=bb: float(np.mean(d[f'mean_{bb}'][idx])), n, b)
            va, eva = jack(lambda idx, bb=bb: float(np.mean(d[f'var_{bb}'][idx])), n, b)
            sk, esk = jack(lambda idx, bb=bb: float(np.mean(d[f'skew_{bb}'][idx])), n, b)
            rows.append(f"   b={bb}: mean={mu:.5f}+-{emu:.5f}  var={va:.6f}+-{eva:.6f}"
                        f"  skew={sk:+.4f}+-{esk:.4f}  (V={V[bb]} blocks)")

        # ---- Gaussian-tail prediction for the same mean/variance
        from math import erfc, sqrt
        for bb in (1, 2):
            if f'mean_{bb}' not in d:
                continue
            mu = float(np.mean(d[f'mean_{bb}']))
            sd = float(np.sqrt(np.mean(d[f'var_{bb}'])))
            pred = []
            for tag, th in (('P<0', 0.0), ('P<P>/2', 0.5 * mu), ('P<0.75P>', 0.75 * mu)):
                z = (mu - th) / (sd * sqrt(2.0)) if sd > 0 else np.inf
                pred.append(f"{tag}: {0.5 * erfc(z):.5f}")
            rows.append(f"   b={bb} Gaussian-tail prediction from the SAME moments: "
                        + ", ".join(pred))

        # ---- defect counts, Poisson errors, c0 proxy
        for bb in (1, 2):
            for tag, lab in (('0', 'P<0'), ('hh', 'P<P>/2'), ('h75', 'P<0.75P>')):
                if f'k_{bb}_{tag}' not in d:
                    continue
                k = int(np.sum(d[f'k_{bb}_{tag}']))
                Ntot = V[bb] * n
                f_ = k / Ntot
                if k > 0:
                    ef = np.sqrt(k) / Ntot
                    c0 = -np.log(f_)
                    ec0 = 1.0 / np.sqrt(k)
                    rows.append(f"   b={bb} {lab}: k={k}  f={f_:.5f}+-{ef:.5f}"
                                f"  c0=-ln f = {c0:.3f} +- {ec0:.3f}")
                else:
                    f95 = 3.0 / Ntot
                    rows.append(f"   b={bb} {lab}: k=0  f<{f95:.5f} (95%)"
                                f"  c0 > {-np.log(f95):.3f} (95% lower bound)")
        out.append('\n'.join(rows))
    txt = '\n\n'.join(out)
    print(txt)
    open(os.path.join(HERE, 'block_summary.txt'), 'w').write(txt + '\n')


if __name__ == '__main__':
    main()
