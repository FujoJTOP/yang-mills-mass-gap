"""Analyse a pilot_numba.py output: cosh fit of the connected correlator,
effective mass, tau_int and effective statistics."""
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def conn(Sm):
    n, L = Sm.shape
    X = Sm - Sm.mean()
    return np.array([np.mean(X * np.roll(X, -t, axis=1)) for t in range(L)])


def fit_m(C, L, tmin=1, tmax=None):
    if tmax is None:
        tmax = min(10, L // 2)
    ts = np.arange(tmin, tmax + 1)
    y = C[ts]
    ms = np.linspace(0.02, 2.0, 6000)
    Z = np.exp(-ms[:, None] * ts[None, :]) + np.exp(-ms[:, None] * (L - ts)[None, :])
    A = (Z @ y) / np.sum(Z * Z, axis=1)
    r = np.sum((y[None, :] - A[:, None] * Z) ** 2, axis=1)
    k = int(np.argmin(r))
    if k == 0 or k == len(ms) - 1:
        return float('nan'), float('nan')
    return float(ms[k]), float(A[k])


def meff(C, L):
    out = []
    for t in range(1, min(L // 2, 15) + 1):
        cm, c0, cp = C[(t - 1) % L], C[t % L], C[(t + 1) % L]
        if min(cm, c0, cp) > 1e-18:
            out.append((t, float(np.arccosh((cm + cp) / (2 * c0)))))
        else:
            out.append((t, np.nan))
    return out


def tau_int(x):
    x = np.asarray(x, float); n = len(x); x = x - x.mean()
    v = np.dot(x, x) / n
    if v <= 0:
        return 0.5
    s = 0.0
    for t in range(1, min(n // 10, 1000)):
        s += np.dot(x[:-t], x[t:]) / (n - t) / v
        if t >= 6 * (0.5 + s):
            break
    return max(0.5 + s, 0.5)


def jk(fn, n, b):
    b = max(b, int(np.ceil(n / 200.0)))
    nb = n // b
    bs = [np.arange(i * b, (i + 1) * b) for i in range(nb)]
    full = fn(np.arange(n))
    vals = np.array([fn(np.concatenate([bs[j] for j in range(nb) if j != i]))
                     for i in range(nb)])
    err = np.sqrt((nb - 1) / nb * np.sum((vals - vals.mean(axis=0)) ** 2, axis=0))
    return full, err


def main():
    files = sys.argv[1:] or sorted(glob.glob(os.path.join(HERE, 'runs', 'pilot_*.npz')))
    for f in files:
        d = dict(np.load(f, allow_pickle=True))
        L = int(d['L']); beta = float(d['beta']); nlev = int(d['nlev'])
        tau = tau_int(d['plaq'])
        nm = len(d['S0'])
        for tag, key in (('n=0', 'S0'), (f'n={nlev}', 'Sn')):
            S = d[key]
            C, eC = jk(lambda idx: conn(S[idx]), nm, max(1, int(np.ceil(2 * tau))))
            m, em = jk(lambda idx: np.array([fit_m(conn(S[idx]), L)[0]]), nm,
                       max(1, int(np.ceil(2 * tau))))
            m = float(m[0]); em = float(em[0])
            C = np.asarray(C)
            print(f"=== {os.path.basename(f)}  L={L} beta={beta:.3f}  {tag}")
            print(f"    tau_int(P)={tau:.2f}  N_meas={nm}  N_eff~{nm/(2*tau):.0f}")
            print(f"    C/C0(t=1..6) = {np.round(C[1:7]/C[0],4)} +- {np.round(eC[1:7]/C[0],4)}")
            print(f"    cosh fit (t=1..{min(10, L//2)}): am={m:.4f}+-{em:.4f}  "
                  f"xi/a={1/m:.3f}+-{em/m**2:.3f}")
            print(f"    m_eff = {[(t, round(v,3)) for t, v in meff(C, L)][:8]}")
            # autocorrelation of the correlator entries (glueball slow mode)
            for tt in (1, 2, 3):
                tv = tau_int(S[:, tt])
                print(f"    tau_int(S[t={tt}]) = {tv:.2f}  -> N_eff there = {nm/(2*tv):.0f}")


if __name__ == '__main__':
    main()
