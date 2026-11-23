"""
Analyse tm_gap_run.py output: connected Euclidean-time correlator of the
0++ operator and its effective-mass plateau (the transfer-matrix gap Delta).

C(t) = <O(t') O(t'+t)>_c   (connected, time-translation + config averaged,
                           periodic time).  C(t) is symmetric, C(t)=C(L-t).

Spectral form on a periodic time direction of length L:
  C(t) = sum_k A_k [ e^{-Delta_k t} + e^{-Delta_k (L-t)} ]
so the lightest state's Delta is read from
  m_eff(t) = arccosh( [C(t-1)+C(t+1)] / [2 C(t)] )  -> Delta   (plateau),
and from a 1-state cosh fit C(t)=A(e^{-m t}+e^{-m(L-t)}) over t=1..L/2.
Errors: bin the configs by 2*tau_int(plaquette) and jackknife.
"""
import glob
import sys

import numpy as np


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


def conn(O):
    """O: (nconf, L).  Connected correlator C(t), t=0..L-1 (periodic)."""
    n, L = O.shape
    X = O - O.mean()
    return np.array([np.mean(X * np.roll(X, -t, axis=1)) for t in range(L)])


def meff(C, L):
    """cosh effective mass, t=1..L//2 (NaN where signal non-positive)."""
    out = {}
    for t in range(1, L // 2 + 1):
        cm, c0, cp = C[(t - 1) % L], C[t % L], C[(t + 1) % L]
        if min(c0, cm + cp) > 1e-14 and (cm + cp) / (2 * c0) >= 1.0:
            out[t] = float(np.arccosh((cm + cp) / (2 * c0)))
        else:
            out[t] = np.nan
    return out


def fit1cosh(C, L, tmin=1, tmax=None):
    """1-state cosh fit C(t)=A[e^{-m t}+e^{-m(L-t)}], grid search on m."""
    if tmax is None:
        tmax = L // 2
    ts = np.arange(tmin, tmax + 1)
    y = C[ts]
    ms = np.linspace(0.02, 5.0, 6000)
    Z = np.exp(-ms[:, None] * ts[None, :]) + np.exp(-ms[:, None] * (L - ts)[None, :])
    A = (Z @ y) / np.sum(Z * Z, axis=1)
    r = np.sum((y[None, :] - A[:, None] * Z) ** 2, axis=1)
    k = int(np.argmin(r))
    m = float(ms[k])
    return (m, float(A[k])) if 0 < k < len(ms) - 1 else (np.nan, np.nan)


def jack(fn, n, b, nbin_max=200):
    """full estimate and jackknife error of fn(indices).
    bin b >= 2*tau_int for decorrelation; also cap the bin count at nbin_max
    so the jackknife stays cheap (bin size grown, still >> tau_int)."""
    b = max(b, 1, n // nbin_max)
    nb = n // b
    if nb < 2:
        return fn(np.arange(n)), np.nan
    bs = [np.arange(i * b, min((i + 1) * b, n)) for i in range(nb)]
    full = fn(np.arange(n))
    vals = [fn(np.concatenate([bs[j] for j in range(nb) if j != i])) for i in range(nb)]
    vals = np.array(vals)
    err = np.sqrt((nb - 1) / nb * np.sum((vals - vals.mean(axis=0)) ** 2, axis=0))
    return full, err


def analyse_file(f):
    d = dict(np.load(f, allow_pickle=True))
    L = int(d['L']); beta = float(d['beta']); levels = [int(x) for x in d['levels']]
    plaq = d['plaq']; n = len(plaq)
    ti = tau_int(plaq)
    b = max(1, int(np.ceil(2 * ti)))
    print(f"\n=== {f}  L={L} beta={beta} N={n}  tau_int={ti:.2f} bin={b}")
    rows = {}
    for lv in levels:
        O = d[f'O{lv}']
        C = conn(O)
        me = meff(C, L)
        # jackknife errors on m_eff(t) and on the 1-cosh fit mass
        def fn_me(idx, t=None):
            return np.array([meff(conn(O[idx]), L)[t] for t in range(1, L // 2 + 1)])
        full_me, err_me = jack(lambda idx: fn_me(idx), n, b)
        def fn_fit(idx):
            return np.array([fit1cosh(conn(O[idx]), L)[0]])
        full_fit, err_fit = jack(fn_fit, n, b)
        tmax = L // 2
        full_fit_long, err_fit_long = jack(
            lambda idx: np.array([fit1cosh(conn(O[idx]), L, tmin=max(1, L // 2 - 1))[0]]),
            n, b)
        full_fit = float(np.atleast_1d(full_fit)[0])
        err_fit = float(np.atleast_1d(err_fit)[0])
        full_fit_long = float(np.atleast_1d(full_fit_long)[0])
        err_fit_long = float(np.atleast_1d(err_fit_long)[0])
        # signal-to-noise of C
        Csn = np.array([C[t] / (abs(C[t]) + 1e-300) for t in range(L)])
        print(f" level n={lv}:")
        print(f"   C/C0 = {np.round(C / C[0], 4)}")
        me_str = ", ".join(f"t={t}:{me[t]:.3f}" for t in me)
        print(f"   m_eff(t) = {me_str}")
        print(f"   m_eff +- {np.round(np.nan_to_num(err_me), 3)}")
        print(f"   1cosh fit m = {full_fit:.4f} +- {err_fit:.4f}   "
              f"(window t>= {max(1, L // 2 - 1)}: {full_fit_long:.4f} +- {err_fit_long:.4f})")
        rows[lv] = dict(C=C, me=me, err_me=np.nan_to_num(err_me),
                        fit_m=float(full_fit), fit_e=float(err_fit),
                        fit_m_long=float(full_fit_long), fit_e_long=float(err_fit_long))
    return dict(L=L, beta=beta, n=n, bin=b, tau_int=ti, rows=rows)


if __name__ == '__main__':
    files = sys.argv[1:]
    if not files:
        files = sorted(glob.glob('tm_gap_data/*.npz'))
    for f in files:
        if '_smoke' in f:
            pass
        analyse_file(f)
