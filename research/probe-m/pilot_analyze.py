"""Analyse pilot_run.py output: connected correlator, cosh fit, jackknife xi/a.

Same fit as lattice-cpb/corr_analyze.py (grid m in [0.05,4], cosh, t=1..L/2),
but with no matplotlib dependency.
"""
import argparse
import glob
import os

import numpy as np


def conn(Sm):
    """Connected correlator C(t) from time-slice sums, shape (n, L)."""
    X = Sm - Sm.mean()
    L = X.shape[1]
    return np.array([np.mean(X * np.roll(X, -t, axis=1)) for t in range(L)])


def fit_xi(Sm, tmin=1, tmax=None, block=None):
    L = Sm.shape[1]
    if tmax is None:
        tmax = L // 2
    ts = np.arange(tmin, tmax + 1)
    C = conn(Sm)
    y = C[ts]
    ms = np.linspace(0.05, 4.0, 4000)
    Z = np.exp(-ms[:, None] * ts[None, :]) + np.exp(-ms[:, None] * (L - ts)[None, :])

    def one(yv):
        A = (Z @ yv) / np.sum(Z * Z, axis=1)
        r = np.sum((yv[None, :] - A[:, None] * Z) ** 2, axis=1)
        k = int(np.argmin(r))
        return ms[k], A[k], r[k]

    m, A, r = one(y)
    # jackknife over config blocks
    n = Sm.shape[0]
    b = block or max(2, n // 40)
    nb = n // b
    vals = []
    for i in range(nb):
        keep = np.ones(n, bool)
        keep[i * b:(i + 1) * b] = False
        mi, _, _ = one(conn(Sm[keep])[ts])
        vals.append(mi)
    vals = np.array(vals)
    err = np.sqrt((nb - 1) / nb * np.sum((vals - vals.mean()) ** 2))
    return m, err, C, (tmin, tmax)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('paths', nargs='+')
    ap.add_argument('--level', default='S8', help='S0, S8, S16, S32, ...')
    ap.add_argument('--tmin', type=int, default=1)
    ap.add_argument('--tmax', type=int, default=0, help='0 = L//2')
    args = ap.parse_args()

    paths = []
    for p in args.paths:
        paths += glob.glob(p)

    print(f"{'beta':>6} {'L':>3} {'level':>5} {'C1/C0':>9} {'C2/C0':>9} "
          f"{'am':>12} {'xi/a':>14} {'P':>9}")
    for p in sorted(paths):
        d = np.load(p)
        S = d[args.level]
        L = int(d['L'])
        m, e, C, _ = fit_xi(S, tmin=args.tmin, tmax=args.tmax or None)
        # jackknife error on C(1) to flag a signal consistent with zero
        n = S.shape[0]
        b = max(2, n // 40)
        nb = n // b
        c1 = []
        for i in range(nb):
            keep = np.ones(n, bool)
            keep[i * b:(i + 1) * b] = False
            Cc = conn(S[keep])
            c1.append(Cc[1] / Cc[0])
        c1 = np.array(c1)
        c1e = np.sqrt((nb - 1) / nb * np.sum((c1 - c1.mean()) ** 2))
        if C[1] <= 2 * c1e:
            print(f"{float(d['beta']):6.2f} {L:3d} {args.level:>5} "
                  f"{C[1]/C[0]:9.4f} {C[2]/C[0]:9.4f} "
                  f"{'unresolved':>12} {'':>14} {np.mean(d['plaq']):9.5f}")
            continue
        xi, xie = 1.0 / m, e / m**2
        print(f"{float(d['beta']):6.2f} {L:3d} {args.level:>5} "
              f"{C[1]/C[0]:9.4f} {C[2]/C[0]:9.4f} "
              f"{m:12.4f} {xi:8.3f}({xie:.3f}) {np.mean(d['plaq']):9.5f}")


if __name__ == '__main__':
    main()
