#!/usr/bin/env python3
"""
na2_t0_fit.py -- T0 cross-check by an INDEPENDENT estimator of X(a0=0):
   instead of the h^2 Richardson extrapolation of a kernel regression, evaluate X(a) at
   a = 0, .25, .5, .75, 1 with a small fixed bandwidth and fit  X(a) = A + B a^2 (+ C a^4);
   the intercept A is an independent estimate of X(0).
Also prints the model-A (N=2 exact) analogue for the same fit, as a control.

Run: py -3.12 -I na2_t0_fit.py [--n 1200000]
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import na2_su3_na as M

def x_at(geo, ct, beta, a0, h):
    r = M.cond_stats(geo, ct, beta, a0, h)
    return r['X'], r['seX'], r['neff']

def fit_X0(geo, ct, beta, h=0.05, aa=(0.0, 0.25, 0.5, 0.75, 1.0), deg=2):
    """fit X(a) = sum_{k<=deg/2} c_k a^{2k}; return c_0 and its error from the point scatters."""
    xs, ys, es = [], [], []
    for a0 in aa:
        X, se, ne = x_at(geo, ct, beta, a0, h)
        xs.append(a0); ys.append(X); es.append(max(se, 1e-12))
    A = np.vstack([np.array(xs) ** (2 * k) for k in range(deg // 2 + 1)]).T
    W = np.diag(1.0 / np.array(es) ** 2)
    cov = np.linalg.inv(A.T @ W @ A)
    c = cov @ (A.T @ W @ np.array(ys))
    chi2 = float(np.sum(((np.array(ys) - A @ c) / np.array(es)) ** 2))
    return c[0], float(np.sqrt(cov[0, 0])), chi2, ys

def main():
    n = 1200000
    for i, a in enumerate(sys.argv):
        if a == '--n':
            n = int(sys.argv[i + 1])
    t0 = time.time()
    print("=" * 104, flush=True)
    print("T0 CROSS-CHECK: X(a0) at a0=0,.25,.5,.75,1 (kernel h=0.05) fitted as A + B a^2 (+ C a^4);", flush=True)
    print("                A is an independent estimate of X(0).", flush=True)
    print("=" * 104, flush=True)
    rng = np.random.default_rng(777)

    # ---- N=2 control (model A): exact X(0)=0 for every theta,beta
    geo2 = M.geometry(2, n, rng)
    print(f"\n[N=2 control, model A: exact X(0)=0, exact X(1)=0.00222 (th=0,beta=2.5)]  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    ct02 = M.action_traces(geo2, [np.diag(np.exp(+1j * np.array([t / 2, -t / 2]))) for t in
                                  (0.0, np.pi / 2, np.pi)])
    for lab, ct in zip(('th=0', 'th=pi/2', 'th=pi'), ct02):
        for beta in (2.5,):
            A, sA, chi2, ys = fit_X0(geo2, ct, beta)
            print(f"     {lab:>8} beta={beta:4.1f}:  X(0) = {A:+.3e} +- {sA:.1e}   chi2/dof={chi2/3:.1f}"
                  f"   X(a) at a=.25,.5,.75,1 = " + " ".join(f"{v:+.4f}" for v in ys[1:]), flush=True)

    # ---- N=3
    print(f"\n[N=3] the real test  [{time.time()-t0:.0f}s]", flush=True)
    geo3 = M.geometry(3, n, rng)
    ths3 = [(0.0, 0.0, 0.0), (0.3, 1.1, -1.4), (np.pi / 2, np.pi / 2, -np.pi),
            (np.pi, np.pi, -2 * np.pi), (2 * np.pi / 3, 2 * np.pi / 3, -4 * np.pi / 3)]
    cts3 = M.action_traces(geo3, [np.diag(np.exp(-1j * np.array(t))) for t in ths3])
    out = []
    print(f"     {'(th1,th2,th3)':>24}{'beta':>6}{'X(0) [fit]':>13}{'err':>9}{'signif':>8}"
          f"{'X(0) [h-ext]':>14}{'chi2/dof':>10}", flush=True)
    for lab, ct in zip([str(np.round(t, 3)) for t in ths3], cts3):
        for beta in (0.5, 1.0, 2.5, 4.0):
            A, sA, chi2, ys = fit_X0(geo3, ct, beta)
            r = M.extrap(geo3, ct, beta, 0.0)
            print(f"     {lab:>24}{beta:6.2f}{A:+13.3e}{sA:9.1e}{A/sA:8.1f}{r['X']:+14.3e}"
                  f"{chi2/3:10.1f}", flush=True)
            out.append(dict(lab=lab, beta=beta, X0_fit=A, err=sA, X0_extrap=r['X'], chi2=chi2))
    with open('na2_t0_fit.json', 'w') as f:
        json.dump(out, f, indent=1, default=float)
    print(f"\n  [{time.time()-t0:.0f}s]  Clay 0%. No QED.", flush=True)

if __name__ == '__main__':
    main()
