#!/usr/bin/env python3
r"""
o12_marg.py -- O12: measure Delta = G''(0) = -d^2/da^2 log rho_v(a)|_{a=0} DIRECTLY from the
marginal density (histogram + local quadratic), no decomposition, no derivative identity.
That sidesteps the Ddyn-sign ambiguity of na2_su3_na.py (reported, not edited).

rho(a) = sum_i w_i delta(a_i - a) ,  w_i = wg_i exp(-S_i)   (Weyl+action weight; kernel-free!)
Fit log rho(a) ~ c0 + c1 (a-a0) + c2 (a-a0)^2 on |a-a0| < H  =>  Delta = -2 c2.

VALIDATION (exact quadrature, o12_su2_theta.log, SU(2) model A):
   theta=0 , beta=2.5 -> +0.72992      theta=pi, beta=2.5 -> +0.13599
   theta=pi, beta=5.0 -> -0.07137      theta=pi, beta=0  -> +0.20512
Run: py -3.12 -u -I o12_marg.py [--n 150000]
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import na2_su3_na as M

def gpp_marg(a, w, a0=0.0, H=0.4, bw=0.08):
    """-d^2/da^2 log rho at a0 via histogram + local quadratic (weighted)."""
    m = np.abs(a - a0) < H + bw
    a, w = a[m], w[m]
    nb = max(int(np.ceil(2 * H / bw)), 6)
    edges = np.linspace(a0 - H, a0 + H, nb + 1)
    idx = np.clip(np.digitize(a, edges) - 1, 0, nb - 1)
    cnt = np.bincount(idx, weights=w, minlength=nb)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    ok = cnt > 1e-12
    if ok.sum() < 5:
        return float('nan'), 0
    x = ctr[ok] - a0
    y = np.log(cnt[ok])
    W = cnt[ok]
    A = np.stack([np.ones_like(x), x, x * x], axis=1)
    Aw = A * W[:, None]
    coef = np.linalg.solve(A.T @ Aw, Aw.T @ y)
    return -2.0 * coef[2], int(ok.sum())

def main():
    n = 150000
    for i, a in enumerate(sys.argv):
        if a == '--n':
            n = int(sys.argv[i + 1])
    t0 = time.time(); rng = np.random.default_rng(20261125)
    print("=" * 118, flush=True)
    print("O12 : Delta = G''(0) from the MARGINAL density (model-free).  a0=0, v=e3.", flush=True)
    print("=" * 118, flush=True)
    Hd = (f"     {'holonomy':>26}{'beta':>6}" + "".join(f"{'H='+str(h):>11}" for h in (0.6, 0.4, 0.25))
          + f"{'Delta(H=.4)':>13}{'sign':>5}{'bins':>6}")

    geo2 = M.geometry(2, n, rng)
    print("\n[V] N=2 vs exact: th=0,b=2.5 +0.72992 | th=pi,b=2.5 +0.13599 | th=pi,b=5 -0.07137", flush=True)
    print(Hd, flush=True)
    for th in (0.0, np.pi):
        Bd = np.diag(np.exp(1j * np.array([th / 2, -th / 2])))
        c0, c1 = M.action_traces(geo2, [Bd])[0]
        a = geo2['a']
        for beta in (2.5, 5.0):
            w = geo2['wg'] * np.exp(-beta * (1.0 - c0 / 2.0))
            vs = [gpp_marg(a, w, 0.0, H) for H in (0.6, 0.4, 0.25)]
            print(f"     {('th=pi' if th else 'th=0'):>26}{beta:6.2f}"
                  + "".join(f"{v[0]:+11.4f}" for v in vs)
                  + f"{vs[1][0]:+13.4f}{('NEG' if vs[1][0] < 0 else 'pos'):>5}{vs[1][1]:6d}", flush=True)

    geo3 = M.geometry(3, n, rng)
    print(f"\n[N=3]  a in [{geo3['a'].min():.3f},{geo3['a'].max():.3f}]  ({time.time()-t0:.0f}s)", flush=True)
    print(Hd, flush=True)
    for t in [(0.0, 0.0, 0.0), (0.3, 1.1, -1.4), (2 * np.pi / 3, 2 * np.pi / 3, -4 * np.pi / 3)]:
        c0, c1 = M.action_traces(geo3, [np.diag(np.exp(-1j * np.array(t)))])[0]
        a = geo3['a']
        for beta in (2.5, 4.0, 6.0, 8.0, 12.0):
            w = geo3['wg'] * np.exp(-beta * (1.0 - c0 / 3.0))
            vs = [gpp_marg(a, w, 0.0, H) for H in (0.6, 0.4, 0.25)]
            print(f"     {str(np.round(t,3)):>26}{beta:6.2f}"
                  + "".join(f"{v[0]:+11.4f}" for v in vs)
                  + f"{vs[1][0]:+13.4f}{('NEG' if vs[1][0] < 0 else 'pos'):>5}{vs[1][1]:6d}", flush=True)
    print(f"\n  [{time.time()-t0:.0f}s]  Clay 0%. No QED.", flush=True)

if __name__ == '__main__':
    main()
