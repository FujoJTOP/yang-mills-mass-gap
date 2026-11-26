#!/usr/bin/env python3
r"""
o12_sun2.py -- O12 candidate 3, WITH A SIGN CORRECTION to the archive's Ddyn estimator.

FOUND BUG in na2_su3_na.py (not edited, reported instead):
  code:  Ddyn = -dEG - X/2 ,  dEG = sum w*s*(G-EG)/Z ,  s=(a-a0)/h^2
  and  s = dK/da0 / K  =>  dEG = +dE[G]/da exactly.
  conditional identity:  dE[G]/da = E[H] - Var(G) - Cov(G,Phi')  =>  Ddyn = E[H]-Var(G)
                                                                    = dEG - X/2 .
  So the code's sign is flipped:  Ddyn_code = -Ddyn_true - X .  (X ~ 0 here, so Ddyn_code ~ -Ddyn_true.)
  Validation below (N=2, exact quadrature o12_su2_theta.log): theta=0,beta=2.5 exact Delta=+0.72992
  with E[H]=0.55306, Var=0, Wgeo=0.176864.

  Delta_true = Ddyn_true + Wgeo + X = dEG - X/2 + Wgeo + X .
Run: py -3.12 -u -I o12_sun2.py [--n 150000]
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import na2_su3_na as M

def stats(geo, ct, beta, a0, h):
    N = geo['N']
    S = beta * (1.0 - ct[0] / N)
    G = -beta * ct[1] / N
    a = geo['a']
    w = geo['wg'] * np.exp(-S) * np.exp(-0.5 * ((a - a0) / h) ** 2)
    Z = w.sum()
    if not np.isfinite(Z) or Z <= 0:
        return None
    wn = w / Z
    E = lambda x: float(np.sum(wn * x))
    EG, EP = E(G), E(geo['P1'])
    VarG = E(G * G) - EG ** 2
    VarP = E(geo['P1'] ** 2) - EP ** 2
    Cov = E(G * geo['P1']) - EG * EP
    W = E(geo['P2']) - VarP
    X = -2 * Cov
    s = (a - a0) / h ** 2
    dEG = float(np.sum(w * s * (G - EG))) / Z         # = dE[G]/da0
    Ddyn = dEG - 0.5 * X                              # corrected identity
    return dict(dEG=dEG, Ddyn=Ddyn, Ddyn_code=-dEG - 0.5 * X,
                Wgeo=W, X=X, Gpp=Ddyn + W + X, Gpp_code=-dEG - 0.5 * X + W + X,
                neff=1.0 / np.sum(wn ** 2), EG=EG)

def ext(geo, ct, beta, a0=0.0, hs=(0.12, 0.045)):
    r1, r2 = stats(geo, ct, beta, a0, hs[0]), stats(geo, ct, beta, a0, hs[1])
    if r1 is None or r2 is None:
        return None
    o = {k: r2[k] + (r2[k] - r1[k]) / 3.0 for k in ('dEG', 'Ddyn', 'Wgeo', 'X', 'Gpp', 'Gpp_code')}
    o['neff'] = r1['neff']
    return o

def main():
    n = 150000
    for i, a in enumerate(sys.argv):
        if a == '--n':
            n = int(sys.argv[i + 1])
    t0 = time.time(); rng = np.random.default_rng(20261125)
    print("=" * 122, flush=True)
    print("O12/CAND 3 (corrected Ddyn) : SU(N) one plaquette, a0=0, v=e3.  Delta = G''(0)", flush=True)
    print("=" * 122, flush=True)
    H = (f"     {'holonomy':>26}{'beta':>6}{'dEG':>10}{'Wgeo':>9}{'X':>10}"
         f"{'Delta_true':>12}{'Delta_code':>12}{'sign':>5}{'neff':>8}")
    print("\n[V] N=2 control vs EXACT quadrature (o12_su2_theta.log):", flush=True)
    print("     (theta=0, b=2.5) exact +0.72992 ; (theta=pi, b=2.5) exact +0.13599 ; (theta=pi,b=5) -0.07137",
          flush=True)
    geo2 = M.geometry(2, n, rng)
    print(H, flush=True)
    rows = []
    for th, lab in ((0.0, 'th=0'), (np.pi, 'th=pi')):
        Bd = np.diag(np.exp(1j * np.array([th / 2, -th / 2])))
        ct = M.action_traces(geo2, [Bd])[0]
        for beta in (2.5, 5.0):
            r = ext(geo2, ct, beta)
            print(f"     {lab:>26}{beta:6.2f}{r['dEG']:10.4f}{r['Wgeo']:9.5f}{r['X']:+10.5f}"
                  f"{r['Gpp']:+12.5f}{r['Gpp_code']:+12.5f}{('NEG' if r['Gpp'] < 0 else 'pos'):>5}"
                  f"{r['neff']:8.0f}", flush=True)

    geo3 = M.geometry(3, n, rng)
    print(f"\n[N=3]  a in [{geo3['a'].min():.3f},{geo3['a'].max():.3f}]  ({time.time()-t0:.0f}s)", flush=True)
    print(H, flush=True)
    ths3 = [(0.0, 0.0, 0.0), (0.3, 1.1, -1.4), (2 * np.pi / 3, 2 * np.pi / 3, -4 * np.pi / 3)]
    for t in ths3:
        ct = M.action_traces(geo3, [np.diag(np.exp(-1j * np.array(t)))])[0]
        for beta in (2.5, 4.0, 5.5, 6.0, 6.5, 8.0, 12.0, 20.0):
            r = ext(geo3, ct, beta)
            print(f"     {str(np.round(t,3)):>26}{beta:6.2f}{r['dEG']:10.4f}{r['Wgeo']:9.5f}"
                  f"{r['X']:+10.5f}{r['Gpp']:+12.5f}{r['Gpp_code']:+12.5f}"
                  f"{('NEG' if r['Gpp'] < 0 else 'pos'):>5}{r['neff']:8.0f}", flush=True)
    print(f"\n  [{time.time()-t0:.0f}s]  Clay 0%. No QED.", flush=True)

if __name__ == '__main__':
    main()
