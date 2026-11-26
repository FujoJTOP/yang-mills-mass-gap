#!/usr/bin/env python3
r"""
o12_sun.py -- O12 candidate 3: does Delta = E[curv|0]-Var[force|0] = G''(0) go NEGATIVE for N>=3
at the physical coupling?  Reuses the validated SU(N) one-plaquette machinery of
na2_su3_na.py (flat chart = principal log, Weyl measure, kernel conditioning a0=0, v=e3).

  Delta = Ddyn + Wgeo + X = G''(0),  Ddyn = -dE[G]/da - X/2 (exact identity),  X = -2Cov(G,Phi').

Controls: N=2 must reproduce the exact quadrature of o12_su2_theta.log
(theta=pi,beta=2.5 -> +0.13599 ; theta=0,beta=2.5 -> +0.72992).
Run: py -3.12 -u -I o12_sun.py [--n 200000]
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import na2_su3_na as M

def run(geo, ths, betas, a0=0.0, hs=(0.12, 0.06)):
    cts = M.action_traces(geo, [np.diag(np.exp(-1j*np.array(t))) for t in ths])
    out = []
    for lab, ct in zip(ths, cts):
        for beta in betas:
            r = M.extrap(geo, ct, beta, a0, hs)
            out.append((lab, beta, r))
            print(f"     {str(np.round(lab,3)):>26}{beta:6.2f}{r['Ddyn']:+11.5f}{r['Wgeo']:+9.5f}"
                  f"{r['X']:+10.5f}{r['Gpp']:+11.5f}{r['Gpp_marg']:+11.5f}"
                  f"{('NEG' if r['Gpp']<0 else 'pos'):>5}{r['neff']:9.0f}{r['seX']:9.4f}", flush=True)
    return out

def main():
    n = 200000
    for i, a in enumerate(sys.argv):
        if a == '--n':
            n = int(sys.argv[i+1])
    t0 = time.time()
    print("="*126, flush=True)
    print("O12/CAND 3 : SU(N) one plaquette, a0=0, v=e3.  Delta = G''(0) = Ddyn+Wgeo+X", flush=True)
    print("="*126, flush=True)
    rng = np.random.default_rng(20261125)
    hdr = (f"     {'holonomy':>26}{'beta':>6}{'Ddyn':>11}{'Wgeo':>9}{'X':>10}{'Delta':>11}"
           f"{'marg':>11}{'sign':>5}{'neff':>9}{'seX':>9}")

    geo2 = M.geometry(2, n, rng)
    print(f"\n[N=2 CONTROL]  exact quadrature says: (pi,2.5)->+0.13599 ; (0,2.5)->+0.72992", flush=True)
    print(hdr, flush=True)
    run(geo2, [(np.pi, -np.pi), (0.0, 0.0)], [2.5], a0=0.0)

    geo3 = M.geometry(3, n, rng)
    print(f"\n[N=3: THE TEST]  a in [{geo3['a'].min():.3f},{geo3['a'].max():.3f}]"
          f"   ({time.time()-t0:.0f}s)", flush=True)
    print(hdr, flush=True)
    ths3 = [(0.0, 0.0, 0.0), (0.3, 1.1, -1.4), (2*np.pi/3, 2*np.pi/3, -4*np.pi/3)]
    run(geo3, ths3, [2.5, 4.0, 5.5, 6.0, 6.5, 8.0, 12.0], a0=0.0)

    print(f"\n  [{time.time()-t0:.0f}s]  Clay 0%. No QED.", flush=True)

if __name__ == '__main__':
    main()
