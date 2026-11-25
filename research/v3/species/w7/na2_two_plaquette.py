#!/usr/bin/env python3
"""
na2_two_plaquette.py -- T1(a): (NA) on TWO plaquettes sharing ONE link.
   Model: shared link U = exp(i xi^a sigma^a/2) (flat chart, Haar measure), coarse field
          a = xi_3, v = e_3.  Two plaquettes  U B1^dag  and  U B2^dag  (the other three links of
          each plaquette are frozen at the staples B1, B2 -- a Dirichlet/mean-field 2-plaquette):
              S = beta [ (1 - 1/2 Re tr(U B1^dag)) + (1 - 1/2 Re tr(U B2^dag)) ] .
          The Jacobian is still J(xi) (radial, rank 1), but G = d_a S and Phi' no longer have a
          parity that forces Cov = 0 -- except at a0 = 0, where Proposition D still applies.
   Reports (model, conditioning, v, a0) + (Ddyn, Wgeo, X, G'') for a grid of staple pairs.

Run: py -3.12 -I na2_two_plaquette.py [--n 800000]
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import na2_su3_na as M

def staple(ax, t):
    """B = exp(i t sigma_ax / 2) for ax in {0,1,2} (Pauli index), complex 2x2."""
    s = M.lambdas(2)
    return np.eye(2, dtype=complex) * np.cos(t / 2) + 1j * np.sin(t / 2) * s[ax]

def model(geo, Bd1, Bd2, beta, a0, h):
    """S = beta[(1-c01/2)+(1-c02/2)];  combine the two action_traces."""
    (c0a, c1a), (c0b, c1b) = M.action_traces(geo, [Bd1, Bd2])
    N = geo['N']
    S = beta * (2.0 - 0.5 * (c0a + c0b))
    G = -0.5 * beta * (c1a + c1b)
    a = geo['a']; P1 = geo['P1']; P2 = geo['P2']
    w = np.exp(-S) * np.exp(-0.5 * ((a - a0) / h) ** 2)
    Z = w.sum()
    wn = w / Z
    E = lambda x: float(np.sum(wn * x))
    EG, EP = E(G), E(P1)
    VarG = E(G * G) - EG ** 2
    VarP = E(P1 * P1) - EP ** 2
    Cov = E(G * P1) - EG * EP
    seC = float(np.sqrt(max(np.sum(wn ** 2 * ((G - EG) * (P1 - EP) - Cov) ** 2), 0.0)))
    s = (a - a0) / h ** 2
    dEG = float(np.sum(w * s * (G - EG))) / Z
    dEP = float(np.sum(w * s * (P1 - EP))) / Z
    W = E(P2) - VarP
    X = -2 * Cov
    return dict(X=X, seX=2 * seC, Wgeo=W, Ddyn=-dEG - 0.5 * X,
                Gpp=-dEG - 0.5 * X + W + X, X_m2=2 * dEP - 2 * W, neff=1.0 / np.sum(wn ** 2))

def extrap(geo, Bd1, Bd2, beta, a0, hs=(0.12, 0.06)):
    r1 = model(geo, Bd1, Bd2, beta, a0, hs[0])
    r2 = model(geo, Bd1, Bd2, beta, a0, hs[1])
    o = {k: r2[k] + (r2[k] - r1[k]) / 3.0 for k in ('X', 'Wgeo', 'Ddyn', 'Gpp', 'X_m2')}
    o['seX'] = max(r1['seX'], r2['seX'])
    o['neff'] = r1['neff']
    return o

def main():
    n = 800000
    for i, a in enumerate(sys.argv):
        if a == '--n':
            n = int(sys.argv[i + 1])
    t0 = time.time()
    print("=" * 112, flush=True)
    print("T1(a): TWO plaquettes sharing ONE link (frozen staples B1,B2).  SU(2) flat chart, a=xi_3, v=e_3.", flush=True)
    print("       S = beta[(1-(1/2)Re tr(UB1^dag)) + (1-(1/2)Re tr(UB2^dag))];  G''=Ddyn+Wgeo+X.", flush=True)
    print("       Proposition D predicts X(a0=0) == 0 for ANY S (rank SU(2)=1): checked below.", flush=True)
    print("=" * 112, flush=True)
    rng = np.random.default_rng(4242)
    geo = M.geometry(2, n, rng)
    print(f"     samples={n}   a in [{geo['a'].min():.2f},{geo['a'].max():.2f}]   ({time.time()-t0:.0f}s)",
          flush=True)

    cases = [('B1=B2=1', staple(2, 0.0), staple(2, 0.0)),
             ('1 & e3(pi)', staple(2, 0.0), staple(2, np.pi)),
             ('1 & e3(pi/2)', staple(2, 0.0), staple(2, np.pi / 2)),
             ('e3(pi) & e3(pi/2)', staple(2, np.pi), staple(2, np.pi / 2)),
             ('e3(pi) & e3(-pi/2)', staple(2, np.pi), staple(2, -np.pi / 2)),
             ('e3(pi) & e1(pi)', staple(2, np.pi), staple(1, np.pi)),
             ('e3(pi) & e1(pi/2)', staple(2, np.pi), staple(1, np.pi / 2)),
             ('e1(pi/2) & e2(pi/2)', staple(1, np.pi / 2), staple(2, np.pi / 2)),
             ('1 & e1(pi/2)', staple(2, 0.0), staple(1, np.pi / 2))]
    print(f"\n     {'staples':>22}{'beta':>6}{'a0':>5}{'X':>11}{'se':>9}{'X_M2':>10}"
          f"{'Ddyn':>10}{'Wgeo':>9}{'Gpp':>10}{'neff':>8}", flush=True)
    rows = []
    for lab, B1, B2 in cases:
        for beta in (1.0, 2.5, 4.0):
            for a0 in (0.0, 0.5, -0.5, 1.0, -1.0):
                r = extrap(geo, B1, B2, beta, a0)
                print(f"     {lab:>22}{beta:6.2f}{a0:5.2f}{r['X']:+11.6f}{r['seX']:9.6f}"
                      f"{r['X_m2']:+10.5f}{r['Ddyn']:+10.5f}{r['Wgeo']:+9.5f}{r['Gpp']:+10.5f}"
                      f"{r['neff']:8.0f}", flush=True)
                rows.append(dict(lab=lab, beta=beta, a0=a0, **r))
    with open('na2_two_plaquette.json', 'w') as f:
        json.dump(rows, f, indent=1, default=float)

    print("\n     [SIGN SUMMARY]  X at a0=0 (all cases must be 0 by Prop. D) and at a0=+-1:", flush=True)
    print(f"     {'staples':>22}{'X(0),beta=2.5':>16}{'X(+1),b=1':>12}{'X(+1),b=2.5':>13}"
          f"{'X(-1),b=2.5':>13}", flush=True)
    for lab, B1, B2 in cases:
        a0_0 = extrap(geo, B1, B2, 2.5, 0.0)
        a1_1 = extrap(geo, B1, B2, 1.0, 1.0)
        a1_25 = extrap(geo, B1, B2, 2.5, 1.0)
        am1_25 = extrap(geo, B1, B2, 2.5, -1.0)
        print(f"     {lab:>22}{a0_0['X']:+16.3e}{a1_1['X']:+12.5f}{a1_25['X']:+13.5f}"
              f"{am1_25['X']:+13.5f}", flush=True)
    print(f"\n  [{time.time()-t0:.0f}s]  Clay 0%. No QED.", flush=True)

if __name__ == '__main__':
    main()
