#!/usr/bin/env python3
"""
na2_chain.py -- T1(b): (NA) on a 1D chain of L SU(2) links, L = 2,3,4 (open and periodic).
   State space SU(2)^L, flat chart per link ( = principal log ), measure prod_i J(xi_i) d^3xi_i .
   Bond action  S = beta sum_{<i,j>} (1 - 1/2 Re tr(U_i U_j^dag))   (nearest neighbours).
   Coarse fields:
      (A) chain average   a = (1/L) sum_i xi_3^{(i)}          <- the ONLY way to have Phi' != 0 at a=0
      (B) endpoint        a = xi_3^{(1)}                       <- control: Prop. D => X(a0=0) == 0
   Derivatives:  d_a = (1/L) sum_i d/dxi_3^{(i)}  (so Phi' = (1/L) sum_i Phi'_i , Phi'' = (1/L^2) sum_i Phi''_i;
   cross terms vanish because Phi is a sum of one-link functions).

WHY (A) IS THE POINT: Phi = sum_i phi(|xi_i|) is radial in EACH link, so d/dxi_3^{(i)} phi is odd in
xi_3^{(i)} and vanishes only when xi_3^{(i)}=0 individually.  The constraint sum_i xi_3^{(i)} = L a does
NOT force each xi_3^{(i)} to vanish at a = 0 => Phi' != 0 on the conditioning surface.  This is the
non-degenerate conditioning the one-plaquette model cannot produce with a single-link coarse field.

Run: py -3.12 -I na2_chain.py [--n 600000]
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import na2_su3_na as M

def link_matrices(geo):
    """U_i and dU_i (d/d xi_3 at fixed other flat-chart coordinates) -- Loewner divided differences."""
    V, Lam, N = geo['V'], geo['Lam'], geo['N']
    lam3 = M.lambdas(N)[2]
    Vd = np.swapaxes(V.conj(), 1, 2)
    e = np.exp(1j * Lam)
    U = (V * e[:, None, :]) @ Vd
    Mt = Vd @ (lam3 / 2.0) @ V
    d = Lam[:, :, None] - Lam[:, None, :]
    same = np.abs(d) < 1e-9
    dd = np.where(same, 1.0, d)
    D1 = 2j * np.exp(0.5j * (Lam[:, :, None] + Lam[:, None, :])) * (np.sin(dd / 2) / dd)
    D1 = np.where(same, (1j * e)[:, :, None], D1)
    return U, (V @ (Mt * D1)) @ Vd

def build(L, m, rng, periodic):
    """L links; returns per-link U, dU, a_i, P1_i, P2_i and the bond list."""
    Us, dUs, ais, P1s, P2s = [], [], [], [], []
    for _ in range(L):
        g = M.geometry(2, m, rng)
        U, dU = link_matrices(g)
        Us.append(U); dUs.append(dU); ais.append(g['a'])
        P1s.append(g['P1']); P2s.append(g['P2'])
    bonds = [(i, i + 1) for i in range(L - 1)] + ([(L - 1, 0)] if periodic and L > 2 else [])
    return dict(U=Us, dU=dUs, a=ais, P1=P1s, P2=P2s, bonds=bonds, L=L, m=m)

def stats(ch, beta, a0, h, mode='avg'):
    L, bonds = ch['L'], ch['bonds']
    if mode == 'avg':
        a = sum(ch['a']) / L
        P1 = sum(ch['P1']) / L
        P2 = sum(ch['P2']) / (L * L)
    else:
        a, P1, P2 = ch['a'][0], ch['P1'][0], ch['P2'][0]
    # action and its a-derivative
    c0 = 0.0
    dS = np.zeros(ch['m'])
    for (i, j) in bonds:
        Wj = np.swapaxes(ch['U'][j].conj(), 1, 2)          # U_j^dag
        Wi = np.swapaxes(ch['U'][i].conj(), 1, 2)
        c0 = c0 + np.real(np.einsum('mab,mba->m', ch['U'][i], Wj))
        t_i = np.real(np.einsum('mab,mba->m', ch['dU'][i], Wj))   # d_i Re tr(U_i U_j^dag)
        t_j = np.real(np.einsum('mab,mba->m', ch['dU'][j], Wi))   # d_j Re tr(U_i U_j^dag)
        dS = dS - 0.5 * beta * (t_i + t_j)
    S = beta * (len(bonds) - 0.5 * c0)
    G = dS / L if mode == 'avg' else np.zeros(ch['m'])
    if mode != 'avg':
        # endpoint coarse field: d_a = d/d xi_3^{(1)} only
        G = np.zeros(ch['m'])
        for (i, j) in bonds:
            Wj = np.swapaxes(ch['U'][j].conj(), 1, 2); Wi = np.swapaxes(ch['U'][i].conj(), 1, 2)
            if i == 0:
                G = G - 0.5 * beta * np.real(np.einsum('mab,mba->m', ch['dU'][0], Wj))
            if j == 0:
                G = G - 0.5 * beta * np.real(np.einsum('mab,mba->m', ch['dU'][0], Wi))
    w = np.exp(-S) * np.exp(-0.5 * ((a - a0) / h) ** 2)
    Z = w.sum()
    if not np.isfinite(Z) or Z <= 0:
        return None
    wn = w / Z
    E = lambda x: float(np.sum(wn * x))
    EG, EP = E(G), E(P1)
    VarG = E(G * G) - EG ** 2; VarP = E(P1 * P1) - EP ** 2
    Cov = E(G * P1) - EG * EP
    seC = float(np.sqrt(max(np.sum(wn ** 2 * ((G - EG) * (P1 - EP) - Cov) ** 2), 0.0)))
    s = (a - a0) / h ** 2
    dEG = float(np.sum(w * s * (G - EG))) / Z
    dEP = float(np.sum(w * s * (P1 - EP))) / Z
    W = E(P2) - VarP; X = -2 * Cov
    return dict(X=X, seX=2 * seC, Wgeo=W, Ddyn=-dEG - 0.5 * X, Gpp=-dEG + 0.5 * X + W,
                X_m2=2 * dEP - 2 * W, neff=1.0 / np.sum(wn ** 2), a2=float(E(a * a)))

def extrap(ch, beta, a0, mode='avg', hs=(0.10, 0.05)):
    r1 = stats(ch, beta, a0, hs[0], mode); r2 = stats(ch, beta, a0, hs[1], mode)
    if r1 is None or r2 is None:
        return None
    o = {k: r2[k] + (r2[k] - r1[k]) / 3.0 for k in ('X', 'Wgeo', 'Ddyn', 'Gpp', 'X_m2')}
    o['seX'] = max(r1['seX'], r2['seX']); o['neff'] = r1['neff']; o['a2'] = r1['a2']
    return o

def main():
    n = 600000
    for i, a in enumerate(sys.argv):
        if a == '--a' or a == '--n':
            if a == '--n':
                n = int(sys.argv[i + 1])
    t0 = time.time()
    print("=" * 112, flush=True)
    print("T1(b): 1D chain, L SU(2) links, flat chart per link, bond action S=beta sum(1-1/2 Re tr U_i U_j^dag)", flush=True)
    print("       coarse field (A) chain average a=(1/L) sum xi_3^i   (B) endpoint xi_3^1 (control, Prop. D)", flush=True)
    print("=" * 112, flush=True)
    rng = np.random.default_rng(99)
    rows = []
    for periodic in (False, True):
        for L in (2, 3, 4):
            if periodic and L == 2:
                continue
            ch = build(L, n, rng, periodic)
            tag = f"L={L}{' ring' if periodic else ' open'}"
            print(f"\n[{tag}]  ({(time.time()-t0):.0f}s)", flush=True)
            for mode in ('avg', 'end'):
                print(f"   coarse field = {'chain average' if mode=='avg' else 'endpoint xi_3^1'}"
                      f"   <a^2>={extrap(ch,1.0,0.0,mode)['a2']:.3f}", flush=True)
                print(f"     {'beta':>5}{'a0':>5}{'X':>11}{'se':>9}{'X_M2':>10}{'Ddyn':>10}"
                      f"{'Wgeo':>9}{'Gpp':>10}{'neff':>8}", flush=True)
                for beta in (0.5, 1.0, 2.5):
                    for a0 in (0.0, 1.0, -1.0):
                        r = extrap(ch, beta, a0, mode)
                        print(f"     {beta:5.2f}{a0:5.2f}{r['X']:+11.6f}{r['seX']:9.6f}{r['X_m2']:+10.5f}"
                              f"{r['Ddyn']:+10.5f}{r['Wgeo']:+9.5f}{r['Gpp']:+10.5f}{r['neff']:8.0f}",
                              flush=True)
                        rows.append(dict(tag=tag, mode=mode, beta=beta, a0=a0, **r))
    with open('na2_chain.json', 'w') as f:
        json.dump(rows, f, indent=1, default=float)
    print(f"\n  [{time.time()-t0:.0f}s]  Clay 0%. No QED.", flush=True)

if __name__ == '__main__':
    main()
