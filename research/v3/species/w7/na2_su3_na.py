#!/usr/bin/env python3
"""
na2_su3_na.py -- T0: (NA) on SU(N) one plaquette at the NATURAL conditioning point a0=0,
                 N=2 (control: reproduces na_plaquette.py) and N=3 (the real test).

WHY (ROUND1-VERDICTS 7.3 item 1 is SU(2)-only):
  Phi = -log J,  J(ksi) = prod_{a>0}( sin(a(ksi)/2)/(a(ksi)/2) )^2 , and the argument is the
  EIGENVALUE DIFFERENCE of A = ksi^a T^a (equivalently of the principal log), NOT the Cartan
  projection:
      A = V diag(Lam_1..Lam_N) V^dag ,  Lam_j = angles of eig(U),  sum_j Lam_j = 0 ,
      J = prod_{j<k}( sin((Lam_j-Lam_k)/2) / ((Lam_j-Lam_k)/2) )^2 .
  N=2 : ONE invariant (Lam_1-Lam_2 = rho = |ksi|) => Phi radial in ksi => d_v Phi odd in ksi_v
        => d_v Phi == 0 pointwise on {ksi_v=0} => X == 0 identically.   (the ROUND1 fact)
  N=3 : TWO invariants => Phi NOT radial => d_v Phi|_{a0=0} != 0 in general => X != 0.
  v = e_3, coarse field a = ksi_3 = tr(A lam^3).

CHART (exact, no domain convention needed): flat chart == principal log.
  dU = const * prod_{j<k} 4 sin^2((Lam_j-Lam_k)/2) dLam dV   (Weyl integration formula), so
  sample V ~ Haar(flag) and Lam with the Weyl density propto prod sin^2(D/2), reweight e^{-S}.
  N=3 sampler: D12,D23 iid propto sin^2(D/2) (exact inversion), residual weight sin^2(D13/2).

MODEL: S = beta(1 - (1/N) Re tr(U B^dag)),  B = diag(e^{i th_1},e^{i th_2},e^{i th_3}), sum th=0.

DECOMPOSITION (Lemma G1): G''(a0) = Ddyn + Wgeo + X ,
  G = d_a S, H = d_a^2 S, Phi = sum_{j<k} g(Lam_j-Lam_k), g(x)=2log(x/2)-2log sin(x/2),
  Ddyn = E[H]-Var(G),  Wgeo = E[Phi'']-Var(Phi'),  X = -2Cov(G,Phi').
d_a U analytic (Loewner divided differences of exp(i.)); d_a Lam by 1st/2nd order perturbation
theory.  H is NOT from a second divided difference (the Hadamard form built on the 2-point
second divided difference is WRONG for exp -- verified, 66% error); instead
     E[H] = -dE[G]/da + Var(G) + Cov(G,Phi')   (exact),   so  Ddyn = -dE[G]/da - X/2 .

ROUTES to X:
  M1 (primary, no differencing)   X = -2(E[G Phi'] - E[G]E[Phi'])
  M2 (exact identity, checked on model A)  X = +2 d/da E[Phi'|a] - 2 Wgeo
  M3 marginal  G'' = -d^2/da^2 log rho(a)   vs  Ddyn+Wgeo+X
Kernel smoothing is a real O(h^2) bias here (checked against model A); we report two
bandwidths and the h^2 Richardson extrapolation.

Run: py -3.12 -I na2_su3_na.py [--n 800000] [--tag run]
"""
import sys, json, time
import numpy as np

# ------------------------------------------------------------------ Lie algebra
def lambdas(N):
    if N == 2:
        return [np.array([[0, 1], [1, 0]], dtype=complex),
                np.array([[0, -1j], [1j, 0]], dtype=complex),
                np.array([[1, 0], [0, -1]], dtype=complex)]
    if N == 3:
        def e(j, k, v):
            M = np.zeros((3, 3), dtype=complex); M[j, k] = v; return M
        return [e(0, 1, 1) + e(1, 0, 1), e(0, 1, -1j) + e(1, 0, 1j),
                np.diag([1, -1, 0]).astype(complex),
                e(0, 2, 1) + e(2, 0, 1), e(0, 2, -1j) + e(2, 0, 1j),
                e(1, 2, 1) + e(2, 1, 1), e(1, 2, -1j) + e(2, 1, 1j),
                np.diag([1, 1, -2]).astype(complex) / np.sqrt(3)]
    raise ValueError(N)

# ------------------------------------------------------------------ samplers
def haar_su(N, M, rng):
    Z = (rng.standard_normal((M, N, N)) + 1j * rng.standard_normal((M, N, N))) / np.sqrt(2.0)
    Q, R = np.linalg.qr(Z)
    d = np.diagonal(R, axis1=1, axis2=2)
    Q = Q * (d / np.abs(d))[:, None, :].conj()
    det = np.linalg.det(Q)
    return Q * np.exp(-1j * np.angle(det) / N)[:, None, None]

def sample_sin2(M, rng):
    """exact sample of (1/pi) sin^2(x/2) on [0,2pi]: Newton inversion of (X-sinX)/(2pi)."""
    u = rng.random(M)
    X = 2 * np.pi * u
    for _ in range(200):
        f = (X - np.sin(X)) / (2 * np.pi) - u
        fp = (1 - np.cos(X)) / (2 * np.pi)
        dX = np.clip(f / np.maximum(fp, 1e-14), -0.5, 0.5)
        X = X - dX
        if np.max(np.abs(dX)) < 1e-15:
            break
    return np.clip(X, 0.0, 2 * np.pi)

def sample_Lam(N, M, rng):
    if N == 2:
        rho = sample_sin2(M, rng)
        return np.stack([rho / 2, -rho / 2], axis=1), np.ones(M)
    d12 = sample_sin2(M, rng); d23 = sample_sin2(M, rng)
    good = (d12 + d23) < 2 * np.pi
    Lam = np.empty((M, 3))
    Lam[:, 0] = (2 * d12 + d23) / 3
    Lam[:, 1] = (-d12 + d23) / 3
    Lam[:, 2] = -(d12 + 2 * d23) / 3
    return Lam, np.where(good, np.sin((d12 + d23) / 2) ** 2, 0.0)

# ------------------------------------------------------------------ geometry
def geometry(N, M, rng, chunk=250000):
    lam3 = lambdas(N)[2]
    acc = {k: [] for k in ('a', 'wg', 'Phi', 'P1', 'P2')}
    Vout, Lout = [], []
    got = 0
    while got < M:
        m = min(chunk, M - got)
        V = haar_su(N, m, rng)
        Lam, wg = sample_Lam(N, m, rng)
        keep = wg > 1e-13
        V, Lam, wg = V[keep], Lam[keep], wg[keep]
        if Lam.shape[0] == 0:
            continue
        m = Lam.shape[0]
        Mt = np.swapaxes(V.conj(), 1, 2) @ (lam3 / 2.0) @ V
        d = Lam[:, :, None] - Lam[:, None, :]
        ix = np.arange(N)
        d[:, ix, ix] = np.inf
        m1 = np.real(np.diagonal(Mt, axis1=1, axis2=2))
        s2 = 2.0 * np.real(np.sum(np.abs(Mt) ** 2 / d, axis=2))
        iu = np.triu_indices(N, 1)
        Dl = Lam[:, iu[0]] - Lam[:, iu[1]]
        g = 2 * np.log(np.abs(Dl) / 2) - 2 * np.log(np.abs(np.sin(Dl / 2)))
        gp = 2.0 / Dl - 1.0 / np.tan(Dl / 2)
        gpp = -2.0 / Dl ** 2 + 0.5 / np.sin(Dl / 2) ** 2
        dPhi = np.zeros((m, N)); d2Phi = np.zeros((m, N, N))
        for idx, (j, k) in enumerate(zip(iu[0], iu[1])):
            dPhi[:, j] += gp[:, idx]; dPhi[:, k] -= gp[:, idx]
            d2Phi[:, j, j] += gpp[:, idx]; d2Phi[:, k, k] += gpp[:, idx]
            d2Phi[:, j, k] -= gpp[:, idx]; d2Phi[:, k, j] -= gpp[:, idx]
        acc['P1'].append(np.einsum('mj,mj->m', dPhi, m1))
        acc['P2'].append(np.einsum('mj,mj->m', dPhi, s2)
                         + np.einsum('mj,mk,mjk->m', m1, m1, d2Phi))
        P = np.abs(V) ** 2
        acc['a'].append(np.sum(Lam * (P[:, 0, :] - P[:, 1, :]), axis=1))
        acc['wg'].append(wg); acc['Phi'].append(g.sum(1))
        Vout.append(V); Lout.append(Lam)
        got += m
    geo = {k: np.concatenate(v) for k, v in acc.items()}
    geo['V'] = np.concatenate(Vout); geo['Lam'] = np.concatenate(Lout); geo['N'] = N
    return geo

def action_traces(geo, Bd_list, chunk=200000):
    """Re tr(U B^dag) and d_a Re tr(U B^dag); Loewner divided differences of exp(i.)."""
    N = geo['N']; M = geo['Lam'].shape[0]
    lam3 = lambdas(N)[2]
    acc = [[] for _ in Bd_list]
    for s in range(0, M, chunk):
        sl = slice(s, min(s + chunk, M))
        V = geo['V'][sl]; Lam = geo['Lam'][sl]
        Vd = np.swapaxes(V.conj(), 1, 2)
        e = np.exp(1j * Lam)
        U = (V * e[:, None, :]) @ Vd
        Mt = Vd @ (lam3 / 2.0) @ V
        d = Lam[:, :, None] - Lam[:, None, :]
        same = np.abs(d) < 1e-9
        dd = np.where(same, 1.0, d)
        D1 = 2j * np.exp(0.5j * (Lam[:, :, None] + Lam[:, None, :])) * (np.sin(dd / 2) / dd)
        D1 = np.where(same, (1j * e)[:, :, None], D1)
        dU = (V @ (Mt * D1)) @ Vd
        for k, Bd in enumerate(Bd_list):
            Bd = np.asarray(Bd, dtype=complex)
            acc[k].append(tuple(np.real(np.einsum('mij,ji->m', X, Bd)) for X in (U, dU)))
    return [tuple(np.concatenate([a[i] for a in acc[k]]) for i in range(2))
            for k in range(len(Bd_list))]

# ------------------------------------------------------------------ conditional stats
def cond_stats(geo, ct, beta, a0, h):
    N = geo['N']
    S = beta * (1.0 - ct[0] / N)
    G = -beta * ct[1] / N
    a = geo['a']
    K = np.exp(-0.5 * ((a - a0) / h) ** 2)
    w = geo['wg'] * np.exp(-S) * K
    Z = w.sum()
    if not np.isfinite(Z) or Z <= 0:
        return None
    wn = w / Z
    E = lambda x: float(np.sum(wn * x))
    def se(x, mu):
        return float(np.sqrt(max(np.sum(wn ** 2 * (x - mu) ** 2), 0.0)))
    P1 = geo['P1']
    EG, EP = E(G), E(P1)
    VarG = E(G * G) - EG ** 2
    VarP = E(P1 * P1) - EP ** 2
    Cov = E(G * P1) - EG * EP
    o = dict(neff=1.0 / np.sum(wn ** 2), rho=float(Z), EG=EG, EP=EP,
             X=-2 * Cov, seX=2 * se((G - EG) * (P1 - EP), Cov),
             Wgeo=E(geo['P2']) - VarP)
    s = (a - a0) / h ** 2
    dEG = float(np.sum(w * s * (G - EG))) / Z
    dEP = float(np.sum(w * s * (P1 - EP))) / Z
    o['X_m2'] = 2 * dEP - 2 * o['Wgeo']
    o['Ddyn'] = -dEG - 0.5 * o['X']
    o['Gpp'] = o['Ddyn'] + o['Wgeo'] + o['X']
    s2 = ((a - a0) ** 2 - h ** 2) / h ** 4
    d1 = float(np.sum(w * s)) / Z; d2 = float(np.sum(w * s2)) / Z
    o['Gpp_marg'] = -(d2 / o['rho'] - (d1 / o['rho']) ** 2)
    return o

def extrap(geo, ct, beta, a0, hs=(0.12, 0.06)):
    """h^2 Richardson extrapolation of the kernel estimator (bias O(h^2))."""
    r1 = cond_stats(geo, ct, beta, a0, hs[0])
    r2 = cond_stats(geo, ct, beta, a0, hs[1])
    if r1 is None or r2 is None:
        return None
    out = {}
    for k in ('X', 'Wgeo', 'Ddyn', 'Gpp', 'X_m2', 'Gpp_marg'):
        out[k] = r2[k] + (r2[k] - r1[k]) / 3.0
        out[k + '_hi'] = r1[k]; out[k + '_lo'] = r2[k]
    out['seX'] = max(r1['seX'], r2['seX'])
    out['neff'] = r1['neff']
    out['rho'] = r1['rho']
    return out

def header():
    print("=" * 118, flush=True)
    print("T0 : (NA) on SU(N) one plaquette, NATURAL point a0=0.   X = -2Cov(G,Phi'|a0);  v = e3", flush=True)
    print("     flat chart = principal log (chart-exact).  G''=Ddyn+Wgeo+X.  kernel h=0.12/0.06, h^2 extrap.",
          flush=True)
    print("=" * 118, flush=True)

def main():
    n, tag = 800000, 'run'
    for i, arg in enumerate(sys.argv):
        if arg == '--n':
            n = int(sys.argv[i + 1])
        if arg == '--tag':
            tag = sys.argv[i + 1]
    t0 = time.time(); header()
    rng = np.random.default_rng(20261124)
    res = {}

    # ---------------- N=2 control
    print("\n[N=2] CONTROL vs na_plaquette.log (exact quadrature, theta convention of model A).", flush=True)
    geo2 = geometry(2, n, rng)
    print(f"     a in [{geo2['a'].min():.3f},{geo2['a'].max():.3f}]   ({time.time()-t0:.0f}s)", flush=True)
    ths = [0.0, np.pi / 2, np.pi]
    cts2 = action_traces(geo2, [np.diag(np.exp(+1j * np.array([t / 2, -t / 2]))) for t in ths])
    ex = {(0.0, 2.5, 2.0): 0.0118, (np.pi, 2.5, 2.0): 0.0612,
          (np.pi, 2.5, 1.0): 0.0316458126, (np.pi, 2.5, -1.0): -0.01617,
          (np.pi / 2, 2.5, -1.0): -0.00179, (0.0, 2.5, 1.0): 0.00222,
          (np.pi, 6.0, 1.0): 0.10750, (0.0, 2.5, 0.0): 0.0, (np.pi, 2.5, 0.0): 0.0}
    print(f"     {'theta':>7}{'beta':>6}{'a0':>6}{'X_extrap':>11}{'se':>9}{'exact':>11}{'ratio':>8}",
          flush=True)
    for t, ct in zip(ths, cts2):
        for a0 in (0.0, 1.0, -1.0, 2.0):
            r = extrap(geo2, ct, 2.5, a0)
            e = ex.get((round(t, 6), 2.5, a0), float('nan'))
            print(f"     {t:7.3f}{2.5:6.2f}{a0:6.2f}{r['X']:+11.6f}{r['seX']:9.6f}{e:11.6f}"
                  f"{r['X']/e if e not in (0.0,) and np.isfinite(e) else float('nan'):8.3f}", flush=True)
    r6 = extrap(geo2, cts2[2], 6.0, 1.0)
    print(f"     [theta=pi, beta=6, a0=+1]  X = {r6['X']:+.5f}  (exact 0.10750)", flush=True)

    # ---------------- N=3 : THE TEST
    print(f"\n[N=3] *** THE TEST: X(a0=0) for SU(3).  ({time.time()-t0:.0f}s)", flush=True)
    geo3 = geometry(3, n, rng)
    print(f"     a in [{geo3['a'].min():.3f},{geo3['a'].max():.3f}]   ({time.time()-t0:.0f}s)", flush=True)
    ths3 = [(0.0, 0.0, 0.0), (0.3, 1.1, -1.4), (np.pi / 2, np.pi / 2, -np.pi),
            (np.pi, np.pi, -2 * np.pi), (2 * np.pi / 3, 2 * np.pi / 3, -4 * np.pi / 3)]
    cts3 = action_traces(geo3, [np.diag(np.exp(-1j * np.array(t))) for t in ths3])
    print(f"     {'(th1,th2,th3)':>24}{'beta':>6}{'a0':>5}{'X':>11}{'se':>9}{'X_M2':>10}"
          f"{'Ddyn':>10}{'Wgeo':>9}{'Gpp':>10}{'Gpp_marg':>10}{'neff':>8}", flush=True)
    rows = []
    for lab, ct in zip([str(np.round(t, 3)) for t in ths3], cts3):
        for beta in (0.5, 1.0, 2.5, 4.0):
            for a0 in (0.0, 0.5, -0.5, 1.0, -1.0):
                r = extrap(geo3, ct, beta, a0)
                print(f"     {lab:>24}{beta:6.2f}{a0:5.2f}{r['X']:+11.6f}{r['seX']:9.6f}"
                      f"{r['X_m2']:+10.5f}{r['Ddyn']:+10.5f}{r['Wgeo']:+9.5f}{r['Gpp']:+10.5f}"
                      f"{r['Gpp_marg']:+10.5f}{r['neff']:8.0f}", flush=True)
                rows.append(dict(lab=lab, beta=beta, a0=a0, **r))
    with open(f'na2_su3_na_{tag}.json', 'w') as f:
        json.dump(rows, f, indent=1, default=float)

    # ---------------- bandwidth convergence of the a0=0 statement
    print("\n[H-SCAN] is X(a0=0)!=0 real, or a residual kernel bias?  beta=2.5, a0=0 and a0=1", flush=True)
    print(f"     {'boundary':>24}{'a0':>5}" + "".join(f"{'h='+str(h):>12}" for h in (0.24, 0.12, 0.06, 0.03))
          + f"{'h->0':>12}{'se':>9}", flush=True)
    for k in (0, 3, 4, 2):
        ct = cts3[k]
        for a0 in (0.0, 1.0):
            vals = []
            for h in (0.24, 0.12, 0.06, 0.03):
                r = cond_stats(geo3, ct, 2.5, a0, h)
                vals.append(r['X'])
            r06 = cond_stats(geo3, ct, 2.5, a0, 0.06)
            r03 = cond_stats(geo3, ct, 2.5, a0, 0.03)
            ext = r03['X'] + (r03['X'] - r06['X']) / 3.0
            print(f"     {str(np.round(ths3[k],3)):>24}{a0:5.1f}"
                  + "".join(f"{v:+12.3e}" for v in vals) + f"{ext:+12.3e}{r03['seX']:9.2e}", flush=True)
    print(f"\n  [{time.time()-t0:.0f}s]  Clay 0%. No QED. Numbers only.", flush=True)

if __name__ == '__main__':
    main()
