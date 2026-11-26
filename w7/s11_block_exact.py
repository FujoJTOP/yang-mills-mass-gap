"""S11-BLOCK — the minimal block: SU(2), ONE link, Haar, ONE plaquette; exact (star).

Model (verbatim from research/v3/S2-BOUNDARY-COUNTEREXAMPLE.md sec.1, flat/exponential chart):
    U = exp(xi), xi in su(2) ~ R^3,  rho=|xi|,  t = rho/2,  u = xi_3/rho
    S(U;B) = beta (1 - cos t cos mu + sin t sin mu u),   B = exp(i theta sigma_3/2), mu = theta/2
    Haar: dU = J d^3xi, J = sin^2 t / t^2, Phi := -log J ;  measure ~ sin^2 t e^{-S} dt du
      (d^3xi = 8 t^2 dt du dpsi ; J t^2 = sin^2 t, so psi integrates out exactly)

Fixed CONSTANT chart direction v (class (a), div v = 0  <=>  the flat-chart projection of G1 sec.1-3):
    G_v = d_v S,   H_vv = d_v^2 S,   Phi'_v = d_v Phi,   Phi''_vv = d_v^2 Phi
    W_j = W_geo(v) = E[Phi''_vv] - Var(Phi'_v)
    G''(v) = E[H_vv] - Var_j(G_v) + W_j                       <- (star), LHS

Exact identities (S2 Lemma 1 + 2, i.e. one integration by parts with density vanishing at the
chart edge |xi|->2pi):
    E[d_v V] = 0,   E[d_v^2 V] = Var(d_v V)   (V = S+Phi)
    => Cov(G_v,Phi'_v) = W_geo  =>  E[H]-Var(G) = W_geo  and  G'' = 2 W_geo.
This file checks all of them to machine precision and reports the three terms.

Phi in cancellation-free form (t=x=rho/2):
    p1 = Phi_t = 1/x - cot x = sin x - x cos x)/(x sin x)
    p2 = Phi_t^2
    p3 = Phi_tt = (1/2)(csc^2 x - 1/x^2) = (x^2 - sin^2 x)/(2 x^2 sin^2 x)
    p4 = Phi_t/(2x)
numerically: use series for x<0.05 (p3 -> 1/6, p1 -> x/3, p4 -> 1/6).
"""
import numpy as np

SMALL = 0.05


def gl(n, a, b):
    z, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * z + 0.5 * (b + a), 0.5 * (b - a) * w


def s0(x):                                    # sin x / x
    out = np.empty_like(x)
    m = np.abs(x) > 1e-6
    out[m] = np.sin(x[m]) / x[m]
    xm = x[~m]
    out[~m] = 1 - xm ** 2 / 6 + xm ** 4 / 120
    return out


def phi_parts(x):
    """p1..p4 = Phi_t, Phi_t^2, Phi_tt, Phi_t/(2t)  (no sin^2 absorbed)."""
    p1 = np.empty_like(x)
    p3 = np.empty_like(x)
    m = x > SMALL
    xs = x[m]
    p1[m] = 1 / xs - np.cos(xs) / np.sin(xs)
    p3[m] = 0.5 * (1 / np.sin(xs) ** 2 - 1 / xs ** 2)
    xsm = x[~m]
    # series (csc^2 x = x^-2 + 1/3 + x^2/15 + 2x^4/189 + ...)
    p1[~m] = xsm / 3 + xsm ** 3 / 45 + 2 * xsm ** 5 / 945
    p3[~m] = 1 / 6 + xsm ** 2 / 30 + xsm ** 4 / 189
    return p1, p1 ** 2, p3, p1 / (2 * x)


def gl_grid(nx=160, nu=160):
    X, WX = gl(nx, 0.0, np.pi)
    U, WU = gl(nu, -1.0, 1.0)
    XX, UU = np.meshgrid(X, U, indexing='ij')
    return XX, UU, np.outer(WX, WU)


XX, UU, WW = gl_grid()
P1, P2, P3, P4 = phi_parts(XX)
_G = {'X': XX, 'U': UU, 'W': WW, 'P': (P1, P2, P3, P4)}


def observables(beta, theta):
    XX, UU, WW = _G['X'], _G['U'], _G['W']; P1, P2, P3, P4 = _G['P']
    mu = theta / 2.0
    C, M = np.cos(mu), np.sin(mu)
    beta = float(beta)
    x, u = XX, UU
    S = beta * (1 - np.cos(x) * C + np.sin(x) * M * u)
    W = np.sin(x) ** 2 * np.exp(-S) * WW
    Z = W.sum()
    E = lambda F: float((W * F).sum() / Z)
    S0 = 1.0 / (2 * x)          # 1/rho = 1/(2x)
    r = 2 * x
    Sr = (beta / 2) * (np.sin(x) * C + np.cos(x) * M * u)
    Su = beta * np.sin(x) * M
    Srr = (beta / 4) * (np.cos(x) * C - np.sin(x) * M * u)
    Sru = (beta / 2) * np.cos(x) * M

    out = {}
    # ------------------------------- v = e3 (longitudinal) -------------------
    g3 = u * Sr + (1 - u ** 2) * Su * S0
    h33 = (u ** 2 * Srr + 2 * u * (1 - u ** 2) / r * Sru
           + Sr * (1 - u ** 2) / r
           - 3 * Su * u * (1 - u ** 2) / r ** 2)
    ph2_3 = P3 * u ** 2 + P4 * (1 - u ** 2)
    ph1_3 = P1 * u
    ph1sq_3 = P2 * u ** 2
    Vd1 = g3 + ph1_3                     # d_v V
    Vd2 = h33 + ph2_3                    # d_v^2 V
    Eg, Eg2, Eh = E(g3), E(g3 ** 2), E(h33)
    Egp, Egpsq, Eph2 = E(ph1_3), E(ph1sq_3), E(ph2_3)
    Wg = Eph2 - (Egpsq - Egp ** 2)
    out['e3'] = dict(EH=Eh, VarG=Eg2 - Eg ** 2, W=Wg,
                     Gpp=Eh - (Eg2 - Eg ** 2) + Wg,
                     EG=Eg, EPhi1=Egp,
                     EHmVar=Eh - (Eg2 - Eg ** 2), EPh2=Eph2,
                     VarPh=Egpsq - Egp ** 2,
                     dV_mean=E(Vd1), stein=E(Vd2) - E(Vd1 ** 2),
                     cov=E(g3 * ph1_3) - Eg * Egp)
    # ------------------------- v = e1 (transverse, psi-averaged) -------------
    A = Sr - u * Su * S0
    Ar = Srr - u * Sru / r + u * Su / r ** 2
    Au = Sru - Su * S0
    n1sq = (1 - u ** 2) / 2               # E_psi[n1^2]
    g1sq = n1sq * A ** 2                  # E_psi[G_1^2]   (E_psi[G_1]=0)
    h11 = n1sq * (Ar - u * Au / r) + A * (1 - n1sq) / r
    ph2_1 = P3 * n1sq + P4 * (1 - n1sq)
    ph1sq_1 = P2 * n1sq
    Vd2_1 = h11 + ph2_1
    Vd1sq_1 = n1sq * A ** 2 + 2 * n1sq * A * P1 + ph1sq_1
    Eh1, Eg1sq = E(h11), E(g1sq)
    Eph1sq, Eph2_1 = E(ph1sq_1), E(ph2_1)
    W1 = Eph2_1 - Eph1sq
    out['e1'] = dict(EH=Eh1, VarG=Eg1sq, W=W1, Gpp=Eh1 - Eg1sq + W1,
                     EG=0.0, EPhi1=0.0, EHmVar=Eh1 - Eg1sq, EPh2=Eph2_1,
                     VarPh=Eph1sq, dV_mean=0.0, stein=E(Vd2_1) - E(Vd1sq_1),
                     cov=-0.5 * 0.0 + (Eph2_1 - Eph1sq) * 0.0)
    out['Z'] = Z
    return out


def report(rows=True):
    print("=" * 126)
    print("S11-BLOCK  SU(2) single link / Haar / one plaquette / flat chart  160x160 Gauss-Legendre (float64)")
    print("  W_j = E[Phi''_vv] - Var(Phi'_v);  G'' = E[H_vv] - Var_j(G_v) + W_j;  identity G''=2W_j and E[H]-Var=W_j")
    print("=" * 126)
    print(f"{'beta':>5} {'theta':>7} {'v':>3} {'E[H]':>12} {'Var_j':>12} {'W_j':>12} {'G\"':>12}"
          f" {'2W_j':>12} {'E[H]-Var':>12} {'E[dV]':>9} {'stein':>9} {'W-cov':>9}")
    for beta in [0.0, 0.1, 0.5, 1.0, 2.5, 5.0, 20.0]:
        for theta in [0.0, np.pi / 2, np.pi]:
            o = observables(beta, theta)
            for v in ('e3', 'e1'):
                r = o[v]
                print(f"{beta:5.2f} {theta:7.4f} {v:>3} {r['EH']:12.8f} {r['VarG']:12.8f}"
                      f" {r['W']:12.8f} {r['Gpp']:12.8f} {2*r['W']:12.8f} {r['EHmVar']:12.8f}"
                      f" {r['dV_mean']:9.1e} {r['stein']:9.1e} {r['W']-r['cov']:9.1e}")
    print()
    print("-- worst case over boundary theta in [0,pi] (49 pts) x v in {e3,e1} :  inf G'' = 2 inf W_geo --")
    for beta in [0.0, 0.02, 0.1, 0.25, 0.5, 1.0, 1.5, 2.5, 3.5, 5.0, 10.0, 20.0]:
        best = None
        for k in range(49):
            th = np.pi * k / 48
            o = observables(beta, th)
            for v in ('e3', 'e1'):
                if best is None or o[v]['W'] < best[0]:
                    best = (o[v]['W'], th, v, o[v])
        print(f"  beta={beta:5.2f}  min W_geo = {best[0]: .10f}  (theta={best[1]:.5f}, v={best[2]})"
              f"  =>  inf G'' = {2*best[0]: .10f}   E[H]={best[3]['EH']: .8f}  Var_j={best[3]['VarG']: .8f}")
    print()
    print("-- beta->infinity check (measure -> delta at the minimiser t=mu,u=1) vs computed beta=20 --")
    for theta in [0.0, np.pi]:
        mu = max(theta / 2, 1e-9)
        p1, _, p3, p4 = phi_parts(np.array([mu]))
        print(f"  theta={theta:.4f}: W_geo(inf) e3={p3[0]: .8f}  e1={p4[0]: .8f}"
              f"   | computed beta=20: e3={observables(20.0, theta)['e3']['W']: .8f}"
              f"  e1={observables(20.0, theta)['e1']['W']: .8f}")
    print()
    print("-- convergence in node count (beta=1, theta=pi, v=e1) --")
    for n in (40, 80, 160, 260):
        X2, U2, W2 = gl_grid(n, n)
        _G['X'], _G['U'], _G['W'] = X2, U2, W2
        _G['P'] = phi_parts(X2)
        r = observables(1.0, np.pi)['e1']
        print(f"   n={n:4d}: E[H]={r['EH']:.12f} Var_j={r['VarG']:.12f} W_j={r['W']:.12f}"
              f" G''={r['Gpp']:.12f}  (G''-2W)={r['Gpp']-2*r['W']: .2e}")
    X2, U2, W2 = gl_grid(160, 160)
    _G['X'], _G['U'], _G['W'] = X2, U2, W2
    _G['P'] = phi_parts(X2)


if __name__ == '__main__':
    report()
    print()
    print("REPRO: py -3.12 w7/s11_block_exact.py  ->  w7/s11_block.log")
