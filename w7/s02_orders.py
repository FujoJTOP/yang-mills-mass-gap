# -*- coding: utf-8 -*-
r"""
s02_orders.py -- S02: the u-orders of E[H_j], Var_j, W_j = W_geo + X (X=-2Cov), delta_W,
and whether (star)  G''_j = E[H_j] - Var_j + W_j >= delta_W > 0  is TRIVIALLY true.

Two independent witnesses:
 (A) single-plaquette SU(2), exact 1D reduction (q11_cert.py conventions), slice xi3=a0:
       w(u) = sin^2 u / u,  u = rho/2 in (a0/2, pi)
       G  = d_a S_W,  H = d_a^2 S_W,  Phi = 2 log(u/sin u),  Phi' = d_a Phi, Phi'' = d_a^2 Phi
     -> prints S=E[Phi''], Var(G), W_geo, X=-2Cov(G,Phi'), W_j=W_geo+X, G'' and the
        direct form E[H+Phi'']-Var(G+Phi').  This tests the exact cancellation OUTSIDE the
        Gaussian/flat-chart class (where G4 Lemma G2 proves X=-2W_geo exactly).
 (B) Gaussian flat-chart block (d=4, random SPD A,B): exact moments give
       E[H]=A, Var(G)=A.S.A, W_geo=B-B.S.B, Cov(G,Phi')=A.S.B
     -> checks the identity Delta_dyn = W_geo = Cov  =>  W_j = -W_geo, G''=0 to 1e-14.

Deterministic.  py -3.12.  Writes s02_orders.log only.
"""
import numpy as np
import mpmath as mp
import os

mp.mp.dps = 30
LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


def quad(f, lo, hi):
    return mp.quad(f, [lo, hi], maxdegree=10)


# ---------------------------------------------------------------- (A) single plaquette
def slice_moments(beta, theta, a0):
    c = mp.cos(theta / 2)
    s = mp.sin(theta / 2)
    lo = a0 / 2 if a0 > 0 else mp.mpf('1e-18')   # 0*(1/u) guards; weight ~ u kills the endpoint

    def SW(u):
        r = 2 * u
        return beta * (1 - mp.cos(u) * c + mp.sin(u) * s * a0 / r)

    def w(u):
        # FULL measure on the slice: Boltzmann x Haar Jacobian
        return mp.e ** (-SW(u)) * mp.sin(u) ** 2 / u

    Z = quad(w, lo, mp.pi)

    def rho(u):
        return 2 * u

    def Gf(u):
        a, r = a0, 2 * u
        return beta * (s * (mp.cos(u) * a ** 2 / (2 * r ** 2)
                            + mp.sin(u) * (1 / r - a ** 2 / r ** 3))
                       + c * mp.sin(u) * a / (2 * r))

    def Hf(u):
        a, r = a0, 2 * u
        up = a / (2 * r)
        A = a ** 2 / (2 * r ** 2)
        Ap = a / r ** 2 - a ** 3 / r ** 4
        B = 1 / r - a ** 2 / r ** 3
        Bp = 3 * a * (a ** 2 - r ** 2) / r ** 5
        C = a / (2 * r)
        Cp = (r ** 2 - a ** 2) / (2 * r ** 3)
        return beta * (s * (-mp.sin(u) * up * A + mp.cos(u) * Ap
                            + mp.cos(u) * up * B + mp.sin(u) * Bp)
                       + c * (mp.cos(u) * up * C + mp.sin(u) * Cp))

    def Phif(u):
        a, r = a0, 2 * u
        return (1 / u - mp.cot(u)) * a / (2 * u)

    def Phippf(u):
        a, r = a0, 2 * u
        Phr = 1 / u - mp.cot(u)
        Phrr = -1 / (2 * u ** 2) + 1 / (2 * mp.sin(u) ** 2)
        return Phrr * a ** 2 / (4 * u ** 2) + Phr * (1 / (2 * u) - a ** 2 / (8 * u ** 3))

    E = lambda f: quad(lambda u: f(u) * w(u), lo, mp.pi) / Z
    EG, EG2 = E(Gf), E(lambda u: Gf(u) ** 2)
    EH = E(Hf)
    Ep, Ep2 = E(Phif), E(lambda u: Phif(u) ** 2)
    Epp = E(Phippf)
    EGp = E(lambda u: Gf(u) * Phif(u))
    VarG = EG2 - EG ** 2
    Varp = Ep2 - Ep ** 2
    Cov = EGp - EG * Ep
    Wgeo = Epp - Varp
    X = -2 * Cov
    Wj = Wgeo + X
    Gpp = EH - VarG + Wj
    Gpp_direct = E(lambda u: Hf(u) + Phippf(u)) - (VarG + Varp + 2 * Cov)
    return dict(Z=Z, S=Epp, EH=EH, VarG=VarG, Wgeo=Wgeo, X=X, Wj=Wj,
                Cov=Cov, Gpp=Gpp, Gpp_direct=Gpp_direct, Varp=Varp)


def sec_A():
    P("=" * 100)
    P("[A] single-plaquette SU(2), slice xi3=a0, weight w=sin^2(u)/u, u in (a0/2,pi)")
    P("    G'' = E[H] - Var(G) + W_j,  W_j = W_geo + X,  W_geo=E[Phi'']-Var(Phi'),  X=-2Cov(G,Phi')")
    P("    'direct' = E[H+Phi'']-Var(G+Phi')  (exact identity, no W_j bookkeeping)")
    P("=" * 100)

    # --- anchors at beta=2.5, theta=pi, a0=0 (q11_cert / h1e_var)
    P("")
    P("[A.0] anchors (a0=0, theta=pi): must match q11_cert / N03")
    for b in (mp.mpf(0), mp.mpf('0.0208333333333333'), mp.mpf('2.5'), mp.mpf('6')):
        r = slice_moments(b, mp.pi, mp.mpf(0))
        P("  beta=%-8s  S=%.12f  Var(G)=%.12f  G''=%.12f" % (mp.nstr(b, 8), float(r['S']), float(r['VarG']), float(r['Gpp'])))
    r25 = slice_moments(mp.mpf('2.5'), mp.pi, mp.mpf(0))
    rr = float(r25['S'] / (r25['VarG'] / mp.mpf('2.5') ** 2))
    P("  beta* = sqrt(S/Var(g)) = %.10f   (anchor 4.30660785)" % (rr ** 0.5))
    P("  G''(2.5) = %.9f (anchor +0.135995) ;  G''(6) = %.9f (anchor -0.193019)"
      % (float(r25['Gpp']), float(slice_moments(mp.mpf(6), mp.pi, mp.mpf(0))['Gpp'])))

    P("")
    P("[A.1] slice table: a0 != 0.  Moments; FD in [A.2b] confirms E[H+Phi'']-Var(G+Phi') = -d^2logZ.")
    P("      Flat-chart (constant-direction) identity is checked in [B].")
    hdr = ("  %5s %6s %5s | %10s %10s %10s %10s %10s %10s %10s"
           % ("beta", "theta", "a0", "E[H]", "Var(G)", "W_geo", "Cov", "X=-2Cov", "W_j", "G''"))
    P(hdr)
    P("  " + "-" * (len(hdr) - 2))
    for b in (mp.mpf('0.02'), mp.mpf('2.5')):
        for th in (mp.mpf(0), mp.pi / 2, mp.pi):
            for a0 in (mp.mpf(0), mp.mpf('0.5'), mp.mpf(1.0)):
                r = slice_moments(b, th, a0)
                P("  %5s %6s %5s | %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f"
                  % (mp.nstr(b, 5), mp.nstr(th, 5), mp.nstr(a0, 4),
                     float(r['EH']), float(r['VarG']), float(r['Wgeo']), float(r['Cov']),
                     float(r['X']), float(r['Wj']), float(r['Gpp'])))

    P("")
    P("[A.2] exactness of the W_j bookkeeping:  |G'' - direct| (should be ~1e-25)")
    worst = mp.mpf(0)
    for b in (mp.mpf('0.02'), mp.mpf('2.5')):
        for th in (mp.mpf(0), mp.pi / 2, mp.pi):
            for a0 in (mp.mpf(0), mp.mpf('0.5'), mp.mpf(1.0)):
                r = slice_moments(b, th, a0)
                worst = max(worst, abs(r['Gpp'] - r['Gpp_direct']))
    P("  worst |G'' - [E[H+Phi'']-Var(G+Phi')]| = %s" % mp.nstr(worst, 6))
    P("")
    P("[A.2b] independent FD check: G'' =? -d^2/da^2 log Z(a), Z(a)=int_{a/2}^{pi} e^{-beta S_W} sin^2u/u du")
    P("  %5s %6s %5s | %14s %14s %10s" % ("beta", "theta", "a0", "-d2logZ", "momentcombo", "gap"))
    for b in (mp.mpf('0.02'), mp.mpf('2.5')):
        for th in (mp.mpf(0), mp.pi):
            for a0 in (mp.mpf('0.5'), mp.mpf('1.0')):
                r = slice_moments(b, th, a0)
                h = mp.mpf('1e-3')

                def logZ(a):
                    cc, ss = mp.cos(th / 2), mp.sin(th / 2)
                    f = lambda u: mp.e ** (-b * (1 - mp.cos(u) * cc + mp.sin(u) * ss * a / (2 * u))) * mp.sin(u) ** 2 / u
                    return mp.log(mp.quad(f, [a / 2, mp.pi], maxdegree=10))

                fpp = (-logZ(a0 + 2 * h) + 16 * logZ(a0 + h) - 30 * logZ(a0)
                       + 16 * logZ(a0 - h) - logZ(a0 - 2 * h)) / (12 * h ** 2)
                fd = -fpp   # G'' = -d^2/da^2 log Z(a)
                P("  %5s %6s %5s | %14.6f %14.6f %10.2e"
                  % (mp.nstr(b, 5), mp.nstr(th, 5), mp.nstr(a0, 4),
                     float(fd), float(r['Gpp']), float(abs(fd - r['Gpp']))))

    P("")
    P("[A.3] u-orders at strong coupling (u=beta/8, N=2): E[H] vs u, Var(G) vs u^2, W_geo vs u^0")
    P("  %6s %6s %12s %12s %12s %12s %10s" % ("beta", "u", "E[H]/u", "Var(G)/u^2", "W_geo", "X", "G''"))
    for b in ('0.002', '0.005', '0.01', '0.02'):
        b = mp.mpf(b)
        u = b / 8
        r = slice_moments(b, mp.pi, mp.mpf('0.5'))
        P("  %6.4f %6.5f %12.6f %12.6f %12.6f %12.6f %10.6f"
          % (float(b), float(u), float(r['EH'] / u), float(r['VarG'] / u ** 2),
             float(r['Wgeo']), float(r['X']), float(r['Gpp'])))


# ---------------------------------------------------------------- (B) Gaussian flat-chart block
def sec_B():
    P("")
    P("=" * 100)
    P("[B] Gaussian flat-chart block, d=4: S=1/2 x'A x, Jacobian Phi=1/2 x'B x, A,B SPD random")
    P("    E[H]=A, Var(G)=A.S.A, W_geo=B-B.S.B, Cov(G,Phi')=A.S.B,  S=(A+B)^-1")
    P("=" * 100)
    rng = np.random.default_rng(20261125)
    worst = [0.0, 0.0, 0.0, 0.0]
    for _ in range(5):
        M = rng.standard_normal((4, 4))
        A = M @ M.T + 4 * np.eye(4)
        M = rng.standard_normal((4, 4))
        B = M @ M.T + 4 * np.eye(4)
        S = np.linalg.inv(A + B)
        EH = A
        VarG = A @ S @ A
        Wgeo = B - B @ S @ B
        Cov = A @ S @ B
        Wj = Wgeo - 2 * Cov
        Gpp = EH - VarG + Wj
        worst[0] = max(worst[0], float(np.max(np.abs((EH - VarG) - Wgeo))))
        worst[1] = max(worst[1], float(np.max(np.abs(Cov - Wgeo))))
        worst[2] = max(worst[2], float(np.max(np.abs(Wj + Wgeo))))
        worst[3] = max(worst[3], float(np.max(np.abs(Gpp))))
    P("  max|(E[H]-Var(G)) - W_geo| = %.3e" % worst[0])
    P("  max|Cov(G,Phi') - W_geo|    = %.3e" % worst[1])
    P("  max|W_j + W_geo|            = %.3e   => W_j = -W_geo exactly" % worst[2])
    P("  max|G''|                    = %.3e   => literal (star) is IDENTICALLY 0" % worst[3])


if __name__ == "__main__":
    sec_A()
    sec_B()
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "s02_orders.log"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    print("\n[log] w7/s02_orders.log written")
