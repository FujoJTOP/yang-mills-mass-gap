#!/usr/bin/env python3
"""
ds_ratio.py -- INDEPENDENT computation of the split-quotient  D/S  for G'' = S - D.

Purpose: decide whether "D/S ~ 2%" can be a coordinate.  Everything here is
recomputed from scratch (own jet-free analytic derivatives), not read from the
corpus's implementations.

MODEL A (one plaquette + fixed boundary holonomy), the corpus's exactly solvable cell
------------------------------------------------------------------------------
  link  U = exp(i xi.sigma/2),  xi in R^3, rho=|xi|,  boundary  B = exp(i th sigma_3/2),
  plaquette  U_p = U B ,   (1/2)Re tr U_p = cos(rho/2)cos(th/2) - sin(rho/2)sin(th/2) u,
  u = xi_3/rho,   S_W = beta (1 - (1/2)Re tr U_p) ,
  flat chart, Haar density J ~ sin^2(rho/2)/rho^2 ,  Phi = - log J ,
  measure  exp(-S_W - Phi) d^3 xi  on |xi| < 2 pi.

  Splits of the exact marginal identity   G''(a) = E[V''|a] - Var(V'|a),  V = S_W + Phi :
     A2   (canonical, no subtraction) :  S = E[V''|a] ,        D = Var(V'|a)
     H2   (Wilson-Hessian split)      :  S = E[H|a] ,          D = Var(V'|a) - E[Phi''|a]
     Schur (SCBL, INVENT-VARIANCE-INEQ):  needs K^{-1} on the conditioned fibre;
                                          UNDEFINED when K is indefinite  (witness [D4])
  Any split is A2 shifted by a function psi:  S_psi = S_A2 - psi, D_psi = D_A2 - psi,
  hence   D_psi/S_psi = 1 - G''/(S_A2 - psi)   -- the RATIO IS A FUNCTION OF THE SPLIT.

Run: py -3.12 -I ds_ratio.py
"""
import numpy as np

PI = np.pi
# ------------------------------------------------------------------ geometry
def onb(v):
    v = np.asarray(v, float); v = v / np.linalg.norm(v)
    a = np.array([1.0, 0.0, 0.0]) if abs(v[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = np.cross(v, a); e1 /= np.linalg.norm(e1)
    return v, e1, np.cross(v, e1)

def slice_nodes(a0, v, ns, nph):
    """nodes on {<v,xi> = a0} inside |xi|<2pi ; weight = s ds dphi (Lebesgue)."""
    v, e1, e2 = onb(v)
    R = np.sqrt(max(4.0 * PI ** 2 - a0 * a0, 1e-300))
    xs, ws = np.polynomial.legendre.leggauss(ns)
    xp, wp = np.polynomial.legendre.leggauss(nph)
    s = 0.5 * R * (xs + 1.0); wsv = 0.5 * R * ws
    ph = PI * (xp + 1.0);      wpv = PI * wp
    S, PH = np.meshgrid(s, ph, indexing='ij')
    WS, WPH = np.meshgrid(wsv, wpv, indexing='ij')
    xi = (a0 * v[None, None, :]
          + (S * np.cos(PH))[..., None] * e1[None, None, :]
          + (S * np.sin(PH))[..., None] * e2[None, None, :])
    return xi.reshape(-1, 3), (WS * S * WPH).reshape(-1), v

def fields(xi, v, beta, theta):
    """V = S_W + Phi and its first/second v-derivatives at the nodes xi."""
    r = np.linalg.norm(xi, axis=1)
    q = (xi @ v) / r                    # = <v,xi>/rho
    u = xi[:, 2] / r                    # xi_3 / rho
    v3 = v[2]
    ct, st = np.cos(theta / 2.0), np.sin(theta / 2.0)
    h = 0.5 * r; c1, s1 = np.cos(h), np.sin(h)
    cth = c1 / s1                       # cot h
    csch2 = 1.0 / s1 ** 2
    # S_W = beta*(1 - ct*c1 + st*s1*u)
    dvh = 0.5 * q
    d2vh = 0.5 * (1.0 - q * q) / r
    dvu = (v3 - u * q) / r
    d2vu = -2.0 * q * (v3 - u * q) / r ** 2 - u * (1.0 - q * q) / r ** 2
    G = beta * (ct * s1 * dvh + st * (c1 * dvh * u + s1 * dvu))
    H = beta * (ct * (c1 * dvh ** 2 + s1 * d2vh)
                + st * (-s1 * dvh ** 2 * u + c1 * d2vh * u
                        + 2.0 * c1 * dvh * dvu + s1 * d2vu))
    # Phi = 2 log rho - 2 log sin h
    P1 = q * (2.0 / r - cth)
    P2 = 2.0 * (1.0 - q * q) / r ** 2 - cth * (1.0 - q * q) / r + 0.5 * q * q * csch2
    return dict(Sw=beta * (1.0 - ct * c1 + st * s1 * u), Phi=2.0 * np.log(r) - 2.0 * np.log(s1),
                G=G, H=H, P1=P1, P2=P2, r=r, q=q)

def slice_moments(a0, v, beta, theta, ns=400, nph=40, want=('rho', 'EH', 'VarG', 'EP2', 'VarP1', 'Cov', 'EG', 'EP1', 'EV2p', 'EVp')):
    xi, w, vu = slice_nodes(a0, v, ns, nph)
    f = fields(xi, vu, beta, theta)
    # FIXED reference (no per-slice max shift) so log rho is a consistent density in a0
    # (the integrand is bounded: w <= s <= 2pi and exp(-Phi) = J <= 1/4)
    lw = np.log(w) - f['Sw'] - f['Phi']
    ww = np.exp(lw); Z = ww.sum()
    E = lambda x: float((ww * x).sum() / Z)
    EG, EP1 = E(f['G']), E(f['P1'])
    out = dict(rho=np.log(Z), EH=E(f['H']), VarG=E(f['G'] ** 2) - EG ** 2,
               EP2=E(f['P2']), VarP1=E(f['P1'] ** 2) - EP1 ** 2,
               Cov=E(f['G'] * f['P1']) - EG * EP1, EG=EG, EP1=EP1,
               EVp=E(f['G'] + f['P1']), EV2p=E((f['G'] + f['P1']) ** 2))
    return {k: out[k] for k in want}

def splits(m):
    """Return the three splits' (S, D) plus G''.  A2 is the canonical one."""
    S_A2 = m['EH'] + m['EP2']
    D_A2 = m['VarG'] + m['VarP1'] + 2.0 * m['Cov']
    S_H2 = m['EH']
    D_H2 = D_A2 - m['EP2']
    Gpp_A2 = S_A2 - D_A2
    Gpp_H2 = S_H2 - D_H2
    return dict(S_A2=S_A2, D_A2=D_A2, Gpp_A2=Gpp_A2, r_A2=D_A2 / S_A2 if S_A2 != 0 else np.nan,
                S_H2=S_H2, D_H2=D_H2, Gpp_H2=Gpp_H2,
                r_H2=D_H2 / S_H2 if S_H2 != 0 else np.nan)

def fd_curv(a0, v, beta, theta, da=2e-4, ns=600, nph=48):
    """G'' = -d^2/da^2 log rho(a) at a0, 5-point -- independent of the identity."""
    l = [slice_moments(a0 + k * da, v, beta, theta, ns, nph)['rho'] for k in (-2, -1, 0, 1, 2)]
    return -(-l[0] + 16 * l[1] - 30 * l[2] + 16 * l[3] - l[4]) / (12 * da * da)

E3 = np.array([0.0, 0.0, 1.0]); E1 = np.array([1.0, 0.0, 0.0])
def P(*a): print(*a, flush=True)

# ------------------------------------------------------------------ [C1]
def c1_anchors():
    P("=" * 100)
    P("[C1] INDEPENDENT closed-form anchors of MODEL A  (Haar 1D quadrature, own code)")
    P("=" * 100)
    x, w = np.polynomial.legendre.leggauss(4000)
    rho = PI * (x + 1.0); wr = PI * w
    J = np.sin(rho / 2) ** 2 / rho ** 2                    # Haar density in flat chart
    # SLICE measure: dm = J d^2xi, with d^2xi = rho drho dphi in the plane  =>  w = sin^2(rho/2)/rho
    # (the naive 3D radial measure sin^2(rho/2) drho is the WRONG weight here; it gives 0.22571/0.011579)
    m = wr * (np.sin(rho / 2) ** 2 / rho)
    P1 = 2.0 / rho - np.cos(rho / 2) / np.sin(rho / 2)       # Phi'(rho)
    P2s = P1 / rho                                          # Phi'' on the slice {<v,xi>=0}
    A = (m * P2s).sum() / m.sum()
    f = np.sin(rho / 2) / rho
    Ef = (m * f).sum() / m.sum(); Vf = (m * f * f).sum() / m.sum() - Ef ** 2
    P(f"   A = E_Haar[Phi''|0] = {A:.10f}      (archived 0.2051152978)")
    P(f"   B = Var_Haar[sin(rho/2)/rho] = {Vf:.10f}   (archived 0.0110592938)")
    P(f"   beta* = sqrt(A/B) = {np.sqrt(A / Vf):.5f}   (archived 4.30661)")
    P(f"   check 1/(2*Cin(2pi)) = {1.0 / (2.0 * _cin(2 * PI)):.10f}")
    return A, Vf

def _cin(x):
    # Cin(x) = int_0^x (1-cos t)/t dt
    t, w = np.polynomial.legendre.leggauss(4000)
    tt = 0.5 * x * (t + 1.0); ww = 0.5 * x * w
    v = np.where(tt > 1e-14, (1.0 - np.cos(tt)) / np.maximum(tt, 1e-14), 0.0)
    return float((ww * v).sum())

# ------------------------------------------------------------------ [C2]
def c2_closed_form(A, B, beta_star):
    P("")
    P("=" * 100)
    P("[C2] EXACT D/S ON MODEL A at (theta=pi, v=e3, a0=0):   D/S = beta^2/beta*^2 ,  G''=A-B beta^2")
    P("     (conditional measure is EXACTLY Haar: S_W == beta on the slice, E[H|0]=0)")
    P("=" * 100)
    P("   %9s%15s%13s%12s%12s%12s" % ('beta', 'D/S closed', 'D/S quad', 'S=E[V2]', 'D=VarV1', 'Gpp=S-D'))
    for beta in (0.0, 0.5, 1.0, 2.0, 2.4, 2.5, beta_star, 6.0, 10.0):
        m = slice_moments(0.0, E3, max(beta, 1e-12), PI)
        s = splits(m)
        cl = (beta / beta_star) ** 2
        P(f"   {beta:9.5f}{cl:15.6f}{s['r_A2']:13.6f}{s['S_A2']:12.6f}{s['D_A2']:12.6f}"
          f"{s['Gpp_A2']:12.6f}")
    P("   -- FD validation of S-D (5-point, da=2e-4):")
    for beta in (0.0, 2.5, 6.0):
        m = slice_moments(0.0, E3, max(beta, 1e-12), PI); s = splits(m)
        g = fd_curv(0.0, E3, max(beta, 1e-12), PI)
        P(f"      beta={beta:5.2f}:  S-D={s['Gpp_A2']:+.10f}   fd={g:+.10f}   |diff|={abs(s['Gpp_A2']-g):.2e}")
    P(f"   ==>  D/S = 1 exactly at beta* = {beta_star:.5f};  D/S > 1 for every beta > beta*.")
    P(f"        at beta=2.5:  D/S = {(2.5/beta_star)**2:.5f}   (NOT 2%)")

# ------------------------------------------------------------------ [C3]
def c3_grid():
    P("")
    P("=" * 100)
    P("[C3] FOUR-TUPLE TABLE (model, conditioning, v, beta) -- canonical split A2")
    P("     S = E[V''|a0],  D = Var(V'|a0),  G'' = S-D ;  ratio = D/S")
    P("=" * 100)
    hdr = (f"   {'theta':>7}{'v':>4}{'a0':>6}{'beta':>8}{'S':>11}{'D':>11}{'Gpp':>11}"
           f"{'D/S':>10}{'H2 split':>11}{'Gpp>0?':>8}")
    P(hdr)
    for theta in (0.0, PI / 2, PI):
        for vn, v in (('e3', E3), ('e1', E1)):
            for a0 in (0.0, 0.5, 1.0, 1.5, 2.0):
                for beta in (1.0, 2.5, 4.30661, 6.0):
                    m = slice_moments(a0, v, beta, theta)
                    s = splits(m)
                    h2 = ('undef' if abs(s['S_H2']) < 0.01 else '%.5f' % s['r_H2'])
                    P(f"   {theta:7.3f}{vn:>4}{a0:6.2f}{beta:8.5f}{s['S_A2']:11.6f}"
                      f"{s['D_A2']:11.6f}{s['Gpp_A2']:11.6f}{s['r_A2']:10.5f}"
                      f"{h2:>11}{'YES' if s['Gpp_A2'] > 0 else 'no':>8}")

# ------------------------------------------------------------------ [C4]
def c4_crossings():
    P("")
    P("=" * 100)
    P("[C4] WHERE DOES D/S CROSS 1 ?  (bisection on G''=0 ; D/S=1 <=> G''=0 iff S>0)")
    P("=" * 100)
    P(f"   {'theta':>7}{'v':>4}{'a0':>6}{'beta_c(G''=0)':>16}{'S at beta_c':>14}"
      f"{'D/S at beta=2.5':>18}{'violating beta':>16}")
    for theta in (0.0, PI / 2, PI):
        for vn, v in (('e3', E3), ('e1', E1)):
            for a0 in (0.0, 0.5, 1.0, 2.0):
                f = lambda b: splits(slice_moments(a0, v, b, theta))['Gpp_A2']
                lo, hi = 1e-3, 30.0
                flo, fhi = f(lo), f(hi)
                bc = np.nan; Sbc = np.nan
                if flo * fhi < 0:
                    for _ in range(32):
                        mid = 0.5 * (lo + hi); fm = f(mid)
                        if flo * fm <= 0: hi, fhi = mid, fm
                        else: lo, flo = mid, fm
                    bc = 0.5 * (lo + hi)
                    Sbc = splits(slice_moments(a0, v, bc, theta))['S_A2']
                r25 = splits(slice_moments(a0, v, 2.5, theta))['r_A2']
                P(f"   {theta:7.3f}{vn:>4}{a0:6.2f}{bc:16.5f}{Sbc:14.6f}{r25:18.5f}"
                  f"{('>%.5f' % bc) if np.isfinite(bc) else 'none<30':>16}")

# ------------------------------------------------------------------ [C5]
def c5_split_freedom(A, B, beta_star):
    P("")
    P("=" * 100)
    P("[C5] THE RATIO IS A FUNCTION OF THE SPLIT, NOT OF THE MEASURE")
    P("     A2 at (theta=pi,e3,a0=0,beta=2.5):  S=%.6f  D=%.6f  G''=%.6f" %
      (A, B * 6.25, A - B * 6.25))
    P("     any split (S_psi, D_psi) = (S_A2 - psi, D_A2 - psi)  keeps  S-D = G''  and gives")
    P("        D_psi/S_psi = 1 - G''/(S_A2 - psi)")
    P("=" * 100)
    S2, D2 = A, B * 6.25
    P(f"   {'psi':>10}{'S_psi':>11}{'D_psi':>11}{'ratio':>11}   note")
    for psi in (-2.0, -1.0, -0.5, -0.2, 0.0, 0.05, 0.1, 0.15, 0.19, 0.2, 0.5, 1.0):
        Sp, Dp = S2 - psi, D2 - psi
        r = Dp / Sp if abs(Sp) > 1e-14 else np.nan
        note = ''
        if abs(r - 0.02) < 0.05: note = '<-- the paper\'s 2%'
        if abs(psi) < 1e-12: note = '<-- canonical A2 (0.337)'
        P(f"   {psi:10.4f}{Sp:11.6f}{Dp:11.6f}{r:11.5f}   {note}")
    P("   psi_admissible for the Schur/SCBL split:  psi = E[<K^{-1}c,c>|a0], which is")
    P("   UNDEFINED once K is indefinite (tci_comp.log column 'elim? UNDEF').  Proof of range:")
    P("   for psi -> S_A2^-  the ratio -> +inf ; for psi -> -inf the ratio -> 1^- .")
    # H2 split at the degeneracy
    m = slice_moments(0.0, E3, 2.5, PI); s = splits(m)
    P(f"   H2 split here:  S_H2 = E[H|0] = {s['S_H2']:+.3e}  (exactly 0)  =>  ratio UNDEFINED (0/0)")
    P("   PROD reading:   D/S := 1 - G''*Var(a).  With Var(a) != 1/G'' this is NOT D/S:")
    P("                   it equals 1 - G''/S only if the marginal is Gaussian (S = 1/Var).")
    for rr in (0.98, 0.5, 0.05):
        P(f"      if G''*Var(a) = {rr}:  'ratio' = {1-rr:.3f} while the A2 defect ratio = {D2/S2:.3f}")

def main():
    A, B = c1_anchors()
    beta_star = float(np.sqrt(A / B))
    c2_closed_form(A, B, beta_star)
    c3_grid()
    c4_crossings()
    c5_split_freedom(A, B, beta_star)
    P("")
    P("NOTE: no lattice data was read in this file; see ds_acid.py for acid tests + npz audit.")
    P("Clay 0%.  No QED.")

if __name__ == '__main__':
    main()

