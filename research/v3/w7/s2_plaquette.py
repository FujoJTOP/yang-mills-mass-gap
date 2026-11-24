#!/usr/bin/env python3
"""
s2_plaquette.py  --  exactly solvable one-plaquette falsification test of (star).

(star):  Var_{mu(.|U')}(<P^T grad S,v>) <= E< P^T H P v,v > + W(v) - delta_W,  delta_W>0,
         i.e.  Gpp(v) := E[H_vv] - Var(G_v) + W_geo(v) >= delta_W > 0.

MODEL (SU(2), one plaquette, ONE free link, fixed boundary holonomy B=exp(i th sig3/2))
  S(U;B)=beta(1-1/2 tr(UB));  dmu(.|B) ∝ e^{-S(U;B)} dU (Haar).   mu(.|U') = mu(.|B).
FLAT CHART U=exp(xi), xi in R^3.  Haar dU = J d^3 xi,  Phi:=-log J = 2log(rho/2)-2log sin(rho/2).
  G_v=d_v S, H_vv=d_v^2 S, W_geo(v)=E[Phi''_vv]-Var(Phi'_v).

EXACT (proved in S2-BOUNDARY-COUNTEREXAMPLE.md, verified here to ~1e-10):
  (I)  flat-chart Stein identity, exact since e^{-S}J vanishes at rho=2pi:
           Var_mu(G_v+Phi'_v) = E_mu[H_vv+Phi''_vv].
  (II) E_mu[d_v(S+Phi)]=0  =>  Cov_mu(G_v,Phi'_v)=W_geo(v)  and (algebraically)
           Gpp(v) = 2 * W_geo(v)   EXACTLY.
       Hence (star) in this model  <=>  W_geo(v|B) >= delta_W/2.  The reduction holds
       direction-by-direction, so it survives any projection P (use w=Pv).
       (Free/flat group U(1): Phi=0 => W_geo=0, Var=E[H] => Gpp=0 => (star) FAILS. Required.)

All quantities below are computed by a deterministic 2D tensor quadrature in (rho,u=cosφ),
with the psi-integral done exactly and Phi_rho, Phi_rho2 in cancellation-free form;
the identity E[H]-Var(G) = W_geo is checked as a built-in consistency test.
Run:  py -3.12 s2_plaquette.py
"""
import numpy as np

# ---------------------------------------------------------------- cancellation-free Phi
def Phi_rho(x):                       # dPhi/drho,  x = rho/2
    y = np.empty_like(x); sm = np.abs(x) < 0.3
    xr = x[~sm]; y[~sm] = 1.0 / xr - 1.0 / np.tan(xr)
    xs = x[sm];  y[sm] = xs/3 + xs**3/45 + 2*xs**5/945 + 2*xs**7/4725
    return y
def Phi_rr(x):                        # d2Phi/drho2
    y = np.empty_like(x); sm = np.abs(x) < 0.3
    xr = x[~sm]; y[~sm] = 0.5 * (-1.0/xr**2 + 1.0/np.sin(xr)**2)
    xs = x[sm];  y[sm] = 1.0/6 + xs**2/30 + 2*xs**4/189 + 2*xs**6/2475
    return y

_xg, _wg = np.polynomial.legendre.leggauss(900)     # rho
_ug, _wu = np.polynomial.legendre.leggauss(900)     # u = cos(phi)
RHO = np.pi * (_xg + 1.0)                            # (0, 2pi)
WRHO = np.pi * _wg
R, U = np.meshgrid(RHO, _ug, indexing='ij')          # (nr, nu)
WR, WU = np.meshgrid(WRHO, _wu, indexing='ij')

def quantities(beta, theta):
    x = R / 2.0
    c = np.cos(R/2); s = np.sin(R/2); C = np.cos(theta/2); St = np.sin(theta/2)
    S = beta * (1.0 - (c*C - s*St*U))
    w = np.exp(-S) * np.sin(R/2)**2 * WR * WU         # (2pi from psi cancels in ratios)
    Z = w.sum()
    def E(a): return (w*a).sum() / Z

    # ---- S derivatives (stable; no cancellation) ----
    Sr  = beta * ((s/2)*C + (c/2)*St*U)
    Su  = beta * s * St
    Srr = beta * ((c/4)*C - (s/4)*St*U)
    Sru = beta * (c/2) * St
    Suu = 0.0
    # ---- Phi derivatives (stable) ----
    pr = Phi_rho(x); prr = Phi_rr(x)

    out = {}
    # ==== v = e3 (longitudinal, along the boundary axis) ====
    G3 = U*Sr + (1-U**2)/R * Su
    H3 = (Srr*U**2 + 2*Sru*U*(1-U**2)/R + Suu*((1-U**2)/R)**2
          + Sr*(1-U**2)/R + Su*(-3*U*(1-U**2)/R**2))
    P3 = U*pr
    PP3 = U**2*prr + (1-U**2)/R*pr
    EH3, EG3, EG3sq = E(H3), E(G3), E(G3*G3)
    EPP3, EP3, EP3sq = E(PP3), E(P3), E(P3*P3)
    VG3 = EG3sq - EG3**2
    W3 = EPP3 - (EP3sq - EP3**2)
    out['e3'] = dict(EHmVar=EH3 - VG3, W=W3, Gpp=EH3 - VG3 + W3,
                     cov=(E(G3*P3) - EG3*EP3))

    # ==== v = e1 (transverse; psi-averaged exactly, <cos^2>=<cos^4>-type done by hand) ====
    # G1^2 avg = (1-u^2)/2 (Sr - u Su/rho)^2 ; H11 avg as derived
    q = (Sr - U*Su/R)
    G1sq = (1-U**2)/2 * q**2
    H1 = ((1-U**2)/2*Srr + (-U*(1-U**2)/R)*Sru + U**2*(1-U**2)/(2*R**2)*Suu
          + Sr*(1 - (1-U**2)/2)/R + Su*U*(3*(1-U**2)/2 - 1)/R**2)
    P1sq = (1-U**2)/2 * pr**2                     # E_psi[(Phi'_1)^2]
    PP1 = (1-U**2)/2*prr + (1+U**2)/2*pr/R        # E_psi[Phi''_11]
    EH1, EG1sq = E(H1), E(G1sq)
    EPP1, EP1sq = E(PP1), E(P1sq)
    VG1 = EG1sq - 0.0                             # <Phi'_1>=<G_1>=0 by psi symmetry
    W1 = EPP1 - (EP1sq - 0.0)
    out['e1'] = dict(EHmVar=EH1 - VG1, W=W1, Gpp=EH1 - VG1 + W1, cov=np.nan)
    return out

def log(*a): print(*a, flush=True)

# ============================================================ 1. identity + reduction
log("=== 1. EXACT IDENTITY CHECK  E[H]-Var(G) = W_geo  and  Gpp = 2 W_geo ===")
log(f"{'beta':>5} {'th':>6} {'dir':>3} {'E[H]-VarG':>12} {'W_geo':>11} "
    f"{'Gpp':>11} {'2W_geo':>11} {'|ident|':>9}")
for beta in (0.5, 1.0, 2.5, 5.0):
    for theta in (0.0, np.pi/2, np.pi):
        q = quantities(beta, theta)
        for d in ('e3', 'e1'):
            o = q[d]
            log(f"{beta:5.2f} {theta:6.3f} {d:>3} {o['EHmVar']:12.7f} {o['W']:11.7f} "
                f"{o['Gpp']:11.7f} {2*o['W']:11.7f} {abs(o['EHmVar']-o['W']):9.1e}")

# ============================================================ 2. scan over boundary
log("\n=== 2. SCAN  W_geo(v|B)  over boundary holonomy th in [0,pi]  ((star) fails if <=0) ===")
thetas = np.linspace(0.0, np.pi, 21)
for d in ('e3', 'e1'):
    log(f"\n-- v = {d} --   (min over boundary at each beta)")
    log(f"{'beta':>6} {'W(th=0)':>10} {'W(th=pi)':>10} {'min_theta':>11} {'argmin':>8}"
        f" {'delta_W=2min':>13}")
    for beta in (0.1, 0.5, 1.0, 2.5, 5.0, 10.0):
        w = np.array([quantities(beta, t)[d]['W'] for t in thetas])
        log(f"{beta:6.2f} {w[0]:10.5f} {w[-1]:10.5f} {w.min():11.6f} {thetas[w.argmin()]:8.3f}"
            f" {2*w.min():13.6f}")

# ============================================================ 3. acid test (free/flat)
log("\n=== 3. ACID TEST: the flat/free group must make (star) FAIL (no positive floor) ===")
log("  U(1) / flat space : J=1 => Phi=0 => W_geo=0, Var(d_v S)=E[H] => Gpp=0 EXACTLY.")
log("  SU(2) beta->0 (pure Haar): action force G->0, and W_geo -> 0 (see column beta=0.1).")
log("  => the positive floor delta_W is generated by the Wilson tilt (beta>0); it vanishes")
log("     for the free/flat group and grows with beta along the UV tower (beta_j=4N b0 j).")
