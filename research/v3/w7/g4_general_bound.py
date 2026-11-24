#!/usr/bin/env python3
"""
g4_general_bound.py  --  general-block verdict on (star): is G'' >= 2 W_geo, and is W_geo>0 at beta>0?

Three distinct objects are separated (the whole point of G4):
  (L) LITERAL (star) as written in C2:  G''_lit := E[H] - Var(G) + W_j,  W_j = W_geo + X,
      E,Var over mu(.|U') (boundary-conditioned block measure, NO coarse-mode conditioning).
      Identity (proved in the .md):  G''_lit = 0  identically (Stein + cross-Stein).  => (star) lit. FALSE.
  (U) S2's "G''_model":  Delta_dyn + W_geo  (UNCONDITIONAL, no X).  Identity: = 2 W_geo.  NOT (star).
  (C) the TRUE coarse curvature (ch05 eq:split/eq:vardecomp with |0):  the marginal curvature along the
      coarse mode v, conditioned at a = <v,X> = 0:
            G''(v) := -d^2/da^2 log rho_v(a)|_{a=0} = E[ d_v^2 V | a=0 ] - Var( d_v V | a=0 ),
      with V = S + Phi the FULL potential.  This is what (star) must mean.

Model: SU(2) one plaquette / one free link, boundary B = exp(i th sig3/2)  (S2's model).
  V(xi) = h(rho) + g(rho) u ,  u = xi_3/rho,  h = beta(1 - cos(rho/2)cos(th/2)) + Phi(rho),
  g = beta sin(rho/2) sin(th/2),  Phi = 2 log(rho/2) - 2 log sin(rho/2)   (Haar, flat chart, ball rho<2pi).

Deterministic Gauss-Legendre quadrature.  Run:  py -3.12 g4_general_bound.py
"""
import numpy as np

# -------- cancellation-free Phi derivatives (as in s2_plaquette.py) --------
def Phi_rho(x):
    y = np.empty_like(x); sm = np.abs(x) < 0.3
    xr = x[~sm]; y[~sm] = 1.0/xr - 1.0/np.tan(xr)
    xs = x[sm];  y[sm] = xs/3 + xs**3/45 + 2*xs**5/945 + 2*xs**7/4725
    return y
def Phi_rr(x):
    y = np.empty_like(x); sm = np.abs(x) < 0.3
    xr = x[~sm]; y[~sm] = 0.5*(-1.0/xr**2 + 1.0/np.sin(xr)**2)
    xs = x[sm];  y[sm] = 1.0/6 + xs**2/30 + 2*xs**4/189 + 2*xs**6/2475
    return y

# ============================================================
# (A)  S2's UNCONDITIONAL quantities  Delta_dyn, W_geo  (3D quadrature in rho,u)
#      reproduces s2_plaquette.py; validates the model and the identity Ddyn+W=2W.
# ============================================================
_xg, _wg = np.polynomial.legendre.leggauss(700)
RHO = np.pi*(_xg+1.0); WRHO = np.pi*_wg          # rho in (0,2pi)
_ug, _wu = np.polynomial.legendre.leggauss(700)  # u = cos(phi) in [-1,1]
R, U = np.meshgrid(RHO, _ug, indexing='ij')
WR, WU = np.meshgrid(WRHO, _wu, indexing='ij')

def uncond(beta, theta):
    x = R/2.0
    c = np.cos(R/2); s = np.sin(R/2); C = np.cos(theta/2); St = np.sin(theta/2)
    S = beta*(1.0 - (c*C - s*St*U))
    w = np.exp(-S)*np.sin(R/2)**2*WR*WU            # e^{-S} J, J ~ sin^2(rho/2)
    Z = w.sum()
    E = lambda a: (w*a).sum()/Z
    Sr = beta*((s/2)*C + (c/2)*St*U); Su = beta*s*St
    Srr = beta*((c/4)*C - (s/4)*St*U); Sru = beta*(c/2)*St; Suu = 0.0
    pr = Phi_rho(x); prr = Phi_rr(x)
    # v=e3
    G3 = U*Sr + (1-U**2)/R*Su
    H3 = (Srr*U**2 + 2*Sru*U*(1-U**2)/R + Suu*((1-U**2)/R)**2
          + Sr*(1-U**2)/R + Su*(-3*U*(1-U**2)/R**2))
    P3 = U*pr; PP3 = U**2*prr + (1-U**2)/R*pr
    VG3 = E(G3*G3)-E(G3)**2; W3 = E(PP3)-(E(P3*P3)-E(P3)**2)
    Ddyn3 = E(H3)-VG3
    # v=e1 (psi-averaged exactly)
    q = (Sr - U*Su/R); G1sq=(1-U**2)/2*q**2
    H1 = ((1-U**2)/2*Srr + (-U*(1-U**2)/R)*Sru + U**2*(1-U**2)/(2*R**2)*Suu
          + Sr*(1-(1-U**2)/2)/R + Su*U*(3*(1-U**2)/2-1)/R**2)
    P1sq=(1-U**2)/2*pr**2; PP1=(1-U**2)/2*prr+(1+U**2)/2*pr/R
    VG1 = E(G1sq)-0.0; W1 = E(PP1)-(E(P1sq)-0.0); Ddyn1=E(H1)-VG1
    return {'e3':(Ddyn3,W3),'e1':(Ddyn1,W1)}

# ============================================================
# (B)  TRUE coarse curvature (C):  G''(v) = -d^2/da^2 log rho_v(a)|_{a=0}
#      computed two ways per direction, cross-checked:
#        (i) analytic  E[d_v^2 V|0] - Var(d_v V|0)  on the slice a=0
#        (ii) finite-difference curvature of the marginal density rho_v(a)
# ============================================================
_rg, _rgw = np.polynomial.legendre.leggauss(1400)
_rr = np.pi*(_rg+1.0); _rrw = np.pi*_rgw          # r,s in (0,2pi)
_ag, _agw = np.polynomial.legendre.leggauss(400)
_aa = np.pi*(_ag+1.0); _aaw = np.pi*_agw          # alpha in (0,2pi)

def hfun(r,beta,theta): return beta*(1.0-np.cos(r/2)*np.cos(theta/2)) + (2*np.log(r/2)-2*np.log(np.sin(r/2)))
def hpr(r,beta,theta):  return beta*0.5*np.sin(r/2)*np.cos(theta/2) + Phi_rho(r/2)
def gfun(r,beta,theta): return beta*np.sin(r/2)*np.sin(theta/2)

def Gpp_e3_analytic(beta,theta):
    r=_rr; h=hfun(r,beta,theta); hp=hpr(r,beta,theta); g=gfun(r,beta,theta)
    w=np.exp(-h)*r                                # 2D cartesian slice measure dxi1 dxi2 = r dr
    Z=w.sum()
    V3=g/r; V33=hp/r
    E3=(w*V3).sum()/Z; E33=(w*V33).sum()/Z; E33sq=(w*V3*V3).sum()/Z
    return E33-(E33sq-E3**2)

def Gpp_e1_analytic(beta,theta):
    s=_rr[:,None]; a=_aa[None,:]
    h=hfun(s,beta,theta); hp=hpr(s,beta,theta); g=gfun(s,beta,theta)
    # d/ds (g/s)
    dgs = (gfun(s+1e-6,beta,theta)-gfun(np.maximum(s-1e-6,1e-9),beta,theta))/2e-6  # numerical g'
    gps = dgs
    d_g_over_s = gps/s - g/s**2
    V11 = hp/s + d_g_over_s*np.sin(a)
    w = np.exp(-h-g*np.sin(a))*s*_rgw[:,None]*_aaw[None,:]
    Z=w.sum()
    return (w*V11).sum()/Z                           # V_1 = 0 on slice => Var=0

def marginal_e3(a0,beta,theta,n=1400):
    # rho(a0) = 2pi int_0^{sqrt(4pi^2-a0^2)} e^{-V(r,u)} r dr ,  V = h(rho)+g(rho) u, rho=sqrt(r^2+a0^2), u=a0/rho
    Rmax=np.sqrt(max(4*np.pi**2-a0*a0,1e-12))
    xg,wg=np.polynomial.legendre.leggauss(n)
    r=0.5*Rmax*(xg+1.0); wr=0.5*Rmax*wg
    rho=np.sqrt(r*r+a0*a0); u=a0/rho
    V=hfun(rho,beta,theta)+gfun(rho,beta,theta)*u
    return 2*np.pi*np.sum(np.exp(-V)*r*wr)

def Gpp_e3_fd(beta,theta,da=1e-3):
    # -d^2/dx^2 log rho(x) at 0 by 5-point stencil
    vals=[marginal_e3(k*da,beta,theta) for k in (-2,-1,0,1,2)]
    logs=np.array([np.log(v) for v in vals])
    # f'' ~ (-f[-2]+16 f[-1]-30 f[0]+16 f[1]-f[2])/(12 da^2)
    fpp=(-logs[0]+16*logs[1]-30*logs[2]+16*logs[3]-logs[4])/(12*da*da)
    return -fpp

def log(*a): print(*a, flush=True)

# ============================================================
log("="*78)
log("(A) S2 UNCONDITIONAL identity:  Delta_dyn = W_geo,  Delta_dyn+W_geo = 2 W_geo  (3D quad)")
log(f"{'beta':>5} {'th':>6} {'dir':>3} {'Delta_dyn':>11} {'W_geo':>11} {'2W_geo':>11} {'|D-W|':>9}")
for beta in (0.5,2.5,5.0,10.0):
    for theta in (0.0,np.pi/2,np.pi):
        q=uncond(beta,theta)
        for d in ('e3','e1'):
            Dd,W=q[d]
            log(f"{beta:5.2f} {theta:6.3f} {d:>3} {Dd:11.6f} {W:11.6f} {2*W:11.6f} {abs(Dd-W):9.1e}")

log("\n"+"="*78)
log("(B) TRUE coarse curvature (C)  G''(v) = -d^2/da^2 log rho|_0   vs  S2's 2W_geo  vs  0")
log(f"{'beta':>5} {'th':>6} {'dir':>3} {'Gpp(analytic)':>13} {'Gpp(fd e3)':>11} {'2W_geo':>11} {'Gpp-2W':>11} {'Gpp<0?':>6}")
for beta in (0.1,0.5,1.0,2.5,5.0,10.0):
    for theta in (0.0,np.pi/2,np.pi):
        q=uncond(beta,theta)
        for d in ('e3','e1'):
            G = Gpp_e3_analytic(beta,theta) if d=='e3' else Gpp_e1_analytic(beta,theta)
            W = q[d][1]; two=2*W
            fd = Gpp_e3_fd(beta,theta) if d=='e3' else np.nan
            log(f"{beta:5.2f} {theta:6.3f} {d:>3} {G:13.6f} {fd:11.6f} {two:11.6f} {G-two:11.6f} {str(G<0):>6}")

log("\n"+"="*78)
log("(C) ACID TESTS")
log("  free / flat group: Phi=0 -> W_geo=0 and Gaussian Var=E[H] -> G''=0 (no false positivity).")
# explicit flat Gaussian (1D analog): V=beta*x^2/2 with no Jacobian -> G''=beta-... check 0
# Use model with Phi forced to 0 (J=1): the measure is exp(-S) dxi (flat), S quadratic? Our S is
# beta(1-cos..), non-quadratic; instead check beta->0 pure-Haar curvature vs W_geo.
log("  beta->0 pure Haar (beta=1e-4): G''(e3) > 0  while  W_geo(uncond) -> 0 :")
q=uncond(1e-4,np.pi)
log(f"     Gpp(e3)={Gpp_e3_analytic(1e-4,np.pi):.6f}  Gpp(e1)={Gpp_e1_analytic(1e-4,np.pi):.6f}  W_geo(e3)={q['e3'][1]:.2e}  2W={2*q['e3'][1]:.2e}")

log("\n"+"="*78)
log("(D) COUNTEREXAMPLE SEARCH:  strong boundary field, is G'' < 0 anywhere?")
log(f"{'beta':>5} {'th':>6} {'dir':>3} {'Gpp':>12}")
worst=(1e9,None)
for beta in (2.5,5.0,10.0,20.0):
    for theta in (np.pi*0.9, np.pi):
        for d in ('e3','e1'):
            G=Gpp_e3_analytic(beta,theta) if d=='e3' else Gpp_e1_analytic(beta,theta)
            log(f"{beta:5.2f} {theta:6.3f} {d:>3} {G:12.6f}")
            if G<worst[0]: worst=(G,(beta,theta,d))
log(f"  MIN observed: G''={worst[0]:.6f} at (beta,theta,dir)={worst[1]}")

# zero crossing along theta=pi, e3 (G'' = E_Haar[Phi'/r] - beta^2 Var_Haar[sin(r/2)/r])
log("\n  beta-scan at theta=pi, v=e3 (the analytic G''=A - beta^2 B):")
for beta in (0.0,1.0,2.0,3.0,4.0,5.0,6.0,8.0,10.0):
    log(f"     beta={beta:5.2f}  Gpp(e3)={Gpp_e3_analytic(max(beta,1e-6),np.pi):12.6f}")
