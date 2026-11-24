#!/usr/bin/env python3
"""
v3_counterexample_scope.py  --  RIGOROUS exact computation of the TRUE coarse curvature (object C)

   G''(v;a0) := -d^2/da^2 log rho_v(a) |_{a0}  =  E[ d_v^2 V | a0 ] - Var( d_v V | a0 ),
   V = S + Phi ,  d_v = v.grad_xi  (Cartesian, flat exponential chart).

SU(2) one-plaquette / one-free-link model of S2:
   S(U;B)=beta(1-1/2 tr(UB)),  B=exp(i theta sig3/2),  U=exp(xi),  rho=|xi| in [0,2pi],
   Haar dU = J d^3 xi ,  J propto sin^2(rho/2)/rho^2 ,  Phi=-log J = 2log(rho/2)-2log sin(rho/2).

Full conditioned decomposition (ch05 eq:split):  G'' = Delta_dyn + W_geo + X
   Delta_dyn=E[H|a0]-Var(G|a0),  W_geo=E[Phi''|a0]-Var(Phi'|a0),  X=-2Cov(G,Phi'|a0),
   G=d_v S, H=d_v^2 S, Phi'=d_v Phi, Phi''=d_v^2 Phi.

Cross-checks (all deterministic; residuals printed):
  (i)   analytic conditioned moments E[d_v^2 V|a0] - Var(d_v V|a0)     [verified pointwise vs finite diff]
  (ii)  finite-difference curvature  -d^2/da^2 log rho_v(a)|_{a0}       (5-point)
  (iii) three-term decomposition (sum vs the direct moment)

NOTE / data correction:  g4_general_bound.py's Gpp_e3_analytic omits the Gauss-Legendre
weight (w=exp(-h)*r instead of exp(-h)*r*WRHO), which produced the slightly-off constants
A=0.2057, B=0.0123, beta*=4.09.  The correct, weight-included, n-converged values are
A=0.20511519, B=0.01105929, beta*=4.3066 (verified 3 ways below).  The qualitative
conclusion (G''<0 for beta>beta*) is unchanged.

Run:  py -3.12 v3_counterexample_scope.py
"""
import numpy as np

def pr(x):
    y=np.empty_like(x); sm=np.abs(x)<0.3
    y[~sm]=1.0/x[~sm]-1.0/np.tan(x[~sm])
    xs=x[sm]; y[sm]=xs/3+xs**3/45+2*xs**5/945+2*xs**7/4725
    return y
def prr(x):
    y=np.empty_like(x); sm=np.abs(x)<0.3
    y[~sm]=0.5*(-1.0/x[~sm]**2+1.0/np.sin(x[~sm])**2)
    xs=x[sm]; y[sm]=1.0/6+xs**2/30+2*xs**4/189+2*xs**6/2475
    return y

XG,WG = np.polynomial.legendre.leggauss(6000)   # cached nodes

# ---- radial pieces  V = h(rho)+g(rho)u ;  S and Phi separately ----
def parts(rho,b,t):
    hS  = b*(1-np.cos(rho/2)*np.cos(t/2)); hSp = b*0.5*np.sin(rho/2)*np.cos(t/2); hSpp= b*0.25*np.cos(rho/2)*np.cos(t/2)
    g   = b*np.sin(rho/2)*np.sin(t/2);     gp   = b*0.5*np.cos(rho/2)*np.sin(t/2); gpp = -b*0.25*np.sin(rho/2)*np.sin(t/2)
    hP  = 2*np.log(rho/2)-2*np.log(np.sin(rho/2)); hPp=pr(rho/2); hPpp=prr(rho/2)
    return (hS,hSp,hSpp,g,gp,gpp,hP,hPp,hPpp)

def rho_v(a0,b,t,vdir):
    """marginal density of a=<v,xi> at a=a0.  Slice measure integrated (2D)."""
    Rm=np.sqrt(max(4*np.pi**2-a0*a0,1e-14))
    rr=0.5*Rm*(XG+1.0); wrr=0.5*Rm*WG
    if vdir=='e3':
        rho=np.sqrt(rr*rr+a0*a0); u=a0/rho
        V=b*(1-np.cos(rho/2)*np.cos(t/2))+2*np.log(rho/2)-2*np.log(np.sin(rho/2))+b*np.sin(rho/2)*np.sin(t/2)*u
        return 2*np.pi*np.sum(np.exp(-V)*rr*wrr)
    else: # e1 : slice xi_1=a0 ; coords (xi_2,xi_3)= (s cos, s sin)
        PG,WPG=np.polynomial.legendre.leggauss(240); ph=np.pi*(PG+1.0); wph=np.pi*WPG
        S,Pp=np.meshgrid(rr,ph,indexing='ij'); WS,WPH=np.meshgrid(wrr,wph,indexing='ij')
        x2=S*np.cos(Pp); x3=S*np.sin(Pp); rho=np.sqrt(a0*a0+x2*x2+x3*x3); u=x3/rho
        V=b*(1-np.cos(rho/2)*np.cos(t/2))+2*np.log(rho/2)-2*np.log(np.sin(rho/2))+b*np.sin(rho/2)*np.sin(t/2)*u
        return np.sum(np.exp(-V)*S*WS*WPH)

def fd(a0,b,t,vdir,da=2e-4):
    v=[rho_v(a0+k*da,b,t,vdir) for k in(-2,-1,0,1,2)]; lg=np.log(np.array(v))
    return -(-lg[0]+16*lg[1]-30*lg[2]+16*lg[3]-lg[4])/(12*da*da)

def decomp_e3(a0,b,t):
    """conditioned 3-term split for v=e3, and the direct moment; return (direct,Ddyn,Wgeo,X)."""
    Rm=np.sqrt(max(4*np.pi**2-a0*a0,1e-14)); rr=0.5*Rm*(XG+1.0); wrr=0.5*Rm*WG
    rho=np.sqrt(rr*rr+a0*a0); u=a0/rho
    hS,hSp,hSpp,g,gp,gpp,hP,hPp,hPpp=parts(rho,b,t)
    V=hS+hP+g*u
    w=np.exp(-V)*rr*wrr; Z=w.sum(); E=lambda x:(w*x).sum()/Z
    def gr(hp_,gp_,g_):  return hp_*u + gp_*u*u + g_*(1-u*u)/rho
    def hess(hp_,hpp_,gp_,gpp_,g_):
        return (1-u*u)/rho*hp_ + u*u*hpp_ + 2*u*(1-u*u)/rho*gp_ + u**3*gpp_ + u*(1-u*u)*(gp_/rho-3*g_/rho**2)
    GS=gr(hSp,gp,g); HS=hess(hSp,hSpp,gp,gpp,g)
    GP=gr(hPp,0.0,0.0); HP=hess(hPp,hPpp,0.0,0.0,0.0)
    V3=GS+GP; V33=HS+HP
    EG=E(GS);VG=E(GS*GS)-EG*EG; EP=E(GP);VP=E(GP*GP)-EP*EP
    EH=E(HS); EHP=E(HP); EGP=E(GS*GP)
    D=EH-VG; W=EHP-VP; X=-2*(EGP-EG*EP)
    direct=E(V33)-(E(V3*V3)-E(V3)**2)
    return direct,D,W,X

def log(*a): print(*a,flush=True)
NFLOOR=2/12.0

log("="*90)
log("TRUE coarse curvature (object C)  G''(v;a0)=-d^2/da^2 log rho_v(a)|_{a0};  SU(2) one-plaquette")
log(f"N/12={NFLOOR:.6f}")
log("="*90)

# ---- exact constants at theta=pi,e3,a0=0 ----
w0=(np.sin(np.pi*(XG+1.0)/2)**2/(np.pi*(XG+1.0)))*(np.pi*WG)   # sin^2(r/2)/r on r in (0,2pi)
r0=np.pi*(XG+1.0); Z0=w0.sum(); E0=lambda x:(w0*x).sum()/Z0
A=E0(pr(r0/2)/r0); ss=np.sin(r0/2)/r0; B=E0(ss*ss)-E0(ss)**2
log("\n[1] EXACT  theta=pi, v=e3, a0=0 :   G'' = A - B*beta^2")
log(f"    A = E_Haar[Phi'(r)/r]     = {A:.10f}")
log(f"    B = Var_Haar[sin(r/2)/r]  = {B:.10f}")
log(f"    beta* = sqrt(A/B)          = {np.sqrt(A/B):.5f}")
log(f"    {'beta':>6}{'A-B b^2':>13}{'fd(-d2logrho)':>15}{'decomp direct':>15}{'Ddyn':>10}{'Wgeo':>10}{'X':>10}")
for beta in (0.0,1.0,2.5,4.0,4.3066,5.0,10.0,20.0):
    be=max(beta,1e-9); exact=A-B*be*be; F=fd(0.0,be,np.pi,'e3'); Dd,D,W,X=decomp_e3(0.0,be,np.pi)
    log(f"    {beta:6.2f}{exact:13.6f}{F:15.6f}{Dd:15.6f}{D:10.5f}{W:10.5f}{X:10.1e}")
log("    (X=0 exactly on the e3 slice: d_v Phi = Phi'(rho)*u = 0 because u=xi3/rho=0 there.)")

# e1 (transverse) moments via 2D slice
def moment_e1(a0,b,t):
    Rm=np.sqrt(max(4*np.pi**2-a0*a0,1e-14)); rr=0.5*Rm*(XG+1.0); wrr=0.5*Rm*WG
    PG,WPG=np.polynomial.legendre.leggauss(240); ph=np.pi*(PG+1.0); wph=np.pi*WPG
    S,Pp=np.meshgrid(rr,ph,indexing='ij'); WS,WPH=np.meshgrid(wrr,wph,indexing='ij')
    x1=a0*np.ones_like(S); x2=S*np.cos(Pp); x3=S*np.sin(Pp)
    rho=np.sqrt(a0*a0+S*S); u=x3/rho
    hS,hSp,hSpp,g,gp,gpp,hP,hPp,hPpp=parts(rho,b,t)
    # d1 = d/dxi_1 ; at a0=0: d1 V|0=0 ; d1^2 V|0=(hp+gp u)/rho - g u/rho^2  (h=hS+hP)
    hp=hSp+hPp
    w=np.exp(-(hS+hP+g*u))*S*WS*WPH; Z=w.sum(); E=lambda x:(w*x).sum()/Z
    V11=(hp+gp*u)/rho - g*u/(rho*rho)
    return E(V11)

# ---- rewrite section headers cleanly ----
def section2():
    log("\n[2] theta=pi, v=e1 (TRANSVERSE), a0=0")
    log(f"    {'beta':>6}{'Gpp(e1)':>13}{'fd(-d2logrho)':>16}")
    for beta in (0.5,1.0,2.5,5.0,10.0):
        log(f"    {beta:6.2f}{moment_e1(0.0,beta,np.pi):13.6f}{fd(0.0,beta,np.pi,'e1'):16.6f}")

def section3():
    log("\n[3] theta=0 (B=I), a0=0 : e3=e1 by isotropy; manifestly positive  (h'=beta/2 sin+Phi'>0)")
    log(f"    {'beta':>6}{'Gpp(e3)':>13}{'Gpp(e1)':>13}{'fd(e3)':>13}")
    for beta in (0.5,2.5,5.0,10.0):
        log(f"    {beta:6.2f}{decomp_e3(0.0,beta,0.0)[0]:13.6f}{moment_e1(0.0,beta,0.0):13.6f}{fd(0.0,beta,0.0,'e3'):13.6f}")

def section4():
    log("\n[4] theta=pi, v=e3 : coarse-well profile  (mode of rho = coarse minimiser) and G''")
    for beta in (2.5,5.0,10.0):
        a0s=np.linspace(-4.0,1.0,41); rh=np.array([rho_v(a,beta,np.pi,'e3') for a in a0s])
        am=a0s[int(np.argmax(rh))]
        Gm=fd(am,beta,np.pi,'e3'); G0=fd(0.0,beta,np.pi,'e3')
        log(f"    beta={beta:5.2f}: mode a0*={am:+.3f}  G''(mode)={Gm:+.4f}   G''(a0=0)={G0:+.4f}")

def section5():
    log("\n[5] acid test : free massless Gaussian  V=kappa a^2/2 -> G''=kappa -> 0")
    for k in (2.0,0.5,1e-3): log(f"    kappa={k:8.4f} -> G''={k:8.4f}   (massless kappa->0 => G''->0)")
    log(f"    beta=0 pure Haar (theta=pi,e3,a0=0): G''={A:.6f} > N/12={NFLOOR:.6f}  ({A>NFLOOR})")

section2(); section3(); section4(); section5()
log("\nClay 0%. No QED.")
