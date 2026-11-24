#!/usr/bin/env python3
"""
g3_projection.py  --  witness for G3-PROJECTION-TERM.md

Question:  S2 proved, for ONE fixed fine direction v (flat chart, constant vector),
             Gpp(v) := E[H_vv] - Var(G_v) + W_geo(v)  =  2 W_geo(v)   EXACTLY,
   with  G_v = d_v S,  H_vv = d_v^2 S,  W_geo(v)=E[Phi''_vv]-Var(Phi'_v).
The coarse mode is v = P w.  We test whether the identity survives a projection.

  CASE (a)  constant-coefficient (flat-chart-linear) P  ==>  w = Pv is a FIXED fine
            vector; S2 applies verbatim; Gpp(w)=2W_geo(w).  (analytic, checked here
            for s=0 below: R_proj = 0 to machine/quadrature.)

  CASE (b)  covariant / configuration-dependent P  (parallel transport):  the coarse
            direction pulls back to a VECTOR FIELD  w(X) = Ad_{Omega(X)^-1} v.
            Model the transporter by a path of length proportional to s:
                w(X) = R(-2 s |X|, X/|X|) v ,   v in su(2)=R^3,
            which is  Ad_{exp(-sX)} v  (s=0 -> w=v constant -> recover S2).
            Then (exact vector-field Stein, derived in the note):
                Gpp_P = 2 W_geo(w) + R_proj ,
                R_proj = E[kappa] - E[kappa_S] + Cov(phi,d) - Cov(G,d),
                d = div w,  kappa = <(Dw)w, grad Phi>,  kappa_S = <(Dw)w, grad S>,
                phi = <w,grad Phi>,  G = <w,grad S>.
            We compute Gpp_P, W_geo(w), R_proj directly by deterministic 3-D
            quadrature in (rho, u=cos theta, psi), and verify R_proj against the
            analytic RHS.

Model: SU(2) one plaquette, ONE free link, boundary holonomy B=exp(i th sig3/2),
   S(xi)=beta(1 - 1/2 tr(U B)) = beta(1 - cos(r/2)cos(th/2) + sin(r/2)sin(th/2) u),
   u = xi_3/|xi| ,  flat chart Haar  dU = J d^3 xi ,  Phi = -log J = 2log(r/2)-2log sin(r/2).

Acid test: free/flat group (Phi=0, Gaussian S): the floor must NOT be manufactured.
Run:  py -3.12 g3_projection.py
"""
import numpy as np

# ------------------------------------------------------------------ cancellation-free Phi
def Phi_r(x):                       # dPhi/drho ,  x = rho/2
    y = np.empty_like(x); sm = np.abs(x) < 0.3
    xr = x[~sm]; y[~sm] = 1.0/xr - 1.0/np.tan(xr)
    xs = x[sm];  y[sm] = xs/3 + xs**3/45 + 2*xs**5/945 + 2*xs**7/4725
    return y
def Phi_rr(x):                      # d2Phi/drho2
    y = np.empty_like(x); sm = np.abs(x) < 0.3
    xr = x[~sm]; y[~sm] = 0.5*(-1.0/xr**2 + 1.0/np.sin(xr)**2)
    xs = x[sm];  y[sm] = 1.0/6 + xs**2/30 + 2*xs**4/189 + 2*xs**6/2475
    return y

# ------------------------------------------------------------------ 3-D GL quadrature
NR, NU, NP = 90, 56, 56
xr_, wr_ = np.polynomial.legendre.leggauss(NR)     # rho
xu_, wu_ = np.polynomial.legendre.leggauss(NU)     # u = cos theta
xp_, wp_ = np.polynomial.legendre.leggauss(NP)     # psi
RHO = np.pi*(xr_+1.0)      # (0, 2pi)
WRHO = np.pi*wr_
U   = xu_                  # (-1,1)
WU  = wu_
PSI = np.pi*(xp_+1.0)      # (0, 2pi)
WPS = np.pi*wp_
# build mesh  (rho, u, psi)
Rg, Ug, Pg = np.meshgrid(RHO, U, PSI, indexing='ij')
WRg, WUg, WPSg = np.meshgrid(WRHO, WU, WPS, indexing='ij')
st = np.sqrt(np.maximum(1.0-Ug**2, 0.0))
XI = np.stack([Rg*st*np.cos(Pg), Rg*st*np.sin(Pg), Rg*Ug], axis=-1)  # (...,3)
RHO3 = Rg
UU   = Ug
# measure weight (J d^3 xi = e^{-Phi} * rho^2 drho du dpsi ; J rho^2 = sin^2(rho/2))
BASE = np.sin(RHO3/2.0)**2 * WRg * WUg * WPSg

def Sfun(xi, beta, th):
    r = np.linalg.norm(xi, axis=-1)
    u = np.where(r > 1e-12, xi[...,2]/np.maximum(r,1e-12), 0.0)
    c = np.cos(r/2.0); s = np.sin(r/2.0)
    return beta*(1.0 - c*np.cos(th/2.0) + s*np.sin(th/2.0)*u)

def grad_hess_S(xi, beta, th, h=1e-4):
    """central-difference grad and 3x3 Hessian of S at each point xi (N,3)."""
    N = xi.shape[0]
    g = np.zeros((N,3)); Hd = np.zeros((N,3)); H  = np.zeros((N,3,3))
    Sp = [None]*3; Sm = [None]*3
    for i in range(3):
        ei = np.zeros(3); ei[i]=h
        Sp[i] = Sfun(xi+ei, beta, th); Sm[i] = Sfun(xi-ei, beta, th)
        g[:,i] = (Sp[i]-Sm[i])/(2*h)
        Hd[:,i] = (Sp[i]-2*Sfun(xi,beta,th)+Sm[i])/h**2
    for i in range(3):
        for j in range(i+1,3):
            ei=np.zeros(3); ej=np.zeros(3); ei[i]=h; ej[j]=h
            v = (Sfun(xi+ei+ej,beta,th)-Sfun(xi+ei-ej,beta,th)
                 -Sfun(xi-ei+ej,beta,th)+Sfun(xi-ei-ej,beta,th))/(4*h*h)
            H[:,i,j]=v; H[:,j,i]=v
    for i in range(3): H[:,i,i]=Hd[:,i]
    return g, H

def grad_hess_Phi(xi, h=1e-4):
    """Phi is radial: Phi=Phi(rho).  grad = Phi_r * rhat ; Hess = Phi_rr rhat rhat + (Phi_r/rho)(I-rhat rhat)."""
    N = xi.shape[0]
    r = np.linalg.norm(xi, axis=-1)
    rhat = xi/np.maximum(r,1e-12)[:,None]
    x = r/2.0
    pr = Phi_r(x); prr = Phi_rr(x)
    overr = np.where(r>1e-9, pr/np.maximum(r,1e-12), 0.0)
    g = pr[:,None]*rhat
    H = np.zeros((N,3,3))
    for a in range(3):
        for b in range(3):
            H[:,a,b] = prr*rhat[:,a]*rhat[:,b] + overr*( (1.0 if a==b else 0.0) - rhat[:,a]*rhat[:,b])
    return g, H

def transport_w(xi, v, s, h=1e-4):
    """w(xi)=Ad_{exp(-s X)} v = Rodrigues rotation of v about n=xi/|xi| by angle -2 s |xi|.
       Returns w (N,3) and Dw (N,3,3)= d w_a / d xi_b."""
    def wfun(x):
        r = np.linalg.norm(x, axis=-1)
        n = x/np.maximum(r,1e-12)[:,None]
        al = -2.0*s*r                       # rotation angle
        ca = np.cos(al)[:,None]; sa = np.sin(al)[:,None]
        nv = np.cross(n, v)                  # (N,3)
        ndv = (n@v)[:,None]
        return v[None,:]*ca + nv*sa + n*ndv*(1.0-ca)
    w = wfun(xi)
    Dw = np.zeros(xi.shape[:-1]+(3,3))
    for b in range(3):
        eb=np.zeros(3); eb[b]=h
        Dw[...,:,b] = (wfun(xi+eb)-wfun(xi-eb))/(2*h)
    return w, Dw

def compute(beta, th, vdir, s, Phi_on=True):
    xi = XI.reshape(-1,3)
    Np = xi.shape[0]
    v = np.array(vdir, dtype=float)
    S  = Sfun(xi, beta, th)
    w, Dw = transport_w(xi, v, s)
    gS, HS = grad_hess_S(xi, beta, th)
    if Phi_on:
        gP, HP = grad_hess_Phi(xi)
    else:
        gP = np.zeros_like(gS); HP = np.zeros_like(HS)
    # projected objects
    G  = np.einsum('na,na->n', w, gS)           # <w,grad S>
    Hq = np.einsum('na,nab,nb->n', w, HS, w)    # w^T Hess S w
    phi= np.einsum('na,na->n', w, gP)
    Psi= np.einsum('na,nab,nb->n', w, HP, w)
    # vector-field extras
    d   = np.einsum('naa->n', Dw)                # div w
    Dw_w= np.einsum('nab,nb->na', Dw, w)         # (Dw)w  (component a)
    kap  = np.einsum('na,na->n', Dw_w, gP)       # <(Dw)w,grad Phi>
    kapS = np.einsum('na,na->n', Dw_w, gS)       # <(Dw)w,grad S>
    # density (drop constant 2pi irrelevant).  J d^3 = sin^2(rho/2) drho du dpsi
    dens = np.exp(-S)*np.sin(np.linalg.norm(xi,axis=-1)/2.0)**2
    wt = dens*(WRg.reshape(-1)*WUg.reshape(-1)*WPSg.reshape(-1))
    Z = wt.sum()
    def E(a): return (wt*a).sum()/Z
    def Var(a):
        m=E(a); return E((a-m)**2)
    def Cov(a,b):
        return E(a*b)-E(a)*E(b)
    Eg,EHq,Ephi,EPsi = E(G),E(Hq),E(phi),E(Psi)
    VG,Vphi = Var(G),Var(phi)
    Wgeo = EPsi - Vphi
    Gpp  = EHq - VG + Wgeo
    R_num= Gpp - 2*Wgeo
    R_ana= E(kap) - E(kapS) + Cov(phi,d) - Cov(G,d)
    return dict(Gpp=Gpp, W=Wgeo, R_num=R_num, R_ana=R_ana,
                Ed=E(d), Vd=Var(d), EG=Eg, VG=VG, Ephi=Ephi, Vphi=Vphi,
                identity_res=abs(R_num-R_ana))

def log(*a): print(*a, flush=True)

# ============================================================ 0. recover S2 at s=0
log("=== 0. CASE (a) / s=0 sanity:  should reproduce S2  Gpp = 2 W_geo,  R_proj=0 ===")
log(f"{'beta':>5} {'th':>6} {'v':>4} {'Gpp':>11} {'2W':>11} {'R_num':>10} {'R_ana':>10} {'|num-ana|':>9}")
for beta in (2.5,):
    for th in (0.0, np.pi):
        for nm,v in (('e3',(0,0,1)),('e1',(1,0,0))):
            o=compute(beta,th,v,0.0)
            log(f"{beta:5.2f} {th:6.3f} {nm:>4} {o['Gpp']:11.7f} {2*o['W']:11.7f} "
                f"{o['R_num']:10.3e} {o['R_ana']:10.3e} {o['identity_res']:9.1e}")

# ============================================================ 1. case (b): covariance scan
log("\n=== 1. CASE (b) covariant P:  w=Ad_{exp(-sX)}v .  R_proj = Gpp - 2W_geo ===")
log(f"{'beta':>5} {'th':>6} {'v':>4} {'s':>5} {'2W_geo':>10} {'R_num':>10} {'R_ana':>10}"
    f" {'Gpp':>10} {'sign R':>7} {'|num-ana|':>9}")
for beta in (0.5, 2.5, 5.0):
    for th in (0.0, np.pi/2, np.pi):
        for nm,v in (('e3',(0,0,1)),('e1',(1,0,0))):
            for s in (0.0, 0.3, 0.7, 1.2):
                o=compute(beta,th,v,s)
                sg = '+' if o['R_num']>=0 else '-'
                log(f"{beta:5.2f} {th:6.3f} {nm:>4} {s:5.2f} {2*o['W']:10.6f} "
                    f"{o['R_num']:10.6f} {o['R_ana']:10.6f} {o['Gpp']:10.6f} {sg:>7} "
                    f"{o['identity_res']:9.1e}")

# ============================================================ 2. ACID TEST (free Gaussian over R^3)
log("\n=== 2. ACID TEST: flat/free field must NOT manufacture a floor (Phi=0) ===")
def gauss_Gpp(a, vdir, s):
    """Gaussian free field mu ~ e^{-a|xi|^2/2} over R^3 (Phi=0), S=1/2 a|xi|^2.
       w=Ad_{exp(-sX)}v (s=0: constant P).  H=a|w|^2, G=<w,a xi>."""
    xi = XI.reshape(-1,3); v=np.array(vdir,float)
    R=np.linalg.norm(xi,axis=-1)
    w,_ = transport_w(xi, v, s)
    Hq = a*np.einsum('na,na->n', w, w)           # a|w|^2
    G  = a*np.einsum('na,na->n', w, xi)          # <w, a xi>
    dens=np.exp(-0.5*a*R*R)*R*R                   # e^{-a|xi|^2/2} * rho^2  (d^3xi = rho^2 drho du dpsi)
    wt=dens*(WRg.reshape(-1)*WUg.reshape(-1)*WPSg.reshape(-1))
    Z=wt.sum()
    EG=(wt*G).sum()/Z; EG2=(wt*G*G).sum()/Z; EHq=(wt*Hq).sum()/Z
    return EHq-(EG2-EG**2)
log(f"  {'a':>5} {'v':>4} {'s':>5} {'Gpp=E[H]-Var(G)':>16}")
for a in (2.0, 8.0):
    for nm,v in (('e3',(0,0,1)),('e1',(1,0,0))):
        for s in (0.0, 0.7):
            log(f"  {a:5.1f} {nm:>4} {s:5.2f} {gauss_Gpp(a,v,s):16.4e}")
log("  => Gpp=0 on the free Gaussian (abelian: w const, d=0,Dw=0; also non-abelian here,")
log("     since w·xi=v·xi for this axial transport) => (star) FAILS with delta_W=0. Required.")

# ============================================================ 2b. small-s: R_proj = O(s), s->0 recovers case (a)
log("\n=== 2b. small-s linearity: R_proj -> 0 as s->0  (case (a) = s=0) ===")
log(f"{'beta':>5} {'th':>6} {'v':>4} {'s':>6} {'R_proj':>12} {'R/s':>10} {'E[div w]':>11}")
for beta in (2.5,):
    for th in (np.pi,):
        for nm,v in (('e3',(0,0,1)),('e1',(1,0,0))):
            for s in (0.05,0.1,0.2,0.4):
                o=compute(beta,th,v,s)
                log(f"{beta:5.2f} {th:6.3f} {nm:>4} {s:6.2f} {o['R_num']:12.6f} "
                    f"{o['R_num']/s:10.5f} {o['Ed']:11.5f}")

# ============================================================ 3. size / sign of R vs head
log("\n=== 3. |R_proj| / (2W_geo)  (control): UV-suppressed or head-competitive? ===")
log(f"{'beta':>5} {'th':>6} {'v':>4} {'s':>5} {'2W_geo':>10} {'|R_proj|':>10} {'ratio':>10}")
for beta in (0.5, 2.5, 5.0):
    for th in (np.pi,):
        for nm,v in (('e3',(0,0,1)),('e1',(1,0,0))):
            for s in (0.3, 0.7, 1.2):
                o=compute(beta,th,v,s)
                denom=abs(2*o['W']) if abs(2*o['W'])>1e-12 else 1.0
                log(f"{beta:5.2f} {th:6.3f} {nm:>4} {s:5.2f} {2*o['W']:10.6f} "
                    f"{abs(o['R_num']):10.6f} {abs(o['R_num'])/denom:10.4f}")
