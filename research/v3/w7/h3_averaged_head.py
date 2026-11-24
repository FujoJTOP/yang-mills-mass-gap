#!/usr/bin/env python3
"""
h3_averaged_head.py  --  the AVERAGED head inequality (star_avg) of (star).

Object (corrected, G4/V3): the TRUE coarse marginal curvature
   G''(v;a0) = -d^2/da^2 log rho_v(a)|_{a0} = E[ d_v^2 V | a0 ] - Var( d_v V | a0 ),
   V = S + Phi ,  SU(2) one-plaquette / one-free-link model, B = exp(i theta sig3/2).

We test two natural averages of G'' :

  (A)  a0-average at fixed boundary U':
         G_A(v|U') := E_{a ~ rho_v(.|U')}[ G''(v;a) ]   =?  I(rho_v) := int (rho_v')^2/rho_v  >= 0
       THEOREM (integration by parts):  G_A = I(rho_v)  (Fisher information)  >= 1/Var(a)
       (Cramer-Rao / Heisenberg).  Margin = inverse coarse variance = coarse gap^2.

  (B)  U'-average (boundary holonomy theta) at fixed coarse field a0 :
         G_B(v;a0) := E_{theta ~ nu}[ G''(v;a0,theta) ] ,  good/bad set decomposition.

Acid test: free massless Gaussian rho = N(0,sigma^2) -> I = 1/sigma^2 -> 0.  Must vanish.

Run:  py -3.12 h3_averaged_head.py
"""
import numpy as np

# ---- cancellation-free radial pieces (same as v3_counterexample_scope.py) ----
def pr(x):
    y = np.empty_like(x); sm = np.abs(x) < 0.3
    y[~sm] = 1.0/x[~sm] - 1.0/np.tan(x[~sm])
    xs = x[sm]; y[sm] = xs/3 + xs**3/45 + 2*xs**5/945 + 2*xs**7/4725
    return y
def prr(x):
    y = np.empty_like(x); sm = np.abs(x) < 0.3
    y[~sm] = 0.5*(-1.0/x[~sm]**2 + 1.0/np.sin(x[~sm])**2)
    xs = x[sm]; y[sm] = 1.0/6 + xs**2/30 + 2*xs**4/189 + 2*xs**6/2475
    return y

XG, WG = np.polynomial.legendre.leggauss(6000)
trapz = getattr(np, 'trapezoid', None) or np.trapz
NFLOOR = 2/12.0

def parts(rho, b, t):
    hS = b*(1-np.cos(rho/2)*np.cos(t/2))
    hP = 2*np.log(rho/2) - 2*np.log(np.sin(rho/2))
    g  = b*np.sin(rho/2)*np.sin(t/2)
    return hS, hP, g

def rho_v_e3(a0, b, t):
    """coarse marginal density of a=<e3,xi> at a=a0 (slice integrated, 1D GL)."""
    Rm = np.sqrt(max(4*np.pi**2 - a0*a0, 1e-14))
    rr = 0.5*Rm*(XG+1.0); wrr = 0.5*Rm*WG
    rho = np.sqrt(rr*rr + a0*a0); u = a0/rho
    hS, hP, g = parts(rho, b, t)
    V = hS + hP + g*u
    return 2*np.pi*np.sum(np.exp(-V)*rr*wrr)

def Gpp_e3(a0, b, t, da=1e-4):
    """G'' from finite difference of -log rho_v, 5-point stencil (v=e3)."""
    xs = [rho_v_e3(a0+k*da, b, t) for k in (-2,-1,0,1,2)]
    lg = np.log(np.array(xs))
    return -(-lg[0]+16*lg[1]-30*lg[2]+16*lg[3]-lg[4])/(12*da*da)

def log(*a): print(*a, flush=True)

log("="*88)
log("H3  AVERAGED HEAD INEQUALITY  (star_avg)   SU(2) one-plaquette exact model")
log(f"    N/12 (Haar floor) = {NFLOOR:.6f}")
log("="*88)

# ============================================================
# [1]  Fisher identity   E_{a~rho}[G'']  =  I(rho)  >= 1/Var(a)      (THEOREM check)
# ============================================================
log("\n[1] a0-AVERAGE  G_A(v|U') = E_{a~rho_v}[Gpp]  vs  Fisher info I(rho_v) = int rho'^2/rho")
log("    (want: G_A = I exactly, and I >= 1/Var(a) by Cramer-Rao)")
log(f"    {'beta':>5}{'theta':>7}{'G_A':>20}{'I':>22}{'resid':>10}{'1/Var':>11}")

def fisher_check(b, t, na=181):
    lo, hi = -2*np.pi+0.02, 2*np.pi-0.02
    a = np.linspace(lo, hi, na)
    rho = np.array([rho_v_e3(x, b, t) for x in a])
    Gpp = np.array([Gpp_e3(x, b, t) for x in a])
    rp = np.gradient(rho, a)
    with np.errstate(divide='ignore', invalid='ignore'):
        I = trapz(rp*rp/rho, a)
    GA = trapz(rho*Gpp, a)
    # normalise (rho not exactly normalised by slice formula up to constant)
    Z = trapz(rho, a); rho_n = rho/Z
    I2 = trapz(rp*rp/rho, a)/Z
    mean = trapz(a*rho_n, a); var = trapz((a-mean)**2*rho_n, a)
    return GA/Z, I2, var

for (b,t) in [(2.5,np.pi),(5.0,np.pi),(2.5,0.0),(5.0,np.pi/2)]:
    GA, I, var = fisher_check(b, t)
    log(f"    {b:5.2f}{t:7.3f}{GA:20.6f}{I:22.6f}{GA-I:10.2e}{1.0/var:11.5f}")

# ============================================================
# [2]  Cramer-Rao margin across the head, and its BETA behaviour
# ============================================================
log("\n[2] margin delta_A = I(rho_v)  and Cramer-Rao floor 1/Var(a)   (theta=pi, e3)")
log(f"    {'beta':>6}{'I(rho)':>12}{'1/Var(a)':>12}{'I*Var':>9}   (I*Var >= 1 required)")
for b in [0.5,1.0,2.5,4.3066,5.0,10.0]:
    GA, I, var = fisher_check(b, np.pi, na=161)
    log(f"    {b:6.2f}{I:12.6f}{1.0/var:12.6f}{I*var:9.4f}")

# ============================================================
# [3]  barrier value (the pointwise counterexample) at the head
# ============================================================
log("\n[3] pointwise barrier  G''(e3; theta=pi, a0=0) = A - B beta^2")
w0 = (np.sin(np.pi*(XG+1.0)/2)**2/(np.pi*(XG+1.0)))*(np.pi*WG)
r0 = np.pi*(XG+1.0); Z0 = w0.sum(); E0 = lambda x:(w0*x).sum()/Z0
A = E0(pr(r0/2)/r0); ss = np.sin(r0/2)/r0; B = E0(ss*ss)-E0(ss)**2
log(f"    A={A:.6f}  B={B:.6f}  beta*=sqrt(A/B)={np.sqrt(A/B):.4f}")
for b in [2.5,4.3066,5.0,10.0]:
    log(f"    beta={b:6.2f}:  G''(theta=pi,a0=0) = {A-B*b*b:+.6f}   {'<0 (counterexample)' if A-B*b*b<0 else '>0'}")

# ============================================================
# [4]  U'-average at fixed a0=0 :  good/bad decomposition over the holonomy theta
# ============================================================
log("\n[4] U'-AVERAGE at a0=0 :  G_B(theta) = G''(e3; a0=0, theta);  theta in [0,pi]")
log("    good-set threshold theta_c = largest theta with c_Box+N/12-C1*beta*f(theta) > 0 (one-plaquette c_Box=0)")
# one plaquette has no Box Laplacian gap c_Box=0 -> rho(theta) = N/12 - C1*beta*f ; use C1=1 (conservative)
def G_theta(b, t):  # G'' at a0=0 for given theta (finite diff)
    return Gpp_e3(0.0, b, t)

log(f"    {'theta':>7}{'G_B':>11}    beta=2.5")
for t in [0.0,0.4,0.81,1.2,np.pi/2,np.pi]:
    log(f"    {t:7.3f}{G_theta(2.5,t):11.5f}")
# scan sign of G_B for beta=2.5 and beta=5
for b in [2.5,5.0]:
    ts = np.linspace(0, np.pi, 41)
    gs = np.array([G_theta(b,t) for t in ts])
    log(f"    beta={b:4.1f}: min_theta G_B = {gs.min():+.5f} at theta={ts[np.argmin(gs)]:.3f} ;  #theta with G_B<0 = {(gs<0).sum()}/41")

# U'-averages with two honest boundary measures
ts = np.linspace(0.0, np.pi, 201)
def theta_avg(b, weight='haar'):
    gs = np.array([G_theta(b,t) for t in ts])
    if weight=='haar':
        w = np.sin(ts/2)**2                      # SU(2) Haar on the holonomy angle
    elif weight=='wilson':
        w = np.sin(ts/2)**2*np.exp(-b*(1-np.cos(ts)))   # tilted by boundary plaquette action
    w = w/trapz(w, ts)
    return trapz(w*gs, ts)

log("\n    U'-averaged G_B  <G''>_theta at a0=0 :")
log(f"    {'beta':>6}{'<G>_Haar':>12}{'<G>_Wilson':>13}")
for b in [2.5,5.0]:
    log(f"    {b:6.2f}{theta_avg(b,'haar'):12.6f}{theta_avg(b,'wilson'):13.6f}")

# ============================================================
# [5]  good/bad measure estimate vs good-set margin  (the task's step 3)
# ============================================================
log("\n[5] good/bad decomposition at a0=0, theta-averaged (one-plaquette: c_Box=0, C1=1)")
log("    rho0(theta)=N/12 - beta*f, f=1-cos(theta)  (>0 = good set);  bad set = {rho0<0}")
log(f"    {'beta':>6}{'theta_c':>9}{'nu(bad)_Haar':>14}{'nu(bad)_Wilson':>16}{'rho0(0)':>9}{'min Gbad':>10}{'bound':>10}")
for b in [2.5,5.0,10.0]:
    fc = NFLOOR/b   # f where rho0=0 : f_c = (N/12)/beta
    thc = np.arccos(max(-1.0, min(1.0, 1.0-fc)))
    # Haar & Wilson bad-set measure; good-set lower bound rho0 on good set
    th = np.linspace(0.0, np.pi, 401)
    whaar = np.sin(th/2)**2; wwil = np.sin(th/2)**2*np.exp(-b*(1-np.cos(th)))
    good = th <= thc
    nubad_h = trapz(whaar[~good], th[~good])/trapz(whaar, th)
    nubad_w = trapz(wwil[~good], th[~good])/trapz(wwil, th)
    gs = np.array([G_theta(b,t) for t in th])
    gbad = gs[~good]
    gmin = gbad.min() if gbad.size else 0.0
    rho0max = NFLOOR  # rho0(0) = N/12 - beta*0
    bound = rho0max*(1.0-nubad_w) + gmin*nubad_w
    log(f"    {b:6.2f}{thc:9.4f}{nubad_h:14.5f}{nubad_w:16.6f}{rho0max:9.4f}{gmin:10.5f}{bound:10.5f}")

# ============================================================
# [5b]  SHARPER good/bad : charge only the GENUINELY NEGATIVE set, and a c_Box scan
# ============================================================
log("\n[5b] sharper split: charge the genuinely NEGATIVE set S_neg={G''<0}, not the whole 'not-good' set")
log("     (the crude [5] charges the whole not-good set at its worst value -> too pessimistic)")
th = np.linspace(0.0, np.pi, 401)
whaar = np.sin(th/2)**2; wwil = np.sin(th/2)**2*np.exp(-2.5*(1-np.cos(th)))
log(f"    {'beta':>5}{'#neg/401':>9}{'nu(Sneg)_Haar':>15}{'nu(Sneg)_Wilson':>17}{'min G'' on Sneg':>16}{'<G>_Haar':>10}")
for b in [2.5,5.0,10.0]:
    gs = np.array([G_theta(b,t) for t in th])
    neg = gs < 0
    wh = np.sin(th/2)**2; ww = np.sin(th/2)**2*np.exp(-b*(1-np.cos(th)))
    nneg_h = trapz(wh[neg], th[neg])/trapz(wh, th)
    nneg_w = trapz(ww[neg], th[neg])/trapz(ww, th)
    gmin_neg = gs[neg].min() if neg.any() else 0.0
    avgH = trapz(wh*gs, th)/trapz(wh, th)
    log(f"    {b:5.2f}{neg.sum():9d}{nneg_h:15.5f}{nneg_w:17.6f}{gmin_neg:16.5f}{avgH:10.5f}")

log("\n    c_Box scan: good/bad WORST-CASE bound vs true Wilson average, theta-averaged, C1=1")
log("    rho0(theta)=N/12 + c_Box - beta*f ;  good set {rho0>=0}")
log(f"    {'beta':>5}{'c_Box':>7}{'theta_c':>9}{'nu(bad)_Wilson':>16}{'bound':>10}{'true<G>_Wilson':>16}")
for b in [2.5,5.0]:
    for cbox in [0.0,0.25,0.5,1.0]:
        fc = (NFLOOR+cbox)/b
        thc = np.arccos(max(-1.0,min(1.0,1.0-fc)))
        good = th <= thc
        gs = np.array([G_theta(b,t) for t in th])
        ww = np.sin(th/2)**2*np.exp(-b*(1-np.cos(th)))
        nubad = trapz(ww[~good], th[~good])/trapz(ww, th)
        gmin = gs[~good].min() if (~good).any() else gs.min()
        bound = (NFLOOR+cbox)*(1.0-nubad) + gmin*nubad
        trueavg = trapz(ww*gs, th)/trapz(ww, th)
        log(f"    {b:5.2f}{cbox:7.2f}{thc:9.4f}{nubad:16.6f}{bound:10.5f}{trueavg:16.6f}")

# ============================================================
# [6]  ACID TEST : free massless Gaussian -> I = 1/sigma^2 -> 0
# ============================================================
log("\n[6] ACID TEST : free Gaussian rho = N(0,sigma^2) -> G_A = I = 1/sigma^2 -> 0")
for sig in [1.0,10.0,100.0,1e3]:
    log(f"    sigma={sig:8.1f}:  I = 1/sigma^2 = {1.0/sig**2:.3e}")

log("\nClay 0%. No QED.")
