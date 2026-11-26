#!/usr/bin/env python3
"""s06_counter_gpp.py -- S06: can (star)  G''_j = E[H_j] - Var_j + W_j >= delta_W > 0  be FALSIFIED?

Model: single plaquette SU(2) / one free link, boundary B = exp(i th sig3/2)  (S2 / G4 / H1A 2.1),
flat chart rho in (0,2pi), u = xi_3/rho :
    V(rho,u) = h(rho) + g(rho) u ,
    h = beta(1 - cos(rho/2) cos(th/2)) + Phi(rho) ,  g = beta sin(rho/2) sin(th/2) ,
    Phi = 2 log(rho/2) - 2 log sin(rho/2)          (Haar; ball chart, e^{-Phi} ~ sin^2(rho/2)).
Conditioned slice <v,xi> = a, v = e3:  measure e^{-h(rho)} r dr (2D plane), rho = sqrt(r^2+a^2), u = a/rho.
    G''_e3(a) = E[V_33 | a] - Var(V_3 | a) = -d_a^2 log rho_v(a)
at a=0 (analytic route):  V_33 = h'(r)/r = H_W(r) + Phi'(r)/r ,  V_3 = g(r)/r ,  E[V_3|0]=0
    => G''_e3(0) = E[H] - Var + W ,  E[H] := beta cos(th/2) E[sin(r/2)/(2r)] ,  W := E[Phi'(r)/r]
All quadrature deterministic Gauss-Legendre.  Run: py -3.12 s06_counter_gpp.py
"""
import sys
import numpy as np

OUT = []

def P(s=""):
    OUT.append(str(s))
    print(s, flush=True)

def Phi_r(x):
    y = np.empty_like(x); sm = np.abs(x) < 0.3
    xr = x[~sm]; y[~sm] = 1.0/xr - 1.0/np.tan(xr)
    xs = x[sm];  y[sm] = xs/3 + xs**3/45 + 2*xs**5/945 + 2*xs**7/4725
    return y

def Phi_rr(x):
    y = np.empty_like(x); sm = np.abs(x) < 0.3
    xr = x[~sm]; y[~sm] = 0.5*(-1.0/xr**2 + 1.0/np.sin(xr)**2)
    xs = x[sm];  y[sm] = 1.0/6 + xs**2/30 + 2*xs**4/189 + 2*xs**6/2475
    return y

# ---- e3 slice at a=0 : r in (0,2pi), measure e^{-h(r)} r dr ----
_xg, _wg = np.polynomial.legendre.leggauss(1200)
RG = np.pi*(_xg+1.0); WGRAM = np.pi*_wg

def slice_e3(beta, th):
    r = RG
    x = r/2.0
    C = np.cos(th/2); St = np.sin(th/2)
    h = beta*(1.0 - np.cos(x)*C) + (2*np.log(x) - 2*np.log(np.sin(x)))
    hp = beta*0.5*np.sin(x)*C + Phi_r(x)          # dh/dr
    g = beta*np.sin(x)*St
    w = np.exp(-h)*r*WGRAM
    Z = w.sum(); E = lambda f: (w*f).sum()/Z
    EH = E(beta*0.5*np.sin(x)*C/r)                # Wilson-Hessian part  E[H]
    W1 = E(Phi_r(x)/r)                            # W = E[Phi'(r)/r]
    W2 = E(Phi_rr(x))                             # W = E[Phi'']  (paper form)
    V3 = g/r
    Var = E(V3*V3) - E(V3)**2
    Gpp = E(hp/r) - Var
    return dict(EH=EH, W1=W1, W2=W2, Var=Var, Gpp=Gpp,
                EHmVar=EH-Var, recompose=EH-Var+W1)

# ---- whole curve: rho_v(a) = 2pi int e^{-h(rho)} r dr , rho=sqrt(r^2+a^2) ; FD on log ----
_xm, _wm = np.polynomial.legendre.leggauss(400)      # nodes computed ONCE

def marg(a, beta, th):
    Rmax = np.sqrt(max(4*np.pi**2 - a*a, 1e-12))
    r = 0.5*Rmax*(_xm+1.0); wr = 0.5*Rmax*_wm
    rho = np.sqrt(r*r + a*a); u = a/rho; x = rho/2.0
    h = beta*(1.0 - np.cos(x)*np.cos(th/2)) + (2*np.log(x) - 2*np.log(np.sin(x)))
    g = beta*np.sin(x)*np.sin(th/2)
    V = h + g*u
    return 2*np.pi*np.sum(np.exp(-V)*r*wr)

def Gpp_curve(a, beta, th, da=1e-2):
    v = np.array([np.log(marg(a+k*da, beta, th)) for k in (-2, -1, 0, 1, 2)])
    fpp = (-v[0]+16*v[1]-30*v[2]+16*v[3]-v[4])/(12*da*da)
    return -fpp

def minA(beta, th, nA=241, da=1e-2):
    ag = np.linspace(0.0, np.pi, nA)
    g = np.array([Gpp_curve(a, beta, th, da) for a in ag])
    return g.min(), ag[g.argmin()], g

# ---- direction j = e1 (transverse) at a=0 : 2D slice xi1=0, (s,alpha), s=sqrt(xi2^2+xi3^2) ----
_sg, _sw = np.polynomial.legendre.leggauss(900)
SG = np.pi*(_sg+1.0); WSG = np.pi*_sw
_ag, _aw = np.polynomial.legendre.leggauss(500)
AG = np.pi*(_ag+1.0); WAG = np.pi*_aw

def slice_e1(beta, th):
    s = SG[:, None]; al = AG[None, :]
    x = s/2.0
    h = beta*(1.0 - np.cos(x)*np.cos(th/2)) + (2*np.log(x) - 2*np.log(np.sin(x)))
    hp = beta*0.5*np.sin(x)*np.cos(th/2) + Phi_r(x)
    g = beta*np.sin(x)*np.sin(th/2)
    gp = beta*0.5*np.cos(x)*np.sin(th/2)                 # dg/ds
    d_g_over_s = gp/s - g/s**2
    V11 = hp/s + d_g_over_s*np.sin(al)
    # dV/dxi1 = h'*xi1/rho + g'u*xi1/rho + g*du/dxi1 with du/dxi1 = -xi3 xi1/rho^3 => 0 on the slice:
    # V_1 == 0 pointwise, so Var(V_1|0) = 0 exactly (no variance term).
    w = np.exp(-(h + g*np.sin(al)))*s*WSG[:, None]*WAG[None, :]
    Z = w.sum()
    return (w*V11).sum()/Z

def main():
    P("#"*100)
    P("# s06_counter_gpp -- counterexample hunt for (star):  G''_j = E[H_j]-Var_j+W_j >= delta_W > 0")
    P("# model: single plaquette SU(2), theta=pi barrier, v=e3 (longitudinal), flat chart rho<2pi")
    P(f"# python {sys.version.split()[0]}  numpy {np.__version__}  GL n=1200 (r) / 900x500 (e1)")
    P("#"*100)

    ref = {0.0: 0.205115297123, 2.5: 0.135994711042, 6.0: -0.193019, }
    P("\n[A] THREE BETAS: E[H], Var, W, G'' at a=0 (theta=pi, v=e3)   [closed form: G''=A-B beta^2]")
    P(f"  {'beta':>5} {'E[H]':>13} {'Var':>13} {'W=E[Phi/r]':>13} {'E[Phi'']':>13} "
      f"{'E[H]-Var':>12} {'G''':>14} {'ref':>12} {'d':>9}")
    rows = {}
    for b in (0.0, 1.0, 2.5, 4.30661, 5.0, 6.0, 10.0):
        q = slice_e3(max(b, 0.0), np.pi)
        q["beta"] = b; rows[b] = q
        r = ref.get(b, float("nan"))
        dd = q["Gpp"] - r if not np.isnan(r) else float("nan")
        P(f"  {b:5.2f} {q['EH']:+13.7f} {q['Var']:13.7f} {q['W1']:13.7f} {q['W2']:13.7f} "
          f"{q['EHmVar']:+12.7f} {q['Gpp']:+14.9f} {r:+12.6f} {dd:+9.1e}")
    A = rows[0.0]["W1"]                                   # E_Haar[Phi'(r)/r] = W_0 (paper anchor)
    B = rows[1.0]["Var"]                                  # beta=1 : Var = Var_Haar(sin(r/2)/r) = B
    P(f"  Haar anchors:  A = E[Phi'(r)/r] = {A:.12f}   (= paper W_0 = 1/(2Cin2pi) = 0.20511530)   "
      f"B = Var_Haar(sin(r/2)/r) = {B:.12f}")
    P(f"  note:  W := E[Phi'(r)/r] (the beta-independent piece of E[V_33]);  "
      f"E[Phi''] = {rows[0.0]['W2']:.9f} is a DIFFERENT number (not the W of (star)).")
    P(f"  closed form A - B beta^2 :  beta*=sqrt(A/B) = {np.sqrt(A/B):.6f}   "
      f"(archive 4.30661);  check vs quadrature: "
      f"max|G''-(A-B beta^2)| = {max(abs(rows[b]['Gpp']-(A-B*b*b)) for b in rows):.2e}")
    P(f"  => G''(beta=6, a=0) = {A-B*36:+.9f}   < 0 ?  {A-B*36 < 0}   (threshold beta*=4.30661)")

    P("\n[B] INDEPENDENT ROUTE to the same a=0 number:  marginal density + 5-pt FD of log rho")
    P(f"  {'beta':>5} {'G''(FD e3,a=0)':>16} {'G''(analytic)':>15} {'diff':>10}")
    for b in (0.0, 2.5, 6.0):
        gf = Gpp_curve(0.0, max(b, 1e-6), np.pi, da=1e-2)
        P(f"  {b:5.2f} {gf:+16.9f} {rows[b]['Gpp']:+15.9f} {gf-rows[b]['Gpp']:+10.1e}")

    P("\n[C] THE DECISIVE CHECK: whole curve G''(a), beta=2.5 (INSIDE the campaign window), theta=pi")
    g, a_, curve = minA(2.5, np.pi, nA=241)
    P(f"  min_a G''(a) = {g:+.9f} at a = {a_:.5f}   (archive o28: -0.079195 @ a=1.822)")
    ag = np.linspace(0.0, np.pi, 241)
    neg = ag[curve < 0]
    P(f"  #grid points with G''<0 : {len(neg)} / {len(ag)} ;  a-range of negativity: "
      f"[{neg.min():.5f}, {neg.max():.5f}]   (archive: [0.7147, 2.9688])")
    P("  a-sample: " + "  ".join(f"a={a:.2f}:{c:+.4f}" for a, c in
                                 zip(ag[::24], curve[::24])))

    P("\n[D] WHERE DOES IT TURN NEGATIVE?  beta_cv = inf{beta : min_a G''(a) < 0}  (theta=pi, e3)")
    lo, hi = 1.0, 3.0
    for _ in range(24):
        mid = 0.5*(lo+hi)
        if minA(mid, np.pi, nA=121)[0] < 0: hi = mid
        else: lo = mid
    P(f"  beta_cv = {0.5*(lo+hi):.6f}   (archive o02_bscan: 1.873291)   "
      f"=> for beta in ({0.5*(lo+hi):.3f}, 2.5] the curve G''(a) takes NEGATIVE values")

    P("\n[E] IS IT j-DEPENDENT?  transverse direction v=e1 at a=0 (theta=pi)")
    P(f"  {'beta':>5} {'G''(e3)':>13} {'G''(e1)':>13} {'sign(e3)':>9} {'sign(e1)':>9}")
    for b in (0.0, 2.5, 5.0, 6.0, 10.0):
        g1 = slice_e1(max(b, 1e-9), np.pi)
        g3 = rows.get(b, slice_e3(max(b, 1e-9), np.pi))["Gpp"] if b in rows else slice_e3(b, np.pi)["Gpp"]
        P(f"  {b:5.2f} {g3:+13.6f} {g1:+13.6f} {str(g3>0):>9} {str(g1>0):>9}")

    P("\n[F] theta-dependence at a=0, beta=2.5 (is the barrier theta=pi special?)")
    for th in (0.0, np.pi/2, np.pi):
        q = slice_e3(2.5, th)
        P(f"  theta={th:5.3f}: E[H]={q['EH']:+9.6f} Var={q['Var']:9.6f} W={q['W1']:9.6f} "
          f"G''={q['Gpp']:+9.6f}")
        P(f"      min_a G''(a) at this theta = {minA(2.5, th, nA=121)[0]:+.6f}")

    P("\n[VERDICT] counterexample EXISTS: (star) as stated (uniform in j,U',a, any beta) is FALSE.")
    P("  witness 1: beta=6, theta=pi, v=e3, a=0 :  G'' = -0.193019  (exact closed form A-B beta^2)")
    P("  witness 2: beta=2.5 (in-window!), theta=pi, v=e3, a=1.822 : G'' = -0.079195")
    P("  witness 3: beta=2.5, any beta>1.8733 : the a-curve dips below 0 => delta_W(a-independent) impossible")
    P("  positive:  beta=0 -> +0.205115 ; beta=2.5, a=0 -> +0.135995 ; theta=0 -> +0.729924")
    P("  mechanism: G'' = W - Var ,  W = 0.205115 (Haar anchor, beta-indep.) ,  Var = B beta^2 -> oo")
    with open("w7/s06_counter_gpp.log", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

main()
