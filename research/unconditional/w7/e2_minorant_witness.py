"""
e2_minorant_witness.py  --  [V-here] witnesses for E2-FLIP-MINORANT.md
Load-bearing DIRECTION claims of the five lower-bound methods, checked on toy spectra.
Pure numpy. No heavy lattice job. Run:  py -3.12 e2_minorant_witness.py
"""
import numpy as np
rng = np.random.default_rng(0)

def hr(t): print("\n" + "="*72 + "\n" + t + "\n" + "="*72)

# ---------------------------------------------------------------- (A) Rayleigh = UPPER
hr("(A) Rayleigh quotient R(psi) is an UPPER bound on m1 = E1-E0  (so RQ is NOT a minorant)")
E = np.array([0.0, 1.0, 2.0])           # vac, first excited (m1=1), second
m1 = E[1] - E[0]
rs = []
for _ in range(200000):                  # trials on (vac)^perp : c1 e1 + c2 e2
    c = rng.normal(size=2); c /= np.linalg.norm(c)
    rs.append(E[1]*c[0]**2 + E[2]*c[1]**2)
rs = np.array(rs)
print(f"   200k random trials on (vac)^perp:  min R = {rs.min():.6f}  (exact inf = E1 = {E[1]:.1f})")
print(f"   all R >= E1 ?  {bool(np.all(rs >= E[1]-1e-12))}   -> RQ bounds m1 FROM ABOVE (Phi>=m1, violates Phi<=m1)")
print(f"   inf over trials = E1 only if psi = exact eigenvector (min-max); not certifiable => no lower bound from RQ.")

# ---------------------------------------------------------------- (B) Temple needs rho <= E2
hr("(B) Temple LOWER bound  E0 >= R - Var/(rho-R)  is sound iff rho <= E2; rho > E2 OVERCLAIMS")
E0, E1, E2 = 0.0, 1.0, 1.5
print(f"   true spectrum E0={E0}, E1={E1}, E2={E2}   (target lower bound on E1=m1)")
for ad in (0.01, 0.05, 0.2, 0.5):
    psi = np.array([0.0, np.sqrt(1-ad**2), ad])          # mostly phi1 + phi2 admixture
    R = float(np.sum(np.array([E0,E1,E2]) * psi**2))
    Var = float(np.sum((np.array([E0,E1,E2]) - R)**2 * psi**2))
    line = f"   admix={ad:4.2f}  R={R:7.4f}  Var={Var:7.5f}"
    for rho, tag in ((E2,"rho=E2  (valid)"), (2*E1,"rho=2E1 (>E2)"), (3.0,"rho=3.0 (>E2)")):
        if R < rho:
            b = R - Var/(rho - R)
            flag = "sound(<=E1)" if (rho <= E2+1e-12) else ("OK" if b <= E1+1e-12 else "**OVERCLAIM b>E1**")
            line += f" | {tag}: b={b:8.4f} {flag}"
    print(line)
print("   => the only sound Temple bound uses rho<=E2. A trial-only 'rho' is either E2 (circular) or unsound.")

# ---------------------------------------------------------------- (C) trace/partition minorant vacuous
hr("(C) partition/trace minorant  m1 >= -(1/t) ln Z_perp(t),  Z_perp=sum_{n>=1} e^{-t E_n}")
print("   single-state (only E1, i.e. E2=+inf): Z_perp=e^{-t E1} => bound=-(1/t)ln e^{-tE1}=E1  EXACT")
print("   ...but that *is* knowing N(E1)=1, i.e. E2, i.e. the gap. With any continuum above E1:")
m1edge = 1.0
Es_cont = np.concatenate([[m1edge], m1edge + np.linspace(0.0, 5.0, 500)**1.4])   # dense continuum
Es_free = 0.0 + np.linspace(0.0, 6.0, 800)                                        # gapless (massless) toy
for tag, Es in (("edge=1 + continuum", Es_cont), ("massless (edge 0)", Es_free)):
    for t in (0.5, 1.0, 2.0, 5.0, 10.0):
        Z = float(np.sum(np.exp(-t*Es[1:] if tag.startswith('edge') else -t*Es)))
        Z = float(np.sum(np.exp(-t*(Es[1:] if tag.startswith('edge') else Es))))
        b = -(1.0/t)*np.log(Z) if Z > 0 else float('inf')
        print(f"   [{tag:18s}] t={t:5.1f}: Z_perp={Z:10.3e}  bound={b:8.4f}  {'VACUOUS(<=0)' if b<=0 else 'positive'}")
print("   => with a continuum the bound is <=0 for all accessible t; positivity needs N(E1)=1 = E2.")

# ---------------------------------------------------------------- (D) genuine strong-coupling minorant
hr("(D) strong-coupling: u(beta)=c_f(beta) in (0,1) => m1 >= 4*ln(1/u)/a > 0  (a REAL [C] minorant, beta<beta1)  〔factor-4 corr. 2026-11-24〕")
# SU(2) plaquette weight e^{(beta/2) Re tr U} = e^{beta cos th};  Haar dU=(2/pi) sin^2 th dth
def c_f_su2(beta, Nth=6000):
    th = np.linspace(1e-6, np.pi-1e-6, Nth)
    w  = (2.0/np.pi)*np.sin(th)**2
    chi_f = 2.0*np.cos(th)                       # j=1/2 character sin(2th)/sin(th)
    integ = np.exp(beta*np.cos(th))*chi_f*w       # = e^{S} * chi_f under Haar
    return float(np.trapezoid(integ, th))
def Z_su2(beta, Nth=6000):
    th = np.linspace(1e-6, np.pi-1e-6, Nth); w=(2.0/np.pi)*np.sin(th)**2
    return float(np.trapezoid(np.exp(beta*np.cos(th))*w, th))
for beta in (0.005, 0.01, 0.02, 0.05, 0.10, 0.20):
    u = c_f_su2(beta)
    m1_lead = 4.0*np.log(1.0/u)   # factor-4 corr. 2026-11-24: P3 flux-tube mass m1 = 4*ln(1/u)/a (lattice a=1)
    print(f"   beta={beta:5.3f}:  u=c_f={u: .6e}  (0<u<1: {0<u<1})   ln(1/u)={np.log(1.0/u):7.4f}>0  => m1=4*ln(1/u)={m1_lead:7.4f}>0 (lattice units, a=1)")
print("   => on beta<beta1 the KP/OS cluster expansion gives an EXPLICIT uniform-in-L lower bound m1>=m0(beta)/a,")
print("      whose leading rate is 4*ln(1/u)/a (>0). This is the one genuine [C] minorant on F_sc;")
print("      it does NOT extend to beta~2.5 (Savvidy: no stable quadratic minorant, MS-BLIND 1.5).")

hr("SUMMARY")
print(" (A) Rayleigh  = UPPER   (B) Temple valid iff rho<=E2   (C) trace minorant vacuous w/ continuum")
print(" (D) strong-coupling minorant is REAL but restricted to beta<beta1.  No uniform-across-crossover Phi constructed.")
