#!/usr/bin/env python3
"""
transport_packing.py  --  Audit of the "dense instanton gas" premise in
research/unconditional/E2-FLIP-TRANSPORT.md (claim: d := n rho^4 ~ O(10) for SU(3)).

Goal: compute, from scratch, the one-loop 't Hooft/CDG instanton size distribution,
the peak, the PACKING FRACTION d = n * rho_bar^4 (n = INTEGRAL of D, not peak value),
and cross-check against ILM / lattice.  Identify exactly where an O(10) could come from.

Units: Lambda = 1 throughout for the one-loop part; conversion 1 fm^-1 = 197.327 MeV.

Run:  py -3.12 transport_packing.py     (single job, numpy/scipy only)
"""
import numpy as np
from math import factorial, exp, log, gamma, pi
from scipy.integrate import quad

FM_INV_MEV = 197.3269804            # 1 fm^-1 in MeV
def fm4_from_MeV4(m): return (m / FM_INV_MEV)**4
def MeV_from_fm_inv(x): return x * FM_INV_MEV

# ---------------------------------------------------------------------------
# 1. The exact one-loop SU(N) instanton size distribution (Schaefer-Shuryak Eq 2.13)
#      D(rho) drho = C_N * (8 pi^2/g^2(rho))^{2N} * exp(-8 pi^2/g^2(rho)) * drho/rho^5
#      C_N = 0.466 * exp(-1.679 N) / [ (N-1)! (N-2)! ]
#      one-loop: 8 pi^2/g^2(rho) = S(rho) = b ln(1/(rho Lambda)),  b = 11N/3  (pure YM)
# ---------------------------------------------------------------------------
def b_of(N):  return 11.0*N/3.0
def C_of(N):  return 0.466*exp(-1.679*N)/(factorial(N-1)*factorial(N-2))
def S_of(rho, N, Lambda=1.0):     # 8 pi^2 / g^2(rho)
    x = rho*Lambda
    if x <= 0: return np.inf
    return b_of(N)*log(1.0/x)
def D_of(rho, N, Lambda=1.0):
    S = S_of(rho, N, Lambda)
    if not np.isfinite(S) or S <= 0: return 0.0
    return C_of(N)*S**(2*N)*exp(-S)/(rho**5)

print("="*78)
print(" 1. ONE-LOOP PREFACTOR C_N   (the note assumed C = O(0.1))")
print("="*78)
for N in (2,3,4,5,6):
    print(f"   SU({N}): b={b_of(N):6.3f}  C_N={C_of(N):.5e}   "
          f"[(N-1)!(N-2)! = {factorial(N-1)*factorial(N-2)}]")
print("   -> for SU(3): C_3 = %.5e   (NOT O(0.1); ~%.0fx smaller)"%(C_of(3), 0.1/C_of(3)))

# ---------------------------------------------------------------------------
# 2. Peak of D(rho)
#    d/drho [S^{2N} e^{-S} rho^{-5}] = 0  =>  S* = 2 N b/(b-5)
# ---------------------------------------------------------------------------
print("\n"+"="*78)
print(" 2. PEAK OF THE ONE-LOOP DISTRIBUTION")
print("="*78)
for N in (2,3):
    b=b_of(N); Spk=2*N*b/(b-5); ln_inv=Spk/b; rpk=exp(-ln_inv)
    G   = Spk**(2*N)*exp(-Spk)          # pure integrand factor S*^{2N} e^{-S*}  (= 29.6 for SU(3))
    Dpk = C_of(N)*G/rpk**5              # D(rho*)
    print(f"   SU({N}): rho*={rpk:.4f}/L  S*=8pi^2/g^2={Spk:.3f}  e^-S*={exp(-Spk):.3e}")
    print(f"          G := S*^(2N) e^-S*        = {G:.4f}   (the note's '(11)^6 e^-11 = 29.6')")
    print(f"          D(rho*) rho*^5 = C_N * G  = {Dpk*rpk**5:.4e}  (real C_N; = {G:.2f}*C_N)")
    print(f"          D(rho*) rho*^4            = {Dpk*rpk**4:.4e}  [dim 1/length -- NOT a packing frac]")

# ---------------------------------------------------------------------------
# 3. THE ACTUAL PACKING FRACTION:  n = INTEGRAL of D(rho) drho ; d = n * rho_bar^4
#    analytic (exact):  n = C_N (2N)! / (b-4)^{2N+1}
#                       rho_bar = (b-4)^{2N+1}/(b-3)^{2N+1}
# ---------------------------------------------------------------------------
print("\n"+"="*78)
print(" 3. FULL-INTEGRAL DENSITY AND PACKING FRACTION  (this is the physical d)")
print("="*78)
def compute_d(N, rho_max=1.0, Lambda=1.0, with_S_factor=True):
    b=b_of(N)
    f=lambda r: D_of(r,N,Lambda) if with_S_factor else C_of(N)*exp(-S_of(r,N,Lambda))/r**5
    n,_=quad(f,1e-9,rho_max,limit=400,epsabs=1e-13,epsrel=1e-10)
    g=lambda r: r*f(r)
    m1,_=quad(g,1e-9,rho_max,limit=400,epsabs=1e-13,epsrel=1e-10)
    rbar=m1/n
    return n, rbar, n*rbar**4

for N in (2,3):
    b=b_of(N)
    n_an = C_of(N)*b**(2*N)*factorial(2*N)/(b-4)**(2*N+1)
    rbar_an = (b-4)**(2*N+1)/(b-3)**(2*N+1)
    n_num, rbar_num, d_num = compute_d(N)
    print(f"   SU({N}):  analytic n = C_N (2N)!/(b-4)^(2N+1) = {n_an:.5e} /L^4  "
          f"(numeric {n_num:.5e})")
    print(f"          analytic rho_bar = ((b-4)/(b-3))^(2N+1) = {rbar_an:.4f}/L   "
          f"(numeric {rbar_num:.4f})")
    print(f"          *** d = n rho_bar^4 = {d_num:.4f}   (analytic {n_an*rbar_an**4:.4f}) ***")
    print(f"          same with C_N replaced by note's C=0.1 :  d = {d_num/C_of(N)*0.1:.3f}")
    print(f"          same with C=0.3                         :  d = {d_num/C_of(N)*0.3:.3f}")

# ---------------------------------------------------------------------------
# 4. Where did O(10) come from?  Reproduce the note's arithmetic EXACTLY.
# ---------------------------------------------------------------------------
print("\n"+"="*78)
print(" 4. REPRODUCING THE NOTE'S O(10)  (diagnosis of the error)")
print("="*78)
N=3; b=b_of(N); Spk=2*N*b/(b-5); rpk=exp(-Spk/b)
factor = Spk**(2*N)*exp(-Spk)     # = (11)^6 e^-11
print(f"   SU(3): (11)^6 e^-11 = {factor:.3f}   (the note's stated 29.6)")
print(f"   note wrote 'n rho*^4 = C * 29.6' with C = O(0.1):")
for C in (0.1,0.2,0.3):
    print(f"        C={C}: C*29.6 = {C*29.6:.1f}   <- this is the 'O(10)'")
print(f"   BUT: (i) 29.6 is  D(rho*) rho*^5 = C_N S^(2N) e^-S*  -- that is the density PER UNIT SIZE")
print(f"          times rho*^5; dimension = 1/length ... it is NOT n*rho^4 (packing fraction).")
print(f"        (ii) the ONE-LOOP coefficient is C_3 = {C_of(3):.4e}, not O(0.1).")
print(f"   Correct peak-size packing fraction (using the true integral and true C_3):")
print(f"        d = {compute_d(3)[2]:.4f}   (dilute/borderline, NOT O(10))")

# ---------------------------------------------------------------------------
# 5. HONEST RANGE: truncate the (invalid-at-large-rho) one-loop at rho_max, vary C_N
# ---------------------------------------------------------------------------
print("\n"+"="*78)
print(" 5. HONEST RANGE OF d FROM THE ONE-LOOP FORMULA (rho_max sweep; C_N sensitivity)")
print("="*78)
for N in (3,):
    print(f"   SU({N}) d vs IR cutoff rho_max/L:")
    for rmax in (0.3,0.5,0.7,1.0):
        n,rb,d = compute_d(N,rho_max=rmax)
        print(f"        rho_max={rmax:.2f}: n={n:.3e}  rho_bar={rb:.3f}  d={d:.4f}")
    print("   C_N multiplier sweep at rho_max=1 (d scales linearly with C_N):")
    d1=compute_d(N)[2]
    for C in (C_of(3), 0.01,0.03,0.1,0.3):
        print(f"        C={C:.4e}: d={d1/C_of(3)*C:.4f}")

# ---------------------------------------------------------------------------
# 6. ILM + LATTICE CROSS-CHECK  (empirical, this is what sets the real d)
# ---------------------------------------------------------------------------
print("\n"+"="*78)
print(" 6. ILM AND LATTICE CROSS-CHECK (empirical d)")
print("="*78)
print("   Schaefer-Shuryak ILM (RMP 70 (1998) 323): n ~ 1 fm^-4, rho ~ 0.33 fm")
for (n_fm4, rho_fm, tag) in [(1.0,0.33,"SS central"),(1.0,0.40,"larger rho"),
                             (1.0,0.50,"cooling-based"),(0.7,0.33,"low chi"),(1.3,0.35,"high chi")]:
    d=n_fm4*rho_fm**4
    volfrac=n_fm4*(pi**2/2)*rho_fm**4     # 4-ball volume = pi^2/2 rho^4
    print(f"        n={n_fm4:.2f} fm^-4, rho={rho_fm:.2f} fm : d=n rho^4 = {d:.4f}   "
          f"(4D vol fraction {volfrac:.3f})   [{tag}]")
print("   Review's own quote: pi^2 * n * rho^4 ~ 0.12  =>  n rho^4 ~ 0.12/pi^2 = %.4f"%(0.12/pi**2))
print("   Lattice chi_top (pure SU(3)) = n (dilute-gas relation chi_top = n):")
for M in (179,187,195,200):
    chi=fm4_from_MeV4(M)
    print(f"        chi_top^1/4 = {M} MeV  ->  chi_top = {chi:.3f} fm^-4  (=n) ; "
          f"d(n=chi,rho=0.33) = {chi*0.33**4:.4f} ; d(.,rho=0.5) = {chi*0.5**4:.4f}")

# ---------------------------------------------------------------------------
# 7. FREE-FIELD ACID TEST: remove running -> scale invariant -> d undefined
# ---------------------------------------------------------------------------
print("\n"+"="*78)
print(" 7. FREE-FIELD (SCALE-INVARIANT) ACID TEST")
print("="*78)
N=3; C=C_of(3); Sfix=8*pi**2/1.0     # fixed coupling, no running
print(f"   no running: S=8pi^2/g^2 fixed; D(rho) = C S^(2N) e^-S / rho^5  ~ 1/rho^5")
print(f"   -> D(lambda rho) = lambda^-5 D(rho): scale invariant; NO preferred size.")
for lo,hi in [(1e-6,1.0),(1e-3,1.0),(1e-2,1.0)]:
    # integral of 1/rho^5 diverges at rho->0 ; truncated it blows up as lo->0
    val,_=quad(lambda r: C*Sfix**(2*N)*exp(-Sfix)/r**5, lo, hi, limit=300)
    print(f"        int_{lo}^{hi} D drho = {val:.3e}  (diverges as lo->0: no finite n)")
print("   => n undefined, rho_bar undefined, d = n rho_bar^4 UNDEFINED.")
print("      The diluteness MASS statement vanishes in the scale-invariant limit. ACID PASSED.")

# ---------------------------------------------------------------------------
# 8. KP criterion sanity: what cluster expansion actually needs
# ---------------------------------------------------------------------------
print("\n"+"="*78)
print(" 8. WHAT KOTECKY-PREISS ACTUALLY NEEDS (vs n rho^4)")
print("="*78)
print("   KP: sum_p |z_p| e^{a|p|} < infinity  -- activity z_p ~ e^{-S_p} and geometric")
print("   size |p| ~ rho of the defect.  n rho^4 is only a proxy for 'dilute'.")
print("   e^{-S} at the physical peak (SU(3), S*=11) = %.3e (small activity) ;" % exp(-11))
print("   but the support/range rho ~ 1/Lambda is LARGE and mutual overlap O(1):")
print("   KP convergence is governed by ACTIVITY x RANGE, not by n rho^4 alone.")
print("   => d ~ 0.01-0.05 is CONSISTENT with dilute, but does NOT by itself prove KP;")
print("      continuum control of the activity at rho ~ 1/Lambda is the open step.")

print("\n"+"="*78)
print(" SUMMARY")
print("="*78)
print(f"  one-loop C_3                      = {C_of(3):.4e}   (note assumed O(0.1))")
print(f"  one-loop full-integral d          = {compute_d(3)[2]:.4f}  (borderline)")
print(f"  ILM d (n=1 fm^-4, rho=0.33 fm)    = {1.0*0.33**4:.4f}  (dilute)")
print(f"  lattice-consistent d range        = 0.01 - 0.06 (rho 0.33-0.5 fm)")
print(f"  note's O(10) reproduced only with C=O(0.1) x D(rho*)rho*^5 -> ARTIFACT")
