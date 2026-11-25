#!/usr/bin/env python3
r"""
uc_beta.py -- UC-COMPLETE witness: the beta-dependence of the SU(2) one-plaquette
activity u(beta), and the exact ceiling beta_1 of every "u-based" convergence route
to charged completeness (u-completeness).

Deterministic, numpy + scipy.special.iv + scipy.optimize.brentq, single process, ~2 s.
Run:  py -3.12 -I uc_beta.py
Out:  uc_beta.log

Maps 1:1 onto research/v3/species/UC-COMPLETE.md:

  [A] exact SU(2) character coefficients   c_n(beta) = I_n(beta) - I_{n+2}(beta),
      n = #boxes = 2j, dim = n+1;  identity  SUM_n (n+1) c_n = e^beta ;  c_n > 0.
      Four exact forms of the activity, verified to agree:
         u(beta) = c_1/c_0 = 2 I_2(beta)/I_1(beta)
                 = 2 d/dbeta ln c_0(beta) = 2 ( I_0/I_1 - 2/beta )
                 = 2 / (4/beta + 1/(6/beta + 1/(8/beta + ...)))     [Gauss J-fraction]
      The J-fraction truncations alternate => two-sided RIGOROUS bounds on u(beta),
      used to bracket beta_tube without trusting scipy.

  [B] EXACT DICTIONARY between the two normalisations used in the corpus:
         T1-POLYMER-RG / P3 / BN-RESIDUAL:   u_P3  = c_1/c_0
         T5-CENTER-VORTEX / T23-TEMPLATE:    u_T5  = I_2/I_1
      From the recurrence I_{n-1}-I_{n+1} = 2n I_n/beta:   u_P3 = 2 u_T5  EXACTLY.
      => T23's "6u = 3.04" and T1's "u(2.5) = 1.01439" are the same sentence.

  [C] MONOTONICITY [T]: u = 2 (ln c_0)'  =>  u'(beta) = 2 Var_beta(cos theta) >= 0,
      strictly > 0 for beta>0.  Hence {beta : u(beta)<u_c} = [0,beta_max(u_c)) is an
      interval: "the largest beta with u<u_c" is well defined and unique.
      Verified: finite difference of u  vs  closed form  1 - 3 I_2/(beta I_1) - (I_2/I_1)^2.

  [D] The four numbers asked for:  u(1/48), u(0.5), u(1), u(2.5)  (both conventions).

  [D2] The KP sum is MONOTONE in the grading a (all c_n>0), so the sharpest member of
      the class sits at a -> 0^+, where it collapses to the single closed form
         Phi(0,beta) = (e^beta - c_0)/c_0 = beta e^beta / (2 I_1(beta)) - 1 .
      Evaluating at a=1 (as T1-POLYMER-RG does) is conservative by Phi(1)/Phi(0).

  [E] THE CEILING THEOREM (the T0(a) answer), made exact:
      every criterion of the class  rho(a,beta) = lambda_a * SUM_{n>=1}(n+1)(c_n/c_0)e^{an} < 1
      (any grading a>0, any animal constant lambda_a>=1) satisfies, by its n=1 term,
            lambda_a * 2 u(beta) e^a  <  1   =>   u(beta) < 1/(2 lambda_a e^a) <= 1/2 ,
      hence its beta-set lies in [0, u^{-1}(1/2)) = [0, beta_KP).
      => beta=2.5 is outside the whole KP class UNCONDITIONALLY (no constant reaches it).
      Hierarchy of exact roots: beta_tube = u^{-1}(1), beta_KP = u^{-1}(1/2),
      beta_f0 (= beta of Phi(0)=1, i.e. lambda_a=1), beta_KP(lambda_a=e*Delta).
      Scan over lambda_a = 1 ... 1e3 shows the root moves by <2x, all far below 2.5.
      The lambda_a that beta=2.5 WOULD need: 1/Phi(0,2.5) = 0.198 < 1 = structural floor.

  [F] Why no repair of "u alone" can work at beta=2.5: the n-profile of the Ursell sum.
      At beta=2.5 the n>=2 terms carry 60% of Phi(0) (90% of Phi(1)); the plaquette
      measure is dense in ALL representations.  So the failure is not "u too big".

  [G] Cross-checks of the corpus: u/(1-u) = -70.49 (BN-RESIDUAL N3), 6 u_T5 = 3.04 (T23),
      m_1 a = 4 ln(1/u) (four-wall tube), free-field acid test u->2.

Nothing here is a proof of the .md's theorems; the geometry (four-wall minimality), the
KP criterion, monotone convergence of the J-fraction and Var>=0 are cited in the .md.
This witness is [V].
"""

import numpy as np
from scipy.special import iv
from scipy.optimize import brentq

# ---------------------------------------------------------------- exact SU(2) data
def Ib(n, b):
    return float(iv(int(n), b))

def c(n, b):
    """exact SU(2) character coefficient; n = #boxes = 2j, dim = n+1"""
    return Ib(n, b) - Ib(n + 2, b)

def c0(b):
    return c(0, b)

def u_p3(b):
    """T1-POLYMER-RG / P3 / BN-RESIDUAL normalisation  u = c_1/c_0"""
    return c(1, b) / c0(b)

def u_t5(b):
    """T5-CENTER-VORTEX / T23-TEMPLATE normalisation  u = I_2/I_1"""
    return Ib(2, b) / Ib(1, b)

def var_cos(b):
    """Var_beta(cos theta) under the plaquette measure e^{b cos} dU / Z, closed form"""
    r = Ib(2, b) / Ib(1, b)
    return 1.0 - 3.0 * r / b - r * r

def Phi(a, b, nmax=300):
    """per-plaquette KP activity at grading a: SUM_{n>=1} (n+1)(c_n/c_0) e^{a n}"""
    s = 0.0
    for n in range(1, nmax + 1):
        t = (c(n, b) / c0(b)) * (n + 1) * np.exp(a * n)
        s += t
        if n > 15 and abs(t) < 1e-18 * max(abs(s), 1e-300):
            break
    return s

def Phi0(b):
    """closed form of Phi(a->0,beta) = (e^beta - c_0)/c_0 = beta e^beta/(2 I_1) - 1"""
    return b * np.exp(b) / (2.0 * Ib(1, b)) - 1.0

def cf_u(b, L):
    """L-th truncation of the Gauss J-fraction for u = 2 I_2/I_1 (a_k = 2(k+2)/b)"""
    tail = 2.0 * (L + 2) / b
    for k in range(L - 1, -1, -1):
        tail = 2.0 * (k + 2) / b + 1.0 / tail
    return 2.0 / tail

def blk(t):
    print("\n" + "=" * 80 + "\n" + t + "\n" + "=" * 80)

BETA_P3 = 1.0 / 48.0            # the proven P3/SZZ point
LAM_CRUDE = np.e * 20.0         # lambda_a = e*Delta, Delta = 4*(2(d-1)-1) = 20 in d=4
LO, HI = 1e-7, 20.0

# =====================================================================  [A]
blk("[A] exact SU(2) coefficients, the completeness identity, and the forms of u")
print(f"  {'beta':>9} {'SUM(n+1)c_n':>14} {'e^beta':>14} {'rel.err':>9} {'min_n c_n':>11}")
for b in [BETA_P3, 0.5, 1.0, 2.5]:
    tot = sum((n + 1) * c(n, b) for n in range(0, 140))
    mn = min(c(n, b) for n in range(0, 60))
    print(f"  {b:9.5f} {tot:14.10f} {np.exp(b):14.10f} "
          f"{abs(tot-np.exp(b))/np.exp(b):9.2e} {mn:11.3e}")
print("  => c_n = I_n - I_{n+2} is exact, and c_n > 0 for all n (SU(2), beta>0): the")
print("     character expansion is a POSITIVE sum over surfaces (no Weingarten sign at")
print("     this level) -- the random-surface / polymer representation is available.")

print("\n  four exact forms of u(beta) = c_1/c_0 and their agreement:")
print(f"  {'beta':>9} {'c1/c0':>16} {'2I2/I1':>16} {'2(I0/I1-2/b)':>16} {'max dev':>10}")
for b in [BETA_P3, 0.5, 1.0, 2.5]:
    a1 = c(1, b) / c0(b)
    a2 = 2.0 * Ib(2, b) / Ib(1, b)
    a3 = 2.0 * (Ib(0, b) / Ib(1, b) - 2.0 / b)
    print(f"  {b:9.5f} {a1:16.13f} {a2:16.13f} {a3:16.13f} "
          f"{max(abs(a1-a2),abs(a1-a3),abs(a2-a3)):10.2e}")

print("\n  Gauss J-fraction (semi-exact; truncations alternate => two-sided rigorous bounds):")
print("     u(beta) = 2 / ( 4/beta + 1/( 6/beta + 1/( 8/beta + ... ) ) )")
print(f"  {'beta':>9} {'L=0':>11} {'L=1':>11} {'L=2':>11} {'L=3':>11} {'L=7':>13} {'exact':>13}")
for b in [BETA_P3, 0.5, 1.0, 2.5]:
    print(f"  {b:9.5f} " + " ".join(f"{cf_u(b,L):11.7f}" for L in (0, 1, 2, 3))
          + f" {cf_u(b,7):13.9f} {u_p3(b):13.9f}")

print("\n  rigorous bracket of beta_tube := u^{-1}(1) from the alternating truncations:")
print(f"  {'L':>3} {'root of u_L(beta)=1':>22}")
roots = []
for L in range(0, 9):
    r = brentq(lambda b: cf_u(b, L) - 1.0, 1e-6, 6.0)
    roots.append(r)
    print(f"  {L:3d} {r:22.10f}")
print(f"  even/odd subsequences converge monotonically to beta_tube; bracket after L=8:")
print(f"     beta_tube in ({min(roots[-2:0:-1]):.9f}, {max(roots[-2:0:-1]):.9f})   [rigorous, no scipy trust]")

# =====================================================================  [B]
blk("[B] exact dictionary between the two corpus normalisations:  u_P3 = 2 u_T5")
print(f"  {'beta':>9} {'u_P3=c1/c0':>14} {'u_T5=I2/I1':>14} {'ratio':>16} {'6 u_T5':>10}")
for b in [BETA_P3, 0.5, 1.0, 2.5]:
    print(f"  {b:9.5f} {u_p3(b):14.9f} {u_t5(b):14.9f} {u_p3(b)/u_t5(b):16.13f} "
          f"{6.0*u_t5(b):10.5f}")
print("  => T23-TEMPLATE's 6u = 3.04 (u=I2/I1=0.5072) and T1-POLYMER-RG's u(2.5)=1.01439")
print("     are THE SAME statement, differing by the exact factor 2.  No contradiction;")
print("     T23 sec.1.3's 1.09-1.40x 'discrepancy' vs T5 is exactly this factor plus rounding.")

# =====================================================================  [C]
blk("[C] monotonicity [T]: u = 2 (ln c_0)'  =>  u' = 2 Var(cos theta) >= 0")
print("  (Cauchy-Schwarz on the plaquette measure => u strictly increasing; hence")
print("   {beta: u < u_c} is an interval and 'the largest beta' is well defined.)")
print(f"  {'beta':>9} {'u':>14} {'u (finite diff)':>18} {'2 Var(cos)':>14} {'dev':>10}")
h = 1e-6
for b in [BETA_P3, 0.5, 1.0, 2.0, 2.5]:
    fd = (u_p3(b + h) - u_p3(b - h)) / (2 * h)
    v2 = 2.0 * var_cos(b)
    print(f"  {b:9.5f} {u_p3(b):14.9f} {fd:18.12f} {v2:14.9f} {abs(fd-v2):10.2e}")

# =====================================================================  [D]
blk("[D] the four requested numbers, both conventions")
print(f"  {'beta':>8} {'u_P3 = c1/c0':>16} {'u_T5 = I2/I1':>16} {'6 u_T5':>10} {'m1 a = 4ln(1/u)':>18}")
for lab, b in [("1/48", BETA_P3), ("0.5", 0.5), ("1", 1.0), ("2.5", 2.5)]:
    uu = u_p3(b)
    print(f"  {lab:>8} {uu:16.9f} {u_t5(b):16.9f} {6.0*u_t5(b):10.5f} "
          f"{4.0*np.log(1.0/uu):18.9f}")
print("  (m1 a < 0 at beta=2.5: the four-wall tube does not decay -- the dilute-tube")
print("   binding picture is INVERTED, not merely slowly convergent.)")

# =====================================================================  [D2]
blk("[D2] the KP sum is monotone in the grading a  =>  its sharpest member is a -> 0+")
print(f"  {'beta':>9} {'Phi(0)':>12} {'Phi(0.25)':>12} {'Phi(0.5)':>12} {'Phi(1)':>12} {'Phi(2)':>12}")
for b in [BETA_P3, 0.5, 1.0, 2.5]:
    print(f"  {b:9.5f} " + " ".join(f"{Phi(a,b):12.6e}" for a in (0.0, 0.25, 0.5, 1.0, 2.0)))
print("  all c_n>0 => Phi strictly increasing in a => the a-optimisation is trivial:")
print("     min_a Phi(a,beta) = Phi(0,beta) = (e^beta - c_0)/c_0 = beta e^beta/(2I_1) - 1 .")
print(f"  {'beta':>9} {'Phi(0) sum':>16} {'Phi(0) closed':>16} {'dev':>9} {'Phi(1)/Phi(0)':>14}")
for b in [BETA_P3, 0.5, 1.0, 2.5]:
    s = sum((n + 1) * (c(n, b) / c0(b)) for n in range(1, 80))
    print(f"  {b:9.5f} {s:16.10e} {Phi0(b):16.10e} {abs(s-Phi0(b)):9.2e} {Phi(1.0,b)/Phi0(b):14.6f}")
print("  => T1-POLYMER-RG's grading a=1 is CONSERVATIVE by the last column; the exact")
print("     (weakest) member of the class is the closed form Phi(0,beta).  Both give the")
print("     same verdict; the exact one moves beta_KP up, not down.")

# =====================================================================  [E]
blk("[E] ceiling theorem: the exact hierarchy of blocking points")
beta_tube = brentq(lambda b: u_p3(b) - 1.0, LO, HI)      # u = 1
beta_half = brentq(lambda b: u_p3(b) - 0.5, LO, HI)      # u = 1/2
beta_f0 = brentq(lambda b: Phi0(b) - 1.0, LO, HI)        # Phi(0) = 1  <=> f0 = 1/2
print(f"  beta_tube = u^-1(1)     = {beta_tube:.9f}   u=1: weight u^(4T) stops decaying")
print(f"  beta_KP   = u^-1(1/2)   = {beta_half:.9f}   ceiling of ANY KP criterion (lam>=1)")
print(f"  beta_f0   = Phi(0)^-1(1)= {beta_f0:.9f}   lambda_a=1, calibration-free (f0=1/2)")
print(f"  u(beta)=1/lambda_a      = {brentq(lambda b: u_p3(b)-1.0/LAM_CRUDE, LO, HI):.9f}"
      f"   crude u<1/lambda_a, lambda_a=e*Delta={LAM_CRUDE:.3f}")
print(f"\n  2.5 - beta_tube = {2.5-beta_tube:.6f}   ({(2.5-beta_tube)/2.5*100:.3f}% of 2.5)")
print(f"  2.5 / beta_KP   = {2.5/beta_half:.3f}      2.5 / beta_f0 = {2.5/beta_f0:.3f}")
print(f"  => sup over u_c in (0,1) of beta_max(u_c) = beta_tube = {beta_tube:.6f} < 2.5 :")
print(f"     NO threshold u_c<1 reaches beta=2.5.  The u-route's ceiling is constant-free.")

print("\n  lambda_a scan of the exact KP root  (root of lambda_a*Phi(0,beta)=1):")
print(f"  {'lambda_a':>10} {'beta_KP(lambda)':>18} {'2.5/beta':>10}")
for lam in [1.0, np.e, 20.0, LAM_CRUDE, 1e2, 1e3, 1e4]:
    r = brentq(lambda b: lam * Phi0(b) - 1.0, LO, HI)
    print(f"  {lam:10.3f} {r:18.9f} {2.5/r:10.2f}")
print(f"  even lambda_a = 1 (an animal constant below the structural floor) gives "
      f"{brentq(lambda b: 1.0*Phi0(b)-1.0, LO, HI):.4f}, still {2.5/brentq(lambda b: Phi0(b)-1.0, LO, HI):.1f}x short.")

print("\n  T1-POLYMER-RG's grading a=1, for continuity with the corpus log:")
for lam, lab in [(LAM_CRUDE, "e*Delta")]:
    r = brentq(lambda b: lam * Phi(1.0, b) - 1.0, LO, HI)
    print(f"    beta_max(KP, a=1, lambda_a={lab}) = {r:.9f}   (corpus log says beta* ~ 0.0067)")

lam_need = 1.0 / Phi0(2.5)
print(f"\n  the animal constant beta=2.5 WOULD need:  lambda_a < 1/Phi(0,2.5) = {lam_need:.6f}")
print(f"  structurally lambda_a >= 1 (at least one polymer meets each point)")
print(f"  => short by >= {1.0/lam_need:.2f}x at the very best grading; "
      f"{LAM_CRUDE*Phi0(2.5):.1f}x at lambda_a = e*Delta.")

# =====================================================================  [F]
blk("[F] the n-profile at beta=2.5: the failure is NOT 'u too big'")
print(f"  grading a=0 (the exact/weakest member).  terms (n+1)c_n/c_0 :")
print(f"  {'beta':>9} " + " ".join(f"{'n='+str(n):>9}" for n in range(1, 7))
      + f" {'sum':>9} {'share n>=2':>11}")
for b in [BETA_P3, 0.5, 1.0, 2.5]:
    ts = [(c(n, b) / c0(b)) * (n + 1) for n in range(1, 7)]
    tot = sum(ts)
    print(f"  {b:9.5f} " + " ".join(f"{t:9.4f}" for t in ts)
          + f" {tot:9.4f} {(tot-ts[0])/tot:11.4f}")
sh = lambda b: 1.0 - 2.0 * u_p3(b) / Phi0(b)
b50 = brentq(lambda b: sh(b) - 0.5, 1e-6, 5.0)
print(f"  share(n>=2) crosses 1/2 at beta = {b50:.6f}   (at a=1: "
      f"{brentq(lambda b: 1.0 - 2*u_p3(b)*np.e/Phi(1.0,b) - 0.5, 1e-6, 5.0):.6f})")
print("  at beta=2.5, 60% of the exact Ursell sum sits in n>=2 (90% at grading a=1).")
print("  => even u -> 0 by fiat would NOT restore convergence at beta=2.5: the plaquette")
print("     measure is dense in ALL representations (f0 = 0.165).  'u too big' is one of")
print("     two independent failures, and not the binding one.")

# =====================================================================  [G]
blk("[G] cross-checks against the corpus logs")
print(f"  {'beta':>9} {'u':>12} {'u/(1-u)':>14}  (BN-RESIDUAL N3 geometric majorant)")
for b in [BETA_P3, 0.5, 1.0, 2.5]:
    uu = u_p3(b)
    print(f"  {b:9.5f} {uu:12.8f} {uu/(1.0-uu):14.4f}")
print(f"  at beta=2.5: u/(1-u) = {u_p3(2.5)/(1-u_p3(2.5)):.4f}  (BN-RESIDUAL: -70.49) =>")
print("  the geometric majorant has NO VALUE: 'summability of epsilon' is not even statable.")
print(f"\n  T23-TEMPLATE  : 6 u_T5(2.5) = {6*u_t5(2.5):.4f}   (log: 3.04)")
print(f"  T1-POLYMER-RG : u(2.5) = {u_p3(2.5):.5f} (log: 1.01439); "
      f"f0(2.5) = {c0(2.5)/np.exp(2.5):.5f} (log: 0.165)")
print(f"  P3 four-wall  : m1 a(2.5) = {4*np.log(1/u_p3(2.5)):.5f}  (log: -0.05715)")

print("\n  free-field acid test  beta -> infinity:")
print(f"  {'beta':>9} {'u_P3':>12} {'Phi(0)':>14} {'rho_Ur (lam=eD)':>17}")
for b in [5.0, 10.0, 50.0]:
    print(f"  {b:9.5f} {u_p3(b):12.8f} {Phi0(b):14.6e} {LAM_CRUDE*Phi0(b):17.6e}")
print("  u -> 2, rho_Ur -> infinity: the diagnostic reports 'no expansion' on the free")
print("  field, where the true answer is f1^2 = 0.  The no-go manufactures no false atom.")

# =====================================================================  summary
blk("SUMMARY (labels: [T]/[C]/[NG]/[O]/[V])")
print(f"[A] c_n = I_n-I_{{n+2}} exact; SUM(n+1)c_n = e^beta; c_n>0                     [T]+[V]")
print(f"[A] u = c1/c0 = 2I2/I1 = 2(ln c0)' = Gauss J-fraction                          [T]+[V]")
print(f"[B] u_P3 = 2 u_T5 EXACTLY  => T23's 6u=3.04 == T1's u=1.01439                   [T]+[V]")
print(f"[C] u' = 2 Var(cos theta) >= 0  =>  {{u<u_c}} is an interval                     [T]+[V]")
print(f"[D] u(1/48)={u_p3(BETA_P3):.6f}  u(0.5)={u_p3(0.5):.6f}  u(1)={u_p3(1.0):.6f}  "
      f"u(2.5)={u_p3(2.5):.6f}  [V]")
print(f"[D2] Phi(0,beta) = beta e^beta/(2I_1) - 1  (a->0 is the sharpest grading)       [T]+[V]")
print(f"[E] beta_tube = {beta_tube:.6f}  <-- ABSOLUTE ceiling of the u-route (u=1)")
print(f"    beta_KP   = {beta_half:.6f}  <-- ceiling of ANY KP criterion (lam>=1, u=1/2)")
print(f"    beta_f0   = {beta_f0:.6f}  ;  beta_KP(lambda_a=e*Delta) = "
      f"{brentq(lambda b: LAM_CRUDE*Phi0(b)-1.0, LO, HI):.6f}")
print(f"    beta=2.5 > {beta_tube:.4f}  =>  u-completeness CANNOT be pushed to 2.5     [NG] class")
print(f"    lambda_a needed at 2.5: < {lam_need:.4f}  (structural floor 1)  => "
      f">= {1.0/lam_need:.2f}x short")
print(f"[F] share(n>=2) at 2.5: 60% (a=0), 90% (a=1); 50% crossing at beta={b50:.4f}       [V]")
print(f"[G] u/(1-u)(2.5) = {u_p3(2.5)/(1-u_p3(2.5)):.2f};  m1 a = "
      f"{4*np.log(1/u_p3(2.5)):.5f} (inverted)   [T]+[V]")

print("\nMD-TABLE u(beta):")
for b in [1.0/48.0, 0.5, 1.0, 2.5]:
    print(f"   beta={b:.10f}  u_P3={u_p3(b):.10f}  u_T5={u_t5(b):.10f}  "
          f"6u_T5={6*u_t5(b):.6f}  m1a={4*np.log(1/u_p3(b)):+.6f}")
print("MD-THRESHOLDS:")
print(f"   beta_tube={beta_tube:.6f}  beta_KP={beta_half:.6f}  beta_f0={beta_f0:.6f}  "
      f"beta_kp_lam_eD={brentq(lambda b: LAM_CRUDE*Phi0(b)-1.0, LO, HI):.6f}  "
      f"beta_kp_a1={brentq(lambda b: LAM_CRUDE*Phi(1.0,b)-1.0, LO, HI):.6f}  "
      f"lam_need={lam_need:.6f}  b50={b50:.6f}")

print("\nNo QED. Clay 0%.")
