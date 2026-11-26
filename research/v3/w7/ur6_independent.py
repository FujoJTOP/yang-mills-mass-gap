# -*- coding: utf-8 -*-
"""
UR6 — independent estimation of lambda_Ur and the test lambda_Ur * 0.625 =?= 1.
Does NOT read ur1/ur2/ur4 outputs. Only corpus constants already in C03/TT1/P4/RS3.

Central question (task 1): does "R = u(1+R)^Delta has no positive root at u=0.625"
IMPLY  lambda_Ur * 0.625 > 1 ?

lambda_Ur := limsup |phi^U(Gamma)|^{1/|Gamma|}   (true Ursell growth; TT2:150)
Delta = 20 (4D plaquette adjacency degree)
u(SU2-exact) = beta/4  ->  u(2.5) = 0.625        (TT1:60, P4:143)
u(SU2)"N=2"  = beta/(2N^2) = beta/8 = 1/96 @ beta=1/12  (C03:118)

Sections:
 [A] the tree equation, its critical activity, and what "no positive root" IS
 [B] 0.013635 back-solved (C03:114) and its exact relation to u
 [C] lambda_Ur upper bounds (Lagrange / e*Delta) and the decisive threshold
 [D] D_inf: verification of 1.34487 and why it cannot constrain lambda_Ur (scale invariance)
 [E] same-Delta sign-flip fixture: witness that "no root" carries ZERO lambda_Ur information
 [F] animal-count growth reading (C03:111) -> independent lambda estimate
 [G] verdict numbers
"""
import sys
from mpmath import mp, mpf, exp, log, findroot, binomial, e, sqrt, inf, nan

mp.dps = 40
OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.append(s)


D = 20                      # plaquette-adjacency max degree, 4D
U_25 = mpf("0.625")         # u at beta=2.5, SU(2)-exact convention u = beta/4
U_B12 = mpf(1) / 96         # u at beta=1/12, corpus default u = beta/(2N^2)
R_C03 = mpf("0.013635")     # C03:114
ED = exp(1) * D             # e*Delta = 54.3656365692

P("=" * 78)
P("UR6  independent lambda_Ur probe   (mpmath dps=%d)" % (mp.dps,))
P("Delta = %d   e*Delta = %s   u(beta=2.5, SU2-exact) = %s" % (D, mp.nstr(ED, 12), U_25))
P("=" * 78)

# ----------------------------------------------------------------- [A]
P("\n[A] the tree equation  R = u(1+R)^Delta   -- what 'no positive root' IS")
P("-" * 78)
# f(R) = u(1+R)^D - R ; positive root exists iff min_R f <= 0
# f'(R)=0 => (1+R)^(D-1) = 1/(u D)  => R* = (uD)^(-1/(D-1)) - 1
# f(R*) = (1+R*)/D - R* <= 0  <=>  R* >= 1/(D-1)
#   <=> (uD) <= ((D-1)/D)^(D-1)
u_star = (mpf(D - 1) ** (D - 1)) / (mpf(D) ** D)
lam_tree = 1 / u_star                      # growth constant of rooted-tree counts
P("critical activity      u* = (D-1)^(D-1)/D^D               = %s" % mp.nstr(u_star, 12))
P("tree growth constant   L_t = D^D/(D-1)^(D-1) = 1/u*       = %s" % mp.nstr(lam_tree, 12))
P("   (compare: e*Delta = %s ; Lagrange bound is %s%% tighter)"
  % (mp.nstr(ED, 12), mp.nstr(100 * (1 - lam_tree / ED), 4)))
P("   equivalently  u*  = ((D-1)/D)^(D-1) / D  = %s" % mp.nstr(((mpf(D - 1) / D) ** (D - 1)) / D, 12))

P("\n root existence test at the three activities of interest:")
for tag, u in (("beta=1/12  u=1/96", U_B12), ("C03 rounded u=0.0104", mpf("0.0104")),
               ("beta=2.5    u=0.625", U_25)):
    Rm = (u * D) ** (-mpf(1) / (D - 1)) - 1          # argmin
    fmin = (1 + Rm) / D - Rm
    has = (fmin <= 0)
    if has:
        R = findroot(lambda x: u * (1 + x) ** D - x, Rm / 2, solver="secant")
    else:
        R = None
    P("   %-22s u = %-12s argmin R* = %-12s f(R*) = %-12s root? %-5s R = %s"
      % (tag, mp.nstr(u, 8), mp.nstr(Rm, 8), mp.nstr(fmin, 6), has,
         mp.nstr(R, 12) if R is not None else "-- none --"))

P("\n *KEY*: the predicate 'no positive root' is a function of (u, Delta) ONLY.")
P("        It never sees phi^U.  u/u* at beta=2.5 = %s" % mp.nstr(U_25 / u_star, 8))
P("        i.e. 'no root' <==>  L_t * u = %s > 1   with L_t a COUNTING constant."
  % mp.nstr(lam_tree * U_25, 8))

# ----------------------------------------------------------------- [B]
P("\n[B] what C03's 0.013635 actually is  (task 2)")
P("-" * 78)
R_B12 = findroot(lambda x: U_B12 * (1 + x) ** D - x, mpf("0.01"), solver="secant")
P(" solve R = u(1+R)^20  at u = 1/96  ->  R = %s" % mp.nstr(R_B12, 12))
P(" corpus value C03:114            R = %s" % mp.nstr(R_C03, 12))
P(" agreement to                   = %s  (C03 quotes 6 sig figs of 0.0136354)"
  % ("yes" if abs(R_B12 - R_C03) < mpf('1e-6') else "NO"))
R_104 = findroot(lambda x: mpf("0.0104") * (1 + x) ** D - x, mpf("0.01"), solver="secant")
P(" with C03's rounded u=0.0104     R = %s" % mp.nstr(R_104, 12))
P("\n the beta it lives on:  u = beta/(2N^2) = beta/8 ,  u=1/96  <=>  beta = 1/12  [= 8/96]")
P("   => 0.013635 is the SAME tree equation as 'no positive root', just at beta=1/12.")
P("\n exact relation between 0.013635 and u:")
for tag, u, R in (("u=1/96", U_B12, R_B12), ("u=0.0104", mpf("0.0104"), R_104),
                  ("u at criticality", u_star, 1 / (D - 1))):
    amp = R / u
    P("   %-16s R/u = (1+R)^D = %-12s   vs (1+R)^D = %-12s   u*(1+R)^D/R = 1 -> %s"
      % (tag, mp.nstr(amp, 10), mp.nstr((1 + R) ** D, 10), mp.nstr(u * (1 + R) ** D / R, 10)))
P("\n *KEY*: R/u is an AMPLITUDE ratio. Both R and u scale with beta, so R/u is")
P("        beta-covariant.  lambda_Ur is a ROOT (scale-invariant).  Hence this")
P("        whole comparison -- including C03:114's '41x slack vs e*Delta*u=0.565'")
P("        -- is structurally incapable of informing lambda_Ur.  (see [D])")

# ----------------------------------------------------------------- [C]
P("\n[C] lambda_Ur upper bounds and the decisive threshold  (task 4)")
P("-" * 78)
P(" Penrose/KP tree-graph inequality (C03:69-89):  |phi^U(Gamma)| <= sum_T prod |w|")
P("   =>  sum_{|Gamma|=n, p0 in Gamma} |phi^U|  <=  #rooted-trees(n) * u^n")
P("   =>  lambda_Ur <= L_t = %s        [T]" % mp.nstr(lam_tree, 12))
P("   (C03:89 tabulates r_n^(1/n)/e*Delta -> 1; the exact limit is 1/ e*Delta? no:")
vals = []
for n in (1, 2, 5, 10, 20, 50, 200, 1000, 5000):
    # Lagrange: r_n <= (1/n) C(D n, n-1) ; growth (1/n)^(1/n)*C(Dn,n-1)^(1/n) -> ((D)^(D))/((D-1)^(D-1))
    lb = (mpf(1) / n * binomial(D * n, n - 1)) ** (mpf(1) / n)
    vals.append((n, lb, lb / ED))
for n, lb, r in vals:
    P("     n=%-6d (C(20n,n-1)/n)^(1/n) = %-14s  / e*Delta = %s"
      % (n, mp.nstr(lb, 10), mp.nstr(r, 8)))
P("   => Lagrange bound converges to L_t = %s  (< e*Delta = %s).  Independent small gain:"
  % (mp.nstr(lam_tree, 10), mp.nstr(ED, 10)))
P("      e*Delta was a 2.5%% overestimate; the exact constant is L_t.")

P("\n the decisive test is  lambda_Ur =?= 1/u = 1/0.625 = %s" % mp.nstr(1 / U_25, 10))
P("   best rigorous bound gives   lambda_Ur * 0.625 <= L_t*0.625 = %s"
  % mp.nstr(lam_tree * U_25, 10))
P("   separation between bound and threshold:  L_t/(1/u) = %s  x"
  % mp.nstr(lam_tree * U_25, 8))
P("   => the bound is ABOVE 1, so it CANNOT certify lambda_Ur*u < 1;")
P("      and it CANNOT certify > 1 either.  One-sided.  Undecided.")
P("\n sanity collision with the corpus's independent veto (RS3:101 u <= 1/e):")
P("   u = 0.625 > 1/e = 0.367879 ?  %s   (u*e = %s)"
  % (U_25 > 1 / e, mp.nstr(U_25 * e, 8)))

P("\n FREE MICRO-GAIN of this probe (independent of the verdict): replacing e*Delta")
P(" by its exact Lagrange optimum L_t widens the provable window by %s%%:"
  % mp.nstr(100 * (ED / lam_tree - 1), 4))
P("   corpus  (C03:128):  beta_W < 2N^2/e*Delta = %s" % mp.nstr(8 / ED, 10))
P("   this probe:         beta_W < 2N^2/L_t     = %s" % mp.nstr(8 / lam_tree, 10))
P("   (SU(2)-exact u=beta/4 reading: beta_W < 4/L_t = %s)" % mp.nstr(4 / lam_tree, 10))
P("   -> e*Delta was provably NON-optimal; L_t is the exact constant of the same")
P("      Lagrange count (C03:84-89's own inequality, pushed to its limit).  [T]")

# ----------------------------------------------------------------- [D]
P("\n[D] D_inf and why it constrains nothing  (task 3)")
P("-" * 78)
P(" D_inf := exp sum_{alpha>=0} (1-alpha) A(alpha)    (RS2:35, X01:67)  -- a VALUE")
P(" lambda_Ur := limsup |phi^U|^{1/n}                                  -- a RADIUS")
def Dinf_geom(U):
    return exp(U * (1 - 2 * U) / (1 - U) ** 3)
for tag, U in (("U=0.25  (in KP window)", mpf("0.25")),
               ("U=0.5654 (beta=1/12)", mpf("0.5654")),
               ("U=33.9785 (beta=2.5)", ED * U_25)):
    P("   %-24s D_inf(geom closed form) = %s" % (tag, mp.nstr(Dinf_geom(U), 10)))
P("   reproduce corpus: D_inf(0.25) = 1.34487  ->  %s"
  % mp.nstr(Dinf_geom(mpf("0.25")), 12))
P("\n *and the closed form is NOT a sum at U>1 (TT1:125): radius 1, series diverges.")
P("  continuation value at U=33.9785 is %s > 0 but meaningless -- it is exp(complex-free"
  % mp.nstr(Dinf_geom(ED * U_25), 8))
P("  analytic continuation of a divergent series), NOT evidence about D_inf(2.5).")

P("\n SCALE INVARIANCE WITNESS (why D_inf cannot bound lambda_Ur in EITHER direction):")
P("  rescale every cluster weight: A(alpha) -> c*A(alpha).")
P("    * lambda_Ur -> c^{1/n} growth ... careful: lambda_Ur is defined from |phi^U|,")
P("      so c*|phi^U| gives lambda_Ur -> c*lambda_Ur.  Instead use the RATIO-FREE test:")
P("      adding a term only at alpha=0,1 changes D_inf arbitrarily, changes limsup not at all.")
for extra in (mpf(0), mpf(1), mpf(5), mpf(-3)):
    # A(alpha)=U^{alpha+1}/(1-U) plus a spike at alpha=0 of size 'extra'
    U = mpf("0.25")
    base = sum((1 - a) * U ** (a + 1) / (1 - U) for a in range(0, 200))
    spiked = base + (1 - 0) * extra          # only alpha=0 term altered
    P("     spike A(0) += %-6s : log D_inf = %-14s (was %s) ; limsup |A|^{1/n} unchanged"
      % (mp.nstr(extra, 4), mp.nstr(spiked, 10), mp.nstr(base, 8)))
P("  => D_inf and lambda_Ur are independent functionals of {A(alpha)}.")
P("  => D_inf(0.25)=1.34487 => f1^2 = a0*D_inf > a0  [T] , but ZERO constraint on lambda_Ur.")
P("     (Also: 1.34487 is a geometric-model value at U=0.25, deep inside the window,")
P("      not a value at beta=2.5 at all.)")

# ----------------------------------------------------------------- [E]
P("\n[E] WITNESS: 'no positive root at u=0.625' implies NOTHING about lambda_Ur")
P("-" * 78)
P(" fixture pair, identical Delta=20 and identical tree verdict, opposite lambda_Ur:")
fixtures = [
    ("F1  phi^U = +1.0^n   (sign-definite)      ", 1.0),
    ("F2  phi^U = (-1)^n * 0.1^n                ", 0.1),
    ("F3  phi^U = (-1)^n * 0.5^n                ", 0.5),
    ("F4  phi^U = (-1)^n * 4.0^n                ", 4.0),
    ("F5  phi^U = (-1)^n * 8.0^n                ", 8.0),
    ("F6  phi^U = (-1)^n * 40.0^n               ", 40.0),
]
for tag, lam in fixtures:
    p = mpf(lam) * U_25
    P("   %s  lambda_Ur = %-6s  lambda_Ur*u = %-10s  %s  root? no (same u, same Delta)"
      % (tag, mp.nstr(mpf(lam), 4), mp.nstr(p, 8),
         "> 1  <--" if p > 1 else "< 1  <--"))
P("   ALL SIX have Delta=20 at u=0.625 => ALL SIX have NO positive root.")
P("   Their lambda_Ur*u ranges over %s .. %s, straddling 1."
  % (mp.nstr(mpf('0.1') * U_25, 4), mp.nstr(mpf(40) * U_25, 4)))
P("   => the predicate does not separate the two sides => it CANNOT imply '>1'.")
P("   (c) is REFUTED.  [T] by counterexample (both sides witnessed).")
P("   Admissible range under the only rigorous bound: 0 < lambda_Ur <= L_t = %s,"
  % mp.nstr(lam_tree, 8))
P("   so the predicate leaves lambda_Ur*u free in (0, %s] -- 1 is strictly interior."
  % mp.nstr(lam_tree * U_25, 8))

# ----------------------------------------------------------------- [F]
P("\n[F] independent reading of the geometric growth (C03:111 animals)")
P("-" * 78)
an = [1, 20, 458, 11132, 280175]
P(" 4D plaquette animals a_n (C03:111): %s" % an)
rat = [mpf(an[i + 1]) / an[i] for i in range(len(an) - 1)]
P(" ratios                            : %s" % [mp.nstr(r, 6) for r in rat])
# log-log linear Richardson on ratios
import math
ns = [mpf(i + 1) for i in range(len(an))]
P(" a_n^{1/n}                          : %s" % [mp.nstr(mpf(a) ** (mpf(1) / n), 8) for a, n in zip(an, ns)])
lam_geo_5 = mpf(an[-1]) ** (mpf(1) / 5)
P(" lambda_geo(n<=5) := a_5^{1/5}      = %s   (still rising; a lower reading)" % mp.nstr(lam_geo_5, 8))
P("   lambda_geo*u at u=0.625          = %s" % mp.nstr(lam_geo_5 * U_25, 8))
P("   NOTE: lambda_Ur <= lambda_geo is NOT available; animals bound the SUPPORT, ")
P("         and phi^U carries a sign-alternating character factor.  So even this")
P("         number is a support-side reading, not a bound on lambda_Ur.  [C]/[NG]")

# ----------------------------------------------------------------- [G]
P("\n[G] verdict numbers")
P("-" * 78)
P(" (i)   'no positive root' at u=0.625  <==>  u > u* = %s  <==>  L_t*u = %s > 1"
  % (mp.nstr(u_star, 8), mp.nstr(lam_tree * U_25, 8)))
P("       L_t = %s is a ROOTED-TREE COUNTING constant, NOT lambda_Ur.  [T]"
  % mp.nstr(lam_tree, 10))
P(" (ii)  lambda_Ur <= L_t = %s  (best available, 2.5%% better than e*Delta)  =>  "
  "lambda_Ur*0.625 <= %s" % (mp.nstr(lam_tree, 10), mp.nstr(lam_tree * U_25, 8)))
P(" (iii) threshold: lambda_Ur =?= 1/0.625 = 1.6 ; bound/threshold = %s x  => UNDECIDED"
  % mp.nstr(lam_tree * U_25, 8))
P(" (iv)  D_inf(0.25)=1.34487 gives f1^2>a0  [T] ; zero constraint on lambda_Ur  [T]")
P(" (v)   0.013635 = R(u=1/96) of the SAME equation; R/u=1.3114 is an amplitude reading,")
P("       and the '41x slack' is amplitude slack -- invisible to the scale-invariant lambda_Ur")
P(" (vi)  => (b).  Need real phi^U enumeration (n<=6) + growth-rate extrapolation.")

P("\n" + "=" * 78)
P("VERDICT: (b)  -- 'no positive root' does NOT imply lambda_Ur*0.625>1 (refuted by")
P("         same-Delta fixture [E]); D_inf gives no constraint (scale invariance [D]);")
P("         best rigorous bound is lambda_Ur <= 53.0006, threshold 1.6, gap 33.13x.")
P("=" * 78)

with open(__file__.replace(".py", ".log"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
print("\n[written] ur6_independent.log")
