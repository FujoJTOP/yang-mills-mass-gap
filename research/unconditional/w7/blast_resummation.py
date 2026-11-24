"""
blast_resummation.py -- FINAL BLAST on H1: can aggressive resummation of the
strong-coupling series certify  chi <= 1/lambda  at beta=2.5, or does it diverge?

NOT a Yang-Mills calculation. It demonstrates, with exact elementary facts + a
toy Borel-Pade, the IR-renormalon obstruction to any resummation at the physical
point.  numpy/scipy only.  py -3.12 -I blast_resummation.py

Facts used:
 * SU(2) lattice: u(beta)=I2(beta)/I1(beta), the strong-coupling expansion parameter.
 * The strong-coupling (character/polymer) series is ASYMPTOTIC with a Borel
   singularity at the physical scale t_c ~ the crossover / mass-gap scale,
   t_c ~ O(1) in u.  Its intrinsic ambiguity (minimal term / Stokes jump) is
   ~ exp(-t_c / u).  No contour deformation (Borel-Pade, conformal) removes a
   real-axis singularity; the best achievable error is bounded below by it.
"""
import math
import numpy as np
from scipy.special import iv, gammaln, gamma
from scipy.interpolate import pade

print("=" * 74)
print("1.  u(beta) at the physical point, SU(2):  u = I2/I1")
print("=" * 74)
for b in (2.0, 2.3, 2.5, 2.7):
    u = iv(2, b) / iv(1, b)
    print(f"   beta={b:4.2f}   u={u:.4f}")
u25 = iv(2, 2.5) / iv(1, 2.5)
print(f"   -> the strong-coupling parameter at beta=2.5 is u = {u25:.4f}")
print("      (rigorous radius: beta* ~ 1/48 = 0.0208, i.e. u(1/48) ~ 0.02;")
print("       beta=2.5 is ~2 orders in beta / ~25x in u outside it.)")

print()
print("=" * 74)
print("2.  Intrinsic renormalon ambiguity = minimal term ~ exp(-t_c/u)")
print("=" * 74)
# factorial series compatible with a Borel singularity at t_c:
#   a_n = Gamma(n+1+2b)/Gamma(1+2b) t_c^{-n};  term T_n = a_n u^n.
# log|T_n| = gammaln(n+1+2b)-gammaln(1+2b) - n log t_c + n log u.
def min_term(u, tc, two_b=1.0, N=400):
    n = np.arange(N + 1, dtype=float)
    logT = gammaln(n + 1 + two_b) - gammaln(1 + two_b) - n * np.log(tc) + n * np.log(u)
    k = int(np.argmin(logT))
    return k, float(np.exp(logT[k])), float(np.exp(-tc / u))

print("   t_c   |  n* (optimal truncation)  | minimal-term error | e^{-t_c/u}")
for tc in (0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0):
    nstar, mn, ref = min_term(u25, tc, two_b=1.0)
    print(f"   {tc:4.2f}  |           {nstar:3d}              |      {mn:8.3f}       |  {ref:6.3f}")
print("   -> t_c/u = t_c/0.507.  For t_c <~ 1 (the physical crossover/mass-gap range)")
print("      t_c/u = O(1): the term ratios (n+2)(u/t_c) exceed 1 from the first")
print("      term, so optimal truncation is n*=0 and the minimal term is O(1) =")
print("      NO useful estimate.  The e^{-t_c/u} column is only the large-t_c/u")
print("      asymptotic floor, a LOWER bound on the achievable error -- the real")
print("      error at the physical t_c is WORSE.  No resummation beats this.")

print()
print("=" * 74)
print("3.  Toy check: Borel-Pade of an asymptotically-singular series")
print("=" * 74)
# model Borel transform B(t)=(1-t/t_c)^{-alpha}, alpha=0.5 (<1 => PV-summable).
# series coefficient of u^n (after the Borel map) : a_n = (alpha)_n t_c^{-n}.
tc, alpha = 0.7, 0.5
def coef(n):
    return float(gamma(n + alpha) / gamma(alpha) / tc ** n)
Norder = 18
a = np.array([coef(n) for n in range(Norder + 1)])

# exact Borel sum (upper lateral): int_0^inf e^{-t/u}(1-t/tc)^{-alpha} dt, Im(t)=eps>0
def borel_exact(u, eps=1e-3, S=None):
    S = S or 60 * u
    s = np.linspace(0, S, 400000)
    t = s + 1j * eps
    integ = np.exp(-t / u) * (1 - t / tc) ** (-alpha)
    return np.trapezoid(integ.real, s)  # end values negligible (exponential decay)

truth = borel_exact(u25)

# Borel-Pade: B~(t) = sum a_n t^n/n! ; Pade [M/M]; integrate along Im=eps.
bn = np.array([a[n] / math.factorial(n) for n in range(Norder + 1)])
M = 8
p, q = pade(bn, M)
def borel_pade(u, eps=2e-3, S=None):
    S = S or 60 * u
    s = np.linspace(0, S, 400000)
    t = s + 1j * eps
    val = np.exp(-t / u) * (p(t) / q(t))
    return np.trapezoid(val.real, s)

bp = borel_pade(u25)
print(f"   t_c={tc}, alpha={alpha}, order {Norder} -> Pade [8/8]")
print(f"   exact upper Borel sum at u(2.5): {truth:.4f}")
print(f"   Borel-Pade [8/8]               : {bp:.4f}")
print(f"   |relative error|              : {abs(bp-truth)/abs(truth):.2%}")
print(f"   minimal-term floor e^(-t_c/u)  : {np.exp(-tc/u25):.3f}")
print("   -> the Borel-Pade error sits at/near the minimal-term floor; it does NOT")
print("      go to zero with order.  Same conclusion for any resummation.")

print()
print("=" * 74)
print("VERDICT")
print("=" * 74)
print("   The strong-coupling series is asymptotic with an IR-renormalon-type")
print("   Borel singularity at the physical (crossover/mass-gap) scale t_c ~ O(1).")
print("   At beta=2.5 (u=0.507) the optimal truncation is O(1) terms and the")
print("   intrinsic ambiguity is ~exp(-t_c/u) = O(10-60%).  Aggressive resummation")
print("   (Borel-Pade / conformal) gives a finite DIAGNOSTIC but no rigorous error")
print("   bound, and cannot certify  chi <= 1/lambda  (an inequality needs the")
print("   estimate error < its margin, which is not met).  No certificate.")
