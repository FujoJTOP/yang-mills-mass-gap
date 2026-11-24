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
print("3.  Toy check: Borel-Pade, Borel-summable vs not")
print("=" * 74)
# model Borel transform B(t)=(1-t/t_c)^{-alpha}, alpha=0.5 (<1 => PV-summable).
# series coefficient of u^n: a_n = (alpha)_n t_c^{-n}; B~(t)=sum a_n t^n/n!.
alpha = 0.5
def build(tc, Norder=18, M=8):
    a = np.array([float(gamma(n + alpha) / gamma(alpha) / tc ** n) for n in range(Norder + 1)])
    bn = np.array([a[n] / math.factorial(n) for n in range(Norder + 1)])
    p, q = pade(bn, M)
    return p, q
def borel_exact(u, tc, eps, S=None):
    S = S or 60 * u
    s = np.linspace(0, S, 400000)
    t = s + 1j * eps
    return float(np.trapezoid((np.exp(-t / u) * (1 - t / tc) ** (-alpha)).real, s))
def borel_pade(u, tc, p, q, eps, S=None):
    S = S or 60 * u
    s = np.linspace(0, S, 400000)
    t = s + 1j * eps
    return float(np.trapezoid((np.exp(-t / u) * (p(t) / q(t))).real, s))

for label, tc in (("Borel-summable   t_c=0.7 > u", 0.7),
                  ("NOT summable     t_c=0.4 < u", 0.4)):
    p, q = build(tc)
    up = borel_exact(u25, tc, +2e-3)
    lo = borel_exact(u25, tc, -2e-3)
    bp = borel_pade(u25, tc, p, q, +2e-3)
    print(f"   {label}:")
    print(f"      upper lateral sum = {up:.4f}   lower = {lo:.4f}   "
          f"lateral difference = {abs(up-lo):.4f}")
    print(f"      Borel-Pade [8/8] (upper) = {bp:.4f}   |BP-upper|/|upper| = "
          f"{abs(bp-up)/abs(up):.2%}")
    print(f"      e^(-t_c/u) = {np.exp(-tc/u25):.3f}")
print("   -> KEY subtlety: the REAL parts of the lateral sums coincide (Borel-Pade")
print("      reproduces the principal value PV to <0.1% even for t_c<u).  So Borel-")
print("      Pade is numerically 'accurate' -- but PV is NOT the physical value:")
print("      the IR renormalon makes the physical answer  PV + sigma*exp(-t_c/u)")
print("      with sigma a FREE transseries coefficient (ROOF-RESURGENCE: no saddle).")
print("      So the physical chi is PV + sigma*e^(-t_c/u), sigma undetermined -> an")
print("      intrinsic O(e^(-t_c/u)) = O(20-45%) ambiguity that NO prescription")
print("      (Borel-Pade, conformal, median) removes.  That is the obstruction.")

print()
print("=" * 74)
print("VERDICT")
print("=" * 74)
print("   The strong-coupling series is asymptotic with an IR-renormalon Borel")
print("   singularity at the physical (crossover/mass-gap) scale t_c ~ O(1).")
print("   At beta=2.5 (u=0.507):")
print("     * plain truncation is useless (optimal n*=0, O(1) error);")
print("     * Borel-Pade is accurate ONLY IF the series is Borel-summable")
print("       (singularity beyond u(2.5), integrable) -- but that continuity is")
print("       UNKNOWN, and it carries a genuine transseries ambiguity sigma*e^{-2s}")
print("       (ROOF-RESURGENCE: no saddle), so even a converged value is one of a")
print("       family, and there is NO rigorous error bound (asymptotic series).")
print("   Hence: aggressive resummation gives a numerical DIAGNOSTIC [V-num]/[C],")
print("   NOT a certificate.  It cannot rigorously establish chi <= 1/lambda at")
print("   beta=2.5, because the estimate's error/ambiguity is not controlled below")
print("   the required margin.  No certificate; Clay 0%.")
