"""
roof_dualsuperconductor.py -- witness for research/unconditional/ROOF-DUALSUPERCONDUCTOR.md
ONE small job, numpy + math only, no heavy MC / no scipy.
Run:  py -3.12 -I research/unconditional/w7/roof_dualsuperconductor.py

Question: does the DUAL-Ginzburg-Landau ('t Hooft-Mandelstam) 0++ glueball mass
escape the u=2 IR renormalon that blocks the gluon-condensate floor?

Pieces:
  [1] dual GL mass spectrum: m_V = g v0, m_H = sqrt(2 lam) v0; kappa_GL = m_H/m_V.
      Algebraic, definite, finite (no series).
  [2] gap equation: dV/d|phi|=0 => |phi|=v0 exactly; m_H^2 = V''(v0) = 2 lam v0^2 > 0.
      A finite algebraic root -> the 0++ mass is a TREE/mean-field gap, no factorial growth.
  [3] contrast: the D=4 gluon condensate <tr F^2> has a u=2 IR renormalon:
      perturbative coefficients a_n ~ (2)^n Gamma(n+1+2b) -> a_n/a_{n-1} grows ~ 2(n+...) linearly
      => factorial divergence => Borel pole at u=2 => ambiguity ~ Lambda^4 => no scheme-indep floor.
  [4] m_H/sqrt(sigma) vs kappa_GL (London limit sigma ~ pi v0^2 ln kappa), BPS anchor at kappa=1/sqrt2.
      Shows the ratio is O(1)-O(4); lattice m_0++/sqrt(sigma) ~ 3.9 needs large kappa.

Labels [T]/[C]/[O]/[NG]. Nothing here is a QCD proof.
"""
import numpy as np
from math import gamma, log, sqrt, pi

np.set_printoptions(precision=4, suppress=True)

print("=" * 78)
print("[1] dual GL mass spectrum (g=1, v0=1): m_V = g v0, m_H = sqrt(2 lam) v0")
print("    kappa_GL := m_H/m_V = sqrt(2 lam)/g")
print("=" * 78)
v0 = 1.0
g = 1.0
print(f"{'kappa_GL':>9} {'lam':>8} {'m_V':>8} {'m_H':>8} {'sigma_London':>13} {'m_H/sqrt(sig)':>14}")
for kappa in [0.7071, 1.0, 2.0, 5.0, 10.0, 20.0]:
    lam = 0.5 * (kappa * g) ** 2
    mH = sqrt(2 * lam) * v0
    mV = g * v0
    if kappa > 1.0:
        sig = pi * v0 ** 2 * log(kappa)          # London-limit vortex tension (large kappa)
        ratio = mH / sqrt(sig)
    else:
        sig = float("nan"); ratio = float("nan")  # type-I: London formula n/a
    print(f"{kappa:9.4f} {lam:8.4f} {mV:8.4f} {mH:8.4f} {sig:13.4f} {ratio:14.4f}")
print("  BPS anchor (kappa=1/sqrt2): sigma = pi v0^2 exactly; m_H/sqrt(sig) = (1/sqrt2)/sqrt(pi) = %.4f"
      % (1/sqrt(2)/sqrt(pi)))
print("  NOTE: lattice m_0++/sqrt(sigma) ~ 3.9 is reached only at kappa ~ 11+ (London) ->")
print("        the plain dual-Higgs identification is LIGHTER than the true 0++ (composite/vortex-ring).")

print("\n" + "=" * 78)
print("[2] gap equation is algebraic: dV/d|phi|=0 at |phi|=v0; m_H^2 = V''(v0) = 2 lam v0^2")
print("=" * 78)
for lam in [0.5, 1.0, 2.0, 4.0]:
    mH2 = 2 * lam * v0 ** 2
    print(f"  lam={lam:4.1f}  m_H^2 = V''(v0) = 2*lam*v0^2 = {mH2:6.3f}  (>0, finite, no series)")

print("\n" + "=" * 78)
print("[3] contrast: D=4 gluon condensate u=2 IR renormalon -> factorial coeffs")
print("    a_n ~ (2)^n Gamma(n+1+2b);  a_n/a_(n-1) ~ 2(n+2b) grows linearly => pole at u=2")
print("=" * 78)
b = 0.5   # illustrative anomalous-dimension exponent (u=2 branch point (1-u/2)^(-1-2b))
def a_coeff(n):
    return (2.0 ** n) * gamma(n + 1 + 2 * b) / gamma(1 + 2 * b)
print(f"{'n':>3} {'a_n':>12} {'a_n/a_(n-1)':>14} {'2(n+2b)':>9}")
for n in [1, 2, 3, 5, 10, 20, 30]:
    an = a_coeff(n)
    if n == 1:
        print(f"{n:3d} {an:12.4g} {'(first)':>14} {'-':>9}")
    else:
        r = an / a_coeff(n - 1)
        print(f"{n:3d} {an:12.4g} {r:14.1f} {2*(n+2*b):9.1f}")
print("  => factorial growth: Borel integral has a singularity at u=2, ambiguity ~ Lambda^4.")

print("\n" + "=" * 78)
print("CONCLUSION")
print("=" * 78)
print("  dual GL 0++ mass: finite ALGEBRAIC gap (mean-field), NO factorial growth")
print("     => structurally renormalon-FREE (different operator: v is dim-2 monopole density,")
print("        natural Borel pole at u=1, and UV-finite on the dual lattice; not the u=2 of <tr F^2>).")
print("  BUT dual GL must ASSUME monopole condensation; deriving it from the QCD measure")
print("     is confinement itself = H1. So: a renormalon-free BASE, not a key.")
print("  Clay 0%; no QED.")
