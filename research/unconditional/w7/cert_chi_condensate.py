"""
cert_chi_condensate.py -- CERT-CHI-CONDENSATE witness  (numpy + scipy only)

Question: does the Savvidy-cell RANK-ONE form give a RENORMALON-FREE uniform
floor on S_bare = <tr F^2> (the bare coarse curvature)?

  G'' = S - D,   D = D_+ + D_Sav,   D_+ <= 0 (Brascamp-Lieb),
                 D_Sav = g^2 Var(phi^2) <= 2 g^2 / lambda4.

Part 1  D_Sav bounded, c-independent: Var(phi^2) -> 2/lambda4 as c=gB grows.   [T]
Part 2  Flooring G'' needs S_bare >= lambda + 2 g^2/lambda4  (from D <= const). [T]
Part 3  S_bare = condensate: u=2 IR renormalon, a_n ~ n!, ambiguity ~ Lambda^4. [T/V-num]
Part 4  rank-one D-bound is a one-sided constant, NOT a subtraction -> no cancel.[O/NG]
Verdict: no renormalon-free floor on S_bare; the constant is the u=2-ambiguous
         condensate = H1.  The rank-one form is a real renormalon-free bound on D,
         not a floor on S.  Clay 0%; no QED.
"""
import numpy as np
from math import gamma, exp, log

print("=" * 74)
print("Part 1  D_Sav = g^2 Var(phi^2): bounded, c-independent  (-> 2/lambda4)")
print("=" * 74)
lam4 = 1.0


def var_phi2(c, lam=lam4, n=800001):
    """Var(phi^2) under V = -(c/2) phi^2 + (lam/4) phi^4, measure ~ exp(-V).
    Wells sit at phi = +-sqrt(c/lam); scale the grid with them."""
    edge = 3.0 * np.sqrt(c / lam) + 6.0
    y = np.linspace(-edge, edge, n)
    V = -(c / 2.0) * y ** 2 + (lam / 4.0) * y ** 4
    V -= V.min()
    w = np.exp(-V)
    w /= np.trapezoid(w, y)
    y2 = y ** 2
    mean = np.trapezoid(w * y2, y)
    m2 = np.trapezoid(w * y2 ** 2, y)
    return m2 - mean ** 2, mean


for c in [10.0, 20.0, 40.0, 80.0, 160.0]:
    v, mp = var_phi2(c)
    print(f"  c=gB={c:6.1f}:  Var(phi^2)={v:.6f}   <phi^2>={mp:.4f}   (2/lambda4={2/lam4:.6f})")
print("  -> Var(phi^2) saturates at 2/lambda4, independent of the negative")
print("     curvature c.  D_Sav = g^2 Var <= 2 g^2/lambda4 : bounded, sign>=0. [T]")

print()
print("=" * 74)
print("Part 2  Flooring G'' requires flooring S_bare (>= lambda + 2 g^2/lambda4)")
print("=" * 74)
lam = 0.1
g = 0.5
Dmax = 2.0 * g ** 2 / lam4
print(f"  lam={lam},  g={g}:  D <= 2 g^2/lambda4 = {Dmax:.4f}")
print(f"  G'' = S - D >= S - {Dmax:.4f} ;  G'' >= lambda  <=>  S >= {lam + Dmax:.4f}")
print("  -> the floor lives entirely on S_bare; D's bound just adds a constant. [T]")

print()
print("=" * 74)
print("Part 3  S_bare = condensate:  u=2 IR renormalon, a_n ~ Gamma(n+1+2b)")
print("=" * 74)
b = 2.0          # u = 2  (dim-4 operator <tr F^2>)
c0 = 1.0         # Borel scale (toy)
prev = None
for n in [2, 3, 5, 10, 20, 30]:
    a_n = (c0 ** n) * gamma(n + 1 + 2 * b)
    ratio = "" if prev is None else f"   a_n/a_(n-1) ~ {a_n/prev:.1f}"
    print(f"  n={n:2d}:  a_n = {a_n:.3e}{ratio}")
    prev = a_n
# Borel transform B(u) ~ (1 - u/2)^{-1-2b}: simple pole at u=2
# lateral sum difference (ambiguity) = 2*pi*Res*e^{-2s},  s = 1/(c0*alpha) ~ a scale
print("  a_n ~ Gamma(n+1+2b) => Borel pole at u=2, ambiguity |S+ - S-| = 2*pi*Res*e^{-2s}.")
print("  Absolute ambiguity ~ Lambda^4 = the condensate's OWN magnitude. [T/V-num]")
print("  (This is exactly ROOF-CONDENSATE-FLOOR: no scheme-independent 'condensate' floor.)")

print()
print("=" * 74)
print("Part 4  the rank-one D-bound is a one-sided CONSTANT, not a subtraction")
print("=" * 74)
print("  G'' = S - D >= S - 2 g^2/lambda4   uses only  D <= const.")
print("  A one-sided inequality with a constant does NOT subtract off S's u=2 pole.")
print("  The EXACT difference G'' = S - D (physical 1/chi) WOULD cancel the u=2")
print("  renormalon, but the rank-one form supplies D <= const, not D exactly;")
print("  exact D = the conditioned fiber measure = H1. [O]")
print()
print("VERDICT: NO renormalon-free uniform floor on S_bare.")
print("  - the DEFECT D is bounded, renormalon-free, c-independent  [T]  (real gain)")
print("  - the CURVATURE floor S_bare is the u=2-ambiguous condensate = H1 [NG]")
print("  Clay 0%; no QED.")
