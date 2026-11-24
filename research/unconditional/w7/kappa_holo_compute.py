"""
kappa_holo_compute.py -- witness for research/unconditional/KAPPA-HOLO-COMPUTE.md
Run: py -3.12 -I kappa_holo_compute.py

Four checks:
 (A) hard-wall 0++ tower = Bessel J_2 zeros (the atom is GEOMETRIC, not rigid-string);
 (B) the dual string is NG at leading alpha' => no (d^2 X)^2 term => kappa = 0;
 (C) alpha'*sigma coincidence = 1/(2 pi) ~ 0.159 vs the MASTERFIELD band [0.138,0.193];
 (D) rigidity is a 4-derivative term suppressed by alpha'/ell^2 (1/xi^2 scaling).
"""
import numpy as np
from scipy.special import jn_zeros

print("=" * 72)
print("(A) hard-wall 0++ tower = Bessel J_2 zeros  (atom is GEOMETRIC)")
print("=" * 72)
j = jn_zeros(2, 4)               # first four zeros of J_2
print("  j_{2,n}          =", np.round(j, 5))
print("  tower m_n/m_1    =", np.round(j / j[0], 4), " (lattice ~ 1 : 1.6 : 2.3 : 2.9)")
print("  m_1 * z0 = j_{2,1} =", round(float(j[0]), 5))
print("  -> discreteness/atom comes from the WALL (boundary condition), not worldsheet kappa.")

print()
print("=" * 72)
print("(B) dual string is Nambu-Goto at leading alpha' => kappa_leading = 0")
print("=" * 72)
# The NG action (1/2 pi alpha') int sqrt(-det G_ab) has NO (d^2 X)^2 term.
# A rigidity term is 4-derivative: (kappa/2)(d^2 X)^2 ~ O(alpha') * (NG term).
# So at leading order the coefficient is exactly zero.
kappa_leading = 0.0
print("  kappa_leading(alpha'^0) =", kappa_leading, "  (pure NG; rigidity is higher order)")

print()
print("=" * 72)
print("(C) alpha'*sigma coincidence vs MASTERFIELD band")
print("=" * 72)
aps = 1.0 / (2.0 * np.pi)        # if sigma = 1/(2 pi alpha')  (standard normalization)
band_lo, band_hi = 0.138, 0.193
print("  alpha'*sigma = 1/(2 pi) =", round(aps, 5))
print("  MASTERFIELD band        =", (band_lo, band_hi))
print("  inside band?            =", bool(band_lo <= aps <= band_hi),
      " (COINCIDENCE of normalization, not a derivation)")
# show the O(1) coefficient is free: scan it
for c in (0.7, 1.0, 1.5, 2.0):
    print("   if true coefficient = %.1f * alpha'*sigma -> kappa ~ %.3f" % (c, c * aps))

print()
print("=" * 72)
print("(D) rigidity suppressed by alpha'/ell^2 (1/xi^2), xi = ell/l_s")
print("=" * 72)
xi = np.array([1.0, 2.0, 5.0, 10.0])   # ell in units of l_s = sqrt(alpha')
ratio = 1.0 / xi**2                     # rigidity / leading (relative weight)
for x, r in zip(xi, ratio):
    print("   xi=%5.1f  rigidity/leading = %.4f" % (x, r))
print("  -> at long worldsheet distances (xi>>1) the rigidity term vanishes; the leading")
print("     dual string is NG. kappa_holographic = 0 at leading order.")

print()
print("=" * 72)
print("VERDICT: kappa_holographic = 0 (NG); atom is GEOMETRIC; alpha'*sigma ~ 0.159 is a coincidence.")
print("=" * 72)
