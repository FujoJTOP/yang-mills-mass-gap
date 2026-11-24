"""
route_new_holo.py -- ROUTE-NEW-HOLO witness.

Hard-wall holographic 0++ spectrum: m_n * z0 = j_{2,n} (zeros of J_2).
Demonstrates (a) the exact scale relation, (b) bounding the gap <=> bounding z0
from above, (c) NEC is SATURATED in pure AdS for every z0 (no scale constraint),
(d) why holography gives the scale-free RATIO but never the absolute scale.

numpy/scipy only, <1 s. No claim of a gap.
"""
import numpy as np
from scipy.special import jn_zeros

print("=" * 72)
print("Hard-wall 0++ spectrum:   m_n * z0 = j_{2,n}   (zeros of Bessel J_2)")
print("=" * 72)
j = jn_zeros(2, 6)
for n in range(5):
    print(f"  n={n+1}: j_2,{n+1} = {j[n]:9.5f}   =>  m_{n+1} = {j[n]:9.5f} / z0")
print(f"  scale-free ratio m(0++*)/m(0++) = j_2,2/j_2,1 = {j[1]/j[0]:.4f}  (lattice ~1.6-1.64)")

print()
print("=" * 72)
print("Bounding the gap  <=>  bounding z0 from above")
print("=" * 72)
delta = 1.0
print(f"  m_0++ = 5.1356/z0 >= delta  <=>  z0 <= {j[0]/delta:.4f}  (IR wall bounded)")
print("  a mass gap forces z0 finite; z0 -> infinity (pure AdS) = no gap = continuum.")

print()
print("=" * 72)
print("NEC in pure AdS is SATURATED for every z0 (constrains SHAPE, not SCALE)")
print("=" * 72)
# T_MN = -Lambda g_MN ; for a NULL vector k, g_MN k^M k^N = 0
# => T_MN k^M k^N = -Lambda * 0 = 0  for ANY Lambda  (i.e. ANY z0 / warp scale).
Lambda = -6.0          # AdS_5: Lambda = -(d-1)(d-2)/(2R^2), d=5 -> -6/R^2
null_norm = 0.0        # by definition of a null vector
print(f"  T_MN k^M k^N = -Lambda * (k^M k_M) = ({-Lambda}) * ({null_norm}) = {-Lambda*null_norm:.1f}")
print("  => NEC holds with EQUALITY for every z0; it cannot bound the confinement scale.")

print()
print("=" * 72)
print("Why holography gives the RATIO but never the absolute scale")
print("=" * 72)
# string tension from the geometry: sigma ~ c / z0^2  =>  m/sqrt(sigma) = j_2,1 / c  (pure number)
c = 1.0
print(f"  sigma ~ c / z0^2   =>   m_0++ / sqrt(sigma) = j_2,1 / c = {j[0]/c:.4f}   (scale-free)")
print("  => holography matched m/sqrt(sigma) ~ 3.9 only by fixing the geometry ratio c;")
print("     the absolute scale z0 (= confinement scale = H1) is never predicted.")

print()
print("=" * 72)
print("VERDICT: holographic bound on m_0++ or chi reduces to ASSUMING the IR wall")
print("scale z0 = the confinement scale = H1 in the bulk-geometry coordinate.")
print("=" * 72)
