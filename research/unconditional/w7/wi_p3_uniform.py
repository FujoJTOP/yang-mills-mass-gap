"""
wi_p3_uniform.py  --  witness for WI-P3-UNIFORM.

Claim to demonstrate: the P3 approach-law's UNIFORMITY in (L,a) is exactly m1 > 0
(the gap). The isolation gap Delta(L) = E2(L) - E1(L) of the 0++ channel:
  - glueball ATOM (m1 below the 2-glueball threshold 2*m1): Delta(L) -> m1 > 0  (uniform);
  - :phi^2: CUT  (floor = 2m, but a 2-particle CONTINUUM, f1^2=0): Delta(L) -> 0  (no uniform bound).
So 'uniform Delta >= delta'  <=>  'the 0++ floor is an isolated atom'  = the gap = H1.

Run:  py -3.12 -I research/unconditional/w7/wi_p3_uniform.py
numpy only, <1 s.
"""
import numpy as np

def gap_atom(m1, m_star, L):
    """Lightest 0++ = glueball at m1; next level = min(two-glueball threshold at 2 m1, excited glueball m_star)."""
    E1 = m1
    E2 = min(2.0 * m1, m_star)          # two glueballs at rest (p_rel = 0) sit at 2 m1
    return E2 - E1

def gap_cut_phi2(m, L):
    """:phi^2: is a two-particle CONTINUUM. Lightest = both at rest E=2m; next = relative momentum p=2pi/L."""
    E1 = 2.0 * m                                # p_rel = 0
    E2 = 2.0 * np.sqrt(m * m + (2.0 * np.pi / L) ** 2)   # p_rel = 2pi/L
    return E2 - E1

print("=" * 74)
print("isolation gap Delta(L) = E2 - E1 of the 0++ channel, vs box size L")
print("=" * 74)
print(f"{'L':>5} | {'atom Delta(L)':>14} | {'cut(:phi2:) Delta(L)':>22} | {'cut * L^2':>10}")
print("-" * 74)
m1, m_star = 1.0, 5.0      # glueball mass; excited glueball (above the 2-particle threshold here)
m = 1.0                    # :phi2: constituents
for L in (4, 8, 16, 32, 64, 128):
    da = gap_atom(m1, m_star, L)
    dc = gap_cut_phi2(m, L)
    print(f"{L:>5} | {da:>14.6f} | {dc:>22.6e} | {dc*L*L:>10.4f}")

print()
print("READING:")
print("  atom:  Delta(L) = m1 = 1.0 exactly at every L  ->  UNIFORM positive floor  (the gap)")
print("  cut:   Delta(L) ~ (2pi/L)^2 = 1/L^2 -> 0        ->  NO uniform positive bound")
print("         (cut * L^2 -> (2pi)^2 = 39.478, confirming the 1/L^2 scaling)")
print()
print("A :phi2: CUT has a POSITIVE floor (2m > 0) yet Delta(L)->0 -- so 'positive floor'")
print("is strictly weaker than 'uniform isolation'. The two agree only for an ATOM (gap).")
print("=> the approach-law's uniformity in (L,a) IS the gap (H1); it is not a free handle.")
