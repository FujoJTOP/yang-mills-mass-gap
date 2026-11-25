"""
z1_holo.py -- witness for research/v3/Z1-HOLO.md.

Question: can a RIGOROUS statement "5D normal mode => 4D mass gap" be extracted
for pure SU(N) YM?  Test on solvable holographic models: hard wall (AdS5) and
BTZ (AdS3).  Acid test: the free massless field must come out GAUPLESS.

Structure of the test:
  (A) HARD WALL, AdS5.  The gauge-invariant 0++ bulk mode solves the Sturm-
      Liouville (SL) equation  phi'' - (3/z)phi' + m^2 phi = 0  (= d_z(z^-3 d_z phi)
      + m^2 z^-3 phi = 0) with IR BC  phi(z0)=0.  Solution phi = z^2 J_2(m z),
      so m_n = j_{2,n}/z0.  Verify the SL residual ~ 0, the Bessel zeros, and
      that the "gap" m1 = j_{2,1}/z0 -> 0 as z0 -> infinity (wall removed = pure
      AdS = free/CFT).  => the gap is EXACTLY the input IR scale z0.
  (B) PURE AdS5 (acid test / free massless dual).  With NO IR boundary condition
      there is no quantization: every m>0 is admitted, inf spectrum = 0, the
      boundary channel is a gapless CFT.  Detector returns NO gap for the free
      massless field: not vacuous.  But (A) shows the only place it returns a
      gap is where a finite z0 was INSERTED.
  (C) BTZ, AdS3.  One scale: horizon r+ (temperature T = r+/(2 pi L^2)).  Scalar
      QNM (Son-Starinets) imaginary parts = (dimensionless O(n+Delta)) x (2 pi T);
      ratios are T-independent, so every frequency ~ r+/L^2 -> 0 as r+ -> 0,
      where BTZ -> pure AdS3 (a gapless 2D CFT).  Same structure as (A): the
      "gap" is the compactification/thermal scale, an input.

numpy/scipy only, deterministic, < 1 s.  Sets NO claim of a YM theorem.
Run: py -3.12 -I research/v3/w7/z1_holo.py
"""
import numpy as np
from scipy.special import jn_zeros, jv

np.set_printoptions(precision=5, suppress=True)
BAR = "=" * 74

# ---------------------------------------------------------------------------
# (A) HARD WALL, AdS5: SL operator, spectrum, gap == input scale z0
# ---------------------------------------------------------------------------
print(BAR)
print("(A) HARD WALL AdS5 -- SL eqn  phi'' - (3/z)phi' + m^2 phi = 0,")
print("    phi = z^2 J_2(m z);  IR BC J_2(m z0)=0  =>  m_n = j_{2,n}/z0")
print(BAR)

# verify the claimed SL solution numerically (the prefactor matters: z^-3, not z^-5)
m = 1.7
zz = np.linspace(0.3, 3.0, 20001)
phi = zz**2 * jv(2, m * zz)
d1 = np.gradient(phi, zz)
d2 = np.gradient(d1, zz)
res = d2[50:-50] - (3.0 / zz[50:-50]) * d1[50:-50] + m * m * phi[50:-50]
print(f"  max |phi'' - 3/z phi' + m^2 phi|  = {np.max(np.abs(res)):.2e}   (SL solution OK)")

j = jn_zeros(2, 5)
print(f"  J_2 zeros  j_{{2,n}}        = {j}")
print(f"  mass ratios m_n/m_1      = {j/j[0]}")
print("  (these ratios are the model's real success: lattice 0++*/0++ ~ 1.6-1.8)")

print("\n  gap m1 = j_{2,1}/z0 vs the IR wall position z0  (z0 is an INPUT):")
for z0 in [0.5, 1.0, 5.0, 50.0, 1e3, 1e6]:
    print(f"    z0 = {z0:9.1f}   m1 = {j[0]/z0:.4e}")
print("  z0 -> infinity (= remove wall = pure AdS = free/CFT)  =>  m1 -> 0.")
print("  => inf over all z0 of the 'gap' is 0: the gap is not produced by the")
print("     SL operator, it is the scale z0 that was inserted.")

# ---------------------------------------------------------------------------
# (B) PURE AdS5: acid test -- free massless field must be gapless
# ---------------------------------------------------------------------------
print()
print(BAR)
print("(B) PURE AdS5 (no IR wall) -- acid test on the free massless field")
print(BAR)
print("  No IR boundary condition => no quantization: the normalizable branch")
print("  phi ~ z^2 J_2(m z) exists for EVERY m>0 (Bessel oscillates for all z),")
print("  so the boundary two-point function has spectral weight from s=0 down.")
print("  inf spec(L) = 0.  The dual is a gapless CFT (free massless scalar).")
print("  => detector returns NO gap for the free massless field: NOT vacuous.")
print("  BUT: (A) shows the ONLY region where it returns a gap is finite z0,")
print("       i.e. where the wall (= the gap) was put in by hand.")
# numeric: confirm J_2 oscillates (infinitely many zeros) => continuum of m admitted
extra = jn_zeros(2, 12)
print(f"  J_2 has infinitely many zeros (first 12): {extra[:6]} ...")
print("  each zero is a wall position/mass pair; with no wall, m is continuous.")

# ---------------------------------------------------------------------------
# (C) BTZ, AdS3: one scale r+ (temperature); gap ~ T -> 0 as r+ -> 0
# ---------------------------------------------------------------------------
print()
print(BAR)
print("(C) BTZ AdS3 -- one scale r+ (T = r+/(2 pi L^2)); scalar QNM")
print(BAR)
print("  AdS3 scalar mass-dimension relation:  Delta(Delta-2) = m^2 L^2.")
print("  Non-rotating BTZ scalar QNM (Cardoso-Lemos; AdS/CFT recipe Son-Starinets):")
print("    Im(omega_n) = -(r+/L^2)(2n + Delta') = -(2 pi T)(2n + Delta'),  Re = +-k")
print("    (Delta' a field-dependent O(1) offset; only linear-in-n, ~T scaling matters).")
T = 1.0
for Dp in [2.0, 3.0]:               # e.g. m^2 L^2 = 0 -> Delta=2 ; Delta'=3 -> m^2 L^2 = 3
    n = np.arange(4)
    im = -(2 * np.pi * T) * (2 * n + Dp)
    print(f"    Delta'={Dp}:  Im(omega)/(2 pi T) = {(im/(2*np.pi*T))}  (T-independent)")
print("  ratios Im(omega)/(2 pi T) are T-INDEPENDENT => the whole tower scales with T.")
print("  r+ -> 0 (T -> 0): every frequency -> 0 and BTZ -> pure AdS3, a gapless 2D CFT.")
print("  => BTZ's discrete 'gap' is the horizon/thermal scale r+ = an input,")
print("     not a dynamically generated zero-temperature mass gap.")

print()
print(BAR)
print("VERDICT (no YM theorem; no QED)")
print(BAR)
print("  [T] within a fixed model: SL theorem + GKPW lift bulk normal modes to")
print("      boundary KL atoms; a bulk spectral gap Delta_5>0 => boundary gap if")
print("      the SL operator is discrete (an IR cap exists).  This is real but is")
print("      a statement about the MODEL's SL operator.")
print("  The '5D normal mode => 4D mass gap' chain FAILS for pure YM because:")
print("    * the IR cap that makes the spectrum discrete IS the confining scale")
print("      (gap by hand) -- pure AdS / pure gravity is gapless (B, acid);")
print("    * the bulk must be dual to 4D pure YM -- no such duality (H1);")
print("    * gap = wall/thermal scale (z0, r+) -> 0 under decompactification (A,C).")
print("  Net: reduces to H1.  Clay 0%.")
