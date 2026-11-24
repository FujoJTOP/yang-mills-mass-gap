"""
construct_holo_fluxtube.py -- witness for CONSTRUCT-HOLO-FLUXTUBE.md
Holographic flux tube -> 0++ glueball (hard-wall AdS5 spectrum + string tension at the IR tip).
Lightweight: scipy only, no lattice MC.
Run:  py -3.12 -I construct_holo_fluxtube.py

What is computed (all [T] within the stated model; no YM claim):
  (1) Hard-wall AdS5, dilaton 0++ (nu=2): m_n = j_{2,n}/z0; discrete zeros => DISCRETE spectrum (atom).
  (2) String tension at the IR wall (Nambu-Goto string at z=z0): sigma = L^2/(2*pi*alpha' * z0^2).
  (3) Ratio  m_0++/sqrt(sigma) = j_{2,1}/sqrt(lam'),  lam' = L^2/(2*pi*alpha').
      -> matches lattice 3.9 ONLY at lam' ~ 1.7 (outside the large-lam SUGRA regime).
  (4) Witten-cigar structural scaling: both m and sqrt(sigma) are set by the SAME IR scale M_KK,
      so the ratio is O(1) (Brower-Mathur-Tan ~3.5-4), independent of lam -- the genuine positive.
"""
import numpy as np
from scipy.special import jn_zeros

print("=" * 74)
print("CONSTRUCT-HOLO-FLUXTUBE | holographic flux tube -> 0++ glueball")
print("=" * 74)

# ---------------------------------------------------------------------------
# (1) Hard-wall AdS5: 0++ = lowest bulk scalar KK mode (dilaton, nu=2).
#     m_n = j_{nu,n}/z0  (j_{2,n} = zeros of the Bessel function J_2)
# ---------------------------------------------------------------------------
nu = 2
j2 = jn_zeros(nu, 4)
print("\n(1) Hard-wall AdS5 scalar spectrum (nu=2, dilaton 0++):  m_n = j_{2,n}/z0")
print("    j_{2,n} =", np.round(j2, 5))
print("    m_n/m_1 =", np.round(j2 / j2[0], 4),
      "  (lattice glueball tower 1 : ~1.6 : ~2.3 : ~2.9)")

# ---------------------------------------------------------------------------
# (2)-(3) String tension at the IR wall and the ratio m_0++/sqrt(sigma).
#     Nambu-Goto string sitting at z=z0: induced metric (L^2/z0^2) diag(-1,1)
#       sigma = L^2/(2*pi*alpha' * z0^2),  lam' := L^2/(2*pi*alpha')
#     m_0++ = j_{2,1}/z0  =>  m/sqrt(sigma) = j_{2,1}/sqrt(lam')
# ---------------------------------------------------------------------------
print("\n(2)-(3) String tension at the IR wall: sigma = lam'/z0^2,  lam' = L^2/(2 pi alpha')")
print("    RATIO  m_0++/sqrt(sigma) = j_{2,1}/sqrt(lam')   (a pure number, set by geometry)")
target = 3.9  # lattice ~3.9
lam_star = (j2[0] / target) ** 2
print(f"    to match lattice {target}:  lam' = (j_2,1/{target})^2 = {lam_star:.3f}")
print("    lambda'  ->   m_0++/sqrt(sigma) = j_2,1/sqrt(lam'):")
for lam in [1.0, 1.734, 3.0, 10.0, 30.0, 100.0]:
    print(f"        lam'={lam:7.3f}  ->  {j2[0]/np.sqrt(lam):7.3f}")
print("    >>> SUGRA requires lam' = L^2/(2 pi alpha') LARGE (>>1).")
print("    >>> In that controlled regime the hard-wall ratio -> 0, NOT O(1).")
print("    >>> So the hard-wall m/sqrt(sigma) matches 3.9 only at lam'~1.7 (uncontrolled) = a FIT, not a prediction.")

# ---------------------------------------------------------------------------
# (4) Witten D4 cigar: structural O(1) ratio (the genuine positive).
#     In the cigar both the glueball mass and the string tension are set by the
#     same IR scale M_KK = (2/3)(u_L/R^3)^{1/2}; hence m/sqrt(sigma) is a pure
#     O(1) number independent of lam.  Brower-Mathur-Tan (hep-th/0003115) find
#     the 0++ lightest with m_0++/sqrt(sigma) ~ 3.5-4, matching lattice ~3.9.
#     We do NOT recompute the full SUGRA reduction here; we record the scaling.
# ---------------------------------------------------------------------------
print("\n(4) Witten-cigar scaling (structural, [C], cited not recomputed):")
print("    m_0++ ~ c * M_KK ,  sqrt(sigma) ~ c' * M_KK   =>   m/sqrt(sigma) = c/c' = O(1).")
print("    Brower-Mathur-Tan (hep-th/0003115): m_0++/sqrt(sigma) ~ 3.5-4  (lattice ~3.9).")
print("    ordering m(0++) < m(2++) < m(1++) reproduced  [T] within the model.")

# ---------------------------------------------------------------------------
# (5) Atom vs continuum diagnostic.
# ---------------------------------------------------------------------------
print("\n(5) ATOM vs CONTINUUM:")
print("    hard-wall / cigar KK spectrum = discrete Bessel/KK zeros => the 0++ is an ISOLATED level")
print("    (an ATOM, a normalizable bulk mode <-> a CFT state), NOT a continuum.")
print("    BUT the discreteness comes from the IR wall / cigar = the dual of CONFINEMENT.")
print("    The atom is the dual of the IR geometry; the IR geometry is imported, not derived = H1.")

print("\n" + "=" * 74)
print("VERDICT: holography gives an ATOM (discrete 0++, right O(1) ratio ~3.5-4) but by")
print("ASSUMING the confining bulk geometry (IR wall / cigar). The atom is the dual of H1.")
print("Clay 0%; no QED.")
print("=" * 74)
