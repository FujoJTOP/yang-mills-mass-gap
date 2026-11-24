"""
roof_holography.py -- witness for ROOF-HOLOGRAPHY (2026-11-24).

Hard-wall + soft-wall glueball spectra: shows the 0++ mass is a NUMBER times
an INPUT IR scale (z0 or dilaton slope c), i.e. the renormalon/dimensional-
transmutation ambiguity is RELOCATED into the IR boundary condition, not
removed. NOT a Yang-Mills calculation.

Run: py -3.12 -I research/unconditional/w7/roof_holography.py
"""
import numpy as np
from scipy.special import jn_zeros

print("=" * 72)
print("HARD WALL (Polchinski-Strassler): 0++ scalar glueball, m_n = j_{2,n}/z0")
print("  wave eqn  d_z(z^-3 d_z phi) + m^2 z^-5 phi = 0  ->  phi = z^2 J_2(m z)")
print("  IR BC  J_2(m z0) = 0  ->  m_n = j_{2,n}/z0")
print("=" * 72)
j = jn_zeros(2, 4)
print("  J_2 zeros j_{2,n}      =", np.round(j, 4))
print("  m_n z0                 =", np.round(j, 4), "   (m_0++ z0 = %.4f)" % j[0])
print("  mass ratios m_n/m_1    =", np.round(j / j[0], 4),
      "   (lattice 0++*/0++ ~ 1.6-1.8)")

print()
print("=" * 72)
print("SOFT WALL (Karch-Katz-Son-Stephanov): dilaton Phi=c z^2, m_n^2=4c(n+1)")
print("=" * 72)
c = 1.0
n = np.arange(4)
m2 = 4 * c * (n + 1)
print("  m_n^2/(4c)             =", (m2 / (4*c)).astype(int))
print("  m_n/sqrt(4c)           =", np.round(np.sqrt(m2 / (4*c)), 4), "  (linear Regge: m_n^2 ~ n)")
print("  m_{n+1}/m_n -> 1       :", np.round(np.sqrt((n[1:-1] + 2) / (n[1:-1] + 1)), 4))

print()
print("=" * 72)
print("SCALE-RELOCATION: the gap is an overall IR scale (an INPUT)")
print("=" * 72)
print("  hard wall : m_0++ = %.4f / z0        (z0 = IR wall position, input)" % j[0])
print("  soft wall : m_0++ = %.4f * sqrt(c)   (c  = dilaton slope, input)" % np.sqrt(8*c))
print("  => m_0++ = (dimensionless number) x (INPUT IR scale).")
print("  => m_0++/sqrt(sigma) is scale-INDEPENDENT (a genuine predicted ratio)")
print("     but the ABSOLUTE gap is NOT predicted:")
print("     dim.-transmutation ambiguity is RELOCATED into z0 / c, not removed.")

print()
print("=" * 72)
print("CONSISTENCY with strong-coupling bridge  m_0++ >= 4 sqrt(sigma)")
print("=" * 72)
print("  A scale-free model predicts only RATIOS. If m_0++ = k/z0 and")
print("  sqrt(sigma)=k'/z0, then m_0++/sqrt(sigma)=k/k' (model-dependent,")
print("  holography ~ O(3-4); lattice ~3.7). The strong-coupling bound >=4")
print("  is a u<=e^-1 statement, NOT contradicted at the crossover.")
print()
print("VERDICT: holography predicts definite FORMULAS and RATIOS, but the")
print("absolute gap is the input scale z0/c. It RELOCATES H1 (dim. transm.)")
print("into the IR boundary condition; it does not derive it. And it is NOT")
print("rigorous: no proven duality between pure YM and these geometries.")
