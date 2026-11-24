"""
impossible_susy_decouple.py -- witness for IMPOSSIBLE-SUSY-DECOUPLE.md

N=1 SYM -> pure YM via gluino decoupling.  Illustrates the STRUCTURAL point:
the exact holomorphic scale matching transfers the SCALE, but the matched
pure-YM IR sits at alpha_s ~ O(1) (crossover), where the SUSY/holomorphic
protection that made SYM exact is absent.  So the exact SYM gap is not
transferred; the pure-YM gap at the matched scale is H1.

Schematic 1-loop running + Konishi/NSVZ shift.  Not a YM calculation.
Run:  py -3.12 -I research/unconditional/w7/impossible_susy_decouple.py
"""
import numpy as np

pi = np.pi
N  = 2
b_SYM = 3 * N           # NSVZ one-loop coefficient (N=1 SYM)
b_YM  = 11 * N / 3.0    # pure-YM one-loop coefficient

Lam_SYM = 1.0
mu0     = 10.0 * Lam_SYM       # high reference scale
m_glu   = 100.0                # heavy gluino (decoupling)
delta   = N / (8 * pi**2)      # Konishi/NSVZ coefficient

# SYM coupling at mu0 (1-loop, illustrative)
g2_SYM = 8 * pi**2 / (b_SYM * np.log(mu0 / Lam_SYM))

# exact holomorphic matching: 1/g_YM^2(mu0) = 1/g_SYM^2(mu0) + (N/8pi^2) ln(m/mu0)
inv_g2_YM = 1.0 / g2_SYM + delta * np.log(m_glu / mu0)

# pure-YM dynamical scale from one-loop definition: Lambda = mu * exp(-8pi^2/(b g^2(mu)))
g2_YM_mu0 = 1.0 / inv_g2_YM
Lam_YM = mu0 * np.exp(-8 * pi**2 / (b_YM * g2_YM_mu0))

print("SU(2), N=1 SYM (b=%.2f) decoupled to pure YM (b=%.3f)" % (b_SYM, b_YM))
print("  Lambda_SYM = %.3f ; gluino mass m = %.1f" % (Lam_SYM, m_glu))
print("  g2_SYM(mu0)               = %.4f" % g2_SYM)
print("  1/g2_YM(mu0) [matched]    = %.4f  (g2_YM=%.4f)" % (inv_g2_YM, 1/inv_g2_YM))
print("  matched pure-YM scale Lambda_YM = %.4f  (exact holomorphic matching)"
      % Lam_YM)
print()
print("  pure-YM coupling near the matched scale (crossover):")
for f in (2.0, 3.0, 5.0):
    mu = f * Lam_YM
    g2 = 8 * pi**2 / (b_YM * np.log(mu / Lam_YM))
    print("    g2_YM(%.1f * Lambda_YM) = %.3f" % (f, g2))
print()
print("=> The SCALE matches exactly (holomorphic).  But at Lambda_YM the")
print("   pure-YM coupling is O(1) -- the crossover -- where holomorphy (the")
print("   source of the SYM exactness) is ABSENT.  The exact SYM gap is NOT")
print("   transferred; the pure-YM gap at the matched scale is H1.")
print()
print("Structural note: the SYM beta is a closed all-orders form (NSVZ);")
print("pure YM has no such closed form -> 'exact' is a SUSY property pure YM lacks.")
