"""
kappa_overdetermined.py -- over-determine kappa (rigid-string) against the lattice glueball spectrum.
Representative SU(3) lattice values (Morningstar-Peardon 1999, hep-lat/9901004), sqrt(sigma)~440 MeV.
NOT a lattice computation. Tests whether the ONE-parameter rigid string fits ALL states.
"""
import numpy as np

# lattice ratios m/sqrt(sigma) (representative literature values; 0++* has large error)
m0  = 3.93   # 0++   ~1730 MeV
m2  = 5.45   # 2++   ~2400 MeV
m0s = 6.07   # 0++*  ~2670 MeV (excited scalar, large error)

print("lattice ratios:  0++ =%.3f  2++ =%.3f  0++* =%.3f  (units sqrt(sigma))"%(m0,m2,m0s))
print("lattice splits:  0++*/0++ =%.3f   2++/0++ =%.3f"%(m0s/m0, m2/m0))

# --- fix kappa from the lightest 0++ (the model's defining relation): m0 = 2*pi*sqrt(2 kappa)
kappa = m0**2/(8*np.pi**2)
print("\nkappa from m_0++ = (m0)^2/(8 pi^2) = %.4f   [ONE datum -> kappa]"%kappa)

# --- radial/breathing mode: m_breath/sqrt(sigma) = 2/sqrt(kappa)
m_breath = 2/np.sqrt(kappa)
print("radial (breathing) 0++: m/sqrt(sigma)=%.3f  ratio/m0=%.3f"%(m_breath, m_breath/m0))

# --- closed-string oscillator tower: E_N = m0 + (N+Nbar-1/6)/R*, R*=sqrt(kappa/2)/sqrt(sigma)
# first level N+Nbar=1: Delta = (5/6)/R* = (5/6) sqrt(2/kappa) sqrt(sigma)
m_osc = m0 + (5.0/6.0)*np.sqrt(2.0/kappa)
print("closed-string osc level (N+Nbar=1): m/sqrt(sigma)=%.3f  ratio/m0=%.3f  [0++* AND 2++ degenerate]"%(m_osc, m_osc/m0))

print("\n--- tensions vs lattice ---")
print("radial:  %.3f vs 0++* %.3f  -> %+.1f%%"%(m_breath/m0, m0s/m0, 100*(m_breath/m0-m0s/m0)/(m0s/m0)))
print("osc:     %.3f vs 0++* %.3f  -> %+.1f%%"%(m_osc/m0, m0s/m0, 100*(m_osc/m0-m0s/m0)/(m0s/m0)))
print("osc degenerate 2++=0++* at %.3f; lattice 2++ = %.3f (split %.3f) -> osc over by %+.1f%%"%(m_osc,m2,m0s-m2,100*(m_osc-m2)/m2))
print("\nVERDICT: one datum fixes kappa; the two natural 'next 0++' predictions (radial 1.15, osc 1.68)")
print("bracket the lattice 1.54 but neither matches, and the osc tower predicts 2++=0++* (lattice splits them).")
