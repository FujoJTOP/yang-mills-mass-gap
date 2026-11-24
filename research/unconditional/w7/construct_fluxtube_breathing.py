"""
construct_fluxtube_breathing.py -- CONSTRUCT-FLUXTUBE-BREATHING witness.

Model the neutral 0++ glueball as a closed flux loop (rigid string:
E(R)=2*pi*sigma*R + pi*kappa/R) and quantize the transveral RADIAL (breathing)
collective mode R(t)=R*+dR cos(wt).  Compute m_breath, compare to the
rigid-string ground state 2*pi*sqrt(2*kappa) (MASTERFIELD) and to lattice
m_0++/sqrt(sigma) ~ 3.9, and locate the transverse-phonon continuum on the loop.

NOT a Yang-Mills proof; a model computation ([C], sigma set to 1).
Run: py -3.12 -I construct_fluxtube_breathing.py
"""
import numpy as np

def analyze(kappa, sigma=1.0):
    R = np.sqrt(kappa / (2 * sigma))                 # interior minimum of E(R)
    m_ground = 2 * np.pi * np.sqrt(2 * sigma * kappa)  # = E(R*), MASTERFIELD 0++
    Epp = 2 * np.pi * kappa / R**3                    # d2E/dR2 at R*
    M_eff = 2 * np.pi * sigma * R                     # ring inertia (tube mass)
    m_breath = np.sqrt(Epp / M_eff)                   # == 2*sqrt(sigma/kappa)
    w1 = 1.0 / R                                      # one-phonon quantum (circle)
    w2 = 2.0 / R                                      # two-phonon threshold
    return dict(kappa=kappa, R=R, m_ground=m_ground, m_breath=m_breath,
                ratio=m_breath / m_ground, w1=w1, w2=w2,
                above_single=m_breath > w1, below_two=m_breath < w2)

print("=" * 84)
print("Closed rigid-string flux loop: ground state, radial BREATHING mode, phonon continuum")
print("=" * 84)
print(f"{'kappa':>7}{'R*':>8}{'m_ground':>10}{'m_breath':>10}{'ratio':>8}"
      f"{'1/R*':>8}{'2/R*':>8}{'>1ph':>6}{'<2ph':>6}")
for kappa in [0.138, 0.16, 0.193, 0.25, 0.5]:
    d = analyze(kappa)
    print(f"{d['kappa']:>7.3f}{d['R']:>8.4f}{d['m_ground']:>10.4f}{d['m_breath']:>10.4f}"
          f"{d['ratio']:>8.3f}{d['w1']:>8.3f}{d['w2']:>8.3f}"
          f"{str(d['above_single']):>6}{str(d['below_two']):>6}")

print()
print("Closed forms (sigma=1):  m_ground = 2*pi*sqrt(2*kappa)   m_breath = 2/sqrt(kappa)")
print("ratio m_breath/m_ground = 1/(pi*sqrt(2*kappa))")
print("continuum on the circle: one-phonon 1/R*, two-phonon threshold 2/R*")
print()
print("Atom bookkeeping at kappa=0.193 (MASTERFIELD fit -> m_ground/sqrt(sigma) = 3.90 = lattice):")
d = analyze(0.193)
print(f"  m_ground/sqrt(sigma) = {d['m_ground']:.3f}   (lattice 0++ = 3.9)")
print(f"  m_breath/sqrt(sigma) = {d['m_breath']:.3f}   (an EXCITED 0++; overshoots 3.9 by {100*(d['m_breath']/3.9-1):.1f}%)")
print(f"  m_breath sits in ({d['w1']:.2f}, {d['w2']:.2f}) = between 1-phonon and 2-phonon threshold")
print(f"  -> cannot decay into TWO phonons (m_breath < 2/R*); single-phonon decay needs J=0 selection")
print()
print("Local (infinite-tube) sausage mode ~ m_H; dual-GL m_H/sqrt(sigma) (CONSTRUCT-DUAL-LATTICE):")
for kGL, mr in [(2, 1.36), (5, 2.22), (10, 3.72)]:
    print(f"  kappa_GL={kGL:>2}: m_H/sqrt(sigma) ~ {mr:.2f}")
