"""
pnt_dual_condensate.py -- witness for PNT-DUAL-CONDENSATE.
Dual Abelian-Higgs ('t Hooft-Mandelstam): can the dual-Higgs mass be LIFTED to the KL atom?

Conventions that reproduce ROOF-DUALSUPERCONDUCTOR's ratios (g~=1):
  m_V = g~ v           (dual photon, 1^{--})
  m_H = sqrt(2 lam) v = kappa_GL * m_V   (dual Higgs magnitude, 0^{++})
  London tension  sigma = pi v^2 ln(kappa_GL)
Perturbative dual-AHM 0^{++} content of the ELECTRIC operator tr F^2 ~ (dual F)^2:
  - dual-photon pair: threshold 2 m_V   (a CUT)
  - dual-Higgs:       m_H              (a POLE only if tr F^2 couples to |phi|^2 -- an extra,
                                        non-minimal coupling the minimal dual AHM does NOT fix)
KL atom f_1^2 = residue at the pole. A CUT is not an atom (that is the :phi^2: failure mode).

Run: py -3.12 -I pnt_dual_condensate.py
"""
import numpy as np

LATTICE_0PP = 3.9   # Morningstar-Peardon: lightest 0++ / sqrt(sigma)

def mV_over_rtsig(kappa):        # with g~=1: m_V/sqrt(sigma) = 1/sqrt(pi ln kappa)
    return 1.0 / np.sqrt(np.pi * np.log(kappa))

print("=" * 80)
print("(1) Perturbative dual-AHM 0++ masses vs string tension and the lattice glueball")
print("=" * 80)
print(f"{'kappa_GL':>9} {'m_H/sqrt(sig)':>14} {'2m_V/sqrt(sig)':>15} {'2m_V/m_H':>9}  verdict")
for kappa in [1.5, 2.0, 5.0, 10.0, 12.0, 20.0, 40.0]:
    r = mV_over_rtsig(kappa)
    mH = kappa * r
    tV = 2.0 * r
    tags = []
    if mH < LATTICE_0PP: tags.append("m_H too light")
    if tV < LATTICE_0PP: tags.append("2m_V too light")
    print(f"{kappa:9.2f} {mH:14.3f} {tV:15.3f} {2.0/kappa:9.3f}  {', '.join(tags)}")
print(f"\n  lattice 0++/sqrt(sigma) ~ {LATTICE_0PP}")
print("  -> for every kappa, at least one perturbative dual-AHM 0++ mode is BELOW the")
print("     glueball; pushing m_H up to 3.9 needs kappa~11, where 2m_V/sqrt(sig)~0.18 is")
print("     far below. So the physical 0++ is NOT a perturbative dual-AHM fluctuation.")

print("\n" + "=" * 80)
print("(2) Atom vs cut: the tr F^2 (electric) 0++ spectral density")
print("=" * 80)
print("""
  rho_trF2(s) = [dual-photon-pair CUT from s=4 m_V^2]                       (always)
              + f_1^2 * delta(s - m_H^2)         (POLE only if tr F^2 <-> |phi|^2 is switched on)

  * The KL atom f_1^2 is the RESIDUE of a POLE.
  * The minimal dual AHM fixes the MASSES m_V, m_H but NOT the coupling tr F^2 <-> |phi|^2,
    hence NOT f_1^2:  any  0 <= f_1^2 < inf  gives the same m_H.
  * A CUT at 2 m_V carries f_1^2 = 0 in the atom sense: gapped, atomless  ==  the :phi^2: adversary.
  * Free massless (v=0): m_V=m_H=0, no pole, f_1^2=0  ->  acid test returns 0 (passes).
""")
print("  => the dual AHM delivers a renormalon-free MASS but not the ATOM.")
print("     Lifting to f_1^2>0 is an independent non-minimal overlap, not fixed here.")

print("\n" + "=" * 80)
print("(3) The likely physical 0++: a vortex RING (collective), atom overlap = open")
print("=" * 80)
for kappa in [2.0, 5.0, 12.0]:
    r = mV_over_rtsig(kappa)
    for c in [0.5, 0.62]:
        print(f"  kappa={kappa:>5.1f}: vortex-ring mass/sqrt(sig) ~ 2*pi*{c} = {2*np.pi*c:.3f}"
              f"   (m_H/sqrt(sig)={kappa*r:.3f}, 2m_V/sqrt(sig)={2*r:.3f})")
print("  -> the ratio-3.9 glueball is a vortex RING. Its tr F^2 overlap <0|tr F^2|ring>")
print("     is a nonperturbative matrix element NOT fixed by dual-AHM semiclassics.")

print("\n" + "=" * 80)
print("BOTTOM LINE: monopole condensation -> renormalon-free MASS [C];")
print("             lifting to the KL ATOM f_1^2>0 is NOT delivered [O]: the atom is a")
print("             nonperturbative tr F^2 <-> vortex-ring overlap. Same wall (H1), dual coordinate.")
print("=" * 80)
