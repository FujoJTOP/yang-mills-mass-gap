"""
roof_defect_cft.py -- witness for research/unconditional/ROOF-DEFECT-CFT.md

Defect CFT (conformal / boundary / defect bootstrap) as a route to the 0++ mass gap.
Three structural checks; cheap (numpy + scipy.special.k0), no lattice, no MC.

 1. REGIME. Defect CFT power-law one-point functions require the BULK to be conformal
    (m=0). At m>0 the decay is exponential with rate m. We illustrate with the 2D
    propagator G(m;r): m=0 gives -log r (power law), m>0 gives K_0(mr) ~ e^{-mr}.
    => "defect CFT" presupposes the gapless regime the mass gap is about (ACID-FAIL
       for a gapped bulk, which is exactly 4D pure YM in the IR).

 2. CHANNEL. The displacement operator of a p-dimensional conformal defect in D dimensions
    has PROTECTED dimension Delta_D = D - p + 1 = codim + 1 (transverse conformal symmetry;
    Billo-Goncalves-Lauria-Meineri, JHEP 04 (2016) 091, arXiv:1601.02883). It is
    renormalon-free and gap-blind (independent of the bulk mass/coupling).

 3. 0-FORM NEUTRALITY. The 't Hooft phase is exp(2*pi*i*lk/N); a local (0-form) operator
    is a point (0-cycle) and links nothing, so lk(point,loop)=0 and the phase is 1. The
    0++ operator tr F^2 is center-neutral. The 't Hooft loop is 1-form charged.

Run:  py -3.12 -I research/unconditional/w7/roof_defect_cft.py
"""
import numpy as np
from scipy.special import k0

def G2d(m, r):
    r = np.asarray(r, float)
    if m == 0.0:
        return -np.log(r) / (2.0 * np.pi)       # 2D massless propagator (up to additive const)
    return k0(m * r) / (2.0 * np.pi)            # 2D massive propagator ~ e^{-mr}/sqrt(r)

print("=" * 72)
print("1. REGIME -- defect CFT power law requires m=0 (bulk conformality)")
print("=" * 72)
r = np.array([4.0, 8.0, 12.0, 16.0, 20.0])
for m in [0.0, 0.5, 1.0]:
    g = np.abs(G2d(m, r))
    m_eff = -np.gradient(np.log(g), r)         # effective exponential rate -d log G / dr
    print(f"  m={m:3.1f}  effective rate -dlogG/dr = {np.round(m_eff, 4)}")
print("  --> m>0: effective rate -> m (exponential decay). m=0: power law (-log r), rate ~0.")
print("      Defect-CFT scaling exists ONLY at m=0, i.e. only in the GAUPLESS regime.")
print("      For the gapped IR of 4D pure YM there is no conformal defect CFT to bootstrap.")

print()
print("=" * 72)
print("2. CHANNEL -- displacement dimension Delta_D = D-p+1 = codim+1 (protected, gap-blind)")
print("=" * 72)
# (defect name, p = worldvolume dim, D = spacetime dim); Delta_D = D - p + 1
rows = [("line", 1, 4), ("surface", 2, 4), ("line", 1, 3), ("surface", 2, 3)]
for name, p, D in rows:
    codim = D - p
    print(f"  {name:8s} D={D}: p={p}, codim={codim}, Delta_D = D-p+1 = {D-p+1}   (PROTECTED)")
print("  --> the 't Hooft loop (line, p=1, D=4) has Delta_D = 4.  defect-CFT data are")
print("      renormalon-free but GAP-BLIND: Delta_D does not depend on the bulk mass or")
print("      coupling, so it cannot carry the 0++ gap.")

print()
print("=" * 72)
print("3. 0-FORM NEUTRALITY -- the 0++ operator is 't Hooft-neutral")
print("=" * 72)
lk_point = 0                                    # a point links no loop
for N in [2, 3]:
    phase = np.exp(2j * np.pi * lk_point / N)
    print(f"  N={N}: exp(2*pi*i*lk(point,loop)/N) = {phase.real:+.3f}{phase.imag:+.3f}j  (lk=0)")
print("  --> a local (0-form) operator has phase 1: it sits in the 1-form-NEUTRAL sector.")
print("      The 't Hooft loop B(C) is 1-FORM charged. In the defect OPE of tr F^2 on B(C),")
print("      the operators that appear are DEFECT operators -- defect dimensions, not the")
print("      BULK 0++ gap.  Channel mismatch persists (same as BASE-ORDER-PARAMETER).")

print()
print("VERDICT: defect CFT is renormalon-free but (a) requires m=0 [= gaplessness, the")
print("negation of the gap] and (b) yields 1-form / defect-local data, not the bulk 0++")
print("gap.  It inherits BOTH the regime and the channel obstruction.  Not the roof.")
