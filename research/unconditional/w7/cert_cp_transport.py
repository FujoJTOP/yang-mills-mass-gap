"""
cert_cp_transport.py -- witness for CERT-CP-TRANSPORT
Certificate attempt on the ATOM half via CP / Luscher-flow residue transport:
    f1^2 = Z_O(t)^2 f1^2_strong
Question: does the flowed residue survive >0 as t->inf, beta->2.5?

Run: py -3.12 -I research/unconditional/w7/cert_cp_transport.py
"""
import numpy as np
np.set_printoptions(precision=4, suppress=True)

print("=" * 76)
print("CERT-CP-TRANSPORT: can CP transport give the physical atom at beta~2.5?")
print("=" * 76)

f1sq = 0.5   # atom residue at mass m1
m1 = 1.0

# ----------------------------------------------------------------------
# (A) Luscher flow: flowed residue decays absolutely; atom only enhances RELATIVELY
#     flow factor per state e^{-t m^2} on the amplitude, so residue e^{-2 t m^2}
# ----------------------------------------------------------------------
print("\n(A) LUSCHER FLOW: flowed residue f1^2(t) = f1^2 * e^{-2 t m^2}")
m2 = np.array([2.0, 2.5, 3.0, 4.0, 5.0])
f2 = np.array([0.20, 0.10, 0.10, 0.06, 0.04])
for t in [0.0, 1.0, 5.0, 10.0]:
    f1t = f1sq * np.exp(-2 * t * m1 * m1)
    tot = f1t + float(np.sum(f2 * np.exp(-2 * t * m2 * m2)))
    print(f"    t={t:5.1f}: f1^2(t)={f1t:.3e}   relative share f1^2(t)/total={f1t/tot:.4f}")
print("    -> ABSOLUTE weight -> 0 as t->inf; RELATIVE share -> 1 (atom enhanced).")
print("       The flow is DOWNSTREAM of f1^2(0): it cannot create the atom.")
print("       'survives t->inf' is convention-dependent (relative yes, absolute no).")

# ----------------------------------------------------------------------
# (B) The multiplicative transport is a normalization identity
# ----------------------------------------------------------------------
print("\n(B) CP TRANSPORT f1^2_phys = Z_O^2 f1^2_strong")
print("    P3: f1^2_strong>0 (beta<1/48). A Hermitian operator has real Z_O, so Z_O^2>=0.")
print("    Given Z_O!=0 and the SAME atom, f1^2_phys>0 follows -- but the transport is a")
print("    finite-(L,a) normalization identity; it is vacuous for the continuum content (see C).")
for ZO in [0.5, 1.0, float(np.e)]:
    print(f"    Z_O={ZO:6.4f}:  Z_O^2 * f1^2_strong = {ZO*ZO*f1sq:.4f}  (positive for any real ZO!=0)")

# ----------------------------------------------------------------------
# (C) Finite-L residue is trivially >0 for atom AND cut; only L->inf distinguishes
# ----------------------------------------------------------------------
print("\n(C) FINITE-L TRIVIALITY: every residue >0 at finite L (compact transfer matrix)")
print("    Two models at fixed finite L:")
print("      ATOM: 2nd level at FIXED gap 0.5   -> atom survives L->inf, f1^2_inf=f1^2>0")
print("      CUT : 2nd level at (2pi/L)^2/m1 -> 0 -> atom dissolves into continuum, f1^2_inf=0")
print(f"    {'L':>4} {'atom gap':>10} {'cut gap':>10}  {'f1^2_L>0 (both)':>16}")
for L in [2, 4, 8, 16, 32, 64]:
    gap_atom = 0.5
    gap_cut = (2 * np.pi / L) ** 2 / m1
    print(f"    {L:>4} {gap_atom:>10.4f} {gap_cut:>10.5f}  {'yes (trivial)':>16}")
print("    -> BOTH have f1^2_L>0 for every finite L => CP transport at finite L cannot")
print("       distinguish atom from cut. Only L->inf does: atom gap stays 0.5; cut gap ->0.")

# ----------------------------------------------------------------------
# (D) Flow / local RG renorm is blind to the L->inf isolation
# ----------------------------------------------------------------------
print("\n(D) FLOW/RENORM IS LOCAL => BLIND TO L->inf ISOLATION")
print("    Z_O(t) (flow) and Z_O(beta) (RG) are LOCAL operator renormalizations (pointwise/")
print("    statewise). They carry no information about whether the lightest level stays")
print("    ISOLATED as L->inf, a->0. That survival = isolation I2 = weak-* discontinuity = A3 = H1.")

print("\nVERDICT: CP transport is a finite-(L,a) normalization identity (vacuous for the continuum);")
print("the flow preserves the level not the weight and is downstream; finite-L f1^2>0 is trivial.")
print("Exact failing step: L->inf, a->0 atom SURVIVAL + ISOLATION = A3 = H1. Clay 0%; no QED.")
