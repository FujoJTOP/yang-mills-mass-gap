#!/usr/bin/env python3
"""
e2_flip_l3_witness.py  --  [V-here] witness for E2-FLIP-L3-BETA-CRITICAL

Tests three claims about a putative beta-critical ("jumps at the Savvidy crossover")
flip functional Phi that is supposed to lower-bound the transfer gap m1 WITHOUT
controlling the nonperturbative crossover transport.

  W1  moment/aggregate functionals are edge-blind (banked moment lemma): two positive
      measures with identical low moments and different support floors => NO function
      of those moments lower-bounds the floor.  (analyticity of that function is irrelevant.)
  W2  beta-criticality (crossover sharpness) is DECOUPLED from the gap magnitude:
      the sharpness max|dm1/dbeta| can be driven to infinity (a true jump) while the
      physical-side gap inf_{beta>beta_c} m1 can be set to ANY value, including 0.
      => "being beta-critical" carries no gap lower-bound information.
  W3  far-side-knowledge lemma: the largest constant Delta such that "Phi == Delta on the
      physical side" is a VALID lower bound (Delta <= m1 for all beta>beta_c) is exactly
      inf_{beta>beta_c} m1 -- i.e. the far-side gap.  A discontinuous minorant must be
      calibrated to the gap it is meant to bound; the seam/jump gives it for free on neither side.

This is a witness of already-banked lemmas (moment lemma GAP-G2-edge / SCHUR-CONVEX;
Theorem 0 / acid test; single-scale direction) plus the L-3-specific decoupling W2.
It makes NO claim about the real SU(N) spectrum; it tests the LOGIC of L-3 exactly.
Pure numpy, single process, cheap.
"""
import numpy as np

def H(t): print("\n" + "="*72 + "\n" + t + "\n" + "="*72)

# ---------------------------------------------------------------------------
H("W1  moment lemma: aggregate functionals cannot lower-bound the spectral edge")
# Two positive measures on [0,inf) with IDENTICAL M0,M1 but different floors
# (the classical delta_0+delta_2 vs 2 delta_1 witness; order-3 upgrade = banked linear-data lemma).
def moments(atoms, kmax=3):
    M = [0.0] * kmax
    for s, w in atoms:
        for k in range(kmax):
            M[k] += w * s**k
    return np.array(M)

A = [(0.0, 1.0), (2.0, 1.0)]      # delta_0 + delta_2   floor = 0
B = [(1.0, 2.0)]                  # 2 delta_1           floor = 1
mA, mB = moments(A), moments(B)
print(f"  A = delta_0 + delta_2 :  M0,M1 = {mA[:2]}  floor=0")
print(f"  B = 2*delta_1         :  M0,M1 = {mB[:2]}  floor=1")
print(f"  identical (M0,M1) to 1e-12?  {np.allclose(mA[:2], mB[:2], atol=1e-12)}")
print("  => any Phi(M0,M1) takes the SAME value on A and B, yet floors differ (0 vs 1).")
print("     No function of the first moments -- smooth OR non-analytic -- lower-bounds the floor.")
print("     (Matching up to ANY finite order is the banked linear-data lemma; floor stays free.)")
# and: M1/M0 is an UPPER bound on the squared floor for any measure supported on [m,inf)
print("\n  (direction) for any measure with floor m:  M1/M0 = int s dnu / int dnu >= m.")
for m_floor in [0.0, 0.5, 1.0]:
    atoms = [(m_floor, 0.3), (m_floor + 3.0, 0.7)]
    M = moments(atoms)
    print(f"    floor={m_floor:4.1f}  M1/M0={M[1]/M[0]:6.3f} >= floor_expect {m_floor:4.1f}  "
          f"sqrt(M1/M0)={np.sqrt(M[1]/M[0]):6.3f}  [UPPER on the floor]")

# ---------------------------------------------------------------------------
H("W2  beta-criticality (crossover sharpness) is decoupled from the gap magnitude")
# edge model: strong-side value 1 for beta<beta_c, relaxes to g_phys for beta>beta_c,
# with crossover sharpness set by kappa (kappa->inf = a genuine jump = 'beta-critical').
beta_c = 2.5
beta = np.linspace(0.0, 8.0, 400001)
print(f"  {'kappa':>8} {'g_phys':>9} {'max|dm1/dbeta|':>16} {'inf_{beta>bc} m1':>18}")
for kappa in [1.0, 10.0, 100.0, 1000.0, 10000.0]:
    for g_phys in [1.0, 0.1, 1e-3, 0.0]:
        sig = 1.0 / (1.0 + np.exp(np.clip(kappa * (beta - beta_c), -700, 700)))
        m1 = g_phys + (2.0 - g_phys) * sig          # strong-side value 2, physical-side g_phys
        sharp = np.max(np.abs(np.gradient(m1, beta)))
        far = m1[beta > beta_c].min()
        print(f"  {kappa:8.0f} {g_phys:9.4g} {sharp:16.3f} {far:18.5f}")
print("  => sharpness grows without bound with kappa (approaching a true jump = 'beta-critical'),")
print("     while the physical-side gap inf_{beta>bc} m1 = g_phys is set INDEPENDENTLY (even 0).")
print("     'Being beta-critical' therefore carries NO lower-bound information about m1.")

# ---------------------------------------------------------------------------
H("W3  far-side-knowledge lemma: the valid jump-constant IS inf_{beta>beta_c} m1")
# A discontinuous candidate Phi(beta) = Delta * 1[beta>beta_c] with Phi<=m1 everywhere.
# The largest admissible Delta on the physical side is  inf_{beta>beta_c} m1.
print(f"  {'g_phys':>10} {'Delta_max = inf m1':>20} {'valid iff Delta<=Delta_max':>30}")
for g_phys in [1.0, 0.5, 0.1, 1e-2, 1e-3, 0.0]:
    b = np.linspace(beta_c, 40.0, 600001)
    m1 = g_phys + (1.0 - g_phys) * np.exp(-(b - beta_c))   # relax to g_phys
    Delta_max = m1.min()
    ok = "Delta_max == g_phys" if abs(Delta_max - g_phys) < 1e-9 else f"{Delta_max:.6f}"
    print(f"  {g_phys:10.4g} {Delta_max:20.6f}   {ok:>30}")
print("  => a discontinuous minorant's admissible constant is the far-side gap itself;")
print("     the seam gives no calibration on either side.  Knowing Delta <= inf_{>bc} m1")
print("     REQUIRES knowledge of m1 on the far (physical) side = the transport target.")

# ---------------------------------------------------------------------------
H("Summary")
print("W1 [moment lemma]      : aggregates -> edge-blind / UPPER bound (analyticity-blind).")
print("W2 [decoupling]        : beta-criticality (sharpness) does NOT control the gap.")
print("W3 [far-side lemma]    : a discontinuous minorant must be <= m1 on the far side,")
print("                         i.e. calibrated to the very gap it should avoid computing.")
print("Together: no beta-critical functional can lower-bound m1 without far-side/transport knowledge.")
