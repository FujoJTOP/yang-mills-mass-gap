"""
star_no_dilaton.py -- STAR-NO-DILATON witness (illustrative, not a Yang-Mills calculation).

Tests the LOGICAL content of "pure YM has no massless dilaton", NOT a dynamics proof:

  (1) ANOMALY: the one-/two-loop beta coefficients give beta(g) < 0 for ALL g>0 in SU(N) pure YM
      (no Banks-Zaks zero).  A non-zero trace anomaly => scale invariance is EXPLICITLY broken =>
      no exact spontaneously-broken scale symmetry => no massless Goldstone dilaton (the would-be
      dilaton is a pseudo-Goldstone, mass^2 proportional to the anomaly, > 0).  Conformal (beta=0,
      mu=0) limit => dilaton mass 0.

  (2) ACID TEST / LOGIC: "no massless dilaton" == "no delta(s) at s=0 in the 0++ spectral measure".
      Demonstrate three spectral measures and the small-s weight fraction w(eps)=mass[0,eps]/total:
        * CONFORMAL continuum (free massless, rho ~ s^{Delta-2}) : w(eps) -> 0   (gapless, NO dilaton)
        * :phi^2: (gapped-atomless, continuum from 4 m^2)       : w(eps) -> 0   (no dilaton)
        * hypothetical DILATON (delta at s=0)                    : w(eps) -> 1  (massless particle)
      => "no delta at 0" holds for the GAUPLESS conformal theory and for the GAPPED :phi^2: alike;
      so the link removes the dilaton MECHANISM but does NOT force "gapped". Necessary, not sufficient.

Run:  py -3.12 -I research/unconditional/w7/star_no_dilaton.py
"""
import numpy as np

print("=" * 74)
print("(1) ANOMALY: beta(g) < 0 for all g>0 in pure SU(N)  => no BZ zero  => scale inv. explicitly broken")
print("=" * 74)
# convention: beta(g) = -b0 g^3/(16 pi^2) - b1 g^5/(16 pi^2)^2 ,  b0 = 11 N/3  (pure Yang-Mills)
for N in (2, 3):
    b0 = 11.0 * N / 3.0
    # two-loop b1 for SU(N) pure YM:  b1 = 34 N^2 / 3
    b1 = 34.0 * N**2 / 3.0
    g = np.linspace(1e-3, 3.0, 6)
    beta = -b0 * g**3 / (16*np.pi**2) - b1 * g**5 / (16*np.pi**2)**2
    print(f"  SU({N}): b0={b0:.2f} b1={b1:.2f}  beta(g) at g=1e-3,0.5,1,2,3:",
          np.round(beta, 6), " all<0:", bool(np.all(beta < 0)))
print("  => theta^mu_mu = (beta/2g) tr F^2  != 0  (anomaly)  => no exact scale symmetry")
print("  Pseudo-dilaton mass^2 ~ anomaly * (scale mu)^2 ; mu=0 (conformal) => m_d = 0 ; mu>0 (YM) => m_d > 0\n")

print("=" * 74)
print("(2) ACID TEST: 'no massless dilaton' == 'no delta(s) at s=0'; conformal-gapless AND :phi^2: both pass")
print("=" * 74)
s = np.linspace(1e-6, 12.0, 400000)
ds = s[1] - s[0]
m = 1.0  # :phi^2: threshold at 4 m^2 = 4

def weight_frac(rho, eps):
    tot = np.sum(rho) * ds
    lo = np.sum(rho[s < eps]) * ds
    return lo / tot if tot > 0 else 0.0

# conformal free-massless: rho(s) ~ s^{Delta-2} for a dimension-Delta operator (Delta=2.5 here) -> ~0 at s=0
rho_conf = np.where(s > 0, s**0.5, 0.0)          # s^{Delta-2}, Delta=2.5
# :phi^2: gapped-atomless: two-particle phase space from s=4 m^2 ; 0 below threshold, NO delta
rho_phi2 = np.where(s > 4*m*m, np.sqrt(np.clip(1 - 4*m*m/s, 0, None)), 0.0)
# hypothetical DILATON: delta at s=0 (massless particle) -- model as a spike at the first grid cells
rho_dil = np.zeros_like(s); rho_dil[s < 0.02] = 1.0

for name, rho in [("CONFORMAL (free massless, gapless)", rho_conf),
                  (":phi^2: (gapped, atomless)", rho_phi2),
                  ("DILATON (massless particle, forbidden)", rho_dil)]:
    w = weight_frac(rho, eps=0.05)
    has_dilaton = w > 0.5
    print(f"  {name:42s}  w(eps=0.05)={w:7.4f}  delta@0? {has_dilaton}")

print("\n  LOGIC:")
print("   * CONFORMAL: gapless, NO dilaton  -> 'no dilaton' is TRUE but it is still GAPLESS.")
print("   * :phi^2:  gapped,  NO dilaton    -> 'no dilaton' TRUE, gap comes from the continuum threshold.")
print("   * DILATON : would have delta@0     -> the ONLY 'no-dilaton' target; the anomaly forbids it in YM.")
print("  => the no-dilaton link is NECESSARY (kills the delta@0 branch) but NOT SUFFICIENT for the gap")
print("     (the unbroken-CFT continuum-from-0 branch remains; that is the OTHER link = no-Banks-Zaks).")
