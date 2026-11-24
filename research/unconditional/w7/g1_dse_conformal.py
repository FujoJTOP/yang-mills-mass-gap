"""
g1_dse_conformal.py -- witness for G1-DSE-CONFORMAL.md

Two cheap, reproducible checks:
 (A) two-loop beta for pure YM (N_f=0): beta(g) < 0 for ALL g>0, NO Banks-Zaks zero;
     contrast N_f=12 (SU(3)) HAS a Banks-Zaks zero => the conformal window is REAL.
 (B) the Landau-gauge scaling-solution IR exponent: a representative leading-order IR balance
     of the ghost+gluon DSE admits a POSITIVE root kappa ~ 0.6 for N_f=0
     => the truncated DSE DOES admit a conformal scaling solution at N_f=0.

Honest note: (B) is a *calibrated representative* of the truncated IR balance; the full
Fischer-Alkofer truncation gives kappa ~ 0.595 (SU(3), 4D). The point is qualitative: a positive
root EXISTS, so the DSE is NOT inconsistent with a scale-covariant IR at N_f=0.  The decisive
content is the physical interpretation (unphysical gauge-fixed propagators, D(0)=0 = positivity
violation = confinement, NOT a massless particle).

Run:  py -3.12 -I g1_dse_conformal.py
"""
import numpy as np
from scipy.special import beta
from scipy.optimize import brentq

print("=" * 72)
print("(A) two-loop beta of pure YM (N_f=0): any Banks-Zaks zero?")
print("=" * 72)

def beta_two_loop(g, Nc=3, Nf=0):
    b0 = (11.0 * Nc - 2.0 * Nf) / (3.0 * (4 * np.pi) ** 2)
    b1 = (34.0 / 3.0 * Nc ** 2 - (10.0 / 3.0) * Nc * Nf - 2.0 * (Nc ** 2 - 1) / Nc * Nf) / (4 * np.pi) ** 4
    return -(b0 * g ** 3 + b1 * g ** 5)

gs = np.linspace(1e-3, 10.0, 200000)
bg = beta_two_loop(gs, 3, 0)
print(f"  N_f=0: max beta on (0,10] = {bg.max():.3e}  (negative everywhere => no BZ zero)")
print(f"  N_f=0: sign changes g>0 = {np.sum(np.diff(np.sign(bg)) != 0)}  (0 => no Banks-Zaks fixed point)")
bg12 = beta_two_loop(gs, 3, 12)
print(f"  [contrast] N_f=12 (SU3): sign changes g>0 = {np.sum(np.diff(np.sign(bg12)) != 0)} "
      f"(>=1 => Banks-Zaks zero => conformal window is REAL)")

print()
print("=" * 72)
print("(B) scaling-solution IR exponent: does the DSE admit a POSITIVE root kappa?")
print("=" * 72)

# Representative calibrated IR balance of the ghost+gluon DSE (angular integration -> Beta kernel):
#   1 = C (1+kappa) B(1+kappa, 1+kappa),  C calibrated so the SU(3),4D root is ~0.595.
#   (Exact coefficient and root are truncation-dependent; Fischer-Alkofer 2003: kappa ~ 0.595.)
C_CAL = 1.90
def ir_balance(kappa):
    return 1.0 - C_CAL * (1.0 + kappa) * beta(1.0 + kappa, 1.0 + kappa)

kappa_root = brentq(ir_balance, 1e-4, 5.0)
print(f"  SU(3), N_f=0: POSITIVE root kappa = {kappa_root:.4f}  (literature ~0.595)")
print(f"    => the truncated DSE ADMITS a conformal scaling solution at N_f=0 (NOT inconsistent)")
print("  f(kappa) sign check:")
for k in (0.0, 0.3, round(kappa_root, 2), 0.9, 1.5):
    print(f"    kappa={k:.2f}: f={ir_balance(k):+.4f}")

print()
print("=" * 72)
print("Physical interpretation (the decisive point):")
print("=" * 72)
print("""  * The scaling solution is about LANDAU-GAUGE (unphysical) gluon/ghost propagators.
  * D_A(0)=0 (gluon suppressed) and G(k)~k^{-kappa} (ghost enhanced) are POSITIVITY-VIOLATING
    / confinement signals, NOT a massless physical particle.
  * The physical 0++ glueball mass is a gauge-invariant correlator the propagator DSE does NOT compute.
  * LATTICE finds the DECOUPLING solution D_A(0)>0 (a 'massive' gluon), not scaling
    (Aguilar-Binosi-Papavassiliou vs Alkofer et al.) -> scaling-vs-decoupling is truncation/scheme dependent.
  => The DSE ADMITS a conformal scaling solution at N_f=0 [C], but it is NOT a physical pure-glue CFT
     and does NOT refute (G1). (G1) stays [O].""")
