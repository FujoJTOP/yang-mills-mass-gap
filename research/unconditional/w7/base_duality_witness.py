"""
base_duality_witness.py  --  BASE-DUALITY witness  (2026-11-24)
One small job, numpy/scipy only.  NOT a Yang-Mills calculation.

Question: does ANY duality coordinate (character/Fourier-on-the-group expansion,
electric-magnetic 't Hooft algebra, deformed small-circle semiclassics) turn the
crossover beta~2.5 into a COMPUTABLE statement with a SMALL PARAMETER?

We measure the one thing that decides this for the character (strong-coupling /
random-surface) dual coordinate:  the character coefficients a_r(beta) and their
ratios to the trivial coefficient a_0.  A controlled dual expansion needs
a_{r!=0}/a_0 << 1.  We evaluate at the crossover.

Run: py -3.12 -I base_duality_witness.py
"""
import numpy as np
from scipy.special import iv

print("="*76)
print("BASE-DUALITY witness  --  is there a dual small parameter at beta~2.5 ?")
print("="*76)

print("""
[t Hooft 1978, rigorous]  The electric (Wilson) and magnetic ('t Hooft) loop
operators satisfy the CENTER ALGEBRA
     W(C) B(C') = B(C') W(C) * exp(2 pi i /N * link(C,C')) .
This is the rigorous content of electric-magnetic duality for pure SU(N):
it fixes the non-local ORDER PARAMETER TYPE (loops/surfaces, not local
observables) and the relevant global structure (the center Z_N).  [T] -- but it
is an ALGEBRA, not a small parameter.

[character / dual coordinate]  e^{beta Re tr U} = sum_r d_r a_r(beta) chi_r(U).
The dual (random-surface) expansion weights a surface by prod a_{r}(beta).
Controlled  <=>  a_{r!=0}/a_0 << 1.
""")

# SU(2)/abelian proxy for the character coefficients: modified Bessel ratios
# I_n(2 beta)/I_0(2 beta).  Robust, convention-independent fact used below:
# for ANY compact group, a_r(beta) -> a_r(infinity) = chi_r(1)/d_r-scaled value,
# hence the RATIOS a_r/a_0 -> 1 as beta -> infinity.  So the dual expansion is a
# small-parameter expansion at SMALL beta and DEGENERATE at LARGE beta.
betas = [0.02, 0.125, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0]
print(f"{'beta':>6} {'a1/a0':>8} {'a2/a0':>8} {'a3/a0':>8}   regime")
print("-"*52)
for b in betas:
    z = 2.0*b
    a0 = iv(0, z)
    r1 = iv(1, z)/a0
    r2 = iv(2, z)/a0
    r3 = iv(3, z)/a0
    if b < 1/48:
        regime = "strong-coupling (gap PROVEN, OS/SZZ)"
    elif b < 1.5:
        regime = "dual expansion marginal"
    else:
        regime = "crossover / weak: NO small parameter"
    print(f"{b:6.3f} {r1:8.4f} {r2:8.4f} {r3:8.4f}   {regime}")

print("""
READ:
  * beta < 1/48  (the PROVEN strong-coupling corner): a_r/a_0 is small; the
    dual surface expansion converges.  This IS a duality with a small parameter
    -- but it is the regime where the gap is already a theorem (OS'78, SZZ'23),
    so the duality is a BACKWARD certificate there, not a key.
  * beta = 2.5 (the crossover): a1/a0 ~ 0.89, a2/a0 ~ 0.68.  The dual series is
    NOT controlled -- the surface weights are O(1).  The character expansion is
    degenerate EXACTLY where the problem lives.
""")

# vortex thickness heuristic: core size in lattice units ~ 1/sqrt(beta); the
# 'dilute center-vortex' picture needs core << spacing.
print("="*76)
print("center-vortex dilution diagnostic  (order parameter = 't Hooft / vortex)")
print("="*76)
print(f"{'beta':>6} {'core~1/sqrt(beta)':>18}   dilute (core<<1)?")
for b in [0.5, 1.0, 2.0, 2.5, 3.0]:
    core = 1.0/np.sqrt(b)
    print(f"{b:6.2f} {core:18.3f}   {'yes' if core < 0.5 else 'NO -- vortices overlap, dilute-gas uncontrolled'}")

print("""
So at beta=2.5 the vortices are fat (core ~ 0.63 lattice spacings): no dilute
gas, no controlled semiclassical/vortex dual picture.  The non-local order
parameter is the RIGHT TYPE (T: 't Hooft algebra), but its area/perimeter law
at the crossover is verbatim the confinement statement = H1 -- the difficulty
does NOT transfer to a controlled dual side; it stays put.

Deformed / small-circle escape [Unsal-Yaffe]:  on R^3 x S^1 with a center-
symmetric deformation the theory becomes weakly coupled (a genuine small
parameter, monopole-instantons/bions, calculable gap).  But that is a DEFORMED
theory at small circle; turning off the deformation and growing the circle
(decompactification) is the transport/decompactification wall = the same H1.
""")
print("DONE (no QED; Clay 0%).")
