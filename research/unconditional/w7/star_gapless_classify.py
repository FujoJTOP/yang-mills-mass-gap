"""
star_gapless_classify.py -- cheap illustrative witness for STAR-GAPLESS-CLASSIFY.md
No heavy machinery; exact/analytic checks of the classification's content:
  (1) a scale-invariant scalar 2-pt function is forced to a pure power law (no scale);
  (2) scale breaking = a gap, in the simplest (free massive) case: gap = m;
  (3) the a-theorem direction: a_UV (massless/free) > a_IR = 0 (gapped).
Run: py -3.12 -I research/unconditional/w7/star_gapless_classify.py
Nothing here is a proof of (star); it illustrates the classification.
"""
import numpy as np

print("="*70)
print("(1) scale invariance forces a power-law 2-pt function (no scale)")
print("="*70)
# A scale-covariant scalar correlator G(x)=c|x|^{-2Delta} under x->lambda x
# scales as lambda^{-2Delta}: the ONLY functional form with no dimensionful
# parameter. Any other form (e.g. a mass term) introduces a scale.
# Check numerically that G(lambda x) = lambda^{-2D} G(x) for the power law:
Delta, c = 0.6, 1.0
for x in (0.3, 1.0, 2.5):
    for lam in (0.5, 2.0):
        lhs = c*abs(lam*x)**(-2*Delta)
        rhs = lam**(-2*Delta)*c*abs(x)**(-2*Delta)
        assert abs(lhs-rhs) < 1e-12
print("  power-law G=c|x|^{-2Delta}: scale-covariant for all x,lambda  [no scale]")

print("\n" + "="*70)
print("(2) scale breaking => a gap, free massive case: E(p)=sqrt(p^2+m^2), gap=m")
print("="*70)
for m in (0.0, 0.2, 1.0):
    p = np.linspace(0, 5, 2000)
    E = np.sqrt(p**2 + m**2)
    print(f"  m={m:.2f}:  min_p E(p) = {E.min():.6f}   (= m; m=0 => gapless)")

print("\n" + "="*70)
print("(3) a-theorem direction: a_UV > a_IR (=0 if gapped)")
print("="*70)
a_scalar, a_weyl, a_vector = 1/90, 11/180, 62/360   # standard 4D anomaly coefficients
for N in (2, 3):
    a_UV = a_vector*(N*N-1)            # pure SU(N) gluons (UV, AF fixed point)
    print(f"  SU({N}): a_UV = 62/360*(N^2-1) = {a_UV:.4f}  >  a_IR = 0 (gapped)  [a-theorem OK]")
print("  (a gapless IR CFT would need a_IR>0 at an *interacting fixed point* = a strong-coupling")
print("   zero of beta; pure YM has no perturbative zero => that zero's absence is the open (G1).)")

print("\nSummary: classification = {gapped, CFT, dilaton(flat dir)}  (D>=3, unitary, relativistic, endpoint).")
print("(star) for pure YM  <=>  'beta has no strong-coupling zero'  (G1, open).  No QED.")
