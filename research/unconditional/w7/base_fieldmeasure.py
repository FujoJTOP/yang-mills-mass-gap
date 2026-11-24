"""
base_fieldmeasure.py -- witness for research/unconditional/BASE-FIELDMEASURE.md

Demonstrates the classification: the measure-inequality family (Poincare/LSI gap,
Bakry-Emery, Wasserstein drift, Brascamp-Lieb defect, FKG defect, entropy production,
transfer entropy, Gaussian-domination) is NOT a kappa-independent SOURCE. Each is
either (a) exactly zero on equilibrium, (b) proportional to kappa (homogeneous), or
(c) sign-indefinite. Hence none can be the additive inhomogeneity S in kappa=F(kappa)+S.
The only inhomogeneous field-measure objects are ANOMALIES (scale breaking).

Run:  py -3.12 -I research/unconditional/w7/base_fieldmeasure.py
"""
import numpy as np

print("="*78)
print("PART A | Gaussian coarse source S=1/Var(a) is LINEAR in kappa (=m^2): homogeneous")
print("="*78)
# 1D periodic free field, action 1/2 sum (d phi)^2 + 1/2 m^2 phi^2. The k=0 coarse
# mode a = (1/L) sum_x phi_x = (1/sqrt(L)) phi_tilde_0, so Var(a) = <phi_tilde_0^2>/L = 1/(L m^2).
# Hence S := 1/Var(a) = L * m^2  -- exactly PROPORTIONAL to kappa = m^2.
L = 64
print(f"  L = {L}; S/kappa should be the constant L, NOT a kappa-independent source.")
print(f"  {'m^2=kappa':>12} {'Var(a)':>14} {'S=1/Var':>14} {'S/kappa':>10}")
for m2 in [0.5, 0.1, 0.01, 1e-3, 1e-4]:
    Var = 1.0/(L*m2)
    S = 1.0/Var
    print(f"  {m2:12.4g} {Var:14.4g} {S:14.4g} {S/m2:10.4g}")
print("  => S -> 0 as kappa -> 0, and S/kappa = L const: S is a HOMOGENEOUS function of kappa.")
print("     A genuine SOURCE must be kappa-INDEPENDENT and nonzero at kappa=0; S is neither.")

print()
print("="*78)
print("PART B | Brascamp-Lieb defect D is SIGN-INDEFINITE (not a sign-definite source)")
print("="*78)
# 2-mode toy (as SLOT-SD-FIXPOINT): V(a,y)=1/2 A a^2 + 1/2 B y^2 + (k/4) y^4 + h a y^2.
# Condition on a=0.  D = Var(f|0) - E[c^2/K|0],  f=h y^2, c=2 h y, K=B+3k y^2.
A, h, k = 1.0, 0.5, 1.0
def quad(f, lo=-14.0, hi=14.0, n=400001):
    y = np.linspace(lo, hi, n); w = np.exp(-f(y)); w /= np.trapezoid(w, y); return y, w
def expect(g, y, w):
    return np.trapezoid(g*w, y)          # E[g] with w a normalized DENSITY (per unit y)
print(f"  A={A}, h={h}, k={k}; vary fiber curvature B.  D = Var(f|0) - E[c^2/K|0], f=h y^2.")
print(f"  {'B':>8} {'E[c^2/K|0]':>14} {'h^2 Var(y^2)':>15} {'D':>12}")
for B in [-3.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0]:
    y, w = quad(lambda y: 0.5*B*y*y + (k/4.0)*y**4)
    c2K   = expect((2*h*y)**2/(B + 3*k*y*y), y, w)
    Ey2   = expect(y*y, y, w); Ey4 = expect(y**4, y, w)
    Vary2 = Ey4 - Ey2*Ey2
    D = h*h*Vary2 - c2K
    print(f"  {B:8.2f} {c2K:14.4g} {h*h*Vary2:15.4g} {D:12.4g}")
print("  => D changes sign with B: NOT sign-definite. Cannot serve as a sign-definite source.")

print()
print("="*78)
print("PART C | Equilibrium functionals vanish identically")
print("="*78)
print("  Entropy production sigma = int J . grad(log err) = 0 for ANY equilibrium measure")
print("     (detailed balance => probability current J = 0). So sigma = 0 on every YM measure.")
print("  Transfer entropy between past/future also = 0 at equilibrium (time-reversal symmetry).")
print("  => these cannot be a positive source at beta=2.5 (they are 0 there).")

print()
print("="*78)
print("PART D | The ONLY inhomogeneous field-measure object: the trace anomaly")
print("="*78)
print("  theta^mu_mu = (beta(g)/2g) tr F^2 = -(d/4) tr F^2 (d=4). The vacuum condensate")
print("  <tr F^2> ~ Lambda^4 is INDEPENDENT of the coarse mass kappa (it is a UV-scale, scale-")
print("  breaking object), so it is the one candidate that can be a kappa-INDEPENDENT source.")
print("  BUT: (i) its value is the gluon condensate; proving <tr F^2> != 0 with a magnitude")
print("       floored uniformly in (L,a) IS H1 (the open step).")
print("  (ii) the OTHER anomaly object, topological susceptibility chi_top, sources the 0^-+")
print("       channel (Witten-Veneziano / eta'), NOT the 0^++ channel Clay needs: channel mismatch.")
print()
print("  VERDICT: no construct. But the BASE classification:")
print("    - measure-inequality family: homogeneous / zero / sign-indefinite => NOT a source [T-class]")
print("    - source must be an ANOMALY (scale-breaking, inhomogeneous)")
print("    - the only 0^++ anomaly is the gluon condensate = H1; uniform a->0 = renormalon wall")
print("  Clay 0%; no QED.")
