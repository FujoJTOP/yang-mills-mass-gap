"""
savvidy_negmode_exact.py -- witness for SAVVIDY-NEGMODE-EXACT.

Exact (nonperturbative) treatment of the Savvidy negative-mode direction.
Toy but faithful model of the unstable mode: a single mode phi with
    V(phi) = - (c/2) phi^2 + (lam/4) phi^4 ,   c = gB > 0 (the instability),
so phi=0 is a MAX of V (curvature V''(0) = -c = the negative mode), and the
true minima are phi = +/- v, v = sqrt(c/lam), with V''(+/-v) = 2c > 0 (STABLE).

Weight e^{-V(phi)} = e^{(c/2) phi^2 - (lam/4) phi^4}, exact 0+0d quadrature
(exponent shifted by V_min to avoid overflow; ratios are shift-invariant).
We compute exact <phi^2>, Var(phi^2) and compare the softening g^2 Var(phi^2)
under the EXACT (true-vacuum) treatment vs the Gaussian-around-0 (wrong saddle).

Run: py -3.12 -I research/unconditional/w7/savvidy_negmode_exact.py
"""
import numpy as np

def stats(c, lam, L=None, n=2000001):
    if L is None:
        L = max(10.0, 8.0*np.sqrt(c/lam) + 6.0)
    phi = np.linspace(-L, L, n)
    V = -0.5*c*phi**2 + 0.25*lam*phi**4
    Vmin = V.min()                      # at +/-v, Vmin = -c^2/(4 lam)
    w = np.exp(-(V - Vmin))             # shift: max weight = 1, no overflow
    Z = np.trapezoid(w, phi)
    m1 = np.trapezoid(phi**2 * w, phi)/Z
    m2 = np.trapezoid(phi**4 * w, phi)/Z
    return m1, m2 - m1**2               # <phi^2>, Var(phi^2)

lam = 1.0
print(f"lambda = {lam};  v=sqrt(c/lam);  V''(minima)=2c>0 (STABLE)")
print(f"{'c':>8} {'v':>8} {'<phi^2>':>10} {'Var(phi^2)':>12} {'2/lam':>7} {'2c':>8}")
for c in [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]:
    m1, var = stats(c, lam)
    print(f"{c:8.3f} {np.sqrt(c/lam):8.3f} {m1:10.5f} {var:12.6f} {2/lam:7.3f} {2*c:8.1f}")

print("\nGaussian-around-0 (WRONG saddle): 1/(2|V''(0)|) = 1/(2c)  -> INF as c->0")
for c in [0.01, 0.1, 1.0, 10.0, 100.0]:
    print(f"  c={c:7.2f}: 1/(2c) = {1/(2*c):.5f}")

print("\nExact asymptote: Var(phi^2) -> 2/lam (bounded, c-INDEPENDENT) for large c")
for c in [10.0, 20.0, 50.0, 100.0]:
    _, var = stats(c, lam)
    print(f"  c={c:7.1f}: Var(phi^2)={var:.6f}   (target 2/lam = {2/lam:.1f})")

print("\nSoftening g^2*Var(phi^2) with g^2=1:  EXACT bounded; WRONG-saddle g^2/(2c)->inf")
for label, c in [("small c", 0.01), ("mid c", 5.0), ("large c", 100.0)]:
    _, var = stats(c, lam)
    print(f"  {label:8s} c={c:7.2f}: exact Var(phi^2)={var:.5f}   wrong-saddle 1/(2c)={1/(2*c):.5f}")
