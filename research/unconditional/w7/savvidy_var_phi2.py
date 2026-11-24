"""
Savvidy negative-mode variance:  Var(phi^2) under  mu(phi) ∝ exp( -(V) )
    V(phi) = -1/2 c phi^2 + 1/4 lam phi^4 ,   c = gB > 0 , lam = lam4 > 0 .

Claim (saddle point, y = phi^2, measure ∝ y^{-1/2} exp(1/2 c y - 1/4 lam y^2)):
    y-distribution peaks at  y* = c/lam ,  Gaussian width  sigma_y^2 = 2/lam ,
    so   Var(phi^2) = Var(y) = 2/lam + 2/c^2 + O(c^{-4})   (bounded in c).
Contrast:  Var(phi) ~ c/lam  GROWS  (two wells) --- but the coupling is to phi^2.

Physically lam = k * g^2  =>   softening  g^2 Var(phi^2) = 2/k + 2/B^2 = O(1).
"""
import numpy as np
from scipy.integrate import quad


def var_phi2(c, lam, M=20.0):
    """Numerically exact Var(phi^2) = <phi^4> - <phi^2>^2 by quadrature in y=phi^2."""
    ys = c / lam                       # peak of y-distribution
    sig = np.sqrt(2.0 / lam)           # Gaussian width in y
    hw = M * sig
    lo = -hw
    if ys + lo <= 0.0:
        lo = -ys + 1e-12               # keep y>0
    def base(d):                       # d = y - ys ; removes the e^{c^2/4lam} magnitude
        y = ys + d
        return y ** (-0.5) * np.exp(-0.25 * lam * d * d)
    I0, _ = quad(base, lo, hw, limit=800, epsabs=1e-15, epsrel=1e-12)
    I1, _ = quad(lambda d: (ys + d) * base(d), lo, hw, limit=800, epsabs=1e-15, epsrel=1e-12)
    I2, _ = quad(lambda d: (ys + d) ** 2 * base(d), lo, hw, limit=800, epsabs=1e-15, epsrel=1e-12)
    m2 = I1 / I0
    m4 = I2 / I0
    return m4 - m2 * m2, m2


print("A) Var(phi^2) vs 2/lam + 2/c^2   (residual*(c^2) should -> 2)")
print("  c     lam     Var(phi2)        2/lam           2/lam+2/c^2     (Var-2/lam)*c^2")
for lam in [0.5, 1.0, 2.0]:
    for c in [15.0, 20.0, 25.0, 30.0, 40.0, 60.0]:
        if c / lam < 1.2 * np.sqrt(2.0 / lam):   # need peak y* >> width sigma
            continue
        V, m2 = var_phi2(c, lam)
        print(f"{c:5.1f} {lam:5.2f}  {V:.8e}  {2/lam:.8e}  {2/lam+2/c**2:.8e}  {(V-2/lam)*c*c:.5f}")

print()
print("B) boundedness in c at FIXED lam=1:  Var -> 2/lam=2 (does not grow)")
for c in [10.0, 20.0, 40.0, 80.0, 160.0]:
    V, m2 = var_phi2(c, 1.0)
    print(f"  c={c:6.1f}  Var(phi^2)={V:.8f}   <phi^2>={m2:.4e}  (y*=c/lam={c:.1f})")

print()
print("C) physical scaling lam = k g^2, c = gB :  softening g^2 Var(phi^2)")
print("   (problem depends on a=c/sqrt(lam)=B/sqrt(k) only => g^2 Var is g-INDEPENDENT = O(1))")
k = 1.0
for B in [1.0, 10.0]:
    print(f"  -- k={k}, B={B}  (a=B/sqrt(k)={B/np.sqrt(k):.3f}, large-a asymptote 2/k+2/B^2={2/k+2/B**2:.4f})")
    print("     g      lam=kg^2    Var(phi2)        g^2*Var(phi2)")
    for g in [1.0, 0.7, 0.5, 0.4, 0.3]:
        lam = k * g * g
        c = g * B
        V, m2 = var_phi2(c, lam)
        print(f"   {g:5.2f}  {lam:.6e}  {V:.6e}   {g*g*V:.6f}")

print()
print("D) Var(phi) (NOT phi^2) grows:  Var(phi) ~ c/lam  -> diverges")
for c in [10.0, 20.0, 40.0]:
    # Var(phi) under the full two-well measure: dominated by +-sqrt(y*)
    ys = c/1.0
    print(f"  c={c:5.1f}:  two wells +-sqrt(y*) = +-{np.sqrt(ys):.3f}  => Var(phi) ~ c/lam = {c:.1f} (grows)")
