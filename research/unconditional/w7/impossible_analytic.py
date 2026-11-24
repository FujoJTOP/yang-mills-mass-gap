"""
impossible_analytic.py -- witness for IMPOSSIBLE-ANALYTIC (cheap, numpy/scipy).
  (1) Dimensional transmutation = ESSENTIAL SINGULARITY at g=0:
      m/mu = exp(-1/(2 b0 g^2)); e^{-1/g^2}/g^{2k} -> 0 for ALL k, so every
      derivative at g=0 vanishes -> the Taylor series at g=0 is IDENTICALLY 0
      -> perturbation theory (analytic in g^2) sees m=0 at every order.
  (2) Strong-coupling (beta->0) is the wrong (infinite-a) limit: m*a=4 ln(1/u)
      grows while a->inf, so physical m = m*a/a -> 0.
  (3) Torus T^3xS^1: m(L) continuous (one phase); extending to L->inf needs a
      UNIFORM lower bound (continuity alone insufficient) = the transport = H1.
"""
import numpy as np
try:
    from scipy.special import ive
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

b0 = 11/(24*np.pi**2)   # SU(2) one-loop

def gap_scale(g2):
    return np.exp(-1.0/(2*b0*g2))

print("="*70)
print("(1) essential singularity: all derivatives of e^{-1/g^2} vanish at g=0")
print("="*70)
print("g        m/mu          e^{-1/g^2}/g^{2k}  for k=1,2,3   (all -> 0)")
for g in [1e-1, 3e-2, 1e-2, 3e-3, 1e-3]:
    g2 = g*g
    row = [gap_scale(g2)/(g2**k) for k in (1,2,3)]
    print(f"{g:.0e}  {gap_scale(g2):.3e}   " + "  ".join(f"{r:.2e}" for r in row))
print("=> every Taylor coefficient at g=0 is 0 -> the weak-coupling expansion")
print("   of the gap is IDENTICALLY ZERO (gap is purely nonperturbative).")
print("   Analytic continuation in g^2 from g^2->0 gives m=0, NOT the g^2~1 value.")

print()
print("="*70)
print("(2) strong-coupling beta->0 is the wrong (infinite-a) limit")
print("="*70)
if HAVE_SCIPY:
    def u(beta):
        return ive(2,beta)/ive(1,beta)
    print("beta   u=I2/I1     m*a=4ln(1/u)  (a->inf as beta->0 => m=m*a/a->0)")
    for beta in [2.5, 1.0, 0.5, 0.2, 0.1]:
        uu = u(beta)
        print(f"{beta:5.2f}   {uu:.3e}   {4*np.log(1/uu):8.3f}")

print()
print("="*70)
print("(3) torus: continuity alone is insufficient; need a uniform lower bound")
print("="*70)
L = np.array([0.05,0.1,0.3,0.7,1.5,3.0,10.0,30.0])
mL = 1.0/(1.0+L) + 0.5   # illustrative continuous positive profile -> 0.5
for l, m in zip(L, mL):
    print(f"L={l:6.2f}  m(L)={m:.4f}")
print("=> m(inf)=lim m(L) requires a UNIFORM lower bound over all L; that")
print("   uniform bound = the transport = H1.")

print()
print("VERDICT: the analytic-continuation base dies at the ESSENTIAL SINGULARITY")
print("(dimensional transmutation). Weak coupling sees 0 (all orders), strong")
print("coupling is the wrong limit, the torus needs a uniform lower bound.")
print("The continuation IS the wall = H1.")
