"""
dse_fixpoint_witness.py  --  SLOT-DSE-FIXPOINT numeric witness  [V-here]

Three schematic witnesses for research/unconditional/SLOT-DSE-FIXPOINT.md.
NOT a Yang-Mills calculation. A reduced 0-dimensional gap equation that
reproduces, faithfully and transparently, the *logical* structure of the
PT-BFM gluon DSE:

  W1  m^2 = lambda * K(m^2) * m^2 ,  x := m^2/scale^2
        (a) REGULAR vertex   K(x) = 1/(1+x)      -> seagull identity regime
        (b) POLE   vertex    K(x) = 1/x + 1/(1+x) -> Schwinger-mechanism regime
      In both:  x = 0 is always a solution.
      Zero branch (x*>0) exists for lambda above a *scheme/vertex-dependent*
      critical coupling; x* -> 0 at the threshold (NON-uniform).

  W2  fixed-point iteration: x=0 is the contraction attractor for lambda<=lambda_c
      (the "provable branch collapses to the gapless solution"), pitchfork above.

  W3  Oehme-Zimmermann: a propagator can have D(0)=1/m_app^2 > 0 AND satisfy
      the superconvergence sum rule  integral rho ds = 0  ONLY because rho is
      sign-indefinite.  Then m_app is NOT a positive spectral atom.  Contrast a
      single positive atom (OZ-inconsistent) and the physical 0++ channel
      (rho>=0, sign-definite -> gives UPPER bounds only).

Run:  py -3.12 research/unconditional/w7/dse_fixpoint_witness.py
"""

import numpy as np

print("=" * 72)
print("W1 - reduced gap equation: regular vertex vs vertex-pole (seagull)")
print("=" * 72)
print("""
Reduced DSE (0-dim, x = m^2 / scale^2):
    x = lambda * K(x) * x          (m^2 = Pi(m^2) structure)
  x = 0 (gapless / free / scaling) is a solution for EVERY K, lambda.
  A nonzero root needs K(x) nonzero as x->0 (divide by x): 1 = lambda K(x).

(a) REGULAR vertex:  K(x) = 1/(1+x),  K(0)=1 finite
      -> nonzero branch x* = lambda - 1, only for lambda > lambda_c = 1.
(b) POLE vertex (Schwinger mechanism): K(x) = 1/x + 1/(1+x)
      -> 1/x term cancels the x prefactor: x = lambda(1 + x/(1+x))
      -> x* = [ (2 lambda -1) + sqrt((2 lambda -1)^2 + 4 lambda) ] / 2 > 0  for all lambda>0.
""")

def xstar_regular(lam):
    # nonzero root of x = lam * x/(1+x)  ->  x = lam - 1
    return lam - 1.0 if lam > 1.0 else 0.0

def xstar_pole(lam):
    # nonzero root of x = lam*(1 + x/(1+x))
    a = 2*lam - 1.0
    return 0.5*(a + np.sqrt(a*a + 4*lam))

print(f"{'lambda':>8} | {'x*(regular)':>12} {'(exists?)':>9} | {'x*(pole)':>12}")
print("-" * 52)
for lam in [0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0]:
    xr = xstar_regular(lam)
    xp = xstar_pole(lam)
    print(f"{lam:8.2f} | {xr:12.6f} {('yes' if xr>0 else 'NO'):>9} | {xp:12.6f}")

print("""
Reading:
  * m^2 = 0 solves BOTH models at every coupling  (gapless branch survives).
  * REGULAR: positive branch appears only for lambda > 1 (=the seagull identity
    regime: below the critical vertex strength the ONLY solution is m=0).
  * POLE:    positive branch exists for every lambda > 0 -- because the 1/q^2
    vertex pole is exactly what the STI/seagull does NOT fix.
  => WHICH branch exists is decided by the vertex, i.e. by the SCHEME, not by
     the equation or the STI.  This is the truncation-dependence, made explicit.
""")

# ---------------------------------------------------------------- W2
print("=" * 72)
print("W2 - x=0 is the contraction attractor below threshold (pitchfork above)")
print("=" * 72)

def iterate(x0, lam, pole=False, n=400):
    x = x0
    for _ in range(n):
        K = (1.0/x + 1.0/(1.0+x)) if pole else (1.0/(1.0+x))
        x = lam * K * x
        if not np.isfinite(x) or abs(x) > 1e12:
            return x
    return x

print("REGULAR vertex, iterate from x0=0.5 (must contract to x=0 if lambda<1):")
for lam in [0.5, 0.9, 1.1, 1.5]:
    xf = iterate(0.5, lam)
    print(f"   lambda={lam:4.2f}: from x0=0.5 -> x_final={xf:.6f}   "
          f"(x*=max(0,lambda-1)={max(0.0,lam-1):.4f})")
print("""
  lambda<1: x=0 is the stable (contraction) fixed point -> the PROVABLE branch
            collapses to the gapless solution (matches MS-GAP-EQUATION).
  lambda>1: x=0 loses stability (pitchfork), x*=lambda-1 attracts.
""")

# ---------------------------------------------------------------- W3
print("=" * 72)
print("W3 - Oehme-Zimmermann: D(0)>0 coexists with integral rho ds = 0")
print("=" * 72)
print("""
Gluon (transverse) spectral function obeys the superconvergence relation
   integral_0^inf rho_A(s) ds = 0      (gluon not a positive-norm asymptotic state)
and  D(k^2) = integral rho_A(s)/(k^2+s) ds .

A single positive atom  rho = delta(s-a):   integral rho ds = 1 != 0  -> OZ-VIOLATING.
   D(0) = 1/a  > 0   ('mass' = sqrt(a), a positive atom).

OZ-CONSISTENT two-atom example  rho = delta(s-a) - delta(s-b),  0<a<b:
   integral rho ds = 1 - 1 = 0        (OZ satisfied)
   D(0) = 1/a - 1/b > 0               (IR-finite, 'mass'-like scale!)
   but the weight at s=b is NEGATIVE  -> the spectrum is sign-INDEFINITE.
   m_app = 1/sqrt(D(0)) is NOT a spectral atom.
""")

def two_atom(a, b, k2=0.0):
    # D(k^2) = 1/(k^2+a) - 1/(k^2+b)
    return 1.0/(k2+a) - 1.0/(k2+b)

for (a, b) in [(1.0, 2.0), (1.0, 10.0), (0.5, 5.0)]:
    D0 = two_atom(a, b, 0.0)
    mapp2 = 1.0/D0 if D0 != 0 else np.inf
    print(f"   a={a:4.1f}, b={b:5.1f}:  int rho ds = 0 ;  D(0)={D0:.4f} ;  "
          f"apparent m^2=1/D(0)={mapp2:.4f} ;  true atoms at s={a} (+1), s={b} (-1)")

print("""
Physical 0++ channel (O = (1/4) tr F^2): in Euclidean signature O >= 0 pointwise
=> <O(x)O(0)> >= 0 => rho_O(s) >= 0 sign-DEFINITE. Hence NO OZ cancellation there;
but a sign-definite rho gives moment/sum-rule UPPER bounds on the atom f_1^2
(IL-i death, KEY-DIRECT-DERIVATION Attack IV).  The two channels therefore give
opposite-direction information:
   gluon (1--, OZ): sign-indefinite rho_A -> 'mass' is unphysical, no transfer;
   0++ (O):         sign-definite  rho_O -> UPPER bounds only.
Neither yields a scheme-independent sign-definite LOWER bound on m_{0++}. [V-here]
""")
print("=" * 72)
print("DONE.  All three witnesses are schematic; no Yang-Mills number is asserted.")
print("=" * 72)
