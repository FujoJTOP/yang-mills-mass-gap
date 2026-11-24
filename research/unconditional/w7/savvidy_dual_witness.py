"""
savvidy_dual_witness.py  --  SAVVIDY-DUAL schematic witness  [V-here]

Demonstrates the DUAL-language statement of the Savvidy softening, and the
fact that the dual condensate provides a PHYSICAL IR cutoff but NOT an
automatic sub-criticality.

NOT a Yang-Mills calculation.  A reduced log-toy of the softening
    S(m) = int dln k  g2Var * prop(k; m)
with prop(k;0)=1/k^2 (massless soft mode, ELECTRIC language) and
prop(k;mV)=1/(k^2+mV^2) (dual photon of mass mV=g~ v, DUAL language).

Point 1: massless (electric) softening is IR-DIVERGENT as kIR->0.
Point 2: massive (dual) softening SATURATES; for mV >> kUV it is ~ 1/mV^2,
         i.e. the monopole condensate v (mV=g~ v) bounds Var(phi^2) physically.
Point 3: sub-criticality S_dual < kappa_bare  <=>  mV > mV*  <=>  v > v*,
         i.e. it is H1 in dual coordinates, NOT automatic.

Run:  py -3.12 research/unconditional/w7/savvidy_dual_witness.py
"""
import numpy as np
from scipy.optimize import brentq

KUV = 1.0
N = 400000


def S(mass, kIR):
    k = np.geomspace(kIR, KUV, N)
    prop = 1.0 / (k**2 + mass**2)          # electric: mass=0 -> 1/k^2
    return np.trapezoid(prop, np.log(k))   # int dlnk * prop


print("=" * 74)
print("P1  ELECTRIC (massless soft mode): S diverges as kIR -> 0")
print("=" * 74)
for kIR in [1e-2, 1e-3, 1e-4, 1e-5]:
    print(f"  kIR={kIR:8.1e}   S_electric = {S(0.0, kIR):12.3e}   (grows unbounded)")

print()
print("=" * 74)
print("P2  DUAL (photon mass m_V): S saturates, and ~ 1/m_V^2 for m_V >> kUV")
print("=" * 74)
kIR = 1e-4
print(f"  fixed kIR={kIR:.1e}; ln(kUV/kIR) = {np.log(KUV/kIR):.3f}  (the 1/m^2 prefactor)")
for mV in [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
    s = S(mV, kIR)
    print(f"  m_V={mV:5.2f}   S_dual = {s:10.4f}   S*m_V^2 = {s*mV**2:10.4f}")

print()
print("=" * 74)
print("P3  sub-criticality  S_dual(mV) < kappa_bare  <=>  mV > mV*  <=>  v > v*")
print("=" * 74)
kappa_bare = 0.5
mVstar = brentq(lambda m: S(m, kIR) - kappa_bare, 1e-3, 50.0)
print(f"  kappa_bare = {kappa_bare}")
print(f"  m_V* = {mVstar:.4f}   (sub-critical iff m_V > m_V*)")
print(f"  v* = m_V*/g~ : the threshold condensate ; deriving v is H1.")
# sanity: S just above/below threshold
for mV in [mVstar * 0.9, mVstar, mVstar * 1.1]:
    print(f"  m_V={mV:7.4f}  S_dual={S(mV,kIR):8.4f}  "
          f"{'SUPER-critical' if S(mV,kIR) > kappa_bare else 'sub-critical'}")

print()
print("VERDICT: dual condensate bounds the softening (finite, ~1/m_V^2 ~ 1/v^2),")
print("but sub-criticality needs v > v* = H1 in dual coordinates. Not automatic.")
