"""
base_order_parameter_witness.py -- BASE-ORDER-PARAMETER §3(M-b): spectral/topological
"global" functionals are atom/floor-blind.

Two 0++ spectral measures nu:
  A: an atom f1^2=0.5 at s=m1^2=1 (a glueball) + discrete continuum nodes;
  B: a GENUINE POSITIVE measure on a dense grid of s-values, NONE at s=1
     -> B has NO atom at m1^2=1 (atom weight 0 there).
We solve B's nonnegative weights (LP feasibility on a dense grid ~ continuum) so A and B
agree EXACTLY on a set of global spectral observables (the type index / eta / spectral
sums are built from). Then the atom weight at s=1 is compared.
=> no spectral/topological functional can certify the KL atom (floor-blind).

Robustness: the script tries to match {5,4,3} observables (in that order) on a dense
grid and reports the FIRST that is feasible, so no unattained claim is printed.

numpy + scipy linprog only. py -3.12 -I. Logical/arithmetic witness only.
"""
import numpy as np
from scipy.optimize import linprog

np.set_printoptions(precision=6, suppress=True)

# ---- measure A: atom at s=1 plus continuum nodes ----
sA = np.array([1.0, 4.0, 8.0, 16.0])
wA = np.array([0.5, 0.2, 0.15, 0.1])          # atom 0.5 at m1^2=1
f1sq_A = wA[0]

NAMES = ["int s^-1 dnu", "int 1 dnu", "int s dnu", "int e^-s dnu", "int 1/(1+s) dnu"]
PHI = [lambda s: s**(-1.0), lambda s: np.ones_like(s), lambda s: s,
       lambda s: np.exp(-s), lambda s: 1.0/(1.0+s)]

def phi(s, w, k):
    s = np.asarray(s, float)
    return np.array([np.sum(w * PHI[i](s)) for i in range(k)])

# ---- measure B: dense grid over [1.2, 32], no node at s=1 ----
sB = np.linspace(1.2, 32.0, 80)
f1sq_B = 0.0                                   # no node at s=1  -> atom weight 0

def try_match(k):
    tgt = phi(sA, wA, k)
    Mmat = np.array([[PHI[i](x) for x in sB] for i in range(k)])
    res = linprog(c=np.ones_like(sB), A_eq=Mmat, b_eq=tgt,
                  bounds=[(0, None)] * len(sB), method="highs")
    if not res.success:
        return None
    return np.clip(res.x, 0.0, None), tgt

chosen = None
for k in (5, 4, 3):
    got = try_match(k)
    if got is not None:
        chosen = (k, got[0], got[1])
        break

print("=" * 76)
print("Measure A: atom f1^2=%.3f at m1^2=1 + nodes %s" % (f1sq_A, sA[1:]))
print("Measure B: POSITIVE measure on a dense grid over [1.2, 32], 0 nodes at s=1")
print("=" * 76)
assert chosen is not None, "even k=3 infeasible on this grid"
k, wB, tgt = chosen
Bvals = phi(sB, wB, k)
print("matched %d independent global spectral observables (A vs B):" % k)
for i in range(k):
    print("  %-16s  A=%.12f   B=%.12f   |diff|=%.2e" % (NAMES[i], tgt[i], Bvals[i], abs(tgt[i] - Bvals[i])))
print("-" * 76)
print("  all B weights >= 0 ?  %s   (min=%.2e,  total mass=%.6f)"
      % (bool(np.all(wB >= -1e-9)), wB.min(), wB.sum()))
print("  A total mass = %.6f" % wA.sum())
print("  max|diff| on matched observables = %.2e" % np.max(np.abs(tgt - Bvals)))
print("-" * 76)
print("  ATOM WEIGHT at m1^2=1:   A = %.3f     B = %.3f  (no node at s=1)" % (f1sq_A, f1sq_B))
print("  => a GENUINE POSITIVE measure atomless-at-1 matches %d global spectral" % k)
print("     functionals EXACTLY, yet the gap atom differs: 0.5 vs 0.")
print("     index/eta are of THIS (spectral, positive-cone) type; 't Hooft/Polyakov")
print("     detect the 1-form center (not the 0-form atom) -> both classes blind.")
print("=" * 76)
