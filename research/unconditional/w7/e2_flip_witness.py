"""
E2-FLIP-NOGO fresh witness  [V-here]
Three tiny checks backing the four-horn dichotomy. No lattice job; < 1 s.

W1  Finite-moment fiber (horn alpha): fix a GAPPED reference measure on [1,3],
    match its first 4 moments with a measure that ALSO carries mass at s=0.
    LP feasibility (scipy.optimize.linprog). If w0>0 feasible -> moments up to
    order 3 do not pin inf supp; no finite-moment functional lower-bounds it.

W2  CM rate on a massless-type correlator (horn gamma / acid test):
    C(t)=(1+t)^{-p}  (power-law decay => support floor 0). r(t)=-(1/t)log(C/C0)
    is STRICTLY POSITIVE at every finite t and only -> 0 as t->inf. So the
    finite-t rate is positive on the free massless field: it can never certify
    a gap (it is positive exactly where it must be 0).

W3  Scale/moment direction (horn delta / scale-detector lemma):
    for any measure on [m,inf), the mean E[S] >= m (moment-ratio upper bound);
    a "scale detector" Psi(M) with Psi(0)=0 can be a LOWER bound on m only if
    Psi(M) <= m for the local measure's scale M  <=>  anti-concentration (H1).
"""
import numpy as np
from scipy.optimize import linprog

print("== W1: finite-moment fiber (mass at 0 with same M0..M3 as a gapped ref) ==")
# gapped reference: discrete measure on [1,3]
xs = np.array([1.0, 1.4, 2.0, 2.7, 3.0])
ws = np.array([0.10, 0.30, 0.25, 0.25, 0.10])
ws = ws / ws.sum()
M = np.array([np.sum(ws * xs**k) for k in range(4)])  # M0..M3
print("  reference atoms at", xs)
print("  reference weights", np.round(ws, 4))
print("  reference moments M0..M3 =", np.round(M, 5))
print("  reference inf supp =", xs.min())

# LP: variables [w0 (mass at s=0), w_j on grid on [1,3]]
grid = np.linspace(1.0, 3.0, 41)
n = 1 + len(grid)
# maximize w0  -> minimize -w0
c = np.zeros(n); c[0] = -1.0
# equality constraints: total mass = 1 ; sum s^k w = M_k for k=1,2,3
A_eq = []
b_eq = []
# mass
row = np.zeros(n); row[0] = 1.0; row[1:] = 1.0
A_eq.append(row); b_eq.append(1.0)
for k in range(1, 4):
    row = np.zeros(n)
    row[0] = 0.0**k if k>0 else 1.0   # s=0 contributes 0 for k>=1
    row[1:] = grid**k
    A_eq.append(row); b_eq.append(M[k])
A_eq = np.array(A_eq); b_eq = np.array(b_eq)
bounds = [(0, None)] * n
res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
print("  LP status:", res.message)
if res.success:
    w0 = res.x[0]
    resid = A_eq @ res.x - b_eq
    print(f"  MAX feasible mass at s=0 with M0..M3 matched: w0 = {w0:.4f}")
    print(f"  moment residual (should be ~0): max|err| = {np.max(np.abs(resid)):.2e}")
    print(f"  => a measure with inf supp = 0 shares M0..M3 = {np.round(M,5)}")
    print(f"     as the gapped ref (inf supp 1.0). Finite moments do NOT pin the edge.")

print()
print("== W2: CM rate r(t) on massless-type C(t)=(1+t)^-p  (acid-test failure) ==")
p = 2.0
def C(t): return (1.0 + t)**(-p)
C0 = C(0.0)
print(f"  p={p}: C(t)=(1+t)^-{p}, C0={C0}")
print(f"  {'t':>8} {'C(t)':>12} {'r(t)=-(1/t)log(C/C0)':>24}")
for t in [0.5, 1.0, 3.0, 10.0, 30.0, 100.0]:
    r = -(1.0/t)*np.log(C(t)/C0)
    print(f"  {t:8.1f} {C(t):12.5f} {r:24.6f}")
print("  r(t)>0 at every finite t, -> 0 only as t->inf (log(t)/t -> 0):")
print("  a gap tool built from finite-t rate rates POSITIVE on free massless => vacuous.")

print()
print("== W3: scale/moment direction: mean E[S] >= floor m  (upper, not lower) ==")
m = 0.7
# two measures on [m, inf) with different means, same floor
for name, s, w in [
    ("atom at floor", np.array([m]), np.array([1.0])),
    ("floor + tail", np.array([m, 2.0, 5.0]), np.array([0.6, 0.3, 0.1])),
]:
    mean = np.sum(w*s)
    print(f"  {name:16s} floor={m:.3f}  E[S]={mean:.4f}  (E[S]>=m: {mean>=m-1e-12})")
print("  => any 'detector' Psi(M)<=m must satisfy Psi(M)<=floor; a moment-based")
print("     Psi is an UPPER bound. A valid LOWER detector is anti-concentration = H1.")
