"""
star_combine.py -- STAR-COMBINE witness (pure standard-library / numpy, runs clean)

Checks the concrete, falsifiable pieces of the (star) chain that CAN be checked:
  - (0) pure YM is asymptotically free: 2-loop beta(g) = -(b0 g^3 + b1 g^5), b0,b1 > 0.
  - (4a) no Banks-Zaks zero: beta(g) < 0 for ALL g > 0 (no g* > 0 with beta(g*) = 0) for SU(N), N=2,3,4.
The nonperturbative links (2),(4b) are NOT numerically checkable here -- they are open
(the report labels them [O]); this witness only plays the beta-function bookkeeping.

Run:  py -3.12 -I star_combine.py
"""
import numpy as np

# 2-loop coefficients for pure SU(N) Yang-Mills (standard: b0 = 11N/3, b1 = 34N^2/3,
# in units of 1/(16 pi^2) absorbed; SIGNS are what matters: both > 0 for pure YM).
def beta_2loop(g, N):
    b0 = 11.0 * N / 3.0
    b1 = 34.0 * N * N / 3.0
    return -(b0 * g**3 + b1 * g**5)

print("=" * 66)
print("(0)+(4a): pure YM AF and NO Banks-Zaks zero")
print("=" * 66)
print(f"{'N':>3} {'beta at g=0.3':>14} {'g=1':>10} {'g=3':>10} {'g=10':>12} {'zero(g>0)?':>12}")
for N in (2, 3, 4):
    gs = np.array([0.3, 1.0, 3.0, 10.0])
    vals = beta_2loop(gs, N)
    # scan for any sign change (zero) on g in (0, 30]
    g = np.linspace(1e-3, 30.0, 300001)
    b = beta_2loop(g, N)
    zero = bool(np.any(np.diff(np.sign(b)) != 0))
    print(f"{N:>3} {vals[0]:>14.4f} {vals[1]:>10.4f} {vals[2]:>10.4f} {vals[3]:>12.4f} {str(zero):>12}")

print()
print("beta(g) = -(b0 g^3 + b1 g^5), b0=11N/3>0, b1=34N^2/3>0  =>  beta<0 for every g>0")
print("==> no g* in (0,inf) with beta(g*)=0  ==> no Banks-Zaks IR fixed point (perturbative).")
print()
print("Chain labels (see STAR-COMBINE.md):")
for i, (ln, lab) in enumerate([
    ("(0) pure YM AF", "[T]"),
    ("(1) AF => IR not free", "[C]/[O] (needs no nonpert fixed point = (4b))"),
    ("(2) gapless interacting 4D relativistic => CFT or dilaton", "[O]  <-- LOAD-BEARING"),
    ("(3) pure YM has no massless dilaton", "[C]"),
    ("(4a) no perturbative IR fixed point (no BZ)", "[T]"),
    ("(4b) no nonperturbative IR fixed point", "[O]  (= the gap)"),
    ("(5) => gapped", "conditional on (2)+(4b)"),
], 1):
    print(f"  {ln:<62} {lab}")
