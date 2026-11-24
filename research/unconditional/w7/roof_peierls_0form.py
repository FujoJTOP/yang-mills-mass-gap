"""
roof_peierls_0form.py -- witness for research/unconditional/ROOF-PEIERLS-0FORM.md

Question: does a 0-form (0++) Peierls/contour argument ESCAPE the circularity that
kills the 1-form (Wilson-loop) Peierls at beta=2.5?

Cheap: single-plaquette character coefficients only (numpy/scipy). No lattice MC.
Run: py -3.12 -I roof_peierls_0form.py
"""
import numpy as np
from scipy.special import iv
from scipy.optimize import brentq

def a_ratio(r, b):
    # SU(2) single-plaquette character coefficient ratio a_r/a_0,
    # convention: plaquette weight exp(b cos theta),  a_r = 2(2r+1) I_{2r+1}(b)/b
    # => a_r/a_0 = (2r+1) I_{2r+1}(b) / I_1(b).
    return (2*r+1)*iv(2*r+1, b)/iv(1, b)

b = 2.5
u1 = a_ratio(1, b)
sigma = -np.log(u1)

print("="*74)
print("SU(2) character coefficients at beta = 2.5   (weight exp(b cos theta))")
print("="*74)
for r in range(0, 5):
    print(f"  a_{r}/a_0 = {a_ratio(r, b):.4f}")
print(f"\n  defect activity  u = a_1/a_0 = {u1:.4f}")
print(f"  surface tension proxy  sigma = -log u = {sigma:.4f}")

# ---- Peierls thresholds: convergence needs sigma > c  (c = log of branching) ----
c1 = np.log(5.0)                          # 1-form 2-surface: each link meets 6 plaquettes -> 5 continuations
c0_lo, c0_hi = np.log(3.0), np.log(7.0)   # 0-form 3-surface: rigorous bounds O(1) in this range
print("\n" + "="*74)
print("Peierls threshold c   (need sigma > c for the contour/polymer sum to converge)")
print("="*74)
print(f"  1-form (2-surface / magnetic string / area law):  c1 = log5      = {c1:.4f}")
print(f"  0-form (3-surface / domain wall / 0++ clustering): c0 in [log3,log7] = [{c0_lo:.4f}, {c0_hi:.4f}]")

print("\n" + "="*74)
print(f"At beta = 2.5:  sigma = {sigma:.4f}")
print("="*74)
print(f"  1-form: sigma - c1 = {sigma - c1:+.4f}   -> {'converges' if sigma > c1 else 'FAILS'}")
print(f"  0-form: sigma = {sigma:.4f} vs c0 >= {c0_lo:.4f}  -> {'converges' if sigma > c0_hi else 'FAILS'}")

# ---- validity window: sigma(beta) > c ----
f = lambda x: -np.log(a_ratio(1, x)) - c1
try:
    bstar = brentq(f, 0.01, 2.4999)
except ValueError:
    bstar = float('nan')

print("\n" + "="*74)
print("Validity window of each Peierls argument  (proves gap for sigma(beta) > c)")
print("="*74)
print("   beta    u=a1/a0    sigma    1form(>log5)   0form(>log7)")
for x in np.linspace(0.02, 2.5, 13):
    u = a_ratio(1, x); s = -np.log(u)
    print(f"  {x:5.3f}   {u:7.4f}   {s:7.4f}     {'ok ' if s > c1 else 'fail':4s}        {'ok ' if s > c0_hi else 'fail':4s}")

print(f"\n  1-form Peierls validity ends at beta* = {bstar:.3f}  (rigorous gap only for beta < beta*)")
print(f"  crossover target beta = 2.5 :  {'INSIDE Peierls regime' if bstar > 2.5 else 'OUTSIDE -> Peierls cannot reach it'}")

print("\n" + "="*74)
print("READOFF")
print("="*74)
print(" * sigma is decreasing in beta (activity grows); Peierls is a small-activity (small-beta) method.")
print(" * 1-form: fails for beta > %.3f; at beta=2.5 grossly outside." % bstar)
print(" * 0-form: uses a DIFFERENT tension (domain-wall, not string), so NOT the same circularity,")
print("   but faces the SAME O(1) activity at beta=2.5 -> the same 'no small parameter' wall.")
print(" * Neither channel's Peierls sum converges at beta=2.5.  [NG] at the crossover; [O] = the exact step.")
