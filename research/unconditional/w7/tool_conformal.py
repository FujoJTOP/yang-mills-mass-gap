"""
tool_conformal.py -- TOOL-CONFORMAL witness  [V-num / NG demo]
Run:  py -3.12 -I research/unconditional/w7/tool_conformal.py

Question: can conformal mapping / extremal length / harmonic measure FORCE a lower
bound on the Kallen-Lehmann atom weight f1^2 > 0?

Measure (solvable):  nu = f1 * delta_{m1} + (1-f1) * delta_{m2},   0 < f1 <= 1.
Resolvent:           m(z) = f1/(z-m1) + (1-f1)/(z-m2)   (Stieltjes/Cauchy transform).

Claims tested:
 (A) SUPPORT geometry (log-capacity, extremal length, harmonic measure) is BLIND to
     f1: it is a function of the SUPPORT only, and a point mass has ZERO capacity /
     ZERO harmonic measure. Two measures with identical support {m1,m2} but different
     atom weight have IDENTICAL support geometry  =>  acid test FAILS.
 (B) The ONLY geometric reader is the conformal 'size' of a level set of the RESOLVENT,
     which equals f1/R near the atom  =  the Stieltjes residue = A3 (read-not-force).
"""
import numpy as np

m1, m2 = 1.0, 4.0

def m_res(z, f1):
    return f1 / (z - m1) + (1.0 - f1) / (z - m2)

def cap_two_points(a, b):
    # logarithmic capacity of two points {a,b}
    return abs(a - b) / 4.0

def preimage_radius(f1, R, tmin=1e-7, tmax=2.0, N=800_000):
    # radius t where |m(z)| = R along z = m1 + i t  (preimage of the level circle)
    ts = np.linspace(tmin, tmax, N)
    vals = np.abs(m_res(m1 + 1j * ts, f1)) - R
    s = np.sign(vals)
    idx = np.where(np.diff(s) != 0)[0]
    return ts[idx[0]] if len(idx) else float("nan")

print("=" * 76)
print("(A) SUPPORT geometry is BLIND to the atom weight f1  (support {m1,m2} fixed).")
print("=" * 76)
print(f"    support = {{{m1}, {m2}}};   log-capacity = |m2-m1|/4 = {cap_two_points(m1,m2):.4f}"
      f"  (CONSTANT in f1)")
print( "    harmonic measure of a single POINT = 0   (a point has zero hm)")
print( "    log-capacity of a single point     = 0   (a lone atom is invisible)")
print()
print("       f1     residue@m1 (=-Res Pi)     cap(support)   <- geometry constant")
print("    " + "-" * 62)
for f1 in [0.1, 0.3, 0.5, 0.9, 1.0]:
    print(f"    {f1:5.2f}        {f1:8.4f}               {cap_two_points(m1, m2):6.4f}")

print()
print("=" * 76)
print("(B) The level-set conformal radius of the RESOLVENT reads f1 -- but = residue.")
print("=" * 76)
print("    near the atom m(z) ~ f1/(z-m1); preimage of |m|=R is a contour of radius")
print("    ~ f1/R  around m1  (its 'conformal size' encodes the residue).")
print()
R = 10.0
print("       f1     preimage radius (|m|=R)      f1/R        ratio")
print("    " + "-" * 62)
for f1 in [0.1, 0.3, 0.5, 0.9, 1.0]:
    r = preimage_radius(f1, R)
    print(f"    {f1:5.2f}        {r:.6f}                {f1/R:.6f}    {r/(f1/R):.4f}")

print()
print("=" * 76)
print("ACID TEST: two measures, SAME support {1,4}, DIFFERENT atom weight")
print("=" * 76)
for name, fa in [("nu_a (f1=0.9)", 0.9), ("nu_b (f1=0.1)", 0.1)]:
    print(f"    {name}:  cap(support)={cap_two_points(m1,m2):.4f}    residue@m1={fa:.4f}")
print("    => support geometry IDENTICAL (0.75) though atom weight differs 9x: BLIND.")
print("       Only the RESOLVENT residue distinguishes them = A3 (read-not-force).")

print()
print("VERDICT: conformal/extremal-length/harmonic-measure on the SUPPORT is blind to")
print("         the atom weight (a point has zero capacity); the only geometric reader")
print("         is the resolvent level-set radius = the Stieltjes residue = A3.")
