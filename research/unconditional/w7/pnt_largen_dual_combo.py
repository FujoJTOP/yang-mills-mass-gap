"""
pnt_largen_dual_combo.py -- witness for PNT-LARGEN-DUAL-COMBO (2026-11-24)

Illustrative only; NO claim about Yang-Mills. Checks four things:
 (1) Borel pole structure: dim-d operator -> leading IR renormalon at u=d/2
     (dim-2 monopole condensate u=1 ; dim-4 gluon condensate u=2).
 (2) dual GL ratio  m_H/sqrt(sigma) = kappa_GL / sqrt(pi ln kappa_GL)  (g~=1).
 (3) large-N 1/N^2 extrapolation from quoted data.
 (4) consistency band and the honesty check: the combined object would
     PREDICT kappa_GL; without the dual master field kappa_GL is free.

Run:  py -3.12 -I research/unconditional/w7/pnt_largen_dual_combo.py
"""
import math

print("=" * 72)
print("(1) Borel structure: dim-d operator => leading IR renormalon at u=d/2")
print("    schematic a_n ~ Gamma(n+1+d/2) (factorial growth => Borel pole at u=d/2)")
print("=" * 72)
def a_n(n, d):
    return math.gamma(n + 1 + d / 2.0)
for d, label in ((2.0, "dim-2 (monopole condensate v^2): pole u=1"),
                 (4.0, "dim-4 (gluon condensate trF^2): pole u=2")):
    print("\n " + label)
    for n in (5, 10, 20, 30):
        r = a_n(n, d) / a_n(n - 1, d)
        print(f"   n={n:2d}: a_n={a_n(n,d):10.3e}  a_n/a_(n-1)={r:6.2f}")
print("\n => dim-2 has a less severe (u=1) Borel pole than dim-4 (u=2).")

print("\n" + "=" * 72)
print("(2) dual GL ratio  m_H/sqrt(sigma) = kappa_GL / sqrt(pi ln kappa_GL)  (g~=1)")
print("=" * 72)
def dual_ratio(k):
    return k / math.sqrt(math.pi * math.log(k))
for k in (2, 5, 10, 15):
    print(f"   kappa_GL={k:2d}: m_H/sqrt(sigma) = {dual_ratio(k):.3f}")

print("\n" + "=" * 72)
print("(3) large-N extrapolation  (quoted: m(SU2)=3.79, m(SU8)=3.33)")
print("=" * 72)
m2, m8 = 3.79, 3.33
c = (m2 - m8) / (0.25 - 1.0 / 64.0)
m_inf = m2 - c * 0.25
print(f"   1/N^2 fit  m(N)=m_inf+c/N^2 :  m_inf={m_inf:.3f}, c={c:.3f}")

print("\n" + "=" * 72)
print("(4) consistency band + honesty check")
print("=" * 72)
print(f"   large-N          m_0++/sqrt(sigma) -> {m_inf:.2f}")
print(f"   lattice          m_0++/sqrt(sigma) ~  3.90")
lo, hi = 2.0, 200.0
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if dual_ratio(mid) < 3.9:
        lo = mid
    else:
        hi = mid
print(f"   dual GL reaches lattice 3.9 at kappa_GL ~ {0.5*(lo+hi):.2f}  (a FIT, not a prediction)")
print("\n   => a combined 'dual master field' would PREDICT kappa_GL (hence m/√σ).")
print("      Without it, kappa_GL is free: the band 3.3-3.9 is consistency, not derivation.")
print("\nNO dual master field constructed. Bridge blocked at large-N Abelian projection. Clay 0%; no QED.")
