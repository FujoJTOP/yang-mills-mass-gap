"""
kappa_2d_exact.py -- witness for KAPPA-2D-EXACT.
Is the rigid-string rigidity kappa universal (transferable from 2D / thin models) or a
microscopic LEC?  Three concrete computations, numpy only, <1 s.

(1) The Luscher term is UNIVERSAL: it is the zeta-regularized zero-point energy of the
    D-2 free transverse bosons, so it depends only on c = D-2 and the worldsheet topology.
      closed (periodic):  E0 = (D-2)*(2*pi/L)*zeta(-1)          = -(D-2)*pi/(6L)
      open   (bc):        E0 = (D-2)*(pi/(2L))*zeta(-1)         = -(D-2)*pi/(24L)
    (zeta(-1) = -1/12).  For D=4 these are -pi/3 (closed) and -pi/12 (open) -- the 4x
    ratio CONSTRUCT-LUSCHER found.  Depends ONLY on D-2 => universal.

(2) 2D pure YM is DEGENERATE for rigidity: D-2 = 0 transverse directions => no transverse
    modes, no Luscher term, and the (d^2 X)^2 rigidity multiplies zero fields.  kappa is
    undefined/vacuous in D=2.

(3) A physical thick tube has rigidity kappa ~ sigma * <r^2>_perp, where <r^2> is the
    second moment of the transverse stress profile.  <r^2> is PROFILE-DEPENDENT: at the
    same nominal width w it differs by O(10) across profiles => the O(1) coefficient of
    kappa is set by the microscopic core shape, NOT by universality.
"""
import numpy as np
from math import pi

zeta_m1 = -1.0 / 12.0   # zeta(-1)

print("=" * 72)
print("(1) Luscher term: zero-point of D-2 transverse bosons -> depends ONLY on c = D-2")
print("=" * 72)
for D in (3, 4, 5, 6):
    L = 1.0
    closed = (D - 2) * (2 * pi / L) * zeta_m1
    open   = (D - 2) * (pi / (2 * L)) * zeta_m1
    print(f"  D={D}: c=D-2={D-2}   E0_closed = {closed:+.6f} (= -(D-2)pi/6L)   "
          f"E0_open = {open:+.6f} (= -(D-2)pi/24L)")
print("  -> coefficient is (D-2)*const: UNIVERSAL (same for any confining string).")

print()
print("=" * 72)
print("(2) 2D pure YM: D-2 = 0 transverse directions -> rigidity DEGENERATE")
print("=" * 72)
D = 2
c = D - 2
print(f"  D=2: c = D-2 = {c}  => # transverse modes = {c}")
print(f"       Luscher E0 = {c}*(2pi/L)*zeta(-1) = {c*(2*pi)*zeta_m1:.3e} (trivial)")
print("       kappa multiplies (d^2 X)^2 of ZERO transverse fields -> kappa undefined.")
print("       (The 2D 'flux tube' is the whole worldsheet; no transverse profile.)")

print()
print("=" * 72)
print("(3) Core rigidity kappa ~ sigma*<r^2>_perp: <r^2> is PROFILE-DEPENDENT")
print("=" * 72)
r = np.linspace(0.0, 60.0, 600001)

def second_moment(rho):
    d2r = 2.0 * pi * r
    num = np.trapezoid(rho * r**2 * d2r, r)
    den = np.trapezoid(rho * d2r, r)
    return num / den

w = 1.0
prof_gauss = np.exp(-r**2 / w**2)          # <r^2> = w^2
prof_disk  = (r <= w).astype(float)        # <r^2> = w^2/2
prof_exp   = np.exp(-r / w)                # <r^2> = 6 w^2
for name, rho, expected in [("Gaussian", prof_gauss, 1.0),
                            ("flat disk", prof_disk, 0.5),
                            ("exponential", prof_exp, 6.0)]:
    m = second_moment(rho) / w**2
    print(f"  {name:12s}: <r^2>/w^2 = {m:.4f}   (expected {expected})")
print("  -> at the SAME nominal width w, <r^2> spans 0.5 .. 6.0 (a 12x spread).")
print("     kappa = sigma*<r^2>*c_profile : the O(1) c_profile is set by the microscopic")
print("     core shape => NON-UNIVERSAL.  (The universal log slope of w^2(R) is kappa-free.)")

print()
print("=" * 72)
print("VERDICT")
print("=" * 72)
print("  Luscher term: UNIVERSAL (only c=D-2).")
print("  Rigidity kappa: NON-UNIVERSAL microscopic LEC (profile-dependent), undefined in 2D.")
print("  => kappa does NOT transfer from exact 2D / thin models to 4D;")
print("     the 4D physical kappa needs the microscopic core = H1.")
