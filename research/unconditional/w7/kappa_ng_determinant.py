"""
kappa_ng_determinant.py  --  KAPPA-NG-DETERMINANT witness  (2026-11-24)

Question: is the rigid-string rigidity kappa computable from the EXACT
Nambu-Goto (Luescher-Weisz) determinant, and does it land in [0.14,0.19]?

Model under attack (CONSTRUCT-MASTERFIELD):
    rigid closed ring:  E(rho) = 2*pi*sigma*rho + pi*kappa/rho,
    interior minimum at rho*=sqrt(kappa/(2 sigma)),  m_0++/sqrt(sigma)=2*pi*sqrt(2*kappa).
Fitting lattice m/sqrt(sigma) in [3.3,3.9] gives kappa in [0.138,0.193].
Pure NG (kappa=0): E = 2*pi*sigma*rho, monotone -> NO interior minimum [T].

We compute what the NG determinant actually supplies:
  (A) tree-level quadratic NG action: free 2D CFT, dispersion linear in k^2,
      NO k^4 term  =>  tree-level rigidity kappa_tree = 0.
  (B) exact 4-derivative NG quartic V_4 (light-cone form) -- the Luescher-Weisz
      quartic.  It is NOT a positive (d^2 X)^2 rigidity; it is the standard
      quartic whose 1-loop effect is the Luescher term.
  (C) the 1-loop (Luescher) correction to a closed ring is
      E_L = -pi (D-2)/(6 L),  L = 2 pi rho  =>  coefficient of 1/rho is
      -(D-2)/(12) = -1/6 (D=4):  NEGATIVE  =>  DESTABILIZING.
      In the MASTERFIELD normalization (stabiliser +pi*kappa/rho) the NG
      "effective kappa" is  kappa_NG = (-1/6)/pi = -1/(6 pi) < 0.
  (D) rigid minimum exists iff the total 1/rho coefficient > 0, i.e.
      pi*kappa_total - 1/6 > 0;  pure NG gives -1/6 < 0 (no minimum, [T]).

CONCLUSION: kappa_NG is negative (wrong sign) and NOT in [0.14,0.19].
The stabilizing kappa ~ 0.19 must be an EXTRA, microscopic (flux-tube core)
extrinsic-curvature coupling -- an irrelevant dim-4 operator (WORLDSHEET-CFT
[T]) set by the UV completion, not by the NG IR expansion.  kappa stays = H1.

Numpy + sympy only.   Run:  py -3.12 -I kappa_ng_determinant.py
"""

import numpy as np
import sympy as sp

BANNER = "=" * 74
print(BANNER); print("KAPPA-NG-DETERMINANT: is kappa derivable from the exact NG determinant?"); print(BANNER)

# ----------------------------------------------------------------------
# (A) Tree-level quadratic NG: free 2D CFT, no k^4.
#     S_quad = (sigma/2) int [(d_t X)^2 + (d_s X)^2]  (Euclidean).
#     For X = exp(i(w t + k s)) the kernel is (w^2 + k^2): linear in k^2,
#     coefficient of k^4 is identically 0  =>  kappa_tree = 0.
# ----------------------------------------------------------------------
w, k = sp.symbols('w k', real=True)
kernel_tree = w**2 + k**2            # quadratic action kernel
kappa_tree = sp.diff(kernel_tree, k, 4)   # coefficient of k^4
print("\n(A) tree-level NG quadratic kernel:", kernel_tree,
      "   d^4/dk^4 =", kappa_tree, " => kappa_tree =", kappa_tree)
assert kappa_tree == 0

# ----------------------------------------------------------------------
# (B) Exact 4-derivative NG quartic.  sqrt(det(eta + dX.dX)), M_ab = d_aX.d_bX.
#     A = M00, B = M11, C = M01.  Expand sqrt((1+A)(1+B)-C^2) to quartic and
#     subtract the quadratic part (1/2)(A+B): the remainder must be
#         V4 = -(1/8)(A-B)^2 - (1/2) C^2 .
# ----------------------------------------------------------------------
A, B, C = sp.symbols('A B C', real=True)
det = (1 + A) * (1 + B) - C**2
# series of sqrt(det) in the expansion parameter: treat A,B,C ~ O(eps^2) so that
# sqrt = 1 + (1/2)(A+B+AB-C^2) - (1/8)(A+B+AB-C^2)^2 + O(eps^6)
eps = sp.symbols('eps')
Ab, Bb, Cb = eps*A, eps*B, eps*C
det_eps = (1 + Ab) * (1 + Bb) - Cb**2
ser = sp.series(sp.sqrt(det_eps), eps, 0, 5).removeO()
quad = (1/2) * (A + B)               # expected O(eps^2) part
ser_c = sp.expand(ser)
# coefficient of eps^4 = V4 (since A,B,C carry eps)
V4_series = sp.expand(ser_c.coeff(eps, 2))
V4_claim = sp.expand(-sp.Rational(1,8)*(A - B)**2 - sp.Rational(1,2)*C**2)
print("\n(B) V4 from expansion :", V4_series)
print("    V4 = -(1/8)(A-B)^2 - (1/2)C^2 :", V4_claim)
print("    difference (must be 0) :", sp.simplify(V4_series - V4_claim))
assert sp.simplify(V4_series - V4_claim) == 0

# light-cone check: A-B = 4R, C = P-Q  =>  V4 = -2 R^2 - (1/2)(P-Q)^2
R, P, Q = sp.symbols('R P Q', real=True)
V4_full_lc = sp.expand(-sp.Rational(1,8)*(4*R)**2 - sp.Rational(1,2)*(P - Q)**2)
print("    light-cone: V4 = -2 R^2 - (1/2)(P-Q)^2 :", V4_full_lc,
      "   (R=d_+X.d_-X, P=d_+X^2, Q=d_-X^2)")

# ----------------------------------------------------------------------
# (C) 1-loop Luescher correction to a closed ring; its 1/rho coefficient.
#     Closed-string Luescher term:  E_L = - pi (D-2)/(6 L),  L = 2 pi rho.
#     => E_L = -(D-2)/(12 rho).   D = 4  =>  coefficient = -2/12 = -1/6.
# ----------------------------------------------------------------------
D = 4
# E_L = -pi (D-2)/(6 L),  L = 2 pi rho  =>  E_L = -(D-2)/(12 rho)  (the pi's cancel)
coeff_luescher = -(D - 2) / 12                  # 1/rho coefficient = -(D-2)/12 = -1/6
print("\n(C) closed-ring Luescher 1/rho coefficient (D=4):",
      coeff_luescher, "  ( = -(D-2)/12 = -1/6 )")
kappa_NG = coeff_luescher / np.pi               # in MASTERFIELD normalisation +pi*kappa/rho
print("    => effective kappa_NG = (-1/6)/pi =", kappa_NG,
      "  < 0  (WRONG SIGN, destabilising)")
print("    needed band from lattice: kappa in [0.138, 0.193];  kappa_NG in band? ",
      (0.138 <= kappa_NG <= 0.193))

# ----------------------------------------------------------------------
# (D) Minimum existence: total 1/rho coefficient = pi*kappa_total - 1/6.
#     E(rho) = 2 pi sigma rho + (pi*kappa_total - 1/6)/rho.
#     Interior minimum iff the 1/rho coefficient > 0.
# ----------------------------------------------------------------------
sigma = 1.0
print("\n(D) minimum existence test  E(rho)=2 pi sigma rho + c/rho, c = pi*kappa - 1/6")
for kappa in [0.0, -1/(6*np.pi), 0.138, 0.193]:
    c = np.pi * kappa - 1/6
    has_min = c > 0
    rho_star = np.sqrt(c / (2*np.pi*sigma)) if has_min else float('nan')
    Emin = 2*np.pi*sigma*rho_star + c/rho_star if has_min else float('nan')
    m_over_sqrt_sigma = Emin / np.sqrt(sigma) if has_min else float('nan')
    print(f"    kappa={kappa:+.5f}: c={c:+.4f}  interior min? {has_min}"
          + (f"  rho*={rho_star:.4f}, m/sqrt(sigma)={m_over_sqrt_sigma:.3f}" if has_min else ""))

# monotonicity of pure NG + Luescher (kappa=0): E'(rho)=2 pi sigma + (1/6)/rho^2 > 0
rho = np.linspace(0.2, 3.0, 7)
c0 = -1/6
E_pure = 2*np.pi*sigma*rho + c0/rho
print("    pure NG+Luescher E(rho):", np.round(E_pure, 4), " strictly increasing? ",
      np.all(np.diff(E_pure) > 0))

print("\n" + BANNER)
print("VERDICT: kappa_NG = -1/(6 pi) = %.5f  (wrong sign, not in [0.14,0.19])." % kappa_NG)
print("NG supplies a NEGATIVE 1/rho (destabilising) term; the +kappa that")
print("stabilises the neutral ring is an EXTRA microscopic (core) coupling,")
print("an irrelevant dim-4 operator => kappa stays = H1.  No QED, Clay 0%.")
print(BANNER)
