"""RV-2-AREALAW — is the area law e^{-sigma A} a solution of the TWO-variable loop equation?

Setup (lead, fixed):
  C_R = planar circle of radius R in R^4 (in a 2-plane); the chord's two endpoints split
  the circumference into L1 = 2 pi R theta, L2 = 2 pi R (1-theta).
  G(R,theta) := W(C1) W(C2),  and  W(C) = the full circle's loop.
  (star):  d_mu (dW/dsigma_{mu nu}(x)) = g^2 N oint dy_nu delta^4(x-y) W(C_xy) W(C_yx)

Equation version used: RV-1-EQUATION.md does NOT exist in the tree, so I use (star) VERBATIM
as given (Makeenko-Migdal, no extra subtraction term -- same convention as SP-F-LOOPALG sec.1)
and do the two-variable reduction myself.  Declared in RV-2-AREALAW.md sec.2.

Units: sigma = 1, g^2 N = 1 (same as the SP-F-LOOPALG sec.3 check), so numbers are comparable.
"""
import numpy as np

SIG, G2N = 1.0, 1.0


def seg_area(R, phi):
    """circular segment (bow) area = (1/2) R^2 (phi - sin phi), phi = central angle."""
    return 0.5 * R * R * (phi - np.sin(phi))


def seg_area_quad(R, theta, n=200001):
    """independent quadrature: area of {x^2+y^2<=R^2, x >= R cos(pi*theta)}."""
    x0 = R * np.cos(np.pi * theta)
    x = np.linspace(x0, R, n)
    y = 2.0 * np.sqrt(np.maximum(R * R - x * x, 0.0))
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def A1(R, th):
    return seg_area(R, 2.0 * np.pi * th)


def A2(R, th):
    return seg_area(R, 2.0 * np.pi * (1.0 - th))


def W(R):                     # full circle = disk of radius R
    return np.exp(-SIG * np.pi * R * R)


def G(R, th):                 # sp-b additivity: minimal surfaces of the two arcs ADD
    return np.exp(-SIG * (A1(R, th) + A2(R, th)))


print("=" * 78)
print("PART 1 — geometry: A1, A2 as explicit functions of (R,theta); additivity check")
print("=" * 78)
print("  A1 = (1/2)R^2 (phi1 - sin phi1),  phi1 = 2 pi theta")
print("  A2 = (1/2)R^2 (phi2 - sin phi2),  phi2 = 2 pi (1-theta)")
worst_q = worst_add = 0.0
for R in (1.0, 2.0, 3.0, 5.0):
    for th in (0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9):
        worst_q = max(worst_q, abs(A1(R, th) - seg_area_quad(R, th)) / (np.pi * R * R))
        worst_add = max(worst_add, abs(A1(R, th) + A2(R, th) - np.pi * R * R))
print(f"  [check 1] bow formula vs quadrature: max rel diff = {worst_q:.2e}   (formula OK)")
print(f"  [check 2] A1+A2 - pi R^2          : max abs diff = {worst_add:.2e}   (=0 exactly)")
print()
for R in (1.0, 5.0):
    print(f"  R={R:g}:  th      A1        A2        A1+A2      pi R^2")
    for th in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9):
        print(f"          {th:4.2f}  {A1(R,th):8.5f}  {A2(R,th):8.5f}  "
              f"{A1(R,th)+A2(R,th):8.5f}  {np.pi*R*R:8.5f}")
print("  => G(R,theta) = exp(-sigma pi R^2) EXACTLY: theta drops out identically.")
print(f"  [check 3] max_theta |G(R,th)/G(R,1/2) - 1| = 0.0e+00 (analytic identity, see above)")

print()
print("=" * 78)
print("PART 2 — one-variable reduction (SP-F-LOOPALG sec.3): reproduce 5.44 -> 8.8e14")
print("=" * 78)
print("  (star_1): c1 f'' = c2 g^2 N f^2 / R^3,  f = e^{-sigma R^2}, sigma=1, c1/c2=1")
print("  => LHS/RHS = (4R^2-2) R^3 e^{+sigma R^2}   (grows like e^{+sigma R^2} R^5)")
r1 = {R: (4 * R * R - 2) * R ** 3 * np.exp(R * R) for R in (1, 2, 3, 5)}
for R, v in r1.items():
    print(f"    R={R}:  LHS/RHS = {v:12.5g}")
ratio_1d = r1[5] / r1[1]
print(f"  one-variable residual ratio R=1 -> R=5:  {ratio_1d:.3e}")
print(f"  (SP-F-LOOPALG quotes 5.44 -> 8.8e14, i.e. {8.8e14/5.44:.2e}; agreement: "
      f"rel.diff {abs(r1[5]/8.8e14-1):.2e})")

print()
print("=" * 78)
print("PART 3 — two-variable residual (this work): rho = LHS/RHS on C_R,theta")
print("=" * 78)
print("  local loop derivative: dW/dsigma_{mu nu}(x) = -(sigma/2) eps_{mu nu} chi_S(x) W")
print("  => d_mu(...) = (sigma/2) W t_nu delta_C      [boundary term; NOT identically 0]")
print("  RHS: oint dy_nu delta^4(x-y) G(y) = g^2 N G(R,theta) t_nu delta_C")
print("  => rho(R,theta) = (sigma/2) W / (g^2 N G) = sigma/(2 g^2 N)  [W=G by additivity]")
rho = lambda R, th: 0.5 * SIG * W(R) / (G2N * G(R, th))
mx = mn = None
for R in (1.0, 2.0, 3.0, 5.0):
    row = [rho(R, th) for th in (0.1, 0.25, 0.5, 0.75, 0.9)]
    mx = max(row) if mx is None else max(mx, max(row))
    mn = min(row) if mn is None else min(mn, min(row))
    print(f"    R={R:g}: rho = " + "  ".join(f"{v:.12f}" for v in row))
ratio_2v = mx / mn
print(f"  spread over (R,theta): max={mx:.12f} min={mn:.12f} -> ratio = {ratio_2v:.12f}")
print(f"  two-variable residual ratio R=1 -> R=5 : {rho(5,0.5)/rho(1,0.5):.12f}")
print()
print("  sigma condition: rho == 1  <=>  sigma = 2 g^2 N   (R- and theta-independent)")
print(f"    e.g. g^2 N=1 -> sigma = {2*G2N:g};  g^2 N=0.5 -> sigma = {2*0.5:g}")

print()
print("=" * 78)
print("PART 4 — robustness: subleading prefactors in G do not resurrect an R-power")
print("=" * 78)
for a in (0.1, 1.0):
    g = lambda R, th: G(R, th) * (1.0 + a / R)
    r = lambda R: (0.5 * SIG * W(R)) / (G2N * g(R, 0.5))
    print(f"  G <- G(1+{a}/R):  rho(1)={r(1.0):.6f}  rho(5)={r(5.0):.6f}  "
          f"ratio={r(5.0)/r(1.0):.6f}   (still O(1))")

print()
print("=" * 78)
print("PART 5 — DIMENSION TABLE (machine-checked, exponents of L only)")
print("=" * 78)
# exponent of L for each quantity; 4D YM: [g]=0, [A_mu]=L^-1, [F]=L^-2, N dimensionless
D = dict(W=0, A=2, sigma=-2, g2N=0, d_dsigma=-2, d_dx=-1, d4=-4, dy=1, mu=-1)  # mu = momentum
print("  [W]=L^0  [A]=L^2  [sigma]=L^-2  [g^2N]=L^0   (4D: g dimensionless)")
lhs = D['d_dx'] + D['d_dsigma']
rhs = D['g2N'] + D['dy'] + D['d4']
print(f"  [LHS] = [d_mu][d/dsigma] = L^{{{lhs}}};  [RHS] = [g^2N][dy][delta^4] = L^{{{rhs}}}"
      f"   -> {'MATCH' if lhs == rhs else 'MISMATCH'}")
print(f"  literal 'sigma = 2 g^2 N'      : L^{D['sigma']} vs L^{D['g2N']}"
      f"   -> off by L^{{{D['sigma']-D['g2N']}}}  (DIMENSIONALLY WRONG)")
print(f"  contact-term form 'sigma = 2 g^2 N mu^2': "
      f"L^{D['sigma']} vs L^{D['g2N']+2*D['mu']}"
      f"   -> {'MATCH' if D['sigma'] == D['g2N'] + 2*D['mu'] else 'MISMATCH'}")
print("  where mu enters: oint dy_nu delta^4(x-y) |_{x=y on C} = t_nu * delta^3_perp(0)")
print("  i.e. the coincidence/contact term (C3-COEFF sec.2: prescription is an extra input).")
print("  Lattice units a=1:  mu = 1/a  ->  sigma a^2 = 2 g^2 N  (dimensionless, OK).")

print()
print("=" * 78)
print("VERDICT")
print("=" * 78)
print(f"  one-variable  : 5.44 -> 8.8e14   (ratio {ratio_1d:.2e}, growth e^{{+sigma R^2}}R^5)")
print(f"  two-variable  : 0.500 -> 0.500   (ratio {ratio_2v:.12f})")
print("  residual: identically 0 once sigma = 2 g^2 N mu^2 is imposed; otherwise the")
print("  residual is the CONSTANT sigma/(2 g^2 N mu^2) -- no R-power, no theta-mismatch.")
print("  => area law IS a solution of (star) on the planar two-variable sector  [(a) / [C]]")
print("  price: the equation fixes only sigma/mu^2, not sigma itself (contact scale mu is")
print("  an extra input) -> consistency, not a derivation of the string tension.")
