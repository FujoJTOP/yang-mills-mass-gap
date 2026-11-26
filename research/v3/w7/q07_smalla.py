#!/usr/bin/env python3
"""
q07_smalla.py -- G(a) = -log rho(a) as an expansion about a=0.
Two independent things:
  (1) EXACT one-plaquette SU(2) marginal rho(a), beta=2.5, theta in {0, pi/2, pi}:
      G', G'', G''', G'''' at a=0 by high-order FD on log rho (+ anchors).
  (2) Formula validation on a smooth toy where the conditional measure is explicit:
      G'''  = E[V'''] - 3Cov(F,kappa) + cum3(F)
      G'''' = 6 m^4 - 12 m^2 u2 + 3 u2^2 - 4 m u3 - u4,  u_n = rho^(n)/rho
              u1=-m, u2=E[F^2-kap], u3=E[3F kap - V''' - F^3],
              u4=E[F^4 + 4F V''' + 3kap^2 - 6F^2 kap - V'''']
  (3) Fisher identity  E_rho[G''] = I(rho) = int (rho')^2/rho  (average, not pointwise).
Run: py -3.12 q07_smalla.py
"""
import numpy as np

GL_N = 700
XG, WG = np.polynomial.legendre.leggauss(GL_N)


def quad(f, lo, hi):
    t = 0.5 * (XG + 1.0) * (hi - lo) + lo
    return 0.5 * (hi - lo) * np.sum(WG * f(t))


# ---------------------------------------------------------------- one plaquette
# dU = J dx, J = 4 sin^2(r/2)/r^2 ,  S_W = beta(1 - (1/2)tr(B U)),
# (1/2)tr(B U) = cos(th/2)cos(r/2) - sin(th/2)sin(r/2)*(a/r),  r = |x|, a = x3
# transverse x1,x2 integrated out:  2 pi r dr  ->   mu(a) = 2pi int_{|a|}^{2pi} r dr e^{-S} J
BETA = 2.5
RMAX = 2.0 * np.pi


def rho_plaq(a, th):
    lo = abs(a)

    def f(r):
        half = np.cos(th / 2) * np.cos(r / 2) - np.sin(th / 2) * np.sin(r / 2) * (a / r)
        S = BETA * (1.0 - half)
        J = 4.0 * np.sin(r / 2) ** 2 / r ** 2
        return 2.0 * np.pi * r * np.exp(-S) * J

    return quad(f, lo, RMAX)


def G_plaq(a, th):
    return -np.log(rho_plaq(a, th))


def fdf(f, a, h, order):
    if order == 1:
        return (-f(a + 2 * h) + 8 * f(a + h) - 8 * f(a - h) + f(a - 2 * h)) / (12 * h)
    if order == 2:
        return (-f(a + 2 * h) + 16 * f(a + h) - 30 * f(a) + 16 * f(a - h) - f(a - 2 * h)) / (12 * h ** 2)
    if order == 3:
        return (-f(a + 3 * h) + 8 * f(a + 2 * h) - 13 * f(a + h) + 13 * f(a - h)
                - 8 * f(a - 2 * h) + f(a - 3 * h)) / (8 * h ** 3)
    if order == 4:
        return (56 * f(a) - 39 * (f(a + h) + f(a - h)) + 12 * (f(a + 2 * h) + f(a - 2 * h))
                - (f(a + 3 * h) + f(a - 3 * h))) / (6 * h ** 4)
    raise ValueError


def derivs(f, a, orders=(1, 2, 3, 4), hs=(0.04, 0.02)):
    out = {}
    for n in orders:
        v1 = fdf(f, a, hs[0], n)
        v2 = fdf(f, a, hs[1], n)
        out[n] = v2 + (v2 - v1) / (2 ** n - 1)          # Richardson
        out[("raw", n)] = (v1, v2)
    return out


# ---------------------------------------------------------------------- toy
def toy_G(a):
    k, g, c, d = 1.3, 0.7, 0.5, 0.2

    def f(y):
        V = 0.5 * k * (y - a) ** 2 + 0.25 * g * y ** 4 + c * a ** 3 / 3.0 + d * a ** 4 / 4.0
        return np.exp(-V)

    return -np.log(quad(f, -14.0, 14.0))


def toy_moments(a):
    """conditional moments at parameter a: nu_a(y) ~ exp(-V(a,y))."""
    k, g, c, d = 1.3, 0.7, 0.5, 0.2
    y = 0.5 * (XG + 1.0) * 28.0 - 14.0
    V = 0.5 * k * (y - a) ** 2 + 0.25 * g * y ** 4 + c * a ** 3 / 3.0 + d * a ** 4 / 4.0
    w = WG * np.exp(-V)
    E = lambda h: np.sum(w * h(y)) / np.sum(w)
    F = -k * (y - a) + c * a ** 2 + d * a ** 3          # dV/da
    kap = k + 2 * c * a + 3 * d * a ** 2                # d2V/da2
    V3 = 2 * c + 6 * d * a
    V4 = 6 * d
    m = E(lambda yy: F)
    A2, A3, A4 = E(lambda yy: kap), E(lambda yy: V3), E(lambda yy: V4)
    return dict(m=m, A2=A2, A3=A3, A4=A4, F=F, kap=kap, V3=V3, V4=V4, E=E, w=w)


def toy_formulas(a):
    t = toy_moments(a)
    E, F, kap, V3, V4 = t["E"], t["F"], t["kap"], t["V3"], t["V4"]
    m, A2, A3, A4 = t["m"], t["A2"], t["A3"], t["A4"]
    E = lambda h: np.sum(t["w"] * h) / np.sum(t["w"])
    u2 = E(F ** 2 - kap)
    u3 = E(3 * F * kap - V3 - F ** 3)
    u4 = E(F ** 4 + 4 * F * V3 + 3 * kap ** 2 - 6 * F ** 2 * kap - V4)
    g1 = m
    g2 = m ** 2 - u2          # G'' = -l'' = -u2 + u1^2 ,  u1 = -m
    g3 = -(u3 - 3 * (-m) * u2 + 2 * (-m) ** 3)
    g4 = -(u4 - 4 * (-m) * u3 - 3 * u2 ** 2 + 12 * m ** 2 * u2 - 6 * m ** 4)
    cum3 = E((F - m) ** 3)
    cov = E((F - m) * kap)
    return dict(G1=g1, G2=g2, G3=g3, G4=g4,
                G3_cum=A3 - 3 * cov + cum3, G4_alt=g4, u=(u2, u3, u4), cum3=cum3, cov=cov)


def main():
    print("=" * 96)
    print("(1) exact one-plaquette SU(2) marginal, beta=2.5,  G(a) = -log rho(a)")
    print("=" * 96)
    print(f"{'th':>6} {'G1(+0)':>13} {'G1(-0)':>13} {'G2(0)':>13} {'G3(+0)':>12} "
          f"{'G4(0)':>12}   stability of G4 (h=0.04/0.02/0.01)")
    for th in (0.0, np.pi / 2, np.pi):
        f = (lambda aa: G_plaq(aa, th))
        d = derivs(f, 0.0)
        dm = derivs(f, -1e-9) if False else None
        # one-sided values (asymmetry check)
        hh = 1e-3
        g1p = (f(hh) - f(0.0)) / hh
        g1m = (f(0.0) - f(-hh)) / hh
        g4a = fdf(f, 0.0, 0.02, 4)
        g4b = fdf(f, 0.0, 0.01, 4)
        print(f"{th:6.3f} {g1p:13.8f} {g1m:13.8f} {d[2]:13.8f} {d[3]:12.8f} "
              f"{d[4]:12.8f}   {g4a:11.6f}/{g4b:11.6f}  G3raw={d[('raw',3)]}")

    print()
    print("anchors (must match O28/P05/W3): th=0: G2=+0.7299241 ; th=pi: G2=+0.1359947")
    for th in (0.0, np.pi):
        f = (lambda aa: G_plaq(aa, th))
        print(f"  th={th:5.3f}  G2(0) = {fdf(f,0.0,0.02,2):+.9f}")

    print()
    print("G(a) profile, th=pi (asymmetry) and th=0:")
    for th in (0.0, np.pi):
        vals = []
        for a in (-0.4, -0.2, -0.05, 0.0, 0.05, 0.2, 0.4):
            vals.append(f"{a:+.2f}:{G_plaq(a, th) - G_plaq(0.0, th):+.6f}")
        print(f"  th={th:5.3f}  G(a)-G(0)  " + "  ".join(vals))

    print()
    print("=" * 96)
    print("(2) toy validation:  V(a,y)=k(y-a)^2/2 + g y^4/4 + c a^3/3 + d a^4/4,  k=1.3 g=0.7 c=0.5 d=0.2")
    print("    a0 = 0.3 : FD-of-G  vs  cumulant/moment formulas (u-recursion)")
    print("=" * 96)
    a0 = 0.3
    fd = derivs(toy_G, a0, hs=(0.02, 0.01))
    fm = toy_formulas(a0)
    print(f"  {'n':>2} {'FD of G':>16} {'formula':>16} {'diff':>12}")
    for n, key in ((1, "G1"), (2, "G2"), (3, "G3"), (4, "G4")):
        print(f"  {n:2d} {fd[n]:16.10f} {fm[key]:16.10f} {fd[n]-fm[key]:12.2e}")
    print(f"  G''' via A3-3Cov+cum3 = {fm['G3_cum']:.10f}  (same route)   u2,u3,u4 = {fm['u']}")
    print(f"  diagnostic: cum3(F)={fm['cum3']:+.6f}  Cov(F,kap)={fm['cov']:+.6f}")

    print()
    print("=" * 96)
    print("(3) Fisher identity  E_rho[G''] = I(rho) = int (rho')^2/rho   (+boundary)")
    print("    [3a] toy over the full line (rho has no support edge)")
    print("=" * 96)
    ys = np.linspace(-30.0, 30.0, 60001)
    for a0 in (0.0, 0.3):
        rv = np.array([toy_G(a0 + 0.0) ])
    aa = np.linspace(-12.0, 12.0, 24001)
    rv = np.array([toy_G(x) for x in aa])
    rv = np.exp(-(rv - rv.min()))
    rv = rv / np.trapezoid(rv, aa)
    G2 = np.array([fdf(toy_G, x, 0.02, 2) for x in aa])
    drho = np.gradient(rv, aa)
    Eavg = np.trapezoid(rv * G2, aa)
    Fish = np.trapezoid(drho ** 2 / rv, aa)
    m1 = np.trapezoid(rv * aa, aa)
    var = np.trapezoid(rv * (aa - m1) ** 2, aa)
    print(f"  E_rho[G2]={Eavg:+.9f}  I(rho)={Fish:+.9f}  diff={Eavg-Fish:+.2e} "
          f"| I={Fish:.6f}  1/Var={1/var:.6f}  I-1/Var={Fish-1/var:+.4e}  (CR needs full line)")
    print()
    print("    [3b] one-plaquette: full support |a|<2pi, boundary term = -(rho'(R)-rho'(-R))")
    for th in (0.0, np.pi):
        N = 4001
        aa = np.linspace(-RMAX + 1e-6, RMAX - 1e-6, N)
        rv = np.array([rho_plaq(x, th) for x in aa])
        rv = rv / np.trapezoid(rv, aa)
        G2 = np.array([fdf(lambda y: G_plaq(y, th), x, 0.02, 2) for x in aa])
        drho = np.gradient(rv, aa)
        Eavg = np.trapezoid(rv * G2, aa)
        Fish = np.trapezoid(drho ** 2 / rv, aa)
        bnd = -(drho[-1] - drho[0])
        m1 = np.trapezoid(rv * aa, aa)
        var = np.trapezoid(rv * (aa - m1) ** 2, aa)
        print(f"  th={th:5.3f}: int rho G2={Eavg:+.6f}  I={Fish:+.6f}  bnd={bnd:+.6f}"
              f"  I+bnd-Eavg={Fish + bnd - Eavg:+.2e}  | I-1/Var={Fish-1/var:+.4e} Var={var:.5f}")


if __name__ == "__main__":
    main()
