# -*- coding: utf-8 -*-
"""
R05-WCF : <w>, Var(w) for
    w(x) = sin x/(2x),  x = rho/2 in [0,pi]
    mu(x) ∝ exp(-beta(1-c cos x)) sin^2(x) / x dx ,  c = cos(theta/2)
Only the product  gam := beta*c  matters:  beta(1-c cos x) = beta(1-c) + gam(1-cos x).
So  <w>(beta,c) = W(gam),  Var(w)(beta,c) = V(gam),  gam = beta c.
[1] exact beta=0 closed forms (Cin)   [2] high-precision quadrature, independent check
[3] 1/gamma asymptotic series (sympy, Watson)   [4] min of <w>/Var over gamma in [0,inf)
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 40
OUT = []


def log(s=""):
    OUT.append(str(s))
    print(s)


# ---------------------------------------------------------------- [1] beta = 0
def Cin(z):
    # Cin(z) = gamma + ln z - Ci(z) , z > 0
    return mp.euler + mp.log(z) - mp.ci(z)


def be0():
    C2 = Cin(2 * mp.pi)
    C1 = Cin(mp.pi)
    C3 = Cin(3 * mp.pi)
    C4 = Cin(4 * mp.pi)
    Z = C2 / 2                      # int_0^pi sin^2 x / x dx
    N1 = 3 * (C3 - C1) / 4          # int_0^pi sin^3 x / x^2 dx
    N2 = C4 - C2                    # int_0^pi sin^4 x / x^3 dx
    ew = N1 / (2 * Z)               # <w> = (1/2) N1/Z
    ew2 = N2 / (4 * Z)              # <w^2> = (1/4) N2/Z
    var = ew2 - ew ** 2
    return dict(Cin1=C1, Cin2=C2, Cin3=C3, Cin4=C4, Z=Z, N1=N1, N2=N2,
                ew=ew, ew2=ew2, var=var, ratio=ew / var, two=2 * var / ew)


# ------------------------------------------------- [2] quadrature for gamma>0
def sinc(x):
    if x == 0:
        return mp.mpf(1)
    return mp.sin(x) / x


def kern(gam, order):
    """e^{gam cos x} * sin^order x / x^{order-1} ; order=2 -> A0, 3 -> A1, 4 -> A2."""
    def f(x):
        s = mp.sin(x)
        w = mp.e ** (gam * mp.cos(x))
        return w * s * sinc(x) ** (order - 1)
    return f


def Amat(gam, split_scale=3.0):
    """A0,A1,A2 by deterministic quadrature (split near the peak at x=0)."""
    if gam == 0:
        x0 = mp.mpf(1) / 4
    else:
        x0 = min(mp.pi, split_scale / mp.sqrt(gam))
    pts = sorted({mp.mpf(0), x0, mp.pi / 2, mp.pi})
    vals = []
    for order in (2, 3, 4):
        f = kern(gam, order)
        vals.append(mp.quad(f, pts, maxdegree=14))
    return vals


def WV(gam):
    A0, A1, A2 = Amat(gam)
    ew = A1 / A0 / 2
    ew2 = A2 / A0 / 4
    var = ew2 - ew ** 2
    return ew, ew2, var


# ------------------------------------------ [3] 1/gamma series via Watson lemma
def series(nord=5):
    z, u = sp.symbols("z u", positive=True)
    y = z * sp.sqrt(u)                      # x = z sqrt(u), u = 1/gamma
    wser = sp.series(z ** 2 / 2 - (1 - sp.cos(y)) / u, u, 0, nord).removeO()
    expfac = sp.series(sp.exp(u * wser), u, 0, nord).removeO()
    res = []
    for order in (2, 3, 4):
        g = sp.series(sp.sin(y) ** order / y ** (order - 1) * sp.sqrt(u), u, 0, nord)
        integ = sp.expand(expfac * g.removeO())
        poly = sp.Poly(integ, u)
        acc = []
        for k in range(nord):
            ck = sp.expand(poly.coeff_monomial(u ** k))
            pk = sp.Poly(ck, z)
            tot = 0
            for (n,), co in pk.terms():
                tot += co * sp.integrate(z ** n * sp.exp(-z ** 2 / 2), (z, 0, sp.oo))
            acc.append(sp.nsimplify(sp.simplify(tot)))
        res.append(acc)          # res[i][k] = coefficient of u^k in A_{order}
    return {2: res[0], 3: res[1], 4: res[2]}


def main():
    b = be0()
    log("=== R05-WCF : <w>, Var(w), w = sin(x)/(2x), x=rho/2 in [0,pi] ===")
    log("mpmath %s , dps=%d" % (mp.__version__, mp.mp.dps))
    log("")
    log("[1] beta=0 closed forms (Cin(z)=int_0^z (1-cos t)/t dt)")
    for k in ("Cin1", "Cin2", "Cin3", "Cin4"):
        log("   %-5s = %s" % (k, mp.nstr(b[k], 25)))
    log("   Z  = Cin(2pi)/2            = %s" % mp.nstr(b["Z"], 25))
    log("   N1 = 3/4(Cin3pi-Cin pi)    = %s" % mp.nstr(b["N1"], 25))
    log("   N2 = Cin4pi - Cin2pi       = %s" % mp.nstr(b["N2"], 25))
    log("   <w>_0   = N1/(2Z)          = %s" % mp.nstr(b["ew"], 25))
    log("   <w^2>_0 = N2/(4Z)          = %s" % mp.nstr(b["ew2"], 25))
    log("   Var_0   = <w^2>-<w>^2      = %s" % mp.nstr(b["var"], 25))
    log("   2Var_0/<w>_0               = %s" % mp.nstr(b["two"], 20))
    log("   <w>_0/Var_0                = %s" % mp.nstr(b["ratio"], 20))
    log("")

    log("[2] independent quadrature at beta=0 -> closed forms (rel. diff)")
    A0, A1, A2 = Amat(mp.mpf(0))
    e1 = A1 / A0 / 2
    e2 = A2 / A0 / 4
    log("   <w> quad  = %s   closed = %s   reldiff = %s"
        % (mp.nstr(e1, 25), mp.nstr(b["ew"], 25), mp.nstr(abs(e1 / b["ew"] - 1), 5)))
    log("   <w2> quad = %s   closed = %s   reldiff = %s"
        % (mp.nstr(e2, 25), mp.nstr(b["ew2"], 25), mp.nstr(abs(e2 / b["ew2"] - 1), 5)))
    v1 = e2 - e1 ** 2
    log("   Var quad  = %s   closed = %s   reldiff = %s"
        % (mp.nstr(v1, 25), mp.nstr(b["var"], 25), mp.nstr(abs(v1 / b["var"] - 1), 5)))
    log("")

    log("[3] 1/gamma (gamma=beta*c) asymptotic series, u = 1/gamma")
    A = series(6)
    u = sp.symbols("u")
    A0s = sum(A[2][k] * u ** k for k in range(6))
    A1s = sum(A[3][k] * u ** k for k in range(6))
    A2s = sum(A[4][k] * u ** k for k in range(6))
    ew_s = sp.series(sp.Rational(1, 2) * A1s / A0s, u, 0, 5).removeO()
    ew2_s = sp.series(sp.Rational(1, 4) * A2s / A0s, u, 0, 5).removeO()
    var_s = sp.series(ew2_s - ew_s ** 2, u, 0, 4).removeO()
    rat_s = sp.series(ew_s / var_s, u, 0, 2).removeO()
    log("   <w>     = %s" % sp.nsimplify(ew_s))
    log("   <w^2>   = %s" % sp.nsimplify(ew2_s))
    log("   Var     = %s" % sp.nsimplify(var_s))
    log("   <w>/Var = %s" % sp.nsimplify(rat_s))
    for kk in range(5):
        log("   u^%d : A0 = %-10s A1 = %-10s A2 = %-10s"
            % (kk, A[2][kk], A[3][kk], A[4][kk]))
    log("")

    log("[3b] series vs quadrature (gamma = 10,50,200,1000)")
    for gam in (10, 50, 200, 1000):
        ew, ew2, var = WV(mp.mpf(gam))
        uu = mp.mpf(1) / gam
        f = sp.lambdify(u, ew_s, "mpmath")
        fv = sp.lambdify(u, var_s, "mpmath")
        log("   gamma=%5d  <w>=%s  ser=%s  rel=%s | Var=%s ser=%s rel=%s"
            % (gam, mp.nstr(ew, 12), mp.nstr(f(uu), 12),
               mp.nstr(abs(f(uu) / ew - 1), 3), mp.nstr(var, 12),
               mp.nstr(fv(uu), 12), mp.nstr(abs(fv(uu) / var - 1), 3)))
    log("")

    log("[4] <w>/Var over gamma in [0,inf)   (grid + golden-section refinement)")
    tab = []
    grid = [mp.mpf(0)] + [mp.mpf(10) ** (mp.mpf(i) / 8) for i in range(-16, 25)]
    for gam in grid:
        ew, ew2, var = WV(gam)
        tab.append((gam, ew, ew2, var, ew / var, 2 * var / ew))
    for row in tab:
        log("   gam=%-10s <w>=%s <w2>=%s Var=%s  <w>/Var=%s  2Var/<w>=%s"
            % (mp.nstr(row[0], 6), mp.nstr(row[1], 12), mp.nstr(row[2], 12),
               mp.nstr(row[3], 12), mp.nstr(row[4], 12), mp.nstr(row[5], 10)))

    # monotonicity / minimum: fine scan on [0,2] plus golden section below it
    def s_of(gam):
        ew, ew2, var = WV(gam)
        return ew / var
    lo, hi = mp.mpf("1e-6"), mp.mpf(2)
    best = (s_of(mp.mpf(0)), mp.mpf(0))
    for i in range(0, 41):
        gam = lo + (hi - lo) * i / 40
        s = s_of(gam)
        if s < best[0]:
            best = (s, gam)
    g0 = best[1]
    a, bb = max(mp.mpf(0), g0 - mp.mpf("0.05")), g0 + mp.mpf("0.05")
    for _ in range(60):
        m1 = a + (bb - a) * mp.mpf("0.382")
        m2 = a + (bb - a) * mp.mpf("0.618")
        if s_of(m1) < s_of(m2):
            bb = m2
        else:
            a = m1
    gm = (a + bb) / 2
    log("")
    log("   min over gam:  <w>/Var = %s at gam = %s   (beta=0 value %s)"
        % (mp.nstr(s_of(gm), 20), mp.nstr(gm, 10), mp.nstr(b["ratio"], 20)))
    log("   => 2Var/<w> max = %s  (beta=0: %s)"
        % (mp.nstr(2 / s_of(gm), 12), mp.nstr(b["two"], 12)))
    log("   s_min/2 = %s   (relaxation factor vs beta t^2 <= c)"
        % mp.nstr(s_of(gm) / 2, 12))

    with open("w7/r05_wcf.log", "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")


main()
