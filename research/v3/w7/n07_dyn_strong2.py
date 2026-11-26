r"""
n07_dyn_strong2.py -- corrected instrument for N07-STRONG.

The W_0 = \int delta(a-<v,xi>) e^{-S} dxi  route has a *boundary term* (the flat slice
cross-section depends on a, and the flat density does not vanish at the chart edge), so
  -d^2/da^2 log W_0  =  D_dyn + BOUNDARY  !=  D_dyn .
Use instead the corpus's exact *conditional* identity (Lemma G1 / ch02):
  G''(v;a)   = E[ d^2V/da^2 | a ] - Var( dV/da | a ) ,  V = S + Phi
  D_dyn(a)   = E[ d^2S/da^2 | a ] - Var( dS/da | a )  =: E[H|a] - Var(G|a)
No boundary term.  In the one-plaquette model S = beta(1-f), so
  D_dyn = -beta <f_aa | a> + beta^2 [ <f_a|a>^2 - <f_a^2|a> ]_Haar + O(beta^3)
   ^O(u) term                  ^ O(u^2) term (<= 0)
=> leading order is O(beta) = O(u); sign = -sign(<f_aa|a>).

<.|a> always denotes the Haar-conditional average on the slice {x3=a, |x|<=pi}.
"""
import mpmath as mp

mp.mp.dps = 30
N = 2
TWO_PI = 2 * mp.pi
EDGE = mp.pi            # SU(2) fundamental domain in the exponential chart


def f_sl(a, rho, th):
    r = mp.sqrt(a * a + rho * rho)
    if r == 0:
        return mp.cos(th / 2)
    return mp.cos(r / 2) * mp.cos(th / 2) - (a / r) * mp.sin(r / 2) * mp.sin(th / 2)


def rmax(a):
    return mp.sqrt(EDGE ** 2 - a * a)


def cond(Xfun, a, th, beta=mp.mpf(0)):
    """Haar-conditional average on the slice, weight e^{-beta(1-f)}."""
    num = lambda rho: rho * mp.sin(mp.sqrt(a * a + rho * rho) / 2) ** 2 / (a * a + rho * rho) \
        * mp.e ** (-beta * (1 - f_sl(a, rho, th))) * Xfun(rho)
    den = lambda rho: rho * mp.sin(mp.sqrt(a * a + rho * rho) / 2) ** 2 / (a * a + rho * rho) \
        * mp.e ** (-beta * (1 - f_sl(a, rho, th)))
    return mp.quad(num, [0, rmax(a)]) / mp.quad(den, [0, rmax(a)])


def f_aa(a, th, h=mp.mpf("1e-12")):
    """d^2f/da^2 for f on the slice: f = cos(r/2)cos(th/2) - (a/r) sin(r/2) sin(th/2)."""
    g = lambda z: f_sl(z, None, th) if False else None
    return None


def faa(rho, a, th, h=mp.mpf("1e-8")):
    return (f_sl(a + h, rho, th) - 2 * f_sl(a, rho, th) + f_sl(a - h, rho, th)) / (h * h)


def fa(rho, a, th, h=mp.mpf("1e-8")):
    return (f_sl(a + h, rho, th) - f_sl(a - h, rho, th)) / (2 * h)


def main():
    out = []

    def P(*xs):
        s = " ".join(str(x) for x in xs)
        out.append(s)
        print(s, flush=True)

    P("# n07_dyn_strong2.log -- D_dyn(a,theta) leading order in u  [exact SU(2) 1-plaquette]")
    P("# N=2, u = beta/(2N^2) = beta/8  [T] (X01),  a0_amp = (2N)^-2 = 1/16 [T] (Z04 s=1)")
    P("# D_dyn = -beta*Faa + beta^2*V2 + O(beta^3),  Faa = <f_aa|a>_Haar, V2 = <f_a|a>^2-<f_a^2|a>")
    P("")

    P("=== (A) Faa(a,theta) = <d^2f/da^2 | a>_Haar   (sign of D_dyn at O(u) is -sign(Faa)) ===")
    P("%9s %7s %20s %20s %20s" % ("theta", "a", "Faa", "-beta*Faa@1/48", "V2"))
    tab = {}
    for th in [mp.mpf(0), mp.pi / 4, mp.pi / 2, 3 * mp.pi / 4, mp.pi]:
        for a in ["0.25", "0.5", "1.0", "1.5", "2.0", "2.5"]:
            a = mp.mpf(a)
            Faa = cond(lambda rho: faa(rho, a, th), a, th)
            Faa2 = cond(lambda rho: fa(rho, a, th) ** 2, a, th)
            m1 = cond(lambda rho: fa(rho, a, th), a, th)
            V2 = m1 ** 2 - Faa2
            b1 = mp.mpf(1) / 48
            tab[(float(th), float(a))] = (Faa, V2)
            P("%9s %7s %20s %20s %20s" % (mp.nstr(th, 6), mp.nstr(a, 4),
                                           mp.nstr(Faa, 10), mp.nstr(-b1 * Faa, 8),
                                           mp.nstr(V2, 10)))
        P("")

    P("=== (B) beta-linearity check of the exact conditional D_dyn at theta=pi, a=1.5 ===")
    P("%9s %14s %18s %18s" % ("beta", "D_dyn exact", "-beta*Faa", "ratio"))
    th, a = mp.pi, mp.mpf("1.5")
    Faa = cond(lambda rho: faa(rho, a, th), a, th)
    for beta in [mp.mpf("0.01"), mp.mpf("0.02"), mp.mpf("0.04"), mp.mpf("0.1")]:
        EH = -beta * cond(lambda rho: faa(rho, a, th), a, th, beta)
        m1 = cond(lambda rho: fa(rho, a, th), a, th, beta)
        m2 = cond(lambda rho: fa(rho, a, th) ** 2, a, th, beta)
        D = EH - beta ** 2 * (m2 - m1 ** 2)
        P("%9s %14s %18s %18s" % (mp.nstr(beta, 4), mp.nstr(D, 8), mp.nstr(-beta * Faa, 8),
                                  mp.nstr(D / (-beta * Faa), 8)))
    P("")

    P("=== (C) parity check: theta=0 slice is even in a => Faa(0) analytic?  and Faa at a=0 ===")
    for th in [mp.mpf(0), mp.pi]:
        Faa0 = cond(lambda rho: faa(rho, mp.mpf("1e-6"), th), mp.mpf("1e-6"), th)
        P("theta=%s  Faa(a->0) = %s" % (mp.nstr(th, 5), mp.nstr(Faa0, 10)))
    P("")

    P("=== (D) window numbers, beta1 = 1/48, N=2, 4D single-link mode => x6 plaquettes ===")
    b1 = mp.mpf(1) / 48
    P("%9s %6s %14s %16s %16s" % ("theta", "a", "Faa", "D_dyn^1plq", "D_dyn^4D(x6)"))
    for th in [mp.mpf(0), mp.pi / 2, mp.pi]:
        for a in [mp.mpf("0.5"), mp.mpf("1.5"), mp.mpf("2.0")]:
            Faa, V2 = tab[(float(th), float(a))]
            P("%9s %6s %14s %16s %16s" % (mp.nstr(th, 6), mp.nstr(a, 3), mp.nstr(Faa, 10),
                                           mp.nstr(-b1 * Faa, 8), mp.nstr(-6 * b1 * Faa, 8)))
    P("")
    P("=== (E) same in u:  D_dyn = -(2N^2 u) Faa = -8u Faa  (N=2);  u(beta1)=%.6e ===" % float(b1 / 8))
    P("===     I1 amplitude a0=(2N)^-2=1/16 ;  u/a0 = 32u  => D_dyn = -(Faa/4)(u/a0)  ===")

    open("research/v3/w7/n07_dyn_strong2.log", "w", encoding="utf-8").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
