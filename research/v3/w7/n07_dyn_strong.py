r"""
n07_dyn_strong.py  --  N07-STRONG:  D_dyn = E[H]-Var(G) in the strong-coupling window.

Exact SU(2) one-plaquette model in the flat exponential chart (the corpus's L=1 anchor,
H1C / D1 / G4 machinery):  U_l = exp(i x.sigma/2),  |x| <= pi  (fundamental domain),
dU = J(r) d^3x,  J(r) = (sin(r/2)/r)^2,  Phi = -log J.

Free link U = U_l, background holonomy B = e^{i theta sigma_3/2} fixed:
    U_P = U B ,      S = beta (1 - (1/2) Re tr(U_P))  =: beta (1 - f) ,   f = (1/2) Re tr(U B)
    f(x;theta) = cos(r/2) cos(theta/2) - (x3/r) sin(r/2) sin(theta/2) ,  r=|x|

Coarse mode: v = e3 on that link, a = x3.  Slice x3 = a  ->  1D polar quadrature over
rho = sqrt(x1^2+x2^2) in [0, sqrt(pi^2-a^2)],  r = sqrt(a^2+rho^2).

Objects (exact, no fitting):
    W_0(a)   = \int dx1 dx2  e^{-beta(1-f)}                 (flat Lebesgue fibre, S only)
    W(a)     = \int dx1 dx2  J(r) e^{-beta(1-f)}            (true Haar fibre)
    D_dyn(a) = -d^2/da^2 log W_0(a)     [exact identity: = E[H] - Var(G)]
    Gpp(a)   = -d^2/da^2 log W(a)       [= D_dyn + W_full]

Leading-order prediction (cumulant on Haar, beta->0):
    D_dyn = -(beta^2/2) d^2/da^2 <f^2 | a>_Haar + O(beta^3)  =: -(beta^2/2) A2 + ...
so the leading power is beta^2 = (2N^2 u)^2 and the sign is -sign(A2).

Tags: [T] identities, [E] numbers produced here.
"""
import sys
import mpmath as mp

mp.mp.dps = 30


def f_slice(a, rho, th):
    """f on the slice x3=a, polar radius rho (x1^2+x2^2=rho^2)."""
    r2 = a * a + rho * rho
    r = mp.sqrt(r2)
    if r == 0:
        return mp.cos(th / 2)
    return mp.cos(r / 2) * mp.cos(th / 2) - (a / r) * mp.sin(r / 2) * mp.sin(th / 2)


def rmax(a):
    return mp.sqrt(mp.pi ** 2 - a * a)


def J(r):
    if r == 0:
        return mp.mpf(1) / 4
    return (mp.sin(r / 2) / r) ** 2


def W0(a, beta, th):
    """flat-Lebesgue fibre integral 2*pi*int rho drho exp(-beta(1-f))"""
    g = lambda rho: rho * mp.e ** (-beta * (1 - f_slice(a, rho, th)))
    return 2 * mp.pi * mp.quad(g, [0, rmax(a)])


def Whaar(a, beta, th):
    g = lambda rho: rho * J(mp.sqrt(a * a + rho * rho)) * mp.e ** (-beta * (1 - f_slice(a, rho, th)))
    return 2 * mp.pi * mp.quad(g, [0, rmax(a)])


def F2(a, th):
    """<f^2 | a>_Haar"""
    num = lambda rho: rho * J(mp.sqrt(a * a + rho * rho)) * f_slice(a, rho, th) ** 2
    den = lambda rho: rho * J(mp.sqrt(a * a + rho * rho))
    return mp.quad(num, [0, rmax(a)]) / mp.quad(den, [0, rmax(a)])


def F1(a, th):
    """<f | a>_Haar  (should be 0 for theta=0 by parity; used as a check)"""
    num = lambda rho: rho * J(mp.sqrt(a * a + rho * rho)) * f_slice(a, rho, th)
    den = lambda rho: rho * J(mp.sqrt(a * a + rho * rho))
    return mp.quad(num, [0, rmax(a)]) / mp.quad(den, [0, rmax(a)])


def cumulant_split(a, beta, th, h=mp.mpf("1e-3")):
    """E[H], Var(G) directly at the given beta (exact derivatives, finite-diff in a)."""
    def lp(W):  # log W(a)
        return mp.log(W(a, beta, th))

    d1 = (lp(W0) - lp(W0)) if False else None
    # d/d of log W0 and second derivative
    Lm = lp(lambda z: W0(z, beta, th))
    return None


def main():
    out = []

    def P(*xs):
        line = " ".join(str(x) for x in xs)
        out.append(line)
        print(line, flush=True)

    N = 2
    a0amp = mp.mpf(1) / (2 * N) ** 2  # I1 amplitude floor a0 = (2N)^-2 = 1/16  [T] (Z04 s=1)

    P("# n07_dyn_strong.log -- D_dyn leading order in u, exact SU(2) 1-plaquette (flat chart)")
    P("# N=2,  I1 amplitude a0=(2N)^-2 =", mp.nstr(a0amp, 8), " ; u = beta/(2N^2) = beta/8  [T] (X01)")
    P("")

    # ---------- (0) sanity: <f|a> = 0 for theta=0 by parity ----------
    for th in [mp.mpf(0), mp.pi / 2, mp.pi]:
        vals = [F1(mp.mpf(a), th) for a in (mp.mpf("0.3"), mp.mpf(1.0), mp.mpf(2.0))]
        P("check <f|a>_Haar  theta=%s :" % mp.nstr(th, 6), [mp.nstr(v, 6) for v in vals])
    P("")

    # ---------- (1) A2 = d^2/da^2 <f^2|a>_Haar ----------
    P("=== (1) A2(a,theta) = d^2/da^2 <f^2|a>_Haar   [exact quadrature, dps=30] ===")
    P("%8s %8s %22s %22s" % ("theta", "a", "<f^2|a>", "A2"))
    A2tab = {}
    for th in [mp.mpf(0), mp.pi / 2, mp.pi]:
        for a in [mp.mpf("0.0"), mp.mpf("0.5"), mp.mpf("1.0"), mp.mpf("1.5"), mp.mpf("2.0")]:
            if a >= mp.pi - mp.mpf("0.05"):
                continue
            val = F2(a, th)
            A2 = mp.diff(lambda z: F2(z, th), a, 2)
            A2tab[(float(th), float(a))] = A2
            P("%8s %8s %22s %22s" % (mp.nstr(th, 6), mp.nstr(a, 4),
                                     mp.nstr(val, 12), mp.nstr(A2, 12)))
    P("")

    # ---------- (2) beta-series: D_dyn / beta^2 -> -A2/2 ? ----------
    P("=== (2) D_dyn(beta,a) = -d^2/da^2 log W_0 ; ratio to -(beta^2/2)A2 ===")
    P("%8s %6s %14s %18s %18s %10s" % ("theta", "a", "beta", "D_dyn", "-(b^2/2)A2", "ratio"))
    for th in [mp.mpf(0), mp.pi]:
        for a in [mp.mpf("0.5"), mp.mpf("1.5")]:
            A2 = mp.diff(lambda z: F2(z, th), a, 2)
            for beta in [mp.mpf("0.01"), mp.mpf("0.02"), mp.mpf("0.05")]:
                D = -mp.diff(lambda z: mp.log(W0(z, beta, th)), a, 2)
                pred = -(beta ** 2 / 2) * A2
                P("%8s %6s %14s %18s %18s %10s" % (mp.nstr(th, 6), mp.nstr(a, 3), mp.nstr(beta, 4),
                                                    mp.nstr(D, 10), mp.nstr(pred, 10),
                                                    mp.nstr(D / pred, 8)))
    P("")

    # ---------- (3) G'' - D_dyn = W_full : O(beta^0)? ----------
    P("=== (3) G''(beta,a) and W_full = G'' - D_dyn  (should -> W_geo(0) = O(1)) ===")
    P("%8s %6s %10s %16s %16s %16s" % ("theta", "a", "beta", "G''", "D_dyn", "W_full"))
    for th in [mp.mpf(0), mp.pi]:
        for a in [mp.mpf("0.5"), mp.mpf("1.5")]:
            for beta in [mp.mpf("0.02"), mp.mpf("0.1"), mp.mpf("0.5")]:
                Gpp = -mp.diff(lambda z: mp.log(Whaar(z, beta, th)), a, 2)
                D = -mp.diff(lambda z: mp.log(W0(z, beta, th)), a, 2)
                P("%8s %6s %10s %16s %16s %16s" % (mp.nstr(th, 6), mp.nstr(a, 3), mp.nstr(beta, 4),
                                                    mp.nstr(Gpp, 10), mp.nstr(D, 10),
                                                    mp.nstr(Gpp - D, 10)))
    P("")

    # ---------- (4) window numbers ----------
    P("=== (4) numbers in the strong-coupling window  beta < beta1 ===")
    beta1 = mp.mpf(1) / 48
    u1 = beta1 / (2 * N ** 2)
    P("beta1 = 1/48 =", mp.nstr(beta1, 8), "   u(beta1)=beta/(2N^2) =", mp.nstr(u1, 8))
    for th in [mp.mpf(0), mp.pi]:
        for a in [mp.mpf("0.0"), mp.mpf("1.0"), mp.mpf("2.0")]:
            A2 = mp.diff(lambda z: F2(z, th), a, 2)
            D = -(beta1 ** 2 / 2) * A2
            P("theta=%s a=%s : A2=%s  D_dyn(1/48)=%s  D_dyn(a0=1/16 as beta? no) " %
              (mp.nstr(th, 5), mp.nstr(a, 3), mp.nstr(A2, 10), mp.nstr(D, 10)))
    P("")
    P("conversion: beta = 2N^2 u = u/(2 a0)  =>  beta^2 = u^2/(4 a0^2)")
    P("            D_dyn = -(A2/8) (u/a0)^2 * [4D plaquette multiplicity 6 for a single-link mode]")
    P("")

    # ---------- (5) exact check of E[H] and Var(G) separately (both O(beta^2)) ----------
    P("=== (5) separate leading orders: E[H] and Var(G) each O(beta^2) ===")
    P("%6s %8s %16s %16s %16s" % ("a", "beta", "E[H]", "Var(G)", "E[H]-Var(G)"))
    th = mp.pi
    for a in [mp.mpf("1.5")]:
        for beta in [mp.mpf("0.01"), mp.mpf("0.02"), mp.mpf("0.04")]:
            def E(delta):
                # E[X] under e^{-beta(1-f)} J ... exact: derivative w.r.t. a of log W
                return None
            # d log W / da  =  < df/da >  (beta-avg); use analytic-free finite diff of the exact integrals
            def avg(Xfun):
                num = lambda rho: rho * J(mp.sqrt(a * a + rho * rho)) * mp.e ** (
                    -beta * (1 - f_slice(a, rho, th))) * Xfun(rho)
                den = lambda rho: rho * J(mp.sqrt(a * a + rho * rho)) * mp.e ** (
                    -beta * (1 - f_slice(a, rho, th)))
                return 2 * mp.pi * mp.quad(num, [0, rmax(a)]) / (2 * mp.pi * mp.quad(den, [0, rmax(a)]))
            # f_a = df/da, f_aa = d^2f/da^2 on the slice
            fa = lambda rho: mp.diff(lambda z: f_slice(z, rho, th), a, 1)
            faa = lambda rho: mp.diff(lambda z: f_slice(z, rho, th), a, 2)
            EH = -beta * avg(faa)
            VarG = beta ** 2 * (avg(lambda rho: fa(rho) ** 2) - avg(fa) ** 2)
            P("%6s %8s %16s %16s %16s" % (mp.nstr(a, 3), mp.nstr(beta, 4),
                                          mp.nstr(EH, 8), mp.nstr(VarG, 8),
                                          mp.nstr(EH - VarG, 8)))

    with open("research/v3/w7/n07_dyn_strong.log", "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
