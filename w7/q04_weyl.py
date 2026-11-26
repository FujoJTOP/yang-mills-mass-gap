"""Q04-WEYL: Weyl integration formula -> explicit 1-D integrals for Var(G_v|a), E[H_vv|a].

Model = single plaquette SU(2), conventions verbatim from H1A-EXACT.md 2.1 / O09-PHICURV 2.2:
    U = exp(xi),  xi in R^3,  rho=|xi| < 2pi,  u = xi_3/rho
    J = sin^2(rho/2)/rho^2,  Phi = -log J = -2 log(sin(rho/2)/rho)
    S_W = beta*(1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2) u)
    v = e_3, conditioning slice xi_3 = a  (u = a/rho)
Exact statements derived here (see Q04-WEYL.md):
    dU = (1/4pi^2) J d^3xi   (Haar),   int_{SU(2)} dU = 1
    a=0 slice:  G(rho) = beta*t*sin(rho/2)/rho,  H(rho) = (beta*c/2)*sin(rho/2)/rho
                => H = (1/2)cot(th/2) * G   pointwise in rho
    any a: sliced measure is 1-D in rho in [|a|,2pi]:
                dmu_a  prop  exp(-beta*S(rho)) * sin^2(rho/2)/rho * drho
Nothing is committed to git.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

TWOPI = 2.0 * np.pi
HALF = np.pi
_GL = {}


def gl(n):
    """cached Gauss-Legendre nodes (leggauss is O(n^2); never call it inside a loop)."""
    if n not in _GL:
        _GL[n] = leggauss(n)
    return _GL[n]


# ---------------------------------------------------------------- Weyl machinery
def cin(x, n=80):
    """Cin(x) = int_0^x (1-cos t)/t dt = sum_{k>=1} (-1)^{k+1} x^{2k}/((2k)! (2k))."""
    s = 0.0
    x2 = x * x
    fac = 1.0          # (2k)!
    pw = 1.0           # x^{2k}
    for k in range(1, n):
        fac *= (2 * k - 1) * (2 * k)
        pw *= x2
        s += (-1.0) ** (k + 1) * pw / (fac * 2 * k)
        if abs(pw / (fac * 2 * k)) < 1e-20:
            break
    return s


def weyl_su2_class(f, n=3000):
    """int_{SU(2)} f dU for class functions, f = f(rho): (1/pi) int_0^{2pi} sin^2(rho/2) f drho."""
    x, wq = gl(n)
    rho = 0.5 * TWOPI * (x + 1.0)
    return np.sum(wq * np.sin(rho / 2) ** 2 * f(rho)) * (0.5 * TWOPI) / np.pi


def weyl_su2_full(f2, n=700):
    """int_{SU(2)} f dU for f = f(rho,n_axis): (1/4pi^2) int drho sin^2(rho/2) int_{S^2} dOmega f."""
    x, wq = gl(n)
    rho = 0.5 * TWOPI * (x + 1.0)
    y, wy = gl(n)
    cz = y                                            # cos of polar angle of the axis
    om = np.sum(wy[None, :] * f2(rho[:, None], cz[None, :]), axis=1) * TWOPI   # int dphi = 2pi
    return np.sum(wq * np.sin(rho / 2) ** 2 * om) * (0.5 * TWOPI) / (4 * np.pi ** 2)


def weyl_sun_torus_check(n=3000):
    """Check the 1/N! and the torus normalisation for N=2:
       (1/2!) (1/2pi) int_0^{2pi} |e^{it}-e^{-it}|^2 dt = 1  (|Delta|^2 = 4 sin^2 t)."""
    x, wq = gl(n)
    t = np.pi * (x + 1.0)
    return np.sum(wq * 4.0 * np.sin(t) ** 2) * np.pi / (2.0 * 2 * np.pi)


# ---------------------------------------------------------------- model functions
def S_W(rho, a, th, beta):
    c, s = np.cos(th / 2), np.sin(th / 2)
    return beta * (1.0 - c * np.cos(rho / 2) + s * np.sin(rho / 2) * (a / rho))


def _w(rho):
    return np.sin(rho / 2) / rho


def _wp(rho):
    return ((rho / 2) * np.cos(rho / 2) - np.sin(rho / 2)) / rho ** 2


def _wpp(rho):
    return -np.sin(rho / 2) / (4 * rho) - np.cos(rho / 2) / rho ** 2 + 2 * np.sin(rho / 2) / rho ** 3


def G_v(rho, a, th, beta):
    """G = d_a S_W at fixed transverse eta (v = e_3), u = a/rho."""
    c, s = np.cos(th / 2), np.sin(th / 2)
    return beta * ((c / 2) * (a / rho) * np.sin(rho / 2) + s * (_w(rho) + (a * a / rho) * _wp(rho)))


def H_vv(rho, a, th, beta):
    """H = d_a^2 S_W (Hess S_W along v)."""
    c, s = np.cos(th / 2), np.sin(th / 2)
    A = np.sin(rho / 2) / rho + (a * a / (2 * rho ** 2)) * np.cos(rho / 2) - (a * a / rho ** 3) * np.sin(rho / 2)
    B = (3 * a / rho) * _wp(rho) - (a ** 3 / rho ** 3) * _wp(rho) + (a ** 3 / rho ** 2) * _wpp(rho)
    return beta * ((c / 2) * A + s * B)


def Phi1(rho, a):
    """Phi' = d_a Phi = Phi'(rho) * a/rho ;  Phi'(rho) = -cot(rho/2) + 2/rho."""
    return (-1.0 / np.tan(rho / 2) + 2.0 / rho) * (a / rho)


def Phi2(rho, a=0.0):
    """slice-second-derivative d_a^2 Phi at fixed transverse part:
         d_a^2 Phi = Phi_{rho rho} (a/rho)^2 + Phi_rho (1/rho)(1-(a/rho)^2)."""
    a = np.asarray(a, dtype=float)
    fr = a / rho
    return (_Phi2r(rho) * fr ** 2 + _Phi1r(rho) * (1.0 - fr ** 2) / rho)


def _Phi1r(rho):
    return -1.0 / np.tan(rho / 2) + 2.0 / rho


def _Phi2r(rho):
    return 0.5 / np.sin(rho / 2) ** 2 - 2.0 / rho ** 2


# ---------------------------------------------------------------- 1-D conditioned integral
def moments(a, th, beta, n=4000):
    """Exact 1-D quadrature of the conditioned slice measure
         dmu_a prop exp(-beta S) sin^2(rho/2)/rho drho on rho in [|a|,2pi]
       written in the transverse radius eta = sqrt(rho^2-a^2) so no sqrt singularity."""
    a = float(abs(a))
    eta_max = np.sqrt(TWOPI ** 2 - a * a)
    x, wq = gl(n)
    eta = 0.5 * eta_max * (x + 1.0)
    rho = np.sqrt(a * a + eta * eta)
    # 2-D plane measure d^2eta = 2pi eta deta, times J = sin^2(rho/2)/rho^2
    jac = 2 * np.pi * eta * np.sin(rho / 2) ** 2 / rho ** 2
    wt = jac * np.exp(-S_W(rho, a, th, beta)) * 0.5 * eta_max * wq
    Z = wt.sum()
    G = G_v(rho, a, th, beta)
    H = H_vv(rho, a, th, beta)
    P1 = Phi1(rho, a)
    P2 = Phi2(rho, a)
    EG = (wt * G).sum() / Z
    VG = (wt * G * G).sum() / Z - EG ** 2
    EH = (wt * H).sum() / Z
    EP1 = (wt * P1).sum() / Z
    VP1 = (wt * P1 * P1).sum() / Z - EP1 ** 2
    EP2 = (wt * P2).sum() / Z
    Cov = (wt * G * P1).sum() / Z - EG * EP1
    VGP1 = VG + VP1 + 2 * Cov
    Gpp = EH + EP2 - VGP1
    wr = _w(rho)
    Ew = (wt * wr).sum() / Z
    Varw = (wt * wr * wr).sum() / Z - Ew ** 2
    return dict(a=a, th=th, beta=beta, EG=EG, VG=VG, EH=EH, EP1=EP1, VP1=VP1,
                EP2=EP2, Cov=Cov, VGP1=VGP1, Gpp=Gpp, Z=Z, EGpP1=EG + EP1,
                Ew=Ew, Varw=Varw)


def slice_haar_moments(n=4000):
    """a=0, c=0 (th=pi): slice measure = sin^2(rho/2)/rho drho / (Cin(2pi)/2): the 'Haar' slice."""
    x, wq = gl(n)
    rho = 0.5 * TWOPI * (x + 1.0)
    w = np.sin(rho / 2) ** 2 / rho
    wt = 0.5 * TWOPI * wq
    Z = (wt * w).sum()
    f = _w(rho)
    Ef = (wt * w * f).sum() / Z
    E2 = (wt * w * f * f).sum() / Z
    EP2 = (wt * w * Phi2(rho, 0.0)).sum() / Z
    return Z, Ef, E2 - Ef ** 2, EP2


def moments_naive_weight(a, th, beta, n=4000):
    """Same 1-D reduction but with the (wrong) 3-D radial weight sin^2(rho/2) - for the discipline check."""
    a = float(abs(a))
    x, wq = gl(n)
    rho = a + (TWOPI - a) * 0.5 * (x + 1.0)
    wt = (TWOPI - a) * 0.5 * wq * np.sin(rho / 2) ** 2 * np.exp(-S_W(rho, a, th, beta))
    Z = wt.sum()
    G, H = G_v(rho, a, th, beta), H_vv(rho, a, th, beta)
    EG = (wt * G).sum() / Z
    return EG, (wt * G * G).sum() / Z - EG ** 2, (wt * H).sum() / Z


def plane_2d_check(a, th, beta, L=None, n=900):
    """Brute-force 2-D Cartesian quadrature on the slice {xi_3=a}: independent of the 1-D reduction."""
    a = float(abs(a))
    eta_max = np.sqrt(TWOPI ** 2 - a * a)
    L = eta_max if L is None else L
    x, wq = gl(n)
    e1 = L * x
    w1 = L * wq
    E1, E2 = np.meshgrid(e1, e1, indexing="ij")
    W1, W2 = np.meshgrid(w1, w1, indexing="ij")
    rho = np.sqrt(a * a + E1 ** 2 + E2 ** 2)
    inside = rho < TWOPI
    rho = np.where(inside, rho, 1.0)
    jac = np.where(inside, np.sin(rho / 2) ** 2 / rho ** 2, 0.0)
    wt = W1 * W2 * jac * np.exp(-S_W(rho, a, th, beta))
    Z = wt.sum()
    G = G_v(rho, a, th, beta)
    H = H_vv(rho, a, th, beta)
    EG = (wt * G).sum() / Z
    return EG, (wt * G * G).sum() / Z - EG ** 2, (wt * H).sum() / Z


def fd_check(a, th, beta, eta=0.37, hs=(1e-3, 5e-4)):
    """Finite-difference verification of G_v and H_vv at fixed transverse part."""
    out = []
    for h in hs:
        def Sfun(aa):
            rho = np.hypot(aa, eta)
            return S_W(rho, aa, th, beta) if False else beta * (1 - np.cos(th / 2) * np.cos(rho / 2)
                                                               + np.sin(th / 2) * np.sin(rho / 2) * (aa / rho))
        g = (Sfun(a + h) - Sfun(a - h)) / (2 * h)
        hh = (Sfun(a + h) - 2 * Sfun(a) + Sfun(a - h)) / h ** 2
        out.append((g, hh))
    rho = np.hypot(a, eta)
    return G_v(rho, a, th, beta), H_vv(rho, a, th, beta), out


def main():
    L = []
    P = L.append
    P("Q04-WEYL  single plaquette SU(2):  Weyl formula -> 1-D integrals for Var(G_v|a), E[H_vv|a]")
    P("=" * 100)

    # ---------------- [1] Weyl normalisation and the J constant
    c2p = cin(TWOPI)
    P("[1] Weyl / Haar constants, SU(2)")
    P("    Cin(2pi)                         = %.10f" % c2p)
    P("    1/(2 Cin(2pi))                   = %.10f   (paper W0 = 0.20511530)" % (1.0 / (2 * c2p)))
    P("    int_{SU(2)} 1 dU  (class formula (1/pi) int_0^2pi sin^2(rho/2) drho) = %.15f"
      % weyl_su2_class(lambda r: np.ones_like(r)))
    P("    int_{SU(2)} 1 dU  (full     (1/4pi^2) int drho sin^2(rho/2) int_S2 dOmega) = %.15f"
      % weyl_su2_full(lambda r, cz: np.ones_like(r * cz)))
    P("    int_T |Delta|^2 dt /(2pi) /2!    = %.15f   (1/N! needed; N=2: |Delta|^2=4sin^2 t)"
      % weyl_sun_torus_check())
    x, wq = gl(4000)
    rho = 0.5 * TWOPI * (x + 1.0)
    P("    int J d^3xi = int_0^2pi sin^2(rho/2)drho * 4pi = %.12f  (= 4 pi^2 = %.12f)"
      % (np.sum(wq * np.sin(rho / 2) ** 2) * 0.5 * TWOPI * 4 * np.pi, 4 * np.pi ** 2))
    P("    => dU = (1/4pi^2) J d^3xi exactly:  the plaquette model IS a Haar integral.")
    P("")

    # ---------------- [2] slice-Haar moments: reproduce H1A anchors
    Z, Ew, Varw, EP2 = slice_haar_moments()
    P("[2] a=0, theta=pi  -> slice measure = 'Haar' slice  sin^2(rho/2)/rho drho")
    P("    Z = Cin(2pi)/2                   = %.12f   vs %.12f" % (Z, 0.5 * c2p))
    P("    E[w]  , w = sin(rho/2)/rho       = %.10f" % Ew)
    P("    V0 = Var[w]                      = %.10f   (H1A archive 0.0110593)" % Varw)
    P("    E[Phi'']                         = %.10f   (O09/H1A archive 0.2051153)  = 1/(2Cin2pi)? %s"
      % (EP2, abs(EP2 - 1.0 / (2 * c2p)) < 1e-9))
    P("    beta^2 V0 at beta=2.5            = %.10f   (H1A archive Var=0.069121)" % (6.25 * Varw))
    Eg_n, Vg_n, Eh_n = moments_naive_weight(0.0, np.pi, 2.5)
    P("    [discipline] same quadrature with the WRONG 3-D radial weight sin^2(rho/2):")
    P("        E[G]=%.6f  Var(G)=%.6f  E[H]=%.6f   -> does NOT reproduce 0.0691210" % (Eg_n, Vg_n, Eh_n))
    P("")

    # ---------------- [3] reproduce H1A theta=0 (independent of any slicing subtlety)
    P("[3] theta=0 check: Var(G)=0 exactly (t=0),  E[H]=(beta/2)E[w]  vs H1A archive")
    for beta, ref in ((0.5, 0.095581), (2.5, 0.553061), (6.0, 1.420010), (10.0, 2.418516)):
        m = moments(0.0, 0.0, beta)
        P("    beta=%6.3f   E[H]=%.6f  (archive %.6f)  Var(G)=%.3e  E[Phi'']=%.7f"
          % (beta, m["EH"], ref, m["VG"], m["EP2"]))
    P("")

    # ---------------- [4] pointwise identity on the a=0 slice
    P("[4] a=0 slice, pointwise in rho:  H(rho)/G(rho) = (1/2)cot(theta/2)  (exact, no inequality)")
    for th in (0.25, 0.5 * np.pi, 0.75 * np.pi, 0.9 * np.pi):
        rr = np.linspace(0.05, TWOPI - 0.05, 777)
        r = H_vv(rr, 0.0, th, 1.0) / G_v(rr, 0.0, th, 1.0)
        P("    theta=%.4f pi   ratio = %.12f   target 0.5*cot(th/2) = %.12f   max dev = %.2e"
          % (th / np.pi, r[0], 0.5 / np.tan(th / 2), np.max(np.abs(r - 0.5 / np.tan(th / 2)))))
    for th in (0.0, np.pi):
        P("    theta=%.1f pi: degenerate (t=0 -> G==0 / c=0 -> H==0)" % (th / np.pi))
    P("")

    # ---------------- [5] reproduce O09 table (beta=2.5, theta=pi, general a): validates the reduction
    P("[5] O09-PHICURV 2.2 table reproduced from the 1-D integral (beta=2.5, theta=pi)")
    ref = {
        0.0: (0.2051, 0.0000, 0.0691, 0.0000, 0.0000, +0.1360),
        0.25: (0.2094, -0.0632, 0.0714, 2.9e-4, -0.0034, +0.0812),
        0.5: (0.2162, -0.1240, 0.0716, 1.4e-3, -0.0072, +0.0336),
        0.75: (0.2262, -0.1816, 0.0697, 3.6e-3, -0.0114, -0.0059),
        1.0: (0.2398, -0.2350, 0.0657, 7.5e-3, -0.0158, -0.0368),
        1.5: (0.2818, -0.3277, 0.0527, 2.30e-2, -0.0243, -0.0730),
        2.0: (0.3516, -0.3990, 0.0361, 5.50e-2, -0.0306, -0.0773),
    }
    P("    a      E[Phi'']  (ref)      E[H]     (ref)      Var(G)   (ref)      Cov     (ref)      G''     (ref)")
    for a, rf in ref.items():
        m = moments(a, np.pi, 2.5)
        P("   %.2f  %8.4f (%7.4f)  %8.4f (%8.4f)  %7.4f (%7.4f)  %8.4f (%8.4f)  %8.4f (%8.4f)"
          % (a, m["EP2"], rf[0], m["EH"], rf[1], m["VG"], rf[2], m["Cov"], rf[4], m["Gpp"], rf[5]))
    P("")

    # ---------------- [6] independent checks of the 1-D reduction
    P("[6] independent checks")
    for a, th in ((0.0, np.pi), (0.5, np.pi), (0.5, 0.5 * np.pi)):
        m1 = moments(a, th, 2.5)
        m2 = plane_2d_check(a, th, 2.5)
        P("    a=%.2f th=%.4fpi : 1-D  E[G]=%+.8f Var(G)=%.8f E[H]=%+.8f" % (a, th / np.pi, m1["EG"], m1["VG"], m1["EH"]))
        P("                      2-D  E[G]=%+.8f Var(G)=%.8f E[H]=%+.8f   max|diff|=%.2e"
          % (m2[0], m2[1], m2[2], max(abs(m1["EG"] - m2[0]), abs(m1["VG"] - m2[1]), abs(m1["EH"] - m2[2]))))
    for a, th in ((0.0, np.pi), (0.4, 0.0), (0.9, 0.6 * np.pi)):
        G0, H0, fd = fd_check(a, th, 2.5)
        P("    a=%.2f th=%.4fpi : G(analytic)=%+.10f  G(FD)=%+.10f (h=1e-3)" % (a, th / np.pi, G0, fd[0][0]))
        P("                    : H(analytic)=%+.10f  H(FD)=%+.10f   rel dev vs h=5e-4: %.2e"
          % (H0, fd[0][1], abs(H0 - fd[1][1]) / max(1e-12, abs(fd[1][1]))))
    P("")

    # ---------------- [7] the two 1-D integrals, and exact values
    P("[7] single plaquette: exact 1-D integrals")
    P("    a=0:  Var(G|0) = beta^2 t^2 Var_{mu_th}(w),   E[H|0] = (beta c/2) E_{mu_th}(w)")
    P("          mu_th prop exp(-beta(1-c cos(rho/2))) sin^2(rho/2)/rho drho,  w=sin(rho/2)/rho")
    P("    a!=0: 1-D in rho in [|a|,2pi], same weight with S=beta(1-c cos(rho/2)+t sin(rho/2) a/rho)")
    P("")
    P("    beta=2.5,  theta dependence at a=0 (t=sin(th/2), c=cos(th/2)):")
    P("      theta      c        E[w]      Var(w)     Var(G)     E[H]     Var/E[H]   E[H]-Var")
    for th in (0.0, 0.1 * np.pi, 0.25 * np.pi, 0.4 * np.pi, 0.5 * np.pi, 0.6 * np.pi,
               0.75 * np.pi, 0.9 * np.pi, 0.99 * np.pi, np.pi):
        m = moments(0.0, th, 2.5)
        rt = m["VG"] / m["EH"] if abs(m["EH"]) > 1e-14 else np.inf
        P("    %.3fpi  %+.5f  %8.5f  %.6f  %9.6f  %+9.6f  %9.4f  %+9.6f"
          % (th / np.pi, np.cos(th / 2), m["Ew"], m["Varw"], m["VG"], m["EH"], rt, m["EH"] - m["VG"]))
    P("")
    P("    beta dependence at theta=pi and theta=pi/2 (a=0):")
    P("      beta   theta      Var(G)      E[H]     Var/E[H]")
    for beta in (0.5, 1.0, 2.5, 4.30660785, 6.0, 10.0):
        for th in (np.pi, 0.5 * np.pi):
            m = moments(0.0, th, beta)
            rt = m["VG"] / m["EH"] if abs(m["EH"]) > 1e-13 else np.inf
            P("     %6.3f   %.2fpi  %10.6f  %+10.6f   %9.4f" % (beta, th / np.pi, m["VG"], m["EH"], rt))
    P("")
    P("    a dependence at beta=2.5, theta=pi (NOTE: E[H]<0 for a>0 => Var<=E[H] fails there):")
    P("      a       Var(G)      E[H]     Var/E[H]   E[H]-Var")
    for a in (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.8, 3.35, 3.9, 4.45):
        m = moments(a, np.pi, 2.5)
        rt = m["VG"] / m["EH"] if abs(m["EH"]) > 1e-13 else np.inf
        P("     %.2f  %10.6f  %+10.6f   %9.4f  %+10.6f" % (a, m["VG"], m["EH"], rt, m["EH"] - m["VG"]))
    P("")

    # ---------------- [8] can Var(G)<=E[H] be seen algebraically?
    P("[8] algebraic verdict")
    P("    On a=0: H = (1/2)cot(th/2) G pointwise, so")
    P("        Var(G) <= E[H]  <=>  2 beta t^2 Var(w) <= c E[w]        (one inequality, one measure)")
    P("    Pointwise bound 0 <= w = sin(rho/2)/rho <= 1/2  =>  Var(w) <= E[w^2] <= (1/2)E[w]")
    P("        => Var(G) <= (beta^2 t^2/2) E[w] <= (beta c/2)E[w] = E[H]   whenever  beta t^2 <= c")
    for beta in (0.5, 1.0, 2.5):
        th_ok = np.linspace(0.0, np.pi, 20001)
        ok = beta * np.sin(th_ok / 2) ** 2 <= np.cos(th_ok / 2)
        th_b = th_ok[ok].max() if ok.any() else 0.0
        P("      beta=%.2f: sufficient condition beta sin^2(th/2)<=cos(th/2) holds for th/pi <= %.4f" % (beta, th_b / np.pi))
    # exact crossing theta_c(beta): Var=E[H]
    P("    exact crossing theta_c(beta) where Var(G)=E[H] (a=0):")
    for beta in (0.5, 1.0, 2.5, 6.0):
        ths = np.linspace(0.02, np.pi - 1e-4, 400)
        vals = np.array([moments(0.0, t, beta, n=1200) for t in ths])
        d = np.array([v["EH"] - v["VG"] for v in vals])
        sgn = np.where(d >= 0)[0]
        i = sgn[-1] if len(sgn) else 0
        if i + 1 < len(ths) and i > 0:
            # linear interpolation in th
            t0, t1 = ths[i], ths[i + 1]
            d0, d1 = d[i], d[i + 1]
            tc = t0 - d0 * (t1 - t0) / (d1 - d0)
        else:
            tc = ths[i]
        P("      beta=%.2f : theta_c/pi = %.4f  (Var<=E[H] for th < theta_c)" % (beta, tc / np.pi))
    P("")

    # ---------------- [9] H1-form check (with W_full), single plaquette
    P("[9] the H1-form inequality Var(G) <= E[H] + W_geo holds iff G''>=0 (single plaquette, a=0)")
    P("      beta    theta   Var(G)    E[H]    E[Phi'']   E[H]+E[Phi'']-Var    G''")
    for beta in (2.5, 4.30660785, 6.0):
        for th in (np.pi, 0.5 * np.pi):
            m = moments(0.0, th, beta)
            P("     %6.3f   %.2fpi  %8.5f  %+8.5f   %8.5f   %+12.6f   %+10.6f"
              % (beta, th / np.pi, m["VG"], m["EH"], m["EP2"], m["EH"] + m["EP2"] - m["VG"], m["Gpp"]))
    P("")
    P("[10] L=2 lattice reference (different object, P02): Var(G)=5.4477, E[H]=6.3315, ratio=%.4f"
      % (5.4477 / 6.3315))
    P("     single plaquette does NOT reproduce these; it is a 4-D 96-plaquette constrained HMC value.")

    txt = "\n".join(L)
    print(txt)
    with open("w7/q04_weyl.log", "w", encoding="utf-8") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()






