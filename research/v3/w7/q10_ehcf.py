"""Q10 — E[H|0] vs Var(G|0) on ONE SU(2) plaquette: closed forms, term by term, and the
beta^2 (theta=pi slice) vs beta^1.57 (L=2, p19) arbitration.

Exact reduction (same conventions as w7/p05_cinf.py and w7/o05_var_scaling.py):
  xi in R^3, |xi|<2pi;  S = beta(1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2)u), u = xi3/rho;
  flat chart dU = J dxi, J = sin^2(rho/2)/rho^2, Phi = -log J;  v = e3;  slice {xi3 = 0} (a0=0).
  On the slice, x := rho/2 in (0,pi):
      w(x) dx  ∝  sin^2 x / x * e^{B cos x} dx ,        B := beta*cos(th/2)      [slice measure]
      G  = d_v S   = beta sin(th/2) sin x / x  =: beta*s*g(x),   g(x) = sin x/(2x)
      H  = d_v^2 S = beta cos(th/2) g(x)       =  beta*c*g(x)
      Phi'' = 2/rho^2 - cot(rho/2)/rho          (slice value; E[Phi''] = W0 at B=0)
  =>  E[H|0]   = (beta c/2) <A0>(B) ,      <A0>(B) := E_{w_B}[g]
      Var(G|0) = beta^2 s^2 V0(B) ,        V0(B)   := Var_{w_B}(g)
      V0(0)    = (N2 Z - N1^2)/(4 Z^2)   [h1e closed form]
      Z = Cin(2pi)/2 ,  N1 = (3/4)[ln3 + Ci(pi) - Ci(3pi)] ,  N2 = ln2 + Ci(2pi) - Ci(4pi)
      V0(B) ~ 1/(36 B^2) as B->inf       [Laplace: y=x^2 ~ Exp(B/2), g = 1/2 - y/12 + ...]

Only this script and its .log are written.  No git commit.
"""
import math

import numpy as np
from scipy.special import sici

PI = math.pi
GAM = float(np.euler_gamma)
LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


def Ci(x):
    return float(sici(x)[1])


def Cin(x):
    return GAM + math.log(x) - Ci(x)


# ------------------------------------------------------------------ quadrature on x in (0,pi)
NX = 3000
_XX, _WW = np.polynomial.legendre.leggauss(NX)
_X = 0.5 * PI * (_XX + 1.0)                      # (0, pi)
_W = 0.5 * PI * _WW
_LOGBASE = 2.0 * np.log(np.sin(_X)) - np.log(_X)  # log(sin^2 x / x)


def wmom(B, f):
    """E_w[f] under w ∝ sin^2 x / x * exp(B cos x); B may be any real (log-sum-exp)."""
    L = _LOGBASE + B * np.cos(_X)
    L -= L.max()
    w = _WW * np.exp(L)
    return float(np.sum(w * f) / np.sum(w))


def gx():
    return np.sin(_X) / (2.0 * _X)


_G = gx()


def A0(B):
    return wmom(B, _G)


def V0(B):
    m = A0(B)
    return wmom(B, _G * _G) - m * m


# ------------------------------------------------------------------ 0. closed-form constants
ZC = 0.5 * Cin(2.0 * PI)
N1C = 0.75 * (math.log(3.0) + Ci(PI) - Ci(3.0 * PI))
N2C = math.log(2.0) + Ci(2.0 * PI) - Ci(4.0 * PI)
V0_ZERO = (N2C * ZC - N1C ** 2) / (4.0 * ZC ** 2)
A0_ZERO = N1C / (2.0 * ZC)
W0 = 1.0 / (2.0 * Cin(2.0 * PI))                  # = E[Phi''|0] at theta=pi (beta-free)
KAPPA = 1.0 / 36.0
BETA_STAR = math.sqrt(W0 / V0_ZERO)

P("=" * 100)
P("q10_ehcf.log — Q10: Var(G|0) closed form vs E[H|0]; beta^2 (theta=pi) vs beta^1.57 (L=2, p19)")
P("=" * 100)
P("  numpy %s / scipy sici OK ; Gauss-Legendre n=%d on x in (0,pi) ; Legendre nodes finite: %s"
  % (np.__version__, NX, bool(np.isfinite(_X).all() and np.isfinite(_W).all())))
P("")

P("[0] h1e closed-form endpoint V0(0) and the two theta=pi anchors")
P("    Ci(pi)=%.12f Ci(2pi)=%.12f Ci(3pi)=%.12f Ci(4pi)=%.12f"
  % (Ci(PI), Ci(2 * PI), Ci(3 * PI), Ci(4 * PI)))
P("    Z=Cin(2pi)/2 = %.15f" % ZC)
P("    N1 = (3/4)[ln3+Ci(pi)-Ci(3pi)] = %.15f" % N1C)
P("    N2 = ln2+Ci(2pi)-Ci(4pi)  = %.15f" % N2C)
P("    V0(0) = (N2 Z - N1^2)/(4 Z^2) = %.15f      [h1e: 0.011059293772973]" % V0_ZERO)
P("    <A0>(0) = N1/(2Z)             = %.15f      [h1e: 0.357411353391]" % A0_ZERO)
P("    W0 = E[Phi''|0](theta=pi) = 1/(2Cin(2pi)) = %.15f   [h1e: 0.205115297123]" % W0)
P("    beta* = sqrt(W0/V0) = %.8f   [h1e: 4.30660785; p19 fit: zero of Gpp at 4+ ... see Sec.3]"
  % BETA_STAR)
P("    quadrature check of the same endpoint: V0(0)=%.15f  <A0>(0)=%.15f"
  % (V0(1e-12), A0(1e-12)))
P("")

P("[1] B -> inf asymptote: the ONLY other exactly known point of the same closed form")
P("    %8s %14s %14s %10s %10s" % ("B", "V0(B)", "B^2*V0(B)", "<A0>(B)", "(1/2-<A0>)B"))
for B in (1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 1e4):
    v = V0(B)
    P("    %8.0f %14.9f %14.9f %10.6f %10.6f" % (B, v, B * B * v, A0(B), (0.5 - A0(B)) * B))
P("    1/36 = %.12f   (p05 [P2]: (1/2-<A0>)B -> 1/6 = %.6f)" % (KAPPA, 1.0 / 6.0))
P("    => V0(B) = 1/(36 B^2) * (1 + o(1));  Var(G) -> sin^2(th/2)/(36 cos^2(th/2)) = tan^2(th/2)/36")
for th in (PI / 2, 3 * PI / 4):
    P("       theta/pi=%.4f : tan^2(th/2)/36 = %.6f   [o05 [P4] sup_beta Var = %.6f]"
      % (th / PI, math.tan(th / 2) ** 2 / 36.0, {PI / 2: 0.027773, 3 * PI / 4: 0.161844}[th]))
P("")

P("[2] theta=pi: Var(G)=beta^2 V0 exactly, E[H]=0 exactly -> the ratio is not beta^2-vs-beta")
P("    %8s %14s %10s %12s %14s %14s" % ("beta", "Var(G)", "E[H]", "E[H+Phi'']", "Var/(E[H+P])", "Var/E[H]"))
for b in (1.0, 2.5, 4.0, BETA_STAR, 6.0, 10.0, 100.0):
    v = b * b * V0_ZERO
    tot = W0
    P("    %8.4f %14.6f %10.1e %12.6f %14.6f %14s" % (b, v, 0.0, tot, v / tot, "inf"))
P("    => Var(G) crosses E[H+Phi''] = W0 (a BETA-FREE number) at beta* = sqrt(W0/V0) = %.6f,"
  % BETA_STAR)
P("       because on the theta=pi slice the measure is beta-free (S == beta identically) while")
P("       E[H|0] = (beta c/2)<A0> dies with c = cos(th/2) = 0.  The loser is Phi'', not H.")
P("")

P("[3] general theta: Var(G)=beta^2 s^2 V0(B), E[H]=(beta c/2)<A0>(B),  B=beta c")
P("    %7s %7s %12s %12s %10s %10s" % ("th/pi", "beta", "Var(G)", "E[H]", "Var/E[H]", "B"))
for t in (0.5, 0.9, 0.95, 0.99, 1.0):
    th = t * PI
    c, s = math.cos(th / 2), math.sin(th / 2)
    for b in (1.0, 2.5, 6.0, 20.0, 100.0, 1000.0):
        B = b * c
        v = b * b * s * s * V0(B)
        e = 0.5 * b * c * A0(B)
        rat = (v / e) if e > 1e-14 else float("inf")
        P("    %7.3f %7.1f %12.6f %12.6f %10.4f %10.3f" % (t, b, v, e, rat, B))
    P("")
P("")

P("[4] same-window exponent test (beta in [1,6], p19's window): dlog Var/dlog beta")
P("    %16s %10s %10s %10s" % ("theta/pi", "k(Var)", "k(E[H])", "k(E[H])-k(Var)"))
BG = np.array([1.0, 2.5, 4.0, 6.0])
for t in (0.5, 0.9, 0.99, 0.999, 1.0):
    th = t * PI
    c, s = math.cos(th / 2), math.sin(th / 2)
    vs = np.array([b * b * s * s * V0(b * c) for b in BG])
    es = np.array([0.5 * b * c * A0(b * c) for b in BG])
    kv = np.polyfit(np.log(BG), np.log(vs), 1)[0]
    ke = np.polyfit(np.log(BG), np.log(es), 1)[0]
    P("    %16.4f %10.4f %10.4f %10.4f" % (t, kv, ke, ke - kv))
P("    p19 (L=2, same beta window, research/unconditional/w7/p19_fit.log l.19):")
P("      k(Var F) = 1.57 ,  k(E[kappa]) = 1.25 ,  k(Gpp) = 0.41 ,  Var/E[kappa] -> 1 from below")
P("    => the beta^2 is the theta=pi ENDPOINT of this family; one step off the slice (th=0.999pi)")
P("       the same window already gives k < 2.  Neither 2 nor 1.57 is a universal exponent.")
P("")

P("[5] criterion: where does Var(G) > E[H|0] actually happen?  (grid over beta, theta)")
bt = np.exp(np.linspace(math.log(0.1), math.log(1e4), 400))
tt = np.linspace(0.0, 1.0, 801) * PI
Cc = np.cos(tt / 2)[None, :]
Ss = np.sin(tt / 2)[None, :]
BB = bt[:, None] * Cc
_bg = np.exp(np.linspace(math.log(1e-8), math.log(1e5), 900))
_vb = np.array([V0(b) for b in _bg])
_ab = np.array([A0(b) for b in _bg])
_lb = np.log(_bg)
_sel = slice(None)                        # B in [1e-8, 1e5] covers the whole grid
Vg = np.interp(np.log(np.maximum(BB, 1e-8)), _lb[_sel], _vb[_sel])
Ag = np.interp(np.log(np.maximum(BB, 1e-8)), _lb[_sel], _ab[_sel])
Var = bt[:, None] ** 2 * Ss ** 2 * Vg
EH = 0.5 * bt[:, None] * Cc * Ag
R = np.where(EH > 1e-300, Var / np.maximum(EH, 1e-300), np.inf)
bad = R > 1.0
P("    violation set {R>1}: %d / %d grid points" % (int(bad.sum()), bad.size))
if bad.any():
    j0 = int(np.where(bad.any(0))[0][0])         # smallest theta that violates somewhere
    P("    theta_crit (first violating theta/pi) = %.4f   at beta = %.3g"
      % (tt[j0] / PI, bt[int(np.where(bad[:, j0])[0][0])]))
    P("    %10s %12s %12s" % ("theta/pi", "beta_lo(R>1)", "beta_hi(R>1)"))
    for t in (tt[j0] / PI, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999, 1.0):
        j = int(round(t * 800))
        idx = np.where(bad[:, j])[0]
        if idx.size:
            P("    %10.4f %12.4g %12.4g" % (t, bt[idx[0]], bt[idx[-1]]))
        else:
            P("    %10.4f %12s %12s" % (t, "-", "-"))
    P("    theta>pi-eps: beta_hi ~ %.3g/eps?  beta_hi(theta) grows as the slice degenerates:"
      % (1.0))
    for t in (0.99, 0.995, 0.999):
        j = int(round(t * 800))
        idx = np.where(bad[:, j])[0]
        eps = PI * (1.0 - t)
        if idx.size:
            P("       eps=%.4f  beta_hi=%10.4g   beta_hi*eps^2=%.4g   beta_hi*eps=%.4g"
              % (eps, bt[idx[-1]], bt[idx[-1]] * eps ** 2, bt[idx[-1]] * eps))
P("    ratio at the moderate slices (no violation):")
for t in (0.5, 0.9):
    j = int(round(t * 800))
    P("       theta/pi=%.3f : max_beta R = %.6f  at beta=%.3g"
      % (t, float(R[:, j].max()), float(bt[int(np.argmax(R[:, j]))])))
P("")

P("[6] verdict quantities")
P("    (i)   closed forms: Var(G|0) = beta^2 sin^2(th/2) V0(beta cos(th/2));  "
  "E[H|0] = (beta/2)cos(th/2)<A0>(beta cos(th/2))")
P("    (ii)  V0(0) = %.12f [T, h1e] ;  B^2 V0 -> 1/36 [C->T] ;  <A0>(0)=%.9f, <A0>(inf)=1/2"
  % (V0_ZERO, A0_ZERO))
P("    (iii) ratio Var/E[H] = 2 beta sin^2(th/2) V0(B) / (c <A0>(B)) ~ sin^2/(18 c^3 <A0> beta)  -> 0: LARGE beta WINS")
P("    (iv)  theta=pi exactly: Var=beta^2 V0, E[H]=0 for ALL beta -> inequality fails for every beta>0,")
P("          NOT at beta*: beta*=%.6f is the root of Var(G)=E[H+Phi'']=W0, the beta-free Jacobian term."
  % BETA_STAR)
P("    (v)   L=2 (p19, beta in [1,6]): Var(F)/E[kappa] = 0.505, 0.787, 0.848, 0.889 -> 1 from below,")
P("          1 - D/S = 0.483 beta^-0.834, Gpp(0) > 0 throughout (0.2835 -> 0.5913).  No crossing in window.")
P("    (vi)  exponent of Var in p19's own window: single plaquette gives k=%.3f (th=0.9pi) ... "
  "k=2.000 (th=pi)" % np.polyfit(np.log(BG), np.log([b * b * math.sin(0.45 * PI) ** 2 * V0(b * math.cos(0.45 * PI)) for b in BG]), 1)[0])
P("")
P("Clay 0%.  No QED.")
P("")

with open("q10_ehcf.log", "w", encoding="utf-8") as f:
    f.write("\n".join(LOG) + "\n")
