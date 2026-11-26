"""R01-HG -- does the q04 pointwise identity H = (1/2)cot(theta/2) G survive a != 0?

Model = H1A-EXACT 2.1 / Q04 (2.4)(2.5) / Q08 (verbatim):  rho=|xi|, xi_3=a, s=|eta|,
  c=cos(th/2), t=sin(th/2), w(rho)=sin(rho/2)/rho,
  S(a,rho) = beta(1 - c cos(rho/2) + t a w(rho)),  rho=sqrt(a^2+s^2)
  G(a,rho) = dS/da  (pointwise force),   H(a,rho) = d^2S/da^2  (pointwise Hessian)

Slice measure (Q04 (1.6)):  dmu_a ~ w(rho) exp(-beta S) drho on rho in [|a|,2pi].
Q08 representation:  Z(a)=e^{-beta} int_|a|^2pi w exp(-beta a phi) drho,  phi=sin(rho/2)/rho.

Sections
  [1] sympy: (2.4)(2.5) vs exact differentiation of S;  Wronskian  Wr = C*B - D*A
  [2] sympy: rho=2pi endpoint values,  H/G at 2pi,  Wr at 2pi
  [3] numeric: pointwise G,H on the rho grid;  spread of H/G   (beta=2.5, theta=pi/2)
  [4] numeric: moments E[H],E[G],Var,Cov; the four candidate constants; R^2
  [5] a=0 control: everything must collapse to K0=(1/2)cot(th/2)
"""
import sys
import numpy as np
import sympy as sp

LOG = open(r"D:\Dev\math\w7\r01_hgident.log", "w", encoding="utf-8")


def say(*x):
    line = " ".join(str(t) for t in x)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


TWOPI = 2.0 * np.pi
_GL = {}


def gl(n):
    if n not in _GL:
        x, wg = np.polynomial.legendre.leggauss(n)
        _GL[n] = (0.5 * (x + 1.0), 0.5 * wg)
    return _GL[n]


def integ(f, lo, hi, n):
    if hi <= lo:
        return 0.0
    tt, wt = gl(n)
    return (hi - lo) * float(np.dot(wt, f(lo + (hi - lo) * tt)))


# ---------------------------------------------------------------- pointwise fields
def wfun(r):
    r = np.asarray(r, dtype=float)
    out = np.empty_like(r)
    m = r != 0.0
    out[m] = np.sin(0.5 * r[m]) / r[m]
    out[~m] = 0.5                      # limit rho->0
    return out


def wpf(r):
    r = np.asarray(r, dtype=float)
    return (0.5 * r * np.cos(0.5 * r) - np.sin(0.5 * r)) / np.where(r == 0.0, 1.0, r) ** 2


def wppf(r):
    r = np.asarray(r, dtype=float)
    rr = np.where(r == 0.0, 1.0, r)
    return -np.sin(0.5 * r) / (4.0 * rr) - np.cos(0.5 * r) / rr ** 2 + 2.0 * np.sin(0.5 * r) / rr ** 3


def fields(a, rho, beta, th):
    """Q04 (2.4)(2.5)."""
    c, t = np.cos(0.5 * th), np.sin(0.5 * th)
    rho = np.asarray(rho, dtype=float)
    u = a / rho
    w, wp, wpp = wfun(rho), wpf(rho), wppf(rho)
    G = beta * ((c / 2.0) * u * np.sin(0.5 * rho) + t * (w + (a * a / rho) * wp))
    H = beta * ((c / 2.0) * (np.sin(0.5 * rho) / rho + (a * a / (2.0 * rho ** 2)) * np.cos(0.5 * rho)
                             - (a * a / rho ** 3) * np.sin(0.5 * rho))
                + t * ((3.0 * a / rho) * wp - (a ** 3 / rho ** 3) * wp + (a ** 3 / rho ** 2) * wpp))
    return G, H


def Sweight(a, rho, beta, th):
    c, t = np.cos(0.5 * th), np.sin(0.5 * th)
    return beta * (1.0 - c * np.cos(0.5 * rho) + t * a * wfun(rho))


# ================================================================ [1] sympy algebra
say("=" * 100)
say("R01-HG  a!=0 generalisation of  H = (1/2)cot(theta/2) G   [single-plaquette SU(2), H1A/Q04/Q08]")
say(f"numpy {np.__version__}   sympy {sp.__version__}   python {sys.version.split()[0]}")
say("=" * 100)
say("")
say("[1] sympy: G,H affine in (c/2, t)  ->  G=(c/2)A+tB, H=(c/2)C+tD ;  Wronskian Wr=C*B-D*A")

r, a, b = sp.symbols("rho a beta", positive=True)
S2, C2 = sp.symbols("S2 C2")                     # sin(rho/2), cos(rho/2)
w = S2 / r
wp = ((r / 2) * C2 - S2) / r ** 2
wpp = -S2 / (4 * r) - C2 / r ** 2 + 2 * S2 / r ** 3

A = b * (a / r) * S2
B = b * (w + a ** 2 * wp / r)
Cc = b * (S2 / r + a ** 2 * C2 / (2 * r ** 2) - a ** 2 * S2 / r ** 3)
Dd = b * ((3 * a / r) * wp - a ** 3 * wp / r ** 3 + a ** 3 * wpp / r ** 2)

Wr = sp.simplify(sp.expand(Cc * B - Dd * A))
say("  Wr(rho,a) = C*B - D*A  =")
say("    " + sp.sstr(sp.factor(sp.simplify(Wr))))
say("")
say("  sanity: A(0)=D(0)=0, B(0)=C(0)=beta*w  ->  rank collapse to span{w}:")
say("    A|a=0 =", sp.simplify(A.subs(a, 0)), "   D|a=0 =", sp.simplify(Dd.subs(a, 0)))
say("    B|a=0 =", sp.simplify(B.subs(a, 0)), "   C|a=0 =", sp.simplify(Cc.subs(a, 0)))

# exact differentiation check of (2.4)(2.5) (symbolic, via rho(a) with rho'=a/rho)
s_ = sp.symbols("s", positive=True)
rho_e = sp.sqrt(a ** 2 + s_ ** 2)
cS, tS = sp.symbols("c t")
Sfull = b * (1 - cS * sp.cos(rho_e / 2) + tS * a * sp.sin(rho_e / 2) / rho_e)
Gd = sp.diff(Sfull, a)
Hd = sp.diff(Sfull, a, 2)
chk_G = sp.simplify(Gd.subs(cS, sp.Symbol("cc")) - (b * (sp.Symbol("cc") / 2) * a / rho_e * sp.sin(rho_e / 2)
                                              + tS * b * (sp.sin(rho_e / 2) / rho_e
                                                          + a ** 2 * ((rho_e / 2) * sp.cos(rho_e / 2) - sp.sin(rho_e / 2)) / rho_e ** 3)))
say("")
say("  exact-diff check  dS/da  - (2.4)  simplifies to :", sp.simplify(chk_G))
say("  (H check is done numerically in [3] below: central differences in a)")

# ================================================================ [2] endpoint rho=2pi
say("")
say("-" * 100)
say("[2] endpoint rho=2pi (w=0 there, but it is where proportionality dies):")
ep = {r: 2 * sp.pi, S2: 0, C2: -1}
Ae, Be, Ce, De = [sp.simplify(x.subs(ep)) for x in (A, B, Cc, Dd)]
say("  A(2pi) =", Ae, "   B(2pi) =", Be)
say("  C(2pi) =", Ce, "   D(2pi) =", sp.simplify(De))
say("  => C(2pi) = B(2pi)  (KILLER: A vanishes, C does not)")
say("  Wr(2pi) = C*B - D*A =", sp.simplify(Wr.subs(ep)))
Kep = sp.simplify((sp.Symbol("cc") / 2 * Ce + sp.Symbol("tt") * De) / (sp.Symbol("tt") * Be))
say("  H/G at rho=2pi  =", sp.simplify(sp.expand(Kep)),
    "  [ = (c/2)/t + 3/a - 3a/(4pi^2) with c=2*cc/2... see below ]")
Kep0 = sp.simplify(Kep.subs({sp.Symbol("cc"): 2 * sp.Symbol("C0"), sp.Symbol("tt"): sp.Symbol("T0")}))
say("  written with c,t:  H/G|2pi = c/(2t) + 3/a - 3a/(4pi^2)   => kill term  3/a - 3a/(4pi^2)")

# ================================================================ [3] numeric pointwise
say("")
say("-" * 100)
BETA, TH = 2.5, np.pi / 2.0
K0 = 0.5 / np.tan(TH / 2.0)
say(f"[3] pointwise  beta={BETA}  theta=pi/2  K0=(1/2)cot(th/2)={K0:.15f}")
say("    a      rho-range        rho*     G(rho*)      H(rho*)     H/G at rho*    min H/G     max H/G    spread     std H/G     max|Wr|/beta^2")
table = []
for av in [0.0, 0.3, 0.75, 1.2]:
    lo = abs(av)
    rho = np.linspace(lo + 1e-9, TWOPI, 200001)
    G, H = fields(av, rho, BETA, TH)
    # Wronskian in beta-free normalisation: C*B - D*A over beta^2
    c_, t_ = np.cos(0.5 * TH), np.sin(0.5 * TH)
    G_, H_ = G / BETA, H / BETA
    # invert: G_ = (c/2)A+tB ; H_ = (c/2)C+tD  -> Wr from the 2x2 determinant (c,t) solve
    det = (c_ / 2.0) * 0.0  # placeholder to keep flake quiet
    # direct Wronskian from the closed forms (written out again, beta-free):
    A_ = (av / rho) * np.sin(0.5 * rho)
    B_ = wfun(rho) + av ** 2 * wpf(rho) / rho
    C_ = np.sin(0.5 * rho) / rho + av ** 2 * np.cos(0.5 * rho) / (2 * rho ** 2) - av ** 2 * np.sin(0.5 * rho) / rho ** 3
    D_ = 3.0 * av * wpf(rho) / rho - av ** 3 * wpf(rho) / rho ** 3 + av ** 3 * wppf(rho) / rho ** 2
    Wrn = C_ * B_ - D_ * A_
    R = np.divide(H, G, out=np.full_like(H, np.nan), where=np.abs(G) > 1e-12)
    fin = np.isfinite(R)
    # interior reference point rho* = max(|a|, 1.0)
    istar = int(np.argmin(np.abs(rho - max(lo + 0.05, 1.0))))
    say(f"  {av:4.2f}   [{lo:6.3f},{TWOPI:6.3f}]   {rho[istar]:6.3f}  {G[istar]:+11.6f}  {H[istar]:+11.6f}"
        f"  {R[istar]:+11.6f}  {np.nanmin(R[fin]):+11.6f} {np.nanmax(R[fin]):+11.6f}"
        f"  {np.nanmax(R[fin])-np.nanmin(R[fin]):10.4f}  {np.nanstd(R[fin]):10.4f}  {np.max(np.abs(Wrn)):11.6f}")
    table.append((av, G, H, rho))

# endpoint numeric spot check
say("")
say("  endpoint check H/G at rho=2pi  (analytic : (1/2)cot(th/2)+3/a-3a/(4pi^2))")
for av in [0.3, 0.75, 1.2]:
    rr = TWOPI - 1e-7
    G, H = fields(av, np.array([rr]), BETA, TH)
    ana = K0 + 3.0 / av - 3.0 * av / (4 * np.pi ** 2)
    say(f"    a={av:4.2f}  numeric {float(H[0]/G[0]):+.10f}   analytic {ana:+.10f}   diff {float(H[0]/G[0])-ana:+.2e}")

# central-difference check of H = dG/da on the slice (independent of (2.5))
say("")
say("  independent check: H vs d/da[G] by central differences in a (fixed transverse s)")
h = 1e-5
for av in [0.0, 0.3, 0.75, 1.2]:
    ss = np.linspace(0.05, 5.0, 40)
    worst = 0.0
    for s0 in ss:
        rp = np.sqrt((av + h) ** 2 + s0 ** 2)
        rm = np.sqrt((av - h) ** 2 + s0 ** 2)
        r0 = np.sqrt(av ** 2 + s0 ** 2)
        Gp, _ = fields(av + h, np.array([rp]), BETA, TH)
        Gm, _ = fields(av - h, np.array([rm]), BETA, TH)
        _, H0 = fields(av, np.array([r0]), BETA, TH)
        num = (float(Gp[0]) - float(Gm[0])) / (2 * h)
        worst = max(worst, abs(num - float(H0[0])) / max(1.0, abs(float(H0[0]))))
    say(f"    a={av:4.2f}   max rel dev |H - dG/da| = {worst:.2e}")

# ================================================================ [4] moments
say("")
say("-" * 100)
say("[4] moments under dmu_a ~ w exp(-beta S) drho  (GL n=4000)")
say("    a     E[G]        E[H]        Var(G)      Var(H)      Cov(G,H)    K0          E[H]/E[G]   Cov/Var(G)  sd(H)/sd(G)  R^2        resid=1-R^2   pointwise K spread/|K0|")
for av in [0.0, 0.3, 0.75, 1.2]:
    lo = abs(av)
    n = 4000
    tt, wt = gl(n)
    rho = lo + (TWOPI - lo) * tt
    wt = (TWOPI - lo) * wt
    G, H = fields(av, rho, BETA, TH)
    dens = wfun(rho) * np.exp(-Sweight(av, rho, BETA, TH))
    Z = np.dot(wt, dens)
    EG = np.dot(wt, dens * G) / Z
    EH = np.dot(wt, dens * H) / Z
    VG = np.dot(wt, dens * G * G) / Z - EG ** 2
    VH = np.dot(wt, dens * H * H) / Z - EH ** 2
    CV = np.dot(wt, dens * G * H) / Z - EG * EH
    Kmean = EH / EG
    Kcov = CV / VG
    Ksd = np.sqrt(VH / VG)
    R2 = CV ** 2 / (VG * VH)
    # pointwise spread of H/G weighted by the measure (interior)
    rg = np.linspace(lo + 1e-6, TWOPI, 400001)
    Gg, Hg = fields(av, rg, BETA, TH)
    Rg = np.divide(Hg, Gg, out=np.full_like(Hg, np.nan), where=np.abs(Gg) > 1e-12)
    fin = np.isfinite(Rg)
    spread = (np.nanmax(Rg[fin]) - np.nanmin(Rg[fin])) / abs(K0)
    say(f"  {av:4.2f}  {EG:+.6f}  {EH:+.6f}  {VG:+.6e}  {VH:+.6e}  {CV:+.6e}  {K0:+.6f}"
        f"  {Kmean:+.6f}  {Kcov:+.6f}  {Ksd:+.6f}  {R2:.6f}  {1-R2:.3e}   {spread:.4f}")

say("")
say("[5] a=0 control: for exact proportionality all of  E[H]/E[G], Cov/Var(G), sd(H)/sd(G)  must equal")
say("    K0 = (1/2)cot(th/2); at a=0 they do to machine precision (row 1 above).")

# ================================================================ [6] verdict block
say("")
say("-" * 100)
say("[6] VERDICT")
say("  * pointwise proportionality H = K*G fails for every a!=0: Wr(2pi) = beta^2 a^4/(64 pi^4) != 0,")
say("    and H/G|rho=2pi = (1/2)cot(th/2) + 3/a - 3a/(4pi^2).")
say("  * at a=0 the four functions collapse: (A,B,C,D)|a=0 = (0,beta*w,beta*w,0) -> rank 1 -> K = (1/2)cot(th/2).")
say("  * numeric spreads above quantify the deviation (column 'spread/|K0|').")
LOG.close()
