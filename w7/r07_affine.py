"""R07-AFFINE: is  H = K(a) G + C(a)  an IDENTITY or an approximation?

Model (H1A §2.1 / Q04 (2.4)(2.5)):
    S_W = beta[ 1 - c cos(rho/2) + t a w ] ,  w(rho)=sin(rho/2)/rho,  rho=sqrt(a^2+s^2)
    G = dS/da ,  H = d^2S/da^2
Part A: sympy, exact.  Part B: numerics with the SAME measure as r01
        (dens = w exp(-beta S), GL n=6000), reproduce r01's K,C and test closed forms.
"""
import numpy as np
import sympy as sp

LOG = open(r"D:\Dev\math\w7\r07_affine.log", "w", encoding="utf-8")


def say(*xs):
    line = " ".join(str(z) for z in xs)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


say("=" * 100)
say("R07-AFFINE  (beta=2.5, theta=pi/2;  model S_W = beta[1 - c cos(rho/2) + t a w])")
say("=" * 100)

# ============================ PART A : sympy exact ============================
x, a = sp.symbols("rho a", positive=True)
c, t = sp.symbols("c t")  # c=cos(th/2), t=sin(th/2), with c^2+t^2=1
U = sp.symbols("U", positive=True)  # u = rho/2

w = sp.sin(x / 2) / x
wp = sp.diff(w, x)
wpp = sp.diff(w, x, 2)

A = a * w
B = w + a ** 2 * wp / x
Co = sp.sin(x / 2) / x + (a ** 2 / (2 * x ** 2)) * sp.cos(x / 2) - (a ** 2 / x ** 3) * sp.sin(x / 2)
D = 3 * a * wp / x - a ** 3 * wp / x ** 3 + a ** 3 * wpp / x ** 2

G = c / 2 * A + t * B
H = c / 2 * Co + t * D

say("")
say("-- A1  structural identities (all should be 0) --")
say("[A1a] r01's A  minus  a*w(ρ)            :", sp.simplify(A - a * w))
say("[A1b] r01's C  minus  B                 :", sp.simplify(Co - B))
say("[A1c] d_aa[-cos(ρ/2)] - B/2             :",
    sp.simplify(sp.diff(-sp.cos(x / 2), a, 2) - B / 2))
say("[A1d] d_a[a*w(ρ)] - B   (=A' vs B)      :",
    sp.simplify(sp.diff(a * sp.sin(sp.sqrt(a ** 2 + sp.symbols('s') ** 2) / 2)
                        / sp.sqrt(a ** 2 + sp.symbols('s') ** 2), a).subs(
        sp.sqrt(a ** 2 + sp.symbols('s') ** 2), x) - B))
say("  =>  {A,B,C,D} 的秩塌缩不是 a=0 现象: 恒有 A = a w,  C = B  (= d_a B ... 见 A1c)")

K0 = c / (2 * t)
Delta = sp.simplify(H - K0 * G)
say("")
say("-- A2  exact decomposition  H = K0 G + beta*Delta,  K0 = c/(2t) --")
say("[A2a] H - K0 G - [ t D - (c^2/4t) A ]   :", sp.simplify(Delta - (t * D - c ** 2 / (4 * t) * A)))
say("[A2b] Delta =", sp.simplify(Delta))
say("[A2c] Delta 小 a 展开:", sp.series(sp.simplify(Delta), a, 0, 3))

# ---- trig-basis decomposition:  f = P(ρ) sin(ρ/2) + Q(ρ) cos(ρ/2) + R(ρ) ----
S, Cv = sp.symbols("S Cv")


def decomp(e):
    e2 = sp.expand(e.subs({sp.sin(x / 2): S, sp.cos(x / 2): Cv}))
    P = sp.expand(e2.coeff(S))
    Q = sp.expand(e2.coeff(Cv))
    R = sp.expand(e2.subs({S: 0, Cv: 0}))
    return sp.simplify(P), sp.simplify(Q), sp.simplify(R)


PG, QG, RG = decomp(G)
PH, QH, RH = decomp(H)
say("")
say("-- A3  trig-basis decomposition  f = P sin(ρ/2) + Q cos(ρ/2) + R --")
say("  P_G =", PG)
say("  Q_G =", QG, "   R_G =", RG)
say("  P_H =", PH)
say("  Q_H =", QH, "   R_H =", RH)

Kc = sp.simplify(QH / QG)
Ks = sp.simplify(PH / PG)
say("  Q_H/Q_G =", Kc)
say("  P_H/P_G =", Ks)
say("  (Q_H/Q_G) - (P_H/P_G)  simplified :", sp.simplify(Kc - Ks))
say("  P_H - (Q_H/Q_G)*P_G    simplified :", sp.simplify(PH - Kc * PG))
say("  K_s - K_c  numerator over common denom, series in a:",
    sp.series(sp.simplify(Ks - Kc), a, 0, 2))

# ---- Wronskians ----
Gx = sp.diff(G, x)
Hx = sp.diff(H, x)
Gxx = sp.diff(Gx, x)
Hxx = sp.diff(Hx, x)
W1 = sp.simplify(G * Hx - Gx * H)              # = 0  <=> H,G proportional
W2 = sp.simplify(Gx * Hxx - Gxx * Hx)          # = 0  <=> H-kG constant (affine!)
say("")
say("-- A4  Wronskians (w.r.t. rho; same zero set as w.r.t. u=rho/2) --")
say("  W1 = W(G,H)      =", W1)
say("  W2 = W(G',H')    =", W2)
say("  W1 series in a (a->0):", sp.series(W1, a, 0, 2))
say("  W2 series in a (a->0):", sp.series(W2, a, 0, 2))

W3 = sp.Matrix([[sp.sin(x / 2), sp.cos(x / 2), 1],
                [sp.diff(sp.sin(x / 2), x), sp.diff(sp.cos(x / 2), x), 0],
                [sp.diff(sp.sin(x / 2), x, 2), sp.diff(sp.cos(x / 2), x, 2), 0]]).det()
say("  W{sin(ρ/2),cos(ρ/2),1} =", sp.simplify(W3), "(nonzero => 1 is NOT in their span)")

# numeric spot check of W1, W2 (c=cos(pi/4), t=sin(pi/4))
subs_num = {c: np.cos(np.pi / 4), t: np.sin(np.pi / 4)}
for av in [0.0, 0.3, 1.2]:
    for xv in [2.5, 4.0]:
        w1 = complex(W1.subs({a: av, x: xv}).subs(subs_num).evalf())
        w2 = complex(W2.subs({a: av, x: xv}).subs(subs_num).evalf())
        say(f"  a={av} rho={xv}:  W1={w1:.6e}   W2={w2:.6e}")

# ---- u-forms (u = rho/2, g(u)=sin u/u) ----
g = sp.Function("g")
gain = sp.simplify(A.subs(x, 2 * U) * 2 * U / sp.sin(U))  # A = a*g/(2*?)  -> report raw
say("")
say("-- A5  u-form (u=rho/2, g=sin(u)/u):  A = a*sin(u)/(2u) = (a/2) g(u) --")
say("  A  in u :", sp.simplify(A.subs(x, 2 * U)))
say("  B  in u :", sp.simplify(B.subs(x, 2 * U)))
say("  Co in u :", sp.simplify(Co.subs(x, 2 * U)))
say("  D  in u :", sp.simplify(D.subs(x, 2 * U)))
say("  G  in u :", sp.simplify(G.subs(x, 2 * U)))
say("  H  in u :", sp.simplify(H.subs(x, 2 * U)))
say("  G(u)|a=0 :", sp.simplify(G.subs({x: 2 * U, a: 0})), "   [= c/2 * sin(u)/(2u)]")
say("  H(u)|a=0 :", sp.simplify(H.subs({x: 2 * U, a: 0})), "   [= (c/2) * sin(u)/(2u)]")

# ============================ PART B : numerics ============================
TWOPI = 2.0 * np.pi
_G = {}


def gl(n):
    if n not in _G:
        xx, wg = np.polynomial.legendre.leggauss(n)
        _G[n] = (0.5 * (xx + 1.0), 0.5 * wg)
    return _G[n]


def wfun(r):
    r = np.asarray(r, float)
    out = np.empty_like(r)
    m = r != 0.0
    out[m] = np.sin(0.5 * r[m]) / r[m]
    out[~m] = 0.5
    return out


def wpf(r):
    return (0.5 * r * np.cos(0.5 * r) - np.sin(0.5 * r)) / np.where(r == 0.0, 1.0, r) ** 2


def wppf(r):
    rr = np.where(r == 0.0, 1.0, r)
    return -np.sin(0.5 * r) / (4.0 * rr) - np.cos(0.5 * r) / rr ** 2 + 2.0 * np.sin(0.5 * r) / rr ** 3


BETA, TH = 2.5, np.pi / 2.0
CC, TT = np.cos(0.5 * TH), np.sin(0.5 * TH)
K0n = 0.5 / np.tan(TH / 2.0)


def fields(av, rho):
    u = av / rho
    w, wp, wpp = wfun(rho), wpf(rho), wppf(rho)
    Gv = BETA * ((CC / 2.0) * u * np.sin(0.5 * rho) + TT * (w + (av * av / rho) * wp))
    Hv = BETA * ((CC / 2.0) * (np.sin(0.5 * rho) / rho + (av * av / (2.0 * rho ** 2)) * np.cos(0.5 * rho)
                               - (av * av / rho ** 3) * np.sin(0.5 * rho))
                 + TT * ((3.0 * av / rho) * wp - (av ** 3 / rho ** 3) * wp + (av ** 3 / rho ** 2) * wpp))
    Aa = av * w
    Dd = 3.0 * av * wp / rho - av ** 3 * wp / rho ** 3 + av ** 3 * wpp / rho ** 2
    Delt = TT * Dd - (CC ** 2 / (4.0 * TT)) * Aa
    return Gv, Hv, Delt, Aa, Dd


def grid(av, N):
    lo = abs(av)
    tt, wt = gl(N)
    rho = lo + (TWOPI - lo) * tt
    wt = (TWOPI - lo) * wt
    dens = wfun(rho) * np.exp(-BETA * (1.0 - CC * np.cos(0.5 * rho) + TT * av * wfun(rho)))
    p = wt * dens
    return rho, p / p.sum()


def moments(y1, y2, p):
    m1 = np.sum(p * y1)
    m2 = np.sum(p * y2)
    v1 = np.sum(p * (y1 - m1) ** 2)
    v2 = np.sum(p * (y2 - m2) ** 2)
    cv = np.sum(p * (y1 - m1) * (y2 - m2))
    return m1, m2, v1, v2, cv


say("")
say("=" * 100)
say("PART B  measure = w exp(-beta S) drho on [|a|,2pi], GL n=6000 (r01 convention)")
say("=" * 100)
say("")
say("  a     E[G]      E[H]      Var(G)     Var(H)     Cov        K_LS      C_LS      "
    "R2(aff)    relres")
rows = {}
for av in [0.0, 0.3, 0.75, 1.2]:
    rho, p = grid(av, 6000)
    Gv, Hv, Delt, Aa, Dd = fields(av, rho)
    EG, EH, VG, VH, CV = moments(Gv, Hv, p)
    K = CV / VG
    C = EH - K * EG
    res = Hv - K * Gv - C
    R2 = 1.0 - np.sum(p * res ** 2) / VH
    relres = np.sqrt(np.sum(p * res ** 2) / np.sum(p * Hv ** 2))
    rows[av] = (EG, EH, VG, VH, CV, K, C, R2, relres, p, Gv, Hv, Delt)
    say(f"{av:5.2f} {EG:+.6f} {EH:+.6f} {VG:.6e} {VH:.6e} {CV:+.6e} "
        f"{K:+.8f} {C:+.8f} {R2:.8f} {relres:.3e}")

say("")
say("-- B1  closed form for r01's K and C :  K_LS = K0 + beta*Cov(G,Delta)/Var(G),")
say("     C_LS = beta*( E[Delta] - (Cov(G,Delta)/Var(G))*E[G] ),  Delta = t D - (c^2/4t) A")
say("  a     K_LS(num)   K0+bCov/Var   diff        C_LS(num)   C_pred      diff")
for av, r in rows.items():
    EG, EH, VG, VH, CV, K, C, R2, relres, p, Gv, Hv, Delt = r
    ED = np.sum(p * Delt)
    cvd = np.sum(p * (Gv - EG) * (Delt - ED))
    Kpred = K0n + BETA * cvd / VG
    Cpred = BETA * (ED - (cvd / VG) * EG)
    say(f"{av:5.2f} {K:+.8f} {Kpred:+.8f} {K-Kpred:+.2e}  {C:+.8f} {Cpred:+.8f} {C-Cpred:+.2e}")

say("")
say("-- B2  exact variance identity:  Var(H) = K_LS^2 Var(G) + beta^2 Var(Delta)(1-corr^2(G,Delta))")
say("  a     Var(H)-K^2Var(G)   beta^2*Var(D)*(1-r^2)    rel./Var(H)")
for av, r in rows.items():
    EG, EH, VG, VH, CV, K, C, R2, relres, p, Gv, Hv, Delt = r
    ED = np.sum(p * Delt)
    VD = np.sum(p * (Delt - ED) ** 2)
    cvd = np.sum(p * (Gv - EG) * (Delt - ED))
    corr2 = cvd ** 2 / (VG * VD) if VD > 0 else 0.0
    lhs = VH - K ** 2 * VG
    rhs = BETA ** 2 * VD * (1.0 - corr2)
    say(f"{av:5.2f} {lhs:+.6e}       {rhs:+.6e}         {abs(lhs-rhs)/VH:.3e}")

say("")
say("-- B3  residual order in a:  relres = rms(H-KG-C)/rms(H)  ;  also sup-norm --")
say("  a       relres      relres/a    relres/a^2   max|H-KG-C|/max|H|   N=6000 vs N=1500")
for av in [0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.2]:
    out = []
    for N in (6000, 1500):
        rho, p = grid(av, N)
        Gv, Hv, Delt, Aa, Dd = fields(av, rho)
        EG, EH, VG, VH, CV = moments(Gv, Hv, p)
        K, C = CV / VG, 0.0
        C = EH - K * EG
        res = Hv - K * Gv - C
        relres = np.sqrt(np.sum(p * res ** 2) / np.sum(p * Hv ** 2))
        sup = np.max(np.abs(res)) / np.max(np.abs(Hv))
        out.append((relres, sup))
    say(f"{av:5.2f}  {out[0][0]:.6e}  {out[0][0]/av:.6e}  {out[0][0]/av**2:.6e}   "
        f"{out[0][1]:.6e}      {out[0][0]:.6e} vs {out[1][0]:.6e}")

say("")
say("-- B4  what happens at G's interior zero rho_z:  H = K0*0 + beta*Delta  =>  H(rho_z)=beta*Delta(rho_z)")
say("  a      rho_z        G(rho_z)      H(rho_z)      beta*Delta(rho_z)   K_c(rho_z)")
from scipy.optimize import brentq  # noqa: E402

for av in [0.3, 0.75, 1.2]:
    lo = abs(av)
    f = lambda r: fields(av, np.array([r]))[0][0]
    rz = brentq(f, lo + 1e-9, TWOPI - 1e-9)
    Gz, Hz, Dz, _, _ = fields(av, np.array([rz]))
    kc = K0n + 3.0 / av - 3.0 * av / rz ** 2
    say(f"{av:5.2f}  {rz:.9f}  {Gz[0]:+.3e}  {Hz[0]:+.8f}  {BETA*Dz[0]:+.8f}    {kc:+.6f}")

say("")
say("-- B5  pointwise 'would-be' K:  K_c(rho)=Q_H/Q_G = K0 + 3/a - 3a/rho^2  (NOT constant)")
say("  a     K_c(2.5)   K_c(4.0)   K_c(6.0)   K_s(2.5)   K_s(4.0)   K_s(6.0)")
for av in [0.3, 0.75, 1.2]:
    kcs, kss = [], []
    for xv in (2.5, 4.0, 6.0):
        sub = {a: av, x: xv}
        kcs.append(float((Kc.subs(sub).subs(subs_num)).evalf()))
        kss.append(float((Ks.subs(sub).subs(subs_num)).evalf()))
    say(f"{av:5.2f}  " + "  ".join(f"{v:+.6f}" for v in kcs) + "  " +
        "  ".join(f"{v:+.6f}" for v in kss))

say("")
say("-- A6  numerical check of A1c (correct version, rho=sqrt(a^2+s^2)):  d_aa[-cos(rho/2)] - B/2 --")
s = sp.symbols("s", positive=True)
rho_s = sp.sqrt(a ** 2 + s ** 2)
A1c = sp.simplify(sp.diff(-sp.cos(rho_s / 2), a, 2)
                  - (B / 2).subs(x, rho_s))
say("   d_aa[-cos(rho/2)] - B/2  =", A1c)
say("")
say("DONE")
LOG.close()
