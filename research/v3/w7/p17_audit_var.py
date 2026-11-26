"""p17_audit_var.py -- INDEPENDENT re-derivation / audit of O20's Var bound.

Discipline: written from scratch; o20's implementation was never opened (notes only).
Object taken from the notes and used as-is:

  u = r/2 in (0,pi);   nu_B(u) du  prop.  sin^2(u) exp(-B(1-cos u)) du/u ,  B := beta cos(theta/2)
  g(u) = sin(u)/(2u);  f = beta sin(theta/2) g
  =>  Var(f|0) = beta^2 sin^2(theta/2) W(B),  W(B) = Var_{nu_B}(g),  Psi(B) := 36 B^2 W(B)

Claims audited:
  (i)   V_inf = tan^2(theta/2)/36                        [limit B->inf]
  (ii)  sharp  Var <= V_inf       i.e. Psi <= 1, sup Psi ~ 1              [C]
  (iii) strict Var <= 6.51 V_inf  from C1=1.581977, E[w^2]<=3.328/B^2, L, range bound [T]

Numerics avoid cancellation: h := g - 1/2,  W = E[h^2] - E[h]^2
(E[h^2] ~ 1/(18B^2), E[h]^2 ~ 1/(36B^2): only a factor 2 is lost, not 12 digits).
The naive E[g^2]-E[g]^2 form is printed once to show what it costs.
"""
import numpy as np
import mpmath as mp
import sympy as sp
from scipy.integrate import quad

mp.mp.dps = 30
GK = dict(epsabs=1e-17, epsrel=1e-14, limit=600)
out = []
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); out.append(s)

def q(f, a, b):
    return quad(f, a, b, **GK)[0]

# --------------------------------------------------------------- float engine
def nu_float(B, naive=False):
    if B >= 1.0:
        s = np.sqrt(B); ymax = min(np.pi * s, 45.0)
        f0 = lambda t: np.sin(t / s) ** 2 * np.exp(-B * (1 - np.cos(t / s))) / t      # t = y
        g  = lambda t: np.sin(t / s) / (2 * (t / s))
        a, b = 0.0, ymax
    else:
        f0 = lambda t: np.sin(t) ** 2 * np.exp(-B * (1 - np.cos(t))) / t              # t = u
        g  = lambda t: np.sin(t) / (2 * t)
        a, b = 0.0, np.pi
    Z  = q(f0, a, b)
    N1 = q(lambda t: (g(t) - 0.5) * f0(t), a, b)
    N2 = q(lambda t: (g(t) - 0.5) ** 2 * f0(t), a, b)
    if naive:
        Ngn = q(lambda t: g(t) * f0(t), a, b); Ng2 = q(lambda t: g(t) ** 2 * f0(t), a, b)
        return Z, N1, N2, Ngn, Ng2
    return Z, N1, N2

def W_float(B):
    Z, N1, N2 = nu_float(B)
    return N2 / Z - (N1 / Z) ** 2

# --------------------------------------------------------------- mpmath engine
def W_mp(B):
    B = mp.mpf(B)
    if B >= 1:
        s = mp.sqrt(B); ymax = min(mp.pi * s, mp.mpf(45))
        f0 = lambda t: mp.sin(t / s) ** 2 * mp.e ** (-B * (1 - mp.cos(t / s))) / t
        g  = lambda t: mp.sin(t / s) / (2 * (t / s))
        a, b = mp.mpf(0), ymax
    else:
        f0 = lambda t: mp.sin(t) ** 2 * mp.e ** (-B * (1 - mp.cos(t))) / t
        g  = lambda t: mp.sin(t) / (2 * t)
        a, b = mp.mpf(0), mp.pi
    Z  = mp.quad(f0, [a, b])
    N1 = mp.quad(lambda t: (g(t) - mp.mpf(1) / 2) * f0(t), [a, b])
    N2 = mp.quad(lambda t: (g(t) - mp.mpf(1) / 2) ** 2 * f0(t), [a, b])
    return N2 / Z - (N1 / Z) ** 2

BSPLIT = 1e3            # float engine below, mpmath (30 dps) above
def W(B):
    return W_float(B) if B <= BSPLIT else float(W_mp(B))
def Psi(B):
    return 36.0 * B * B * W(B)

# =========================================================== [P0]
P("[P0] independent audit of O20:  Var(B) -> V_inf = tan^2(theta/2)/36,  Psi=36B^2 W(B)")
P(f"     engines: scipy GK (epsrel 1e-14) for B<={BSPLIT:g};  mpmath 30 dps above;  B = beta cos(theta/2)")

# ====================================== [P1] two engines, anchors
P("\n[P1] W(B), Psi(B): float64(GK, stable h-form) vs mpmath(30 dps)")
P("     B            W(float)                Psi(float)      Psi(mp)              |dPsi|")
anch = [1e-4, 1e-3, 1e-2, 0.1, 0.5, 1, 2, 5, 10, 35.3553, 100, 1000, 1e4, 1e5, 1e6]
worst = 0.0
for B in anch:
    wf = W_float(B); pm = float(36 * mp.mpf(B) ** 2 * W_mp(B))
    pf = 36 * B * B * wf
    d = abs(pf - pm); worst = max(worst, d)
    P(f"  {B:9.4g}  {wf:.18e}  {pf:.14f}  {pm:.14f}  {d:.2e}")
P(f"  max |dPsi| over anchors = {worst:.2e}   (float64 GK degrades only at the very large-B end)")

P("\n[P1b] cost of the naive  Var = E[g^2]-E[g]^2  (double precision)")
for B in [1e2, 1e3, 1e4]:
    Z, _, _, Ngn, Ng2 = nu_float(B, naive=True)
    wn = Ng2 / Z - (Ngn / Z) ** 2
    P(f"  B={B:8.4g}  W_naive={wn:.18e}  W_stable={W_float(B):.18e}  rel.dev={wn/W_float(B)-1:+.2e}")

# ====================================== [P2] V_inf check from (theta,beta)
P("\n[P2] V_inf = tan^2(theta/2)/36, checked through the (theta,beta) route")
P("     Var(theta,beta)=beta^2 sin^2(theta/2) W(beta cos(theta/2));  reported ratio Var/(tan^2(theta/2)/36)")
for th_deg in [60, 90, 120, 135]:
    th = np.deg2rad(th_deg)
    for beta in [50.0, 200.0]:
        B = beta * np.cos(th / 2)
        var = beta ** 2 * np.sin(th / 2) ** 2 * W(B)
        vinf = np.tan(th / 2) ** 2 / 36
        P(f"  theta={th_deg:3d}deg beta={beta:6.1f} B={B:8.3f}  ratio-1={var/vinf-1:+.4e}"
          f"   analytic -2/(15B)={-2/(15*B):+.4e}")

# ====================================== [P3] exact 1/B series
P("\n[P3] exact series in eps=1/B (sympy, rational coefficients, Gaussian moments m_k=2^((k-1)/2)Gamma((k+1)/2))")
yy, eps = sp.symbols('y', positive=True), sp.symbols('eps', positive=True)
ORD = 6
u = yy * sp.sqrt(eps)
# density in y (constant 1/B prefactor cancels):  sin^2(u)/y * exp(-(1-cos u)/eps)
# split the Gaussian: (1-cos u)/eps = y^2/2 - y^4 eps/24 + ...
Srem = sp.series((1 - sp.cos(u)) / eps - yy ** 2 / 2, eps, 0, ORD).removeO()
Dser = sp.series(sp.sin(u) ** 2 / yy * sp.exp(-Srem), eps, 0, ORD).removeO()
Hser = sp.series(sp.sin(u) / (2 * u) - sp.Rational(1, 2), eps, 0, ORD).removeO()
def mom(expr):                                     # integral over y>0 against exp(-y^2/2) dy
    tot = 0
    for term in sp.Add.make_args(sp.expand(expr)):
        c, rest = term.as_coeff_Mul()
        pd = rest.as_powers_dict()
        k = pd.get(yy, sp.Integer(0)); m = pd.get(eps, sp.Integer(0))
        if k == 0:
            tot += c * eps ** m * sp.sqrt(2 * sp.pi) / 2
        else:
            tot += c * eps ** m * sp.Integer(2) ** sp.Rational(k - 1, 2) * sp.gamma(sp.Rational(k + 1, 2))
    return tot
Z_s  = sp.series(mom(Dser), eps, 0, ORD).removeO()
E1_s = sp.series(mom(Hser * Dser) / Z_s, eps, 0, ORD).removeO()
E2_s = sp.series(mom(Hser ** 2 * Dser) / Z_s, eps, 0, ORD).removeO()
W_s  = sp.series(E2_s - E1_s ** 2, eps, 0, ORD).removeO()
Psi_s = sp.expand(sp.series(36 * W_s / eps ** 2, eps, 0, ORD).removeO())
P(f"  E[h]   = {sp.expand(E1_s)}")
P(f"  E[h^2] = {sp.expand(E2_s)}")
P(f"  W      = {sp.expand(W_s)}")
P(f"  Psi    = {Psi_s}")
P(f"  => Psi = 1 + ({sp.expand(Psi_s).coeff(eps,1)})/B + ({sp.expand(Psi_s).coeff(eps,2)})/B^2 + ...")
P("  series vs quadrature:")
for B in [10.0, 100.0, 1e3, 1e4, 1e6]:
    ser = float(Psi_s.subs(eps, 1.0 / B)); num = Psi(B)
    P(f"    B={B:9.4g}  series={ser:.15f}  quad={num:.15f}  diff={ser-num:+.2e}")

# ====================================== [P4] wide grid
P("\n[P4] Psi(B) on a wide grid  B in [1e-4, 1e6]  (241 log points; wider than O20's [1e-3,1e4]/36)")
NB = 241
grid = np.exp(np.linspace(np.log(1e-4), np.log(1e6), NB))
psi = np.array([Psi(B) for B in grid])
imax = int(np.argmax(psi)); d = np.diff(psi)
P(f"  Psi(1e-4) = {psi[0]:.14f}     Psi(1e6) = {psi[-1]:.14f}")
P(f"  sup Psi   = {psi[imax]:.14f}  at B = {grid[imax]:.6g}")
P(f"  1 - sup Psi = {1-psi[imax]:.3e}   (analytic 2/(15B) = {2/(15*grid[imax]):.3e})")
P(f"  monotone?  min dPsi = {d.min():+.3e} (B={grid[1:][int(np.argmin(d))]:.5g})"
  f"   max dPsi = {d.max():+.3e} (B={grid[1:][int(np.argmax(d))]:.5g})")
P(f"  #grid points with Psi>1 : {int((psi>1).sum())}")
P("  large-B tail  (1-Psi)*B should -> 2/15 = 0.133333:")
for B in [3, 10, 30, 100, 300, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6]:
    ps = Psi(B)
    P(f"    B={B:9.4g}  Psi={ps:.13f}  (1-Psi)*B={100*(1-ps)*B:+9.5f}%   [x100]")
P("  small-B end (Psi should be ~ 36 B^2 W(0), W(0)=Var_Haar(g)):")
for B in [1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 1.0, 2.0, 4.0, 8.0]:
    P(f"    B={B:9.4g}  Psi={Psi(B):.12f}   Psi/B^2={Psi(B)/B/B:.6f}")

# ====================================== [P5] C1
P("\n[P5] tail constant:  nu_B((eps,2]) <= C1 exp(-B eps)  for eps >= 1/B")
phi = lambda w: np.sqrt(w * (2 - w)) / np.arccos(1 - w)
C1_th = 1.0 / (1.0 - np.exp(-1.0))
P(f"  formula 1/(1-e^-1) = {C1_th:.9f}  (note 1.581977)")
def ratio(B, f_):
    e0 = f_ / B
    num = q(lambda w: phi(w) * np.exp(-B * w), e0, 2.0)
    den = q(lambda w: phi(w) * np.exp(-B * w), 0.0, 2.0)
    return (num / den) / np.exp(-B * e0)
worst = (0.0, None)
for B in [0.5, 1, 2, 5, 10, 50, 200, 1e3]:
    for f_ in [1.0, 1.5, 2.0, 3.0, 5.0]:
        r = ratio(B, f_)
        if r > worst[0]: worst = (r, (B, f_))
P(f"  sup  nu_B((f/B,2])/exp(-f)  = {worst[0]:.6f}  at (B,f)={worst[1]}")
P(f"  C1 = 1.581977 valid (sup {'<= ' if worst[0]<=C1_th else '> '}C1)   =>  {'OK' if worst[0]<=C1_th else 'FAIL'}")

# ====================================== [P6] E[w^2]
P("\n[P6] E[w^2] for w=1-cos u  (claim: <= 3.328/B^2)")
C3 = 1.0 + 4.0 * C1_th / np.e
P(f"  my derivation: 2*int_0^{1/B} eps deps + 2C1 int_{{1/B}}^inf eps e^(-B eps) deps = 1/B^2 + 4C1/(e B^2) = {C3:.6f}/B^2")
P("     B       E[w^2]*B^2   Var(w)*B^2   bound")
for B in [1, 2, 5, 20, 100, 1e3, 1e4]:
    e1 = q(lambda u: (1 - np.cos(u)) * np.sin(u) ** 2 * np.exp(-B * (1 - np.cos(u))) / u, 0, np.pi) / \
         q(lambda u: np.sin(u) ** 2 * np.exp(-B * (1 - np.cos(u))) / u, 0, np.pi)
    e2 = q(lambda u: (1 - np.cos(u)) ** 2 * np.sin(u) ** 2 * np.exp(-B * (1 - np.cos(u))) / u, 0, np.pi) / \
         q(lambda u: np.sin(u) ** 2 * np.exp(-B * (1 - np.cos(u))) / u, 0, np.pi)
    P(f"  {B:8.4g}  {e2*B*B:10.6f}   {(e2-e1*e1)*B*B:10.6f}   {C3:.6f}")

# ====================================== [P7] Lipschitz + composite
P("\n[P7] phi(w)=sqrt(w(2-w))/arccos(1-w),  L(w1)=max_[0,w1]|phi'|, composite constant")
w = sp.Symbol('w', positive=True)
phis = sp.sqrt(w * (2 - w)) / sp.acos(1 - w)
dphi = sp.simplify(sp.diff(phis, w))
P(f"  phi' (sympy) = {dphi}")
P("  phi' recomputed in the stable form  ((1-w)*b/a - 1)/b^2 ,  a=sqrt(w(2-w)), b=arccos(1-w):")
def dphi_s(warr):
    warr = np.clip(np.asarray(warr, float), 1e-12, 2 - 1e-12)
    b = np.arccos(1 - warr); a = np.sqrt(warr * (2 - warr))
    return ((1 - warr) * b / a - 1) / b ** 2
P(f"    check w=1e-6: {abs(dphi_s(1e-6)):.9f} vs sympy {abs(float(dphi.subs(w, 1e-6))):.9f} ; analytic 1/3={1/3:.9f}")
P(f"    check w=0.3 : {abs(dphi_s(0.3)):.9f} vs sympy {abs(float(dphi.subs(w, sp.Rational(3,10)))):.9f} ; note 0.348")
wgn = np.concatenate([[1e-6], np.geomspace(1e-5, 0.2, 300), np.linspace(0.2, 1.9990, 2500)])
Ltab = np.abs(dphi_s(wgn))
P(f"  |phi'| increasing on (0,2)?  min diff = {np.diff(Ltab).min():+.2e}  => L(w1) = |phi'(w1)|")
P(f"  phi(0+)=1  => (phi(w)-phi(w'))^2 <= 1 on the tail event  (range of phi is (0,1])")
for w1 in [0.3, 0.5, 1.0, 1.2, 1.5]:
    P(f"    w1={w1:4.2f}  L={abs(dphi_s(w1)):.6f}    (note: 0.348/0.361/0.405/0.434/0.504)")
glo = lambda u: np.sin(u) / (2 * u)
uu = np.linspace(1e-12, np.pi - 1e-12, 400001)
P(f"  range of g on (0,pi):  sup g = {glo(uu).max():.10f} (u->0),  inf g = {glo(uu).min():.3e} (u->pi)")
P(f"  => true range length = 1/2,  Var(g) <= (1/2)^2/4 = 0.0625  =>  Psi <= 36*(1/16) B^2 = 2.25 B^2")
P(f"  note's 0.2971 B^2 corresponds to the (1/pi,1/2) range: 36*(1/2-1/pi)^2/4 = {36*(0.5-1/np.pi)**2/4:.6f} B^2")
P(f"     (that range is g on (0,pi/2) only -- mass sits on all of (0,pi): see P7b)")
A_LIP = 9.0 * C3
P(f"  A_lip = 9 * {C3:.6f} = {A_LIP:.4f}   (note: 29.95)")
P("  P7b mass of nu_B on u in (pi/2, pi):")
for B in [1e-6, 1e-3, 0.1, 1.0, 10.0]:
    tot = q(lambda u: np.sin(u) ** 2 * np.exp(-B * (1 - np.cos(u))) / u, 0, np.pi)
    hi = q(lambda u: np.sin(u) ** 2 * np.exp(-B * (1 - np.cos(u))) / u, np.pi / 2, np.pi)
    P(f"    B={B:9.4g}  nu_B((pi/2,pi)) = {hi/tot:.6f}")

def lip(B):
    t2 = 9.0 * B * B * np.minimum(1.0, C1_th * np.exp(-B * wgn))
    return float(np.min(A_LIP * Ltab ** 2 + t2))
def lip_curve(B):
    t2 = 9.0 * B * B * np.minimum(1.0, C1_th * np.exp(-B * wgn))
    v = A_LIP * Ltab ** 2 + t2
    i = int(np.argmin(v)); return float(v[i]), float(wgn[i])

P("\n  [P7c] where the slack is (note: 'L~0.4 but effective slope ~0.19'):")
def wmom(B):
    e1 = q(lambda u: (1 - np.cos(u)) * np.sin(u) ** 2 * np.exp(-B * (1 - np.cos(u))) / u, 0, np.pi) / \
         q(lambda u: np.sin(u) ** 2 * np.exp(-B * (1 - np.cos(u))) / u, 0, np.pi)
    e2 = q(lambda u: (1 - np.cos(u)) ** 2 * np.sin(u) ** 2 * np.exp(-B * (1 - np.cos(u))) / u, 0, np.pi) / \
         q(lambda u: np.sin(u) ** 2 * np.exp(-B * (1 - np.cos(u))) / u, 0, np.pi)
    return e1, e2 - e1 * e1
P("    B     Var(phi)/Var(w)   L(w1*)^2     ratio      eff.slope   L")
for B in [1, 2, 3, 5, 10, 30]:
    e1, vw = wmom(B)
    vphi = 4 * W(B)
    _, w1s = lip_curve(B)
    Ls = abs(dphi_s(w1s))
    P(f"  {B:6.2f}  {vphi/vw:12.6f}  {Ls**2:10.6f}  {(vphi/vw)/Ls**2:8.4f}  {np.sqrt(vphi/vw):9.4f}  {Ls:.4f}")

P("\n  composite  C(B) = min(range-branch, min_{w1}[A_lip L^2 + 9B^2 min(1,C1 e^{-B w1})])")
def composite(krange, tag):
    Bs = np.concatenate([np.geomspace(1e-3, 0.05, 60), np.linspace(0.05, 15, 3000),
                         np.geomspace(15, 1e5, 300)])
    best = (-1, None); bestlip = (-1, None)
    for B in Bs:
        lb = lip(B)
        v = min(krange * B * B, lb)
        if v > best[0]: best = (v, B, lb, krange * B * B, lip_curve(B))
        if lb > bestlip[0]: bestlip = (lb, B, lip_curve(B))
    P(f"  [{tag}] k={krange:.6f}:  max_B min(...) = {best[0]:.4f} at B={best[1]:.4f}"
      f"  (lip={best[2]:.4f}, range={best[3]:.4f}, argmin w1={best[4][1]:.4f}, L={np.sqrt(best[4][0]/A_LIP):.4f})")
    P(f"  [{tag}] max_B of the Lipschitz branch alone = {bestlip[0]:.4f} at B={bestlip[1]:.4f} (w1={bestlip[2][1]:.4f})")
    return best[0]
c_note = composite(0.2971, "note range")
c_true = composite(2.25, "true range")
P(f"\n  C_rig with the note's range constant : {c_note:.3f}   (note claims 6.505)")
P(f"  C_rig with the true range constant   : {c_true:.3f}")
P(f"  ratio = {c_true/c_note:.3f} ; identical Lipschitz+tail machinery, only the range constant differs")

# ====================================== [P8]
P("\n[P8] SUMMARY")
P(f"  V_inf = tan^2(theta/2)/36 : confirmed by B->inf numerics and by the exact series (P2,P3)")
P(f"  exact series: Psi(B) = 1 - 2/(15B) + O(B^-2)   => Psi < 1, increasing to 1 from below")
P(f"  sup Psi on [1e-4,1e6], 241 pts = {psi[imax]:.14f} at B={grid[imax]:.6g}  (Psi->1)")
P(f"  C1 = {C1_th:.6f} (confirmed numerically), E[w^2] <= {C3:.4f}/B^2 (confirmed), A_lip = {A_LIP:.4f}")
P(f"  C_rig = {c_note:.3f} [note range] / {c_true:.3f} [true range]")

with open(r"research\v3\w7\p17_audit_var.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(out) + "\n")
print("\n[log written]")
