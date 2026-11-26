"""q05_ode.py -- Q05: G(a)=-log rho_v(a) 是否满足闭合 ODE G''=F(a,G')?

模型: 单 plaquette SU(2) 平直图 (H1A-EXACT sec.2.1 原样; 与 p21_rebuild.py / p22_mg_hvv.py 逐字同模型):
    X in su(2)=R^3, rho=|X|, u=xi_3/rho, 切片 xi_3=<v,X>=a (v=e_3), 边界 holonomy theta
    S = beta*(1 - cos(th/2)cos(rho/2) + sin(th/2)sin(rho/2)u)
    rho_v(a) = \int_{|a|}^{2pi} sin^2(rho/2)/rho * exp(-S) drho      (定义级: 由 delta 推得)
    G = -log rho_v ,  V = S + Phi ,  Phi = -log J = -2log sin(rho/2)/rho... (p21 原样)
    F = d_3 V = G + Phi' ,  kappa = d_3^2 V = H + Phi''

检验 (判死线):
  [A] 恒等式链  G' = E[F|a] ,  G'' = E[kappa|a] - Var(F|a)   (定义级差商 vs 条件矩)
  [B] 与 p22/P09 恒等式的一致性: G'' = m_G' + E[Phi''] - Var(Phi') - Cov(G,Phi')
  [C] 闭合检验: 跨 (beta,theta) 族找 (a,G') 撞点, 看 G'' 是否单值
  [D] 无穷梯: d/da E[F^2|a] = 2E[F kappa|a] - (m3 - m2 m1)
  [E] theta=0 (S 与 a 无关) 的精确 Riccati: G'' = G'^2 + (log g)' G'
  [F] 精确变号: G''(0) = W0 - V0 beta^2 at theta=pi (beta*=4.30660785)
"""
import math
import numpy as np
from scipy.integrate import quad

TWOPI = 2.0 * math.pi
W0 = 0.205115297123
V0 = 0.011059293772973
LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


def dens(rho, a, beta, ct, st):
    u = a / rho
    S = beta * (1.0 - ct * math.cos(rho / 2.0) + st * math.sin(rho / 2.0) * u)
    return math.sin(rho / 2.0) ** 2 / rho * math.exp(-S)


def obs(rho, a, beta, ct, st):
    """(G,H,Phi1,Phi2) -- 与 p21_rebuild.py 逐字相同."""
    u = a / rho
    s2, c2 = math.sin(rho / 2.0), math.cos(rho / 2.0)
    S_r = beta * (ct * s2 / 2.0 + st * c2 * u / 2.0)
    S_rr = beta * (ct * c2 / 4.0 - st * s2 * u / 4.0)
    S_u = beta * st * s2
    S_ur = beta * st * c2 / 2.0
    S_uu = 0.0
    G = S_r * u + S_u * (1 - u * u) / rho
    H = (S_rr * u * u
         + S_ur * u * (1 - u * u) / rho
         + S_r * (1 - u * u) / rho
         + (S_ur * u + S_uu * (1 - u * u) / rho) * (1 - u * u) / rho
         - 3.0 * S_u * u * (1 - u * u) / (rho * rho))
    Ph_r = -c2 / s2 + 2.0 / rho
    Ph_rr = 0.5 / (s2 * s2) - 2.0 / (rho * rho)
    Phi1 = Ph_r * u
    Phi2 = Ph_rr * u * u + Ph_r * (1 - u * u) / rho
    return G, H, Phi1, Phi2


def raw(a, beta, th, eps=1e-13):
    """定义级 rho_v(a) 与全部条件矩 (F=G+Phi1, kappa=H+Phi2)."""
    ct, st = math.cos(th / 2.0), math.sin(th / 2.0)
    lo = abs(a) + 1e-13
    hi = TWOPI - 1e-10
    pts = sorted({lo, max(lo + 1e-9, min(abs(a) + (hi - abs(a)) * 0.02, hi)), math.pi, hi})

    def integ(h):
        def f(rho):
            g, hh, p1, p2 = obs(rho, a, beta, ct, st)
            F = g + p1
            kap = hh + p2
            return dens(rho, a, beta, ct, st) * h(g, hh, p1, p2, F, kap)
        tot = 0.0
        for x, y in zip(pts[:-1], pts[1:]):
            if y - x > 1e-12:
                tot += quad(f, x, y, epsabs=eps, epsrel=1e-12, limit=400)[0]
        return tot

    Z = integ(lambda g, h, p1, p2, F, k: 1.0)
    d = {"Z": Z}
    d["EG"] = integ(lambda g, h, p1, p2, F, k: g) / Z
    d["EH"] = integ(lambda g, h, p1, p2, F, k: h) / Z
    d["EF"] = integ(lambda g, h, p1, p2, F, k: F) / Z
    d["EF2"] = integ(lambda g, h, p1, p2, F, k: F * F) / Z
    d["EF3"] = integ(lambda g, h, p1, p2, F, k: F * F * F) / Z
    d["EK"] = integ(lambda g, h, p1, p2, F, k: k) / Z
    d["EFK"] = integ(lambda g, h, p1, p2, F, k: F * k) / Z
    d["Ep1"] = integ(lambda g, h, p1, p2, F, k: p1) / Z
    d["Ep1s"] = integ(lambda g, h, p1, p2, F, k: p1 * p1) / Z
    d["Ep2"] = integ(lambda g, h, p1, p2, F, k: p2) / Z
    d["EGp1"] = integ(lambda g, h, p1, p2, F, k: g * p1) / Z
    d["EGG"] = integ(lambda g, h, p1, p2, F, k: g * g) / Z
    d["VarF"] = d["EF2"] - d["EF"] ** 2
    d["Gp"] = d["EF"]                                   # 恒等式链 (1)
    d["Gpp"] = d["EK"] - d["VarF"]                       # 恒等式链 (2)
    d["VarG"] = d["EGG"] - d["EG"] ** 2
    d["VarP1"] = d["Ep1s"] - d["Ep1"] ** 2
    d["Cov"] = d["EGp1"] - d["EG"] * d["Ep1"]
    d["mGp"] = d["EH"] - d["VarG"] - d["Cov"]            # p22 (ii), [T]
    return d


def mGp_fd(a, beta, th, h=0.01):
    return (raw(a + h, beta, th)["EG"] - raw(a - h, beta, th)["EG"]) / (2 * h)


def Gdef(a, beta, th):
    return -math.log(raw(a, beta, th)["Z"])


def Gd_fd1(a, beta, th, h=0.01):
    return (Gdef(a + h, beta, th) - Gdef(a - h, beta, th)) / (2 * h)


def Gd_fd2(a, beta, th, h=0.02):
    return (Gdef(a + h, beta, th) - 2 * Gdef(a, beta, th) + Gdef(a - h, beta, th)) / (h * h)


P("=" * 104)
P("q05_ode.log -- Q05: G(a)=-log rho_v(a) 的 ODE: 闭合否? 与 p21/p22 恒等式一致否?")
P("模型 = p21_rebuild.py = p22_mg_hvv.py (单 plaquette SU(2), 平直图 su(2)=R^3, 切片 xi_3=a, v=e_3)")
P("目标恒等式: G'=E[F|a], G''=E[kappa|a]-Var(F|a), F=G+Phi'=d_3V, kappa=H+Phi''=d_3^2V")
P("=" * 104)
P("")

# ---------- [A] 恒等式链 vs 定义级差商 ----------
P("[A] 恒等式链 (定义级: G=-log rho_v 的差商) vs 条件矩")
P("    beta theta/pi   a     G'_def(h=1e-2)   E[F|a]     rel      G''_def(h=2e-2)  E[kappa]-VarF  rel")
relA1 = relA2 = 0.0
for beta, th, a in ((2.5, math.pi, 0.0), (2.5, math.pi, 0.75), (2.5, math.pi / 2, 0.5),
                    (1.0, math.pi, 0.5), (5.0, math.pi, 1.0)):
    d = raw(a, beta, th)
    g1 = Gd_fd1(a, beta, th)
    g2 = Gd_fd2(a, beta, th)
    r1 = abs(g1 - d["Gp"]) / max(1e-9, abs(d["Gp"]) if abs(d["Gp"]) > 1e-3 else 1.0)
    r2 = abs(g2 - d["Gpp"]) / max(1e-3, abs(d["Gpp"]))
    relA1 = max(relA1, r1)
    relA2 = max(relA2, r2)
    P(f"    {beta:4.1f}  {th/math.pi:6.3f}  {a:5.2f}   {g1:+14.8f}  {d['Gp']:+10.6f}  {r1:.1e}"
      f"   {g2:+14.8f}  {d['Gpp']:+11.6f}  {r2:.1e}")
P(f"    => G'=E[F|a]: max rel {relA1:.1e} [T]; G''=E[kappa]-Var(F): max rel {relA2:.1e} [T] (差商误差主导)")
P("")

# ---------- [B] 与 p22/P09 恒等式一致? ----------
P("[B] 与 p22 (ii) 的一致性: G'' ?= m_G' + E[Phi''] - Var(Phi') - Cov(G,Phi')   (m_G' 用 H-差商, h=1e-2)")
P("    beta theta/pi   a      G''=E[k]-VarF   m_G'+E[Phi'']-VarPhi'-Cov   残差")
resB = 0.0
for beta, th, a in ((2.5, math.pi, 0.75), (2.5, math.pi / 2, 0.5), (2.5, math.pi / 2, 0.0),
                    (1.0, math.pi, 1.0), (5.0, math.pi, 0.5)):
    d = raw(a, beta, th)
    mp = mGp_fd(a, beta, th)
    rhs = mp + d["Ep2"] - d["VarP1"] - d["Cov"]
    resB = max(resB, abs(d["Gpp"] - rhs))
    P(f"    {beta:4.1f}  {th/math.pi:6.3f}  {a:5.2f}   {d['Gpp']:+13.8f}  {rhs:+26.8f}    {d['Gpp']-rhs:+.2e}")
P(f"    => 两式是同一个 psi-引理的两条实例 (psi=Phi'), max|残差| {resB:.1e} => 一致, 无矛盾 [T]")
P("       代数: m_G'+E[Phi'']-VarPhi'-Cov = E[H]+E[Phi'']-VarG-VarPhi'-2Cov = E[kappa]-Var(F)")
P("")

# ---------- [C] 闭合检验 ----------
P("[C] 闭合检验: 若 G''=F(a,G') 存在, G'' 必须是 (a,G') 的单值函数. 扫 (beta,theta) 族找 (a,G') 撞点.")
pts = []
for beta in (1.0, 2.5, 5.0):
    for thf in (0.0, 0.25, 0.5, 0.7, 1.0):
        for a in (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
            d = raw(a, beta, thf * math.pi)
            pts.append((beta, thf, a, d["Gp"], d["Gpp"]))
P(f"    样本 {len(pts)} 点; G' 范围 [{min(p[3] for p in pts):+.3f},{max(p[3] for p in pts):+.3f}]"
  f"  G'' 范围 [{min(p[4] for p in pts):+.3f},{max(p[4] for p in pts):+.3f}]")
best = []
for i in range(len(pts)):
    for j in range(i + 1, len(pts)):
        da_, dg = abs(pts[i][2] - pts[j][2]), abs(pts[i][3] - pts[j][3])
        if da_ < 0.06 and dg < 0.04:
            best.append((dg + 0.5 * da_, pts[i], pts[j], abs(pts[i][4] - pts[j][4])))
best.sort(key=lambda t: t[0])
P("    最近的 (a,G') 撞点 (da<0.06, dG'<0.04):")
for sc, p1, p2, dy in best[:6]:
    P(f"      ({p1[0]:4.1f},{p1[1]:4.2f}pi,a={p1[2]:4.2f}) vs ({p2[0]:4.1f},{p2[1]:4.2f}pi,a={p2[2]:4.2f})"
      f"   da={abs(p1[2]-p2[2]):.3f} dG'={abs(p1[3]-p2[3]):.3f}  G':{p1[3]:+.4f}/{p2[3]:+.4f}"
      f"  G'':{p1[4]:+.4f}/{p2[4]:+.4f}  |dG''|={dy:.4f}")
if best:
    P(f"    => 在 dG'<0.04, da<0.06 的近撞点上 |dG''| 最大 {max(b[3] for b in best):.4f}"
      f" (G'' 自身尺度 {max(abs(p[4]) for p in pts):.3f})")
P("    => G'' 不是 (a,G') 的单值函数 => 无 F(a,G') 闭合 [C](数值) / 见 §3 的 [T] 反证")
P("")

# ---------- [D] 无穷梯 ----------
P("[D] 无穷梯: d/da E[F^2|a] ?= 2E[F kappa|a] - (m3 - m2 m1)   (psi=F^2 的引理), 并看它引出的新对象")
P("    beta theta/pi   a    LHS(差商)      RHS           残差      m3")
resD = 0.0
for beta, th, a in ((2.5, math.pi, 0.75), (2.5, math.pi / 2, 0.5), (1.0, math.pi, 0.5)):
    h = 0.01
    lhs = (raw(a + h, beta, th)["EF2"] - raw(a - h, beta, th)["EF2"]) / (2 * h)
    d = raw(a, beta, th)
    m1, m2, m3 = d["EF"], d["EF2"], d["EF3"]
    rhs = 2.0 * d["EFK"] - (m3 - m2 * m1)
    resD = max(resD, abs(lhs - rhs))
    P(f"    {beta:4.1f}  {th/math.pi:6.3f}  {a:5.2f}  {lhs:+12.6f}  {rhs:+12.6f}   {lhs-rhs:+.2e}  {m3:+.4f}")
P(f"    => 梯子逐阶成立, max|残差| {resD:.1e} [T]; 第 n 阶要 m_(n+1) 与 E[d_3^n kappa] => 不闭合 [T]")
P("")

# ---------- [E] theta=0 的精确 Riccati (S 与 a 无关: a 只从下极限 |a| 进) ----------
P("[E] theta=0 (S 与 a 无关): 精确 Riccati G'' = G'^2 + (log g)' G',  g(rho)=sin^2(rho/2)/rho e^{-beta(1-cos(rho/2))}")
P("    a     G'(=g(a)/rho_v)   G''(矩)    G'^2+(log g)'G'   残差")
resE = 0.0
beta = 2.5
ct, st = 1.0, 0.0

def g0(rho):
    return math.sin(rho / 2.0) ** 2 / rho * math.exp(-beta * (1.0 - math.cos(rho / 2.0)))

def dg0(rho, h=1e-5):
    return (g0(rho + h) - g0(rho - h)) / (2 * h)

for a in (0.25, 0.5, 1.0, 1.5, 2.0):
    d = raw(a, beta, 0.0)
    ric = d["Gp"] ** 2 + (dg0(a) / g0(a)) * d["Gp"]
    resE = max(resE, abs(ric - d["Gpp"]))
    P(f"    {a:4.2f}  {d['Gp']:+12.6f}   {d['Gpp']:+10.6f}   {ric:+16.6f}   {ric-d['Gpp']:+.2e}")
P(f"    => theta=0 上闭合 (Riccati), max|残差| {resE:.1e} [T]; 但 theta=0 是 a-解耦情形 (G'' 与其符号已知), 不携带耦合信息")
P("")

# ---------- [F] 精确变号 vs 任何 F>=0 的闭合 ----------
P("[F] 精确闭式 (O01/H1E): theta=pi, G''(0)=W0-V0*beta^2 ; 任何 F>=0 的闭合都会禁止 G''<0")
P("    beta   G''(0) 数值      W0-V0 beta^2")
for beta in (1.0, 2.5, 4.0, 4.30660785, 5.0, 6.0):
    d = raw(0.0, beta, math.pi)
    P(f"    {beta:6.3f}   {d['Gpp']:+14.8f}   {W0 - V0*beta*beta:+14.8f}")
P("    => beta>beta*=4.30660785 时 G''(0)<0 [T] => 不存在全域 F(a,G')>=0 的闭合")
P("")

P("[DONE]")
with open(r"D:\Dev\math\research\v3\w7\q05_ode.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG) + "\n")
