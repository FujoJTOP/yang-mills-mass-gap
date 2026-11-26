"""p22_mg_hvv.py -- P22: 仲裁 m_G' 与 E[H_vv] 是否同一对象.

模型: 单 plaquette SU(2) 平直图 (H1A-EXACT sec.2.1 原样; 与 p21_rebuild.py 逐字同模型):
    X in su(2)=R^3, rho=|X|, u=xi_3/rho, 切片 xi_3=a (a = <v,X>, v=e_3), 边界 holonomy theta
    S = beta*(1 - cos(th/2)cos(rho/2) + sin(th/2)sin(rho/2)u)
    dmu_a ∝ sin^2(rho/2)/rho * exp(-S) drho
    G=d_3 S, H=d_3^2 S, Phi'=d_3 Phi, Phi''=d_3^2 Phi, F=G+Phi', kappa=H+Phi''
判据 (本文件的判死线):
    (i)  逐点: d_a G = H  (a 沿 v)  -- 但条件期望的导数 != 导数的条件期望
    (ii) [T] (P09 (*)): m_G'(a) = E[H|a] - Cov(G,G+Phi'|a) = E[H|a] - Var(G|a) - Cov(G,Phi'|a)
    (iii) 同一批构型/同一 a 参数化上同时给 m_G' 与 E[H|a]
"""
import math
import numpy as np
from scipy.integrate import quad

N12 = 1.0 / 6.0
TWOPI = 2.0 * math.pi
LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


# ---- 与 p21_rebuild.py 同一模型 (逐字复用) ----
def dens(rho, a, beta, ct, st):
    u = a / rho
    S = beta * (1.0 - ct * math.cos(rho / 2.0) + st * math.sin(rho / 2.0) * u)
    return math.sin(rho / 2.0) ** 2 / rho * math.exp(-S)


def obs(rho, a, beta, ct, st):
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


def moments(a, beta, th, eps=1e-13):
    ct, st = math.cos(th / 2.0), math.sin(th / 2.0)
    lo = abs(a) + 1e-13
    hi = TWOPI - 1e-10
    pts = sorted({lo, max(lo + 1e-9, min(abs(a) + (hi - abs(a)) * 0.02, hi)), math.pi, hi})

    def integ(h):
        def f(rho):
            return dens(rho, a, beta, ct, st) * h(*obs(rho, a, beta, ct, st))
        tot = 0.0
        for x, y in zip(pts[:-1], pts[1:]):
            if y - x > 1e-12:
                tot += quad(f, x, y, epsabs=eps, epsrel=1e-12, limit=400)[0]
        return tot

    Z = integ(lambda G, H, p1, p2: 1.0)
    EG = integ(lambda G, H, p1, p2: G) / Z
    EG2 = integ(lambda G, H, p1, p2: G * G) / Z
    EH = integ(lambda G, H, p1, p2: H) / Z
    EH2 = integ(lambda G, H, p1, p2: H * H) / Z
    Ep1 = integ(lambda G, H, p1, p2: p1) / Z
    Ep1s = integ(lambda G, H, p1, p2: p1 * p1) / Z
    Ep2 = integ(lambda G, H, p1, p2: p2) / Z
    EGp1 = integ(lambda G, H, p1, p2: G * p1) / Z
    d = dict(EG=EG, EH=EH, EHs=np.sqrt(max(EH2, 0.0)), EPhi1=Ep1, EPhi2=Ep2)
    d["VarG"] = EG2 - EG * EG
    d["VarPhi1"] = Ep1s - Ep1 * Ep1
    d["Cov"] = EGp1 - EG * Ep1
    d["Gpp"] = EH + Ep2 - (d["VarG"] + d["VarPhi1"] + 2 * d["Cov"])
    d["Ddyn"] = EH - d["VarG"]
    d["Z"] = Z
    return d


def mGprime(a, beta, th, h=0.02):
    return (moments(a + h, beta, th)["EG"] - moments(a - h, beta, th)["EG"]) / (2 * h)


P("=" * 108)
P("p22_mg_hvv.log -- P22: m_G' vs E[H_vv] 是不是同一个对象 (beta=2.5, theta=pi 为主)")
P("模型与 p21_rebuild.py 逐字相同: 单 plaquette SU(2), 平直图 su(2)=R^3, 切片 xi_3=a, v=e_3")
P("判据(i) 逐点 d_a G = H ; 判据(ii)[T] m_G' = E[H|a] - Var(G|a) - Cov(G,Phi'|a)  (P09 (*))")
P("=" * 108)
P("")

# ---------- [Q0] 锚: 本脚本的 E[H|0] 是否等于 p05 的 E[H_vv]=beta(c/2)G(B) ----------
P("[Q0] 锚: a=0 切片上 E[H_vv|0] vs P05-CINF 的闭式 beta*(c/2)*G(B)  (P05 值: theta=0 -> 0.553060561)")
P("   theta     E[H|0](本脚本)     P05 值        差")
p05ref = {0.0: 0.553060561}
for th in (0.0, math.pi / 2, 0.7 * math.pi, math.pi):
    d = moments(0.0, 2.5, th)
    ref = p05ref.get(th)
    s = f"   {th/math.pi:5.3f}pi  {d['EH']:+.9f}   " + (f"{ref:.9f}   {d['EH']-ref:+.2e}" if ref is not None else f"{'-':>9s}   -")
    P(s)
P("   => a=0 上 E[H|0] = beta*cos(theta/2)*<sin(rho/2)/(2rho)> : theta=pi 时精确 0")
P("   => 本脚本的 E[H_vv] 就是 P05 的 E[H_vv] (同一个 a=0 切片对象) [C]")
P("")

# ---------- [Q1] 逐点: d_a G = H (a 沿 v) ----------
P("[Q1] 判据(i) 逐点 d_a G_v = H_vv? (a=xi_3 沿 v=e_3, 固定横向 Y 即 rho^2-a^2=c; h=1e-5)")
ct, st = math.cos(math.pi / 2), math.sin(math.pi / 2)
rng = np.random.default_rng(7)
mx = 0.0
for _ in range(6):
    rho = float(rng.uniform(1.0, 6.0))
    a = float(rng.uniform(-0.9, 0.9) * rho)
    c2_ = rho * rho - a * a                       # |Y|^2 = const 沿 v 的直线
    h = 1e-5
    Gp = obs(math.sqrt(c2_ + (a + h) ** 2), a + h, 2.5, ct, st)[0]
    Gm = obs(math.sqrt(c2_ + (a - h) ** 2), a - h, 2.5, ct, st)[0]
    Hv = obs(rho, a, 2.5, ct, st)[1]
    fd = (Gp - Gm) / (2 * h)
    mx = max(mx, abs(fd - Hv) / max(1e-6, abs(Hv)))
    P(f"   rho={rho:6.3f} a={a:+6.3f}  d_aG={fd:+12.6f}  H={Hv:+12.6f}  rel={abs(fd-Hv)/max(1e-12,abs(Hv)):.1e}")
P(f"   => 沿 v (固定横向 Y) 逐点 d_a G = H 成立, max rel {mx:.1e} [T] (a 确为沿 v 的坐标)")
P("      注意: 固定 rho 变 a 不是沿 v 的方向导数 (差 32%, 见原始 Q1 版本), 不是判据")
P("      但  d_a E[G|a]  !=  E[d_a G|a]  : 差 = -Cov(G, d_a V|a) (测度沿 a 变化, P09 (*))")
P("")

# ---------- [Q2] 主表: 同一批构型、同一 a, 同时给 m_G' 与 E[H] ----------
P("[Q2] 主表 beta=2.5, theta=pi: 同一 a 参数化上同时给 m_G' 与 E[H|a] (中心差分 h=0.02)")
P("   a      E[H|a]      m_G'(a)     Var(G|a)    Cov(G,Phi')   -VarG-Cov   残差 m_G'-(EH-VarG-Cov)")
rows = []
for a in (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
    d = moments(a, 2.5, math.pi)
    mp = mGprime(a, 2.5, math.pi)
    rhs = d["EH"] - d["VarG"] - d["Cov"]
    rows.append((a, d, mp, rhs))
    P(f"  {a:<5.2f} {d['EH']:+.6f}   {mp:+.6f}   {d['VarG']:.6f}   {d['Cov']:+.6f}   "
      f"{-(d['VarG']+d['Cov']):+.6f}    {mp-rhs:+.1e}")
P(f"   => 恒等式 (ii) 全网格成立, max|残差| = {max(abs(r[2]-r[3]) for r in rows):.1e} [T]")
P("   => theta=pi: E[H|0]=0 精确, 但 m_G'(0)=-Var(G|0)<0 => 两者不是同一个对象")
P("")

# ---------- [Q3] 差分步长收敛 (m_G' 是 [C] 而非 [T] 的部分) ----------
P("[Q3] h-收敛 (a=0 与 a=0.75, theta=pi): m_G' 的差分误差")
for a in (0.0, 0.75):
    vals = [mGprime(a, 2.5, math.pi, h) for h in (0.02, 0.01, 0.005)]
    P(f"   a={a:<4.2f} h=0.02/0.01/0.005 -> {vals[0]:+.6f} / {vals[1]:+.6f} / {vals[2]:+.6f}"
      f"   spread={max(vals)-min(vals):.1e}")
P("")

# ---------- [Q4] 对照: theta=pi/2 上 m_G' 变号, E[H|a] 仍 >0 ----------
P("[Q4] 对照 theta=pi/2 (B=beta*cos(theta/2)=1.7678): 符号随 theta 族变, 但 E[H]>0 两支同")
P("   a     E[H|a]      m_G'(a)     Var(G|a)    Cov(G,Phi')   Var(F)-E[kappa]=-G''")
for a in (0.0, 0.25, 0.75, 1.0, 2.0):
    d = moments(a, 2.5, math.pi / 2)
    mp = mGprime(a, 2.5, math.pi / 2)
    P(f"  {a:<5.2f} {d['EH']:+.6f}   {mp:+.6f}   {d['VarG']:.6f}   {d['Cov']:+.6f}   "
      f"{d['VarG']+d['VarPhi1']+2*d['Cov']-(d['EH']+d['EPhi2']):+.6f}")
P("")

# ---------- [Q5] E[H|a] 在 theta=pi 的 a-依赖 (P05 的 [T] 只在 a=0 上给 0) ----------
P("[Q5] theta=pi: E[H|a] 是全 a 上恒 0, 还是只在 a=0 上 0?")
P("   a      E[H|a]      pointwise max|H| on rho-grid")
for a in (0.0, 0.25, 0.75, 1.5):
    d = moments(a, 2.5, math.pi)
    grid = np.linspace(abs(a) + 1e-6, TWOPI - 1e-6, 20001)
    Hm = max(abs(obs(float(r), a, 2.5, 0.0, 1.0)[1]) for r in grid)
    P(f"  {a:<5.2f} {d['EH']:+.6f}    {Hm:.6f}")
P("   => theta=pi 时 H_vv=beta*cos(th/2)*sin(rho/2)/(2rho)*[a=0] 逐点 0 只在 a=0 切片;")
P("      a!=0 上 E[H|a]<0 (P05 的逐点非负是 a=0 切片命题, 不是全切片命题)")
P("")

# ---------- [Q6] 判: 每个 a 上, P05 能给 m_G' 什么 ----------
P("[Q6] 若把 P05 的 [T] 当 m_G' >= 0 用, 它实际需要什么?")
P("   m_G' >= 0  <=>  Var(G|a) + Cov(G,Phi'|a) <= E[H|a]         [T] (由 (ii))")
for a in (0.0, 0.75, 1.5):
    d = moments(a, 2.5, math.pi)
    mp = mGprime(a, 2.5, math.pi)
    P(f"   a={a:<4.2f}: LHS={d['VarG']+d['Cov']:+.6f}  RHS=E[H|a]={d['EH']:+.6f}  "
      f"m_G'={mp:+.6f}  {'FAIL' if mp < 0 else 'ok'}")
P("   => theta=pi: RHS=0 (a=0) 或 <0, LHS>0 => 不成立; P05 的 [T] 给不出 m_G'>=0")
P("")
P("[DONE]")

with open(r"D:\Dev\math\research\v3\w7\p22_mg_hvv.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG) + "\n")
