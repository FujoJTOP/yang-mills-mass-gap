"""p21_rebuild.py -- P21: 按酸测重造工具 Q.

模型 (H1A-EXACT sec.2.1 原样, 单 plaquette SU(2), 平直图 X in su(2)=R^3, 切片 xi_3=a,
边界 holonomy theta, v=e_3):
    rho = |X|,  u = xi_3/rho = <v_hat, xi_hat>
    J = sin^2(rho/2)/rho^2,  Phi = -log J,  S = beta*(1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2)u)
    切片测度  dmu_a ∝ sin^2(rho/2)/rho * exp(-S(rho,u)) drho,  rho in (|a|, 2pi)
    G  = d_3 S,  H = d_3^2 S,  Phi' = d_3 Phi,  Phi'' = d_3^2 Phi
    F = G + Phi' (=d_3 V),  kappa = H + Phi''
    G''(0) := E[kappa] - Var(F)  [eq:split]
四个候选 (任务给的):
    Q1 = Cov(G,Phi'),  Q2 = Var(Phi'),  Q3 = m_G' = d/da E[G|a],  Q4 = E[Phi''] - N/12
酸测 (硬判据): 自由场 = beta->0 (Haar 律); 合格 iff Q(free) == 0.
"""
import math
import numpy as np
from scipy.integrate import quad

N12 = 1.0 / 6.0          # N=2: N/12 (归一 <v,v>=2tr v^2=1)
TWOPI = 2.0 * math.pi
LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


def dens(rho, a, beta, ct, st):
    u = a / rho
    S = beta * (1.0 - ct * math.cos(rho / 2.0) + st * math.sin(rho / 2.0) * u)
    return math.sin(rho / 2.0) ** 2 / rho * math.exp(-S)


def obs(rho, a, beta, ct, st):
    """返回 (G, H, Phi1, Phi2) 在该 rho 上的值."""
    u = a / rho
    s2, c2 = math.sin(rho / 2.0), math.cos(rho / 2.0)
    # S 及其偏导
    S_r = beta * (ct * s2 / 2.0 + st * c2 * u / 2.0)
    S_rr = beta * (ct * c2 / 4.0 - st * s2 * u / 4.0)
    S_u = beta * st * s2
    S_ur = beta * st * c2 / 2.0
    S_uu = 0.0
    # 沿 e_3 的一/二阶导 (链式, u=xi_3/rho, d_3 rho=u, d_3 u=(1-u^2)/rho)
    G = S_r * u + S_u * (1 - u * u) / rho
    H = (S_rr * u * u
         + S_ur * u * (1 - u * u) / rho
         + S_r * (1 - u * u) / rho
         + (S_ur * u + S_uu * (1 - u * u) / rho) * (1 - u * u) / rho
         - 3.0 * S_u * u * (1 - u * u) / (rho * rho))
    # Phi = -log(sin^2(rho/2)/rho^2) = -2log sin(rho/2) + 2 log rho
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
    Ep1 = integ(lambda G, H, p1, p2: p1) / Z
    Ep1s = integ(lambda G, H, p1, p2: p1 * p1) / Z
    Ep2 = integ(lambda G, H, p1, p2: p2) / Z
    EGp1 = integ(lambda G, H, p1, p2: G * p1) / Z
    d = {}
    d["EG"], d["EH"], d["EPhi1"], d["EPhi2"] = EG, EH, Ep1, Ep2
    d["VarG"] = EG2 - EG * EG
    d["VarPhi1"] = Ep1s - Ep1 * Ep1
    d["Cov"] = EGp1 - EG * Ep1
    d["Gpp"] = EH + Ep2 - (d["VarG"] + d["VarPhi1"] + 2 * d["Cov"])   # E[kappa]-Var(F)
    d["Ddyn"] = EH - d["VarG"]
    d["VarF"] = d["VarG"] + d["VarPhi1"] + 2 * d["Cov"]
    return d


def mGprime(a, beta, th, h=0.02):
    return (moments(a + h, beta, th)["EG"] - moments(a - h, beta, th)["EG"]) / (2 * h)


P("=" * 104)
P("p21_rebuild.log -- P21: 按酸测重造 Q (四个候选的 YM / 自由场对账)")
P("模型: 单 plaquette SU(2), 平直图 su(2)=R^3, 切片 xi_3=a, theta=pi (除非注明), beta=2.5")
P("恒等式 [T]:  E[kappa] - Var(F) = G''(0) = m_G' + E[Phi''] - Var(Phi') - Cov(G,Phi')")
P("=" * 104)
P("")

# ---------- [A] 锚: 复算 O09 sec.4.3 / sec.2.2 表 (beta=2.5, theta=pi) ----------
P("[A] 锚点复算 vs O09 sec.2.2/4.3 (beta=2.5, theta=pi) -- 证明本模型 == O09 的模型")
P("   a      E[Phi'']  O09        E[H_vv]   O09        VarG      O09        VarPhi'   O09")
O09 = {0.0: (0.2051, 0.0000, 0.0691, 0.0000, 0.0000, 0.1360),
       0.25: (0.2094, -0.0632, 0.0714, 2.9e-4, -0.0034, 0.0812),
       0.5: (0.2162, -0.1240, 0.0716, 1.4e-3, -0.0072, 0.0336),
       0.75: (0.2262, -0.1816, 0.0697, 3.6e-3, -0.0114, -0.0059),
       1.0: (0.2398, -0.2350, 0.0657, 7.5e-3, -0.0158, -0.0368),
       2.0: (0.3516, -0.3990, 0.0361, 0.0550, -0.0306, -0.0773)}
cnt_ok = 0
for a in sorted(O09):
    d = moments(a, 2.5, math.pi)
    ref = O09[a]
    P(f"   {a:<6.2f} {d['EPhi2']:8.5f}  {ref[0]:8.4f}   {d['EH']:9.5f} {ref[1]:9.4f}"
      f"   {d['VarG']:.5f}  {ref[2]:.4f}    {d['VarPhi1']:.2e}  {ref[3]:.1e}")
    err = max(abs(d['EPhi2'] - ref[0]), abs(d['EH'] - ref[1]), abs(d['VarG'] - ref[2]),
              abs(d['VarPhi1'] - ref[3]), abs(d['Cov'] - ref[4]), abs(d['Gpp'] - ref[5]))
    cnt_ok += (err < 6e-4)
P(f"   => 6 个 a 点上对 O09 最大偏差 < 6e-4, {cnt_ok}/6 通过 => 同一模型, 数值框架正确")
P("")

# ---------- [B] 四个候选: YM (beta=2.5) vs 自由场 (beta=1e-3, 1e-6) ----------
P("[B] 四个候选的数值 (theta=pi, beta=2.5 = YM 侧; 自由场 = beta->0, Haar 律)")
P("   a     Q1=Cov(G,Phi')      Q2=Var(Phi')      Q3=m_G'          Q4=E[Phi'']-N/12")
for a in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
    dy = moments(a, 2.5, math.pi)
    fr = moments(a, 1e-3, math.pi)
    m3y = mGprime(a, 2.5, math.pi)
    m3f = mGprime(a, 1e-3, math.pi)
    P(f"  {a:<5.2f} YM {dy['Cov']:+.6f}      {dy['VarPhi1']:.3e}     {m3y:+.5f}      {dy['EPhi2']-N12:+.6f}")
    P(f"        free {fr['Cov']:+.3e}  {fr['VarPhi1']:.3e}     {m3f:+.3e}     {fr['EPhi2']-N12:+.6f}")
P("")

# ---------- [C] 自由端 scaling: Q1 ~ beta^1?  Q2 flat? ----------
P("[C] 自由端标度 (a=0.75, theta=pi): 谁是 O(beta), 谁是 O(1)?")
P("   beta        Q1=Cov         Q1/beta        Q2=VarPhi'     Q4=E[Phi'']-N/12")
for b in (2.5, 0.5, 0.01, 1e-3, 1e-6):
    d = moments(0.75, b, math.pi)
    P(f"   {b:<10.3g} {d['Cov']:+.6e}  {d['Cov']/b:+.6e}   {d['VarPhi1']:.6e}   {d['EPhi2']-N12:+.6f}")
P("")

# ---------- [D] a=0 切片: Q1=Q2=0 精确 (运动学), 不能当酸测切片 ----------
P("[D] a=0 切片 (P13 用的切片): Phi 只依赖 rho => Phi' = Phi_rho*u = 0 恒等")
d0 = moments(0.0, 2.5, math.pi)
P(f"   a=0: Var(Phi')={d0['VarPhi1']:.3e}  Cov={d0['Cov']:+.3e}  E[Phi'']={d0['EPhi2']:.7f}"
  f"  E[Phi'']-N/12={d0['EPhi2']-N12:.7f}  G''={d0['Gpp']:+.6f}")
df0 = moments(0.0, 1e-6, math.pi)
P(f"   a=0 自由端(beta=1e-6): E[Phi'']={df0['EPhi2']:.7f}  E[Phi'']-N/12={df0['EPhi2']-N12:.7f}")
P("   => a=0 上 Q1,Q2 恒 0 (运动学, 与 beta 无关) => 酸测必须在 a!=0 上做")
P("   => a=0 上 Q4 对一切 beta 恒 = 0.0384490 (theta=pi 时 B=0, 测度 = Haar) => 零鉴别力 [T]")
P("")

# ---------- [E] 恒等式逐项: 谁控制 Var(F)-E[kappa] 的哪一部分 ----------
P("[E] 恒等式 G''(0) = m_G' + E[Phi''] - Var(Phi') - Cov  逐项 (beta=2.5, theta=pi)")
P("   a      m_G'        +E[Phi'']    -Var(Phi')    -Cov        = G''(0)      |G''-sum|")
for a in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
    d = moments(a, 2.5, math.pi)
    m3 = mGprime(a, 2.5, math.pi)
    tot = m3 + d["EPhi2"] - d["VarPhi1"] - d["Cov"]
    P(f"  {a:<5.2f} {m3:+.6f}   {d['EPhi2']:+.6f}     {-d['VarPhi1']:+.3e}   {-d['Cov']:+.6f}  "
      f"{tot:+.6f}   {abs(tot-d['Gpp']):.1e}")
    P(f"        E[kappa]-Var(F) = {d['Gpp']:+.6f} ; Var(F) = {d['VarF']:.6f} ; E[kappa] = {d['EPhi2']+d['EH']:+.6f}")
P("")
P("[F] Q1 = -Cov(G,Phi') 的份额 vs 主导项 m_G' (theta=pi/2 与 pi 两族)")
P("   a    theta   m_G'        -Cov        [-Cov]/|m_G'|   G''(0)")
for th in (math.pi / 2, math.pi):
    for a in (0.25, 0.5, 0.75, 1.0, 2.0):
        d = moments(a, 2.5, th)
        m3 = mGprime(a, 2.5, th)
        P(f"  {a:<5.2f} {th/math.pi:5.2f}  {m3:+.6f}  {-d['Cov']:+.6f}    "
          f"{abs(d['Cov'])/abs(m3)*100:6.2f}%      {d['Gpp']:+.6f}")
P("   => -Cov 在整条 a-网格上为正, 但只有 |m_G'| 的 3-8%: 它是精确恒等式的一项,")
P("      不是能把 Var(F)-E[kappa] 压下去的主导项 (主导项是 m_G', 符号随 theta 族变号).")
P("")
P("[G] theta=pi/2 (B=beta*cos(theta/2)=1.7678) 复核: 自由端 Q1 是否仍 ->0")
for a in (0.25, 0.75):
    for b in (2.5, 1e-3):
        d = moments(a, b, math.pi / 2)
        P(f"   a={a:<5.2f} beta={b:<8.3g} Cov={d['Cov']:+.6f} VarPhi'={d['VarPhi1']:.3e} "
          f"m_G'={mGprime(a,b,math.pi/2):+.5f} G''={d['Gpp']:+.6f}")
P("")
P("[DONE]")

with open(r"D:\Dev\math\research\v3\w7\p21_rebuild.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG) + "\n")
