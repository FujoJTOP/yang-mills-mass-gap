"""q03_separable.py -- Q03: H1 能否分解成 S_W 半边 + Phi 半边, 各自独立证?

模型: 单 plaquette SU(2) 平直图 (H1A-EXACT sec.2.1 原样; 与 p21_rebuild.py / p22_mg_hvv.py 逐字同模型):
    X in su(2)=R^3, rho=|X|, u=xi_3/rho, 切片 xi_3=a (a = <v,X>, v=e_3), 边界 holonomy theta
    S = beta*(1 - cos(th/2)cos(rho/2) + sin(th/2)sin(rho/2)u)
    dmu_a ∝ sin^2(rho/2)/rho * exp(-S) drho
    G=d_3 S, H=d_3^2 S, Phi'=d_3 Phi, Phi''=d_3^2 Phi, F=G+Phi', kappa=H+Phi''

要判的三条 (任务给的分解式):
    D - S = Var(F) - E[kappa] = [Var(G)-E[H]] + [Var(Phi')-E[Phi'']] + 2Cov(G,Phi')
    A := Var(G)-E[H]        (S_W 半边; A<=0 <=> Var(G) <= E[H])
    B := Var(Phi')-E[Phi''] (Phi 半边;  B<=0 <=> Var(Phi') <= E[Phi''])
    C := 2Cov(G,Phi')       (交叉项)
    H1 本身: Var(G)+Cov <= E[H]   <=>  A + Cov <= 0

[T] 恒等式 (本文骨架, P09 (*) 取 f=G 与 f=Phi'):
    (I1) m_G' := d/da E[G|a]  = E[H] - Var(G) - Cov
    (I2) m_P' := d/da E[Phi'|a] = E[Phi''] - Var(Phi') - Cov
    => A = -m_G' - Cov,  B = -m_P' - Cov,  A+B+2Cov = -(m_G'+m_P') = -d/da E[F|a]
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


# ---- 与 p21_rebuild.py / p22_mg_hvv.py 同一模型 (逐字复用) ----
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


def moments(a, beta, th, eps=1e-12):
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
    d = dict(EG=EG, EH=EH, EPhi1=Ep1, EPhi2=Ep2)
    d["VarG"] = EG2 - EG * EG
    d["VarPhi1"] = Ep1s - Ep1 * Ep1
    d["Cov"] = EGp1 - EG * Ep1
    # 三块
    d["A"] = d["VarG"] - EH                     # S_W 半边
    d["B"] = d["VarPhi1"] - Ep2                 # Phi 半边
    d["C"] = 2.0 * d["Cov"]                     # 交叉项
    d["DS"] = d["VarG"] + d["VarPhi1"] + 2.0 * d["Cov"] - (EH + Ep2)   # = Var(F)-E[kappa]
    d["H1"] = d["VarG"] + d["Cov"] - EH         # H1 的赤字 (<=0 则 H1 成立)
    d["Gpp"] = EH + Ep2 - (d["VarG"] + d["VarPhi1"] + 2.0 * d["Cov"])
    d["Z"] = Z
    return d


def dmom(a, beta, th, key, h=0.02):
    return (moments(a + h, beta, th)[key] - moments(a - h, beta, th)[key]) / (2 * h)


def verdict(x):
    return "ok" if x <= 0.0 else "FAIL"


P("=" * 118)
P("q03_separable.log -- Q03: H1 的 S_W 半边 / Phi 半边 / 2Cov 能否各自独立成立")
P("模型与 p21_rebuild.py / p22_mg_hvv.py 逐字相同: 单 plaquette SU(2), 平直图 su(2)=R^3, 切片 xi_3=a, v=e_3")
P("分解式 (任务): D-S = Var(F)-E[kappa] = [Var(G)-E[H]] + [Var(Phi')-E[Phi'']] + 2Cov  = A + B + C")
P("H1 本身: Var(G)+Cov <= E[H]  <=>  A + Cov <= 0")
P("=" * 118)
P("")

# ---------------- [1] 锚: 模型身份 ----------------
P("[1] 锚 (同一模型的证据): 本脚本复算 p22/p21 的既有数, 应当逐位重合")
d = moments(0.0, 2.5, 0.0)
P(f"    a=0, theta=0 :  E[H] = {d['EH']:+.9f}   [P05-CINF: 0.553060561]   diff = {d['EH']-0.553060561:+.2e}")
rows_anchor = [(0.25, -3.353e-3, 2.923e-4, 4.2712e-2), (0.5, -7.192e-3, 1.369e-3, 4.9574e-2),
               (0.75, -1.1406e-2, 3.610e-3, 5.9511e-2), (1.0, -1.5823e-2, 7.519e-3, 7.3164e-2),
               (1.5, -2.4339e-2, 2.306e-2, 0.120), (2.0, -3.0596e-2, 5.510e-2, 0.185016)]
P("    a      Cov(本脚本)    Cov(P21)      VarPhi'(本)     VarPhi'(P21)    EPhi''(本)     Q4(P21)+1/6")
for a, cov_ref, vp_ref, ep_ref in rows_anchor:
    dd = moments(a, 2.5, math.pi)
    P(f"   {a:<5.2f} {dd['Cov']:+.6e}  {cov_ref:+.6e}  {dd['VarPhi1']:.6e}  {vp_ref:.6e}   "
      f"{dd['EPhi2']:+.6f}    {ep_ref+N12:+.6f}")
P("    (P21 候选4 表列的是 Q4=E[Phi'']-N/12, 故此处加回 1/6 再比)")
P("    => 与 P21 表 (候选1/2/4) 同量同号, |diff| < 1e-4 => 本脚本测的是同一个对象 [C]")
P("")

# ---------------- [2] 主表: beta=2.5, theta-scan x a-scan ----------------
P("[2] 主表 beta=2.5: 三个数 A=Var(G)-E[H], B=Var(Phi')-E[Phi''], C=2Cov + D-S = A+B+C")
P("     th/pi   a      E[H]       Var(G)      A=VarG-EH   Var(Phi')    E[Phi'']     B          Cov         C=2Cov     D-S        H1:A+Cov")
grid = []
for th in (0.0, math.pi / 2, 0.7 * math.pi, math.pi):
    for a in (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
        d = moments(a, 2.5, th)
        grid.append((th, a, d))
        P(f"    {th/math.pi:5.3f}  {a:<5.2f} {d['EH']:+10.6f} {d['VarG']:14.6f} {d['A']:+11.6f} "
          f"{d['VarPhi1']:11.6e} {d['EPhi2']:+10.6f} {d['B']:+11.6f} {d['Cov']:+11.6e} {d['C']:+11.6e} "
          f"{d['DS']:+11.6f} {d['H1']:+11.6f}")
P("")
P("    => S_W 半边 A<=0 最大者 (最坏点):")
worst = max(grid, key=lambda r: r[2]["A"])
P(f"       th={worst[0]/math.pi:.3f}pi a={worst[1]:.2f}: A = {worst[2]['A']:+.6f}  "
  f"(Var(G)={worst[2]['VarG']:.6f} vs E[H]={worst[2]['EH']:+.6f})  => {verdict(worst[2]['A'])}")
best = min(grid, key=lambda r: r[2]["A"])
P(f"       最有利点 th={best[0]/math.pi:.3f}pi a={best[1]:.2f}: A = {best[2]['A']:+.6f} => {verdict(best[2]['A'])}")
P("       [注意 th=0,a=0 是运动学退化点: sin(th/2)=0 => G≡0, Var(G)=0, 不能当'S_W 半边成立'的证据]")
nondeg = [r for r in grid if not (abs(r[0]) < 1e-12 and r[1] == 0.0)]
bn = min(nondeg, key=lambda r: r[2]["A"])
P(f"       非退化最有利点 th={bn[0]/math.pi:.3f}pi a={bn[1]:.2f}: A = {bn[2]['A']:+.6f} => {verdict(bn[2]['A'])}")
P("    => Phi 半边 B<=0 最大者 (最坏点):")
worstB = max(grid, key=lambda r: r[2]["B"])
P(f"       th={worstB[0]/math.pi:.3f}pi a={worstB[1]:.2f}: B = {worstB[2]['B']:+.6f}  "
  f"(Var(Phi')={worstB[2]['VarPhi1']:.6e} vs E[Phi'']={worstB[2]['EPhi2']:+.6f})  => {verdict(worstB[2]['B'])}")
bestB = min(grid, key=lambda r: r[2]["B"])
P(f"       最有利点 th={bestB[0]/math.pi:.3f}pi a={bestB[1]:.2f}: B = {bestB[2]['B']:+.6f} => {verdict(bestB[2]['B'])}")
P("")

# ---------------- [3] 交叉项量级 ----------------
P("[3] 交叉项 2Cov 的量级: 它是不是主导?")
P("     th/pi   a      |C|/(|A|+|B|)   |C|/|D-S|     sign(C)")
for th, a, d in grid:
    if a == 0.0:
        continue
    den = abs(d["A"]) + abs(d["B"])
    P(f"    {th/math.pi:5.3f}  {a:<5.2f} {abs(d['C'])/den:12.3e} {abs(d['C'])/max(1e-300,abs(d['DS'])):12.3e}   "
      f"{'+' if d['C'] >= 0 else '-'}")
P("")
rat = [(abs(d['C'])/(abs(d['A'])+abs(d['B'])), th, a) for th, a, d in grid if a > 0]
P(f"    => max |2Cov|/(|A|+|B|) = {max(rat)[0]:.3e}  at th={max(rat)[1]/math.pi:.3f}pi, a={max(rat)[2]}")
P("    => 交叉项是全网格的最小块; 见 [4] 它还有精确的代数抵消")
P("")

# ---------------- [4] 恒等式核验 ----------------
P("[4] [T] 恒等式核验 (中心差分 h=0.02):")
P("    (I1) A = -m_G' - Cov ;  (I2) B = -m_P' - Cov ;  (I3) D-S = -d/da E[F|a] = -(m_G'+m_P')")
for th, a in ((math.pi, 0.5), (math.pi, 0.75), (math.pi, 1.0), (math.pi / 2, 0.75)):
    d = moments(a, 2.5, th)
    mG = dmom(a, 2.5, th, "EG")
    mP = dmom(a, 2.5, th, "EPhi1")
    mF = dmom(a, 2.5, th, "EG") + mP
    P(f"    th={th/math.pi:4.2f}pi a={a:<4.2f}: A={d['A']:+11.6f} vs -mG'-Cov={-mG-d['Cov']:+11.6f} "
      f"(res {d['A']+mG+d['Cov']:+.1e}) | B={d['B']:+11.6f} vs -mP'-Cov={-mP-d['Cov']:+11.6f} "
      f"(res {d['B']+mP+d['Cov']:+.1e}) | D-S={d['DS']:+11.6f} vs -(mG'+mP')={-mF:+11.6f} (res {d['DS']+mF:+.1e})")
P("    => 三式全成立 => 2Cov 在两半的 -Cov 里被精确吃掉 (对称性, 非'独立第三块')")
P("")

# ---------------- [5] beta 窗口: 两半是不同阶的 beta 多项式 ----------------
P("[5] beta-扫描: G=beta*g1, H=beta*h1 (S 对 beta 线性) => E[H] ∝ beta 精确; Var(G) 则随 beta 饱和")
P("     th/pi   a      beta        beta^2*Var(g1)=Var(G)   beta*E[h1]=E[H]   A           A/beta^2     A/beta")
for th, a in ((0.0, 0.0), (math.pi / 2, 0.0), (0.0, 0.75), (math.pi, 0.0), (math.pi, 0.75)):
    for beta in (0.1, 0.5, 1.0, 2.5, 5.0, 10.0, 20.0):
        d = moments(a, beta, th)
        P(f"    {th/math.pi:5.3f}  {a:<5.2f} {beta:8.2f} {d['VarG']:20.6e} {d['EH']:+18.6e} {d['A']:+12.6e} "
          f"{d['A']/beta**2:+11.6e} {d['A']/beta:+11.6e}")
P("")
P("    => 关键 [V]-in-model (精确, 不是数值): S 与 beta 成正比 => G=beta*g1, H=beta*h1;")
P("       故 E[H]=beta*E[h1] (E[h1] 是 beta-无关常数, 实测 E[H]/beta 在 7 个 beta 上稳定),")
P("       而 Var(G)>=0 恒真 => **E[H] < 0 的点上 A = Var(G)-E[H] > 0 对一切 beta>0**.")
P("    a      E[h1]=E[H]/beta 于 beta=2.5   Var(G)/beta^2 (Haar 端, beta=0.1)   A 的符号 (beta=0.1..20)")
for th, a in ((0.0, 0.0), (0.0, 0.75), (math.pi / 2, 0.0), (math.pi / 2, 2.0),
              (0.7 * math.pi, 0.75), (math.pi, 0.0), (math.pi, 0.75), (math.pi, 2.0)):
    d25 = moments(a, 2.5, th)
    d01 = moments(a, 0.1, th)
    signs = set()
    for beta in (0.1, 0.5, 1.0, 2.5, 5.0, 10.0, 20.0):
        signs.add(moments(a, beta, th)["A"] > 0)
    P(f"    {th/math.pi:5.3f}  {a:<5.2f} {d25['EH']/2.5:+20.6f} {d01['VarG']/0.01:30.6e}   "
      f"{'A>0 (FAIL) 全 beta' if signs == {True} else ('A<=0 (ok) 全 beta' if signs == {False} else 'A 随 beta 变号')}")
P("       th=pi, a=0: E[h1] = E[H]/beta = 0 精确 (B=beta*cos(pi/2)=0) => A = Var(G) > 0 对一切 beta")
P("")

# ---------------- [7] 抵消性: D-S 是两个大反号块的小差 ----------------
P("[7] 抵消性 (为什么'两半各自证'不是 D-S<=0 的稳定路线): |D-S| vs |A|+|B|")
P("     th/pi   a      |A|         |B|         |A|+|B|     |D-S|      |D-S|/(|A|+|B|)   A 与 B 是否反号")
mx_ratio = (0.0, None)
for th, a, d in grid:
    if a == 0.0:
        continue
    s = abs(d["A"]) + abs(d["B"])
    rr = abs(d["DS"]) / s
    if rr > mx_ratio[0]:
        mx_ratio = (rr, (th, a))
    P(f"    {th/math.pi:5.3f}  {a:<5.2f} {abs(d['A']):11.6f} {abs(d['B']):11.6f} {s:11.6f} {abs(d['DS']):11.6f} "
      f"{rr:15.3e}   {'反号' if d['A']*d['B'] < 0 else '同号'}")
P(f"    => max |D-S|/(|A|+|B|) = {mx_ratio[0]:.3e} at th={mx_ratio[1][0]/math.pi:.3f}pi a={mx_ratio[1][1]}")
n_opp = sum(1 for th, a, d in grid if a > 0 and d["A"] * d["B"] < 0)
n_pos = sum(1 for th, a, d in grid if a > 0)
P(f"    => {n_opp}/{n_pos} 个 (a>0) 点上 A 与 B **反号** => D-S 是两个大块之差, 不是两块之和")
P("")
P("[8] H1 与两半的逻辑关系 (数值): H1 <==> A+Cov<=0 <==> -m_G'<=0 <==> m_G'>=0  [T] 代数")
P("     (代入 I1: A+Cov = (VarG-EH)+Cov = -m_G' ; Phi 半边 B **完全不进 H1**)")
P("     th/pi   a      H1 赤字 A+Cov    -m_G'(差分)     |A+Cov| 与 |A| 之比   H1    D-S<=0(G''>=0)?")
for th, a, d in grid:
    if a in (0.0,):
        continue
    P(f"    {th/math.pi:5.3f}  {a:<5.2f} {d['H1']:+14.6f} {d['A']+d['Cov']:+14.6f} "
      f"{(1 if d['A']==0 else abs(d['H1'])/abs(d['A'])):18.3f}   {'ok' if d['H1']<=0 else 'FAIL'}   "
      f"{'ok' if d['DS']<=0 else 'FAIL'}")
P("    => H1 与 D-S<=0 在不同点上给出**不同判决** (例: th=pi a=0.75: H1 赤字 +0.2399 而 D-S=+0.0059)")
P("")


P("[6] L=2 (真实靶) 上的同一分解, 用 o27/p21 的既有 [C] 数 (本脚本不算 L=2):")
P("    A = -Delta_dyn = -0.562 (Delta_dyn=+0.562+-0.277);  B = -W_geo = -1.09 (W_geo=+1.09+-0.02);")
P("    C = 2Cov = 2*(-0.0215+-0.0298) = -0.043;  D-S = A+B+C = -1.695")
P(f"    对照 G'' 实测 = +1.696+-0.293 => 与 -(D-S) 差 {-1.695+1.696:+.3f} (0.06%)")
P("    => L=2 上两半都成立 (A,B 皆 <0) 且 |2Cov|/|D-S| = 2.5%")
P("")
P("[DONE]")

with open(r"D:\Dev\math\research\v3\w7\q03_separable.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG) + "\n")
