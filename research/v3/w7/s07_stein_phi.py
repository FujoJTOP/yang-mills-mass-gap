"""s07_stein_phi.py -- S07: 把 q15 的 Stein 恒等式推广到带 Phi 的测度 V=S+Phi.

模型: 单 plaquette SU(2) 平直图 (p21/p22 同一模型):
    X in su(2)=R^3, rho=|X|=2*psi, u=xi_3/rho, a=xi_3 (沿 v=e_3 的坐标), 边界 holonomy theta
    S = beta*(1 - cos(th/2)cos(rho/2) + sin(th/2)sin(rho/2)u)
    dxi = 16 pi psi^2 dpsi du ; dHaar = J dxi, J = sin^2(psi)/(4 psi^2), Phi = -log J
    => dmu = e^{-S} dHaar ~ sin^2(psi) e^{-S} dpsi du          (全测度)
    => 切片 a=const: dmu_a ~ sin^2(psi)/psi e^{-S} dpsi
两条独立求积: (i) 全测度 2D Gauss-Legendre (psi,u); (ii) a-切片 1D GL + a-marginal.

题设:
 [S1] 0 = int d_v(f e^{-V}) dxi  =>  E[d_v f] = E[f d_v V],  f = V' = G+Phi'
      => E[H+Phi''] = E[(G+Phi')^2]                                  (main.tex eq:stein2)
 [S2] 展开 = E[G^2] + 2E[G Phi'] + E[Phi'^2] + 修正项 ; 修正项 =?
 [S3] E[kappa]-Var(F) = ?  (kappa=H+Phi'', F=G+Phi'); cross-Stein(f=1): E[F]=0
 [S4] 条件化 G''(a) = E[kappa|a]-Var(F|a); 总方差: E_a[G''] = Var(m_F)+E[m_F]^2
 [S5] 边界项/衰减: psi->pi 时 J->0 且 |f| e^{-V} -> 0
"""
import math
import numpy as np
from numpy.polynomial.legendre import leggauss

B = 2.5
PI = math.pi


def kern(rho, u, beta, ct, st):
    s2, c2 = np.sin(rho / 2.0), np.cos(rho / 2.0)
    S = beta * (1.0 - ct * c2 + st * s2 * u)
    S_r = beta * (ct * s2 / 2.0 + st * c2 * u / 2.0)
    S_rr = beta * (ct * c2 / 4.0 - st * s2 * u / 4.0)
    S_u = beta * st * s2
    S_ur = beta * st * c2 / 2.0
    G = S_r * u + S_u * (1 - u * u) / rho
    H = (S_rr * u * u + S_ur * u * (1 - u * u) / rho + S_r * (1 - u * u) / rho
         + S_ur * u * (1 - u * u) / rho - 3.0 * S_u * u * (1 - u * u) / (rho * rho))
    Ph_r = -c2 / s2 + 2.0 / rho
    Ph_rr = 0.5 / (s2 * s2) - 2.0 / (rho * rho)
    P1 = Ph_r * u
    P2 = Ph_rr * u * u + Ph_r * (1 - u * u) / rho
    return S, G, H, P1, P2


def full(beta, th, N=2400, M=700):
    """全测度: (psi,u) 上的 2D Gauss-Legendre, 权重 sin^2(psi) e^{-S}."""
    ct, st = math.cos(th / 2.0), math.sin(th / 2.0)
    x, wx = leggauss(N)
    y, wy = leggauss(M)
    psi = 0.5 * PI * (x + 1.0)
    wpsi = 0.5 * PI * wx
    PS, UU = np.meshgrid(psi, y, indexing="ij")
    W = np.outer(wpsi, wy)
    S, G, H, P1, P2 = kern(2.0 * PS, UU, beta, ct, st)
    w = W * np.sin(PS) ** 2 * np.exp(-S)
    Z = w.sum()
    E = lambda q: (w * q).sum() / Z
    F = G + P1
    d = dict(Z=Z, EG=E(G), EG2=E(G * G), EH=E(H), Ep1=E(P1), Ep1s=E(P1 * P1),
             Ep2=E(P2), EGp1=E(G * P1), EF=E(F), EF2=E(F * F))
    d["Ek"] = d["EH"] + d["Ep2"]
    d["VarG"] = d["EG2"] - d["EG"] ** 2
    d["VarP1"] = d["Ep1s"] - d["Ep1"] ** 2
    d["Cov"] = d["EGp1"] - d["EG"] * d["Ep1"]
    d["VarF"] = d["EF2"] - d["EF"] ** 2
    d["Wgeo"] = d["Ep2"] - d["VarP1"]
    return d


def slices(a, beta, th, npsi=400):
    """a-切片条件矩: 权重 sin^2(psi)/psi e^{-S} dpsi, psi in (|a|/2, pi)."""
    ct, st = math.cos(th / 2.0), math.sin(th / 2.0)
    a = np.asarray(a, dtype=float)[:, None]
    x, wx = leggauss(npsi)
    t = 0.5 * (x + 1.0)
    wt = 0.5 * wx
    lo = np.abs(a) / 2.0
    hi = PI - 1e-15
    psi = lo + (hi - lo) * t[None, :]
    dpsi = (hi - lo) * wt[None, :]
    rho = 2.0 * psi
    u = np.divide(a, rho, out=np.zeros_like(rho), where=rho > 0)
    S, G, H, P1, P2 = kern(rho, u, beta, ct, st)
    w = np.sin(psi) ** 2 / psi * np.exp(-S) * dpsi
    Z = w.sum(axis=1)
    E = lambda q: (w * q).sum(axis=1) / Z
    F = G + P1
    d = dict(Z=Z, EG=E(G), EH=E(H), Ep1=E(P1), Ep2=E(P2))
    d["EF"] = d["EG"] + d["Ep1"]
    d["EF2"] = E(F * F)
    d["Ek"] = d["EH"] + d["Ep2"]
    d["VarF"] = d["EF2"] - d["EF"] ** 2
    d["Gpp"] = d["Ek"] - d["VarF"]
    return d


LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


P("=" * 118)
P("s07_stein_phi.log -- S07: 带 Phi 的 Stein 恒等式 (q15 -> V=S+Phi)   beta=2.5, theta=pi/2 主")
P("模型: 单 plaquette SU(2), 平直图 su(2)=R^3, 切片 xi_3=a (<v,X>, v=e_3), V=S+Phi, Phi=-log J")
P("求积: (i) 全测度 2D Gauss-Legendre (psi,u); (ii) a-切片 1D GL + a-marginal (a=2pi sin phi)")
P("=" * 118)
P("")

# ---------------- [S0] 锚: a=0 切片 E[H|0] = P05 闭式 ----------------
P("[S0] 锚 (与 P05/p22 同一对象): a=0 切片 E[H_vv|0] vs P05 的 beta*(c/2)*G(B)")
d0 = slices([0.0], B, 0.0)
P(f"   theta=0: 本脚本 E[H|0] = {d0['EH'][0]:.9f}   P05 值 0.553060561   "
  f"差 {d0['EH'][0]-0.553060561:+.2e}")
for th in (0.0, PI / 2, PI):
    dd = slices([0.0], B, th)
    P(f"   theta={th/PI:.3f}pi: E[H|0]={dd['EH'][0]:+.6f}  E[Phi''|0]={dd['Ep2'][0]:+.6f}"
      f"  E[F|0]={dd['EF'][0]:+.2e}  G''(0)={dd['Gpp'][0]:+.6f}")
P("")

# ---------------- [S1] 主恒等式 ----------------
P("[S1] 全测度 Stein:  E[H+Phi''] ?= E[F^2] = Var(F) + E[F]^2 ,  F=G+Phi'")
P("   theta      E[H]       E[Phi'']    E[H+Phi'']   E[F^2]       Var(F)       E[F]       残差")
for th in (0.0, PI / 2, 0.7 * PI, PI):
    d = full(B, th)
    P(f"   {th/PI:5.3f}pi  {d['EH']:+.6f}  {d['Ep2']:+.6f}  {d['Ek']:+.9f}  {d['EF2']:+.9f}  "
      f"{d['VarF']:+.9f}  {d['EF']:+.2e}  {d['Ek']-d['EF2']:+.1e}")
P("   => E[H+Phi''] = E[(G+Phi')^2] 全 theta 成立 (残差=机器精度) [T](代数)+[C](数值)")
P("   => 且 cross-Stein(f=1): E[F]=E[G]+E[Phi']=E[d_vV]=0  =>  E[H+Phi''] = Var(F) (逐字 eq:stein2)")
P("")

# ---------------- [S2] 展开与修正项 ----------------
P("[S2] 展开 (theta=pi/2): E[H+Phi''] = E[G^2] + 2E[G Phi'] + E[Phi'^2] + 修正项")
d = full(B, PI / 2)
three = d["EG2"] + 2 * d["EGp1"] + d["Ep1s"]
P(f"   E[G^2]           = {d['EG2']:+.9f}")
P(f"   2E[G Phi']       = {2*d['EGp1']:+.9f}")
P(f"   E[Phi'^2]        = {d['Ep1s']:+.9f}")
P(f"   三项和            = {three:+.9f}")
P(f"   E[H+Phi'']       = {d['Ek']:+.9f}")
P(f"   修正项 = 差       = {d['Ek']-three:+.3e}")
P(f"   -E[F]^2          = {-d['EF']**2:+.3e}   (代数恒等: 修正项 == -E[F]^2)")
P("   => 修正项 = 0 (机器精度). Phi 的 v-导数结构**不产生**额外项: 平直图里 v 是常系数场,")
P("      散度算子 d_v 与 Hessian 的 d_v^2 是同一个算子 (无 transport/Christoffel 项).")
P("      只有在非无散帧 (div_mu v != 0) 才会出现 修正项 = -E[(div_mu v) d_vV].")
P("")

# ---------------- [S3] 判: 平方是不是真的非负内容 ----------------
P("[S3] 判: E[kappa]-Var(F) = E[F]^2 ?  (kappa=H+Phi'', F=G+Phi')")
for th in (0.0, PI / 2, 0.7 * PI, PI):
    d = full(B, th)
    P(f"   {th/PI:5.3f}pi: E[kappa]-Var(F) = {d['Ek']-d['VarF']:+.3e}"
      f"   E[F]^2 = {d['EF']**2:.2e}   E[kappa]=Var(F)={d['Ek']:.9f}")
P("   => 平方项 = E[F]^2, 而 E[F]=E[d_vV]=0 (score 零均值) 是**恒等**的,")
P("      不是靠对称性; 故 平方 == 0, Stein 只给**等式** E[kappa]=Var(F), 没有余量.")
P("")
P("[S3b] 对照: 无 Phi 的 q15 式 (E[H]-Var(G)=E[G]^2) 在 V-测度下**不成立**")
for th in (0.0, PI / 2, PI):
    d = full(B, th)
    P(f"   {th/PI:5.3f}pi: E[H]-Var(G) = {d['EH']-d['VarG']:+.6f}   E[G]^2 = {d['EG']**2:+.6f}"
      f"   差 {d['EH']-d['VarG']-d['EG']**2:+.6f}   W_geo=E[Phi'']-Var(Phi')={d['Wgeo']:+.6f}"
      f"   Cov(G,Phi')={d['Cov']:+.6f}")
P("   => V-测度下正确的关系是 E[H]-Var(G) = W_geo = Cov(G,Phi') (cross-Stein), 不是 E[G]^2 [T]+[C]")
P("")

# ---------------- [S4] 条件化层 ----------------
P("[S4] 条件化 (a-切片): G''(a)=E[kappa|a]-Var(F|a);  总方差分解 => E_a[G''] = Var_a(m_F)+E[m_F]^2")
x, wx = leggauss(121)
phi = 0.5 * PI * x
aw = 0.5 * PI * wx
a_grid = 2.0 * PI * np.sin(phi)
amarg_w = np.cos(phi) * aw                      # da = 2pi cos(phi) dphi
ds = slices(a_grid, B, PI / 2)
wei = amarg_w * ds["Z"]
wei = wei / wei.sum()
mean_mF = float((wei * ds["EF"]).sum())
var_mF = float((wei * (ds["EF"] - mean_mF) ** 2).sum())
avg_Gpp = float((wei * ds["Gpp"]).sum())
P("   a       m_F(a)=E[F|a]   E[kappa|a]    Var(F|a)    G''(a)      rho(a)权重")
for i in (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120):
    P(f"  {a_grid[i]:+6.3f}   {ds['EF'][i]:+9.5f}   {ds['Ek'][i]:+9.5f}   {ds['VarF'][i]:+.5f}"
      f"    {ds['Gpp'][i]:+9.5f}   {wei[i]:.4f}")
imin, imax = int(np.argmin(ds["Gpp"])), int(np.argmax(ds["Gpp"]))
P(f"   min G''(a) = {ds['Gpp'][imin]:+.6f} (a={a_grid[imin]:+.3f});  "
  f"max = {ds['Gpp'][imax]:+.6f} (a={a_grid[imax]:+.3f})")
P(f"   E_a[G''] = {avg_Gpp:+.6f}    Var_a(m_F) = {var_mF:+.6f}    E[m_F]^2 = {mean_mF**2:+.2e}"
  f"    差 {avg_Gpp-var_mF-mean_mF**2:+.1e}")
P("   => 条件层只有 G''(a) = m_F'(a) (p22/P09 (*)), 无平方; 但总方差分解给")
P("      E_a[G''] = Var_a(m_F) + E[m_F]^2 = Var_a(E[F|a]) >= 0  (**只 a-平均后**非负)")
P("   => 逐 a 的 G''(a) 有正有负 (见 min) => 点态 G''>=0 (H1) 不由此得到")
P("")

# ---------------- [S5] 前提: 边界项 ----------------
P("[S5] 前提核验: 分部积分的边界项 = |f| e^{-V} (psi->pi, 即 rho->2pi 的 chart 边)")
P("   psi       rho        J=sin^2(psi)/(4psi^2)   Phi'        |f| e^{-V}*const")
for psi in (PI - 1e-2, PI - 1e-3, PI - 1e-4, PI - 1e-6):
    rho = 2.0 * psi
    J = math.sin(psi) ** 2 / (4.0 * psi ** 2)
    S, G, H, P1, P2 = kern(rho, 0.3, B, math.cos(PI / 4), math.sin(PI / 4))
    P(f"   {psi:.6f}  {rho:8.6f}   {float(J):.3e}          {float(P1):+.3e}   "
      f"{float((abs(G)+abs(P1))*math.exp(-S)*J*16*PI*psi**2*1e-3):.3e}")
P("   => J ~ sin^2(psi) -> 0 二次, Phi' ~ 1/sin(psi) 一次 => |f| e^{-V} -> 0 : 边界项=0 [T](模型内)")
P("      psi->0 (rho->0): J->1/4, Phi'->0 (Ph_r = -cot(psi)+1/psi = O(psi)) => 无奇异 [T]")
P("      另: f=G+Phi' 在 rho->2pi 上按 1/sin(psi) 发散但被 J 杀掉, 故 V' 的边界项 0.")
P("")

# ---------------- [S6] 点态 G''(a) 窗口扫描 ----------------
P("[S6] 点态 G''(a) 的符号 (窗口 |a|<=5.5; |a|->2pi 是 cut-locus 伪影: Phi 发散, 权重->0)")
P("   theta     min G''      位置 a      max G''      位置 a     #负点/总点数")
for th in (PI / 2, PI, 0.75 * PI):
    ag = np.linspace(-5.5, 5.5, 221)
    dg = slices(ag, B, th)
    g = dg["Gpp"]
    P(f"   {th/PI:5.3f}pi  {g.min():+10.6f}  {ag[int(np.argmin(g))]:+8.3f}   "
      f"{g.max():+10.6f}  {ag[int(np.argmax(g))]:+8.3f}     {int((g<0).sum())}/{g.size}")
P("   => theta=pi/2: 窗口内 G''(a) 全正 (最小 +4e-3) => 与 H1 数值相容, 但**不是**由 Stein 推出")
P("")
P("[DONE]")

with open(r"D:\Dev\math\research\v3\w7\s07_stein_phi.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG) + "\n")
