# -*- coding: utf-8 -*-
r"""
v06_decide.py -- V06: 裁决 G''(0) = 0 (s07) vs 0.5429 (r03)

独立实现 (从 su(2) 参数化重新推导, 不抄 s07/r03 的核函数):
    X in su(2) = R^3, rho=|X|=2psi, u=X_3/rho, a=X_3
    S = beta(1 - ct*cos(psi) + st*sin(psi)*u),  ct=cos(th/2), st=sin(th/2)
    dHaar = 4 pi sin^2(psi) dpsi du          (J dxi, J=sin^2(psi)/(4psi^2), dxi=16pi psi^2 dpsi du)
    Phi = -log J = 2log(2psi) - 2log sin psi

    d_a rho = u,  d_a u = (1-u^2)/rho
    d_a^2 f = u^2 f_rr + 2u(1-u^2)/rho f_ru + (1-u^2)/rho f_r
              - 3u(1-u^2)/rho^2 f_u + (1-u^2)^2/rho^2 f_uu     (f_uu = 0 for S)

两套测度都算:
  [F] 全测度 (2D: psi,u),  w = sin^2(psi)
  [S] a=0 切片 (1D: psi), w = sin^2(psi)/psi      <- r03 的 nu, 也是 s07 的 slices(a=0)
额外: [L] S-测度 (Lebesgue dxi, 无 J) 用于核 q15 的 E[H]=E[G^2]
"""
import math
import numpy as np
from numpy.polynomial.legendre import leggauss

B = 2.5
TH = math.pi / 2
CT, ST = math.cos(TH / 2), math.sin(TH / 2)


def kern(psi, u):
    rho = 2.0 * psi
    sp, cp = np.sin(psi), np.cos(psi)
    S = B * (1.0 - CT * cp + ST * sp * u)
    Sr = B * (CT * sp / 2.0 + ST * cp * u / 2.0)
    Srr = B * (CT * cp / 4.0 - ST * sp * u / 4.0)
    Sru = B * ST * cp / 2.0
    Su = B * ST * sp
    G = u * Sr + (1.0 - u * u) / rho * Su
    H = (u * u * Srr + 2.0 * u * (1.0 - u * u) / rho * Sru
         + (1.0 - u * u) / rho * Sr - 3.0 * u * (1.0 - u * u) / rho**2 * Su)
    Pr = 2.0 / rho - cp / sp                      # Phi_rho
    Prr = -2.0 / rho**2 + 0.5 / sp**2             # Phi_rho rho
    P1 = u * Pr                                   # Phi'
    P2 = u * u * Prr + (1.0 - u * u) / rho * Pr   # Phi''
    return S, G, H, P1, P2


def moments(w, G, H, P1, P2):
    Z = w.sum()
    E = lambda q: (w * q).sum() / Z
    d = {}
    d["EH"] = E(H)
    d["EP2"] = E(P2)
    d["EG"] = E(G)
    d["EG2"] = E(G * G)
    d["EP1"] = E(P1)
    d["EP1s"] = E(P1 * P1)
    d["EGP1"] = E(G * P1)
    F = G + P1
    d["EF"] = E(F)
    d["EF2"] = E(F * F)
    d["VarG"] = d["EG2"] - d["EG"] ** 2
    d["VarP1"] = d["EP1s"] - d["EP1"] ** 2
    d["Cov"] = d["EGP1"] - d["EG"] * d["EP1"]
    d["S"] = d["EH"] + d["EP2"]
    d["D"] = d["VarG"] + d["VarP1"] + 2.0 * d["Cov"]
    d["D2"] = d["EF2"] - d["EF"] ** 2
    return d


def full(N=2000, M=400):
    x, wx = leggauss(N)
    y, wy = leggauss(M)
    psi = 0.5 * math.pi * (x + 1.0)
    wpsi = 0.5 * math.pi * wx
    PS, UU = np.meshgrid(psi, y, indexing="ij")
    S, G, H, P1, P2 = kern(PS, UU)
    w = np.outer(wpsi, wy) * np.sin(PS) ** 2 * np.exp(-S)
    return moments(w, G, H, P1, P2)


def full_leb(N=2000, M=400):
    """S-测度: 权 dxi = 16 pi psi^2 dpsi du (无 J), 即 w = psi^2."""
    x, wx = leggauss(N)
    y, wy = leggauss(M)
    psi = 0.5 * math.pi * (x + 1.0)
    wpsi = 0.5 * math.pi * wx
    PS, UU = np.meshgrid(psi, y, indexing="ij")
    S, G, H, P1, P2 = kern(PS, UU)
    w = np.outer(wpsi, wy) * PS ** 2 * np.exp(-S)
    return moments(w, G, H, P1, P2)


def slice_a0(N=4000):
    """a=0 切片: psi in (0,pi), 权 e^{-S} sin^2(psi)/psi, u=0."""
    x, wx = leggauss(N)
    psi = 0.5 * math.pi * (x + 1.0)
    u = np.zeros_like(psi)
    S, G, H, P1, P2 = kern(psi, u)
    w = (0.5 * math.pi * wx) * np.sin(psi) ** 2 / psi * np.exp(-S)
    return moments(w, G, H, P1, P2)


LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


P("=" * 100)
P("v06_decide.log -- V06 裁决: G''(0)=0 (s07) vs 0.5429 (r03)   beta=2.5, theta=pi/2, 单 plaquette SU(2)")
P("独立实现: S = beta(1-ct cos psi + st sin psi u); d_a^2 用一般链式算子 (rho,u 坐标)")
P("=" * 100)

f1 = full(2000, 400)
f2 = full(1200, 200)
s1 = slice_a0(4000)
s2 = slice_a0(2000)
leb = full_leb(2000, 400)

P("")
P("[F] 全测度 dHaar = 4 pi sin^2(psi) dpsi du  (s07 的测度)")
P("    量            GL(2000,400)        GL(1200,200)        |diff|")
for k, lab in (("EH", "E[H_vv]"), ("EP2", "E[Phi'']"), ("VarG", "Var(G)"),
               ("VarP1", "Var(Phi')"), ("Cov", "Cov(G,Phi')")):
    P("    %-12s %+.15f   %+.15f   %.1e" % (lab, f1[k], f2[k], abs(f1[k] - f2[k])))
P("    %-12s %+.15f   %+.15f   %.1e" % ("E[G]", f1["EG"], f2["EG"], abs(f1["EG"] - f2["EG"])))
P("    %-12s %+.15f   %+.15f   %.1e" % ("E[Phi']", f1["EP1"], f2["EP1"], abs(f1["EP1"] - f2["EP1"])))
P("    %-12s %+.15f   %+.15f   %.1e" % ("E[G+Phi']", f1["EF"], f2["EF"], abs(f1["EF"] - f2["EF"])))
P("    %-12s %+.15f   %+.15f   %.1e" % ("E[G^2]", f1["EG2"], f2["EG2"], abs(f1["EG2"] - f2["EG2"])))
P("    %-12s %+.15f   %+.15f   %.1e" % ("E[F^2]", f1["EF2"], f2["EF2"], abs(f1["EF2"] - f2["EF2"])))
P("")
P("    S = E[H]+E[Phi'']          = %+.15f" % f1["S"])
P("    D = Var(G)+Var(Phi')+2Cov  = %+.15f" % f1["D"])
P("    D' = Var(G+Phi')           = %+.15f" % f1["D2"])
P("    S - D                      = %+.3e" % (f1["S"] - f1["D"]))
P("    S - D'                     = %+.3e" % (f1["S"] - f1["D2"]))
P("    E[H] - E[G^2]              = %+.15f   ( != 0 => q15 的 E[H]=E[G^2] 在 V-测度下不成立 )"
  % (f1["EH"] - f1["EG2"]))
P("    E[G]*E[Phi']               = %+.6f   ( E[H]-E[G^2] = E[G Phi'] ? )" % (f1["EG"] * f1["EP1"]))
P("")
P("[L] S-测度 dxi = 16 pi psi^2 dpsi du (无 J, Lebesgue) -- 核 q15 的 E[H]=E[G^2]")
P("    E_L[H] = %+.15f   E_L[G^2] = %+.15f   diff = %+.3e"
  % (leb["EH"], leb["EG2"], leb["EH"] - leb["EG2"]))
P("    E_L[Phi''] = %+.15f  (同一测度下 Phi 只是被积分项, 不参与恒等式)" % leb["EP2"])
P("")
P("[S] a=0 切片 dmu_a = 2 pi sin^2(psi)/psi dpsi du   ( = r03 的 nu, = s07 的 slices(a=0) )")
P("    量            GL(4000)            GL(2000)            |diff|")
for k, lab in (("EH", "E[H|0]"), ("EP2", "E[Phi''|0]"), ("VarG", "Var(G|0)"),
               ("VarP1", "Var(Phi'|0)"), ("Cov", "Cov(G,Phi'|0)")):
    P("    %-12s %+.15f   %+.15f   %.1e" % (lab, s1[k], s2[k], abs(s1[k] - s2[k])))
P("    %-12s %+.15f   %+.15f   %.1e" % ("E[G+Phi'|0]", s1["EF"], s2["EF"], abs(s1["EF"] - s2["EF"])))
P("    E[G|0]      = %+.15f   ( != 0: 切片上 score 恒等式不成立 )" % s1["EG"])
P("    E[G^2|0]    = %+.15f" % s1["EG2"])
P("    S_slice     = %+.15f" % s1["S"])
P("    D_slice     = %+.15f" % s1["D"])
P("    S_slice - D_slice = %.15f" % (s1["S"] - s1["D"]))
P("    (r03:  S=0.558183642724587...,  D=0.015239720370897...,  gap=+0.542943922353690)")
P("")
P("[A] a-切片族: G''(a) = E[H+Phi''|a] - Var(F|a),  边际权 da ∝ Z(a)da")
x, wx = leggauss(81)
phi = 0.5 * math.pi * x
aw = 0.5 * math.pi * wx
a_grid = 2.0 * math.pi * np.sin(phi)
amarg = np.cos(phi) * aw
npsi = 800
xg, wg = leggauss(npsi)
t = 0.5 * (xg + 1.0)
wt = 0.5 * wg
res = []
for a in a_grid:
    lo = abs(a) / 2.0
    psi = lo + (math.pi - 1e-15 - lo) * t
    dpsi = (math.pi - 1e-15 - lo) * wt
    rho = 2.0 * psi
    u = np.where(rho > 0, a / rho, 0.0)
    S, G, H, P1, P2 = kern(psi, u)
    w = np.sin(psi) ** 2 / psi * np.exp(-S) * dpsi
    Z = w.sum()
    E = lambda q: (w * q).sum() / Z
    F = G + P1
    res.append((Z, E(G), E(H), E(P1), E(P2), E(F), E(F * F)))
res = np.array(res)
wei = amarg * res[:, 0]
wei = wei / wei.sum()
mF = res[:, 5]
Gpp = (res[:, 2] + res[:, 4]) - (res[:, 6] - res[:, 5] ** 2)
mean_mF = float((wei * mF).sum())
var_mF = float((wei * (mF - mean_mF) ** 2).sum())
avg_Gpp = float((wei * Gpp).sum())
i0 = int(np.argmin(np.abs(a_grid)))
mask = a_grid <= 5.5
P("    G''(0) (a=0 切片, 独立路径) = %+.9f" % Gpp[i0])
P("    min_a G''(a) = %+.6f @ a=%+.3f ;  max(|a|<=5.5) = %+.6f"
  % (Gpp[mask].min(), a_grid[mask][int(np.argmin(Gpp[mask]))], Gpp[mask].max()))
P("    E_a[G''] = %+.6f    Var_a(m_F) = %+.6f    E[m_F]^2 = %+.2e    diff = %.1e"
  % (avg_Gpp, var_mF, mean_mF ** 2, avg_Gpp - var_mF - mean_mF ** 2))
P("    E_a[Var(F|a)] = Var(F) - Var_a(m_F) = %.6f - %.6f = %+.6f"
  % (f1["D2"], var_mF, f1["D2"] - var_mF))
P("    (s07 [S4] 报 E_a[G'']=+0.646732, Var_a(m_F)=+0.646732)")
P("")
P("[V] 判:  s07 的 0 与 r03 的 0.5429 是**不同测度**上的两个量")
P("    全测度:  S = %+.9f,  D = %+.9f,  S-D = %+.1e   <- Stein 等式 (E[V']=0)" % (f1["S"], f1["D"], f1["S"] - f1["D"]))
P("    a=0 切片: S = %+.9f,  D = %+.9f,  S-D = %+.9f   <- m_F'(0) = G''(0)" % (s1["S"], s1["D"], s1["S"] - s1["D"]))
P("    差: S_full - S_slice = %+.9f ;  D_full - D_slice = %+.9f" % (f1["S"] - s1["S"], f1["D"] - s1["D"]))

with open("v06_decide.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG) + "\n")
