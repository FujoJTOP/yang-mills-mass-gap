# -*- coding: utf-8 -*-
r"""
q11_theta.py -- Q11 配套 [C] 数值: theta-表 + beta-窗口 + 与严格区间互证
   单 plaquette SU(2), 切片 xi3=0, v=e3, 密度 ∝ sin^2(u)/u du  (u=rho/2)
   一般 theta 约化 (同 q11_cert.py 的推导):
       c=cos(th/2), s=sin(th/2), B = beta c
       nu_B ∝ sin^2(u)/u * exp(B cos u) du
       G = beta s g(u),  H = (beta c/2) g(u),  Phi'≡0 (切片切于等 rho 球面),  Phi''=1/(2u^2)-cot(u)/(2u)
       S = (beta c/2) E_nu[g] + E_nu[Phi''] ,   D = beta^2 s^2 Var_nu(g)
theta=pi 时 c=0,s=1: S=E_nu[Phi'']=1/(4Z) 与 beta 无关, D=beta^2 Var(g).
所有被积函数用整函数形式 (cot 已消):  Sf = sin^2u/(2u^3) - sin(2u)/(4u^2).
"""
import math
import numpy as np

NG = 400
XP, WP = np.polynomial.legendre.leggauss(NG)
NPAN = 600


def gl(f, a, b, npan=NPAN):
    e = np.linspace(a, b, npan + 1)
    mid = 0.5 * (e[:-1] + e[1:]); half = 0.5 * (e[1:] - e[:-1])
    t = (mid[:, None] + half[:, None] * XP[None, :]).ravel()
    w = (half[:, None] * WP[None, :]).ravel()
    return float(np.sum(w * f(t)))


def g(u):
    return np.where(u < 1e-12, 0.5, np.sin(u) / (2.0 * np.where(u < 1e-12, 1.0, u)))


def phipp(u):
    """Phi'' = 1/(2u^2) - cot(u)/(2u); 小 u 用级数 1/6+u^2/90+u^4/945+u^6/9450 (避免 1/u^2 抵消)."""
    small = u < 1e-2
    us = np.where(small, 1.0, u)
    ex = 0.5 / us ** 2 - np.cos(us) / (2.0 * us * np.sin(us))
    ser = 1.0 / 6.0 + u ** 2 / 90.0 + u ** 4 / 945.0 + u ** 6 / 9450.0
    return np.where(small, ser, ex)


def Sf(u):
    """Phi'' * sin^2 u/u 的整函数形式 (= -(1/4) d/du[sin^2u/u^2]); 小 u 用级数 u/6-2u^3/45+u^5/210."""
    small = u < 1e-2
    us = np.where(small, 1.0, u)
    ex = np.sin(us) ** 2 / (2 * us ** 3) - np.sin(2 * us) / (4 * us ** 2)
    ser = u / 6.0 - 2.0 * u ** 3 / 45.0 + u ** 5 / 210.0
    return np.where(small, ser, ex)


LOG = []


def P(s=""):
    print(s); LOG.append(s)


P("=" * 96)
P("q11_theta.log -- [C] theta/beta 表 + 与严格区间 q11_cert.log 互证")
P("=" * 96)
P("")
P("[C0] 独立求积 (numpy Gauss-Legendre, %d 面板 x %d 点) vs 严格区间中点" % (NPAN, NG))
Z = gl(lambda u: np.where(u < 1e-12, u, np.sin(u) ** 2 / np.where(u < 1e-12, 1.0, u)), 0.0, math.pi)
NA = gl(lambda u: g(u) * np.sin(u) ** 2 / np.where(u < 1e-12, 1.0, u), 0.0, math.pi)
NB = gl(lambda u: g(u) ** 2 * np.sin(u) ** 2 / np.where(u < 1e-12, 1.0, u), 0.0, math.pi)
NS = gl(Sf, 0.0, math.pi)
RIG = {"Z": 1.2188266965286122059053551707, "NA": 0.4356224991559818334675841019,
       "NB": 0.1691757894863797036654118569, "NS": 0.25}
for nm, v in (("Z", Z), ("NA", NA), ("NB", NB), ("NS", NS)):
    P("   %-3s GL = %.16f   rigid(mid) = %.16f   |d| = %.2e" % (nm, v, RIG[nm], abs(v - RIG[nm])))
Eg, Eg2 = NA / Z, NB / Z
VG = Eg2 - Eg * Eg
P("   E[g]=%.16f  E[g^2]=%.16f  Var(g)=%.16f  (严格区间 0.0110592937729735899114)" % (Eg, Eg2, VG))
P("   S(NS/Z)=%.16f  (严格区间 0.2051152971230731479833 ; 语料锚 1/(2Cin2pi)=0.205115297123)"
  % (NS / Z))
P("")


def moments(B, beta, theta):
    """nu_B ∝ sin^2 u/u e^{B cos u};  返回 (S, D, E[g], E[Phi''], Var(g))."""
    c, s = math.cos(0.5 * theta), math.sin(0.5 * theta)
    w = lambda u: np.sin(u) ** 2 / np.where(u < 1e-12, 1.0, u) * np.exp(B * (np.cos(u) - 1.0))
    Z0 = gl(w, 0.0, math.pi)
    E1 = gl(lambda u: w(u) * g(u), 0.0, math.pi) / Z0
    E2 = gl(lambda u: w(u) * g(u) ** 2, 0.0, math.pi) / Z0
    F = gl(lambda u: w(u) * phipp(u), 0.0, math.pi) / Z0
    V = E2 - E1 * E1
    S = 0.5 * beta * c * E1 + F
    D = beta * beta * s * s * V
    return S, D, E1, F, V


P("[C1] beta=2.5, 一般 theta:  S - D  (S=(bc/2)E[g]+E[Phi''], D=b^2 s^2 Var(g))")
P("     theta/pi     c        s       E[g]       E[Phi'']     Var(g)       S          D        S-D      sign")
for t in (0.0, 0.25, 0.5, 0.75, 0.9, 0.95, 1.0):
    th = t * math.pi
    c, s = math.cos(0.5 * th), math.sin(0.5 * th)
    S, D, E1, F, V = moments(2.5 * c, 2.5, th)
    P("     %6.4f  %+8.5f %+8.5f  %.8f  %.8f  %.9f  %.8f  %.8f  %+.8f   %s"
      % (t, c, s, E1, F, V, S, D, S - D, "+" if S - D > 0 else "-"))
P("   注意 theta=0 时 s=0 ⇒ D=0 (v 方向无涨落), theta=pi 时 c=0 ⇒ H≡0 (唯一失去 beta-线性正项的点)")
P("")
P("[C2] theta=pi 的 beta-窗口 (S 与 beta 无关, D=beta^2*Var(g)):")
P("     beta      S          D          S-D       sign")
for b in (0.5, 1.0, 2.5, 3.0, 4.0, 4.2, 4.3, 4.3066, 4.4, 5.0, 6.0):
    S, D, *_ = moments(0.0, b, math.pi)
    P("     %6.3f  %.9f  %.9f  %+.9f   %s" % (b, S, D, S - D, "+" if S - D > 0 else "-"))
bst = math.sqrt((NS / Z) / VG)
P("     beta* = sqrt(S/Var(g)) = %.10f   (h1e_var.log:20 = 4.30660785)" % bst)
P("     => theta=pi 上 S-D>0 <=> beta < beta* ; beta=2.5 有 %.1f%% 余量" % (100 * (1 - 2.5 / bst)))
P("")
P("[C3] L=2: 两个 plaquette 的 SU(2) 格点. 不可化. 理由:")
P("  - L=2 (1x1x1x1 环面, 16 link) 的 Z 与边缘是 SU(2)^{16} 上的 16 重积分, 只能做 character 展开;")
P("    beta=2.5 上展开系数 ~ (beta/4)^n 衰减极慢 (n 需 >>10), 无闭式截断; 精确值不存在有限维表示.")
P("  - 单 plaquette 的可积性来自: 只有一个 link 变量 + 一个约束 => 2D 切片上的 1D 条带 (u 一条), 因此可化 1D.")
P("    L=2 的条件化切片 (固定 <v,xi_1>=0) 仍是 15 维流形, 该 1D 约化不复存在.")
P("  - 结论: L=2 只能 MC (或 character 展开的逐阶界). 本文件不做.")
P("")
P("done.")

with open("q11_theta.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG) + "\n")
