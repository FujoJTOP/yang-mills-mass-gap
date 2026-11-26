#!/usr/bin/env python3
"""v0_sign.py -- V0(N) 的符号：闭式(N=2)是否可延到 N>=3；若不可，从 MC 三点拟合。

只读 betastar_n34.log 的已发布数字（不重跑 MC）。
运行: py -3.12 research/unconditional/w7/v0_sign.py  >  ...v0_sign.log
"""
import numpy as np
from scipy.special import sici
from numpy.polynomial.legendre import leggauss

np.set_printoptions(precision=15, suppress=False)
EULER = 0.5772156649015328606
SEP = "=" * 78


def Ci(x):
    return sici(float(x))[0]


def Cin(x):
    """Cin x = int_0^x (1-cos t)/t dt = gamma + ln x - Ci(x)."""
    return EULER + np.log(float(x)) - Ci(x)


def gl(f, a, b, n=4000):
    """Gauss-Legendre on [a,b] (被积函数在端点解析)."""
    x, w = leggauss(n)
    xs = 0.5 * (b - a) * x + 0.5 * (b + a)
    return 0.5 * (b - a) * np.sum(w * f(xs))


def fit(x, y, s, tag):
    """加权最小二乘 y = a + c*x，x = beta^2。返回 V0eff = -c。"""
    x = np.asarray(x, float); y = np.asarray(y, float); s = np.asarray(s, float)
    w = 1.0 / s ** 2
    X = np.vstack([np.ones_like(x), x]).T
    A = (X.T * w) @ X
    rhs = (X.T * w) @ y
    cov = np.linalg.inv(A)
    p = cov @ rhs
    r = y - X @ p
    chi2 = float(np.sum(w * r ** 2)); dof = len(x) - 2
    # 独立复算：numpy.polyfit 加权
    p2 = np.polyfit(x, y, 1, w=np.sqrt(w))
    print(f"  [{tag}] n={len(x)} dof={dof}")
    print(f"    a  = {p[0]:+.6f} +- {np.sqrt(cov[0,0]):.6f}   (polyfit {p2[1]:+.6f})")
    print(f"    c  = {p[1]:+.6e} +- {np.sqrt(cov[1,1]):.3e}  (polyfit {p2[0]:+.6e})")
    print(f"    V0eff = -c = {-p[1]:+.6e} +- {np.sqrt(cov[1,1]):.3e}")
    print(f"    chi2/dof = {chi2/max(dof,1):.3f}   sigma(V0eff) 意义: |V0eff|/sig = {abs(p[1])/np.sqrt(cov[1,1]):.2f}")
    sig = np.sqrt(cov[1,1])
    if p[1] > 0:
        print(f"    => V0eff < 0 （delta 随 beta^2 上升），{abs(p[1])/sig:.1f} sigma")
    else:
        print(f"    => V0eff > 0，{abs(p[1])/sig:.1f} sigma")
    return -p[1], sig, chi2 / max(dof, 1)


# ----------------------------------------------------------------------------
print(SEP)
print("[0] 定义（先固定记号）")
print("  Cin x = int_0^x (1-cos t)/t dt = gamma + ln x - Ci(x)      （余弦积分 Ci）")
print("  Ci x  = gamma + ln x + int_0^x (cos t - 1)/t dt")
print(SEP)

print("[1] Cin / Ci 数值 + 独立求积核对（5 位核准）")
checks = {
    "Cin(2pi)": (Cin(2 * np.pi), gl(lambda t: (1 - np.cos(t)) / t, 1e-12, 2 * np.pi, 20000)),
    "Cin(pi)":  (Cin(np.pi),     gl(lambda t: (1 - np.cos(t)) / t, 1e-12, np.pi, 20000)),
}
for k, (cf, qu) in checks.items():
    print(f"  {k:9s} 闭式 {cf:.15f}   求积 {qu:.15f}   差 {cf-qu:+.2e}")
print(f"  Cin(2pi) 5 位: {Cin(2*np.pi):.5f}    （H1E 存档 2.43765）")
for m in (1, 2, 3, 4):
    print(f"  Ci({m}pi) = {Ci(m*np.pi):+.12f}")

# ----------------------------------------------------------------------------
print(SEP)
print("[2] N=2 闭式：V0, W0, beta*（只依赖 Z, N1, N2）")
Z_cf = 0.5 * Cin(2 * np.pi)
N1_cf = 0.75 * (np.log(3.0) + Ci(np.pi) - Ci(3 * np.pi))
N2_cf = np.log(2.0) + Ci(2 * np.pi) - Ci(4 * np.pi)
Z_q = gl(lambda x: np.sin(x) ** 2 / x, 0.0, np.pi)
N1_q = gl(lambda x: np.sin(x) ** 3 / x ** 2, 0.0, np.pi)
N2_q = gl(lambda x: np.sin(x) ** 4 / x ** 3, 0.0, np.pi)
for nm, cf, qu in (("Z", Z_cf, Z_q), ("N1", N1_cf, N1_q), ("N2", N2_cf, N2_q)):
    print(f"  {nm:3s} 闭式 {cf:.15f}   求积 {qu:.15f}   差 {cf-qu:+.2e}")
V0_2 = N2_cf / (4 * Z_cf) - (N1_cf / (2 * Z_cf)) ** 2
W0_2 = 1.0 / (2 * Cin(2 * np.pi))
Ef, Ef2 = N1_cf / (2 * Z_cf), N2_cf / (4 * Z_cf)
print(f"  E[f]  = N1/(2Z) = {Ef:.15f}")
print(f"  E[f2] = N2/(4Z) = {Ef2:.15f}   (= E[G^2]/beta^2, H1E 0.138802)")
print(f"  V0(2) = {V0_2:+.15f}   [T]")
print(f"  W0(2) = {W0_2:+.15f}   [T]")
print(f"  beta*(2) = sqrt(W0/V0) = {np.sqrt(W0_2/V0_2):.8f}   (H1A 4.30661)")
print(f"  1/(N2-4W0N1^2) = {1.0/(N2_cf - 4*W0_2*N1_cf**2):.8f}")

# ----------------------------------------------------------------------------
print(SEP)
print("[3] N-依赖逐项追查：Z, N1, N2 里有 N 吗？")
print("  模型: g = (1/N) Re tr(U B),  S = beta(1 - g),  G = <v,grad S> = -beta * d_v g")
print("  N 只能从下面 4 个槽进入，逐个点名:")
print("   (i) 归一化 1/N  : 显式, 进 G 的整体标度")
print("   (ii) 维数 d=N^2-1 : 进测度/Haar 权（N=2 时 d=3 -> 2D 切片 -> 权 sin^2x/x）")
print("   (iii) tr(T_a T_b)=delta/2, tr(1)=N : 进算符界 |grad g| <= 1/sqrt(2N)")
print("   (iv) 切片结构: theta=pi 的双侧平分赤道 Re tr(U B)=0")
print()
print("  Z, N1, N2 是 x in (0,pi) 上的一维积分:")
print("   Z  = int sin^2 x / x      dx      -> 无 N")
print("   N1 = int sin^3 x / x^2    dx      -> 无 N")
print("   N2 = int sin^4 x / x^3    dx      -> 无 N")
print("  => 闭式里没有 N 槽。它不是 'V0(N) 的 N=2 值'，是 'N=2 这一格的积分值'。")
print("  且 (iv) 是 SU(2) 独有: SU(N>=3) 的 Re tr(U B) 不能恒为 0 且割迹非线性切片。")
print("  [T] 结论: 闭式不可延到 N>=3。H1E §5 已标 (c)，此处逐项确认。")
print()
print("  算符界（任何 N）: V0(N) <= 1/(2N); 但 N=2 实测 E[G^2]/beta^2 = %.6f = 界的 %.1f%%"
      % (Ef2, 100 * Ef2 / 0.25))
for n in (2, 3, 4, 5):
    print(f"    N={n}: bound 1/(2N) = {1/(2*n):.6f}")

# ----------------------------------------------------------------------------
print(SEP)
print("[4] 闭式硬套到 N=3 的证伪（不是'不可延'的口头判断，是数据反证）")
b3 = np.array([0.0, 6.0, 10.0])
d3 = np.array([0.34368, 0.63393, 0.95757])
s3 = np.array([0.0077, 0.0676, 0.1446])
pred3 = W0_2 - V0_2 * b3 ** 2
for b, d, s, pr in zip(b3, d3, s3, pred3):
    print(f"  beta={b:5.2f}: MC delta={d:+.5f}+-{s:.4f} | 闭式(N=2)={pr:+.5f} | 差 {(d-pr)/s:+.1f} sigma")
print("  => beta=6 一处就差 %.1f sigma。闭式在 N=3 被数据直接排除。[T]/[E]" % ((d3[1]-pred3[1])/s3[1]))

# ----------------------------------------------------------------------------
print(SEP)
print("[5] (B) 从 MC 拟合 V0(3): delta = a + c beta^2, V0eff = -c")
f3 = fit(b3 ** 2, d3, s3, "N=3 全 3 点")
print("  留一法（判符号是否靠单点）:")
for i in range(3):
    m = np.ones(3, bool); m[i] = False
    v, sg, _ = fit(b3[m] ** 2, d3[m], s3[m], f"N=3 去掉 beta={b3[i]:.0f}")
print("  区间斜率（判 delta 是否在 beta^2 上凹）:")
for i in range(2):
    sl = (d3[i + 1] - d3[i]) / (b3[i + 1] ** 2 - b3[i] ** 2)
    ss = np.hypot(s3[i], s3[i + 1]) / (b3[i + 1] ** 2 - b3[i] ** 2)
    print(f"    [{b3[i]:.0f},{b3[i+1]:.0f}]: c = {sl:+.5e} +- {ss:.2e}  (V0eff={-sl:+.5e})")
sl1 = (d3[1] - d3[0]) / (b3[1] ** 2 - b3[0] ** 2)
sl2 = (d3[2] - d3[1]) / (b3[2] ** 2 - b3[1] ** 2)
ss1 = np.hypot(s3[0], s3[1]) / (b3[1] ** 2 - b3[0] ** 2)
ss2 = np.hypot(s3[1], s3[2]) / (b3[2] ** 2 - b3[1] ** 2)
print(f"    Delta c = {sl2-sl1:+.3e} +- {np.hypot(ss1,ss2):.2e}  -> 凹性 {(sl2-sl1)/np.hypot(ss1,ss2):+.1f} sigma（不显著）")
turn = b3[2] ** 2 - sl2 / ((sl1 - sl2) / (b3[2] ** 2 - b3[1] ** 2))
print(f"    [NG] 若凹性为真：c->0 外推 beta^2 ~ {turn:.0f} => beta ~ {np.sqrt(abs(turn)):.1f}（在实测窗外，不可判）")

# ----------------------------------------------------------------------------
print(SEP)
print("[6] N=4 同法（MC 6 点，误差大）")
b4 = np.array([0.0, 4.0, 8.0, 11.15, 14.0, 18.0])
d4 = np.array([0.48678, 0.50754, 0.57527, 0.66871, 0.78610, 0.99508])
s4 = np.array([0.0400, 0.0687, 0.1581, 0.3794, 0.8509, 2.0740])
f4all = fit(b4 ** 2, d4, s4, "N=4 全 6 点")
m8 = b4 <= 8.0
f4lo = fit(b4[m8] ** 2, d4[m8], s4[m8], "N=4 beta<=8")
print("  => N=4 符号未定（|V0eff|/sigma < 1..2）；N=5 无 MC 数据，不给数。")

# ----------------------------------------------------------------------------
print(SEP)
print("[7] 分解：符号为什么翻？(log 列 H, V=varG, Phi2 都在 betastar_n34.log 里)")
print("  betastar_scaling.py:130  δ = E[H] - varG + W_geo + X,  H = -β ∂²f,  G = -β ∂f")
print("  N=2 在赤道上 E[H]=0 精确、测度 β-无关 ⇒ δ = W0 - V0 β²（严格）")
print("  N>=3 无赤道 ⇒ E[H] != 0，且随 β 涨：")
for nm, bb, dd, ss, HH, VV, P2 in (
    ("N=3", b3, d3, s3, np.array([0.0, 0.3064, 0.6478]), np.array([0.0, 0.03239, 0.07483]),
     np.array([0.3444, 0.3592, 0.3816])),
    ("N=4", b4, d4, s4, np.array([0.0, 0.0257, 0.1082, 0.2211, 0.3595, 0.5971]),
     np.array([0.0, 0.00777, 0.03375, 0.06884, 0.10880, 0.16788]),
     np.array([0.4908, 0.4942, 0.5070, 0.5236, 0.5417, 0.5670])),
):
    print(f"  {nm}:")
    print(f"    {'beta':>6s} {'delta':>10s} {'dH/db2':>11s} {'V/b2 (=V0var)':>15s} {'dPhi2/db2':>11s}")
    for i in range(1, len(bb)):
        x2 = bb[i] ** 2
        print(f"    {bb[i]:6.2f} {dd[i]:+10.5f} {HH[i]/x2:+11.5f} {VV[i]/x2:+15.6e} {(P2[i]-P2[0])/x2:+11.5f}")
    print(f"    残差检查 delta-(H-varG+Phi2) = "
          f"{np.max(np.abs(dd-(HH-VV+P2))):.2e}  (未测的 varD+2cov, 应为 O(1e-3))")
    print(f"    => 有效斜率的符号 = dH/dβ² + dPhi2/dβ² - V0var;  "
          f"第一项压过第三项 ⇒ V0eff < 0")

# ----------------------------------------------------------------------------
print(SEP)
print("[8] 统计量需求（判 N=4 符号 / 判 N=3 是否在窗外翻转）")
v4, s4v, _ = f4all
print(f"  N=4: 现 V0eff = {v4:+.2e} +- {s4v:.2e}  ({(abs(v4)/s4v):.1f} sigma)")
need = (3.0 * s4v / max(abs(v4), 1e-12)) ** 2 if abs(v4) > 0 else np.inf
print(f"       要 3 sigma: sigma 需缩 {3*s4v/abs(v4):.1f}x => 统计量 x{need:.0f}"
      f" => 每 β 有效样本 {1e5*need/1e6:.2f}e6, 约 {420.1*need/60:.0f} min CPU (accept 7.8e-4)")
print(f"  N=3: 现 V0eff = {-f3[0]:+.2e} +- {f3[1]:.2e} ({abs(f3[0])/f3[1]:.1f} sigma) -> 符号已定")
print("       要在 beta in [10,20] 判是否翻转: 需 beta=14,20 两点, sigma(delta)<=0.1 各点 (~2e5 有效样本)")
print(SEP)
