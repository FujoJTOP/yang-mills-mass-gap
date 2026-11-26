"""Q17 — p19 的 Var(F) ~ beta^1.57: 真幂, 还是 beta in [1,6] 的窗口效应?

对照 (单 plaquette, theta=pi 切片: S == beta 严格 ⇒ 无高阶项): Var = beta^2 * V0,
局域指数恒 == 2 (h1e_var.log:12/46-48 精确, o05_var_scaling.log:65-72 k=2.0000 全区间)。

数据源: research/unconditional/w7/p19_fit.py:7   Var(F) = 0.2888 1.5232 2.8421 4.7329
        research/unconditional/w7/p19_b{1,2.5,4,6}.log:22  jackknife / tau 误差
只读本脚本, 不重跑 MC。
"""
import numpy as np

B = np.array([1.0, 2.5, 4.0, 6.0])
V = np.array([0.2888, 1.5232, 2.8421, 4.7329])      # p19_fit.py:7
SJ = np.array([0.0058, 0.0426, 0.0948, 0.1716])     # jackknife
ST = np.array([0.0096, 0.1018, 0.2184, 0.4634])     # tau-corrected (Neff~25-1000)
V0 = 0.011059293772974                             # h1e_var.log:12 (theta=pi, 单 plaquette)
L = np.log(B)


def slope(v, i, j, sig=None):
    e = np.log(v[j] / v[i]) / (L[j] - L[i])
    s = np.nan if sig is None else np.hypot(sig[i] / v[i], sig[j] / v[j]) / (L[j] - L[i])
    return e, s


def wfit(p, q):
    """相对加权最小二乘 Var ~ a*b^p + c*b^q (min sum ((model/V)-1)^2)."""
    X = np.stack([B ** p, B ** q], 1) / V[:, None]
    ab = np.linalg.lstsq(X, np.ones(4), rcond=None)[0]
    r = (ab[0] * B ** p + ab[1] * B ** q) / V - 1.0
    return ab, r


def lexp(ab, p, q, b):
    """该模型在 beta=b 的局域指数."""
    a1, a2 = ab[0] * b ** p, ab[1] * b ** q
    return (p * a1 + q * a2) / (a1 + a2)


print("=" * 96)
print("Q17  p19 (L=2, SU(2)) Var(F) 的 beta-幂审计   对照: 单 plaquette Var=beta^2*V0 (指数严格 2)")
print("  数据: p19_fit.py:7   beta =", B, "\n                       Var  =", V)
print("  题面给的四舍五入值 0.289/1.523/2.842/4.733 与全精度偏差 < 0.1% (不影响任何结论)")

print("\n[1] 四点局域对数导数 dlnVar/dlnbeta")
rules = [(0, 1, "FWD  [1,2.5]"), (0, 2, "CEN  [1,4]"), (1, 3, "CEN  [2.5,6]"), (2, 3, "BWD  [4,6]")]
EP = []
for i, j, lab in rules:
    e, sj = slope(V, i, j, SJ)
    _, st = slope(V, i, j, ST)
    EP.append((e, sj, st))
    print(f"   beta={B[i]:>4g}  {lab:<14s}  exponent = {e:5.3f}   +-{sj:.3f}(jk)  +-{st:.3f}(tau)")
print("   逐区间坡度 (相邻两点):", "  ".join(f"[{B[i]:g},{B[i+1]:g}]={slope(V,i,i+1)[0]:.3f}"
                                      for i in range(3)))
print(f"   V/beta^2 = {np.round(V/B**2,5)}  ->  单调下降, 没有任何 beta^2 主导的迹象")

print("\n[2] 单 plaquette 对照 (Var = beta^2*V0, V0 = %.12f)" % V0)
VC = B ** 2 * V0
print("   Var   =", np.round(VC, 6))
print("   四个局域指数 =", np.round([slope(VC, i, j)[0] for i, j, _ in rules], 12),
      "  (恒等 2, o05_var_scaling.log:65-72 独立复算 k=2.0000 全 7 个区间)")

print("\n[3] 窗口敏感性: 纯幂拟合 Var~beta^alpha 换窗口后 alpha 怎么变")
for ii, jj, lab in [(0, 1, "[1,2.5]"), (0, 2, "[1,4]"), (0, 3, "[1,6]"), (1, 3, "[2.5,6]"),
                    (2, 3, "[4,6]")]:
    a = slope(V, ii, jj)[0]
    print(f"   窗口 {lab:<10s} alpha = {a:.3f}")
a_full = np.polyfit(L, np.log(V), 1)[0]
r_full = np.log(V) - np.polyval(np.polyfit(L, np.log(V), 1), L)
print(f"   [1,6] 最小二乘 alpha = {a_full:.3f} (= p19_fit.log:19 的 1.57), 对数残差 = {np.round(r_full,4)}")
print("   -> 残差符号 (-- - +) 是系统性曲率, 不是噪声: 单幂律在 [1,6] 已不成立")

print("\n[4] 模型比较 (相对加权, 残差 = 拟合值/V - 1)")
models = []
ab, r = wfit(2.0, 0.0)
models.append(("Var = a*beta^2      (p19 若为 beta^2)", r, (ab, 2.0, 0.0)))
ab, r = wfit(2.0, 1.0)
models.append(("Var = a*beta^2+b*beta (含 beta^2 项)", r, (ab, 2.0, 1.0)))
ab, r = wfit(1.0, 0.0)
models.append(("Var = a*beta + c    (渐近幂 1)", r, (ab, 1.0, 0.0)))
ab, r = wfit(1.25, 0.0)
models.append(("Var = a*beta^1.25   (5/4)", r, (ab, 1.25, 0.0)))
ab, r = wfit(1.5714, 0.0)
models.append(("Var = a*beta^1.571  (p19 的拟合)", r, (ab, 1.5714, 0.0)))
for lab, r, _ in models:
    print(f"   {lab:<34s} rms={np.sqrt((r**2).mean()):6.3f}  max|r|={np.abs(r).max():6.3f}  r={np.round(r,3)}")
bq, bqq = None, None
for q in np.linspace(0.0, 1.95, 196):
    ab, r = wfit(2.0, q)
    nr = np.sqrt((r ** 2).mean())
    if bq is None or nr < bq[0]:
        bq = (nr, q, ab, r)
print(f"   {'Var = a*beta^2+b*beta^q (q 自由)':<34s} rms={bq[0]:6.3f}  q={bq[1]:.2f}  "
      f"b={bq[2][1]:+.4f}  r={np.round(bq[3],3)}")

best, fam = None, 0
for p in np.linspace(1.0, 2.6, 161):
    for q in np.linspace(0.3, 2.2, 191):
        if q > p - 0.05:
            continue
        ab, r = wfit(p, q)
        nr = np.sqrt((r ** 2).mean())
        fam += nr < 0.01
        if best is None or nr < best[0]:
            best = (nr, p, q, ab, r)
print(f"   最优双幂 a*b^p+b*b^q: p={best[1]:.2f} q={best[2]:.2f} rms={best[0]:.3f} "
      f"a={best[3][0]:+.3f} b={best[3][1]:+.3f}")
print(f"   rms<0.01 的 (p,q) 组合数 = {fam} / ~3e4  ->  四点定不住双幂 (参数退化)")
print("   各模型外推的局域指数:")
for lab, r, (ab, p, q) in models:
    print(f"     {lab:<34s} beta=8:{lexp(ab,p,q,8):5.2f}  10:{lexp(ab,p,q,10):5.2f}  20:{lexp(ab,p,q,20):5.2f}")

print("\n[5] 判决")
e = np.array([x[0] for x in EP])
s = np.array([x[2] for x in EP])
print(f"   四点局域指数 = {np.round(e,3)}  (tau 误差 {np.round(s,3)})")
print(f"   单调下降: {bool(np.all(np.diff(e) < 0))}   总降幅 = {e[0]-e[-1]:.3f} "
      f"({(e[0]-e[-1])/np.hypot(s[0],s[-1]):.1f} sigma)")
print(f"   高 beta 端 e(6) = {e[-1]:.3f} +- {s[-1]:.3f};  距 2 还有 {(2-e[-1])/s[-1]:.1f} sigma 且方向相反")
print("   => 不漂向 2: (a) 否;  1.57 是 [1,6] 上单调下降曲线的弦平均 (换窗口 1.26~1.81),")
print("      故 1.57 是窗口/弦效应, 但渐近幂也不是 2 — beta^2 只在单 plaquette 严格成立。")

print("\n[6] 可判别的下一步 (同类 MC, 每点 ~60 s: ns=5000, G=4 链, 见 p19_b6.log:2)")
for lab, f in [("H(2)  beta^2 ", lambda b: 0.1843 * b ** 2),
               ("H(5/4) beta^1.25", lambda b: 0.4970 * b ** 1.25),
               ("H(1)  beta^1 ", lambda b: 0.7497 * b)]:
    print(f"   {lab:<16s} Var(8)={f(8):5.2f}  Var(10)={f(10):5.2f}  Var(12)={f(12):5.2f}")
print("   (锚点: 0.1843 = mean(V/b^2, b>=2.5), 0.4970 = mean(V/b^1.25, b>=2.5), 0.7497 = mean(V/b, b>=4))")
print("   beta=10 一点即可分开 H(2)=13.1  vs  H(5/4)=8.8  vs  H(1)=7.5; beta in {8,10,12,16} 分辨 5/4 与 1。")

# ---------------- self-check ----------------
assert abs(a_full - 1.571) < 0.005, a_full
assert abs(e[0] - 1.814) < 0.005 and abs(e[1] - 1.649) < 0.005
assert abs(e[2] - 1.295) < 0.005 and abs(e[3] - 1.258) < 0.005
assert np.all(np.diff(e) < 0)
assert np.allclose([slope(VC, i, j)[0] for i, j, _ in rules], 2.0, atol=1e-12)
ab2, r2 = wfit(2.0, 1.0)
assert np.abs(r2).max() > 0.1, "beta^2 项应当拟合失败"
lv = V[1:] / B[1:] ** 1.25
assert (lv.max() - lv.min()) / lv.mean() < 0.05, lv
assert np.abs(V / np.array([0.289, 1.523, 2.842, 4.733]) - 1).max() < 1e-3
print("\n[SELF-CHECK] 全部断言通过")
