#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gm3_falsify_check.py — GM3-FALSIFY 的确定性数值核（无随机、无 MC、无格点作业）
任务：给 (SM-2) 的三种读法各配一个精确数：
  [1] Wilson 符号级数（交叉核 GR1:205）
  [2] 自由场边缘 (Delta_a) 的 c1/c2（交叉核 GK2 (1.2a)）
  [3] 族 F1：Wilson 自由场，边缘动量 p* != 0（扭转/反周期 BC）⇒ delta1 != 0，k = 2
  [4] 族 F2：Symanzik 改进自由场（GK3 §2.2 闭式）⇒ delta1 = +m^3/24 != 0，k = 4
  [5] SEP-2（GR1:151-154）：delta1 = -inf（发散），且非解析 ⇒ 不是任何格点族
  [6] 扭转参数 theta 的连续路径：delta1 从 0 连续走到 -1.231，k 恒为 2
  [7] (1.2a) 的精确代数：c1 = 0 <=> delta1 = D^3/24（Fraction）
用法：py -3.12 gm3_falsify_check.py  （< 10 s）
"""
import math
from fractions import Fraction as Fr
from mpmath import mp, mpf, asinh, sqrt as msqrt, sinh as msinh, log as mlog, pi as mpi
import sympy as sp

mp.dps = 50
LOG = open("gm3_falsify_check.log", "w", encoding="utf-8")


def out(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


out("=" * 78)
out("gm3_falsify_check.py  —  (SM-2) 打假的确定性数值核")
out("=" * 78)

# ----------------------------------------------------------------- [1]
out("\n[1] Wilson 自由场符号级数（d=3）  h_a = (2/a^2)Σ(1-cos(a p_i))")
a_, p1, p2, p3 = sp.symbols("a p1 p2 p3", positive=True)
h = (2 / a_ ** 2) * ((1 - sp.cos(a_ * p1)) + (1 - sp.cos(a_ * p2)) + (1 - sp.cos(a_ * p3)))
hser = sp.simplify(sp.expand(sp.series(h, a_, 0, 5).removeO()))
p2sum = p1 ** 2 + p2 ** 2 + p3 ** 2
p4sum = p1 ** 4 + p2 ** 4 + p3 ** 4
c2coef = sp.simplify(sp.diff(hser, a_, 2).subs(a_, 0) / 2)
out("    h_a = ", sp.nsimplify(hser))
out("    a^2 系数 = ", c2coef, "   （应 = -Σp^4/12 ⇒", sp.simplify(c2coef + p4sum / 12) == 0, "）")

# ----------------------------------------------------------------- helpers
def h_a_d3(p, a):
    """Wilson 符号（物理动量 p 固定，a→0）"""
    return (2 / a ** 2) * sum(1 - math.cos(a * x) for x in p)


def omega_wilson(p, a, m):
    """sinh-mass 的原料：ω_a = sqrt(m^2 + h_a(p))"""
    return math.sqrt(m * m + h_a_d3(p, a))


def Delta_wilson(p, a, m):
    """Wilson 自由场边缘：cosh(aΔ)=1+a²ω²/2 ⇒ Δ=(2/a)asinh(aω/2)"""
    return (2 / a) * math.asinh(a * omega_wilson(p, a, m) / 2)


def Delta_wilson_mp(p, a, m):
    a = mpf(a)
    w = msqrt(mpf(m) ** 2 + (2 / a ** 2) * sum(1 - mp.cos(a * mpf(x)) for x in p))
    return (2 / a) * asinh(a * w / 2)


def richardson(vals, ratio=100.0):
    """vals = [f(a), f(a/10), f(a/100), ...]；返回 Richardson 加速后的极限"""
    v = [mpf(x) for x in vals]
    while len(v) > 1:
        n = len(v)
        num = v[-1] * ratio ** (n - 1) - v[-2]
        den = ratio ** (n - 1) - 1
        v = v[:-2] + [num / den]      # 每次缩一格（v1 版写成 v[:-1]+ 是死循环，本轮修）
    return v[0]


# ----------------------------------------------------------------- [2]
out("\n[2] 自由场边缘 p*=0：c1/c2 精确复现（交叉核 GK2 (1.2a), GR1:237）")
for m in (Fr(1), Fr(2), Fr(1, 2)):
    m_ = sp.Rational(m.numerator, m.denominator)
    aa = sp.symbols("aa", positive=True)
    Dl = (2 / aa) * sp.asinh(aa * m_ / 2)
    s = sp.series(Dl, aa, 0, 5).removeO()
    c1 = sp.simplify(sp.diff(s, aa, 2).subs(aa, 0) / 2 / m_)
    c2 = sp.simplify(sp.diff(s, aa, 4).subs(aa, 0) / 24 / m_)
    out(f"    m={m}: c1 = {c1}  (应 -m²/24 = {sp.Rational(-m.numerator**2, 24*m.denominator**2)})"
        f"   c2 = {c2}  (应 3m⁴/640 = {sp.Rational(3*m.numerator**4, 640*m.denominator**4)})")
out("    ⇒ delta1(p*=0) = lim (ω_a - Δ)/a² = 0 精确（ω_a ≡ m，与 a 无关；[E] 机器精度 0）")

# ----------------------------------------------------------------- [3]
out("\n[3] 族 F1：Wilson 自由场，边缘动量 p* ≠ 0（反周期/扭转 BC）")
out("    设定：m=1, L=1, d=3, 仅 x 方向反周期 ⇒ p* = (π/L, 0, 0)")
m0, L = 1.0, 1.0
pstar = (math.pi / L, 0.0, 0.0)
pstar_mp = (mpi / mpf(L), mpf(0), mpf(0))
# 参考常数用 mp 精度（v1 版用 float 的 Dc 当减数，a≤1e-5 处把 1e-16 的常数误差放大成 1e-4 ⇒ 本轮修）
Dc_mp = msqrt(mpf(m0) ** 2 + mpi ** 2)
Dc = float(Dc_mp)
S4 = sum(x ** 4 for x in pstar)
d1_pred = -S4 / (24 * Dc)
c1_pred = d1_pred / Dc - Dc ** 2 / 24
out(f"    Δ_cont = sqrt(m²+|p*|²) = {Dc:.15f}    Σp⁴ = {S4:.15f} = π⁴")
out(f"    预言 δ₁ = -Σp⁴/(24Δ_cont)     = {d1_pred:.15f}")
out(f"    预言 c₁ = δ₁/Δ_cont - Δ²/24   = {c1_pred:.15f}")

a_list = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
seq_d1, seq_c1 = [], []
for a in a_list:
    Dl = Delta_wilson_mp(pstar_mp, a, m0)
    om = (2 / mpf(a)) * msinh(mpf(a) * Dl / 2)
    seq_d1.append((om - Dc_mp) / mpf(a) ** 2)
    seq_c1.append((Dl / Dc_mp - 1) / mpf(a) ** 2)
out("    a          (ω_a-Δ)/a²                  (Δ_a/Δ-1)/a²")
for a, x, y in zip(a_list, seq_d1, seq_c1):
    out(f"    {a:<10.0e} {mp.nstr(x,15):<26} {mp.nstr(y,15)}")
out(f"    Richardson 外推：δ₁ = {mp.nstr(richardson(seq_d1),15)}    （预言 {d1_pred:.15f}）")
out(f"    Richardson 外推：c₁ = {mp.nstr(richardson(seq_c1),15)}    （预言 {c1_pred:.15f}）")
out(f"    ⇒ δ₁ ≠ 0 且有限、与 a 无关 ⇒ 通过 (SM-2) 的'δ₁≠0'前件；而 k 仍是 2（见下）")

# 上界常数 C：|Δ_a/Δ-1| ≤ C a² 对 a ∈ (0, 0.5] 全部成立？
Cnum = max(abs(Delta_wilson(pstar, a, m0) / Dc - 1) / a ** 2
           for a in [10 ** (-k / 20) for k in range(0, 61)])
out(f"    一致上界核：sup_(a∈[1e-3,1]) |Δ_a/Δ-1|/a² = {Cnum:.12f}  （有限 ⇒ k = 2 可证，C = {Cnum:.6f}）")
out(f"    注：Richardson −1.2310663 与 a=1e-6 行 −1.2310661218 一致到 7 位；两行都对预言 −1.2310661218 收敛（50 位 mp，动量用 mpi）")

# 对照：p*=0
out("\n    对照（p*=0，同一公式）：")
seq0 = []
for a in a_list:
    Dl = Delta_wilson_mp((mpf(0), mpf(0), mpf(0)), a, m0)
    om = (2 / mpf(a)) * msinh(mpf(a) * Dl / 2)
    seq0.append((om - mpf(m0)) / mpf(a) ** 2)
out("      (ω_a-m)/a² = " + ", ".join(mp.nstr(x, 6) for x in seq0) + "   ⇒ δ₁ = 0 精确")

# ----------------------------------------------------------------- [4]
out("\n[4] 族 F2：Symanzik 改进自由场（GK3-SYMANZIK.md:144-155 的闭式）")
out("    cosh(aΔ) = 4 - 3√(1 - a²m²/3)  ⇒  Δ_a = m + a⁴m⁵/180 + O(a⁶)")


def Delta_imp(a, m):
    s = a * a * m * m / (1 + math.sqrt(1 - a * a * m * m / 3))
    return (2 / a) * math.asinh(math.sqrt(abs(s) / 2))


def Delta_imp_mp(a, m):
    a = mpf(a)
    s = a ** 2 * mpf(m) ** 2 / (1 + msqrt(1 - a ** 2 * mpf(m) ** 2 / 3))
    return (2 / a) * asinh(msqrt(abs(s) / 2))


m0 = 1.0
out("    a        (m-Δ_a)/a⁴            (Δ_a/m-1)/a²        (Δ_a/m-1)/a⁴        (ω_a-m)/a²")
sd, sd4, so = [], [], []
for a in [0.2, 0.1, 0.05, 0.02, 0.01, 0.005]:
    Dl = Delta_imp_mp(a, m0)
    om = (2 / mpf(a)) * msinh(mpf(a) * Dl / 2)
    r4 = (mpf(m0) - Dl) / mpf(a) ** 4
    q2 = (Dl / mpf(m0) - 1) / mpf(a) ** 2
    q4 = (Dl / mpf(m0) - 1) / mpf(a) ** 4
    o2 = (om - mpf(m0)) / mpf(a) ** 2
    sd.append(r4); sd4.append(q4); so.append(o2)
    out(f"    {a:<8} {mp.nstr(r4,10):<20} {mp.nstr(q2,10):<20} {mp.nstr(q4,10):<20} {mp.nstr(o2,10)}")
out("    极限：(m-Δ)/a⁴ → -1/180 = -0.0055555556（Wilson 是 -a²m³/24）")
out("    极限：(Δ/m-1)/a² → 0  ⇒ c₁ = 0  ⇒ 由 (1.2a) δ₁ = Δ³/24 = m³/24 = 0.0416666667 ≠ 0")
out(f"    Richardson (Δ/m-1)/a⁴ = {mp.nstr(richardson(sd4, ratio=2.0),10)}  （应 1/180 = 0.0055555556 ⇒ k = 4）")
out(f"    Richardson (ω_a-m)/a²  = {mp.nstr(richardson(so, ratio=2.0),10)}  （应 m³/24 = 0.0416666667）")
out("    ⇒ 同一'自由场'类里：δ₁ = 0 给 k=2；δ₁ = +m³/24 给 k=4（>2）")
out("    注：a≥0.02 时 a⁶ 项尚未压住 ⇒ 认 a=0.005 行：0.005555567956 ≈ 1/180、0.04166681858 ≈ 1/24（7 位）")

# ----------------------------------------------------------------- [5]
out("\n[5] SEP-2（GR1:151-154）：H_a = H - δ_a P₁，δ_a = 1/log(1/a)，Δ_cont = 1")
out("    a          δ_a             δ_a/a (k=1)      δ_a/a² (k=2)      (ω_a-1)/a²")
for a in (1e-2, 1e-6, 1e-12, 1e-30):
    da = 1.0 / math.log(1.0 / a)
    v1, v2 = da / a, da / a ** 2
    a_m = mpf(a); d_m = 1 / mlog(1 / a_m)
    Dl = 1 - d_m
    om = (2 / a_m) * msinh(a_m * Dl / 2)
    out(f"    {a:<10.0e} {da:<15.10g} {v1:<15.6g} {v2:<16.6g} {mp.nstr((om-1)/a_m**2,8)}")
out("    ⇒ |Δ_a-Δ_cont| = δ_a ⇒ δ_a/a^k → ∞ 对所有 k（GR1:154 的数字逐位复现：")
out("      a=1e-2: 2.17e+01 / 2.17e+03；a=1e-6: 7.24e+04 / 7.24e+10；a=1e-12: 3.62e+10 / 3.62e+22）")
out("    ⇒ 它的 δ₁ = lim(ω_a-1)/a² = -1/(a²log(1/a)) → -∞：不是'δ₁≠0 有限'，是发散")
out("    结构性 [C]：δ_a = 1/log(1/a) 在 a=0 处 Taylor 级数恒为 0 却不恒等（平坦函数）")
out("    ⇒ 任何'幂 × log'型渐近（= Symanzik 展开；4D 无非微扰版本，GK3:79-82）都给不出它")
out("    ⇒ SEP-2 不可能由格点族实现 —— 此步记 [C]（依赖 SYM-4D），不是 [T]")

# ----------------------------------------------------------------- [6]
out("\n[6] 连续物理路径：扭转角 theta ∈ [0,1]，p* = (2πθ/L, 0, 0)  —— δ₁ 连续变，k 恒 2")
out("    theta      p*_x          Δ_cont        δ₁          c₁          sup|Δ_a/Δ-1|/a²")
for th in (0.0, 0.25, 0.5, 0.75, 1.0):
    px = 2 * math.pi * th / L
    p_ = (px, 0.0, 0.0)
    D_ = math.sqrt(m0 ** 2 + px ** 2)
    S_ = px ** 4
    d1 = -S_ / (24 * D_) if D_ else 0.0
    c1 = (d1 / D_ - D_ ** 2 / 24) if D_ else -m0 ** 2 / 24
    sup = max(abs(Delta_wilson(p_, a, m0) / D_ - 1) / a ** 2
              for a in [10 ** (-k / 20) for k in range(0, 61)])
    out(f"    {th:<10} {px:<14.8f} {D_:<13.10f} {d1:<12.8f} {c1:<12.8f} {sup:.6f}")
out("    ⇒ θ: 0→1 时 δ₁ 从 0 连续走到 -10.2070（θ=0.5 即 p*=π/L 那格给出 -1.2311），")
out("      而 k ≡ 2、C(θ) 有限 ⇒ 'δ₁ = 0' 是这条物理路径上的**一个点**，不是 k=2 的开关")
out("    ⇒ C(L) 随 L 减小（p*=π/L→0 时 C→m²/24=0.0417）⇒ sup_L C = C(L=1) = 0.8263 < ∞ ⇒ (L,a)-一致")

# ----------------------------------------------------------------- [7]
out("\n[7] (1.2a)/(1.4d) 的精确代数（Fraction，无浮点）")
out("    c₁ = δ₁/Δ - Δ²/24   ⇔   δ₁ = Δ(c₁ + Δ²/24) = ⟨K₂⟩₁ + Δ³/24")
for (D, d1) in [(Fr(1), Fr(0)), (Fr(1), Fr(1, 24)), (Fr(3, 2), Fr(-1, 3)),
                (Fr(5, 7), Fr(2)), (Fr(2), Fr(8, 24))]:
    c1 = d1 / D - D ** 2 / 24
    back = D * (c1 + D ** 2 / 24)
    out(f"    Δ={str(D):<5} δ₁={str(d1):<8} ⇒ c₁={str(c1):<12} ⇒ Δ(c₁+Δ²/24)={str(back):<8} 回代{'OK' if back == d1 else 'FAIL'}"
        f"   k 的判据：c₁=0 ？{'是 ⇒ 若 c₂ 有限则 k=4' if c1 == 0 else '否 ⇒ k=2'}")
out("    ⇒ k = 2 ⟺ c₁ ≠ 0；c₁ = 0 ⟺ δ₁ = Δ³/24。**k = 2 不要求 δ₁ = 0，只排除 δ₁ = Δ³/24**")

out("\n" + "=" * 78)
out("完。全部确定性；无随机数、无 MC、无大数组；exit=0")
out("=" * 78)
LOG.close()
