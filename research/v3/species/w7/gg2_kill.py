#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gg2_kill.py -- GG2-KILL 见证脚本
靶: 判 GF2 的"形状盲代价"是否杀死 GF3 的 mu_eff <= e 这个 bit
内容:
  P1  sigma = beta(1-cos(2pi/N)) 逐值复算 (beta=2.5)
  P2  独立 DP 重算 a_n(D) / a_n(S1) 精确整数 (n<=24), 与 GF3-ENTROPY 表逐位比对
  P3  x_c 与 mu 的精确/高精读数 (闭式 + 二分 + 比值/Richardson)
  P4  加权增长率恒等式: b_n = a_n * y^n,  y = e^{-sigma}  =>  lim b_n/b_{n-1} = y*mu
      并给出错式 sigma*mu_class 的数值 (对照)
  P5  闸门表: 真判据 sigma >= kappa + ln(1+2PR) + gamma_class  (P=R=1, kappa=0.5)
      对每个候选 sigma 给 gamma_max / mu_max, 并逐类判过/不过
  P6  判决窗口宽度: sigma_req(S1) - sigma_req(D) = ln(mu(S1)/mu(D))
无随机、无 MC; 纯整数 DP + 浮点仅作读数打印。
"""
import math

print("=" * 84)
print("P1  sigma = beta*(1-cos(2*pi/N))  @ beta = 2.5   (GE2 的界面率, ge2_defdim_witness.py:32)")
print("=" * 84)
BETA = 2.5
for N in (2, 3, 4, 5, 6):
    s = BETA * (1.0 - math.cos(2.0 * math.pi / N))
    print(f"  N={N}:  1-cos(2pi/N) = {1-math.cos(2*math.pi/N):.6f}   sigma = {s:.6f}")
print("  精确: N=2 -> 2beta = 5.00 ; N=3 -> 1.5beta = 3.75 ; N=4 -> beta = 2.50 ;")
print("        N=6 -> 0.5beta = 1.25 ; N=5 -> 2.5*(1-cos72deg) = 1.7274577259...")

print()
print("=" * 84)
print("P2  独立 DP: a_n 精确整数 (轮廓自动机, 约束 h_j >= max(1, d_j+1, 1-d_{j-1}))")
print("=" * 84)

def class_counts(deltas, nmax):
    """a[n] = 类内 n 片簇的个数 (模平移)。deltas = 允许的界面位移集合 (含 0)。"""
    D = sorted(deltas)
    F = {}                       # F[(n,d)] : 已放 >=1 列, 末位移 = d, 前 j 列约束全满足
    for d in D:                  # 第一列: h1 >= max(1, d+1)
        for h in range(max(1, d + 1), nmax + 1):
            F[(h, d)] = F.get((h, d), 0) + 1
    for n in range(1, nmax):     # 升序扩散 (每次转移 n 严格增大)
        for d in D:
            c = F.get((n, d), 0)
            if not c:
                continue
            for d2 in D:
                L = max(1, d2 + 1, 1 - d)
                for h in range(L, nmax - n + 1):
                    F[(n + h, d2)] = F.get((n + h, d2), 0) + c
    a = [0] * (nmax + 1)
    for n in range(1, nmax + 1):
        a[n] += 1                                   # m = 1 (单列 = 高一列的条)
    for (n, d), c in F.items():                     # 收尾列: h_m >= max(1, 1-d)
        for h in range(max(1, 1 - d), nmax - n + 1):
            a[n + h] += c
    return a

NMAX = 24
DD = list(range(0, -(NMAX + 1), -1))     # 有向: delta <= 0, 无界 (截断到 n<=24 足够)
S1 = [-1, 0, 1]                          # 双侧 |delta| <= 1
aD = class_counts(DD, NMAX)
aS = class_counts(S1, NMAX)

print("  n        :" + "".join(f"{n:>10d}" for n in range(1, 15)))
print("  a_n(D)   :" + "".join(f"{aD[n]:>10d}" for n in range(1, 15)))
print("  a_n(S1)  :" + "".join(f"{aS[n]:>10d}" for n in range(1, 15)))
fib = [0, 1]
while len(fib) <= 2 * NMAX + 2:
    fib.append(fib[-1] + fib[-2])
okD = all(aD[n] == fib[2 * n - 1] for n in range(1, NMAX + 1))
print(f"  D  vs F_(2n-1) 逐位相等 (n<=24): {okD}")
# S1 的独立校验用"分母给出的递推"而不是抄表: P=1-4x+3x^2+x^3-x^5  =>  a_n = 4a_{n-1}-3a_{n-2}-a_{n-3}+a_{n-5}
print(f"  S1: a_n = 4a_(n-1)-3a_(n-2)-a_(n-3)+a_(n-5) 对 n>=6 全成立:",
      all(aS[n] == 4 * aS[n - 1] - 3 * aS[n - 2] - aS[n - 3] + aS[n - 5] for n in range(6, NMAX + 1)))
print(f"  S1: 与 GF3-ENTROPY §2.3 已公布的 a_8=1090, a_14=561346 对齐: {aS[8]==1090 and aS[14]==561346}")
print(f"  D : a_n = 3a_(n-1)-a_(n-2) (x^2-3x+1) 对 n>=3 全成立:",
      all(aD[n] == 3 * aD[n - 1] - aD[n - 2] for n in range(3, NMAX + 1)))

print()
print("=" * 84)
print("P3  x_c 与 mu")
print("=" * 84)
phi2 = (3.0 + math.sqrt(5.0)) / 2.0
xcD = (3.0 - math.sqrt(5.0)) / 2.0
print(f"  D : x_c = (3-sqrt5)/2 = {xcD:.15f} ; mu = phi^2 = {phi2:.15f}  (闭式, GF3 §2.2)")

def poly_val(c, x):
    return sum(ci * x ** i for i, ci in enumerate(c))

P_S1 = [1, -4, 3, 1, 0, -1]          # 1 - 4x + 3x^2 + x^3 - x^5   (GF3 §0④)
lo, hi = 0.0, 1.0
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if poly_val(P_S1, mid) * poly_val(P_S1, lo) > 0:
        lo = mid
    else:
        hi = mid
xcS = 0.5 * (lo + hi)
muS = 1.0 / xcS
print(f"  S1: P = 1-4x+3x^2+x^3-x^5 ; 最小正根 x_c = {xcS:.15f} ; mu = {muS:.15f}")

def ratio_limit(a, k):
    """a_n/a_{n-1} 的 Richardson (2 步) 外推"""
    r = [a[n] / a[n - 1] for n in range(2, len(a))]
    r1 = r[-1]; r2 = r[-2]; r3 = r[-3]
    # 假定 r_n = L + c/n  =>  L ~ 2 r1 - r2 ; 再一阶
    L1 = 2 * r1 - r2
    L2 = 2 * r2 - r3
    return r1, (2 * L1 - L2)
for name, a, mu in (("D", aD, phi2), ("S1", aS, muS)):
    r1, ext = ratio_limit(a, NMAX)
    print(f"  {name:2s}: a_N/a_(N-1) = {r1:.9f} ; Richardson = {ext:.9f} ; 精确 mu = {mu:.9f} ; 差 {ext-mu:+.2e}")

print()
print("=" * 84)
print("P4  加权增长率:  b_n = a_n * y^n ,  y = e^{-sigma}")
print("=" * 84)
print("  恒等式(换元 x -> x/y):  G(y x) 的收敛半径 = x_c/y  =>  增长率 = y * mu_class   [T]")
print("  错式(任务书 rho = sigma*mu_class) 与正解的数值对照:")
sig_vals = [("2.598612 (=0.5+1+ln3, gamma_ent=1 窗)", 2.598612288668),
            ("5.00 (=2beta, beta=2.5 SU(2) Z_2 twist)", 5.0)]
for lab, sig in sig_vals:
    y = math.exp(-sig)
    print(f"\n  sigma = {lab}")
    print(f"    y = e^(-sigma) = {y:.9f}")
    print(f"    mu_w(D)  = y*phi^2      = {y*phi2:.9f}     (错式 sigma*phi^2 = {sig*phi2:.6f})")
    print(f"    mu_w(S1) = y*mu(S1)     = {y*muS:.9f}      (错式 sigma*2.8312 = {sig*muS:.6f})")
    print(f"    mu_w/mu_class = y = {y:.9f}  (两类同一个因子 => 排序不变)  [T]")
print("\n  排序不变性的数值核: mu(S1)/mu(D) 是常数:")
print(f"    mu(S1)/mu(D) = {muS/phi2:.9f}  (与 sigma 无关)  [T]")

print()
print("=" * 84)
print("P5  真判据: sigma >= kappa + ln(1+2PR) + gamma_class   (P=R=1 => ln3=1.098612, kappa=0.5)")
print("      <=>  r := e^{kappa - sigma} * mu_class <= 1/3   [BianchiLemma.md:109-111]")
print("=" * 84)
KAPPA = 0.5
LN3 = math.log(3.0)
C0 = KAPPA + LN3               # 1.598612
print(f"  kappa = {KAPPA}, ln3 = {LN3:.9f}, kappa+ln3 = {C0:.9f}")
print(f"  检: 0.5+1.0+ln3 = {KAPPA+1.0+LN3:.12f}   (GD2-WSTAR.md:265 的 2.598612)")
classes = [
    ("bargraph D0      ", 2.0),
    ("D (有向, GF3 类) ", phi2),
    ("S1 (双侧 1 步)   ", muS),
    ("S4 (列凸下界)    ", 3.190001268135846),
    ("平面 polyomino 硬界", 2.927334013684),
    ("平面 polyomino 文献", 4.0626),
    ("动物 a_n 下界   ", 13.902),
    ("动物 B_n 下界   ", 15.58),
    ("动物 中位 B     ", 31.773),
]
print("\n  sigma_req(class) = kappa + ln3 + ln mu_class :")
for lab, mu in classes:
    g = math.log(mu)
    print(f"    {lab:22s} mu = {mu:9.6f}  gamma = {g:.6f}  sigma_req = {g + C0:.6f}")
cands = [
    ("实测张力锚 c0*q   (c0=1/3.405^2, q=6.3)", 0.086251 * 6.3),
    ("语料工作值 c0*q   (c0=1, q=1)", 1.0),
    ("实测锚 二次式 c0*q^2", 0.086251 * 6.3 ** 2),
    ("SU(3) 1.5beta     ", 3.75),
    ("SU(2) 2beta       ", 5.0),
    ("gamma_ent=1 窗 2.598612", 0.5 + 1.0 + LN3),
]
print("\n  每个候选 sigma 的闸门 (gamma_max = sigma - (kappa+ln3), mu_max = e^gamma_max):")
for lab, sig in cands:
    gm = sig - C0
    mm = math.exp(gm)
    passing = [l.strip() for l, m in classes if math.log(m) <= gm]
    print(f"    sigma = {sig:9.6f}  [{lab}]")
    print(f"      gamma_max = {gm:+9.6f} ;  mu_max = {mm:10.6f} ;  r(D) = {math.exp(KAPPA-sig)*phi2:.6f} (<=1/3 ? {math.exp(KAPPA-sig)*phi2 <= 1/3})")
    print(f"      过闸的类 ({len(passing)}/{len(classes)}): {', '.join(passing) if passing else '(空)'}")

print()
print("=" * 84)
print("P6  判决窗口 (bit 唯一有判别力的 sigma 区间)")
print("=" * 84)
sD = math.log(phi2) + C0
sS = math.log(muS) + C0
print(f"  sigma_req(D)  = {sD:.9f}")
print(f"  sigma_req(S1) = {sS:.9f}")
print(f"  Delta_sigma   = ln(mu(S1)/mu(D)) = {sS - sD:.9f}   (窗口宽)")
print(f"  Delta_sigma / sigma_req(D) = {(sS-sD)/sD:.4%} ;  Delta_sigma / (kappa+ln3) = {(sS-sD)/C0:.4%}")
print(f"  2.598612 落在此窗口内? {sD < 0.5+1.0+LN3 < sS}")
for lab, sig in cands:
    inside = sD < sig < sS
    print(f"    候选 sigma = {sig:9.6f}  在窗口内? {inside}   距窗口中心 {sig - 0.5*(sD+sS):+9.6f}")
print("\n  同一窗口的另一种读数 (gamma 上):")
print(f"    过/不过 的翻转在 gamma_class = sigma - (kappa+ln3):")
for lab, sig in cands:
    print(f"      sigma = {sig:9.6f}  ->  gamma_class 边界 = {sig - C0:+9.6f}  (mu = {math.exp(sig-C0):10.4f})")
print("\nPASS: 全部为确定性计算 (无随机).")

print()
print("=" * 84)
print("P7  'mu_eff <= e' 的两种读法 (判决点)")
print("    读法(i)  mu_eff = y*mu_class (raw 加权率)  =>  阈值 sigma_raw = ln(mu_class/e)")
print("    读法(ii) mu_eff = mu_class                =>  只在 sigma = 2.598612 成立")
print("=" * 84)
print("  读法(i): 每个类掉到 e 以下所需的 sigma_raw = ln(mu/e):")
for lab, mu in classes:
    print(f"    {lab:22s} mu = {mu:9.6f}  sigma_raw = {math.log(mu/math.e):+9.6f}"
          f"   (在 sigma=2.598612 处过? {math.exp(-0.5-1.0-LN3)*mu <= math.e})   (在 sigma=5.0 处过? {math.exp(-5.0)*mu <= math.e})")
print("\n  => 读法(i) 在 sigma >= 2.4576 (=ln(31.773/e)) 时对**全部 9 类**成立 => 判据恒真、无内容")
print("  => 读法(ii) 是 sigma=2.598612 的切片 => 依赖 gamma_ent=1 (GC2-H23.md:20-28 已作废该估计)")
print("\n  真判据 (BL:109-111) 的余量读数: r = e^{kappa-sigma}*mu_class <= 1/3")
for lab, sig in cands:
    rD = math.exp(KAPPA - sig) * phi2
    rS = math.exp(KAPPA - sig) * muS
    print(f"    sigma = {sig:9.6f} :  r(D) = {rD:9.6f}   r(S1) = {rS:9.6f}   (r<=1/3 => {rD<=1/3}, {rS<=1/3})")
print("\nDONE.")

