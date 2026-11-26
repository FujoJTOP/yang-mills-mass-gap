# -*- coding: utf-8 -*-
"""TN2 — lambda_T 的逐字定义核验 + 2.4 的结构意义 + 动物常数真值探针。

只读语料常数（P3-FULLREP-KP:89/116/121、ADJ2、KP6），全部算术本脚本独立复算。
py -3.12 + mpmath / decimal(prec=50)。
"""
import itertools, sys, time
from decimal import Decimal, getcontext
import mpmath as mp

getcontext().prec = 50
mp.mp.dps = 50

def H(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)

# ---------------------------------------------------------------- 1. 常数
H("1 · lambda_T 的常数：e*Delta 与 48 的两种写法")
E = mp.e
d = mp.mpf(4)
Delta = 4 * (2 * d - 3)          # P3-FULLREP-KP:116  Δ = 4(6-1) = 20, 6 = 2(d-1)
lam_T = E * Delta                # P3-FULLREP-KP:116  λ_a = eΔ = 20e
SZZ_den = 16 * (d - 1)           # SZZ: |β| < 1/(16(d-1))  ⇒ d=4: 1/48
print(f"Delta = 4(2d-3)          = {mp.nstr(Delta,20)}   (d=4)")
print(f"lambda_T = e*Delta       = {mp.nstr(lam_T,25)}")
print(f"16(d-1) = 48             = {mp.nstr(SZZ_den,20)}")
print(f"48/Delta                 = {mp.nstr(SZZ_den/Delta,25)}")
print(f"e*Delta/48               = {mp.nstr(lam_T/SZZ_den,25)}")
print(f"e/(48/Delta)             = {mp.nstr(E/(SZZ_den/Delta),25)}   <- e/2.4")
print(f"1/lambda_T               = {mp.nstr(1/lam_T,25)}  = beta_KP")
print(f"1/48                     = {mp.nstr(1/SZZ_den,25)}  = beta_SZZ")
print(f"beta_SZZ/beta_KP         = {mp.nstr((1/SZZ_den)/(1/lam_T),25)}")
print(f"lambda_T - 48            = {mp.nstr(lam_T-48,20)}   (>0  ⇒ lambda_T<48 假)")

H("2 · 2.4 的 d-结构：48/Delta = 4(d-1)/(2d-3)，12/5 只是 d=4 的值")
print(" d | Delta=4(2d-3) | 16(d-1) | 48/Delta | e*Delta/48 (缺口)")
for dd in range(2, 9):
    D = 4 * (2 * dd - 3); S = 16 * (dd - 1)
    print(f" {dd} | {D:12d} | {S:7d} | {mp.nstr(mp.mpf(S)/D,10):>8} | {mp.nstr(E*D/S,12)}")  # noqa
# 缺口 = 1 的实根： e(2d-3) = 4(d-1)
dstar = (3 * E - 4) / (2 * E - 4)
print(f"\ng(d)=e*Delta/48=e(2d-3)/(4(d-1));  g(d*)=1 ⇒ d* = (3e-4)/(2e-4) = {mp.nstr(dstar,20)}")
print(f"g'(d) ∝ 1/(d-1)^2 > 0 ⇒ 缺口随 d 单调升；整数 d≥3 全部 >1")
for dd in (3, 4):
    D = 4 * (2 * dd - 3); S = 16 * (dd - 1)
    print(f"  d={dd}: 缺口 = {mp.nstr(E*D/S,12)}  (需把 lambda_T 压 "
          f"{mp.nstr((1 - S/(E*D))*100, 6)}%)")

H("3 · 12 的两个来源（去混淆）：12/5 的 12 与 '立方体 12 边' 的 12")
print("(i)  12/5 的 12 = 48/4  ← gcd(48,20)=4，48=16(d-1)，20=4(2d-3) ⇒ 4(d-1)/ (2d-3)")
print("(ii) 语料的 12 = 立方体盒 K_{2,2,2} 的 12 条边（UR3、KP6 §1，进 Ĉ=2^{12/6}=4）")
print("(iii) 12 还是 P3:116 里被更正掉的旧错值 Δ_old=12（'degree ≤12' → 错，Δ=20）")
print(f"用旧 Δ=12 时的比值： 48/12 = {mp.nstr(mp.mpf(48)/12,10)}  (≠2.4)；Δ=20 才给 2.4")
print("⇒ 三处 12 数值相撞，但在算式里分别属于 Ĉ 槽 / 被撤的 Δ 槽 / gcd 约分，互不连通")

H("4 · 动物常数真值探针：T^4_L 的 plaquette 邻接图上过定点的连通动物数")
def build(L):
    plaq, link = [], {}
    for mu in range(4):
        for nu in range(mu + 1, 4):
            for x in itertools.product(range(L), repeat=4):
                lid = len(plaq); plaq.append([])
                e = [0, 0, 0, 0]; e[mu] = 1
                emu = tuple((x[i] + e[i]) % L for i in range(4))
                e = [0, 0, 0, 0]; e[nu] = 1
                enu = tuple((x[i] + e[i]) % L for i in range(4))
                for key in ((mu, x), (nu, x), (mu, enu), (nu, emu)):
                    link.setdefault(key, []).append(lid)
    adj = [set() for _ in plaq]
    for ids in link.values():
        for a, b in itertools.combinations(ids, 2):
            adj[a].add(b); adj[b].add(a)
    return plaq, adj

for L in (6,):
    t0 = time.time()
    plaq, adj = build(L)
    deg = sorted({len(s) for s in adj})
    print(f"L={L}: plaquettes={len(plaq)}  邻接度集合={deg}  (P3:116 期望 {{20}})")
    root = 0
    levels = [[frozenset([root])]]
    counts = [1]
    n_max = 5
    for n in range(2, n_max + 1):
        seen = set()
        for S in levels[-1]:
            for v in S:
                for w in adj[v]:
                    if w not in S:
                        T = frozenset(S | {w})
                        seen.add(T)
        levels.append(list(seen)); counts.append(len(seen))
        print(f"  a_{n} = {len(seen):>9d}    ratio a_n/a_(n-1) = "
              f"{len(seen)/counts[-2]:.4f}    a_n^(1/n) = {len(seen)**(1.0/n):.4f}")
    print(f"  (ADJ2 语料 n≤4: 1, 20, 458, 11132   →  本脚本 {counts[:4]})")
    print(f"  运行 {time.time()-t0:.1f}s")
    if len(counts) >= 5:
        r = [counts[i] / counts[i - 1] for i in range(1, len(counts))]
        print(f"  相邻比 {['%.3f' % v for v in r]}  升 ⇒ 生长常数 λ_* ≳ {r[-1]:.1f}")
        print(f"  λ_*/eΔ ≈ {r[-1]/float(lam_T):.3f}   ⇒ 若 λ_*≈{r[-1]:.0f}，则 1/λ_* = {1/r[-1]:.5f}"
              f" > 1/48 = {1/48:.5f}")
        print(f"  严格包络（P3:89）: a_n ≤ (eΔ)^n  —— 这是**上界**，非真值：a_n^(1/n)/eΔ at n=5 "
              f"= {counts[-1]**(1.0/5)/float(lam_T):.4f}")

H("5 · 判死线落点")
print(f"lambda_T := λ_a := eΔ = {mp.nstr(lam_T,12)}  （定义，非真值；P3:116/89 逐字）")
print(f"lambda_T < 48 ?  {lam_T < 48}   缺口 = {mp.nstr(lam_T/48,12)} = e/(48/Δ)")
print("真值层：eΔ 是 Penrose 型**包络**（P3:89 'one may take λ_a = eΔ'、P3:121 'precise animal")
print("        constant is immaterial and only rescales β₁'、P3:128 'Not optimized'）")
print("       ⇒ (a) 的方向成立（eΔ 只是上界），但 λ_* ≤ 48 目前在语料内**无定理**")
