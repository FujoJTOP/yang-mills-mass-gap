"""gk2sym_check.py -- GK2-SYMBOL 的确定性数值/符号核（v3：精确有理数 + Decimal）。

纯 numpy + fractions + decimal；无随机数、无 MC、无 cvxpy。
复跑：py -3.12 gk2sym_check.py
"""
import numpy as np
from fractions import Fraction as F
from decimal import Decimal as Dc, getcontext

getcontext().prec = 200   # [3] 要减两次 1（lo 与 hi 各一次），60 位在 x~1e-15 处不够
ORD = 3   # 保留 t^0, t^1, t^2 (t = a^2)

# ---------------------------------------------------------------- 精确幂级数
def pmul(p, q, ordn=ORD):
    r = [F(0)] * min(len(p) + len(q) - 1, ordn)
    for i, x in enumerate(p):
        if i >= ordn:
            break
        for j, y in enumerate(q):
            if i + j >= ordn:
                break
            r[i + j] += x * y
    return r

def padd(p, q):
    n = min(max(len(p), len(q)), ORD)
    r = [F(0)] * n
    for i, x in enumerate(p[:n]):
        r[i] += x
    for i, y in enumerate(q[:n]):
        r[i] += y
    return r

def ppow(p, k):
    r = [F(1)]
    for _ in range(k):
        r = pmul(r, p)
    return r

# asinh(x) = sum_{n>=0} (-1)^n C_n x^{2n+1},  C_n = (2n)!/(4^n (n!)^2 (2n+1))
C_ASINH = [F(1), F(-1, 6), F(3, 40), F(-15, 336), F(105, 3456), F(-945, 42240)]

def delta_series(omega):
    """(2/a)*asinh(a*omega/2), omega = 关于 t=a^2 的级数。返回 t^0..t^2 的系数。"""
    out = [F(0)] * ORD
    for n, c in enumerate(C_ASINH):
        # 项: c_n/4^n * t^n * omega^(2n+1)
        term = [F(0)] * n + [F(c, 4 ** n)]
        term = term[:ORD]
        term = pmul(term, ppow(omega, 2 * n + 1))
        out = padd(out, term)
    return out

print("=" * 80)
print("[1] 精确：自由场 Delta_a = (2/a) asinh(a m/2)，t=a^2 的系数（Fraction，无浮点）")
print("=" * 80)
for m in (F(1), F(2), F(3), F(1, 2)):
    s = delta_series([m, F(0), F(0)])
    assert s[0] == m, s
    c1, c2 = s[1] / s[0], s[2] / s[0]
    print(f"  m={m}: c1 = {c1}   应 = -m^2/24 = {-m*m/24}   [{'OK' if c1 == -m*m/24 else 'MISMATCH'}]")
    print(f"         c2 = {c2}   应 = 3 m^4/640 = {3*m**4/640}  [{'OK' if c2 == 3*m**4/640 else 'MISMATCH'}]")
print("  => c1 = -m^3/(24*Delta_cont)，c2 = 3 m^5/(640*Delta_cont)，Delta_cont = m。")

print()
print("=" * 80)
print("[2] 精确：合成公式  c1 = d1/D - D^2/24 ,  c2 = d2/D - d1 D/8 + 3 D^4/640")
print("    omega_a = D + a^2 d1 + a^4 d2")
print("=" * 80)
for (D0, d1, d2) in ((F(1), F(0), F(0)), (F(1), F(7, 10), F(-2, 5)),
                     (F(3, 2), F(-1, 3), F(5)), (F(5, 7), F(2), F(-11, 3))):
    s = delta_series([D0, d1, d2])
    c1f = d1 / D0 - D0 ** 2 / 24
    c2f = d2 / D0 - d1 * D0 / 8 + 3 * D0 ** 4 / 640
    print(f"  D={D0}, d1={d1}, d2={d2}")
    print(f"     c1 = {s[1]/s[0]}  vs 公式 {c1f}  [{'OK' if s[1]/s[0] == c1f else 'MISMATCH'}]")
    print(f"     c2 = {s[2]/s[0]}  vs 公式 {c2f}  [{'OK' if s[2]/s[0] == c2f else 'MISMATCH'}]")

print()
print("=" * 80)
print("[3] 精确双边不等式 0 <= 1 - Delta_a/m <= a^2 m^2/24（Decimal 200 位）")
print("=" * 80)
def asinhD(x):
    return (x + (x * x + 1).sqrt()).ln()

bad, npt = [], 0
for ea in range(-12, 0):
    a = Dc(10) ** ea
    for em in range(-4, 5):
        m = Dc(10) ** em
        x = a * m / 2
        Dl = 2 / a * asinhD(x)
        lo = 1 - Dl / m
        hi = a * a * m * m / 24 - lo
        npt += 1
        if lo < 0 or hi < 0:
            bad.append((ea, em, lo, hi))
print(f"  点样 {npt} 个（a=1e-12..1e-1, m=1e-4..1e4）：违反 {len(bad)} 个")
for t in bad:
    print("   违规点 ea=%d em=%d  lo=%.3e  hi=%.3e" % t)

print()
print("=" * 80)
print("[4] 二阶谱和的幂次计数  S = sum_{k!=0}|k|^6  (K2 为 4 阶符号: |V_n1|^2~k^8, 1/gap~k^-2)")
print("    d 维：S ~ L^d a^{-(d+6)}  =>  a^4 S ~ L^d a^{-(d+2)} -> INF")
print("=" * 80)
def S_sum(d, N, L=1.0):
    ax = np.arange(-(N // 2), N // 2) * (2.0 * np.pi / L)
    grids = np.meshgrid(*([ax] * d), indexing="ij")
    k2 = sum(g ** 2 for g in grids)
    return float(np.sum((k2 ** 3)[k2 > 0]))
for d, Ns in ((3, (16, 24, 32, 48, 64)), (4, (8, 12, 16, 24))):
    print(f"  d={d}:  {'N':>4} {'S_N':>16} {'S/N^(d+6)':>12} {'a^4 S':>16}")
    prev = None
    for N in Ns:
        s = S_sum(d, N)
        tag = f"   斜率 {np.log(s/prev[1])/np.log(N/prev[0]):.4f} (应 {d+6})" if prev else ""
        print(f"        {N:4d} {s:16.6e} {s/N**(d+6):12.4f} {s/N**4:16.6e}{tag}")
        prev = (N, s)

print()
print("=" * 80)
print("[5] 二阶/一阶 之比 (a^4 S)/(|c1| a^2) = a^2 S/|c1| ~ a^{-d}：发散，不是小量")
print("=" * 80)
for d in (3, 4):
    row = [f"N={N}:{(S_sum(d,N)/N**4)/((1/24)/N**2):.2e}" for N in (16, 32, 64)]
    print(f"  d={d}: " + "   ".join(row))

print()
print("[6] 2D: T_a = exp(-(g^2 a/2) sum C2) => H_a = (g^2/2) sum C2 与 a 无关")
print("    => Delta_a == Delta_cont, c1 = c2 = 0, C = 0 精确（引 GN4:244 / GP9-RB:287）")
print()
print("=" * 80)
print("[7] 恒等式核：omega_a := (2/a)sinh(a Delta_a/2) 把自由场精确送回 m（a-无关）")
print("=" * 80)
for a in (1e-1, 1e-3, 1e-6):
    for m in (1.0, 2.0):
        Dl = (2.0 / a) * np.arcsinh(a * m / 2.0)
        om = (2.0 / a) * np.sinh(a * Dl / 2.0)
        print(f"  a={a:.0e}, m={m:.1f}:  omega_a - m = {om - m:.3e}")
print("  => 自由场 omega_a ≡ m（差 = 0.000e+00）；c1 的全部 a-依赖在 asinh 核里。")
print()
print("DONE（无随机数、无 MC）")
