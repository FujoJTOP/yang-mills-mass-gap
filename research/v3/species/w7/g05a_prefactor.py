# -*- coding: utf-8 -*-
# g05a_prefactor.py -- 核 G05-BETAOFa.md §1.2/1.3 的两圈前因子：符号与 51/121 来源
# 约定：β(g) = -b0 g^3 - b1 g^5 ,  b0 = 11N/(48π²) = 11N/(3·16π²) ,  b1 = 34N²/(3·(16π²)²)
#        β_W = 2N/g² ,  β_A = 2/(g²N) = β_W/N²
# 两写法：std  aΛ = exp(-1/(2b0g²))·(b0g²)^{-51/121}   |   g05  aΛ = exp(-1/(2b0g²))·(g²/16π²)^{+51/121}
from fractions import Fraction
from mpmath import mp, mpf, exp, log, pi, nstr

mp.dps = 40
PI2 = pi ** 2
EXP = mpf(51) / 121          # b1/(2 b0^2)，待证 N-无关


def b0(N):
    return mpf(11) * N / (mpf(48) * PI2)


def b1(N):
    return mpf(34) * N ** 2 / (mpf(3) * (mpf(16) * PI2) ** 2)


def aL(g2, N, form):
    e = exp(-1 / (2 * b0(N) * g2))
    if form == 'std':
        return e * (b0(N) * g2) ** (-EXP)
    return e * (g2 / (16 * PI2)) ** (+EXP)


def g2_of(beta, N, conv):
    betaW = N * N * beta if conv == 'A' else beta
    return 2 * N / betaW


# ---------------------------------------------------------------- A. 系数与指数
print("=" * 78)
print("A. b0/b1 逐字 + 指数 51/121 的精确来源（有理数，不含 π）")
print("=" * 78)
print("  b0 = 11N/(3·16π²) = 11N/(48π²)      = (11/3)·N·(16π²)^{-1}")
print("  b1 = 34N²/(3·(16π²)²)               = (34/3)·N²·(16π²)^{-2}")
r = Fraction(34, 3) / (2 * (Fraction(11, 3)) ** 2)
print(f"  b1/(2b0²) = (34/3)/(2·(11/3)²) = {Fraction(34,3)}/{2*(Fraction(11,3))**2}"
      f" = 306/726 = {r} = {r.numerator}/{r.denominator}")
print("  逐 N 数值验证（40 位）：")
for N in (1, 2, 3, 4, 5, 8, 10):
    print(f"    N={N:>2}: b0={nstr(b0(N),12)}  b1={nstr(b1(N),12)}  "
          f"b1/(2b0²)={nstr(b1(N)/(2*b0(N)**2),20)}")
print(f"  51/121 = {nstr(EXP,20)}")
print("  => [T] 指数对纯 SU(N) 与 N 无关（b0∝N, b1∝N²）。仅对纯 YM 成立；加费米子后 b1 变。")

# ---------------------------------------------------------------- B. 符号：RG 不变性
print()
print("=" * 78)
print("B. 符号判：把 Λ=μ·exp(-Φ) 沿两圈 RG 流求 dlnΛ/dlnμ（闭合形式）")
print("=" * 78)
print("  由 ∫dg/β(g)： lnΛ = lnμ - 1/(2b0g²) - (b1/(2b0²))·ln(g²) + const")
print("  即 Φ(g) = 1/(2b0g²) + s·EXP·ln(g²)；s=+1 <=> 前因子 (g²)^{-EXP} 即 (b0g²)^{-51/121}（std）")
print("                                              s=-1 <=> 前因子 (g²)^{+EXP}（g05）")
print("  dlnΛ/dlnμ = 1 - Φ'(g)·β(g)")
for s, name in ((+1, 'std  (b0g²)^{-51/121}'), (-1, 'g05  (g²/16π²)^{+51/121}')):
    print(f"\n  {name}")
    for N in (2, 3):
        g2v = mpf('0.5')
        g = mp.sqrt(g2v)
        Phi_p = -1 / (b0(N) * g ** 3) + 2 * s * EXP / g
        beta = -b0(N) * g ** 3 - b1(N) * g ** 5
        resid = 1 - Phi_p * beta
        pred = (b1(N) / b0(N)) ** 2 * g ** 4 if s > 0 else \
               -2 * (b1(N) / b0(N)) * g ** 2 - (b1(N) / b0(N)) ** 2 * g ** 4
        print(f"    N={N}, g²=0.5: resid={nstr(resid,10):>16}   closed-form={nstr(pred,10):>16}")
print("  [T] s=+1(std) 的残差是 +(b1/b0)²g⁴ = O(g⁴)=三圈效应（两圈解应有的精度）；")
print("      s=-1(g05) 的残差是 -2(b1/b0)g² -(b1/b0)²g⁴ = O(g²)，与已保留的一圈项同阶")
print("      => 该式在两圈阶**不是** RG 不变式（符号错在代数上就被判死，不依赖数值标定）。")

# ---------------------------------------------------------------- C. ODE 实测漂移
print()
print("=" * 78)
print("C. 独立核验：数值积分两圈 ODE，看两写法的 Λ 是否真的恒定（N=3）")
print("=" * 78)


def beta(g, N):
    return -b0(N) * g ** 3 - b1(N) * g ** 5


def flow(N, g_start, g_end, h=mpf('1e-4')):
    """RK4 解 dt/dg = 1/β(g)（沿 RG 流单调），返回 [(g, μ)]"""
    g, t = mpf(g_start), mpf(0)
    out = [(g, exp(t))]
    sign = 1 if g_end > g_start else -1
    while sign * (g_end - g) > 0:
        step = min(h, abs(g_end - g)) * sign
        k1 = 1 / beta(g, N)
        k2 = 1 / beta(g + step / 2, N)
        k3 = 1 / beta(g + step / 2, N)
        k4 = 1 / beta(g + step, N)
        t += step / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        g += step
        out.append((g, exp(t)))
    return out


for (gs, ge) in ((mpf('1.0'), mpf('0.3')), (mpf('1.0'), mpf('0.1'))):
    traj = flow(3, gs, ge, mpf('2e-4'))
    print(f"\n  N=3, g: {nstr(gs,4)} -> {nstr(ge,4)}   (lnμ 跨度 {nstr(log(traj[-1][1]),8)})")
    for form in ('std', 'g05'):
        vals = [mu * aL(g * g, 3, form) for (g, mu) in traj[::max(1, len(traj)//8)]]
        print(f"    {form:4s}: Λ = {nstr(vals[0],8):>14} ... {nstr(vals[-1],8):>22}"
              f"   max/min = {nstr(max(vals)/min(vals),8)}")
print("  [T] std 漂移 = O(g⁴)（数值上是常量）；g05 的 Λ 随 μ 单调漂移，非 RG 不变。")

# ---------------------------------------------------------------- D. β_W=6 标定（复现 ta04）
print()
print("=" * 78)
print("D. β_W=6.00 SU(3) 标定（复现 ta04_ratio.log:52-55）")
print("=" * 78)
g2 = g2_of(mpf(6), 3, 'W')
print(f"  g²=2N/β_W = 6/6 = {nstr(g2,6)}   b0g² = {nstr(b0(3)*g2,12)}")
aLtrue = mpf('0.263') / mpf('28.809') / mpf('2.91')
print(f"  语料真值 aΛ_lat = 263MeV/28.809/2.91GeV = {nstr(aLtrue,8)}")
for form, tag in (('std', '两圈/真 = 0.748'), ('g05', '两圈/真 = 0.029')):
    v = aL(g2, 3, form)
    print(f"  {form:4s}: aΛ = {nstr(v,8):>14}   两圈/真 = {nstr(v/aLtrue,6):>10}   "
          f"(ta04 log 标 {tag})")
print(f"  两写法之比 std/g05 = {nstr(aL(g2,3,'std')/aL(g2,3,'g05'),8)}"
      f"   = [(16π²)/(b0g⁴)]^{{51/121}} = "
      f"{nstr((16*PI2/(b0(3)*g2**2))**EXP,8)}")
print(f"  0.7479/0.0288 = {nstr(mpf('0.7479')/mpf('0.0288'),8)}  <= 与上式一致 => 0.748 与 0.029 的分歧")
print("  100% 来自前因子符号（指数项两种写法完全相同）。")

# ---------------------------------------------------------------- E. β_A=1/48 表
print()
print("=" * 78)
print("E. β_A=1/48（SZZ 半径）：正确式给 aΛ 与 N 无关，g05 式伪造出 N 分裂")
print("=" * 78)
print(f"{'N':>3} {'β_W':>10} {'g²':>8} {'b0g²':>12} {'std aΛ':>14} {'g05 aΛ':>14} {'g05表值':>10}")
tab = {2: '0.4837', 3: '0.4077'}
for N in (2, 3):
    g2v = g2_of(mpf(1)/48, N, 'A')
    print(f"{N:>3} {nstr(mpf(N*N)/48,6):>10} {nstr(g2v,6):>8} {nstr(b0(N)*g2v,10):>12} "
          f"{nstr(aL(g2v,N,'std'),8):>14} {nstr(aL(g2v,N,'g05'),8):>14} {tab[N]:>10}")
print("  b0g² = 11N/(48π²)·96/N = 22/π² = 2.229066 与 N 无关 => 前因子在 β_A=1/48 处 [T] 严格 N-无关；")
print("  正确 aΛ(1/48) = 0.57000（SU2 = SU3），g05 §1.3 的 0.4837/0.4077 与它的 19% N 分裂是符号伪影。")

# ---------------------------------------------------------------- F. 反解
print()
print("=" * 78)
print("F. 反解 β_W(aΛ=0.01)：正确式的「外推窗」在哪")
print("=" * 78)


def solve_betaW(target, N, form='std'):
    lo, hi = mpf('0.2'), mpf('50')
    for _ in range(200):
        mid = (lo + hi) / 2
        if aL(2 * N / mid, N, form) > target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


for target in (mpf('0.01'), mpf('0.003')):
    for N in (2, 3):
        print(f"  aΛ={nstr(target,4)}  N={N}: β_W(std)={nstr(solve_betaW(target,N,'std'),6):>8}"
              f"   β_W(g05)={nstr(solve_betaW(target,N,'g05'),6):>8}")
print("  (g05 §1.3 逐字给 β_W≈1.12(SU2)/2.39(SU3) —— 复现见右列 g05 值；左列是正确式)")

# 交叉点：指数项与幂次项的对数导数相等处（β_W = 22·(102/121)/(48π²)·... 见下）
print()
print("  「指数压过幂次」的判据：dln(aΛ)/dln(g²) = 1/(2b0g²) - 51/121，翻号处 b0g² = 121/102")
for N in (2, 3):
    g2c = (mpf(121) / 102) / b0(N)
    print(f"    N={N}: 翻号于 g²={nstr(g2c,6)}（β_W={nstr(2*N/g2c,6)}）；"
          f"β_A=1/48 处 dln(aΛ)/dln g² = {nstr(1/(2*b0(N)*g2_of(mpf(1)/48,N,'A')) - EXP,6)}")
print("  => β_A=1/48 处幂次项已压过指数项（|−0.42| > 0.22），g05 §1.4 '指数压过任何幂次修正' 在该点不成立；")
print("     弱耦合区（g²≲1）指数项 7.18 ≫ 0.42，单调性 [T] 仍成立。")

print()
print("=" * 78)
print("H. 副产物：正确符号下 aΛ(β_W) 在强耦合端**非单调**（§1.4 的 [T] 降级）")
print("=" * 78)
print(f"{'β_W':>9} {'g²':>9} {'std aΛ':>12} {'g05 aΛ':>12}  N=3")
for bW in ('0.05', '0.10', '0.1875', '0.30', '0.35', '0.50', '1.0', '2.5', '6.0'):
    g2v = 2 * 3 / mpf(bW)
    print(f"{bW:>9} {nstr(g2v,6):>9} {nstr(aL(g2v,3,'std'),8):>12} "
          f"{nstr(aL(g2v,3,'g05'),8):>12}")
print("  正确式：aΛ 从 β_W→0 的 0 升到 β_W≈0.35 的峰 0.607，再随 β_W 降到 6 处的 2.3e-3。")
print("  => {β<1/48}={β_W<0.1875} 整段落在**上升支**：截断式内 a(β) 不是降的。")
print("  => g05 §1.4 的 '[T] 严格降' 依赖加号（其值随 β 全程单调降）；改对符号后该 [T] 只是 [C]")
print("     （真实 a(β) 的单调性 = AF 的物理内容，不是两圈式给的）。§3.2 的 no-go 不用此条，不受影响。")
print()
print("=" * 78)
print("G. 结论（标签）")
print("=" * 78)

print("  [T] 51/121 = b1/(2b0²)，纯 SU(N) 与 N 无关。")
print("  [T] 前因子必为 (b0g²)^{-51/121}（负号）：与 exp(-1/(2b0g²)) 由同一个 ∫dg/β 给出，符号锁定。")
print("  [T] g05 的 (g²/16π²)^{+51/121} 是 a^{-1}/Λ 那一支的写法（numerics/REPORT.md:61 的 F(g²)），")
print("      与它自己写的 aΛ 方向（exp 取负号）内部矛盾 —— 是符号错，不是口径差。")
print("  [E] β_W=6, SU(3)：std 0.7479×（≈25% 的两圈截断误差，合理），g05 0.0288×（34×，任何 O(g²) 常数都补不上）。")



