# b02_largen.py — B02-LARGEN: `N=inf` + `1/N` 路线的定量核验
# 运行: py -3.12 b02_largen.py > b02_largen.log
# 纪律: 精确算术优先（Fraction）；无拟合；标签 [T]/[M]
from fractions import Fraction as F
import numpy as np
import mpmath as mp

mp.mp.dps = 60

def C2f(N):   return F(N * N - 1, 2 * N)   # 基本表示 Casimir
def C2adj(N): return F(N)                  # 伴随（N>=3 时为第二低；N=2 时 adj=fund 的 3 倍）

hdr = lambda s: print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)

# ---------------------------------------------------------------- PART A
hdr("PART A  2D SU(N) Casimir 塔: 精确 rho^2(N)=C2(adj)/C2(fund)=2N^2/(N^2-1) 与 1/N 展开")
print(f"{'N':>6} {'rho2 (exact)':>22} {'rho2-2':>14} {'N(rho2-2)':>12} {'N^2(rho2-2)':>13} {'rel.err':>9}")
rows = [2, 3, 4, 5, 10, 20, 50, 100, 1000, 10**6]
for N in rows:
    r2 = F(C2adj(N), 1) / C2f(N)
    d = r2 - 2
    print(f"{N:>6} {float(r2):>22.12f} {float(d):>14.3e} {float(N*d):>12.6f} {float(N*N*d):>13.9f} {float(d/2):>8.4%}")
print("\n[T] 恒等式: rho^2 - 2 = 2/(N^2-1)  (精确);  => c1 = lim N(rho2-2) = 0,  c2 = lim N^2(rho2-2) = 2")
print("     => 1/N 展开无 1/N 项, 首修正是 1/N^2 且系数 2; 剩余项精确有界 |rho2-2| = 2/(N^2-1) <= 6/N^2 (N>=2)")
print("     => 一致收敛(几何型): 是. 但注意这不是 4D; 这是唯一有闭式的规范理论.")
print("     !! 关键: rho^2(N) 的收敛半径大, 但物理 N=2,3 上的首项修正 = 33.3% / 12.5% (见 rel.err 列)")

# ---------------------------------------------------------------- PART B
hdr("PART B  't Hooft 定标 (lambda=g^2 N 固定) 下 2D 的质量: m_R = (lambda/2) C2(R)/N")
lam = F(1)
print(f"{'N':>6} {'m_fund':>16} {'m_adj':>12} {'|m_f-m_f(inf)|':>18} {'= lam/(4N^2)':>14}")
for N in [2, 3, 10, 100, 10**6]:
    mf = lam * (N * N - 1) / (4 * N * N)
    ma = lam * C2adj(N) / (2 * N)
    rem = mf - lam / 4
    print(f"{N:>6} {float(mf):>16.12f} {float(ma):>12.9f} {float(abs(rem)):>18.3e} {float(lam/(4*N*N)):>14.3e}")
print("\n[T] 精确: m_fund(N) = (lam/4)(1 - 1/N^2),  m_adj(N) = lam/2 对一切 N (完全 N-无关)")
print("    => 2D 里 m(inf) = lam/4 或 lam/2 > 0, 余项恰为 lam/(4N^2) 几何收敛")
print("    => 逻辑 m(inf)>0 + |rem|<=C/N^2 => N0 = ceil(2C/m(inf)) = 1  => 一切 N 的 m>0  [T](2D 内, 且 2D 早已 [T])")
print("    => 即: 这条路的『最好情形』在 2D 里是定理, 而它给的正是 2D 已知的东西 => 零新地.")

# ---------------------------------------------------------------- PART C
hdr("PART C  锚点 Lambda* = 2 m1^2 的 1/N 缺陷: Lambda* - 2m1^2 = N - (N^2-1)/N = 1/N (精确)")
print(f"{'N':>6} {'m1^2=C2(f)':>18} {'Lam*=C2(adj)':>14} {'Lam*-2m1^2':>14} {'相对 1/N^2':>12}")
for N in [2, 3, 4, 10, 100, 10**6]:
    m1 = C2f(N); Ls = C2adj(N); d = Ls - 2 * m1
    print(f"{N:>6} {float(m1):>18.9f} {float(Ls):>14.9f} {float(d):>14.3e} {float(d/Ls):>11.4%}")
print("\n[T] v05 的 anchor `Lambda* = 2 m1^2` 只在 N=inf 成立, 缺陷恰为 1/N (相对 1/N^2).")
print("    => 把 N=inf 的 anchor 搬到物理 N: N=2 相对错 1/N^2 = 25%, N=3 错 11.1%  [T]")

# ---------------------------------------------------------------- PART D
hdr("PART D  自由场酸测: 有限环面 T^4 上的自由无质量场 rho^2 = ?")
vals = sorted({a*a + b*b + c*c + d*d for a in range(-4, 5) for b in range(-4, 5) for c in range(-4, 5) for d in range(-4, 5)})
print("  环面本征值 (2pi/L)^2 * n, 不同的 n:", vals[:8])
e1, e2 = vals[1], vals[2]
print(f"  e1={e1} e2={e2} => rho^2 = e2/e1 = {e2/e1:.15f}   (n=0 是真空, 必须扣掉)")
print("  !! 与 PART A 的 N=inf 值 2 完全同值 => N=inf 的 rho^2 与自由场不可区分 (RB-9 的『地板盲』, 本轮独立复算精确 2)")
print("\n  加质量 mu (单位 (2pi/L)^2=1):  rho^2 = (mu^2+2)/(mu^2+1)")
for mu in [0, 0.5, 1, 2, 5, 10, 100]:
    r2 = (mu * mu + 2) / (mu * mu + 1)
    print(f"    mu={mu:>6}: rho^2 = {r2:.12f}   (I2 判据 rho^2>1: {r2>1})")
print("  => rho^2>1 对任意手放大质量都真, 且 mu->inf 时 rho^2->1; 判据只反映『两个数的比』, 完全不看是否有原子")

# ---------------------------------------------------------------- PART E
hdr("PART E  1/N 级数的两种尾: 几何 (2D 实况) vs 阶乘 (1-matrix/Penner 属型 [M])")
print("  E1 几何尾: a_k = 2/(N^{2k}) 型 => 部分和单调收敛, 误差 <= 首项被略项")
for N in [2, 3, 10]:
    s = sum(F(2, N**(2 * k)) for k in range(1, 200))
    print(f"    N={N:>3}: sum_(k>=1) 2/N^2k = {float(s):.15f}  (精确 2/(N^2-1) = {float(F(2, N*N-1)):.15f})")
print("\n  E2 阶乘尾 [M] (Penner/Gaussian 型: F_g = B_2g/(2g(2g-2))): 部分和不收敛, 有最优截断阶 g*")
for N in [2, 3, 5, 10]:
    terms = []
    for g in range(2, 400):   # g=1 时分母 2g(2g-2)=0, 从 g=2 起
        B = abs(mp.bernoulli(2 * g))
        terms.append((mp.mpf(B) / (2 * g * (2 * g - 2)) / mp.mpf(N) ** (2 * g), g))
    best = min(terms)
    # 首次开始上升的阶
    gstar = best[1]
    print(f"    N={N:>3}: 最优截断 g*={gstar:>3}, 最小项 ~ {mp.nstr(best[0], 5)}, e^(-2pi N) = {mp.nstr(mp.e**(-2*mp.pi*N), 5)}")
print("  => 阶乘尾 => 固定阶 k 的界存在(对 k 一致? 否) => 『一致收敛』在 4D 不成立 [C];")
print("     但最优截断误差 ~ e^(-cN) 仍 -> 0 => 收敛性不是唯一障碍, 真正的障碍是 Borel 可和性 [O]")

# ---------------------------------------------------------------- PART F
hdr("PART F  四模型试金石 (N=inf 型量的直接取值)")
print(f"{'model':<26}{'m':>10}{'G\" / f1^2':>16}{'N=inf 型量 (Lam*, rho^2)':>34}{'应给':>8}")
print(f"{'free field (massless)':<26}{'0':>10}{'0 / 0':>16}{'rho^2=2 ; Lam*=2e1 (torus)':>34}{'0':>8}")
print(f"{'2D YM (N=inf)':<26}{'>0':>10}{'0 / =0':>16}{'rho^2=2 ; Lam*=2m1^2+1/N':>34}{'>0':>8}")
print(f"{'M6 (osc (x) free)':<26}{'0':>10}{'2w>0 / -':>16}{'Lam*-m1^2=0 (0=0) ; rho^2=0/0':>34}{'0':>8}")
print(f"{'single plaquette':<26}{'undef':>10}{'+0.136 / -':>16}{'无 transfer => 无定义':>34}{'undef':>8}")
print("\n  净判: 自由场一行 与 2D YM N=inf 一行 的 rho^2 都是 2 => 该量不能区分『有原子』与『无原子』[T]+[M]")
print("        M6 一行: 若读 rho^2 => 0/0 无定义(不给 0); 若读 Lam*-m1^2 => 0=0 (真但空)")
print("        单 plaquette: 无定义 => 唯一干净的『合规』")
print("\n  => 通过试金石的唯一量是 m(inf) 自身(= 目标在 N=inf 的改写); 代理量 rho^2 / Lam*=2m1^2 在自由场上硬失败.")
print("\n[END]")
