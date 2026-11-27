# GK2-BIRKHOFF 酸测/精确算术（本轮唯一数值作业，单进程，无 MC，<=300 s）
# 目的：[T1] BH 不等式与"Δ 下界"两种写法严格等价；[T2] d=2 锐性；
#       [T3] 对称正核上 λ2/λ1 <= tanh(Δ/4) 的方向核验（含随机核）；
#       [T4] 单 link toy 的 Δ = 4β 精确值与 β=2.5 处 BH 地板的定量；
#       [T5] 自由场（去掉零模）Δ 的下界随 L 增长；[T6] 地板永不消失；
#       [T7] β=2.5 由实测 m·a=1.707 反解的 Δ 下界；[T8] 3D 片 L=1 的精确下界。
import math, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from decimal import Decimal as D, getcontext
import numpy as np
from scipy.special import iv

getcontext().prec = 50
out = []
def P(*a):
    s = " ".join(str(x) for x in a)
    out.append(s)
    print(s)

def dtanh(x):                      # Decimal tanh
    e = (2 * D(x)).exp()
    return (e - 1) / (e + 1)
def dcoth(x):
    e = (2 * D(x)).exp()
    return (e + 1) / (e - 1)

P("=" * 78)
P("[T1] 两写法的严格等价： g >= -log tanh(D/4)  <=>  D >= 2 log coth(g/2)")
P("=" * 78)
n = 0
worst = D(0)
for i in range(1, 41):
    for j in range(1, 41):
        Dv = D(i) / 2                 # 0.5 .. 20
        g = D(j) / 20                 # 0.05 .. 2.0
        lhs = -(dtanh(Dv / 4)).ln()
        rhs = 2 * (dcoth(g / 2)).ln()
        b1 = (g >= lhs); b2 = (Dv >= rhs)
        assert b1 == b2, (Dv, g)
        gb = -(dtanh(Dv / 4)).ln()          # 边界：g 取地板值
        Db = 2 * (dcoth(gb / 2)).ln()
        worst = max(worst, abs(Db - Dv))
        n += 1
P(f"  网格 {n} 点：两条不等式真假完全一致；边界回代 max|D - 2log coth(gb/2)| = {worst:.3E}")
P(f"  D=10 : -log tanh(2.5) = {float(-(dtanh(D(10)/4)).ln()):.12f}")

P("")
P("=" * 78)
P("[T2] d=2 锐性：K=[[1,b],[b,1]]，Δ=|log(K11 K22/(K12 K21))|，比 = tanh(Δ/4)")
P("=" * 78)
for b in ["0.9", "0.5", "0.1", "0.01"]:
    bd = D(b)
    Dv = (1 / bd ** 2).ln()
    ratio = (1 - bd) / (1 + bd)
    P(f"  b={b:>6}  Δ={float(Dv):>16.10f}  λ2/λ1={float(ratio):>18.14f}"
      f"  tanh(Δ/4)={float(dtanh(Dv/4)):>18.14f}  diff={abs(float(ratio-dtanh(Dv/4))):.2E}")

P("")
P("=" * 78)
P("[T3] 方向核验：随机对称正核，λ2/λ1 <= tanh(Δ/4) ?")
P("=" * 78)
rng = np.random.default_rng(20261126)
worst_slack = 1e9; worst_case = None; tested = 0
for trial in range(3000):
    m = int(rng.integers(2, 7))
    A = rng.random((m, m)) ** rng.uniform(0.3, 3.0)
    K = A + A.T + 1e-12
    ev = np.sort(np.linalg.eigvalsh(K))[::-1]
    if ev[0] <= 0 or ev[1] <= 0:
        continue
    ratio = ev[1] / ev[0]
    Lg = np.log(K)
    Dm = 0.0
    for i in range(m):
        for j in range(m):
            for k in range(m):
                for l in range(m):
                    v = abs(Lg[i, j] + Lg[k, l] - Lg[i, l] - Lg[k, j])
                    if v > Dm:
                        Dm = v
    slack = math.tanh(Dm / 4) - ratio
    tested += 1
    if slack < worst_slack:
        worst_slack = slack; worst_case = (m, ratio, math.tanh(Dm / 4), Dm)
P(f"  测试 {tested} 组：min(tanh(Δ/4) - λ2/λ1) = {worst_slack:.6E}   (>=0 ⇒ 方向确认)")
P(f"  最紧一组: m={worst_case[0]}  λ2/λ1={worst_case[1]:.12f}  tanh(Δ/4)={worst_case[2]:.12f}  Δ={worst_case[3]:.6f}")

P("")
P("=" * 78)
P("[T4] 单 link toy：核 K(U,V)=e^{β cos(UV†)}，Δ_toy = 4β 精确")
P("=" * 78)
P("  U=U'=1, V=V'=-1：四个 cos=(+1,+1,-1,-1) ⇒ 比=e^{4β}；且 |Σ±cos|<=4 ⇒ Δ_toy = 4β")
P("  真隙 λ2/λ1 = I2(β)/I1(β)，g_toy = -log(I2/I1)（D02:110）")
P(f"  {'β':>5} {'Δ=4β':>8} {'BH地板-ln tanh(β)':>20} {'真 g_toy':>16} {'地板/真值':>10}")
for beta in [1.0, 2.0, 2.5, 3.0, 5.0, 10.0]:
    Dv = 4 * beta
    floor = -math.log(math.tanh(Dv / 4))
    g = -math.log(iv(2, beta) / iv(1, beta))
    P(f"  {beta:>5} {Dv:>8.3f} {floor:>20.9f} {g:>16.9f} {floor/g:>10.5f}")

P("")
P("=" * 78)
P("[T5] 自由场（去零模、无质量）：g_L = 2 asinh(sin(π/L))，Δ >= 2 log coth(g/2) → ∞")
P("=" * 78)
P(f"  {'L':>6} {'g_L':>14} {'Δ 下界':>14} {'2log(L/π)':>14} {'BH地板':>14}")
for L in [4, 8, 16, 64, 256, 1024, 4096]:
    g = 2 * math.asinh(math.sin(math.pi / L))
    Dlow = 2 * math.log(1 / math.tanh(g / 2))
    P(f"  {L:>6} {g:>14.8f} {Dlow:>14.6f} {2*math.log(L/math.pi):>14.6f} "
      f"{-math.log(math.tanh(Dlow/4)):>14.8f}")

P("")
P("=" * 78)
P("[T6] 地板永不消失：-log tanh(Δ/4) 单调递减、恒 > 0")
P("      大 Δ 用精确恒等式 -log tanh(x) = -log1p(-2/(e^{2x}+1))，避免浮点下溢成 0")
P("=" * 78)
for Dv in [0.7324, 4.0, 8.0, 10.0, 20.0, 40.0, 100.0, 200.0]:
    x = Dv / 4.0
    val = -math.log1p(-2.0 / (math.exp(min(2 * x, 700.0)) + 1.0))
    P(f"  Δ={Dv:>8.4f}   -log tanh(Δ/4) = {val:.6e}   （tanh 直接算给 {math.tanh(x):.17g}）")

P("")
P("=" * 78)
P("[T7] β=2.5：实测 m·a=1.707（B04:98，RUN3 [C]）反解 Δ 下界")
P("=" * 78)
gmeas = 1.707
Dlow = 2 * math.log(1 / math.tanh(gmeas / 2))
P(f"  Δ >= 2 log coth(1.707/2) = {Dlow:.6f}    （远小于 D02:101 预期 8）")
P(f"  若 Δ=8：BH 地板 = {-math.log(math.tanh(2.0)):.6f}，占真值 1.707 的 {(-math.log(math.tanh(2.0)))/1.707*100:.2f}%")

P("")
P("=" * 78)
P("[T8] 3D 片 L=1（奇 L，-1 链与 1 不同轨）：Δ >= 3βN 精确")
P("=" * 78)
for N in [2, 3]:
    P(f"  SU({N}) β=2.5: Δ(⊔L=1) >= 3βN = {3*2.5*N}")
P("  （L 偶时 -1 链是纯规范 ⇒ 被商掉；故该界只对奇 L=1 有内容）")
P("")
P("=" * 78)
P("[T9] 关键：Δ 不是谱量 —— 块对角核，δ→0 时 Δ→∞ 而 λ2/λ1 固定")
P("=" * 78)
A = np.array([[2.0, 1.0], [1.0, 2.0]])          # 谱 {3,1}
B = np.array([[1.5, 0.5], [0.5, 1.5]])          # 谱 {2,1}
J = np.ones((2, 2))
P("  A=[[2,1],[1,2]] (谱 3,1)，B=[[1.5,.5],[.5,1.5]] (谱 2,1)；K_δ=[[A,δJ],[δJ,B]]")
P("  预期：λ2/λ1 → 2/3（谱不随 δ 变），Δ → ∞（跨块比 A_ij B_kl/δ²）")
P(f"  {'δ':>10} {'λ2/λ1':>14} {'g=-log':>10} {'Δ':>12} {'tanh(Δ/4)':>14} {'slack':>12}")
for d in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8]:
    K = np.block([[A, d * J], [d * J.T, B]])
    ev = np.sort(np.linalg.eigvalsh(K))[::-1]
    ratio = ev[1] / ev[0]
    Lg = np.log(K)
    Dm = 0.0
    n4 = 4
    for i in range(n4):
        for j in range(n4):
            for k in range(n4):
                for l in range(n4):
                    v = abs(Lg[i, j] + Lg[k, l] - Lg[i, l] - Lg[k, j])
                    Dm = max(Dm, v)
    P(f"  {d:>10.0e} {ratio:>14.9f} {-math.log(ratio):>10.6f} {Dm:>12.6f} "
      f"{math.tanh(Dm/4):>14.9f} {math.tanh(Dm/4)-ratio:>12.3e}")
P("  ⇒ 谱（及 g）对 δ 不敏感，Δ 随 δ→0 无界；BH 的 slack = tanh(Δ/4)−λ2/λ1 → 1−2/3")
P("  ⇒ 『Δ 的上界』是独立于谱的额外输入，BH 自身不含它；且可任意失效。")

P("")
P("DONE")
open("research/v3/species/gk2birk_check.log", "w", encoding="utf-8").write("\n".join(out) + "\n")

