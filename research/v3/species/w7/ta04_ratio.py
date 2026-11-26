# ta04_ratio.py — A04-RATIO2: m_lat/(a*Lambda) 在四个 beta 上的比值；是常数吗？
# 只读已知量：两圈 a*Lambda_lat，m1a=4ln(1/u)=18.2638 (X02:131)，m/Lambda_MS=6.28/6.51 (L2)，
# Lambda_MS/Lambda_lat = 28.809(SU3)/19.82(SU2) (A10-SCHEME tab)
import math

PI2 = math.pi ** 2
EXPO = 51.0 / 121.0          # b1/(2 b0^2)，N-无关

def b0(N):
    return 11.0 * N / (48 * PI2)

def aL(beta, N, conv, form='std'):
    """a*Lambda_lat 两圈。conv='A': beta=beta_A=2/(g^2 N)；conv='W': beta=beta_W=2N/g^2。
       form='std'  : (b0 g^2)^{-EXPO} exp(-1/(2 b0 g^2))        <- 标准
       form='g05'  : (g^2/16pi^2)^{+EXPO} exp(-1/(2 b0 g^2))    <- G05 1.2/1.3 的写法"""
    betaW = N * N * beta if conv == 'A' else beta
    g2 = 2.0 * N / betaW
    e = math.exp(-1.0 / (2 * b0(N) * g2))
    if form == 'std':
        return e * (b0(N) * g2) ** (-EXPO)
    return e * (g2 / (16 * PI2)) ** (+EXPO)

def table(betas, N, conv, form):
    C_ms = 6.28 if N == 3 else 6.51
    C_lat = C_ms * (28.809 if N == 3 else 19.82)
    print(f"\n== SU({N}), beta 读作 {conv}, 两圈式={form} ==")
    print(f"{'beta':>9} {'beta_W':>9} {'g^2':>8} {'a*L_lat':>11} {'A: 18.26/(aL)':>14} "
          f"{'B: C(L_lat)':>12} {'B: m/L_MS':>10}")
    for b in betas:
        bW = N * N * b if conv == 'A' else b
        g2 = 2.0 * N / bW
        a_ = aL(b, N, conv, form)
        print(f"{b:>9.6g} {bW:>9.5g} {g2:>8.4g} {a_:>11.5g} {18.2638/a_:>14.4g} "
              f"{C_lat:>12.5g} {C_ms:>10.4g}")

BETAS = [1/48, 0.1, 1.0, 2.5]
for f in ('std', 'g05'):
    for N in (3, 2):
        table(BETAS, N, 'A', f)
print("\n-- beta 读作 beta_W（2.5 的物理点就在这一套）--")
for f in ('std',):
    for N in (3, 2):
        table(BETAS, N, 'W', f)

print("\n== 恒等式（本题核心）==")
print("  m_lat := a*m_phys, a*Lambda := a*Lambda_phys => m_lat/(a*Lambda) = m_phys/Lambda  [T]")
print("  => 比值对 beta、对 a 都不显含；d/dbeta = 0 **由定义给出**，不是动力学结论。")
print("  A 支（分子 = a-无关纯数 r=18.2638）：比值 = r/(a*Lambda)，全部 beta-依赖在分母。")

print("\n== 两圈式的口径误差标定（beta_W=6.00 SU(3)，真值 = 263MeV/28.809/2.91GeV）==")
aLtrue = 0.263 / 28.809 / 2.91
for f in ('std', 'g05'):
    v = aL(6.00, 3, 'W', f)
    print(f"  {f:4s}: a*L_lat = {v:.4e}   两圈/真 = {v/aLtrue:.4f}")
print(f"  真值 : {aLtrue:.4e}   (RUNBOOK:995 的 3.137e-3 同源)")
print(f"  校准后 C(SU3, Lambda_lat) = 6.28*28.809 = {6.28*28.809:.2f}；"
      f"数值直算 0.568/{aLtrue:.4e} = {0.568/aLtrue:.2f}  ← 两条独立路径对到 0.1%")

print("\n== A 支爆炸幅度（SU(3)，beta 读作 beta_A，标准两圈）==")
v = [18.2638 / aL(b, 3, 'A') for b in BETAS]
print("  " + " -> ".join(f"{x:.4g}" for x in v) + f"   (极差 {max(v)/min(v):.3g}x)")
print("  scaling violation O((a*Lambda)^2): beta=1/48 -> %.2f ; beta=2.5 -> %.1e"
      % (aL(1/48, 3, 'A') ** 2, aL(2.5, 3, 'A') ** 2))

print("\n== 格点 m_lat 数据盘点（判死线 c）==")
print("  [T]: beta=2.5, L=1, a=0, SU(2) 的 G''>0 (T10)；L=2 只有 [E] (D/S=0.164)，L=2 的 [T] 被维数判死。")
print("  [E]: 4D m_lat 的实测只有 beta_W=6.00 一点：语料 ma=1.073(88)，真值 0.568 (L2 6.1)。")
print("  => 四个 beta 的 m_lat 格点栅格**不存在**；唯一四点栅格是两圈 aLambda（纯微扰）。")
print("  SU(2) 强耦合 u: beta_A=1/48 u=0.0104 (m1a=18.26)；beta_W=2.5 u=1.014>1 => 闭式失效。")
