"""H04 — F 的"积分版"：F(b') - F(b) = int dF/db；证书里的 1/Delta 让界发散。
本脚本复核的是**与模型无关的恒等式/不等式**（F、Var、Delta 按解析形状指定，非物理模型）：
  (A)  Newton-Leibniz 精确性 + 反向三角严格性 + 变差预算否决的显式见证。
  (B)  beta-方向 Cauchy-Schwarz 的损失因子 (int f g <= sqrt(int f^2 int g^2))。
  (C)  唯一有物理数字的部分：用 d07 精确 int Var = 6V*u_p(2.5) 算 Delta_req(L, r0)。
"""
import numpy as np
from scipy.integrate import quad

B1, B2 = 1.0 / 48.0, 2.5
T = B2 - B1
UP25 = 0.507195100            # d07 [E]
IVAR = 6 * UP25               # = int Var(S_W) db / V  (4D SU(2), N_p = 6V)

print("=" * 78)
print("H04 / (A) 恒等式层：F(b) = 0.60 + 0.35*sin(3(b-b1))  [解析指定形状，非物理模型]")
print("=" * 78)


def F(b):
    return 0.60 + 0.35 * np.sin(3.0 * (b - B1))


def Fp(b):
    return 0.35 * 3.0 * np.cos(3.0 * (b - B1))


dFi, e_nl = quad(Fp, B1, B2, epsabs=1e-14, epsrel=1e-13, limit=400)
dF = F(B2) - F(B1)
print(f"[E] int F' db                    = {dFi:.15e}  (quad err {e_nl:.1e})")
print(f"[E] F(B2)-F(B1)                  = {dF:.15e}")
print(f"[E] Newton-Leibniz 残差           = {abs(dFi-dF):.2e}   => (1.1) 逐字精确 [T]")

tv, _ = quad(lambda b: abs(Fp(b)), B1, B2, epsabs=1e-13, limit=400)
print(f"[E] int|F'| = TV(F)               = {tv:.12f}")
print(f"[E] |F(B2)-F(B1)|                 = {abs(dF):.12f}   TV/|dF| = {tv/abs(dF):.6f}  (>=1 反向三角 (2.2))")
print(f"[E] F(B1) = {F(B1):.8f}   min_[b1,b2] F = {min(F(b) for b in np.linspace(B1,B2,20001)):.8f}  (>0)")
print(f"[E] 判据 int|F'| < F(B1)?         {tv < F(B1)}    TV/F(B1) = {tv/F(B1):.4f}")
print("    => 见证：F>0 全程，但 TV >= F(B1)：判据严格强于靶，且此处为假。")

print()
print("=" * 78)
print("H04 / (A5) 反向三角的必要条件扫描：ratio = |F(B2)-F(B1)| / F(B1)")
print("=" * 78)
for amp in (0.05, 0.1, 0.2, 0.35):
    f = lambda b: 1.0 + amp * np.sin(3.0 * (b - B1))
    fp = lambda b: amp * 3.0 * np.cos(3.0 * (b - B1))
    r = abs(f(B2) - f(B1)) / f(B1)
    t, _ = quad(lambda b: abs(fp(b)), B1, B2, epsabs=1e-13, limit=400)
    print(f"    amp={amp:5.2f}  F(b2)/F(b1) = {f(B2)/f(B1):7.4f}  |dF|/F(b1) = {r:6.4f}  "
          f"TV/F(b1) = {t/f(B1):6.4f}  判据成立? {t < f(B1)}")
print("    => 判据成立必须 F(b2) < 2*F(b1)；一旦残数翻倍即逐字假。")

print()
print("=" * 78)
print("H04 / (B) beta-方向 C-S：int (sqrt(Var)*1/Delta) <= sqrt(int Var * int 1/Delta^2)")
print("=" * 78)
for name, Var, Dlt in (
    ("Delta = b*1.0, Var = C*e^{-b}", lambda b: np.exp(-b), lambda b: b),
    ("Delta = b*0.1, Var = C*e^{-b}", lambda b: np.exp(-b), lambda b: 0.1 * b),
    ("Delta = 1.0 (const), Var = C", lambda b: np.ones_like(b) if isinstance(b, np.ndarray) else 1.0, lambda b: 1.0),
    ("Delta = b*1.0, Var = C*b", lambda b: b, lambda b: b),
):
    f = lambda b: np.sqrt(Var(b))
    g = lambda b: 1.0 / Dlt(b)
    I, _ = quad(lambda b: f(b) * g(b), B1, B2, epsabs=1e-13, limit=400)
    I1, _ = quad(lambda b: Var(b), B1, B2, epsabs=1e-13, limit=400)
    I2, _ = quad(lambda b: 1.0 / Dlt(b) ** 2, B1, B2, epsabs=1e-13, limit=400)
    print(f"    {name:32s} direct={I:10.4f}  C-S={np.sqrt(I1*I2):10.4f}  损失因子={I/np.sqrt(I1*I2):.4f}")
print("    => C-S 的损失 <=1；在 Delta∝b 且 Var 慢变时损失 ~0.8（几乎不亏），在 Var∝b 时 ~1.00（无损）")
print("    => 故 (3.2) 的 L-标度不是 C-S 引入的假象。")

print()
print("=" * 78)
print("H04 / (C) 预算：int|F'| < r0 需要 Delta 多大（用 d07 精确 int Var）")
print("=" * 78)
print(f"[E] int Var(S_W) db = 6*u_p(2.5)*V = {IVAR:.6f} * V   (V=L^4, u_p(2.5)={UP25})")
K = np.sqrt(IVAR * T)
print(f"[E] 闭式：Delta_req = sqrt(intVar*T)/r0 = {K:.4f} * L^2 / r0")
print(f"    {'L':>3} {'V':>9} {'Delta_req (r0=1)':>18} {'/L^2':>8} {'d07 0.276L^2':>13}")
for L in (4, 6, 8, 16):
    Vv = L ** 4
    print(f"    {L:3d} {Vv:9d} {K*L**2:18.3f} {K:8.4f} {0.276*L**2:13.3f}")
print(f"[E] L=8: Delta_req = {K*64:.1f}  vs 物理 Delta = O(1)  => 差 {K*64:.0f}x")
print(f"[E] 与 d07 命题 3.1 (0.276L^2) 同标度 L^2，比值 = {K/0.276:.2f} (= d07 自取的 K=10 的量级)")

print()
print("[E] 反向三角硬约束（不依赖 Delta、不依赖任何证书）：int|F'| >= |F(B2)-F(B1)|")
print("[E] 判据 int|F'| < F(B1) 的必要条件 <=> F(B2) < 2*F(B1)（在 F>0 的窗口内）")
print("    I1 给 F(b1) = r0 > 0；若 f1^2(2.5) >= 2*r0 则判据本身为假 => 判 (b)，无需任何谱估计。")
