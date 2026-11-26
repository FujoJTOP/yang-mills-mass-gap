# w7/m05_ratio.py -- M05-RATIO: 新靶 liminf_{a->0} m_lat(a)/(a*Lambda) >= c0 的可判化
# 纯 stdlib（ponytail: 一行能算就不装库）。输出 = w7/m05_ratio.log
# 四件事: (1) 两圈 a*Lambda 复算 g05 的数; (2) 常数地板 vs 比值地板 两分支;
#          (3) 四模型比值; (4) liminf != 极限 的自造序列。
import math

PI2 = math.pi ** 2


def aL_two_loop(beta_A, N):
    """a*Lambda_lat = exp(-1/(2 b0 g^2)) * (g^2/16pi^2)^(51/121), beta_A = 2/(g^2 N)."""
    beta_W = N * N * beta_A
    g2 = 2.0 * N / beta_W
    b0 = 11.0 * N / (48.0 * PI2)
    expo = 1.0 / (2.0 * b0 * g2)
    return math.exp(-expo) * (g2 / (16.0 * PI2)) ** (51.0 / 121.0), beta_W, g2, expo


print("=" * 78)
print("[1] 两圈 a*Lambda_lat（Wilson beta=2N/g^2, b1/2b0^2=51/121）-- 复核 G05 1.3")
print("=" * 78)
print(f"{'beta_A':>9} {'N':>2} {'beta_W':>8} {'g^2':>7} {'1/(2b0g2)':>10} {'aL':>8}")
for beta_A in (1.0 / 48.0, 1.0 / 72.0):
    for N in (2, 3):
        aL, bW, g2, ex = aL_two_loop(beta_A, N)
        print(f"{beta_A:9.5f} {N:2d} {bW:8.4f} {g2:7.3f} {ex:10.5f} {aL:8.4f}")
print("  G05 靶值: beta_A=1/48 -> 0.4837(SU2) / 0.4077(SU3)")

# 反解: aL<=0.01 需要多大 beta_W
print("\n  反读: 求 beta_W 使 aL=0.01 / 0.003 (两圈, 消去 N 后按 beta_W 解)")
for target in (0.01, 0.003):
    lo, hi = 0.05, 200.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        g2 = 6.0 / mid                      # SU(3): g^2=6/beta_W, beta_A=beta_W/9
        aL, *_ = aL_two_loop(mid / 9.0, 3)
        if aL > target:
            lo = mid
        else:
            hi = mid
    print(f"    aL={target:6.3f}  ->  beta_W ~ {0.5*(lo+hi):7.3f}  (SU3)")

print()
print("=" * 78)
print("[2] 两分支: 常数地板(自我否定) vs 比值地板")
print("=" * 78)
print("  A 支: m_lat = r = 18.26 (I2 闭式, 纯数)  => ratio = r/(aL), m_phys/Lambda = 同值")
print("  B 支: m_lat = c0*(aL)                   => ratio = c0 (常数), m_phys = c0*Lambda")
r_const = 4.0 * math.log(1.0 / 0.0104)        # 4 ln(1/u), u=0.0104
print(f"  r = 4ln(1/0.0104) = {r_const:.4f}")
c0 = 6.3                                      # L2-GLUEBALL-LAMBDA: m_0++/Lambda_MS = 6.3(3) [C]
x_phys = 1.0                                  # 固定物理间距 x = 1/Lambda
print(f"\n{'aL':>10} {'A: ratio':>12} {'A: corr e^-m x':>16} {'B: ratio':>9} {'B: corr':>10}")
for aL in (0.4837, 0.10, 0.01, 1e-3, 1e-4):
    rA = r_const / aL
    rB = c0
    corrA = math.exp(-rA * x_phys)             # <O(x)O(0)> ~ e^{-m_phys x}, x=1/Lambda
    corrB = math.exp(-rB * x_phys)
    print(f"{aL:10.4f} {rA:12.3f} {corrA:16.3e} {rB:9.3f} {corrB:10.4f}")
print("  读法: A 支 ratio->inf 且 corr->0 (连续极限平凡) = 自我否定; B 支两者都有限.")

print()
print("=" * 78)
print("[3] 四模型: 新靶 = liminf m_lat/(a*尺度) >= c0 ?")
print("=" * 78)
print(f"{'模型':<22} {'m(物理)':<28} {'比值':<22} {'新靶'}")
print(f"{'自由场':<20} {'0 (零模)':<28} {'0':<22} {'假 ✓(酸测)'}")
print(f"{'2D YM (SU2)':<18} {'(g^2/2)C2 = 3g^2/8':<26} {'3/8 = 0.375 精确, a-无关':<20} {'真 (T)'}")
L_list = (2.0, 10.0, 40.0, 1000.0)
om = 1.0
m6 = [min(om, 2.0 * PI2 / (L * L)) for L in L_list]
print(f"{'M6 (omega=1)':<18} {'m_L=min(w,2pi^2/L^2)':<27} "
      f"{'/'.join(f'{m/om:.4f}' for m in m6)}", " L->inf 得 0")
print(f"{'单 plaquette':<17} {'无 transfer':<28} {'无定义':<22} {'不可判 (无定义!=0)'}")
print("  签名 (0,>0,0,undef) 与 m 的签名逐位相同 => 四模型 2-bit, 无鉴别力")

print()
print("=" * 78)
print("[4] liminf != 极限: 目标不要求极限存在")
print("=" * 78)
tail = [6.3 + 3.0 * ((-1) ** n) for n in range(12)]     # 振荡: 9.3 / 3.3 交替
print("  序列 C_n = 6.3 + 3(-1)^n:", " ".join(f"{v:.1f}" for v in tail))
print(f"  liminf = {min(tail):.2f}, limsup = {max(tail):.2f}, 极限不存在")
print("  => liminf >= 3.3 成立而极限不存在: 陈述层不需要 L3/W3 (存在性), 只证明层需要")

print()
print("=" * 78)
print("[5] G05 1.3 的 22x 不自洽 复核")
print("=" * 78)
aL = 0.4837
for m_over_L in (1.7, 6.3):
    m_lat = m_over_L * aL
    print(f"  m/Lambda={m_over_L:5.2f} -> m_lat = {m_lat:.3f};  I2 闭式 18.26;  比 = {18.26/m_lat:.1f}x")
print("  注: c0 依赖 scheme 常数(Λ_lat/Λ_MS), 故 c0=O(10)-O(100) 而非 O(1)")
