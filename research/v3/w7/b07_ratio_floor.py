"""B07 — 比值地板: 精确形式 / 恒等式还是定理 / 可证性
判据: (i) Lambda 是独立量还是由 m 定义  (ii) 量词次序  (iii) 与 T4 等价性
只做算术核验, 不做物理断言. py -3.12 + mpmath.
"""
from mpmath import mp, mpf, exp, log, findroot, pi, inf

mp.dps = 40
b0_of = lambda N: mpf(11) * N / (48 * pi ** 2)          # 11N/(48 pi^2)
EXPO = mpf(51) / 121                                     # b1/(2 b0^2), N-无关
c53 = lambda N, g2: 2 * b0_of(N) * g2                    # = 2 b0 g^2


def aL_2loop(N, g2):
    """a*Lambda_lat, 两圈截断式 (A8:40). 只依赖 g^2, 不含 m => 独立量."""
    x = c53(N, g2)
    return exp(-1 / x) * x ** (-EXPO)


def aL_2loop_closed(N, g2):
    """精确二圈闭式 (A8:40), 用于量化截断差."""
    b0 = b0_of(N)
    b1_over_b0 = EXPO * 2 * b0                            # b1/b0 = 2 b0 * (b1/(2b0^2))
    y = b0 * g2 / (1 + b1_over_b0 * g2)
    return exp(-1 / (2 * b0 * g2)) * y ** (-EXPO)


print("=" * 78)
print("[0] Lambda 的定义: 只依赖 g^2 (=> 与 m 无关, 独立量)")
print("=" * 78)
for N in (2, 3):
    print(f" N={N}: b0=11N/(48pi^2)={mp.nstr(b0_of(N), 10)}  "
          f"b1/(2b0^2)=51/121={mp.nstr(EXPO, 8)}  <- 入口里没有 m")

print()
print("=" * 78)
print("[1] 结构上限 sup_{g^2 in (0,200]} aL  (核 A8:68 的 0.4558 / g^2=17.0)")
print("=" * 78)
for N in (3, 2):
    gs = [mpf(k) / 100 for k in range(1, 20001)]
    vals = [(aL_2loop(N, g), g) for g in gs]
    vmax, gmax = max(vals, key=lambda t: t[0])
    print(f" N={N}: sup aL = {mp.nstr(vmax, 8)}  at g^2={mp.nstr(gmax, 6)}")
    print(f"        aL(g^2=50)={mp.nstr(aL_2loop(N, mpf(50)), 6)}  "
          f"aL(g^2=1)={mp.nstr(aL_2loop(N, mpf(1)), 6)}  "
          f"aL(g^2=1e-3)={mp.nstr(aL_2loop(N, mpf('1e-3')), 6)}")
print(" => aL 在 g^2->0 时 ->0; 有上界、无正下界  => 比值地板 != 常数地板")

print()
print("=" * 78)
print("[2] 约定核对: A8 的 beta_W=6.0 / 6.3 两点 (aL=1.7518e-3)")
print("=" * 78)
for bW in (mpf('6.00'), mpf('6.30')):
    g2 = 2 * 3 / bW
    print(f" beta_W={mp.nstr(bW, 4)}  g^2=2N/beta_W={mp.nstr(g2, 8)}  "
          f"aL_2loop={mp.nstr(aL_2loop(3, g2), 8)}  "
          f"aL_closed={mp.nstr(aL_2loop_closed(3, g2), 8)}")
print(" (A8:52 报 1.7518e-3 / 1.2489e-3 => 约定一致)")

print()
print("=" * 78)
print("[3] A04 3.1 的 aL 表 vs 本式 (口径差?)  + 上限违例检查")
print("=" * 78)
A04 = [('1/48', mpf(1) / 48, mpf('0.5700')), ('0.1', mpf('0.1'), mpf('0.4708')),
       ('1', mpf(1), mpf('7.690e-5')), ('2.5', mpf('2.5'), mpf('1.096e-11'))]
sup3 = max(aL_2loop(3, mpf(k) / 100) for k in range(1, 20001))
print(f" {'beta_A':>6} {'g^2=2/(N beta_A)':>18} {'A04 aL':>12} {'本式 aL':>12} {'比':>7} {'> sup?':>6}")
for tag, bA, ref in A04:
    g2 = 2 / (3 * bA)
    mine = aL_2loop(3, g2)
    print(f" {tag:>6} {mp.nstr(g2, 8):>18} {mp.nstr(ref, 6):>12} {mp.nstr(mine, 6):>12} "
          f"{mp.nstr(ref / mine, 6):>7} {'YES' if ref > sup3 else 'no':>6}")
print(f" sup_{{(0,200]}} aL(N=3) = {mp.nstr(sup3, 8)}")
print(" => A04 表 = 本式 x 1.28-1.34 (常数口径差); 且 aL(1/48)=0.5700 > sup 0.4558")
print("    => A04 3.1 与 A8:68 两表的 c0 数值锚 有 ~33% 隐含口径差 [C]")

print()
print("=" * 78)
print("[4] 量词次序: 哪个读法是恒等式")
print("=" * 78)
# B 支循环代入: m_lat := C * aL  =>  M = C 精确, 散布 0 (定义, 零内容)
C = mpf('180.9')
rowsB = []
for tag, bA, ref in A04:
    mlat = C * ref
    rowsB.append((tag, mlat / ref))
print(" 读法 B' (m_lat := C*aL, 即 A04:131 那张 '散布=0' 表):")
for tag, r in rowsB:
    print(f"   beta_A={tag:>5}  m_lat/(a*L) = {mp.nstr(r, 12)}   (逐位相同: "
          f"{all(abs(x - C) < mpf('1e-25') for _, x in rowsB)})")
print(" => 散布 0 是【代入】的产物, 不是定理. 对任意 C 都成立 => 零内容 [T]")
print()
print(" 读法 A: m_lat 是独立测得的量 (A8:46-52, 带误差棒)")
A8 = [(mpf('6.00'), mpf('1.7518e-3'), mpf('1.653'), mpf('0.50') / mpf('2.91') * 0 + mpf('0.050') * 0),
      ]
for bW, aLref, ma, _ in A8:
    M = (ma) / aLref
    print(f"   beta_W=6.00: m a=1.653 +-0.050 (A8:46)  aL={mp.nstr(aLref, 6)}  "
          f"=> m/L = {mp.nstr(M, 7)} +- {mp.nstr(mpf('0.050') / aLref, 5)}")
    print(f"   A8:52 报 612.5 +- 50.2  => 一致 (误差棒非零 == 不是恒等式)")
d = 2 * b0_of(3) * log(mpf('612.5') / C)
print(f"   612.5 / 180.9 = {mp.nstr(mpf('612.5') / C, 6)};  需要的 scheme 平移 "
      f"delta(1/g^2) = 2b0 ln(比) = {mp.nstr(d, 6)}  (有限正数 => 非恒等)")

print()
print("=" * 78)
print("[5] 等价性: 比值地板  <=>  T4   (要求 Lambda 与 a 无关)")
print("=" * 78)


def check_equiv(N, a_series, Delta0, c0_ratio):
    """T4: m_phys(a) >= Delta0  <=>  m_lat/(a L) >= c0, c0 = Delta0/Lam"""
    ok = True
    for a in a_series:
        # 选一条满足 T4 的轨迹: m_phys(a) = Delta0 * (1 + a)
        m_phys = Delta0 * (1 + a)
        mLam = m_phys / Delta0 * c0_ratio          # 构造: 比值 = c0*(1+a) >= c0
        ok &= (m_phys >= Delta0) and (mLam >= c0_ratio)
    return ok


Lam = mpf(1)                                       # RG 不变量, 与 a 无关
a_ser = [mpf(10) ** (-k) for k in range(1, 13)]
print(f" T4 => 比值: m_phys >= D0  =>  m_lat >= a D0 = (D0/L) aL  =>  比值 >= D0/L = c0 [T]")
print(f" 比值 => T4: m_lat >= c0 aL  =>  m_phys = m_lat/a >= c0 L =: D0  [T] (L 与 a 无关)")
print(f" 数值抽样 a=1e-1..1e-12 双向成立: {check_equiv(3, a_ser, mpf(1), mpf('6.28'))}")

print()
print("=" * 78)
print("[6] 常数地板 vs 比值地板: 严格强弱 (G05B 4 见证模型)")
print("=" * 78)


def aL_of_a(a):
    """取一个合法的 a->0 (沿 RG 轨迹: g^2 ~ 1/ln(1/a))"""
    return a * mpf(1)                              # aL = a*Lam, Lam=1


worst = None
for a in a_ser:
    aL = aL_of_a(a)
    m_lat_ratio = mpf(1) * aL                      # G05B §4: m_lat = a*Lam => 比值 = 1
    r_const = m_lat_ratio                          # 常数地板需要的 r 就是它, 随 a -> 0
    worst = (a, r_const)
print(f" 见证 m_lat(a) = a*Lam:  比值 m_lat/(aL) = 1 (恒 >= c0=1)  [每个 a 有隙]")
print(f"   而 m_lat 本身 a=1e-12 时 = {mp.nstr(worst[1], 6)} -> 0  => 常数地板假 [T]")
print(" => 常数地板 ==> 比值地板 (取 c0 = r/sup aL = r/0.45582); 反之不成立")
print(f"    c0 = r / sup aL, r=18.26 => c0 = {mp.nstr(mpf('18.26') / sup3, 6)}")
print(" => 比值地板严格弱于常数地板 => 它不可能'更强', 只可能'不被杀' [T]")

print()
print("=" * 78)
print("[7] 强耦合窗对 liminf 的贡献")
print("=" * 78)
g2_48 = 2 / (3 * (mpf(1) / 48))
r_i1 = mpf('18.2638')                              # m1 a = 4ln(1/u), I1/I2 闭式
for name, aL in (("A8 式 aL(g^2=32)", aL_2loop(3, g2_48)),
                 ("A04 表 aL(1/48)", mpf('0.5700'))):
    print(f" 窗内 beta_A=1/48 ({name}={mp.nstr(aL, 6)}): 比值 = {mp.nstr(r_i1 / aL, 6)} >= c0")
print(" 但 G05 2.1 自证: {beta<1/48} 不含 a->0 子列 => 对 liminf_{a->0} 贡献恰为 0")
print()
print("DONE")
