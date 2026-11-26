"""
h01_osschro.py -- D3 OS-SCHRO 判定性数值（规格逐字取 research/v3/species/D06-SPEC.md §3.3, L179-183）

规格原文：
    U : H_phys -> L^2(A/G, dmu_W)，  U H U^{-1} = -Delta_{A/G} + V
    T = (-Delta_{A/G}) 的正平方根
    T^2 - m T >= 0  <=>  spec(T) ∩ (0,m) = ∅      （不需 E_0，因为 T >= 0 由构造）

本脚本只判两件事，不生产任何关于 4D SU(N) 的说法：
 (I)  在 "T = e^{-aH} 转移矩阵" 读法下，T^2 - mT >= 0 到底是不是自动的、以及它给什么；
 (II) 在规格的 "T = (-Delta_{A/G})^{1/2}" 读法下，同一不等式是什么、四模型上是否自动。
纯 numpy，py -3.12。
"""
import numpy as np

np.set_printoptions(linewidth=150, precision=6, suppress=True)


def chk(T2mT_min):
    return "OK(>=0)" if T2mT_min >= -1e-12 else "FAIL(<0)"


print("=" * 78)
print("PART I -- 读法 A: T = e^{-aH}（转移矩阵/收缩算子，0 < lambda <= 1）")
print("  问题: T^2 - mT >= 0 是否自动? m 取 PF 的 lambda_1 还是别的?")
print("=" * 78)

# ---------- 四模型的本征值塔 E_n（只取能写成精确式的那几族） ----------
models = {}

# 1) 自由场（无质量，环上，环长 L，c=1）
L, a = 20.0, 0.25
E_free = np.array([2 * np.pi * n / L for n in range(1, 41)])
models["free(m=0)"] = E_free

# 2) 2D YM（SU(2)，圆上，H = (g^2/2) C_2(j)，精确解）
g, a2 = 1.0, 1.0
j = np.arange(0, 7)
E_2dym = 0.5 * g**2 * j * (j + 1.0)
models["2D YM"] = E_2dym

# 3) M6 = 谐振子(n*omega) (x) 环上自由模( (2pi k/L)^2 /2 )，见 U09-M6.md
om, Lm = 1.0, 10.0
E_m6 = np.concatenate([np.arange(0, 6) * om,
                       0.5 * (2 * np.pi * np.arange(0, 21) / Lm) ** 2])
models["M6"] = np.sort(E_m6)

# 4) 单 plaquette: 无时间方向 => 无 T（纪律：无定义 != 0）

rows = []
for name, E in models.items():
    aa = a if name == "free(m=0)" else (a2 if name == "2D YM" else 1.0)
    lam = np.exp(-aa * np.asarray(E, float))
    lam = np.clip(lam, 0.0, None)
    lam1, lammin = lam.max(), lam.min()
    nA = (lam**2 - lam1 * lam).min()        # T^2 - lambda_1 T
    nB = (lam1 * lam - lam**2).min()        # lambda_1 T - T^2   <- 自动的那一个
    nC = (lam**2 - lammin * lam).min()      # T^2 - lambda_min T
    nD = (lam**2 - 0.0 * lam).min()         # m = 0
    rows.append((name, len(E), lam1, lammin, nA, nB, nC, nD))

print(f"{'model':12s} {'N':>4s} {'lam1':>10s} {'lam_min':>11s} "
      f"{'T^2-lam1*T':>12s} {'lam1*T-T^2':>12s} {'T^2-lammin*T':>13s} {'T^2':>9s}")
for (name, N, l1, lm, nA, nB, nC, nD) in rows:
    print(f"{name:12s} {N:4d} {l1:10.6f} {lm:11.3e} {nA:12.3e} {nB:12.3e} "
          f"{nC:13.3e} {nD:9.3e}   {chk(nA)}/{chk(nB)}")
print()
print("读法 A 的三条判据:")
print("  (A1) T^2 - lambda_1 T >= 0 ?  -> 全部 FAIL（lam_n < lam_1 处为负）; 且 lam1=1 时")
print("       (T^2 - T) = diag(lam_n(lam_n-1)) <= 0，即差一个整体符号。")
print("  (A2) lambda_1 T - T^2 >= 0 ? -> 全部 OK，且它是 PF 的直接推论（lam_n <= lam_1）。")
print("  (A3) T^2 - lambda_min T >= 0 ? -> 全部 OK，但取等在本征值 lam_min 上；")
print("       而 lam_min = exp(-a*E_max)，塔往 UV 延长就 -> 0。")
print()

# A3 的 "空" 量化：lam_min 随塔长度/盒长如何走
print("  A3 的空性量化: m* := lam_min 随截断/体积的走向（'隙' 读法 E_1 >= -ln(m)/a 的右端）")
for name, E in models.items():
    aa = a if name == "free(m=0)" else (a2 if name == "2D YM" else 1.0)
    lam = np.exp(-aa * np.asarray(E, float))
    mstar = lam.min()
    rhs = -np.log(max(mstar, 1e-300)) / aa
    print(f"    {name:12s}: m* = {mstar:.3e},  -ln(m*)/a = {rhs:8.3f}  (= E_max = {np.max(E):8.3f})")
print("    => 'm <= lam_min' 的右端是 E_max（UV 端），不是 E_1（IR 端）。塔越全，m* 越小。")
print()

print("  2x2 显式见证（闭合形式，不靠塔）: T = diag(1, r), 0 < r < 1")
for r in (0.4, 0.9):
    T = np.diag([1.0, r])
    A = T @ T - 1.0 * T
    B = 1.0 * T - T @ T
    print(f"    r={r}: eig(T^2 - 1*T) = {np.linalg.eigvalsh(A)},  eig(1*T - T^2) = {np.linalg.eigvalsh(B)}")
print()

print("  (A4) 'lam_2 <= lam_1^2' 是否自动? 也不是（且在 lam_1<=1 下它与 'lam_2<=lam_1' 不等价:")
print("       lam_1<=1 => lam_1^2 <= lam_1，方向是 'lam_2<=lam_1^2 更强'）")
for (l1, l2) in ((0.8, 0.7), (0.9, 0.99 * 0.9), (0.5, 0.2)):
    print(f"       lam=({l1},{l2}): lam_2<=lam_1 -> {l2 <= l1};  lam_2<=lam_1^2 -> {l2 <= l1*l1}")
print()

print("=" * 78)
print("PART II -- 读法 B（规格原读法）: T = (-Delta_{A/G})^{1/2} >= 0，无界")
print("  等价: T^2 - mT >= 0  <=>  T >= m  <=>  spec(T) ∩ (0,m) = ∅  <=>  -Delta >= m^2 1_(ker)^perp")
print("  => m* := inf{ t : t in spec(T), t > 0 }，且不等式只给 m <= m*（上界）")
print("=" * 78)


def mstar(tvals):
    t = np.asarray(tvals, float)
    t = t[t > 1e-14]
    return t.min() if t.size else 0.0


def form_min(tvals, m):
    t = np.asarray(tvals, float)
    return (t**2 - m * t).min()


# (B1) 2D YM: A/G = SU(2)/Ad = 类函数，-Delta 的本征值 = Casimir C_2(j)=j(j+1)
#      本征函数 = 特征标 chi_j（j = 0,1/2,1,...）。A/G 上 -Delta 的谱是离散的 => 有隙。
jj = np.arange(0, 8)                     # 只取整数 j 的塔，够判 inf>0
t_2dym = np.sqrt(jj * (jj + 1.0))
m2 = mstar(t_2dym)
print(f"(B1) 2D YM: spec(-Delta) = {{j(j+1)}} = {jj*(jj+1.0)}")
print(f"     spec(T) = sqrt(j(j+1)) = {t_2dym}")
print(f"     m* = inf spec(T) = {m2:.6f} = sqrt(2)   [离散谱 => 不自动为零]")
for mm in (0.5 * m2, m2, 1.2 * m2):
    print(f"     m={mm:8.5f}: min(t^2-m t) = {form_min(t_2dym, mm):+.3e}  -> "
          f"{'HOLDS' if form_min(t_2dym, mm) >= -1e-12 else 'FAILS'}")
print()

# (B2) 自由场: A/G 的 -Delta 是 (2pi/L)^2 n^2，spec 触 0 => 无 m>0
print("(B2) 自由场: -Delta 在盒上的谱 (2pi n/L)^2 触 0 => 任何 m>0 都 FAIL（酸测：自由场给 0）")
for LB in (10.0, 20.0, 40.0, 80.0):
    t = np.array([2 * np.pi * n / LB for n in range(1, 200)])
    print(f"     L={LB:5.1f}: m* = {mstar(t):.6f} = 2pi/L = {2*np.pi/LB:.6f},  "
          f"min(t^2-0.5 t) = {form_min(t, 0.5):+.3e}")
print()

# (B3) M6: 振子 x 平坦因子 => 平坦因子那一支把 m* 拉到 0（L -> oo）
print("(B3) M6: osc (x) 环上自由 => spec(T) 含 2pi k/L 一支 => m* = 2pi/L -> 0")
for LB in (10.0, 20.0, 40.0, 80.0):
    t = np.sqrt(np.concatenate([np.arange(1, 6) * 1.0, [0.5 * (2 * np.pi * k / LB) ** 2 for k in range(1, 60)]]))
    print(f"     L={LB:5.1f}: m* = {mstar(t):.6f}, min(t^2-0.1 t) = {form_min(t, 0.1):+.3e}")
print()

# (B4) 单 plaquette: 无空间切片（0 维），T 不存在
print("(B4) 单 plaquette: 无空间切片 => T 无定义（不是 0、不是 >0）。")
print()

# (B5) 等价式的直接数值核对：随机自伴 T，比较三种写法
rng = np.random.default_rng(7)
X = rng.normal(size=(6, 6))
T = (X + X.T) / 2
T = T @ T.T                                # T >= 0
t = np.linalg.eigvalsh(T)
print("(B5) 等价式核对（T>=0 随机 6x6）: spec(T) =", np.sort(t))
for mm in (0.02, 0.2, 1.0, 3.0):
    a1 = np.linalg.eigvalsh(T @ T - mm * T).min()
    a2 = np.linalg.eigvalsh(T - mm * np.eye(6)).min()
    a3 = np.linalg.eigvalsh(T @ T - mm * mm * np.eye(6)).min()   # -Delta >= m^2 形式
    sg = "all>=0" if min(a1, a2, a3) >= -1e-10 else "all<0"
    print(f"     m={mm:4.2f}: min eig(T^2-mT) = {a1:+.3e} | min eig(T-m) = {a2:+.3e} "
          f"| min eig(T^2-m^2) = {a3:+.3e}   => {sg}")
print("     三个 min eig 逐 m 同号 => 'T^2-mT>=0'、'T>=m'、'-Delta>=m^2' 在 T>=0 下等价（符号层）。")
print()

# (B6) 量纲/耦合核对：2D YM 上 m* 与真隙 3g^2/8 是否同一个数
print("(B6) 2D YM 上 m*（-Delta 的谱底）与真隙 3g^2/8 的量纲核对（同一 g，单位 = 同一个）")
print(f"     {'g':>6s} {'m*(T-floor)=g*sqrt2':>22s} {'true gap=3g^2/8':>16s} {'m*/gap':>9s}")
for gg in (0.5, 1.0, 2.0, 8 * np.sqrt(2) / 3, 5.0):
    mstar_g = gg * np.sqrt(2.0)
    gap = 3.0 * gg**2 / 8.0
    print(f"     {gg:6.3f} {mstar_g:22.6f} {gap:16.6f} {mstar_g/gap:9.4f}")
print("     => m* ∝ g, 真隙 ∝ g^2：两者只在 g = 8*sqrt2/3 = 3.771 相交。")
print("        m* > 真隙 时，'m <= m*' 对 certify 真隙毫无约束（它是上界，方向反）。")
print()

print("=" * 78)
print("收束（本脚本支持的三句，不含任何 4D 断言）")
print("  [T] 读法 A：T^2 - lambda_1 T >= 0 为假；自动的是 lambda_1 T - T^2 >= 0；")
print("      T^2 - lambda_min T >= 0 为真但取等在 UV 端（m* -> 0，空）。")
print("  [T] 读法 B：T^2 - mT >= 0  <=>  T >= m  <=> -Delta >= m^2 1_(ker)^perp，只给 m <= inf spec T。")
print("  [C] 2D YM 上 inf spec T = sqrt2（g-无关的谱底）!= 真隙 3g^2/8；自由场/M6 上谱底 -> 0；")
print("      单 plaquette 无定义。=> '自动成立' 在两个读法里都不成立；读法 A 是反号，读法 B 是定理。")
print("=" * 78)
