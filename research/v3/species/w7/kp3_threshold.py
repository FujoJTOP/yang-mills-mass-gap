# kp3_threshold.py — K̂=8 第三问：门槛 7.0633 的精确值、它约束的是谁、(1.7)^6 的关系
# 只算数，不下结论。dps=50。整数精确枚举 phi^U 双路径。
from mpmath import mp, mpf, e as mp_e, log, power
mp.dps = 50

N = 2
Delta = 20
eDelta = Delta * mp_e                      # eΔ = e·Δ = 54.3656365691809047...
beta = mpf(1) / 48
u = beta / (2 * N * N)                     # u = β/(2N²) = β/8
cf = beta / (2 * N)                        # c_f = β/(2N) = 2u
T = 1 / (eDelta * u)                       # 门槛：判据 X·eΔ·u < 1 ⟺ X < T
T2 = 2 * N * N * 48 / eDelta               # 题目给的另一写法

print("=" * 78)
print("§0 门槛的精确值")
print("=" * 78)
for k, v in [
    ("eΔ = 20e", eDelta),
    ("u(1/48) = β/8", u),
    ("c_f(1/48) = β/4", cf),
    ("1/(eΔ·u)      ", T),
    ("2N²·48/eΔ     ", T2),
    ("比值           ", T / T2),
    ("8·48/eΔ        ", 8 * 48 / eDelta),
    ("T/4 (逐片率门槛)", T / 4),
    ("T/2            ", T / 2),
    ("eΔ/48          ", eDelta / 48),
    ("8/T = eΔ/48 ?  ", 8 / T - eDelta / 48),
]:
    print(f"  {k:22s} = {mp.nstr(v, 25)}")

print()
print("=" * 78)
print("§1 (1.7)^6 与 T 的关系")
print("=" * 78)
q = mpf("1.7")
q6 = q ** 6
rate_star = T / 4                          # 精确临界率 = 1.7658213...
star6 = rate_star ** 6                     # 精确临界 |φ^U|
print(f"  (1.7)^6                = {mp.nstr(q6, 25)}")
print(f"  精确临界率 T/4         = {mp.nstr(rate_star, 25)}")
print(f"  (T/4)^6                = {mp.nstr(star6, 25)}")
print(f"  (T/4)^6 / (1.7)^6      = {mp.nstr(star6 / q6, 25)}")
print(f"  ((T/4)/1.7)^6          = {mp.nstr((rate_star / q) ** 6, 25)}")
print(f"  4·1.7                  = {mp.nstr(4 * q, 25)}   (对比 T = {mp.nstr(T, 12)})")
print(f"  T - 4·1.7              = {mp.nstr(T - 4 * q, 25)}")
print(f"  (1.7)^6 / T            = {mp.nstr(q6 / T, 25)}")
print(f"  (1.7)^6 / 2            = {mp.nstr(q6 / 2, 25)}")
print(f"  (1.7)^6 / (1/48)       = {mp.nstr(q6 * 48, 25)}")
print(f"  T^6                    = {mp.nstr(T ** 6, 25)}")
print(f"  (1.7)^6 是 T 的哪一次幂？  log(q6)/log(T) = {mp.nstr(log(q6) / log(T), 20)}")
print(f"  (1.7)^6^(1/6)·4        = {mp.nstr(4 * power(q6, mpf(1) / 6), 25)}")

print()
print("=" * 78)
print("§2 M07 表内两个 1.7658 尾数的自查")
print("=" * 78)
print(f"  T/4                = {mp.nstr(T / 4, 25)}")
print(f"  8/(4eΔ)            = {mp.nstr(8 / (4 * eDelta), 25)}")
print(f"  M07 表写 1.76582025642，差 = {mp.nstr(T / 4 - mpf('1.76582025642'), 10)}")
print(f"  1.766/48           = {mp.nstr(mpf('1.766') / 48, 25)}")
print(f"  (T/4)/48           = {mp.nstr(T / 4 / 48, 25)}")

print()
print("=" * 78)
print("§3 phi^U 精确枚举（整数，双路径）")
print("=" * 78)


def components(n, edge_mask, edges):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    cnt = n
    for i, (a, b) in enumerate(edges):
        if edge_mask >> i & 1:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
                cnt -= 1
    return cnt


def phi_two_paths(name, n, edges):
    m = len(edges)
    total_bfs = 0
    comp_poly = [0] * (n + 1)          # [y^c] 累加 (-1)^|A|
    for mask in range(1 << m):
        c = components(n, mask, edges)
        comp_poly[c] += -1 if bin(mask).count("1") % 2 else 1
        if c == 1:
            total_bfs += -1 if bin(mask).count("1") % 2 else 1
    return name, n, m, total_bfs, comp_poly[1]


# 立方体 6 面的黏附图 = 八面体 K_{2,2,2}：6 点 12 边
oct_edges = [(i, j) for i in range(6) for j in range(i + 1, 6)
             if not ((i, j) in [(0, 1), (2, 3), (4, 5)])]
K6_edges = [(i, j) for i in range(6) for j in range(i + 1, 6)]

rows = []
for nm, n, ed in [("K_{2,2,2} (立方体簇/八面体)", 6, oct_edges),
                  ("K_6 (一个 link 上 6 片)", 6, K6_edges)]:
    rows.append(phi_two_paths(nm, n, ed))

for nm, n, m, p1, p2 in rows:
    ok = "一致" if p1 == p2 else "不一致!"
    print(f"  {nm:30s} n={n} |E|={m:2d}  2^|E|={2**m:6d}  phi^U 路径①={p1:6d}  [y^1]={p2:6d}  {ok}")

phi_cube = rows[0][3]
print()
print("=" * 78)
print("§4 三条 X 与同一个门槛")
print("=" * 78)
Chat = 4                                   # sup_gamma (∏R)^{1/|γ|}（立方体取到）
Cbar = cf / u                              # = N，权重槽
tau_cube = 384
X_route1 = Chat * Cbar                     # C06 (1.1)
X_route2 = (4096 * tau_cube) ** (1 / 6)    # C06 (5.1) K̂_τ，树界版
X_route3 = (4096 * abs(phi_cube)) ** (1 / 6)   # ur3 的 Λ_Ur，真 phi^U 版
for nm, X in [("X1 = Ĉ·C̄", X_route1), ("X2 = K̂_τ (sup ∏R·τ)", X_route2),
              ("X3 = Λ_Ur (sup ∏R·|φ^U|)", X_route3)]:
    Xm = mpf(X)
    print(f"  {nm:26s} X = {mp.nstr(Xm, 10):>10s}   X/T = {mp.nstr(Xm / T, 8)}   "
          f"β_max = 8/(X·eΔ) = {mp.nstr(8 / (Xm * eDelta), 10)}  = {mp.nstr(48 * 8 / (Xm * eDelta), 8)}×(1/48)")
print()
print(f"  |φ^U(立方体)|^{{1/6}} = {mp.nstr(power(mpf(abs(phi_cube)), mpf(1) / 6), 20)}   vs  C̄ = N = {mp.nstr(Cbar, 10)}")
print(f"  4·|φ^U|^{{1/6}}      = {mp.nstr(4 * power(mpf(abs(phi_cube)), mpf(1) / 6), 20)}")
print(f"  τ 版: 4·τ^{{1/6}}     = {mp.nstr(4 * mpf(tau_cube) ** (1 / 6), 20)}")

print()
print("=" * 78)
print("§5 判据在 β=1/48 处是否成立（三读）")
print("=" * 78)
for nm, X in [("X1=8", X_route1), ("X2=10.784", X_route2), ("X3=8", X_route3)]:
    val = mpf(X) * eDelta * u
    print(f"  {nm:12s} X·eΔ·u = {mp.nstr(val, 12)}  {'<1 ✔' if val < 1 else '>1 ✘'}")
print(f"  eΔ·β = {mp.nstr(eDelta * beta, 12)}   (= X1·eΔ·u，因 8u = β)")
