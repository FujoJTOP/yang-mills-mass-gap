"""KP1-b —— SU(2) 逐置换类 / 键多重度：R(kappa)、Wg(sigma) 类求和、以及与整数 phi^U 的归一化核对。

三件事：
 1) mult(1 in f^k) 用 SU(2) CG 精确迭代（整数）⇒ Catalan；
 2) 每 link 的因子化权重 R(kappa) = mult(1 in f^{kappa+2}) / mult(1 in f^{kappa})（精确 Fraction）
    ⇒ R(2)=2, R(4)=5/2, R(6)=14/5（与 C06 §2.3 逐字对照）；
    R(2) 再用"非交叉配对计数"独立算一次（Temperley–Lieb / Catalan）；
 3) U(2) 的 Wg(sigma) 类表（k<=2）并核恒等式 sum_sigma Wg(sigma) N^{c(sigma)} = 1：
    归一化里的 1/N 在这里**恰好消掉**成整数 1；N<k 处该表退化（如实记），故 k>=3 用 TL/Catalan。
"""
from fractions import Fraction
from itertools import permutations
import sympy as sp

N = 2  # SU(2)

# ---------- 1) mult(1 in f^k) ----------
def mult_triv(kmax):
    reps = {Fraction(0): 1}   # 真空
    out = [reps.get(Fraction(0), 0)]
    for _ in range(kmax):
        new = {}
        for j, m in reps.items():
            for jj in (j + Fraction(1, 2), j - Fraction(1, 2)):
                if jj >= 0:
                    new[jj] = new.get(jj, 0) + m
        reps = new
        out.append(reps.get(Fraction(0), 0))
    return out


m = mult_triv(8)
print("=" * 74)
print("1) mult(1 in f^k)（SU(2) CG 精确，k=0..8）")
print("=" * 74)
print("k   :", list(range(9)))
print("mult:", m)
print("Catalan(k/2) 应为 1,1,2,5,14（偶数 k）—— ", m[::2])
assert m[::2] == [1, 1, 2, 5, 14]

# ---------- 2) 每 link 因子化权重 R(kappa) ----------
print()
print("=" * 74)
print("2) R(kappa) = mult(1 in f^{kappa+2}) / mult(1 in f^{kappa})   [C06 §2.3 口径]")
print("=" * 74)
R = {}
for kap in (2, 4, 6):
    R[kap] = Fraction(m[kap + 2], m[kap])
    print(f"  R({kap}) = {m[kap+2]}/{m[kap]} = {R[kap]}   (= {float(R[kap]):.6f})")
assert R[2] == 2 and R[4] == Fraction(5, 2) and float(R[6]) == 2.8
print("  语料 C06: R(2)=2, R(4)=5/2, R(6)=2.8  ⇒ 逐位复现 [E]")


# ---------- 2b) 非交叉配对（TL）独立算 R(2) ----------
def noncross_pairings(slots):
    """非交叉完美配对计数；slots 为方向串（+1=f, -1=fbar），要求配对面方向相反。"""
    slots = tuple(slots)

    def rec(pos, alive):
        if not alive:
            return 1
        i = min(alive)
        tot = 0
        for j in alive:
            if j <= i:
                continue
            if slots[i] * slots[j] != -1:
                continue
            # 非交叉：i 与 j 之间的存活点数为偶数，且 i,j 之间内部配对
            between = [x for x in alive if i < x < j]
            if len(between) % 2:
                continue
            inside = rec(i, set(between))
            outside = rec(i, set(x for x in alive if x > j))
            tot += inside * outside
        return tot

    return rec(0, set(range(len(slots))))


print()
print("   非交叉配对（Temperley–Lieb）:")
for kk in (2, 4, 6):
    slots = [1 if i % 2 == 0 else -1 for i in range(kk)]
    c = noncross_pairings(slots)
    print(f"     (f fbar)^{kk//2}: 非交叉配对数 = {c}   mult(1 in f^{kk}) = {m[kk]}")
    assert c == m[kk]
print("   ⇒ R(2) = 2 = 两条非交叉配对各自权 1，无 1/N 残留 [T]/[E]")

# ---------- 3) U(2) Wg(sigma) 类表 ----------
print()
print("=" * 74)
print("3) U(2) Weingarten 类表（k<=2）与恒等式 sum_sigma Wg(sigma) N^{c(sigma)} = dim Inv(f^k)")
print("=" * 74)


def cycles(sig):
    seen, c = set(), 0
    for i in range(len(sig)):
        if i in seen:
            continue
        c += 1
        j = i
        while j not in seen:
            seen.add(j)
            j = sig[j]
    return c


# U(N) Weingarten 类值（k=1: Wg=1/N；k=2: id -> 1/(N^2-1), transposition -> -1/(N(N^2-1))），
# 这两条由矩恒等式本身定死（下两行逐一核）：
WG = {
    1: {1: Fraction(1, N)},
    2: {2: Fraction(1, N * N - 1), 1: Fraction(-1, N * (N * N - 1))},
}
for k in (1, 2):
    print(f"  k={k}: 按循环数 c 的类值 -> " +
          ", ".join(f"c={c}: Wg={WG[k][c]}" for c in sorted(WG[k], reverse=True)))
    # 矩恒等式（这是 Wg 的定义性恒等式）: sum_sigma Wg(sigma) N^{c(sigma^{-1} tau)} = delta_{tau,id}
    for tau in permutations(range(k)):
        tot2 = 0
        for sig in permutations(range(k)):
            comp = tuple(sig[tau[i]] for i in range(k))
            tot2 += WG[k][cycles(comp)] * N ** cycles(sig)
        assert tot2 == (1 if tau == tuple(range(k)) else 0), (k, tau, tot2)
    print("       矩恒等式 sum_sigma Wg(sigma)N^{c(sigma^-1 tau)} = delta_{tau,id}：逐 tau 全过")

# 只在与矩恒等式匹配的对象上读"重数"：k=2 时 dim Hom(1, f (x) fbar) = sum_sigma Wg N^{c(sigma)}
tot = sum(WG[2][cycles(sig)] * N ** cycles(sig) for sig in permutations(range(2)))
print(f"  dim Hom(1, f (x) fbar) = sum_(sigma in S_2) Wg(sigma) N^c = {tot}  (= multiplicities 1) ")
assert tot == 1

print("  注：k>=3 且 N=2 时 sum 公式出现 s_lambda(1^2)=0 的项 ⇒ 该表退化，")
print("      故 per-link 多重度一律用 TL/Catalan（上面 §2），不硬套 Wg。")

# ---------- 4) 立方体簇装配 ----------
print()
print("=" * 74)
print("4) 立方体簇装配：全 12 条 link 皆 kappa=2")
print("=" * 74)
prodR = Fraction(1)
for _ in range(12):
    prodR *= R[2]
phi = -64
print(f"  prodR = R(2)^12 = {prodR} = 2^12")
print(f"  |phi^U| (整数图系数，本文件不重算，见 kp1_phi.log) = {abs(phi)}")
print(f"  prodR*|phi^U| = {prodR*abs(phi)} = 2^18")
print(f"  每片权重若取 u:      (2^18)^(1/6) = 8")
print(f"  每片权重若取 c_f=Nu: 再乘 N=2  ⇒ 16")
print("  ⇒ 归一化只在 'u vs c_f=Nu' 这一处差 N=2；语料 C̄=c_f/u=N 已把该 N 计入，故报 8（见 KP1.md §3）")
print()
print("done")
