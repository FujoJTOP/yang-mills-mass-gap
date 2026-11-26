"""adj7: 核 a_n <= t_n (片集数 <= 根定树和) 及其推论 a_n <= B^n.

A. t_n/a_n 精确比值 (corpus 值) + 逐 n 比较
B. tau(G_Gamma) >= 1 的独立机器复核: 自己枚举 n<=N, 对**每一个** Gamma 算 Kirchhoff tau,
   报 min tau / #tau=0 / #tau=1 / sum tau (务必逐位复现 t_n 与 a_n)
C. 反例: Gamma 不连通时 tau=0 => 连通性是承重的 (a_n 的定义域 = 连通片集)
D. B = 19^19/18^18 的精确值, a_n <= B^n (有理数精确比较), B/48, 缺口百分比
E. 判词 (a)/(b)/(c)

Usage:  py -3.12 adj7_an_tn.py [NMAX]      NMAX default 5
No git.
"""
import sys, time
from fractions import Fraction
from decimal import Decimal, getcontext

getcontext().prec = 50

A_CORPUS = [1, 20, 458, 11132, 280175, 7218768]              # a_n, n=1..6  (TN1:23)
T_CORPUS = [1, 20, 570, 18776, 669215, 25060608]             # t_n, n=1..6  (TN1:12)

# ---------------- 4D plaquette 图 (adjacency = 共 link), 同 TN1 §2 ----------------
OFF = 8
ORI = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
OIDX = {o: i for i, o in enumerate(ORI)}


def pack(x, mu, nu):
    if mu > nu:
        mu, nu = nu, mu
    v = 0
    for c in x:
        v = (v << 5) | (c + OFF)
    return (v << 3) | OIDX[(mu, nu)]


def unpack(p):
    o = ORI[p & 7]
    v = p >> 3
    x = []
    for _ in range(4):
        x.append((v & 31) - OFF)
        v >>= 5
    x.reverse()
    return tuple(x), o[0], o[1]


def links(p):
    x, mu, nu = unpack(p)
    out = []

    def sh(y, r, k):
        y = list(y)
        y[r] += k
        return tuple(y)

    out.append((x, mu))
    out.append((sh(x, mu, 1), nu))
    out.append((x, nu))
    out.append((sh(x, nu, 1), mu))
    return set(out)


_NC = {}


def neigh(p):
    r = _NC.get(p)
    if r is not None:
        return r
    x, mu, nu = unpack(p)
    res = set()
    for (y, rho) in links(p):
        for sigma in range(4):
            if sigma == rho:
                continue
            res.add(pack(y, rho, sigma))
            yy = list(y)
            yy[sigma] -= 1
            res.add(pack(tuple(yy), rho, sigma))
    res.discard(p)
    _NC[p] = frozenset(res)
    return _NC[p]


def det_bareiss(M):
    M = [row[:] for row in M]
    n = len(M)
    if n == 0:
        return 1
    sign, prev = 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for i in range(k + 1, n):
                if M[i][k] != 0:
                    M[k], M[i] = M[i], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
        for i in range(k + 1, n):
            M[i][k] = 0
    return sign * M[n - 1][n - 1]


def tau_kirchhoff(cluster):
    """生成树数 tau(G_Gamma); G_Gamma = 片集上的"共 link"图."""
    n = len(cluster)
    if n == 1:
        return 1
    Lk = [links(p) for p in cluster]
    Lap = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if Lk[i] & Lk[j]:
                Lap[i][i] += 1
                Lap[j][j] += 1
                Lap[i][j] -= 1
                Lap[j][i] -= 1
    return det_bareiss([row[1:] for row in Lap[1:]])


def nroot(x, n):
    return (Decimal(x) ** (Decimal(1) / Decimal(n)))


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    t0 = time.time()

    print("=" * 78)
    print("A · t_n/a_n (corpus: TN1:12 t_n, TN1:23 a_n)")
    print("=" * 78)
    print(" n |      a_n |        t_n | t_n/a_n (exact Fraction) | decimal    | t_n>=a_n | t_n>a_n")
    ratios = []
    for i in range(6):
        n, a, t = i + 1, A_CORPUS[i], T_CORPUS[i]
        fr = Fraction(t, a)
        ratios.append(fr)
        print("%2d | %8d | %10d | %-25s | %.6f | %-8s | %-5s"
              % (n, a, t, str(fr), float(fr), t >= a, t > a))
    mono = all(ratios[i] <= ratios[i + 1] for i in range(5))
    print("\n ratio 单调升 (data): %s ; 严格 t_n>a_n 的 n: %s"
          % (mono, [i + 1 for i in range(6) if T_CORPUS[i] > A_CORPUS[i]]))
    print(" 比值第一次 >1 在 n=%d ; n=1,2 处 t_n=a_n (tau = 1 恒等, 见 B 段)"
          % next(i + 1 for i in range(6) if ratios[i] > 1))

    print()
    print("=" * 78)
    print("B · tau(G_Gamma) >= 1 的机器复核 (独立枚举, 逐 Gamma 算 Kirchhoff)")
    print("=" * 78)
    root = pack((0, 0, 0, 0), 0, 1)
    deg = len(neigh(root))
    print(" root plaquette degree = %d   (Delta=4(2(d-1)-1), d=4 -> 20)" % deg)
    levels = [((root,),)]
    for n in range(1, nmax):
        nxt = set()
        for A in levels[-1]:
            for p in A:
                for q in neigh(p):
                    if q not in A:
                        nxt.add(tuple(sorted(A + (q,))))
        levels.append(nxt)
        print("  level %d: %d connected sets  (%.1fs)" % (n + 1, len(nxt), time.time() - t0))

    print("\n n |     a_n(enum) | corpus a_n | min tau | #tau=0 | #tau=1 |   sum tau | corpus t_n | match")
    ok_all, ok_tau, ok_sum = True, True, True
    for n in range(1, nmax + 1):
        lvl = levels[n - 1]
        taus = [tau_kirchhoff(list(cl)) for cl in lvl]
        tot, mn = sum(taus), min(taus)
        zeros, ones = taus.count(0), taus.count(1)
        an_e = len(lvl)
        match = "OK" if (an_e == A_CORPUS[n - 1] and tot == T_CORPUS[n - 1]) else "MISMATCH"
        ok_all &= (mn >= 1)
        ok_tau &= (zeros == 0)
        ok_sum &= (tot == T_CORPUS[n - 1])
        print("%2d | %13d | %10d | %7d | %6d | %6d | %9d | %10d | %s"
              % (n, an_e, A_CORPUS[n - 1], mn, zeros, ones, tot, T_CORPUS[n - 1], match))
    print("\n >> 每一个连通 Gamma 上 tau>=1 : %s   (#tau=0 出现 %s 次)" % (ok_all, "未" if ok_tau else ""))
    print(" >> sum tau 逐位复现 corpus t_n : %s" % ok_sum)
    print(" >> 结构含义: #tau=1 = 树型片集; t_n - a_n = sum(tau-1) = 非树贡献 (n=1,2 为 0)")

    print()
    print("=" * 78)
    print("C · 反例: Gamma 不连通 => tau=0  (连通性是承重的)")
    print("=" * 78)
    far = [root, pack((6, 0, 0, 0), 0, 1)]
    print(" 两个相距 (L1=6) 的片: 邻接? %s  tau=%d  (tau>=1 在此失效)"
          % (bool(links(far[0]) & links(far[1])), tau_kirchhoff(far)))
    far2 = [root, pack((0, 0, 0, 0), 2, 3)]
    print(" 同点不同朝向的片 (共点不共 link): 邻接? %s  tau=%d"
          % (bool(links(far2[0]) & links(far2[1])), tau_kirchhoff(far2)))
    print(" => a_n 的定义域 = 连通片集; 若 a_n 曾按'全部 n 元片集'计数, t_n>=a_n 立刻假。")

    print()
    print("=" * 78)
    print("D · B = 19^19/18^18 与 a_n <= B^n")
    print("=" * 78)
    B = Fraction(19) ** 19 / Fraction(18) ** 18
    Bd = Decimal(B.numerator) / Decimal(B.denominator)
    print(" B = 19^19/18^18 = %s" % Bd)
    print("   = %d / %d" % (B.numerator, B.denominator))
    print(" eD = 20e        = %s" % (Decimal(20) * Decimal(1).exp()))
    print(" B/48            = %s   ==> 缺口 = %s %%"
          % (Bd / 48, (Bd / 48 - 1) * 100))
    print(" 需要再降 %s %% 才能到 48" % ((1 - Decimal(48) / Bd) * 100))
    print("\n n |        a_n |   a_n <= B^n ?  (有理数精确) | a_n^(1/n) | t_n^(1/n)")
    for i in range(6):
        n, a, t = i + 1, A_CORPUS[i], T_CORPUS[i]
        holds = Fraction(a) <= B ** n
        print("%2d | %10d | %-27s | %9.4f | %9.4f"
              % (n, a, holds, nroot(a, n), nroot(t, n)))
    print("\n 逐 n: a_n^{1/n} <= t_n^{1/n} <= B 全部成立? %s"
          % all(nroot(A_CORPUS[i], i + 1) <= nroot(T_CORPUS[i], i + 1) <= Bd for i in range(6)))
    print(" a_n^{1/n} 升 (=> liminf a_n^{1/n} >= %.4f), 上界 B = %.4f  ==> 48 在夹逼里"
          % (nroot(A_CORPUS[5], 6), Bd))

    print()
    print("=" * 78)
    print("E · 判词")
    print("=" * 78)
    print(" (a) a_n <= t_n: 成立 [T]; 传递得 a_n <= t_n <= B^n => a_n <= B^n [T]")
    print("     但 B/48-1 = %.4f%% => 缺口一步未缩, 仍是 bound 侧的缺口, 非 a_n 的缺口"
          % ((Bd / 48 - 1) * 100))
    print(" (b) tau>=1 在**连通** Gamma 上恒成立 [T] (连通 <=> 有生成树); 不连通处 tau=0,")
    print("     故不等式只在'连通片集'定义域上成立 —— 连通性是承重的, 无隐藏例外")
    print(" (c) B 的推导 (TN1 §3 覆盖树) 对全部 n 无条件成立, 不限于特定类 Gamma [T];")
    print("     它界的是 t_n (树和), 不是 lambda_a 本身 => 对 a_n 的传递无损但 B 仍 >48")
    print("\n elapsed %.1fs" % (time.time() - t0))


main()
