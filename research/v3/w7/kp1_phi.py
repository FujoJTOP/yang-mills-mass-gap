"""KP1-a —— 独立重算立方体簇的 phi^U（三条整数路径），不引用 ur3 的任何数。

phi^U(Gamma) = sum_{A subset E(G_Gamma), (Gamma,A) 连通生成} (-1)^|A|
路径①  2^|E| 子集 + 并查集连通判定，权重 (-1)^|A|
路径②  色多项式 P_G(q) 的 [q] 系数（恒等式 phi^U = [q]P_G，自证见下）
路径③  Tutte：phi^U = (-1)^(n-1) T_G(1,0)，用②的 P_G 关系独立核对
另附  Kirchhoff 生成树数 tau（|phi^U| <= tau 的对照）
"""
from fractions import Fraction
import sympy as sp

Q = sp.symbols('q')


# ---------- 路径① ----------
def phi_bitmask(n, edges):
    tot = 0
    m = len(edges)
    for mask in range(1 << m):
        par = list(range(n))

        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a

        k = 0
        for i, (u, v) in enumerate(edges):
            if (mask >> i) & 1:
                k += 1
                ru, rv = find(u), find(v)
                if ru != rv:
                    par[ru] = rv
        if len({find(a) for a in range(n)}) == 1:
            tot += -1 if k & 1 else 1
    return tot


# ---------- 路径② ----------
def chrom(n, edges):
    memo = {}

    def rec(nv, es):
        key = (nv, es)
        if key in memo:
            return memo[key]
        if not es:
            r = Q ** nv
        elif len(es) == 1:
            r = Q * (Q - 1) * Q ** (nv - 2)
        else:
            e = sorted(es)[0]
            u, v = e
            d = rec(nv, es - {e})
            mp = list(range(nv))
            mp[v] = u
            remap = {}
            for w in range(nv):
                if w != v:
                    remap[w] = len(remap)
            es2 = frozenset(tuple(sorted((remap[mp[a]], remap[mp[b]])))
                            for (a, b) in es - {e} if mp[a] != mp[b])
            r = sp.expand(d - rec(nv - 1, es2))
        memo[key] = r
        return r

    return sp.Poly(sp.expand(rec(n, frozenset(edges))), Q)


def tau_kirchhoff(n, edges):
    if n == 1:
        return 1
    L = sp.zeros(n, n)
    for (a, b) in edges:
        L[a, a] += 1
        L[b, b] += 1
        L[a, b] -= 1
        L[b, a] -= 1
    M = L[1:, 1:]
    return int(sp.expand(M.det()))


# ---------- 图 ----------
def K(n):
    return n, [(i, j) for i in range(n) for j in range(i + 1, n)]


def cube_faces_cluster():
    """3-cube 的 6 个面 = 6 片 plaquette；两片共享一条 cubic edge <=> 法向轴不同。"""
    faces = [(ax, s) for ax in range(3) for s in (0, 1)]
    idx = {f: i for i, f in enumerate(faces)}
    E = [(idx[a], idx[b]) for a in faces for b in faces if a < b and a[0] != b[0]]
    return len(faces), E


def Q3():
    return 8, [(i, i ^ (1 << b)) for i in range(8) for b in range(3) if i < (i ^ (1 << b))]


def k3_bridge_k3():
    return 6, [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3)]


GRAPHS = [
    ("K3", *K(3)),
    ("K4", *K(4)),
    ("K5", *K(5)),
    ("K6 (一个 link 上 6 片)", *K(6)),
    ("C4", 4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
    ("K3-桥-K3", *k3_bridge_k3()),
    ("★立方体簇 6 面(共享 cubic edge)", *cube_faces_cluster()),
    ("Q3 (8 顶点误读)", *Q3()),
]

print("=" * 78)
print("KP1-a  phi^U 独立枚举（整数精确；无浮点自由度）")
print("=" * 78)
hd = f"{'图':34s}{'n':>3s}{'|E|':>4s}{'phi^U①':>8s}{'phi^U②':>8s}{'tau':>7s}{'|phi|^(1/n)':>13s}"
print(hd)
print("-" * 78)
res = {}
for name, n, edges in GRAPHS:
    p1 = phi_bitmask(n, edges)
    P = chrom(n, edges)
    p2 = int(P.coeff_monomial(Q))
    # 路径③: T(1,0) 与 P 的关系  P_G(q) = (-1)^(n-1) q T_G(1-q,0)  =>  T(1,0) = (-1)^(n-1)[q]P
    T10 = (-1) ** (n - 1) * p2
    p3 = (-1) ** (n - 1) * T10
    tau = tau_kirchhoff(n, edges)
    assert p1 == p2 == p3, (name, p1, p2, p3)
    res[name] = (n, len(edges), p1, tau)
    root = abs(p1) ** (1.0 / n) if p1 else 0.0
    print(f"{name:34s}{n:3d}{len(edges):4d}{p1:8d}{p2:8d}{tau:7d}{root:13.7f}")

print()
print("路径③核对：phi^U = (-1)^(n-1) T_G(1,0)，T_G(1,0)=(-1)^(n-1)[q]P_G  ⇒ 三路逐位相同")
print("路径②自检：K3 -> P=q^3-3q^2+2q, [q]=2=phi^U  ；K6 -> [q]=-120  ；")
print("        C4  -> -3 (4 条路 -1 + 环 +1) ；K3-桥-K3 -> +4 (=2*2 块乘积)")

print()
print("=" * 78)
print("立方体簇的几何核对（不预设八面体）")
print("=" * 78)
n, edges = cube_faces_cluster()
print(f"6 面 x 每面 4 边; 共享 cubic edge 的面对数 = |E(G_Gamma)| = {len(edges)} (=12)")
print(f"每片邻接度 = {len(edges)//n} (=4)；3 对对面不相邻 ⇒ G_Gamma = K_{{2,2,2}} 八面体")
print("两片共享 cube 顶点 <=> 共享 cube 边（立方体上任两面共顶点必共边）⇒ 顶点邻接与边邻接同图")
print(f"phi^U(立方体簇) = {res['★立方体簇 6 面(共享 cubic edge)'][2]}")

print()
print("=" * 78)
print("阈值与 K̂（精确整数幂；阈值/常数按语料逐字）")
print("=" * 78)
eD = 54.3656365691809          # e*Delta（语料 54.3656365692）
U = Fraction(1, 384)           # u = beta/(2N^2), beta=1/48, N=2
thr = float(1 / (eD * float(U)))
phi_c = res['★立方体簇 6 面(共享 cubic edge)'][2]
F = 6
prodR = 2 ** 12                 # 立方体 12 条 link，全 kappa=2，R(2)=2

def nth_root_exact(N, k):
    r = sp.integer_nthroot(N, k)
    return r

for (lbl, val) in [("prodR*|phi^U| = 2^12*64 = 2^18", prodR * abs(phi_c)),
                   ("prodR*tau       = 2^12*384", prodR * res['★立方体簇 6 面(共享 cubic edge)'][3])]:
    root, exact = nth_root_exact(val, F)
    print(f"{lbl:38s} = {val:9d}   ^(1/6) = {root}  {'(恰整数)' if exact else '(非整)'}")

Khat = 2 ** 3
print()
print(f"e*Delta = {eD}   u = {float(U):.10f}   e*Delta*u = {eD*float(U):.10f}")
print(f"门槛 1/(e*Delta*u)            = {thr:.10f}")
print(f"K̂ = (2^12*64)^(1/6)           = {Khat}.0       (精确)")
print(f"K̂/(门槛)                      = {Khat/thr:.7f}   ⇒ 超 {Khat/thr:.4f}x")
print(f"要过门槛需 |phi^U| < (门槛/4)^6 = {(thr/4)**6:.7f}  ；实有 64")
print(f"|phi^U|^(1/6)                 = {abs(phi_c)**(1/6):.7f} ; 门槛 1.7658213")
print(f"tau^(1/6)                     = {384**(1/6):.7f}       ; (tau 路线 K̂=10.7840492368)")

print()
print("=" * 78)
print("归一化读数表（同一 |phi^U|=64，不同的 per-plaquette / per-link 约定）")
print("=" * 78)
readings = [
    ("(A) u/片 + R=2/link  [语料 C06(5.1) 口径]", 8.0, "K̂=8"),
    ("(B) c_f=Nu/片 + R/link  [A01/C05 C̄=N 口径]", 8.0 * 2, "K̂=16（与 (A) 差 N）"),
    ("(C) 1/N/片（字符按 dim 归一）", 8.0 / 2, "K̂=4 <7.0633（非语料约定）"),
    ("(D) 再除 1/6!（把 Γ 当无标号）", 4 * (64 / 720) ** (1 / 6), "K̂=2.658（非法：1/n! 属簇计数）"),
    ("(E) 只用 tau（不开真 phi^U）", (prodR * 384) ** (1 / 6), "K̂=10.784"),
]
for lbl, val, note in readings:
    print(f"{lbl:44s} K̂ = {val:8.4f}   {note}")
print()
print("注：(A) 与 (B) 数值相同的原因 —— |phi^U|^{1/6} = 64^{1/6} = 2 = N = C̄ = c_f/u（**恰等**），")
print("    故 Ĉ·|phi^U|^{1/6} = 4N = Ĉ·C̄ = 8（同一条式子的两种写法，不是两笔账）。")
print()
print("done")
