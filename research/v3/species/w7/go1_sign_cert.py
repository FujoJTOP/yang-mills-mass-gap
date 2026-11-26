"""
GO1 — 号定理的精确证书（不跑 MC，纯整数算术）

靶：phi^U(Gamma) = sum_{A subset E(G_Gamma), (Gamma,A) connected} (-1)^|A|
                 = mu(G_Gamma) = (-1)^(n-1) * #{NBC spanning trees of G_Gamma}
    其中 NBC = 不含破圈（broken circuit）的生成树，边序 = 排序后的边表次序。

三段独立算法（互为对照，任一不一致即报）：
  (1) mu_bf   : 2^m 子集枚举 + 并查集            （定义直算，m<=20）
  (2) mu_dp   : 3^n 子集 DP（与 ur1_enum.conn_sum 同式，独立重写）
  (3) nbc_txt : 生成树全枚举 + 破圈判据          （Whitney 侧，只对小 m）
另加 tau       : Kirchhoff 整数行列式，给 |mu| <= tau。

检验对象：
  A. 4D plaquette 图（Delta=20）：n<=4 全枚举，n=5..10 抽样
  B. 2D plaquette 图（4 邻居）：n<=6 全枚举（酸测）
  C. 任意连通图（K_n / 圈 / 树 / 随机图）：n<=7 全对照，n=8..10 只查号
  D. 约定 B（f_e=u）：P(u)=sum_A(-1)^|A| u^|A| 在 (0,1] 的符号 + 首个正根

用法： py -3.12 go1_sign_cert.py [--json OUT]
"""
import sys, os, json, time, random, itertools
from fractions import Fraction
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ur1_enum as UR                      # 只借用 4D 邻接（neighbors/pc），不借用其统计

FAIL = []
FACT = [1] + [__import__("math").factorial(k) for k in range(1, 25)]


# ---------------------------------------------------------------- 通用图工具
def edges_of(adj, n):
    return [(i, j) for i in range(n) for j in range(i + 1, n) if adj[i] >> j & 1]


def mu_bf(n, edges):
    """mu = sum_{A connected spanning} (-1)^|A|   （定义直算）"""
    m = len(edges)
    tot = 0
    for mask in range(1 << m):
        par = list(range(n))

        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a

        k, pc = n, 0
        for i in range(m):
            if mask >> i & 1:
                pc += 1
                u, v = edges[i]
                ru, rv = find(u), find(v)
                if ru != rv:
                    par[ru] = rv
                    k -= 1
        if k == 1:
            tot += -1 if pc & 1 else 1
    return tot


def mu_dp(adj, n):
    """3^n 子集 DP；与 ur1_enum.conn_sum(x=1) 同式，独立重写。"""
    M = 1 << n
    ec = [0] * M
    for T in range(1, M):
        low = T & -T
        i = low.bit_length() - 1
        R = T ^ low
        ec[T] = ec[R] + bin(adj[i] & R).count('1')
    g = [1 if ec[T] == 0 else 0 for T in range(M)]        # (1-1)^e(T)
    f = [0] * M
    for T in range(1, M):
        low = T & -T
        rest = T ^ low
        tot = g[T]
        U = rest
        while True:
            if U != rest:
                P = low | U
                tot -= f[P] * g[T ^ P]
            if U == 0:
                break
            U = (U - 1) & rest
        f[T] = tot
    return f[M - 1]


def tau_kirchhoff(n, edges):
    """生成树数（Matrix-Tree，整数/有理精确）"""
    if n == 1:
        return 1
    L = [[0] * n for _ in range(n)]
    for (u, v) in edges:
        L[u][u] += 1
        L[v][v] += 1
        L[u][v] -= 1
        L[v][u] -= 1
    d = n - 1
    A = [[Fraction(L[i][j]) for j in range(d)] for i in range(d)]
    det = Fraction(1)
    for c in range(d):
        p = next((r for r in range(c, d) if A[r][c] != 0), None)
        if p is None:
            return 0
        if p != c:
            A[c], A[p] = A[p], A[c]
            det = -det
        det *= A[c][c]
        inv = A[c][c]
        for r in range(c + 1, d):
            if A[r][c] != 0:
                f = A[r][c] / inv
                for k in range(c, d):
                    A[r][k] -= f * A[c][k]
    return int(det)


def _tree_path_edges(adjT, n, u0):
    par = {u0: (None, None)}
    stack = [u0]
    while stack:
        x = stack.pop()
        for (y, ei) in adjT[x]:
            if y not in par:
                par[y] = (x, ei)
                stack.append(y)
    return par


def nbc_tree_count(n, edges):
    """生成树全枚举 + 破圈判据（边序 = edges 列表次序，已排序）"""
    m = len(edges)
    if n == 1:
        return 1
    cnt = 0
    for comb in itertools.combinations(range(m), n - 1):
        par = list(range(n))

        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a

        ok = True
        for ei in comb:
            u, v = edges[ei]
            ru, rv = find(u), find(v)
            if ru == rv:
                ok = False
                break
            par[ru] = rv
        if not ok:
            continue
        T = set(comb)
        adjT = [[] for _ in range(n)]
        for ei in T:
            u, v = edges[ei]
            adjT[u].append((v, ei))
            adjT[v].append((u, ei))
        nbc = True
        for e, (u0, v0) in enumerate(edges):
            if e in T:
                continue
            pd = _tree_path_edges(adjT, n, u0)
            path = []
            x = v0
            while x != u0:
                px, pi = pd[x]
                path.append(pi)
                x = px
            if min([e] + path) == e:        # e 是其基本圈的最小元 => 破圈在 T 里
                nbc = False
                break
        if nbc:
            cnt += 1
    return cnt


def poly_conn(n, edges):
    """约定 B 的多项式系数：P(u) = sum_{A connected} (-1)^|A| u^|A|"""
    m = len(edges)
    co = [0] * (m + 1)
    for mask in range(1 << m):
        par = list(range(n))

        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a

        k, pc = n, 0
        for i in range(m):
            if mask >> i & 1:
                pc += 1
                u, v = edges[i]
                ru, rv = find(u), find(v)
                if ru != rv:
                    par[ru] = rv
                    k -= 1
        if k == 1:
            co[pc] += -1 if pc & 1 else 1
    return co


def first_pos_root(co, hi=50.0):
    deg = len(co) - 1
    while deg > 0 and co[deg] == 0:
        deg -= 1
    co = co[:deg + 1]
    if deg == 0:
        return None
    r = np.roots(list(reversed(co)))
    pos = sorted(x.real for x in r if abs(x.imag) < 1e-9 and x.real > 1e-12)
    return pos[0] if pos else None


def poly_val(co, u):
    return sum(c * u ** k for k, c in enumerate(co))


# ---------------------------------------------------------------- 检验主体
def check_graph(name, n, adj, tag, want_nbc, want_poly, u, res):
    edges = edges_of(adj, n)
    m = len(edges)
    mu = mu_dp(adj, n)
    sgn = 1 if mu > 0 else (-1 if mu < 0 else 0)
    rec = dict(name=name, n=n, tag=tag, m=m, mu=mu, sign_ok=(sgn == ((-1) ** (n - 1))))
    if not rec["sign_ok"]:
        FAIL.append(f"{name} n={n} mu={mu}")
    if m <= 20:
        mb = mu_bf(n, edges)
        rec["mu_bf"] = mb
        if mb != mu:
            FAIL.append(f"DP!=BF {name} n={n} {mu} vs {mb}")
    tt = tau_kirchhoff(n, edges)
    rec["tau"] = tt
    if abs(mu) > tt:
        FAIL.append(f"|mu|>tau {name} n={n} {mu} {tt}")
    if want_nbc and n <= 7:
        nc = nbc_tree_count(n, edges)
        rec["nbc"] = nc
        rec["nbc_ok"] = (abs(mu) == nc)
        if not rec["nbc_ok"]:
            FAIL.append(f"|mu|!=NBC {name} n={n} {mu} vs {nc}")
    if want_poly and m <= 16:
        co = poly_conn(n, edges)
        rec["poly"] = co
        rec["P_u"] = poly_val(co, u)
        rec["P_1"] = poly_val(co, 1.0)
        s_u = 1 if rec["P_u"] > 0 else (-1 if rec["P_u"] < 0 else 0)
        rec["sign_u_ok"] = (s_u == ((-1) ** (n - 1)))
        if not rec["sign_u_ok"]:
            FAIL.append(f"u-conv sign {name} n={n} P={rec['P_u']}")
        rr = first_pos_root(co)
        rec["root"] = rr
    res.append(rec)
    return rec


def main():
    out_json = None
    if "--json" in sys.argv:
        out_json = sys.argv[sys.argv.index("--json") + 1]
    t0 = time.time()
    rng = random.Random(20261126)
    u_phys = 0.625

    # ---------- C. 具名图（定理的通用性） ----------
    print("=== C. 具名/随机图：mu vs NBC vs tau ===")
    gres = []
    for k in range(1, 8):
        adj = [0] * k
        for i in range(k):
            for j in range(k):
                if i != j:
                    adj[i] |= 1 << j
        check_graph(f"K{k}", k, adj, "complete", True, False, u_phys, gres)
        r = gres[-1]
        print(f"  K{k}: mu={r['mu']:>8d} NBC={r.get('nbc')} tau={r['tau']:>7d} "
              f"| mu/(n-1)! = {abs(r['mu'])}/{np.prod(range(1,k)) if k>1 else 1}")
    # 圈
    for k in range(3, 9):
        adj = [0] * k
        for i in range(k):
            adj[i] |= 1 << ((i + 1) % k)
            adj[(i + 1) % k] |= 1 << i
        check_graph(f"C{k}", k, adj, "cycle", True, False, u_phys, gres)
        print(f"  C{k}: mu={gres[-1]['mu']} NBC={gres[-1].get('nbc')} (expect (-1)^(n-1)*(n-1)={(-1)**(k-1)*(k-1)})")
    # 树（路径）
    for k in range(2, 9):
        adj = [0] * k
        for i in range(k - 1):
            adj[i] |= 1 << (i + 1)
            adj[i + 1] |= 1 << i
        check_graph(f"P{k}", k, adj, "tree", True, False, u_phys, gres)
    print(f"  路径 P2..P8: mu = {[r['mu'] for r in gres if r['tag']=='tree']} (expect (-1)^(n-1))")

    # 随机连通图
    nr = 0
    for n in range(5, 11):
        for trial in range(25):
            while True:
                adj = [0] * n
                for i in range(n):
                    for j in range(i + 1, n):
                        if rng.random() < 0.45:
                            adj[i] |= 1 << j
                            adj[j] |= 1 << i
                e = edges_of(adj, n)
                if len(e) <= 20 and mu_dp_is_conn(adj, n):
                    break
            check_graph(f"rand{n}_{trial}", n, adj, "random",
                        n <= 7, n <= 6, u_phys, gres)
            nr += 1
    print(f"  随机图 {nr} 个：全部 sign_ok={all(r['sign_ok'] for r in gres)}，"
          f"NBC 对照 {sum(1 for r in gres if 'nbc' in r)} 个，"
          f"不一致 {sum(1 for r in gres if r.get('nbc_ok') is False)}")

    # ---------- E. 证明两步的逐例执行（删-缩递归） ----------
    print("\n=== E. 删-缩递归：mu(G)=mu(G-e)-mu(G/e)，N(G)=N(G-e)+N(G/e)（e=最大边）===")

    def is_conn_edges(n, el):
        par = list(range(n))

        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a
        k = n
        for (u, v) in el:
            ru, rv = find(u), find(v)
            if ru != rv:
                par[ru] = rv
                k -= 1
        return k == 1

    def contract(n, edges, ei):
        u, v = edges[ei]
        keep = [x for i, x in enumerate(edges) if i != ei]
        mp = lambda w: u if w == v else w
        ne = []
        for (a, b) in keep:
            a, b = mp(a), mp(b)
            if a == b:
                continue
            ne.append((a, b))
        verts = sorted({x for ab in ne for x in ab})
        if len(verts) != n - 1:
            return None
        rel = {w: i for i, w in enumerate(verts)}
        out = sorted(tuple(sorted((rel[a], rel[b]))) for (a, b) in ne)
        if len(set(out)) != len(out):
            return None                       # 出现重边 => 跳过（本测试只走简单图）
        return out

    nE = 0
    for n in range(4, 9):
        for trial in range(40):
            while True:
                adj = [0] * n
                for i in range(n):
                    for j in range(i + 1, n):
                        if rng.random() < 0.5:
                            adj[i] |= 1 << j
                            adj[j] |= 1 << i
                el = edges_of(adj, n)
                if 3 <= len(el) <= 16 and mu_dp_is_conn(adj, n):
                    break
            mG = mu_bf(n, el)
            # 最大边（edge 列表已按 (i,j) 升序，取末位）
            ei = len(el) - 1
            Ge = el[:ei] + el[ei + 1:]
            if not is_conn_edges(n, Ge):
                continue                      # e 是桥：走块积，不在本测试
            con = contract(n, el, ei)
            if con is None:
                continue
            m1 = mu_bf(n, Ge)
            m2 = mu_bf(n - 1, con)
            ok1 = (mG == m1 - m2)
            N_all = nbc_tree_count(n, el) if n <= 7 else None
            N1 = nbc_tree_count(n, Ge) if n <= 7 else None
            N2 = nbc_tree_count(n - 1, con) if n <= 7 else None
            ok2 = (N_all == N1 + N2) if n <= 7 else None
            nE += 1
            if not ok1 or (ok2 is False):
                FAIL.append(f"recursion n={n}: mu {mG} vs {m1}-{m2}; N {N_all} vs {N1}+{N2}")
    print(f"  {nE} 个 (图,最大非桥边) 实例：mu 递归全对 = {not any('recursion' in f for f in FAIL)}；"
          f"N 递归（n<=7）= {not any('recursion' in f for f in FAIL)}")

    # ---------- A. 4D plaquette 图 ----------
    print("\n=== A. 4D plaquette 图（Delta=20）===")
    ares = []
    levels = UR.animals(4)                      # n<=4 全枚举（1+20+458+11132）
    for n in range(1, 5):
        t1 = time.time()
        for S in levels[n - 1]:
            adj = UR.graph_data(S)
            check_graph("4D", n, adj, "exhaustive", True, True, u_phys, ares)
        sub = [r for r in ares if r["n"] == n]
        print(f"  n={n} 全枚举 {len(sub)} 簇: sign-viol={sum(1 for r in sub if not r['sign_ok'])} "
              f"|mu|!=NBC {sum(1 for r in sub if r.get('nbc_ok') is False)} "
              f"max|mu|={max(abs(r['mu']) for r in sub)} [... {time.time()-t1:.1f}s]")

    # 抽样 n=5..10
    root = UR.pc((0, 0, 0, 0), 0, 1)
    print("  抽样（从根随机生长，去重）：")
    for n in range(5, 11):
        t1 = time.time()
        seen, got = set(), 0
        want = 120 if n <= 8 else 60
        while got < want:
            S = {root}
            bnd = set(UR.neighbors(root))
            while len(S) < n and bnd:
                q = rng.choice(sorted(bnd))
                S.add(q)
                bnd.discard(q)
                for rr in UR.neighbors(q):
                    if rr not in S:
                        bnd.add(rr)
            if len(S) != n:
                continue
            key = frozenset(S)
            if key in seen:
                continue
            seen.add(key)
            got += 1
            adj = UR.graph_data(tuple(sorted(key)))
            check_graph("4D", n, adj, "sampled", n <= 7, n <= 8, u_phys, ares)
        sub = [r for r in ares if r["n"] == n and r["tag"] == "sampled"]
        print(f"  n={n:>2} 抽样 {len(sub):>3} 簇: sign-viol={sum(1 for r in sub if not r['sign_ok'])} "
              f"|mu| max={max(abs(r['mu']) for r in sub):>6d} "
              f"|mu| 中位={int(np.median([abs(r['mu']) for r in sub])):>4d} "
              f"NBC对照={sum(1 for r in sub if 'nbc' in r)} "
              f"min root={min([r['root'] for r in sub if r.get('root')] or [float('nan')]):.3f} "
              f"[{time.time()-t1:.1f}s]")

    # 两个具名大簇（复现 C02:92-96 / GV1 §2.2）
    print("\n  具名簇复现：")
    # (a) K6：过同一条 link 的 6 片
    def plq(x, mu, nu):
        return UR.pc(x, mu, nu)
    K6 = frozenset([plq((0, 0, 0, 0), 0, 1), plq((0, 0, 0, 0), 0, 2), plq((0, 0, 0, 0), 0, 3),
                    plq((0, -1, 0, 0), 0, 1), plq((0, 0, -1, 0), 0, 2), plq((0, 0, 0, -1), 0, 3)])
    check_graph("K6_cluster", 6, UR.graph_data(tuple(sorted(K6))), "named", True, False, u_phys, ares)
    r = ares[-1]
    print(f"    K6（一 link 上 6 片）: |E|={r['m']} mu={r['mu']} NBC={r.get('nbc')} tau={r['tau']}"
          f"  [GV1:120, tau=1296]")
    # (b) 3-cube 的 6 个面
    CU = frozenset([plq((0, 0, 0, 0), 1, 2), plq((0, 0, 0, 1), 1, 2),     # 面 ⟂ e3
                    plq((0, 0, 0, 0), 1, 3), plq((0, 0, 1, 0), 1, 3),     # 面 ⟂ e2
                    plq((0, 0, 0, 0), 2, 3), plq((0, 1, 0, 0), 2, 3)])    # 面 ⟂ e1
    check_graph("cube6", 6, UR.graph_data(tuple(sorted(CU))), "named", True, False, u_phys, ares)
    r = ares[-1]
    print(f"    3-cube 六面: |E|={r['m']} mu={r['mu']} NBC={r.get('nbc')} tau={r['tau']}"
          f"  [C02:96: mu=-64, tau=384]")

    # ---------- B. 2D plaquette 图（酸测） ----------
    print("\n=== B. 2D plaquette 图（4 邻居）— 酸测 ===")
    bres = []

    def nb2(p):
        x, y = p
        return [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]

    def animals2(maxn):
        lv = [{frozenset({(0, 0)})}]
        for _ in range(1, maxn):
            nxt = set()
            for S in lv[-1]:
                for p in S:
                    for q in nb2(p):
                        if q not in S:
                            nxt.add(frozenset(S | {q}))
            lv.append(nxt)
        return [sorted(L) for L in lv]

    def adj2(S):
        S = list(S)
        idx = {p: i for i, p in enumerate(S)}
        adj = [0] * len(S)
        for p, i in idx.items():
            for q in nb2(p):
                if q in idx:
                    adj[i] |= 1 << idx[q]
        return adj

    lv2 = animals2(6)
    for n in range(1, 7):
        for S in lv2[n - 1]:
            check_graph("2D", n, adj2(S), "exhaustive", True, True, u_phys, bres)
        sub = [r for r in bres if r["n"] == n]
        print(f"  n={n} 全枚举 {len(sub):>4} 个含原点多联骨牌（每个平移类出现 n 次）: sign-viol="
              f"{sum(1 for r in sub if not r['sign_ok'])} "
              f"|mu|!=NBC {sum(1 for r in sub if r.get('nbc_ok') is False)} "
              f"|mu| 取值={sorted({abs(r['mu']) for r in sub})[:8]}")
    print("  ⇒ 2D 的号 = (-1)^(n-1)，与 4D 同；逐项全正在**特征标/曲面**表象，不在此表象。")

    # ---------- D. 约定 B ----------
    print("\n=== D. 约定 B：P(u)=sum_A(-1)^|A|u^|A| 的符号窗口 ===")
    roots = [r["root"] for r in ares if r.get("root")]
    roots2 = [r["root"] for r in bres if r.get("root")]
    print(f"  4D n<=10 抽样+全枚举 {len(roots)} 簇：min 首个正根 = {min(roots):.6f}；"
          f"P(0.625) 与 P(1) 符号全对 = "
          f"{all(r.get('sign_u_ok', True) for r in ares)}")
    print(f"  2D {len(roots2)} 簇：min 首个正根 = {min(roots2):.6f}")
    print(f"  理论对照：C3 在 u=3 翻号（P=3u^2-u^3）⇒ u-约定**不是**定理，只是窗口事实")

    # ---------- F. 路由 1 的恒等式 + 多线性 log Z 端到端 ----------
    print("\n=== F. 路由 1（分拆格/多线性对数）：phi^U(Γ) = Σ_π (-1)^(|π|-1)(|π|-1)! [各块独立] ===")

    def partitions_of(lst):
        if not lst:
            yield []
            return
        first, rest = lst[0], lst[1:]
        for part in partitions_of(rest):
            yield [[first]] + part
            for i in range(len(part)):
                yield part[:i] + [[first] + part[i]] + part[i + 1:]

    def route1_mu(adj, n):
        tot = 0
        for part in partitions_of(list(range(n))):
            ok = True
            for B in part:
                m = 0
                for v in B:
                    m |= 1 << v
                for v in B:
                    if adj[v] & m & ~(1 << v):
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                k = len(part)
                tot += (-1) ** (k - 1) * FACT[k - 1]
        return tot

    def sc_prod(f, g):
        out = {}
        for a, ca in f.items():
            for b, cb in g.items():
                if a & b:
                    continue
                m = a | b
                out[m] = out.get(m, 0) + ca * cb
        return out

    def mv_log_check(name, adj):
        """Z(w) = Σ_{independent S} ∏ w_p；log Z 的多线性系数 ?= φ^U(Γ)"""
        nv = len(adj)
        Z = {mask: 1 for mask in range(1, 1 << nv)
             if all(not (adj[i] & mask) for i in range(nv) if mask >> i & 1)}
        N = dict(Z)
        logZ = {}
        P = dict(N)
        for k in range(1, nv + 1):
            coef = Fraction((-1) ** (k + 1), k)
            for m, c in P.items():
                logZ[m] = logZ.get(m, Fraction(0)) + coef * c
            P = sc_prod(P, N)
        bad = 0
        for mask in range(1, 1 << nv):
            n = bin(mask).count('1')
            vs = [i for i in range(nv) if mask >> i & 1]
            idx = {v: i for i, v in enumerate(vs)}
            a = [0] * n
            for v in vs:
                for w in vs:
                    if v != w and (adj[v] >> w & 1):
                        a[idx[v]] |= 1 << idx[w]
            conn = mu_dp_is_conn(a, n)
            mu = mu_dp(a, n) if conn else 0
            got = logZ.get(mask, Fraction(0))
            if got != mu:
                bad += 1
                FAIL.append(f"logZ {name} mask={mask}: {got} vs {mu}")
        print(f"  {name}: nv={nv}，2^{nv}-1 个 Γ 的多线性系数逐个比对：不一致 {bad}")
        return bad == 0

    # 2D 3x3 torus / 4D 1^4 torus
    def torus2(L):
        verts = [(x, y) for x in range(L) for y in range(L)]
        idx = {p: i for i, p in enumerate(verts)}
        adj = [0] * len(verts)
        for p, i in idx.items():
            for d in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
                q = ((p[0] + d[0]) % L, (p[1] + d[1]) % L)
                adj[i] |= 1 << idx[q]
        return adj

    adj2_3 = torus2(3)
    PL = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    adj4 = [0] * 6
    for i in range(6):
        for j in range(6):
            if i != j and set(PL[i]) & set(PL[j]):
                adj4[i] |= 1 << j
    mv_log_check("2D 3x3 torus", adj2_3)
    mv_log_check("4D 1^4 torus", adj4)

    # 路由 1 恒等式：4D n<=4 全枚举 + 2D n<=4 全枚举 + 抽样
    n1 = 0
    for n in range(1, 5):
        for S in levels[n - 1]:
            adj = UR.graph_data(S)
            if route1_mu(adj, n) != mu_dp(adj, n):
                FAIL.append(f"route1 4D n={n}")
            n1 += 1
    for n in range(1, 5):
        for S in lv2[n - 1]:
            adj = adj2(S)
            if route1_mu(adj, n) != mu_dp(adj, n):
                FAIL.append(f"route1 2D n={n}")
            n1 += 1
    print(f"  路由 1 恒等式在 {n1} 个簇（4D+2D，n<=4 全枚举）上逐簇相等 = "
          f"{not any('route1' in f for f in FAIL)}")
    print("  注：log Z 的**单变量**系数不是簇和（差在重复聚合物项）——"
          "必须是**多线性**系数；这一条是本轮实测到的记账陷阱。")

    # ---------- 汇总 ----------
    dt = time.time() - t0
    print(f"\n=== 汇总 ===  用时 {dt:.1f}s")
    print(f"FAIL 计数 = {len(FAIL)}")
    for f in FAIL:
        print("  !", f)
    if out_json:
        with open(out_json, "w") as fh:
            json.dump(dict(fail=FAIL, graphs=gres, f4d=ares, f2d=bres,
                           elapsed=dt), fh, indent=1)
        print(f"json -> {out_json}")


def mu_dp_is_conn(adj, n):
    seen = {0}
    st = [0]
    while st:
        x = st.pop()
        for y in range(n):
            if adj[x] >> y & 1 and y not in seen:
                seen.add(y)
                st.append(y)
    return len(seen) == n


if __name__ == "__main__":
    main()
