"""
GV1 — exact Ursell shell masses on the 4D SU(2) plaquette graph (beta=2.5, u=0.625)

Definition (GU4-NECESSARY.md §4.5, 口径钉死):
    c_Gamma := phi^U(Gamma)/u^{|Gamma|} ,   lambda_Ur := limsup_Gamma |c_Gamma|^{1/|Gamma|}
    m_n     := sum_{Gamma: p0 in Gamma, |Gamma|=n} |c_Gamma|

Model actually computed here (= the corpus's own object, C02-URSELL (1.1) block product,
= leading term of the SU(2) character expansion: every plaquette fundamental, every
shared link contracted pairwise):
    phi^U(Gamma) = u^{|Gamma|} * mu(G_Gamma) ,   G_Gamma = plaquette-adjacency graph
    mu(G) = sum_{A subset E(G), (V,A) connected & spanning} (-1)^{|A|}   (integer)
=>  c_Gamma = mu(G_Gamma)  (u-independent)   and   m_n = A_n := sum |mu|

So two readings are both reported:
    lambda^sum_n := A_n^{1/n}                  (per-site sum: this is what bounds the
                                                radius r <= lambda_Ur * u of A(alpha))
    lambda^sup_n := max_{|Gamma|=n} |mu|^{1/n}  (per-cluster sup: the literal GU4 limsup,
                                                the one for which |c_G| <= tau holds =>
                                                lambda^sup <= L_t = 53.0006865)

Usage:  py -3.12 gv1_enum.py [maxn] [--json OUT]
"""
import sys, json, time
from fractions import Fraction

DIRS = ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1))
PLANES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
PMAP = {p: i for i, p in enumerate(PLANES)}
OFF = 8
U = 0.625                     # = beta/4 at beta_W = 2.5   (GU4 §2.2, UR6:U1)
EDELTA = 20 * 2.718281828459045235360287471352662497757


# ---------------------------------------------------------------- plaquette graph
def pc(x, a, b):
    if a > b:
        a, b = b, a
    c = (((x[0] + OFF) * 32 + (x[1] + OFF)) * 32 + (x[2] + OFF)) * 32 + (x[3] + OFF)
    return c * 8 + PMAP[(a, b)]


def dec(c):
    pl = c & 7
    s = c >> 3
    d = [0, 0, 0, 0]
    for k in (3, 2, 1, 0):          # x0 is the most significant digit
        d[k] = s % 32 - OFF
        s //= 32
    return tuple(d), PLANES[pl]


def add(x, d):
    return (x[0] + d[0], x[1] + d[1], x[2] + d[2], x[3] + d[3])


def sub(x, d):
    return (x[0] - d[0], x[1] - d[1], x[2] - d[2], x[3] - d[3])


def links_of(p):
    x, (a, b) = dec(p)
    return ((x, a), (x, b), (add(x, DIRS[a]), b), (add(x, DIRS[b]), a))


_NB = {}


def neighbors(c):
    r = _NB.get(c)
    if r is not None:
        return r
    out = set()
    for (y, d) in links_of(c):
        for k in range(4):
            if k == d:
                continue
            out.add(pc(y, d, k))
            out.add(pc(sub(y, DIRS[k]), k, d))
    out.discard(c)
    r = tuple(sorted(out))
    _NB[c] = r
    return r


def animals(maxn, root):
    """levels[i] = set of connected (i+1)-plaquette sets containing root"""
    cur = {(root,)}
    levels = [cur]
    for _ in range(1, maxn):
        nxt = set()
        for S in cur:
            base = set(S)
            last = S[-1]
            for p in S:
                for q in neighbors(p):
                    if q not in base:
                        nxt.add(tuple(sorted(S + (q,))))
        cur = nxt
        levels.append(cur)
    return levels


def adj_of(S):
    n = len(S)
    groups = {}
    for i, p in enumerate(S):
        for lk in links_of(p):
            g = groups.get(lk)
            if g is None:
                groups[lk] = [i]
            else:
                g.append(i)
    adj = [0] * n
    m = 0
    for g in groups.values():
        for a in range(len(g)):
            for b in range(a + 1, len(g)):
                i, j = g[a], g[b]
                if not (adj[i] >> j) & 1:
                    adj[i] |= 1 << j
                    adj[j] |= 1 << i
                    m += 1
    return adj, m


def mu_dp(adj, n):
    """mu(G) = sum_{A subset E, (V,A) connected spanning} (-1)^|A|   (exact integer)"""
    m = 1 << n
    ec = [0] * m
    for T in range(1, m):
        low = T & -T
        i = low.bit_length() - 1
        R = T ^ low
        ec[T] = ec[R] + bin(adj[i] & R).count("1")
    g = [1 if k == 0 else 0 for k in ec]        # (1-1)^e
    f = [0] * m
    for T in range(1, m):
        low = T & -T
        rest = T ^ low
        tot = g[T]
        # iterate over proper subsets of rest (including empty)
        sub_ = rest
        while True:
            if sub_ != rest:
                if g[T ^ (low | sub_)]:
                    tot -= f[low | sub_]
            if sub_ == 0:
                break
            sub_ = (sub_ - 1) & rest
        f[T] = tot
    return f[m - 1]


def mu_fast(S):
    """exact mu with tree / unicyclic fast paths, DP otherwise"""
    n = len(S)
    adj, m = adj_of(S)
    sgn = (-1) ** (n - 1)
    if m == n - 1:
        return sgn, m, "tree"
    if m == n:
        deg = [bin(adj[i]).count("1") for i in range(n)]
        alive = [True] * n
        stack = [i for i in range(n) if deg[i] == 1]
        rem = 0
        while stack:
            i = stack.pop()
            if not alive[i] or deg[i] != 1:
                continue
            alive[i] = False
            rem += 1
            for j in range(n):
                if alive[j] and (adj[i] >> j) & 1:
                    deg[j] -= 1
                    if deg[j] == 1:
                        stack.append(j)
        L = n - rem
        return sgn * (L - 1), m, "unicyclic"
    return mu_dp(adj, n), m, "dp"


# ---------------------------------------------------------------- abstract-graph tests
def g_from_edges(n, E):
    adj = [0] * n
    for (i, j) in E:
        adj[i] |= 1 << j
        adj[j] |= 1 << i
    return adj


def tau_kirchhoff(n, E):
    if n == 1:
        return 1
    L = [[Fraction(0)] * n for _ in range(n)]
    for (i, j) in E:
        L[i][i] += 1
        L[j][j] += 1
        L[i][j] -= 1
        L[j][i] -= 1
    A = [[L[i][j] for j in range(n - 1)] for i in range(n - 1)]
    det = Fraction(1)
    for c in range(n - 1):
        piv = None
        for r in range(c, n - 1):
            if A[r][c] != 0:
                piv = r
                break
        if piv is None:
            return 0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            det = -det
        det *= A[c][c]
        inv = A[c][c]
        for r in range(c + 1, n - 1):
            if A[r][c]:
                f = A[r][c] / inv
                for k in range(c, n - 1):
                    A[r][k] -= f * A[c][k]
    return det


def main():
    maxn = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    out_json = None
    if "--json" in sys.argv:
        out_json = sys.argv[sys.argv.index("--json") + 1]
    t00 = time.time()
    print("=" * 78)
    print("GV1 [A] abstract-graph unit tests on mu (exact integers)")
    print("=" * 78)
    tests = []
    # K_k -> (k-1)! * (-1)^(k-1)   (C02:14 判词 2, exact)
    for k in (2, 3, 4, 5, 6):
        E = [(i, j) for i in range(k) for j in range(i + 1, k)]
        mu = mu_dp(g_from_edges(k, E), k)
        exp = (-1) ** (k - 1)
        for t in range(1, k):
            exp *= t
        tests.append((f"K_{k}", mu, exp))
    # C_4 -> -3 ; C_5 -> 4 ; C_6 -> -5   (C02:14  "C_n -> n-1")
    for k in (4, 5, 6):
        E = [(i, (i + 1) % k) for i in range(k)]
        tests.append((f"C_{k}", mu_dp(g_from_edges(k, E), k), (-1) ** (k - 1) * (k - 1)))
    # path on 4 -> -1 (tree: mu = (-1)^(n-1))
    tests.append(("P_4 (tree)", mu_dp(g_from_edges(4, [(0, 1), (1, 2), (2, 3)]), 4), -1))
    # K_3 --bridge-- K_3 : mu = 2 * 2 * (-1) = -4  (block product, C02 (1.1))
    tests.append(("K_3-br-K_3", mu_dp(g_from_edges(6, [(0, 1), (1, 2), (0, 2),
                                                     (3, 4), (3, 5), (4, 5), (2, 3)]), 6), -4))
    ok = True
    for name, got, exp in tests:
        flag = "OK " if got == exp else "FAIL"
        ok &= (got == exp)
        print(f"  {name:12s} mu = {got:>8d}   expect {exp:>8d}   {flag}")
    print(f"  all abstract tests: {'PASS' if ok else 'FAIL'}")

    # geometric tests: the 6 plaquettes on one link = K_6 ; the 6 faces of a 4D cube
    root = pc((0, 0, 0, 0), 0, 1)
    lk0 = links_of(root)[0]
    K6 = set()
    y, d = lk0
    for k in range(4):
        if k == d:
            continue
        K6.add(pc(y, d, k))
        K6.add(pc(sub(y, DIRS[k]), k, d))
    K6 = tuple(sorted(K6))
    mu6, mE6, tag = mu_fast(K6)
    adj6, _ = adj_of(K6)
    E6 = [(i, j) for i in range(6) for j in range(i + 1, 6) if (adj6[i] >> j) & 1]
    print(f"  K_6 on link {lk0}: n={len(K6)} |E|={mE6} mu={mu6} (expect -120) "
          f"tau={tau_kirchhoff(6, E6)} (expect 1296) "
          f"{'OK' if mu6 == -120 and len(E6) == 15 and tau_kirchhoff(6, E6) == 1296 else 'FAIL'}")
    cube = tuple(sorted(pc((0, 0, 0, 0), a, b) for (a, b) in PLANES))
    muc, mEc, _ = mu_fast(cube)
    adjc, _ = adj_of(cube)
    Ec = [(i, j) for i in range(6) for j in range(i + 1, 6) if (adjc[i] >> j) & 1]
    tauc = tau_kirchhoff(6, Ec)
    print(f"  4D-cube 6 faces:   n=6 |E|={mEc} mu={muc} (expect -64) tau={tauc} (expect 384) "
          f"{'OK' if muc == -64 and mEc == 12 and tauc == 384 else 'FAIL'}")

    print()
    print("=" * 78)
    print(f"GV1 [B] exact enumeration, 4D plaquette graph (Delta=20), root = {dec(root)}")
    print("=" * 78)
    levels = animals(maxn, root)
    print(f"  enumeration to n={maxn} in {time.time()-t00:.1f}s")
    rows = []
    for n in range(1, maxn + 1):
        t1 = time.time()
        S_list = levels[n - 1]
        A = 0
        bmax = 0
        ntree = nuni = ndp = 0
        hist = {}
        nchk = 0
        bad_chk = 0
        for idx, S in enumerate(S_list):
            mu, mE, tag = mu_fast(S)
            if idx % 4999 == 0:                      # fast-path vs DP cross-check
                nchk += 1
                if tag != "dp":
                    adj, _ = adj_of(S)
                    if mu_dp(adj, n) != mu:
                        bad_chk += 1
            a = abs(mu)
            A += a
            if tag == "tree":
                ntree += 1
            elif tag == "unicyclic":
                nuni += 1
            else:
                ndp += 1
            if a > bmax:
                bmax = a
            hist[mu] = hist.get(mu, 0) + 1
        cnt = len(S_list)
        lam_sum = A ** (1.0 / n)
        lam_sup = bmax ** (1.0 / n)
        rows.append(dict(n=n, count=cnt, trees=ntree, unicyclic=nuni, dp=ndp,
                         A=A, max_abs_mu=bmax, lam_sum=lam_sum, lam_sup=lam_sup,
                         avg_abs_mu=A / cnt, chk=nchk, chk_bad=bad_chk,
                         hist={str(k): v for k, v in sorted(hist.items())}))
        print(f"  n={n}: sets={cnt:>9d} trees={ntree:>9d} unicyc={nuni:>8d} dp={ndp:>7d} "
              f"max|mu|={bmax:>7d}  A_n=sum|mu|={A:<12d} avg={A/cnt:8.4f}")
        print(f"        lambda^sum_n = A_n^(1/n) = {lam_sum:10.6f}   "
              f"lambda^sup_n = max|mu|^(1/n) = {lam_sup:8.6f}   "
              f"sum*u = {lam_sum*U:9.6f}   sup*u = {lam_sup*U:8.6f}   [{time.time()-t1:.1f}s]")
        print(f"        fast-path vs DP cross-check: {nchk} sampled, {bad_chk} mismatch")
        levels[n - 1] = None
        del S_list

    print()
    print("=" * 78)
    print("GV1 [C] ratio / growth diagnostics")
    print("=" * 78)
    print("   n   A_n            A_n/A_(n-1)   lambda^sum_n   lambda^sup_n")
    prev = None
    for r in rows:
        rat = (r["A"] / prev["A"]) if prev else float("nan")
        print(f"  {r['n']:2d}  {r['A']:<14d} {rat:12.6f}   {r['lam_sum']:12.6f}   {r['lam_sup']:11.6f}")
        prev = r

    print()
    print("=" * 78)
    print("GV1 [D] sup-version probe beyond n=6: clusters containing a K_6 (exact mu)")
    print("=" * 78)
    base = set(K6)
    ball = set()
    for p in K6:
        ball |= set(neighbors(p))
    ball -= base
    ball = sorted(ball)
    print(f"  K_6 has {len(ball)} plaquette-neighbours outside itself")
    best = []
    m0 = mu_dp(adj_of(K6)[0], 6)
    best.append((6, abs(m0), "K_6"))
    for nadd in (1, 2):
        tb = 0
        arg = None
        from itertools import combinations
        for extra in combinations(ball, nadd):
            S = tuple(sorted(base | set(extra)))
            mu, mE, tag = mu_fast(S)
            if abs(mu) > tb:
                tb = abs(mu)
                arg = (extra, mE, mu)
        nn = 6 + nadd
        if arg:
            best.append((nn, tb, f"K_6+{nadd} (|E|={arg[1]}, mu={arg[2]})"))
            print(f"  n={nn}: best |mu| over K_6+{nadd} = {tb}   -> |mu|^(1/n) = {tb**(1.0/nn):.6f}"
                  f"   [arg: {arg[0]}, |E|={arg[1]}]")
    print("  sup-version floor table:")
    for (nn, m, tag) in best:
        print(f"    n={nn:2d}: |mu|={m:<10d} |mu|^(1/n)={m**(1.0/nn):.6f}  "
              f"(lambda^sup*u = {m**(1.0/nn)*U:.6f})   {tag}")

    if out_json:
        json.dump(dict(u=U, eDelta=EDELTA, rows=rows,
                       sup_probe=[dict(n=a, abs_mu=b, tag=c) for (a, b, c) in best]),
                  open(out_json, "w"), indent=1)
        print(f"\n  json -> {out_json}")
    print(f"\ntotal {time.time()-t00:.1f}s")


if __name__ == "__main__":
    main()
