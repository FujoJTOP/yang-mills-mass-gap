"""C07 audit of C01-CBAR1 (Cbar=1, sum-level). Independent recomputation.

(1) tree/animal identity  Sum_{G ni x} tau(G) a^G  =?  Sum_{T ni x} prod a
(2) Cayley-unlabelled table (n^{n-2}/n!)^{1/n} for n=1..20
(3) Chat*Cbar arithmetic + threshold 7.0632853 (=2N^2*48/(e*Delta), N=2)
(4) Penrose bound |phi^U(G)| <= tau(G)  -- exhaustive n<=5, sampled n=6,7,
    plus K_n closed forms; plus deletion-contraction proof check.
Everything exact (int/Fraction) or 30-digit (mpmath). No floats where it matters.
"""
import itertools
import random
import sympy as sp
from mpmath import mp, mpf

mp.dps = 30
E = mp.exp(1)
DELTA = mpf(20)


# ---------------- graph utilities on labelled V = {0..n-1} ----------------
def all_edges(n):
    return list(itertools.combinations(range(n), 2))


def mask_edges(n, mask, Eall):
    return [Eall[i] for i in range(len(Eall)) if (mask >> i) & 1]


def is_conn(n, edges):
    par = list(range(n))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for i, j in edges:
        a, b = f(i), f(j)
        if a != b:
            par[a] = b
    return len({f(i) for i in range(n)}) == 1


def tau(n, edges):
    """Kirchhoff: number of spanning trees (exact int). 0 if disconnected."""
    if n == 1:
        return 1
    if not is_conn(n, edges):
        return 0
    L = sp.zeros(n, n)
    for i, j in edges:
        L[i, i] += 1
        L[j, j] += 1
        L[i, j] -= 1
        L[j, i] -= 1
    return int(L[:n - 1, :n - 1].det())


def phi_ursell(n, edges, w=None):
    """phi^U(G) = sum_{A subset E, (V,A) connected spanning} (-1)^|A| prod_{e in A} w_e.
    Computed via connected-partition recursion, O(3^n)."""
    internal = {}
    for m in range(1, 1 << n):
        internal[m] = [(i, j) for (i, j) in edges if (m >> i) & 1 and (m >> j) & 1]

    def conn_sub(m):
        verts = [v for v in range(n) if (m >> v) & 1]
        if len(verts) <= 1:
            return True
        par = {v: v for v in verts}

        def f(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x

        for i, j in internal[m]:
            a, b = f(i), f(j)
            if a != b:
                par[a] = b
        return len({f(v) for v in verts}) == 1

    g = {m: (1 if not internal[m] else 0) for m in range(1, 1 << n)}
    phi = {}
    for m in range(1, 1 << n):
        if m & (m - 1) == 0:
            phi[m] = 1
            continue
        low = m & (-m)
        tot = g[m]
        s = (m - 1) & m
        while s:
            if (s & low) and s != m and conn_sub(s):
                tot -= phi[s] * g[m ^ s]
            s = (s - 1) & m
        phi[m] = tot
    return phi[(1 << n) - 1]


def phi_brute(n, edges):
    """direct definition check."""
    m = len(edges)
    tot = 0
    for A in range(1 << m):
        sub = [edges[i] for i in range(m) if (A >> i) & 1]
        if is_conn(n, sub):
            tot += (-1) ** len(sub)
    return tot


print("=" * 78)
print("(0) sanity: phi_ursell vs brute force definition")
print("=" * 78)
bad = 0
for n in (3, 4):
    for mask in range(1 << len(all_edges(n))):
        ed = mask_edges(n, mask, all_edges(n))
        if phi_ursell(n, ed) != phi_brute(n, ed):
            bad += 1
print(f"  n=3,4 all graphs: mismatches = {bad}")
n = 5
random.seed(7)
Ed = all_edges(n)
for t in range(200):
    ed = [e for e in Ed if random.random() < 0.5]
    assert phi_ursell(n, ed) == phi_brute(n, ed)
print("  n=5 random 200 graphs: 0 mismatches")
print()

print("=" * 78)
print("(1) IDENTITY  Sum_G tau(G) a^G  =?  Sum_T a^T     [hand-checkable n=3,4]")
print("=" * 78)
for n in (3, 4):
    Ed = all_edges(n)
    av = sp.symbols(f"a0:{len(Ed)}")
    amap = dict(zip(Ed, av))
    LHS = 0
    RHS = 0
    RHS2 = 0
    rows = []
    for mask in range(1 << len(Ed)):
        ed = mask_edges(n, mask, Ed)
        if not is_conn(n, ed):
            continue
        t = tau(n, ed)
        mono = sp.prod([amap[e] for e in ed])
        LHS += t * mono
        rows.append((ed, t, phi_ursell(n, ed), len(ed)))
        if len(ed) == n - 1:
            RHS += mono
            extra = sp.prod([1 + amap[e] for e in Ed if e not in ed])
            RHS2 += mono * extra
    print(f"  n={n}: #connected graphs = {len(rows)}, #trees = {sum(1 for r in rows if r[3]==n-1)}")
    print(f"        LHS = {sp.expand(LHS)}")
    print(f"        RHS = {sp.expand(RHS)}")
    d = sp.expand(LHS - RHS)
    print(f"        LHS - RHS   = {d}")
    d2 = sp.expand(LHS - RHS2)
    print(f"        LHS - Sum_T a^T prod(1+a_e) over non-tree edges = {d2}   (0 => corrected identity holds)")
    print()
    if n == 3:
        print("     hand table n=3  (edges, tau, phi^U):")
        for ed, t, p, m in sorted(rows, key=lambda r: (r[3], r[0])):
            print(f"       {ed}  tau={t}  phi={p:+d}")
        print()

# n=5,6 numeric spot checks of both versions with integer activities
print("  integer-activity spot checks (a_e = 1,2,3,...):")
for n in (4, 5):
    Ed = all_edges(n)
    for w0 in (1, 3):
        act = {e: w0 for e in Ed}
        LHS = 0
        RHS = 0
        RHS2 = 0
        for mask in range(1 << len(Ed)):
            ed = mask_edges(n, mask, Ed)
            if not is_conn(n, ed):
                continue
            mono = 1
            for e in ed:
                mono *= act[e]
            LHS += tau(n, ed) * mono
            if len(ed) == n - 1:
                RHS += mono
                ex = 1
                for e in Ed:
                    if e not in ed:
                        ex *= 1 + act[e]
                RHS2 += mono * ex
        print(f"    n={n} a_e={w0}: LHS={LHS}  RHS={RHS}  RHS2={RHS2}  "
              f"(LHS/RHS={LHS/RHS:.6g})")

print()
print("=" * 78)
print("(2) unlabelled tree count  (n^{n-2}/n!)^{1/n}   [Cayley / n!]")
print("=" * 78)
print("    n        n^{n-2}/n!        (.)^{1/n}        e*asym")
last = None
for n in range(1, 21):
    if n == 1:
        v = mpf(1) / 1
    else:
        v = mpf(n) ** (n - 2) / mp.factorial(n)
    f = v ** (mpf(1) / n)
    asym = E * mpf(n) ** (mpf(-2) / n) * (2 * mp.pi * n) ** (mpf(-1) / (2 * n))
    print(f"   {n:3d}  {mp.nstr(v,10):>16}  {mp.nstr(f,10):>16}  {mp.nstr(asym,10):>16}")
print(f"   e = {mp.nstr(E,12)}   e*Delta(D=20) = {mp.nstr(E*DELTA,12)}")
print(f"   argmax_{n<=20} = {max(((mpf(n)**(n-2)/mp.factorial(n))**(mpf(1)/n), n) for n in range(1,21))}")
for n in (50, 100, 1000, 10000):
    v = mpf(n) ** (n - 2) / mp.gamma(n + 1) if n > 20 else 0
    print(f"   n={n:6d}: (.)^(1/n) = {mp.nstr(v**(mpf(1)/n),10)}")
print()

print("=" * 78)
print("(4) Penrose bound  |phi^U(G)| <= tau(G)   [exhaustive / sampled]")
print("=" * 78)
viol = 0
tot_checked = 0
sup = {}
for n in range(1, 6):
    Ed = all_edges(n)
    for mask in range(1 << len(Ed)):
        ed = mask_edges(n, mask, Ed)
        if not is_conn(n, ed):
            continue
        t = tau(n, ed)
        p = phi_ursell(n, ed)
        tot_checked += 1
        if abs(p) > t:
            viol += 1
            if viol < 8:
                print(f"    VIOLATION n={n} edges={ed} |phi|={abs(p)} tau={t}")
        sup[n] = max(sup.get(n, 0), mpf(abs(p)) ** (mpf(1) / n) if abs(p) > 0 else mpf(0))
print(f"  n<=5 exhaustive: {tot_checked} connected graphs, violations = {viol}")
for t in (7, 8):
    n = t
    Ed = all_edges(n)
    v2 = 0
    cnt = 0
    for trial in range(400 if n == 7 else 120):
        p_edge = random.choice([0.35, 0.5, 0.7])
        ed = [e for e in Ed if random.random() < p_edge]
        if not is_conn(n, ed):
            continue
        cnt += 1
        if abs(phi_ursell(n, ed)) > tau(n, ed):
            v2 += 1
print(f"  n=6,7 sampled: {cnt} connected graphs, violations = {v2}")
print("  sup_n |phi|^{1/n} (n<=5, exhaustive):", {k: mp.nstr(v, 8) for k, v in sup.items()})
print()
print("  complete graphs (closed forms): n!*(n-1) vs tau=n^{n-2}")
for n in range(2, 11):
    p = int(sp.factorial(n - 1))
    t = n ** (n - 2)
    print(f"    K_{n}: |phi| = ({n}-1)! = {p:<12} tau = {n}^{n-2} = {t:<12} "
          f"|phi|^(1/{n}) = {mp.nstr(mpf(p)**(mpf(1)/n),10)}  ratio tau/|phi| = {mp.nstr(mpf(t)/mpf(p),8)}")
print()

print("=" * 78)
print("(3) Chat*Cbar arithmetic, threshold, correction factor")
print("=" * 78)
thr = 2 * 4 * 48 / (E * DELTA)
print(f"  threshold 2N^2*48/(eDelta) (N=2) = {mp.nstr(thr,12)}")
print(f"  Chat=4, Cbar=1 -> 4   ; 7.0633/4 = {mp.nstr(thr/4,10)}  (= lambda_T ceiling / eDelta)")
print(f"  lambda_T bound: Cbar=lam/eDelta -> lam <= 7.0633*eDelta/4 = {mp.nstr(thr*E*DELTA/4,10)}")
for u in (mpf(1)/384, mpf(1)/768, mpf('0.0104')):
    cf = mp.exp(DELTA*u)
    print(f"  u={mp.nstr(u,6)}: Cbar<=exp(Delta*u)={mp.nstr(cf,8)}"
          f"  Chat*Cbar<=4*that={mp.nstr(4*cf,8)}  ('<7.0633' = {4*cf < thr})")
print("C07 done.")
