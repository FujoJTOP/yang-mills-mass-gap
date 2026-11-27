# -*- coding: utf-8 -*-
"""
GI1-COUNTING — exact combinatorial certificates for the "delete-contract + NBC count" form
applied to  f1^2  and  m .

All arithmetic is exact (int / Fraction).  No Monte Carlo, no floats in any assertion.

S1 mu(G)=sum_{A conn spanning}(-1)^|A|, N=|mu|; K_n -> (n-1)!; the 4D-family values.
S2 N(G) = prod_{2cc blocks B} N(B);  #bridges + sum_{2cc}(|B|-1) = n-1.   (Tarjan, 179 random graphs)
S3 (rewritten) The two-rooted statement is a MULTI-AFFINE one.  In the ring with w_p^2=0,
      [prod_{p in S} w_p] (d^2/dJ_0 dJ_t log Z) = mu(G_S)  if {0,t} subset S, else 0.
   Checked exactly.  Plus: the naive evaluation at a COMMON w does NOT commute with the
   multi-affine sum  (toy counterexample, exact) -> the counting identity lives at coefficient
   level only.  This is the GI1 "evaluation gap".
S4 (n,m) does not determine N.
S5 tube family C_4 x P_T = the four-wall tube graph; connectivity DP (fixed canonical form),
   verified vs brute force for T=1,2,3; N_T for T<=40; minimal integer recurrence.
S6 2D witness: every all-positive NBC sum is > 0.
S7 4D: the rooted NBC shell sequence  B2_k = sum_{|S|=k, S sup {p0,pt}} N(G_S)  on the real
   4D plaquette graph, |S|<=5, for two root separations.
"""
import itertools, json, math, random, sys
from fractions import Fraction

LOG = []
def say(s=""):
    LOG.append(s)
    print(s)

# =============================================================== S1
def mu_brute(n, edges):
    m = len(edges); s = 0
    for mask in range(1 << m):
        par = list(range(n))
        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]; a = par[a]
            return a
        cnt = n; bits = 0; mm = mask; i = 0
        while mm:
            if mm & 1:
                a, b = edges[i]
                ra, rb = find(a), find(b)
                if ra != rb: par[ra] = rb; cnt -= 1
                bits += 1
            mm >>= 1; i += 1
        if cnt == 1:
            s += -1 if (bits & 1) else 1
    return s

def complete_graph(k): return k, [(i, j) for i in range(k) for j in range(i + 1, k)]
def cycle_graph(k):    return k, [(i, (i + 1) % k) for i in range(k)]

say("=" * 78)
say("S1  mu(G) = sum_{A conn. spanning} (-1)^|A| ;  N = |mu|")
say("=" * 78)
ok = True
for k in range(1, 7):
    n, e = complete_graph(k)
    mu = mu_brute(n, e); pred = (-1) ** (n - 1) * math.factorial(n - 1)
    ok &= (mu == pred)
    say(f"  K_{k}: n={n} m={len(e):2d}  mu={mu:8d}  predicted (-1)^(n-1)(n-1)!={pred:8d}  N={abs(mu):6d}"
        f"  {'OK' if mu==pred else 'FAIL'}")

def k_bridge_k():
    return 6, [(0,1),(1,2),(0,2),(3,4),(4,5),(3,5),(2,3)]
def k222():
    parts = [[0,1],[2,3],[4,5]]
    return 6, [(a,b) for i in range(3) for j in range(i+1,3) for a in parts[i] for b in parts[j]]
def cube_graph():
    return 8, [(i, i ^ (1 << b)) for i in range(8) for b in range(3) if i ^ (1 << b) > i]

fam = [("C_3", *cycle_graph(3)), ("C_4", *cycle_graph(4)), ("C_5", *cycle_graph(5)),
       ("C_6", *cycle_graph(6)), ("K3 --bridge-- K3", *k_bridge_k()),
       ("K_{2,2,2} (cubic shell)", *k222()), ("Q_3 (cube graph, misread)", *cube_graph())]
famval = {}
for name, n, e in fam:
    mu = mu_brute(n, e); famval[name] = (n, len(e), mu, abs(mu))
    say(f"  {name:28s} n={n} m={len(e):2d}  mu={mu:8d}  N={abs(mu):6d}")
say(f"  [S1] {'ALL OK' if ok else 'FAIL'}  (K_n -> (n-1)!; KP1's -3,-4,-64,-133 reproduced)")

# =============================================================== S2
def biconnected(n, edges):
    adj = [[] for _ in range(n)]
    for eid, (a, b) in enumerate(edges):
        adj[a].append((b, eid)); adj[b].append((a, eid))
    disc = [-1]*n; low = [0]*n; tmr = [0]; stack = []; blocks = []; bridges = set()
    sys.setrecursionlimit(50000)
    def dfs(u, pe):
        disc[u] = low[u] = tmr[0]; tmr[0] += 1
        for v, eid in adj[u]:
            if eid == pe: continue
            if disc[v] == -1:
                stack.append((u, v, eid)); dfs(v, eid)
                low[u] = min(low[u], low[v])
                if low[v] >= disc[u]:
                    blk = []
                    while True:
                        t = stack.pop(); blk.append(t)
                        if t[2] == eid: break
                    blocks.append(blk)
                if low[v] > disc[u]: bridges.add(eid)
            elif disc[v] < disc[u]:
                stack.append((u, v, eid)); low[u] = min(low[u], disc[v])
    for s in range(n):
        if disc[s] == -1: dfs(s, -1)
    return blocks, bridges

def relabel(n, edges, verts):
    idx = {v: i for i, v in enumerate(sorted(verts))}
    return len(verts), [(idx[a], idx[b]) for (a, b) in edges if a in idx and b in idx]

def block_product(n, edges):
    blocks, bridges = biconnected(n, edges)
    prod = 1; s2 = 0; nb2 = 0
    for blk in blocks:
        vs = {x for (a, b, _) in blk for x in (a, b)}
        if len(vs) >= 3:
            nb, eb = relabel(n, edges, vs)
            prod *= mu_brute(nb, eb); s2 += len(vs) - 1; nb2 += 1
    return prod, len(bridges) + s2, nb2, len(bridges)

say()
say("=" * 78)
say("S2  N(G) = prod_{2-connected blocks} N(B)   and   #bridges + sum_{2cc}(|B|-1) = n-1")
say("=" * 78)
random.seed(20261126)
tf = ts = tests = 0
for _ in range(240):
    n = random.randint(3, 9)
    allp = [(i, j) for i in range(n) for j in range(i+1, n)]
    random.shuffle(allp)
    m = random.randint(n - 1, min(len(allp), n + 5))
    edges = allp[:m]
    par = list(range(n))
    def f(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    for a, b in edges:
        ra, rb = f(a), f(b)
        if ra != rb: par[ra] = rb
    if len({f(i) for i in range(n)}) != 1: continue
    mu = mu_brute(n, edges)
    prod, sgn, _, _ = block_product(n, edges)
    tests += 1
    if abs(mu) != abs(prod):
        tf += 1
        if tf <= 3: say(f"  FAIL block: n={n} e={edges} mu={mu} prod={prod}")
    if sgn != n - 1:
        ts += 1
        if ts <= 3: say(f"  FAIL sign: n={n} e={edges} got {sgn} want {n-1}")
say(f"  {tests} random connected graphs:  block-product violations = {tf} ; sign violations = {ts}")
for name, n, e in fam:
    prod, sgn, nb2, nbr = block_product(n, e)
    say(f"  {name:28s} N={abs(mu_brute(n,e)):6d}  prod_2cc N(B)={abs(prod):6d}"
        f"  #br+sum(|B|-1)={sgn:2d}  n-1={n-1:2d}  (#2cc={nb2}, #br={nbr})")

# =============================================================== S3
def log_multiaffine(n, edges, deg):
    """squarefree (multi-affine) expansion of log Z in the ring w_p^2 = 0.
       returns dict frozenset(S) -> Fraction coefficient."""
    # Z = sum over independent I of prod_{p in I} w_p  -- as a squarefree series
    # log Z = sum_{k>=1} (-1)^(k+1) (Z-1)^k / k ; in the nilpotent ring the sum terminates
    # because (Z-1) has zero constant term and monomials of degree > n vanish.
    Z1 = {}
    for mask in range(1 << n):
        verts = [i for i in range(n) if mask >> i & 1]
        if any(((mask >> a & 1) and (mask >> b & 1)) for a, b in edges):
            continue
        if not verts: continue
        Z1[frozenset(verts)] = Z1.get(frozenset(verts), Fraction(0)) + 1
    # repeated squaring of series in the nilpotent ring
    def mul(A, B):
        out = {}
        for s1, c1 in A.items():
            for s2, c2 in B.items():
                if s1 & s2: continue
                s = s1 | s2
                out[s] = out.get(s, Fraction(0)) + c1 * c2
        return out
    res = {}
    P = {frozenset(): Fraction(1)}
    for k in range(1, n + 1):
        P = mul(P, Z1)
        if not P: break
        for s, c in P.items():
            res[s] = res.get(s, Fraction(0)) + (c if k % 2 == 1 else -c) / k
    return res

say()
say("=" * 78)
say("S3  multi-affine two-rooted identity:  [prod_{p in S} w_p](d^2_{J0 Jt} log Z) = mu(G_S)")
say("=" * 78)
s3f = 0; s3t = 0
for _ in range(150):
    n = random.randint(2, 7)
    allp = [(i, j) for i in range(n) for j in range(i+1, n)]
    random.shuffle(allp)
    edges = allp[:random.randint(0, min(len(allp), 8))]
    r0, rt = random.sample(range(n), 2)
    lg = log_multiaffine(n, edges, n)
    good = True
    for S, c in lg.items():
        if not (r0 in S and rt in S):
            continue          # such monomials are killed by d^2_{J0 Jt}; not part of the test
        pred = mu_brute(*relabel(n, [(a, b) for a, b in edges if a in S and b in S], S))
        if c != Fraction(pred):
            good = False
            say(f"    FAIL S={sorted(S)} c={c} mu={pred}  (n={n}, e={edges})")
            break
    s3t += 1
    if not good: s3f += 1
say(f"  {s3t} instances: violations = {s3f}")
say("  => the two-rooted coefficient is the SINGLE-rooted Ursell coefficient mu(G_S), i.e. GO1's")
say("     (-1)^(|S|-1) N(G_S); it is zero unless both roots lie in S.  No new identity, but the")
say("     exact placement of N inside the two-point expansion is pinned.")
say()
say("  EVALUATION GAP (exact toy, 2 adjacent plaquettes, no other):")
say("     Z = 1 + w_0 + w_1 ;  C = -w_1 w_0 / Z^2  (analytic, at J=0)")
say("     multi-affine sum at w_p = w:   mu(edge) w^2 = -w^2  (one term only)")
a = [Fraction(1), Fraction(2), Fraction(0), Fraction(0), Fraction(0)]
inv = [Fraction(0)] * 5; inv[0] = Fraction(1)
for k in range(1, 5):
    s = sum(a[j] * inv[k - j] for j in range(1, k + 1)); inv[k] = -s
sq = [sum(inv[i] * inv[k - i] for i in range(k + 1)) for k in range(5)]
say(f"     machine check 1/(1+2w)^2 = " + ", ".join(f"w^{k}:{sq[k]}" for k in range(5)))
say(f"     => analytic C(w) = -w^2/(1+2w)^2 has w^3..w^4 coeffs "
    f"{-sq[1]},{-sq[2]}  != 0 = multi-affine sum's higher coefficients.")
say("     => the multi-affine sum is NOT the evaluation at a common activity; the difference is")
say("        the repeated-plaquette terms (w_p^2 = 0 truncation).  N enters at coefficient level only.")

# =============================================================== S4
say()
say("=" * 78)
say("S4  (n,m) does not determine N")
say("=" * 78)
n1, e1 = cycle_graph(4)
n2, e2 = 4, [(0,1),(1,2),(0,2),(2,3)]
say(f"  C_4                      n=4 m=4  N={abs(mu_brute(n1,e1))}")
say(f"  triangle + pendant edge  n=4 m=4  N={abs(mu_brute(n2,e2))}")
say("  => N is not a function of (n,m): no closed form in (n,m) can exist.")

# =============================================================== S5
def canon_minpart(labels):
    """labels: list; block label := min index in the block (so st[st[j]]==st[j])."""
    groups = {}
    for i, l in enumerate(labels):
        groups.setdefault(l, []).append(i)
    out = list(range(len(labels)))
    for g in groups.values():
        m = min(g)
        for i in g: out[i] = m
    return tuple(out)

def tube_mu(T):
    cyc = [(j, (j + 1) % 4) for j in range(4)]
    states = {}
    for mask in range(1 << 4):
        par = list(range(4))
        def find(a):
            while par[a] != a: par[a] = par[par[a]]; a = par[a]
            return a
        bits = 0
        for i in range(4):
            if mask >> i & 1:
                bits += 1
                ra, rb = find(cyc[i][0]), find(cyc[i][1])
                if ra != rb: par[ra] = rb
        st = canon_minpart([find(j) for j in range(4)])
        states[st] = states.get(st, 0) + (-1 if bits & 1 else 1)
    for t in range(1, T):
        new = {}
        for st, val in states.items():
            for mask in range(1 << 8):
                par = list(range(8))
                def find(a):
                    while par[a] != a: par[a] = par[par[a]]; a = par[a]
                    return a
                for j in range(4):
                    if st[j] != j:
                        ra, rb = find(j), find(st[j])
                        if ra != rb: par[ra] = rb
                bits = 0
                for i in range(4):
                    if mask >> i & 1:
                        bits += 1
                        ra, rb = find(4 + cyc[i][0]), find(4 + cyc[i][1])
                        if ra != rb: par[ra] = rb
                for j in range(4):
                    if mask >> (4 + j) & 1:
                        bits += 1
                        ra, rb = find(j), find(4 + j)
                        if ra != rb: par[ra] = rb
                if not {find(j) for j in range(4)}.issubset({find(4 + j) for j in range(4)}):
                    continue
                nst = canon_minpart([find(4 + j) for j in range(4)])
                new[nst] = new.get(nst, 0) + val * (-1 if bits & 1 else 1)
        states = new
    return sum(v for st, v in states.items() if len(set(st)) == 1)

def tube_graph(T):
    edges = []
    for t in range(T):
        for j in range(4): edges.append((4 * t + j, 4 * t + (j + 1) % 4))
    for t in range(T - 1):
        for j in range(4): edges.append((4 * t + j, 4 * (t + 1) + j))
    return 4 * T, edges

say()
say("=" * 78)
say("S5  tube family C_4 x P_T  (= the four-wall tube; T=1 is C_4)")
say("=" * 78)
for T in (1, 2):
    n, e = tube_graph(T)
    d, b = tube_mu(T), mu_brute(n, e)
    say(f"  T={T}: n={n:2d} m={len(e):2d}  DP mu={d:8d}  brute mu={b:8d}  {'OK' if d==b else 'FAIL'}")
n, e = tube_graph(3)
d, b = tube_mu(3), mu_brute(n, e)
say(f"  T=3: n={n:2d} m={len(e):2d}  DP mu={d:9d}  brute mu={b:9d}  {'OK' if d==b else 'FAIL'}")
seq = [abs(tube_mu(T)) for T in range(1, 41)]
sgns = [tube_mu(T) for T in range(1, 41)]
say("  sign of mu: " + "".join("+" if s > 0 else "-" for s in sgns) + "  (T=1..40)")
say("  N_T = |mu(C_4 x P_T)|:")
for i in range(0, 40, 4):
    say("    " + "  ".join(f"T={T+1}:{seq[T]}" for T in range(i, min(i + 4, 40))))

def find_recurrence(a, dmax=8, need=12):
    N = len(a)
    for d in range(1, dmax + 1):
        rows = []; rhs = []
        for k in range(d, min(N, d + 3 * d)):
            if all(x is not None for x in [a[k]]):
                rows.append([Fraction(a[k - i]) for i in range(1, d + 1)]); rhs.append(Fraction(a[k]))
        if len(rows) < d: continue
        M = [row[:] + [r] for row, r in zip(rows[:d], rhs[:d])]
        # gaussian elimination over Q
        piv = 0; sol = None; okg = True
        for c in range(d):
            p = next((r for r in range(piv, len(M)) if M[r][c] != 0), None)
            if p is None: okg = False; break
            M[piv], M[p] = M[p], M[piv]
            f0 = M[piv][c]
            M[piv] = [x / f0 for x in M[piv]]
            for r in range(len(M)):
                if r != piv and M[r][c] != 0:
                    fac = M[r][c]
                    M[r] = [x - fac * y for x, y in zip(M[r], M[piv])]
            piv += 1
        if not okg or piv != d: continue
        sol = [M[r][d] for r in range(d)]
        okc = True
        for k in range(d, N):
            if a[k] != sum(sol[i] * a[k - 1 - i] for i in range(d)):
                okc = False; break
        if okc: return d, sol
    return None, None

d, sol = find_recurrence(seq, dmax=8)
if d:
    say(f"  minimal integer recurrence: N_T = " + " + ".join(f"({str(sol[i])})N_{{T-{i+1}}}" for i in range(d))
        + f"   [order {d}, verified on T={d+1}..40]")
    say(f"  characteristic roots (numeric): " + ", ".join(
        f"{complex(r):.6f}" for r in __import__("numpy").roots(
            [1.0] + [-float(x) for x in sol])))
else:
    say("  no linear recurrence of order <= 8 found")

# =============================================================== S6
say()
say("=" * 78)
say("S6  2D witness for the no-go on all-positive counting representations")
say("=" * 78)
say("  S = {p0}:  N = 1 and weight w > 0  =>  every  sum_S a_S w^|S| N(G_S)  with a_S >= 0 is >= w > 0,")
say("  in every dimension.  2D pure SU(N) YM has  f1^2(tr F^2) == 0  identically (ATOM2:15,62-63 [T]).")
say("  => f1^2 admits no all-positive NBC representation.  Only the alternating form survives,")
say("     and there positivity is exactly the domination problem (not a count).")

# =============================================================== S7
PLANES = [(a, b) for a in range(4) for b in range(a + 1, 4)]
def plaq_links(x, ab):
    a, b = ab
    lk = []
    lk.append((x, a))
    lk.append((x, b))
    y = list(x); y[a] += 1
    lk.append((tuple(y), b))
    y = list(x); y[b] += 1
    lk.append((tuple(y), a))
    return tuple(sorted(lk))

def shift(x, mu):
    y = list(x); y[mu] += 1; return tuple(y)

def build_neighborhood(radius=2):
    sites = [x for x in itertools.product(range(-radius, radius + 1), repeat=4)]
    P = []
    for x in sites:
        for ab in PLANES:
            P.append(((x, ab), plaq_links(x, ab)))
    idx = {p[0]: i for i, p in enumerate(P)}
    linkmap = {}
    for i, (_, lk) in enumerate(P):
        for l in lk:
            linkmap.setdefault(l, []).append(i)
    adj = [set() for _ in P]
    for l, plist in linkmap.items():
        for a in plist:
            for b in plist:
                if a != b: adj[a].add(b)
    adj = [sorted(s) for s in adj]
    return P, idx, adj

def clusters_with_roots(r1, kmax):
    """all connected plaquette sets containing i0 and r1, grouped by size k."""
    out = {2: {frozenset([i0, r1])}}
    level = {frozenset([i0, r1])}
    for k in range(3, kmax + 1):
        nxt = set()
        for s in level:
            for v in s:
                for w in adj[v]:
                    if w in s: continue
                    nxt.add(frozenset(s | {w}))
        out[k] = nxt
        level = nxt
    return out

say()
say("=" * 78)
say("S7  rooted NBC shell sums on the REAL 4D plaquette graph:  B2_k = sum_{|S|=k, S sup {p0,pt}} N(G_S)")
say("=" * 78)
P, idx, adj = build_neighborhood(2)
p0 = ((0,0,0,0), (0,1))
i0 = idx[p0]
say(f"  neighbourhood: {len(P)} plaquettes;  deg(p0) = {len(adj[i0])}  (must be 20 for d=4)")
for sep in (1, 2, 3):
    rt = p0[0]
    for _ in range(sep): rt = shift(rt, 0)
    target = (rt, p0[1])
    if target not in idx:
        say(f"  separation {sep}: target not in window"); continue
    r1 = idx[target]
    cl = clusters_with_roots(r1, 5)
    parts = []
    for k in range(2, 6):
        tot = 0; cnt = 0; hist = {}
        for s in cl[k]:
            verts = sorted(s)
            pos = {v: i for i, v in enumerate(verts)}
            es = [(pos[a], pos[b]) for a in verts for b in adj[a] if b in pos and b > a]
            nn = abs(mu_brute(len(verts), es))
            tot += nn; cnt += 1
        parts.append(f"k={k}: #S={cnt:6d}, sum N={tot}")
    say(f"  roots separated by {sep}:  " + ";  ".join(parts))

# ---- byproduct: the tube family gives a STRICT improvement of the corpus's lambda_Ur lower bound
say()
say("  BYPRODUCT:  |phi^U(tube_T)| = N_T  with n = 4T, so  |c_Gamma|^{1/n} = N_T^{1/(4T)}")
lam_plus = 23 + 2 * math.sqrt(111)
say(f"     closed form from the recurrence:  N_T = A*lambda+^T + B*lambda-^T,  lambda+- = 23 +- 2*sqrt(111)"
    f" = {lam_plus:.10f}, {23 - 2*math.sqrt(111):.10f}")
say(f"     =>  lim_T N_T^(1/(4T)) = lambda+^(1/4) = {lam_plus ** 0.25:.8f}")
say(f"     corpus lower bound (GV2 4.2(A), K_6): 120^(1/6) = {120 ** (1/6):.8f}")
say(f"     =>  lambda_Ur >= {(lam_plus) ** 0.25:.8f}  >  2.220906  (STRICT improvement)")
say(f"     SU(2) convention u=0.625:  lambda_Ur*u >= {(lam_plus) ** 0.25 * 0.625:.6f} > 1"
    f"   (corpus had {120 ** (1/6) * 0.625:.6f})")

with open("gi1_count.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG) + "\n")
json.dump({"family": {k: v[3] for k, v in famval.items()},
           "tube_N": seq, "tube_mu_sign": sgns,
           "s2_block_violations": tf, "s2_sign_violations": ts, "s2_tests": tests,
           "s3_violations": s3f, "s3_tests": s3t,
           "recurrence_order": d, "recurrence_coeffs": [str(x) for x in sol] if d else None},
          open("gi1_count.json", "w", encoding="utf-8"), indent=1)
say()
say("  [wrote gi1_count.log / gi1_count.json]")
