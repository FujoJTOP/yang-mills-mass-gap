# C01 — is  Cbar = 1?  Ursell coefficient exact values, tree vs animal entropy, one-bill test.
# py -3.12 ; stdlib only (itertools, math, random, fractions-free integer Bareiss).
import math, random, itertools, sys

# ---------- exact Ursell (Mayer) coefficient ----------
def phi_ursell(n, edges):
    """phi^U(G) = sum_{A subset E, (V,A) connected spanning} (-1)^{|A|}."""
    m = len(edges)
    tot = 0
    for mask in range(1 << m):
        par = list(range(n))
        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]; a = par[a]
            return a
        k = 0
        for i in range(m):
            if mask >> i & 1:
                a, b = find(edges[i][0]), find(edges[i][1])
                if a != b:
                    par[a] = b; k += 1
        if k == n - 1:                       # connected spanning
            tot += -1 if (bin(mask).count("1") & 1) else 1
    return tot

def tau_span(n, edges):
    """integer spanning-tree count, Bareiss on reduced Laplacian."""
    if n == 1:
        return 1
    L = [[0] * n for _ in range(n)]
    for a, b in edges:
        L[a][a] += 1; L[b][b] += 1; L[a][b] -= 1; L[b][a] -= 1
    M = [[L[i][j] for j in range(1, n)] for i in range(1, n)]
    s = n - 1; sign = 1; prev = 1
    for i in range(s - 1):
        if M[i][i] == 0:
            for r in range(i + 1, s):
                if M[r][i]:
                    M[i], M[r] = M[r], M[i]; sign = -sign; break
            else:
                return 0
        for r in range(i + 1, s):
            for c in range(i + 1, s):
                M[r][c] = (M[r][c] * M[i][i] - M[r][i] * M[i][c]) // prev
        prev = M[i][i]
        for r in range(i + 1, s):
            M[r][i] = 0
    return sign * M[s - 1][s - 1]

def conn(n, edges):
    par = list(range(n))
    def find(a):
        while par[a] != a:
            par[a] = par[par[a]]; a = par[a]
        return a
    for a, b in edges:
        x, y = find(a), find(b)
        if x != y: par[x] = y
    return len({find(i) for i in range(n)}) == 1

# ---------- 1. per-Gamma witness: cliques ----------
print("== S1  phi^U on cliques K_n  (witness for the pointwise sup) ==")
print(" n   |phi^U(K_n)|=(n-1)!   ^(1/n)    tau(K_n)=n^(n-2)")
Cs = {}
for n in range(2, 11):
    ph = math.factorial(n - 1)
    Cs[n] = ph ** (1.0 / n)
    print(f"{n:3d} {ph:20d} {Cs[n]:10.6f} {n**(n-2):20d}")
print(" K_6 = the 6 plaquettes at one link (all pairwise share that link):")
print(f" |phi^U| = {math.factorial(5)}  ->  Cbar >= {math.factorial(5)**(1/6):.6f} > 1.7658  => pointwise Cbar=1 is FALSE")

# ---------- 2. Penrose/tree-graph inequality:  sum |phi| over connected sets <= tree sum ----------
def rooted_tree_weight(adj, x, n, w):
    """sum over trees rooted at x with vertex set of size n, weight prod w, in graph adj (list of sets)."""
    total = 0
    def rec(used, frontier):
        nonlocal total
        if len(used) == n:
            total += math.prod(w[v] for v in used); return
        for v in list(frontier):
            if v in used: continue
            nb = (frontier | adj[v]) - used - {v}
            rec(used | {v}, nb)
    rec({x}, set(adj[x]))
    return total

def all_connected_sets(adj, x, nmax, cap=400000):
    """all connected vertex sets containing x, size<=nmax (as frozensets)."""
    out = []
    stack = [frozenset([x])]
    seen = {frozenset([x])}
    while stack:
        S = stack.pop(); out.append(S)
        if len(S) >= nmax: continue
        nb = set()
        for v in S: nb |= adj[v]
        for v in nb - S:
            T = S | {v}
            if T not in seen:
                seen.add(T); stack.append(T)
                if len(seen) > cap: raise RuntimeError("cap")
    return out

print("\n== S2  tree-graph inequality  LHS=sum_{Gamma connected,ni x}|phi^U|  <=  RHS=#rooted trees ni x ==")
random.seed(7)
bad = 0
for n in range(2, 7):
    trials = 1 if n <= 4 else 12
    for t in range(trials):
        for wmode in ("one", "rand"):
            if n <= 4:
                es = [(i, j) for i in range(n) for j in range(i + 1, n)]
                # exhaustive over all graphs for n<=4
                graphs = [g for g in range(1 << len(es))]
            else:
                graphs = [None] * trials
            gl = graphs if n <= 4 else [graphs[t]]
            for g in gl:
                if n <= 4:
                    edges = [es[i] for i in range(len(es)) if g >> i & 1]
                else:
                    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if random.random() < 0.55]
                if not conn(n, edges): continue
                adj = [set() for _ in range(n)]
                for a, b in edges: adj[a].add(b); adj[b].add(a)
                if max(len(a) for a in adj) > 12: continue
                w = [1.0] * n if wmode == "one" else [random.uniform(0.05, 1.0) for _ in range(n)]
                L = 0.0
                for S in all_connected_sets(adj, 0, n):
                    vl = sorted(S); pos = {v: i for i, v in enumerate(vl)}
                    es = [(pos[u], pos[v]) for u in vl for v in adj[u] if v in pos and u < v]
                    L += abs(phi_ursell(len(vl), es)) * math.prod(w[v] for v in S)
                # LHS restricted to vertex sets containing 0, any size <= n  -> compare with tree sum
                R = sum(rooted_tree_weight(adj, 0, k, w) for k in range(1, n + 1))
                if L > R + 1e-9:
                    bad += 1
                    if bad <= 5: print(f"  COUNTEREXAMPLE n={n} edges={edges} L={L} R={R}")
print(f"  counterexamples: {bad}   (0 => inequality holds on all sampled graphs, n<=6)")

# ---------- 3. Cayley with and without the 1/n! (unlabelling) ----------
print("\n== S3  tree entropy: labelled vs unlabelled (positioned) count ==")
print(" n   (n^(n-2))^(1/n)   (n^(n-2)/n!)^(1/n)   -> limit e")
for n in (2, 3, 5, 10, 20, 50, 100, 200, 500, 1000):
    a = math.exp((n - 2) * math.log(n) / n) if n > 1 else 1.0
    b = math.exp(((n - 2) * math.log(n) - math.lgamma(n + 1)) / n)
    print(f"{n:5d} {a:14.6f} {b:16.6f}")
print(" e =", round(math.e, 6), " ; 1/n! division turns n^(n-2) divergence into entropy e per vertex")

# ---------- 4. 4D plaquette adjacency, exact phi on small connected sets ----------
def lattice4(L):
    sites = list(itertools.product(range(L), repeat=4))
    def nxt(s, mu):
        t = list(s); t[mu] = (t[mu] + 1) % L; return tuple(t)
    plaq = []
    for s in sites:
        for mu in range(4):
            for nu in range(mu + 1, 4):
                links = frozenset({(s, mu), (s, nu), (nxt(s, mu), nu), (nxt(s, nu), mu)})
                plaq.append((s, mu, nu, links))
    return plaq

L = 6
plaq = lattice4(L)
idx = {}
for i, p in enumerate(plaq):
    for lk in p[3]: idx.setdefault(lk, []).append(i)
nbr = [set() for _ in plaq]
for lk, lst in idx.items():
    for a in lst:
        for b in lst:
            if a != b: nbr[a].add(b)
print(f"\n== S4  4D torus L={L}: 6L^4={6*L**4} plaquettes, deg = {len(nbr[0])} (Delta=20 expected) ==")
x = 0
sets = all_connected_sets(nbr, x, 5)          # sizes 1..5
bysize = {}
for S in sets:
    bysize.setdefault(len(S), []).append(S)
print(" n   #sets(ni x)   sum|phi^U|      max|phi^U|   (sum)^(1/n)   (max)^(1/n)   viol |phi|>tau")
Cbar_emp = 0.0
viol = 0
for n in sorted(bysize):
    Ss = bysize[n]
    vals = []
    for S in Ss:
        vl = sorted(S); pos = {v: i for i, v in enumerate(vl)}
        es = [(pos[u], pos[v]) for i, u in enumerate(vl) for v in nbr[u] if v in pos and u < v]
        ph = abs(phi_ursell(n, es))
        ta = tau_span(n, es)
        if ph > ta + 1e-9: viol += 1
        vals.append(ph)
    s = sum(vals); mx = max(vals)
    Cbar_emp = max(Cbar_emp, mx ** (1.0 / n))
    print(f"{n:3d} {len(Ss):12d} {s:14.0f} {mx:12.0f} {s**(1.0/n):14.6f} {mx**(1.0/n):12.6f} {viol:8d}")
print(f" empirical per-Gamma Cbar (n<=5, exact) = {Cbar_emp:.6f}")
print(f" |phi^U(S)| <= tau(G[S]) violations: {viol}  (0 => Penrose pointwise bound holds)")

# ---------- 5. the 6-plaquette clique at a link, exact ----------
lk = next(iter(plaq[x][3]))
clique = idx[lk]
print(f"\n== S5  plaquettes through one link: {len(clique)}  -> adjacency among them = K_{len(clique)} ==")
vl = clique; pos = {v: i for i, v in enumerate(vl)}
es = [(pos[u], pos[v]) for u in vl for v in nbr[u] if v in pos and u < v]
print(f" graph: n={len(vl)}, |E|={len(es)} ; phi^U = {phi_ursell(len(vl), es)} ; tau = {tau_span(len(vl), es)}")

# ---------- 6. window arithmetic ----------
print("\n== S6  window:  Cbar=1 (one bill) vs Cbar=20 (a01 product bound) ==")
N = 2; Dl = 20; lam = math.e * Dl
thr = 2 * N * N * 48 / lam
print(f" threshold  Chat*Cbar < 2N^2*48/(eDelta) = {thr:.8f}")
for Chat in (2, 4):
    for Cb, tag in ((1.0, "one-bill Cbar=1"), (20.0, "a01 product"), (math.factorial(5)**(1/6), "exact per-Gamma sup (K_6)")):
        ok = Chat * Cb < thr
        b1 = 2 * N * N / (Chat * lam * Cb)
        print(f"  Chat={Chat} Cbar={Cb:8.4f} [{tag:28s}] Chat*Cbar={Chat*Cb:9.4f} {'OK ' if ok else 'NG '} beta1={b1:.6e} = {(b1*48):.4f}*(1/48)")
print(f" e*Delta = {lam:.6f}")
