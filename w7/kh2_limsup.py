"""kh2_limsup.py -- KH2: sup vs limsup vs average for K_hat.

Q1  What is limsup_n |phi^U_n|^{1/n}?  (sup was 2.9072 at n=11 and still rising)
Q2  Does the KP criterion need sup (for all gamma) -> yes, that is the criterion.
Q3  Is there an averaged (geometric-mean) version?  C_hat already *is* a geometric mean.
Q4  Can limsup or averaging push K_hat below 7.0632853?

Numbers produced here
  [A] fast exact phi^U (independent-set / cluster DP) validated against the O(3^n) DP
  [B] the multi-beam family: k link-bundles at one site, pairwise sharing one plaquette
        k=1 n=6  -> 120        (K_6)
        k=2 n=11 -> 125400     (cross-check: KP4 / kp4_chain.py)
        k=3 n=15 -> ?          (new: the limsup lower-bound point)
  [C] block-product reference rates (what a *tree of blocks* could give)
  [D] the averaged readings: m_n^{1/n} (arithmetic sum over clusters) vs sup
  [E] verdict numbers: K_hat = 4*limsup-layer, threshold 7.0632853
"""
import itertools, time, math

D = 4
T0 = time.time()
def log(*a): print(*a, flush=True)
def unit(i):
    v = [0]*D; v[i] = 1; return tuple(v)
def add(u, v): return tuple(u[k]+v[k] for k in range(D))
DIRS = [(a, b) for a in range(D) for b in range(a+1, D)]

THR = 96/(5*math.e)          # 7.06328527049169

# ---------- exact phi^U, two independent implementations ----------
def phi_3n(adj, n):
    """O(3^n) reference DP (kp4_higher.phi_int)."""
    size = 1 << n; H = bytearray(size); H[0] = 1
    for S in range(1, size):
        low = S & -S; i = low.bit_length()-1; rest = S ^ low
        H[S] = 1 if (adj[i] & rest) == 0 and H[rest] else 0
    F = [0]*size
    for S in range(1, size):
        low = S & -S; rest = S ^ low; tot = H[S]; t = rest
        while True:
            blk = low | t
            if blk != S: tot -= F[blk]*H[S ^ blk]
            if t == 0: break
            t = (t-1) & rest
        F[S] = tot
    return F[size-1]

def phi_fast(adj, n):
    """F(S) = H(S) - sum_{empty!=J subset S, min(S) not in J, J independent} F(S\\J).
    Only independent J survive, so the cost is 2^n * #indep, not 3^n."""
    size = 1 << n; H = bytearray(size); H[0] = 1
    for S in range(1, size):
        low = S & -S; i = low.bit_length()-1; rest = S ^ low
        H[S] = 1 if (adj[i] & rest) == 0 and H[rest] else 0
    indep = [J for J in range(1, size) if H[J]]
    log(f"    [phi_fast] n={n} 2^n={size} #indep sets={len(indep)}")
    full = size-1; F = [0]*size
    for S in range(1, size):
        low = S & -S; tot = H[S]
        for J in indep:
            if J & low: continue
            if J & ~S: continue
            tot -= F[S ^ J]
        F[S] = tot
    return F[full]

# ---------- small lattice pool ----------
def build_box(RMAX):
    sites = [x for x in itertools.product(range(-RMAX, RMAX+1), repeat=D)]
    plaq = [(x, a, b) for x in sites for (a, b) in DIRS]
    idx = {p: i for i, p in enumerate(plaq)}
    linkmap = {}
    for i, (x, a, b) in enumerate(plaq):
        for L in ((x, a), (add(x, unit(a)), b), (x, b), (add(x, unit(b)), a)):
            linkmap.setdefault(L, []).append(i)
    nb = [set() for _ in plaq]
    for L, ps in linkmap.items():
        for i in range(len(ps)):
            for j in range(i+1, len(ps)):
                nb[ps[i]].add(ps[j]); nb[ps[j]].add(ps[i])
    return plaq, idx, linkmap, nb

def rows_of(vs, nb):
    n = len(vs)
    rows = [0]*n
    for i, v in enumerate(vs):
        r = 0
        for j, w in enumerate(vs):
            if j != i and w in nb[v]: r |= 1 << j
        rows[i] = r
    return rows

def edge_count(n, rows):
    return sum(bin(r).count("1") for r in rows)//2

def prodR_root(vs, linkmap, n, activated=False):
    """(prod_l R(l))^{1/n}, R(kappa)=2^{kappa/2}  (Z12-FIXC:206).
    sum_l kappa_l = 4n ALWAYS (each plaquette owns 4 links) => root = 4 exactly.
    activated=True: links with kappa=1 get R=1 instead of 2^{1/2}."""
    kap = {}
    for v in vs:
        for L, ps in linkmap.items():
            if v in ps:
                kap[L] = kap.get(L, 0) + 1
    tot = sum((0.0 if (activated and c == 1) else c/2.0) for c in kap.values())
    hist = {}
    for c in kap.values():
        hist[c] = hist.get(c, 0) + 1
    return 2.0 ** (tot/n), sum(kap.values()), hist

log("="*78)
log("[0] threshold and the two readings")
log(f"    threshold T = 96/(5e) = {THR:.13f}")
log(f"    K_hat < T  <=>  sup_n|phi^U|^{{1/n}} < T/4 = {THR/4:.13f}   (C_hat=4, (prodR)^{{1/n}}=4)")

# ---------- [A] validation ----------
log("\n[A] validation of phi_fast against phi_3n on named graphs")
def from_edges(n, E):
    rows = [0]*n
    for (i, j) in E:
        rows[i] |= 1 << j; rows[j] |= 1 << i
    return rows
def K(k):
    return from_edges(k, [(i, j) for i in range(k) for j in range(i+1, k)])
for k in (3, 4, 5, 6):
    r = K(k)
    a, b = phi_fast(r, k), phi_3n(r, k)
    log(f"    K_{k}: phi_fast={a} phi_3n={b}  {'OK' if a == b else 'MISMATCH'}  |phi|^{{1/{k}}}={abs(a)**(1/k):.6f}")

# ---------- [B] multi-beam family ----------
log("\n[B] k link-bundles at one site (pairwise sharing exactly one plaquette)")
plaq, idx, linkmap, nb = build_box(2)
log(f"    box R=2: {len(plaq)} plaquettes, links carry {sorted({len(v) for v in linkmap.values()})}")
base = (0, 0, 0, 0)
res = {}
for k in (1, 2, 3):
    links = [(base, a) for a in range(k)]
    vs = set()
    for L in links:
        vs |= set(linkmap[L])
    vs = sorted(vs)
    n = len(vs)
    rows = rows_of(vs, nb)
    t1 = time.time()
    ph = phi_fast(rows, n)
    ph3 = phi_3n(rows, n) if n <= 12 else None
    rk, sumkap, hist = prodR_root(vs, linkmap, n, False)
    rk_act, _, _ = prodR_root(vs, linkmap, n, True)
    res[k] = (n, abs(ph), rk)
    log(f"    k={k}: links={links} n={n} |E|={edge_count(n,rows)} phi^U={ph} |phi|={abs(ph)} "
        f"|phi|^{{1/n}}={abs(ph)**(1/n):.6f} (prodR)^{{1/n}}={rk:.6f} "
        f"K_hat>={rk*abs(ph)**(1/n):.6f} [{time.time()-t1:.1f}s]"
        + (f"  3n-DP={ph3} {'OK' if ph3 == ph else 'MISMATCH'}" if ph3 is not None else ""))
    log(f"          link kappa histogram {hist}  sum(kappa)={sumkap} (=4n={4*n})  "
        f"(prodR)^{{1/n}} activated={rk_act:.6f}")
    # shared plaquettes
    sp = []
    for a in range(len(links)):
        for b in range(a+1, len(links)):
            sp.append(len(set(linkmap[links[a]]) & set(linkmap[links[b]])))
    log(f"          shared plaquette counts per pair: {sp}   (block product would be 120^{k})")

log("\n[B2] the KP4 table (sup reading) vs the new n=15 point")
kp4 = {5: 24, 6: 120, 7: 240, 8: 924, 9: 5712, 10: 25080, 11: 125400}
for n, m in kp4.items():
    log(f"    n={n:>2}: max|phi^U| >= {m:>7}  ^{1}/{n} = {m**(1/n):.6f}")
log(f"    n={res[3][0]:>2}: max|phi^U| >= {res[3][1]:>7}  ^{1}/{res[3][0]} = {res[3][1]**(1/res[3][0]):.6f}   <- NEW")

log("\n[B3] limsup witness: r replicas of the k=3 cluster, joined by single bridge plaquettes")
V, n3 = res[3][1], res[3][0]
log(f"    replica value V={V} (n3={n3});  bridge factorisation is the block rule phi^U=prod(blocks)")
log(f"    (checked at n=13 in kp4_higher.py: 120*120=14400 for two K_6 + one bridge plaquette [E])")
for r in (1, 2, 3, 5, 10, 100, 1000):
    n = n3*r + (r-1)
    rate = V ** (r/n)
    log(f"    r={r:>4}: n={n:>6}  (V^r)^(1/n) = {rate:.6f}")
log(f"    r->inf limit: V^(1/{n3+1}) = {V**(1.0/(n3+1)):.6f}   (each bridge adds exactly 1 plaquette)")
log(f"    => limsup_n >= {V**(1.0/(n3+1)):.6f}   [C, uses cut-vertex factorisation]")
log(f"    => sup_n    >= {V**(1.0/n3):.6f}        [E, the n={n3} cluster itself]")

# ---------- [C] reference rates ----------
log("\n[C] reference rates for a *tree of blocks* (exact block multiplicativity)")
log(f"    120^(1/6) = {120**(1/6):.6f}   (K_6, 6 plaquettes)")
log(f"    120^(1/7) = {120**(1/7):.6f}   (UR2 bridge chain, 7k-1 law)")
log(f"    120^(1/5) = {120**(1/5):.6f}   (shared-plaquette block product, 6+5k law)")
log(f"    125400^(1/11) = {125400**(1/11):.6f}  (n=11 measured, block product 14400)")
for (k, (n, m, rk)) in res.items():
    log(f"    k={k}: measured rate {m**(1/n):.6f} vs name-plate block rate {120**(1/(n if k==1 else 5*n/(1+0))):.6f}"
        if False else f"    k={k}: measured rate {m**(1/n):.6f}")

# ---------- [D] averaged readings ----------
log("\n[D] the averaged readings (arithmetic sum over clusters) -- KP4 sec.1 m_n")
m_n = {1: 1, 2: 20, 3: 514, 4: 14748, 5: 450255}
for n, m in m_n.items():
    log(f"    n={n}: m_n=sum|phi^U|={m}  m_n^{{1/n}}={m**(1/n):.6f}   (sup over clusters at fixed n <= m_n)")
log(f"    => the summed/averaged layer rate m_5^(1/5) = {450255**(1/5):.6f} is "
    f"{450255**(1/5)/2.90724:.1f}x the measured sup 2.9072")

# ---------- [E] verdict numbers ----------
log("\n[E] verdict numbers")
best = max(m**(1/n) for (n, m, rk) in res.values())
log(f"    sup_n |phi^U|^{{1/n}} >= max over computed witnesses = {best:.6f}")
log(f"    K_hat >= 4*{best:.6f} = {4*best:.6f}   vs threshold {THR:.6f}  "
    f"ratio {4*best/THR:.4f}")
log(f"    n=11 point: 125400^(1/11) = {125400**(1/11):.6f}  -> K_hat >= {4*125400**(1/11):.6f}")
log(f"    T/4 = {THR/4:.13f} ;  |phi|^{{1/n}} must fall below this for K_hat<T")
log(f"\nelapsed {time.time()-T0:.1f}s")
