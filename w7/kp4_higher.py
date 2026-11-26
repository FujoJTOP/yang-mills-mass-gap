"""kp4_higher.py -- KP4: where does sup_n |phi^U_n|^{1/n} sit, and hence what is K_hat?

Q: max |phi^U| over connected plaquette clusters of size n -- which n wins?
   K_hat := sup_n (prod_l R(l) * |phi^U|)^{1/n} = 4 * sup_n |phi^U|^{1/n}   (prod R)^{1/n}=4 on the cube family.

Exact integer phi^U via the vertex-subset DP (Mayer/Ursell, Fernandez-Procacci (2.4)):
   H(S)=1 iff E(S)=empty,  F(S)=H(S)-sum_{T proper, min(S) in T} F(T)H(S\T),   phi^U=F(V).
O(3^n); n<=13 fine.

Parts
  A  lattice pool + verify Delta=20, m_n = 1,20,514,14748,450255 (n<=5)  [E]
  B  named graphs: K_k, cube K_{2,2,2}, K_6+q~r
  C  every 7-cluster containing a fixed link-bundle K_6 (+1 plaquette)   [E] exhaustive in that family
  D  every 8-cluster containing a fixed link-bundle K_6 (+2 plaquettes)  [E] exhaustive in that family
  E  two link-bundles sharing one plaquette (n=11) and UR2 bridge chain (n=13) exact
  F  beam search n=7,8 unrestricted (no K_6 assumption) -- lower bound check
  G  table + K_hat
"""
import itertools, sys, time, math

D = 4
T0 = time.time()


def log(*a):
    print(*a, flush=True)


def unit(i):
    v = [0] * D
    v[i] = 1
    return tuple(v)


def add(u, v):
    return tuple(u[k] + v[k] for k in range(D))


def sub_(u, v):
    return tuple(u[k] - v[k] for k in range(D))


DIRS = [(a, b) for a in range(D) for b in range(a + 1, D)]


# ---------------- exact phi^U ----------------
def phi_int(adj, n):
    """adj[i] = bitmask of neighbours of i.  Returns phi^U (signed integer)."""
    size = 1 << n
    full = size - 1
    H = bytearray(size)
    H[0] = 1
    for S in range(1, size):
        low = S & -S
        i = low.bit_length() - 1
        rest = S ^ low
        H[S] = 1 if (adj[i] & rest) == 0 and H[rest] else 0
    F = [0] * size
    for S in range(1, size):
        low = S & -S
        rest = S ^ low
        tot = H[S]
        t = rest
        while True:
            blk = low | t
            if blk != S:
                tot -= F[blk] * H[S ^ blk]
            if t == 0:
                break
            t = (t - 1) & rest
        F[S] = tot
    return F[full]


def phi_rows(rows, n):
    return phi_int(rows, n)


# ---------------- lattice ----------------
def build_pool(R):
    sites = sorted(x for x in itertools.product(range(-R, R + 1), repeat=D)
                   if sum(map(abs, x)) <= R)
    plaq = [(x, a, b) for x in sites for (a, b) in DIRS]
    idx = {p: i for i, p in enumerate(plaq)}
    P = len(plaq)
    linkmap = {}
    for i, (x, a, b) in enumerate(plaq):
        for L in ((x, a), (add(x, unit(a)), b), (x, b), (add(x, unit(b)), a)):
            linkmap.setdefault(L, []).append(i)
    nb = [set() for _ in range(P)]
    sizes = set()
    for L, ps in linkmap.items():
        sizes.add(len(ps))
        for i in range(len(ps)):
            for j in range(i + 1, len(ps)):
                nb[ps[i]].add(ps[j])
                nb[ps[j]].add(ps[i])
    return plaq, idx, linkmap, nb, sorted(sizes), P


def rows_of(vs, idxset, nb):
    n = len(vs)
    rows = [0] * n
    for i, v in enumerate(vs):
        r = 0
        for j, w in enumerate(vs):
            if j != i and w in nb[v]:
                r |= 1 << j
        rows[i] = r
    return rows


# ---------------- A: verify ----------------
R = 2 * (8 - 1)   # n<=8 needs L1 reach 2*(n-1)=14
plaq, idx, linkmap, nb, lsizes, P = build_pool(R)
log(f"[A] pool R={R}: {P} plaquettes, links carry {lsizes}, links={len(linkmap)}")
p0 = idx[((0, 0, 0, 0), 0, 1)]
log(f"[A] deg(p0) = {len(nb[p0])}   (expect 20)")


def shell_enum(nmax):
    """connected sets containing p0, as frozensets of pool ids, level by level"""
    levels = {1: [frozenset([p0])]}
    for n in range(2, nmax + 1):
        t1 = time.time()
        cur = set()
        for s in levels[n - 1]:
            cand = set()
            for v in s:
                cand |= nb[v]
            cand -= s
            for c in cand:
                cur.add(s | {c})
        levels[n] = list(cur)
        log(f"    n={n}: {len(cur)} sets ({time.time()-t1:.1f}s)")
    return levels


log("[A] shell enumeration for m_n check (n<=5):")
lv = shell_enum(5)
for n in range(1, 6):
    tot = 0
    mx = 0
    for s in lv[n]:
        vs = sorted(s)
        ph = abs(phi_int(rows_of(vs, idx, nb), n))
        tot += ph
        mx = max(mx, ph)
    log(f"    n={n}: sets={len(lv[n])} m_n={tot} m_n^(1/n)={tot**(1/n):.6f} max|phi|={mx} "
        f"max^(1/n)={mx**(1/n):.6f}")

# ---------------- B: named graphs ----------------
log("\n[B] named graphs (exact integer phi^U)")


def from_edges(n, E):
    rows = [0] * n
    for (i, j) in E:
        rows[i] |= 1 << j
        rows[j] |= 1 << i
    return rows


def report(name, rows, n, extra=""):
    ph = phi_rows(rows, n)
    log(f"    {name:<22} n={n} |E|={sum(bin(r).count('1') for r in rows)//2:>2}  "
        f"phi^U={ph:>9}  |phi|={abs(ph):>9}  ^{1}/{n}={abs(ph)**(1/n):.6f} {extra}")
    return abs(ph)


for k in range(1, 7):
    E = [(i, j) for i in range(k) for j in range(i + 1, k)]
    report(f"K_{k}", from_edges(k, E), k)
# cube K_{2,2,2}: 8 vertices? no -- the 6 cube faces: cube's face adjacency graph = K_{2,2,2} with 6 vertices
E = []
# cube faces: pairs of opposite? use K_{2,2,2} on parts {0,1},{2,3},{4,5}
for i in range(2):
    for j in range(2, 4):
        E.append((i, j))
for i in range(2):
    for j in range(4, 6):
        E.append((i, j))
for i in range(2, 4):
    for j in range(4, 6):
        E.append((i, j))
report("K_{2,2,2} (cube)", from_edges(6, E), 6, "<- C06/UR3: phi=-64")
# K_6 + one vertex joined to r of the six
for r in (1, 2, 3):
    E = [(i, j) for i in range(6) for j in range(i + 1, 6)] + [(6, 0), (6, 1), (6, 2)][:r]
    report(f"K_6 + q~{r}", from_edges(7, E), 7)
# K_6 minus an edge
E = [(i, j) for i in range(6) for j in range(i + 1, 6) if (i, j) != (0, 1)]
report("K_6 - e", from_edges(6, E), 6)
# C_4
report("C_4", from_edges(4, [(0, 1), (1, 2), (2, 3), (3, 0)]), 4)

# ---------------- C: all 7-clusters containing one link-bundle ----------------
log("\n[C] all 7-clusters = (6 plaquettes of link (0,0,0,0),dir1) + 1 adjacent plaquette")
L0 = ((0, 0, 0, 0), 0)
bundle = sorted(linkmap[L0])
assert len(bundle) == 6
rowsB = rows_of(bundle, idx, nb)
assert all(bin(r).count("1") == 5 for r in rowsB), "bundle must be K_6"
log(f"    bundle = 6 plaquettes at one link, K_6 confirmed, phi^U={phi_rows(rowsB,6)}")
adj = set()
for v in bundle:
    adj |= nb[v]
adj -= set(bundle)
log(f"    plaquettes adjacent to the bundle: {len(adj)}")
best7 = (0, None)
for q in sorted(adj):
    vs = bundle + [q]
    ph = abs(phi_rows(rows_of(vs, idx, nb), 7))
    if ph > best7[0]:
        best7 = (ph, q)
    log(f"      +{plaq[q]}  |phi^U|={ph}")
log(f"    ==> max 7-cluster in this family: |phi^U|={best7[0]}  ({best7[0]**(1/7):.6f})")

# ---------------- D: all 8-clusters containing the same bundle ----------------
log("\n[D] all 8-clusters = bundle + 2 extra plaquettes (any 2 st. connected)")
adj2 = set(adj)
for v in adj:
    adj2 |= nb[v]
adj2 -= set(bundle)
cand = sorted(adj | (adj2 - adj))
log(f"    candidate extra plaquettes within distance 2: {len(cand)}")
seen = set()
best8 = (0, None)
t1 = time.time()
for a in range(len(cand)):
    for b in range(a + 1, len(cand)):
        q1, q2 = cand[a], cand[b]
        if not (q2 in nb[q1] or q1 in adj or q2 in adj):
            continue
        vs = tuple(sorted(bundle + [q1, q2]))
        if vs in seen:
            continue
        seen.add(vs)
        ph = abs(phi_rows(rows_of(list(vs), idx, nb), 8))
        if ph > best8[0]:
            best8 = (ph, vs)
log(f"    tested {len(seen)} distinct 8-sets ({time.time()-t1:.1f}s)")
log(f"    ==> max 8-cluster in this family: |phi^U|={best8[0]}  ({best8[0]**(1/8):.6f})")
log(f"        witness: {[plaq[v] for v in best8[1]]}" if best8[1] else "")

# ---------------- E: two bundles sharing a plaquette; UR2 bridge chain ----------------
log("\n[E] multi-bundle configurations")
L0 = ((0, 0, 0, 0), 0)          # link at origin, dir 0
L1 = ((0, 0, 0, 0), 1)          # link at origin, dir 1  -> shares plaquette ((0,0,0,0),0,1) with L0
B0 = sorted(linkmap[L0])
B1 = sorted(linkmap[L1])
sh = set(B0) & set(B1)
log(f"    bundles at the two links share {len(sh)} plaquette(s): {[plaq[x] for x in sh]}")
vs = sorted(set(B0) | set(B1))
rows = rows_of(vs, idx, nb)
ph = phi_rows(rows, 11)
log(f"    union (n={len(vs)}, |E|={sum(bin(r).count('1') for r in rows)//2}): phi^U={ph} "
    f"|phi|={abs(ph)} ^{1}/11={abs(ph)**(1/11):.6f}   (block-product would be 120^2=14400)")

# UR2 construction: bundle at ((0,0,0,0),0) and at ((0,0,0,3),0)?? rebuild from UR2 text:
# links (0,(0,0,0,0)) and (0,(0,0,3,0)), bridging single plaquette P(0,2,(0,0,1,0)) -> n=13, phi=14400
# links are (site, dir); UR2 wrote (0,dir) meaning site 0 dir 'dir'.
try:
    La = ((0, 0, 0, 0), 0)
    Lb = ((0, 0, 0, 3), 0)
    Ba = sorted(linkmap[La])
    Bb = sorted(linkmap[Lb])
    # bridging plaquette with base 0, dirs 0,2
    bp = idx[((0, 0, 0, 0), 0, 2)]
    vs2 = sorted(set(Ba) | set(Bb) | {bp})
    rows2 = rows_of(vs2, idx, nb)
    ph2 = phi_rows(rows2, len(vs2))
    log(f"    UR2 two-bundle+bridge (n={len(vs2)}): phi^U={ph2} |phi|={abs(ph2)} "
        f"^{1}/{len(vs2)}={abs(ph2)**(1/len(vs2)):.6f}   (UR2 claimed n=13, 14400)")
except KeyError as e:
    log(f"    UR2 rebuild: KeyError {e}")

# ---------------- F: unrestricted beam search n=7,8 ----------------
log("\n[F] unrestricted beam search (no K_6 assumption)")


def beam(nmax, width):
    cur = {frozenset([p0])}
    best = {}
    for n in range(2, nmax + 1):
        nxt = set()
        for s in cur:
            cand = set()
            for v in s:
                cand |= nb[v]
            cand -= s
            for c in cand:
                nxt.add(s | {c})
        scored = []
        for s in nxt:
            vs = sorted(s)
            ph = abs(phi_int(rows_of(vs, idx, nb), n))
            scored.append((ph, s))
            if n not in best or ph > best[n][0]:
                best[n] = (ph, s)
        scored.sort(key=lambda t: -t[0])
        cur = {s for _, s in scored[:width]}
        log(f"    n={n}: expanded {len(nxt)} -> kept {len(cur)}; best |phi|={best[n][0]} "
            f"({best[n][0]**(1/n):.6f})")
    return best


best = beam(8, 400)
for n in range(2, 9):
    if n in best:
        log(f"    beam n={n}: max|phi^U|>={best[n][0]}  ^{1}/{n}={best[n][0]**(1/n):.6f}")

# ---------------- G: table ----------------
log("\n[G] table")
log(f"{'n':>2} {'max|phi^U|':>12} {'^{1}/n':>10}   source")
rows_tab = [(1, 1, "K_1"), (2, 2, "K_2"), (3, 6, "K_3"), (4, 6, "K_4 [E]"),
            (5, 24, "K_5 [E]"), (6, 120, "K_6 [T/E]"),
            (7, best7[0], "K_6+q (all 7-clusters w/ K_6) [E]"),
            (8, best8[0], "K_6+2q (all 8-clusters w/ K_6) [E]")]
for (n, m, src) in rows_tab:
    log(f"{n:>2} {m:>12} {m**(1/n):>10.6f}   {src}")
sup = max(m ** (1 / n) for (n, m, _) in rows_tab)
log(f"\n    sup_n |phi^U|^{1/n} = {sup:.6f}  (n=6, 120^{1/6})")
log(f"    K_hat = 4 * sup = {4*sup:.6f}      (threshold 7.0632853; K_hat/7.0632853 = {4*sup/7.0632853:.4f})")
log(f"    120^(1/7)={120**(1/7):.6f}  144^(1/7)={144**(1/7):.6f}  24576^(1/8)={24576**(1/8):.6f}")
log(f"\nelapsed {time.time()-T0:.1f}s")
