"""ur4_verify.py -- independent checks of the ur4_enum pipeline on random connected sets.

For a random 6-set through p0:
   (a) mycanon (WL-restricted min) -> rebuild graph -> bruteforce canonical form (min over all 6! perms)
   (b) bruteforce canonical form of the original labelled graph
   require (a) == (b)   [canonicalisation is faithful]
   (c) phi^U computed from the rebuilt graph == phi^U computed directly from the labelled graph
   (d) phi^U depends only on the canonical class: two relabelled copies give the same phi
Also checks the exact K_6 star and the exact cube Q_3 (6 faces of a unit 3-cube) values.
"""
import itertools, random, sys

D = 4
DIRS = [(a, b) for a in range(D) for b in range(a + 1, D)]


def unit(i):
    v = [0] * D
    v[i] = 1
    return tuple(v)


def add(u, v):
    return tuple(u[k] + v[k] for k in range(D))


def links_of(x, a, b):
    return frozenset({(x, a), (add(x, unit(a)), b), (x, b), (add(x, unit(b)), a)})


def build_graph(plaqs):
    L = [links_of(*p) for p in plaqs]
    n = len(plaqs)
    rows = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if not L[i].isdisjoint(L[j]):
                rows[i] |= 1 << j
                rows[j] |= 1 << i
    return rows


def phi_int(rows, n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if (rows[i] >> j) & 1]
    m = len(edges)
    tot = 0
    for mask in range(1 << m):
        par = list(range(n))

        def find(u):
            while par[u] != u:
                par[u] = par[par[u]]
                u = par[u]
            return u
        k = 0
        bits = 0
        for t in range(m):
            if (mask >> t) & 1:
                bits += 1
                a, b = edges[t]
                ra, rb = find(a), find(b)
                if ra != rb:
                    par[ra] = rb
                    k += 1
        if k == n - 1:
            tot += -1 if (bits & 1) else 1
    return tot


def brute_canon(rows, n):
    best = None
    for perm in itertools.permutations(range(n)):
        k = 0
        for i in range(n):
            oi = perm[i]
            for j in range(i + 1, n):
                oj = perm[j]
                if ((rows[oi] >> oj) & 1) or ((rows[oj] >> oi) & 1):
                    k |= 1 << (i * n + j)
        if best is None or k < best:
            best = k
    return best


def wl_colors(rows, n):
    c = [bin(rows[i]).count("1") for i in range(n)]
    for _ in range(4):
        sig = []
        for i in range(n):
            x = rows[i]
            lst = []
            while x:
                b = x & -x
                lst.append(c[b.bit_length() - 1])
                x ^= b
            lst.sort()
            sig.append((c[i], tuple(lst)))
        order = {s: q for q, s in enumerate(sorted(set(sig)))}
        nc = [order[s] for s in sig]
        if nc == c:
            break
        c = nc
    return c


def mycanon(rows, n):
    c = wl_colors(rows, n)
    groups = {}
    for i, ci in enumerate(c):
        groups.setdefault(ci, []).append(i)
    best = None
    gl = [groups[q] for q in sorted(groups)]
    for perm in itertools.product(*[itertools.permutations(g) for g in gl]):
        order = [v for t in perm for v in t]
        k = 0
        for i in range(n):
            oi = order[i]
            for j in range(i + 1, n):
                oj = order[j]
                if ((rows[oi] >> oj) & 1) or ((rows[oj] >> oi) & 1):
                    k |= 1 << (i * n + j)
        if best is None or k < best:
            best = k
    return best


def rows_from_key(key, n):
    rows = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if (key >> (i * n + j)) & 1:
                rows[i] |= 1 << j
                rows[j] |= 1 << i
    return rows


def main():
    random.seed(20261125)
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
    p0 = ((0, 0, 0, 0), 0, 1)
    bad = 0
    ncalls = 0
    seen = {}
    # random 6-sets by random growth inside a fixed local window
    sites = [x for x in itertools.product(range(-10, 11), repeat=D) if sum(map(abs, x)) <= 10]
    plaqs = [(x, a, b) for x in sites for (a, b) in DIRS]
    lk = {p: links_of(*p) for p in plaqs}
    nbr = {p: set() for p in plaqs}
    linkmap = {}
    for p in plaqs:
        for L in lk[p]:
            linkmap.setdefault(L, []).append(p)
    for L, ps in linkmap.items():
        for i in range(len(ps)):
            for j in range(i + 1, len(ps)):
                nbr[ps[i]].add(ps[j])
                nbr[ps[j]].add(ps[i])

    for t in range(trials):
        cur = [p0]
        S = {p0}
        while len(cur) < 6:
            cands = set()
            for p in cur:
                cands |= nbr[p]
            cands -= S
            cands = [c for c in cands if c in nbr]
            if not cands:
                break
            q = random.choice(sorted(cands))
            cur.append(q)
            S.add(q)
        if len(cur) < 6:
            continue
        rows = build_graph(cur)
        k = mycanon(rows, 6)
        j = brute_canon(rows, 6)
        rec = rows_from_key(k, 6)
        ok = (brute_canon(rec, 6) == j) and (phi_int(rec, 6) == phi_int(rows, 6))
        # relabel invariance
        perm = list(range(6))
        random.shuffle(perm)
        rows2 = [0] * 6
        for i in range(6):
            for jj in range(6):
                if (rows[perm[i]] >> perm[jj]) & 1:
                    rows2[i] |= 1 << jj
        ok = ok and (mycanon(rows2, 6) == k) and (phi_int(rows2, 6) == phi_int(rows, 6))
        if not ok:
            bad += 1
            print("FAIL at trial", t, "phi=", phi_int(rows, 6), "phi(rec)=", phi_int(rec, 6))
        ncalls += 1
        seen[abs(phi_int(rows, 6))] = seen.get(abs(phi_int(rows, 6)), 0) + 1
    print(f"checked {ncalls} random 6-sets: failures = {bad}")
    print("|phi| values seen:", dict(sorted(seen.items())))

    # exact anchors
    for L in sorted(lk[p0])[:1]:
        s = [p0] + sorted([p for p in linkmap[L] if p != p0])
        if len(s) == 6:
            print("K_6 on one link of p0: |phi^U| =", abs(phi_int(build_graph(s), 6)), "(expect 120)")
    # cube Q_3 on directions 0,1,2: 6 faces, p0 = face at c=0
    faces = []
    dirs = [0, 1, 2]
    for (a, b) in [(0, 1), (0, 2), (1, 2)]:
        faces.append(((0, 0, 0, 0), a, b))
    for (a, b) in [(0, 1), (0, 2), (1, 2)]:
        x = [0, 0, 0, 0]
        c = 3 - a - b  # third direction among 0,1,2
        x[c] = 1
        faces.append((tuple(x), a, b))
    print("cube Q_3 (6 faces): |phi^U| =", abs(phi_int(build_graph(faces), 6)), "(F07:64 says 64)")


if __name__ == "__main__":
    main()
