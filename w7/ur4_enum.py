"""ur4_enum.py -- m_n = sum_{|Gamma|=n, p0 in Gamma} |phi^U(Gamma)|  (4D plaquette animals through a fixed plaquette).

phi^U(Gamma) = sum_{A subset E(G_Gamma), (V,A) connected spanning} (-1)^{|A|}   (Fernandez-Procacci (2.4))
sign depends only on n:  phi^U = (-1)^{n-1}|phi^U|  (F-P (2.8)); |phi^U| is a positive integer, N-independent
(no Weingarten enters: the cluster coefficients are pure graph combinatorics).

Geometry / pool.  A plaquette is (base site x in Z^4, direction pair a<b).  It carries 4 links
(x,a),(x+e_a,b),(x,b),(x+e_b,a).  Two plaquettes are adjacent iff they share a link.  Sharing a link
can move the base site by L1 = 1 or (when the shared link is the "far" side of both) by L1 = 2, so a
connected set of n plaquettes through p0 has every base site within L1 <= 2(n-1).
Pool = all plaquettes with base site in the L1 ball of radius R = 2*(NMAX-1).

Adjacency is built link-wise (each link lies in exactly 2(d-1)=6 plaquettes), O(P) not O(P^2).
"""
import itertools, json, sys, time

D = 4
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
OUT = sys.argv[2] if len(sys.argv) > 2 else "ur4_enum.json"
R = 2 * (NMAX - 1)

DIRS = [(a, b) for a in range(D) for b in range(a + 1, D)]


def unit(i):
    v = [0] * D
    v[i] = 1
    return tuple(v)


def add(u, v):
    return tuple(u[k] + v[k] for k in range(D))


def build():
    sites = [x for x in itertools.product(range(-R, R + 1), repeat=D)
             if sum(map(abs, x)) <= R]
    sites.sort()
    plaq = [(x, a, b) for x in sites for (a, b) in DIRS]
    P = len(plaq)
    idx = {p: i for i, p in enumerate(plaq)}
    W = max(1, (P - 1).bit_length())
    # link -> list of plaquette ids
    linkmap = {}
    for i, (x, a, b) in enumerate(plaq):
        for L in ((x, a), (add(x, unit(a)), b), (x, b), (add(x, unit(b)), a)):
            linkmap.setdefault(L, []).append(i)
    nbrset = [set() for _ in range(P)]
    sizes = set()
    for L, ps in linkmap.items():
        sizes.add(len(ps))
        for i in range(len(ps)):
            for j in range(i + 1, len(ps)):
                nbrset[ps[i]].add(ps[j])
                nbrset[ps[j]].add(ps[i])
    nbr = [0] * P
    for i in range(P):
        m = 0
        for j in nbrset[i]:
            m |= 1 << j
        nbr[i] = m
    p0 = idx[((0, 0, 0, 0), 0, 1)]
    return plaq, P, W, p0, nbr, sorted(sizes)


def main():
    t0 = time.time()
    plaq, P, W, p0, nbr, lsizes = build()
    print(f"pool: {P} plaquettes (R={R}), key width {W} bits/field, "
          f"links carry {lsizes} plaquettes, deg(p0)={bin(nbr[p0]).count('1')}", flush=True)
    MASK = (1 << W) - 1

    def decode(k, n):
        return [(k >> (W * t)) & MASK for t in range(n)]

    def encode_insert(k, n, c):
        out = []
        placed = False
        for t in range(n):
            v = (k >> (W * t)) & MASK
            if not placed and c < v:
                out.append(c)
                placed = True
            out.append(v)
        if not placed:
            out.append(c)
        r = 0
        for t, v in enumerate(out):
            r |= v << (W * t)
        return r

    # ---- enumerate connected sets containing p0 ----
    levels = {1: {p0}}
    for n in range(2, NMAX + 1):
        prev = levels[n - 1]
        cur = set()
        t1 = time.time()
        for k in prev:
            lst = decode(k, n - 1)
            m = 0
            sm = 0
            for i in lst:
                m |= nbr[i]
                sm |= 1 << i
            cand = m & ~sm
            while cand:
                b = cand & -cand
                cand ^= b
                cur.add(encode_insert(k, n - 1, b.bit_length() - 1))
        levels[n] = cur
        print(f"n={n}: {len(cur)} sets  ({time.time()-t1:.1f}s)", flush=True)

    # ---- phi^U per set via canonical-class cache ----
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

    def pack(rows, order, n):
        k = 0
        for i in range(n):
            oi = order[i]
            for j in range(i + 1, n):
                oj = order[j]
                bit = ((rows[oi] >> oj) & 1) or ((rows[oj] >> oi) & 1)
                if bit:
                    k |= 1 << (i * n + j)
        return k

    def canon(rows, n):
        c = wl_colors(rows, n)
        groups = {}
        for i, ci in enumerate(c):
            groups.setdefault(ci, []).append(i)
        if len(groups) == n:
            return pack(rows, sorted(range(n), key=lambda i: c[i]), n)
        gl = [groups[q] for q in sorted(groups)]
        best = None
        for perm in itertools.product(*[itertools.permutations(g) for g in gl]):
            order = [v for t in perm for v in t]
            key = pack(rows, order, n)
            if best is None or key < best:
                best = key
        return best

    def phi_int(rows, n):
        edges = [(i, j) for i in range(n) for j in range(i + 1, n) if (rows[i] >> j) & 1]
        m = len(edges)
        tot = 0
        for mask in range(1 << m):
            par = list(range(n))
            k = 0
            bits = 0
            for t in range(m):
                if (mask >> t) & 1:
                    bits += 1
                    a, b = edges[t]
                    ra = a
                    while par[ra] != ra:
                        par[ra] = par[par[ra]]
                        ra = par[ra]
                    rb = b
                    while par[rb] != rb:
                        par[rb] = par[par[rb]]
                        rb = par[rb]
                    if ra != rb:
                        par[ra] = rb
                        k += 1
            if k == n - 1:
                tot += -1 if (bits & 1) else 1
        return tot

    def rows_from_key(key, n):
        rows = [0] * n
        for i in range(n):
            for j in range(i + 1, n):
                if (key >> (i * n + j)) & 1:
                    rows[i] |= 1 << j
                    rows[j] |= 1 << i
        return rows

    res = {}
    classes = {}
    for n in range(1, NMAX + 1):
        t1 = time.time()
        cache = {}
        dist = {}
        tot = 0
        cnt = 0
        for k in levels[n]:
            lst = decode(k, n)
            rows = [0] * n
            for i in range(n):
                ni = nbr[lst[i]]
                ri = 0
                for j in range(n):
                    if j != i and ((ni >> lst[j]) & 1):
                        ri |= 1 << j
                rows[i] = ri
            ck = canon(rows, n)
            hit = cache.get(ck)
            if hit is None:
                ph = phi_int(rows_from_key(ck, n), n)
                hit = (abs(ph), ph, ck)
                cache[ck] = hit
            a = hit[0]
            tot += a
            cnt += 1
            dist[a] = dist.get(a, 0) + 1
        res[n] = dict(sets=cnt, m_n=tot, dist={str(a): b for a, b in sorted(dist.items())},
                      nclasses=len(cache), maxabs=max(dist),
                      mroot=tot ** (1.0 / n), maxroot=max(dist) ** (1.0 / n))
        classes[n] = {hex(k): v[1] for k, v in cache.items()}
        print(f"n={n}: sets={cnt} m_n={tot} m_n^(1/n)={res[n]['mroot']:.6f} "
              f"max|phi|={max(dist)} classes={len(cache)}  ({time.time()-t1:.1f}s)", flush=True)

    with open(OUT, "w") as f:
        json.dump(dict(NMAX=NMAX, R=R, pool=P, res=res, classes=classes,
                       elapsed=time.time() - t0), f, indent=1)

    print("\n=== table ===", flush=True)
    print(f"{'n':>2} {'#sets':>10} {'m_n':>13} {'m_n^(1/n)':>11} {'ratio':>9} "
          f"{'max|phi|':>9} {'max^(1/n)':>10} {'#classes':>9}")
    prevm = None
    for n in range(1, NMAX + 1):
        r = res[n]
        ratio = (r["m_n"] / prevm) if prevm else float("nan")
        print(f"{n:>2} {r['sets']:>10} {r['m_n']:>13} {r['mroot']:>11.6f} {ratio:>9.4f} "
              f"{r['maxabs']:>9} {r['maxroot']:>10.6f} {r['nclasses']:>9}")
        prevm = r["m_n"]
    print("\n=== |phi^U| distribution ===", flush=True)
    for n in range(1, NMAX + 1):
        print(f"n={n}: " + "  ".join(f"{a}:{b}" for a, b in res[n]["dist"].items()), flush=True)
    print(f"\nelapsed {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
