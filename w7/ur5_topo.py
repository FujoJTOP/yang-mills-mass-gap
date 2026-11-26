"""ur5_topo.py -- topology classification of connected plaquette sets Gamma subset of Z^4, |Gamma| = n <= 6.

Conventions identical to ur4_enum.py:
  plaquette = (base site x in Z^4, plane a<b);  its 4 links are
      (x,a) , (x+e_a,b) , (x,b) , (x+e_b,a) ;  link = (vertex, direction).
  adjacency  = share a link  (each link lies in exactly 2(d-1)=6 plaquettes, d=4);
  site-adjacency = share >=1 vertex (reported separately, weaker notion).
  canonical form = min over the 384-element hypercubic group B4 (signed coordinate
  permutations) after quotienting translations -> "free" (= symmetry-class) count.

Exact cross-check: the number of *fixed* sets through p0=(0,0,0,0;0,1) is
      Nfix(n) = sum over free shapes S of n * |Orb(S)| / 6
( |{g in B4 : g(q)=p0}| = |Stab(p0)| = 64 for each q, so #placements containing p0
  = n*64/|Stab(S)| = n*|Orb(S)|/6 ).  Compare with ur4_enum.log: 1, 20, 458, 11132, 280175.

Output: counts per n, per-shape invariants (e = #shared links, link-degree multiset, V, chi,
closedness), family labels, and the n=6 prediction for ur4.
"""
import itertools, json, math, sys, time
import numpy as np

D = 4
PLANES = [(a, b) for a in range(D) for b in range(a + 1, D)]


def unit(i):
    v = [0] * D
    v[i] = 1
    return tuple(v)


def add(u, v):
    return tuple(u[k] + v[k] for k in range(D))


def sub(u, v):
    return tuple(u[k] - v[k] for k in range(D))


def corners(pl):
    base, a, b = pl
    return (base, add(base, unit(a)), add(base, unit(b)),
            add(add(base, unit(a)), unit(b)))


def links(pl):
    base, a, b = pl
    return ((base, a), (add(base, unit(a)), b), (base, b), (add(base, unit(b)), a))


def plaqs_through_link(v, a):
    out = []
    for b in range(D):
        if b == a:
            continue
        lo, hi = (a, b) if a < b else (b, a)
        out.append((v, lo, hi))
        out.append((sub(v, unit(b)), lo, hi))
    return out


# ---------------- B4 group, vectorised ----------------
def make_group():
    perms, signs = [], []
    for pm in itertools.permutations(range(D)):
        for sg in itertools.product((1, -1), repeat=D):
            perms.append(pm)
            signs.append(sg)
    return np.array(perms, dtype=np.int64), np.array(signs, dtype=np.int64)


PERM, SIGN = make_group()
NG = len(PERM)


def canon_np(S):
    """S: list of (base,a,b). returns (key tuple of packed int codes, orbit size, representative S)"""
    n = len(S)
    X = np.array([corners(pl) for pl in S], dtype=np.int64)     # (n,4,4)
    Y = np.empty((NG, n, 4, D), dtype=np.int64)
    for c in range(D):                                          # output coordinate c
        Y[:, :, c, :] = SIGN[:, c, None, None] * np.moveaxis(np.take(X, PERM[:, c], axis=1), 1, 0)
    mn = Y.min(axis=2).min(axis=1)                              # (NG,D) global min corner
    Z = Y - mn[:, None, None, :]                                # all coords >= 0
    M = Z.min(axis=2)                                           # (NG,n,D) per-plaquette min corner
    off = Z - M[:, :, None, :]
    L1 = np.abs(off).sum(-1)                                    # (NG,n,4)
    flag = (L1 == 1)
    w = np.where(flag[..., None], np.abs(off), 0)
    dvec = w.argmax(-1)                                         # (NG,n,4)
    order = np.argsort(~flag, axis=-1, kind='stable')           # flagged corners first
    d1 = np.take_along_axis(dvec, order[..., :1], axis=-1)[..., 0]
    d2 = np.take_along_axis(dvec, order[..., 1:2], axis=-1)[..., 0]
    plane = np.minimum(d1, d2) * 4 + np.maximum(d1, d2)
    B = 16
    code = (((M[..., 0] * B + M[..., 1]) * B + M[..., 2]) * B + M[..., 3]) * 32 + plane
    Sc = np.sort(code, axis=1)                                  # (NG,n)
    order_rows = np.lexsort(Sc.T[::-1])
    gi = order_rows[0]
    key = tuple(int(v) for v in Sc[gi])
    orbit = len(np.unique(Sc, axis=0))
    # rebuild the representative geometry from group element gi
    pm, sg = PERM[gi], SIGN[gi]
    rep = []
    for pl in S:
        cs = [tuple(int(sg[k] * c[pm[k]]) for k in range(D)) for c in corners(pl)]
        m0 = tuple(min(c[k] for c in cs) for k in range(D))
        ds = sorted(i for i in range(D) for c in [None] if True)  # placeholder
        dirs = []
        for c in cs:
            dd = tuple(c[k] - m0[k] for k in range(D))
            if sum(abs(x) for x in dd) == 1:
                dirs.append(dd.index(1) if 1 in dd else dd.index(-1))
        dirs.sort()
        rep.append((m0, dirs[0], dirs[1]))
    return key, orbit, tuple(sorted(rep))


def shift_to_origin(S):
    cs = [c for pl in S for c in corners(pl)]
    m0 = tuple(min(c[k] for c in cs) for k in range(D))
    return tuple(sorted((sub(pl[0], m0), pl[1], pl[2]) for pl in S))


def invariants(S):
    """S: tuple of (base,a,b)"""
    n = len(S)
    lk = {}
    for pl in S:
        for L in links(pl):
            lk.setdefault(L, []).append(pl)
    degs = sorted(len(v) for v in lk.values())
    shared = sum(1 for v in lk.values() if len(v) >= 2)
    verts = set()
    for pl in S:
        for c in corners(pl):
            verts.add(c)
    V = len(verts)
    E = len(lk)
    chi = V - E + n
    closed = (2 * shared == 4 * n)
    # family-ish labels
    planes = set((pl[1], pl[2]) for pl in S)
    sig = (shared, tuple(degs), V, chi, closed)
    return dict(n=n, e=shared, degs=degs, V=V, E=E, chi=chi, closed=closed,
                nplanes=len(planes), sig=sig)


def enumerate_shapes(nmax, tcap=900.0, log=print):
    t0 = time.time()
    start = ((0, 0, 0, 0), 0, 1)
    level = {canon_np([start])[0]: (1, shift_to_origin([start]))}   # key -> (orbit, rep)
    res = {}
    for n in range(1, nmax + 1):
        res[n] = dict(level)
        log("  n=%d : %d free shapes  (n*orbit/6 sum = %d)  [%.1fs]" %
            (n, len(level), sum(n * o for o, _ in level.values()), time.time() - t0))
        if n == nmax:
            break
        nxt = {}
        calls = 0
        for key, (orb, S) in level.items():
            Sset = set(S)
            bd = set()
            for pl in S:
                for L in links(pl):
                    for q in plaqs_through_link(*L):
                        if q not in Sset:
                            bd.add(q)
            for q in sorted(bd):
                T = tuple(sorted(Sset | {q}))
                k2, o2, rep2 = canon_np(T)
                calls += 1
                if k2 not in nxt:
                    nxt[k2] = (o2, rep2)
            if time.time() - t0 > tcap:
                log("  !! time cap hit at n=%d after %d canon calls" % (n + 1, calls))
                return res
        level = nxt
    return res


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    print("=== connected plaquette sets in Z^4, |Gamma| = n <= %d ===" % nmax, flush=True)
    print("adjacency = share a link;  quotient = translations x B4 (|B4| = 384)", flush=True)
    res = enumerate_shapes(nmax)
    t0 = time.time()
    print("\n=== counts ===")
    ur4 = {1: 1, 2: 20, 3: 458, 4: 11132, 5: 280175}
    tot = {}
    for n in sorted(res):
        shapes = res[n]
        Nfix = sum(n * o for o, _ in shapes.values())
        tot[n] = Nfix
        chk = ""
        if n in ur4:
            chk = "  ur4=%d  %s" % (ur4[n], "MATCH" if ur4[n] == Nfix else "MISMATCH")
        print("  n=%d : %7d free shapes   Nfix(=sets through p0)=%9d%s" % (n, len(shapes), Nfix, chk))

    print("\n=== family/invariant census (free shapes) ===")
    for n in sorted(res):
        buckets = {}
        for key, (orb, S) in res[n].items():
            iv = invariants(S)
            buckets.setdefault((iv['e'], iv['closed'], iv['chi'], iv['nplanes']), 0)
            buckets[(iv['e'], iv['closed'], iv['chi'], iv['nplanes'])] += 1
        print("  n=%d  (e=shared links, closed, chi, #planes used) -> #shapes" % n)
        for k in sorted(buckets):
            print("       e=%2d closed=%-5s chi=%3d planes=%d : %d" % (k[0], k[1], k[2], k[3], buckets[k]))

    print("\n=== named special shapes ===")
    for n in sorted(res):
        for key, (orb, S) in res[n].items():
            iv = invariants(S)
            labs = []
            if iv['nplanes'] == 1:
                labs.append("planar patch(one 2-plane)")
            if iv['closed']:
                labs.append("CLOSED surface (every link shared twice)")
            if iv['e'] == n - 1:
                labs.append("tree-like (e=n-1, max #free indices)")
            if iv['e'] == 2 * n:
                labs.append("max sharing")
            if labs:
                print("  n=%d %-46s e=%2d chi=%3d V=%2d : %s" %
                      (n, str(S[0]) + ("..." if n > 1 else ""), iv['e'], iv['chi'], iv['V'], ", ".join(labs)))

    # dump n=6 shapes that are closed or maximal-sharing, plus all n<=4 for inspection
    out = {}
    for n in sorted(res):
        rows = []
        for key, (orb, S) in res[n].items():
            iv = invariants(S)
            rows.append(dict(rep=[list(map(list, [pl[0]])) + [pl[1], pl[2]] for pl in S],
                             orbit=orb, e=iv['e'], degs=iv['degs'], V=iv['V'], chi=iv['chi'],
                             closed=iv['closed'], nplanes=iv['nplanes']))
        out[n] = rows
    with open("ur5_topo.json", "w") as f:
        json.dump({str(k): v for k, v in out.items()}, f)
    print("\nwrote ur5_topo.json ; elapsed %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
