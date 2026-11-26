"""
P3B — can a d=4 plaquette animal be decomposed into d=2 layers?  What bound does that give?

Object (same as lc4/ur1): 4D hypercubic lattice, vertices = plaquettes (2-cells),
edge = shared 1-cell (link).  Degree Delta(d) = 4(2d-3); d=4 -> 20.
Rooted animals A_n = connected n-plaquette sets containing p0 = ((0,0,0,0), plane (0,1)).

A "layer" := (plane type {i,j}, values of the two complementary coordinates).
Within one layer the cells live on a Z^2 spanned by e_i, e_j, and the in-plane
adjacency is exactly the d=2 plaquette graph (= square lattice = polyominoes).

Sections
 S0  Delta(d) = 4(2d-3) by explicit face-graph degrees, d=2..5; in-plane/cross-plane split
 S1  d=2 face graph == Z^2; rooted polyomino counts c_n = n*A_n, lambda_2 reference
 S2  layer test on ALL rooted 4D animals n<=5:
       (i)  is every (plane,layer) piece in-plane connected?  (needed for "layers are 2D animals")
       (ii) is every plane-family's layer-index set connected?  (needed for a product bound)
 S3  arithmetic: lambda_2^k, k-window for lambda_2^k in (lambda_* probe, 48], slack kappa

Usage: py -3.12 p3b_layered.py [maxn]
"""
import sys, time
from decimal import Decimal, getcontext
import ur1_enum as U1

getcontext().prec = 40

D4 = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]


def add(x, d):
    return tuple(a + b for a, b in zip(x, d))


def sub(x, d):
    return tuple(a - b for a, b in zip(x, d))


def pl(i, j):
    return (i, j) if i < j else (j, i)


def e(dim, k):
    v = [0] * dim
    v[k] = 1
    return tuple(v)


def nb(dim, x, plane):
    """all plaquettes of Z^dim sharing exactly one link with plaquette (x, plane)"""
    i, j = plane
    links = [(x, i), (add(x, e(dim, i)), j), (add(x, e(dim, j)), i), (x, j)]
    out = set()
    for y, dd in links:
        for k in range(dim):
            if k == dd:
                continue
            out.add((y, pl(dd, k)))
            out.add((sub(y, e(dim, k)), pl(k, dd)))
    out.discard((x, plane))
    return out


def conn(bases, i, j):
    """is the cell set `bases` connected under in-plane (e_i,e_j) adjacency?"""
    s = set(bases)
    st = [next(iter(s))]
    seen = set(st)
    while st:
        v = st.pop()
        for d in (D4[i], D4[j]):
            for w in (add(v, d), sub(v, d)):
                if w in s and w not in seen:
                    seen.add(w)
                    st.append(w)
    return len(seen) == len(s)


def layers_conn(lts):
    """is the layer-index set (2-tuples of the two complement coords) connected
       under +-1 steps in either coordinate?"""
    s = set(lts)
    st = [next(iter(s))]
    seen = set(st)
    while st:
        v = st.pop()
        for c in (0, 1):
            for sgn in (1, -1):
                w = list(v)
                w[c] += sgn
                w = tuple(w)
                if w in s and w not in seen:
                    seen.add(w)
                    st.append(w)
    return len(seen) == len(s)


def s0():
    print("P3B S0 — Delta(d) = 4(2d-3) from the face graph, and the in-plane split")
    print(f"  {'d':>2} {'deg':>4} {'4(2d-3)':>8}  {'in-plane':>8} {'cross':>6}  check")
    for dim in (2, 3, 4, 5):
        x = (0,) * dim
        tot = nb(dim, x, (0, 1))
        inpl = {q for q in tot if q[1] == (0, 1)}
        print(f"  {dim:>2} {len(tot):>4} {4*(2*dim-3):>8}  {len(inpl):>8} {len(tot)-len(inpl):>6}"
              f"  {'OK' if len(tot) == 4*(2*dim-3) else 'MISMATCH'}")
    print("  cross-check vs ur1_enum.neighbors (d=4): ", end="")
    root = U1.pc((0, 0, 0, 0), 0, 1)
    x, plane = U1.dec(root)
    print(f"ur1 = {len(U1.neighbors(root))}, generic = {len(nb(4, x, plane))}",
          "OK" if len(U1.neighbors(root)) == len(nb(4, x, plane)) else "MISMATCH")
    print("  => d=2: each square has 4 neighbours, all in-plane (cross-plane = 0)")
    print("     d=4: 4 in-plane + 16 cross-plane = 20  (matches lc4_dbg: 4 coplanar / 16 skew)\n")


def s1():
    print("P3B S1 — the d=2 plaquette graph is exactly Z^2 (= square lattice = polyominoes)")
    # 2D: all 4 neighbours of a square are in-plane translates; graph = Z^2
    seen = {(0, 0)}
    lvl = {(0, 0)}
    cnts = [1]
    for n in range(2, 10):
        nxt = set()
        for c in lvl:
            for d in D4[:2]:
                for w in (add(c, d), sub(c, d)):
                    if w not in seen:
                        nxt.add(w)
        # rooted animal enumeration (sets containing (0,0))
        pass
    # rooted animal counts: BFS over sets
    lvl = {((0, 0),)}
    c = []
    for n in range(1, 10):
        c.append(len(lvl))
        if n == 9:
            break
        nxt = set()
        for S in lvl:
            base = set(S)
            for p in S:
                for d in D4[:2]:
                    for q in (add(p, d), sub(p, d)):
                        if q not in base:
                            nxt.add(tuple(sorted(S + (q,))))
        lvl = nxt
    print("  rooted 2D animal counts c_n (root = fixed square):", c)
    print("  literature: fixed polyominoes A_n = 1,2,6,19,63,216,760,2725,9910 -> c_n = n*A_n")
    print("              = 1,4,18,76,315,1296,5320,21800,89190")
    print("  match:", "OK" if c == [1, 4, 18, 76, 315, 1296, 5320, 21800, 89190] else "MISMATCH")
    print(f"  ratios c_n/c_(n-1): {[round(c[i]/c[i-1],4) for i in range(1,len(c))]}")
    print("  lambda_2 (literature, fixed polyominoes) = 4.062570...  [E]\n")
    return c


def s2(maxn):
    print(f"P3B S2 — layer test over ALL rooted 4D animals, n <= {maxn}")
    t0 = time.time()
    levels = U1.animals(maxn)
    a = [len(L) for L in levels]
    print("  a_n =", a, f"[{time.time()-t0:.1f}s]")
    print(f"  {'n':>2} {'a_n':>8} {'bad-piece':>10} {'bad-layerset':>13} {'max#pieces':>10}"
          f" {'1st witness'}")
    first = None
    for n in range(1, maxn + 1):
        badp = badl = 0
        mpiece = 0
        for S in levels[n - 1]:
            pieces = {}
            fam = {}
            for code in S:
                xx, (mu, nu) = U1.dec(code)
                rest = tuple(xx[k] for k in range(4) if k not in (mu, nu))
                pieces.setdefault((mu, nu, rest), []).append(xx)
                fam.setdefault((mu, nu), set()).add(rest)
            mpiece = max(mpiece, len(pieces))
            okp = all(conn(v, mu, nu) for (mu, nu, _), v in pieces.items())
            okl = True
            for (mu, nu), lts in fam.items():
                if not layers_conn(lts):
                    okl = False
                    break
            if not okp:
                badp += 1
                if first is None:
                    first = (n, S)
            if not okl:
                badl += 1
        print(f"  {n:>2} {len(levels[n-1]):>8} {badp:>10} {badl:>13} {mpiece:>10}"
              f"  {'' if first is None or first[0]!=n else 'yes'}")
    if first is not None:
        n, S = first
        print(f"\n  minimal entanglement witness at n={n} ({len(S)} plaquettes):")
        pl_by = {}
        for code in sorted(S):
            xx, (mu, nu) = U1.dec(code)
            rest = tuple(xx[k] for k in range(4) if k not in (mu, nu))
            print(f"    plaquette plane ({mu},{nu}) base {xx}  -> layer {rest}")
            pl_by.setdefault((mu, nu, rest), []).append(xx)
        for key, v in pl_by.items():
            print(f"    layer {key}: bases {v}  in-plane connected? {conn(v, key[0], key[1])}")
        print("  => at least one layer piece is NOT a 2D animal: the animal leaves the layer")
        print("     and comes back (connected through another plane type).")
    else:
        print("  no entanglement witness up to n = %d" % maxn)
    print(f"  [{time.time()-t0:.1f}s]\n")
    return first


def s3():
    print("P3B S3 — what bound could a layer decomposition give?")
    l2 = Decimal("4.0626")                     # literature value, tagged [E]
    print(f"  lambda_2                       = {l2}   [E] (fixed polyominoes, 4.062570...)")
    for k in range(1, 7):
        v = l2 ** k
        print(f"  lambda_2^{k}                     = {v:.6f}")
    lt20 = Decimal(20) ** 20 / Decimal(19) ** 19
    lt18 = Decimal(18) ** 18 / Decimal(17) ** 17
    probe = Decimal("25.17")                   # LC4 ratio probe a6/a5, tagged [E]
    print(f"  probe lambda_* >= a_6/a_5      = {probe}   [E] (LC4.md; ratio probe, not a theorem)")
    print(f"  lambda_* >= a_6^(1/6)          = {Decimal(7218768) ** (Decimal(1)/6):.4f}"
          f"   [T] (rooted counts are supermultiplicative)")
    print(f"  48                             = 48")
    print(f"  lambda_2^2                     = {l2**2:.6f}   (< 48 -> would close the gap)")
    print(f"  lambda_2^3                     = {l2**3:.6f}   (> 54 -> (b))")
    print(f"  lambda_2^(Delta4/Delta2) = ^5  = {l2**5:.4f}   (Delta-ratio reading)")
    ln = Decimal(48).ln() / l2.ln()
    ln2 = probe.ln() / l2.ln()
    print(f"  k window: lambda_2^k in [{probe}, 48]  <=>  k in [{ln2:.4f}, {ln:.4f}]"
          f"   (contains no integer; truth needs k >= {ln2:.3f}, target needs k <= {ln:.3f})")
    print(f"  slack form:  lambda_* <= kappa * lambda_2^2  closes iff kappa <= 48/lambda_2^2"
          f" = {Decimal(48)/(l2**2):.4f}")
    print(f"               observed kappa (probe) = {probe/(l2**2):.4f}")
    print(f"  L_t(20) = {lt20:.6f}   L_t(18) = {lt18:.6f}   (lc4 numbers, 48 between them)")
    print("\n  slice-set count (why lambda_2 cannot be used for a slice):")
    print(f"  {'k':>3} {'#k-subsets of Z^2 up to transl. >= C(k^2,k-1)':>46} {'k-th root':>10}")
    from math import comb
    for k in range(2, 13, 2):
        q = comb(k * k, k - 1)
        print(f"  {k:>3} {q:>46} {q ** (1.0 / k):>10.3f}")
    print("  (grows like e*k -> unbounded per site: disconnected slice contents are NOT")
    print("   controlled by lambda_2; only connected slices would be.)\n")


if __name__ == "__main__":
    maxn = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    s0()
    s1()
    s2(maxn)
    s3()
