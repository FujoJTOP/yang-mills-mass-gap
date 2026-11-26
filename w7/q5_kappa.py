"""Q5: pin kappa = log mu for edge-connected plaquette animals in Z^4.

(a) brute-force: degree of one plaquette (edge-adjacency) in Z^d, d=2..5;
    lattice coordination 2d; check the formula 4*(2d-3).
(b) exact enumeration of plaquette animals in Z^4, A <= AMAX (translation classes),
    report a_A and a_A^(1/A)  [E];  cross-check the single-plane sub-family = Z^2
    site animals (1,2,6,19,63,216).
(c) kappa bounds:  kmin = log 4.06 (q2, Z^2 sub-family), kmax = log(4*(Delta-1))
    with Delta = plaquette degree; also log(4*Delta) for q2's own encoding.
Only prints; writes nothing but the redirected log.
"""
import sys, time
from itertools import combinations, product

D = 4
AX = tuple(range(D))
PLANES = list(combinations(AX, 2))          # 6 coordinate 2-planes


def verts(pl):
    """vertex set (4 lattice vertices) of plaquette (origin, plane-index)."""
    v, p, = pl
    a, b = PLANES[p]
    out = [v]
    for s in (0, 1):
        for t in (0, 1):
            w = list(v); w[a] += s; w[b] += t
            out.append(tuple(w))
    return frozenset(out)


def canon(pls):
    """translate a set of plaquettes so that its lex-min VERTEX sits at origin."""
    lo = min(min(v for v in verts(p)) for p in pls)
    d = tuple(-c for c in lo)
    return frozenset((tuple(x + y for x, y in zip(p[0], d)), p[1]) for p in pls)


def neighbours(p):
    """all plaquettes sharing an EDGE (=2 vertices) with p; brute force, no formula."""
    V = verts(p)
    out = []
    for qi in range(len(PLANES)):
        for off in product((-1, 0, 1), repeat=D):
            w = tuple(p[0][i] + off[i] for i in range(D))
            cand = (w, qi)
            if cand == p:
                continue
            if len(V & verts(cand)) == 2:
                out.append(cand)
    return out


def main():
    t0 = time.time()
    print("== (1) exact adjacency counts, brute force ==")
    for d in range(2, 6):
        global D, AX, PLANES
        Dsave, AXsave, PLANEsave = D, AX, PLANES
        D = d; AX = tuple(range(D)); PLANES = list(combinations(AX, 2))
        base = (tuple(0 for _ in range(D)), 0)
        nb = set(neighbours(base))
        print(f"  d={d}: lattice coord 2d={2*d:2d} | plaquettes sharing an edge with one plaquette"
              f" = {len(nb):2d} | 4*(2d-3)={4*(2*d-3):2d} | match={len(nb)==4*(2*d-3)}"
              f" | neighbours per edge (incl. self) = {len(nb)//4+1}")
        D, AX, PLANES = Dsave, AXsave, PLANEsave
    print("  -> 3D value is 12, 4D value is 20 (q2's Delta=12 is the d=3 number)")

    print("\n== (2) plaquette animals in Z^4, exact, translation classes ==")
    AMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    cap = float(sys.argv[2]) if len(sys.argv) > 2 else 150.0
    root = frozenset({(tuple(0 for _ in range(D)), 0)})
    cur = {canon(root)}
    print(f"  A=1  a_A=1  a_A^(1/A)=1.000000  ({time.time()-t0:.1f}s)")
    for A in range(2, AMAX + 1):
        nxt = set()
        for S in cur:
            for p in S:
                for q in neighbours(p):
                    if q not in S:
                        nxt.add(canon(S | {q}))
        cur = nxt
        a = len(cur)
        print(f"  A={A}  a_A={a}  a_A^(1/A)={a**(1.0/A):.6f}  ({time.time()-t0:.1f}s)")
        if time.time() - t0 > cap:
            print(f"  [time cap {cap}s hit, stop]")
            break

    print("\n== (3) cross-check: single 2-plane sub-family == Z^2 site animals ==")
    # plaquettes in plane (0,1) with all other coords 0: adjacency within plane = Z^2 site adj.
    known = {1: 1, 2: 2, 3: 6, 4: 19, 5: 63, 6: 216, 7: 760}
    curz = {canon(root)}
    for A in range(1, 8):
        if A > 1:
            nxt = set()
            for S in curz:
                for p in S:
                    for q in neighbours(p):
                        if q not in S and q[1] == 0:
                            nxt.add(canon(S | {q}))
            curz = nxt
        if A in known:
            print(f"  n={A}: enumerated in-plane={len(curz)}  Z^2 site animals={known[A]}"
                  f"  match={len(curz)==known[A]}  n^(1/n)={known[A]**(1.0/A):.4f}")

    print("\n== (4) kappa bounds ==")
    import math
    Dlt = 4 * (2 * D - 3)
    print(f"  Delta(plaquette degree, d=4) = {Dlt}")
    print(f"  lower: mu >= 4.06      -> kappa >= {math.log(4.06):.6f}")
    print(f"  upper(q2 form 4D)   : mu <= 4*Delta     = {4*Dlt:5d} -> kappa <= {math.log(4*Dlt):.6f}")
    print(f"  upper(refined, no parent) : mu <= 4*(Delta-1) = {4*(Dlt-1):5d} -> kappa <= {math.log(4*(Dlt-1)):.6f}")
    print(f"  q2 as printed (Delta=12)  : mu <= {4*12:5d} -> kappa <= {math.log(48):.6f}")
    for name, v in (("4*100", 400),):
        pass
    print(f"  widths: q2={math.log(48)-math.log(4.06):.3f}  corrected={math.log(4*(Dlt-1))-math.log(4.06):.3f}")


main()
