#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GE2-DEFDIM witness (single job, pure python).

Part 1  BL-3 slab on T^4_L with exact Z_N centre-valued plaquettes:
        costs, and the closure of {S, upper layer, lower layer, D = defect, full-volume stack}.
Part 2  codim-2 thin centre vortex ('t Hooft) membrane.
Part 3  fixed polyomino (cell-animal) enumeration -> planar class growth constant (hard bound).
Part 4  bargraphs (flat-bottom height sequences) -> mu = 2 exactly.
Part 5  window arithmetic 0.5 + gamma + ln3  vs  2.598612; per-plaquette rates vs 3.73-3.84.
"""
import math

def plaq_exps(L, N, twist):
    def k(x, mu):
        return twist.get((x[0] % L, x[1] % L, x[2] % L, x[3] % L, mu), 0)
    E = {}
    for x0 in range(L):
        for x1 in range(L):
            for x2 in range(L):
                for x3 in range(L):
                    x = (x0, x1, x2, x3)
                    for mu in range(4):
                        for nu in range(mu + 1, 4):
                            a = list(x); a[mu] = (a[mu] + 1) % L
                            b = list(x); b[nu] = (b[nu] + 1) % L
                            E[(mu, nu, x)] = (k(x, mu) + k(tuple(a), nu)
                                              - k(tuple(b), mu) - k(x, nu)) % N
    return E

def cu(exp, N):
    return 1.0 - math.cos(2.0 * math.pi * exp / N)

def boundary(chain, L):
    """oriented boundary of the 2-chain {(mu,nu,base) with mu<nu}, all coefficient +1."""
    b = {}
    def add(x, idx, s):
        key = ((x[0] % L, x[1] % L, x[2] % L, x[3] % L), idx)
        b[key] = b.get(key, 0) + s
    for (mu, nu, x) in chain:
        m = list(x); m[mu] += 1
        n = list(x); n[nu] += 1
        add(x, mu, +1); add(tuple(m), nu, +1); add(tuple(n), mu, -1); add(x, nu, -1)
    return {e: v for e, v in b.items() if v != 0}

def cost_of(chain, E, N):
    return sum(cu(E[p], N) for p in chain)

print("=" * 82)
print("PART 1 -- BL-3 slab: U_l = sigma0 on the 3 spatial dirs of the slice x3=0, on T^4_L")
print("=" * 82)
for L in (2, 3, 4):
    for N in (2, 3):
        tw = {}
        for a in range(L):
            for b in range(L):
                for c in range(L):
                    for mu in (0, 1, 2):
                        tw[(a, b, c, 0, mu)] = 1
        E = plaq_exps(L, N, tw)
        total = sum(cu(s, N) for s in E.values())
        site = lambda a, b, c, d: (a % L, b % L, c % L, d % L)
        S = [(mu, nu, site(a, b, c, 0)) for (mu, nu) in ((0, 1), (0, 2), (1, 2))
             for a in range(L) for b in range(L) for c in range(L)]
        up = [(mu, 3, site(a, b, c, 0)) for mu in (0, 1, 2)
              for a in range(L) for b in range(L) for c in range(L)]
        lo = [(mu, 3, site(a, b, c, L - 1)) for mu in (0, 1, 2)
              for a in range(L) for b in range(L) for c in range(L)]
        D = [p for p in E if E[p] != 0]
        print(f"\n L={L} N={N}  beta=2.5   S_W = {2.5*total:.4f}   W/beta = {total:.4f}"
              f"   (rate 1-cos(2pi/N) = {1-math.cos(2*math.pi/N):.4f})")
        print(f"   |S|={len(S)} (3L^3={3*L**3})   cost(S)/beta = {cost_of(S,E,N):.6f}")
        print(f"   closure:  S {len(boundary(S,L))==0} | upper-layer {len(boundary(up,L))==0}"
              f" | lower-layer {len(boundary(lo,L))==0} | D=upp+low {len(boundary(D,L))==0}"
              f"  (D bad edges: {len(boundary(D,L))})")
        print(f"   |D| = {len(D)} (6L^3 = {6*L**3}) ; cost(D)/beta = {cost_of(D,E,N):.6f}"
              f" ; per-plaquette = {cost_of(D,E,N)/len(D):.6f}"
              f"   ; S_W/|S| = {total/len(S):.6f}")

    # twist on *every* spatial link: pure-gauge check
    tw = {}
    for a in range(L):
        for b in range(L):
            for c in range(L):
                for d in range(L):
                    for mu in (0, 1, 2):
                        tw[(a, b, c, d, mu)] = 1
    E = plaq_exps(L, N, tw)
    D = [p for p in E if E[p] != 0]
    tot = sum(cu(s, N) for s in E.values())
    print(f"   [uniform twist, L={L} N={N}] |D| = {len(D)}  S_W/beta = {tot:.1f}"
          f"  => uniform twist on spatial links is PURE GAUGE (all plaquettes trivial)")

print()
print("=" * 82)
print("PART 2 -- thin codim-2 centre vortex: U_3 = sigma0 on the links with x0=x1=0")
print("=" * 82)
for L in (2, 3):
    for N in (2, 3):
        tw = {}
        for c in range(L):
            for d in range(L):
                tw[(0, 0, c, d, 3)] = 1
        E = plaq_exps(L, N, tw)
        D = [p for p in E if E[p] != 0]
        tot = sum(cu(s, N) for s in E.values())
        print(f" L={L} N={N}: |D| = {len(D)} (4L^2 = {4*L*L})  per-plaquette = "
              f"{cost_of(D,E,N)/len(D):.6f}  (=1-cos(2pi/N) {1-math.cos(2*math.pi/N):.6f})  "
              f"S_W = {2.5*tot:.3f}  closed {len(boundary(D,L))==0}")

print()
print("   comparison at beta=2.5, both rate = 1-cos(2pi/N):")
for N in (2, 3, 4, 5, 6):
    r = 1 - math.cos(2 * math.pi / N)
    print(f"     N={N}:  rate = {2.5*r:.6f} per plaquette;  slab |C| = 6L^3, membrane |C| = 4L^2 ->"
          f" cost ratio = (6L^3)/(4L^2) = 1.5L")

print()
print("=" * 82)
print("PART 3 -- fixed polyominoes (planar class)")
print("=" * 82)
def poly_counts(nmax):
    cur = {((0, 0),)}
    out = [0] * (nmax + 1)
    out[1] = 1
    for n in range(2, nmax + 1):
        nxt = set()
        for Pt in cur:
            P = set(Pt)
            for (x, y) in P:
                for (dx, dy) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    c = (x + dx, y + dy)
                    if c in P:
                        continue
                    Q = P | {c}
                    mx = min(q[0] for q in Q); my = min(q[1] for q in Q)
                    nxt.add(tuple(sorted((q[0] - mx, q[1] - my) for q in Q)))
        cur = nxt
        out[n] = len(cur)
    return out

NMAX = 11
A = poly_counts(NMAX)
print("   n       : " + "".join(f"{n:>8d}" for n in range(1, NMAX + 1)))
print("   A_n     : " + "".join(f"{A[n]:>8d}" for n in range(1, NMAX + 1)))
print("   lnA/n   : " + "".join(f"{math.log(A[n])/n:>8.4f}" for n in range(1, NMAX + 1)))
print("   ratio   : " + "".join(f"{A[n]/A[n-1]:>8.4f}" for n in range(2, NMAX + 1)))
g = math.log(A[NMAX]) / NMAX
print(f"   hard bound (Fekete, supermultiplicative): gamma_surf >= {g:.6f} -> mu_surf >= {math.exp(g):.4f}")

print()
print("=" * 82)
print("PART 4 -- bargraphs: height seq (h_1..h_k), h_i>=1, sum n  => 2^(n-1) shapes")
print("=" * 82)
print("   n       : " + "".join(f"{n:>8d}" for n in range(1, 13)))
print("   B_n     : " + "".join(f"{2**(n-1):>8d}" for n in range(1, 13)))
print("   lnB/n   : " + "".join(f"{math.log(2**(n-1))/n:>8.4f}" for n in range(1, 13)))
print(f"   limit: gamma = ln2 = {math.log(2):.6f};  mu = 2 exactly")

print()
print("=" * 82)
print("PART 5 -- window 0.5+gamma+ln3 <= 2.598612 ?   and rate vs 3.73-3.84")
print("=" * 82)
TGT = 0.5 + 1.0 + math.log(3.0)
print(f"   target (gamma = 1) = {TGT:.6f}")
def root3(c0, c1, c2, c3):
    lo, hi = 1e-6, 1.0
    f = lambda x: c0 + c1 * x + c2 * x * x + c3 * x ** 3
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)
mu_rc = 1.0 / root3(1.0, -5.0, 7.0, -4.0)
rows = [("animals a_n rate 2.632032", math.exp(2.632032)),
        ("animals B_n rate 3.458607", math.exp(3.458607)),
        ("planar polyominoes n<=11 (hard)", math.exp(g)),
        ("planar polyominoes lit. 4.0626 [M]", 4.0626),
        ("row-convex (GF denom.) [M]", mu_rc),
        ("bargraphs (exact)", 2.0),
        ("fixed shapes (translations)", 1.0)]
for name, mu in rows:
    gg = math.log(mu); w = 0.5 + gg + math.log(3.0)
    v = "PASS" if w <= TGT + 1e-12 else f"FAIL {w-TGT:+.4f}"
    print(f"   {name:34s} mu={mu:10.5f} gamma={gg:8.5f} window={w:8.6f}  {v:12s}"
          f" mu<=e? {'yes' if mu <= math.e else 'no x%.2f' % (mu / math.e)}")
print()
print("   thinnest codim-2 centre vortex vs the required per-plaquette rate 3.73-3.84 (beta=2.5):")
for N in (2, 3, 4, 5, 6):
    c = 2.5 * (1 - math.cos(2 * math.pi / N))
    tag = "passes" if c >= 3.84 else ("MARGINAL" if c >= 3.73 else "fails")
    print(f"     SU({N}): {c:6.4f}   {tag}")
print()
print("done.")
