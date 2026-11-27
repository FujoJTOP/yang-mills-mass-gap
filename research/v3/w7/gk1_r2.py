# gk1_r2.py -- GK1-R2: constant family of R2 (P3A:154): tree layer has no K3 => each link carries <=2 plaquettes
# py -3.12 -I gk1_r2.py     (ASCII-only output on purpose: console codepage mangling)
import sys, time
from fractions import Fraction
from math import comb

T0 = time.time()
D = 4
L = 12                      # torus side; for n<=5 coordinate span <=4  => no wrap-around fake adjacency
PLANES = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
PIDX = {pl:i for i,pl in enumerate(PLANES)}
P6 = 6
NP = P6 * L**4

def pid(x, pl):
    return (((x[0]*L + x[1])*L + x[2])*L + x[3])*P6 + PIDX[pl]

def unpid(p):
    pl = PLANES[p % P6]; t = p // P6
    x3 = t % L; t //= L
    x2 = t % L; t //= L
    x1 = t % L; t //= L
    return (t, x1, x2, x3), pl

def ad(x, v): return tuple((x[i]+v[i]) % L for i in range(4))
E4 = [tuple(1 if k == i else 0 for k in range(4)) for i in range(4)]

def links_of(x, pl):
    mu, nu = pl
    return [(x, mu), (ad(x, E4[mu]), nu), (ad(x, E4[nu]), mu), (x, nu)]

from collections import defaultdict
l2p = defaultdict(list)
for i in range(NP):
    x, pl = unpid(i)
    for lk in links_of(x, pl):
        l2p[lk].append(i)

print("== A local structure ==")
print("A1 #plaquettes             =", NP)
print("A2 #links                  =", len(l2p))
print("A3 plaquettes-per-link set =", sorted(set(len(v) for v in l2p.values())))
nbr = {}
for i in range(NP):
    x, pl = unpid(i)
    s = set()
    for lk in links_of(x, pl):
        for q in l2p[lk]:
            if q != i: s.add(q)
    nbr[i] = sorted(s)
print("A4 plaquette-degree set    =", sorted(set(len(v) for v in nbr.values())))

p0 = pid((0,0,0,0), (0,1))
x0, pl0 = unpid(p0)
groups = [sorted(q for q in l2p[lk] if q != p0) for lk in links_of(x0, pl0)]
print("A5 link-groups of p0 sizes =", [len(g) for g in groups])
allg = [q for g in groups for q in g]
print("A6 groups partition N(p0)  =", len(set(allg)) == 20 and sorted(allg) == nbr[p0])

def shared_links(a, b):
    xa, pla = unpid(a); xb, plb = unpid(b)
    return len(set(links_of(xa, pla)) & set(links_of(xb, plb)))
print("A7 #shared links p0-neigh  =", sorted(set(shared_links(p0, q) for q in nbr[p0])))
print("A8 group members pairwise adj (K5, +p0 = K6) =",
      all(shared_links(u, v) == 1 for g in groups for i1,u in enumerate(g) for v in g[i1+1:]))
t3 = [p0] + groups[0][:2]
print("A9 3 plaquettes on one link form K3 =",
      all(shared_links(u, v) >= 1 for i1,u in enumerate(t3) for v in t3[i1+1:]))

print("\n== B full enumeration n<=5 ==")
nbr_set = {i: set(nbr[i]) for i in range(NP)}
layers = [None, {(p0,)}]
for n in range(1, 5):
    nxt = set()
    for A in layers[n]:
        sA = set(A); cand = set()
        for p in A: cand |= (nbr_set[p] - sA)
        for q in cand: nxt.add(tuple(sorted(A + (q,))))
    layers.append(nxt)
    print("B1 |layer n=%d|             = %d   (%.1f s)" % (n+1, len(nxt), time.time()-T0))

def stats(A):
    cnt = defaultdict(int)
    for p in A:
        x, pl = unpid(p)
        for lk in links_of(x, pl): cnt[lk] += 1
    Eadj = sum(c*(c-1)//2 for c in cnt.values())
    return Eadj, max(cnt.values())

print("\nB2 k-layer (k=|E|-n+1) and link multiplicity m")
for n in range(1, 6):
    per = defaultdict(int); mper = defaultdict(int); lab = defaultdict(int)
    for A in layers[n]:
        Eadj, m = stats(A); k = Eadj - n + 1
        per[k] += 1; mper[m] += 1; lab[(k == 0, m <= 2)] += 1
    print("  n=%d a_n=%d" % (n, len(layers[n])))
    print("      k dist :", dict(sorted(per.items())))
    print("      m dist :", dict(sorted(mper.items())))
    print("      (k==0,m<=2):", dict(lab))

viol = 0; mm2_k1 = 0; tot_k0 = 0; tot_mm2 = 0
for n in range(1, 6):
    for A in layers[n]:
        Eadj, m = stats(A); k = Eadj - n + 1
        if k == 0:
            tot_k0 += 1
            if m > 2: viol += 1
        if m <= 2:
            tot_mm2 += 1
            if k >= 1: mm2_k1 += 1
print("\nB3a trees n<=5 = %d ; with m>2 (violations of link<=2) = %d" % (tot_k0, viol))
print("B3b m<=2 animals n<=5 = %d ; of them k>=1 = %d" % (tot_mm2, mm2_k1))

print("\nB4 per-layer: k=0 count vs m<=2 count")
for n in range(1, 6):
    c0 = 0; cm2 = 0; cg3 = 0
    for A in layers[n]:
        Eadj, m = stats(A)
        if Eadj - n + 1 == 0: c0 += 1
        if m <= 2: cm2 += 1
        if m >= 3: cg3 += 1
    print("  n=%d  k=0: %7d  m<=2: %7d (%6.3f%%)  m>=3: %7d" %
          (n, c0, cm2, 100.0*cm2/len(layers[n]), cg3))

print("\nB5 local ratios")
tot = [0] + [len(layers[n]) for n in range(1, 6)]
for n in range(2, 6):
    print("  n=%d->%d  total ratio %8.4f   a_n^(1/n) %8.4f" % (n-1, n, tot[n]/tot[n-1], tot[n]**(1.0/n)))

print("\n== C R2 constant family: two-type Lagrange ==")
def Gn(n):
    return Fraction(comb(3*n, n-1) * 5**(n-1), n) if n >= 1 else 0
Gc = [0] + [Gn(n) for n in range(1, 9)]
Fc = [0]*9
for n in range(1, 9):
    s = Gc[n]
    for i in range(1, n): s += 5 * Gc[i] * Gc[n-i]
    Fc[n] = s
def Fclosed(n):
    t = comb(3*n, n-1) * 5**(n-1)
    if n >= 2: t += 10 * comb(3*n, n-2) * 5**(n-2)
    return Fraction(t, n)
print("C1 [x^n]G =", [int(Gc[n]) for n in range(1, 7)])
print("C2 [x^n]F =", [int(Fc[n]) for n in range(1, 7)])
print("C3 closed form == G+5G^2 :", all(Fc[n] == Fclosed(n) for n in range(1, 9)))
a0 = {1:1, 2:20, 3:402, 4:8328, 5:177295}
print("C4 a_{n,0} vs F_n :")
for n in range(1, 6):
    print("     n=%d  a0=%7d  F_n=%7d  ratio=%.4f" % (n, a0[n], Fc[n], a0[n]/Fc[n]))
u = Fraction(1, 10)
xstar = u / (1 + 5*u)**3
print("C5 P(u)-uP'(u)=0 at u*=1/10 :", (1+5*u)**3 - u*3*5*(1+5*u)**2 == 0)
print("C6 x* = u*/(1+5u*)^3 =", xstar, " => constant 135/4 =", 1/xstar, "=", float(1/xstar))
print("C7 47.6270568122 / 33.75 =", 47.627056812155747/float(1/xstar))

print("\n== D three readings of the same constraint ==")
import mpmath as mp
mp.mp.dps = 30
def crit(P, dP):
    f = lambda t: P(t) - t*dP(t)
    lo, hi = mp.mpf('1e-7'), mp.mpf('0.9')
    for _ in range(300):
        mid = (lo+hi)/2
        if f(mid) > 0: lo = mid
        else: hi = mid
    uu = (lo+hi)/2
    return uu, P(uu)/uu
u1, k1 = crit(lambda t: (1+5*t)**3, lambda t: 15*(1+5*t)**2)
print("D1 group 3 links x 5 slots : u*=%s  const=%s" % (mp.nstr(u1, 12), mp.nstr(k1, 12)))
u2, k2 = crit(lambda t: (1+5*t)**4, lambda t: 20*(1+5*t)**3)
print("D2 group 4 links (no parent slot removed): u*=%s  const=%s" % (mp.nstr(u2, 12), mp.nstr(k2, 12)))
u3, k3 = crit(lambda t: 1+19*t+171*t**2+969*t**3, lambda t: 19+342*t+2907*t**2)
print("D3 19 independent slots, out-degree<=3   : u*=%s  const=%s" % (mp.nstr(u3, 12), mp.nstr(k3, 12)))
print("D4 L_t(19)=19^19/18^18 =", mp.nstr(mp.mpf(19)**19/mp.mpf(18)**18, 20))
print("D5 L_t(20)=20^20/19^19 =", mp.nstr(mp.mpf(20)**20/mp.mpf(19)**19, 20))
print("D6 threshold 47.627056812155747 (GJ1-UC4:203)")
for nm, v in [("D1", k1), ("D2", k2), ("D3", k3)]:
    print("   %s vs 47.627 : %s (margin %s x)" % (nm, "PASS" if v < mp.mpf('47.6270568122') else "FAIL",
                                                 mp.nstr(mp.mpf('47.6270568122')/v, 8)))
print("   L_t(19)/47.627 = %s ; L_t(19)/48 = %s" %
      (mp.nstr(mp.mpf(19)**19/mp.mpf(18)**18/mp.mpf('47.6270568122'), 8),
       mp.nstr(mp.mpf(19)**19/mp.mpf(18)**18/48, 8)))

print("\n== E dimension family kappa(d) = 6.75*(2d-3) ==")
for d in (2, 3, 4, 5):
    c = 2*(d-1) - 1
    uu = Fraction(1, c*2)
    kd = (1 + c*uu)**3 / uu
    print("  d=%d  c=2(d-1)-1=%d  u*=1/%d  kappa=%s = %.4f" % (d, c, c*2, kd, float(kd)))
print("  d=2 truth 4.0626 -> slack %.4f x" % (float(Fraction(27,4))/4.0626))
print("\ntotal %.1f s" % (time.time()-T0))
