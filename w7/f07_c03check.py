"""f07 -- adjudicate c03's Cbar-free bound.

Checks (ASCII output only):
  A  identity (L1-c.1) + geometric absorption, exact Fraction
  B  root-tree recursion / Lagrange r_n=(1/n)C(nD,n-1) vs (eD)^n
  C  lattice facts: plaquettes/link=6, adjacency degree D=20
  D  hand-check n<=5 (exhaustive) and n=6 (cube exact + random sample):
     |phi(Gamma)| <= tau(G_Gamma) <= (eD)^n ; A(alpha) vs bound
  E  weight ratio |w|/u  (cube / 4-wall tube: link incidence 2, opposite)
  F  windows: x = eD*u  for u in {c_f/N, c_f, Phi(0)}, at beta=1/12, N=2
"""
import sys, itertools, math, random
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 30
D = 20
eD = mp.e * D                      # 54.3656365692...
N = 2
beta = mp.mpf(1) / 12
cf = beta / (2 * N)                # c_f = beta/(2N)  = 0.0208333
u_norm = cf / N                    # X01 convention u := c_f/N = beta/(2N^2)
Phi0 = 2 * N**2 * u_norm           # = beta

L = []
def p(*a):
    s = " ".join(str(x) for x in a)
    print(s); L.append(s)

# ---------------------------------------------------------------- A
p("== A. identity (L1-c.1) and geometric absorption ==")
A = [Fraction(7, 10)**(k + 1) for k in range(6)]
phi = sum(A)
logD = sum((1 - k) * A[k] for k in range(6))
res = []
for T in range(6):
    Lft = sum((T - k + 1) * A[k] for k in range(T + 1))
    eps = sum((k - 1 - T) * A[k] for k in range(T + 1, 6))
    res.append(Lft - (phi * T + logD + eps))
p("  residuals T=0..5 :", [str(x) for x in res])
for x in (mp.mpf('0.57'), mp.mpf('0.566309')):
    worst = 0
    for T in range(9):
        tail = mp.nsum(lambda k: (k - 1 - T) * x**(T + k), [1, mp.inf]) if False else \
               sum((k - 1) * x**(T + k) for k in range(1, 400))
        worst = max(worst, abs(tail - x**(T + 2) / (1 - x)**2))
    p("  x=%s  max|sum_{a>T}(a-1-T)x^a - x^(T+2)/(1-x)^2| = %s" % (mp.nstr(x, 8), mp.nstr(worst, 5)))
p("  model A(a)=x^(a+1)/(1-x) => |eps_T|/x^T = x^3/(1-x)^3 ; x=0.5663 -> %s"
  % mp.nstr(mp.mpf('0.566309')**3 / (1 - mp.mpf('0.566309'))**3, 8))

# ---------------------------------------------------------------- B
p("== B. root-tree recursion R<=zu(1+R)^D, Lagrange r_n=(1/n)C(nD,n-1) ==")
p("   n   r_n            r_n^(1/n)/eD")
for n in (1, 2, 3, 4, 5, 6, 10, 20, 50, 200, 1000):
    rn = Fraction(1, n) * math.comb(D * n, n - 1)
    rat = (mp.mpf(rn.numerator) / rn.denominator) ** (mp.mpf(1) / n) / eD
    p("  %4d  %-14s %s" % (n, rn, mp.nstr(rat, 5)))
R = mp.findroot(lambda z: z - u_norm * (1 + z)**D, 0.01)
p("  fixed point R=u(1+R)^D at u=%s : R=%s  (u*eD=%s)"
  % (mp.nstr(u_norm, 6), mp.nstr(R, 6), mp.nstr(u_norm * eD, 6)))
p("  eD = %s ; Delta = %d" % (mp.nstr(eD, 12), D))

# ---------------------------------------------------------------- C
def addsite(x, mu, k=1):
    y = list(x); y[mu] += k; return tuple(y)
def norm(x, mu, nu):
    return (x, mu, nu) if mu < nu else (x, nu, mu)
def plaqs_at_link(z, mu):
    out = []
    for nu in range(4):
        if nu == mu:
            continue
        out.append(norm(z, mu, nu))
        out.append(norm(addsite(z, nu, -1), nu, mu))
    return out
def nbrs(pp):
    x, mu, nu = pp
    out = set()
    for (z, l) in ((x, mu), (x, nu), (addsite(x, mu), nu), (addsite(x, nu), mu)):
        for q in plaqs_at_link(z, l):
            if q != pp:
                out.add(q)
    return out

p("== C. lattice facts (open 4D, interior plaquette) ==")
p0 = ((0, 0, 0, 0), 0, 1)
nb = nbrs(p0)
p("  p0 = %s  degree = %d  (Delta=20 expected)" % (p0, len(nb)))
cnt = {}
for (z, l) in ((p0[0], p0[1]), (p0[0], p0[2]),
               (addsite(p0[0], p0[1]), p0[2]), (addsite(p0[0], p0[2]), p0[1])):
    for q in plaqs_at_link(z, l):
        cnt[q] = cnt.get(q, 0) + 1
p("  plaquettes at one link = %d (c05 [A])" % len(plaqs_at_link((0, 0, 0, 0), 0)))

# ---------------------------------------------------------------- helpers
def partitions(n):
    """all set partitions of range(n) as tuples of bitmasks"""
    out = []
    def rec(i, blocks):
        if i == n:
            out.append(tuple(blocks)); return
        for j in range(len(blocks)):
            blocks[j] |= 1 << i; rec(i + 1, blocks); blocks[j] ^= 1 << i
        blocks.append(1 << i); rec(i + 1, blocks); blocks.pop()
    rec(0, [])
    return out
PARTS = {n: partitions(n) for n in range(1, 7)}

def ursell(adjm, n):
    """phi(Gamma) with all activities = 1 :
       phi = sum_{A subset E} (-1)^|A| * 1[V(A) spanning and connected]
       (Fernandez-Procacci (2.4) form; brute force, |E|<=12 here)"""
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if (adjm[i] >> j) & 1]
    m = len(edges)
    tot = 0
    for mask in range(1 << m):
        par = list(range(n))
        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]; a = par[a]
            return a
        bits = 0
        for k, (i, j) in enumerate(edges):
            if (mask >> k) & 1:
                ra, rb = find(i), find(j)
                if ra != rb:
                    par[ra] = rb
        if len({find(v) for v in range(n)}) == 1:
            tot += (-1) ** bin(mask).count("1")
    return tot

def nspanning(adjm, n):
    """tau via Kirchhoff, reduced Laplacian determinant (exact int)"""
    if n == 1:
        return 1
    M = [[0] * (n - 1) for _ in range(n - 1)]
    for i in range(n - 1):
        deg = bin(adjm[i]).count("1")
        M[i][i] = deg
        for j in range(n - 1):
            if (adjm[i] >> j) & 1:
                M[i][j] -= 1
    det = 1
    for c in range(n - 1):
        piv = next((r for r in range(c, n - 1) if M[r][c]), None)
        if piv is None:
            return 0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]; det = -det
        det *= M[c][c]
        for r in range(c + 1, n - 1):
            f = Fraction(M[r][c], M[c][c])
            for k in range(c, n - 1):
                M[r][k] -= f * M[c][k]
    return abs(int(det))

def lo(p_):
    return p_[0][0]
def hi(p_):
    x, mu, nu = p_
    return x[0] + (1 if (mu == 0 or nu == 0) else 0)

def metrics(animal):
    vs = sorted(animal)
    n = len(vs)
    pos = {v: i for i, v in enumerate(vs)}
    adjm = [0] * n
    for v in vs:
        for w in nbrs(v):
            if w in pos:
                adjm[pos[v]] |= 1 << pos[w]
    return ursell(adjm, n), nspanning(adjm, n)

# ---------------------------------------------------------------- D
p("== D. hand-check: |phi| <= tau <= (eD)^n, and A(alpha) vs bound ==")
def enum(nmax):
    lvl = {frozenset([p0])}
    out = [lvl]
    for n in range(1, nmax):
        nxt = set()
        for A_ in lvl:
            bd = set()
            for q in A_:
                bd |= nbrs(q)
            bd -= A_
            for q in bd:
                nxt.add(A_ | {q})
        out.append(nxt); lvl = nxt
    return out

levels = enum(5)
p("  pinned animal counts n=1..5 :", [len(s) for s in levels],
  "(c03 [C]/c05 [B] log: 1,20,458,11132,280175)")

S = {}
TAU = {}
bad = 0
for n in range(1, 5):
    s = t = 0
    for A_ in levels[n - 1]:
        ph, ta = metrics(A_)
        if abs(ph) > ta:
            bad += 1
        s += abs(ph); t += ta
    S[n], TAU[n] = s, t
random.seed(7)
lv5 = list(levels[4]); random.shuffle(lv5)
s = t = 0
for A_ in lv5[:4000]:
    ph, ta = metrics(A_)
    if abs(ph) > ta:
        bad += 1
    s += abs(ph); t += ta
S[5], TAU[5] = s * len(lv5) / 4000, t * len(lv5) / 4000
p("   |phi|>tau violations : %d  (all levels incl. 4000-sample of n=5)" % bad)
p("   n | sum|phi| (pinned) | sum tau (=#trees) | (eD)^n | trees/(eD)^n")
for n in range(1, 6):
    tag = "" if n < 5 else "  [n=5 scaled from 4000 sample]"
    p("   %d | %-16d | %-16d | %-10s | %s%s"
      % (n, round(S[n]), round(TAU[n]), mp.nstr(eD**n, 8),
         mp.nstr(TAU[n] / mp.mpf(eD)**n, 5), tag))

p("  -- n=6 by hand (unit spatial cube, 6 faces) --")
x0 = (0, 0, 0, 0)
cube = frozenset([norm(x0, 1, 2), norm(x0, 1, 3), norm(x0, 2, 3),
                  norm(addsite(x0, 1), 2, 3), norm(addsite(x0, 2), 1, 3),
                  norm(addsite(x0, 3), 1, 2)])
ph6, ta6 = metrics(cube)
p("   cube: |phi|=%d  tau=%d  (tau(Q3)=384)  |phi|<=tau : %s"
  % (abs(ph6), ta6, abs(ph6) <= ta6))
# link incidence on the cube: each link in 2 plaquettes, opposite orientation
inc = {}
for q in cube:
    x, mu, nu = q
    for (z, l) in ((x, mu), (x, nu), (addsite(x, mu), nu), (addsite(x, nu), mu)):
        inc[(z, l)] = inc.get((z, l), 0) + 1
p("   cube link incidences : %s  (all 2 -> closed surface)"
  % (sorted(set(inc.values())),))

p("  -- A(alpha) vs c03 bound, pinned sum, u=%s --" % mp.nstr(u_norm, 6))
x0u = eD * u_norm
X_alpha = {}
for n in range(1, 5):
    for A_ in levels[n - 1]:
        ph, ta = metrics(A_)
        a = hi(max(A_, key=lambda q: hi(q))) - lo(min(A_, key=lambda q: lo(q)))
        X_alpha[a] = X_alpha.get(a, 0) + abs(ph) * u_norm**n
p("   alpha : |A_pinned(alpha)|  vs  eD*u^alpha/(1-eD*u)  (c03 bound /24)")
for a in sorted(X_alpha):
    bd = x0u**a / (1 - x0u)
    p("     %2d   %-12s  %-12s  ratio=%s"
      % (a, mp.nstr(X_alpha[a], 6), mp.nstr(bd, 6),
         mp.nstr(X_alpha[a] / bd, 4)))

# ---------------------------------------------------------------- E
p("== E. weight ratio |w|/u (why Cbar is not 1) ==")
p("   per plaquette: |w| <= c_f  (Z06-HOLEA:83, PROOF lem.2'); exact for cube/tube")
p("   c_f = %s ; u_norm=c_f/N = %s ; c_f/u_norm = %s = N"
  % (mp.nstr(cf, 8), mp.nstr(u_norm, 8), mp.nstr(cf / u_norm, 8)))
p("   Phi(0) = sum_r c_r dim r = beta = %s ; Phi(0)/u_norm = 2N^2 = %s"
  % (mp.nstr(Phi0, 8), mp.nstr(Phi0 / u_norm, 8)))
uu = mp.besseli(2, beta) / mp.besseli(1, beta)
p("   SU(2) exact u = I2/I1 = %s  (c05 [E]: 0.020827)" % mp.nstr(uu, 8))

# ---------------------------------------------------------------- F
p("== F. windows (beta < 2N^2/(eD*Cbar)), N=2, eD=%s ==" % mp.nstr(eD, 12))
rows = [("c03 : x=eD*u_norm        ", eD * u_norm, mp.mpf(1)),
        ("f07 : x=eD*c_f (exact |w|)", eD * cf, mp.mpf(N)),
        ("    x=eD*Phi(0)=eD*beta  ", eD * Phi0, mp.mpf(2 * N**2)),
        ("c05 : x=4*eD*Phi(0)      ", 4 * eD * Phi0, mp.mpf(4 * 2 * N**2)),
        ("a01 : x=eD*20*u_norm     ", eD * 20 * u_norm, mp.mpf(20))]
for name, x, ce in rows:
    win = 2 * N**2 * u_norm / x        # = beta at which U=1
    p("  %s x=%-10s win_beta=%-12s = %-8s x(1/48)   beta=1/12 in window: %s"
      % (name, mp.nstr(x, 7), mp.nstr(win, 7),
         mp.nstr(win * 48, 6), win > beta))
p("  B(1/12)=24*x^2/(1-x)^3 at c03 x : %s ; 24 is N_piv=4 walls*6 plaq/link"
  % mp.nstr(24 * (eD * u_norm)**2 / (1 - eD * u_norm)**3, 8))

# convention-free window: mean plaquette < 1/eD
I1b = lambda b: (mp.besseli(0, b) + mp.besseli(2, b)) / (2 * mp.besseli(1, b)) - 1 / b
bst = mp.findroot(lambda b: I1b(b) - 1 / eD, mp.mpf('0.05'))
p("  exact SU(2) window (convention-free): <(1/N)tr U_p> < 1/eD  <=>  beta < %s"
  % mp.nstr(bst, 8))
p("     at beta=1/12 : <(1/N)tr U_p> = %s , *eD = %s  ->  U<1 : %s"
  % (mp.nstr(I1b(beta), 8), mp.nstr(I1b(beta) * eD, 8), I1b(beta) * eD < 1))
p("     beta*=%.8f = %s x (1/48) ; beta=1/12 = 0.0833333 > beta* : %s"
  % (float(bst), mp.nstr(bst * 48, 6), bst < beta))

open("w7/f07_c03check.log", "w", encoding="utf-8").write("\n".join(L) + "\n")
