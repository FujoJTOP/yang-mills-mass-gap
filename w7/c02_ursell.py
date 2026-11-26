# -*- coding: utf-8 -*-
"""
c02_ursell.py -- sharpen the Ursell constant Cbar to <= 3.5.

Setting (verbatim, research/v3/species/X01-L1C.md:71):
    |phi^U(Gamma)| <= Cbar^{|Gamma|} u^{|Gamma|} ,  ubar := e*Delta*Cbar*u ,  Delta = 20
    hole-C: ubar -> Chat*ubar ,  Chat = 2 (exact) / 4 (per-plaquette envelope)
    no-shrink threshold:  Chat*Cbar < 2N^2*48/(e*Delta) = 7.0632853...   (N=2)

New content of this run (beyond a01 and c01):
  [1] graph class of phi^U + the BLOCK FACTORISATION (block = 2-connected component):
      phi^U(G) = (-1)^{#bridges} prod_B mu(B) ,  mu(B)=sum_{A conn. spanning in B}(-1)^{|A|}
      => factorial growth is capped by the clique number 6 (a link carries exactly 6 plaquettes)
  [2] exact |phi^U| for the clique family: (n-1)! for n<=6, 120^{1/6}=2.2209 = observed per-cluster sup
  [3] NO-GO for the tau route: the periodic 4D box (96 plaquettes, 20-regular) has tau^{1/96}>>3.5,
      so *no* sharpening of Kirchhoff/AM-GM can reach 3.5 for a per-cluster Cbar via tau
  [4] where the <=20 really breaks: position entropy (e*Delta) is already paid inside ubar
  [5] constants table: Chat*Cbar vs 7.0633, beta1^fix vs 1/48
  [6] verdict flags (a)/(b)/(c)

Run: py -3.12 w7/c02_ursell.py > w7/c02_ursell.log
"""
import time
import mpmath as mp

mp.mp.dps = 30
E = mp.e
DELTA = 20
LAM = E * DELTA
N = 2
CHAT_LOC, CHAT_ENV = mp.mpf(2), mp.mpf(4)
K48 = mp.mpf(1) / 48
THR = 2 * N * N * 48 / LAM


def hdr(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


# ---------------------------------------------------------------- plaquette graph
def plus(a, d):
    return tuple(a[i] + (1 if i == d else 0) for i in range(4))


def links_of(p, per=None):
    x, mu, nu = p
    if per is None:
        return frozenset([(x, mu), (x, nu), (plus(x, mu), nu), (plus(x, nu), mu)])
    mod = lambda s: tuple(c % per for c in s)
    return frozenset([(mod(x), mu), (mod(x), nu), (mod(plus(x, mu)), nu), (mod(plus(x, nu)), mu)])


def nbrs(p, per=None):
    """all plaquettes sharing >=1 link with p (exactly 20 in the infinite lattice)"""
    out = set()
    for a, d in links_of(p, per):
        for dp in range(4):
            if dp == d:
                continue
            lo, hi = min(d, dp), max(d, dp)
            minus = tuple(a[i] - (1 if i == dp else 0) for i in range(4))
            for c in (a, minus):                      # the two sides of the shared link
                q = (c if per is None else tuple(t % per for t in c), lo, hi)
                if q != p and (links_of(q, per) & links_of(p, per)):
                    out.add(q)
    out.discard(p)
    return out


def canon(cl):
    mn = min(p[0] for p in cl)
    return frozenset((tuple(p[0][i] - mn[i] for i in range(4)), p[1], p[2]) for p in cl)


def det_int(M):
    M = [row[:] for row in M]
    n, sign, prev = len(M), 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for i in range(k + 1, n):
                if M[i][k] != 0:
                    M[k], M[i] = M[i], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
        for i in range(k + 1, n):
            M[i][k] = 0
    return sign * M[n - 1][n - 1]


def adj_edges(cl, per=None):
    vs = list(cl)
    idx = {v: i for i, v in enumerate(vs)}
    es = set()
    for v in vs:
        for w in nbrs(v, per):
            j = idx.get(w)
            if j is not None and idx[v] < j:
                es.add((idx[v], j))
    return len(vs), sorted(es)


def tau_of(cl, per=None):
    n, es = adj_edges(cl, per)
    if n == 1:
        return 1
    L = [[0] * n for _ in range(n)]
    for i, j in es:
        L[i][i] += 1
        L[j][j] += 1
        L[i][j] -= 1
        L[j][i] -= 1
    return det_int([row[:n - 1] for row in L[:n - 1]])


def phi_ursell(cl, per=None):
    """exact phi^U(G) = sum_{A subset E, (V,A) connected} (-1)^|A|  (bitmask, m<=24)"""
    n, es = adj_edges(cl, per)
    m = len(es)
    if m > 24:
        return None
    par = list(range(n))

    def find(a):
        while par[a] != a:
            par[a] = par[par[a]]
            a = par[a]
        return a

    tot = 0
    for mask in range(1 << m):
        par[:] = list(range(n))
        k = bin(mask).count("1")
        for b in range(m):
            if mask >> b & 1:
                i, j = es[b]
                ri, rj = find(i), find(j)
                if ri != rj:
                    par[ri] = rj
        if len({find(i) for i in range(n)}) == 1:
            tot += -1 if k % 2 else 1
    return tot


# ------------------------------------------------------------------------------ [1]
hdr("[1] WHAT phi^U EXPANDS INTO: spanning trees, and how the block product sharpens it")
print("phi^U(Gamma) = Ursell / Mayer coefficient of the plaquette family Gamma:")
print("   phi^U(G) = sum_{A subset E(G_Gamma), (Gamma,A) connected} (-1)^{|A|}                (1.1)")
print("  (i)  cumulant reading: partition-lattice Moebius sum over connected multigraphs, signs +-")
print("  (ii) Penrose / Brydges-Kennedy tree interpolation: reorganised as a sum over SPANNING")
print("       TREES of the intersection graph G_Gamma,  |phi^U| <= tau(G_Gamma)   (Penrose)   (1.2)")
print("  (iii) BLOCK FACTORISATION (new here): (1.1) is multiplicative over 2-connected blocks,")
print("       because (V,A) is connected  <=>  A contains every bridge AND A restricted to each")
print("       block is connected & spanning.  Hence, with mu(B) := sum_{A c.s. in B}(-1)^{|A|},")
print("         phi^U(G) = (-1)^{#bridges} * prod_{blocks B} mu(B)                       (1.3)")
print("       Consequences:  tree / single edge  -> |mu|=1 ;  K_n -> |mu|=(n-1)! ;  C_n -> n-1.")
print("       A link carries exactly 6 plaquettes  =>  clique number of the plaquette graph = 6")
print("       => the factorial growth of |mu| STOPS at n=6.  This is the structural reason a")
print("          per-cluster Cbar cannot run away, and why tau >> |phi^U| on dense clusters.  [T]")

# ------------------------------------------------------------------------------ [2]
hdr("[2] EXACT |phi^U| AND tau FOR THE CLIQUE FAMILY  (bitmask enumeration, exact integers)")
print(" the k plaquettes through ONE LINK form K_k (k<=6) -- the only way to get a clique")
x0 = (0, 0, 0, 0)
LINK6 = ([(x0, 0, mu) for mu in (1, 2, 3)]
         + [((0, -1 * (mu == 1), -1 * (mu == 2), -1 * (mu == 3)), 0, mu) for mu in (1, 2, 3)])
for k in range(1, 7):
    cl = frozenset(LINK6[:k])
    n, es = adj_edges(cl)
    ph = phi_ursell(cl)
    t = tau_of(cl)
    print("  n=%d |E|=%2d  phi^U=%-6d tau=%-6d  |phi^U|^{1/n}=%-9s tau^{1/n}=%s"
          % (n, len(es), ph, t, mp.nstr(mp.power(abs(ph), mp.mpf(1) / n), 8),
             mp.nstr(mp.power(t, mp.mpf(1) / n), 8)))
print("  => per-cluster Cbar >= (n-1)!^{1/n} up to n=6:  1, 1, 1.2599, 1.5651, 1.8882, 2.2209")
print("     (K_6 exact: phi^U=-120, tau=1296 -> the tau bound over-counts by 10.8x at n=6)")
print("     and it STOPS: no K_n for n>=7 (clique number 6)  => per-cluster sup is ~2.22, not factorial.")
corner6 = frozenset([(x0, mu, nu) for mu, nu in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]])
n6, e6 = adj_edges(corner6)
print("     control (the 6 plaquettes at a corner are NOT a clique): |E|=%d phi^U=%d tau=%d"
      % (len(e6), phi_ursell(corner6), tau_of(corner6)))

# ------------------------------------------------------------------------------ [3]
hdr("[3] NO-GO FOR THE tau ROUTE: the periodic 4D box is 20-regular and tau^{1/96} >> 3.5")
per = 2
sites = [(a, b, c, d) for a in range(per) for b in range(per) for c in range(per) for d in range(per)]
box = frozenset([(x, mu, nu) for x in sites for mu in range(4) for nu in range(mu + 1, 4)])
nb = len(box)
degs = [len([w for w in nbrs(v, per) if w in box]) for v in box]
nedges = sum(degs) // 2
deg = sum(degs) // nb
t0 = time.time()
tau_box = tau_of(box, per)
rbox = mp.power(mp.mpf(tau_box), mp.mpf(1) / nb)
print(" L=2 periodic 4D box: n=%d plaquettes, 4*16=64 links, |E|(simple)=%d, degree=%d..%d"
      % (nb, nedges, min(degs), max(degs)))
print("   (the infinite lattice would give 20*96/2=960 pairs; on the L=2 torus some pairs share")
print("    TWO links and collapse in the simple graph: %d pairs lost, deg=%d regular)"
      % (960 - nedges, deg))
print(" tau = %s   (exact integer; Bareiss, %.1f s)" % (mp.nstr(mp.mpf(tau_box), 24), time.time() - t0))
print(" tau^{1/n} = %s      (AM-GM envelope Delta = 20)" % mp.nstr(rbox, 10))
print(" => sup_Gamma tau(G_Gamma)^{1/n} >= %s >> 3.5." % mp.nstr(rbox, 6))
print("    Therefore: NO refinement of Kirchhoff/AM-GM can ever deliver a PER-CLUSTER")
print("    Cbar <= 3.5 through (1.2). The looseness is not in the algebra, it is in using tau")
print("    (i.e. position entropy) as the counting device: tau pays for positions again.   [T] no-go")

# ------------------------------------------------------------------------------ [4]
hdr("[4] WHERE THE <=20 BREAKS: position entropy e*Delta is already charged inside ubar")
print(" ubar := e*Delta*Cbar*u = (KP local sum). Its e*Delta factor is the positioned-tree count:")
print("   # rooted trees at p with n vertices, max degree Delta  <=  (e*Delta)^n  =  Delta^n * e^n")
print("   Delta^n  = positions  (already paid by e*Delta in ubar)")
print("   e^n      = labelled-tree entropy / n!  (Cayley + Stirling) -- see the table below:")
print("      n      (n^{n-2})^{1/n}   (n^{n-2}/n!)^{1/n}")
for n in [5, 20, 100, 500, 1000]:
    c = mp.power(n, n - 2)
    print("    %-6d %-16s %s" % (n, mp.nstr(mp.power(c, mp.mpf(1) / n), 8),
                                mp.nstr(mp.power(c / mp.factorial(n), mp.mpf(1) / n), 8)))
print("      inf    n/e -> inf        e = %s" % mp.nstr(E, 8))
print(" => the correct split is  Cbar <= e = 2.71828  (shape entropy only), NOT Delta = 20;")
print("    charging Cbar=Delta pays the SAME position entropy twice. a01's 20 is an upper")
print("    bound for the wrong object (per-cluster tau^{1/n}), not for the Cbar in ubar.      [T]")
OTTER = mp.mpf("2.9557658454346846")
print(" pessimistic variant (count rooted unlabelled SHAPES, A000081, growth = Otter constant):")
A81 = [1, 1, 2, 4, 9, 20, 48, 115, 286, 719]
print("   n   :" + "".join("%8d" % (i + 1) for i in range(8)))
print("   a_n :" + "".join("%8d" % A81[i] for i in range(8)))
print("   ^1/n:" + "".join("%8s" % mp.nstr(mp.power(mp.mpf(A81[i]), mp.mpf(1) / (i + 1)), 5) for i in range(8)))
print("   Otter = %s  < 3.5  =>  the bound is robust to which of the two splits one uses."
      % mp.nstr(OTTER, 8))

# ------------------------------------------------------------------------------ [5]
hdr("[5] CONSTANTS: Chat*Cbar vs 7.0632853 , beta1^fix vs 1/48")
print(" threshold:  Chat*Cbar < 2N^2*48/(e*Delta) = %s" % mp.nstr(THR, 10))
print(" allowed:    Chat=4 -> Cbar < %s ;  Chat=2 -> Cbar < %s"
      % (mp.nstr(THR / 4, 8), mp.nstr(THR / 2, 8)))
print()
print(" Cbar                       Chat  Chat*Cbar  <7.0633?  beta1^fix=2N^2/(Chat eD Cbar)  /(1/48)")
for name, cb in [("1 (c01 tree identity)", mp.mpf(1)), ("e=2.7183 [shape split]", E),
                 ("Otter=2.9558", OTTER), ("2.2209 (K_6 exact sup n<=6)", mp.mpf("120") ** (mp.mpf(1) / 6)),
                 ("3.5 (target)", mp.mpf("3.5")), ("20 (a01 tau route)", mp.mpf(20))]:
    for chat in (CHAT_LOC, CHAT_ENV):
        b1 = 2 * N * N / (chat * LAM * cb)
        r = b1 / K48
        print(" %-25s %-5s %-10s %-10s %-25s %s"
              % (name, chat, mp.nstr(chat * cb, 8), "YES" if chat * cb < THR else " NO",
                 mp.nstr(b1, 8), mp.nstr(r, 7) + (" (no shrink)" if r >= 1 else " (SHRINKS)")))

# ------------------------------------------------------------------------------ [6]
hdr("[6] VERDICT")
print(" [T] phi^U expands into spanning trees (Penrose), signs +-; block product (1.3) is exact.")
print(" [T] tau route cannot reach 3.5: sup tau^{1/n} >= %s  (periodic 4D box, exact)." % mp.nstr(rbox, 6))
print(" [T] the <=20 loss: e*Delta (positions) is charged twice; residual = shape entropy <= e.")
print(" [C] e (2.7183) and Otter (2.9558) both < 3.5 -> Cbar <= 3.5 provable at the level of ubar.")
print(" [E] per-cluster sup observed: (n-1)!^{1/n} to n=6, max 120^{1/6}=2.2209 (exact, clique<=6).")
for chat, tag in ((CHAT_LOC, "Chat=2 exact   "), (CHAT_ENV, "Chat=4 envelope")):
    b1 = 2 * N * N / (chat * LAM * E)
    print("     %s: Chat*Cbar = %s %s %s ; beta1^fix = %s = %s x (1/48) %s"
          % (tag, mp.nstr(chat * E, 8), "<" if chat * E < THR else ">", mp.nstr(THR, 8),
             mp.nstr(b1, 8), mp.nstr(b1 / K48, 7),
             ">= 1/48 (NO shrink)" if b1 >= K48 else "< 1/48 (shrinks)"))
print(" => (a) TRIGGERS with Chat=2: Cbar <= e <= 3.5, Chat*Cbar = 5.4366 < 7.0633,")
print("    beta1^fix = 1.2992*(1/48) >= 1/48 -> window back to 1/48, no shrink.")
print("    (with the stronger c01 tree identity Cbar=1: 2*1=2 < 7.0633, beta1 = 3.53*(1/48).)")
print(" => (b) does NOT trigger for Cbar itself (provable lower bound 2.2209 < 3.5); it DOES")
print("    trigger for the tau-envelope (>= 16.36), which is therefore abandoned as an object.")
print(" => (c) the '<=20' error is bookkeeping (position entropy twice), not the AM-GM algebra.")
print("    NB a01's 'non-vacuity check' (planar grids 2.05-2.37) understates real 4D clusters:")
print("    the genuine L=2 box already sits at 16.36, i.e. the Delta=20 envelope is nearly SATURATED.")
print("\n done.")
