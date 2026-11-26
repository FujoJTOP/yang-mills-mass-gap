"""kh2_prodR.py -- KH2 part 2: is the leading factor 4 in K_hat = 4*sup correct?

K_hat := sup_gamma ( prod_l R(l) * tau_gamma )^{1/|gamma|}          (C06-WINDOW.md:118, (5.1))
R(kappa): kappa_l = # plaquettes of the cluster on link l;  R(1)=1 ("kappa>=2 才激活"),
          exact R(2)=2, R(4)=5/2, R(6)=2.8 (Z12-FIXC:70, C06:107), proven envelope R(kappa)<=2^{kappa/2}.
tau_gamma = spanning trees of the cluster's adjacency graph (Kirchhoff).

For a *closed surface* cluster (every link kappa=2) sum kappa = 4n gives (prodR)^{1/n}=4 exactly.
For the multi-beam clusters the boundary links have kappa=1 -> they are INACTIVE, so
(prodR)^{1/n} is much smaller than 4.  This script measures it.
"""
import itertools, math
from fractions import Fraction

D = 4
def log(*a): print(*a, flush=True)
def unit(i):
    v = [0]*D; v[i] = 1; return tuple(v)
def add(u, v): return tuple(u[k]+v[k] for k in range(D))
DIRS = [(a, b) for a in range(D) for b in range(a+1, D)]

Rexact = {1: 1.0, 2: 2.0, 4: 2.5, 6: 2.8}     # tabulated
def Rbound(k): return 1.0 if k == 1 else 2.0**(k/2.0)

def build_box(RMAX):
    sites = [x for x in itertools.product(range(-RMAX, RMAX+1), repeat=D)]
    plaq = [(x, a, b) for x in sites for (a, b) in DIRS]
    idx = {p: i for i, p in enumerate(plaq)}
    linkmap = {}
    for i, (x, a, b) in enumerate(plaq):
        for L in ((x, a), (add(x, unit(a)), b), (x, b), (add(x, unit(b)), a)):
            linkmap.setdefault(L, []).append(i)
    nb = [set() for _ in plaq]
    for L, ps in linkmap.items():
        for i in range(len(ps)):
            for j in range(i+1, len(ps)):
                nb[ps[i]].add(ps[j]); nb[ps[j]].add(ps[i])
    return plaq, idx, linkmap, nb

def kappa_hist(vs, linkmap):
    h = {}
    for L, ps in linkmap.items():
        c = sum(1 for p in ps if p in vs)
        if c:
            h[c] = h.get(c, 0) + 1
    return h

def spanning_trees(n, rows):
    """Kirchhoff: any cofactor of the Laplacian (exact integer determinant)."""
    L = [[0]*n for _ in range(n)]
    for i in range(n):
        d = bin(rows[i]).count("1")
        L[i][i] = d
        for j in range(n):
            if i != j and (rows[i] >> j) & 1:
                L[i][j] = -1
    A = [[Fraction(L[i][j]) for j in range(n-1)] for i in range(n-1)]
    det = Fraction(1)
    for c in range(n-1):
        p = next((r for r in range(c, n-1) if A[r][c] != 0), None)
        if p is None: return 0
        if p != c:
            A[c], A[p] = A[p], A[c]; det = -det
        det *= A[c][c]
        inv = A[c][c]
        for r in range(c+1, n-1):
            f = A[r][c]/inv
            if f:
                for k2 in range(c, n-1):
                    A[r][k2] -= f*A[c][k2]
    return int(det)

def rows_of(vs, nb):
    n = len(vs); rows = [0]*n
    for i, v in enumerate(vs):
        r = 0
        for j, w in enumerate(vs):
            if j != i and w in nb[v]: r |= 1 << j
        rows[i] = r
    return rows

plaq, idx, linkmap, nb = build_box(2)
base = (0, 0, 0, 0)
THR = 96/(5*math.e)
PHI = {1: 120, 2: 125400, 3: 155244960}     # kh2_limsup.py, exact

log(f"threshold 7.0632852705;  K_hat < T  <=>  K_hat/4 < 1.7658213 (only if (prodR)^{{1/n}}=4)")
log(f"{'cluster':<22}{'n':>3}{'kappa hist':>26}{'prodR(bd)':>11}{'root_bd':>9}{'tau':>12}{'K_bd':>9}{'root_ex':>9}{'K_ex':>9}")
rows_out = []
cases = {}
for k in (1, 2, 3):
    vs = set()
    for a in range(k):
        vs |= set(linkmap[(base, a)])
    vs = sorted(vs); n = len(vs)
    h = kappa_hist(vs, linkmap)
    rb = 1.0
    for c, m in h.items(): rb *= Rbound(c)**m
    rex = 1.0
    ex_known = True
    for c, m in h.items():
        if c in Rexact: rex *= Rexact[c]**m
        else: ex_known = False; rex *= Rbound(c)**m
    rr = rows_of(vs, nb)
    tau = spanning_trees(n, rr)
    ph = PHI[k]
    Kbd = (rb*ph)**(1/n)
    Kex = (rex*ph)**(1/n)
    cases[k] = (n, ph, rb**(1/n), tau / ph, Kbd, Kex)
    log(f"beam x{k:<16}{n:>3}{str(h):>26}{rb:>11.3f}{rb**(1/n):>9.4f}{tau:>12}"
        f"{Kbd:>9.4f}{(rex**(1/n)):>9.4f}{Kex:>9.4f}"
        + ("" if ex_known else " (R3/R5 by bound)"))

# cube witness (all kappa = 2): the C06 5.3 number
n = 6
rows = [0]*6
E = [(i, j) for i in range(2) for j in range(2, 4)] + \
    [(i, j) for i in range(2) for j in range(4, 6)] + \
    [(i, j) for i in range(2, 4) for j in range(4, 6)]
for (i, j) in E:
    rows[i] |= 1 << j; rows[j] |= 1 << i
tau_c = spanning_trees(n, rows)
log(f"\ncube K_2,2,2: tau = {tau_c} (C06 5.3 quotes 384) ; prodR=2^12=4096 ; "
    f"(prodR*tau)^{{1/6}} = {(4096*tau_c)**(1/6):.10f}")

log("\n[verdict] per-witness K_hat (phi^U route, exact integer phi, per-cluster prodR):")
for k, (n, ph, rr, ratio, Kbd, Kex) in cases.items():
    log(f"  beam x{k}: n={n} K_hat <= {Kex:.4f} (exact-R where tabulated) / {Kbd:.4f} (envelope R<=2^(k/2))")
log(f"  cube:     n=6 K_hat = {(4096*tau_c)**(1/6):.4f}   > T = {THR:.4f}   [T]+[E] (C06 5.3: 10.784)")
log(f"\n  => sup_gamma is still attained on the CLOSED-SURFACE (all kappa=2) family;")
log(f"     the phi^U giants (n=11,15) carry (prodR)^{{1/n}} = {cases[2][2]:.4f}, {cases[3][2]:.4f}, NOT 4.")
log(f"     so 'K_hat >= 4*sup' OVER-counts those witnesses by 4/{cases[3][2]:.4f} = {4/cases[3][2]:.2f}x.")
log(f"     but K_hat >= {(4096*tau_c)**(1/6):.4f} > {THR:.4f} still holds (ratio "
    f"{(4096*tau_c)**(1/6)/THR:.3f}).")
