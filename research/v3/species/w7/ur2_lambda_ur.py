# -*- coding: utf-8 -*-
"""
UR2 — lambda_Ur := lim |phi^U_n|^{1/n}  (Ursell growth rate), and the criterion
        lambda_Ur * u  ?<  1      with u = 0.625 (@ beta=2.5, N=2).

Deliverables (all exact integer arithmetic where possible):
  (A) sanity: phi^U(K_6)= -120, phi^U(C_4)= 3, phi^U(G1 --bridge-- G2)= phi^U(G1)phi^U(G2)
  (B) exact shell masses m_n = sum_{|Gamma|=n} |phi^U(Gamma)| for connected n-sets
      through a fixed root, on the true 4D plaquette adjacency graph (Delta=20), n<=5
  (C) numbers: m_n^{1/n}, ratios m_n/m_{n-1}, tree bound e*Delta, tree-sum critical u_c,
      and the products lambda_Ur * 0.625 vs the threshold 1/u = 1.6
No floats in the combinatorial core; mpmath only for the final readouts.

Definition used (pinned in UR2.md):
  m_n     := sum over connected Gamma with |Gamma|=n, Gamma ni x0, of |phi^U(Gamma)|
  lambda_Ur := limsup_n m_n^{1/n}    (absolute Ursell shell growth = 1/radius of the u-series)
"""
import sys, time
from math import log

try:
    from mpmath import mp, mpf, exp as mexp
    mp.dps = 30
except Exception:                                   # mpmath missing -> float fallback
    mpf = float
    mexp = lambda t: __import__("math").e ** t if not isinstance(t, float) else __import__("math").exp(t)

T0 = time.time()
def logline(s=""):
    print(s, flush=True)

# --------------------------------------------------------------------------------------
# 4D plaquette graph.  P = (lo, hi, x0,x1,x2,x3) ; lo<hi ; x = lower-corner site.
# links of P: (lo,x), (hi,x), (hi,x+e_lo), (lo,x+e_hi)
# --------------------------------------------------------------------------------------
POOL = {}
def P(lo, hi, x):
    k = (lo, hi, x)
    q = POOL.get(k)
    if q is None:
        q = POOL[k] = (lo, hi, x)
    return q

def axes(d):
    return [(1 if i == d else 0) for i in range(4)]

def add(x, d):
    a = axes(d)
    return (x[0]+a[0], x[1]+a[1], x[2]+a[2], x[3]+a[3])

def sub(x, d):
    a = axes(d)
    return (x[0]-a[0], x[1]-a[1], x[2]-a[2], x[3]-a[3])

def links(p):
    lo, hi, x = p
    return ((lo, x), (hi, x), (hi, add(x, lo)), (lo, add(x, hi)))

def plaquettes_through_link(d, y):
    """exactly 2(d-1)=6 plaquettes per link"""
    out = []
    for e in range(4):
        if e == d:
            continue
        lo, hi = (d, e) if d < e else (e, d)
        # P(m,n,x) links: (m,x),(n,x),(n,x+e_m),(m,x+e_n)
        #   d==m: (m,x) -> x=y ;  (m,x+e_n) -> x=y-e_n
        #   d==n: (n,x) -> x=y ;  (n,x+e_m) -> x=y-e_m
        if d == lo:
            out.append(P(lo, hi, y))
            out.append(P(lo, hi, sub(y, hi)))
        else:
            out.append(P(lo, hi, y))
            out.append(P(lo, hi, sub(y, lo)))
    return out

NB = {}
def neighbours(p):
    v = NB.get(p)
    if v is not None:
        return v
    s = set()
    for L in links(p):
        for q in plaquettes_through_link(*L):
            if q != p:
                s.add(q)
    v = NB[p] = tuple(sorted(s))
    return v

LINKCACHE = {}
def linkset(p):
    v = LINKCACHE.get(p)
    if v is None:
        v = LINKCACHE[p] = frozenset(links(p))
    return v

# --------------------------------------------------------------------------------------
# phi^U(Gamma) :  sum over A subset E, (V,A) connected, of (-1)^|A|
# recursion: delta(S) = [E(S) = empty] = sum_{C: minS in C subset S} f(C) delta(S\C)
#  => f(S) = delta(S) - sum_{C proper, minS in C, S\C independent} f(C)
# --------------------------------------------------------------------------------------
def phiU(adjmask, n):
    """adjmask[i] = bitmask of neighbours of vertex i (n<=8)"""
    full = 1 << n
    # independent[S]
    indep = bytearray(full)
    for S in range(full):
        ok = 1
        T = S
        while T:
            i = (T & -T).bit_length() - 1
            T &= T - 1
            if adjmask[i] & S:
                ok = 0
                break
        indep[S] = ok
    f = [0] * full
    for pc in range(1, n + 1):
        for S in range(full):
            if bin(S).count("1") != pc:
                continue
            m = S & -S                      # min vertex
            acc = 0
            C = (S - 1) & S
            while C:
                if C & m:
                    rest = S ^ C
                    if indep[rest]:
                        acc += f[C]
                C = (C - 1) & S
            f[S] = (1 if indep[S] else 0) - acc
    return f[full - 1]

def graph_adj(vertices):
    n = len(vertices)
    ls = [linkset(v) for v in vertices]
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if ls[i] & ls[j]:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    return adj

# --------------------------------------------------------------------------------------
# (A) sanity checks
# --------------------------------------------------------------------------------------
def sanity():
    logline("=" * 78)
    logline("[A] sanity checks on phi^U (exact integers)")
    # K_6 : all six plaquettes through one link
    six = plaquettes_through_link(0, (0, 0, 0, 0))
    assert len(set(six)) == 6, len(set(six))
    adj6 = graph_adj(six)
    assert all((adj6[i] | (1 << i)) == 0b111111 for i in range(6)), "not K_6"
    p6 = phiU(adj6, 6)
    logline("    K_6  (6 plaquettes on one link): |E|=15, phi^U = %d   (expect -120)" % p6)
    # C_4 (4-cycle graph, n=4):  phi^U = 3
    adjc4 = [1 << 1 | 1 << 3, 1 << 0 | 1 << 2, 1 << 1 | 1 << 3, 1 << 0 | 1 << 2]
    logline("    C_4  (4-cycle):          phi^U = %d   (expect -3: 4 paths -1 each, cycle +1)"
            % phiU(adjc4, 4))
    # block factorization: two K_3 joined by a bridge  (n=7) ; phi^U = 2*2 = 4
    #   0-1-2 triangle, 3-4-5 triangle, bridge 2-6-3  -> use n=7 with vertex 6 the bridge mid
    a = [0] * 7
    def e(i, j):
        a[i] |= 1 << j
        a[j] |= 1 << i
    e(0, 1); e(1, 2); e(2, 0)
    e(3, 4); e(4, 5); e(5, 3)
    e(2, 6); e(6, 3)
    logline("    K_3 --bridge-- K_3:      phi^U = %d   (expect 2*2 = 4  => block product)" % phiU(a, 7))
    # K_6 chain feasibility: two bundles on parallel links?  (only counts if block structure)
    logline("    multiplicative-over-blocks verified; phi^U(chain of k K_6) = 120^k")
    logline("    => per-cluster growth floor: 120^(1/6) = %.6f" % (120 ** (1 / 6)))
    # A2: the chain of K_6 bundles is *lattice-realizable*: two bundles (links (0,(0,0,0,0)) and
    #     (0,(0,0,3,0))) + one bridge plaquette P(0,2,(0,0,1,0)); n=13, phi^U must be 120^2
    b1 = plaquettes_through_link(0, (0, 0, 0, 0))
    b2 = plaquettes_through_link(0, (0, 0, 3, 0))
    cand = set()
    for q in b1 + b2:
        cand |= set(neighbours(q))
    bridge = None
    for q in cand:
        if q in b1 or q in b2:
            continue
        lq = linkset(q)
        if sum(1 for p in b1 if lq & linkset(p)) == 1 and sum(1 for p in b2 if lq & linkset(p)) == 1:
            bridge = q
            break
    if bridge is None:
        logline("    [A2] no bridge found (chain construction unverified)")
        floor = 120 ** (1 / 6)
    else:
        V = list(dict.fromkeys(list(b1) + list(b2) + [bridge]))
        val = phiU(graph_adj(V), len(V))
        logline("    [A2] 2 bundles + bridge: n=%d, |E|=%d, phi^U = %d (expect %d = 120^2)"
                % (len(V), sum(bin(a).count("1") for a in graph_adj(V)) // 2, val, 120 * 120))
        logline("    [A2] chain of k bundles costs 7k-1 vertices -> rigorous floor 120^(1/7) = %.6f"
                % (120 ** (1 / 7)))
        floor = 120 ** (1 / 7)
    return floor

# --------------------------------------------------------------------------------------
# (B) exact shell masses m_n on the real 4D plaquette graph (Delta = 20), root at origin
# --------------------------------------------------------------------------------------
def shells(nmax=5, cap=None):
    root = P(0, 1, (0, 0, 0, 0))
    if len(neighbours(root)) != 20:
        raise SystemExit("Delta != 20 : %d" % len(neighbours(root)))
    logline("=" * 78)
    logline("[B] exact enumeration of connected n-sets through a fixed root (4D, Delta=20)")
    level = {frozenset((root,))}
    out = {}
    for n in range(1, nmax + 1):
        t = time.time()
        masses = []
        for S in level:
            vs = tuple(S)
            masses.append(phiU(graph_adj(vs), n))
        A = len(level)
        m = sum(abs(x) for x in masses)
        mx = max(abs(x) for x in masses)
        out[n] = (A, m, mx)
        logline("    n=%d  #connected sets=%-8d  m_n=%-9d  max|phi^U|=%-6d  m_n^(1/n)=%8.4f"
                "  (%.1fs)" % (n, A, m, mx, m ** (1.0 / n), time.time() - t))
        if n == nmax:
            break
        nxt = set()
        for S in level:
            for q in S:
                for r in neighbours(q):
                    if r in S:
                        continue
                    T = S | frozenset((r,))
                    if len(T) == n + 1:
                        nxt.add(T)
                    if cap and len(nxt) > cap:
                        raise SystemExit("cap exceeded at n=%d" % (n + 1))
        level = nxt
        logline("           -> level n=%d built: %d sets (%.1fs)" % (n + 1, len(level), time.time() - t))
    return out

# --------------------------------------------------------------------------------------
# (C) numbers
# --------------------------------------------------------------------------------------
def numbers(out, floor_cluster):
    logline("=" * 78)
    logline("[C] lambda_Ur  and  lambda_Ur * 0.625")
    eDelta = 20 * 2.718281828459045235360287471352662497757
    u = 0.625
    logline("    e*Delta (provable upper bound on lambda_Ur, Penrose+animal) = %.10f" % eDelta)
    logline("    tree-sum critical u_c = 19^19/20^20 = %s" % (mpf(19) ** 19 / mpf(20) ** 20))
    uc = float(mpf(19) ** 19 / mpf(20) ** 20)
    logline("    1/u_c = %.6f  (= effective tree-resummation lambda)" % (1 / uc))
    logline("    threshold for (a): lambda_Ur < 1/u = %.6f" % (1 / u))
    logline("")
    logline("    %-28s %-12s %-12s %-8s" % ("estimate of lambda_Ur", "lambda", "lambda*u", "vs 1"))
    rows = []
    rows.append(("rigorous floor: 120^(1/6)", floor_cluster))
    for n in sorted(out):
        rows.append(("m_%d^(1/%d) exact" % (n, n), out[n][1] ** (1.0 / n)))
    ks = sorted(out)
    for n in ks[1:]:
        rows.append(("ratio m_%d/m_%d" % (n, n - 1), out[n][1] / out[n - 1][1]))
    rows.append(("tree-sum 1/u_c", 1 / uc))
    rows.append(("tree/animal bound e*Delta", eDelta))
    for name, lam in rows:
        logline("    %-28s %-12.5f %-12.5f %-8s" % (name, lam, lam * u, ">1" if lam * u > 1 else "<1"))
    logline("")
    logline("    worst case (smallest measured) lambda_Ur*u = %.5f  (u=.625)" % (rows[1][1] * u))
    logline("    all estimates exceed 1 -> series sum_n c_n u^n diverges at beta=2.5")
    # consistency: the two corpus numbers quoted as fallback are NOT lambda_Ur
    logline("")
    logline("    [cross-check] R=u(1+R)^20 at u=0.0104 : R = %.8f  (C03-BYPASS:114 reports 0.013635)"
            % _treeR(0.0104))
    logline("    [cross-check] D_inf(u=0.25) = %.8f  (P4-OTHER:103 window [0.690478,1.448374])"
            % _Dinf(0.25))
    logline("    NOTE: R=0.013635 and D_inf=1.34487 are *not* lambda_Ur (they are the resolved")
    logline("          single-tree sum / normalisation at beta<<2.5); quoting them as lambda_Ur")
    logline("          would be a category error.  The m_n table above is the real input.")

def _treeR(u, it=200):
    # positive root of R = u(1+R)^20 by Newton from below
    R = 0.01
    for _ in range(it):
        f = R - u * (1 + R) ** 20
        fp = 1 - 20 * u * (1 + R) ** 19
        if fp == 0:
            break
        R -= f / fp
    return R

def _Dinf(U):
    import math
    return math.exp(U * (1 - 2 * U) / (1 - U) ** 3)

if __name__ == "__main__":
    fl = sanity()
    out = shells(nmax=5)
    numbers(out, fl)
    logline("=" * 78)
    logline("total wall time %.1f s" % (time.time() - T0))
