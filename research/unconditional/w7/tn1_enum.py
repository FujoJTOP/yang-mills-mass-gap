"""tn1: exact enumeration of rooted plaquette clusters in the 4D hypercubic
plaquette graph (adjacency = sharing a link), and of the rooted TREE-SUM

    t_n := #{ (Gamma,T) : Gamma ni p0, |Gamma| = n, Gamma connected,
                        T spanning tree of G_Gamma }  =  sum_Gamma tau(G_Gamma)

Also: a_n = #animals (connected Gamma), m_n = sum |phi^U(Gamma)| (Ursell shell
mass, for cross-checking the corpus values 1,20,514,14748,450255).

No git. Pure python3 + fractions. Usage:  py -3.12 tn1_enum.py N
"""
import sys, time, itertools
from fractions import Fraction

OFF = 8
ORI = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
OIDX = {o: i for i, o in enumerate(ORI)}


def pack(x, mu, nu):
    if mu > nu:
        mu, nu = nu, mu
    v = 0
    for c in x:
        v = (v << 5) | (c + OFF)
    return (v << 3) | OIDX[(mu, nu)]


def unpack(p):
    o = ORI[p & 7]
    v = p >> 3
    x = []
    for _ in range(4):
        x.append((v & 31) - OFF)
        v >>= 5
    x.reverse()
    return tuple(x), o[0], o[1]


def links(p):
    x, mu, nu = unpack(p)
    out = []

    def sh(y, r, k):
        y = list(y)
        y[r] += k
        return tuple(y)

    out.append((x, mu))
    out.append((sh(x, mu, 1), nu))
    out.append((x, nu))
    out.append((sh(x, nu, 1), mu))
    return set(out)


_NC = {}


def neigh(p):
    r = _NC.get(p)
    if r is not None:
        return r
    x, mu, nu = unpack(p)
    res = set()
    for (y, rho) in links(p):
        for sigma in range(4):
            if sigma == rho:
                continue
            res.add(pack(y, rho, sigma))
            yy = list(y)
            yy[sigma] -= 1
            res.add(pack(tuple(yy), rho, sigma))
    res.discard(p)
    _NC[p] = frozenset(res)
    return _NC[p]


# ---------------- exact integer determinant (Bareiss) ----------------
def det_bareiss(M):
    M = [row[:] for row in M]
    n = len(M)
    if n == 0:
        return 1
    sign = 1
    prev = 1
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


def tau_kirchhoff(cluster):
    """number of spanning trees of the induced plaquette graph (link-sharing)."""
    n = len(cluster)
    if n == 1:
        return 1
    Lk = {}
    for p in cluster:
        Lk[p] = links(p)
    Lap = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if Lk[cluster[i]] & Lk[cluster[j]]:
                Lap[i][i] += 1
                Lap[j][j] += 1
                Lap[i][j] -= 1
                Lap[j][i] -= 1
    M = [row[1:] for row in Lap[1:]]
    return det_bareiss(M)


def phi_ursell(cluster):
    """phi^U(Gamma) = sum_{A subset E, (Gamma,A) connected} (-1)^|A|."""
    n = len(cluster)
    if n == 1:
        return 1
    Lk = [links(p) for p in cluster]
    E = []
    for i in range(n):
        for j in range(i + 1, n):
            if Lk[i] & Lk[j]:
                E.append((i, j))
    s = 0
    for mask in range(1 << len(E)):
        par = list(range(n))

        def f(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a

        bits = 0
        for b, (i, j) in enumerate(E):
            if mask >> b & 1:
                bits += 1
                ri, rj = f(i), f(j)
                if ri != rj:
                    par[ri] = rj
        if len({f(i) for i in range(n)}) == 1:
            s += -1 if bits & 1 else 1
    return s


def main():
    maxn = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    urell_max = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    t0 = time.time()
    root = pack((0, 0, 0, 0), 0, 1)
    print("root plaquette degree =", len(neigh(root)), flush=True)
    levels = [((root,),)]
    for n in range(1, maxn):
        nxt = set()
        for A in levels[-1]:
            for p in A:
                for q in neigh(p):
                    if q not in A:
                        nxt.add(tuple(sorted(A + (q,))))
        levels.append(nxt)
        print("  level %d: %d animals   (%.1fs)" % (n + 1, len(nxt), time.time() - t0), flush=True)

    animals = [list(lv) if isinstance(lv, set) else list(lv) for lv in levels]
    print("\n n |     a_n     |   t_n=sum tau  |  tau_max | distinct tau | t_n^(1/n) | a_n^(1/n)")
    tn_list = []
    for n in range(1, maxn + 1):
        lvl = animals[n - 1]
        tot = 0
        taus = []
        for cl in lvl:
            tt = tau_kirchhoff(list(cl))
            tot += tt
            taus.append(tt)
        tn_list.append(tot)
        print("%2d | %11d | %13d | %8d | %3d | %.6f | %.6f"
              % (n, len(lvl), tot, max(taus), len(set(taus)), tot ** (1.0 / n), len(lvl) ** (1.0 / n)), flush=True)

    # Ursell shell mass cross-check (corpus: 1,20,514,14748,450255)
    print("\n n |  m_n = sum |phi^U| ")
    for n in range(1, min(urell_max, maxn) + 1):
        lvl = animals[n - 1]
        m = sum(abs(phi_ursell(list(cl))) for cl in lvl)
        print("%2d | %12d" % (n, m), flush=True)

    print("\nelapsed %.1fs" % (time.time() - t0))
    print("t_n =", tn_list)


main()
