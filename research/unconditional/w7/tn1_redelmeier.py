"""tn1 part 2: Redelmeier-style exhaustive enumeration (each animal once, O(depth)
memory) giving exact a_n and t_n = sum_Gamma tau(G_Gamma) for n <= N.

Cross-checks (must reproduce): a_n = 1,20,458,11132,280175 ; t_n = 1,20,570,18776,669215 (n<=5)
Usage: py -3.12 tn1_redelmeier.py N
"""
import sys, time
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
    return frozenset(out)


_LK = {}
_NC = {}


def L(p):
    r = _LK.get(p)
    if r is None:
        r = _LK[p] = links(p)
    return r


def neigh(p):
    r = _NC.get(p)
    if r is not None:
        return r
    x, mu, nu = unpack(p)
    res = set()
    for (y, rho) in L(p):
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
        pv = M[k][k]
        for i in range(k + 1, n):
            Mik = M[i][k]
            if Mik:
                for j in range(k + 1, n):
                    M[i][j] = (M[i][j] * pv - Mik * M[k][j]) // prev
            for j in range(k + 1, n):
                pass
        prev = pv
        for i in range(k + 1, n):
            M[i][k] = 0
    return sign * M[n - 1][n - 1]


def tau_cluster(cl):
    n = len(cl)
    if n == 1:
        return 1
    Lap = [[0] * n for _ in range(n)]
    for i in range(n):
        Li = L(cl[i])
        for j in range(i + 1, n):
            if Li & L(cl[j]):
                Lap[i][i] += 1
                Lap[j][j] += 1
                Lap[i][j] -= 1
                Lap[j][i] -= 1
    return det_bareiss([row[1:] for row in Lap[1:]])


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    root = pack((0, 0, 0, 0), 0, 1)
    counts = [0] * (N + 1)
    tsum = [0] * (N + 1)
    clusters = []          # stack
    clset = set()
    t0 = time.time()

    # Redelmeier: sites popped along the current path are "used" forever on that
    # path (excl); a new frontier site must be free (not in cluster, not excl).
    def rec(untried, depth, excl):
        while untried:
            p = untried.pop()
            clset.add(p)
            clusters.append(p)
            n = len(clusters)
            counts[n] += 1
            tsum[n] += tau_cluster(clusters)
            if n < N:
                rec(untried | (neigh(p) - clset - excl), depth + 1, excl | {p})
            clusters.pop()
            clset.discard(p)

    clset.add(root)
    clusters.append(root)
    counts[1] = 1
    tsum[1] = 1
    rec(set(neigh(root)), 1, {root})

    print(" N |     a_n (animals) |        t_n = sum tau | t_n^(1/n) | a_n^(1/n) | ratio a_n/a_n-1 | ratio t_n/t_n-1")
    prev_a = prev_t = None
    for n in range(1, N + 1):
        ra = ("%.4f" % (counts[n] / prev_a)) if prev_a else "-"
        rt = ("%.4f" % (tsum[n] / prev_t)) if prev_t else "-"
        print("%2d | %17d | %20d | %.6f | %.6f | %8s | %8s"
              % (n, counts[n], tsum[n], tsum[n] ** (1.0 / n), counts[n] ** (1.0 / n), ra, rt))
        prev_a, prev_t = counts[n], tsum[n]
    print("elapsed %.1fs" % (time.time() - t0))
    print("a_n =", counts[1:])
    print("t_n =", tsum[1:])
    # exact rationals for the per-step ratios
    print("t_n/n! style check, sqrt etc. not needed here")


main()
