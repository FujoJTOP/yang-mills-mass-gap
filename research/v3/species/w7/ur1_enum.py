"""
UR1 — exact enumeration of the Ursell shell masses

    m_n := sum_{Gamma : p0 in Gamma, |Gamma| = n} |phi^U(Gamma)|

on the 4D plaquette graph (vertices = plaquettes, edge = shared link, Delta = 20),
p0 = fixed plaquette at the origin in the (0,1) plane   (convention C03-BYPASS:77).

phi^U(Gamma) = sum_{A subset E(G_Gamma), (Gamma,A) connected} (-1)^|A| prod_{e in A} f_e

Two bookkeeping conventions (both reported):
  (A) edge weight 1, one activity u per plaquette  => phi^U = u^n * mu(G_Gamma),
      mu(G) = sum_{A connected spanning} (-1)^|A|   (integer),  |phi^U| <= tau(G) u^n  (C03)
  (B) edge weight f_e = u as well                  => phi^U = u^n * C_Gamma(u),
      C_Gamma(x) = sum_{A connected spanning} x^|A|

SU(2), beta = 2.5, u = beta/4 = 0.625, eDelta = e*20 = 54.3656365691809047.

Usage:  py -3.12 ur1_enum.py [max_n] [--json OUT]
"""
import sys, time, json
from math import comb

DIRS = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
PLANES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
PMAP = {p: i for i, p in enumerate(PLANES)}
OFF = 20
U = 0.625
EDELTA = 20 * 2.718281828459045235360287       # e*20


def sc(x):
    return (((x[0] + OFF) << 15) | ((x[1] + OFF) << 10)
            | ((x[2] + OFF) << 5) | (x[3] + OFF))


def pc(x, mu, nu):
    if mu > nu:
        mu, nu = nu, mu
    return (sc(x) << 3) | PMAP[(mu, nu)]


def dec(c):
    s, pl = c >> 3, c & 7
    x = (((s >> 15) & 31) - OFF, ((s >> 10) & 31) - OFF,
         ((s >> 5) & 31) - OFF, (s & 31) - OFF)
    return x, PLANES[pl]


def add(x, d):
    return (x[0] + d[0], x[1] + d[1], x[2] + d[2], x[3] + d[3])


def sub(x, d):
    return (x[0] - d[0], x[1] - d[1], x[2] - d[2], x[3] - d[3])


_NB = {}


def neighbors(c):
    r = _NB.get(c)
    if r is not None:
        return r
    x, (mu, nu) = dec(c)
    links = [(x, mu), (add(x, DIRS[mu]), nu), (add(x, DIRS[nu]), mu), (x, nu)]
    out = set()
    for (y, d) in links:
        for k in range(4):
            if k == d:
                continue
            out.add(pc(y, d, k))                       # plaquette at y in plane (d,k)
            out.add(pc(sub(y, DIRS[k]), k, d))         # plaquette at y-k in plane (k,d)
    out.discard(c)
    _NB[c] = tuple(sorted(out))
    return _NB[c]


def animals(maxn, root=None):
    """connected n-plaquette sets containing root, n = 1..maxn (dedup by sorted tuple)"""
    if root is None:
        root = pc((0, 0, 0, 0), 0, 1)
    levels = [{(root,)}]
    for _ in range(1, maxn):
        prev, nxt = levels[-1], set()
        for S in prev:
            base = set(S)
            for p in S:
                for q in neighbors(p):
                    if q not in base:
                        nxt.add(tuple(sorted(S + (q,))))
        levels.append(nxt)
        print(f"  level n={len(levels)}: {len(nxt)} animals", flush=True)
    return levels


def graph_data(S):
    n = len(S)
    idx = {p: i for i, p in enumerate(S)}
    adj = [0] * n
    nbrsets = [set(neighbors(p)) for p in S]
    for i in range(n):
        for j in range(i + 1, n):
            if S[j] in nbrsets[i]:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    return adj


def conn_sum(adj, n, x):
    """sum over A subset E connected and spanning of x^|A|  (x=1 -> exact int)"""
    m = 1 << n
    ec = [0] * m
    for T in range(1, m):
        low = T & -T
        i = low.bit_length() - 1
        R = T ^ low
        ec[T] = ec[R] + bin(adj[i] & R).count('1')
    mE = ec[m - 1]
    one = 1 if isinstance(x, int) else 1.0
    # P(x) := sum_{A conn spanning} (-x)^|A| ;  g[T] = sum_{A subset E(T)} (-x)^|A| = (1-x)^e(T)
    g = [((one - x) ** k) for k in ec]
    f = [0] * m
    for T in range(1, m):
        low = T & -T
        rest = T ^ low
        tot = g[T]
        Uu = rest
        while True:
            if Uu != rest:
                P = low | Uu
                tot -= f[P] * g[T ^ P]
            if Uu == 0:
                break
            Uu = (Uu - 1) & rest
        f[T] = tot
    return f[m - 1], mE


def main():
    maxn = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    out_json = None
    if "--json" in sys.argv:
        out_json = sys.argv[sys.argv.index("--json") + 1]
    t0 = time.time()
    print(f"UR1: 4D plaquette graph (Delta=20), u={U}, eDelta={EDELTA:.16g}, maxn={maxn}")
    levels = animals(maxn)
    print(f"enumeration done in {time.time()-t0:.1f}s\n")

    rows = []
    for n in range(1, maxn + 1):
        S_list = levels[n - 1]
        t1 = time.time()
        sum_abs_mu = 0          # convention A: sum |mu|
        sum_abs_Cu = 0.0        # convention B: sum |C_Gamma(u)|
        sign_viol = 0           # count of sign(mu) != (-1)^(n-1)
        n_trees = 0
        max_mu = 0
        mu_hist = {}
        for S in S_list:
            adj = graph_data(S)
            mu, mE = conn_sum(adj, n, 1)
            if mE == n - 1:
                n_trees += 1
            cu, _ = conn_sum(adj, n, U)
            sum_abs_mu += abs(mu)
            sum_abs_Cu += abs(cu)
            max_mu = max(max_mu, abs(mu))
            mu_hist[mu] = mu_hist.get(mu, 0) + 1
            sgn = 1 if mu > 0 else (-1 if mu < 0 else 0)
            if sgn != ((-1) ** (n - 1)) and mu != 0:
                sign_viol += 1
        A = sum_abs_mu
        m_A = (U ** n) * A
        m_B = (U ** n) * sum_abs_Cu
        rows.append(dict(n=n, count=len(S_list), trees=n_trees,
                         sum_abs_mu=A, sum_abs_Cu=sum_abs_Cu, max_mu=max_mu,
                         sign_viol=sign_viol, m_A=m_A, m_B=m_B,
                         mu_hist={str(k): v for k, v in sorted(mu_hist.items())}))
        lam = A ** (1.0 / n) if A > 0 else 0.0
        print(f"n={n}: count={len(S_list):>8d} trees={n_trees:>8d} max|mu|={max_mu:>6d} "
              f"sign-viol={sign_viol:>7d} | sum|mu|={A:<12d} lam=|mu|^(1/n)={lam:.6f} "
              f"| m_A={m_A:.6e} m_B={m_B:.6e}  [{time.time()-t1:.1f}s]", flush=True)

    print("\n--- growth test: m_n ?= c * ubar^n,  ubar = eDelta*u = %.10f ---" % (EDELTA * U))
    for r in rows:
        n = r["n"]
        A = r["sum_abs_mu"]
        lam = A ** (1.0 / n)
        print(f"n={n}: m_n^(1/n)={r['m_A']**(1.0/n):.6f}  (u*lam={U*lam:.6f})  "
              f"lam/u^(? ) : lam^(1/n) ratio to eDelta = {lam/EDELTA:.6f}  "
              f"m_n/ubar^n={r['m_A']/(EDELTA*U)**n:.6e}")

    if out_json:
        with open(out_json, "w") as fh:
            json.dump(dict(u=U, eDelta=EDELTA, rows=rows), fh, indent=1)
        print(f"\njson -> {out_json}")
    print(f"total {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
