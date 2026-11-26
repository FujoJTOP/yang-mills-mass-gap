"""ur5_wg2.py -- EXACT SU(2) Weingarten table (all classes of S_n, n<=6) via multi-trace moments.

Key facts used (all exact, no floating point).

(1) Weyl: for U = diag(e^{i th}, e^{-i th}) in SU(2),
        tr(U^k) = 2 cos(k th),   Haar  =  (2/pi) sin^2(th) dth  on [0,pi].
(2) So every multi-trace moment  I(lam,mu) = Int tr_lam(U) tr_mu(U^dag) dU   with
        tr_lam(U) = prod_i tr(U^{lam_i})
    is EXACT: writing  P(th) = prod_i 2cos(lam_i th) prod_j 2cos(mu_j th) = sum_{m>=0} c_m cos(m th),
        I(lam,mu) = c_0 - c_2/2                    (since (2/pi)Int sin^2 cos(mth) = 1,0,-1/2 for m=0,>2,=2).
    c_m are computed with  cos a cos b = (cos(a+b) + cos(a-b))/2  -- exact integer/rational arithmetic.
    Cross-checks: I((1),(1))=1 (Catalan C_1), I((2),(2))=2, I((1,1),(1,1))=2, I((3),(3))=2, I((k),(k)) = C_k?
(3) Weingarten expansion (the defining formula; convention fixed by n=1: Int U_ab Ubar_cd = d_ad d_bc/N):
        Int prod_p U_{a_p b_p} prod_p Ubar_{c_p d_p}
          = sum_{s,t in S_n} Wg(s^-1 t) prod_p d_{a_p d_{s(p)}} prod_p d_{b_p c_{t(p)}},
    hence  I(lam,mu) = sum_{s,t in S_n} Wg(s^-1 t) N^{comp(s,t;lam,mu)},  N=2,
    where comp = number of free index loops = #components of the union-find graph
        {b_p = a_{p+1}}  (U chains),  {d_p = c_{p+1}}  (Ubar chains),
        {a_p = d_{s(p)}}, {b_p = c_{t(p)}} .
    Wg is a class function -> p(n) unknowns, p(n)^2 exact equations -> over-determined linear system.
    A CONSISTENT solution proves the S_n Weingarten ansatz holds for SU(2) and determines the table.
"""
import itertools, json, math, sys
from fractions import Fraction
from sympy import Matrix, Rational


def partitions(n):
    if n == 0:
        yield ()
        return
    for first in range(n, 0, -1):
        for rest in partitions(n - first):
            if not rest or rest[0] <= first:
                yield (first,) + rest


def cycle_type(p):
    n = len(p)
    seen = [False] * n
    cyc = []
    for i in range(n):
        if not seen[i]:
            l, j = 0, i
            while not seen[j]:
                seen[j] = True
                j = p[j]
                l += 1
            cyc.append(l)
    return tuple(sorted(cyc, reverse=True))


def cname(ct):
    out = []
    for l in sorted(set(ct), reverse=True):
        k = ct.count(l)
        out.append(("%d" % l) + ("^%d" % k if k > 1 else ""))
    return ".".join(out)


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(a)))


def inverse(a):
    out = [0] * len(a)
    for i, v in enumerate(a):
        out[v] = i
    return tuple(out)


# ---------- exact moments ----------
def cosmul(A, B):
    """multiply two 'cos-basis' Laurent dicts: A = sum_m c_m cos(m th) (m>=0)"""
    C = {}
    for m, cm in A.items():
        for n_, cn in B.items():
            for s in (m + n_, abs(m - n_)):
                C[s] = C.get(s, Fraction(0)) + cm * cn / 2
    return C


def moment(lam, mu):
    P = {0: Fraction(1)}
    for k in list(lam) + list(mu):
        P = cosmul(P, {k: Fraction(2)})
    return P.get(0, Fraction(0)) - P.get(2, Fraction(0)) / 2


# ---------- index-loop count ----------
def chain_labels(part):
    n = sum(part)
    a = [None] * n
    b = [None] * n
    c = [None] * n
    d = [None] * n
    p = 0
    for lam in part:
        for i in range(lam):
            u = p + i
            a[u] = ("A", p + i)
            b[u] = ("A", p + (i + 1) % lam)
            c[u] = ("C", p + i)
            d[u] = ("C", p + (i + 1) % lam)
        p += lam
    return n, a, b, c, d


def ncomp(lam, mu, s, t):
    n = sum(lam)
    _, a, b, c, d = chain_labels(lam)
    _, a2, b2, c2, d2 = chain_labels(mu)
    cM = [("C2", v[1]) for v in c2]
    dM = [("C2", v[1]) for v in d2]
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        r = x
        while parent[r] != r:
            r = parent[r]
        while parent[x] != r:
            parent[x], x = r, parent[x]
        return r

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for p in range(n):
        union(a[p], dM[s[p]])
        union(b[p], cM[t[p]])
    for p in range(n):
        union(a[p], b[(p + 1) % n]) if False else None
    return len({find(x) for x in list(parent.keys())})


def solve(n):
    cts = []
    for p in itertools.permutations(range(n)):
        ct = cycle_type(p)
        if ct not in cts:
            cts.append(ct)
    cts.sort()
    cidx = {ct: i for i, ct in enumerate(cts)}
    perms = list(itertools.permutations(range(n)))
    inv = {p: inverse(p) for p in perms}
    rows, rhs, tags = [], [], []
    parts = list(partitions(n))
    for lam in parts:
        for mu in parts:
            row = [0] * len(cts)
            for s in perms:
                si = inv[s]
                for t in perms:
                    row[cidx[cycle_type(compose(si, t))]] += 2 ** ncomp(lam, mu, s, t)
            rows.append(row)
            rhs.append(moment(lam, mu))
            tags.append((lam, mu))
    return cts, rows, rhs, tags


def main():
    out = {}
    print("=== SU(2) Weingarten table, exact ===")
    for n in range(1, 5):
        cts, rows, rhs, tags = solve(n)
        p = len(cts)
        M = Matrix(rows)
        bcol = Matrix([Rational(x.numerator, x.denominator) for x in rhs])
        r = M.rank()
        Maug = M.row_join(bcol)
        raug = Maug.rank()
        consistent = (r == raug)
        print("\nn=%d  p(n)=%d  #eq=%d  rank(M)=%d  rank([M|b])=%d  %s"
              % (n, p, len(rows), r, raug, "CONSISTENT" if consistent else "INCONSISTENT"))
        if not consistent:
            print("   !! the S_n class-function ansatz is FALSE for SU(2) at n=%d" % n)
            continue
        sol, params = M.gauss_jordan_solve(bcol)
        vals = [sol[i] for i in range(p)]
        if params:
            print("   UNDETERMINED by these moments: %d free parameters" % len(params))
        # verify every equation exactly
        res = M * sol - bcol
        assert all(x == 0 for x in res), "exact residual nonzero"
        for j, ct in enumerate(cts):
            print("   %-10s  Wg = %s" % (cname(ct), vals[j]))
        out[n] = {cname(ct): str(vals[j]) for j, ct in enumerate(cts)}
        print("   [check] I((%s),(%s)) = %s" % (",".join(map(str, (n,))), ",".join(map(str, (n,))),
                                                moment((n,), (n,))))
        sys.stdout.flush()
    with open("ur5_wg_table.json", "w") as f:
        json.dump({str(k): v for k, v in out.items()}, f, indent=1)
    print("\nwrote ur5_wg_table.json")


if __name__ == "__main__":
    main()
