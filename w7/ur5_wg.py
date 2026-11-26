"""ur5_wg.py -- exact SU(2) Weingarten table Wg(sigma) for all classes of S_n, n<=6.

METHOD [T] (independent of published SU(N) formulas, hence a check on them):
  SU(2) = S^3:  U = a0*1 + i (a1 s1 + a2 s2 + a3 s3),  a uniform on S^3  (= Haar measure).
  For 2m linear forms L_k in a,
      Int_{S^3} prod_{k=1}^{2m} L_k = 1/(2^m (m+1)!) * sum_{perfect matchings P} prod_{ij in P} (v_i . v_j)
  with v_k in C^4 the coefficient vector of L_k in the basis (a0,a1,a2,a3).
  (S^3 moment formula: E[prod a_{i_1}..a_{i_2m}] = (2^m (m+1)!)^{-1} sum over pairings of prod deltas.)
  Hence every SU(2) matrix-element integral is an exact hafnian -> Gaussian rational.

  Then  Int prod_p U_{a_p b_p} prod_p Ubar_{c_p d_p}
      = sum_{s,t in S_n} Wg(s^-1 t) prod_p d_{a_p d_{s(p)}} prod_p d_{b_p c_{t(p)}} ,
  convention fixed by n=1:  Int U_ab Ubar_cd = d_ad d_bc / N.
  Wg is a class function, so a few dozen index configurations (indices in {0,1}, N=2) give an
  overdetermined rational linear system for the p(n) class values; solved exactly (Gauss-Jordan).

CHECKS  (i) n=1 -> 1/N = 1/2 ; n=2 -> 1/(N^2-1)=1/3, -1/(N(N^2-1))=-1/6  [the given values]
       (ii) Catalan moments <|tr U|^{2k}> = C_k for k<=6, computed *from the solved table*
      (iii) Monte-Carlo on raw index configurations (table not used)
       (iv) the naive "restricted lambda-sum" 1/n!^2 sum_{l(lam)<=N} f^2 chi/s_lam(1^N)
            printed alongside: it FAILS for n>=3 -> the det/epsilon-tensor correction is real.
"""
import itertools, math, sys, random, time
from fractions import Fraction as F
from sympy import Rational

N = 2
P2 = [[[1, 0], [0, 1]], [[0, 1], [1, 0]], [[0, 0], [0, 0]], [[1, 0], [0, -1]]]
# sigma2 = [[0,-i],[i,0]];  i*sigma2 = [[0,1],[-1,0]] (real) -- handled explicitly below
SIG = [None, None, None, None]


def _s1():
    return [[0, 1], [1, 0]]


# coefficient vectors: entry (i,j) of U = a0*d + i*s1*a1 + i*s2*a2 + i*s3*a3
#   i*s1 = [[0,i],[i,0]] ; i*s2 = [[0,1],[-1,0]] ; i*s3 = [[i,0],[0,-i]]
C1 = [[(0, 0), (0, 1)], [(0, 1), (0, 0)]]          # i*s1  -> (re,im)
C2 = [[(0, 0), (1, 0)], [(-1, 0), (0, 0)]]         # i*s2
C3 = [[(0, 1), (0, 0)], [(0, 0), (0, -1)]]         # i*s3


def mvec(i, j, conj=False):
    """(a0,a1,a2,a3) coefficients of U_ij (conj=False) or Ubar_ij, as complex pairs."""
    d = (1, 0) if i == j else (0, 0)
    out = [(d, None)]
    v = [d]
    for M in (C1, C2, C3):
        c = M[i][j]
        v.append((c[0], -c[1]) if conj else c)
    return tuple(v)


def cmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def cadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def dot(u, v):
    s = (F(0), F(0))
    for a, b in zip(u, v):
        s = cadd(s, cmul(a, b))
    return s


def haf(vecs, idx):
    if not idx:
        return (F(1), F(0))
    i0 = idx[0]
    rest = idx[1:]
    tot = (F(0), F(0))
    for k in range(len(rest)):
        j = rest[k]
        t = haf(vecs, rest[:k] + rest[k + 1:])
        tot = cadd(tot, cmul(dot(vecs[i0], vecs[j]), t))
    return tot


def integ(vecs):
    m = len(vecs) // 2
    h = haf(vecs, tuple(range(len(vecs))))
    den = 2 ** m * math.factorial(m + 1)
    return Rational(h[0], den) + Rational(h[1], den) * 1j * 0  # imaginary part must vanish


def integ_exact(vecs):
    m = len(vecs) // 2
    h = haf(vecs, tuple(range(len(vecs))))
    den = 2 ** m * math.factorial(m + 1)
    assert h[1] == 0, "imaginary part nonzero -- index pattern inconsistent"
    return Rational(h[0], den)


# ---------------- class bookkeeping ----------------
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


def classes_of(n):
    out = {}
    for p in itertools.permutations(range(n)):
        out.setdefault(cycle_type(p), 0)
        out[cycle_type(p)] += 1
    return sorted(out), out


def solve_wg(n, ncfg=None, seed=7):
    perms = list(itertools.permutations(range(n)))
    ctlist, _ = classes_of(n)
    cidx = {ct: i for i, ct in enumerate(ctlist)}
    inv = {p: inverse(p) for p in perms}
    cache = {}

    def classrow(s):
        if s not in cache:
            si = inv[s]
            cache[s] = [cidx[cycle_type(compose(si, t))] for t in perms]
        return cache[s]

    rng = random.Random(seed)
    if ncfg is None:
        ncfg = max(16, 4 * len(ctlist))
    eqs, rhs = [], []
    for _ in range(ncfg):
        ab = [rng.randint(0, 1) for _ in range(n)]
        cd = [rng.randint(0, 1) for _ in range(n)]
        ab2 = [rng.randint(0, 1) for _ in range(n)]
        cd2 = [rng.randint(0, 1) for _ in range(n)]
        vecs = [mvec(ab[p], cd[p]) for p in range(n)] + [mvec(ab2[p], cd2[p], True) for p in range(n)]
        lhs = integ_exact(vecs)
        row = [0] * len(ctlist)
        for s in perms:
            bad = False
            for p in range(n):
                if ab[p] != ab2[s[p]]:
                    bad = True
                    break
            if bad:
                continue
            cr = classrow(s)
            for ti, t in enumerate(perms):
                ok = True
                for p in range(n):
                    if cd[p] != cd2[t[p]]:
                        ok = False
                        break
                if ok:
                    row[cr[ti]] += 1
        if any(row):
            eqs.append(row)
            rhs.append(lhs)
    from sympy import Matrix
    M = Matrix(eqs)
    b = Matrix(rhs)
    sol, params = M.gauss_jordan_solve(b)
    if params:
        raise SystemExit("rank deficient n=%d: %d free params" % (n, len(params)))
    res = M * sol - b
    assert all(x == 0 for x in res)
    return ctlist, [sol[i] for i in range(len(ctlist))], len(eqs)


# ---------------- naive restricted lambda-sum (diagnostic) ----------------
def partitions(n):
    if n == 0:
        yield ()
        return
    for first in range(n, 0, -1):
        for rest in partitions(n - first):
            if not rest or rest[0] <= first:
                yield (first,) + rest


def hook_dim(lam):
    n = sum(lam)
    conj = [sum(1 for r in lam if r > j) for j in range(lam[0])] if lam else []
    prod = 1
    for i, r in enumerate(lam):
        for j in range(r):
            prod *= (r - j) + (conj[j] - i) - 1
    return math.factorial(n) // prod


def schur_1N(lam):
    conj = [sum(1 for r in lam if r > j) for j in range(lam[0])] if lam else []
    num, den = 1, 1
    for i, r in enumerate(lam):
        for j in range(r):
            num *= (N + (j + 1) - (i + 1))
            den *= ((r - j) + (conj[j] - i) - 1)
    return Rational(num, den)


def chi_frob(lam, ct):
    """chi^lam(class ct) via Frobenius: coefficient of x^(lam+delta) in Delta * p_ct, 2 variables."""
    m = 2
    lam = tuple(lam) + (0,) * (m - len(lam))
    # polynomial in x1,x2 as dict {(e1,e2): coeff}
    def mul(A, B):
        C = {}
        for k1, v1 in A.items():
            for k2, v2 in B.items():
                k = (k1[0] + k2[0], k1[1] + k2[1])
                C[k] = C.get(k, 0) + v1 * v2
        return C
    p = {(0, 0): 1}
    for l in ct:
        # p_l(x1,x2) = x1^l + x2^l
        p = mul(p, {(l, 0): 1, (0, l): 1})
    p = mul(p, {(0, 0): 1})  # noop
    # multiply by Delta = x1 - x2
    D = {}
    for k, v in p.items():
        D[(k[0] + 1, k[1])] = D.get((k[0] + 1, k[1]), 0) + v
        D[(k[0], k[1] + 1)] = D.get((k[0], k[1] + 1), 0) - v
    key = (lam[0] + 1, lam[1])          # delta = (1,0)
    return D.get(key, 0)


def naive_table(n):
    ctlist, _ = classes_of(n)
    vals = {}
    for lam in partitions(n):
        if len(lam) > N:
            continue
        f = hook_dim(lam)
        d = schur_1N(lam)
        if d == 0:
            continue
        for ct in ctlist:
            ch = chi_frob(lam, ct)
            if ch:
                vals[ct] = vals.get(ct, Rational(0)) + Rational(f * f) * Rational(ch) / (Rational(math.factorial(n)) ** 2 * d)
    return vals


def main():
    t0 = time.time()
    print("=== SU(2) Weingarten table, exact, n<=6  (S^3 hafnian method) ===", flush=True)
    table = {}
    for n in range(1, 7):
        ctlist, vals, neq = solve_wg(n)
        table[n] = dict(zip(ctlist, vals))
        print("\nn = %d : %d classes, n! = %d, %d independent-looking equations" %
              (n, len(ctlist), math.factorial(n), neq))
        for ct, v in zip(ctlist, vals):
            print("  %-10s  Wg = %-16s = %+.10f" % (cname(ct), v, float(v)))
        print("  [%.1fs]" % (time.time() - t0), flush=True)
    print("\ngiven-value check: Wg(2 at n=2) = %s  (must be -1/6)" % table[2][(2,)])
    print("given-value check: Wg(1^2)        = %s  (must be 1/3)" % table[2][(1, 1)])

    print("\n=== check (ii): Catalan moments from the table ===")
    cat = [1, 1, 2, 5, 14, 42, 132]
    for k in range(1, 7):
        ctlist, mult = classes_of(k)
        s = Rational(0)
        for ct in ctlist:
            s += mult[ct] * table[k][ct] * N ** len(ct)
        val = math.factorial(k) * s
        print("  k=%d  <|trU|^%2d> = %-6s (%.6f)   C_%d=%d  %s" %
              (k, 2 * k, val, float(val), k, cat[k], "OK" if val == cat[k] else "MISMATCH"))

    print("\n=== check (iii): Monte-Carlo (2e6 samples; table not used) ===")
    import numpy as np
    rng = np.random.default_rng(11)
    MS = 2_000_000
    A = rng.normal(size=(4, MS))
    A /= np.linalg.norm(A, axis=0)
    a0, a1, a2, a3 = A
    U = np.empty((MS, 2, 2), dtype=complex)
    U[:, 0, 0] = a0 + 1j * a3
    U[:, 0, 1] = a2 + 1j * a1
    U[:, 1, 0] = -a2 + 1j * a1
    U[:, 1, 1] = a0 - 1j * a3
    Ub = U.conj()
    for n in (2, 3, 4):
        ab = [0, 1, 0, 1, 0, 1][:n]
        cd = [1, 0, 1, 0, 1, 0][:n]
        ab2 = [1, 0, 1, 0, 1, 0][:n]
        cd2 = [0, 1, 0, 1, 0, 1][:n]
        prod = np.ones(MS, dtype=complex)
        for p in range(n):
            prod *= U[:, ab[p], cd[p]]
        for p in range(n):
            prod *= Ub[:, ab2[p], cd2[p]]
        mc = prod.mean()
        err = prod.std() / math.sqrt(MS)
        ex = integ_exact([mvec(ab[p], cd[p]) for p in range(n)] +
                         [mvec(ab2[p], cd2[p], True) for p in range(n)])
        print("  n=%d  MC = %+.5f +- %.5f   exact = %s (%.5f)  %s" %
              (n, mc.real, err.real, ex, float(ex), "OK" if abs(mc.real - float(ex)) < 5 * err.real else "CHECK"))

    print("\n=== check (iv): naive sum_{l(lam)<=N} f^2 chi/s_lam(1^N) vs exact ===")
    for n in range(1, 7):
        nv = naive_table(n)
        bad = [ct for ct in table[n] if nv.get(ct) != table[n][ct]]
        print("  n=%d : %d/%d classes differ" % (n, len(bad), len(table[n])))
        for ct in bad[:2]:
            print("        %-8s naive %-14s exact %s" % (cname(ct), nv.get(ct), table[n][ct]))
    print("\nelapsed %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
