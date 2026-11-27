#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GF3-ENTROPY witness: entropy constant mu of defect-cluster film classes via
transfer matrices + Perron-Frobenius; exact algebraic values; Sturm decision
of "mu <= e" with no enumeration.

Classes (2D films of plaquettes, counted up to TRANSLATION):
  D_s  : directed (one-sided / monotone base) no-overhang films, drop <= s   (s=0 = bargraph)
  D_inf: same, unbounded drop  (closed form; mu = phi^2)
  S_s  : two-sided no-overhang films, |base step| <= s        (column-convex truncation)
  CV   : column-convex (S_s with s >= n)                      [subclass of planar polyominoes]
  PLA  : planar films = all fixed polyominoes (brute force n<=11)

Encoding:
  a film = column sequence (h_1,d_1,h_2,d_2,...,d_{m-1},h_m), h_j>=1, d_j = b_{j+1}-b_j
  overlap  <=>  h_j >= max(1, d_j+1, 1-d_{j-1})  =: L(d_{j-1}, d_j)
  height summation is exact:  sum_{h>=L} x^h = x^L/(1-x)
  => finite automaton on d, matrix F(x) = E(x)/(1-x), E[d][d'] = x^{L(d,d')}
  => mu = 1/x_c,  x_c = smallest positive root of  det((1-x) I - E(x))  (integer poly)

All decisions use exact integer / Fraction arithmetic.  Floats appear only as
printed cross-checks (numpy spectral radius).
"""
import time
import math
from fractions import Fraction
from math import gcd
import decimal

t0 = time.time()
OUT = []


def say(*a):
    s = ' '.join(str(x) for x in a)
    OUT.append(s)
    print(s, flush=True)


# ----------------------------------------------------------------- poly over Z
ZERO = [0]
ONE = [1]
X1M = [1, -1]          # 1 - x


def ptrim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def padd(a, b):
    n = max(len(a), len(b))
    r = [0] * n
    for i in range(len(a)):
        r[i] += a[i]
    for i in range(len(b)):
        r[i] += b[i]
    return ptrim(r)


def psub(a, b):
    n = max(len(a), len(b))
    r = [0] * n
    for i in range(len(a)):
        r[i] += a[i]
    for i in range(len(b)):
        r[i] -= b[i]
    return ptrim(r)


def pmul(a, b):
    if a == ZERO or b == ZERO:
        return ZERO[:]
    r = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    r[i + j] += ai * bj
    return ptrim(r)


def pscal(a, c):
    if c == 0:
        return ZERO[:]
    return ptrim([c * x for x in a])


def mono(k):
    return ptrim([0] * k + [1])          # x^k


def peval(p, x):
    s = Fraction(0)
    for c in reversed(p):
        s = s * x + c
    return s


def pdiv_exact(a, b):
    a = ptrim(a)
    b = ptrim(b)
    if b == ZERO:
        raise ZeroDivisionError('poly div by 0')
    if len(a) < len(b):
        return ZERO[:]
    rem = [Fraction(x) for x in a]
    db = len(b) - 1
    q = [Fraction(0)] * (len(a) - len(b) + 1)
    for i in range(len(a) - len(b), -1, -1):
        c = rem[i + db] / Fraction(b[db])
        q[i] = c
        if c:
            for j in range(db + 1):
                rem[i + j] -= c * b[j]
    for r in rem:
        if r != 0:
            raise ZeroDivisionError('inexact poly division')
    out = []
    for c in q:
        if c.denominator != 1:
            raise ZeroDivisionError('non-integer quotient')
        out.append(int(c))
    return ptrim(out)


def det_poly(M):
    n = len(M)
    if n == 0:
        return ONE[:]
    A = [[ptrim(list(e)) for e in row] for row in M]
    prev = ONE[:]
    sign = 1
    for k in range(n - 1):
        if A[k][k] == ZERO:
            piv = -1
            for i in range(k + 1, n):
                if A[i][k] != ZERO:
                    piv = i
                    break
            if piv < 0:
                return ZERO[:]
            A[k], A[piv] = A[piv], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = pdiv_exact(
                    psub(pmul(A[k][k], A[i][j]), pmul(A[i][k], A[k][j])), prev)
        prev = A[k][k]
    return pscal(A[n - 1][n - 1], sign)


def det_perm(M):
    from itertools import permutations
    n = len(M)
    tot = ZERO[:]
    for perm in permutations(range(n)):
        sgn = 1
        p = list(perm)
        for i in range(n):
            for j in range(i + 1, n):
                if p[i] > p[j]:
                    sgn = -sgn
        term = ONE[:]
        for i in range(n):
            term = pmul(term, M[i][perm[i]])
        tot = padd(tot, pscal(term, sgn))
    return ptrim(tot)


# ------------------------------------------------------------------ Sturm
def prem(a, b):
    a = [Fraction(x) for x in ptrim(a)]
    b = ptrim(b)
    db = len(b) - 1
    while a != [Fraction(0)] and len(a) - 1 >= db:
        d = len(a) - 1 - db
        c = a[-1] / b[-1]
        for j in range(db + 1):
            a[d + j] -= c * b[j]
        while a and a[-1] == 0:
            a.pop()
        if not a:
            a = [Fraction(0)]
            break
    L = 1
    for x in a:
        L = L * x.denominator // gcd(L, x.denominator)
    return ptrim([int(x * L) for x in a])


def pderiv(p):
    p = ptrim(p)
    if len(p) <= 1:
        return ZERO[:]
    return ptrim([p[i] * i for i in range(1, len(p))])


def sturm_chain(p):
    p = ptrim(p)
    ch = [p, pderiv(p)]
    if ch[1] == ZERO:
        return [p]
    while True:
        r = prem(ch[-2], ch[-1])
        if r == ZERO:
            break
        ch.append(pscal(r, -1))
    return ch


def sign_variations(chain, x):
    v = 0
    last = 0
    for p in chain:
        s = peval(p, x)
        if s == 0:
            continue
        sg = 1 if s > 0 else -1
        if last != 0 and sg != last:
            v += 1
        last = sg
    return v


def count_roots(chain, a, b):
    return sign_variations(chain, a) - sign_variations(chain, b)


def isolate_smallest_positive(p, bits=64):
    ch = sturm_chain(p)
    hi = Fraction(1)
    while count_roots(ch, Fraction(0), hi) == 0:
        hi *= 2
        if hi > 4:
            return None, None, 0
    lo = Fraction(0)
    n0 = count_roots(ch, lo, hi)
    for _ in range(bits):
        mid = (lo + hi) / 2
        if count_roots(ch, lo, mid) >= 1:
            hi = mid
        else:
            lo = mid
    return lo, hi, n0


def inv_e_bounds(m=30):
    """1/e = sum_{k>=0} (-1)^k/k! with alternating partial sums => rigorous bounds"""
    lo = Fraction(0)
    hi = Fraction(1)
    term = Fraction(1)
    s = Fraction(0)
    for k in range(m + 1):
        if k > 0:
            term = term / k
        s += term if (k % 2 == 0) else -term
        if k % 2 == 0:
            hi = s
        else:
            lo = s
    return lo, hi


# ------------------------------------------------------------ class matrices
def E_matrix(D, ycost=None):
    n = len(D)
    smax = max(abs(z) for z in D)
    M = [[None] * n for _ in range(n)]
    for i, d in enumerate(D):
        for j, dp in enumerate(D):
            L = max(1, dp + 1, 1 - d)
            e = mono(L)
            if ycost is not None:
                p, q = ycost
                e = pscal(e, (p ** abs(dp)) * (q ** (smax - abs(dp))))
            M[i][j] = e
    return M


def charpoly(D, ycost=None):
    """det of  (1-x)I - E(x)  (integer poly).

    ycost=(p,q) attaches the Boltzmann weight y^|d'| per transition with
    y = p/q (y<1 suppresses interface slope).  The WHOLE matrix is multiplied
    by the common denominator `den` (a positive integer), which multiplies the
    determinant by den^n and leaves its roots unchanged.
    """
    n = len(D)
    Lm = [[max(1, dp + 1, 1 - d) for dp in D] for d in D]
    if ycost is None:
        M = [[psub(X1M if i == j else ZERO, mono(Lm[i][j])) for j in range(n)]
             for i in range(n)]
        return det_poly(M)
    p, q = ycost
    S = [[Fraction(p, q) ** abs(D[j]) for j in range(n)] for i in range(n)]
    den = 1
    for i in range(n):
        for j in range(n):
            den = den * S[i][j].denominator // gcd(den, S[i][j].denominator)
    M = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            c = int(S[i][j] * den)
            M[i][j] = psub(pscal(X1M, den) if i == j else ZERO,
                           pscal(mono(Lm[i][j]), c))
    return det_poly(M)


# ------------------------------------------------------- exact DP counts
def dp_counts(D, N):
    """a_1..a_N : number of class members of size n, counted up to translation.

    state = current column height h only: the incoming offset constraint is
    enforced when the column is created (h >= max(1,1-d_prev)); the outgoing
    offset then only constrains the NEXT height (h_next >= max(1,1-d)) and the
    current one (h >= d+1).  Hence: state (h), transitions over d in D.
    Every prefix is a complete shape (its last column is the last column), so
    a_n = sum_h states_n[h].  Aggregated => no path explosion.
    """
    st = [dict() for _ in range(N + 1)]
    for h in range(1, N + 1):
        st[h][h] = 1                      # single column of height h
    for s in range(1, N + 1):
        if not st[s]:
            continue
        for h, c in list(st[s].items()):
            for d in D:
                if h < d + 1:             # overlap with the next column fails
                    continue
                L = max(1, 1 - d)
                for hp in range(L, N - s + 1):
                    tgt = st[s + hp]
                    tgt[hp] = tgt.get(hp, 0) + c
    a = [0] * (N + 1)
    for n in range(1, N + 1):
        a[n] = sum(st[n].values())
    return a


def gf_series(D, N):
    """expand u^T (I-F)^-1 w + x/(1-x) to x^N, exact integers"""
    n = len(D)

    def ser(L):
        return [0] * L + [1] * (N - L + 1)

    F = [[ser(max(1, dp + 1, 1 - d)) for dp in D] for d in D]
    u = [ser(max(1, d + 1)) for d in D]
    w = [ser(max(1, 1 - d)) for d in D]

    def mulser(A, B):
        r = [0] * (N + 1)
        for i, ai in enumerate(A):
            if ai:
                m = min(len(B) - 1, N - i)
                for j in range(m + 1):
                    if B[j]:
                        r[i + j] += ai * B[j]
        return r

    acc = [0] * (N + 1)
    v = [row[:] for row in w]
    for _ in range(N + 1):
        for i in range(n):
            acc = [x + y for x, y in zip(acc, mulser(u[i], v[i]))]
        vnew = [[0] * (N + 1) for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if any(F[i][j]):
                    vnew[i] = [x + y for x, y in zip(vnew[i], mulser(F[i][j], v[j]))]
        v = vnew
    out = acc[:]
    for k in range(1, N + 1):
        out[k] += 1
    return out


def rat_series(coef_num, coef_den, N):
    """expand a rational GF num(x)/den(x) (integer coeff lists) to x^N"""
    def pediv(a, b, k):
        pass
    d0 = coef_den[0]
    assert d0 != 0
    out = [0] * (N + 1)
    num = [Fraction(c, d0) for c in coef_num]
    den = [Fraction(c, d0) for c in coef_den]
    for k in range(N + 1):
        s = num[k] if k < len(num) else Fraction(0)
        for j in range(1, min(k, len(den) - 1) + 1):
            s -= den[j] * out[k - j]
        out[k] = s
    return [int(x) if x.denominator == 1 else x for x in out]


# ---------------------------------------------------- polyomino brute force
def enum_polyominoes(N):
    shapes = {1: {((0, 0),)}}
    cur = {((0, 0),)}
    for n in range(2, N + 1):
        nxt = set()
        for sh in cur:
            cells = set(sh)
            cand = set()
            for (x, y) in cells:
                for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    p = (x + dx, y + dy)
                    if p not in cells:
                        cand.add(p)
            for p in cand:
                s2 = cells | {p}
                mnx = min(c[0] for c in s2)
                mny = min(c[1] for c in s2)
                nxt.add(tuple(sorted((c[0] - mnx, c[1] - mny) for c in s2)))
        cur = nxt
        shapes[n] = cur
    return shapes


def cols_of(sh):
    cols = {}
    for (x, y) in sh:
        cols.setdefault(x, []).append(y)
    return cols


def is_colconvex(cols):
    for x in cols:
        ys = sorted(cols[x])
        if ys[-1] - ys[0] + 1 != len(ys):
            return False
    return True


def is_directed_down(cols):
    """column-convex + base (min row) non-increasing left->right"""
    if not is_colconvex(cols):
        return False
    xs = sorted(cols)
    b = [min(cols[x]) for x in xs]
    return all(b[i + 1] <= b[i] for i in range(len(b) - 1))


# ======================================================================= main
def main():
    N = 14
    say('=' * 78)
    say('GF3-ENTROPY  |  transfer matrix + Perron-Frobenius for the entropy constant mu')
    say('=' * 78)

    say('')
    say('[A] poly machinery self-test (Bareiss vs permutation expansion)')
    tests = [
        [[mono(1), psub(X1M, mono(2))], [psub(X1M, ZERO), mono(1)]],
        [[mono(1), mono(2), X1M],
         [mono(3), psub(X1M, mono(1)), ZERO],
         [mono(0), mono(2), mono(1)]],
        [[mono(1), psub(X1M, mono(1)), mono(2), ZERO],
         [mono(0), mono(2), psub(X1M, ZERO), mono(3)],
         [mono(4), ZERO, mono(1), mono(1)],
         [psub(X1M, mono(2)), mono(1), ZERO, mono(0)]],
    ]
    allok = True
    for M in tests:
        d1 = det_poly(M)
        d2 = det_perm(M)
        allok = allok and (d1 == d2)
    say('    det_poly == det_perm : %s' % ('PASS' if allok else 'FAIL'))

    say('')
    say('    t=%.1fs' % (time.time() - t0))
    say('[B] classes and exact DP counts (n <= %d)' % N)
    classes = {
        'D_0 bargraph': [0],
        'D_1 directed drop<=1': [0, -1],
        'D_2': [0, -1, -2],
        'D_3': [0, -1, -2, -3],
        'D_4': [0, -1, -2, -3, -4],
        'S_1 two-sided |st|<=1': [-1, 0, 1],
        'S_2': [-2, -1, 0, 1, 2],
        'S_3': [-3, -2, -1, 0, 1, 2, 3],
        'S_4': list(range(-4, 5)),
        'D_inf (unbounded drop)': list(range(-N, 1)),
        'CV column-convex': list(range(-N, N + 1)),
    }
    counts = {}
    for name, D in classes.items():
        counts[name] = dp_counts(D, N)
        say('    %-24s states=%2d  a_1..a_%d = %s' %
            (name, len(D), N, counts[name][1:]))

    say('')
    say('    t=%.1fs' % (time.time() - t0))
    say('[C] cross-check: transfer-matrix GF series == DP counts (exact ints)')
    for name, D in classes.items():
        if len(D) > 9:
            continue
        g = gf_series(D, N)
        same = (g[1:] == counts[name][1:])
        say('    %-24s %s' % (name, 'MATCH' if same else 'MISMATCH'))

    say('')
    say('[D] closed form check: a_n(D_inf) = Fibonacci(2n-1) ?')
    fib = [0, 1]
    while len(fib) < 2 * N + 3:
        fib.append(fib[-1] + fib[-2])
    pred = [fib[2 * n - 1] for n in range(1, N + 1)]
    got = counts['D_inf (unbounded drop)'][1:]
    say('    F_(2n-1)      : %s' % pred)
    say('    DP D_inf      : %s' % got)
    say('    equality      : %s' % ('PASS' if pred == got else 'FAIL'))
    # rational GF  x(1-x)/(1-3x+x^2)
    rs = rat_series([0, 1, -1], [1, -3, 1], N)
    say('    GF x(1-x)/(1-3x+x^2) series: %s' % rs[1:])
    say('    matches DP    : %s' % ('PASS' if [int(x) for x in rs[1:]] == got else 'FAIL'))

    say('')
    say('    t=%.1fs' % (time.time() - t0))
    say('[E] brute-force polyominoes (independent enumeration)')
    MAXB = 11
    shapes = enum_polyominoes(MAXB)
    plan, cc, dd = {}, {}, {}
    for n in range(1, MAXB + 1):
        plan[n] = len(shapes[n])
        cc[n] = 0
        dd[n] = 0
        for sh in shapes[n]:
            cols = cols_of(sh)
            if is_colconvex(cols):
                cc[n] += 1
                if is_directed_down(cols):
                    dd[n] += 1
    rng = range(1, MAXB + 1)
    say('    n            : ' + ' '.join('%6d' % n for n in rng))
    say('    planar A_n   : ' + ' '.join('%6d' % plan[n] for n in rng))
    say('    col-convex   : ' + ' '.join('%6d' % cc[n] for n in rng))
    say('    directed(dn) : ' + ' '.join('%6d' % dd[n] for n in rng))
    say('    DP CV  == brute col-convex : %s' %
        ('PASS' if [counts['CV column-convex'][n] for n in rng] == [cc[n] for n in rng] else 'FAIL'))
    say('    DP Dinf== brute directed   : %s' %
        ('PASS' if [counts['D_inf (unbounded drop)'][n] for n in rng] == [dd[n] for n in rng] else 'FAIL'))

    say('')
    say('    t=%.1fs' % (time.time() - t0))
    say('[F] char polys P(x) = det((1-x)I - E(x)) and exact smallest positive roots')
    order = ['D_0 bargraph', 'D_1 directed drop<=1', 'D_2', 'D_3', 'D_4',
             'S_1 two-sided |st|<=1', 'S_2', 'S_3', 'S_4']
    table = []
    for name in order:
        D = classes[name]
        P = charpoly(D)
        lo, hi, nr = isolate_smallest_positive(P)
        table.append((name, len(D), P, lo, hi))
        say('    %-24s states=%d deg=%2d  P = %s' % (name, len(D), len(P) - 1, P))
        say('        x_c in (%.20g, %.20g]   mu = 1/x_c in [%.15g, %.15g)' %
            (float(lo), float(hi), float(1 / hi), float(1 / lo)))

    say('')
    say('[G] decision mu <= e (Sturm brackets of x_c vs rational brackets of 1/e)')
    elo, ehi = inv_e_bounds(30)
    decimal.getcontext().prec = 60
    say('    1/e in (%.25f, %.25f)  width=%.2e' % (float(elo), float(ehi), float(ehi - elo)))
    say('    %-24s %-20s %-20s %s' % ('class', 'mu lower', 'mu upper', 'mu <= e ?'))
    for (name, nst, P, lo, hi) in table:
        if lo > ehi:
            v = 'YES'
        elif hi < elo:
            v = 'NO'
        else:
            v = 'UNDECIDED'
        say('    %-24s %-20.15f %-20.15f %s' % (name, float(1 / hi), float(1 / lo), v))
    xc_inf = (3 - decimal.Decimal(5).sqrt()) / 2
    mu_inf = 1 / xc_inf
    say('    %-24s %-20s %-20s %s' % ('D_inf', 'phi^2 = ' + str(mu_inf)[:18], '(exact)',
                                    'YES' if Fraction(str(xc_inf)) > ehi else 'NO'))
    say('    phi^2  = %s' % str(mu_inf)[:40])
    say('    e      = %s' % str(decimal.Decimal(1).exp())[:40])

    say('')
    say('[H] spectral-radius cross-check (float): rho(F(x_c)) must equal 1')
    import numpy as np
    for (name, nst, P, lo, hi) in table:
        D = classes[name]
        n = len(D)
        x = float((lo + hi) / 2)
        F = np.zeros((n, n))
        for i, d in enumerate(D):
            for j, dp in enumerate(D):
                F[i, j] = x ** max(1, dp + 1, 1 - d)
        F = F / (1 - x)
        rho = float(np.max(np.abs(np.linalg.eigvals(F))))
        say('    %-24s x_c=%.15f  rho(F)=%.15f  |rho-1|=%.2e' %
            (name, x, rho, abs(rho - 1)))

    say('')
    say('[I] one-sided sequence D_s -> phi^2 (all must satisfy mu <= e)')
    for s in [0, 1, 2, 3, 4, 6, 8, 10, 12, 14]:
        D = list(range(-s, 0 + 1))
        P = charpoly(D)
        lo, hi, nr = isolate_smallest_positive(P, bits=48)
        say('    s=%2d states=%2d  mu in [%.12f, %.12f)   limit phi^2 = %.12f' %
            (s, len(D), float(1 / hi), float(1 / lo), float(mu_inf)))

    say('')
    say('[J] two-sided S_1 as lower bound for the planar / column-convex class')
    for name in ['S_1 two-sided |st|<=1', 'S_2', 'S_3', 'S_4']:
        t = [z for z in table if z[0] == name][0]
        mu_lo = 1 / t[4]
        say('    %-24s states=%d  mu >= %.12f   > e ? %s' %
            (name, t[1], float(mu_lo), 'YES' if mu_lo > ehi else 'NO'))
    say('    planar mu >= A_10^(1/10) = %.12f' % (plan[10] ** (1 / 10)))
    say('    planar mu >= A_11^(1/11) = %.12f' % (plan[11] ** (1 / 11)))

    say('')
    say('[K] interface-slope cost: transition d carries weight y^|d| (y<1 suppresses)')
    say('    mu(y) <= e  <=>  y <= y*   ;  cost per unit slope kappa* = -ln y*')
    for s in [1, 2, 3]:
        D = list(range(-s, s + 1))

        def ok(y):
            P = charpoly(D, ycost=(y.numerator, y.denominator))
            lo, hi, nr = isolate_smallest_positive(P, bits=40)
            if lo is None:
                return None
            if lo > ehi:
                return True         # x_c > 1/e  =>  mu < e
            if hi < elo:
                return False
            return None

        f1 = ok(Fraction(1))
        f0 = ok(Fraction(1, 16))
        say('    s=%d : mu(1)<=e ? %s ; mu(1/16)<=e ? %s' % (s, f1, f0))
        if f1 is False and f0 is True:
            loY, hiY = Fraction(1, 16), Fraction(1)
            for _ in range(22):
                mid = (loY + hiY) / 2
                r = ok(mid)
                if r is True:
                    loY = mid
                else:
                    hiY = mid
            say('    s=%d : y* in [%.10f, %.10f]  =>  kappa* = -ln y* in [%.6f, %.6f]' %
                (s, float(loY), float(hiY), -math.log(float(hiY)), -math.log(float(loY))))

    say('')
    say('[M] closed form check (one-sided family): P_s == (1-x)^s * (1-2x-x^2-...-x^{s+1})')
    cfok = True
    for s in range(0, 7):
        D = list(range(-s, 1))
        P = charpoly(D)
        Q = ptrim([1, -2] + [-1] * s)
        R = ONE[:]
        for _ in range(s):
            R = pmul(R, X1M)
        R = pmul(R, Q)
        same = (P == R) or (P == pscal(R, -1))
        cfok = cfok and same
        say('    s=%d P=%s' % (s, P))
        say('        (1-x)^%d*(1-2x-...-x^%d)=%s  %s' % (s, s + 1, R, 'OK' if same else 'FAIL'))
    say('    s=0..6 all %s ; hence x_s solves x+x^2+...+x^{s+1}=1-x and mu(D_s)=1/x_s' %
        ('PASS' if cfok else 'FAIL'))

    say('')
    say('[L] runtime %.1f s' % (time.time() - t0))
    with open('gf3_entropy.log', 'w', encoding='utf-8') as f:
        f.write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
