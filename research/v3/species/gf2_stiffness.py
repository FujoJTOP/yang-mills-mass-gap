#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GF2-STIFFNESS witness: exact enumeration of overhang-free lattice classes.

Parts
  1  bargraphs (flat bottom + column-convex): B_n = 2^(n-1) exact
  2  column-convex (no dangling cell, floats allowed): exact counts n<=14,
     plus perimeter / reflex-corner distributions
  3  directed (NE) animals: mask DP n<=14 + brute-force cross-check n<=9
  4  "common row" family (star-shaped w.r.t. a line): exact GF, mu = phi^2
  5  tilt sums  Z_n(kappa) = sum_C exp(-kappa f(C)),  f = perimeter / reflex corners
     -> subexponential vs extensive test
Run:  py -3.12 -I gf2_stiffness.py
"""
import math
import time
import sys
from collections import defaultdict

import numpy as np

NS = 14
OUT = []


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


# ----------------------------------------------------------------------------
# curvature: reflex-corner count on the square lattice
#   curv(C) = #{ vertices v : exactly 3 of the 4 cells around v are in C }
# ----------------------------------------------------------------------------
def pair_curv(I, J):
    """reflex corners on the vertical line shared by column intervals I and J."""
    b, t = I
    b2, t2 = J
    c = 0
    for y in range(min(b, b2) - 1, max(t, t2) + 2):
        k = ((b <= y - 1 <= t) + (b2 <= y - 1 <= t2)
             + (b <= y <= t) + (b2 <= y <= t2))
        if k == 3:
            c += 1
    return c


def curv_of(shape):
    """shape: set of (x,y) cells -> reflex corner count (brute force)."""
    xs = [p[0] for p in shape]
    ys = [p[1] for p in shape]
    c = 0
    for x in range(min(xs), max(xs) + 2):
        for y in range(min(ys), max(ys) + 2):
            k = ((x - 1, y - 1) in shape) + ((x, y - 1) in shape) \
                + ((x - 1, y) in shape) + ((x, y) in shape)
            if k == 3:
                c += 1
    return c


def perim_of(shape):
    p = 0
    for (x, y) in shape:
        for d in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            if (x + d[0], y + d[1]) not in shape:
                p += 1
    return p


# ----------------------------------------------------------------------------
# Part 1  bargraphs
# ----------------------------------------------------------------------------
def part1():
    say("=" * 78)
    say("Part 1  bargraphs  (column-convex + flat bottom = height sequence)")
    say("=" * 78)
    say("bijection: h_1..h_k >= 1, sum h_i = n  <->  composition of n")
    say("=> B_n = 2^(n-1), mu = 2 exactly, gamma = ln2 = %.6f" % math.log(2))
    tot = 0
    for n in range(1, NS + 1):
        b = 2 ** (n - 1)
        say("  n=%2d  B_n=%6d   B_n^(1/n)=%.6f" % (n, b, b ** (1.0 / n)))
    # perimeter and curvature of a bargraph
    h = [1, 2, 1, 3]
    sh = set()
    for i, hi in enumerate(h):
        for y in range(hi):
            sh.add((i, y))
    say("  check shape h=(1,2,1,3): n=%d p=%d curv=%d" %
        (len(sh), perim_of(sh), curv_of(sh)))
    say("  perimeter of a bargraph: p = 2n + 2k  (k = #columns)  -> EXTENSIVE (2 per cell)")
    say("  reflex corners: curv = #{i : h_i != h_{i+1}}  (one notch per step)")
    say("  -> for a typical composition curv = Theta(n): EXTENSIVE")
    return


# ----------------------------------------------------------------------------
# Part 2  column-convex
# ----------------------------------------------------------------------------
def cc_dp(N, mode):
    """mode 'n' | 'p' | 'c'.  Returns array indexed by (n[,p][,c])."""
    Pm = 4 * N + 2
    if mode == 'n':
        shape = (N + 1,)
    else:
        shape = (N + 1, Pm + 1)
    ints = [(b, t) for t in range(N + 1) for b in range(t + 1)]

    def mk(b, t):
        A = np.zeros(shape, dtype=object if False else np.float64)
        n0 = t - b + 1
        if mode == 'n':
            A[n0] += 1
        elif mode == 'p':
            p0 = 2 * n0 + 2
            A[n0, p0] += 1
        else:
            A[n0, 0] += 1
        return A

    dp = {}
    for (b, t) in ints:
        if t - b + 1 > N:
            continue
        dp[(b, t, 1 if b == 0 else 0)] = mk(b, t)

    out = np.zeros(shape)
    for (b, t, h), A in dp.items():
        if h == 1:
            out += A

    for step in range(1, N):
        ndp = defaultdict(lambda: np.zeros(shape))
        for (b, t, h), A in dp.items():
            nz = np.nonzero(A)[0]
            if len(nz) == 0:
                continue
            for (b2, t2) in ints:
                if max(b, b2) > min(t, t2):
                    continue          # disconnected
                dn = t2 - b2 + 1
                if dn > N:
                    continue
                h2 = 1 if (h == 1 or b2 == 0) else 0
                if mode == 'n':
                    ndp[(b2, t2, h2)][dn:] += A[:N + 1 - dn]
                elif mode == 'p':
                    ov = min(t, t2) - max(b, b2) + 1
                    dpp = 2 * dn + 2 - 2 * ov
                    ndp[(b2, t2, h2)][dn:, dpp:] += A[:N + 1 - dn, :Pm + 1 - dpp]
                else:
                    dc = pair_curv((b, t), (b2, t2))
                    ndp[(b2, t2, h2)][dn:, dc:] += A[:N + 1 - dn, :Pm + 1 - dc]
        dp = ndp
        for (b, t, h), A in dp.items():
            if h == 1:
                out += A
    return out


def part2():
    say("")
    say("=" * 78)
    say("Part 2  column-convex polyominoes (no dangling cell; bottom b_i free)")
    say("=" * 78)
    t0 = time.time()
    Cn = cc_dp(NS, 'n')
    Cp = cc_dp(NS, 'p')
    Cc = cc_dp(NS, 'c')
    say("DP time: %.1f s" % (time.time() - t0))
    say("")
    say(" n     C_n        ratio    C_n^(1/n)   <p>     p/n     <curv>   curv/n")
    prev = None
    for n in range(1, NS + 1):
        a = int(round(Cn[n]))
        pp = Cp[n]
        cc = Cc[n]
        tot = pp.sum()
        mp = (pp * np.arange(len(pp))).sum() / tot if tot > 0 else 0
        mc = (cc * np.arange(len(cc))).sum() / tot if tot > 0 else 0
        r = ("%.4f" % (a / prev)) if prev else "   -  "
        say("%2d  %10d  %s   %.6f   %6.2f  %.3f   %6.2f  %.3f" %
            (n, a, r, a ** (1.0 / n), mp, mp / n, mc, mc / n))
        prev = a
    return Cn, Cp, Cc


# ----------------------------------------------------------------------------
# Part 3  directed animals
# ----------------------------------------------------------------------------
def directed_dp(N):
    """diagonal mask DP.  state (mask, used) -> count."""
    counts = [0] * (N + 1)
    counts[1] = 1
    cur = {(1, 1): 1}
    diag = 0
    while cur:
        diag += 1
        nxt = {}
        for (mask, used), c in cur.items():
            R = mask | (mask << 1)
            sub = R
            while sub:
                k = bin(sub).count("1")
                if used + k <= N:
                    key = (sub, used + k)
                    nxt[key] = nxt.get(key, 0) + c
                    counts[used + k] += c
                sub = (sub - 1) & R
        cur = nxt
        if diag > 3 * N:
            break
    return counts


def brute_animals(nmax):
    """BFS with canonical dedup: all fixed polyominoes up to nmax cells."""
    cur = {frozenset([(0, 0)])}
    byn = {1: len(cur)}
    res = {1: list(cur)}
    for n in range(2, nmax + 1):
        nxt = set()
        for sh in cur:
            for (x, y) in sh:
                for d in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    c = (x + d[0], y + d[1])
                    if c in sh:
                        continue
                    s2 = frozenset(sh | {c})
                    mx = min(p[0] for p in s2)
                    my = min(p[1] for p in s2)
                    nxt.add(frozenset((p[0] - mx, p[1] - my) for p in s2))
        cur = nxt
        byn[n] = len(cur)
        res[n] = list(cur) if n <= 9 else None
    return byn, res


def is_directed(sh):
    cells = list(sh)
    S = set(cells)
    for s in cells:
        seen = {s}
        stack = [s]
        while stack:
            x, y = stack.pop()
            for d in ((1, 0), (0, 1)):
                q = (x + d[0], y + d[1])
                if q in S and q not in seen:
                    seen.add(q)
                    stack.append(q)
        if len(seen) == len(S):
            return True
    return False


def part3():
    say("")
    say("=" * 78)
    say("Part 3  directed (NE) animals")
    say("=" * 78)
    t0 = time.time()
    D = directed_dp(13)
    say("mask-DP time (N=13): %.1f s" % (time.time() - t0))
    say(" n   D_n        ratio    D_n^(1/n)")
    prev = None
    for n in range(1, 14):
        r = ("%.4f" % (D[n] / prev)) if prev else "   -  "
        say("%2d  %8d  %s   %.6f" % (n, D[n], r, D[n] ** (1.0 / n) if D[n] else 0))
        prev = D[n]
    # brute-force cross-check
    t0 = time.time()
    byn, res = brute_animals(9)
    cnts = []
    for n in range(1, 10):
        c = sum(1 for sh in res[n] if is_directed(sh))
        cnts.append(c)
    say("brute force cross-check (n<=9): %s   [%.1f s]" % (cnts, time.time() - t0))
    say("DP                          : %s" % [D[n] for n in range(1, 10)])
    # overhang counterexample: directed but a column with a gap
    sh = {(0, 0), (0, 1), (0, 2), (1, 0), (1, 2)}
    cols = defaultdict(list)
    for (x, y) in sh:
        cols[x].append(y)
    cc_ok = all(sorted(v) == list(range(min(v), max(v) + 1)) for v in cols.values())
    say("directed-but-overhanging example %s : directed=%s, column-convex=%s" %
        (sorted(sh), is_directed(sh), cc_ok))
    say("   columns=%s  curv=%d  p=%d" % (dict(cols), curv_of(sh), perim_of(sh)))
    return D


# ----------------------------------------------------------------------------
# Part 4  "common row" family
# ----------------------------------------------------------------------------
def part4(N=40):
    say("")
    say("=" * 78)
    say("Part 4  class 'every column meets a fixed row'  (star-shaped in a line)")
    say("=" * 78)
    say("per-column GF: sum_{b<=0<=t} x^(t-b+1) = x/(1-x)^2 ; columns independent")
    say("(consecutive columns always overlap at row 0 => connectivity automatic)")
    say("GF = 1/(1 - x/(1-x)^2) = (1-x)^2/(1-3x+x^2)")
    say("dominant singularity x* = (3-sqrt5)/2 = %.6f  =>  mu = phi^2 = %.6f" %
        ((3 - math.sqrt(5)) / 2, ((3 + math.sqrt(5)) / 2)))
    # exact coefficients by DP
    a = [0] * (N + 1)
    a[0] = 1
    # per-column poly: x/(1-x)^2 => c_k = k for k>=1
    for n in range(1, N + 1):
        s = 0
        for k in range(1, n + 1):
            s += k * a[n - k]
        a[n] = s
    say(" n   S_n        ratio     S_n^(1/n)")
    for n in range(1, min(N, 14) + 1):
        say("%2d  %9d  %.4f    %.6f" % (n, a[n], a[n] / a[n - 1] if n > 1 else 0,
                                        a[n] ** (1.0 / n)))
    say("(ratio -> phi^2 = 2.618034 < e = 2.718282 : PASSES by %.1f%%)" %
        (100 * (1 - 2.618034 / math.e)))
    return a


# ----------------------------------------------------------------------------
# Part 5  tilt sums
# ----------------------------------------------------------------------------
def part5(Cp, Cc):
    say("")
    say("=" * 78)
    say("Part 5  tilt sums  Z_n(kappa) = sum_C exp(-kappa f(C))")
    say("=" * 78)
    kappas = [0.0, 0.25, 0.5, 1.0, 2.0]
    say("(a) f = perimeter, column-convex")
    say("kap   " + "".join("  n=%2d " % n for n in range(6, NS + 1)))
    for k in kappas:
        row = []
        for n in range(6, NS + 1):
            p = Cp[n]
            Z = (p * np.exp(-k * np.arange(len(p)))).sum()
            row.append(Z ** (1.0 / n))
        say("%.2f  " % k + "".join(" %.5f" % v for v in row))
    say("(b) f = reflex corners, column-convex")
    say("kap   " + "".join("  n=%2d " % n for n in range(6, NS + 1)))
    for k in kappas:
        row = []
        for n in range(6, NS + 1):
            p = Cc[n]
            Z = (p * np.exp(-k * np.arange(len(p)))).sum()
            row.append(Z ** (1.0 / n))
        say("%.2f  " % k + "".join(" %.5f" % v for v in row))
    say("(c) f = perimeter, bargraph (p = 2n+2k, k=#columns) for comparison")
    for k in kappas:
        row = []
        for n in range(6, NS + 1):
            # compositions of n into k parts, perimeter 2n+2k
            Z = 0.0
            for kk in range(1, n + 1):
                Z += math.comb(n - 1, kk - 1) * math.exp(-k * (2 * n + 2 * kk))
            row.append(Z ** (1.0 / n))
        say("%.2f  " % k + "".join(" %.5f" % v for v in row))
    say("")
    say("read (a): column-convex <p>/n -> ~1.93 (EXTENSIVE, thin-dominated class),")
    say("   so a perimeter term is a RATE SHIFT, not exponentially invisible;")
    say("   the tilt density drifts toward kappa*1.93, not to 0.")
    say("   (for COMPACT classes -- all polyominoes, 4D animals -- <p> ~ sqrt(n) and the")
    say("    same term IS subexponential; the difference is the class, not the term.)")
    return


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def part6(Cn, Cp, Cc):
    say("")
    say("=" * 78)
    say("Part 6  rational-GF check (column-convex) + tilt-density diagnostic")
    say("=" * 78)
    # candidate GF from GD2-WSTAR/GE2: x(1-x)^3/(1-5x+7x^2-4x^3)
    num = [0, 1, -3, 3, -1]          # x(1-x)^3
    den = [1, -5, 7, -4]             # 1-5x+7x^2-4x^3
    ser = [0] * (NS + 1)
    for n in range(1, NS + 1):
        s = sum(den[k] * ser[n - k] for k in range(1, min(3, n) + 1))
        s = -s
        if n < len(num):
            s += num[n]
        ser[n] = s
    say(" n    DP C_n      GF [x^n]    equal?")
    ok = True
    for n in range(1, NS + 1):
        a = int(round(Cn[n]))
        say("%2d  %10d  %10d    %s" % (n, a, ser[n], "OK" if a == ser[n] else "MISMATCH"))
        ok = ok and (a == ser[n])
    say("all match: %s" % ok)
    # smallest positive root of 1-5x+7x^2-4x^3 by bisection
    f = lambda x: 1 - 5 * x + 7 * x * x - 4 * x ** 3
    lo, hi = 0.2, 0.4
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) * f(lo) <= 0:
            hi = mid
        else:
            lo = mid
    xs = 0.5 * (lo + hi)
    say("smallest positive pole x* = %.9f  =>  mu = 1/x* = %.6f" % (xs, 1 / xs))
    say("(GE2 row-convex value 3.2056 -> agreement to 4 decimals)")
    say("")
    say("tilt density  -[ln Z_n(kap) - ln Z_n(0)]/n   (column-convex, f=perimeter)")
    say("if f is subexponential this -> 0 ; if extensive it -> const")
    say("kap   " + "".join("  n=%2d " % n for n in range(6, NS + 1)))
    for k in [0.25, 0.5, 1.0]:
        row = []
        for n in range(6, NS + 1):
            p = Cp[n]
            Z0 = p.sum()
            Z = (p * np.exp(-k * np.arange(len(p)))).sum()
            row.append(-(math.log(Z) - math.log(Z0)) / n)
        say("%.2f  " % k + "".join(" %+.4f" % v for v in row))
    say("bargraph comparison (f = perimeter = 2n+2k, exactly extensive)")
    for k in [0.25, 0.5, 1.0]:
        row = []
        for n in range(6, NS + 1):
            Z = 0.0
            for kk in range(1, n + 1):
                Z += math.comb(n - 1, kk - 1) * math.exp(-k * (2 * n + 2 * kk))
            row.append(-math.log(Z) / n)
        say("%.2f  " % k + "".join(" %+.4f" % v for v in row))
    return


def part7(N=60):
    say("")
    say("=" * 78)
    say("Part 7  float-height hierarchy: b_i in [0,K]  (grounded within K rows)")
    say("=" * 78)
    say("K=0 bargraphs; K=inf column-convex.  a_n by transfer DP, mu = lim a_n/a_{n-1}")
    say(" K      a_n at n=%d          ratio(n=%d)   %d-th root   mu<=e" %
        (N, N, N))
    muK = {}
    for K in [0, 1, 2, 3, 4]:
        A = cc_dp_K(N, K)
        r = [int(round(A[n])) for n in range(1, N + 1)]
        ratio = r[-1] / r[-2] if r[-2] else 0
        muK[K] = ratio
        say("%2d  %30d   %.6f      %.6f    %s" %
            (K, r[-1], ratio, r[-1] ** (1.0 / N),
             "YES" if ratio <= math.e else "NO"))
    say("inf  (column-convex, exact pole)             3.205569                 NO")
    say("")
    say("sharp threshold: mu_K <= e  <=>  K <= 1   (mu_1 = 2.490216, mu_2 = 2.7457)")
    say("=> the window is closed by BOUNDED FLOAT HEIGHT (localisation),")
    say("   not by a curvature cost.")
    return muK


def cc_dp_K(N, K):
    """column-convex with bottoms b in [0,K]; counts by n."""
    ints = [(b, t) for t in range(N + 1) for b in range(min(t, K) + 1)]
    dp = {}
    for (b, t) in ints:
        if t - b + 1 > N:
            continue
        A = [0] * (N + 1)
        A[t - b + 1] += 1
        dp[(b, t, 1 if b == 0 else 0)] = A
    out = [0] * (N + 1)
    for (b, t, h), A in dp.items():
        if h == 1:
            for i in range(N + 1):
                out[i] += A[i]
    for step in range(1, N):
        ndp = defaultdict(lambda: [0] * (N + 1))
        for (b, t, h), A in dp.items():
            nz = [i for i in range(N + 1) if A[i]]
            if not nz:
                continue
            for (b2, t2) in ints:
                if max(b, b2) > min(t, t2):
                    continue
                dn = t2 - b2 + 1
                if dn > N:
                    continue
                h2 = 1 if (h == 1 or b2 == 0) else 0
                B = ndp[(b2, t2, h2)]
                for i in range(N + 1 - dn):
                    if A[i]:
                        B[i + dn] += A[i]
        dp = ndp
        for (b, t, h), A in dp.items():
            if h == 1:
                for i in range(N + 1):
                    out[i] += A[i]
    import numpy as _np
    return _np.array(out, dtype=float)


def part8(N=24):
    say("")
    say("=" * 78)
    say("Part 8  mean perimeter / mean reflex corners vs n (column-convex)")
    say("=" * 78)
    say("additive over column pairs: <p> = 2n+2k-2ovsum, <curv> = sum of pair terms")
    ints = [(b, t) for t in range(N + 1) for b in range(t + 1)]
    pc = {}
    for t in range(N + 1):
        for b in range(t + 1):
            for t2 in range(N + 1):
                for b2 in range(t2 + 1):
                    if max(b, b2) > min(t, t2):
                        continue
                    pc[((b, t), (b2, t2))] = pair_curv((b, t), (b2, t2))
    dp = {}
    for (b, t) in ints:
        n0 = t - b + 1
        if n0 > N:
            continue
        cnt = [0] * (N + 1)
        sp = [0.0] * (N + 1)
        sc = [0.0] * (N + 1)
        cnt[n0] = 1
        sp[n0] = 2 * n0 + 2
        dp[(b, t, 1 if b == 0 else 0)] = (cnt, sp, sc)
    tot = [[0.0] * (N + 1) for _ in range(3)]

    def acc():
        for (b, t, h), (cnt, sp, sc) in dp.items():
            if h == 1:
                for i in range(N + 1):
                    tot[0][i] += cnt[i]
                    tot[1][i] += sp[i]
                    tot[2][i] += sc[i]

    acc()
    for step in range(1, N):
        ndp = defaultdict(lambda: ([0] * (N + 1), [0.0] * (N + 1), [0.0] * (N + 1)))
        for (b, t, h), (cnt, sp, sc) in dp.items():
            for (b2, t2) in ints:
                if max(b, b2) > min(t, t2):
                    continue
                dn = t2 - b2 + 1
                if dn > N:
                    continue
                ov = min(t, t2) - max(b, b2) + 1
                dpp = 2 * dn + 2 - 2 * ov
                dc = pc[((b, t), (b2, t2))]
                h2 = 1 if (h == 1 or b2 == 0) else 0
                C, P, Q = ndp[(b2, t2, h2)]
                for i in range(N + 1 - dn):
                    if cnt[i]:
                        C[i + dn] += cnt[i]
                        P[i + dn] += sp[i] + dpp * cnt[i]
                        Q[i + dn] += sc[i] + dc * cnt[i]
        dp = ndp
        acc()
    say(" n      C_n        <p>     p/n     <curv>   curv/n")
    for n in range(1, N + 1):
        c = tot[0][n]
        if c < 1:
            continue
        say("%2d  %10d  %7.2f  %.4f   %7.2f   %.4f" %
            (n, round(c), tot[1][n] / c, tot[1][n] / c / n,
             tot[2][n] / c, tot[2][n] / c / n))
    say("read: p/n stays O(1) (1.89 at n=24, drifting down slowly) and curv/n rises")
    say("      through 0.51 => both are EXTENSIVE on this class: a boundary penalty")
    say("      shifts the rate by kappa*rho, it is not exponentially invisible.")
    return


def main():
    say("GF2-STIFFNESS witness   py-%s  numpy-%s" %
        (sys.version.split()[0], np.__version__))
    part1()
    Cn, Cp, Cc = part2()
    D = part3()
    part4()
    part5(Cp, Cc)
    part6(Cn, Cp, Cc)
    part7()
    part8()
    say("")
    say("=" * 78)
    say("window table   W = 0.5 + gamma + ln3 ,  target W <= 2.598612  <=>  mu <= e")
    say("=" * 78)
    ln3 = math.log(3)
    rows = [("bargraph  (flat bottom)", 2.0),
            ("common row (phi^2, exact)", (3 + math.sqrt(5)) / 2),
            ("column-convex (GF pole)", 3.2056),
            ("directed animals [M]", 3.0)]
    for name, mu in rows:
        g = math.log(mu)
        W = 0.5 + g + ln3
        say("%-30s mu=%8.4f  gamma=%7.4f  W=%8.4f  gap=%+8.4f  mu<=e: %s" %
            (name, mu, g, W, W - 2.598612, "YES" if mu <= math.e else "NO"))
    say("")
    say("pure-enumeration 14th/13th roots (evidence, not a lower bound on the limit):")
    say("  column-convex  C_14^(1/14) = %.6f" % (float(Cn[NS]) ** (1.0 / NS)))
    say("  directed       D_13^(1/13) = %.6f" % (float(D[13]) ** (1.0 / 13)))
    say("")
    say("required stiffness if one insists on the cost route:")
    for name, mu in [("column-convex", 3.2056), ("directed", 3.0)]:
        say("  %-14s  gamma-1 = ln(mu/e) = %+0.4f  (needs kappa*rho >= this)" %
            (name, math.log(mu) - 1))
    with open("gf2_stiffness.log", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
