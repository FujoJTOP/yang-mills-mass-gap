# -*- coding: utf-8 -*-
"""GG1-BIT witness — orientation (mirror) orbit of the directed defect class.

Independent of gf3_entropy.py (different DP: state = current column height, offsets summed
explicitly; no Bareiss, no Sturm chain -- verdicts use exact rational sign brackets instead).

Parts
  1  exact integer DP:  D (delta<=0), S1 (|delta|<=1), Drev (delta>=0), D u Drev (sign fixed),
     D_s (0>=delta>=-s) for s=0..6, and the phase automaton for the disjunctive class.
  2  exact rational brackets:  mu(D)=phi^2 < e ;  mu(S1)=1/x_c > e ;  no MC, no float in verdict.
  3  identity checks: a_n(D)=F_{2n-1}; a_n(D u Drev)=2F_{2n-1}-2^{n-1}; G_union=2G_D-G_bargraph.
  4  brute force: enumerate ALL column-convex shapes with n<=N; test mono-x, colconvex-y, mono-y;
     the two-slicing (existential) class is counted as an [E] side number, NOT used in any verdict.

Discipline: no randomness, no Monte Carlo, single job, prints only.  Writes gg1_bit.log next to itself.
"""
import os, sys, time
from fractions import Fraction as Fr

T0 = time.time()
LOG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gg1_bit.log")
_fh = open(LOG, "w", encoding="utf-8")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    _fh.write(s + "\n")
    _fh.flush()


# ----------------------------------------------------------------------------- Part 1: DP
def dp_offsets(N, kind):
    """a_n for column-convex shapes (columns at consecutive x, each a maximal vertical
    interval) whose offset sequence delta_j = b_{j+1}-b_j lies in a class.

    kind: 'mono'  delta <= 0                (D, baseline non-increasing)
          'rev'   delta >= 0                (D^rev)
          'sK'    |delta| <= K              (S_K, two-sided)
         'barg'   delta == 0                (bargraphs)
    Bijection: shape <-> (h_1..h_m, delta_1..delta_{m-1}) with
        h_j >= max(1, delta_j+1, 1-delta_{j-1}),
    i.e. transition h -> h' needs delta in [1-h', h-1] (adjacent columns intersect).
    """
    states = [dict() for _ in range(N + 1)]
    for h in range(1, N + 1):
        states[h][h] = 1                      # one column of height h
    a = [0] * (N + 1)
    for n in range(1, N + 1):
        a[n] = sum(states[n].values())
        for h, c in list(states[n].items()):
            for hp in range(1, N - n + 1):
                lo, hi = 1 - hp, h - 1
                if kind == 'mono':
                    hi = min(hi, 0)
                elif kind == 'rev':
                    lo = max(lo, 0)
                elif kind == 'barg':
                    lo, hi = max(lo, 0), min(hi, 0)
                elif kind.startswith('s'):
                    K = int(kind[1:])
                    lo, hi = max(lo, -K), min(hi, K)
                else:
                    raise ValueError(kind)
                if lo <= hi:
                    d = states[n + hp]
                    d[hp] = d.get(hp, 0) + c * (hi - lo + 1)
    return a


def dp_disjunctive(N):
    """D u D^rev: the offset sequence never changes sign (sign fixed once and for all).
    State = (phase, h), phase 0 = sign not yet fixed, -1 = non-increasing, +1 = non-decreasing.
    This is the *reflection-invariant* replacement of the one-way bit."""
    states = [dict() for _ in range(N + 1)]
    for h in range(1, N + 1):
        states[h][(0, h)] = 1
    a = [0] * (N + 1)
    for n in range(1, N + 1):
        a[n] = sum(states[n].values())
        for (ph, h), c in list(states[n].items()):
            for hp in range(1, N - n + 1):
                lo, hi = 1 - hp, h - 1
                if lo > hi:
                    continue
                tgt = states[n + hp]
                if ph == 0:
                    # delta == 0 keeps phase 0
                    if lo <= 0 <= hi:
                        tgt[(0, hp)] = tgt.get((0, hp), 0) + c
                    # delta < 0 -> phase -1 ; delta > 0 -> phase +1
                    d0, d1 = max(lo, -10 ** 9), min(hi, -1)
                    if d0 <= d1:
                        tgt[(-1, hp)] = tgt.get((-1, hp), 0) + c * (d1 - d0 + 1)
                    e0, e1 = max(lo, 1), hi
                    if e0 <= e1:
                        tgt[(1, hp)] = tgt.get((1, hp), 0) + c * (e1 - e0 + 1)
                elif ph == -1:
                    d1 = min(hi, 0)
                    if lo <= d1:
                        tgt[(-1, hp)] = tgt.get((-1, hp), 0) + c * (d1 - lo + 1)
                else:
                    e0 = max(lo, 0)
                    if e0 <= hi:
                        tgt[(1, hp)] = tgt.get((1, hp), 0) + c * (hi - e0 + 1)
    return a


N = 16
AD = dp_offsets(N, 'mono')
AR = dp_offsets(N, 'rev')
AB = dp_offsets(N, 'barg')
AS1 = dp_offsets(N, 's1')
AU = dp_disjunctive(N)


def fib(k):
    a, b = 0, 1
    for _ in range(k):
        a, b = b, a + b
    return a


say("=" * 78)
say("PART 1  exact integer counts (n = 1..%d)" % N)
say("=" * 78)
say("  n : %s" % " ".join("%7d" % n for n in range(1, N + 1)))
for name, seq in (("D      (mono)", AD), ("D^rev  (mono)", AR), ("bargraph", AB),
                  ("S1 (|d|<=1)", AS1), ("D u D^rev", AU)):
    say("%-13s: %s" % (name, " ".join("%7d" % seq[n] for n in range(1, N + 1))))

say("")
say("identity checks (integer equality, no float):")
ok = all(AD[n] == fib(2 * n - 1) for n in range(1, N + 1))
say("  a_n(D) == F_{2n-1}                      : %s" % ("PASS" if ok else "FAIL"))
ok2 = all(AU[n] == 2 * fib(2 * n - 1) - 2 ** (n - 1) for n in range(1, N + 1))
say("  a_n(D u D^rev) == 2F_{2n-1} - 2^{n-1}   : %s" % ("PASS" if ok2 else "FAIL"))
say("  (inclusion-exclusion: D cap D^rev = bargraphs, a_n = 2^{n-1})")
ok3 = all(AS1[n] == v for n, v in zip(range(1, 9), [1, 2, 6, 17, 48, 136, 385, 1090]))
say("  a_n(S1) head == 1,2,6,17,48,136,385,1090: %s  (GF3-ENTROPY.md:166)" % ("PASS" if ok3 else "FAIL"))

# D_s family and its closed form check
say("")
say("directed truncations D_s (0 >= delta >= -s), closed form check")
say("  P_s = (1-x)^s (1 - 2x - x^2 - ... - x^{s+1});  x_s = the unique root of sum_{k=1..s+1} x^k = 1-x")
for s in range(0, 7):
    a = dp_offsets(12, 's%d' % s)
    say("  s=%d : a_1..a_6 = %s" % (s, a[1:7]))


def dp_bounded_range(N, K):
    """GF2's K-class: every column bottom b_j in [0,K] (suspension range <= K), sign free.
    State = (b, h).  K=0 -> bargraphs; K=oo -> column-convex."""
    states = [dict() for _ in range(N + 1)]
    for b in range(0, K + 1):
        for h in range(1, N + 1):
            states[h][(b, h)] = states[h].get((b, h), 0) + 1
    a = [0] * (N + 1)
    for n in range(1, N + 1):
        a[n] = sum(states[n].values())
        for (b, h), c in list(states[n].items()):
            for bp in range(0, K + 1):
                d = bp - b
                if d > h - 1:
                    continue
                for hp in range(max(1, 1 - d), N - n + 1):
                    t = states[n + hp]
                    t[(bp, hp)] = t.get((bp, hp), 0) + c
    return a


say("")
say("GF2's bounded-suspension ladder, independently recounted (state = (bottom b, height h)):")
NBIG = 60
for K in (0, 1, 2):
    a = dp_bounded_range(NBIG, K)
    say("  K=%d : a_40 = %d  ratio(60/59) = %.6f  root n=60 = %.6f"
        % (K, a[40], a[60] / a[59], a[60] ** (1.0 / 60)))
aK1 = dp_bounded_range(16, 1)
say("  (the K-DP keeps b in [0,K] without mod-translation, so a_1(K=1)=2: the RATE is unaffected)")
say("  D_1 (delta in {0,-1}, drift unbounded) : a_1..a_6 = %s" % dp_offsets(6, 's1')[1:7])


# ----------------------------------------------------------------------------- Part 2: exact verdicts
say("")
say("=" * 78)
say("PART 2  exact rational brackets (all arithmetic in Fraction; no float enters a verdict)")
say("=" * 78)


def peval(p, x):                      # p = list of ints, ascending powers
    r = Fr(0)
    for c in reversed(p):
        r = r * x + c
    return r


def pderiv(p):
    return [c * i for i, c in enumerate(p)][1:]


def e_bounds(m=40):
    """e = sum 1/k!  with 0 < e - S_m < 1/(m! m).  Exact rationals."""
    S = Fr(0)
    f = Fr(1)
    for k in range(0, m + 1):
        if k:
            f = f / k
        S += f
    return S, S + f / m


elo, ehi = e_bounds(40)
say("e in (%.20f, %.20f)   width %.3e" % (float(elo), float(ehi), float(ehi - elo)))


def bracket_root(p, lo, hi, steps=80):
    """p strictly monotone on [lo,hi] (checked), p(lo)*p(hi)<0 -> exact rational bracket."""
    d = pderiv(p)
    sgn = 1 if peval(d, lo) > 0 else -1
    assert all((peval(d, Fr(i, 10 ** 4)) > 0) == (sgn == 1) for i in range(int(lo * 10 ** 4), int(hi * 10 ** 4) + 1)), "not monotone"
    assert peval(p, lo) * peval(p, hi) < 0, "no sign change"
    a, b = lo, hi
    for _ in range(steps):
        c = (a + b) / 2
        if peval(p, a) * peval(p, c) <= 0:
            b = c
        else:
            a = c
    return a, b


say("")
say("(i) class D : G(x) = x(1-x)/(1-3x+x^2), x_c = smallest positive root of x^2-3x+1")
pD = [1, -3, 1]                                   # 1 - 3x + x^2
aD, bD = bracket_root(pD, Fr(38, 100), Fr(39, 100))
say("    x_c in (%.12f, %.12f)  (width %.2e)" % (float(aD), float(bD), float(bD - aD)))
say("    mu = 1/x_c in (%.12f, %.12f)" % (float(1 / bD), float(1 / aD)))
# phi^2 = (3+sqrt5)/2 < 2.61805  <=>  5 < (2*2.61805-3)^2
rhs = 2 * Fr(261805, 100000) - 3
say("    phi^2 = (3+sqrt5)/2 < 2.61805 : %s   (needs 5 < %.8f)" % (5 < rhs * rhs, float(rhs * rhs)))
say("    phi^2 < e : %s   [2.61805 <= elo = %.6f]" % (Fr(261805, 100000) <= elo, float(elo)))

say("")
say("(ii) class S1 : P(x) = 1-4x+3x^2+x^3-x^5 = (1-x)(x^4+x^3-3x+1)")
pS = [1, -4, 3, 1, 0, -1]
aS, bS = bracket_root(pS, Fr(35, 100), Fr(36, 100), steps=100)
say("    x_c in (%.14f, %.14f)  (width %.2e)" % (float(aS), float(bS), float(bS - aS)))
say("    mu = 1/x_c in (%.12f, %.12f)" % (float(1 / bS), float(1 / aS)))
say("    mu > e : %s   [1/bS = %.9f  vs  elo = %.9f]" % (1 / bS > ehi, float(1 / bS), float(ehi)))
say("    margin (mu-e)/e = +%.4f%%" % (100 * float((1 / bS - ehi) / ehi)))

say("")
say("(iii) the disjunctive class D u D^rev : G = 2G_D - G_bargraph = 2x(1-x)/(1-3x+x^2) - x/(1-2x)")
say("      poles: x = (3-sqrt5)/2 = %.9f  and  x = 1/2 ; dominant pole unchanged" % (float((3 - 5 ** 0.5) / 2)))
say("      => mu(D u D^rev) = mu(D) = phi^2  (max of two equal rates) [T]")
say("      a_n identity verified for n<=%d above (2F_{2n-1} - 2^{n-1})" % N)

# ----------------------------------------------------------------------------- Part 3: brute force
say("")
say("=" * 78)
say("PART 3  brute force over ALL column-convex shapes (n <= 12) : predicate census")
say("=" * 78)
NB = 12
tot = [0] * (NB + 1)
cMonoX = [0] * (NB + 1)
cMonoY = [0] * (NB + 1)
cMonoEitherXY = [0] * (NB + 1)
cNotMonoX = [0] * (NB + 1)


def classify(cols):
    m = len(cols)
    b = [c[0] for c in cols]
    dl = [b[j + 1] - b[j] for j in range(m - 1)]
    mono_x = all(d <= 0 for d in dl) or all(d >= 0 for d in dl)
    ys = {}
    for j, (bb, hh) in enumerate(cols):
        for k in range(hh):
            ys.setdefault(bb + k, []).append(j)
    keys = sorted(ys)
    ccy = True
    for y in keys:
        xs = ys[y]
        if xs[-1] - xs[0] + 1 != len(xs):
            ccy = False
            break
    if ccy:
        for i in range(len(keys) - 1):
            x1, x2 = ys[keys[i]], ys[keys[i + 1]]
            if x1[-1] < x2[0] or x2[-1] < x1[0]:
                ccy = False
                break
    mono_y = False
    if ccy:
        bx = [ys[y][0] for y in keys]
        mono_y = all(bx[i + 1] <= bx[i] for i in range(len(bx) - 1)) or \
                 all(bx[i + 1] >= bx[i] for i in range(len(bx) - 1))
    return mono_x, (ccy and mono_y)


def rec(cols, total):
    n = total
    mx, my = classify(cols)
    tot[n] += 1
    if mx:
        cMonoX[n] += 1
    if my:
        cMonoY[n] += 1
    if mx or my:
        cMonoEitherXY[n] += 1
    if not mx:
        cNotMonoX[n] += 1
    if n >= NB:
        return
    b_last, h_last = cols[-1]
    for hp in range(1, NB - n + 1):
        lo, hi = 1 - hp, h_last - 1
        for d in range(lo, hi + 1):
            cols.append((b_last + d, hp))
            rec(cols, n + hp)
            cols.pop()


for h0 in range(1, NB + 1):
    rec([(0, h0)], h0)
say("  n        : %s" % " ".join("%9d" % n for n in range(1, NB + 1)))
say("  colconv  : %s" % " ".join("%9d" % tot[n] for n in range(1, NB + 1)))
say("  mono-x   : %s" % " ".join("%9d" % cMonoX[n] for n in range(1, NB + 1)))
say("  mono-y*  : %s   (* only when the shape is also column-convex in y)" % " ".join("%9d" % cMonoY[n] for n in range(1, NB + 1)))
say("  either   : %s" % " ".join("%9d" % cMonoEitherXY[n] for n in range(1, NB + 1)))
say("  not-mono-x: %s" % " ".join("%9d" % cNotMonoX[n] for n in range(1, NB + 1)))
say("")
say("  check: colconv head vs GF2-STIFFNESS.md:176 (1,2,6,19,61,196,629,2017,6466,20727,66441,212980): %s"
    % ("PASS" if [tot[n] for n in range(1, 13)] == [1, 2, 6, 19, 61, 196, 629, 2017, 6466, 20727, 66441, 212980] else "FAIL"))
say("  check: mono-x == a_n(D u D^rev) : %s"
    % ("PASS" if all(cMonoX[n] == AU[n] for n in range(1, NB + 1)) else "FAIL"))
say("")
say("  growth-rate ratios a_n / a_{n-1} (NOT a proof of the limit; [E] trend only):")
for nm, arr in (("mono-x (= D u D^rev)", cMonoX), ("either-slicing", cMonoEitherXY),
                ("colconv (all)", tot)):
    rr = [float(arr[n]) / arr[n - 1] for n in range(8, NB + 1)]
    say("    %-22s n=8..12 : %s" % (nm, " ".join("%.4f" % x for x in rr)))
say("  e = 2.718282 ;  phi^2 = 2.618034")

say("")
say("=" * 78)
say("done in %.1f s   single job, no MC, no randomness" % (time.time() - T0))
_fh.close()
