# -*- coding: utf-8 -*-
"""
GG3-CLASS witness:  is the CPB defect-cluster class a polyomino class?

Pure exact integer enumeration. No MC, no randomness, no floating point in the
counting layer (floats only for logs).

Objects
-------
plaquette  p = (x0,x1,x2,x3, mu, nu),  0 <= mu < nu <= 3      (2-face of Z^4)
adjacency  = share exactly one link (a d=4 plaquette graph, degree Delta=20)

a_n  := # connected n-plaquette clusters containing the root plaquette
        p0 = ((0,0,0,0), 0, 1)            <- CPB / H3 class  (GV2:116-125)
m_n  := same, restricted to clusters all of whose plaquettes carry the SAME
        label (mu,nu)                     <- = d=2 face graph = fixed polyominoes
        (P3B.md:30 [T]; predicted m_n = n * A_n)
A_n  := # fixed polyominoes with n cells (translation classes)

Verdict target:  is  CPB class  ==  polyomino class?   strict inclusion?  disjoint?
"""
import sys, math
from collections import defaultdict

D = 4
ROOT = (tuple([0] * D), 0, 1)


def links(p):
    """the 4 links of plaquette p=(x,mu,nu).

    For direction d in {mu,nu}: the d-links sit at x and at x + e_other,
    other = the OTHER index of the pair (that is what closes the square).
    """
    x, mu, nu = p
    out = []
    for d in (mu, nu):
        other = nu if d == mu else mu
        out.append((x, d))                 # d-link at x
        y = list(x)
        y[other] += 1
        out.append((tuple(y), d))          # d-link at x + e_other
    return out


def neighbors(p):
    """all plaquettes q != p sharing exactly one link with p"""
    x, mu, nu = p
    nb = set()
    for (lx, d) in links(p):
        # 6 plaquettes contain link (lx,d): for each e != d, at position lx and lx-e
        for e in range(D):
            if e == d:
                continue
            a, b = (d, e) if d < e else (e, d)
            nb.add((lx, a, b))
            y = list(lx)
            y[e] -= 1
            nb.add((tuple(y), a, b))
    nb.discard(p)
    assert len(nb) == 20, (p, len(nb))
    return nb


NBR = {}


def nb(p):
    if p not in NBR:
        NBR[p] = tuple(sorted(neighbors(p)))
    return NBR[p]


def enumerate_clusters(nmax):
    """level-by-level exact enumeration of rooted connected clusters (sets)."""
    levels = []
    cur = {frozenset([ROOT])}
    levels.append(cur)
    for _ in range(nmax - 1):
        nxt = set()
        for C in cur:
            for p in C:
                for q in nb(p):
                    if q not in C:
                        nxt.add(C | {q})
        levels.append(nxt)
        cur = nxt
    return levels


def mono_label(C):
    labs = {(mu, nu) for (x, mu, nu) in C}
    return len(labs) == 1


def polyominoes(nmax):
    """fixed polyominoes (translation classes), brute force on Z^2, canonical by translation."""
    A = [0] * (nmax + 1)
    if nmax >= 1:
        A[1] = 1
    cur = {frozenset([(0, 0)])}
    for n in range(2, nmax + 1):
        nxt = set()
        for S in cur:
            for (a, b) in S:
                for (da, db) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    q = (a + da, b + db)
                    if q in S:
                        continue
                    T = S | {q}
                    mx = min(t[0] for t in T)
                    my = min(t[1] for t in T)
                    nxt.add(frozenset((t[0] - mx, t[1] - my) for t in T))
        A[n] = len(nxt)
        cur = nxt
    return A


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    say = print

    say('=' * 78)
    say('GG3-CLASS witness  |  is the CPB defect class a polyomino class?')
    say('=' * 78)

    # ---------------- S0: graph, degree, two orbits -------------------------
    say('\n[S0] 4D plaquette graph: vertices = plaquettes, edges = shared link')
    R = (tuple([0] * 4), 0, 1)
    nbs = neighbors(R)
    same = [q for q in nbs if (q[1], q[2]) == (0, 1)]
    diff = [q for q in nbs if (q[1], q[2]) != (0, 1)]
    say('  degree of root plaquette          : %d' % len(nbs))
    say('  same-label (in-plane) neighbours  : %d   %s' % (len(same), sorted(same)))
    say('  other-label neighbours            : %d' % len(diff))
    labs = defaultdict(int)
    for q in diff:
        labs[(q[1], q[2])] += 1
    say('  label histogram of the 16 oblique : %s' % dict(sorted(labs.items())))
    say('  => Delta = 4*(2d-3) = %d  (d=4)   [matches P3B.md:26]' % (4 * (2 * 4 - 3)))

    # ---------------- S1: full class enumeration ---------------------------
    say('\n[S1] full CPB class  a_n  (rooted connected plaquette clusters)')
    levels = enumerate_clusters(NMAX)
    a = [0] + [len(L) for L in levels]
    say('  n        : ' + ' '.join('%8d' % n for n in range(1, NMAX + 1)))
    say('  a_n      : ' + ' '.join('%8d' % a[n] for n in range(1, NMAX + 1)))
    ref = {1: 1, 2: 20, 3: 458, 4: 11132, 5: 280175, 6: 7218768}
    ok = all(a[n] == ref[n] for n in range(1, min(NMAX, 6) + 1))
    say('  ref(GC2-H23:54 / GV2:118 / P3B:42): %s  => %s'
        % ({n: ref[n] for n in range(1, min(NMAX, 6) + 1)}, 'MATCH' if ok else 'MISMATCH'))
    say('  ratios   : ' + ' '.join('   %6.4f' % (a[n] / a[n - 1]) for n in range(2, NMAX + 1)))
    say('  (1/n)lna : ' + ' '.join('%8.6f' % (math.log(a[n]) / n) for n in range(1, NMAX + 1)))

    # ---------------- S2: planar (mono-label) sub-class ---------------------
    say('\n[S2] polyomino sub-class  m_n  = clusters with ONE plaquette label')
    m = [0] * (NMAX + 1)
    for n in range(1, NMAX + 1):
        m[n] = sum(1 for C in levels[n - 1] if mono_label(C))
    A = polyominoes(NMAX)
    say('  n        : ' + ' '.join('%8d' % n for n in range(1, NMAX + 1)))
    say('  m_n      : ' + ' '.join('%8d' % m[n] for n in range(1, NMAX + 1)))
    say('  A_n      : ' + ' '.join('%8d' % A[n] for n in range(1, NMAX + 1)))
    say('  n*A_n    : ' + ' '.join('%8d' % (n * A[n]) for n in range(1, NMAX + 1)))
    say('  m_n == n*A_n ?  %s' % ('YES  => d=2 face graph IS Z^2, sub-class IS polyominoes'
                                 if all(m[n] == n * A[n] for n in range(1, NMAX + 1)) else 'NO'))

    # ---------------- S3: inclusion is STRICT and the gap explodes ----------
    say('\n[S3] a_n / m_n   (how much bigger the CPB class is than the polyomino sub-class)')
    say('  n        : ' + ' '.join('%12d' % n for n in range(1, NMAX + 1)))
    say('  a_n/m_n  : ' + ' '.join('%12.4f' % (a[n] / m[n]) for n in range(1, NMAX + 1)))
    say('  fraction : ' + ' '.join('%11.4f%%' % (100.0 * m[n] / a[n]) for n in range(1, NMAX + 1)))
    say('  n=2 witness of non-planarity: %d of %d two-clusters are NOT in one plane'
        % (a[2] - m[2], a[2]))
    say('  explicit: {(12)@0, (13)@0} shares the 1-link at 0 but spans e1,e2 and e1,e3')
    say('            => connected in the CPB class, lying in NO single 2-plane.')
    say('  => CPB class  STRICTLY CONTAINS the polyomino class (not equal, not disjoint).')

    # ---------------- S4: the monotone chain of growth constants -----------
    say('\n[S4] monotone chain  (all numbers are LOWER bounds for mu_CPB)')
    e = math.e
    chain = [('phi^2 = mu(D)      (directed col-convex)  GF3:24', 2.6180339887498948),
             ('1+sqrt2 = mu(D_1)  (directed col-convex)  GF3:135', 1 + math.sqrt(2)),
             ('mu(S_1)=2.831177   (two-sided |d|<=1)     GF3:140', 2.831177207208337),
             ('mu(S_4)>=3.190001  (two-sided |d|<=4)     GF3:143', 3.190001268135846),
             ('mu(row-convex)     [M]                    GE2:250', 3.2056),
             ('mu(planar)=4.0626  fixed polyominoes      [M]', 4.062570),
             ('a_6^(1/6)=13.9020  CPB class, hard        P3B:83', 13.9020),
             ('e^(3.458607)=31.77 CPB shell-mass median  GD2:259', math.exp(3.458607))]
    for name, v in chain:
        say('  %-52s %12.6f   %s' % (name, v, 'ABOVE e' if v > e else 'below e'))
    say('  e = %.6f  sits BETWEEN mu(D)=phi^2 (below) and every super-class (above).' % e)
    lo, hi = 2.6180339887498948, 2.831177207208337
    say('  => the "mu <= e" PASS holds ONLY on the SMALLEST class in the chain.')
    say('  => flipping the single "directed" bit moves mu by %.4f (%.2f%%), vs margin %.4f (%.2f%%)'
        % (hi - lo, 100 * (hi - lo) / lo, e - lo, 100 * (e - lo) / lo))

    # ---------------- S5: gaps to e ---------------------------------------
    say('\n[S5] the criterion  mu <= e  evaluated on each class')
    say('  %-34s %10s %10s %9s' % ('class', 'mu', 'mu/e', 'verdict'))
    for name, v in [('(a) bargraph        D_0', 2.0),
                    ('(c) D  directed col-convex', 2.6180339887498948),
                    ('(c\') S_1 two-sided |d|<=1', 2.831177207208337),
                    ('(b) S_4 two-sided |d|<=4', 3.190001268135846),
                    ('(b) planar polyomino (hard)', 2.927334013684),
                    ('(b) planar polyomino (lit.)', 4.062570),
                    ('(d) CPB class a_6 hard', 13.9020),
                    ('(d) CPB class extrap median', math.exp(3.249025))]:
        say('  %-34s %10.4f %10.4f %9s' % (name, v, v / e, 'YES' if v <= e else 'NO'))
    say('  => the passing rows are exactly the classes STRICTLY CONTAINED in the CPB class.')

    # ---------------- S6: what c0 the CPB class would really need ----------
    say('\n[S6] CPB-S window with the TRUE class entropy   (BianchiLemma:109-111)')
    say('  window:  c0 >= 0.5 + gamma_ent + ln3')
    for tag, g in [('gamma = phi^2-class   ln2.618=0.962424', math.log(2.6180339887498948)),
                   ('gamma = planar hard   ln2.927=1.074092', 1.074092),
                   ('gamma = planar lit.   ln4.0626=1.401823', 1.401823),
                   ('gamma = CPB hard      (1/6)ln a_6 = 2.632032', 2.632032),
                   ('gamma = CPB shell     (1/6)ln B_6 = 2.746029', 2.746029),
                   ('gamma = CPB median    (M3-M8)  = 3.249025', 3.249025)]:
        w = 0.5 + g + math.log(3)
        say('  %-46s c0 >= %8.6f   (vs 2.598612: %+.6f)' % (tag, w, w - 2.598612))
    say('  => phi^2 would report c0 >= 2.561036 <= 2.598612  "PASS" (margin 1.45%)')
    say('     the CPB class needs   c0 >= 4.230644          -> understate by 1.669608')

    say('\nexit 0')


if __name__ == '__main__':
    main()
