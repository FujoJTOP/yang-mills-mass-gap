# -*- coding: ascii -*-
# w7/m07_window2.py -- M07: recompute the window from Cbar == N; decide C_hat in {2,4}
# py -3.12 + mpmath (dps=40).  Sections:
#   S0 four constants; S1 provenance of 7.0633; S2 windows; S3 C_hat=4 (Weyl + cube + random closed
#   complexes); S4 consistency checks.  ASCII output only (log is ASCII).
from mpmath import mp, mpf, e as ME, log, nstr
import random, itertools

mp.dps = 40
N = 2
Delta = 20
eD = ME * Delta                      # lambda_a = e*Delta (P3-FULLREP-KP:89, Delta=20)
one48 = mpf(1) / 48

def line(s=""):
    print(s)

line("S0 -- four constants, dps=40")
vals = [("1/48", one48),
        ("1.766/48", mpf("1.766") / 48),
        ("3.532/48", mpf("3.532") / 48),
        ("7.0633/2/48", mpf("7.0633") / 2 / 48),
        ("7.06328527049/2/48 (exact thr)", mpf("7.06328527049") / 2 / 48),
        ("1.76582025642/48 (exact 8/(4eD)/48)", mpf("1.76582025642") / 48)]
for name, v in vals:
    line("  %-36s = %s" % (name, nstr(v, 20)))

line("")
line("S1 -- where 7.0633 comes from")
u48 = mpf(1) / 384                   # u = beta/(2N^2) at beta=1/48, N=2 -> 1/384
cf48 = one48 / (2 * N)               # c_f = beta/(2N) -> 1/192
T = 1 / (eD * u48)
line("  e*Delta                       = %s" % nstr(eD, 20))
line("  u(1/48) = beta/(2N^2)         = %s" % nstr(u48, 20))
line("  cf(1/48) = beta/(2N)          = %s" % nstr(cf48, 20))
line("  1/(e*Delta*u(1/48))           = %s" % nstr(T, 20))
line("  2*N^2*48/(e*Delta)            = %s" % nstr(2 * N * N * 48 / eD, 20))
line("  ratio of the two (must be 1)  = %s" % nstr(T / (2 * N * N * 48 / eD), 20))
line("  cf/u at beta=1/48             = %s   (= N = %d)" % (nstr(cf48 / u48, 20), N))
line("  1/(e*Delta*cf(1/48))          = %s   (= thr / N)" % nstr(1 / (eD * cf48), 20))
line("  7.0633/2 (4-digit)            = %s" % nstr(mpf("7.0633") / 2, 20))
line("  exact T/2                     = %s" % nstr(T / 2, 20))
line("  3.532 vs T/2  rel.diff        = %s" % nstr((mpf("3.532") - T / 2) / mpf("3.532"), 8))

line("")
line("S2 -- windows  beta < 8/(C_hat*Cbar*e*Delta)   (u = beta/8 at N=2)")
line("  C_hat  Cbar      beta_max            x(1/48)        1/48 inside?")
for ch in (1, 2, 4):
    for cb in (1, 2, 8, 20):
        bm = 8 / (mpf(ch) * mpf(cb) * eD)
        line("  %-6d %-9d %-19s %-14s %s" % (ch, cb, nstr(bm, 12), nstr(48 * bm, 12),
                                             "YES" if bm > one48 else "NO"))

line("")
line("S3 -- C_hat: per-gamma sup of (prod_l R(l))^(1/F)")
def spin_mult(kappa):
    """multiplicities of total spin (in units of 2j) in (1/2)^{otimes kappa}."""
    d = {1: 1}
    for _ in range(kappa - 1):
        nd = {}
        for two_j, m in d.items():
            for s in (two_j - 1, two_j + 1):
                if s >= 0:
                    nd[s] = nd.get(s, 0) + m
        d = nd
    return d

line("  kappa   mult(1 in f^k)   mult(adj in f^k)   R=1+m_adj/m_1")
Rtab = {}
for k in range(1, 9):
    d = spin_mult(k)
    m1, madj = d.get(0, 0), d.get(2, 0)
    line("  %-7d %-16d %-18d %s" % (k, m1, madj, ("-" if m1 == 0 else nstr(1 + mpf(madj) / m1, 10))))
    if m1:
        Rtab[k] = 1 + mpf(madj) / m1
line("  check mult(1 in f^(2m)) = Catalan(m): %s" % [spin_mult(2 * m).get(0, 0) for m in range(1, 6)])

line("  4D cube: F=6, E2=12 -> prod R = 2^12 = %d, (prod R)^(1/6) = %s"
     % (2 ** 12, nstr(mpf(2) ** (mpf(12) / 6), 12)))
line("  bound  prod R <= 2^{E2} <= 2^{2F} = 4^F  (E2<=2F for closed complexes): 4^6 = %d" % 4 ** 6)

random.seed(20261125)
L = 4
dirs = list(range(4))
def coords(x):
    return [(x // L ** i) % L for i in range(4)]
def idx(c):
    return sum(c[i] * L ** i for i in range(4))
def sadd(x, d):
    c = coords(x); c[d] = (c[d] + 1) % L
    return idx(c)

def cube_faces(x, S):
    """dict: frozenset(4 corners) -> tuple of its 4 link ids; 6 faces of the 3-cube.

    Faces are keyed by their corner SET (not the minimal corner) so that the two
    plaquettes in a 4-periodic direction that share a minimal corner stay distinct.
    """
    v = {(): x}
    for d in S:
        for k in list(v.keys()):
            v[k + (d,)] = sadd(v[k], d)
    def lid(u, w):
        return (min(u, w), max(u, w))
    out = {}
    for pair in itertools.combinations(S, 2):
        a, b = pair
        for other in S:
            if other in pair:
                continue
            for base in (v[()], v[(other,)]):
                pa, pb, pab = sadd(base, a), sadd(base, b), sadd(sadd(base, a), b)
                corners = (base, pa, pb, pab)
                out[frozenset(corners)] = (lid(base, pa), lid(base, pb), lid(pa, pab), lid(pb, pab))
    return out

def evaluate(faces):
    """faces: list of 4-link tuples; returns (F, khist, handshake_residual, gm)."""
    kmult = {}
    for edges in faces:
        for e in edges:
            kmult[e] = kmult.get(e, 0) + 1
    khist = {}
    for kap in kmult.values():
        khist[kap] = khist.get(kap, 0) + 1
    F = len(faces)
    hs = 4 * F - sum(k * n for k, n in khist.items())
    known = all(k in Rtab for k in khist)
    gm = mpf(2) ** (sum(n * log(Rtab[k], 2) for k, n in khist.items()) / F) if known else None
    return F, dict(sorted(khist.items())), hs, gm

cube_sets = [cube_faces(x, S) for x in range(L ** 4) for S in itertools.combinations(dirs, 3)]
best = mpf(0); bestinfo = None; seen = 0; cnt4 = 0; hsbad = 0; kmax = 0; odd = 0
for trial in range(4000):
    pick = random.sample(range(len(cube_sets)), random.randint(1, 6))
    occ = {}
    for i in pick:
        for f, edges in cube_sets[i].items():
            if f in occ:
                del occ[f]
            else:
                occ[f] = edges
    faces = list(occ.values())
    if not faces:
        continue
    seen += 1
    F, khist, hs, gm = evaluate(faces)
    if hs != 0:
        hsbad += 1
    if any(k % 2 for k in khist):
        odd += 1
    kmax = max(kmax, max(khist))
    if gm is not None and gm > best:
        best, bestinfo = gm, (F, khist)
    if gm is not None and abs(gm - 4) < mpf("1e-25"):
        cnt4 += 1
line("  random closed complexes: %d nonempty, max (prod R)^(1/F) = %s  [F=%d, hist=%s]"
     % (seen, nstr(best, 14), bestinfo[0], bestinfo[1]))
line("  samples attaining exactly 4: %d ; max kappa seen: %d ; odd-kappa samples: %d ; handshake violations: %d"
     % (cnt4, kmax, odd, hsbad))
line("  all gammas in sample have gm <= 4 : %s" % (best <= 4))

dmax = mpf(0); dhit = 0; dbest = None
for fs in cube_sets:
    F, khist, hs, gm = evaluate(list(fs.values()))
    if gm > dmax:
        dmax, dbest = gm, (F, khist)
    if abs(gm - 4) < mpf("1e-25"):
        dhit += 1
line("  all %d single 3-cubes: max gm = %s [F=%d, hist=%s], number equal to 4: %d"
     % (len(cube_sets), nstr(dmax, 14), dbest[0], dbest[1], dhit))

line("")
line("S4 -- consistency")
line("  e*Delta*beta at beta=1/48      = %s   (>1 means KP fails)" % nstr(eD * one48, 12))
line("  e*Delta*beta at beta=1/12      = %s" % nstr(eD / 12, 12))
line("  beta at e*Delta*beta=1         = %s" % nstr(1 / eD, 16))
line("  (1/e*Delta)/(1/48)             = %s" % nstr(48 / eD, 12))
line("  C_hat*Cbar (=8) vs threshold   = %s   (= e*Delta/48 ?)" % nstr(8 / T, 12))
line("  e*Delta/48                     = %s" % nstr(eD / 48, 12))
line("  8*u*eD / (eD*beta)             = %s   (8u = beta)" % nstr(8 * u48 / one48, 12))
line("  Cbar bound at C_hat=4          : Cbar < %s" % nstr(T / 4, 12))
line("  Cbar bound at C_hat=2          : Cbar < %s" % nstr(T / 2, 12))
line("  window ratio c06(Cbar=20) vs exact(Cbar=2): %s" % nstr(20 / 2, 6))
