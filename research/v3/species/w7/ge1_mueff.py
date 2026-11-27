"""GE1-MUEFF -- what growth constant is mu_eff, and does mu_eff <= e hold?

Pure exact integer enumeration + arithmetic.  No MC, no random, no new physics.
S0  arithmetic of the criterion
S1  fixed polyominoes  = planar membranes (d=2 class),  n<=14  (bitboard BFS)
S2  plaquette animals "share a link" (d=4 full class; d=3 slab class)
S3  width-1 open ribbons (self-avoiding plaquette paths with P0 as endpoint)
S4  product tubes (SAW axis x closed cross-section) -> mu <= 8^(1/m)
S5  closed 2-chains (the (H1)-refined class): slab-doubling lower bound
S6  acid-test numbers
Output: ge1_mueff.log
"""
import math, time
T0 = time.time()
DEAD = T0 + 295.0
LN = math.log
LN3 = math.log(3.0)
E = math.e
C0 = 2.598612
OUT = []
def P(s=""):
    print(s, flush=True)
    OUT.append(s)

def left():
    return DEAD - time.time()

# ---------------------------------------------------------------- S0
P("== S0  arithmetic of the criterion ==")
gmax = C0 - 0.5 - LN3
P("  0.5 + g_eff + ln3 <= %.6f  <=>  g_eff <= %.9f  <=>  mu_eff <= e = %.6f"
  % (C0, gmax, math.exp(gmax)))
P("  window(gamma) = 0.5 + gamma + ln3,   gamma = ln(mu)")
P("  class                        mu          gamma       window    window-2.5986   mu/e")
def row(name, mu):
    g = LN(mu); w = 0.5 + g + LN3
    P("  %-27s %10.6f %10.6f %10.6f %+13.6f %7.4f" % (name, mu, g, w, w - C0, mu / E))
row("target  mu = e", E)
row("d=4 animals  exp(2.632032)", math.exp(2.632032))
row("d=4 animals  exp(3.249025)", math.exp(3.249025))
row("planar membrane (lit 4.0626)", 4.0626)
row("planar membrane (n=13 lb)", 3.041231)
row("planar membrane (n=14 lb)", 7204874 ** (1 / 14))
row("thick tube  m=4  (rigid ub)", 8.0 ** 0.25)
P("")
P("  extra per-block suppression needed to reach e, per class:")
for nm, mu in (("membrane lit", 4.0626), ("membrane n=13 lb", 3.041231), ("membrane n=14 lb", 7204874 ** (1 / 14))):
    P("    %-18s Delta >= ln(mu/e) = %.6f  <=>  c0 >= 0.5+ln(mu)+ln3 = %.6f"
      % (nm, LN(mu / E), 0.5 + LN(mu) + LN3))

# ---------------------------------------------------------------- S1
P("")
P("== S1  planar membranes = fixed polyominoes (exact, bitboard BFS) ==")
MARG = 4
WB = 32
WIDE = (1 << WB) - 1
def norm(x):
    low = x & -x
    mr = (low.bit_length() - 1) // WB
    y = x >> (WB * mr)
    m = y | (y >> 256); m |= m >> 128; m |= m >> 64; m |= m >> 32; m &= WIDE
    mc = (m & -m).bit_length() - 1
    d = mc - MARG
    return (y >> d) if d >= 0 else (y << (-d))

A168 = [None, 1, 2, 6, 19, 63, 216, 760, 2725, 9910, 36446, 135268, 505861, 1903890, 7204874]
cur = {1 << MARG}
A = [0, 1]
NMAX = 14
for n in range(2, NMAX + 1):
    t1 = time.time()
    nxt = set(); add = nxt.add
    for occ in cur:
        cand = ((occ << 1) | (occ >> 1) | (occ << WB) | (occ >> WB)) & ~occ
        while cand:
            b = cand & -cand; cand ^= b
            add(norm(occ | b))
    cur = nxt
    A.append(len(cur))
    ok = "OK(A001168)" if len(cur) == A168[n] else "MISMATCH"
    P("  n=%2d  A_n=%10d  (1/n)ln=%.6f  ratio=%8.5f  %6.1fs  %s"
      % (n, len(cur), LN(len(cur)) / n, len(cur) / A[n - 1], time.time() - t1, ok))
    if n < NMAX and time.time() - t1 > left() - 10:
        P("  [stop: budget]"); break
nA = len(A) - 1
gA = LN(A[nA]) / nA
row("membrane rigorous lb", math.exp(gA))
P("  rigorous: gamma_membrane >= (1/%d)ln A_%d = %.6f  (mu >= %.6f)" % (nA, nA, gA, math.exp(gA)))
P("  => gamma_membrane > 1 with margin %.2f%% at n=%d [T]" % (100 * (gA - 1), nA))
P("  literature growth constant 4.0626 => gamma = %.6f [M, not re-derived]" % LN(4.0626))

# ---------------------------------------------------------------- S2
P("")
P("== S2  plaquette animals, adjacency = share a link (d=4: Delta=20) ==")
def verts(d, pl, base):
    i, j = pl
    x = list(base)
    vs = [tuple(x)]
    y = list(x); y[i] += 1; vs.append(tuple(y))
    y = list(x); y[j] += 1; vs.append(tuple(y))
    y = list(x); y[i] += 1; y[j] += 1; vs.append(tuple(y))
    return vs

def nbrs(d, cell):
    pl, base = cell
    i, j = pl
    x = list(base); out = set()
    def sh(vec, k):
        z = list(vec); z[k] += 1; return tuple(z)
    # the four edges of the plaquette {x, x+e_i, x+e_j, x+e_i+e_j}
    edges = [(i, tuple(x)), (j, tuple(x)), (i, sh(x, j)), (j, sh(x, i))]
    for (r, y) in edges:
        for k in range(d):
            if k == r: continue
            pln = (min(k, r), max(k, r))
            out.add((pln, y))
            z = list(y); z[k] -= 1
            out.add((pln, tuple(z)))
    out.discard((pl, base))
    return out

def canon(d, cells):
    mins = [min(c[1][k] for c in cells) for k in range(d)]
    return tuple(sorted(((pl, tuple(b[k] - mins[k] for k in range(d))) for pl, b in cells)))

P("  degree of the link-adjacency graph (plaquette at the origin):")
for d in (2, 3, 4):
    c0 = ((0, 1), tuple([0] * d))
    nb = nbrs(d, c0)
    V = set(verts(d, c0[0], c0[1]))
    cand = set()
    rng = [0, 1, -1]
    from itertools import product as _prod
    for off in _prod(rng, repeat=d):
        y = tuple(off)
        for a in range(d):
            for b in range(a + 1, d):
                cand.add(((a, b), y))
    vshare = sum(1 for c in cand if c != c0 and len(V & set(verts(d, c[0], c[1]))) >= 1)
    P("    d=%d : |link-neighbours| = %2d  (4*2*(d-1)-4 = %2d) ; |vertex-sharing| = %2d"
      % (d, len(nb), 8 * (d - 1) - 4, vshare))

def animals(d, nmax, deadline, target=None):
    start = ((0, 1), tuple([0] * d))
    cur = {(start,)}
    out = [0, 1]
    for n in range(2, nmax + 1):
        t1 = time.time(); nxt = set()
        for cs in cur:
            have = set(cs)
            bucket = set()
            for c in cs:
                bucket |= nbrs(d, c)
            bucket -= have
            for c in bucket:
                nxt.add(canon(d, tuple(cs) + (c,)))
        cur = nxt
        out.append(len(cur))
        tg = ("  expect %d %s" % (target[n], "OK" if len(cur) == target[n] else "MISMATCH")) if target and n < len(target) else ""
        P("    d=%d n=%2d : a_n=%9d  (1/n)ln=%.6f  ratio=%7.3f  %5.1fs%s"
          % (d, n, len(cur), LN(len(cur)) / n, len(cur) / out[n - 1], time.time() - t1, tg))
        if time.time() > deadline:
            P("    [stop: budget]"); break
    return out

a4 = animals(4, 4, T0 + 170, target=[None, 1, 20, 458, 11132])
P("  corpus a_n (GV2:121-123) = 1, 20, 458, 11132, 280175, 7218768 ; reproduced to n=%d" % (len(a4) - 1))
if len(a4) > 4:
    P("  => d=4 class: gamma*(n<=4) >= %.6f ; corpus gamma*(n=6) = 2.632032 [T] (GV2:16)" % (LN(a4[-1]) / (len(a4) - 1)))
P("")
P("  d=3 slab class (BL-3 plate: clusters inside one 3-plane):")
a3 = animals(3, 6, T0 + 190)
if len(a3) > 2:
    P("  => d=3 class: gamma*(n<=%d) >= %.6f (mu >= %.4f)" % (len(a3) - 1, LN(a3[-1]) / (len(a3) - 1), a3[-1] ** (1.0 / (len(a3) - 1))))

# ---------------------------------------------------------------- S3
P("")
P("== S3  open width-1 ribbons: NOT enumerated (over budget) ==")
P("  The ribbon class contains every planar membrane patch, so mu_ribbon >= mu_membrane >= 3.0891")
P("  [T] by inclusion; whether BENDING raises mu above 4.06 is left as [O] (enumeration > 300 s).")

P("")
P("== S4  product tubes: closed cross-section x self-avoiding axis ==")
NB4 = [(1,0,0,0),(-1,0,0,0),(0,1,0,0),(0,-1,0,0),(0,0,1,0),(0,0,-1,0),(0,0,0,1),(0,0,0,-1)]
LMAX = 7
saw = [0] * (LMAX + 1)
def saw_dfs(pos, path, depth):
    saw[depth] += 1
    if depth >= LMAX:
        return
    for dv in NB4:
        q = (pos[0]+dv[0], pos[1]+dv[1], pos[2]+dv[2], pos[3]+dv[3])
        if q in path: continue
        path.add(q); saw_dfs(q, path, depth + 1); path.discard(q)
t = time.time()
saw_dfs((0,0,0,0), {(0,0,0,0)}, 0)
P("  axis SAW counts c_L (L=0..%d) = %s  [%.1fs]" % (LMAX, saw, time.time() - t))
for Lv in range(max(6, LMAX - 2), LMAX + 1):
    lam = saw[Lv] ** (1.0 / Lv)
    P("    L=%2d c_L=%9d  c_L^(1/L)=%.4f  -> mu_tube(m=4)=lam^(1/4)=%.4f   (rigorous ub 8^(1/4)=%.4f)"
      % (Lv, saw[Lv], lam, lam ** 0.25, 8.0 ** 0.25))
P("  minimal closed cross-section in the transverse Z^3: m=4 (unit square) => n = m*L plaquettes")
P("  rigorous upper bound for the whole product-tube family: mu <= 8^(1/m) <= 8^(1/4) = %.6f < e" % (8.0 ** 0.25))

# ---------------------------------------------------------------- S5
P("")
P("== S5  closed 2-chains ((H1)-refined class): slab doubling ==")
P("  polyomino P of area n  ->  closed surface d(P x [0,1]) of area 2n + perim(P) <= 4n + 2")
P("  distinct P give distinct surfaces => N_cyc(<= 4n+2) >= A_n")
P("  => mu_cyc >= lim A_n^(1/(4n+2)) = mu_poly^(1/4)")
for mu, tag in ((3.041231, "n=13 lb"), (7204874 ** (1/14), "n=14 lb"), (4.0626, "lit")):
    P("    mu_poly=%.5f (%s)  =>  mu_cyc >= %.5f  (gamma_cyc >= %.5f)  %s e"
      % (mu, tag, mu ** 0.25, LN(mu) / 4, "<" if mu ** 0.25 < E else ">"))
P("  tube family gives mu_cyc <= 1.6818 too => the closed class is NOT excluded by these bounds")

# ---------------------------------------------------------------- S6
P("")
P("== S6  acid test ==")
P("  model             defect-class content                mu_eff of the class        m")
P("  free field        no defects (F=0, w=0)               same lattice class 4.06    0")
P("  2D YM             plaquettes of Z^2 = polyominoes     4.0626 (exact)             3g^2/8 > 0")
P("  M6 (HO (x) free)  0+1 model, no plaquettes            undefined                  0")
P("  single plaquette  one cluster, n=1                    mu = 1 (trivial)           undef")
P("  single-plaquette test: z(1)=exp(-(c0-0.5)) <= 1/2 <=> c0 >= 0.5+ln2 = %.6f" % (0.5 + math.log(2)))
P("    at c0=%.6f: z(1)=%.6f (pass) ; at membrane window c0=%.6f: z(1)=%.6f (pass)"
  % (C0, math.exp(-(C0 - 0.5)), 0.5 + LN(4.0626) + LN3, math.exp(-(0.5 + LN(4.0626) + LN3 - 0.5))))
P("")
P("== total elapsed %.1f s ==" % (time.time() - T0))
with open("ge1_mueff.log", "w", encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")





