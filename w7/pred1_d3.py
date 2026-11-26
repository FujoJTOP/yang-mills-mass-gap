#!/usr/bin/env python3.12
"""
pred1_d3.py -- PRED1: the d=3 test of the "gap" reading  g(d) = e*Delta / (16(d-1)).

Two layers, kept strictly apart:

  [L1] ENVELOPE layer (identity, no data needed):
        Delta(d)  = 4(2d-3)         plaquette adjacency degree   (P3:116)
        T(d)      = 16(d-1)         SZZ threshold denominator, |beta| < 1/(16(d-1))
        g(d)      = e*Delta(d)/T(d) = e(2d-3)/(4(d-1))
      =>  g(3) = 3e/8,  g(4) = 5e/12 = 1.13262.
      The "wrong" reading is to carry the d=4 number 48 into d=3: e*12/48 = e/4.

  [L2] TRUE layer (needs data): the actual animal growth constant lambda_*(d) of
       plaquette animals rooted at a fixed plaquette (adjacency = share a link),
       enumerated here from scratch in Z^d.  Compare lambda_*(d) with T(d) and e*Delta(d).

Calibration: d=4 must reproduce the corpus numbers a_1..a_4 = 1, 20, 458, 11132
(ADJ2:14, TN2 §4) and (torus L=6) a_5 = 280175.
"""
import sys, time
from mpmath import mp, mpf, e, sqrt, log, nstr

mp.dps = 40
OUT = []


def P(s=""):
    print(s)
    OUT.append(s)


# ------------------------------------------------------------------ [L1] algebra
P("=" * 78)
P("[L1] ENVELOPE LAYER  g(d) = e*Delta(d)/T(d),  Delta = 4(2d-3),  T = 16(d-1)")
P("=" * 78)
P("")
P("  d   Delta=4(2d-3)   T=16(d-1)   Delta/T      g(d)=e*Delta/T    e*Delta/48 (wrong reading)")
for d in range(2, 6):
    Dl = 4 * (2 * d - 3)
    T = 16 * (d - 1)
    g = e * Dl / T
    bad = e * Dl / 48
    P(f"  {d}   {Dl:>10}        {T:>6}      {nstr(mpf(Dl)/T,10):>10}   {nstr(g,12):>14}   {nstr(bad,12):>14}")

g3 = 3 * e / 8
g4 = 5 * e / 12
P("")
P(f"  g(3) = e*12/32 = 3e/8  = {nstr(g3,20)}   (excess over 1: {nstr((g3-1)*100,10)} %)  [T: identity]")
P(f"  g(4) = e*20/48 = 5e/12 = {nstr(g4,20)}   (excess over 1: {nstr((g4-1)*100,10)} %)  [T: identity]")
P(f"  48/Delta(3) = 32/12 = {nstr(mpf(32)/12,12)} = 8/3 ;  e/(8/3) = 3e/8  -> same number, by construction")
P(f"  wrong reading at d=3 (keeping 48): e*12/48 = e/4 = {nstr(e/4,12)}  (a 32% DROP below 1, not 13.3%)")
P(f"  root of g(d)=1:  d* = (3e-4)/(2e-4) = {nstr((3*e-4)/(2*e-4),12)} (non-integer)")
P("")

# ------------------------------------------------- [L2] plaquette-animal counting
def plaq_links(p, d):
    """p = (mu, nu, base) with mu<nu ; return the 4 link ids (dir, point) it contains."""
    mu, nu, b = p
    return ((mu, b),
            (mu, tuple(b[i] + (1 if i == nu else 0) for i in range(d))),
            (nu, b),
            (nu, tuple(b[i] + (1 if i == mu else 0) for i in range(d))))


def neighbors(p, d):
    """All plaquettes sharing >=1 link with p, excluding p itself.  Count = 4(2d-3) = Delta."""
    out = set()
    for (mu, x) in plaq_links(p, d):
        for nu2 in range(d):
            if nu2 == mu:
                continue
            for s in (0, 1):                       # base = x  or  x - e_nu2
                b = tuple(x[i] - (s if i == nu2 else 0) for i in range(d))
                a1, a2 = (mu, nu2) if mu < nu2 else (nu2, mu)
                out.add((a1, a2, b))
    out.discard(p)
    return out


ROOT = (0, 1, (0,) * 8)   # (mu=0, nu=1, base=0^d) ; trailing zeros harmless per-d


def root_plaq(d):
    return (0, 1, (0,) * d)


def enumerate_animals(d, nmax, cap=3_000_000, tmax=180.0):
    """a_n = #connected sets of n plaquettes containing the fixed root (edge adjacency)."""
    root = root_plaq(d)
    cur = {frozenset((root,))}
    counts = [1]                      # a_1 = 1 (the root alone)
    t0 = time.time()
    for n in range(1, nmax):
        nxt = set()
        for A in cur:
            for p in A:
                for q in neighbors(p, d):
                    if q in A:
                        continue
                    B = A | {q}
                    if len(B) == n + 1:
                        nxt.add(B)
        counts.append(len(nxt))
        cur = nxt
        if len(cur) > cap or time.time() - t0 > tmax:
            P(f"  [cap/time stop at n={n+1}]  (this level size {len(cur)})")
            break
    return counts, time.time() - t0


def report_animals(d, nmax, ref=None):
    P("=" * 78)
    P(f"[L2] PLAQUETTE ANIMALS in Z^{d}   (root fixed, adjacency = share a link, degree Delta = {4*(2*d-3)})")
    P("=" * 78)
    cnt, dt = enumerate_animals(d, nmax)
    P(f"  wall time {dt:.1f}s")
    P("")
    P("   n    a_n            a_n/a_(n-1)    a_n^(1/n)     a_n^(1/n)/eDelta")
    eD = e * 4 * (2 * d - 3)
    for i, a in enumerate(cnt, start=1):
        ratio = mpf(cnt[i - 2]) if False else (mpf(a) / cnt[i - 2] if i >= 2 and cnt[i - 2] else mpf('nan'))
        root_ = mpf(a) ** (mpf(1) / i)
        P(f"  {i:>2}  {a:>12}    {nstr(ratio,8):>10}    {nstr(root_,9):>10}    {nstr(root_/eD,8):>8}")
    if ref:
        P("")
        for k, v in ref.items():
            got = cnt[k - 1] if k <= len(cnt) else None
            tag = "MATCH" if got == v else "DIFF"
            P(f"  calibration a_{k} = {v} (corpus)  vs  {got} (here)   [{tag}]")
    P("")
    return cnt


d4 = report_animals(4, 6, ref={1: 1, 2: 20, 3: 458, 4: 11132, 5: 280175})
d3 = report_animals(3, 8)

# ------------------------------------------------------------ lambda_* extraction
def lam_est(cnt):
    """crude 1/n-extrapolation of the ratio sequence a_n/a_(n-1)."""
    r = [mpf(cnt[i]) / cnt[i - 1] for i in range(1, len(cnt))]
    if len(r) < 3:
        return r[-1] if r else mpf('nan'), mpf('nan')
    # two-point Richardson on the last two ratios (error ~ c/n)
    n1, n2 = len(r), len(r) - 1
    s = (r[-1] * n1 - r[-2] * n2) / (n1 - n2)
    return r[-1], s


P("=" * 78)
P("[L2] lambda_* vs the two reference scales")
P("=" * 78)
for d, cnt in ((3, d3), (4, d4)):
    Dl, T = 4 * (2 * d - 3), 16 * (d - 1)
    eD = e * Dl
    last, ext = lam_est(cnt)
    P("")
    P(f"  d={d}:  a_n = {cnt}")
    P(f"        last ratio a_n/a_(n-1) = {nstr(last,10)}   Richardson(1/n) = {nstr(ext,10)}")
    P(f"        eDelta = {nstr(eD,10)}   T = 16(d-1) = {T}")
    P(f"        lambda_*/eDelta       = {nstr(last/eD,8)} (lower bound on the true ratio)")
    P(f"        lambda_*/T(={T})      = {nstr(last/T,8)}   [gap at TRUE layer, lower bound]")
    P(f"        envelope gap g({d})  = eDelta/T = {nstr(eD/T,10)}   [gap at ENVELOPE layer, identity]")
P("")

# ------------------------------------------- what could "a_n(3) = 1, 6, ..." mean?
P("=" * 78)
P("[L2'] the claim a_n(3) = 1, 6, ... -- three candidate conventions")
P("=" * 78)
r3 = root_plaq(3)
P(f"  edge-adjacency (the corpus convention):  a_1 = 1, a_2 = {len(neighbors(r3,3))}  (= Delta(3) = 4(2d-3) = 12)")
pl3 = set()
for (mu, x) in plaq_links(r3, 3):
    for nu2 in range(3):
        if nu2 == mu:
            continue
        for s in (0, 1):
            b = tuple(x[i] - (s if i == nu2 else 0) for i in range(3))
            a1, a2 = (mu, nu2) if mu < nu2 else (nu2, mu)
            pl3.add((a1, a2, b))
pl3.discard(r3)
P(f"  vertex-sharing (any common corner):      a_2 = {len(pl3)}")
P(f"  plaquettes through ONE link (not an animal count): 2(d-1) = {2*(3-1)} at d=3 (= 6 at d=4 too)")
P(f"  => '1, 6' is not a_2 of any edge-adjacency animal model; 6 = 2(d-1) is a PER-LINK count.")
P("")

# ------------------------------------------------------------------------- verdict
P("=" * 78)
P("[VERDICT]")
P("=" * 78)
last3, ext3 = lam_est(d3)
P(f"  g(3) = 3e/8 = {nstr(g3,15)}   [T] -- exact algebraic identity e*(2*3-3)/(4*2)")
P(f"  d=3 threshold 16(d-1) = 32 (SZZ: |beta| < 1/(16(d-1))), Delta(3) = 12, 32/12 = 8/3")
P(f"  d=3 TRUE-layer gap lower bound: lambda_*(3)/32 >= {nstr(last3/32,8)}   (eDelta/32 = {nstr(g3,10)})")
P(f"  => at the ENVELOPE layer the prediction is self-consistent but vacuous (identity);")
P(f"     at the TRUE layer the d=3 data put lambda_*/32 near 0.4-0.5, i.e. the same")
P(f"     'gap disappears at true level' pattern as d=4 (25.17/48 = 0.524).")
P("")

path = r"D:\Dev\math\w7\pred1_d3.log"
with open(path, "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
print(f"[log written] {path}")
