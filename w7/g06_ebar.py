"""g06 -- is `e-bar` (= u-bar/u = e*D*C-bar) a-independent, i.e. is
delta_A = ln(1/u-bar)/a EXACTLY proportional to 1/a ?

Target (f05's only real danger): delta_A = kappa/4 - ln(e-bar)/a is a
DIFFERENCE OF TWO 1/a TERMS; if e-bar carried an O(a) piece, the 1/a
would amplify it into an O(1) term and flip the sign of delta_A.

Checks (ASCII output only):
  S1  e-bar := u-bar/u = e*D*C-bar ; numeric agreement with 8/beta_1
  S2  lattice fact: plaquette-adjacency degree D = 20 (pure integer)
  S3  animal counts a_n = 1,20,458,11132,... (D enters only here)
  S4  exactness: delta_A(a)*a independent of a to 50 dps ; identity residual
  S5  danger quantification: e-bar(1+c1 a) -> O(1) shift -c1 ; flip threshold
  S6  other (non-1/a) amplifiers: O(u) corrections in c_f, and C-hat/C-bar
  S7  verdict table
"""
import sys, math
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 50
N = 2
D = 20                                   # 4 links/plaquette x (6-1) others = 20
eD = mp.e * D                            # = 54.3656365691809...
Cbar_sum = mp.mpf(1)                     # sum-level value (C01 sec.4 / C07)
Cbar_urcell_hi = mp.exp(mp.mpf(1))       # C02 pessimistic reading C-bar<=e
e_bar = eD * Cbar_sum
beta = mp.mpf(1) / 12
cf = beta / (2 * N)
u = cf / N                               # u := c_f/N = beta/(2N^2) = 1/96
u_round = mp.mpf('0.0104')               # the truncated u used in f05's [E] line

L = []
def p(*a):
    s = " ".join(str(x) for x in a)
    print(s); L.append(s)

# ------------------------------------------------------------------ S1
p("== S1. definition of e-bar and the 8/e-bar round trip ==")
p("  definition (Z07:84,100):  e-bar := ub/u = e*D*C-bar ,  D=20, C-bar=1 (sum level)")
p("  e*D                = %s" % mp.nstr(eD, 20))
p("  corpus number (Z12:97, Y04:117, C03:126): 54.3656365692")
p("  |e*D - corpus|     = %s" % mp.nstr(abs(eD - mp.mpf('54.3656365692')), 6))
p("  beta_1 = 8/e-bar   = %s   (Z12:97 gives 1.4715178e-1)" % mp.nstr(8 / e_bar, 12))
p("  8/0.147152         = %s" % mp.nstr(8 / mp.mpf('0.147152'), 12))
p("  8/0.14715          = %s   (f05 quoted the 5-digit window 0.14715)" % mp.nstr(8 / mp.mpf('0.14715'), 12))
p("  round trip 8/(8/e-bar) - e-bar = %s" % mp.nstr(8 / (8 / e_bar) - e_bar, 6))
p("  => e-bar(from 0.147152) = 54.3656  vs  e*D = %s : agree to 6 s.f.;" % mp.nstr(eD, 8))
p("     the 54.366 from 8/0.14715 is pure 5-digit rounding of beta_1 (err %s)."
  % mp.nstr(8 / mp.mpf('0.14715') - eD, 6))

# ------------------------------------------------------------------ S2
p("== S2. D = 20 is a combinatorial integer (no a) ==")
Lt = 8                                    # torus size; cluster spread <= 4 < L
plaq = []
for x in range(Lt):
    for y in range(Lt):
        for z in range(Lt):
            for w in range(Lt):
                X = (x, y, z, w)
                for mu in range(4):
                    for nu in range(mu + 1, 4):
                        plaq.append((X, mu, nu))
linkof = {}
def links_pl(q):
    X, mu, nu = q
    out = []
    e_mu = [0, 0, 0, 0]; e_mu[mu] = 1
    e_nu = [0, 0, 0, 0]; e_nu[nu] = 1
    A = tuple(X); B = tuple((X[k] + e_mu[k]) % Lt for k in range(4))
    C = tuple((X[k] + e_nu[k]) % Lt for k in range(4))
    out.append((A, mu)); out.append((B, nu)); out.append((A, nu)); out.append((C, mu))
    return sorted(out)
for i, q in enumerate(plaq):
    for lk in links_pl(q):
        linkof.setdefault(lk, []).append(i)
deg = set()
for lk, ps in linkof.items():
    deg.add(len(ps))
p("  #plaquettes on L=%d torus: %d ; #links: %d" % (Lt, len(plaq), len(linkof)))
p("  plaquettes per link: %s   (= 2(d-1) = 6)" % sorted(deg))
root = 0
nb = set()
for lk in links_pl(plaq[root]):
    nb |= set(linkof[lk])
nb.discard(root)
p("  adjacency degree of a root plaquette: %d   (= 4*(6-1) = 20)" % len(nb))
p("  => D = %d, integer, independent of the lattice spacing a  [T]" % len(nb))
adj = {}
for i, q in enumerate(plaq):
    s = set()
    for lk in links_pl(q):
        s |= set(linkof[lk])
    s.discard(i)
    adj[i] = s
p("  (cross-check: min/max degree over all plaquettes = %d / %d)"
  % (min(len(v) for v in adj.values()), max(len(v) for v in adj.values())))

# ------------------------------------------------------------------ S3
p("== S3. animal counts through a fixed plaquette (D enters only here) ==")
cur = {frozenset([root])}
counts = [1]
for n in range(2, 6):
    nxt = set()
    for S in cur:
        for q in S:
            for r in adj[q]:
                if r not in S:
                    nxt.add(frozenset(tuple(sorted(S)) + (r,)))
    cur = nxt
    counts.append(len(cur))
p("  a_n (n=1..5): %s" % counts)
p("  corpus (C03:111): 1, 20, 458, 11132, 280175")
p("  match: %s" % (counts == [1, 20, 458, 11132, 280175]))
r = [mp.mpf(c) ** (mp.mpf(1) / (i + 1)) for i, c in enumerate(counts)]
p("  a_n^(1/n): %s   -> e*D = %s is a valid envelope" % ([mp.nstr(x, 6) for x in r], mp.nstr(eD, 6)))

# ------------------------------------------------------------------ S4
p("== S4. EXACTNESS: delta_A(a)*a independent of a (50 dps) ==")
def delta_A(a, ebar=e_bar, uu=u):
    return mp.log(1 / uu) / a - mp.log(ebar) / a
ref = None
p("     a            delta_A(a)                delta_A(a)*a          dev")
for a in (mp.mpf(1), mp.mpf('0.5'), mp.mpf('0.1'), mp.mpf('0.01'),
          mp.mpf('1e-3'), mp.mpf('1e-6'), mp.mpf('1e-12')):
    v = delta_A(a) * a
    if ref is None:
        ref = v
    p("  %-12s %-25s %-21s %s" % (mp.nstr(a, 6), mp.nstr(delta_A(a), 12),
                                  mp.nstr(v, 20), mp.nstr(v - ref, 3)))
p("  ub = e-bar*u = %s ;  ln(1/ub) = %s" % (mp.nstr(e_bar * u, 10), mp.nstr(mp.log(1 / (e_bar * u)), 20)))
p("  identity residual  delta_A(a)*a - ln(1/ub)  at a=1: %s"
  % mp.nstr(delta_A(mp.mpf(1)) - mp.log(1 / (e_bar * u)), 5))
p("  kappa/4*a = ln(1/u) = %s   (u = 1/96 exactly)" % mp.nstr(mp.log(1 / u), 20))
p("  same with the truncated u=0.0104 (f05's [E] line): kappa/4*a = %s , delta_A*a = %s"
  % (mp.nstr(mp.log(1 / u_round), 12), mp.nstr(delta_A(mp.mpf(1), uu=u_round), 12)))

# ------------------------------------------------------------------ S5
p("== S5. what the danger would look like: e-bar -> e-bar(1+c1*a) ==")
Lnum = mp.log(1 / (e_bar * u))
p("  L := ln(1/ub_0) = %s  (the pure numerator)" % mp.nstr(Lnum, 12))
p("  e-bar(a) = e-bar_0*(1+c1*a)  =>  delta_A(a)*a = L - ln(1+c1*a):")
p("  NOT a constant any more; the extra piece is -ln(1+c1*a) = O(1) at a=1.")
f = mp.exp(Lnum) - 1
p("  exact flip threshold:  c1*a > e^L - 1 = 1/ub_0 - 1 = %s  (so c1 > %s at a=1)"
  % (mp.nstr(f, 10), mp.nstr(f, 10)))
p("     a        c1*a (c1=1)   delta_A*a (c1=1)   sign")
for a in (mp.mpf(1), mp.mpf('0.7'), mp.mpf('0.5'), mp.mpf('0.1'), mp.mpf('1e-2')):
    v = Lnum - mp.log(1 + mp.mpf(1) * a)
    p("  %-8s %-13s %-19s %s" % (mp.nstr(a, 4), mp.nstr(a, 8), mp.nstr(v, 12),
                                 "+" if v > 0 else "-"))
for c1 in (mp.mpf('0.5'), mp.mpf('0.766'), mp.mpf(1), mp.mpf(3)):
    p("  c1=%-6s -> sign flips for a > %s   (a=1: delta_A = %s)"
      % (mp.nstr(c1, 4), mp.nstr(f / c1, 8),
         mp.nstr(delta_A(mp.mpf(1), ebar=e_bar * (1 + c1)), 12)))
p("  => the danger, IF e-bar had an O(a) piece, is LIVE AT a=1 (where beta_1=0.147")
p("     and delta_A*a=0.5686 are quoted): c1 > 0.7658 (a 77% lattice artifact in")
p("     the animal constant) flips the sign.  c1 must be 0 EXACTLY, not 'small'.")
p("     same threshold if the O(a) sits in u instead: numerator = L - ln(1+c*a).")

# ------------------------------------------------------------------ S6
p("== S6. the OTHER amplifiers (no 1/a, hence harmless) ==")
p("  (i) c_f = beta/(2N) + O(beta^3) => u = beta/(2N^2) + O(beta^3):")
p("      this shifts the NUMERATOR by O(u), relative size u/L = %s ; margin L/u = %s"
  % (mp.nstr(u / Lnum, 6), mp.nstr(Lnum / u, 8)))
p("      no 1/a factor => cannot flip the sign (numerator stays ~0.57 >> u=0.0104).")
p("  (ii) C-hat (Z12=4) and C-bar (>=1):  numerator becomes ln(1/(C-hat*e-bar*C-bar))")
for btag, uu in (('beta=1/12 (u=1/96)', u), ('beta=1/48 (u=1/384)', mp.mpf(1) / 384)):
    for Ch, Cb, tag in ((mp.mpf(1), mp.mpf(1), 'C-hat=1, C-bar=1  (C03 route)'),
                        (mp.mpf(4), mp.mpf(1), 'C-hat=4, C-bar=1  (Z12 window)'),
                        (mp.mpf(4), mp.exp(mp.mpf(1)), 'C-hat=4, C-bar=e  (C02 pessim.)')):
        num = mp.log(1 / uu) - mp.log(Ch * eD * Cb)
        win = 8 / (Ch * eD * Cb)
        p("      %-18s %-28s numerator=%-13s win beta<%-10s sign=%s"
          % (btag, tag, mp.nstr(num, 8), mp.nstr(win, 6), "+" if num > 0 else "-"))
p("      => the numerator IS a pure number in every reading, but its SIGN is a")
p("         beta-window statement, not an a-statement: beta=1/12 sits OUTSIDE the")
p("         C-hat=4 window (beta_1^fix=3.679e-2 < 0.0833, Z12:161-162), that is why")
p("         the C-hat=4 row at beta=1/12 is negative.  a-dependence: none anywhere.")
p("         => delta_A = (pure number)/a in every in-window reading; no O(1) remainder.")

# ------------------------------------------------------------------ S7
p("== S7. verdict ==")
p("  (1) e-bar = u-bar/u = e*D*C-bar = e*20*C-bar ; = %s with C-bar=1;" % mp.nstr(e_bar, 12))
p("      8/e-bar = %s = beta_1 ; 8/0.147152 = 54.3656 : consistent [E]" % mp.nstr(8 / e_bar, 10))
p("  (2) D=20 combinatorial integer, C-bar pure (tree=animal identity), e constant")
p("      => e-bar is a pure number, a-INDEPENDENT [T]")
p("  (3) delta_A = (1/a)[ln(1/u) - ln e-bar] = ln(1/ub)/a , numerator pure [T]")
p("      => delta_A is EXACTLY proportional to 1/a ; no O(1) remainder (S4: dev 0 to 50 dps)")
p("  (4) f05's danger does NOT exist ; but it is not vacuous: c1 > 0.7658 at a=1 flips (S5).")
p("  (5) residual a-dependence is only: 1/a^p conversions, L/a counting, e^{-cL/a}.")

open('g06_ebar.log', 'w', encoding='utf-8').write("\n".join(L) + "\n")
p("  [log written: g06_ebar.log, %d lines]" % len(L))
