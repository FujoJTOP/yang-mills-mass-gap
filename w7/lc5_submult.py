# lc5_submult.py -- LC5: does a_n (connected plaquette clusters through a fixed
# plaquette, 4D) satisfy submultiplicativity?  Fekete bounds?  Can lambda* <= 48
# be proved from finitely many terms?
# a_n corpus: w7/ur4_enum_n6.log:19-21  (n=1..6) = C03-BYPASS:111 for n<=5
import math, itertools, sys

A = [None, 1, 20, 458, 11132, 280175, 7218768]      # 1-indexed, n=1..6
N = 6
out = []
def p(s=""):
    out.append(s); print(s)

p("LC5  submultiplicativity / Fekete / lambda* <= 48")
p("corpus a_n (n=1..6): %s   [ur4_enum_n6.log:19-21]" % A[1:])

# ---- 1. roots and ratios -------------------------------------------------
p("\n[1] n   a_n          a_n^(1/n)     a_n/a_(n-1)")
roots, ratios = {}, {}
for n in range(1, N + 1):
    r = A[n] ** (1.0 / n); roots[n] = r
    q = (A[n] / A[n - 1]) if n >= 2 else float('nan'); ratios[n] = q
    p("    %d  %12d  %12.6f  %14.6f" % (n, A[n], r, q))
p("    roots monotone increasing: %s" % all(roots[n] < roots[n+1] for n in range(1, N)))
p("    ratios monotone increasing: %s" % all(ratios[n] < ratios[n+1] for n in range(2, N)))

# cross-check the 15.58 claim
p("\n[1b] user claims n=6 gives 15.58 -- CHECK")
p("     a_6^(1/6)            = %.6f   [tn1_enum_n6_final.log:14 = 13.901997]" % roots[6])
p("     15.58^6              = %.6e   vs a_6 = %d" % (15.58 ** 6, A[6]))
p("     m_6^(1/6), m_n=sum|phi| = 15.580645  [ur4_enum_n6.log:12,21]  -> MISATTRIBUTED")
p("     => 15.58 is the root of m_n (sum|phi|), NOT of a_n")

# ---- 2. sub- vs super-multiplicativity on all pairs m+n <= 6 -------------
p("\n[2] pair test a_(m+n) vs a_m*a_n  (all pairs, m+n<=%d)" % N)
sub_fail, sup_fail, worst_sub, worst_sup = [], [], (0.0, None), (float('inf'), None)
for m, n in itertools.combinations_with_replacement(range(1, N), 2):
    s = m + n
    if s > N: continue
    t = A[m] * A[n]
    sub_ratio = A[s] / t          # submultiplicative  <=> sub_ratio <= 1
    p("    (%d,%d) sum=%d  a_%d=%-9d  a_%d*a_%d=%-9d  ratio=%.4f  %s"
      % (m, n, s, s, A[s], m, n, t, sub_ratio, "SUB-OK" if sub_ratio <= 1 else "SUB-FAIL"))
    if sub_ratio > 1: sub_fail.append((m, n))
    if sub_ratio < 1: sup_fail.append((m, n))
    if sub_ratio > worst_sub[0]: worst_sub = (sub_ratio, (m, n))
    if sub_ratio < worst_sup[0]: worst_sup = (sub_ratio, (m, n))
npair = len([1 for m, n in itertools.combinations_with_replacement(range(1, N), 2) if m + n <= N])
p("    pairs tested = %d ; SUB-FAIL = %d ; SUPER-FAIL = %d"
  % (npair, len(sub_fail), len(sup_fail)))
p("    worst submultiplicativity violation: a_%d/(a_%d*a_%d) = %.6f"
  % (worst_sub[1][0] + worst_sub[1][1], worst_sub[1][0], worst_sub[1][1], worst_sub[0]))
p("    smallest supermultiplicativity margin: %.6f  at (m,n)=%s"
  % (worst_sup[0], worst_sup[1]))
p("    VERDICT: a_(m+n) <= a_m*a_n is FALSE (%d/%d pairs fail); "
  "a_(m+n) >= a_m*a_n holds on %d/%d pairs  [E]" % (len(sub_fail), npair, npair - len(sup_fail), npair))

# a_1 = 1 alone kills exact submultiplicativity
p("    trivial killer: a_1 = 1  =>  sub would force a_n <= 1 for all n, but a_2 = 20")

# ---- 3. best constant C:  a_(m+n) <= C a_m a_n  needs C >= ? -------------
p("\n[3] constant-C form  a_(m+n) <= C*a_m*a_n :  observed lower bound on C")
p("    C >= max_observed ratio = %.6f   (only covers m+n<=6 -> NOT a theorem)" % worst_sub[0])

def fk_upper(C, n):        # Fekete upper bound IF uniform C holds: min_k (C a_k)^(1/k)
    return min((C * A[k]) ** (1.0 / k) for k in range(1, n + 1))

p("    hypothetical: if C were uniformly = %.4f then lambda* = inf_n (C a_n)^(1/n):" % worst_sub[0])
for n in range(1, N + 1):
    p("      n=%d  (C a_n)^(1/n) = %.6f" % (n, (worst_sub[0] * A[n]) ** (1.0 / n)))
p("      inf over n<=6 = %.6f   (< 48, but C is NOT proven)" % fk_upper(worst_sub[0], N))

# ---- 4. what C would suffice for lambda* <= 48 --------------------------
p("\n[4] SUFFICIENT (and weak) lemma:  C_48(n) = 48^n / a_n")
p("    if a_(m+n) <= C a_m a_n with C <= max_n C_48(n), then lambda* <= 48")
best = (0, 0)
for n in range(1, N + 1):
    c = 48 ** n / A[n]
    p("      n=%d  C_48 = %.3f   ->  (C_48 a_n)^(1/n) = %.6f" % (n, c, (c * A[n]) ** (1.0 / n)))
    if c > best[0]: best = (c, n)
p("    max C_48 = %.3f at n=%d  => uniform C <= %.0f closes lambda*<=48" % (best[0], best[1], best[0]))
p("    (with n=4 alone: C <= %.3f suffices)" % (48 ** 4 / A[4]))

# ---- 5. unconditional tree envelope (repo route) ------------------------
p("\n[5] unconditional ceiling in repo: a_n <= r_n = C(n*Delta, n-1)/n, Delta=20")
D = 20
for n in range(2, N + 1):
    r = math.comb(n * D, n - 1) / n
    p("    n=%d  r_n^(1/n) = %.6f   a_n^(1/n) = %.6f   a_n<=r_n: %s"
      % (n, r ** (1.0 / n), roots[n], A[n] <= r))
p("    r_n^(1/n) -> e*Delta = %.6f  [T]  =>  lambda* <= 54.3656 (known ceiling, > 48)"
  % (math.e * D))
p("    gap to close: 54.3656 / 48 = %.4f" % (math.e * D / 48.0))

# ---- 5b. local ratio route (the only route that can hit 48) -------------
p("\n[5b] local ratio lemma  a_(n+1) <= B a_n  (B explicit)  =>  lambda* <= B  [no Fekete]")
p("     observed max ratio = %.6f ; target B = 48 ; slack 48/%.3f = %.4f"
  % (ratios[6], ratios[6], 48.0 / ratios[6]))
p("     evidence ratios: %s  (all < 48 already at n=2)"
  % ", ".join("%.2f" % ratios[n] for n in range(2, N + 1)))

# ---- 6. verdict --------------------------------------------------------
p("\n[6] VERDICT")
p("    (i)   exact submultiplicativity: REFUTED [T]  (a_1=1; worst violation %.2f)" % worst_sub[0])
p("    (ii)  sequence is empirically SUPERmultiplicative on all %d pairs tested [E]" % npair)
p("    (iii) Fekete(super) => lambda* = sup_n a_n^(1/n) = lim exists;")
p("          finite terms give LOWER bounds only: lambda* >= a_6^(1/6) = %.6f" % roots[6])
p("    (iv)  Fekete(sub or C-form) would give UPPER bounds; neither hypothesis holds:")
p("          sub is false; uniform C is unproven (and observed C >= %.2f is rising)." % worst_sub[0])
p("    (v)   => (b): finitely many terms + Fekete CANNOT prove lambda* <= 48.")
p("          lambda*<=48 needs a uniform-glue lemma  a_(m+n) <= %.0f a_m a_n  [C]," % best[0])
p("          or a local ratio lemma a_(n+1) <= B a_n with B<=48 [T]-level bookkeeping.")
p("    (vi)  more terms (n=7..12) would only raise the LOWER bound (to ~26-30);")
p("          they can never yield an upper bound. No finite n suffices.")

with open("w7/lc5_run.log", "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")
print("\n[written] w7/lc5_run.log")
