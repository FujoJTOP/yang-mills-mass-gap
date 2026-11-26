# ENV1 — is L_t = Delta^Delta/(Delta-1)^(Delta-1) the "optimized lambda_a"?
#   relation L_t/(e Delta), legality of the substitution, gap to 48, best constant.
# py -3.12 + mpmath.  writes nothing; prints to stdout (tee'd into env1_optimized_lambda.log)
import mpmath as mp

mp.mp.dps = 30
D = mp.mpf(20)

eD = mp.e * D                                   # e*Delta  = "unoptimized animal constant"
Lt = D**D / (D - 1)**(D - 1)                    # rooted-tree growth constant
Lt2 = D * (1 + 1 / (D - 1))**(D - 1)            # same thing, binomial form


def p(s=""):
    print(s)


p("=" * 78)
p("ENV1  Delta = 20   (4D plaquette adjacency degree)")
p("=" * 78)
p("  e*Delta        = %s   (corpus 54.3656365692)" % mp.nstr(eD, 20))
p("  L_t            = %s   (corpus 53.0006865328)" % mp.nstr(Lt, 20))
p("  L_t alt form   = %s   (Delta*(1+1/(Delta-1))^(Delta-1))" % mp.nstr(Lt2, 20))
p("  |L_t - alt|    = %s" % mp.nstr(abs(Lt - Lt2), 3))

p("\n== S1  the relation  L_t/(e*Delta) ==")
r = Lt / eD
p("  L_t/(e*Delta)  = %s" % mp.nstr(r, 20))
p("  => L_t = e*Delta*(1 - %.4f%%)" % (100 * (1 - r)))
p("  => e*Delta = L_t*(1 + %.4f%%)" % (100 * (1 / r - 1)))
p("  first-order asymptotic 1/(2(Delta-1)) = %s  (= the '2.6%%' folklore)" % mp.nstr(1 / (2 * (D - 1)), 8))
p("  exact reduction                        = %s" % mp.nstr(1 - r, 8))
p("  >>> cb6's 'Delta*e*(1-2.6%%)' = %s   <-- that is L_t itself, NOT the ratio"
  % mp.nstr(eD * (1 - 1 / (2 * (D - 1))), 12))
p("  >>> ratio slot: (1 + 1/19)^19/e = %s  = (Delta/(Delta-1))^(Delta-1)/e  [cb6 CB6:102]"
  % mp.nstr((D / (D - 1))**(D - 1) / mp.e, 12))

p("\n== S2  legality: (1+1/m)^m < e for every finite m  =>  L_t < e*Delta strictly ==")
for m in (1, 2, 5, 19, 100):
    p("  m=%-4d (1+1/m)^m = %s   vs e = %s" % (m, mp.nstr((1 + mp.mpf(1) / m)**m, 10), mp.nstr(mp.e, 10)))
p("  same route, same object:  rooted-tree Lagrange count")
p("  r_n := (1/n)*C(n*Delta, n-1) = #rooted trees of size n on degree-Delta graph;  r_n^(1/n) -> L_t")
ladder = [(2, 4.472135955), (5, 15.09698946), (10, 25.52820276), (20, 34.93181899),
          (50, 43.64702383), (200, 49.96763161), (1000, 52.25342124), (5000, 52.82487289)]
for n, v in ladder:
    p("    n=%-5d (corpus) %s" % (n, v))
for n in (10000, 40000, 160000):
    lv = (mp.loggamma(D * n + 1) - mp.loggamma(n) - mp.loggamma((D - 1) * n + 2) - mp.log(n)) / n
    p("    n=%-6d  (recomputed) r_n^(1/n) = %s   L_t - this = %s"
      % (n, mp.nstr(mp.e**lv, 12), mp.nstr(Lt - mp.e**lv, 6)))
p("  >> monotone increasing to L_t with no slack left: the Lagrange/Lagrange-type count is TERMINAL at L_t")

p("\n== S3  gap to 48 ==")
p("  e*Delta/48      = %s" % mp.nstr(eD / 48, 12))
p("  L_t/48          = %s" % mp.nstr(Lt / 48, 12))
p("  (L_t/48)/(e*Delta/48) = %s" % mp.nstr(r, 12))
p("  needed further reduction from L_t : %s" % mp.nstr(1 - 48 / Lt, 8))
p("  needed further reduction from e*Delta: %s" % mp.nstr(1 - 48 / eD, 8))
p("  beta-window 2N^2/lambda (N=2, Cbar=1): e*Delta -> %s ; L_t -> %s ; 48 -> %s (=1/12)"
  % (mp.nstr(8 / eD, 10), mp.nstr(8 / Lt, 10), mp.nstr(8 / 48, 10)))
p("  beta-window 4/lambda (Cbar=N=2 reading): e*Delta -> %s ; L_t -> %s"
  % (mp.nstr(4 / eD, 10), mp.nstr(4 / Lt, 10)))

p("\n== S4  what would be needed: Delta* <= with L_t(Delta*) = 48 ==")
f = lambda x: x * (1 + 1 / (x - 1))**(x - 1) - 48
Ds = mp.findroot(f, 18.2)
p("  Delta*(continuous) = %s   (Delta = 20 is 1.10x too big)" % mp.nstr(Ds, 12))
for d in (18, 19, 20):
    d = mp.mpf(d)
    p("    Delta=%-3s L_t = %s" % (d, mp.nstr(d * (1 + 1 / (d - 1))**(d - 1), 12)))

p("\n== S5  the third number 4.8572: which slot? ==")
phi = (1 + mp.sqrt(5)) / 2
T2 = phi**2 * mp.e**(1 / phi)
p("  phi^2*e^(1/phi) = %s   (corpus 4.857178)" % mp.nstr(T2, 12))
h = lambda a: (a / (1 + a)) * mp.e**(-a)          # d/da: -1/a + 1/(1+a) + 1 = 0 <=> a^2+a-1 = 0
al = 1 / phi
best = max(h(mp.mpf(k) / 2000) for k in range(1, 8000))   # dense real scan: no complex branch
p("  alpha* = %s = 1/phi  (a^2+a-1=0) ;  scan max h = %s ;  h(1/phi) = %s"
  % (mp.nstr(al, 12), mp.nstr(best, 12), mp.nstr(h(al), 12)))
p("  C := max_alpha [alpha/(1+alpha)]e^-alpha = %s = 1/%s"
  % (mp.nstr(1 / T2, 10), mp.nstr(T2, 10)))
p("  => the criterion reads  lambda_a*K*u <= 1/%s  (a LAMBDA-K-U condition,"
  % mp.nstr(T2, 8))
p("     NOT a bound on lambda_a).  Reading it as 'lambda_a <= 4.8572' gives")
p("     beta-window 8/4.857178 = %s > 1/12 = %s  -> would falsely 'close the gap'"
  % (mp.nstr(8 / T2, 6), mp.nstr(mp.mpf(1) / 12, 6)))
p("  >> 4.8572 is the reciprocal of an alpha-optimized coefficient; different slot [T]")

p("\n== S6  verdict numbers ==")
p("  lambda_a := e*Delta = %s   (definition, P3-FULLREP-KP:89 'one may take')" % mp.nstr(eD, 12))
p("  best rigorous upper bound on it = L_t = %s   [T]" % mp.nstr(Lt, 12))
p("  48/53.000687 = %s  -> still short of 48 by 10.4%% ; verdict (b)"
  % mp.nstr(Lt / 48, 10))
p("  empirical probe (a_2..a_5 = 20,458,11132,280175): ratios 20,22.9,24.3,25.17 => lambda_* >~ 25.2")
p("  => L_t is 2.10x above the probed truth: legal but far from tight; 48 is out of reach of this route")
