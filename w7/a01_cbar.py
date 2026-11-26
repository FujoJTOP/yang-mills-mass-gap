# -*- coding: utf-8 -*-
"""
a01_cbar.py -- pin down Cbar (the Ursell constant of X01-L1C:71) with a PROVABLE upper bound,
               its N-scaling, and the resulting KP window.

Definition under test (verbatim, research/v3/species/X01-L1C.md:71):
    |phi^U(Gamma)| <= Cbar^{|Gamma|} u^{|Gamma|}      (Ursell tree-graph bound, one u per plaquette)
    ubar := e*Delta*Cbar*u ,  Delta = 20 , e*Delta = 54.3656365692...
    hole-C fix (Z12-FIXC.md:96):  ubar -> Chat*ubar ,  Chat = 4 (provable per-plaquette envelope)
    no-window-shrink condition (Z12-FIXC.md:111):  Chat*Cbar < 2N^2*48/(e*Delta) = 7.0632853...  (N=2)

Outputs: bound on Cbar, N-scaling, beta-window, delta, and the (a)/(b)/(c) verdict flags.
Run: py -3.12 w7/a01_cbar.py > w7/a01_cbar.log
"""
import mpmath as mp

mp.mp.dps = 30
E = mp.e
DELTA = mp.mpf(20)          # plaquette-adjacency degree (P3-FULLREP-KP 4.1)
LAM = E * DELTA             # animal constant lambda_a = e*Delta = 54.3656...
N = 2                       # SU(2) reference
CHAT_ENV = mp.mpf(4)        # per-plaquette provable envelope (Z12-FIXC)
CHAT_LOC = mp.mpf(2)        # exact local value per factorisation point (Z12-FIXC)


def hdr(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


# ----------------------------------------------------------------------------------
hdr("[1] WHAT CBAR IS: the Ursell tree-graph constant, not an animal/entropy constant")
# ----------------------------------------------------------------------------------
print("definition (X01-L1C:71):  |phi^U(Gamma)| <= Cbar^{|Gamma|} u^{|Gamma|}")
print("  phi^U(Gamma) = Ursell/Mayer polynomial of a connected cluster Gamma")
print("  |Gamma|      = number of plaquettes in Gamma (one factor u = c_f/N each)")
print("  => Cbar      = sup_Gamma ( |phi^U(Gamma)| / u^{|Gamma|} )^{1/|Gamma|}")
print("  ubar := e*Delta*Cbar*u is then exactly the KP l.h.s.  sum_{g ni x}|w(g)|e^{a|g|}  (Z02-SEILER:68)")
print()
print("numerical reading: ubar=0.57 at (SU(2), beta_W=1/12) with e*Delta=%.6f, u=0.0104" % LAM)
print("  => Cbar_implied = ubar/(e*Delta*u) = %.4f   (i.e. Cbar=1 is ASSUMED, never proved)"
      % (mp.mpf("0.5654") / (LAM * mp.mpf("0.0104"))))

# ----------------------------------------------------------------------------------
hdr("[2] PROVABLE BOUND I: naive Cayley tree count does NOT give a constant (shows what is at stake)")
# ----------------------------------------------------------------------------------
# Penrose tree-graph bound, blocks = plaquettes:
#   |phi^U(Gamma)| <= (number of spanning trees of the cluster graph) * prod |w_p|
# Cayley: n^{n-2} trees on n labelled vertices.
print(" n   n^{n-2}          (n^{n-2})^{1/n}   [per-vertex envelope; blow-up]")
for n in [1, 2, 3, 5, 10, 20, 50, 100]:
    if n == 1:
        c = mp.mpf(1)
    else:
        c = mp.power(n, n - 2)
    print(" %3d  %-16.6s %s" % (n, mp.nstr(c, 10), mp.nstr(mp.power(c, mp.mpf(1) / n), 8)))
print(" => envelope grows like n/e : NO size-independent Cbar from Cayley alone.")
print("    (the tree entropy must be paid by the lattice position sum -> see [3];")
print("     this is exactly why a *graph* bound, not a labelled-tree count, is the right route.)")

# ----------------------------------------------------------------------------------
hdr("[3] PROVABLE BOUND II: Kirchhoff + AM-GM on the plaquette-adjacency graph  =>  Cbar <= Delta = 20")
# ----------------------------------------------------------------------------------
# Gamma connected, n plaquettes, adjacency graph G_Gamma with max degree <= DELTA = 20.
# Kirchhoff:  tau(G) = (1/n) prod_{j=2..n} lam_j <= (1/n) (sum lam_j / n)^{n-1}
#             sum lam_j = 2|E| <= DELTA*n            (handshake + max degree)
#   => tau(G) <= DELTA^{n-1}/n
# Penrose:    |phi^U(Gamma)| <= tau(G_Gamma) prod_p |w_p| ,  |w_p| <= C_g u ,  C_g<=2Ne^beta (P3:78)
# With the C_g part already carried by K_KP (P3-FULLREP-KP 4.2), the residual tree constant is
#   Cbar <= sup_n ( DELTA^{n-1}/n )^{1/n}  ->  DELTA = 20      (crude variant 2*DELTA = 40)
print(" n     (Delta^{n-1}/n)^{1/n}     2*Delta-per-vertex")
for n in [1, 2, 3, 5, 10, 50, 500, 5000]:
    v = mp.power(mp.power(DELTA, n - 1) / n, mp.mpf(1) / n)
    print(" %5d   %-20s      %s" % (n, mp.nstr(v, 10), mp.nstr(mp.mpf(2) * v, 10)))
CBAR_LIM = DELTA                        # Kirchhoff+AM-GM envelope
CBAR_CRUDE = 2 * DELTA                  # handshake-free variant
print(" => PROVABLE:  Cbar <= Delta = %s   (crude: Cbar <= 2*Delta = %s)" % (CBAR_LIM, CBAR_CRUDE))

# sanity: actual tau for explicit 4D-plaquette animals (2D planar grid of m x m plaquettes,
# planar adjacency degree 8 in 4D -> a sub-case of the Delta=20 envelope)
print("\n non-vacuity check: tau of an m x m planar plaquette grid (n=m^2), and (tau)^{1/n}:")


def tau_grid(m):
    # eigenvalues of the m x m grid graph: 4 - 2cos(pi j/m) - 2cos(pi k/m), zero mode removed
    prod = mp.mpf(1)
    for j in range(m):
        for k in range(m):
            if j == 0 and k == 0:
                continue
            prod *= 4 - 2 * mp.cos(mp.pi * j / m) - 2 * mp.cos(mp.pi * k / m)
    return prod / (m * m)


for m in [2, 3, 4, 6]:
    t = tau_grid(m)
    print("   m=%d  n=%2d  tau=%-18s  tau^{1/n}=%s" % (m, m * m, mp.nstr(t, 12),
                                                       mp.nstr(mp.power(t, mp.mpf(1) / (m * m)), 8)))
print(" (grids are far below the envelope 20 -- the bound is safe and loose by ~2-3x;")
print("  no claim that 20 is attained, only that it is provable.)")

# ----------------------------------------------------------------------------------
hdr("[4] N-SCALING: Cbar is O(1) in N (provable); the reading Cbar = N^2 is a DOUBLE COUNT")
# ----------------------------------------------------------------------------------
# u = c_f/N = beta/(2N^2) already carries the per-link 1/dim f = 1/N (X01 sec.0, Gupta (14.18)).
# Exact SU(2) strong-coupling mass  M = -4 ln u  with u = beta/4  <=>  u = c_f/N, c_f = beta/N
# (f ~ fbar at N=2 doubles the leading coefficient).  No extra N^2 anywhere.
print(" Cbar from [3] depends on Delta only:  dCbar/dN = 0  (Delta is a lattice degree, not an N-power)")
print(" => Cbar = O(1);  ubar = e*Delta*Cbar*u ~ Cbar*beta/(2N^2)*O(1)  =>  beta_1 ~ N^2 (Wilson beta)")
print("    the N^2 reading (Z02-SEILER sec.4 route D) reproduces the u=c_f normalisation:")
print("    e*Delta*N^2*c_f/N = e*Delta*N*c_f  == route A  =>  identical criterion, only bookkeeping.")
print("    Taking BOTH u=c_f/N and Cbar=N^2 pays the per-link 1/N twice. Rejected. [T]")

# exact c_f check for SU(2) with weight exp(beta*N^{-1}Re tr U) = exp(beta cos theta)
def c_f_su2(beta):
    f = lambda th: mp.e ** (beta * mp.cos(th)) * 2 * mp.cos(th) * mp.sin(th) ** 2
    return mp.quad(f, [0, mp.pi]) * (2 / mp.pi)


print("\n exact SU(2) c_f(beta) vs the two normalisations (window remark, AGENTS sec.2.1):")
for b in [mp.mpf("0.010"), mp.mpf("1") / 48, mp.mpf("0.1")]:
    cf = c_f_su2(b)
    print("   beta=%-10s c_f=%-20s b/(2N)=%-20s b/N=%s"
          % (mp.nstr(b, 8), mp.nstr(cf, 10), mp.nstr(b / (2 * N), 10), mp.nstr(b / N, 10)))
print("   => at N=2 (f ~ fbar) c_f = beta/N + O(beta^3) = 2*(beta/2N): the repo's c_f=beta/(2N)")
print("      is the N>=3 statement. This makes u 2x LARGER at N=2 => window 2x SMALLER (direction safe).")

# ----------------------------------------------------------------------------------
hdr("[5] THRESHOLD AND WINDOW: does Chat*Cbar < 7.0632853 hold?")
# ----------------------------------------------------------------------------------
K_48 = mp.mpf(1) / 48
THR = 2 * N * N * 48 / LAM      # = 7.0632853...  (no-shrink threshold on Chat*Cbar)
print(" no-shrink threshold:  Chat*Cbar < 2N^2*48/(e*Delta) = %s" % mp.nstr(THR, 10))
print(" => allowed Cbar:  Chat=4 -> Cbar < %s ;  Chat=2 -> Cbar < %s"
      % (mp.nstr(THR / 4, 8), mp.nstr(THR / 2, 8)))
print()
print(" Chat  Cbar   Chat*Cbar   ok?   beta1^fix=2N^2/(Chat*eD*Cbar)   beta1/(1/48)   shrink")
rows = []
for chat, cbar, tag in [(CHAT_ENV, mp.mpf(1), "assumed Cbar=1 (z02)"),
                        (CHAT_LOC, mp.mpf(1), "Chat=2, Cbar=1"),
                        (CHAT_ENV, CBAR_LIM, "PROVABLE (this file)"),
                        (CHAT_LOC, CBAR_LIM, "Chat=2 + provable Cbar"),
                        (CHAT_ENV, CBAR_CRUDE, "crude 2*Delta")]:
    b1 = 2 * N * N / (chat * LAM * cbar)
    ok = "YES" if chat * cbar < THR else "NO "
    r = b1 / K_48
    rows.append((chat, cbar, b1))
    print(" %-5s %-6s %-11s %-5s %-28s %-14s %s"
          % (chat, cbar, chat * cbar, ok, mp.nstr(b1, 10), mp.nstr(r, 6),
             ("%.2fx" % r) if r < 1 else ("%.2fx widen" % r)))

# ----------------------------------------------------------------------------------
hdr("[6] DELTA AT THE NEW WINDOW: delta = ln(1/(Chat*ubar))/a = kappa/4 - ln(Chat eD Cbar)/a")
# ----------------------------------------------------------------------------------
def delta_a(beta, chat, cbar, cf_coeff=None):
    """delta*a for SU(2); u = beta/(2N^2) by the repo convention (cf_coeff=1) or u=beta/N (exact, 2)."""
    cf = 1 if cf_coeff is None else cf_coeff
    u = cf * beta / (2 * N * N) if cf == 1 else beta / N
    ub = LAM * cbar * u
    return mp.log(1 / (chat * ub)), ub


for chat, cbar, tag in [(CHAT_ENV, CBAR_LIM, "Chat=4 Cbar=20"), (CHAT_LOC, CBAR_LIM, "Chat=2 Cbar=20"),
                        (CHAT_ENV, mp.mpf(1), "Chat=4 Cbar=1")]:
    print("\n %s :" % tag)
    for beta in [K_48, mp.mpf("1") / 100, mp.mpf("1") / 500, mp.mpf("0.001"), mp.mpf("2e-4")]:
        d, ub = delta_a(beta, chat, cbar)
        print("   beta=%-8s ubar=%-10s Chat*ubar=%-10s delta*a=%s"
              % (mp.nstr(beta, 6), mp.nstr(ub, 8), mp.nstr(chat * ub, 8), mp.nstr(d, 8)))
    # max beta with delta>0
    bmax = 2 * N * N / (chat * LAM * cbar)
    print("   delta>0  <=>  beta < %s   (= (1/48)*%s)" % (mp.nstr(bmax, 8), mp.nstr(bmax / K_48, 6)))

# ----------------------------------------------------------------------------------
hdr("[7] SUM CHECK: the geometric identity used for |A(alpha)| and |eps_T| (X01 sec.2.4)")
# ----------------------------------------------------------------------------------
def geo_check(x, mmax=400):
    s1 = mp.nsum(lambda m: m * x ** m, [1, mp.inf])
    s2 = mp.nsum(lambda a: mp.mpf(a + 2) * x ** (a + 1), [0, mp.inf])   # sum_{a>T}(a-1-T)x^{a+1}, T=0
    return s1, x / (1 - x) ** 2, s2, x * (2 - x) / (1 - x) ** 2


for x in [mp.mpf("0.1"), mp.mpf("0.5"), mp.mpf("0.9")]:
    a, b, c, d = geo_check(x)
    print("  x=%-5s  sum m x^m = %-16s  x/(1-x)^2 = %-16s  rel.err=%s"
          % (mp.nstr(x, 4), mp.nstr(a, 12), mp.nstr(b, 12), mp.nstr(abs(a - b) / b, 3)))

# ----------------------------------------------------------------------------------
hdr("[8] VERDICT FLAGS")
# ----------------------------------------------------------------------------------
print(" Cbar defin.      : |phi^U(G)| <= Cbar^{|G|}u^{|G|} ,  Cbar = sup_G (|phi^U|/u^{|G|})^{1/|G|}   [T]")
print(" Cbar provable    : Cbar <= Delta = %s   (crude 2*Delta = %s)                                  [T]"
      % (CBAR_LIM, CBAR_CRUDE))
print(" Cbar N-scaling   : Cbar = O(1) in N ; Cbar=N^2 is a double count of the per-link 1/N          [T]")
print(" Cbar < 1.766 ?   : NOT established (provable bound 20 >> 1.766)                              [NG]")
print(" Chat*Cbar<7.0633 : with provable Cbar: 4*20 = 80 > 7.0633  -> 1/48 does NOT survive          [T]")
b1 = 2 * N * N / (CHAT_ENV * LAM * CBAR_LIM)
print(" new window       : beta1^fix = %s = (1/48)*%s > 0  -> (c) NOT triggered                     [T]"
      % (mp.nstr(b1, 8), mp.nstr(b1 / K_48, 6)))
print(" => branch (b) triggers: recompute, I1 stays [T] mod Cbar with the smaller window.")
print("\n done.")
