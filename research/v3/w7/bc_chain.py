#!/usr/bin/env python3
"""bc_chain.py -- BC-CHAIN constant-chain arithmetic (N3 route, D:\\Dev\\math).

WHAT THIS IS
------------
Pure arithmetic on the constants that ALREADY APPEAR in
  research/v3/C2-BALABAN-UNIFORM.md   (Lemma B, c*, t*, rho(c), one-loop flow)
  docs/paper-v2/ch13-machine.tex      (lambda_tH in [1,7.3]; "Balaban reaches lambda~1")
  preprint/CPB-preprint.md  §3.3      ("with the coupling as expansion parameter; stops at lambda~1")
  research/v3/.../W1-STRONG-ALIGN.md  (u = c_f/N; K <= 80 e^{2+beta} N, P3-FULLREP-KP)

NOTHING here is a theorem about YM.  It is:
  (i)  re-derivation / verification of the published-in-repo constants,
  (ii) the unit conversion between the THREE "lambda"s that the corpus calls "lambda",
  (iii) the multiplicative slack M of the constant chain, link by link,
  (iv) the explicit value of the small-field expansion parameter at lambda_tH = 1.

Run:  py -3.12 research/v3/w7/bc_chain.py
"""
import math

# ---------------------------------------------------------------- conventions
# Wilson lattice action: S = beta sum_p (1 - (1/N) Re tr U_p),  beta = 2N/g^2.
# 't Hooft coordinate used by ch13 / CTL-CLOSEST:  lambda_tH := beta * N = 2 N^2 / g^2.
# one-loop (campaign normalisation, G1-L3-TOWER-POINCARE): beta0 = 11 N / (48 pi^2).

def beta0(N):
    return 11.0 * N / (48.0 * math.pi ** 2)

def lam_tH_from_g2(N, g2):
    return 2.0 * N * N / g2

def g2_from_lam_tH(N, lam):
    return 2.0 * N * N / lam

def beta_from_g2(N, g2):
    return 2.0 * N / g2

# ---------------------------------------------------------------- the chain
def c_star_union(d, L, N):
    """Lemma B (C2 §2.2, crude block-by-block union bound): 2 c beta0 > d ln L."""
    return d * math.log(L) / (2.0 * beta0(N))

def c_star_KP(d, N):
    """Kotecky-Preiss refinement: count per polymer, entropy per block <= ln(2d).

    A fixed connected polymer of n blocks at one scale has at most (2d)^n
    realisations (each new block attaches to an existing one; 2d = 8 in d=4),
    and carries suppression exp(-c n / g_j^2).  The per-block KP sum converges
    iff  ln(2d) < c / g_j^2 , worst at the head shell (j=1, g_1^2 = 1/(2 beta0)).
    =>  c** = ln(2d) / (2 beta0).
    """
    return math.log(2 * d) / (2.0 * beta0(N))

def q_ratio(c, d, L, N):
    """Geometric ratio of the per-unit-volume large-field series at suppression c."""
    return math.exp(d * math.log(L) - 2.0 * c * beta0(N))

def rho(c, d, L, N):
    """Per-unit-volume large-field block density (C2 §2.2 / witness c2_largefield_count.py)."""
    q = q_ratio(c, d, L, N)
    if q >= 1.0:
        return float('inf')
    return a0_floor * q / (1.0 - q)          # a0^-d = 1 in lattice units a0 = 1

a0_floor = 1.0

# ---------------------------------------------------------------- T3 quantity
def eta_B(N, lam_tH):
    """Effective small-field expansion parameter at a shell whose coupling is g^2 = 2N^2/lam_tH.

    Why this product:
      * Lemma B forces the small-field radius to satisfy t_j > t* (else the large-field
        entropy does not converge), i.e. the 'small' field is |F| <= t* g_j;
      * the per-shell expansion is in (coupling) x (amplitude in coupling units)
        ~ g_j * t* .
    t* is the radius in units of g_j (C2 §2.3 table: c* = t*^2).
    """
    t_star = math.sqrt(c_star_union(4, 2, N))
    g = math.sqrt(g2_from_lam_tH(N, lam_tH))
    return t_star * g, t_star, g

def lam_at_eta_one(N):
    """lambda_tH at which eta_B = 1 (i.e. g^2 = 1/t*^2)."""
    cs = c_star_union(4, 2, N)
    g2 = 1.0 / cs
    return lam_tH_from_g2(N, g2)

def main():
    d, L = 4, 2
    print("=" * 78)
    print("P0  VERIFY the in-repo constants (C2-BALABAN-UNIFORM.md §2.2/§2.3)")
    print("=" * 78)
    for N in (2, 3):
        cs = c_star_union(d, L, N)
        print(f"  N={N}: beta0={beta0(N):.5f}  c*=d lnL/(2 beta0)={cs:.2f}  t*=sqrt(c*)={math.sqrt(cs):.2f}")
    for N in (2, 3):
        r = rho(40.0, d, L, N)
        print(f"  N={N}, c=40: q={q_ratio(40.0,d,L,N):.4f}  rho(40)={r:.4f}   "
              f"(doc: 0.638 for N=2)")
    print("  L=3 (why L=2 is forced):")
    for N in (2, 3):
        print(f"    N={N}: c*(L=3)={c_star_union(d,3,N):.2f}  vs c*(L=2)={c_star_union(d,2,N):.2f}"
              f"   ratio={c_star_union(d,3,N)/c_star_union(d,2,N):.3f}")

    print()
    print("=" * 78)
    print("P1  LINK 2 (polymer): the ONLY large-field slack, union bound -> Kotecky-Preiss")
    print("=" * 78)
    for N in (2, 3):
        cu, ck = c_star_union(d, L, N), c_star_KP(d, N)
        print(f"  N={N}: c*_union={cu:.2f} (t*={math.sqrt(cu):.2f})   "
              f"c**_KP={ck:.2f} (t**={math.sqrt(ck):.2f})")
        print(f"        improvement factor c*/c** = {cu/ck:.4f}  "
              f"(= d lnL/ln(2d) at L=2)")
        print(f"        in t*: {math.sqrt(cu/ck):.4f}   in admissible g^2: {cu/ck:.4f}")
    print("        NOTE: c**_KP = ln(2d)/(2 beta0) is L-INDEPENDENT; c*_union = d lnL/(2 beta0)")
    print("              grows with L  =>  L=2 is the optimal integer blocking factor:")
    for Lx in (2, 3, 4):
        print(f"              L={Lx}: c*_union={d*math.log(Lx)/(2*beta0(3)):6.2f}   "
              f"c*/c**_KP = {d*math.log(Lx)/math.log(2*d):5.3f}")

    print()
    print("=" * 78)
    print("P2  UNIT CONVERSION: the corpus has THREE objects called 'lambda'")
    print("=" * 78)
    print("   (R1) lambda = g^2            [CPB-preprint §3.3: 'with the coupling as expansion parameter']")
    print("   (R2) lambda = lambda_tH=betaN [ch13-machine.tex: crossover at lambda_tH in [1,7.3]]")
    print("   (R3) lambda = beta            [lattice]")
    print()
    print(f"   {'point':>28} | {'beta':>7} | {'g^2':>8} | {'lambda_tH':>10}")
    for N in (2, 3):
        print(f"   --- N={N} ---")
        rows = [
            ("SZZ radius  beta=1/48", 1.0 / 48.0),
            ("crossover LOW edge lam_tH=1", 1.0 / N),
            ("physical point beta=2.4", 2.4),
            ("physical point beta=2.5", 2.5),
            ("lambda_tH = 7.3 (ch13)", 7.3 / N),
            ("g^2 = 1  (Balaban marginal)", 2.0 * N),
        ]
        for name, b in rows:
            g2 = 2.0 * N / b
            print(f"   {name:>28} | {b:7.4f} | {g2:8.3f} | {b*N:10.3f}")
    print()
    print("   => 'lambda ~ 1' (R1, g^2 ~ 1) is lambda_tH = 2N^2 = "
          f"{2*2*2} (N=2), {2*3*3} (N=3).")
    print("   => 'lambda_tH = 1' is g^2 = 2N^2 = "
          f"{2*2*2} (N=2), {2*3*3} (N=3).")
    print("      The two '1's differ by exactly 2N^2.  They are NOT the same number.")

    print()
    print("=" * 78)
    print("P3  LINK 6 (the lambda threshold) + T3: the small-field expansion parameter")
    print("=" * 78)
    print("   eta_B := t* * g  with g^2 = 2N^2/lambda_tH   (t* from Lemma B, radius in units of g)")
    print(f"   {'lambda_tH':>10} | {'g^2':>9} {'g':>7} | {'eta_B(N=2)':>11} {'eta_B(N=3)':>11}")
    for lam in (1.0, 2.0, 7.3, 18.0, 50.0, 100.0, 358.0):
        e2, _, _ = eta_B(2, lam)
        e3, _, _ = eta_B(3, lam)
        g2 = g2_from_lam_tH(3, lam)
        print(f"   {lam:10.1f} | {g2:9.4f} {math.sqrt(g2):7.3f} | {e2:11.3f} {e3:11.3f}")
    print()
    for N in (2, 3):
        e, t, g = eta_B(N, 1.0)
        print(f"   N={N}: at lambda_tH=1 -> g^2={g*g:.3f}, g={g:.3f}, t*={t:.3f}, eta_B={e:.2f}")
    print()
    print("   eta_B = 1  (the formal 'O(1)' line) is reached only at")
    for N in (2, 3):
        print(f"     N={N}: lambda_tH = 2 N^2 c* = {lam_at_eta_one(N):.1f}   "
              f"(g^2 = 1/c* = {1.0/c_star_union(d,L,N):.4f})")

    print()
    print("=" * 78)
    print("P4  head of the one-loop tower (no new physics; the repo's own flow)")
    print("=" * 78)
    print("   g_j^2 = 1/(2 beta0 j),  beta_j = 4 N beta0 j ,  lambda_tH(j) = beta_j N")
    for N in (2, 3):
        b0 = beta0(N)
        print(f"   --- N={N} ---")
        for j in (1, 2, 5, 10, 50):
            g2 = 1.0 / (2.0 * b0 * j)
            print(f"     j={j:3d}: g_j^2={g2:8.4f}  beta_j={4*N*b0*j:8.3f}  "
                  f"lambda_tH(j)={4*N*b0*j*N:9.3f}")
        j1 = 1.0 / (2.0 * b0)       # shell where g^2 = 1  (tower depth to the head)
        print(f"     shells to reach g^2=1 : j = {j1:.2f}   (K = log_2(1/a) = {j1:.2f})")
        print(f"     head-shell (j=1) lambda_tH = {4*N*b0*N:.3f}  -> "
              f"{'INSIDE' if 1 <= 4*N*b0*N <= 7.3 else 'outside'} [1,7.3]")

    print()
    print("=" * 78)
    print("P5  T2: the multiplicative slack M of the constant chain")
    print("=" * 78)
    print("   per-link factor = how much the ADMISSIBLE g^2 could rise if that link")
    print("   were replaced by the finest estimate known (or computable) today.")
    print()
    fac = {
        1: ("small-field decomposition  (radius t_j = kappa g_j)", 1.00,
            "kappa is *set by* link 2 (Lemma B); no independent slack"),
        2: ("polymer / cluster estimate (union bound -> KP)", None,
            "c*/c** = d lnL/ln(2d) = 1.333 (exact, L-independent)"),
        3: ("K-quantity  (polymer counting constant)", 10.0,
            "crude K=80 e^{2+beta} N (P3-FULLREP-KP); 'the true K is smaller' -- NOT quantified -> <=10"),
        4: ("gauge fixing", 1.00,
            "tower Poincare constant c_L = 4 pi^2/L^2 is the exact torus spectral gap; sharp"),
        5: ("uv-stability exponent", 1.00,
            "superrenormalizable type = fixed by 4D power counting; dimension counting is exact"),
        6: ("the lambda threshold itself", 1.00,
            "NOT a bound: it is the hypothesis 'g in (0,gamma], gamma small'; gamma unpublished"),
    }
    fac[2] = (fac[2][0], c_star_union(d, L, 3) / c_star_KP(d, 3), fac[2][2])
    lo, hi = 1.0, 1.0
    for k in sorted(fac):
        name, f, why = fac[k]
        print(f"   L{k}: {name:52s} factor {f:6.3f}")
        print(f"       {why}")
        if k == 3:
            lo *= 1.0; hi *= f
        else:
            lo *= f; hi *= f
    print()
    print(f"   M in [{lo:.2f}, {hi:.2f}]   (geometric mid ~ {math.sqrt(lo*hi):.2f}; order 10^0.5)")
    print()
    print("   TARGETS (factor by which the admissible g^2 must rise):")
    for N in (2, 3):
        g2_marg = 1.0                       # Balaban's marginal g^2 (R1: lambda~1)
        tgt_phys = 2.0 * N * N / 7.3        # physical point in lambda_tH=7.3
        tgt_edge = 2.0 * N * N / 1.0        # crossover lower edge lambda_tH=1
        print(f"     N={N}: to physical point (lam_tH=7.3): x{tgt_phys/g2_marg:6.2f}"
              f"   |  to crossover lower edge (lam_tH=1): x{tgt_edge/g2_marg:6.2f}")
    print()
    print("   VERDICT (arithmetic only): M_max ~ 10^1 is the same ORDER as the physical-point")
    print("   target 2.5-3 but BELOW the crossover-edge target 8-18; and the deciding link")
    print("   (L6) contributes exactly 1.00 because it is a parameter, not a bound.")

    print()
    print("=" * 78)
    print("P6  the counterfactual: what M would have to be, in each coordinate")
    print("=" * 78)
    for N in (2, 3):
        print(f"   N={N}:  lambda_tH=1 needs g^2 rise x{2*N*N:.0f};  "
              f"equivalently a SMALL PARAMETER at g^2={2*N*N} does not exist unless")
        print(f"         a hidden large parameter is found (N, q=m a, or 1/L).")

if __name__ == "__main__":
    main()
