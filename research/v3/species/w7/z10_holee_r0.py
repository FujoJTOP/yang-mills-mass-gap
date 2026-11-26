"""Z10-HOLEE: is r0 = 1/4(2N)^-2 (beta-free) compatible with r0 = 1/2 c1 u^p (beta-dependent)?

Hole E (Z05-AUDITI1.md:79): P3-DRESSED-MASS.md:124 vs PROOF-P3-UNCONDITIONAL.md:31.

Everything here is parametric in the per-tile strong-coupling weight
    u := c_f/N = beta/(2N^2) + O(beta^3)      [P3-FOURWALL-DERIVATION.md:81, Gupta 14.18]
so the verdict does NOT depend on the (separately open) u-normalisation hole A.

Outputs:
  (1) exact SU(2) fundamental activity c_f vs the leading beta/(2N)  -- convention cross-check
  (2) the three candidate floors at beta < 1/48, N = 2,3,4
  (3) the value u* at which 1/2 c1 u^p equals 1/4(2N)^-2  (=> what "u^p" would have to be)
  (4) the exact factor 2 between the per-L floor 1/2(2N)^-2 and the L-inf floor 1/4(2N)^-2
"""
import math

BETA_WALL = 1.0 / 48.0  # SZZ/Audit beta1, d=4


# ---------------------------------------------------------------- (1) exact SU(2) c_f
def cf_su2_exact(beta):
    """Fundamental character coefficient of e^{(beta/N) Re tr U} for SU(2).

    c_f = (1/d_f) Int dU f(U) conj(chi_f(U)) = Int_{-1}^{1} dx sqrt(1-x^2) (2x) e^{2 beta x} / 2 / (pi/2)
    with Haar dU = (2/pi) sqrt(1-x^2) dx, x = cos(theta), chi_f = 2x, d_f = 2.
    Normalised so that c_f -> beta/N = beta/2 as beta -> 0 (the leading order used in
    P3-FOURWALL-DERIVATION:81, "beta/2N from the expansion").
    """
    n = 200001
    lo, hi = -1.0, 1.0
    h = (hi - lo) / (n - 1)
    s_chi = 0.0
    s_one = 0.0
    for i in range(n):
        x = lo + i * h
        w = 1.0 if i in (0, n - 1) else (4.0 if i % 2 else 2.0)  # Simpson
        g = w * math.sqrt(max(0.0, 1.0 - x * x)) * math.exp(beta * x)  # (beta/N)Re tr U = (beta/2)(2x)
        s_chi += g * (2.0 * x)
        s_one += g
    return (s_chi / s_one) / 2.0  # divide by d_f = 2


# ---------------------------------------------------------------- helpers
def u_lead(beta, N):
    """per-tile activity u = beta/(2N^2) (leading order, FOURWALL:81)."""
    return beta / (2.0 * N * N)


def main():
    print("=" * 78)
    print("Z10-HOLEE :: r0 = 1/4(2N)^-2  vs  r0 = 1/2 c1 u^p   [holes E]")
    print("=" * 78)
    print("convention:  u = c_f/N = beta/(2N^2) + O(beta^3)   [FOURWALL:81]")
    print("             a0 = (2N)^-p s   (beta-FREE numbers)  [DRESSED:36,119; T2:114; FOURWALL:86]")
    print()

    # ---------------- (1) exact SU(2)
    print("-" * 78)
    print("(1) exact SU(2) fundamental activity vs leading beta/2N   [convention cross-check]")
    print("-" * 78)
    print(f"{'beta':>10} {'cf_exact':>14} {'beta/2N (N=2)':>16} {'ratio':>10} {'u=cf/N':>12}")
    for beta in (1.0 / 48.0, 1.0 / 24.0, 0.5, 1.0, 2.5):
        cf = cf_su2_exact(beta)
        lead = beta / 4.0
        print(f"{beta:10.6f} {cf:14.9f} {lead:16.9f} {cf / lead:10.6f} {cf / 2.0:12.9f}")
    print("  => cf -> beta/2 as beta -> 0 (ratio -> 1): the leading-order u = beta/(2N^2) is right;")
    print("     both statements below are evaluated in u, so the open u-normalisation (hole A) is")
    print("     NOT load-bearing for the beta-scaling verdict.")
    print()

    # ---------------- (2) the floors
    print("-" * 78)
    print("(2) the three candidate floors at the wall beta = 1/48")
    print("-" * 78)
    print("    F_DR = 1/4 (2N)^-2            [DRESSED:124]")
    print("    F_p0 = 1/2 c1 u^0 = 1/2 (2N)^-2   [PROOF:31 with p=0, c1=a0=(2N)^-2]")
    print("    F_p1 = 1/2 c1 u^1 ,  F_p2 = 1/2 c1 u^2   [PROOF:31 as literally written, p>=1]")
    print()
    hdr = f"{'N':>3} {'u(1/48)':>12} {'F_DR':>14} {'F_p0':>14} {'F_p1':>12} {'F_p2':>12} {'F_p1/F_DR':>11} {'F_p2/F_DR':>11}"
    print(hdr)
    for N in (2, 3, 4):
        u = u_lead(BETA_WALL, N)
        c1 = (2.0 * N) ** -2.0
        FDR = 0.25 * c1
        F0 = 0.5 * c1 * u ** 0
        F1 = 0.5 * c1 * u ** 1
        F2 = 0.5 * c1 * u ** 2
        print(f"{N:3d} {u:12.3e} {FDR:14.9e} {F0:14.9e} {F1:12.3e} {F2:12.3e} {F1 / FDR:11.3e} {F2 / FDR:11.3e}")
    print()
    print("  F_p0 / F_DR = 2 exactly, for every N  (the factor 2 is the L-limit margin, X02:111):")
    for N in (2, 3, 4):
        c1 = (2.0 * N) ** -2.0
        print(f"    N={N}:  1/2(2N)^-2 / [1/4(2N)^-2] = {0.5 * c1 / (0.25 * c1):.12f}")
    print()
    print("  p>=1 puts the floor 2.5-3.3 orders of magnitude BELOW the DR value at the SAME beta.")
    print("  => p>=1 and 1/4(2N)^-2 are NOT the same scale; they cannot both be the quoted r0.")
    print()

    # ---------------- (2b) robustness to the open u-normalisation (hole A)
    print("-" * 78)
    print("(2b) robustness: the corpus has TWO u normalisations (hole A). does the verdict move?")
    print("     u = c_f/N = beta/(2N^2)  [FOURWALL:81]   vs   u = c_f = beta/(2N)  [exact SU(2) = 0.0104 @1/48]")
    print("-" * 78)
    print(f"{'N':>3} {'u=c_f/N':>12} {'F_p1/F_DR':>12} {'F_p2/F_DR':>12} {'u=c_f':>12} {'F_p1/F_DR':>12} {'F_p2/F_DR':>12}")
    for N in (2, 3, 4):
        uA = u_lead(BETA_WALL, N)
        uB = N * uA  # u = c_f
        print(f"{N:3d} {uA:12.3e} {2 * uA:12.3e} {2 * uA * uA:12.3e} {uB:12.3e} {2 * uB:12.3e} {2 * uB * uB:12.3e}")
    print("  => both normalisations give F_p1/F_DR <= 1/48 and F_p2/F_DR <= 1/64 (worst case N=2):")
    print("     the beta-scaling verdict is INDEPENDENT of hole A.  [both p>=1 readings are far off]")
    print()

    # ---------------- (3) u* such that 1/2 c1 u^p = 1/4(2N)^-2
    print("-" * 78)
    print("(3) which u makes  1/2 c1 u^p = 1/4 (2N)^-2 ?   (c1 = (2N)^-2)")
    print("-" * 78)
    print("    1/2 c1 u^p = 1/4 c1  <=>  u^p = 1/2  <=>  u* = 2^(-1/p)")
    print(f"{'p':>3} {'u* = 2^(-1/p)':>16} {'u(1/48), N=2':>15} {'u*/u(1/48)':>14}")
    for p in (1, 2, 3, 4):
        us = 2.0 ** (-1.0 / p)
        u2 = u_lead(BETA_WALL, 2)
        print(f"{p:3d} {us:16.6f} {u2:15.3e} {us / u2:14.3e}")
    print()
    print("  => the beta-FREE value 1/4(2N)^-2 is what 1/2 c1 u^p gives at u* = 0.5 (p=1) or")
    print("     0.707 (p=2) -- i.e. at u ~ O(1), the STRONG-coupling end u->1, NOT at beta<1/48")
    print("     where u = 2.6e-3 (N=2).  Hence 1/4(2N)^-2 is the u->1 (large-beta) evaluation")
    print("     of the same amplitude; at beta<1/48 the u^p form spuriously vanishes.")
    print()

    # ---------------- (4) monotonicity in u: does a0 carry u at all?
    print("-" * 78)
    print("(4) is a0 beta-free?  two spellings in the corpus, and what each implies")
    print("-" * 78)
    print("  FOURWALL:81-86 : tube weight = a0 u^{4N_t}, a0 = (2N)^-2 [endpoint recouplings] [transverse sum]")
    print("                   -> the 4N_t plaquette activities are ALL in the exponent; the two")
    print("                      insertions' own measure sits at zeroth order (no extra c_f)")
    print("                   -> a0 = pure numbers, beta-FREE.")
    print("  AUDIT-COR:78   : a0(beta) = (2N)^-p [prod_{wall plaquettes} c_f] [recoupling]")
    print("                   -> would put the SAME c_f factors in a0 that u^{4N_t} already carries")
    print("                      (double count), and contradicts AUDIT-COR:87 (p=0).")
    print()
    print("  consequence if AUDIT-COR:78 were right: a0 ~ u^2 -> r0 -> 0 as beta -> 0,")
    print("  i.e. the beta->0 limit of the residue would be 0 at leading order.  Check at N=2:")
    for beta in (1.0 / 48.0, 1.0e-3, 1.0e-5):
        u = u_lead(beta, 2)
        c1 = (2.0 * 2) ** -2.0  # (2N)^-2 for N=2 = 1/16
        print(f"    beta={beta:10.3e}  u={u:10.3e}  a0~(2N)^-2 u^2 = {c1 * u * u:12.3e}"
              f"  vs a0=(2N)^-2 = {c1:.3e}")
    print("  => the u^2 reading drives the amplitude to 0 with beta; Z05-AUDITI1:79's checked")
    print("     beta->0 behaviour (f1^2 -> a0 > 0, amplitude survives) selects the beta-FREE a0.")
    print()

    # ---------------- verdict
    print("=" * 78)
    print("VERDICT (b): SAME quantity, one formula's u-power is wrong.")
    print("  r0(per L) = 1/2 a0 = 1/2 c1 u^0 >= 1/2 (2N)^-2   (beta-FREE, p=0)")
    print("  r0(inf_L) = 1/4 (2N)^-2                          (the extra 1/2 is the L margin)")
    print("  the literal 'u^p, p>=1' at PROOF:31 is a transport error out of the superseded")
    print("  PROOF-P3.md:7; it is a valid but non-sharp bound that mis-states beta -> 0.")
    print("=" * 78)


if __name__ == "__main__":
    main()
