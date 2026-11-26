# d02_dimcount.py -- D02-DIRECT: exact dimension accounting for the SU(2) transfer matrix at L=2
# Distinguishes THREE different objects that the corpus conflates:
#   (A) 4D full torus config space (r04's 143)          -- NOT the transfer matrix
#   (B) one time slice, spatial links only (standard Wilson T)
#   (C) one time slice including temporal links (the "4*L^3 = 32 link" object)
# Plus the matrix-element integral dimension and the truncation cost.
# Only integers printed; everything is a counting identity.
import sys, math
from fractions import Fraction

def main():
    print("="*78)
    print("D02 dimension accounting, SU(2), dim(SU(2))=3 (real)")
    print("="*78)
    print("\n[A] 4D FULL TORUS  (r04's object: the S-D / Hessian configuration space)")
    print(f"{'L':>3} {'sites=L^4':>10} {'links=4L^4':>11} {'comps=3*links':>14} {'plaq=6L^4':>10} {'gauge=3*sites':>14} {'phys=comps-gauge':>17} {'cond=phys-1':>11}")
    A = {}
    for L in (1,2,3,4,6,8):
        sites = L**4; links = 4*sites; comps = 3*links; plaq = 6*sites
        gauge = 3*sites; phys = comps - gauge
        A[L] = (sites, links, comps, plaq, gauge, phys, phys-1)
        print(f"{L:>3} {sites:>10} {links:>11} {comps:>14} {plaq:>10} {gauge:>14} {phys:>17} {phys-1:>11}")
    print("  r04 check: L=2 -> 192 raw, 144 gauge-fixed, 143 after codim-1 conditioning:",
          A[2][2] == 192, A[2][5] == 144, A[2][6] == 143)

    print("\n[B] ONE TIME SLICE, SPATIAL LINKS ONLY  (carrier of the Wilson transfer matrix T)")
    print(f"{'L':>3} {'spatial links=3L^3':>18} {'manifold dim=9L^3':>18} {'gauge params=3L^3':>18} {'orbit dim':>10} {'physical':>9}")
    B = {}
    for L in (2,3,4,6,8):
        links = 3*L**3; mdim = 3*links; gp = 3*L**3
        orbit = min(mdim, gp)      # generic stabiliser trivial when gp <= mdim
        phys = mdim - orbit
        B[L] = (links, mdim, gp, orbit, phys)
        print(f"{L:>3} {links:>18} {mdim:>18} {gp:>18} {orbit:>10} {phys:>9}")
    print("  L=2 -> 24 links, 72-dim manifold, 24-dim gauge orbit, 48 physical")

    print("\n[C] ONE TIME SLICE INCLUDING TEMPORAL LINKS  (the '4*L^3 = 32 link' reading)")
    print(f"{'L':>3} {'links=4L^3':>11} {'manifold=12L^3':>15} {'gauge=3L^3':>11} {'physical':>9}")
    C = {}
    for L in (2,3,4,6,8):
        links = 4*L**3; mdim = 3*links; gp = 3*L**3
        C[L] = (links, mdim, gp, mdim-gp)
        print(f"{L:>3} {links:>11} {mdim:>15} {gp:>11} {mdim-gp:>9}")
    print("  L=2 -> 32 links, 96-dim manifold, 72 physical   (user's 4*L^3=32 confirmed)")

    print("\n[D] MATRIX-ELEMENT INTEGRAL DIMENSION (the hard sub-object of T)")
    print("    T[U,U'] = e^{-S_sp(U)/2} e^{-S_sp(U')/2} W(U,U'),")
    print("    W(U,U') = int prod_x dV(x) exp((beta/2) sum_{x,i} Re tr[U_i(x) V(x+i) U'_i(x)^dag V(x)^dag])")
    print("    -> one integral over L^3 temporal links = 3 L^3 real dimensions.")
    print(f"{'L':>3} {'integral dim=3L^3':>18} {'bits at eps=1e-6':>18} {'ratio to 143':>13}")
    for L in (2,4,6,8):
        d = 3*L**3
        print(f"{L:>3} {d:>18} {d*math.log2(1e6):>18.1f} {143.0/d:>13.2f}")
    print("  r04's feasibility threshold was '<=20 dims'; the transfer route needs 24 at L=2")

    print("\n[E] TRUNCATION COST: full SU(2) character (Peter-Weyl) basis per link, j<=J")
    print("    basis size b_J = sum_{m=1..2J+1} m^2  (half-integer j included: m=1,2,3,...)")
    # b_J = sum_{m=1}^{2J+1} m^2 = (2J+1)(2J+2)(4J+3)/6
    for J in (2,3,4):
        M = 2*J+1
        b = sum(m*m for m in range(1, M+1))
        print(f"  J={J}: M={M}, basis per link b={b}", end="")
        for L in (2,4):
            n = 3*L**3
            print(f" | L={L}: states=b^{n}=1e{ n*math.log10(b):.1f}", end="")
        print()
    print("  (J=3 is not a free choice: see d02_cert.log tail bound)")

    print("\n[F] SANITY: the three counts are distinct objects at L=2")
    print("    full 4D torus raw comps      :", A[2][2])
    print("    slice spatial-only manifold  :", B[2][1])
    print("    slice incl. temporal manifold:", C[2][1])
    print("    r04 conditional dim          :", A[2][6])
    print("    -> 96 != 192 ; 143 is not a transfer-matrix dimension (it is 4*L^4*3 - 3*L^4 - 1)")
    return 0

if __name__ == "__main__":
    sys.exit(main())
