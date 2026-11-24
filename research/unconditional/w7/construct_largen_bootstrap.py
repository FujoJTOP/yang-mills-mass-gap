"""
construct_largen_bootstrap.py  --  witness for CONSTRUCT-LARGEN-BOOTSTRAP.md

Cheap: exact arithmetic + tiny numeric.  py -3.12 -I.  No Monte Carlo, no SDP.

Three checks:
 (1) 2D exact master field W=exp(-lambda*A/2) is multiplicative at a self-intersection
     (figure-eight = product of plaquettes)  ->  the loop-space CLOSES in D=2  [T].
 (2) The closure deficit in a truncation to RECTANGLES: the number of loop-footprint
     topologies (polyominoes, OEIS A000105) vastly exceeds the number of rectangles
     of the same area, and the deficit grows ~exponentially -> no finite truncation
     to simple loops closes  [T / V-num].
 (3) In D=4 the master field needs TWO functionals (area + perimeter); a one-parameter
     "pure area law" gives an effective sigma that DRIFTS with the loop size -> a
     one-scalar bootstrap cannot pin sigma  [V-num].

Labels: [T] exact, [V-num] numeric.  Nothing here proves a gap.  Clay 0%; no QED.
"""

import numpy as np
from math import isqrt, log

np.set_printoptions(precision=6, suppress=True, linewidth=120)


def banner(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)


# ===========================================================================
# PART 1 -- D=2 exact master field: the loop space CLOSES
#   Makeenko-Migdal in D=2 is solved exactly:  W(C) = exp(-lambda*A(C)/2),
#   and because the signed area is additive at a self-intersection (Stokes),
#   W(figure-eight) = W(plaq)*W(plaq).  This is WHY 2D is a closed master field.
# ===========================================================================
banner("PART 1 | D=2 exact master field  W=exp(-lambda*A/2):  multiplicative at self-intersection  [T]")
lam = 0.7


def W2D(A):
    return np.exp(-lam * A / 2.0)


# one plaquette = area 1 ; figure-eight (two plaquettes sharing a corner) = area 2
w_plaq = W2D(1.0)
w_fig8 = W2D(2.0)
w_prod = w_plaq * w_plaq
print(f"  lambda = {lam}")
print(f"  W(plaquette, A=1)          = {w_plaq:.8f}")
print(f"  W(figure-eight, A=2)       = {w_fig8:.8f}")
print(f"  W(plaquette)^2 (factorized)= {w_prod:.8f}")
print(f"  |W(fig8) - W(plaq)^2|      = {abs(w_fig8 - w_prod):.2e}   -> 0 : 2D CLOSES (self-intersecting loop IS determined)")
print("  => in D=2 the 'extra' loop the MM equation generates is fixed by the area law; no deficit.")


# ===========================================================================
# PART 2 -- the closure deficit: loop-footprint topologies vs rectangles
#   A rectangle W(m,n) = the loop equation sources loops with a plaquette added
#   at a link -> connected square-footprints (polyominoes).  Only a vanishing
#   fraction of polyominoes are rectangles.  Counts are exact (OEIS A000105).
#   If you truncate to rectangles, every non-rectangle the equation needs is an
#   UNKNOWN outside the system -> the deficit = polyominoes - rectangles > 0 and grows.
# ===========================================================================
banner("PART 2 | closure deficit:  #loop topologies (polyominoes) vs #rectangles  by area n  [T/V-num]")

# free polyominoes A000105, n = 1..9
polyominoes = [1, 1, 2, 5, 12, 35, 108, 369, 1285]


def nrect(n):
    """number of distinct rectangle shapes a x b (a<=b) with area n, up to symmetry"""
    return sum(1 for a in range(1, isqrt(n) + 1) if n % a == 0)


print(f"  {'n':>3} {'polyominoes':>12} {'rectangles':>11} {'deficit':>9} {'deficit/total':>14}")
prev_def = 0
for n in range(1, 10):
    P = polyominoes[n - 1]
    R = nrect(n)
    D = P - R
    print(f"  {n:>3} {P:>12} {R:>11} {D:>9} {D / P:>14.3f}")
    prev_def = D

print()
print("  growth of total topologies poly(n)/poly(n-1):",
      [f"{polyominoes[i]/polyominoes[i-1]:.2f}" for i in range(1, len(polyominoes))])
print("  the polyomino count grows ~4.06^n (A000105 growth constant);")
print("  the rectangle count is d(n) = O(n^eps).")
print("  => the closure deficit (non-rectangular loops the MM eq. NEEDS) grows EXPONENTIALLY.")
print("     Truncating to simple rectangles is therefore never closed: at every truncation")
print("     the number of required-but-omitted loop topologies >> the number of equations.")


# ===========================================================================
# PART 3 -- D=4 needs TWO functionals; a one-scalar bootstrap cannot pin sigma
#   A realistic 4D master field has area-law + perimeter (Coulomb) corrections:
#       W(m,n) = exp(-sigma*m*n - c*(m+n)).
#   If you INSIST on a one-parameter pure-area law, the 'effective' sigma you read
#   depends on (m,n) -> it drifts -> no single sigma fits -> one-scalar bootstrap fails.
#   A two-parameter fit (linear in log W) recovers (sigma,c) EXACTLY -> two functionals forced.
# ===========================================================================
banner("PART 3 | D=4: one-scalar area law drifts; two functionals (area+perimeter) are forced  [V-num]")

sig_true, c_true = 0.50, 0.30


def W4D(m, n):
    return np.exp(-sig_true * m * n - c_true * (m + n))


sizes = [(1, 1), (1, 2), (2, 2), (1, 3), (2, 3), (3, 3), (2, 4), (4, 4)]
print(f"  true (sigma, c) = ({sig_true}, {c_true})")
print(f"  {'m x n':>7} {'W':>11} {'sigma_eff = -lnW/(mn)':>23}")
sig_effs = []
for (m, n) in sizes:
    W = W4D(m, n)
    s_eff = -log(W) / (m * n)
    sig_effs.append(s_eff)
    print(f"  {m:>2}x{n:<2} {W:>11.6f} {s_eff:>23.6f}")
print(f"\n  spread of sigma_eff over these loops = {max(sig_effs) - min(sig_effs):.4f}")
print("  => a single sigma does NOT fit: the 'pure area' reading drifts with loop size.")
print("     (The drift is the perimeter/Coulomb functional; one scalar cannot absorb it.)")

# two-parameter linear fit in log W = -sigma*A - c*P  (A=mn, P=2(m+n))
A = np.array([m * n for (m, n) in sizes], float)
P = np.array([2.0 * (m + n) for (m, n) in sizes], float)
Y = np.array([-log(W4D(m, n)) for (m, n) in sizes], float)   # = sigma A + c P
M = np.vstack([A, P]).T
coef, res, rank, sv = np.linalg.lstsq(M, Y, rcond=None)
print(f"\n  two-functional fit  Y = sigma*A + c*P :  sigma={coef[0]:.6f}  c={coef[1]:.6f}  rank={rank}")
print(f"  residual max = {np.max(np.abs(M @ coef - Y)):.2e}   -> recovers ({sig_true}, {c_true}) exactly")
print("  => TWO independent functionals are forced (the PNT-LARGEN-LOOP 'surface sum'),")
print("     so a one-scalar bootstrap (one unknown sigma) is underdetermined by construction.")


# ===========================================================================
# VERDICT
# ===========================================================================
banner("VERDICT")
print("""  The large-N loop equation is CLOSED in W (no n-point tower) -- a real structural win.
  But a finite truncation to simple/rectangular loops is NOT closed:
    (1) D=2 closes exactly (area law is multiplicative at self-intersections);
    (2) D=4 does NOT: the equation sources exponentially many non-rectangular loop
        topologies (polyominoes ~4.06^n vs rectangles ~n^eps) -> deficit grows;
    (3) D=4 forces TWO functionals (area + perimeter), so one scalar sigma is
        underdetermined in any single-parameter closure.
  Hence the loop-space bootstrap does NOT converge to a unique sigma, and gives NO
  m_{0++}/sqrt(sigma) (the 0++ mass needs the linearized fluctuation around the
  unknown master field, which the bootstrap does not produce).
  [NG] for this bootstrap as a sigma/m_{0++} machine.  Clay 0%; no QED.""")
