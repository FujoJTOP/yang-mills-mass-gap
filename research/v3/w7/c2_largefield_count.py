#!/usr/bin/env python3
"""
c2_largefield_count.py -- illustrative counting estimate for the (a,Lambda)-uniform
Balaban-type RG step (attack surface C, second angle; the L3 node).

WHAT THIS IS
------------
A *counting* witness, not a proof. It computes the entropy of large-field blocks in the
UV tower under the crude (block-by-block union) bound, and checks whether the count is
*per-unit-volume* (volume-independent) or *total* (grows like Lambda^d).

The point it illustrates (Lemma B in C2-BALABAN-UNIFORM.md):
  the crude large-field entropy channel is VOLUME-INDEPENDENT per unit volume; whether it
  converges is a CONSTANT inequality (large-field suppression vs block entropy), NOT a
  statement about the number of blocks. Hence "the volume growth is intrinsic" cannot be
  sourced from large-field entropy.

CONVENTIONS (campaign normalisation, fixed by the campaign files; NOT rederived here)
  - d = 4, blocking factor L (integer), torus volume Lambda^d (Lambda = side, lattice units).
  - shell index j: larger j = deeper UV, block spacing a_j = a_0 L^{-j}.
  - one-loop: g_j^2 = 1/(2 beta0 j), beta_j = 2N/g_j^2 = 4 N beta0 j   (TOWER-POINCARE).
  - beta0 = 11 N / (48 pi^2)                                          (SU(N) one-loop).
  - large-field weight at scale j (field above a fixed threshold t): p_j = exp(-t^2/g_j^2)
      = exp(-2 beta0 t^2 j).  Put c := t^2 (the suppression coefficient).
  - number of blocks at scale j in volume Lambda^d: (Lambda/a_j)^d = (Lambda/a_0)^d L^{d j}.
  => per-scale density (per unit physical volume, in units a_0^{-d}) at scale j:
        n_j(c) = L^{d j} exp(-2 c beta0 j) = exp( j (d ln L - 2 c beta0) ).
  => per-unit-volume total: rho(c) = sum_{j>=1} n_j(c), finite iff  2 c beta0 > d ln L.
  => threshold c* = d ln L / (2 beta0).
This is a standard geometric/probability estimate; nothing here is a theorem about YM.
"""

import math

def beta0(N):
    return 11.0 * N / (48.0 * math.pi**2)

def critical_c(d, L, N):
    # 2 c beta0 = d ln L  =>  c* = d ln L / (2 beta0)
    return d * math.log(L) / (2.0 * beta0(N))

def rho(c, d, L, N, Jmax):
    # per-unit-volume large-field block density, crude union bound, truncated at Jmax shells
    b = d * math.log(L) - 2.0 * c * beta0(N)
    # geometric series sum_{j=1}^{Jmax} exp(b j), guarded
    if abs(b) < 1e-12:
        return float(Jmax) * math.exp(b)
    return math.exp(b) * (math.exp(b * Jmax) - 1.0) / (math.exp(b) - 1.0)

def main():
    d = 4
    print("=" * 78)
    print("Large-field block counting in the UV tower (crude union bound)")
    print("=" * 78)
    for L in (2, 3):
        for N in (2, 3):
            b0 = beta0(N)
            cstar = critical_c(d, L, N)
            print(f"\nL={L}  N={N}  beta0={b0:.5f}   d lnL={d*math.log(L):.4f}"
                  f"   2 beta0={2*b0:.5f}   c* = t*^2 = {cstar:.2f}  (t* = {math.sqrt(cstar):.2f})")

    print("\n" + "-" * 78)
    print("Per-unit-volume density rho(c)  [Jmax = 60 shells, ~ a = L^-60]")
    print("Volume-INDEPENDENT by construction (no number of spatial blocks appears).")
    print("-" * 78)
    d = 4
    L = 2
    for N in (2, 3):
        cstar = critical_c(d, L, N)
        print(f"\nN={N}, L={L}, c*={cstar:.2f}")
        print(f"{'c':>6} {'rho(c)':>16} {'converges?':>12}")
        for c in (1.0, 5.0, 10.0, 20.0, math.ceil(cstar) + 1.0, 40.0, 80.0):
            r = rho(c, d, L, N, 60)
            conv = "YES" if (2*c*beta0(N) > d*math.log(L)) else "NO (diverges)"
            print(f"{c:6.1f} {r:16.4g} {conv:>12}")

    print("\n" + "-" * 78)
    print("Total expected large-field blocks = rho * (number of spatial coarse cells)")
    print("This DOES grow like (Lambda/a_0)^d -- but that is extensivity, not a wall.")
    print("-" * 78)
    N = 2; L = 2; c = 40.0  # comfortably above c* ~ 30 for the crude bound
    r = rho(c, d, L, N, 60)
    print(f"\nN={N}, L={L}, c={c}: per-unit-volume rho = {r:.5g} (FIXED)")
    print(f"{'Lambda/a_0':>12} {'# coarse cells':>16} {'total LF blocks':>18}")
    for M in (4, 8, 16, 32, 64):
        cells = M**d
        print(f"{M:12d} {cells:16d} {r*cells:18.4g}")

    print("\n" + "=" * 78)
    print("READING:")
    print(" * per-unit-volume rho is independent of the number of blocks; its finiteness is a")
    print("   CONSTANT inequality 2 c beta0 > d ln L  (c > c*), not a volume-growth statement.")
    print(" * total LF blocks grow like Lambda^d: that is the extensivity of a density, not an")
    print("   obstruction -- the density, not the total, is what an action-density bound needs.")
    print(" * NO INTRINSIC VOLUME GROWTH in this channel: 'intrinsic volume growth' is not")
    print("   sourced by large-field entropy.  The residual is the analytic per-shell (star).")
    print(" * Caveat: the crude union bound needs a LARGE suppression constant c* ~ 20-30;")
    print("   a Kotecky-Preiss connected-polymer bound lowers this threshold (counts per")
    print("   polymer, not per block).  Either way the criterion is const-level, volume-free.")
    print("=" * 78)

if __name__ == "__main__":
    main()
