#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""B03 — block-spin decimation at fixed b as a mass-gap route.

Question bank (task B03):
  Q1  one-step change of the effective coupling after blocking a -> b a ;
  Q2  is beta_eff = beta * b^x with x>0  (i.e. beta_eff > beta) ?
  Q3  iterate from beta=1/48 up to 2.5 : beta_k = ?, how many steps ?
  Q4  four-model test (free field / 2D YM / M6 / single plaquette).

Conventions (pinned, cf. research/v3/species/BETA-CONVENTION.md):
  Wilson   beta_W = 2N/g^2                       (dimensionless in d=4)
  't Hooft  beta_A = 2/(g^2 N) = beta_W / N^2    (SZZ; the wall 1/48 lives here)
  => beta_A = 1/48  <=>  beta_W = N^2/48  (SU(2): 1/12 = 0.08333, SU(3): 3/16)
  The target 2.5 is a *Wilson* value for SU(2)  (SU(3) crossover: 5.69).
  Mixing the two (2.5*48 = 120x) is the known P0 error; the like-for-like
  distance is ~30x.  This script never mixes them: every beta is tagged _W or _A.

One-loop input (standard, textbook):
  dg/dln mu = -b0 g^3,   b0 = 11N/(48 pi^2)      (pure YM, d=4)
  => d(1/g^2)/dln mu = 2 b0 = 11N/(24 pi^2)
  => d beta_W / dln mu = 2N * 2 b0 = 11 N^2/(12 pi^2) =: C_N   (0.3715 for N=2)

Coarse-graining is mu -> mu/b (Delta ln mu = -ln b), hence
  beta_W,eff = beta_W - C_N ln b ,     C_N = 11 N^2/(12 pi^2) > 0.
Nothing here is fitted; every number below is a closed form.
"""
import math
from fractions import Fraction

PI2 = math.pi ** 2


def C(N):
    """d beta_W / d ln(mu) at one loop."""
    return 11.0 * N * N / (12.0 * PI2)


def dbeta(N, b):
    """One-step change of beta_W under a -> b a (coarse-graining direction)."""
    return -C(N) * math.log(b)


def beta_eff(N, b, beta):
    return beta - C(N) * math.log(b)


def g2(N, beta_W):
    return 2.0 * N / beta_W


def x_exponent(N, b, beta):
    """If one insists on beta_eff = beta * b^x, what is x ?  (x<0 <=> beta_eff<beta)"""
    r = beta_eff(N, b, beta) / beta
    return math.log(r) / math.log(b) if r > 0 else float("nan")


WALL_A = Fraction(1, 48)                      # SZZ wall (beta_A)
def wall_W(N):
    return float(WALL_A) * N * N              # same wall in Wilson units


TARGET = {"SU(2)": (2, 2.5), "SU(3)": (3, 5.69)}   # crossover, Wilson units

print("=" * 78)
print("B03 · block-spin decimation (fixed b) — single-step beta and iteration")
print("=" * 78)

# ---------------------------------------------------------------- 0 conventions
print("\n[0] conventions")
print("  wall  beta_A = 1/48 = %.6f   ->  beta_W = N^2/48" % float(WALL_A))
for name, (N, bc) in TARGET.items():
    print("  %s: beta_W(wall) = %.5f ,  beta_W(crossover) = %.3f ,  ratio = %.2fx"
          % (name, wall_W(N), bc, bc / wall_W(N)))
print("  (2.5*48 = 120x is the convention-mixing error; 30x is like-for-like)")

# ------------------------------------------------------- 1 single-step d beta
print("\n[1] single step:  beta_W,eff = beta_W - C_N ln b ,  C_N = 11N^2/(12 pi^2)")
print("    %-6s %-4s %10s %14s %14s" % ("N", "b", "C_N", "d beta_W / step", "ratio at beta=1"))
for N in (2, 3):
    for b in (2, 3):
        print("    %-6d %-4d %10.5f %14.5f %14.5f"
              % (N, b, C(N), dbeta(N, b), 1.0 - C(N) * math.log(b) / 1.0))

print("\n    sign check: d beta_W / step < 0 for every N>=1, b>1  ->  beta_eff < beta.")
print("    (asymptotic freedom: g grows toward the IR, so beta_W = 2N/g^2 shrinks)")

print("\n[1b] ratio beta_eff/beta at the two ends of the window (b=2)")
print("     %-6s %12s %12s %12s %12s" % ("N", "beta_W", "beta_eff", "ratio", "x=log_b(ratio)"))
for name, (N, bc) in TARGET.items():
    for tag, beta in (("cross", bc), ("wall", wall_W(N))):
        be = beta_eff(N, 2, beta)
        r = be / beta
        print("     %-6s %12.5f %12.5f %12.5f %12s"
              % ("%s/%s" % (name, tag), beta, be, r,
                 ("%12.5f" % x_exponent(N, 2, beta)) if r > 0 else "   INVALID"))
print("     => at its own crossover the one-step ratio is 0.897 (N=2) / 0.898 (N=3): WEAKER,")
print("        not stronger (nearly N-independent because beta_c ~ N^2).  Note: 2.5 is the")
print("        SU(2) crossover; SU(3) would give 0.768 at beta_W=2.5 -- wrong pairing, unused.")
print("     => at the wall the one-loop formula returns ratio<0: ONE step leaves perturbation theory.")
print("     => no x>0 exists; the d=4 power ansatz beta*b^x is the wrong ansatz (see [5]).")

# ------------------------------------------------------------ 2 iteration
print("\n[2] iteration in the *physical* (coarse-graining) direction, from the crossover")
print("    (one-loop, constant step: valid only while g^2 = 2N/beta_W < 1)")
for name, (N, bc) in TARGET.items():
    beta = bc
    print("    %s (N=%d, b=2):" % (name, N))
    for k in range(4):
        print("      k=%d  beta_W=%8.4f  g^2=%6.3f  %s"
              % (k, beta, g2(N, beta), "one-loop OK" if g2(N, beta) < 1 else "OUT of one-loop domain"))
        beta = beta_eff(N, 2, beta)
    print("      steps until beta_W<=0 (formal): %.2f" % (bc / C(N) / math.log(2)))

print("\n[3] iteration in the *UV* (inverse-blocking) direction, wall -> crossover")
print("    k* = (beta_c - beta_wall) / (C_N ln b) ;  scale factor 2^k*")
print("    %-6s %10s %10s %8s %10s" % ("N", "beta_wall", "beta_c", "k*", "2^k*"))
for name, (N, bc) in TARGET.items():
    step = C(N) * math.log(2)
    kstar = (bc - wall_W(N)) / step
    print("    %-6d %10.5f %10.3f %8.3f %10.1f" % (N, wall_W(N), bc, kstar, 2 ** kstar))
print("    => ~9.4-9.5 steps for BOTH N: k* = (beta_c-beta_wall)/(C_N ln b) is N-independent")
print("       because beta_c ~ N^2 and C_N ~ N^2 cancel.  Ceiling: 10 steps.")
print("    BUT this direction is *adding* degrees of freedom (mu increases): it is not")
print("    the pushforward of a measure, it needs a UV completion.  No theorem; no map.")

# --------------------------------------------------- 4 what a power law would need
print("\n[4] why a power law cannot help at d=4 (engineering dimension)")
print("    beta_W = 2N/(g^2 a^(4-d)) ; at fixed g, a -> b a gives beta_eff = b^(4-d) beta")
for d in (2, 3, 4):
    print("      d=%d : exponent 4-d = %+d  ->  beta_eff/beta = b^%+d = %.4f (b=2)"
          % (d, 4 - d, 4 - d, 2.0 ** (4 - d)))
print("    d=4 is exactly marginal: exponent 0, so the leading term is beta_eff = beta.")
print("    The whole flow is logarithmic and its SIGN is the one-loop coefficient (AF).")
print("    d<4 is the opposite case: exponent 4-d>0, beta grows, IR weak coupling")
print("    (super-renormalisable) -- which is why only d=4 has the hard end at the IR.")
print("    d=2 'beta_eff = b^2 beta' is thus a dimensional artifact of the normalisation,")
print("    not a flow statement: in 2D the physical gap m = C2 g^2/2 does not move.")

# --------------------------------------------- 5 cross-check vs the repo anchor
print("\n[5] cross-check against the stored activity-matched step (T1-POLYMER-RG §4, §6)")
u, u_eff, beta_anchor, beta_c2 = 1.014, 0.988, 2.40, 2.5
print("    activity-matched:  u(2.5)=%.3f -> u_eff=%.3f  =>  beta_eff^act=%.2f, ratio=%.3f"
      % (u, u_eff, beta_anchor, beta_anchor / beta_c2))
print("    coupling one-loop:                             beta_eff    =%.3f, ratio=%.3f"
      % (beta_eff(2, 2, beta_c2), beta_eff(2, 2, beta_c2) / beta_c2))
print("    two independent definitions of the same step -> both <1 (0.96 / 0.90).")
print("    sign is scheme-independent; the SIZE is scheme-dependent. No scheme gives >1.")

# ------------------------------------------------------------ 6 four-model test
print("\n[6] four-model test of a block-spin-type quantity (Q_BS := one-step beta flow)")
rows = [
    ("free field", "0", "beta=inf, g=0 exactly: no flow, Q_BS=0", "[T]"),
    ("2D YM", ">0", "exact Migdal blocking, Q_BS!=0 and m=C2 g^2/2>0", "[T]"),
    ("M6", "NOT 0", "Gaussian marginal is an exact RG fixed point: the a-marginal of M6",
     "[T]"),
    ("M6 (cont.)", "", "equals that of the gapped HO (sup|drho|=2.3e-14, U09 §3)", ""),
    ("single plaquette", "undefined", "no thermodynamic limit, no tower", "[T]"),
]
for r in rows:
    print("    %-18s %-10s %-58s %s" % r)
print("    => two horns, both fatal: (i) if Q_BS sees only the blocked (a) marginal, then")
print("       Q_BS(M6) = Q_BS(HO) exactly (2.3e-14) => it cannot read 0 on M6 and >0 on HO;")
print("       (ii) if it sees the free (b) direction, it must be min-type => it is the gap")
print("       itself (U05) => circular.  4-model: FAIL on both horns.")

# --------------------------------------------------------------- 7 verdict
print("\n[7] verdict")
print("    (a) beta_eff > beta  with a theorem : NO  -> (a) not triggered")
print("    (b) fixed-b block-spin theorems exist for phi^4_4 (GK) and 4D U(1) (Gopfert-Mack,")
print("        strong-coupling/dipole route); for SU(N) at fixed b: none  -> (b) triggered")
print("    (c) beta_eff < beta : YES, one loop exactly, ratio 0.90 (N=2) at the crossover;")
print("        the only beta-raising direction is inverse blocking (UV), 10 steps, no map")
print("        -> (c) TRIGGERED: the flow points away from 2.5.")

# ----------------------------------------------------------------- 8 assertions
assert all(dbeta(N, b) < 0 for N in (1, 2, 3) for b in (2, 3, 4))
assert all(beta_eff(N, 2, 2.5) < 2.5 for N in (2, 3))
assert abs(x_exponent(2, 2, 2.5) + 0.15685) < 5e-4
k2 = (2.5 - wall_W(2)) / (C(2) * math.log(2))
k3 = (5.69 - wall_W(3)) / (C(3) * math.log(2))
assert 9.0 < k2 < 10.0 and 9.0 < k3 < 10.0 and abs(k2 - k3) < 0.2
assert 2.0 ** (4 - 4) == 1.0
print("\n    all self-checks passed")
print("=" * 78)
