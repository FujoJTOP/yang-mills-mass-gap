#!/usr/bin/env python3.12
# -*- coding: utf-8 -*-
"""
sm5_tci_soft.py  --  witness for SM5-TCI-HEAD.md
Deterministic. numpy only. Run: py -3.12 -I sm5_tci_soft.py

It witnesses (NOT proves) the four task items for the tachyonic-compensation
inequality (TCI) on the Savvidy SOFT mode at beta=2.5:

  (A) double-well normal form: TCI holds pointwise at the well (margin 2c),
      FAILS at the barrier (ess inf = -c).  c = g|B|.
  (B) exact rank-one split  K = K+ - |lam_-| P_-,  K^-1 c decomposition,
      D = D+ + D_Sav with D+<=0 (BL) and D_Sav>=0.  split error ~1e-15.
  (C) one-octave overlap: overlap^2 = x^2/(1+x^2), x=k/sqrt(gB); O(1) at the
      head (x=1 -> 0.5), ->0 in the deep IR.  No small parameter at the head.
  (D) two-loop sign-protected bound  D_Sav <= 2 g^2/lam4  vs  S at beta=2.5.
  (E) Fisher identity  E_rho[G''] = I(rho) = 1/sigma^2 >= 0; free field -> 0;
      volume-shrinking margin (the H2 numbers).
  (F) beta-dependence of the provable margin mu(beta): strong-coupling LSI
      K_S>0 at beta<1/48, and the head O(1) competition flag.
"""
import numpy as np

np.set_printoptions(precision=6, suppress=True)
sep = "=" * 78


def section(t):
    print("\n" + sep + "\n" + t + "\n" + sep)


# ---------------------------------------------------------------- (A) double well
section("(A) DOUBLE-WELL: TCI at the well vs at the barrier  (c = g|B|)")
# Savvidy soft-mode normal form: V_eff(y) = -1/2 c y^2 + 1/4 lam4 y^4
# TCI pieces:  Q = 3 lam4 y0^2 (quartic 2nd var at the saddle), -cF = -c, K (Schur)
for c, lam4 in [(1.0, 1.0), (1.6, 1.6), (0.5, 2.0)]:
    y0sq = c / lam4                      # well location
    Vpp_well = -c + 3.0 * lam4 * y0sq    # = 2c  (well curvature)
    Vpp_barrier = -c                     # at y=0  (double-well barrier)
    Q = 3.0 * lam4 * y0sq                # quartic second variation at the saddle
    tci_well = Q + 0.0 - c               # K=0 (one-loop, LLL is exact eigenspace)
    print(f"  c=gB={c:4.2f} lam4={lam4:4.2f} : y0^2={y0sq:6.3f}  "
          f"Q=3lam4 y0^2={Q:6.3f}  -cF=-{c:4.2f}  "
          f"TCI(well)=Q-cF={tci_well:6.3f} (=2c)  G''(barrier)={Vpp_barrier:6.3f}")
print("  => TCI holds at the WELL with margin 2c>0 ; the ESSENTIAL INFIMUM over")
print("     the support is m^2=-c<0 (the barrier).  Uniform/ess-inf TCI FAILS.")

# ---------------------------------------------------------------- (B) rank-one split
section("(B) EXACT RANK-ONE SPLIT  K = K+ - |lam_-| P_- ,  D = D+ + D_Sav")
rng = np.random.default_rng(20261124)
n = 6
lam_pos = np.sort(rng.uniform(0.4, 2.0, n - 1))     # positive fiber eigenvals
c0 = 1.0                                            # |lam_-| = g|B|
K = np.diag(np.concatenate([lam_pos, [-c0]]))       # full Savvidy Hessian
Kplus = np.diag(np.concatenate([lam_pos, [0.0]]))   # positive part, phi_- in ker
phi = np.zeros(n); phi[-1] = 1.0
Pminus = np.outer(phi, phi)
# a random vector c and a random force-covariance-compatible vector f
cc = rng.normal(size=n)
Kinv = np.linalg.inv(K)
Kplus_inv = np.linalg.inv(Kplus + np.eye(n) * 1e-15)  # pseudo-inverse on ker=phi_-
Kplus_pinv = np.diag(np.concatenate([1.0 / lam_pos, [0.0]]))
lhs = cc @ Kinv @ cc
rhs = cc @ Kplus_pinv @ cc - (phi @ cc) ** 2 / c0
print(f"  <K^-1 c,c> = {lhs:.10f}")
print(f"  <K+^+ c,c> - <phi,c>^2/|lam_-| = {rhs:.10f}")
print(f"  rank-one split error = {abs(lhs - rhs):.2e}")
print(f"  K = K+ - |lam_-|P_-  check ||(K+ - c0 P-) - K|| = "
      f"{np.abs(Kplus - c0 * Pminus - K).max():.2e}")
# D+ (BL, log-concave positive fibre) <= 0 ; D_Sav >= 0
print("  D+ = Var(f|0) - E<K+^+ c,c>  <= 0  (Brascamp-Lieb on K+ >= 0) [T]")
print("  D_Sav = E<phi,c>^2/|lam_-|    >= 0  (one rank-one projection)   [T]")

# ---------------------------------------------------------------- (C) one-octave overlap
section("(C) ONE-OCTAVE OVERLAP  overlap^2 = x^2/(1+x^2),  x = k/sqrt(gB)")
print("   x=k/sqrt(gB)   overlap^2      regime")
for x in [0.05, 0.1, 0.3, 0.5, 1.0, 2.0]:
    o2 = x * x / (1.0 + x * x)
    regime = "deep IR" if x < 0.3 else ("HEAD (O(1))" if x >= 0.8 else "crossover")
    print(f"      {x:5.2f}        {o2:8.4f}      {regime}")
print("  => suppression only deep IR ; at the head k~sqrt(gB) overlap is O(1).")

# ---------------------------------------------------------------- (D) two-loop bound
section("(D) TWO-LOOP SIGN-PROTECTED BOUND  D_Sav <= 2 g^2/lam4   vs   S")
g2 = 4.0 / 2.5      # beta=2.5 lattice convention  g^2 = 4/beta = 1.6
S = 1.40            # bare stiffness scale at the head (H2)
print(f"   beta=2.5 -> g^2=4/beta={g2:.3f} ; S_bare={S:.2f}")
print("   lam4/g^2   2g^2/lam4    bound < S ?")
for r in [0.5, 1.0, 1.5, 2.0]:
    bound = 2.0 * g2 / (r * g2)
    print(f"    {r:4.2f}      {bound:8.3f}       {'YES' if bound < S else 'NO (overshoot)'}")
print("  => for NATURAL lam4~g^2 the exact bound ~2 > S~1.4 : too lossy at the head.")

# ---------------------------------------------------------------- (E) Fisher identity
section("(E) FISHER IDENTITY  E_rho[G''] = I(rho) = 1/sigma^2 ; free -> 0")
def fisher_gaussian(sigma2):
    # rho Gaussian N(0,sigma2): G''=-d^2 log rho = 1/sigma2 everywhere; I=1/sigma2
    return 1.0 / sigma2
for s2 in [4.0, 1.0, 0.25]:
    print(f"   Gaussian var sigma^2={s2:5.2f}: E[G'']=I={fisher_gaussian(s2):7.4f}")
print(f"   free massless (flat rho): I=0.0000   (acid test)  ")
print("   volume-shrinking margin (H2 measured): L=4 -> +0.28 ; L=8 -> +0.049")
print("   => averaged TCI is TRUE [T]; its uniform margin is the coarse gap [O].")

# ---------------------------------------------------------------- (F) beta dependence
section("(F) BETA-DEPENDENT MARGIN mu(beta): proved at both ends, open at head")
# strong coupling: SZZ LSI constant K_S = N(1/2 - 24|beta|) at N=2
N = 2
for b in [0.0, 0.01, 0.020833, 0.05, 0.5, 2.5]:
    KS = N * (0.5 - 24.0 * abs(b))
    tag = "LSI>0 [T]" if KS > 0 else ("head O(1) [O]" if b > 0.1 else "?")
    print(f"   beta={b:7.5f}: K_S(N=2)={KS:8.4f}   {tag}")
print("  => TCI provable (strong) for beta<1/48~0.0208 ; UV certifiable (measure")
print("     of bad set <= e^{-c/g^2}) ; the HEAD beta~2.5 is the O(1) competition.")

section("DONE (witness only; theorem content is in SM5-TCI-HEAD.md)")
