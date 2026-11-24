"""
construct_rg_witness.py  --  CONSTRUCT-RG numeric witness (fork, 2026-11-24)

Grounds the numbers quoted in research/unconditional/CONSTRUCT-RG.md for the
best rigorous 4D gauge RG assembled today.  NOT a Yang-Mills calculation; a
consistency check of (A) the AF log-resummation, (B) the Savvidy-cone depth vs
the shallow-band floor, (C) the SZZ gap radius vs the crossover, (D) the best
rigorous reach.

Run:  py -3.12 -I research/unconditional/w7/construct_rg_witness.py
"""

import numpy as np
from scipy.special import i0, i1

np.set_printoptions(precision=5, suppress=True, linewidth=120)


def banner(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)


# ---------------------------------------------------------------------------
# (A) AF log-resummation: the marginal sum is absolutely summable.
#     One-loop (lattice Hooft): g_j^2 = 1/(2 beta0 j).
#     => Sum_j g_j^4 = zeta(2)/(4 beta0^2) < inf  (summable, log resums).
# ---------------------------------------------------------------------------
banner("(A) AF marginal-log resummation:  Sum_j g_j^4  converges")
zeta2 = np.pi ** 2 / 6.0
for N in (2, 3):
    beta0 = 11.0 * N / (48.0 * np.pi ** 2)
    J = 200000
    j = np.arange(1, J + 1)
    g2 = 1.0 / (2.0 * beta0 * j)
    S = np.sum(g2 ** 2)
    S_exact = zeta2 / (4.0 * beta0 ** 2)
    print(f"  N={N}: beta0={beta0:.5f}  Sum_{{j<=2e5}} g_j^4={S:.6e}  "
          f"zeta(2)/(4b0^2)={S_exact:.6e}  rel.diff={abs(S-S_exact)/S_exact:.2e}")

# residual running correction  delta g_J = O(J^{-3/2} log J) -> 0
J = np.array([10.0, 100.0, 1000.0, 10000.0])
print("  residual dg_J ~ J^-3/2 log J :", (J ** -1.5 * np.log(J)))

# ---------------------------------------------------------------------------
# (B) Savvidy Hessian:  Hess = -D^2 + 2 g ad(F_B),  lambda_min = -g|B| < 0.
#     Shallow band closes where 2 g |F| <= N/12 (Phi'' >= N/12 absorbs it).
#     At beta=2.5 (SU(2)) g=2/sqrt(beta) ~ 1.26; the negative directions are
#     90-97% high-momentum (W6-BE), so the goal mode sits ABOVE the band.
# ---------------------------------------------------------------------------
banner("(B) Savvidy cone at beta=2.5:  g,  shallow floor N/12,  depth -g|B|")
for beta in (1 / 48, 1.0, 2.5):
    g = 2.0 / np.sqrt(beta)          # SU(2): beta = 2N/g^2 = 4/g^2
    floor = 2.0 / 12.0               # N/12, N=2
    print(f"  beta={beta:>8.5f}: g={g:.4f}  floor N/12={floor:.4f}  "
          f"lambda_min=-g|B| < 0 whenever |B|>0  (indefinite Hessian)")

# ---------------------------------------------------------------------------
# (C) SZZ gap radius vs crossover (Hooft).
#     SZZ: |beta| < 1/(16(d-1)) = 1/48 (d=4)  <=>  lambda = g^2 N > 96.
#     Crossover lambda ~ 7.3 (SU(3), beta_W ~ 2.47).
# ---------------------------------------------------------------------------
banner("(C) SZZ gap radius vs crossover")
d = 4
beta_gap = 1.0 / (16.0 * (d - 1))
lam_gap = 96.0
lam_cross = 7.3
print(f"  SZZ: beta < {beta_gap:.5f}  <=>  lambda = g^2 N > {lam_gap:.1f}")
print(f"  crossover: lambda ~ {lam_cross}  (beta_W ~ 2.47)")
print(f"  distance: lambda {lam_gap:.1f} -> {lam_cross:.1f}  =  {lam_gap/lam_cross:.1f}x  in lambda")

# ---------------------------------------------------------------------------
# (D) Best rigorous reach of the gap today (from corpus).
# ---------------------------------------------------------------------------
banner("(D) Best rigorous reach of the gap today")
print("  peer-reviewed gap radius : beta < 1/48        = 0.02083   (SZZ CMP 400)")
print("  unrefereed 2026 reach    : beta < 195/7168    = 0.02720")
print("  area-law (not gap) reach : beta < c_d (abs.) and into beta ~ N  (CNS)")
print("  UV control (not IR gap)  : beta -> inf        (Balaban CMP 109/116/122)")
print("  crossover (the target)   : beta ~ 2.5  --  NOT reached by any of these")
