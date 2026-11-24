"""
base_transport_scalar.py -- witness for ROOF-TRANSPORT-SCALAR.md  (2026-11-24)

Four tiny checks (numpy only, <1 s). NOT a Yang-Mills calculation -- it isolates
the LOGICAL structure of "transport the base's single scalar across the crossover."

A. NAIVE SCALAR IS WRONG.  The lattice-unit inverse susceptibility 1/sqrt(chi)
   = m_lattice.  Holding the PHYSICAL gap fixed and a(beta)=exp(-c beta)
   (asymptotic freedom), m_lattice = a*m_phys -> 0 as beta->inf.  So "1/sqrt(chi)
   stays >0 across [1/48,inf)" is unachievable for a SCALING reason, not a gap loss.

B. BETA-CONTINUITY IS ACID-DIRTY.  On a finite lattice the free field has box gap
   4 sin^2(pi/n) > 0 for every finite n (compactness => "positive at fixed L"),
   but it -> 0 as n->inf.  So positivity-by-continuity returns the FREE-FIELD box
   gap = vacuous (the BASE-PF acid failure), and gives nothing at L->inf.

C. THE TRANSPORT REDUCES TO ONE MONOTONICITY.  Write R_vac = exp(-DeltaF(beta)).
   If DeltaF is NON-INCREASING (R_vac non-decreasing), inf over [b0,2.5] = R_vac(b0)
   = the strong-coupling value -> the floor TRANSPORTS and the roof closes.
   If DeltaF is NON-DECREASING, inf sits at 2.5 -> needs a weak-coupling IR anchor
   that does not rigorously exist.  Show both infima.

D. LATTICE-MEASURED f1^2 IS NOT MONOTONE.  A hump (rises then "dissolves into the
   cutoff") -- the empirical beta=1.0 dissolution.  Distinguish lattice-measured
   (non-monotone, resolution-limited) from physical R_vac (could be monotone).
"""
import numpy as np

print("=" * 74); print("A. naive scalar 1/sqrt(chi)=m_lattice -> 0 (wrong transport variable)"); print("=" * 74)
c = 0.20
for beta in [1/48, 0.1, 0.5, 1.0, 2.5, 5.0, 10.0]:
    a = np.exp(-c * beta)           # a(beta) ~ exp(-c beta)
    print(f"  beta={beta:8.4f}   a={a:.3e}   m_lattice=a*m_phys(1)={a:.3e}")
print("  => fixed physical gap, m_lattice -> 0. 'stays >0' fails by SCALING, not by a gap loss.")

print("\n" + "=" * 74); print("B. box-gap acid: free field Z_n gap = 4 sin^2(pi/n) > 0 at finite n, -> 0 as n->inf"); print("=" * 74)
for n in [4, 8, 16, 32, 64, 128]:
    g = 4.0 * np.sin(np.pi / n) ** 2
    print(f"  n={n:4d}   gap={g:.6e}")
print("  => 'positive at fixed L' is TRIVIALLY true (box gap) and ACID-DIRTY (free field has it).")

print("\n" + "=" * 74); print("C. transport = monotonicity of R_vac = exp(-DeltaF(beta))"); print("=" * 74)
b0 = 1/48
beta = np.linspace(b0, 2.5, 200001)
D0, A = 0.8, 1.5
# case 1: DeltaF NON-INCREASING in beta => R_vac NON-DECREASING (inf at b0)
D1 = D0 + A * (2.5 - beta) / (2.5 - b0)
Rvac1 = np.exp(-D1)
# case 2: DeltaF NON-DECREASING in beta => R_vac NON-INCREASING (inf at 2.5)
D2 = D0 + A * (beta - b0) / (2.5 - b0)
Rvac2 = np.exp(-D2)
print(f"  beta in [{b0:.5f}, 2.5]")
print(f"  case1 R_vac non-decreasing: R_vac(b0)={Rvac1[0]:.5f}  R_vac(2.5)={Rvac1[-1]:.5f}  inf={Rvac1.min():.5f} (at b0)  -> TRANSPORTS")
print(f"  case2 R_vac non-increasing: R_vac(b0)={Rvac2[0]:.5f}  R_vac(2.5)={Rvac2[-1]:.5f}  inf={Rvac2.min():.5f} (at 2.5) -> needs RIGHT anchor")
print(f"  both R_vac in (0,1]: case1 [{Rvac1.min():.3f},{Rvac1.max():.3f}]  case2 [{Rvac2.min():.3f},{Rvac2.max():.3f}]")
print("  => if R_vac (=f1^2/W_tube) is NON-DECREASING in beta, the strong-coupling floor extends to [b0,inf).")
print("     The whole transport is this ONE monotonicity statement.")

print("\n" + "=" * 74); print("D. lattice-measured f1^2 is NON-monotone (hump + dissolution)"); print("=" * 74)
bb = np.linspace(0.05, 3.0, 7)
f1 = 0.0 + 0.45 * (1 - np.exp(-3.0 * bb)) * np.exp(-0.90 * np.maximum(bb - 1.0, 0))  # rises then dissolves
for x, y in zip(bb, f1):
    print(f"  beta={x:5.2f}   f1^2_meas={y:.4f}")
print("  => not monotone: an apparent dissolution near beta~1 (a LATTICE-RESOLUTION artifact,")
print("     atom leaking into the cutoff continuum). Physical R_vac may still be monotone --")
print("     that gap between measurement and proof is exactly the open step.")
