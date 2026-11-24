"""
star_polchinski.py -- witness for research/unconditional/STAR-POLCHINSKI.md

Two cheap checks (numpy only):

(A) Pure YM 2-loop beta function has NO zero for g>0 (no Banks-Zaks fixed point):
        beta(g) = -(b0 g^3 + b1 g^5)/(16 pi^2),
    b0 = 11 N/3, b1 = 34 N^2/3  (both > 0 for pure SU(N)).
    Show beta(g) < 0 for all g in (0, gmax], for N=2,3,4.

(B) The trace anomaly T^mu_mu = (beta/2g) tr F^2 + (Weyl^2, Euler).
    IR conformality requires the ANOMALY OPERATOR to vanish in the IR, i.e. beta(g_IR)=0.
    Since beta(g)<0 for all g>0 (no fixed point), there is NO g_IR>0 at which the anomaly vanishes
    perturbatively -> the "no conformal IR" step (b) is not delivered by the perturbative beta alone.

This is an illustrative/schematic witness for the trace-anomaly logic; NOT a Yang-Mills theorem.
"""
import numpy as np

# ---------------------------------------------------------------- beta function
def beta_2loop(g, N):
    b0 = 11.0 * N / 3.0
    b1 = 34.0 * N * N / 3.0
    return -(b0 * g**3 + b1 * g**5) / (16.0 * np.pi**2)

print("=" * 70)
print("(A) Pure YM 2-loop beta: is there a zero at g>0? (Banks-Zaks would need beta(g*)=0)")
print("=" * 70)
g = np.linspace(0.0, 6.0, 60001)
for N in (2, 3, 4):
    b = beta_2loop(g, N)
    has_zero_at_gpos = np.any(np.abs(b[1:] < 0) == False) and np.any(b[1:] >= 0.0)
    # more directly: is beta ever >= 0 for g>0 ?
    any_nonneg = np.any(b[1:] >= 0.0)
    bmin, bmax = b[1:].min(), b[1:].max()
    print(f"  N={N}: beta(g)<0 for all g in (0,6]? {not any_nonneg}   "
          f"range on g>0: [{bmin:.3e}, {bmax:.3e}]")
    print(f"       => perturbative IR fixed point (beta=0, g*>0) exists? {any_nonneg}")

# ---------------------------------------------------------------- anomaly logic
print()
print("=" * 70)
print("(B) Trace anomaly  T^mu_mu = (beta/2g) tr F^2 + (Weyl^2, Euler).")
print("    IR conformality  <=>  anomaly operator vanishes in IR  <=>  beta(g_IR)=0.")
print("=" * 70)
# At a putative IR fixed point g_IR, the anomaly coefficient beta(g_IR)/(2 g_IR) must be 0.
# Scan for any g_IR>0 with beta(g_IR) ~ 0 (perturbatively), for N=3:
N = 3
ghat = np.linspace(0.2, 6.0, 58001)
bg = beta_2loop(ghat, N)
idx = np.argmin(np.abs(bg))
g_closest = ghat[idx]
print(f"  N={N}: g>0 minimizing |beta(g)| is g={g_closest:.4f} with beta={bg[idx]:.4e} "
      f"(anomaly coeff beta/2g = {bg[idx]/(2*g_closest):.4e}).")
print("  beta never reaches 0 at g>0  => no perturbative g_IR with vanishing anomaly.")
print("  => step (b) ['no conformal IR'] is NOT closed by the perturbative beta;")
print("     a strongly-coupled CFT IR would need a nonperturbative argument.")
print()
print("  NOTE: the anomaly also carries the Weyl^2 / Euler (curvature) terms; in flat space")
print("  the obstruction to T^mu_mu=0 is precisely the nonzero operator (beta/2g) tr F^2.")
