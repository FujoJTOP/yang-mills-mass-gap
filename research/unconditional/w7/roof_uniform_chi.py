"""
roof_uniform_chi.py -- witness for research/unconditional/ROOF-UNIFORM-CHI.md
ONE tiny job (numpy only, <1s). Exact 1D facts underlying the roof.

Three checks:
 (1) CHEEGER FRACTION IS EXACTLY 1/(2pi). For the Gaussian invariant density,
     h = inf_x p(x)/min(F(x),1-F(x)); verify h^2/4 = gamma/(2pi), gamma=1.
 (2) BOUNDEDNESS IS VACUOUS (acid). Uniform on [-A,A]: gamma = pi^2/(2A)^2 -> 0
     and 1/chi = 3/A^2 -> 0 as A->inf. A bounded state space gives a KINEMATIC
     lattice gap that vanishes in the continuum; it is NOT a physical floor.
 (3) Cheeger inequality gamma >= h^2/4 holds in both examples.

No Yang-Mills claim here.  Illustrative identities only.  [T]/[V-num].
"""
import numpy as np

def cheeger_const(p, F, mfloor=1e-3):
    # h = inf_x p(x) / min(F(x), 1-F(x));  F = CDF.
    # Restrict to the interior (m>mfloor) to avoid the finite-window residual
    # artifact at the edges; the true inf of a symmetric unimodal density is at
    # the mode (m=1/2), so this restriction does not change h.
    m = np.minimum(F, 1.0 - F)
    ok = m > mfloor
    return np.min(p[ok] / m[ok])

def banner(s):
    print("\n" + "=" * 70); print(s); print("=" * 70)

# ---------------------------------------------------------------- (1) Gaussian
banner("(1) Gaussian: cheeger fraction = 1/(2pi) exactly")
for sigma in (0.5, 1.0, 2.0):
    x = np.linspace(-12 * sigma, 12 * sigma, 2_000_001)
    p = np.exp(-x**2 / (2 * sigma**2)) / (sigma * np.sqrt(2 * np.pi))
    dx = x[1] - x[0]
    F = np.cumsum(p) * dx
    F = F - F[0] + p[0] * dx  # left tail negligible on this window
    F = np.clip(F, 0, 1)
    h = cheeger_const(p, F)
    gamma = 1.0 / sigma**2                 # OU spectral gap for N(0,sigma^2)
    chi = sigma**2
    print(f"  sigma={sigma:4.2f}  h={h:.6f}  h*sqrt(chi)={h*np.sqrt(chi):.6f} "
          f"(sqrt(2/pi)={np.sqrt(2/np.pi):.6f})")
    print(f"           gamma={gamma:.6f}  (h^2/4)/gamma={(h**2/4)/gamma:.6f} "
          f"(1/2pi={1/(2*np.pi):.6f})   1/chi={1/chi:.6f}")

# ---------------------------------------------------------------- (2) Uniform
banner("(2) Uniform on [-A,A]: kinematic gap & 1/chi both vanish as A->inf")
print("   A      gamma=pi^2/(2A)^2    1/chi=3/A^2    h=1/A    h^2/4 <= gamma?")
for A in (1.0, 10.0, 100.0, 1e3):
    x = np.linspace(-A, A, 2_000_001)
    p = np.full_like(x, 1.0 / (2 * A))
    dx = x[1] - x[0]
    F = np.clip(np.cumsum(p) * dx, 0, 1)
    h = cheeger_const(p, F)
    gamma = np.pi**2 / (4 * A**2)          # Neumann gap, interval length 2A
    chi = A**2 / 3.0
    print(f"  {A:7.1f}   {gamma:12.3e}   {1/chi:12.3e}   {h:9.3e}   "
          f"{h**2/4 <= gamma + 1e-12}")

# ---------------------------------------------------------------- (3) summary
banner("(3) reading")
print("""  (1) On the Gaussian the Cheeger fraction h^2/4 = gamma/(2pi) is exact at
      every sigma -> the base's 1/(2pi) is the true shape constant, scale-free.
  (2) On a bounded interval (compact group / bounded logU) gamma~1/A^2 -> 0 and
      1/chi~1/A^2 -> 0 together as A->inf (A ~ 1/a_ell).  Boundedness alone gives
      a lattice gap that VANISHES in the continuum: it can never supply a physical
      floor 1/chi >= lambda.  Confirms R1 is acid-vacuous.
  => The roof's only missing strut is the MAGNITUDE of 1/chi_phys at beta=2.5,
     i.e. chi_phys <= 1/lambda uniformly = H1 verbatim.  No QED.
""")
