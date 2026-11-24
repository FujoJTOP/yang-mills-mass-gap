"""
falsify_entropy.py -- FALSIFY-ENTROPY witness.
Does any NONLINEAR information functional (differential entropy, Fisher information,
relative entropy, Renyi) LOWER-BOUND the atom weight f1^2 from FINITE data, escaping
both YM-structure=H1 and J->infinity=H1?

Witness: an atom f2*delta_{s0} + uniform continuum is smoothed to a Gaussian bump of
width eps at s0.  As eps -> 0 the measure -> (atom + continuum).  We show:
  (1) FINITE moments M_0..M_3 are essentially eps-insensitive => any finite-moment
      functional is blind (OPE-RHO / Delta_min=0);
  (2) differential entropy H -> -inf and Fisher J -> +inf as eps -> 0: these are
      PRESENCE detectors via divergence, requiring the full density (eps->0 = infinite
      data = J->infty = H1), NOT finite-data forcing of the weight.
Run:  py -3.12 -I research/unconditional/w7/falsify_entropy.py
"""
import numpy as np

N = 8001
s = np.linspace(0.5, 5.0, N)
ds = s[1] - s[0]
f2, lo, hi = 0.5, 1.0, 4.0


def density(eps):
    uni = np.where((s >= lo) & (s <= hi), (1 - f2) / (hi - lo), 0.0)
    bump = f2 * np.exp(-0.5 * ((s - lo) / eps) ** 2) / (eps * np.sqrt(2 * np.pi))
    return uni + bump


def moments(rho):
    return np.array([np.trapezoid(rho * s ** j, s) for j in (0, 1, 2, 3)])


def entropy(rho):
    m = rho > 1e-300
    return -np.trapezoid(np.where(m, rho * np.log(np.where(m, rho, 1.0)), 0.0), s)


def fisher_analytic(eps):
    # uniform part has zero Fisher; a normalized Gaussian bump of width eps
    # has Fisher 1/eps^2, so the mass-f2 bump contributes f2/eps^2.
    return f2 / eps ** 2


print("=" * 78)
print("eps      M0       M1       M2       M3        H(dens)    Fisher(analytic)")
print("=" * 78)
ref = None
rows = []
for eps in (0.20, 0.10, 0.05, 0.02):
    rho = density(eps)
    M = moments(rho)
    if ref is None:
        ref = M
    H = entropy(rho)
    J = fisher_analytic(eps)
    relM = np.max(np.abs(M - ref) / np.abs(ref))
    rows.append((eps, M, H, J, relM))
    print(f"{eps:5.3f}  " + "  ".join(f"{x:8.4f}" for x in M) + f"   {H:9.4f}  {J:12.1f}")

print()
print("Max relative change of the FINITE moments M0..M3 from eps=0.20 to eps=0.02:")
print(f"  {rows[-1][4]:.3e}   (=> finite-moment data is essentially blind to the atom)")
print()
print("As eps -> 0, H -> -inf and Fisher -> +inf (presence via divergence),")
print("but that needs the full density (eps->0 = infinite data = J->infty = H1).")
print()
print("Analytic Fisher of a Gaussian bump width eps ~ 1/eps^2:")
for eps in (0.20, 0.10, 0.05, 0.02):
    print(f"  eps={eps:5.3f}: 1/eps^2 = {1/eps**2:10.1f}")
