"""
cert_chi_lattice.py -- CERT-CHI-LATTICE witness (2026-11-24)
The a->0 extrapolation is NOT a separate tool; it IS dimensional transmutation = H1.
Lightweight, numpy only. Not a Yang-Mills computation -- it makes the SCALING structure exact.
"""
import numpy as np

BOUND = (2 * np.pi) ** 2  # compact-group diameter bound, BLAST-LATTICE-CERT


def rel_err(x, y):
    return abs(x - y) / max(abs(y), 1e-300)


print("=" * 72)
print("CERT-CHI-LATTICE: the a->0 wall = chi_lattice = O(a^2) = dimensional transmutation = H1")
print("=" * 72)

print(f"\n[1] COMPACTNESS (BLAST-LATTICE-CERT): 0 <= chi_lattice <= (2*pi)^2 = {BOUND:.3f}, a-INDEPENDENT.")
print("    This is a flat ceiling: it does NOT shrink as a -> 0.")

print("\n[2] PHYSICAL FINITENESS: with a mass scale m, dimensional analysis gives")
print("    chi_lattice(a) = c * (a*m)^2  (c an O(1) number),  so chi_lattice -> 0 like a^2 as a -> 0.")
print("    Free massive field (exact): chi_phys = 1/m^2  =>  chi_lattice(a) = a^2/m^2.")
m = 1.0
print("      a       chi_lattice=a^2/m^2   chi_lattice/a^2 (=1/m^2)")
for a in [0.5, 0.2, 0.1, 0.05, 0.02]:
    cl = a ** 2 / m ** 2
    print(f"      {a:.3f}   {cl:.3e}          {cl / a**2:.6f}")

print("\n[3] THE WALL. Physical demand: chi_phys <= 1/lam  <=>  chi_lattice <= a^2/lam.")
print("    Compactness gives O(1); the demand is O(a^2). The gap between them IS the a^2 coefficient.")
lam = 1.0
print("      a     compactness(4pi^2)   demand(a^2/lam)   demand/compactness")
for a in [0.5, 0.2, 0.1, 0.05, 0.02]:
    dem = a ** 2 / lam
    print(f"      {a:.3f}  {BOUND:12.3f}        {dem:.3e}       {dem / BOUND:.3e}")
print("    -> the demand shrinks like a^2; compactness is flat. The a->0 gap = the a^2 coefficient")
print("       = chi_phys = the inverse mass scale = the gap itself.")

print("\n[4] STRONG-COUPLING SIDE (proven region, beta<1/48): chi_lattice = O(1), a-INDEPENDENT.")
print("    There is NO small-a regime in the rigorously controlled region.")
print("    Any a->0 extrapolation MUST cross the crossover, where the O(a^2) coefficient (the IR mass")
print("    scale) is exactly the nonperturbative input -> H1.")

print("\n[5] THE THREE NAMED ROUTES ALL FAIL:")
print("    (i) RENORMALIZATION CONDITION: fixes a UV coupling at scale mu; chi_phys depends on the IR")
print("        mass scale, which no UV condition supplies.")
print("    (ii) STEP-SCALING: the RG flow gives the ANOMALOUS DIMENSION of trF^2 (perturbative, UV), not")
print("         the CONSTANT (the IR mass scale).")
print("    (iii) MONOTONE QUANTITY IN a: the only candidates (effective curvature kappa(a), R_vac(a)) ARE the")
print("          existing H1 attempts (TRANSPORT-SCALAR: R_vac monotone = H1); they read the gap, not force it.")

# Numeric sanity: show that the 'a->0' extrapolation is underdetermined without the gap.
# Fit candidate chi_lattice(a) = C*a^p to 'data' at strong coupling is impossible: there is no a->0 data.
print("\n[6] UNDETERMINED WITHOUT THE SCALE: the exponent p in chi_lattice ~ a^p is set by the")
print("    physical mass (p=2 for a mass gap; p=0 at strong coupling). Both are compatible with")
print("    compactness 0<=chi_lattice<=4pi^2. The exponent IS the gap; compactness cannot fix it.")

print("\n" + "=" * 72)
print("VERDICT: the a->0 extrapolation is dimensional transmutation = H1, not a separate certificate.")
print("  compactness: chi_lattice <= 4pi^2 (flat);  demand: chi_lattice <= a^2/lam (a^2-small).")
print("  O(a^2) scaling  <=>  chi_phys finite  <=>  the gap.  No renorm-cond / step-scaling / monotone")
print("  supplies the a->0 certificate.  Clay 0%; no QED.")
print("=" * 72)
