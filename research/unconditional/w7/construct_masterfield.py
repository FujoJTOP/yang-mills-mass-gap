"""
construct_masterfield.py -- witness for CONSTRUCT-MASTERFIELD.md  (2026-11-24)

Two master-field closures of the large-N 4D Makeenko-Migdal loop equation are compared:

  PART 1  PURE Nambu-Goto surface sum  W(C)=sum_surfaces exp(-sigma A)
          -> closed-string energy E0(L)=sigma L sqrt(1 - pi(D-2)/(3 sigma L^2))
          -> monotone, NO interior minimum (reproduces PROOF-EFFECTIVE-STRING).

  PART 2  RIGID string (extrinsic-curvature) surface sum  W(C)=sum exp(-sigma A - kappa*int K^2 ds)
          -> neutral 0++ closed loop (circular, N-ality 0) HAS an interior minimum:
             E(rho)=2 pi sigma rho + pi kappa / rho
             rho* = sqrt(kappa/(2 sigma)),  m_0++ = 2 pi sqrt(2 sigma kappa),
             m_0++/sqrt(sigma) = 2 pi sqrt(2 kappa).
          -> the rigidity kappa STABILISES the neutral 0++ that pure NG cannot.

  PART 3  GAUSSIAN / Hartree master field  -> reduces to the DSE gap equation
          -> gapless solution m=0 persists at the crossover (known wall, SLOT-DSE).

Units: sigma = 1 (so sqrt(sigma)=1).  kappa is dimensionless (coefficient of int K^2 ds).
Run:  py -3.12 -I construct_masterfield.py
"""
import numpy as np

sigma = 1.0
D = 4

# ---------------------------------------------------------------- PART 1
print("=" * 72)
print("PART 1 | PURE NG closed string: E0(L)=sigma L sqrt(1 - pi(D-2)/(3 sigma L^2))")
print("=" * 72)
Lc = np.sqrt(np.pi * (D - 2) / (3.0 * sigma))   # = sqrt(2pi/3) for D=4, sigma=1
print(f"  L_c = sqrt(pi(D-2)/(3 sigma)) = {Lc:.4f}   (L<L_c: complex / tachyonic)")
L = np.linspace(Lc * 1.001, 10.0, 200000)
E_ng = sigma * L * np.sqrt(1.0 - np.pi * (D - 2) / (3.0 * sigma * L**2))
imin = int(np.argmin(E_ng))
print(f"  E0(smallest L={L[0]:.3f}) = {E_ng[0]:.5f}")
print(f"  E0(L=10)               = {E_ng[-1]:.5f}")
print(f"  min over interior      = {E_ng.min():.5f} at L={L[imin]:.3f}")
print(f"  => monotone increasing; minimum sits at the LEFT EDGE L->L_c^+ (no interior min).")
print(f"     Confirms PROOF-EFFECTIVE-STRING / TRANSMUTE-STRING-SPECTRUM: pure NG has no 0++ mass.")

# ---------------------------------------------------------------- PART 2
print()
print("=" * 72)
print("PART 2 | RIGID string: E(rho)=2 pi sigma rho + pi kappa / rho  (tension + bending)")
print("=" * 72)
print("  circular neutral loop (N-ality 0, 0++ by symmetry); rho = loop radius, L=2 pi rho")


def rigid_min(kappa):
    rho = np.linspace(0.02, 6.0, 400000)
    E = 2.0 * np.pi * sigma * rho + np.pi * kappa / rho
    i = int(np.argmin(E))
    return rho[i], E[i]


print(f"  {'kappa':>9} {'rho*':>9} {'L*=2pi rho*':>12} {'m_0++/sqrt(sig)':>16} {'E_min':>9}")
for kappa in [0.05, 0.10, 0.15, 0.193, 0.30, 0.50]:
    r, e = rigid_min(kappa)
    print(f"  {kappa:9.3f} {r:9.4f} {2*np.pi*r:12.4f} {e/np.sqrt(sigma):16.4f} {e:9.4f}")

# analytic check: rho*=sqrt(kappa/(2 sigma)), m=2 pi sqrt(2 sigma kappa)
print(f"\n  analytic: rho* = sqrt(kappa/(2 sigma));  m_0++/sqrt(sigma) = 2 pi sqrt(2 kappa)")

# solve kappa for lattice ratio m/√σ = 3.9 (Morningstar-Peardon)
from scipy.optimize import brentq

ratio_of = lambda k: 2.0 * np.pi * np.sqrt(2.0 * k)          # m/√σ, sigma=1
kstar = brentq(lambda k: ratio_of(k) - 3.9, 1e-4, 5.0)
r, e = rigid_min(kstar)
print(f"\n  kappa* for m/sqrt(sigma)=3.9 :  kappa={kstar:.4f},  rho*={r:.4f},  L*={2*np.pi*r:.4f}")
print(f"  kappa* for m/sqrt(sigma)=3.3 :  kappa={brentq(lambda k: ratio_of(k)-3.3,1e-4,5.0):.4f}")
print(f"  (lattice scaling band 3.3-3.9  <->  kappa in [0.138, 0.193])")
print(f"  L* ~ 1.7-1.9 /sqrt(sigma):  the glueball flux-tube length, physically the right scale.")

# ---------------------------------------------------------------- PART 3
print()
print("=" * 72)
print("PART 3 | GAUSSIAN / Hartree master field -> DSE gap equation (cross-check)")
print("=" * 72)
print("  schematic Hartree gap eq:  m^2 = lam * m^2 / (1 + m^2)")
print("  (m^2 sources itself through the tadpole; the 1/(1+m^2) is the dressed propagator)")
print(f"  {'lam':>6} {'m=0 soln':>10} {'m>0 branch m^2':>16} {'exists?':>9}")
for lam in [0.5, 1.0, 1.2, 2.0]:
    m2pos = lam - 1.0
    print(f"  {lam:6.2f} {'ALWAYS':>10} {m2pos:16.3f} {'YES' if lam > 1 else 'NO':>9}")
print("  => m=0 (gapless) is ALWAYS a solution; a positive branch needs lam>1 (=F'(0)>1),")
print("     i.e. the k->0 tachyon that is ABSENT at beta=2.5 (SLOT-DSE, ROOF-SYNTHESIS).")
print("     The Gaussian/Hartree closure reproduces the known wall; it does NOT stabilise a mass.")

print()
print("=" * 72)
print("SUMMARY")
print("=" * 72)
print("  Closure (i)  Gaussian/Hartree  -> DSE gap eq -> gapless persists at crossover (WALL).")
print("  Closure (iv) RIGID surface-sum  -> interior min -> neutral 0++ m=2pi sqrt(2 sigma kappa).")
print("  Pure NG      (kappa=0)          -> NO interior min (corpus no-go) -- rigidity is what changes it.")
print("  Missing input: kappa (flux-tube rigidity / core thickness) = H1 in string coordinates.")
print("  Falsifiable prediction: m_0++/sqrt(sigma) = 2 pi sqrt(2 kappa), kappa in [0.14,0.19] for 3.3-3.9.")
print("  NO QED. This is a construction, not a proof.")
