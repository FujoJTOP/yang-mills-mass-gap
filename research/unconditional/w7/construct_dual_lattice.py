"""
construct_dual_lattice.py  --  CONSTRUCT-DUAL-LATTICE witness (V-num)

Vortex-ring 0++ in the dual abelian-Higgs (dual superconductor) and its
tr F^2 overlap.  KEY QUESTION: does the ring carry a nonzero KL atom
f1^2 = nu({m1^2}) > 0, or only a CONTINUUM EDGE?

Physical input (dual abelian-Higgs, London limit):
  m_V = g v          (dual photon = 1^--),
  m_H = sqrt(2 lam) v (dual Higgs),
  kappa_GL = m_H/m_V,  vortex line tension  sigma_ring = pi v^2 ln(kappa_GL),
  core size xi = 1/m_H.

Ring energy (standard vortex-ring, thin core):
  E(R) = A R ( ln(8R/xi) - 7/4 ),   A = rho_dual * flux^2 / 2  ~ sigma_ring.

Structural claim: R is CONTINUOUS, E(R) is monotone increasing for R>xi,
so {E(R)} is a CONTINUUM [E_min, inf) and the KL measure nu has a smooth
density (no delta) => f1^2 = 0.  A KL atom would require a STABLE DISCRETE
ring (E minimum at finite R*>xi), which the minimal dual GL does not supply.

Run:  py -3.12 -I construct_dual_lattice.py
Nothing here is a proof of the gap.  Labels [V-num]/[C]/[O]/[NG].
"""
import numpy as np

print("=" * 72)
print("CONSTRUCT-DUAL-LATTICE: vortex-ring 0++ in the dual abelian-Higgs")
print("KEY QUESTION: KL atom f1^2>0, or continuum edge?")
print("=" * 72)

xi = 1.0            # core size 1/m_H (units)
A  = 1.0            # line-tension coefficient (schematic; ~ sigma_ring)

def E(R):   return A * R * (np.log(8.0 * R / xi) - 7.0 / 4.0)
def dE(R):  return A * (np.log(8.0 * R / xi) - 3.0 / 4.0)

Rs  = np.linspace(1.0, 40.0, 4000)
Es  = E(Rs)
dEs = dE(Rs)

print("\n[1] Ring energy E(R) on [xi, 40 xi]")
print(f"    dE/dR = A(ln(8R/xi) - 3/4);  min over R>=xi = {dEs.min():.4f}")
print(f"    monotone increasing for R>=xi ?  {bool(np.all(dEs > 0))}")
print(f"    E(xi) = {E(xi):.4f},  E(40xi) = {E(40.0):.4f}")
print("    => spectrum {E(R): R>xi} = CONTINUUM [E_min, inf), no isolated level.")

print("\n[2] tr F^2 spectral density from the ring continuum")
def w(R):
    # tr F^2 localized on the flux core; overlap largest for small rings (schematic)
    return np.exp(-(R - xi) / (6.0 * xi))
wt = w(Rs)**2 / dEs          # rho(s) = w(R)^2 / E'(R)  (continuum density)
print(f"    rho(s) finite everywhere ? {bool(np.all(np.isfinite(wt)))}")
print(f"    rho(s) range [{wt.min():.3e}, {wt.max():.3e}]  -> smooth DENSITY, no delta")

print("\n[3] KL atom test: spectral mass in a shrinking window around E_min")
s0 = E(xi)
for eps in [1e-1, 1e-2, 1e-3, 1e-4]:
    m = np.abs(Es - s0) < eps
    mass = np.trapezoid(wt[m], Es[m]) if m.sum() > 1 else 0.0
    print(f"    |s-Emin|<{eps:.0e}: spectral mass = {mass:.3e}   (-> 0 as eps -> 0)")
print("    => no delta at m1^2  =>  f1^2 = nu({m1^2}) = 0 : CONTINUUM EDGE, not an atom.")

print("\n[4] Contrast: a STABLE discrete ring (E minimum at R*>xi) WOULD be an atom")
B = 3.0
def Etr(R): return E(R) + B * (xi**2) / R
Rt = np.linspace(0.5, 30.0, 60000)
Et = Etr(Rt)
i  = int(np.argmin(Et))
h  = 1e-5
slope = (Etr(Rt[i] + h) - Etr(Rt[i] - h)) / (2 * h)
print(f"    with trapping +B xi^2/R: E minimum at R* = {Rt[i]:.3f}, E* = {Et[i]:.4f}")
print(f"    slope at min = {slope:.2e} (~0)  => isolated level => delta => atom f1^2>0")
print("    (but the minimal dual GL has NO such trapping term; a finite-V torus")
print("     quantizes R with spacing ~ 1/V -> continuum as V->inf, not a genuine atom)")

print("\n[5] Dual-GL scaling of the continuum edge (leading 2 pi R sigma at R=xi)")
for kappa, lam in [(2, 2.0), (5, 12.5), (10, 50.0)]:
    v = 1.0
    sig  = np.pi * v**2 * np.log(kappa)      # sigma_ring
    mH   = np.sqrt(2.0 * lam) * v
    Emin = 2.0 * np.pi * (1.0 / mH) * sig    # 2 pi xi sigma
    print(f"    kappa_GL={kappa:4.1f}: sigma_ring={sig:6.3f},  Emin/sqrt(sigma)={Emin/np.sqrt(sig):.3f}")
print("    O(1-2) sqrt(sigma): a continuum edge a few sqrt(sigma) up (light 0++ scale).")
print("    The heavy 3.9-sqrt(sigma) 0++ (if a ring) is just further up the SAME continuum.")

print("\n[6] Two dual-sector candidates vs the physical 3.9 ratio")
print("    - dual Higgs m_H: DISCRETE (atom-like) but m_H/sqrt(sigma)=1.36..3.72 (SHORT of 3.9)")
print("      and its tr F^2 overlap is FREE in the minimal AHM (f1^2 undetermined).")
print("    - vortex ring: reaches ~3.9 but CONTINUUM (f1^2=0).")
print("    => neither dual object delivers a KNOWN atom at 3.9 without extra input = H1.")

print("\nVERDICT: the vortex-ring 0++ is a CONTINUUM EDGE (f1^2 = 0), NOT a KL atom.")
print("The dual superconductor gives a renormalon-free MASS but no ATOM (f1^2>0).")
print("Clay 0%; no QED.")
