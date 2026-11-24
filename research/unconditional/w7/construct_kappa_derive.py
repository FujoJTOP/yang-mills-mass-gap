"""
construct_kappa_derive.py -- witness for CONSTRUCT-KAPPA-DERIVE.

Question: can the rigid-string rigidity kappa (in m_{0++}/sqrt(sigma)=2*pi*sqrt(2*kappa))
be DERIVED from (a) the flux-tube width / Luscher term, (b) PK-EFT RG, (c) microscopic string?

Pieces (cheap, numpy/scipy):
 (1) Flux-tube width w^2(R) for a Gaussian string with rigidity:
       w^2(R) = D_perp/(2*pi) * int_{1/R}^{1/a} dk / (sigma k + kappa k^3)
     Fit the ln(R) slope for kappa=0 and kappa>0.  Claim: slope is UNIVERSAL
     (D_perp/(2*pi*sigma)), independent of kappa; kappa only shifts the additive constant
     (= the effective core radius r0).  => the Luscher/width route CANNOT fix kappa.
 (2) Circularity: invert m_{0++}=2*pi*sqrt(2*kappa) => kappa = m^2/(8*pi^2*sigma).
     So "derive kappa" == "derive m/sqrt(sigma)" == the gap (not a softer handle).
 (3) Dimensional analysis: kappa = c * sigma * xi^2 (dimensionless), xi=core radius.
     At the glueball scale xi ~ 1/sqrt(sigma) => kappa ~ c (O(1)), consistent with the
     fitted band [0.138,0.193] but c is the microscopic constant (nonperturbative).

No claim about rate of convergence / no theorem.  Labels inline.
Run: py -3.12 -I construct_kappa_derive.py
"""
import numpy as np
from scipy.integrate import quad

# ---------- (1) width and the universal log slope ----------
def width2(R, sigma, kappa, Dperp=2, a=0.04):
    """Coincident transverse Green function of the rigid worldsheet.
       G(k)=1/(sigma k^2 + kappa k^4); w^2 = Dperp ∫d^2k/(2π)^2 G(k)."""
    lo, hi = 1.0 / R, 1.0 / a
    val, _ = quad(lambda k: 1.0 / (sigma * k + kappa * k**3), lo, hi, limit=500)
    return Dperp / (2.0 * np.pi) * val

sigma = 1.0
Rs = np.array([40, 80, 160, 320, 640, 1280], float)
print("=" * 70)
print("(1) flux-tube width: is the ln(R) slope kappa-dependent?")
print("=" * 70)
pred = 2 / (2 * np.pi * sigma)   # D_perp/(2*pi*sigma), D_perp=2
print(f"  predicted universal slope D_perp/(2*pi*sigma) = {pred:.6f}\n")
slopes = {}
for kappa in [0.0, 0.138, 0.193]:
    ws = np.array([width2(R, sigma, kappa) for R in Rs])
    slope, intercept = np.polyfit(np.log(Rs), ws, 1)
    slopes[kappa] = slope
    print(f"  kappa={kappa:.3f}: slope={slope:.6f}  intercept={intercept:.6f}")
print(f"\n  slope spread over kappa in {{0,0.138,0.193}} = "
      f"{max(slopes.values())-min(slopes.values()):.2e}  (== 0 up to quadrature error)")
print("  => the Luscher/width log-slope is UNIVERSAL; kappa only moves the additive")
print("     constant (the effective core radius r0).  The width cannot derive kappa.\n")

# show the additive constant (intercept) DOES depend on kappa (it is the core datum)
ints = {}
for kappa in [0.0, 0.138, 0.193]:
    ws = np.array([width2(R, sigma, kappa) for R in Rs])
    _, intercept = np.polyfit(np.log(Rs), ws, 1)
    ints[kappa] = intercept
print(f"  intercept shift kappa=0 -> 0.193: {ints[0.193]-ints[0.0]:+.6f}  "
      f"(= log of the core-radius ratio; a UV/microscopic datum)\n")

# ---------- (2) circularity: kappa <-> m ----------
print("=" * 70)
print("(2) circularity:  m_{0++} = 2*pi*sqrt(2*kappa)  =>  kappa = m^2/(8*pi^2*sigma)")
print("=" * 70)
for ratio in [3.3, 3.9]:
    m = ratio * np.sqrt(sigma)
    k = m**2 / (8 * np.pi**2 * sigma)
    print(f"  m/sqrt(sigma)={ratio}  =>  kappa = {k:.4f}")
for kappa in [0.138, 0.193]:
    m = 2 * np.pi * np.sqrt(2 * kappa)
    print(f"  kappa={kappa:.3f}  =>  m/sqrt(sigma) = {m:.4f}")
print("  => kappa and m/sqrt(sigma) are the SAME unknown re-expressed;")
print("     deriving kappa is deriving the gap (not a softer handle).\n")

# ---------- (3) dimensional analysis ----------
print("=" * 70)
print("(3) dimensional analysis: kappa = c * sigma * xi^2,  xi = core radius")
print("=" * 70)
print("  [sigma]=mass^2, [xi]=1/mass => sigma*xi^2 dimensionless, c = O(1) coefficient.")
print("  At the glueball scale xi ~ 1/sqrt(sigma):  kappa = c * (sigma * 1/sigma) = c.")
print("  Fitted band kappa in [0.138,0.193] => c in [0.138,0.193]: an O(1) number,")
print("  set by the microscopic flux-tube structure (nonperturbative), not by sigma alone.")
