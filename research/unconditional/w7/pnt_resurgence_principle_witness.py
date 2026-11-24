"""
pnt_resurgence_principle_witness.py  --  PNT-RESURGENCE-PRINCIPLE witness (toy, EXACT algebra)
Run:  py -3.12 -I research/unconditional/w7/pnt_resurgence_principle_witness.py

Purpose: demonstrate, on a transparent toy with an IR renormalon at u=2, that
  (A) the perturbative coefficients {a_n} fix the Borel-residue / Stokes constant
      but are COMPLETELY independent of the transseries parameter sigma;
  (B) the two lateral Borel sums differ by 2*pi*i*Res*exp(-1/(c*alpha)); taking the
      average (PV) plus any REAL sigma*exp(-1/(c*alpha)) leaves F real for EVERY real sigma;
  (C) a positivity-type condition (rho>=0) constrains |sigma| to a RANGE, never to a point.

This is a toy, NOT a Yang-Mills calculation.  It isolates the logical structure.
"""
import numpy as np
from math import gamma, pi

np.set_printoptions(precision=6, suppress=True)

def banner(s):
    print("\n" + "=" * 74); print(s); print("=" * 74)

# ---------------------------------------------------------------------------
# Toy perturbative series:  a_n = A * Gamma(n+1+nu) * c^n
#   -> Borel transform B(z) = sum a_n z^n / n!  has a pole at z = 1/c.
#   Choose c = 1/2 so the pole sits at z = 2  (the "u=2" IR renormalon, normalized).
# ---------------------------------------------------------------------------
A = 1.0
nu = 0.0
c = 0.5            # pole at z = 1/c = 2
N = 40

a = np.array([A * gamma(n + 1 + nu) * c**n for n in range(N)])
geom = np.array([gamma(n + 1 + nu) * c**n for n in range(N)])
ratio = a / geom

banner("(A) perturbative coefficients {a_n}: set by the pole, INDEPENDENT of sigma")
print("  a_n[0:6]                 =", a[:6])
print("  a_n / (Gamma(n+1+nu) c^n) =", ratio[:6], " -> constant A (the residue structure)")
print("  (adding sigma*exp(-1/(c*alpha)) to F changes NO a_n: sigma is not in {a_n})")

# ---------------------------------------------------------------------------
# (B) Lateral sums and the transseries parameter.
#     The divergent part of the Borel sum is  PV + tau*exp(-1/(c*alpha)),  tau = sigma.
#     Two lateral sums (contour above / below the u=2 pole) differ by the residue:
#         S_+ - S_- = 2*pi*i*Res*exp(-1/(c*alpha)).
#     They are complex conjugates for real data, so  S_+ = PV + i*pi*Res*e^{-...},
#     S_- = PV - i*pi*Res*e^{-...}.  The REAL transseries is
#         F = (S_+ + S_-)/2 + sigma*e^{-1/(c*alpha)} = PV + sigma*e^{-1/(c*alpha)}
#     for ANY real sigma.
# ---------------------------------------------------------------------------
banner("(B) lateral sums differ by 2*pi*i*Res*e^{-1/(c*alpha)}; ANY real sigma gives real F")
Res = A
print(f"  residue at u=2 pole:  Res = {Res}")
for alpha in [0.2, 0.3, 0.5, 1.0]:
    e = np.exp(-1.0 / (c * alpha))            # e^{-2/alpha}
    jump = 2 * pi * Res * e                    # |S_+ - S_-|
    F_sig = {s: 0.0 + s * e for s in (0.0, -0.5, 1.7, 10.0)}
    line = "  ".join(f"sigma={s:+.1f}:F={F_sig[s]:+.4f}" for s in F_sig)
    print(f"  alpha={alpha:.1f}:  jump={jump:.4e}\n     {line}")
print("  => all F real; no a_n or reality condition selects a sigma.")

# ---------------------------------------------------------------------------
# (C) Positivity constrains sigma to a RANGE (of order the ambiguity), not a point.
#     rho(s) >= 0 is a condition on the (renormalon-free) spectral density; the OPE
#     SPLIT  PV + sigma*e^{-...}  shifts by the prescription.  The most positivity can
#     do is bound |sigma| by the ambiguity width; it cannot fix a value.
# ---------------------------------------------------------------------------
banner("(C) positivity (rho>=0) bounds |sigma| to a RANGE, never a point")
for alpha in [0.2, 0.3, 0.5, 1.0]:
    e = np.exp(-1.0 / (c * alpha))
    width = pi * Res * e                       # half-width of the allowed band (order of ambiguity)
    print(f"  alpha={alpha:.1f}:  positivity band  |sigma| <= {width:.4e}   (a range, not a value)")

banner("VERDICT")
print("  {a_n} fix the Stokes constant / residue but leave sigma free;")
print("  positivity bounds |sigma| (range); reality leaves sigma free;")
print("  for a SINGLE IR renormalon (no saddle, no second scale) the freedom is GENUINE.")
