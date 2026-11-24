"""construct_width_kappa_lattice.py -- CONSTRUCT-WIDTH-KAPPA-LATTICE witness.

Question: can the Polyakov-Kleinert rigidity alpha (and/or the MASTERFIELD kappa)
be extracted from the known flux-tube width w^2(R) = (D-2)/(2 pi sigma) ln(R/r0) + w0^2 ?
Honest answer produced below: NO -- the log coefficient is universal (alpha-free),
and the constant w0^2 fixes only the combination alpha*r0^2 (degenerate with the UV
core scale r0).  PK alpha is dimensionful and r0-degenerate; MASTERFIELD kappa is
dimensionless -- no direct map.  numpy only; run: py -3.12 -I
"""
import numpy as np

D = 4
# physical units: sqrt(sigma) ~ 440 MeV ; 1 GeV = 5.068 fm^-1
sqrt_sigma = 0.44 * 5.068          # fm^-1
sigma = sqrt_sigma**2              # fm^-2
print("sqrt(sigma) = %.3f fm^-1 , sigma = %.3f fm^-2" % (sqrt_sigma, sigma))
print("D = %d" % D)
print("Luscher per-direction log coeff 1/(2 pi sigma) = %.4f fm^2" % (1/(2*np.pi*sigma)))
print("(D-2)/(2 pi sigma) (summed)   = %.4f fm^2" % ((D-2)/(2*np.pi*sigma)))
print("(D-2)/(4 pi sigma) (const term)= %.4f fm^2" % ((D-2)/(4*np.pi*sigma)))

# PK constant term:  w0^2 = (D-2)/(4 pi sigma) ln(alpha sigma r0^2)
# => alpha*r0^2 = exp( 4 pi sigma w0^2/(D-2) ) / sigma
print("\n-- invert PK width constant for alpha*r0^2 (representative lattice w0^2) --")
for w0sq in [0.10, 0.15, 0.20]:
    L = 4*np.pi*sigma*w0sq/(D-2)
    alpha_r0sq = np.exp(L)/sigma
    print("w0^2=%.2f fm^2 : ln(alpha sigma r0^2)=%.3f , alpha*r0^2=%.3f fm^2" %
          (w0sq, L, alpha_r0sq))
    for r0 in [0.15, 0.20, 0.30]:
        print("      r0=%.2f fm -> alpha=%.2f  (varies as 1/r0^2 -> NOT determined)" %
              (r0, alpha_r0sq/r0**2))

# MASTERFIELD kappa from m/sqrt(sigma) = 2 pi sqrt(2 kappa)  => kappa=(m/sqrt(sigma))^2/(8 pi^2)
print("\n-- MASTERFIELD kappa from the lattice band (a FIT, not derived) --")
for ratio in [3.3, 3.7, 3.9]:
    kappa = (ratio**2)/(2*(2*np.pi)**2)
    print("m/sqrt(sigma)=%.2f -> kappa=%.4f" % (ratio, kappa))
print("\nPK alpha: dim fm^2, r0-degenerate. MASTERFIELD kappa: dimensionless.")
print("=> no direct map; width data does NOT determine the MASTERFIELD kappa.")
