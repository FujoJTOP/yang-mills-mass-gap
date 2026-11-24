"""
blast_exclude_bump.py -- witness for research/unconditional/BLAST-EXCLUDE-BUMP.md
Can a YM-specific property EXCLUDE the atomless flat bump (same mass, floor, width h->0)?
This script checks the ONLY-EXIT adversary's moments and the candidates for exclusion.
Run: py -3.12 -I research/unconditional/w7/blast_exclude_bump.py
"""
import numpy as np

m1sq = 1.0
f2 = 0.5
cont_lo, cont_hi, cont_mass = 4.0, 6.0, 1.0 - f2

def mom_atom():
    # atom f2 at m1sq + uniform continuum [cont_lo,cont_hi]
    M0 = f2 + cont_mass
    chi = f2 / m1sq + cont_mass * np.log(cont_hi / cont_lo) / (cont_hi - cont_lo)
    A = f2 / m1sq**2 + cont_mass * (1/cont_lo - 1/cont_hi) / (cont_hi - cont_lo)
    return M0, chi, A

def mom_bump(h, n=200000):
    # uniform bump of total mass f2 on [m1sq, m1sq+h] + same continuum
    s = np.linspace(m1sq, m1sq + h, n)
    w = (f2 / h) * np.gradient(s)  # density f2/h
    M0b = np.sum(w)
    chib = np.sum(w / s)
    Ab = np.sum(w / s**2)
    # continuum
    sc = np.linspace(cont_lo, cont_hi, n)
    wc = (cont_mass / (cont_hi - cont_lo)) * np.gradient(sc)
    M0c = np.sum(wc); chic = np.sum(wc / sc); Ac = np.sum(wc / sc**2)
    return M0b + M0c, chib + chic, Ab + Ac

Ma, Xa, Aa = mom_atom()
print("atom :  M0=%.8f  chi=%.8f  A=%.8f   f1^2=%.4f" % (Ma, Xa, Aa, f2))
print("squeeze atom: f1^2/m1^2=%.4f <= chi=%.4f <= M0/m1^2=%.4f" % (f2/m1sq, Xa, Ma/m1sq))
print()
for h in [1e-1, 1e-2, 1e-3, 1e-4]:
    Mb, Xb, Ab = mom_bump(h)
    dM0 = abs(Mb - Ma); dchi = abs(Xb - Xa); dA = abs(Ab - Aa)
    # bump's atom weight at the floor is ZERO by construction
    print("h=%.0e  |dM0|=%.2e  |dchi|=%.2e  |dA|=%.2e   f1^2_bump=0.0" % (h, dM0, dchi, dA))
print()
print("=> moments (matching M0, chi, A = all FINITE integral data) converge as O(h);")
print("   the bump has f1^2=0 while the atom has f1^2=0.5, yet ALL finite moments agree.")
print()
# Is any positive nu a valid 0++ KL measure? YES -- generalized free / :phi^2: yields
# an atomless gapped 0++ measure (f1^2=0). So OS-positivity alone does NOT exclude the bump.
print("Cross-channel candidate: theta = (beta/2g) tr F^2 ties 0++ to tr T_mu^mu; 2++ (traceless T)")
print("is a SEPARATE measure; the dilatation Ward relates their spectral densities -> needs the")
print("2++ spectrum / gap = H1. No finite 0++ property separates the bump from the atom.")
