"""
roof_stochastic.py -- witness for ROOF-STOCHASTIC.md
Cheap and exact. numpy only. Run: py -3.12 -I research/unconditional/w7/roof_stochastic.py

(A) Free-field Fokker-Planck gap == mass: one-particle energies omega_k = sqrt(khat^2 + m^2);
    the gap is the lowest omega over k. For the free field this IS the exact FP spectrum, so
    gap == m (mass), and gap -> 0 as m -> 0  (the acid: the massless field is gapless).
(B) Coupling dimension [g^2] = 4 - d: super-renormalizable for d<4, MARGINAL/critical at d=4.
    Regularity structures need subcriticality (d<4); d=4 sits on the critical line.
(C) Perturbative coefficients grow like beta0^n n! (IR renormalon): the 4D counterterm series is
    asymptotic/non-summable -- the same u=2 obstruction as the condensate floor.

No Yang-Mills computation here; these are exact structural facts. Labels inline.
"""
import numpy as np
import math

def khat2(n, L):
    # lattice momentum squared for k = 2*pi*n/L
    k = 2.0 * np.pi * n / L
    return 4.0 * np.sin(k / 2.0) ** 2

print("=" * 74)
print("(A) free-field FP gap == mass  (one-particle spectrum omega_k = sqrt(khat^2 + m^2))")
print("=" * 74)
L = 64
ms = [0.0, 0.01, 0.1, 0.5, 1.0]
ns = np.arange(L)
kk = np.array([khat2(n, L) for n in ns])
for m in ms:
    omega = np.sqrt(kk + m * m)
    gap = omega.min()
    print(f"  m = {m:5.2f}   min_k omega_k = {gap:.6f}   (= m ? {abs(gap-m) < 1e-12})")
print("  => gap(FP) = m exactly; m -> 0 gives gap -> 0 (acid: massless is gapless).")
print("  => a rigorous gap lower bound for FP IS the mass gap (the identity is a tautology).\n")

print("=" * 74)
print("(B) coupling dimension [g^2] = 4 - d  (super-renorm vs MARGINAL/critical)")
print("=" * 74)
for d in [2, 3, 4]:
    dim_g2 = 4 - d
    cls = "super-renormalizable" if dim_g2 > 0 else ("MARGINAL (critical)" if dim_g2 == 0 else "non-renormalizable")
    ok = "RS applies (finitely many counterterms)" if dim_g2 > 0 else "RS NOT applicable (all-orders counterterms)"
    print(f"  d = {d}:  [g^2] = {dim_g2:+d}   {cls:22s}  -> {ok}")
print("  => 4D Yang-Mills is the CRITICAL dimension: g^2 dimensionless, asymptotically free.\n")

print("=" * 74)
print("(C) perturbative counterterm series is asymptotic (IR renormalon at u = d/2 = 2)")
print("=" * 74)
beta0 = 11.0 / (48.0 * np.pi**2)   # SU(2), n_f = 0, schematic normalization
print("  a_n = beta0^n * n!   (factorial growth; ratio a_n/a_{n-1} = beta0 * n -> infinity)")
prev = 1.0
for n in range(1, 9):
    a_n = beta0**n * math.factorial(n)
    print(f"    n = {n}:  a_n = {a_n:.3e}   ratio a_n/a_(n-1) = {a_n/prev:8.4f}   (n*beta0 = {n*beta0:8.4f})")
    prev = a_n
print("  => non-summable (Borel pole at u=2); the 4D SPDE counterterms = the renormalon = H1.")
print("\nVerdict: gap(FP)=mass [T]; d<4 RS supplies the strut [T]; d=4 MARGINAL -> renormalon wall [NG]/[O].")
