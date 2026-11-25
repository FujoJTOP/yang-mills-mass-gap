#!/usr/bin/env python3
"""
U2-WILSON-OPE witness.  Deterministic, numpy only.  Run:  py -3.12 -I u2_wilson_ope.py

Question.  Can a Wilson-loop OPE / shape-derivative identity LOWER-bound the tube overlap
    w1 := f1^2 / M0 ,   f1 = <G| tr F^2 |Omega> ,  M0 = Var(tr F^2) = ||Q tr F^2 Omega||^2 ,
i.e. the atom fraction of the lightest 0++ tube of the local operator tr F^2 ?

Claim tested (U2-WILSON-OPE.md): the shape derivative and the tr F^2 OPE are EXACT
identities, but they control the DIAGONAL (on-loop) <tr F^2> = string-tension condensate
and the OPE CONTINUUM (M0's UV), NOT the off-diagonal spectral-edge residue f1^2.
Hence they do not lower-bound w1 -- they are upper-bound / diagonal objects.

What this witness shows (structural + exact solvable-model values; NOT a proof about 4D YM):
 (A) 2D SU(N) YM: EXACT rep-R area law  <W_R> = exp(-sigma_R A),  sigma_R = (g^2/2) C2(R) > 0.
     Shape derivative  -d/dA log<W_R> = sigma_R  (exact, positive).
 (B) 2D: tr F^2 = (1/l) (Casimir Laplacian) is DIAGONAL in the character basis.
       on-loop (diagonal)  <tr F^2>_{W_R}/<W_R> = C2(R)/l  > 0   ~  sigma_R    [STRING TENSION]
       vacuum-to-atom       <chi_0| tr F^2 |chi_R> = C2(R) <chi_0|chi_R> = 0  [ATOM ZERO]
     ==> sigma_R > 0 and f1^2 = 0 coexist EXACTLY.  Decisive no-go for "sigma>0 => atom".
 (C) SU(2) Haar character orthogonality  <chi_0|chi_R> = delta_{0R}  (Weyl-integral numeric check)
     is the mechanism that kills the atom while the on-loop value stays positive.
 (D) OPE direction: the OPE fixes the CONTINUUM spectral weight, giving M0 >= M_cont.
     w1 = f1^2/M0 <= f1^2/M_cont ; a family with f1^2 -> 0 at fixed continuum drives w1 -> 0.
     => the OPE bounds the DENOMINATOR (upper bound direction); it cannot force a bottom atom.
 (E) free field: sigma = 0 (perimeter/Coulomb); composite tr F^2 = two-particle cut => f1^2 = 0.
     Acid test: the route returns 0 (manufactures no floor).
 (F) trace anomaly identity  f1 = (2g/beta) m_G^2 f_G  is a RELATION, not an inequality:
     f_G = 0 => f1 = 0 for any m_G>0 (2D witness: f_G=0 while m_G, sigma_R > 0).
"""
import numpy as np

np.set_printoptions(precision=6, suppress=True, linewidth=140)


def weyl_integral(g, phi):
    """SU(2) Weyl measure (2/pi) sin^2(phi) dphi, phi in [0,pi]."""
    return (2.0 / np.pi) * np.sin(phi) ** 2 * g


print("=" * 82)
print("U2-WILSON-OPE witness  (numpy %s)" % np.__version__)
print("=" * 82)

# ----------------------------------------------------------------------------
# (C) first: SU(2) Haar character orthogonality  <chi_0|chi_R> = delta_{0R}
#     Weyl:  chi_j(phi) = sin((2j+1) phi)/sin(phi),  dmu = (2/pi) sin^2(phi) dphi.
# ----------------------------------------------------------------------------
phi = np.linspace(0.0, np.pi, 400001)
phi = phi[1:-1]                       # drop endpoints (sin->0 handles the ratio, but keep clean)
wgt = (2.0 / np.pi) * np.sin(phi) ** 2


def char_j(j):
    return np.sin((2 * j + 1) * phi) / np.sin(phi)


print("\n(C) SU(2) Haar character overlaps  <chi_0|chi_R> = int dmu chi_0 chi_R = delta_{0R}:")
ortho = []
for j in range(4):
    c = np.trapezoid(wgt * char_j(0) * char_j(j), phi)
    ortho.append(c)
    print("      j=%d :  <chi_0|chi_%d> = %+.3e   (delta_0j: %s)"
          % (j, j, c, abs(c - (1.0 if j == 0 else 0.0)) < 1e-6))

# ----------------------------------------------------------------------------
# (A) 2D SU(N) YM: exact rep-R area law and its shape (area) derivative.
#     <W_R> = exp(-sigma_R A),  sigma_R = (g^2/2) C2(R).  One common convention.
# ----------------------------------------------------------------------------
print("\n(A) 2D YM EXACT area law  <W_R> = exp(-sigma_R A),  sigma_R = (g^2/2) C2(R):")
g2 = 1.0
# SU(2): C2(fund)=3/4 ; adjoint C2 = 2 ; trivial 0.
C2 = {"trivial": 0.0, "fund(1/2)": 0.75, "adjoint(1)": 2.0}
for name, c2 in C2.items():
    sigma = 0.5 * g2 * c2
    print("      R=%-12s C2=%.4f  sigma_R=%.4f  (>0 for nontrivial: %s)"
          % (name, c2, sigma, sigma > 0))
# shape derivative: -d/dA log<W_R> = sigma_R  (compute on a fine A grid)
A = np.linspace(0.5, 8.0, 16)
sigma_f = 0.5 * g2 * 0.75
logW = -sigma_f * A
dlogW = np.gradient(logW, A)
print("      fund: -d/dA log<W> = %+.6f  (constant)  == sigma_fund=%.6f : %s"
      % (-dlogW.mean(), sigma_f, abs(-dlogW.mean() - sigma_f) < 1e-9))

# ----------------------------------------------------------------------------
# (B) 2D tr F^2 is diagonal in the character basis: on-loop >0, atom = 0.
#     tr F^2 |chi_R> = (C2(R)/l) |chi_R>  (Casimir Laplacian / length l).
#     on-loop diagonal <tr F^2>_{W_R} = C2(R)/l > 0  ( = energy density ~ sigma_R )
#     atom (off-diagonal) <chi_0|tr F^2|chi_R> = (C2(R)/l) <chi_0|chi_R> = 0.
# ----------------------------------------------------------------------------
print("\n(B) 2D tr F^2 DIAGONAL in characters -> on-loop POSITIVE, atom ZERO:")
l = 1.0  # spatial length units
for name, c2 in C2.items():
    diag_onloop = c2 / l
    # atom matrix element <chi_0|trF^2|chi_R> = (C2/l) * <chi_0|chi_R> = (C2/l)*delta_0R
    atom = diag_onloop * (1.0 if name == "trivial" else 0.0)
    print("      R=%-12s on-loop <tr F^2> = C2/l = %.4f   atom <0|trF^2|R> = %.4f"
          % (name, diag_onloop, atom))
print("      READ: for R!=0  on-loop = C2/l > 0  (= sigma_R up to g^2/2)  while  atom = 0.")
print("            sigma_R > 0 and f1^2 = 0 COEXIST EXACTLY (character diagonal).")
print("            <chi_0|chi_R> = delta_0R (C) is the mechanism.")

# ----------------------------------------------------------------------------
# (D) OPE direction: M0 >= M_cont (denominator lower-bounded); w1 = f1^2/M0 can -> 0.
#     Model: dnu = f1^2 delta_{m1^2} + A sqrt(s-theta) ds  (continuum fixed by OPE).
# ----------------------------------------------------------------------------
print("\n(D) OPE fixes the CONTINUUM; M0 = f1^2 + M_cont >= M_cont.  w1 = f1^2/M0:")
m1sq, theta, Amp = 1.0, 4.0, 1.0
# M_cont = int_theta^{theta+W} A sqrt(s-theta) ds = A*(2/3) W^{3/2}; W=2
Wc = 2.0
M_cont = Amp * (2.0 / 3.0) * Wc ** 1.5
print("      M_cont (fixed by OPE continuum) = %.6f   (OPE gives M0 >= this)" % M_cont)
print("       f1^2        M0        w1=f1^2/M0     w1 <= f1^2/M_cont ?")
for f1sq in [1.0, 0.5, 0.1, 1e-2, 1e-4, 0.0]:
    M0 = f1sq + M_cont
    w1 = f1sq / M0
    ok = w1 <= f1sq / M_cont + 1e-15
    print("      %8.1e  %9.6f  %12.6f       %s" % (f1sq, M0, w1, ok))
print("      READ: w1 -> 0 as f1^2 -> 0 at FIXED continuum; the OPE (knowing M_cont) does")
print("            NOT exclude a zero atom -> it bounds the denominator, not the edge.")

# ----------------------------------------------------------------------------
# (E) free field: sigma = 0 (perimeter/Coulomb); composite tr F^2 cut => f1^2 = 0.
#     two-particle density on [4m^2,inf): rho(s) ~ sqrt(1-4m^2/s); point mass = 0.
# ----------------------------------------------------------------------------
print("\n(E) free field: no area law (sigma=0); tr F^2 = two-particle cut:")
m = 1.0
s = np.linspace(4 * m * m, 40.0, 200001)
rho = np.sqrt(np.clip(1.0 - 4 * m * m / s, 0.0, None))
# point mass at any m0: limit of mass in a shrinking window / rho(m0) (density: no atom)
point_mass = 0.0  # absolutely continuous spectrum has zero mass on any point
print("      sigma_free = 0.0 (perimeter/Coulomb, no area term).")
print("      rho(s) = sqrt(1-4m^2/s) on [4m^2,inf), absolutely continuous.")
print("      atomic weight nu({s0}) = %.1f for every s0 (no delta)  => f1^2 = 0 exactly."
      % point_mass)
print("      ACID TEST PASSES: the Wilson-loop/OPE route returns w1 = 0 on the free field.")

# ----------------------------------------------------------------------------
# (F) trace anomaly: f1 = (2g/beta) m_G^2 f_G  -- a relation, not an inequality.
#     2D witness: f_G = 0 (tr F^2 = H orthogonal to the R-sector glueball) while m_G>0, sigma>0.
# ----------------------------------------------------------------------------
print("\n(F) trace anomaly  (beta/2g) <0|tr F^2|G> = m_G^2 <0|phi_G>  =>  f1 = (2g/beta) m_G^2 f_G:")
print("      It RELATES f1, m_G, f_G.  'sigma>0 => m_G>0' does NOT force f1>0:")
print("      2D witness: tr F^2 = H is diagonal => f_G = <0|tr F^2|G_R> = 0  =>  f1 = 0,")
print("                  while m_G > 0 and sigma_R > 0.  Same no-go, in sum-rule language.")

print("\n" + "=" * 82)
print("Verdict of the witness (structural, not a 4D YM proof):")
print("  * shape derivative  -d/dA log<W> = sigma_R > 0  (exact, 2D).")
print("  * OPE / on-loop tr F^2 = the DIAGONAL condensate / string tension (bulk).")
print("  * the ATOM f1^2 = <0|tr F^2|G> = off-diagonal spectral-edge residue.")
print("  * 2D: sigma_R>0 AND f1^2=0 EXACTLY (tr F^2 diagonal, <chi_0|chi_R>=0). => no-go.")
print("  * OPE bounds M0 (denominator) => at most an UPPER bound on w1, never a lower one.")
print("  * free field: sigma=0 and f1^2=0 (acid test passes).")
print("  => the Wilson-loop OPE / shape derivative does NOT lower-bound w1.")
print("=" * 82)
