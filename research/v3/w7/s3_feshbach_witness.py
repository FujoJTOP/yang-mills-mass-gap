#!/usr/bin/env python3
"""
s3_feshbach_witness.py -- witness for S3-UV-FESHBACH.md (the two-block/Feshbach
decomposition of the UV shell inequality (star): G''_j >= delta_W).

WHAT THIS IS / IS NOT
---------------------
An *illustrative, deterministic* numerical witness, NOT a proof.  It checks, on finite
deterministic models, the EXACT algebra of the Schur-variance identity
        G'' = E[Sigma] - D ,   Sigma = A - <K^{-1} c, c>  (Schur complement),
        D = Var(f | a=0) - E[<K^{-1} c, c> | a=0]         (Brascamp-Lieb defect),
and the two structural facts the note turns on:

  (1) GAUSSIAN: the coarse marginal curvature G'' equals the Schur complement Sigma
      exactly (D = 0), and Sigma = A - B C^{-1} B^T  <=  A  (eliminating a POSITIVE fibre
      HURTS the coarse curvature).  Sign of the Schur correction when C > 0.

  (2) SIGN FLIP (Lemma NS): if the eliminated complement is NEGATIVE (C < 0, the Savvidy
      LLL, the D-BBS / FORGE-FESHBACH convention), then Sigma = A + B|C|^{-1}B^T  >=  A:
      the coupling HELPS.  This is the *opposite* convention from TCI/(star).

  (3) NON-GAUSSIAN, 2D by quadrature: the identity G'' = E[Sigma] - D holds exactly, and
      D <= 0 whenever the fibre K >= 0 (log-concave fibre; Prékopa + Brascamp-Lieb 1976),
      hence G'' >= E[Sigma].  This is the UV reduction (star) <= [fibre gap]+[Schur bound].

  (4) ACID TEST: the massless free field returns G'' = 0 and Sigma = 0 -- the identity
      manufactures no false positivity (it must NOT be >0 on the free field).

  (5) W_j ROLE: the action-only coarse curvature E[Sigma_S] can fail to clear delta, while
      adding a convex Jacobian Phi'' = N/12 > 0 (the PROOF-JACOBIAN-SAVES floor) raises
      E[Sigma_V] above delta.  W_j's provable content is exactly this Phi'' floor.

Conventions: the flat-chart measure mu ~ exp(-V) dX, V = S + Phi (Phi = -log J Haar),
coarse mode a = <P^T X, v>, condition a = 0; A = v^T Hess V v, K = fibre Hessian,
c = coupling, f = force.  These are the objects of INVENT-VARIANCE-INEQ / PROOF-JACOBIAN-SAVES.

Run:  python s3_feshbach_witness.py
"""

import numpy as np

np.set_printoptions(precision=6, suppress=True)


def block_matrix(A, B, C):
    n, m = A.shape[0], C.shape[0]
    M = np.zeros((n + m, n + m))
    M[:n, :n] = A
    M[:n, n:] = B
    M[n:, :n] = B.T
    M[n:, n:] = C
    return M


# ----------------------------------------------------------------------------------
# (1) Gaussian:  G'' == Schur complement;  Sigma = A - B C^{-1} B^T <= A   (C > 0)
# ----------------------------------------------------------------------------------
def gaussian_schur():
    rng = np.random.default_rng(0)
    d = 2
    X = rng.standard_normal((d, d))
    A = X @ X.T + 1.5 * np.eye(d)          # SPD coarse block
    X = rng.standard_normal((d, d))
    C = X @ X.T + 1.5 * np.eye(d)          # SPD fibre (positive!) block
    B = 0.7 * rng.standard_normal((d, d))

    M = block_matrix(A, B, C)
    # coarse marginal covariance = (M^{-1})[:d,:d]; coarse curvature = its inverse
    Minv = np.linalg.inv(M)
    Gpp = np.linalg.inv(Minv[:d, :d])      # = (Sigma) exactly for a Gaussian
    Sigma = A - B @ np.linalg.inv(C) @ B.T
    err = np.max(np.abs(Gpp - Sigma))
    # Sigma - A = -B C^{-1} B^T  is NEGATIVE semidefinite  (coupling hurts)
    diff = Sigma - A
    lam_diff = np.linalg.eigvalsh(diff)
    print("(1) GAUSSIAN, C > 0 (fibre positive):")
    print(f"    || G'' - Schur Sigma ||_max        = {err:.3e}   (should be ~0)")
    print(f"    eig(Sigma - A)                     = {np.round(lam_diff,6)}")
    print(f"    Sigma <= A ?  (all eig(Sigma-A)<=0) = {bool(np.all(lam_diff <= 1e-10))}")
    print(f"    => eliminating a POSITIVE fibre HURTS: Sigma = A - BC^-1B^T <= A\n")


# ----------------------------------------------------------------------------------
# (2) Sign flip (Lemma NS):  C < 0  =>  Sigma = A + B|C|^{-1}B^T >= A   (coupling helps)
# ----------------------------------------------------------------------------------
def lemma_NS():
    rng = np.random.default_rng(1)
    d = 2
    X = rng.standard_normal((d, d))
    A = X @ X.T + 1.0 * np.eye(d)
    X = rng.standard_normal((d, d))
    C = -(X @ X.T + 1.0 * np.eye(d))       # NEGATIVE definite (Savvidy LLL toy)
    B = 0.7 * rng.standard_normal((d, d))

    Sigma = A - B @ np.linalg.inv(C) @ B.T          # = A + B|C|^{-1}B^T
    diff = Sigma - A
    lam_diff = np.linalg.eigvalsh(diff)
    print("(2) SIGN FLIP, C < 0 (eliminated block negative; D-BBS / FORGE-FESHBACH):")
    print(f"    eig(Sigma - A)                     = {np.round(lam_diff,6)}")
    print(f"    Sigma >= A ?  (all eig(Sigma-A)>=0) = {bool(np.all(lam_diff >= -1e-10))}")
    print(f"    => eliminating a NEGATIVE block HELPS (Lemma NS): Sigma >= A.")
    print(f"    (TCI/(star) use the OPPOSITE convention, C>0 => Sigma<=A.)\n")


# ----------------------------------------------------------------------------------
# 2D non-Gaussian potential by direct quadrature.
#   V(x,y) = 1/2 a x^2 + 1/4 lam y^4 + 1/2 c y^2 + 1/2 m (x-y)^2   [+ 1/2 phi (x^2+y^2)]
#   at x=0:  h_n = a+m+phi, f = -m y, K = 3 lam y^2 + c + m + phi, c_coupling = -m
#   G'' = E[h_n] - Var(f);  E[Sigma] = E[h_n - c_c^2/K];  D = Var(f) - E[c^2/K]
# ----------------------------------------------------------------------------------
def fibre_stats(a, lam, c, m, phi, Y=9.0, Ny=400001):
    y = np.linspace(-Y, Y, Ny)
    V0 = 0.25 * lam * y**4 + 0.5 * c * y**2 + 0.5 * m * y**2 + 0.5 * phi * y**2
    V0 -= V0.min()
    w = np.exp(-V0)
    Z = np.trapezoid(w, y)
    pdf = w / Z
    h_n = a + m + phi                       # d^2/dx^2 V at x=0 (constant)
    f = -m * y                              # d/dx V at x=0
    K = 3.0 * lam * y**2 + c + m + phi      # d^2/dy^2 V
    cc = -m                                 # d^2/dxdy V
    Ef = np.trapezoid(pdf * f, y)
    Ef2 = np.trapezoid(pdf * f * f, y)
    Var_f = Ef2 - Ef**2
    # Schur correction E[cc^2 / K]; if the coupling vanishes it is 0 even when K=0 (free field).
    if abs(cc) < 1e-15:
        E_schur = 0.0
    else:
        E_schur = cc**2 * np.trapezoid(pdf / K, y)
    E_Sigma = h_n - E_schur
    Gpp = h_n - Var_f
    D = Var_f - E_schur
    return dict(Gpp=Gpp, E_Sigma=E_Sigma, D=D, Kmin=K.min(), Var_f=Var_f,
                E_schur=E_schur, h_n=h_n)


def non_gaussian_identity():
    print("(3) NON-GAUSSIAN 2D by quadrature (log-concave fibre, K>0):")
    a, lam, c, m = 0.05, 0.10, 0.05, 1.0
    s = fibre_stats(a, lam, c, m, phi=0.0)
    print(f"    params a={a} lam={lam} c={c} m={m}")
    print(f"    K_min = {s['Kmin']:.4f} (>0 log-concave)")
    print(f"    G''      = {s['Gpp']:.8f}")
    print(f"    E[Sigma]= {s['E_Sigma']:.8f}")
    print(f"    D        = {s['D']:.3e}   (<=0 ? {s['D'] <= 1e-9})")
    ident = s['Gpp'] - (s['E_Sigma'] - s['D'])
    print(f"    identity G'' - (E[Sigma]-D) = {ident:.2e}   (should be ~0)")
    print(f"    G'' >= E[Sigma] (since D<=0) = {s['Gpp'] >= s['E_Sigma'] - 1e-9}\n")

    # concave-fibre contrast: lam<0, K changes sign somewhere -> D can be > 0
    a2, lam2, c2, m2 = 0.05, -0.5, 1.2, 0.3
    s2 = fibre_stats(a2, lam2, c2, m2, phi=0.0)
    print(f"    contrast (concave fibre, lam={lam2}, K_min={s2['Kmin']:.3f}):")
    print(f"      D = {s2['D']:.3e}  (>0 possible when fibre not log-concave; here "
          f"{'>0' if s2['D'] > 0 else '<=0'})")
    print(f"      G''- (E[Sigma]-D) = {s2['Gpp'] - (s2['E_Sigma'] - s2['D']):.2e}\n")


# ----------------------------------------------------------------------------------
# (4) Acid test: massless free field -> G'' = 0, Sigma = 0.
# ----------------------------------------------------------------------------------
def acid_free():
    print("(4) ACID TEST, massless free field (a=m=lam=c=phi=0):")
    s = fibre_stats(0.0, 0.0, 0.0, 0.0, 0.0, Y=30.0, Ny=600001)
    print(f"    G'' = {s['Gpp']:.3e}   E[Sigma] = {s['E_Sigma']:.3e}   D = {s['D']:.3e}")
    print(f"    NO false positivity: G''==0 and Sigma==0 (max |.|={max(abs(s['Gpp']),abs(s['E_Sigma'])):.2e})\n")


# ----------------------------------------------------------------------------------
# (5) W_j role: action-only E[Sigma_S] fails delta; adding convex Jacobian phi=N/12 fixes it.
# ----------------------------------------------------------------------------------
def W_role():
    print("(5) W_j ROLE: Jacobian floor Phi'' = N/12 raises E[Sigma] above delta.")
    a, lam, c, m = 0.01, 0.10, 0.01, 5.0
    delta = 0.10
    c_phi = 1.0 / 6.0            # N/12 for SU(2) (PROOF-JACOBIAN-SAVES Lemma B)
    s0 = fibre_stats(a, lam, c, m, phi=0.0)
    s1 = fibre_stats(a, lam, c, m, phi=c_phi)
    print(f"    params a={a} lam={lam} c={c} m={m} (strong coupling B), delta={delta}")
    print(f"    action only  (phi=0):      E[Sigma_S] = {s0['E_Sigma']:.4f}   "
          f"G'' = {s0['Gpp']:.4f}   clears delta? {s0['E_Sigma'] >= delta}")
    print(f"    + Jacobian   (phi=N/12):   E[Sigma_V] = {s1['E_Sigma']:.4f}   "
          f"G'' = {s1['Gpp']:.4f}   clears delta? {s1['E_Sigma'] >= delta}")
    print(f"    => W_j must supply the Phi''>0 floor to pay the Schur loss <K^-1 c,c>.\n")


def main():
    print("=" * 82)
    print("S3 witness: two-block/Schur decomposition of the UV shell inequality (star)")
    print("=" * 82 + "\n")
    gaussian_schur()
    lemma_NS()
    non_gaussian_identity()
    acid_free()
    W_role()
    print("=" * 82)
    print("READING:")
    print(" * (1),(2): the Schur correction's SIGN is set by the eliminated block: C>0 subtracts")
    print("   (hurts; TCI/(star) convention), C<0 adds (helps; Lemma NS / D-BBS convention).")
    print(" * (3): the exact identity G''=E[Sigma]-D holds; log-concave fibre => D<=0 => G''>=E[Sigma].")
    print(" * (4): free massless field -> G''=Sigma=0 (no manufactured positivity).")
    print(" * (5): the Jacobian Phi''>=N/12 floor is the provable content of W_j; it is what")
    print("   lifts E[Sigma] over delta_W.  Bounding the Schur loss itself remains OPEN.")
    print(" NOT a proof.  Witness only.")
    print("=" * 82)


if __name__ == "__main__":
    main()
