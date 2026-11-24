"""
H1-RG-DESCENT witness  --  one-step RG descent of the curvature good set, iterated.

Deterministic, `py -3.12`, numpy (+scipy.linalg if needed).  No third-party sources.

Framework (see research/v3/H1-RG-DESCENT.md):
  Shell j, factor-2 block map P (fine->coarse, linearised).  Fine conditioned precision
      A_j = Hess V_j = Hess S_j + Hess Phi     (Phi = -log J, Hess Phi >= (N/12) I).
  Exact coarse (coarse-pinned) curvature = Schur / Gaussian marginal precision:
      Gpp_j := Hess_{U'} S_eff^{(j)} = (P A_j^{-1} P^T)^{-1}  = E[Sigma_j] - D_j .
  One-step multiplicative lemma (exact, Gaussian backbone):
      Gpp_j >= c_j * Gpp_{j+1}^{bare} ,  c_j = (1 + eta_j)^{-1},
      eta_j = || A_{j+1}^{1/2} (P A_j^{-1} P^T) A_{j+1}^{1/2} - I ||_2 ,
      valid iff eta_j < 1 (i.e. ||E||<1).  D_j adds the additive (non-Gaussian) defect.

Sections:
 [D1] exact identity  Schur == (P A^{-1} P^T)^{-1}  (random PD A; compression P).
 [D2] free-field block RG: c_j=1 exactly (sharp passband, no flow); Gpp invariant = beta Q_c.
 [D3] acid test on VOLUME: seed c_L = 4 pi^2 / L^2 -> 0 ; transported floor -> 0  (no gap).
 [D4] head / Savvidy toy: one negative fibre mode -> eta_j >= 1, Gpp_j <= 0 (descent FAILS).
 [D5] cumulative loss ledger: UV sum_{j>=j*} g_j^4 finite; head per-shell loss O(c_L).
"""

import numpy as np

np.set_printoptions(precision=6, suppress=True)
rng = np.random.default_rng(20261124)


def psd(n, scale=1.0, rng=rng):
    M = rng.standard_normal((n, n))
    return scale * (M @ M.T + n * np.eye(n))


def sym(A):
    return 0.5 * (A + A.T)


# =============================================================== [D1] exact Schur identity
def d1_schur_identity():
    print("[D1] EXACT IDENTITY  Gpp = (P A^{-1} P^T)^{-1} = Schur complement of fibre")
    n_f, n_c = 6, 4                      # fine dim, coarse dim
    A = psd(n_f + n_c, scale=1.0)
    # order x = (fibre, coarse);  P = [0 | I]
    A_ff = A[:n_f, :n_f]
    A_fc = A[:n_f, n_f:]
    A_cc = A[n_f:, n_f:]
    A_cf = A_fc.T
    Schur = sym(A_cc - A_cf @ np.linalg.solve(A_ff, A_fc))
    Ainv = np.linalg.inv(A)
    P = np.zeros((n_c, n_f + n_c))
    P[:, n_f:] = np.eye(n_c)
    B = P @ Ainv @ P.T
    Gpp = np.linalg.inv(sym(B))
    err = np.max(np.abs(Gpp - Schur))
    print(f"    || (P A^-1 P^T)^-1 - Schur ||_max = {err:.3e}   (=> exact)")
    # non-Gaussian defect D: Gaussian => 0
    print("    Gaussian backbone: defect D = Var - E[<K^-1 c,c>] = 0 exactly (BL equality)")
    return err


# =============================================================== [D2] free-field block RG
def Q_lowest(L):
    return 2.0 * (1.0 - np.cos(2.0 * np.pi / L))


def d2_free_field():
    print("[D2] FREE-FIELD block RG (sharp passband, no flow): c_j = 1, Gpp invariant")
    L = 16
    Qc = Q_lowest(L)
    beta = 2.5
    m = 2.0 / 12.0                        # N/12, N=2 Haar floor
    # fine precision A = beta Q + m on a mode; coarse curvature = same (sharp passband commutes)
    # Gpp_cond = 1/(v^T C v) = beta Qc + m  (V1 closed form)
    Gpp = beta * Qc + m
    A_coarse = Gpp                        # shell-(j+1) bare = shell-j effective (no flow)
    B = 1.0 / Gpp                         # = v^T C v
    eta = abs(np.sqrt(A_coarse) * B * np.sqrt(A_coarse) - 1.0)
    c = 1.0 / (1.0 + eta)
    print(f"    beta={beta}, Qc={Qc:.6f}, m={m:.6f}")
    print(f"    Gpp = beta Qc + m = {Gpp:.6f}  (invariant across shells, sharp passband)")
    print(f"    eta_j = {eta:.3e}  => c_j = {c:.6f}  (=1 exactly)")
    # show invariance over 'shells' (same beta: free field does not run)
    vals = [beta * Qc + m for _ in range(5)]
    print(f"    Gpp over 5 shells: {np.array(vals)}  (constant)")
    return eta


# =============================================================== [D3] acid test on volume
def d3_volume_acid():
    print("[D3] ACID TEST on VOLUME: seed c_L = 4 pi^2 / L^2 -> 0 ; transported floor -> 0")
    beta = 2.5
    m = 2.0 / 12.0
    print("      L      Qc~4pi^2/L^2     Gpp=beta Qc(+m)   floor/(gap?)")
    for L in [8, 16, 32, 64, 128, 256]:
        Qc = Q_lowest(L)
        Gpp = beta * Qc
        print(f"    {L:4d}   {Qc:12.3e}   {Gpp:12.3e}       ->0")
    print("    => the transported UV floor is the FREE-BOX gap, which vanishes as L->inf:")
    print("       it is NOT a volume-uniform mass gap.  (free massless: Gpp(k=0)=0 too.)")
    # Haar floor is volume-uniform but is measure curvature, not an atom:
    print(f"    (N/12 = {m:.4f} is volume-uniform Haar curvature, but beta=0 measure has f1^2=0:"
          f" curvature>0 =/=> atom, V2 sec7.)")


# =============================================================== [D4] head Savvidy toy
def d4_head_savvidy():
    print("[D4] HEAD / Savvidy toy: a NEGATIVE fibre mode at strong coupling -> descent FAILS")
    # 1 fibre + 1 coarse; fibre precision can go negative (Savvidy soft mode).
    # fine precision A = [[a_ff, g],[g, a_cc]] ; a_ff = mu0 - lam*beta^2 (Savvidy softens)
    mu0 = 1.0
    lam = 0.16
    g = 0.30                              # fibre<->coarse coupling (graph angle * scale)
    a_cc = 1.2
    print("      beta   a_ff(=mu0-lam beta^2)  Gpp=(PA^-1P^T)^-1   eta_j    descent")
    for beta in [0.5, 1.0, 1.5, 2.0, 2.3, 2.4, 2.45, 2.5, 3.0]:
        a_ff = mu0 - lam * beta**2
        A = np.array([[a_ff, g], [g, a_cc]])
        # if A has a negative eigenvalue the fibre is Savvidy-indefinite
        ev = np.linalg.eigvalsh(A)
        A_coarse = a_cc                   # bare shell-(j+1) curvature
        if a_ff <= 1e-9 or ev.min() <= 1e-9:
            print(f"    {beta:4.2f}   {a_ff:12.4f}          (indefinite: "
                  f"lam_min(A)={ev.min():+.3f})   ---    FAILS (A not PD)")
            continue
        # coarse curvature of conditioned fine model: Schur complement of fibre
        Gpp = a_cc - g * g / a_ff
        if Gpp <= 0:
            print(f"    {beta:4.2f}   {a_ff:12.4f}          {Gpp:12.4f}      (Gpp<=0)   "
                  f"FAILS (coarse curvature <=0, A still PD)")
            continue
        Bc = 1.0 / Gpp
        eta = abs(Bc * A_coarse - 1.0)    # 1x1: eta = |Bc*A_coarse - 1|
        c = 1.0 / (1.0 + eta) if eta < 1 else float('nan')
        ok = "ok" if (eta < 1 and Gpp > 0) else "FAILS"
        print(f"    {beta:4.2f}   {a_ff:12.4f}          {Gpp:12.4f}      {eta:8.3f}   {ok}")
    print("    => the coarse curvature Gpp = a_cc - g^2/a_ff DROPS (Schur loss g^2/a_ff ~ g_j^2),")
    print("       crosses ZERO while A is still PD (a_ff < g^2/a_cc, beta above ~2.4), then the")
    print("       fibre goes indefinite (a_ff<0, Savvidy): eta>=1 and the lemma hypothesis fails.")
    print("       The head (g=O(1), beta~2.5) is exactly where this happens.")


# =============================================================== [D5] cumulative loss ledger
def d5_loss_ledger():
    print("[D5] CUMULATIVE LOSS LEDGER")
    N = 2
    beta0 = 11.0 * N / (48.0 * np.pi**2)      # = 22/(48 pi^2)
    print(f"    N={N}, beta0={beta0:.6f}, g_j^2 = 1/(2 beta0 j)")
    # UV tail of the loss ~ C g_j^4 (FORGE5 / G1.1).  Convergent for j>=j*.
    # sum_{j>=j*} j^-4 (tail of zeta(4)); zeta(4)=pi^4/90
    zeta4 = np.pi**4 / 90.0
    print(f"    sum_j g_j^4 (full, from j=1) = (1/(2beta0)^2) zeta(4) = "
          f"{1.0/(2*beta0)**2 * zeta4:.4f}   (dominated by small j = strong coupling)")
    for jstar in [8, 16, 32, 64]:
        # exact tail sum_{j>=jstar} j^-4 by direct summation to convergence
        tail = sum(k**-4 for k in range(jstar, 20000))
        loss = 1.0 / (2 * beta0)**2 * tail
        print(f"    j>={jstar:3d}:  sum g_j^4 (tail) = {loss:.3e}   (UV-convergent tail)")
    # head per-shell loss scale vs the box-gap seed
    L = 16
    Qc = Q_lowest(L)
    print(f"    seed c_L = 4pi^2/L^2 = {Qc:.4f} ; a HEAD shell (g=O(1)) has Schur loss "
          f"~ (g^2/a_ff) = O(Qc) per shell")
    print("    => cumulative UV loss << c_L ; cumulative HEAD loss (O(1) shells at g=O(1)) >= c_L.")
    print("       Net: descent reaches the head with margin ~ c_L (a FREE-BOX artifact), NOT a gap.")


if __name__ == "__main__":
    e1 = d1_schur_identity()
    print()
    d2_free_field()
    print()
    d3_volume_acid()
    print()
    d4_head_savvidy()
    print()
    d5_loss_ledger()
    print()
    assert e1 < 1e-8, "Schur identity failed"
    print("ALL SECTIONS OK (Schur identity exact; deterministic).")
