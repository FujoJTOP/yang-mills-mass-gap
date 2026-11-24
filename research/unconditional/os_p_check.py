"""
OS-P free-field check -- exact, tiny, single process (numpy only).

Tests the three structural claims about the spatial coarse projection P_V = E[.|V]
(V = spatial block-average on a fixed time slice) on the free Gaussian scalar
on a periodic L_t x L_s lattice with site (through-site) time reflection t -> -t:

  (A) E[.|A_0] (conditional expectation onto the time-0 slice) is OS-ORTHOGONAL:
      the residual phi_1 - E[phi_1|phi_0] has zero OS inner product with phi_0.
      This is the free-field witness that P_V descends to H_phys as an OS-orthogonal
      projection (uses reflection symmetry C((-1,i),(0,j)) = C((1,i),(0,j)) ).

  (B) [P_V, T] != 0: the one-step transfer does NOT factor through the coarse
      variable.  Witness: E[V_1 | phi_0] != E[V_1 | V_0] as functions of phi_0,
      i.e. the operator  R = B X C0^{-1} - B X B^T (B C0 B^T)^{-1} B  is nonzero.
      (P_F T P_V != 0 -> coarse sector not autonomous -> P_V and T do not commute.)

  (C) vacuum + soft mode retained: the constant (k=0) mode is block-constant,
      so P_V Omega = Omega and the global soft mode lies in range(P_V).

No untrusted data; pure exact linear algebra.  N <= a few hundred.
"""
import numpy as np
np.set_printoptions(precision=6, suppress=True, linewidth=150)


def build_M(Lt, Ls, m2):
    """Full spacetime Hessian: -Delta_t - Delta_s + m^2, periodic, index (t,x)."""
    N = Lt * Ls

    def idx(t, x):
        return t * Ls + x

    M = np.zeros((N, N))
    for t in range(Lt):
        for x in range(Ls):
            i = idx(t, x)
            M[i, i] += 4.0 + m2
            for dt, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                M[i, idx((t + dt) % Lt, (x + dx) % Ls)] += -1.0
    return M


def block_map(Ls, b):
    assert Ls % b == 0
    n = Ls // b
    B = np.zeros((n, Ls))
    for I in range(n):
        for k in range(b):
            B[I, I * b + k] = 1.0 / b
    return B


def run(Lt, Ls, b, m2, tag):
    print(f"\n================  {tag}:  L_t={Lt} L_s={Ls} b={b} m^2={m2}  ================")
    N = Lt * Ls

    def idx(t, x):
        return t * Ls + x

    M = build_M(Lt, Ls, m2)
    C = np.linalg.inv(M)

    # ---------- (A) OS orthogonality of E[.|A_0] ----------
    # time-0 slice block, and cross t=1 | t=0
    C0 = np.array([[C[idx(0, x), idx(0, xp)] for xp in range(Ls)] for x in range(Ls)])
    X1 = np.array([[C[idx(1, x), idx(0, xp)] for xp in range(Ls)] for x in range(Ls)])
    Xm1 = np.array([[C[idx(Lt - 1, x), idx(0, xp)] for xp in range(Ls)] for x in range(Ls)])  # t=-1
    # reflection symmetry check: C((-1,i),(0,j)) should equal C((1,i),(0,j))
    refl_defect = np.linalg.norm(Xm1 - X1)
    # OS tangent: <phi_{1,i} - (X1 C0^{-1} phi_0)_i , phi_{0,j}>_OS
    #            = C((-1,i),(0,j)) - (X1 C0^{-1} C0)_{ij} = Xm1 - X1
    E1_given_0 = X1 @ np.linalg.inv(C0)          # E[phi_1 | phi_0] = this @ phi_0
    res_tangent = Xm1 - (E1_given_0 @ C0)        # should be 0 (== Xm1 - X1)
    print(f"  (A) reflection defect  ||C(-1,0) - C(1,0)||           = {refl_defect:.3e}")
    print(f"  (A) OS-tangent defect  ||<phi_1-E[phi_1|phi_0], phi_0>_OS|| = "
          f"{np.linalg.norm(res_tangent):.3e}   (0 => E[.|A_0] OS-orthogonal)")

    # ---------- (B) [P_V, T] != 0 : transfer does not factor through coarse ----------
    B = block_map(Ls, b)                          # n x Ls
    BC0Bt = B @ C0 @ B.T                          # coarse marginal covariance (n x n)
    # full predictor of V_1 = B phi_1 given phi_0:  B X1 C0^{-1} phi_0
    P_full = (B @ X1) @ np.linalg.inv(C0)         # n x Ls
    # coarse predictor given only V_0 = B phi_0:  B X1 B^T (B C0 B^T)^{-1} B phi_0
    P_coarse = (B @ X1 @ B.T) @ np.linalg.inv(BC0Bt) @ B   # n x Ls
    R = P_full - P_coarse                         # n x Ls ; nonzero => not autonomous
    print(f"  (B) ||P_full|| = ||B X C0^-1||                        = {np.linalg.norm(P_full):.6f}")
    print(f"  (B) ||R|| = ||E[V_1|phi_0] - E[V_1|V_0]|| operator    = {np.linalg.norm(R):.6f}")
    print(f"  (B) relative defect ||R||/||P_full||                  = "
          f"{np.linalg.norm(R)/np.linalg.norm(P_full):.6f}   (>0 => [P_V,T]!=0)")
    # how much predictability the fiber carries: variance explained lost
    # covariance of V_1 conditioned on phi_0 vs on V_0
    CovV1 = B @ np.array([[C[idx(1, x), idx(1, xp)] for xp in range(Ls)] for x in range(Ls)]) @ B.T
    var_cond_phi0 = CovV1 - P_full @ (C0 @ P_full.T)
    var_cond_V0   = CovV1 - P_coarse @ (C0 @ P_coarse.T)
    extra = np.trace(var_cond_V0) - np.trace(var_cond_phi0)
    print(f"  (B) extra Var(V_1) unexplained by V_0 vs phi_0 (trace) = {extra:.6f}  (>0)")

    # ---------- (C) vacuum + soft mode in coarse ----------
    e0 = np.ones(Ls) / np.sqrt(Ls)                # constant spatial mode at t=0
    # project onto row space of B (block-constant subspace)
    Pc = B.T @ np.linalg.inv(B @ B.T) @ B         # orthonormal proj onto row(B)
    resid_soft = np.linalg.norm(Pc @ e0 - e0)
    # constant 1 is block-constant => trivially in range(P_V); report the soft-mode overlap
    print(f"  (C) ||e0 - P_{'{row B}'} e0||  (0 => k=0 mode in V)       = {resid_soft:.3e}")

    # ---------- (B') time-shift: V_0 and V_1 spaces are transverse ----------
    # canonical correlation between coarse fields at t=0 and t=1
    X01 = B @ np.array([[C[idx(0, x), idx(1, xp)] for xp in range(Ls)] for x in range(Ls)]) @ B.T
    # correlation matrix in the coarse marginal metric
    Wc, Vc = np.linalg.eigh(BC0Bt)
    isqrt = Vc @ np.diag(1.0 / np.sqrt(Wc)) @ Vc.T
    Kcorr = isqrt @ X01 @ isqrt
    kappas = np.sort(np.linalg.eigvalsh(Kcorr))
    print(f"  (B') max canonical correlation |corr(V_0,V_1)|        = {kappas[-1]:.6f}"
          f"   ( =1 => frozen/deterministic; <1 => T moves the coarse space )")

    return dict(refl=refl_defect, resA=np.linalg.norm(res_tangent),
                Rrel=np.linalg.norm(R) / np.linalg.norm(P_full),
                extra=extra, soft=resid_soft, corr=kappas[-1])


if __name__ == "__main__":
    print("############ OS-P FREE-FIELD STRUCTURAL CHECK ############")
    r1 = run(Lt=8,  Ls=8, b=2, m2=1.0,  tag="massive")
    r2 = run(Lt=8,  Ls=8, b=2, m2=0.25, tag="light")
    r3 = run(Lt=10, Ls=12, b=3, m2=0.5, tag="b=3")
    print("\n############ SUMMARY ############")
    for nm, r in (("massive", r1), ("light", r2), ("b=3", r3)):
        print(f"  {nm:8s}: (A) OS-tangent defect={r['resA']:.2e}   "
              f"(B) relative [P_V,T] defect={r['Rrel']:.4f}   corr(V0,V1)={r['corr']:.4f}   "
              f"(C) soft mode in V: {r['soft']:.1e}")
