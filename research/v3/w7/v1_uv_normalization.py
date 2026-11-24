"""
V1-UV-NORMALIZATION witness  --  the free-block limit of (★) at the UV end of the tower.

Deterministic, `py -3.12`, numpy only.  No third-party sources.

Objects (flat chart X_l = log U_l in su(N); Wilson action S = beta_p sum_p (1 - Re tr/N);
Haar dU = J dX, Phi = -log J; measure  exp(-S - Phi) dX):

  H(v)      = v^T Hess S v                       (bare force Hessian, bilinear)
  G(v)      = <v, grad S> = v.grad S             (force)
  Phi'(v)   = <v, grad Phi>,  Phi''(v) = v^T Hess Phi v
  W_geo(v)  = E[Phi''] - Var(Phi')               (Haar geometric curvature)
  Gpp_lit   = E[H] - Var(G) + W_geo              (ch05 eq:split, NO cross term X, boundary-conditioned)
  Gpp_model = same as Gpp_lit  (S2 "2 W_geo" object)
  Gpp_cond  = -d^2/da^2 log rho_v(a)|_{a=0} = 1/(v^T C v)   (coarse-PINNED marginal curvature; G4 (C))
  defect    = E[H] - Var(G)
  Delta_dyn = defect (= W_geo unconditionally by Stein+cross-Stein, G4 Lemma G2)

Free block = quadratic Wilson action around the flat connection (Gaussian fluctuations),
PLUS the exact Haar measure.  Then the fiber is exactly Gaussian in the flat chart:
   C = (beta_p Q + m I)^{-1},  m = N/12  (Hess of the Haar potential at the identity).
For a unit coarse eigenmode v of Q (eigenvalue Qc) all moments are closed form.

Tower (C2 / S5 conventions): g_j^2 = 1/(2 beta0 j), beta0 = 11N/(48 pi^2),
   beta_j = 2N/g_j^2 = 4N beta0 j  -> infinity (weak coupling),  a_j = a0 L^{-j}.
4D periodic lattice, side L sites: lowest nonzero Laplacian eigenvalue
   Qc = 2 (1 - cos(2 pi / L))  ->  4 pi^2 / L^2  as L -> infinity  (== c_L box gap).
"""

import numpy as np

np.set_printoptions(precision=6, suppress=True)


# ----------------------------------------------------------------------------- lattice mode
def Q_coarse(L):
    """Lowest nonzero eigenvalue of Q = 2 sum_mu (1 - cos k_mu) on a periodic 4D lattice, side L.
    Exciting a single direction gives the box gap;  -> 4 pi^2/L^2."""
    return 2.0 * (1.0 - np.cos(2.0 * np.pi / L))


# ----------------------------------------------------------------------------- closed forms
def free_block_closed(beta, Qc, m):
    """Closed-form free-block moments for a unit coarse eigenmode of Q (eigenvalue Qc),
    fiber covariance C=(beta Q + m I)^{-1}.  Returns a dict."""
    H = beta * Qc                              # E[H], since H = beta Q is constant
    Cvv = 1.0 / (beta * Qc + m)                # (C)_{vv}
    VarG = beta**2 * Qc**2 * Cvv               # Var(G) = beta^2 (Q C Q)_{vv}, v eigenvector
    EPhi2 = m                                  # E[Phi''] -> Hess Phi = m I  (per unit mode)
    VarPhi1 = m**2 * Cvv                       # Var(Phi') = m^2 (C)_{vv}
    Wgeo = EPhi2 - VarPhi1
    defect = H - VarG
    Gpp_model = defect + Wgeo
    Gpp_lit = Gpp_model                        # X = -2Cov, Cov = Wgeo (G4 G2) -> lit = model + X
    Gpp_lit = Gpp_model - 2.0 * Wgeo           # = 0 exactly shewn below
    Gpp_cond = 1.0 / Cvv                        # marginal curvature of a = <v,xi> at 0
    return dict(H=H, VarG=VarG, defect=defect, Wgeo=Wgeo,
                Gpp_model=Gpp_model, Gpp_lit=Gpp_lit,
                Gpp_cond=Gpp_cond, Cvv=Cvv, VarPhi1=VarPhi1, EPhi2=EPhi2)


def tower(N, jmax, L):
    beta0 = 11.0 * N / (48.0 * np.pi**2)
    Qc = Q_coarse(L)
    m = N / 12.0
    rows = []
    for j in range(1, jmax + 1):
        g2 = 1.0 / (2.0 * beta0 * j)
        beta = 2.0 * N / g2
        d = free_block_closed(beta, Qc, m)
        d.update(j=j, beta=beta, g2=g2, Qc=Qc, m=m)
        rows.append(d)
    return rows


def fmt_table(rows, N, L):
    print(f"\n=== FREE BLOCK, N={N}, 4D periodic L={L} (Qc={rows[0]['Qc']:.6f}), m=N/12={rows[0]['m']:.6f} ===")
    hdr = ("   j     beta_j     g_j^2     E[H]    Var(G)   defect=E[H]-Var    W_geo"
           "    Gpp_model   Gpp_lit   Gpp_cond   Gpp_cond/beta")
    print(hdr)
    for d in rows:
        print(f"{d['j']:4d} {d['beta']:9.4f} {d['g2']:9.5f} {d['H']:9.4f} {d['VarG']:9.4f}"
              f" {d['defect']:14.6f} {d['Wgeo']:9.6f} {d['Gpp_model']:10.6f} {d['Gpp_lit']:9.2e}"
              f" {d['Gpp_cond']:10.4f} {d['Gpp_cond']/d['beta']:11.6f}")


# ----------------------------------------------------------------------------- explicit dense check
def dense_check(N, beta, Qc, m, seed=0):
    """Build an explicit small Gaussian block to confirm every closed form and identity.
    Use an orthonormal basis {v, u1..u_{d-1}} with v the coarse eigenvector of Q (eig Qc),
    the rest higher eigenvalues Qh, and compare E[H], Var(G), W_geo, G''_lit=0, G''_model=2W_geo,
    G''_cond=1/(v^T C v)."""
    rng = np.random.default_rng(seed)
    d = 6
    # random orthogonal basis, v = e0 (the coarse eigenvector)
    A = rng.standard_normal((d, d))
    Qb, _ = np.linalg.qr(A)
    e0 = np.zeros(d); e0[0] = 1.0
    # Q in the basis: eigenvalues [Qc, Qh*(1..d-1)]  (higher modes)
    Qh = np.linspace(Qc * 4.0, Qc * 12.0, d - 1)
    evals = np.concatenate([[Qc], Qh])
    Q = Qb @ np.diag(evals) @ Qb.T
    Q = 0.5 * (Q + Q.T)
    H = beta * Q
    C = np.linalg.inv(beta * Q + m * np.eye(d))
    v = e0  # unit mode = coarse eigenvector

    def Emat(M):       # E[v^T M v] under N(0,C)
        return v @ M @ v
    # E[H_vv] = beta (v,Qv) (H constant)
    EH = v @ H @ v
    # Var(G_v) = Var(v^T H xi) = v^T H C H v
    VarG = v @ H @ C @ H @ v
    # Phi = (m/2)|xi|^2 -> Phi' = m v.xi, Phi'' = m |v|^2 = m
    EPhi2 = m * (v @ v)
    VarPhi1 = m**2 * (v @ C @ v)
    Wgeo = EPhi2 - VarPhi1
    # Cov(G,Phi') = Cov(v^T H xi, m v.xi) = m v^T H C v
    CovGPhi = m * (v @ H @ C @ v)
    Gpp_model = EH - VarG + Wgeo
    Gpp_lit = EH - VarG + Wgeo - 2.0 * CovGPhi
    Gpp_cond = 1.0 / (v @ C @ v)
    return dict(EH=EH, VarG=VarG, Wgeo=Wgeo, Cov=CovGPhi, Gpp_model=Gpp_model,
                Gpp_lit=Gpp_lit, Gpp_cond=Gpp_cond,
                twoW=2.0 * Wgeo, defect=EH - VarG)


# ----------------------------------------------------------------------------- main
def main():
    print("=" * 108)
    print("V1-UV-NORMALIZATION : free-block limit of (★) at the UV end (beta_j -> infinity)")
    print("=" * 108)

    # --- Part A: closed-form scaling tables
    fmt_table(tower(N=2, jmax=40, L=16), N=2, L=16)

    rows = tower(N=2, jmax=40, L=16)
    d1, dN = rows[0], rows[-1]

    print("\n--- limits (N=2) ---")
    print(f"  defect  = E[H]-Var  : {d1['defect']:.6f} (j=1)  ->  {dN['defect']:.6f} (j=40)  "
          f"=> N/12 = {2/12:.6f}  [O(beta) cancels, Haar floor survives]")
    print(f"  W_geo               : {d1['Wgeo']:.6f}         ->  {dN['Wgeo']:.6f}        "
          f"=> N/12 = {2/12:.6f}  (coupling-free)")
    print(f"  Gpp_model=2 W_geo   : {d1['Gpp_model']:.6f}         ->  {dN['Gpp_model']:.6f}        "
          f"=> N/6 = {2/6:.6f}   (O(1), j-uniform, in LATTICE units)")
    print(f"  Gpp_lit             : {d1['Gpp_lit']:.2e}       ->  {dN['Gpp_lit']:.2e}        "
          f"=> 0  (X=-2Cov cancels the model)  [G4 G2]")
    print(f"  Gpp_cond            : {d1['Gpp_cond']:.4f}       ->  {dN['Gpp_cond']:.4f}        "
          f"=> beta_j * Qc  (O(beta_j), grows)")
    print(f"  Gpp_cond/beta_j     : {d1['Gpp_cond']/d1['beta']:.6f}      ->  "
          f"{dN['Gpp_cond']/dN['beta']:.6f}       => Qc = {d1['Qc']:.6f}  (O(1), j-uniform positive)")

    # --- Part B: dense explicit check of every identity
    print("\n=== dense explicit Gaussian block (d=6), confirm identities ===")
    for beta, tagn in [(2.5, "head-ish"), (40.0, "UV"), (400.0, "deep UV")]:
        r = dense_check(N=2, beta=beta, Qc=Q_coarse(16), m=2/12, seed=1)
        print(f"  beta={beta:7.2f} ({tagn:9s}): "
              f"|lit|={abs(r['Gpp_lit']):.2e}  |model-2W|={abs(r['Gpp_model']-r['twoW']):.2e}  "
              f"|defect-W|={abs(r['defect']-r['Wgeo']):.2e}  |lit_full|={abs(r['Gpp_lit']):.2e}  "
              f"cond={r['Gpp_cond']:.4f}")

    # --- Part C: acid test, pure free field Phi=0
    print("\n=== acid test: pure free field (Phi=0, abelian) ===")
    beta, Qc, m = 50.0, Q_coarse(16), 0.0
    d = free_block_closed(beta, Qc, m)
    print(f"  beta={beta}:  W_geo={d['Wgeo']:.3e}  Gpp_model={d['Gpp_model']:.3e}  "
          f"Gpp_lit={d['Gpp_lit']:.3e}  defect={d['defect']:.3e}   => G''=0  ((★) fails, delta_W=0)")

    # --- Part D: scaling law confirmation
    print("\n=== scaling law: bare delta_W must track beta_j  (NOT a power of a_j) ===")
    print("  a_j = a0 L^{-j};  beta_j = 4N beta0 j = 2N/g_j^2.  j = log_L(a0/a_j), so")
    print("     beta_j  ~  g_j^{-2}  ~  log(1/a_j)        (asymptotic freedom: LOGARITHMIC in a_j)")
    print(f"  N=2,L=16:  Gpp_cond(j=40)/Gpp_cond(j=1) = {dN['Gpp_cond']/d1['Gpp_cond']:.4f} ; "
          f"beta_40/beta_1 = {dN['beta']/d1['beta']:.4f}  (ratio = 1 => delta_W ∝ beta_j)")
    print(f"  Z_j^{-1} delta_W  ->  Qc = {d1['Qc']:.6f} = c_L (free-box gap), O(1), j-uniform.")
    print("=" * 108)


if __name__ == "__main__":
    main()
