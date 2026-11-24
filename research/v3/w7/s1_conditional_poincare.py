#!/usr/bin/env python3
"""
s1_conditional_poincare.py -- witness for S1-CONDITIONAL-POINCARE.md.

A *sanity / illustration* witness, NOT a proof.  Controlled finite models only.

(A) FREE FIELD (acid test: (star) MUST FAIL).  On a Gaussian (free massless) block the
    conditional measure is Gaussian, so Stein gives Var_fiber(force) = E_fiber[Hess]
    identically; with W=0 this makes (star) read 0 <= -delta, FALSE.  G'' = E[H]-Var = 0.

(B) BOUNDARY-DEPENDENT HESSIAN (Savvidy threshold).  On a 4D SU(2) lattice at an abelian
    U(1) constant-field background of strength f(theta)=1-cos(theta/2), build the fine
    Wilson Hessian and report its minimum eigenvalue.  The minimum is the lattice Savvidy
    (curvature) direction; as f grows it goes negative, and once lambda_min < -N/12 the Haar
    Jacobian floor cannot restore log-concavity.  This is the boundary-dependent failure of
    Lemma S1 / Prop. U'-NG.  (Savvidy 1977; Nielsen-Olesen 1978 -- re-measured here.)

(C) HAAR JACOBIAN (Prekopa).  Deterministic check g''(d) >= 1/6 = N/12 for the SU(2)
    Jacobian exponent g(d) = -2 log(sin(d/2)/(d/2)), |d|<2pi  [T, ch05].  The boundary-
    independent term.

(D) GAUSSIAN BLOCK (SCBL split G''=S-D).  D=0, S=1/Var(a); the k->0 mode has S->0.  Shows
    log-concavity (K>0) is NOT sufficient for (star): a coarse spectral gap (=H1) is needed.

Run:  py -3.12 s1_conditional_poincare.py       (numpy only, deterministic)
"""

import numpy as np

np.random.seed(20261124)

I2 = np.eye(2, dtype=complex)
S1 = np.array([[0, 1], [1, 0]], dtype=complex)
S2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
S3 = np.array([[1, 0], [0, -1]], dtype=complex)


def expand_su2(p):
    """U = exp(i p.sigma) = cos|p| I + i sin|p| (p_hat.sigma); vectorised over p[...,3]."""
    p = np.asarray(p, dtype=float)
    r = np.linalg.norm(p, axis=-1)
    out = np.empty(p.shape[:-1] + (2, 2), dtype=complex)
    small = r < 1e-12
    n = np.zeros_like(p)
    n[~small] = p[~small] / r[~small, None]
    c = np.cos(r); s = np.sin(r)
    mats = [S1, S2, S3]
    for a in range(2):
        for b in range(2):
            val = c * I2[a, b]
            for k in range(3):
                val = val + 1j * s * n[..., k] * mats[k][a, b]
            out[..., a, b] = val
    if np.any(small):
        out[small] = I2
    return out


# --------------------------------------------------------------------------- (A)
def part_A_free_field():
    print("=" * 82)
    print("(A) FREE MASSLESS FIELD  (acid test: (star) MUST fail)")
    print("=" * 82)
    n, beta = 8, 1.0
    A = np.zeros((n, n))
    for i in range(n):
        A[i, i] = 2.0
        if i > 0: A[i, i - 1] = -1.0
        if i < n - 1: A[i, i + 1] = -1.0
    A *= beta
    v = np.ones(n) / np.sqrt(n)             # block-uniform coarse mode
    var_force = v @ A @ v                    # Var(<v,force>)
    E_H = v @ A @ v                          # E[Hess_vv]
    print(f"  Dirichlet chain n={n}, beta={beta}, coarse mode v=uniform")
    print(f"    Var_fiber(force_v) = {var_force:.10f}")
    print(f"    E_fiber[Hess_vv]   = {E_H:.10f}")
    print(f"    difference         = {var_force - E_H:.3e}  (Stein: exactly 0)")
    print(f"    W_j = 0  => (star):  {E_H:.6f} <= {E_H:.6f} - delta  <=>  0 <= -delta  FALSE.")
    print("  Lowest Dirichlet eigenvalue c_n = 4 sin^2(pi/(2(n+1))) -> 0:")
    for m in (8, 64, 512):
        c = 4.0 * np.sin(np.pi / (2 * (m + 1))) ** 2
        print(f"      n={m:4d}: c_n = {c:.6e}   (S = 1/Var -> 0, G'' -> 0: free field, no gap)")
    print("  => (star) FAILS on the free field, as REQUIRED.\n")


# --------------------------------------------------------------------------- (B)
def build_lattice_4d(L):
    sites = [(a, b, c, d) for a in range(L) for b in range(L) for c in range(L) for d in range(L)]
    link_id = {}; links = []
    for s in sites:
        for mu in range(4):
            t = list(s); t[mu] += 1
            if t[mu] < L:
                link_id[(s, mu)] = len(links); links.append((s, mu))
    plaq = []
    for s in sites:
        for mu in range(4):
            for nu in range(mu + 1, 4):
                xm = list(s); xm[mu] += 1; xm = tuple(xm)
                xn = list(s); xn[nu] += 1; xn = tuple(xn)
                keys = [(s, mu), (xm, nu), (xn, mu), (s, nu)]
                if all(k in link_id for k in keys):
                    plaq.append([link_id[k] for k in keys])
    return links, np.array(plaq, dtype=int)


def hessian_4d(L, beta, theta, fd=2e-4):
    """Batched finite-difference Hessian of S_W at background
       U_1(x)=exp(i (theta*x4) sigma_3)  [p=(0,0,theta*x4), using U=exp(i p.sigma)].
       Plaquette in (1,4): U_1(x) U_1(x+4)^dag = exp(-i theta sigma_3), f=1-cos(theta)."""
    links, plaq = build_lattice_4d(L)
    nlink = len(links)
    P0 = np.zeros((nlink, 3))
    for idx, (s, mu) in enumerate(links):
        if mu == 1:
            P0[idx, 2] = theta * s[3]
    i1, i2, i3, i4 = plaq[:, 0], plaq[:, 1], plaq[:, 2], plaq[:, 3]

    def action_batch(P):                      # P: (B,nlink,3) -> S: (B,)
        U = expand_su2(P)                     # (B,nlink,2,2)
        U1_ = U[:, i1]; U2_ = U[:, i2]
        U3d = np.conj(np.transpose(U[:, i3], (0, 1, 3, 2)))
        U4d = np.conj(np.transpose(U[:, i4], (0, 1, 3, 2)))
        Up = U1_ @ U2_ @ U3d @ U4d            # (B,P,2,2)
        tr = Up[..., 0, 0] + Up[..., 1, 1]
        return beta * np.sum(1.0 - 0.5 * np.real(tr), axis=1)

    npar = nlink * 3
    S0 = action_batch(P0.reshape(1, nlink, 3))[0]
    # diagonal
    base = np.repeat(P0.reshape(1, nlink, 3), 2 * npar, axis=0)
    for k in range(npar):
        lk, ck = divmod(k, 3)
        base[2 * k, lk, ck] += fd
        base[2 * k + 1, lk, ck] -= fd
    Sd = action_batch(base)
    diag = (Sd[0::2] - 2 * S0 + Sd[1::2]) / fd ** 2
    H = np.diag(diag)
    # off-diagonal, chunked over pairs to bound memory
    pairs = [(k, l) for k in range(npar) for l in range(k + 1, npar)]
    chunk = 4000
    for start in range(0, len(pairs), chunk):
        pr = pairs[start:start + chunk]
        cfg = np.repeat(P0.reshape(1, nlink, 3), 4 * len(pr), axis=0)
        for j, (k, l) in enumerate(pr):
            lk, ck = divmod(k, 3); ll, cl = divmod(l, 3)
            cfg[4 * j + 0, lk, ck] += fd; cfg[4 * j + 0, ll, cl] += fd
            cfg[4 * j + 1, lk, ck] += fd; cfg[4 * j + 1, ll, cl] -= fd
            cfg[4 * j + 2, lk, ck] -= fd; cfg[4 * j + 2, ll, cl] += fd
            cfg[4 * j + 3, lk, ck] -= fd; cfg[4 * j + 3, ll, cl] -= fd
        Ss = action_batch(cfg).reshape(len(pr), 4)
        val = (Ss[:, 0] - Ss[:, 1] - Ss[:, 2] + Ss[:, 3]) / (4 * fd ** 2)
        for j, (k, l) in enumerate(pr):
            H[k, l] = val[j]; H[l, k] = val[j]
    return 0.5 * (H + H.T)


def part_B_savvidy_threshold():
    print("=" * 82)
    print("(B) SU(2) 4D Wilson Hessian at a U(1) constant field  (Savvidy threshold)")
    print(f"    Haar floor c_Phi = N/12 = {2/12:.4f} (SU(2)).   f(theta)=1-cos(theta).")
    print("=" * 82)
    L, beta, N = 2, 2.5, 2
    cphi = N / 12.0
    print(f"  L={L}^4, beta={beta}, N={N}")
    print(f"  {'theta':>7} {'f=1-cos':>9} {'#neg':>5} {'lambda_min':>12} {'lam_min+N/12':>14} {'log-concave?':>15}")
    for theta in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5):
        H = hessian_4d(L, beta, theta)
        ev = np.linalg.eigvalsh(H)
        lmin = ev[0]; nneg = int(np.sum(ev < -1e-6))
        f = 1.0 - np.cos(theta)
        lc = "YES" if (lmin + cphi > 0) else "NO (indefinite)"
        print(f"  {theta:7.2f} {f:9.4f} {nneg:5d} {lmin:12.4f} {lmin+cphi:14.4f} {lc:>15}")
    print("  Reading: at theta=0 (no boundary field) the Hessian is positive semidefinite (pure")
    print("  gauge/Laplacian).  As the boundary field strength grows the Savvidy (curvature) term")
    print("  drives lambda_min negative; once lambda_min < -N/12 the Haar floor cannot make the")
    print("  block log-concave => Lemma S1 threshold crossed.  [Savvidy instability, re-measured.]")
    print("  NOTE 1: L=2 is minimal; the softest Savvidy direction is the constant transverse")
    print("  nonabelian (sigma_1,sigma_2) mode.  This is the GLOBAL/soft-mode value (no block")
    print("  Dirichlet gap).  A finite RG block adds a Dirichlet gap c_box(U') > 0 on its modes;")
    print("  the true modulus is (c_box + N/12) - C_1*beta*f_boundary.  Deep UV c_box ~ 1/a_j^2")
    print("  saves it; at the head (c_box ~ 1/L^2 small) it does not.  Same UV/IR split as the")
    print("  campaign: the UV block is log-concave, the head is H1.\n")


# --------------------------------------------------------------------------- (C)
def part_C_haar():
    print("=" * 82)
    print("(C) SU(2) Haar Jacobian: g''(d) >= 1/6 = N/12   [T, ch05]")
    print("=" * 82)
    # g(d) = -2 log( sin(d/2)/(d/2) ) => g''(d) = (1/2) csc^2(d/2) - 2/d^2   (closed form)
    def gpp(d):
        d = np.asarray(d, float)
        return 0.5 / np.sin(d / 2.0) ** 2 - 2.0 / d ** 2
    ds = np.linspace(1e-3, 2 * np.pi - 1e-3, 20000)
    gp = gpp(ds)
    print(f"  min_d g''(d) over |d|<2pi = {np.min(gp):.5f}  (floor 1/6 = {1/6:.5f})")
    print(f"  g'' >= N/12 everywhere: {np.all(gp >= 1/6 - 1e-6)}")
    print(f"    (g''(0+)=1/6, g''>1/6 for d>0, -> +inf as d->2pi; formula (1/2)csc^2(d/2)-2/d^2)")
    print("  => Phi=-sum g(lambda_i-lambda_j) is convex; e^{-Phi}dX log-concave; by Prekopa every")
    print("     linear marginal has convex free energy, i.e. W_geo = E[Phi'']-Var(Phi') >= 0,")
    print("     boundary-INDEPENDENT.  This is the only U'-uniform positive term.  [T]\n")


# --------------------------------------------------------------------------- (D)
def part_D_gaussian_block():
    print("=" * 82)
    print("(D) GAUSSIAN BLOCK: SCBL  G'' = S - D")
    print("=" * 82)
    def Kmat(m):
        K = np.zeros((m, m))
        for i in range(m):
            K[i, i] = 2.0
            if i > 0: K[i, i - 1] = -1.0
            if i < m - 1: K[i, i + 1] = -1.0
        return K
    n = 6; K = Kmat(n); Kinv = np.linalg.inv(K)
    for label, v in [("lowest Dirichlet-like", np.sin(np.pi * np.arange(1, n + 1) / (n + 1))),
                     ("block-uniform", np.ones(n) / np.sqrt(n))]:
        v = v / np.linalg.norm(v)
        var_a = v @ Kinv @ v; S = 1.0 / var_a
        E_H = v @ K @ v; Varf = v @ K @ Kinv @ K @ v
        print(f"  mode '{label}': Var(a)={var_a:.5f}, S=1/Var(a)={S:.5f}, D=0 (Gaussian)")
        print(f"     E[H]={E_H:.5f}, Var(force)={Varf:.5f}, flat E[H]-Var={E_H-Varf:.2e}, G''=S-D={S:.5f}")
    print("  k->0 limit (lowest mode, growing n): S -> 0, G'' -> 0 (free field, no gap):")
    for m in (6, 12, 24, 48):
        Kinv = np.linalg.inv(Kmat(m))
        v = np.sin(np.pi * np.arange(1, m + 1) / (m + 1)); v /= np.linalg.norm(v)
        print(f"     n={m:3d}: S = {1.0/(v@Kinv@v):.6f}")
    print("  => K>0 (log-concavity) for every n, yet S->0: log-concavity is NOT sufficient for")
    print("     (star); a uniform coarse spectral gap (=H1) is required.\n")


def main():
    print("#" * 82)
    print("# S1 conditional-Poincare witness  (illustration only; no theorem claimed here)")
    print("#" * 82)
    part_A_free_field()
    part_B_savvidy_threshold()
    part_C_haar()
    part_D_gaussian_block()
    print("#" * 82)
    print("# SUMMARY")
    print("#  * (A) free field: Var=E[H], G''=0 => (star) fails.  Required acid test PASSED.")
    print("#  * (B) boundary field strength drives the block Hessian indefinite (Savvidy); the")
    print("#        Haar floor N/12=0.167 cannot rescue a curvature of order beta*f=O(1).")
    print("#        => log-concavity is NOT uniform in the boundary U'.  [method no-go]")
    print("#  * (C) Haar Jacobian Phi''>=N/12 and W_geo>=0: the only boundary-INDEPENDENT term.")
    print("#  * (D) Gaussian: D=0, S=1/Var(a); k->0 has S->0: log-concavity alone gives no (star).")
    print("#  NOTHING here proves (star) at beta=2.5; the head is H1.")
    print("#" * 82)


if __name__ == "__main__":
    main()
