"""
ICA witness: is the Irreversible-Constraint-Accumulation distance D_n a NEW monotone,
or a known one (conditional entropy / relative entropy of recovery / free energy / RG monotone)?

Exact models only.  We work with the EXACT, finite object
    I_n := I(x ; y_n)          mutual information between the fine field x and its n-fold
                              block-average y_n = T^n x,
and the accumulated irreversibility
    D_n := I_0 - I_n           ( I_0 = H(x) ),   so  D_{n+1} >= D_n  exactly  (DPI).
This D_n is identical to H(x | y_n) (differential, constant-mode fixed); it IS the data-
processing deficit = conditional entropy = relative entropy of recovery (Gaussian/Petz case).

Models:
  (A) Gaussian free scalar, 1D (L=64) and 2D (16x16) periodic; m=0 (massless) and m>0.
      Constant mode fixed to 0 -> everything finite and exact.
  (B) exact 1D Z_N gauge on a circle: no local DOF -> D_n = 0 identically (gap but ICA=0).
  (C) exact 2D Z_N gauge (independent plaquettes): D_block = b^2 H(q) - H(q^{*b^2}),
      an entropy / free-energy difference of merged plaquettes (known thermodynamic quantity).

Run:  py -3.12 -I ica_witness.py
"""
import numpy as np

np.set_printoptions(precision=6, suppress=True)


# ---------- centered (zero-mean) subspace basis ----------
def centered_basis(n):
    """Orthonormal basis (n x (n-1)) of the orthogonal complement of 1/sqrt(n)."""
    Q, _ = np.linalg.qr(np.ones((n, 1)), mode="complete")
    return Q[:, 1:]


# ---------- Gaussian free field, 1D ----------
def lap_eigs_1d(L):
    k = np.arange(L)
    return 4.0 * np.sin(np.pi * k / L) ** 2


def cov_1d(L, m):
    lam = lap_eigs_1d(L) + m * m
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = 1.0 / lam
    inv[0] = 0.0                      # fix the constant mode
    corr = np.fft.ifft(inv).real      # C[i,j] = corr[(i-j) mod L]
    d = np.arange(L)
    return corr[(d[:, None] - d[None, :]) % L]


def block_matrix_1d(L, b):
    M = L // b
    A = np.zeros((M, L))
    for I in range(M):
        for j in range(b):
            A[I, I * b + j] = 1.0 / b
    return A


def _conditional_cov(C, A):
    Lf = C.shape[0]
    Mc = A.shape[0]
    Bf = centered_basis(Lf)
    Bc = centered_basis(Mc)
    Cr = Bf.T @ C @ Bf
    Ar = Bc.T @ A @ Bf
    Ccr = Ar @ Cr @ Ar.T
    Cxcr = Cr @ Ar.T
    Conr = Cr - Cxcr @ np.linalg.solve(Ccr, Cxcr.T)
    Conr = 0.5 * (Conr + Conr.T)
    return Cr, Conr


def mutual_info_centered(C, A):
    """I(x;y) for zero-mean jointly-Gaussian (x,y=A x), constant mode removed."""
    Cr, Conr = _conditional_cov(C, A)
    _, lr = np.linalg.slogdet(Cr)
    so, lo = np.linalg.slogdet(Conr)
    if so <= 0:                        # numerical PSD safeguard
        w = np.linalg.eigvalsh(Conr)
        lo = np.sum(np.log(np.clip(w, 1e-14, None)))
    return 0.5 * (lr - lo)


def cond_entropy_centered(C, A):
    """H(x | y=A x) (differential), constant mode removed = the info LOST after the step.
       Equals D_n (accumulated irreversibility) since D_0 = 0."""
    _, Conr = _conditional_cov(C, A)
    w = np.linalg.eigvalsh(Conr)
    w = w[w > 1e-12]
    if w.size == 0:
        return 0.0
    return 0.5 * float(np.sum(np.log(2 * np.pi * np.e * w)))


def part_A_1d():
    print("== (A) Gaussian free scalar, 1D periodic L=512, block b=16 ==")
    print("   D_1 = H(x | y_1) = conditional (differential) entropy of the detail one block")
    print("   step discards  == the one-step irreversibility distance.  Ratchet: D_n >= D_{n-1}.")
    L, b = 512, 16
    A = block_matrix_1d(L, b)
    C0 = cov_1d(L, 0.0)
    D0 = cond_entropy_centered(C0, A)
    print(f"   {'m':>6}  {'D_1':>10}  {'D_1(0)-D_1(m)':>15}  (massless: LARGEST lost info)")
    for m in [0.0, 0.02, 0.05, 0.1, 0.3, 1.0]:
        C = cov_1d(L, m)
        d = cond_entropy_centered(C, A)
        print(f"   {m:>6.2f}  {d:>10.4f}  {D0 - d:>15.4f}")
    print(f"   leading term D_1 ~ {D0:.0f} nats = UV mode count (gap-blind);")
    print("   the entire mass dependence is the O(1) column on the right.")
    print("   => SIGN: m=0 (no gap) gives the biggest D; a finite mass makes the field")
    print("      locally recoverable and SHRINKS D.  A gap implies SMALL D -> wrong direction.")
    print()
    print("   ratchet check at m=0 (D_n non-decreasing; never saturates -> no gap):")
    C = cov_1d(L, 0.0)
    prev = 0.0
    for n in range(1, 5):
        A = block_matrix_1d(L, b * (2 ** (n - 1)))
        d = cond_entropy_centered(C, A)
        print(f"     n={n}: D_n = {d:10.4f}   increment = {d - prev:8.4f}")
        prev = d
    print()


def part_A_2d():
    print("== (A') Gaussian free scalar, 2D periodic 16x16, block 2x2 ==")
    L2 = 16
    ky = np.arange(L2)
    lam1d = 4.0 * np.sin(np.pi * ky / L2) ** 2
    base = lam1d[:, None] + lam1d[None, :]      # (L2,L2) Laplacian eigenvalues
    for m in [0.0, 0.3, 1.0]:
        with np.errstate(divide="ignore", invalid="ignore"):
            inv2 = 1.0 / (base + m * m)
        inv2[0, 0] = 0.0                           # fix constant mode
        corr2 = np.fft.ifft2(inv2).real            # C(dy,dx) = (1/L^2) sum_k inv_k e^{...}
        idx = np.arange(L2 * L2)
        yi, xi = np.divmod(idx, L2)
        Dy = (yi[:, None] - yi[None, :]) % L2
        Dx = (xi[:, None] - xi[None, :]) % L2
        C2 = 0.5 * (corr2[Dy, Dx] + corr2[Dy, Dx].T)
        M = L2 // 2
        A2 = np.zeros((M * M, L2 * L2))
        for Y in range(M):
            for X in range(M):
                row = Y * M + X
                for dy in range(2):
                    for dx in range(2):
                        A2[row, (2 * Y + dy) * L2 + (2 * X + dx)] = 0.25
        D1 = cond_entropy_centered(C2, A2)
        I1 = mutual_info_centered(C2, A2)
        print(f"   m={m:>4.2f}:  D_1 = H(x|y_1) = {D1:8.4f} nats  (I(x;y_1) = {I1:8.4f})"
              + ("   <- free massless: nonzero -> FAIL" if m == 0.0 else ""))
    print()


# ---------- exact Z_N gauge ----------
def zn_single_plaquette_q(beta, N):
    s = np.arange(N)
    w = np.exp(beta * np.cos(2 * np.pi * s / N))
    return w / w.sum()


def zn_fold_convolve(q, k, N):
    p = q.copy()
    for _ in range(k - 1):
        raw = np.convolve(p, q)
        p2 = np.zeros(N)
        for i, val in enumerate(raw):
            p2[i % N] += val
        p = p2
    return p


def H(p):
    p = p[p > 0]
    return -np.sum(p * np.log(p))


def part_B_C():
    print("== (B) exact 1D Z_N gauge on a circle ==")
    print("   1D has NO local gauge-invariant DOF; only the global holonomy survives blocking")
    print("   => D_n = 0 identically for every N, beta.   [transfer gap but ICA = 0]")
    print()

    print("== (C) exact 2D Z_N gauge (independent plaquettes), block b x b ==")
    print("   D_block = b^2 H(q) - H(q^{*b^2}) = entropy / free-energy difference of merged plaquettes")
    for N in [2, 3, 4]:
        print(f"   N={N}:")
        print(f"     {'beta':>5} {'H(q)':>9} {'D_b=2':>9} {'D_b=3':>9} {'D_b=4':>9}")
        for beta in [0.5, 1.0, 2.0]:
            q = zn_single_plaquette_q(beta, N)
            hq = H(q)
            row = [b * b * hq - H(zn_fold_convolve(q, b * b, N)) for b in [2, 3, 4]]
            print(f"     {beta:>5.1f} {hq:>9.4f} " + " ".join(f"{v:>9.4f}" for v in row))
        print("     (all D_block > 0: coarse holonomy does NOT recover the merged plaquettes;")
        print("      nonzero in a theory with sigma>0 and NO 0++ atom -> free energy, not the atom)")
    print()


def main():
    print("ICA witness  (numpy only, deterministic).  No QED; verdict is a KILL candidate.\n")
    part_A_1d()
    part_A_2d()
    part_B_C()
    print("SUMMARY")
    print(" A: D_n = conditional entropy / data-processing deficit -- a KNOWN monotone;")
    print("    its non-decrease IS the data-processing inequality; >0 on free massless -> ACID FAIL.")
    print(" B: 1D exact -> D_n = 0 with a transfer gap -> GAP but ICA = 0.")
    print(" C: 2D exact -> D_block = free-energy/entropy difference of merged plaquettes -> KNOWN,")
    print("    nonzero WITHOUT a 0++ atom, and extensive (gap-blind leading term).")


if __name__ == "__main__":
    main()
