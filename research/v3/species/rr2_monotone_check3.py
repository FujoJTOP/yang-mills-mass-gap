"""RR2-MONOTONE v3 -- cross-checks, derivative identity, and the multi-block (additivity) test.

Conventions as v2.  Layer form for degenerate first level (QQ8 0.3):
    e = sum_{n in L1} |e_0n|^2 / sum_{n>=1} |e_0n|^2 ,  L1 = lowest excited layer.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal


def ladder(N):
    return np.array([k * (k + 2) / 4.0 for k in range(N)], float)


def e_layer(E, W, XEmat, rtol=1e-8):
    """e in the eigenbasis W of H; e_0n = <E_0|X_E|E_n> = (W^T XE W)[0,n]; layer-sum over L1."""
    ep = W.T @ XEmat @ W
    row = ep[0, 1:]
    Eg = E[1:]
    L1 = np.abs(Eg - Eg[0]) <= rtol * max(1.0, abs(Eg[0]))
    num = np.sum(row[L1] ** 2)
    den = np.sum(row ** 2)
    return num / den


def one_block(N, t, lam=1.0, basis="elec"):
    d, off = ladder(N), np.ones(N - 1)
    if basis == "elec":
        E, V = eigh_tridiagonal(lam * d, t * off)
        return e_layer(E, V, np.diag(d))
    else:
        b, Vb = eigh_tridiagonal(np.zeros(N), off)
        XEb = Vb.T @ np.diag(lam * d) @ Vb
        E, W = np.linalg.eigh(t * np.diag(b) + XEb)
        return e_layer(E, W, np.diag(b))


def r_inf_1block(N, lam=1.0):
    """static r = sum_{L1}|E0n|^2 / |E01|^2  in the X_B eigenbasis (simple spectrum)."""
    d, off = ladder(N), np.ones(N - 1)
    b, Vb = eigh_tridiagonal(np.zeros(N), off)
    XE = Vb.T @ np.diag(lam * d) @ Vb
    E0 = XE[0, 1:]
    return float(np.sum(E0[1:] ** 2) / E0[0] ** 2)


def block_chain(M, N, lam=1.0):
    """M blocks (links) in tensor product: X_B = sum_i P_i (commuting), X_E = sum_i D_i."""
    d = ladder(N)
    off = np.ones(N - 1)
    b, Vb = eigh_tridiagonal(np.zeros(N), off)          # single-block X_B eigenbasis
    eps = np.diag(ladder(N))
    DE = Vb.T @ eps @ Vb                                 # single-block X_E in b-basis (dense)
    # build M-fold Kronecker sums
    bsum = np.zeros(1)
    DEsum = np.zeros((1, 1))
    for _ in range(M):
        bsum = np.add.outer(bsum, b).ravel()
        DEsum = np.kron(DEsum, np.eye(N)) + np.kron(np.eye(DEsum.shape[0]), DE)
    return bsum, DEsum


def r_inf_M(M, N, lam=1.0, rtol=1e-8):
    bsum, DEsum = block_chain(M, N, lam)
    dim = len(bsum)
    if dim > 4096:
        return np.nan, np.nan
    E = bsum                                    # X_B is diagonal in this basis (eigenvalues of H at t=inf)
    # X_E in the same basis is DEsum; the "atoms" are bsum, weights are DEsum[0,n]^2
    row = DEsum[0, :]
    order = np.argsort(E)
    E, row = E[order], row[order]
    L1 = np.abs(E[1:] - E[0]) <= rtol * max(1.0, abs(E[0]))
    num = np.sum(row[1:][L1] ** 2)
    den = np.sum(row[1:] ** 2)
    return float(np.sum(row[1:][~L1] ** 2) / num), float(num / den)


def main():
    print("=== [A] basis cross-check (N=8): electric vs X_B-diagonal ===")
    for t in (0.5, 0.78125, 1.0, 2.0, 3.0, 6.0, 20.0):
        a = one_block(8, t, basis="elec")
        bb = one_block(8, t, basis="b")
        print(f"  t={t:8.5f}  e_elec={a:.12f}  e_b={bb:.12f}  diff={abs(a-bb):.2e}")

    print("\n=== [B] analytic dm_n/dt vs central finite difference (N=6, t=1.3) ===")
    N, t = 6, 1.3
    d, off = ladder(N), np.ones(N - 1)
    P = np.diag(off, 1) + np.diag(off, -1)
    E, V = eigh_tridiagonal(d, t * off)
    m = V.T @ (P @ V[:, 0])                      # m_n = <E_0|X_B|E_n>
    M = V.T @ P @ V
    dm = np.zeros(N)
    for n in range(1, N):
        s = 0.0
        for k in range(1, N):
            if k != n:
                s += m[k] * M[k, n] * (1.0 / (E[n] - E[k]) - 1.0 / (E[k] - E[0]))
        s += m[n] * m[0] / (E[n] - E[0])
        dm[n] = s
    h = 1e-6
    Ep, Vp = eigh_tridiagonal(d, (t + h) * off)
    Em, Vm = eigh_tridiagonal(d, (t - h) * off)
    mp = Vp.T @ (P @ Vp[:, 0])
    mm = Vm.T @ (P @ Vm[:, 0])
    fd = (mp - mm) / (2 * h)
    print(f"  max|analytic - FD| / max|FD| = {np.max(np.abs(dm[1:]-fd[1:]))/np.max(np.abs(fd[1:])):.2e}")
    br = np.array([2 * E[k] - E[0] - E[n] for k in range(1, N) if k != n])
    print("  bracket numerators 2E_k-E_0-E_n (k=1..):", np.round(np.array([2*E[k]-E[0]-E[n] for k in range(1,N)]),4))

    print("\n=== [C] local extrema of e(t), N=8, window t in [0.2,30] ===")
    tg = np.logspace(np.log10(0.2), np.log10(30), 4000)
    v = np.array([one_block(8, t) for t in tg])
    dv = np.diff(v)
    sgn = np.sign(dv)
    flips = np.where(np.diff(sgn) != 0)[0]
    print(f"  e(0.2)={v[0]:.6f}  e(t_A)={one_block(8,0.78125):.6f}  e(3)={one_block(8,3.0):.6f}  "
          f"e(30)={v[-1]:.6f}  e_inf={1/(1+r_inf_1block(8)):.6f}")
    for f in flips[:6]:
        print(f"    turning point near t={tg[f+1]:.4f}: e={v[f+1]:.8f}  (prev slope {dv[f]:+.2e})")
    tv = np.sum(np.abs(np.diff(np.log(v))))
    print(f"  TV(ln e) over [0.2,30] = {tv:.5f}")

    print("\n=== [D] r_inf: single block vs M-block chain (additivity test) ===")
    for M in (1, 2, 3):
        for N in (3, 4, 5, 6, 8):
            if N ** M > 5000:
                continue
            r, ein = r_inf_M(M, N) if M > 1 else (r_inf_1block(N), 1 / (1 + r_inf_1block(N)))
            print(f"  M={M} N={N:>2} dim={N**M:>5}: r_inf={r:12.6f}  e_inf={ein:.6f}")
    print("\n  single-block r_inf(N):", [round(r_inf_1block(N), 6) for N in (3, 4, 5, 6, 8, 12, 16)])
    print(f"  63/r_inf(120) = {63/r_inf_1block(120):.2f} ; 21.920736/r_inf(120) = {21.920736/r_inf_1block(120):.2f}")


if __name__ == "__main__":
    main()
