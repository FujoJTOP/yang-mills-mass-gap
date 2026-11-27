"""RR2-MONOTONE v4 -- fixes: layer mask relative to first excited level; sign-aligned FD check."""
import numpy as np
from scipy.linalg import eigh_tridiagonal


def ladder(N):
    return np.array([k * (k + 2) / 4.0 for k in range(N)], float)


def align(V, ref):
    for j in range(V.shape[1]):
        if V[:, j] @ ref[:, j] < 0:
            V[:, j] *= -1.0
    return V


def layer_ratio(E, row, rtol=1e-9):
    """row_n = <E_0|X_E|E_n>, n>=1 ; L1 = first excited layer (degenerate with E[1])."""
    Eg, rg = E[1:], row[1:]
    if len(Eg) == 0:
        return np.nan
    tol = rtol * max(1.0, abs(Eg[0]))
    L1 = np.abs(Eg - Eg[0]) <= tol
    num = np.sum(rg[L1] ** 2)
    den = np.sum(rg ** 2)
    return num / den, int(np.sum(L1))


def one_block(N, t, lam=1.0, basis="elec"):
    d, off = ladder(N), np.ones(N - 1)
    if basis == "elec":
        E, V = eigh_tridiagonal(lam * d, t * off)
        ep = V.T @ np.diag(d) @ V
        return layer_ratio(E, ep[0, :])
    b, Vb = eigh_tridiagonal(np.zeros(N), off)
    XEb = Vb.T @ np.diag(lam * d) @ Vb
    E, W = np.linalg.eigh(t * np.diag(b) + XEb)
    ep = W.T @ np.diag(b) @ W
    return layer_ratio(E, ep[0, :])


def r_inf_1block(N, lam=1.0):
    d, off = ladder(N), np.ones(N - 1)
    b, Vb = eigh_tridiagonal(np.zeros(N), off)
    XE = Vb.T @ np.diag(lam * d) @ Vb
    E0 = XE[0, 1:]
    return float(np.sum(E0[1:] ** 2) / E0[0] ** 2)


def block_static(M, N, lam=1.0, rtol=1e-9):
    """M independent blocks: X_B = sum_i P_i (diag in product b-basis), X_E = sum_i D_i."""
    off = np.ones(N - 1)
    b, Vb = eigh_tridiagonal(np.zeros(N), off)
    DE = Vb.T @ np.diag(ladder(N)) @ Vb
    bsum = np.zeros(1)
    DEsum = np.zeros((1, 1))
    I = np.eye(N)
    for _ in range(M):
        bsum = np.add.outer(bsum, b).ravel()
        DEsum = np.kron(DEsum, I) + np.kron(np.eye(DEsum.shape[0]), DE)
    order = np.argsort(bsum)
    E = bsum[order]
    row = DEsum[0, :][order]
    return layer_ratio(E, row, rtol)


def main():
    print("=== [B'] dm_n/dt: analytic vs sign-aligned central FD (N=6, t=1.3) ===")
    N, t = 6, 1.3
    d, off = ladder(N), np.ones(N - 1)
    P = np.diag(off, 1) + np.diag(off, -1)
    E, V = eigh_tridiagonal(d, t * off)
    m = V.T @ (P @ V[:, 0])
    M = V.T @ P @ V
    dm = np.zeros(N)
    for n in range(1, N):
        s = 0.0
        for k in range(1, N):
            if k != n:
                s += m[k] * M[k, n] * (1.0 / (E[n] - E[k]) - 1.0 / (E[k] - E[0]))
        s += m[n] * m[0] / (E[n] - E[0])
        dm[n] = s
    h = 1e-5
    Ep, Vp = eigh_tridiagonal(d, (t + h) * off)
    Em, Vm = eigh_tridiagonal(d, (t - h) * off)
    Vp, Vm = align(Vp, V), align(Vm, V)
    mp = Vp.T @ (P @ Vp[:, 0])
    mm = Vm.T @ (P @ Vm[:, 0])
    fd = (mp - mm) / (2 * h)
    rel = np.max(np.abs(dm[1:] - fd[1:])) / np.max(np.abs(fd[1:]))
    print(f"  max|analytic-FD|/max|FD| = {rel:.3e}   (analytic {np.round(dm[1:],5)})")
    print(f"                                  (FD       {np.round(fd[1:],5)})")

    print("\n=== [D'] r_inf and e_inf for M blocks x N levels (layer form) ===")
    print(f"{'M':>2} {'N':>3} {'dim':>5} {'#L1':>4} {'e_inf':>12} {'r_inf':>12} {'r_inf/M':>9}")
    for M in (1, 2, 3, 4, 5):
        for N in (3, 4, 5, 6, 8):
            if N ** M > 20000:
                continue
            e, nl = block_static(M, N)
            r = 1 / e - 1
            print(f"{M:>2} {N:>3} {N**M:>5} {nl:>4} {e:12.8f} {r:12.6f} {r/M:9.5f}")

    print("\n=== [D''] same, larger N, M=1 (converged value) ===")
    for N in (20, 40, 80, 160):
        r = r_inf_1block(N)
        print(f"  M=1 N={N:>4}: r_inf={r:.6f}  e_inf={1/(1+r):.8f}")

    print("\n=== [C'] canonical turning points and TV, N=8 (fine) ===")
    tg = np.logspace(np.log10(0.2), np.log10(30), 4000)
    v = np.array([one_block(8, t)[0] for t in tg])
    dv = np.diff(v)
    flips = np.where(np.diff(np.sign(dv)) != 0)[0]
    for f in flips[:8]:
        print(f"    turning point t={tg[f+1]:.4f}: e={v[f+1]:.9f}")
    print(f"  TV(ln e) [0.2,30] = {np.sum(np.abs(np.diff(np.log(v)))):.5f} ; "
          f"net = {abs(np.log(v[-1]/v[0])):.5f}")
    print(f"  e(0.2)={v[0]:.6f} e(t_A)={one_block(8,0.78125)[0]:.6f}  max={v.max():.6f}@{tg[int(np.argmax(v))]:.3f}"
          f"  e(30)={v[-1]:.6f}")


if __name__ == "__main__":
    main()
