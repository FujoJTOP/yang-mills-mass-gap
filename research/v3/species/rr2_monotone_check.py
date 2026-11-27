"""RR2-MONOTONE -- final consolidated check (deterministic eigh only; zero MC, zero sampling).

Q (round 245): can the O-sector Jacobi (tridiagonal) structure sign d_t e(t) ?
Canonical pair (SU(2), single block, chi_{1/2} fusion = path; QQ3-CR (P1)):
    H(t) = t*P + diag(eps) ,  P = path adjacency (off-diag 1) ,  eps_k = k(k+2)/4 = C2(j=k/2)
    e(t) = |e_01|^2 / sum_{n>=1}|e_0n|^2 ,  e_0n = <E_0|X_E|E_n>   (layer form, QQ8 (R1))
    exact: <E_0|H|E_n> = 0  =>  t*m_n + e_0n = 0 , m_n = <E_0|X_B|E_n>
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal

TA = 2.5 ** 2 / 8.0
W0 = 1.0 / 64.0
B_BUDGET = np.log(21.46752)          # QQ8: B = 3.066541 nats (eps(2.5) = 21.46752 w0)


def ladder(N):
    return np.array([k * (k + 2) / 4.0 for k in range(N)], float)


def layer(E, row, rtol=1e-9):
    Eg, rg = E[1:], row[1:]
    if Eg.size == 0:
        return np.nan, 0
    tol = rtol * max(1.0, abs(Eg[0]))
    m = np.abs(Eg - Eg[0]) <= tol
    return float(np.sum(rg[m] ** 2) / np.sum(rg ** 2)), int(m.sum())


def e_pair(eps, off, t):
    E, V = eigh_tridiagonal(eps, t * off)
    return layer(E, (V.T @ np.diag(eps) @ V)[0, :])[0]


def e_t(t, N, shift=0.0, lam=1.0):
    eps, off = lam * (ladder(N) + shift), np.ones(N - 1)
    return layer(eigh_tridiagonal(eps, t * off)[:1] + (0,), 0) if False else e_pair(eps, off, t)


def e_bend(t, N):
    """same e in the X_B-diagonal basis (independent code path)"""
    eps, off = ladder(N), np.ones(N - 1)
    b, Vb = eigh_tridiagonal(np.zeros(N), off)
    XEb = Vb.T @ np.diag(eps) @ Vb
    E, W = np.linalg.eigh(t * np.diag(b) + XEb)
    return layer(E, (W.T @ np.diag(b) @ W)[0, :])[0]


def r_inf(N, shift=0.0, lam=1.0):
    eps, off = lam * (ladder(N) + shift), np.ones(N - 1)
    b, Vb = eigh_tridiagonal(np.zeros(N), off)
    XE = Vb.T @ np.diag(eps) @ Vb
    E0 = XE[0, 1:]
    return float(np.sum(E0[1:] ** 2) / E0[0] ** 2)


def align(V, ref):
    for j in range(V.shape[1]):
        if V[:, j] @ ref[:, j] < 0:
            V[:, j] *= -1
    return V


def main():
    print(f"t_A={TA:.6f}  w0={W0:.8f}  B={B_BUDGET:.6f}")

    print("\n[S1] exact identity |t*m_n + e_0n| (N=6)")
    d0, off0 = ladder(6), np.ones(5)
    P6 = np.diag(off0, 1) + np.diag(off0, -1)
    for t in (0.5, 0.78125, 3.0, 50.0, 1000.0):
        E, V = eigh_tridiagonal(d0, t * off0)
        u, m = V.T @ (d0 * V[:, 0]), V.T @ (P6 @ V[:, 0])
        print(f"   t={t:9.5f}  max|t*m+e| = {np.max(np.abs(t * m[1:] + u[1:])):.2e}")

    print("\n[S2] derivative formula  dm_n/dt  vs sign-aligned central FD (N=6, t=1.3)")
    N, t = 6, 1.3
    d, off = ladder(N), np.ones(N - 1)
    P = np.diag(off, 1) + np.diag(off, -1)
    E, V = eigh_tridiagonal(d, t * off)
    m, M = V.T @ (P @ V[:, 0]), V.T @ P @ V
    dm = np.zeros(N)
    for n in range(1, N):
        s = sum(m[k] * M[k, n] * (1 / (E[n] - E[k]) - 1 / (E[k] - E[0])) for k in range(1, N) if k != n)
        dm[n] = s + m[n] * (M[0, 0] - M[n, n]) / (E[n] - E[0])
    h = 1e-5
    Vp = align(eigh_tridiagonal(d, (t + h) * off)[1], V)
    Vm = align(eigh_tridiagonal(d, (t - h) * off)[1], V)
    fd = (Vp.T @ (P @ Vp[:, 0]) - Vm.T @ (P @ Vm[:, 0])) / (2 * h)
    print(f"   rel.err = {np.max(np.abs(dm[1:] - fd[1:])) / np.max(np.abs(fd[1:])):.2e}")
    print("   brackets 1/(E_n-E_k)-1/(E_k-E_0), k=2..:", np.round([1 / (E[1] - E[k]) - 1 / (E[k] - E[0])
                                                                  for k in range(2, N)], 4))

    print("\n[S3] basis cross-check (N=8): electric vs X_B-diagonal")
    for t in (0.78125, 1.0, 3.0, 20.0):
        print(f"   t={t:8.5f}: elec={e_t(t, 8):.12f}  b={e_bend(t, 8):.12f}  diff={abs(e_t(t,8)-e_bend(t,8)):.1e}")

    print("\n[S4] canonical ladder: shape on [0.2,30] (turning points) and e_inf")
    for N in (4, 6, 8, 16):
        tg = np.logspace(np.log10(0.2), np.log10(30), 3000)
        v = np.array([e_t(x, N) for x in tg])
        fl = np.where(np.diff(np.sign(np.diff(v))) != 0)[0]
        ei = 1 / (1 + r_inf(N))
        s = ", ".join(f"t={tg[f+1]:.4f}/e={v[f+1]:.8f}" for f in fl[:3])
        print(f"   N={N:>2}: e(0.2)={v[0]:.7f} e(t_A)={e_t(TA,N):.7f} e(30)={v[-1]:.7f} "
              f"e_inf={ei:.7f} #turn={len(fl)} [{s}] TV={np.sum(np.abs(np.diff(np.log(v)))):.5f}")

    print("\n[S5] budget on [t_A,1e7]")
    tg = np.unique(np.concatenate([np.logspace(np.log10(TA), 7, 1500), np.logspace(-1, 1.5, 800)]))
    for N in (4, 8, 16):
        v = np.array([e_t(x, N) for x in tg])
        dl = np.diff(np.log(v))
        tv = np.sum(np.abs(dl))
        print(f"   N={N:>2}: +var={np.sum(dl[dl>0]):.5f} -var={-np.sum(dl[dl<0]):.5f} TV={tv:.5f} "
              f"B/TV={B_BUDGET/tv:.2f} e_min={v.min():.8f}@t={tg[int(np.argmin(v))]:.3g} "
              f"e_inf={1/(1+r_inf(N)):.8f}")

    print("\n[S6] value layer (static): r_inf vs N  -- independent repro of QQ3 0.3")
    for N in (3, 4, 8, 16, 64, 128, 256):
        print(f"   N={N:>4}: r_inf={r_inf(N):.9f}  e_inf={1/(1+r_inf(N)):.9f}")
    r = r_inf(256)
    r25 = 1 / e_t(TA, 16) - 1
    print(f"   r_inf(256)={r:.9f} = {63/r:.2f}x inside 63 ; 21.920736/r = {21.920736/r:.2f} ; "
          f"e_inf = {64/(1+r):.5f} w0")
    print(f"   r(2.5) in the pair = {r25:.6f} ; 63/that = {63/r25:.1f}")

    print("\n[S7] invariance 1: X_E -> X_E + c*I   (atoms move, e must not)")
    for c in (0.0, 0.7, 5.0):
        print(f"   c={c:4.1f}: e(t_A)={e_t(TA, 8, shift=c):.12f}  e(3.0)={e_t(3.0, 8, shift=c):.12f}")
    for sh in (0.0, 0.25, 3.0):
        print(f"   shift={sh:5.2f}: r_inf(16)={r_inf(16, shift=sh):.9f}")

    print("\n[S8] invariance 2: independent-block tensor extension (QQ3 (T2), layer form)")
    N = 4
    b, Vb = eigh_tridiagonal(np.zeros(N), np.ones(N - 1))
    DE = Vb.T @ np.diag(ladder(N)) @ Vb
    I = np.eye(N)
    bsum, DEs = np.zeros(1), np.zeros((1, 1))
    for M in (1, 2, 3, 4):
        bsum = np.add.outer(bsum, b).ravel()
        DEs = np.kron(DEs, I) + np.kron(np.eye(DEs.shape[0]), DE)
        o = np.argsort(bsum)
        val, nl = layer(bsum[o], DEs[0, :][o])
        print(f"   M={M}: dim={N**M:>4} |L1|={nl} e_inf={val:.9f} r_inf={1/val-1:.9f}")

    print("\n[S9] generic Jacobi class: monotone? (60 draws, seed 20261127)")
    rng = np.random.default_rng(20261127)
    tg = np.unique(np.concatenate([np.logspace(np.log10(TA), 6, 500), np.logspace(-1, 1.5, 400)]))
    for N in (4, 6, 8):
        bad, mx, tvs, rs = 0, 0.0, [], []
        for _ in range(60):
            dd = np.sort(rng.uniform(0, 12, N)); dd[0] = 0.0
            of = rng.uniform(0.2, 1.8, N - 1)
            v = np.array([e_pair(dd, of, x) for x in tg])
            dl = np.diff(np.log(v))
            mx = max(mx, float(np.max(dl)))
            tvs.append(float(np.sum(np.abs(dl))))
            bad += int(np.max(dl) > 1e-6)
            bb, Vb2 = eigh_tridiagonal(np.zeros(N), of)
            XE = Vb2.T @ np.diag(dd) @ Vb2
            E0 = XE[0, 1:]
            rs.append(np.sum(E0[1:] ** 2) / E0[0] ** 2)
        print(f"   N={N}: non-monotone {bad}/60 ; max dln e = {mx:.2e} ; TV med {np.median(tvs):.3f} "
              f"max {max(tvs):.2f} ; r_inf in [{min(rs):.3g},{max(rs):.3g}] med {np.median(rs):.3g}")

    print("\n[S10] QQ8 0.10 vs QQ3 0.3 : end convention excluded?  (diagonal element must be dropped)")
    for N in (8, 64, 256):
        eps, off = ladder(N), np.ones(N - 1)
        b, Vb = eigh_tridiagonal(np.zeros(N), off)
        XE = Vb.T @ np.diag(eps) @ Vb
        rl = np.sum(XE[1:, 0] ** 2) / XE[1, 0] ** 2
        rh = np.sum(XE[:-1, -1] ** 2) / XE[-2, -1] ** 2
        print(f"   N={N:>3}: r(b_0=lowest)={rl:.9f} r(b_0=highest)={rh:.9f} ratio={rh/rl:.6f}")


if __name__ == "__main__":
    main()
