"""RR2-MONOTONE check v2 -- does the Jacobi (tridiagonal) O-sector structure sign d_t e(t)?

Canonical O-sector Jacobi model (SU(2) single-block KS truncation; FOSP-A: fusion graph = path):
    H(t) = t*P + lam*diag(eps),  eps_k = k(k+2)/4 = C2(j=k/2),  P = path adjacency (off-diag 1)
    e(t)  = |u_1|^2 / sum_{n>=1}|u_n|^2,   u_n = <E_0|X_E|E_n> = sum_k eps_k V[k,0]V[k,n]
    exact identity  <E_0|H|E_n>=0  =>  t*m_n + u_n = 0  (m_n = <E_0|X_B|E_n>)
Deterministic (numpy/scipy eigh) only. No MC.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal

TA = 2.5 ** 2 / 8.0  # 0.78125
W0 = 1.0 / 64.0


def ladder(N):
    return np.array([k * (k + 2) / 4.0 for k in range(N)], float)


def e_of_t(t, d, off, lam=1.0, want=False):
    E, V = eigh_tridiagonal(lam * d, t * off)
    v0 = V[:, 0]
    u = V.T @ (d * v0)          # u_n = <E_n|D|E_0>
    exc = u[1:]
    val = exc[0] ** 2 / np.sum(exc ** 2)
    return (val, E, V, u) if want else val


def pos_neg_split(grid, vals):
    dl = np.diff(np.log(vals))
    return float(np.sum(dl[dl > 0])), float(-np.sum(dl[dl < 0])), float(np.max(dl))


def r_inf(N, lam=1.0):
    """static ratio r = sum_{n>=2}|E0n|^2/|E01|^2 in the X_B (path) eigenbasis."""
    d = ladder(N)
    off = np.ones(N - 1)
    Eb, Vb = eigh_tridiagonal(np.zeros(N), off)
    XE = Vb.T @ np.diag(lam * d) @ Vb
    E0 = XE[0, 1:]
    return float(np.sum(E0[1:] ** 2) / E0[0] ** 2)


def main():
    print(f"t_A={TA:.6f}  w0={W0:.8f}")
    print("\n=== [1] canonical ladder, lam=1: TV split on [t_A, inf) ===")
    print(f"{'N':>3} {'e(t_A)':>10} {'e_min':>10} {'t_min':>8} {'e_inf':>10} {'r_inf':>9} "
          f"{'+var':>8} {'-var':>8} {'max dln+':>9}")
    for N in (3, 4, 5, 6, 8, 10, 12, 16):
        d, off = ladder(N), np.ones(N - 1)
        tg = np.unique(np.concatenate([np.logspace(np.log10(TA), 7, 1400),
                                       np.logspace(np.log10(0.5), np.log10(60), 2000)]))
        v = np.array([e_of_t(t, d, off) for t in tg])
        pv, nv, mx = pos_neg_split(tg, v)
        i = int(np.argmin(v))
        print(f"{N:>3} {e_of_t(TA,d,off):10.6f} {v[i]:10.6f} {tg[i]:8.3g} "
              f"{1/(1+r_inf(N)):10.6f} {r_inf(N):9.5f} {pv:8.4f} {nv:8.4f} {mx:9.2e}")

    print("\n=== [2] lam-sensitivity (N=8): normalization of X_E vs X_B ===")
    for lam in (0.1, 0.25, 1.0, 4.0, 10.0):
        d, off = ladder(8), np.ones(7)
        tg = np.unique(np.concatenate([np.logspace(np.log10(TA), 7, 1200),
                                       np.logspace(np.log10(0.5), np.log10(60), 1500)]))
        v = np.array([e_of_t(t, d, off, lam) for t in tg])
        pv, nv, mx = pos_neg_split(tg, v)
        print(f"  lam={lam:5.2f}: e(t_A)={v[0]:.6f} e_min={v.min():.6f}@{tg[int(np.argmin(v))]:.3g} "
              f"+var={pv:.4f} -var={nv:.4f} max dln+={mx:.2e}  r_inf={r_inf(8,lam):.5f}")

    print("\n=== [3] r_inf vs N (static, exact; convergence) ===")
    for N in (3, 4, 5, 6, 8, 10, 12, 16, 20, 30, 40, 60, 80, 120):
        print(f"  N={N:>4}: r_inf={r_inf(N):.6f}   e_inf={1/(1+r_inf(N)):.6f}")

    print("\n=== [4] N=3 shape of q(t)=u_2/u_1 (is e=1/(1+q^2) monotone?) ===")
    d, off = ladder(3), np.ones(2)
    tg = np.logspace(np.log10(0.05), 4, 4000)
    qs = []
    for t in tg:
        _, E, V, u = e_of_t(t, d, off, want=True)
        qs.append(u[2] / u[1])
    qs = np.array(qs)
    dq = np.diff(qs)
    print(f"  q: {qs[0]:.6e} -> {qs[-1]:.6e};  #dq>0={int(np.sum(dq>0))}/{len(dq)}; "
          f"max dq={dq.max():.3e}; min dq={dq.min():.3e};  sign-flips={int(np.sum(np.diff(np.sign(dq))!=0))}")

    print("\n=== [5] random Jacobi class: monotone? r_inf range ===")
    rng = np.random.default_rng(20261127)
    tg = np.unique(np.concatenate([np.logspace(np.log10(TA), 6, 800),
                                   np.logspace(np.log10(0.5), np.log10(60), 900)]))
    for N in (4, 5, 6, 8):
        bad = 0
        rs, tvs, mxall = [], [], []
        for _ in range(60):
            d = np.sort(rng.uniform(0.0, 12.0, N)); d[0] = 0.0
            off = rng.uniform(0.2, 1.8, N - 1)
            v = np.array([e_of_t(t, d, off) for t in tg])
            pv, nv, mx = pos_neg_split(tg, v)
            tvs.append(pv + nv); mxall.append(mx)
            if mx > 1e-6:
                bad += 1
            Eb, Vb = eigh_tridiagonal(np.zeros(N), off)
            XE = Vb.T @ np.diag(d) @ Vb
            E0 = XE[0, 1:]
            rs.append(np.sum(E0[1:] ** 2) / E0[0] ** 2)
        print(f"  N={N}: non-monotone {bad}/60; TV med {np.median(tvs):.3f} max {max(tvs):.2f}; "
              f"max step {max(mxall):.2e}; r_inf min {min(rs):.3g} med {np.median(rs):.3g} max {max(rs):.3g}")

    print("\n=== [6] invariance: e is exactly invariant under D -> D + c*I (spectral shift) ===")
    d, off = ladder(8), np.ones(7)
    for c in (0.0, 0.7, 5.0):
        v = e_of_t(0.78125, d + c, off)
        w = e_of_t(3.0, d + c, off)
        print(f"  c={c:5.2f}: e(t_A)={v:.12f}  e(3.0)={w:.12f}")


if __name__ == "__main__":
    main()
