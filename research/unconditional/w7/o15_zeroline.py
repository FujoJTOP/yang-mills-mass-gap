"""O15b — zero set of E[H|0] in (N,theta,beta) + BETASTAR §4 magnitude cross-check."""
import numpy as np
from o15_etheta import build, stats, wmean, slices, Hd
from betastar_scaling import centre_B


def theta_c(pool, N, beta, lo=0.0, hi=np.pi):
    f = lambda t: stats(pool, N, t, beta)['EH']
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


print("[Z1] zero crossing of E[H|0](theta) between 3pi/4 and pi  (bisection, 40 iters)")
pools = {N: build(N, n) for N, n in ((2, 120000), (3, 120000), (4, 12000))}
for N in (2, 3, 4):
    for beta in (2.5, 6.0):
        tc = theta_c(pools[N], N, beta)
        print(f"   N={N} beta={beta:4.1f}: theta_c/pi = {tc / np.pi:.6f}   "
              f"E[H](theta_c)={stats(pools[N], N, tc, beta)['EH']:+.2e}")

print("[Z2] BETASTAR §4 check: Z_N centre, N=3, beta=6 -> their E[H]=+0.306, Var(G)=+0.0317")
for N, nw in ((2, 40000), (3, 40000)):
    B = centre_B(N)
    off, Tv, P, V = slices(N)
    iu = np.triu_indices(N, 1)
    rng = np.random.default_rng(9)
    hs, gs, ws, got = [], [], [], 0
    while got < nw:
        xi = rng.uniform(-2 * np.pi, 2 * np.pi, size=(40000, N * (N - 1)))
        X = np.zeros((40000, N, N), complex)
        for a, T in enumerate(off):
            X += xi[:, a, None, None] * T[None]
        mu = np.linalg.eigvalsh(X)
        d = (mu[:, :, None] - mu[:, None, :])[:, iu[0], iu[1]]
        keep = np.all(np.abs(d) <= 2 * np.pi, axis=1)
        if not keep.any():
            continue
        X, d = X[keep], d[keep]
        hf = d / 2
        J = np.prod((np.sin(hf) / hf) ** 2, axis=1)

        def fB(Xs):
            mu2, V2 = np.linalg.eigh(Xs)
            U = (V2 * np.exp(1j * mu2)[:, None, :]) @ V2.conj().transpose(0, 2, 1)
            return np.real(np.trace(U @ B, axis1=1, axis2=2)) / N

        W = np.exp(-6.0 * (1 - fB(X))) * J
        f0, fp, fm = fB(X), fB(X + Hd * Tv), fB(X - Hd * Tv)
        hh = -6.0 / N * (fp - 2 * f0 + fm) / Hd ** 2          # H = -(beta/N) d2 f
        hh_noN = -6.0 * (fp - 2 * f0 + fm) / Hd ** 2          # no 1/N (BETASTAR hypothesis)
        gg = -6.0 * (fp - fm) / (2 * Hd)                      # G = -beta d_v f
        m = min(nw - got, len(J))
        hs.append(hh[:m]); gs.append(gg[:m]); ws.append(W[:m]); got += m
    hh, gg, W = np.concatenate(hs), np.concatenate(gs), np.concatenate(ws)
    mh = wmean(hh, W)[0]
    mv = wmean(gg ** 2, W)[0] - wmean(gg, W)[0] ** 2
    print(f"   N={N}: <H>_(with 1/N)={mh:+.4f}   <H>_(no 1/N)={N * mh:+.4f}   Var(G)={mv:+.4f}")
    print(f"        BETASTAR-N §4 table (N=3, beta=6): E[H]=+0.3060, Var(G)=+0.0317")
