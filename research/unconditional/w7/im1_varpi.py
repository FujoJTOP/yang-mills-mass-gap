"""
IM1: measure

    varpi(bhat) = bhat^4 * Var_mu( (1/N_B) sum_{x in B} q_x ),
    q_x = 1 - (1/N) Re tr U_square(x),   N = 2 (SU(2)),   N_B = bhat^4

4D SU(2) Wilson at beta=2.5 (S = beta sum_p (1 - 1/2 Re tr U_p)).
bhat=1 gives Var(q) itself (control).

Two q_x families (the literal formula has ONE plaquette per site; which
orientation is a convention, so both are reported):
  p01  : q_x = 1 - (1/2) Re tr U_{01}(x)        <- single plaquette, literal
  pall : q_x = 1 - (1/6) sum_{mu<nu} (1/2)Re tr U_{mu nu}(x)   <- orientation avg

Memory: no config is stored; only per-config per-block sums and sums of
squares (bhat=1 folds 4096 sites/config immediately). One job at a time.
Reuses the tested SU(2) heat bath (su2_heatbath + staple_sum2) from
research/lattice-cpb/w4_blocked.py, kappa = beta/2.
"""
import argparse
import sys
import time

import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
from w4_blocked import adj, random_su2, staple_sum2, su2_heatbath  # noqa: E402

FAMS = ('p01', 'pall')


def run(L, beta, ntherm, ncfg, gap, seed, bhats, log=print):
    rng = np.random.default_rng(seed)
    U = random_su2(rng, (4, L, L, L, L))
    ns = L ** 4
    g = np.meshgrid(*[np.arange(L)] * 4, indexing='ij')
    par = ((g[0] + g[1] + g[2] + g[3]) % 2).reshape(-1)
    idxp = [np.nonzero(par == c)[0] for c in (0, 1)]
    kappa = beta / 2.0
    Uf = U.reshape(4, ns, 2, 2)

    def sweep():
        for c in (0, 1):
            ids = idxp[c]
            for mu in range(4):
                S = staple_sum2(U, mu).reshape(ns, 2, 2)[ids]
                A = Uf[mu, ids]
                Uf[mu, ids] = su2_heatbath(kappa, A @ S, 0, 1, rng) @ A

    def qfields():
        """dict fam -> per-site q_x of shape (L,L,L,L)"""
        ps = []
        for mu in range(4):
            for nu in range(mu + 1, 4):
                P = (U[mu] @ np.roll(U[nu], -1, axis=mu)
                     @ adj(np.roll(U[mu], -1, axis=nu)) @ adj(U[nu]))
                ps.append(0.5 * np.real(np.trace(P, axis1=-2, axis2=-1)))
        return {'p01': 1.0 - ps[0], 'pall': 1.0 - sum(ps) / 6.0}

    def coarsen(q, b):
        nb = L // b
        return q.reshape(nb, b, nb, b, nb, b, nb, b).mean(axis=(1, 3, 5, 7))

    bhats = [b for b in bhats if L % b == 0]
    nblk = {b: (L // b) ** 4 for b in bhats}
    q1 = {(f, b): np.empty(ncfg) for f in FAMS for b in bhats}
    q2 = {(f, b): np.empty(ncfg) for f in FAMS for b in bhats}
    plaq = np.empty(ncfg)

    t0 = time.time()
    for _ in range(ntherm):
        sweep()
    log(f"# L={L} beta={beta} N={ncfg} gap={gap} therm={ntherm} seed={seed} "
        f"bhats={bhats}  fams={list(FAMS)}  therm_wall={time.time()-t0:.1f}s")

    t0 = time.time()
    for k in range(ncfg):
        for _ in range(gap):
            sweep()
        qs = qfields()
        plaq[k] = 1.0 - qs['pall'].mean()
        for f in FAMS:
            q = qs[f]
            for b in bhats:
                bm = coarsen(q, b).ravel()
                q1[f, b][k] = bm.sum()
                q2[f, b][k] = np.dot(bm, bm)
        if (k + 1) % 500 == 0:
            el = time.time() - t0
            log(f"  {k+1}/{ncfg}  wall={el:.0f}s  eta={el/(k+1)*(ncfg-k-1):.0f}s"
                f"  plaq_so_far={plaq[:k+1].mean():.4f}")
    wall = time.time() - t0

    res = {}
    for f in FAMS:
        for b in bhats:
            n = nblk[b]
            S1 = q1[f, b].sum()
            S2 = q2[f, b].sum()
            tot = ncfg * n
            theta = (S2 - S1 * S1 / tot) / tot
            th = np.empty(ncfg)
            for i in range(ncfg):
                s1 = S1 - q1[f, b][i]
                s2 = S2 - q2[f, b][i]
                m = (ncfg - 1) * n
                th[i] = (s2 - s1 * s1 / m) / m
            err = np.sqrt((ncfg - 1) / ncfg * np.sum((th - th.mean()) ** 2))
            # autocorrelation-robust check: 20 contiguous bins, spread of bin means
            nb = 20
            edges = np.linspace(0, ncfg, nb + 1).astype(int)
            tb = []
            for lo, hi in zip(edges[:-1], edges[1:]):
                s1 = q1[f, b][lo:hi].sum()
                s2 = q2[f, b][lo:hi].sum()
                m = (hi - lo) * n
                tb.append((s2 - s1 * s1 / m) / m)
            err_bin = np.std(tb, ddof=1) / np.sqrt(nb)
            res[f, b] = (b ** 4 * theta, b ** 4 * err, n, b ** 4 * err_bin)

    log(f"# wall={wall:.1f}s ({wall/ncfg:.2f}s/cfg)  plaq={plaq.mean():.4f}"
        f" +- {plaq.std()/np.sqrt(ncfg):.4f}  (target 0.652)")
    for f in FAMS:
        for b in bhats:
            v, e, n, eb = res[f, b]
            log(f"FAM {f} BHAT {b}  blocks={n}  varpi={v:.6f}  err_jk={e:.6f}"
                f"  err_bin={eb:.6f}")
    return plaq, res, wall


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=8)
    ap.add_argument('--beta', type=float, default=2.5)
    ap.add_argument('--N', type=int, default=2000)
    ap.add_argument('--therm', type=int, default=300)
    ap.add_argument('--gap', type=int, default=2)
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--bhats', default='1,2,4')
    a = ap.parse_args()
    run(a.L, a.beta, a.therm, a.N, a.gap, a.seed,
        [int(x) for x in a.bhats.split(',')])
