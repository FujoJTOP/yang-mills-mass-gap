"""P19 — fixed L=2, SU(2): beta-scan of G''(0) = E[kappa|0] - Var(F|0).

Same model / sampler / slice discipline as o24_lat.py (flat chart on the SU(2) fundamental
domain |x_l| < 2*pi, pair-transfer proposal preserving <v,X>=0 exactly), reusing its tables and
vectorised plaquette helpers.  Added here: the exact split

    G''(0) = kappa.mean() - F.var() = Delta_dyn + W_geo + X
    Delta_dyn = E[kappa_W] - Var(F_W)      (Wilson bracket,  V = beta*sum(1-1/2ReTr U_P))
    W_geo     = E[kappa_Phi] - Var(F_Phi)  (Jacobian block,  Phi = -log J)
    X         = -2*Cov(F_W, F_Phi)

so one pass gives G''(0), Delta_dyn and Var(F) (= Var(F_W+F_Phi)) on the SAME configurations.

usage: py -3.12 p19_sweep.py <beta> <nsweep> <eps> [epsg-ignored] <seed0> <burn> <tag>
writes %TEMP%/p19_<tag>.npz  (repo gains only *.py / *.log)  and prints the per-chain table.
"""
import os
import sys
import time

import numpy as np

from o24_lat import C0, C1, C2, C3, DOM2, NL, NPQ, ROWS, ht, phi_r2, quat, geyer_tau


def split_V(x, beta):
    """(S_W, Phi) separately, from flat coordinates."""
    w, a = quat(x)
    tr = ht(w, a, C0, C1, C2, C3)
    return beta * (NPQ - tr.sum(-1)), phi_r2((x * x).sum(-1)).sum(-1)


def run(beta, v, G, ns, eps, rng, h=1e-3, burn=800):
    assert np.allclose(v, v[0]), "pair-transfer assumes one common v per link"
    x = rng.normal(0, 0.55, (G, NL, 3))
    x -= np.einsum('gla,la->g', x, v)[:, None, None] * v[None]
    w, a = quat(x)
    tr = ht(w, a, C0, C1, C2, C3)
    phS = phi_r2((x * x).sum(-1)).sum(-1)
    ar = np.arange(G)
    nacc = np.zeros(G)
    ndom = 0
    cols = np.empty((ns, G, 6))          # kappa, F, kappa_W, F_W, kappa_Phi, F_Phi
    for sw in range(ns + burn):
        for _ in range(NL):
            l = rng.integers(NL, size=G)
            l2 = (l + rng.integers(1, NL, size=G)) % NL
            wv = rng.uniform(-eps, eps, (G, 3))
            xn = x.copy()
            xn[ar, l] += wv
            xn[ar, l2] -= wv
            r2n = (xn * xn).sum(-1)
            phn = phi_r2(r2n).sum(-1)
            wn, an = quat(xn)
            rows = np.concatenate([ROWS[l], ROWS[l2]], axis=1)
            new = ht(wn, an, C0[rows], C1[rows], C2[rows], C3[rows])
            old = tr[ar[:, None], rows]
            dS = beta * (old - new).sum(-1) + (phn - phS)
            m = (dS <= 0) | (rng.random(G) < np.exp(-np.maximum(dS, 0.0)))
            dom = (r2n < DOM2).all(-1)
            ndom += int((~dom).sum())
            m &= dom
            nacc += m
            if m.any():
                x[m] = xn[m]
                w[m] = wn[m]
                a[m] = an[m]
                phS[m] = phn[m]
                trn = tr.copy()
                trn[ar[:, None], rows] = new
                tr[m] = trn[m]
        if sw < burn:
            continue
        j = sw - burn
        SW0 = beta * (NPQ - tr.sum(-1))
        PH0 = phS
        SWp, PHp = split_V(x + h * v[None], beta)
        SWm, PHm = split_V(x - h * v[None], beta)
        V0, Vp, Vm = SW0 + PH0, SWp + PHp, SWm + PHm
        cols[j, :, 0] = (Vp - 2 * V0 + Vm) / h ** 2
        cols[j, :, 1] = (Vp - Vm) / (2 * h)
        cols[j, :, 2] = (SWp - 2 * SW0 + SWm) / h ** 2
        cols[j, :, 3] = (SWp - SWm) / (2 * h)
        cols[j, :, 4] = (PHp - 2 * PH0 + PHm) / h ** 2
        cols[j, :, 5] = (PHp - PHm) / (2 * h)
    return cols, nacc / ((ns + burn) * NL), ndom / (G * (ns + burn) * NL), x


if __name__ == "__main__":
    beta = float(sys.argv[1])
    ns = int(sys.argv[2])
    eps = float(sys.argv[3])
    seed0 = int(sys.argv[5]) if len(sys.argv) > 5 else 400
    burn = int(sys.argv[6]) if len(sys.argv) > 6 else 800
    tag = sys.argv[7] if len(sys.argv) > 7 else 'b'
    G = 4
    v = np.zeros((NL, 3))
    v[:, 2] = 1.0 / np.sqrt(NL)
    seeds = [seed0 + i for i in range(G)]
    t0 = time.time()
    cols, acc, domrej, x = run(beta, v, G, ns, eps, np.random.default_rng(seeds), burn=burn)
    el = time.time() - t0
    out = os.path.join(os.environ.get('TEMP', '.'), f"p19_{tag}.npz")
    np.savez(out, cols=cols, seeds=np.array(seeds), beta=beta, ns=ns, eps=eps,
             acc=acc, domrej=domrej, el=el)
    r = np.linalg.norm(x, axis=-1)
    print(f"[P19] beta={beta:g} tag={tag} nsw={ns} burn={burn} eps={eps} G={G} seeds={seeds} "
          f"wall={el:.1f}s ({1e3*el/(ns+burn):.2f} ms/sweep for {G} chains) -> {out}")
    print(f"[P19] acc_link={np.round(acc,3)} domain-reject={domrej:.4f} "
          f"final r: mean={r.mean():.3f} max={r.max():.3f} (2pi={2*np.pi:.3f})")
    print(f"   {'g':>2} {'E[kappa]':>9} {'Var(F)':>8} {'Gpp(0)':>9} {'Delta_dyn':>10} "
          f"{'W_geo':>8} {'X':>8} {'sum-check':>10} {'E[F]':>8} {'tauK':>6} {'tauF':>6} {'Neff':>6}")
    for g in range(G):
        k, f, kW, fW, kP, fP = (cols[:, g, i] for i in range(6))
        dd = kW.mean() - fW.var()
        wg = kP.mean() - fP.var()
        xx = -2.0 * float(np.cov(fW, fP)[0, 1])
        gpp = k.mean() - f.var()
        tk, tf = geyer_tau(k), geyer_tau(f)
        print(f"   {g:>2} {k.mean():9.4f} {f.var():8.4f} {gpp:+9.4f} {dd:+10.4f} "
              f"{wg:+8.4f} {xx:+8.4f} {dd+wg+xx-gpp:+10.1e} {f.mean():+8.4f} "
              f"{tk:6.1f} {tf:6.1f} {ns/max(tk,tf):6.0f}")
