"""Q13 — push G''(0) at L=2, SU(2), beta=6 from 1.3 sigma to >=5 sigma.

SAME model / sampler / slice discipline as o24_lat.py + p19_sweep.py (model A: one common
v = e_z/sqrt(NL) on every link, flat chart |x_l| < 2*pi on the SU(2) fundamental domain,
pair-transfer proposal that conserves <v,X> = 0 exactly, hard rejection r2n < DOM2),
only more statistics:  G=8 independent chains x ns=20000 sweeps  (p19 used 4 x 5000).

Observables per sweep (identical definitions to p19):
    kappa = d^2/da^2 V(X + a v) |_{a=0},   F = d/da V(X + a v) |_{a=0}
    G''(0) = E[kappa] - Var(F) = Delta_dyn + W_geo + X     (exact split, [T])

usage: py -3.12 q13_beta6.py [ns] [eps] [seed0] [burn] [G] [tag]
writes %TEMP%/q13_<tag>.npz  (repo gains only *.py / *.log)
"""
import os
import sys
import time

import numpy as np

from o24_lat import C0, C1, C2, C3, DOM2, NL, NPQ, ROWS, geyer_tau, ht, phi_r2, quat


def split_V(x, beta):
    w, a = quat(x)
    tr = ht(w, a, C0, C1, C2, C3)
    return beta * (NPQ - tr.sum(-1)), phi_r2((x * x).sum(-1)).sum(-1)


def run(beta, v, G, ns, eps, rng, h=1e-3, burn=800, prog=2000):
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
    t0 = time.time()
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
        if prog and sw % prog == 0:
            r = np.linalg.norm(x, axis=-1)
            print(f"   [sweep {sw:6d}/{ns + burn}] wall={time.time() - t0:7.1f}s "
                  f"acc={nacc.mean() / ((sw + 1) * NL):.3f} domrej={ndom / (G * (sw + 1) * NL):.4f} "
                  f"r_mean={r.mean():.3f}", flush=True)
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
    return cols, nacc / ((ns + burn) * NL), ndom / (G * (ns + burn) * NL), x, time.time() - t0


if __name__ == "__main__":
    ns = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    eps = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5809
    seed0 = int(sys.argv[3]) if len(sys.argv) > 3 else 700
    burn = int(sys.argv[4]) if len(sys.argv) > 4 else 2000
    G = int(sys.argv[5]) if len(sys.argv) > 5 else 8
    tag = sys.argv[6] if len(sys.argv) > 6 else 'b6'
    beta = 6.0
    v = np.zeros((NL, 3))
    v[:, 2] = 1.0 / np.sqrt(NL)        # mode A, |v|^2 = 1
    seeds = [seed0 + i for i in range(G)]
    print(f"[Q13] beta={beta:g} mode=A G={G} ns={ns} burn={burn} eps={eps} seeds={seeds}", flush=True)
    cols, acc, domrej, x, el = run(beta, v, G, ns, eps, np.random.default_rng(seeds), burn=burn)
    out = os.path.join(os.environ.get('TEMP', '.'), f"q13_{tag}.npz")
    np.savez(out, cols=cols, seeds=np.array(seeds), beta=beta, ns=ns, eps=eps, burn=burn,
             acc=acc, domrej=domrej, el=el)
    r = np.linalg.norm(x, axis=-1)
    print(f"[Q13] beta={beta:g} tag={tag} nsw={ns} burn={burn} eps={eps} G={G} seeds={seeds} "
          f"wall={el:.1f}s ({1e3 * el / (ns + burn):.2f} ms/sweep for {G} chains) -> {out}")
    print(f"[Q13] acc_link={np.round(acc, 3)} domain-reject={domrej:.4f} "
          f"final r: mean={r.mean():.3f} max={r.max():.3f} (2pi={2 * np.pi:.3f})")
    print(f"   {'g':>2} {'E[kappa]':>9} {'Var(F)':>8} {'Gpp(0)':>9} {'Delta_dyn':>10} "
          f"{'W_geo':>8} {'X':>8} {'sum-check':>10} {'E[F]':>8} {'tauK':>7} {'tauF':>6} {'Neff':>6}")
    for g in range(G):
        k, f, kW, fW, kP, fP = (cols[:, g, i] for i in range(6))
        dd = kW.mean() - fW.var()
        wg = kP.mean() - fP.var()
        xx = -2.0 * float(np.cov(fW, fP)[0, 1])
        gpp = k.mean() - f.var()
        tk, tf = geyer_tau(k, 2000), geyer_tau(f, 2000)
        print(f"   {g:>2} {k.mean():9.4f} {f.var():8.4f} {gpp:+9.4f} {dd:+10.4f} "
              f"{wg:+8.4f} {xx:+8.4f} {dd + wg + xx - gpp:+10.1e} {f.mean():+8.4f} "
              f"{tk:7.1f} {tf:6.1f} {ns / max(tk, tf):6.0f}")
