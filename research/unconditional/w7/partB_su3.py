"""PART B on SU(3): exact marginal curvature of the b=2 block-uniform softest
blocked mode (cond_var.py's SU(3) object; proxy gave G''_soft ~ -0.165).
Hbar = E[P^T H P] built with matvec_F (SU(3) gradient), softest eigenvector v*,
exact V_eff''(0) along Vfull=P v* from the equilibrium histogram.
"""
import sys
import time
import itertools
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
import su3 as S3          # noqa: E402
import w4_blocked as W    # noqa: E402
import gpp_core as C      # noqa: E402
import gpp_core3 as G3    # noqa: E402


def su3_lift(L, b=2):
    D = L // b
    rows, cols = [], []
    Nc = 4 * D ** 4 * 8
    cidx = lambda mu, xp: (((mu*D+xp[0])*D+xp[1])*D+xp[2])*D+xp[3]
    for mu in range(4):
        for xp in itertools.product(range(D), repeat=4):
            c0 = cidx(mu, xp) * 8
            for a in range(8):
                for d in itertools.product(range(b), repeat=4):
                    x = tuple(xp[k]*b + d[k] for k in range(4))
                    fi = (((mu*L+x[0])*L+x[1])*L*L + x[2]*L + x[3]) * 8 + a
                    rows.append(fi); cols.append(c0 + a)
    vals = np.ones(len(rows)) / np.sqrt(b ** 4)
    return sp.coo_matrix((vals, (rows, cols)), shape=(4 * L ** 4 * 8, Nc)).tocsr()


def run(L=4, beta=5.7, N=1500, nH=12, seed=505):
    print(f"# PART B SU(3) L={L} beta={beta} N={N} nH={nH}", flush=True)
    P = su3_lift(L, 2); Nc = P.shape[1]
    cfgs = G3.sample_cfgs(L, beta, 200, N, 3, seed=seed)
    Hbar = np.zeros((Nc, Nc))
    t0 = time.time()
    for c in range(min(nH, len(cfgs))):
        U = cfgs[c]
        FP = np.empty((P.shape[0], Nc))
        for j in range(Nc):
            vf = np.asarray(P[:, j].todense()).ravel()
            FP[:, j] = W.matvec_F(U, beta, vf, S3.T_GEN, S3.action_gradient, h=1e-5)
        Hb = P.T @ FP
        Hbar += 0.5 * (Hb + Hb.T)
    Hbar /= min(nH, len(cfgs))
    ev, evec = np.linalg.eigh(Hbar)
    print(f"  Hbar dim={Nc} in {time.time()-t0:.1f}s: min={ev[0]:+.4f} "
          f"lowest5={np.round(ev[:5],4)}", flush=True)
    vc = evec[:, 0]
    Vfull = np.asarray(P @ vc).reshape(4, L, L, L, L, 8)
    a, H, G = [], [], []
    for U in cfgs:
        av, hv, gv = G3.mode_obs(U, beta, [Vfull], h=1e-3)
        a.append(av[0]); H.append(hv[0]); G.append(gv[0])
    a = np.array(a); H = np.array(H); G = np.array(G)
    rng = np.random.default_rng(2026)
    kap, se, lam, _, inv = C.veff_curvature(a, rng=rng)
    Ev, Ese, Ev0, EH, VG, nb = C.ev_estimator(a, H, G, rng=rng)
    print(f"  EXACT softest blocked (SU3): Gpp = {kap:+.4f} +- {se:.4f} "
          f"(1/Var={inv:+.4f})  E[H]={EH:.4f} VarG={VG:.4f} Ev0={Ev0:+.4f} "
          f"lam4={lam:+.4f}", flush=True)
    print(f"  [cond_var proxy -= -0.165; bare E-part lam_E(Hbar) = {ev[0]:+.4f}]"
          , flush=True)
    # branch-sensitivity: recompute dropping |a| outliers beyond 3 sd
    m = np.abs(a) < 3 * a.std()
    k2, s2, *_ = C.veff_curvature(a[m], rng=rng)
    print(f"  branch/outlier-trimmed (|a|<3sd): Gpp={k2:+.4f}+-{s2:.4f} "
          f"({m.sum()}/{N})", flush=True)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=4)
    ap.add_argument('--beta', type=float, default=5.7)
    ap.add_argument('--N', type=int, default=1500)
    ap.add_argument('--nH', type=int, default=12)
    ap.add_argument('--seed', type=int, default=505)
    a = ap.parse_args()
    run(a.L, a.beta, a.N, a.nH, a.seed)
