"""
W7 driver: exact marginal curvature G''(0), SU(2) and SU(3).

PART A: panel of transverse plane-wave gluon modes at the lowest lattice momenta
        (k = 2*pi*n/L), exact V_eff''(0) per mode from the equilibrium histogram.
PART B: the b=2 block-uniform softest mode (cond_var.py's object, whose Var was
        only proxied -> G''_soft ~ -0.165), now done EXACTLY: Hbar=P^T H P,
        softest eigenvector v*, exact V_eff''(0) along v*.
"""
import sys
import time
import itertools
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
import w4_blocked as W   # noqa: E402
import gpp_core as C     # noqa: E402

T2 = W.T2


# ======================================================================
# PART A : plane-wave panel, SU(2)
# ======================================================================
def plane_panel(L, nmax=2):
    modes = {}
    for axis in range(4):
        for pol in range(4):
            if pol == axis:
                continue
            for n in range(1, nmax + 1):
                for kind in ('cos', 'sin'):
                    # representative colours: one Cartan (2), one non-Cartan (0)
                    for color in (2, 0):
                        key = f"ax{axis}_p{pol}_n{n}_{kind}_c{color}"
                        modes[key] = C.plane_mode(L, axis, pol, color, n=n, kind=kind)
    return modes


def run_partA_su2(L, beta, N, seed=101, gapevery=3, ntherm=400, nmax=1):
    print(f"# PART A SU(2) L={L} beta={beta} N={N} nmax={nmax}", flush=True)
    t0 = time.time()
    cfgs = C.sample_configs(L, beta, ntherm, N, gapevery, seed=seed)
    print(f"  sampled {N} cfgs in {time.time()-t0:.1f}s", flush=True)
    modes = plane_panel(L, nmax=nmax)
    names = list(modes); ml = [modes[n] for n in names]
    acc = {n: ([], [], []) for n in names}
    t0 = time.time()
    for i, U in enumerate(cfgs):
        a, H, G = C.mode_obs(U, beta, ml, h=1e-3)
        for n, av, hv, gv in zip(names, a, H, G):
            acc[n][0].append(av); acc[n][1].append(hv); acc[n][2].append(gv)
        if (i + 1) % 1000 == 0:
            print(f"    {i+1}/{N} ({time.time()-t0:.0f}s)", flush=True)
    rng = np.random.default_rng(2024)
    rows = []
    print(f"\n  {'mode':26s} {'Gpp=kap':>9s} {'+-SE':>8s} {'1/Var':>7s} "
          f"{'E[H]':>7s} {'VarG':>7s} {'Ev0':>7s} {'lam4':>7s} {'sd(a)':>6s}")
    # save all raw per-mode data for re-analysis / width-scan
    np.savez(f'partA_su2_N{N}.npz',
             **{f'{n}_{k}': np.array(acc[n][i]) for n in names
                for i, k in enumerate(('a', 'H', 'G'))})
    for n in names:
        a, H, G = map(np.array, acc[n])
        kap, se, lam, mb, inv = C.veff_curvature(a, rng=rng)
        Ev, Ese, Ev0, EH, VG, nb = C.ev_estimator(a, H, G, rng=rng)
        rows.append((n, kap, se, inv, EH, VG, Ev0, lam, a.std()))
        print(f"  {n:26s} {kap:+9.4f} {se:8.4f} {inv:7.3f} "
              f"{EH:7.3f} {VG:7.3f} {Ev0:+7.3f} {lam:+7.3f} {a.std():6.3f}")
    return rows, acc


# ======================================================================
# PART B : block-uniform softest mode, SU(2)
# ======================================================================
def su2_lift(L, b=2):
    D = L // b
    rows, cols = [], []
    Nc = 4 * D ** 4 * 3
    cidx = lambda mu, xp: (((mu*D+xp[0])*D+xp[1])*D+xp[2])*D+xp[3]
    for mu in range(4):
        for xp in itertools.product(range(D), repeat=4):
            c0 = cidx(mu, xp) * 3
            for a in range(3):
                for d in itertools.product(range(b), repeat=4):
                    x = tuple(xp[k]*b + d[k] for k in range(4))
                    fi = (((mu*L+x[0])*L+x[1])*L*L + x[2]*L + x[3]) * 3 + a
                    rows.append(fi); cols.append(c0 + a)
    vals = np.ones(len(rows)) / np.sqrt(b ** 4)
    return sp.coo_matrix((vals, (rows, cols)),
                         shape=(4 * L ** 4 * 3, Nc)).tocsr()


def run_partB_su2(L, beta, N, nH=12, seed=77, gapevery=3, ntherm=400):
    print(f"\n# PART B SU(2) L={L} beta={beta} N={N} nH={nH}", flush=True)
    P = su2_lift(L, 2)
    Nc = P.shape[1]
    rng = np.random.default_rng(seed)
    lat = W.random_su2(rng, (4, L, L, L, L))
    # thermalize
    cfgs = W.sample2(L, beta, rng, ntherm, N, gapevery)
    # Hbar = E[P^T H P] from nH configs via matvec_F
    Hbar = np.zeros((Nc, Nc))
    t0 = time.time()
    for c in range(min(nH, len(cfgs))):
        U = cfgs[c]
        FP = np.empty((P.shape[0], Nc))
        for j in range(Nc):
            vfull = np.asarray(P[:, j].todense()).ravel()
            FP[:, j] = W.matvec_F(U, beta, vfull, T2, W.grad2, h=1e-5)
        Hb = (P.T @ FP)
        Hbar += 0.5 * (Hb + Hb.T)
    Hbar /= min(nH, len(cfgs))
    ev, evec = np.linalg.eigh(Hbar)
    print(f"  Hbar dim={Nc} built in {time.time()-t0:.1f}s: "
          f"min={ev[0]:+.4f}  lowest5={np.round(ev[:5],4)}", flush=True)
    vc = evec[:, 0]                       # softest coarse direction
    Vfull = np.asarray(P @ vc).reshape(4, L, L, L, L, 3)
    print(f"  ||Vfull||={np.linalg.norm(Vfull):.4f} (should be 1)", flush=True)
    # exact marginal along Vfull
    a, H, G = [], [], []
    for U in cfgs:
        av, hv, gv = C.mode_obs(U, beta, [Vfull], h=1e-3)
        a.append(av[0]); H.append(hv[0]); G.append(gv[0])
    a = np.array(a); H = np.array(H); G = np.array(G)
    rng2 = np.random.default_rng(2025)
    kap, se, lam, mb, inv = C.veff_curvature(a, rng=rng2)
    Ev, Ese, Ev0, EH, VG, nb = C.ev_estimator(a, H, G, rng=rng2)
    print(f"  EXACT softest blocked mode: Gpp_soft = {kap:+.4f} +- {se:.4f} "
          f"(1/Var={inv:+.4f})  E[H]={EH:.4f} VarG={VG:.4f} "
          f"Ev0={Ev0:+.4f} lam4={lam:+.4f}", flush=True)
    print(f"  [cond_var proxy was -0.165; bare E-part lam_E(Hbar) = {ev[0]:+.4f}]",
          flush=True)
    return dict(ev=ev, vc=vc, kap=kap, se=se, inv=inv, EH=EH, VG=VG, Ev0=Ev0,
                lam=lam, a=a, H=H, G=G)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--part', default='A')
    ap.add_argument('--L', type=int, default=4)
    ap.add_argument('--beta', type=float, default=2.5)
    ap.add_argument('--N', type=int, default=4000)
    ap.add_argument('--nH', type=int, default=12)
    ap.add_argument('--nmax', type=int, default=1)
    ap.add_argument('--seed', type=int, default=101)
    args = ap.parse_args()
    if args.part in ('A', 'all'):
        run_partA_su2(args.L, args.beta, args.N, seed=args.seed, nmax=args.nmax)
    if args.part in ('B', 'all'):
        run_partB_su2(args.L, args.beta, args.N, nH=args.nH, seed=args.seed + 7)
