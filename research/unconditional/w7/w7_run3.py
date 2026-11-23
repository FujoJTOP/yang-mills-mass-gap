"""SU(3) driver: exact marginal curvature G''(0) for transverse plane modes.

Branch caveat: the SU(3) principal log has a cut at eigenvalue angle pi.  We
record, per config, max|eigangle| and report the curvature BOTH on all configs
and on configs with max|angle| < pi - 0.15 (branch-safe), to test sensitivity.
The force H_vv, G_v are computed from the action along the ray (log-free).
"""
import sys
import time
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
import su3 as S3    # noqa: E402
import gpp_core3 as G3   # noqa: E402
import gpp_core as C     # noqa: E402


def max_eigangle(U):
    ang = np.abs(np.angle(np.linalg.eigvals(U)))
    return float(ang.max())


def run(L=4, beta=5.7, N=1500, seed=101, gapevery=3, ntherm=200,
        colors=(0, 3, 7), nmax=1):
    print(f"# SU(3) L={L} beta={beta} N={N}", flush=True)
    t0 = time.time()
    cfgs = G3.sample_cfgs(L, beta, ntherm, N, gapevery, seed=seed)
    print(f"  sampled {N} cfgs in {time.time()-t0:.1f}s; "
          f"mean plaquette={np.mean([S3.mean_plaquette(U) for U in cfgs[:50]]):.4f}",
          flush=True)
    modes = G3.panel(L, colors=colors, nmax=nmax)
    names = list(modes); ml = [modes[n] for n in names]
    acc = {n: ([], [], []) for n in names}
    ang = []
    t0 = time.time()
    for i, U in enumerate(cfgs):
        ang.append(max_eigangle(U))
        a, H, G = G3.mode_obs(U, beta, ml, h=1e-3)
        for n, av, hv, gv in zip(names, a, H, G):
            acc[n][0].append(av); acc[n][1].append(hv); acc[n][2].append(gv)
        if (i + 1) % 500 == 0:
            print(f"    {i+1}/{N} ({time.time()-t0:.0f}s)", flush=True)
    ang = np.array(ang)
    safe = ang < (np.pi - 0.15)
    print(f"  branch-safe configs: {safe.sum()}/{N} "
          f"(dropped {N-safe.sum()} with eigangle>pi-0.15)", flush=True)
    np.savez(f'partA_su3_N{N}.npz',
             **{f'{n}_{k}': np.array(acc[n][i]) for n in names
                for i, k in enumerate(('a', 'H', 'G'))}, ang=ang)
    rng = np.random.default_rng(4242)
    print(f"\n  {'mode':26s} {'kap_all':>9s} {'+-SE':>7s} {'kap_safe':>9s} "
          f"{'1/Var':>7s} {'E[H]':>7s} {'VarG':>7s}")
    rows = []
    for n in names:
        a, H, G = map(np.array, acc[n])
        kap, se, lam, mb, inv = C.veff_curvature(a, rng=rng)
        ks, ses, lams, _, invs = C.veff_curvature(a[safe], rng=rng)
        rows.append((n, kap, se, ks, inv, H.mean(), G.var()))
        print(f"  {n:26s} {kap:+9.4f} {se:7.4f} {ks:+9.4f} "
              f"{inv:7.3f} {H.mean():7.3f} {G.var():7.3f}")
    return rows, acc, safe


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=4)
    ap.add_argument('--beta', type=float, default=5.7)
    ap.add_argument('--N', type=int, default=1500)
    ap.add_argument('--seed', type=int, default=101)
    ap.add_argument('--nmax', type=int, default=1)
    args = ap.parse_args()
    run(args.L, args.beta, args.N, seed=args.seed, nmax=args.nmax)
