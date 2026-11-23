"""
KZERO driver: exact marginal curvature G''(0) at the SOFTEST nonzero mode
k = 2*pi*n/L (n=1) on larger volumes, SU(2) Wilson at beta=2.5.

Reuses rerun_run.measure / rerun_core (flat exponential chart) exactly:
  G''(0) = ( E[H_vv | a=0] - Var(G_v | a=0) ) + W
  W      = E[Phi''] - Var(Phi'|0) - 2 Cov(G_v,Phi'|0)
  gold-standard kappa_kde = curvature of -log P(a) at a=0.

Saves per-mode a,H,G,P1,P2 so the analysis can be redone offline.
"""
import sys
import time
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
import w4_blocked as W   # noqa: E402
import gpp_core as C     # noqa: E402
import rerun_core as RC  # noqa: E402
import rerun_run as R    # noqa: E402


def mode_spec(nmax=1, extra_colors=False):
    """single-axis standing waves, n=nmax, transverse polarisations.
    All modes have |k| = 2*pi*nmax/L along one axis (softest nonzero)."""
    spec = [(1, 2, 2), (1, 3, 2), (0, 2, 2), (2, 3, 2), (3, 1, 2), (1, 2, 0)]
    if extra_colors:
        spec += [(0, 1, 2), (2, 1, 2), (3, 2, 2), (0, 3, 2), (1, 0, 2), (3, 0, 0)]
    return spec


def run(L, beta, N, therm, gap, seed, nmax=1, extra_colors=False,
        out=None, log=print):
    nlist = nmax if isinstance(nmax, (list, tuple)) else [nmax]
    spec = mode_spec(1, extra_colors)
    modes, names, ntags, kvals = [], [], [], []
    for n in nlist:
        for (ax, pl, c) in spec:
            modes.append(C.plane_mode(L, ax, pl, c, n=n, kind='cos'))
            names.append(f'n{n}_ax{ax}_p{pl}_c{c}')
            ntags.append(n)
            kvals.append(2.0 * np.pi * n / L)
    ntags = np.array(ntags); kvals = np.array(kvals)
    log(f"# KZERO SU2 L={L} beta={beta} N={N} therm={therm} gap={gap} "
        f"nmax={nlist} k={np.round(kvals,3)} modes={len(modes)}")
    t0 = time.time()
    cfgs = R.sample_su2(L, beta, therm, N, gap, seed)
    log(f"  sampled {len(cfgs)} cfgs in {time.time()-t0:.1f}s "
        f"(plaq={np.mean([W.action2(U,beta) for U in cfgs[:20]]):.4f})")
    acc = {k: [] for k in ('a', 'H', 'G', 'P1', 'P2')}
    t0 = time.time()
    for i, U in enumerate(cfgs):
        r = R.measure(U, beta, modes, 'su2')
        for k in acc:
            acc[k].append(r[k])
        if (i + 1) % 200 == 0:
            log(f"    {i+1}/{N} ({time.time()-t0:.0f}s)")
    for k in acc:
        acc[k] = np.array(acc[k])            # (Ncfg, Nmode)
    if out is None:
        out = f'kzero_su2_L{L}_n{"_".join(map(str, nlist))}.npz'
    np.savez(out, **acc, names=np.array(names), ntags=ntags, kvals=kvals,
             L=L, beta=beta, N=N, gap=gap, therm=therm, seed=seed)
    log(f"  saved {out}  shape a={acc['a'].shape}")
    return out


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=8)
    ap.add_argument('--beta', type=float, default=2.5)
    ap.add_argument('--N', type=int, default=1200)
    ap.add_argument('--therm', type=int, default=500)
    ap.add_argument('--gap', type=int, default=4)
    ap.add_argument('--seed', type=int, default=101)
    ap.add_argument('--nmax', default='1')
    ap.add_argument('--extra', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    nmax = [int(x) for x in str(a.nmax).split(',')]
    run(a.L, a.beta, a.N, a.therm, a.gap, a.seed, nmax, a.extra, a.out)
