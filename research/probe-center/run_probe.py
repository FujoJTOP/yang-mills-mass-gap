"""Production probe: MCG fix + Z(3) projection + vortex + sector observables.

Usage:
  python run_probe.py --npz <config.npz> --label b2.6_L6 --rmax 3 [--nmax N]
Saves probe_<label>.npz with per-snapshot aggregates and prints a summary.
"""
import argparse
import sys
import time

import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/probe-center')
import mcg


def random_su3(L, rng):
    g = np.zeros((L, L, L, L, 3, 3), np.complex128)
    a = rng.normal(size=(L, L, L, L, 3, 3)) + 1j * rng.normal(size=(L, L, L, L, 3, 3))
    q, r = np.linalg.qr(a)
    d = np.einsum('...ii->...i', r)
    d = d / np.abs(d)
    q = q * d[..., None, :]                       # uniform-ish SU(3)
    det = np.linalg.det(q)
    q = q / (det ** (1.0 / 3.0))[..., None, None]  # det -> 1
    g[...] = q
    return g


def gauge_fix(U, L, nsw, maxit, rng):
    """Two restarts (identity + random); keep the one with larger R."""
    best = None
    for kind in ('id', 'rand'):
        if kind == 'id':
            g = np.zeros((L, L, L, L, 3, 3), np.complex128)
            for i in range(3):
                g[..., i, i] = 1.0
        else:
            g = random_su3(L, rng)
        for sw in range(nsw):
            mcg.mcg_sweep(U, g, L, maxit)
        V = mcg.apply_gauge(U, g, L)
        R = mcg.center_functional(V, L)
        if best is None or R > best[0]:
            best = (R, V, g)
    return best


def analyze_config(U, L, rmax, nsw, maxit, rng):
    R, V, g = gauge_fix(U, L, nsw, maxit, rng)
    n, amb = mcg.center_project(V, L)
    flux = mcg.plaquette_flux(n, L)
    plaqs, labels, sizes = mcg.vortex_clusters(flux, L)
    ncl = len(sizes)
    k1 = sum(1 for p in plaqs if p[3] == 1)
    k2 = sum(1 for p in plaqs if p[3] == 2)
    # per-cluster: N1-N2 and size
    cq = np.zeros(ncl, np.int64)
    for p, (x, mu, nu, k) in enumerate(plaqs):
        cq[labels[p]] += (1 if k == 1 else -1)
    sizes_sorted = np.sort(sizes)[::-1]

    # fundamental / adjoint Wilson loops
    Tf, Tf2 = mcg.wilson_loop_stats(U, L, rmax)
    Uz = mcg.expand_center_links(n, L)
    Tp, Tp2 = mcg.wilson_loop_stats(Uz, L, rmax)

    # Polyakov correlators
    pff = mcg.polyakov_traces(U, L)
    pfp = mcg.polyakov_traces(Uz, L)
    rp = min(rmax, L - 1)
    Cf = mcg.polyakov_corr(pff, L, rp)          # fundamental (charged)
    CA = mcg.adjoint_polyakov_corr(pff, L, rp)  # adjoint (neutral)
    Cfp = mcg.polyakov_corr(pfp, L, rp)
    CAp = mcg.adjoint_polyakov_corr(pfp, L, rp)

    return dict(
        R=R, amb=amb, nvortex=len(plaqs), nclust=ncl, k1=k1, k2=k2,
        sizes=sizes_sorted, sizes_raw=sizes, cq=cq,
        density=len(plaqs) / (6.0 * L ** 4),
        Tf=Tf, Tf2=Tf2, Tp=Tp, Tp2=Tp2,
        Cf=Cf, CA=CA, Cfp=Cfp, CAp=CAp,
        poly_mean_full=float(np.abs(pff).mean() / 3.0),
        poly_mean_proj=float(np.abs(pfp).mean() / 3.0),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--npz', required=True)
    ap.add_argument('--label', required=True)
    ap.add_argument('--rmax', type=int, default=3)
    ap.add_argument('--nsw', type=int, default=30)
    ap.add_argument('--maxit', type=int, default=6)
    ap.add_argument('--nmax', type=int, default=0)
    ap.add_argument('--seed', type=int, default=1234)
    args = ap.parse_args()

    d = np.load(args.npz, allow_pickle=True)
    Uall = d['snapshots']
    L = int(d['L'])
    beta = float(d['beta'])
    ns = Uall.shape[0] if args.nmax <= 0 else min(args.nmax, Uall.shape[0])
    print(f"{args.label}: L={L} beta={beta} nsnap={ns}", flush=True)

    rng = np.random.default_rng(args.seed)
    recs = []
    t0 = time.time()
    for i in range(ns):
        rec = analyze_config(Uall[i].copy(), L, args.rmax, args.nsw, args.maxit, rng)
        recs.append(rec)
        if (i + 1) % 5 == 0 or i == 0:
            print(f"  [{i+1}/{ns}] R={rec['R']:.3f} dens={rec['density']:.3f} "
                  f"ncl={rec['nclust']} k1={rec['k1']} k2={rec['k2']} "
                  f"Wf11={rec['Tf'][1,1].real:.4f} Wp11={rec['Tp'][1,1].real:.4f} "
                  f"t={time.time()-t0:.0f}s", flush=True)

    def stack(key):
        return np.array([r[key] for r in recs])

    # cluster-size histogram across all configs (finite clusters)
    allsizes = np.concatenate([r['sizes'] for r in recs]) if recs else np.array([])

    np.savez(f"probe_{args.label}.npz",
             L=L, beta=beta, nsnap=ns, rmax=args.rmax,
             R=stack('R'), amb=stack('amb'), nvortex=stack('nvortex'),
             nclust=stack('nclust'), k1=stack('k1'), k2=stack('k2'),
             density=stack('density'),
             Tf=np.array([r['Tf'] for r in recs]),
             Tf2=np.array([r['Tf2'] for r in recs]),
             Tp=np.array([r['Tp'] for r in recs]),
             Tp2=np.array([r['Tp2'] for r in recs]),
             Cf=np.array([r['Cf'] for r in recs]),
             CA=np.array([r['CA'] for r in recs]),
             Cfp=np.array([r['Cfp'] for r in recs]),
             CAp=np.array([r['CAp'] for r in recs]),
             poly_mean_full=stack('poly_mean_full'),
             poly_mean_proj=stack('poly_mean_proj'),
             allsizes=allsizes,
             cq_all=np.concatenate([r['cq'] for r in recs]) if recs else np.array([]),
             cl_size=np.concatenate([r['sizes_raw'] for r in recs]) if recs else np.array([]),
             )
    print(f"saved probe_{args.label}.npz  ({time.time()-t0:.0f}s)")


if __name__ == '__main__':
    main()
