"""Production driver: one (L, beta) chain, heat-bath, measurements.

Usage:
  python run.py --L 6 --beta 1.7 --sweeps 20000 --therm 2000 --meas 2 --seed 1
"""
import argparse
import json
import os
import time

import numpy as np

import su3 as S


def block_mean(P, b):
    L = P.shape[0]
    nb = L // b
    return P.reshape(nb, b, nb, b, nb, b, nb, b).mean(axis=(1, 3, 5, 7))


def polyakov_corr(P, rmax):
    """P: (L,L,L) complex field.  Returns C(r) averaged over x and 3 directions."""
    L = P.shape[0]
    out = np.zeros(rmax + 1)
    npow = 0
    for r in range(rmax + 1):
        acc = 0.0
        n = 0
        for d in range(3):
            acc += np.real(P * np.roll(np.conj(P), -r, axis=d)).sum()
            n += P.size
        out[r] = acc / n
    return out


def measure(lat, want_loops, rmax, block_sizes, thr):
    U = lat.U
    tr = S.plaquette_traces(U)
    allp = np.stack([a.reshape(-1) for a in tr.values()])
    rec = {}
    rec['plaq'] = float(allp.mean())
    dens = S.magnetic_density(U)
    rec['mag_slice'] = dens.sum(axis=(1, 2, 3))
    P = np.mean(np.stack([a for a in tr.values()]), axis=0)
    for b in block_sizes:
        O = block_mean(P, b).reshape(-1)
        rec[f'blk{b}'] = np.array([O.sum(), (O ** 2).sum(), (O ** 3).sum(),
                                   (O ** 4).sum(), O.size, O.min(), O.max()])
    rec['defect'] = np.array([np.mean(allp < t) for t in thr])
    rec['plaq_min'] = float(allp.min())
    rec['plaq_max'] = float(allp.max())
    pol = S.polyakov_field(U)
    rec['poly_corr'] = polyakov_corr(pol, rmax)
    rec['poly_mean'] = np.array([pol.real.mean(), pol.imag.mean()])
    if want_loops:
        L = U.shape[1]
        loops = np.zeros((rmax, rmax))
        for R in range(1, rmax + 1):
            for T in range(1, rmax + 1):
                loops[R - 1, T - 1] = S.wilson_loop(U, R, T)
        rec['loops'] = loops
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=6)
    ap.add_argument('--beta', type=float, required=True)
    ap.add_argument('--sweeps', type=int, default=5000)
    ap.add_argument('--therm', type=int, default=1000)
    ap.add_argument('--meas', type=int, default=1)
    ap.add_argument('--loop_every', type=int, default=10)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--cold', type=int, default=1)
    ap.add_argument('--out', type=str, required=True)
    ap.add_argument('--snap_every', type=int, default=0)
    ap.add_argument('--snap_max', type=int, default=32)
    ap.add_argument('--rmax', type=int, default=0)
    args = ap.parse_args()

    L = args.L
    rmax = args.rmax if args.rmax else min(3, L // 2)
    block_sizes = [b for b in (1, 2, 3, 4) if L % b == 0]
    thr = np.array([0.0, 0.1, 0.2, 1.0 / 3.0])

    rng = np.random.default_rng(args.seed)
    lat = S.Lattice(L, args.beta, rng, cold=bool(args.cold))

    t0 = time.time()
    plaq_hist = []
    for k in range(args.therm):
        lat.sweep()
        if k % max(1, args.therm // 10) == 0:
            plaq_hist.append((k, float(S.mean_plaquette(lat.U))))
    t_therm = time.time() - t0

    recs = []
    snaps = []
    t0 = time.time()
    nmeas = 0
    for k in range(args.sweeps):
        lat.sweep()
        if (k + 1) % args.meas == 0:
            loops = (k + 1) % args.loop_every == 0
            rec = measure(lat, loops, rmax, block_sizes, thr)
            rec['sweep'] = k + 1
            recs.append(rec)
            nmeas += 1
        if args.snap_every and (k + 1) % args.snap_every == 0 and len(snaps) < args.snap_max:
            snaps.append(lat.U.copy())
        if (k + 1) % max(1, args.sweeps // 20) == 0:
            el = time.time() - t0
            print(f"  [{k+1}/{args.sweeps}] {el:.1f}s  "
                  f"P={S.mean_plaquette(lat.U):.5f}  "
                  f"eff={lat.sampler_efficiency():.3f}", flush=True)
    t_prod = time.time() - t0

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out = {
        'L': L, 'beta': args.beta, 'seed': args.seed, 'sweeps': args.sweeps,
        'therm': args.therm, 'meas': args.meas, 'loop_every': args.loop_every,
        'rmax': rmax, 'block_sizes': np.array(block_sizes), 'thr': thr,
        'plaq': np.array([r['plaq'] for r in recs]),
        'mag_slice': np.array([r['mag_slice'] for r in recs]),
        'poly_corr': np.array([r['poly_corr'] for r in recs]),
        'poly_mean': np.array([r['poly_mean'] for r in recs]),
        'defect': np.array([r['defect'] for r in recs]),
        'plaq_min': np.array([r['plaq_min'] for r in recs]),
        'plaq_max': np.array([r['plaq_max'] for r in recs]),
        'sweep': np.array([r['sweep'] for r in recs]),
        'loop_sweep': np.array([r['sweep'] for r in recs if 'loops' in r]),
        'loops': np.array([r['loops'] for r in recs if 'loops' in r]),
        'therm_trace': np.array(plaq_hist),
        't_therm': t_therm, 't_prod': t_prod,
        'sampler_eff': lat.sampler_efficiency(),
    }
    for b in block_sizes:
        out[f'blk{b}'] = np.array([r[f'blk{b}'] for r in recs])
    if snaps:
        out['snapshots'] = np.array(snaps)
    np.savez(args.out, **out)
    print(f"saved {args.out}: {nmeas} measurements, therm {t_therm:.1f}s, "
          f"prod {t_prod:.1f}s, sampler eff {lat.sampler_efficiency():.3f}")


if __name__ == '__main__':
    main()
