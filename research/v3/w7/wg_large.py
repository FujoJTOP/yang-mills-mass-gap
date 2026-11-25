#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
wg_large.py -- extend the WG gauge-invariant estimator (WG-INVARIANT-ESTIMATOR.md) to L=12
(and L=16 if the measured L=12 sweep rate projects < 20 min).

Object (identical estimator, same measure, same protocol as wg_invariant.py):
    O_s(x) = sum_{mu<nu} ( 1 - (1/2) Re tr P_{mu nu}^s(x) ),   P^s = plaquette of APE-smeared links
    q      = V^{-1/2} sum_x ( O_s(x) - <O_s> ),   V = L^4
    chi_O  = Var(q) = V * Var_c( Obar ),          Obar = (1/V) sum_x O_s(x)
Total-normalized (non-self-averaging), exactly gauge-invariant, APE ns=2 alpha=0.5 (fixed lattice radius,
SAME as L=4,6,8 so the points are directly comparable).  Reports chi_O, contact C_s(0), chi_phys,
chi_bare (ns=0 reference), chi_flat and chi_P references, jackknife SEs.

Protocol matches wg_invariant.py: therm=500, gap=8, N in 300..500, one process (charter s4: never run
two heavy SU(N) L>=6 jobs at once).  L=16 is run only if the L=12 measured s/sweep projects the full
L=16 job to < 1200 s (20 min); otherwise it is recorded as a projection and skipped.

Run:  py -3.12 -I research/v3/w7/wg_large.py            # L=12 N=400 -> maybe L=16 N=300
      py -3.12 -I research/v3/w7/wg_large.py 300 0      # N12=300, max16_min=0 (skip L=16)
Pure numpy/numba; no downloads; ONE process.  Labels [V]/[O].  No QED.  Clay 0%.
"""
import sys, os, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CPB = 'D:/Dev/math/research/lattice-cpb'
W7U = 'D:/Dev/math/research/unconditional/w7'
for p in (HERE, CPB, W7U):
    if p not in sys.path:
        sys.path.insert(0, p)

import w4_blocked as W            # ape_all2, random_su2, adj
import m1_fast as F               # Numba heat-bath sweep, random_start, chi_flat
import wg_invariant as G          # scalar_density, polyakov_spatial, jk_var_se, tau_int

RESULT = os.path.join(HERE, 'wg_large_results.json')

# Existing L=4,6,8 anchors (wg_invariant_results.json) -- for the cross-volume trend.
ANCHORS = {
    'L4': dict(chi_O=0.32765014859553365, chi_O_se=0.0222733807527846,
               chi_phys=0.2974008604510313, N=800),
    'L6': dict(chi_O=0.37743611373983316, chi_O_se=0.02133943765327901,
               chi_phys=0.3431821313604035, N=400),
    'L8': dict(chi_O=0.4104591281817326, chi_O_se=0.057637361481588344,
               chi_phys=0.3754055557778887, N=160),
}


def _dump(out):
    with open(RESULT, 'w') as fh:
        json.dump(out, fh, indent=1,
                  default=lambda o: float(o) if isinstance(o, (np.floating, np.ndarray)) else str(o))


def run_volume(L, beta, N, therm, gap, seed, ns_ape=2, ckpt_every=100, out=None, key=None):
    """One volume, one process.  Times therm and measurement separately for the L=16 projection.
    Checkpoints the raw per-config series to `out` every ckpt_every configs (incremental JSON)."""
    np.random.seed(seed)
    U = F.random_start(L, seed)
    F.sweep(U, beta, L)                       # JIT warm-up (cache=True), not timed

    t0 = time.time()
    for _ in range(therm):
        F.sweep(U, beta, L)
    t_therm = time.time() - t0
    per_sweep = t_therm / max(1, therm)       # measured s/sweep at this L

    Obar = np.empty(N); Osq = np.empty(N); Obar_bare = np.empty(N)
    chif = np.empty(N); Pbar = np.empty(N)
    meas_time = 0.0
    for i in range(N):
        for _ in range(gap):
            F.sweep(U, beta, L)
        tm0 = time.time()
        Us = W.ape_all2(U, alpha=0.5, nsteps=ns_ape)
        Os = G.scalar_density(Us)
        Obar[i] = Os.mean(); Osq[i] = (Os * Os).mean()
        Obar_bare[i] = G.scalar_density(U).mean()
        chif[i] = F.chi_flat(U, L)
        Pbar[i] = G.polyakov_spatial(U).mean()
        meas_time += time.time() - tm0
        if out is not None and ckpt_every and ((i + 1) % ckpt_every == 0):
            out.setdefault('partial', {})[key] = dict(
                L=L, beta=beta, N=N, done=i + 1, therm=therm, gap=gap, ns_ape=ns_ape,
                per_sweep=per_sweep, per_meas=meas_time / (i + 1),
                Obar=Obar[:i + 1].tolist(), Osq=Osq[:i + 1].tolist(),
                Obar_bare=Obar_bare[:i + 1].tolist(), chif=chif[:i + 1].tolist(),
                Pbar=Pbar[:i + 1].tolist())
            _dump(out)
    per_meas = meas_time / max(1, N)

    V = L ** 4
    chiO, chiO_se = G.jk_var_se(Obar); chiO *= V; chiO_se *= V
    C0 = float(Osq.mean() - Obar.mean() ** 2)
    chi_phys = chiO - C0
    chi_bare, _ = G.jk_var_se(Obar_bare); chi_bare *= V
    chifm = float(chif.mean()); chifse = float(chif.std(ddof=1) / np.sqrt(N))
    chiP, chiP_se = G.jk_var_se(Pbar); chiP *= L ** 3
    wall = t_therm + N * gap * per_sweep + meas_time

    return dict(L=L, beta=beta, N=N, therm=therm, gap=gap, ns_ape=ns_ape,
                per_sweep=per_sweep, per_meas=per_meas,
                wall_s=wall, wall_therm_s=t_therm, wall_meas_s=meas_time,
                chi_O=chiO, chi_O_se=chiO_se, Var_mean=chiO / V,
                contact_C0=C0, chi_phys=chi_phys, chi_bare=chi_bare,
                chi_flat=chifm, chi_flat_se=chifse, tau_flat=G.tau_int(chif),
                chi_P=chiP, chi_P_se=chiP_se,
                Obar_mean=float(Obar.mean()), Obar_sd=float(Obar.std(ddof=1)),
                chi_O_series=(V * (Obar - Obar.mean()) ** 2).tolist(),
                Obar_series=Obar.tolist())


def project(L_from, per_sweep, per_meas, L_to, N2, therm2, gap2):
    """Extrapolate the wall-clock of an L_to job using the measured L_from rates (cost ~ V = L^4)."""
    scale = (L_to / L_from) ** 4.0
    sweeps = therm2 + N2 * gap2
    return (sweeps * per_sweep + N2 * per_meas) * scale


def main():
    N12 = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    max16_min = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0
    print("# WG-L12-L16: gauge-invariant total-normalized APE-smeared tr F^2 zero mode (beta=2.5 SU(2))")
    print(f"# estimator = wg_invariant.py (identical); APE ns=2 alpha=0.5; one process; beta=2.5")
    out = dict(protocol=dict(beta=2.5, ns_ape=2, alpha=0.5, N12=N12,
                             max16_min=max16_min, anchors=ANCHORS),
               volumes={}, projection={})

    # ---------------- L=12 ----------------
    print(f"\n=== L=12  N={N12} therm=500 gap=8 ===", flush=True)
    r12 = run_volume(12, 2.5, N12, 500, 8, seed=1200, ns_ape=2,
                     ckpt_every=100, out=out, key='L12')
    out['volumes']['L12'] = r12
    out.setdefault('partial', {}).pop('L12', None)
    _dump(out)
    print(f"  L=12 done: per_sweep={r12['per_sweep']*1000:.1f} ms  per_meas={r12['per_meas']*1000:.1f} ms"
          f"  wall={r12['wall_s']:.0f}s")
    print(f"  chi_O   = {r12['chi_O']:.4f} +- {r12['chi_O_se']:.4f}   (total-norm, gauge-inv)")
    print(f"  C_s(0)  = {r12['contact_C0']:.4f}  ->  chi_phys = {r12['chi_phys']:.4f}")
    print(f"  chi_bare= {r12['chi_bare']:.4f}   chi_flat={r12['chi_flat']:.4f}"
          f"   chi_P={r12['chi_P']:.4f}", flush=True)

    # ---------------- project / run L=16 ----------------
    proj = {}
    for N2, therm2, gap2 in [(300, 400, 8), (300, 400, 6), (300, 400, 5), (300, 350, 4)]:
        proj[(N2, therm2, gap2)] = project(12, r12['per_sweep'], r12['per_meas'],
                                           16, N2, therm2, gap2)
    chosen = None
    for proto, sec in proj.items():
        if sec < max16_min * 60.0:
            chosen = (proto, sec); break
    out['projection'] = dict(
        per_sweep_L12=r12['per_sweep'], per_meas_L12=r12['per_meas'],
        candidates={f"N{p[0]}_therm{p[1]}_gap{p[2]}": s for p, s in proj.items()},
        max16_s=max16_min * 60.0,
        chosen=(f"N{chosen[0][0]}_therm{chosen[0][1]}_gap{chosen[0][2]}" if chosen else None),
        projected_s=(chosen[1] if chosen else None))
    _dump(out)
    print("\n=== L=16 projection (cost ~ L^4 from measured L=12 rates) ===")
    for p, s in proj.items():
        print(f"  N={p[0]} therm={p[1]} gap={p[2]}  ->  {s/60:.1f} min"
              + ("   <-- chosen" if chosen and (p, s) == chosen else ""))
    if chosen is None:
        print(f"  no candidate under {max16_min:.0f} min -> SKIP L=16 (recorded as projection)")
    else:
        (N2, therm2, gap2), sec = chosen
        print(f"  projection {sec/60:.1f} min < {max16_min:.0f} min -> RUN L=16 "
              f"N={N2} therm={therm2} gap={gap2}", flush=True)
        r16 = run_volume(16, 2.5, N2, therm2, gap2, seed=1600, ns_ape=2,
                         ckpt_every=50, out=out, key='L16')
        out['volumes']['L16'] = r16
        out.setdefault('partial', {}).pop('L16', None)
        _dump(out)
        print(f"  L=16 done: per_sweep={r16['per_sweep']*1000:.1f} ms  wall={r16['wall_s']:.0f}s "
              f"(actual {r16['wall_s']/60:.1f} min)")
        print(f"  chi_O   = {r16['chi_O']:.4f} +- {r16['chi_O_se']:.4f}   (total-norm, gauge-inv)")
        print(f"  C_s(0)  = {r16['contact_C0']:.4f}  ->  chi_phys = {r16['chi_phys']:.4f}")
        print(f"  chi_bare= {r16['chi_bare']:.4f}   chi_flat={r16['chi_flat']:.4f}"
              f"   chi_P={r16['chi_P']:.4f}", flush=True)

    # ---------------- trend / ratios vs free (L2/L1)^2 ----------------
    seq = [('L4', ANCHORS['L4']), ('L6', ANCHORS['L6']), ('L8', ANCHORS['L8'])]
    for k in ('L12', 'L16'):
        if k in out['volumes']:
            r = out['volumes'][k]
            seq.append((k, dict(chi_O=r['chi_O'], chi_O_se=r['chi_O_se'],
                                chi_phys=r['chi_phys'], N=r['N'])))
    trend = []
    for (ka, a), (kb, b) in zip(seq, seq[1:]):
        La = int(ka[1:]); Lb = int(kb[1:])
        ratio = b['chi_O'] / a['chi_O']
        rel = np.hypot(a['chi_O_se'] / a['chi_O'], b['chi_O_se'] / b['chi_O'])
        rp = b['chi_phys'] / a['chi_phys']
        trend.append(dict(pair=f"{ka}->{kb}", free=(Lb / La) ** 2,
                          chi_O_ratio=ratio, chi_O_ratio_se=ratio * rel,
                          chi_phys_ratio=rp,
                          chi_O_a=a['chi_O'], chi_O_b=b['chi_O']))
    # log-log slope of chi_O vs L (free massless = 2, plateau = 0); unweighted and inverse-variance weighted
    Ls = np.array([int(k[1:]) for k, _ in seq], float)
    ys = np.array([v['chi_O'] for _, v in seq], float)
    ys_se = np.array([v['chi_O_se'] for _, v in seq], float)
    slope = float(np.polyfit(np.log(Ls), np.log(ys), 1)[0]) if len(seq) >= 2 else None
    wslope = wse = None
    if len(seq) >= 2:
        x = np.log(Ls); y = np.log(ys); w = (ys / ys_se) ** 2
        A = np.vstack([x, np.ones_like(x)]).T
        Wm = np.diag(w)
        cov = np.linalg.inv(A.T @ Wm @ A)
        beta = cov @ (A.T @ Wm @ y)
        wslope = float(beta[0]); wse = float(np.sqrt(cov[0, 0]))
    out['trend'] = dict(sequence=[(k, v['chi_O'], v['chi_O_se']) for k, v in seq],
                        steps=trend, loglog_slope=slope,
                        loglog_slope_weighted=wslope, loglog_slope_se=wse,
                        free_slope=2.0, plateau_slope=0.0)
    _dump(out)

    print("\n=== L-trend: chi_O and per-step ratios vs free (L2/L1)^2 ===")
    for k, v, se in out['trend']['sequence']:
        print(f"  {k:>4}: chi_O = {v:.4f} +- {se:.4f}  (N={dict(seq)[k]['N']})")
    for t in trend:
        print(f"  {t['pair']}: ratio = {t['chi_O_ratio']:.3f} +- {t['chi_O_ratio_se']:.3f}"
              f"   vs free {t['free']:.2f}   (plateau = 1)")
    print(f"  log-log slope d ln chi_O / d ln L = {slope:+.3f} unweighted; "
          f"{wslope:+.3f} +- {wse:.3f} weighted   (free massless 2.0; plateau 0.0)")
    print(f"\n# wrote {RESULT}")
    print("Labels: [V] the SU(2) numbers; [O] continuum/uniform-in-L H1.  No QED.  Clay 0%.")


if __name__ == '__main__':
    main()
