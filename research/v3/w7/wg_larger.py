#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
wg_larger.py -- extend the WG gauge-invariant estimator (WG-INVARIANT-ESTIMATOR.md) to L=24
(and L=32 only if the measured L=24 sweep+measure rates project < 60 min).

Object (identical estimator, same measure, same protocol as wg_invariant.py / wg_large.py):
    O_s(x) = sum_{mu<nu} ( 1 - (1/2) Re tr P_{mu nu}^s(x) ),   P^s = plaquette of APE-smeared links
    q      = V^{-1/2} sum_x ( O_s(x) - <O_s> ),   V = L^4
    chi_O  = Var(q) = V * Var_c( Obar ),          Obar = (1/V) sum_x O_s(x)
Total-normalized (non-self-averaging), exactly gauge-invariant, APE ns=2 alpha=0.5 (fixed lattice
radius, SAME as L=4,6,8,12,16 so all points are directly comparable).  Reports chi_O, contact
C_s(0), chi_phys, chi_bare (ns=0 reference), chi_flat and chi_P references, jackknife SEs, tau_int.

Protocol: L=24 N=100 therm=200 gap=3 (suggested in the task); L=32 N=50 therm=150 gap=2 is run only
if the L=24 measured rates project the full job to < 3600 s (60 min).  ONE process (charter s4:
never run two heavy SU(N) L>=6 jobs at once).  Raw per-config series checkpointed to JSON.

Run:  py -3.12 -I research/v3/w7/wg_larger.py            # L=24 N=100 -> project/run L=32 N=50
      py -3.12 -I research/v3/w7/wg_larger.py 60 0       # N24=60, max32_min=0 (skip L=32)
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

RESULT = os.path.join(HERE, 'wg_larger_results.json')

# Existing anchors L=4..16 -- for the cross-volume trend.  L12,L16 from wg_large_results.json.
ANCHORS = {
    'L4':  dict(chi_O=0.32765014859553365,  chi_O_se=0.0222733807527846,  chi_phys=0.2974008604510313,  N=800),
    'L6':  dict(chi_O=0.37743611373983316,  chi_O_se=0.02133943765327901, chi_phys=0.3431821313604035,  N=400),
    'L8':  dict(chi_O=0.4104591281817326,   chi_O_se=0.057637361481588344, chi_phys=0.3754055557778887, N=160),
    'L12': dict(chi_O=0.468797, chi_O_se=0.027880, chi_phys=0.431825, N=400),
    'L16': dict(chi_O=0.586532, chi_O_se=0.054732, chi_phys=0.549158, N=300),
}


def _dump(out):
    with open(RESULT, 'w') as fh:
        json.dump(out, fh, indent=1,
                  default=lambda o: float(o) if isinstance(o, (np.floating, np.ndarray)) else str(o))


def run_volume(L, beta, N, therm, gap, seed, ns_ape=2, ckpt_every=10, out=None, key=None):
    """One volume, one process.  Times therm and measurement separately for the L=32 projection.
    Checkpoints the raw per-config series to `out` every ckpt_every configs (incremental JSON),
    together with a running chi_O estimate so an interrupted run still yields a reading."""
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
    V = L ** 4
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
            chi_now, se_now = G.jk_var_se(Obar[:i + 1]); chi_now *= V; se_now *= V
            out.setdefault('partial', {})[key] = dict(
                L=L, beta=beta, N=N, done=i + 1, therm=therm, gap=gap, ns_ape=ns_ape,
                per_sweep=per_sweep, per_meas=meas_time / (i + 1),
                chi_O_running=chi_now, chi_O_running_se=se_now,
                Obar=Obar[:i + 1].tolist(), Osq=Osq[:i + 1].tolist(),
                Obar_bare=Obar_bare[:i + 1].tolist(), chif=chif[:i + 1].tolist(),
                Pbar=Pbar[:i + 1].tolist())
            _dump(out)
            print(f"    [{key}] {i+1}/{N} configs  chi_O(run)={chi_now:.4f}+-{se_now:.4f}"
                  f"  per_meas={meas_time/(i+1):.2f}s", flush=True)
    per_meas = meas_time / max(1, N)

    chiO, chiO_se = G.jk_var_se(Obar); chiO *= V; chiO_se *= V
    C0 = float(Osq.mean() - Obar.mean() ** 2)
    chi_phys = chiO - C0
    chi_bare, _ = G.jk_var_se(Obar_bare); chi_bare *= V
    chifm = float(chif.mean()); chifse = float(chif.std(ddof=1) / np.sqrt(N))
    chiP, chiP_se = G.jk_var_se(Pbar); chiP *= L ** 3
    # integrated autocorrelation time of the chi_O per-config series (in config units)
    chi_series = V * (Obar - Obar.mean()) ** 2
    wall = t_therm + N * gap * per_sweep + meas_time

    return dict(L=L, beta=beta, N=N, therm=therm, gap=gap, ns_ape=ns_ape,
                per_sweep=per_sweep, per_meas=per_meas,
                wall_s=wall, wall_therm_s=t_therm, wall_meas_s=meas_time,
                chi_O=chiO, chi_O_se=chiO_se, Var_mean=chiO / V,
                contact_C0=C0, chi_phys=chi_phys, chi_bare=chi_bare,
                chi_flat=chifm, chi_flat_se=chifse, tau_flat=G.tau_int(chif),
                chi_P=chiP, chi_P_se=chiP_se,
                tau_chiO=G.tau_int(chi_series),
                Obar_mean=float(Obar.mean()), Obar_sd=float(Obar.std(ddof=1)),
                chi_O_series=chi_series.tolist(),
                Obar_series=Obar.tolist())


def project(L_from, per_sweep, per_meas, L_to, N2, therm2, gap2):
    """Extrapolate the wall-clock of an L_to job using the measured L_from rates (cost ~ V = L^4)."""
    scale = (L_to / L_from) ** 4.0
    sweeps = therm2 + N2 * gap2
    return (sweeps * per_sweep + N2 * per_meas) * scale


def saturating_fit(seq):
    """Fit chi_O(L) = chi_inf - c / L^p on the given (L, chi_O) points (unweighted, [C] only).
    Returns (chi_inf, c, p, residual rms) or None."""
    Ls = np.array([int(k[1:]) for k, _ in seq], float)
    ys = np.array([v['chi_O'] for _, v in seq], float)
    best = None
    for p in np.linspace(0.1, 8.0, 160):
        A = np.vstack([np.ones_like(Ls), -1.0 / Ls ** p]).T
        try:
            (ci, c), *_ = np.linalg.lstsq(A, ys, rcond=None)
        except Exception:
            continue
        r = A @ np.array([ci, c]) - ys
        s = float(np.sqrt(np.mean(r * r)))
        if best is None or s < best[3]:
            best = (float(ci), float(c), float(p), s)
    return best


def main():
    N24 = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    max32_min = float(sys.argv[2]) if len(sys.argv) > 2 else 60.0
    print("# WG-L24: gauge-invariant total-normalized APE-smeared tr F^2 zero mode (beta=2.5 SU(2))")
    print("# estimator = wg_invariant.py (identical); APE ns=2 alpha=0.5; one process; beta=2.5")
    out = dict(protocol=dict(beta=2.5, ns_ape=2, alpha=0.5, N24=N24,
                             max32_min=max32_min, anchors=ANCHORS),
               volumes={}, projection={})

    # ---------------- L=24 ----------------
    print(f"\n=== L=24  N={N24} therm=200 gap=3 ===", flush=True)
    r24 = run_volume(24, 2.5, N24, 200, 3, seed=2400, ns_ape=2,
                     ckpt_every=10, out=out, key='L24')
    out['volumes']['L24'] = r24
    out.setdefault('partial', {}).pop('L24', None)
    _dump(out)
    print(f"  L=24 done: per_sweep={r24['per_sweep']*1000:.1f} ms  per_meas={r24['per_meas']*1000:.1f} ms"
          f"  wall={r24['wall_s']:.0f}s ({r24['wall_s']/60:.1f} min)")
    print(f"  chi_O   = {r24['chi_O']:.4f} +- {r24['chi_O_se']:.4f}   (total-norm, gauge-inv)")
    print(f"  C_s(0)  = {r24['contact_C0']:.4f}  ->  chi_phys = {r24['chi_phys']:.4f}")
    print(f"  chi_bare= {r24['chi_bare']:.4f}   chi_flat={r24['chi_flat']:.4f}"
          f"   chi_P={r24['chi_P']:.4f}   tau_chiO={r24['tau_chiO']:.2f}", flush=True)

    # ---------------- project / run L=32 ----------------
    proj = {}
    for N2, therm2, gap2 in [(50, 150, 2), (60, 150, 2), (50, 150, 3)]:
        proj[(N2, therm2, gap2)] = project(24, r24['per_sweep'], r24['per_meas'],
                                           32, N2, therm2, gap2)
    chosen = None
    for proto, sec in proj.items():
        if sec < max32_min * 60.0:
            chosen = (proto, sec); break
    out['projection'] = dict(
        per_sweep_L24=r24['per_sweep'], per_meas_L24=r24['per_meas'],
        candidates={f"N{p[0]}_therm{p[1]}_gap{p[2]}": s for p, s in proj.items()},
        max32_s=max32_min * 60.0,
        chosen=(f"N{chosen[0][0]}_therm{chosen[0][1]}_gap{chosen[0][2]}" if chosen else None),
        projected_s=(chosen[1] if chosen else None))
    _dump(out)
    print("\n=== L=32 projection (cost ~ L^4 from measured L=24 rates) ===")
    for p, s in proj.items():
        print(f"  N={p[0]} therm={p[1]} gap={p[2]}  ->  {s/60:.1f} min"
              + ("   <-- chosen" if chosen and (p, s) == chosen else ""))
    if chosen is None:
        print(f"  no candidate under {max32_min:.0f} min -> SKIP L=32 (recorded as projection)")
    else:
        (N2, therm2, gap2), sec = chosen
        print(f"  projection {sec/60:.1f} min < {max32_min:.0f} min -> RUN L=32 "
              f"N={N2} therm={therm2} gap={gap2}", flush=True)
        r32 = run_volume(32, 2.5, N2, therm2, gap2, seed=3200, ns_ape=2,
                         ckpt_every=5, out=out, key='L32')
        out['volumes']['L32'] = r32
        out.setdefault('partial', {}).pop('L32', None)
        _dump(out)
        print(f"  L=32 done: per_sweep={r32['per_sweep']*1000:.1f} ms  wall={r32['wall_s']:.0f}s "
              f"(actual {r32['wall_s']/60:.1f} min)")
        print(f"  chi_O   = {r32['chi_O']:.4f} +- {r32['chi_O_se']:.4f}   (total-norm, gauge-inv)")
        print(f"  C_s(0)  = {r32['contact_C0']:.4f}  ->  chi_phys = {r32['chi_phys']:.4f}")
        print(f"  chi_bare= {r32['chi_bare']:.4f}   chi_flat={r32['chi_flat']:.4f}"
              f"   chi_P={r32['chi_P']:.4f}   tau_chiO={r32['tau_chiO']:.2f}", flush=True)

    # ---------------- trend / ratios vs free (L2/L1)^2 ----------------
    seq = [('L4', ANCHORS['L4']), ('L6', ANCHORS['L6']), ('L8', ANCHORS['L8']),
           ('L12', ANCHORS['L12']), ('L16', ANCHORS['L16'])]
    for k in ('L24', 'L32'):
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

    def loglog(pts):
        Ls = np.array([int(k[1:]) for k, _ in pts], float)
        ys = np.array([v['chi_O'] for _, v in pts], float)
        ys_se = np.array([v['chi_O_se'] for _, v in pts], float)
        slope = float(np.polyfit(np.log(Ls), np.log(ys), 1)[0]) if len(pts) >= 2 else None
        wslope = wse = None
        if len(pts) >= 2:
            x = np.log(Ls); y = np.log(ys); w = (ys / ys_se) ** 2
            A = np.vstack([x, np.ones_like(x)]).T
            Wm = np.diag(w)
            cov = np.linalg.inv(A.T @ Wm @ A)
            beta = cov @ (A.T @ Wm @ y)
            wslope = float(beta[0]); wse = float(np.sqrt(cov[0, 0]))
        return slope, wslope, wse

    slope_all, wslope_all, wse_all = loglog(seq)
    seq_big = [s for s in seq if int(s[0][1:]) >= 12]
    _, wslope_big, wse_big = loglog(seq_big)
    sat = saturating_fit(seq)
    out['trend'] = dict(sequence=[(k, v['chi_O'], v['chi_O_se']) for k, v in seq],
                        steps=trend,
                        loglog_slope=slope_all,
                        loglog_slope_weighted=wslope_all, loglog_slope_se=wse_all,
                        loglog_slope_Lge12=wslope_big, loglog_slope_Lge12_se=wse_big,
                        saturating_fit=dict(chi_inf=sat[0], c=sat[1], p=sat[2], rms=sat[3])
                        if sat else None,
                        free_slope=2.0, plateau_slope=0.0)
    _dump(out)

    print("\n=== L-trend: chi_O and per-step ratios vs free (L2/L1)^2 ===")
    for k, v, se in out['trend']['sequence']:
        print(f"  {k:>4}: chi_O = {v:.4f} +- {se:.4f}  (N={dict(seq)[k]['N']})")
    for t in trend:
        print(f"  {t['pair']}: ratio = {t['chi_O_ratio']:.3f} +- {t['chi_O_ratio_se']:.3f}"
              f"   vs free {t['free']:.2f}   (plateau = 1)")
    print(f"  log-log slope d ln chi_O / d ln L = {slope_all:+.3f} unweighted; "
          f"{wslope_all:+.3f} +- {wse_all:.3f} weighted   (all L; free 2.0, plateau 0.0)")
    print(f"  log-log slope L>=12: {wslope_big:+.3f} +- {wse_big:.3f} weighted")
    if sat:
        print(f"  saturating ansatz chi_inf - c/L^p: chi_inf={sat[0]:.4f}, c={sat[1]:.4f}, "
              f"p={sat[2]:.2f}, rms={sat[3]:.4f}   ([C], not a theorem)")
    print(f"\n# wrote {RESULT}")
    print("Labels: [V] the SU(2) numbers; [O] continuum/uniform-in-L H1.  No QED.  Clay 0%.")


if __name__ == '__main__':
    main()
