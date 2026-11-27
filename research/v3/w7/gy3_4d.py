#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gy3_4d.py -- GY3 tasks (c)+(d): 4D SU(2) Wilson, chi_L L-scaling at beta=2.5 and strong
coupling.  Reuses the M1 witness code (m1_4d_run.py, same normalization/protocol:
chi_flat = (1/3)|A|^2, A = N^{-1/2} sum_l log U_l, k=0 flat-chart mode; heatbath,
therm/gap; space-summed-plaquette 0++ correlator; Polyakov correlator).  No Landau.

usage:  py -3.12 -I gy3_4d.py <beta> <L> <N> <therm> <gap> [--noape] [--seed S]
Budget: one (beta,L) pair per invocation, <= 5 min real run.
Output: gy3_4d_results.json (key b<beta>_L<L>), log lines to stdout.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import m1_4d_run as M          # reuses make_stream, k0_vector, estimators (M1 protocol)
import w4_blocked as W

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, 'gy3_4d_results.json')


def run(L, beta, N, therm, gap, seed, ape=True):
    U, sweep = M.make_stream(L, beta, seed, cold=True)
    t0 = time.time()
    for _ in range(therm):
        sweep()
    therm_s = time.time() - t0

    Aflat = np.empty((N, 3))
    chF = np.empty(N); plq = np.empty(N); Pavg = np.empty(N)
    Cp = np.zeros(L)
    C_raw0 = np.zeros(L); sum_o0 = np.zeros(L)
    C_raw2 = np.zeros(L); sum_o2 = np.zeros(L)
    n_corr = 0
    for i in range(N):
        for _ in range(gap):
            sweep()
        Av = M.k0_vector(U)
        Aflat[i] = Av
        chF[i] = float(Av @ Av / 3.0)
        plq[i] = M.plaq(U, beta)
        P = M.polyakov_field(U)
        Pavg[i] = P.mean()
        Cp += np.array([float(np.mean(P * np.roll(P, -r, axis=0))) for r in range(L)])
        o0 = M.plaq_operator_t(U)
        for t0i in range(L):
            oc = np.roll(o0, -t0i)
            C_raw0 += o0 * oc
        sum_o0 += o0
        if ape:
            U2 = W.ape_all2(U, alpha=0.5, nsteps=2)
            o2 = M.plaq_operator_t(U2)
            for t0i in range(L):
                oc = np.roll(o2, -t0i)
                C_raw2 += o2 * oc
            sum_o2 += o2
        n_corr += 1

    Cp /= N
    Pm = float(Pavg.mean())
    C_P = Cp - Pm ** 2
    Om0 = float(sum_o0.mean()) / n_corr
    C0 = C_raw0 / n_corr - Om0 ** 2
    fit0 = M.fit_glueball(C0, L)
    out = dict(L=L, beta=beta, N=N, therm=therm, gap=gap, seed=seed, ape=ape,
               wall_therm_s=therm_s, wall_total_s=time.time() - t0,
               chi_flat=float(chF.mean()), chi_flat_sd=float(chF.std()),
               plaq=float(plq.mean()), Pavg=Pm, chi_P=float(L ** 3 * Pavg.var()),
               C_P=[float(c) for c in C_P], xi_P=M.fit_corr_length(C_P, L))
    ms, se, mc = M.jk_ms_var(Aflat)
    mF, sF = M.jk_mean_se(chF)
    tF = M.tau_int(chF)
    out.update(chi_flat_jk=mF, chi_flat_se=sF, chi_flat_tau=tF, chi_flat_neff=N / tF,
               chi_ms=ms, chi_ms_se=se, chi_mean_contrib=mc)
    Om2 = float(sum_o2.mean()) / n_corr if ape else 0.0
    C2 = C_raw2 / n_corr - Om2 ** 2 if ape else [0.0] * L
    out['glueball'] = dict(O_unsmeared=[float(c) for c in C0], Om0=Om0,
                           fit_unsmeared=fit0, meff_unsmeared=M.eff_mass(C0, L),
                           meff_ape2=M.eff_mass(list(C2), L) if ape else None,
                           fit_ape2=M.fit_glueball(list(C2), L) if ape else None)
    return out


def main():
    args = sys.argv[1:]
    beta = float(args[0]); L = int(args[1])
    N = int(args[2]); therm = int(args[3]); gap = int(args[4])
    ape = '--noape' not in args
    seed = 7000 + int(beta * 100) * 10 + L
    if '--seed' in args:
        seed = int(args[args.index('--seed') + 1])
    print(f"# gy3_4d: beta={beta} L={L} N={N} therm={therm} gap={gap} ape={ape} seed={seed}",
          flush=True)
    r = run(L, beta, N, therm, gap, seed, ape)
    out = {}
    if os.path.exists(RES):
        try:
            out = json.load(open(RES))
        except Exception:
            out = {}
    out[f'b{beta}_L{L}'] = r
    with open(RES, 'w') as fh:
        json.dump(out, fh, indent=1, default=lambda o: int(o) if isinstance(o, np.integer)
                  else (float(o) if isinstance(o, (np.floating, np.ndarray)) else str(o)))
    print(f"  chi_flat = {r['chi_flat']:.4f} +- {r['chi_flat_se']:.4f}  "
          f"(tau={r['chi_flat_tau']:.2f} Neff={r['chi_flat_neff']:.0f} sd={r['chi_flat_sd']:.2f} "
          f"ms={r['chi_ms']:.4f}+-{r['chi_ms_se']:.4f} meanc={r['chi_mean_contrib']:.3f})", flush=True)
    print(f"  plaq={r['plaq']:.4f}  chi_P={r['chi_P']:.3f}  xi_P={r['xi_P']:.3f}", flush=True)
    gb = r['glueball']
    f = gb['fit_unsmeared']
    print(f"  0++ unsmeared: m1={f['m1']:.4f} A={f['A']:.4g} rms={f.get('logfit_rms', float('nan')):.3f}"
          f"  m_eff={np.round(gb['meff_unsmeared'], 3)}", flush=True)
    if gb['fit_ape2']:
        f2 = gb['fit_ape2']
        print(f"  0++ ape2     : m1={f2['m1']:.4f} A={f2['A']:.4g} rms={f2.get('logfit_rms', float('nan')):.3f}"
              f"  m_eff={np.round(gb['meff_ape2'], 3)}", flush=True)
    print(f"  wall: therm {r['wall_therm_s']:.0f}s total {r['wall_total_s']:.0f}s -> {RES}",
          flush=True)


if __name__ == '__main__':
    main()
