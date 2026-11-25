#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
m1_4d_analyze.py -- read m1_4d_run_results.json and produce clean summary numbers:
  chi_flat (and mean-subtracted Var), plateau ratios with propagated SE,
  Polyakov correlation length via a periodic cosh fit + nearest-neighbour contrast,
  0^{++} correlator effective masses (resolved or not), Landau chi + Gribov spread.

py -3.12 -I ; no third-party sources.
"""
import json, os
import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
PATH = os.path.join(HERE, 'm1_4d_run_results.json')


def cosh_xi(Cp, L):
    """Fit Cp(r)=A[exp(-r/xi)+exp(-(L-r)/xi)] over r=0..L-1 (periodic), return (xi, xi_se_boot)."""
    r = np.arange(L)
    y = np.asarray(Cp, float)
    def model(p):
        A, xi = np.exp(p[0]), np.exp(p[1])
        return A * (np.exp(-r / xi) + np.exp(-(L - r) / xi))
    def resid(p):
        return model(p) - y
    try:
        sol = least_squares(resid, x0=[np.log(max(y[0], 1e-6)), np.log(1.0)],
                            bounds=([-60, np.log(0.05)], [60, np.log(1e3)]), max_nfev=4000)
        return float(np.exp(sol.x[1])), float(np.sqrt(np.mean(sol.fun ** 2)))
    except Exception:
        return float('nan'), float('nan')


def nn_xi(Cp):
    """Nearest-neighbour correlation length from r=0->1 contrast: xi = -1/ln(Cp1/Cp0)."""
    if len(Cp) < 2 or Cp[1] <= 0 or Cp[0] <= 0:
        return float('nan')
    r = Cp[1] / Cp[0]
    if not (0 < r < 1):
        return float('nan')
    return float(-1.0 / np.log(r))


def meff(Cc, L):
    out = []
    for t in range(1, L // 2):
        if Cc[t] > 0 and Cc[t + 1] > 0:
            out.append(round(float(np.log(Cc[t] / Cc[t + 1])), 4))
        else:
            out.append(float('nan'))
    return out


def main():
    d = json.load(open(PATH))
    print("# M1-4D analysis\n")
    print(f"protocol: {d.get('meta',{}).get('protocol')}\n")
    have = [L for L in (4, 6, 8, 12, 16) if f'L{L}' in d]
    print(f"{'L':>3} {'chi_flat':>16} {'chi_ms':>14} {'plaq':>7} {'C_P(0)':>8} {'xi_cosh':>8} {'xi_nn':>6} "
          f"{'chi_Land':>15} {'Grib_%':>7}")
    for L in have:
        r = d[f'L{L}']
        cp = r['C_P']
        xi_c, rms = cosh_xi(cp, L)
        xi_n = nn_xi(cp)
        ld = r.get('landau', {})
        gr = r.get('gribov', {})
        print(f"{L:>3} {r['chi_flat']:>10.4f}+-{r['chi_flat_se']:.4f} "
              f"{r.get('chi_flat_ms',float('nan')):>7.3f}+-{r.get('chi_flat_ms_se',float('nan')):.3f} "
              f"{r['plaq']:>7.4f} {cp[0]:>8.4f} {xi_c:>8.3f} {xi_n:>6.3f} "
              f"{ld.get('chi_ms',float('nan')):>9.3f}+-{ld.get('chi_ms_se',float('nan')):.2f} "
              f"{100*gr.get('mean_spread_rel',float('nan')):>6.1f}")
    print("\n## C_P(r) tables (gauge-invariant Polyakov spatial correlator)")
    for L in have:
        print(f"  L={L}: {np.round(d[f'L{L}']['C_P'],4).tolist()}")
    print("\n## 0++ plaquette correlator effective mass m_eff(t)=ln C(t)/C(t+1)  (flat => unresolved)")
    for L in have:
        g = d[f'L{L}']['glueball']
        print(f"  L={L} unsmeared: {meff(g['O_unsmeared'],L)}")
        print(f"  L={L} ape2     : {meff(g['O_ape2'],L)}")
    print("\n## plateau ratios chi_flat (propagated SE)")
    for a, b in zip(have, have[1:]):
        ra, rb = d[f'L{a}'], d[f'L{b}']
        ratio = rb['chi_flat'] / ra['chi_flat']
        rel = np.hypot(ra['chi_flat_se'] / ra['chi_flat'], rb['chi_flat_se'] / rb['chi_flat'])
        free = (b / a) ** 2
        print(f"  L={a}->{b}: {ra['chi_flat']:.3f}->{rb['chi_flat']:.3f}  ratio={ratio:.3f}+-{ratio*rel:.3f}"
              f"   (free-massless L^2 = {free:.2f}; free L^4(norm.) = {(b/a)**4:.1f})")
    print("\n## Landau chi (mean-subtracted Var) trend")
    for L in have:
        ld = d[f'L{L}'].get('landau', {})
        if ld:
            print(f"  L={L}: chi_Landau(Var)={ld.get('chi_ms',float('nan')):.2f}+-{ld.get('chi_ms_se',float('nan')):.2f}"
                  f"  raw={ld.get('chi_mean',float('nan')):.2f}  mean_contrib={ld.get('mean_contrib',float('nan')):.3f}"
                  f"  Gribov rel mean/max={100*d[f'L{L}']['gribov']['mean_spread_rel']:.1f}%/"
                  f"{100*d[f'L{L}']['gribov']['max_spread_rel']:.1f}%")


if __name__ == '__main__':
    main()
