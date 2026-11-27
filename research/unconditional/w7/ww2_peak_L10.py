"""
WW2: three-point scan of <plaq>(beta) at L=10, beta in {2.40, 2.45, 2.50}.

Purpose: ww1-peak established chi_q non-increasing on [2.4,2.5] at L=8
(12.8 sigma) and beta_peak(8) <= 2.40.  The only residual was the L-uniformity:
    beta_peak(L) < 2.5   for all L
This run adds L=10 to the L=4/6/8 chain, i.e. it tests whether beta_peak drifts
upward with L.

Also reports, per beta:
    plaq, chi_q = -(1/6) d<qbar>/dbeta, and the (EP) test
    plaq'(2.5) <= plaq(2.5)/2.5

Reuses the tested SU(2) heat bath via im1_varpi.run.
One beta at a time (memory discipline).
"""
import json
import sys
import time

import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
from im1_varpi import run  # noqa: E402

BETAS = (2.40, 2.45, 2.50)
L = 10
NCFG = 4000
THERM = 300
GAP = 2
SEED = 20261129


def main():
    print("#" * 64)
    print(f"# WW2 three-point scan   L={L}  N={NCFG}  therm={THERM}  gap={GAP}")
    print("#" * 64)
    out = {}
    for b in BETAS:
        t0 = time.time()
        plaq, res, wall = run(L, b, THERM, NCFG, GAP, SEED, [1],
                              log=lambda *a: None)
        v = float(plaq.mean())
        e = float(plaq.std(ddof=1) / np.sqrt(len(plaq)))
        out[b] = [v, e]
        print(f"L{L} BETA {b:.3f}  plaq={v:.6f} +- {e:.6f}"
              f"  wall={wall:.0f}s  total={time.time()-t0:.0f}s")
        sys.stdout.flush()

    p = {b: out[b][0] for b in BETAS}
    e = {b: out[b][1] for b in BETAS}
    h = 0.05
    print()
    print("# ---- summary L=%d ----" % L)
    for b in BETAS:
        print(f"SUM BETA {b:.3f}  plaq={p[b]:.6f} +- {e[b]:.6f}"
              f"  plaq/beta={p[b]/b:.6f}  qbar={1-p[b]:.6f}")

    # chi_q = -(1/6) d<qbar>/dbeta = +(1/6) dplaq/dbeta
    chi = {}
    for b in (2.45,):
        dp = (p[2.50] - p[2.40]) / (2 * h)
        edp = float(np.hypot(e[2.50], e[2.40]) / (2 * h))
        chi[b] = [dp / 6.0, edp / 6.0]
        print(f"CHI_Q chi_q(2.45, central) = {dp/6.0:.6f} +- {edp/6.0:.6f}")

    # one-sided slope at 2.5
    dp25 = (p[2.50] - p[2.45]) / h
    edp25 = float(np.hypot(e[2.50], e[2.45]) / h)
    print(f"DERIV plaq'(2.5) = {dp25:.6f} +- {edp25:.6f}")
    print(f"CHI_Q chi_q(2.5) = {dp25/6.0:.6f} +- {edp25/6.0:.6f}")

    thr = p[2.50] / 2.5
    print(f"THRESH plaq(2.5)/2.5 = {thr:.6f}")
    print(f"EP_TEST thr - plaq'(2.5) = {thr - dp25:.6f} +- {edp25:.6f}"
          f"   ({'PASS' if dp25 <= thr else 'FAIL'})")

    # concavity (centred at 2.45)
    d2 = (p[2.50] - 2 * p[2.45] + p[2.40]) / h ** 2
    print(f"CONCAV plaq''(2.45) = {d2:.6f}   ({'CONCAVE' if d2 <= 0 else 'CONVEX'})")

    # monotonicity of chi_q over the window (midpoint finite differences)
    m1 = (p[2.45] - p[2.40]) / h / 6.0
    m2 = (p[2.50] - p[2.45]) / h / 6.0
    print(f"MONO chi_q midpoints: {m1:.6f} -> {m2:.6f}  diff={m1-m2:+.6f}"
          f"   ({'NONINC' if m2 <= m1 else 'INC'})")

    # peak from the grid
    am = max(BETAS, key=lambda b: p[b] / b)
    print(f"PEAK beta_peak grid argmax(plaq/beta) = {am:.2f}  (L={L})")

    json.dump({'L': L, 'N': NCFG, 'betas': list(BETAS), 'plaq': out,
               'chi_q_245': chi[2.45], 'plaq_p_25': dp25, 'plaq_p_25_err': edp25,
               'ep_pass': bool(dp25 <= thr), 'concav': d2,
               'peak_grid': am},
              open(f'ww2_peak_L{L}.json', 'w'), indent=2)
    print(f"# written ww2_peak_L{L}.json")


if __name__ == '__main__':
    main()
