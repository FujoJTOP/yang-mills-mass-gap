"""
WW1: three-point scan of <plaq>(beta) at L=6, beta in {2.4, 2.45, 2.5},
to decide the sign of  d/dbeta ( plaq/beta )  at beta=2.5.

Target statement (A cell, minimal form):
    plaq'(2.5) <= plaq(2.5)/2.5     <=>   d/dbeta(plaq/beta)|_{2.5} <= 0

Also gives chi_q = -(1/6) d<qbar>/dbeta  (qbar = 1 - plaq)  and its slope,
i.e. the sign of <Q^3>_c  since   d_beta chi_q = -6 L^{-4} <Q^3>_c .

Reuses the tested SU(2) heat bath via im1_varpi.run.
One beta at a time (memory discipline).
"""
import sys
import time

import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
from im1_varpi import run  # noqa: E402

BETAS = (2.40, 2.45, 2.50)
L = 6
NCFG = 5000
THERM = 300
GAP = 2
SEEDS = (20261127, 20261128)


def jk_err(x):
    """jackknife error of the mean"""
    n = len(x)
    tot = x.sum()
    th = (tot - x) / (n - 1)
    m = th.mean()
    return np.sqrt((n - 1) / n * np.sum((th - m) ** 2))


def main():
    print("#" * 60)
    print(f"# WW1 three-point scan   L={L}  N={NCFG}  therm={THERM}  gap={GAP}")
    print("#" * 60)
    out = {}
    for b in BETAS:
        acc = []
        errs = []
        for s in SEEDS:
            t0 = time.time()
            plaq, res, wall = run(L, b, THERM, NCFG, GAP, s, [1],
                                  log=lambda *a: None)
            acc.append(plaq.mean())
            errs.append(plaq.std(ddof=1) / np.sqrt(len(plaq)))
            print(f"BETA {b:.3f} seed={s}  <plaq>={plaq.mean():.6f} "
                  f"+- {errs[-1]:.6f}  wall={wall:.0f}s")
        out[b] = (float(np.mean(acc)),
                  float(np.sqrt(np.mean(np.square(errs)) / len(errs))))
    print()
    print("# ---- summary ----")
    for b in BETAS:
        v, e = out[b]
        print(f"SUM BETA {b:.3f}  plaq={v:.6f} +- {e:.6f}  plaq/beta={v/b:.6f}")

    # finite differences
    b1, b2, b3 = BETAS
    p1, p2, p3 = (out[b][0] for b in BETAS)
    e1, e2, e3 = (out[b][1] for b in BETAS)
    h = b2 - b1
    # first derivative at 2.5: one-sided (p3-p2)/h
    dp = (p3 - p2) / h
    edp = float(np.hypot(e3, e2) / h)
    print(f"DERIV plaq'(2.5) = {dp:.6f} +- {edp:.6f}")
    # chi_q = -(1/6) d<qbar>/dbeta = +(1/6) dplaq/dbeta
    chi = dp / 6.0
    echi = edp / 6.0
    print(f"CHI_Q chi_q(2.5) = {chi:.6f} +- {echi:.6f}")
    # statement: plaq'(2.5) <= plaq(2.5)/2.5 ?
    thr = p3 / 2.5
    print(f"THRESH plaq(2.5)/2.5 = {thr:.6f}")
    print(f"GAP    thr - plaq' = {thr - dp:.6f} +- {edp:.6f}"
          f"   ({'PASS' if dp <= thr else 'FAIL'})")
    # d/dbeta(plaq/beta) = plaq'/beta - plaq/beta^2
    dpp = dp / 2.5 - p3 / 2.5 ** 2
    print(f"DERIV2 d/dbeta(plaq/beta)|2.5 = {dpp:.6f}"
          f"   ({'NONINC' if dpp <= 0 else 'INC'})")
    # rho = plaq / (15 C_q(0));  C_q(0) = Var(q) ~ 0.0720 at beta=2.5 (IM1)
    #   rho* comparison uses measured C_q(0); we only report plaq-side numbers.
    print()
    print("# NOTE: rho = plaq/(15 C_q(0)); use C_q(0)=0.0720 (IM1, L=8)")
    cq0 = 0.0720
    print(f"RHO_INFER rho(2.5) ~ {p3 / (15 * cq0):.4f}   (rho* = 3.066877)")
    # second difference for concavity at 2.5 via three points
    d2 = (p3 - 2 * p2 + p1) / h ** 2
    print(f"CONCAV plaq''(approx, centred at 2.45) = {d2:.6f}"
          f"   ({'CONCAVE' if d2 <= 0 else 'CONVEX'})")


if __name__ == '__main__':
    main()
