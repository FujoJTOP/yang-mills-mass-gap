"""B01 — w1(raw) at beta = 3, 4, 6, 10.  Same operator / same lattice as run1+run3.

Definition (unchanged from run3, so numbers are comparable to w1 = 0.3427(119) at beta=2.5):
    S(t) = (1/3) sum_{space} sum_{3 spatial planes} (1/2) Re tr U_p   (zero momentum, slice t)
    C(t) = (1/Lt) sum_{t0} <S(t0+t) S(t0)>_c
    w1   = A / C(0),  A = sum_{t=1..3} g_t C(t) / sum g_t^2,  g_t = e^{-m t} + e^{-m(Lt-t)}
    m = 1.546 fixed  -> headline w1 (exactly run1/run3's estimator)
    w1_fit           -> free-mass variant (m by least squares on t=1..6), 0.4058 at beta=2.5
Errors: delete-one 20-block jackknife over configurations.

Physics core (heat bath, link update, plaquette, slice operator) is IMPORTED from
run3_apebase.py, with nlev=1 so the T-APE levels are simply never built (raw only).
Sequential by design: one beta at a time, single process, numba threads = 8.

Run:  py -3.12 b01_scan.py --betas 3 4 6 10 --n 24000 --tag scan
"""
import os
import sys
import time
import json
import argparse

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from run3_apebase import _chains, index_tables, hot, correlators, fit_A, A_at  # noqa: E402

M1546 = 1.546
OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def tau_int(x):
    """Integrated autocorrelation time, in configurations (Sokal window)."""
    x = np.asarray(x, float)
    n = len(x)
    y = x - x.mean()
    v = float(np.dot(y, y)) / n
    if v <= 0:
        return 0.5
    t = 0.5
    for lag in range(1, min(n // 4, 300)):
        c = float(np.dot(y[:-lag], y[lag:])) / (n - lag) / v
        if c <= 0:
            break
        t += c
    return t


def jack(fn, S, nb=20):
    """Delete-one block jackknife over the config index; fn takes a correlation vector.
    NOTE: correlators() sums over axis=1, so it needs the 2-D (N, Lt) slice, not (N,1,Lt)."""
    idx = np.array_split(np.arange(S.shape[0]), nb)
    vals = np.empty(nb)
    for b in range(nb):
        keep = np.concatenate([idx[j] for j in range(nb) if j != b])
        Sm = S[keep]
        vals[b] = fn(correlators(Sm[:, 0, :], Sm.shape[2]))
    m = np.nanmean(vals)
    e = np.sqrt((nb - 1) / nb * np.nansum((vals - m) ** 2))
    return float(m), float(e)


def point(beta, n, gap, therm, chains, threads, seed):
    import numba
    numba.set_num_threads(threads)
    Ls, Lt = 8, 12
    per = max(1, n // chains)
    N = per * chains
    plus, minus = index_tables(Ls, Lt)
    U = hot(chains, Ls, Lt, seed)
    rng = np.empty((chains, 2), np.uint64)
    ss = np.uint64(seed)
    for c in range(chains):
        ss = ss * np.uint64(6364136223846793005) + np.uint64(1442695040888963407)
        rng[c, 0] = ss | np.uint64(1)
        s1 = ss ^ np.uint64(0x9E3779B97F4A7C15)
        rng[c, 1] = s1 if s1 != np.uint64(0) else np.uint64(0x2545F4914F6CDD1D)
    S = np.zeros((N, 1, Lt))
    plaq = np.zeros(N)

    t0 = time.time()
    _chains(U[:1], rng[:1], plus, minus, beta / 2.0, 2, gap, 2, 0.5, 1, Lt,
            np.zeros((2, 1, Lt)), np.zeros(2))            # compile / warm-up (1 chain)
    t_comp = time.time() - t0

    t1 = time.time()
    _chains(U, rng, plus, minus, beta / 2.0, per, gap, therm, 0.5, 1, Lt, S, plaq)
    dtw = time.time() - t1

    Sb = S[:, 0, :].reshape(chains, per, Lt)
    pb = plaq.reshape(chains, per)
    tau_p = float(np.mean([tau_int(pb[c]) for c in range(chains)]))
    tau_s = float(np.mean([tau_int(Sb[c, :, 0]) for c in range(chains)]))

    C = correlators(S[:, 0, :], Lt)
    C0 = float(C[0])
    m_fit, A_fit, edge = fit_A(C, Lt, ts=np.arange(1, 7))
    f_fix = lambda Cc: A_at(Cc, Lt, M1546) / Cc[0]
    f_ft = lambda Cc: fit_A(Cc, Lt, ts=np.arange(1, 7))[1] / Cc[0]
    w_fix, w_fix_e = jack(f_fix, S)
    w_ft, w_ft_e = jack(f_ft, S)
    _, s1 = jack(lambda Cc: Cc[1], S)
    _, s2 = jack(lambda Cc: Cc[2], S)
    snr2 = C[2] / s2 if s2 > 0 else float('inf')

    res = dict(beta=float(beta), N=N, chains=chains, per=per, Lt=Lt, gap=gap, therm=therm,
               seed=seed, sec=dtw, sec_compile=t_comp,
               plaq=float(plaq.mean()), plaq_err=float(plaq.std(ddof=1) / np.sqrt(N)),
               C=[float(x) for x in C], C1o0=float(C[1] / C0), C2o0=float(C[2] / C0),
               C1e=float(s1 / C0), C2e=float(abs(s2) / C0),
               m_fit=float(m_fit), m_edge=bool(edge),
               w1=float(w_fix), w1_err=float(w_fix_e),
               w1_fit=float(A_fit / C0), w1_fit_err=float(w_ft_e),
               snr2=float(snr2), snr_w=float(abs(w_fix) / w_fix_e if w_fix_e > 0 else np.inf),
               tau_plaq=tau_p, tau_S0=tau_s,
               positive=bool(w_fix - 3 * w_fix_e > 0))
    return res, S, plaq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--betas', type=float, nargs='+', default=[3.0, 4.0, 6.0, 10.0])
    ap.add_argument('--n', type=int, default=24000)
    ap.add_argument('--chains', type=int, default=8)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--budget', type=float, default=900.0)
    ap.add_argument('--seed', type=int, default=20261125)
    ap.add_argument('--tag', type=str, default='scan')
    ap.add_argument('--gap', type=int, default=0)
    ap.add_argument('--therm', type=int, default=0)
    a = ap.parse_args()

    import numba
    P("=" * 96)
    P("B01  SU(2) Wilson  w1(raw) vs beta   Ls=8^3 x Lt=12   operator = run1/run3 lev0 (unchanged)")
    P(f"m = {M1546} fixed ; w1 = A/C0 ; errors = 20-block delete-one jackknife ; "
      f"numba={numba.__version__} threads={a.threads}")
    P(f"target N = {a.n}/point ; gap,therm scale with beta ; strictly sequential (one process)")
    P("=" * 96)

    plan = {3.0: (5, 900), 4.0: (6, 1200), 6.0: (9, 2000), 10.0: (14, 3000)}
    T0 = time.time()
    rows = []
    for i, beta in enumerate(a.betas):
        gap, therm = plan.get(float(beta), (6, 1200))
        gap = a.gap or gap
        therm = a.therm or therm
        left = a.budget - (time.time() - T0)
        npts = len(a.betas) - i
        n = max(a.chains, (a.n // a.chains) * a.chains)
        if rows:
            spr = min(r['N'] * r['gap'] / r['sec'] for r in rows)     # sweeps/s, worst case
            est = n * gap / spr
            allow = max(1.0, left / npts)
            if est > allow:                     # cut N (never beta) if the wall clock forces it
                n = int(max(2 * a.chains, (allow * spr / gap) // a.chains * a.chains))
                P(f"[budget] beta={beta:g}: projected {est:.0f}s > {allow:.0f}s -> N={n}")
        P("")
        P(f"--- beta = {beta:g}  N = {n} ({a.chains} x {n//a.chains})  gap = {gap}  "
          f"therm = {therm}   (t+{time.time()-T0:.0f}s)")
        r, S, plaq = point(float(beta), n, gap, therm, a.chains, a.threads,
                           a.seed + int(1000 * beta))
        rows.append(r)
        P(f"[E] {r['N']} cfg in {r['sec']:.1f}s ({r['N']/r['sec']:.1f} cfg/s, "
          f"{r['N']*gap/r['sec']:.0f} sweeps/s)  compile {r['sec_compile']:.1f}s  "
          f"tau_int: plaq {r['tau_plaq']:.1f}, S(0) {r['tau_S0']:.1f} cfg")
        P(f"[E] <1/2 Re tr U_p> = {r['plaq']:.5f} +- {r['plaq_err']:.5f}")
        P(f"[E] C1/C0 = {r['C1o0']:.5f} +- {r['C1e']:.5f} ; C2/C0 = {r['C2o0']:.5f} "
          f"+- {r['C2e']:.5f} ; S/N(2) = {r['snr2']:.1f}")
        P(f"[E] m_fit = {r['m_fit']:.3f}{' (edge)' if r['m_edge'] else ''} ; "
          f"w1(m={M1546}) = {r['w1']:.4f} +- {r['w1_err']:.4f}  S/N = {r['snr_w']:.1f} ; "
          f"w1_fit = {r['w1_fit']:.4f} +- {r['w1_fit_err']:.4f}")
        np.savez(os.path.join(HERE, f"b01_b{beta:g}.npz"), S=S, plaq=plaq,
                 beta=np.float64(beta), N=np.int64(r['N']), Lt=np.int64(r['Lt']),
                 gap=np.int64(gap), therm=np.int64(therm), chains=np.int64(r['chains']),
                 seed=np.int64(r['seed']), w1=np.float64(r['w1']),
                 w1_err=np.float64(r['w1_err']))
        json.dump(rows, open(os.path.join(HERE, f"b01_{a.tag}.json"), 'w'), indent=1)
        del S, plaq

    P("")
    P("=" * 96)
    P(f"{'beta':>5} {'N':>7} {'<plaq>':>9} {'tau':>5} {'C1/C0':>9} {'C2/C0':>9} "
      f"{'m_fit':>6} {'w1':>8} {'+-':>7} {'S/N':>5} {'>0?':>4} {'s':>6}")
    for r in rows:
        P(f"{r['beta']:>5.1f} {r['N']:>7d} {r['plaq']:>9.5f} {r['tau_plaq']:>5.1f} "
          f"{r['C1o0']:>9.5f} {r['C2o0']:>9.5f} {r['m_fit']:>6.3f} {r['w1']:>8.4f} "
          f"{r['w1_err']:>7.4f} {r['snr_w']:>5.1f} {'yes' if r['positive'] else 'NO':>4} "
          f"{r['sec']:>6.1f}")
    P("")
    ok = all(rows[j]['plaq'] < rows[j + 1]['plaq'] -
             3 * np.hypot(rows[j]['plaq_err'], rows[j + 1]['plaq_err'])
             for j in range(len(rows) - 1))
    P(f"[check] <plaq> strictly increasing in beta (3 sigma step): {'PASS' if ok else 'FAIL'}")
    P(f"[wall] total {time.time()-T0:.0f}s")
    json.dump(rows, open(os.path.join(HERE, f"b01_{a.tag}.json"), 'w'), indent=1)
    open(os.path.join(HERE, f"b01_{a.tag}.log"), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
