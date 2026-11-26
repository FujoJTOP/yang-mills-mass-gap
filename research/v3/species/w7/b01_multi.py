"""B01 — w1(raw) at beta = 3, 4, 6, 10 IN ONE PASS.

All four betas are independent, so they go into a single prange launch:
    32 chains = 4 betas x 8 chains, one chain per core (cpu_count = 32).
Per-chain kappa / gap / therm come from arrays, so one compiled kernel covers all 4 betas.
Wall time ~= the cost of ONE chain (1500 cfg x gap sweeps), not 4 x that.

Operator (unchanged from run1/run3, lev0 "raw"):
    S(t) = (1/3) sum_{space} sum_{3 spatial planes} (1/2) Re tr U_p   (zero momentum, slice t)
    C(t) = (1/Lt) sum_{t0} <S(t0+t) S(t0)>_c
    w1   = A / C(0),  A = sum_{t=1..3} g_t C(t) / sum g_t^2,  g_t = e^{-m t} + e^{-m(Lt-t)}
    m = 1.546 fixed  (= run3's headline estimator, 0.3427(119) at beta=2.5)
    w1_fit: free-mass variant (m by least squares on t=1..6)
Errors: delete-one 20-block jackknife over configurations.

BUG FIX vs the first draft: the jackknife must receive a 2-D (N, Lt) array; run3's
correlators() sums over axis=1, so a 3-D (N,1,Lt) argument silently correlates the
size-1 level axis instead of time (it returned a constant vector -> zero error bars).

Run:  py -3.12 b01_multi.py --n 12000
"""
import os
import sys
import time
import json
import argparse

import numpy as np
from numba import njit, prange

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from run3_apebase import _sweep, _plaq_all, _slices, index_tables, hot, correlators, fit_A, A_at  # noqa: E402

M1546 = 1.546
OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


@njit(cache=True, fastmath=True, parallel=True)
def _chains_multi(U, rng, plus, minus, kappa, gaps, therms, per, Lt, S, plaq):
    """One chain per index c; kappa/gap/therm are per-chain (so betas can be mixed)."""
    nch = U.shape[0]
    for c in prange(nch):
        _sweep(U[c], plus, minus, kappa[c], therms[c], rng, c)
        for k in range(per):
            j = c * per + k
            _sweep(U[c], plus, minus, kappa[c], gaps[c], rng, c)
            plaq[j] = _plaq_all(U[c], plus, Lt)
            _slices(U[c], plus, Lt, S[j])


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


def jack(fn, Sm, nb=20):
    """Delete-one block jackknife.  Sm MUST be 2-D (N, Lt)."""
    assert Sm.ndim == 2, Sm.shape
    idx = np.array_split(np.arange(Sm.shape[0]), nb)
    vals = np.empty(nb)
    for b in range(nb):
        keep = np.concatenate([idx[j] for j in range(nb) if j != b])
        vals[b] = fn(correlators(Sm[keep], Sm.shape[1]))
    m = float(np.nanmean(vals))
    e = float(np.sqrt((nb - 1) / nb * np.nansum((vals - m) ** 2)))
    return m, e


def analyse(Sm, plaq, Lt, beta, n_per, nch, gap, therm, sec):
    C = correlators(Sm, Lt)
    C0 = float(C[0])
    m_fit, A_fit, edge = fit_A(C, Lt, ts=np.arange(1, 7))
    f_fix = lambda Cc: A_at(Cc, Lt, M1546) / Cc[0]
    f_ft = lambda Cc: fit_A(Cc, Lt, ts=np.arange(1, 7))[1] / Cc[0]
    w_fix, w_e = jack(f_fix, Sm)
    w_ft, w_ft_e = jack(f_ft, Sm)
    _, s1 = jack(lambda Cc: Cc[1], Sm)
    _, s2 = jack(lambda Cc: Cc[2], Sm)
    pb = plaq.reshape(nch, n_per)
    tau_p = float(np.mean([tau_int(pb[c]) for c in range(nch)]))
    tau_s = float(np.mean([tau_int(Sm.reshape(nch, n_per, Lt)[c, :, 0]) for c in range(nch)]))
    snr2 = float(C[2] / s2) if s2 > 0 else float('inf')
    return dict(beta=float(beta), N=int(Sm.shape[0]), per=int(n_per), chains=int(nch),
                Lt=int(Lt), gap=int(gap), therm=int(therm), sec=float(sec),
                plaq=float(plaq.mean()),
                plaq_err=float(plaq.std(ddof=1) / np.sqrt(Sm.shape[0])),
                C=[float(x) for x in C], C1o0=float(C[1] / C0), C1e=float(s1 / C0),
                C2o0=float(C[2] / C0), C2e=float(abs(s2) / C0),
                m_fit=float(m_fit), m_edge=bool(edge),
                w1=float(w_fix), w1_err=float(w_e),
                w1_fit=float(A_fit / C0), w1_fit_err=float(w_ft_e),
                snr2=snr2, snr_w=float(abs(w_fix) / w_e) if w_e > 0 else float('inf'),
                tau_plaq=tau_p, tau_S0=tau_s,
                positive=bool(w_fix - 3 * w_e > 0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--betas', type=float, nargs='+', default=[3.0, 4.0, 6.0, 10.0])
    ap.add_argument('--n', type=int, default=12000)
    ap.add_argument('--chains', type=int, default=8)
    ap.add_argument('--threads', type=int, default=0)
    ap.add_argument('--seed', type=int, default=20261125)
    ap.add_argument('--gap', type=int, default=0)
    ap.add_argument('--therm', type=int, default=0)
    ap.add_argument('--tag', type=str, default='multi')
    a = ap.parse_args()

    import numba
    ncpu = os.cpu_count() or 8
    threads = a.threads or min(ncpu, len(a.betas) * a.chains)
    numba.set_num_threads(threads)

    Ls, Lt = 8, 12
    per = max(1, a.n // a.chains)
    nch_b = a.chains
    nch = len(a.betas) * nch_b
    N = per * nch_b
    P("=" * 100)
    P("B01  SU(2) Wilson  w1(raw) vs beta   Ls=8^3 x Lt=12   operator = run1/run3 lev0 (unchanged)")
    P(f"m = {M1546} fixed ; w1 = A/C0 ; errors = 20-block delete-one jackknife")
    P(f"cpu_count={ncpu}  numba threads={numba.get_num_threads()} ({numba.__version__})")
    P(f"ONE prange pass: {len(a.betas)} betas x {nch_b} chains = {nch} independent chains, "
      f"N={N}/beta  ({per} cfg/chain)")
    P("=" * 100)

    plan = {3.0: (4, 600), 4.0: (4, 800), 6.0: (6, 1200), 10.0: (6, 1500)}
    gaps = np.array([plan.get(float(b), (4, 600))[0] for b in a.betas], np.int64)
    thr = np.array([plan.get(float(b), (4, 600))[1] for b in a.betas], np.int64)
    if a.gap > 0:
        gaps = np.full(len(a.betas), a.gap, np.int64)
    if a.therm > 0:
        thr = np.full(len(a.betas), a.therm, np.int64)
    kap = np.array([b / 2.0 for b in a.betas], np.float64)

    plus, minus = index_tables(Ls, Lt)
    U = hot(nch, Ls, Lt, a.seed)
    rng = np.empty((nch, 2), np.uint64)
    ss = np.uint64(a.seed)
    for c in range(nch):
        ss = ss * np.uint64(6364136223846793005) + np.uint64(1442695040888963407)
        rng[c, 0] = ss | np.uint64(1)
        s1 = ss ^ np.uint64(0x9E3779B97F4A7C15)
        rng[c, 1] = s1 if s1 != np.uint64(0) else np.uint64(0x2545F4914F6CDD1D)
    S = np.zeros((nch * per, Lt))
    plaq = np.zeros(nch * per)
    P(f"[T] sites={Ls**3*Lt} links={4*Ls**3*Lt}  U={U.nbytes/1e6:.1f} MB total  "
      f"S={S.nbytes/1e6:.2f} MB   gap={list(gaps)}  therm={list(thr)}")
    for i, b in enumerate(a.betas):
        P(f"    beta={b:g}: kappa={kap[i]:g} gap={gaps[i]} therm={thr[i]} "
          f"chains {i*nch_b}..{(i+1)*nch_b-1}")

    t0 = time.time()
    _chains_multi(U[:1], rng[:1], plus, minus, kap[:1], gaps[:1], np.array([2], np.int64),
                  2, Lt, np.zeros((2, Lt)), np.zeros(2))
    P(f"[T] compile+smoke: {time.time()-t0:.1f}s")

    t1 = time.time()
    _chains_multi(U, rng, plus, minus, kap, gaps, thr, per, Lt, S, plaq)
    dtw = time.time() - t1
    P(f"[E] production: {nch} chains x {per} cfg x gap in {dtw:.1f}s "
      f"({nch*per/dtw:.1f} cfg/s wall, {N/dtw:.0f} cfg/s per beta)  "
      f"per-chain {per*int(gaps.max())+int(thr.max())} sweeps")

    rows = []
    for i, b in enumerate(a.betas):
        sl = slice(i * nch_b * per, (i + 1) * nch_b * per)
        r = analyse(S[sl], plaq[sl], Lt, float(b), per, nch_b, int(gaps[i]), int(thr[i]), dtw)
        rows.append(r)
        P("")
        P(f"--- beta = {b:g}  N = {r['N']} ({nch_b} x {per})  gap = {r['gap']}  "
          f"therm = {r['therm']}")
        P(f"[E] <1/2 Re tr U_p> = {r['plaq']:.5f} +- {r['plaq_err']:.5f}")
        P(f"[E] C/C0 t=1,2,3 = {r['C'][1]/r['C'][0]:.5f} {r['C'][2]/r['C'][0]:.5f} "
          f"{r['C'][3]/r['C'][0]:.5f}  (+- {r['C1e']:.5f} {r['C2e']:.5f})  "
          f"S/N(2)={r['snr2']:.1f}")
        P(f"[E] m_fit = {r['m_fit']:.3f}{' (edge)' if r['m_edge'] else ''} ; "
          f"w1(m={M1546}) = {r['w1']:.4f} +- {r['w1_err']:.4f}  S/N = {r['snr_w']:.1f} ; "
          f"w1_fit = {r['w1_fit']:.4f} +- {r['w1_fit_err']:.4f}")
        P(f"[E] tau_int: plaq {r['tau_plaq']:.1f} cfg, S(0) {r['tau_S0']:.1f} cfg "
          f"(block = {per//20} cfg, so bins >> tau_int)")

    P("")
    P("=" * 100)
    P(f"{'beta':>5} {'N':>6} {'<plaq>':>9} {'+-':>8} {'tau':>5} {'C1/C0':>8} {'C2/C0':>8} "
      f"{'m_fit':>6} {'w1':>8} {'+-':>7} {'S/N':>5} {'>0?':>4} {'s':>6}")
    for r in rows:
        P(f"{r['beta']:>5.1f} {r['N']:>6d} {r['plaq']:>9.5f} {r['plaq_err']:>8.5f} "
          f"{r['tau_plaq']:>5.1f} {r['C1o0']:>8.5f} {r['C2o0']:>8.5f} {r['m_fit']:>6.3f} "
          f"{r['w1']:>8.4f} {r['w1_err']:>7.4f} {r['snr_w']:>5.1f} "
          f"{'yes' if r['positive'] else 'NO':>4} {r['sec']:>6.0f}")
    P("")
    # correctness ladder: <plaq> must rise with beta and head to 1
    bs = [r['beta'] for r in rows]
    pl = [r['plaq'] for r in rows]
    pe = [r['plaq_err'] for r in rows]
    mono = all(pl[j] < pl[j + 1] - 3 * np.hypot(pe[j], pe[j + 1]) for j in range(len(pl) - 1))
    P(f"[check] <plaq> monotone increasing in beta, 3 sigma step: {'PASS' if mono else 'FAIL'}")
    P(f"[check] beta=2.5 reference 0.652 -> " +
      " ".join(f"{b:g}:{p:.3f}" for b, p in zip(bs, pl)))
    P(f"[wall] total {time.time()-t0:.0f}s (compile+production+analysis)")

    json.dump(rows, open(os.path.join(HERE, f"b01_{a.tag}.json"), 'w'), indent=1)
    np.savez(os.path.join(HERE, f"b01_{a.tag}.npz"), S=S, plaq=plaq, betas=np.array(a.betas),
             gaps=gaps, therms=thr, per=np.int64(per), chains=np.int64(nch_b),
             Lt=np.int64(Lt), seed=np.int64(a.seed))
    open(os.path.join(HERE, f"b01_{a.tag}.log"), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
