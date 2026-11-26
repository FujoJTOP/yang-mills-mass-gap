"""GH2 -- L_t artefact scan: is the paper's 8.3% "volume dependence" an L_t = 12 artefact?

Zero new physics, zero new operator.  Everything below reuses the *unmodified* kernel and
estimator of run3_apebase:

    S(t)      = (1/3) sum_{space slice t} sum_{i<j<=3} (1/2) Re tr U_p       (raw, nlev = 1)
    C(t)      = (1/Lt) sum_{t0} <S(t0+t) S(t0)>_c                            (periodic, t0-averaged)
    Rbar      = -ln C(1)/C(0)
    w1        = A/C(0),  A = sum_t g_t C_t / sum_t g_t^2 ,  g_t = e^{-m t} + e^{-m (Lt-t)}
                window t = 1,2,3 ,  m = 1.546 pinned
    m_eff(t)  = ln C(t)/C(t+1)                                               (energy background)
    errors    = 20-block delete-one jackknife (same block layout as run3/gw7)

The only thing this file varies that run1/run3/b01/gw7 did not is (Ls, Lt) *independently* and
*at fixed beta*.  gw7_fixeda.py already parametrised (Ls, Lt) but keyed both off beta along the
fixed-physical-volume line; it also names its output by Ls only, so it cannot hold two Lt at one
Ls.  Hence a new thin driver -- run3_apebase.py is imported, never edited.

Run:  py -3.12 gh2_ltscan.py --points 2.5:6:12:24000 2.5:8:12:24000 2.5:6:20:24000 2.5:8:20:24000
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
NB = 20
OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def tau_int(x):
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


# ------------------------------------------------------------------ estimators
# All of these take a correlator vector Cc of length Lt and are applied both to the full
# sample (central value) and to the 20 delete-one blocks (error).
def e_Rbar(Cc, Lt):
    return -np.log(max(Cc[1] / Cc[0], 1e-12))


def e_w1(Cc, Lt):
    return A_at(Cc, Lt, M1546) / Cc[0]


def e_kappa(Cc, Lt):
    """w1 / (e^m C1/C0): GW2-NUMERICAL-PRIOR.md:78's 'ceiling saturation' ratio."""
    return (A_at(Cc, Lt, M1546) / Cc[0]) / (np.exp(M1546) * Cc[1] / Cc[0])


def e_w1_nowrap(Cc, Lt):
    """Same fit with the anti-periodic image term dropped: isolates the wrap term's weight."""
    ts = np.arange(1, 4)
    g = np.exp(-M1546 * ts)
    return float(np.sum(g * Cc[ts]) / np.sum(g * g)) / Cc[0]


def e_meff(Cc, Lt, t):
    return np.log(max(Cc[t] / Cc[t + 1], 1e-300))


def e_r(Cc, Lt, t):
    return Cc[t] / Cc[0]


def jack_all(Sm, Lt, nb=NB):
    """Evaluate correlators once per delete-one block; derive every quantity from those."""
    N = Sm.shape[0]
    idx = np.array_split(np.arange(N), nb)
    Cs = np.empty((nb, Lt))
    for b in range(nb):
        keep = np.concatenate([idx[j] for j in range(nb) if j != b])
        Cs[b] = correlators(Sm[keep], Lt)
    return correlators(Sm, Lt), Cs


def jkm(fn, Cs, Lt):
    v = np.array([fn(Cs[b], Lt) for b in range(Cs.shape[0])], float)
    m = float(np.nanmean(v))
    e = float(np.sqrt((Cs.shape[0] - 1) / Cs.shape[0] * np.nansum((v - m) ** 2)))
    return m, e


# ------------------------------------------------------------------ one point
def point(beta, Ls, Lt, n, gap, therm, chains, threads, seed):
    import numba
    numba.set_num_threads(threads)
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
            np.zeros((2, 1, Lt)), np.zeros(2))
    t_comp = time.time() - t0

    t1 = time.time()
    _chains(U, rng, plus, minus, beta / 2.0, per, gap, therm, 0.5, 1, Lt, S, plaq)
    dtw = time.time() - t1

    Sb = S[:, 0, :].reshape(chains, per, Lt)
    pb = plaq.reshape(chains, per)
    tau_p = float(np.mean([tau_int(pb[c]) for c in range(chains)]))
    tau_s = float(np.mean([tau_int(Sb[c, :, 0]) for c in range(chains)]))

    C, Cs = jack_all(S[:, 0, :], Lt)
    C0 = float(C[0])

    r_t, r_t_e, meff, meff_e = {}, {}, {}, {}
    for t in range(1, Lt // 2 + 1):
        r_t[t], r_t_e[t] = jkm(lambda Cc, L, t=t: e_r(Cc, L, t), Cs, Lt)
    for t in range(0, Lt // 2):
        meff[t], meff_e[t] = jkm(lambda Cc, L, t=t: e_meff(Cc, L, t), Cs, Lt)

    w1, w1_e = jkm(e_w1, Cs, Lt)
    w1_full = float(e_w1(C, Lt))
    Rbar, Rbar_e = jkm(e_Rbar, Cs, Lt)
    Rbar_full = float(e_Rbar(C, Lt))
    kap, kap_e = jkm(e_kappa, Cs, Lt)
    w1nw, w1nw_e = jkm(e_w1_nowrap, Cs, Lt)
    c1, c1_e = jkm(lambda Cc, L: Cc[1] / Cc[0], Cs, Lt)
    c2, c2_e = jkm(lambda Cc, L: Cc[2] / Cc[0], Cs, Lt)
    _, s2 = jkm(lambda Cc, L: Cc[2], Cs, Lt)
    snr2 = C[2] / s2 if s2 > 0 else float('inf')
    m_fit, A_fit, edge = fit_A(C, Lt, ts=np.arange(1, min(7, Lt)))

    res = dict(beta=float(beta), Ls=int(Ls), Lt=int(Lt), N=int(N), chains=chains, per=per,
               gap=gap, therm=therm, seed=int(seed), sec=float(dtw), sec_compile=float(t_comp),
               plaq=float(plaq.mean()), plaq_err=float(plaq.std(ddof=1) / np.sqrt(N)),
               C=[float(x) for x in C], C0=C0,
               r={str(k): v for k, v in r_t.items()}, r_e={str(k): v for k, v in r_t_e.items()},
               meff={str(k): v for k, v in meff.items()},
               meff_e={str(k): v for k, v in meff_e.items()},
               C1o0=float(c1), C1o0_err=float(c1_e), C2o0=float(c2), C2o0_err=float(c2_e),
               Rbar=float(Rbar), Rbar_err=float(Rbar_e), Rbar_full=Rbar_full,
               w1=float(w1), w1_err=float(w1_e), w1_full=w1_full,
               kappa=float(kap), kappa_err=float(kap_e),
               w1_nowrap=float(w1nw), w1_nowrap_err=float(w1nw_e),
               m_fit=float(m_fit), m_edge=bool(edge), snr2=float(snr2),
               tau_plaq=tau_p, tau_S0=tau_s)
    return res, S, plaq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--points', type=str, nargs='+',
                    default=['2.5:6:12:24000', '2.5:8:12:24000', '2.5:6:20:24000', '2.5:8:20:24000'],
                    help='beta:Ls:Lt:N ...')
    ap.add_argument('--chains', type=int, default=8)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--gap', type=int, default=3)
    ap.add_argument('--therm', type=int, default=400)
    ap.add_argument('--budget', type=float, default=540.0)
    ap.add_argument('--seed', type=int, default=20261125)
    ap.add_argument('--tag', type=str, default='ltscan')
    a = ap.parse_args()

    pts = []
    for p in a.points:
        q = p.split(':')
        pts.append((float(q[0]), int(q[1]), int(q[2]), int(q[3])))

    import numba
    P("=" * 100)
    P("GH2  L_t-ARTEFACT SCAN   SU(2) Wilson  beta fixed, (Ls, Lt) varied independently")
    P("kernel/operator/estimator = run3_apebase.py VERBATIM (_chains, correlators, fit_A, A_at); "
      "nlev=1 raw")
    P(f"m = {M1546} pinned ; fit window t=1,2,3 ; errors = {NB}-block delete-one jackknife ; "
      f"numba={numba.__version__} threads={a.threads}")
    P(f"gap={a.gap} therm={a.therm} chains={a.chains} ; strictly sequential, one process ; "
      f"budget={a.budget:.0f}s")
    P("=" * 100)

    T0 = time.time()
    rows = []
    for beta, Ls, Lt, n in pts:
        seed = a.seed + int(1000 * beta) + 101 * Ls + 1009 * Lt
        P("")
        P(f"--- beta = {beta:g}  Ls={Ls}  Lt={Lt}  N = {n} ({a.chains} x {n//a.chains})  "
          f"gap = {a.gap}  therm = {a.therm}   (t+{time.time()-T0:.0f}s)")
        r, S, plaq = point(beta, Ls, Lt, n, a.gap, a.therm, a.chains, a.threads, seed)
        rows.append(r)
        P(f"[E] {r['N']} cfg in {r['sec']:.1f}s ({r['N']/r['sec']:.1f} cfg/s, "
          f"{r['N']*a.gap/r['sec']:.0f} sweeps/s)  compile {r['sec_compile']:.1f}s  "
          f"tau_int: plaq {r['tau_plaq']:.1f}, S(0) {r['tau_S0']:.1f} cfg")
        P(f"[E] <1/2 Re tr U_p> = {r['plaq']:.5f} +- {r['plaq_err']:.5f}")
        P(f"[E] C(0) = {r['C0']:.4f}   C(t)/C(0), t=0..4 : "
          + "  ".join(f"1.00000" if t == 0 else
                      (f"{r['r'][str(t)]:.5f}+-{r['r_e'][str(t)]:.5f}" if str(t) in r['r'] else "n/a")
                      for t in range(0, 5)))
        P(f"[E] C1/C0 = {r['C1o0']:.5f} +- {r['C1o0_err']:.5f} ; C2/C0 = {r['C2o0']:.5f} ; "
          f"S/N(2) = {r['snr2']:.1f}")
        P(f"[E] Rbar = -ln C1/C0 = {r['Rbar']:.4f} +- {r['Rbar_err']:.4f} (jack) ; "
          f"full-sample {r['Rbar_full']:.4f}")
        P(f"[E] w1(1.546) = {r['w1']:.5f} +- {r['w1_err']:.5f} (jack) ; full-sample "
          f"{r['w1_full']:.5f} ; w1_nowrap = {r['w1_nowrap']:.5f} +- {r['w1_nowrap_err']:.5f}")
        P(f"[E] kappa = w1/(e^m C1/C0) = {r['kappa']:.4f} +- {r['kappa_err']:.4f}")
        P(f"[E] m_fit(t=1..{min(6, Lt-1)}) = {r['m_fit']:.3f}{' (EDGE)' if r['m_edge'] else ''} ; "
          f"m_eff(t)=ln C(t)/C(t+1): "
          + "  ".join(f"{r['meff'][str(t)]:.3f}" for t in range(0, min(Lt // 2, 10))))
        np.savez(os.path.join(HERE, f"gh2_b{beta:g}_Ls{Ls}_Lt{Lt}.npz"), S=S, plaq=plaq,
                 beta=np.float64(beta), Ls=np.int64(Ls), Lt=np.int64(Lt),
                 N=np.int64(r['N']), gap=np.int64(a.gap), therm=np.int64(a.therm),
                 seed=np.int64(seed))
        del S, plaq
        json.dump(rows, open(os.path.join(HERE, f"gh2_{a.tag}.json"), 'w'), indent=1)
        P(f"[wall] elapsed {time.time()-T0:.0f}s / budget {a.budget:.0f}s")

    P("")
    P("=" * 100)
    P(f"{'beta':>5} {'Ls':>3} {'Lt':>3} {'N':>6} {'<plaq>':>9} {'C1/C0':>9} {'Rbar':>7} {'+-':>6} "
      f"{'w1':>8} {'+-':>7} {'kap':>6} {'S/N2':>5} {'s':>5}")
    for r in rows:
        P(f"{r['beta']:>5.3f} {r['Ls']:>3d} {r['Lt']:>3d} {r['N']:>6d} {r['plaq']:>9.5f} "
          f"{r['C1o0']:>9.5f} {r['Rbar']:>7.3f} {r['Rbar_err']:>6.3f} {r['w1']:>8.4f} "
          f"{r['w1_err']:>7.4f} {r['kappa']:>6.3f} {r['snr2']:>5.1f} {r['sec']:>5.1f}")

    # ---- the actual verdict arithmetic: L_s gap at each fixed L_t, and L_t drift at each L_s
    P("")
    P("--- 8.3% recomputed at each available L_t  (w1, same estimator, same m) ---")
    for Lt in sorted(set(r['Lt'] for r in rows)):
        sub = {r['Ls']: r for r in rows if r['Lt'] == Lt}
        if 6 in sub and 8 in sub:
            a6, a8 = sub[6], sub[8]
            d = a6['w1'] - a8['w1']
            e = np.hypot(a6['w1_err'], a8['w1_err'])
            P(f"[J] Lt={Lt:2d}: w1(6)-w1(8) = {d:+.4f} +- {e:.4f}  = "
              f"{100*d/(a6['w1']+a8['w1'])*2:+.2f}% (of mean) / {100*d/a6['w1']:+.2f}% (of Ls=6)")
            P(f"           Rbar(8)-Rbar(6) = {a8['Rbar']-a6['Rbar']:+.4f} +- "
              f"{np.hypot(a6['Rbar_err'], a8['Rbar_err']):.4f}")
    P("")
    P("--- L_t drift at each L_s ---")
    for Ls in sorted(set(r['Ls'] for r in rows)):
        sub = sorted([r for r in rows if r['Ls'] == Ls], key=lambda r: r['Lt'])
        for i in range(len(sub) - 1):
            p, q = sub[i], sub[i + 1]
            dw = q['w1'] - p['w1']
            ew = np.hypot(p['w1_err'], q['w1_err'])
            dR = q['Rbar'] - p['Rbar']
            eR = np.hypot(p['Rbar_err'], q['Rbar_err'])
            P(f"[J] Ls={Ls}: Lt {p['Lt']:2d}->{q['Lt']:2d}  dw1 = {dw:+.4f} +- {ew:.4f} "
              f"({100*dw/p['w1']:+.2f}%)   dRbar = {dR:+.4f} +- {eR:.4f} "
              f"({100*dR/p['Rbar']:+.2f}%)")

    P("")
    P(f"[wall] total MC+analysis {time.time()-T0:.0f}s")
    json.dump(rows, open(os.path.join(HERE, f"gh2_{a.tag}.json"), 'w'), indent=1)
    open(os.path.join(HERE, f"gh2_{a.tag}.log"), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
