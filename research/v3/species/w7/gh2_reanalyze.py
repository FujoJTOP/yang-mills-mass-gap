"""GH2 re-analysis of ALREADY-STORED configurations (no Markov chain, no new sampling).

Purpose: put w1 and Rbar on the same footing before comparing the paper's 8.3% (a w1 ratio,
RUN1.md:47,52 / RUN3.md:51) with GW7's 0.1% (an Rbar difference, GW7-FIXEDA-SCAN.md:34).

Nothing here is a new measurement: it reads the archived S(t) arrays and re-derives
C(t)/C(0), w1(m=1.546), kappa = w1/(e^m C1/C0), Rbar with the *same* estimator
(run3_apebase.correlators / A_at, 20-block delete-one jackknife).

Run:  py -3.12 gh2_reanalyze.py
"""
import os
import sys
import json
import glob

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from run3_apebase import correlators, A_at, fit_A  # noqa: E402

M = 1.546
NB = 20
FILES = ['run1_b25_l6.npz', 'run1_b2.5_Ls8Lt12.npz', 'run3_L8_b2.5_N12k.npz',
         'gw7_b2.59_Ls10.npz', 'gw7_b2.66_Ls12.npz']
OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def one(path):
    d = np.load(os.path.join(HERE, path))
    Lt = int(d['Lt'])
    S = d['S']
    Sm = S[:, 0, :] if S.ndim == 3 else S
    N = Sm.shape[0]
    C = correlators(Sm, Lt)
    idx = np.array_split(np.arange(N), NB)
    Cs = np.array([correlators(Sm[np.concatenate([idx[j] for j in range(NB) if j != b])], Lt)
                   for b in range(NB)])

    def jk(fn):
        v = np.array([fn(Cs[b]) for b in range(NB)], float)
        m = float(np.nanmean(v))
        e = float(np.sqrt((NB - 1) / NB * np.nansum((v - m) ** 2)))
        return m, e

    f_w = lambda c: A_at(c, Lt, M) / c[0]
    f_R = lambda c: -np.log(max(c[1] / c[0], 1e-12))
    f_k = lambda c: (A_at(c, Lt, M) / c[0]) / (np.exp(M) * c[1] / c[0])
    w1, w1e = jk(f_w)
    Rb, Rbe = jk(f_R)
    ka, kae = jk(f_k)
    r = [C[t] / C[0] for t in range(Lt)]
    re_ = {}
    for t in range(1, Lt // 2 + 1):
        m, e = jk(lambda c, t=t: c[t] / c[0])
        re_[t] = (m, e)
    ceil = float(np.exp(M) * C[1] / C[0])
    mfit, _, edge = fit_A(C, Lt, ts=np.arange(1, min(7, Lt)))
    beta = float(d['beta']) if 'beta' in d.files else float('nan')
    Ls = int(d['Ls']) if 'Ls' in d.files else -1
    return dict(file=path, beta=beta, Ls=Ls, Lt=Lt, N=N, C0=float(C[0]),
                r=[float(x) for x in r], r_err={str(k): v[1] for k, v in re_.items()},
                w1=w1, w1e=w1e, Rbar=Rb, Rbar_e=Rbe, kappa=ka, kappa_e=kae,
                ceiling=ceil, m_fit=mfit, edge=edge)


def main():
    P("=" * 104)
    P("GH2 / S1   archived-configuration re-analysis   (same estimator: run3_apebase.correlators+A_at)")
    P(f"m = {M} pinned, fit window t=1,2,3, {NB}-block delete-one jackknife.  NO new sampling.")
    P("=" * 104)
    rows = []
    for f in FILES:
        if not os.path.exists(os.path.join(HERE, f)):
            P(f"[TOOL-FAIL] missing {f}")
            continue
        r = one(f)
        rows.append(r)
        P("")
        P(f"--- {f}   beta={r['beta']:g} Ls={r['Ls']} Lt={r['Lt']} N={r['N']}")
        P(f"    C(0) = {r['C0']:.4f} ; m_fit = {r['m_fit']:.3f}{' EDGE' if r['edge'] else ''}")
        P("    C(t)/C(0), t=0..%d : %s" % (min(r['Lt'] - 1, 11),
          "  ".join(f"{r['r'][t]:.5f}" for t in range(min(r['Lt'], 12)))))
        P(f"    Rbar = {r['Rbar']:.4f} +- {r['Rbar_e']:.4f} ; w1 = {r['w1']:.5f} +- {r['w1e']:.5f}")
        P(f"    ceiling e^m*C1/C0 = {r['ceiling']:.5f} ; kappa = w1/ceiling = "
          f"{r['kappa']:.4f} +- {r['kappa_e']:.4f}")

    P("")
    P("=" * 104)
    P(f"{'file':>26} {'beta':>5} {'Ls':>3} {'Lt':>3} {'N':>6} {'Rbar':>7} {'w1':>8} {'kap':>6} "
      f"{'w1/ceil':>7}")
    for r in rows:
        P(f"{r['file']:>26} {r['beta']:>5.3f} {r['Ls']:>3d} {r['Lt']:>3d} {r['N']:>6d} "
          f"{r['Rbar']:>7.3f} {r['w1']:>8.4f} {r['kappa']:>6.3f} {r['kappa']:>7.3f}")

    # ---- dimension conversion: w1 = kappa * exp(m - Rbar)  =>  dln w1 = dln kappa - dRbar
    P("")
    P("--- dimension conversion  w1 = kappa*exp(m)*C1/C0 = kappa*exp(m-Rbar) ---")
    byf = {r['file']: r for r in rows}
    pairs = [('run1_b25_l6.npz', 'run1_b2.5_Ls8Lt12.npz',
              'the 8.3% pair as the paper computes it (both rows from RUN1)'),
             ('run1_b25_l6.npz', 'run3_L8_b2.5_N12k.npz',
              'L_s=6 from RUN1 vs the N=12000 L_s=8 row the paper headlines'),
             ('run3_L8_b2.5_N12k.npz', 'gw7_b2.59_Ls10.npz', 'GW7 style: 8x12 -> 10x15'),
             ('run3_L8_b2.5_N12k.npz', 'gw7_b2.66_Ls12.npz', 'GW7 style: 8x12 -> 12x18')]
    for f1, f2, lab in pairs:
        A, B = byf.get(f1), byf.get(f2)
        if A is None or B is None:
            continue
        dlw = np.log(A['w1'] / B['w1'])
        dR = B['Rbar'] - A['Rbar']
        dlk = np.log(A['kappa'] / B['kappa'])
        P(f"[T] {lab}")
        P(f"    ln(w1_A/w1_B) = {dlw:+.4f}  =  Rbar_B-Rbar_A ({dR:+.4f})  +  ln(kappa_A/kappa_B) "
          f"({dlk:+.4f})   residual {dlw - dR - dlk:+.2e}")
        P(f"    => w1 ratio {100*(A['w1']/B['w1']-1):+.2f}%  <->  Delta Rbar {dR:+.4f} "
          f"= {100*dR/A['Rbar']:+.2f}% of Rbar ; kappa carries {100*dlk:+.2f}% of the log")

    json.dump(rows, open(os.path.join(HERE, 'gh2_reanalyze.json'), 'w'), indent=1)
    open(os.path.join(HERE, 'gh2_reanalyze.log'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
