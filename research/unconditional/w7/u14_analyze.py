"""U14 analysis — w1(L) with the CORPUS estimator (m fixed by the smeared level).

Why this file exists: a *free* cosh fit on the unsmeared (s=0) correlator runs to the
m-grid edge (m=6.0) and returns garbage; the archive p8 numbers were made with m FIXED
to the plateau mass of the most-smeared level (`_f1sq_su2_analyze.analyse_C`).  This
script reproduces that estimator for every L, so the L-trend is estimator-consistent.

  m_G   := mean of cosh m_eff(t>=2) of the highest APE level that has a positive signal
  A_s   := <w|C_s>/<w|w>,  w(t)=e^{-m_G t}+e^{-m_G(T-t)},  t=tmin..T/2
  M0_s  := C_s(0)
  w1_s  := A_s/M0_s        (Z_O-invariant: O -> lam O cancels)
  w1hat := e^{m_G t} C_s(t)/C_s(0)  >=  w1   (T2 sec 2.3, decreasing in t)
Errors: block bootstrap over configurations (400 resamples).
"""
import sys

import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/unconditional')
import _f1sq_su2_analyze as A     # noqa: E402


def conn(S):
    T = S.shape[1]
    mu = S.mean(axis=1)
    P = np.array([np.mean(S * np.roll(S, -t, axis=1), axis=1) for t in range(T)]).T
    m = mu.mean()
    return P.mean(axis=0) - m * m


def meff(C, T):
    out = {}
    for t in range(1, T // 2 + 1):
        if t + 1 <= T - 1 and min(C[t - 1], C[t], C[t + 1]) > 0:
            r = (C[t - 1] + C[t + 1]) / (2.0 * C[t])
            out[t] = float(np.arccosh(r)) if r >= 1.0 else np.nan
    return out


def point(Cs, levels):
    """m_G from the most-smeared level, then fixed-m amplitudes for every level."""
    T = Cs[levels[0]].shape[0]
    mG = np.nan
    for s in sorted(levels, reverse=True):
        v = [meff(Cs[s], T)[t] for t in sorted(meff(Cs[s], T))
             if t >= 2 and np.isfinite(meff(Cs[s], T)[t])]
        if v:
            mG = float(np.mean(v)); break
    if not np.isfinite(mG):
        for s in levels:
            v = [x for x in meff(Cs[s], T).values() if np.isfinite(x)]
            if v:
                mG = float(v[-1]); break
    tmin = 2 if T >= 6 else 1
    ts = np.arange(tmin, T // 2 + 1)
    w = np.exp(-mG * ts) + np.exp(-mG * (T - ts)) if np.isfinite(mG) else None
    out = {}
    for s in levels:
        C = Cs[s]
        A_s = float((w @ C[ts]) / (w @ w)) if w is not None else np.nan
        out[s] = dict(A=A_s, w1=A_s / C[0] if C[0] != 0 else np.nan,
                      what={int(t): float(np.exp(mG * t) * C[t] / C[0])
                            for t in range(1, T // 2 + 1)} if np.isfinite(mG) else {})
    return mG, out


def run(S, levels, nboot=400, seed=7):
    T = S[levels[0]].shape[1]
    N = S[levels[0]].shape[0]
    Cs = {s: conn(S[s]) for s in levels}
    mG, cen = point(Cs, levels)
    rng = np.random.default_rng(seed)
    block = max(1, N // 20)
    bs = {s: [] for s in levels}
    bm = []
    for idx in A.block_boot_indices(N, nboot, block, rng):
        Cb = {s: conn(S[s][idx]) for s in levels}
        if any(Cb[s][0] == 0 for s in levels):
            continue
        mb, ob = point(Cb, levels)
        bm.append(mb)
        for s in levels:
            bs[s].append(ob[s]['w1'])
    return dict(T=T, N=N, mG=mG, mG_se=float(np.nanstd(bm)), cen=cen,
                se={s: float(np.nanstd(bs[s])) for s in levels},
                C0={s: float(Cs[s][0]) for s in levels},
                C1C0={s: float(Cs[s][1] / Cs[s][0]) for s in levels},
                snap=sorted(bs[levels[-1]])[:5])


def show(tag, beta, L, S, levels):
    r = run(S, levels)
    print(f'\n=== {tag}  beta={beta}  L={L}  N={r["N"]}  T={r["T"]}  '
          f'm_G={r["mG"]:.3f}+-{r["mG_se"]:.3f} ===')
    print(f'  {"lvl":>3s} {"M0=C(0)":>9s} {"C1/C0":>7s} {"A=f1^2":>9s} {"w1=A/M0":>15s} '
          f'{"w1hat(t)":>20s}')
    for s in levels:
        c = r['cen'][s]
        wh = ' '.join('%.2f' % c['what'][t] for t in sorted(c['what']))
        print(f'  {s:3d} {r["C0"][s]:9.4g} {r["C1C0"][s]:7.4f} {c["A"]:9.4g} '
              f'{c["w1"]:8.3f}+-{r["se"][s]:.3f} {wh:>20s}')
    return r


def main():
    d6 = np.load('D:/Dev/math/research/unconditional/f1sq_b2.2_L6.npz')
    lv6 = [int(x) for x in d6['levels']]
    r6 = show('archive(N=4000)', 2.2, 6, {s: d6[f'S{s}'] for s in lv6}, lv6)

    res = {(2.2, 6): r6}
    for L in (4, 8):
        f = f'D:/Dev/math/research/unconditional/w7/u14_b2.2_L{L}_b2.2.npz'
        try:
            d = np.load(f)
        except FileNotFoundError:
            print(f'\n[MISSING] {f} -- fresh run did not finish'); continue
        lv = [int(x) for x in d['levels']]
        res[(2.2, L)] = show(f'fresh(N={int(d["N"])})', 2.2, L,
                             {s: d[f'S{s}'] for s in lv}, lv)

    d4b = 'D:/Dev/math/research/unconditional/w7/u14_b2.5_L4_b2.5.npz'
    try:
        d = np.load(d4b)
        lv = [int(x) for x in d['levels']]
        res[(2.5, 4)] = show('fresh (DECONFINED side)', 2.5, 4,
                             {s: d[f'S{s}'] for s in lv}, lv)
    except FileNotFoundError:
        print(f'\n[MISSING] {d4b}')

    p = np.load('D:/Dev/math/research/lattice-cpb/pilot_b2.5.npz')
    lvp = [int(x) for x in p['smears']]
    res[(2.5, 8)] = show('pilot(N=40)', 2.5, 8,
                         {s: p[f'O_{s}'] for s in lvp}, lvp)

    print('\n=============== U14 SUMMARY: w1(L), level s=0 ===============')
    print(f'  {"beta":>5s} {"L":>2s} {"N":>5s} {"m_G":>5s} {"M0":>9s} {"w1":>7s} '
          f'{"+-":>6s} {"w1hat(1)":>8s} {"w1hat(2)":>8s}')
    for (b, L) in sorted(res):
        r = res[(b, L)]
        if 0 not in r['cen']:
            continue
        c = r['cen'][0]
        print(f'  {b:5.1f} {L:2d} {r["N"]:5d} {r["mG"]:5.2f} {r["C0"][0]:9.4g} '
              f'{c["w1"]:7.3f} {r["se"][0]:6.3f} {c["what"].get(1, np.nan):8.3f} '
              f'{c["what"].get(2, np.nan):8.3f}')

    print('\n--- same, level s=2 (smeared; least contact pollution) ---')
    for (b, L) in sorted(res):
        r = res[(b, L)]
        if 2 not in r['cen']:
            continue
        c = r['cen'][2]
        print(f'  beta={b:4.1f} L={L:2d} N={r["N"]:4d} m_G={r["mG"]:5.2f} '
              f'w1={c["w1"]:7.3f}+-{r["se"][2]:.3f}  w1hat=' +
              ' '.join('%.2f' % c['what'][t] for t in sorted(c['what'])))

    print('\n--- M0 and f1^2 L-scaling at beta=2.2 (s=0) ---')
    for L in (4, 6, 8):
        r = res.get((2.2, L))
        if r:
            print(f'  L={L}: M0={r["C0"][0]:.4g}  A={r["cen"][0]["A"]:.4g}')


if __name__ == '__main__':
    main()
