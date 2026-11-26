"""U14 — m-FIXED pass: one common mass for every (beta,L) so the L-trend is not
contaminated by a per-volume mass fit (at L=8, N=250 the free m_G runs to 0.43 and
the resulting w1 collapses; that is an estimator failure, not physics).

  m fixed at 1.50 (corpus free-m values: L=4 -> 1.33, L=6 -> 1.53; paper route 1.51)
  A_s = <w|C_s>/<w|w>, w(t)=e^{-m t}+e^{-m(T-t)}, t=tmin..T/2   (tmin=1 if T<=4)
  w1  = A_s / C_s(0) ;  w1hat(t) = e^{m t} C_s(t)/C_s(0)  >= w1
"""
import sys

import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/unconditional')
import _f1sq_su2_analyze as A     # noqa: E402

MFIX = 1.50


def conn(S):
    T = S.shape[1]
    mu = S.mean(axis=1)
    P = np.array([np.mean(S * np.roll(S, -t, axis=1), axis=1) for t in range(T)]).T
    return P.mean(axis=0) - mu.mean() ** 2


def w1_of(S, m=MFIX, tmin=None):
    C = conn(S)
    T = S.shape[1]
    tmin = (1 if T <= 4 else 2) if tmin is None else tmin
    ts = np.arange(tmin, T // 2 + 1)
    w = np.exp(-m * ts) + np.exp(-m * (T - ts))
    Aamp = float((w @ C[ts]) / (w @ w))
    return Aamp / C[0], Aamp, C[0], {int(t): float(np.exp(m * t) * C[t] / C[0])
                                     for t in range(1, T // 2 + 1)}


def table(tag, beta, L, S, levels, nboot=400):
    print(f'\n=== {tag}  beta={beta} L={L} N={S[levels[0]].shape[0]} T={S[levels[0]].shape[1]} '
          f'(m fixed {MFIX}) ===')
    print(f'  {"lvl":>3s} {"M0":>9s} {"w1":>7s} {"+-":>6s} {"w1hat(1)":>8s} {"w1hat(2)":>8s}')
    out = {}
    for s in levels:
        w1, Aamp, M0, wh = w1_of(S[s])
        rng = np.random.default_rng(9)
        N = S[s].shape[0]
        bl = max(1, N // 20)
        bs = []
        for idx in A.block_boot_indices(N, nboot, bl, rng):
            v, _, _, _ = w1_of(S[s][idx])
            bs.append(v)
        se = float(np.nanstd(bs))
        out[s] = (w1, se, wh, M0)
        print(f'  {s:3d} {M0:9.4g} {w1:7.3f} {se:6.3f} {wh.get(1, np.nan):8.3f} '
              f'{wh.get(2, np.nan):8.3f}')
    return out


def main():
    res = {}
    d6 = np.load('D:/Dev/math/research/unconditional/f1sq_b2.2_L6.npz')
    lv = [int(x) for x in d6['levels']]
    res[(2.2, 6)] = table('archive', 2.2, 6, {s: d6[f'S{s}'] for s in lv}, lv)
    for (b, L, tag) in ((2.2, 4, 'fresh N=2000'), (2.2, 8, 'fresh N=250'),
                        (2.5, 4, 'fresh N=1500 DECONF'), (2.5, 8, 'pilot N=40')):
        f = (f'D:/Dev/math/research/unconditional/w7/u14_b2.2_L{L}_b2.2.npz' if b == 2.2
             else f'D:/Dev/math/research/unconditional/w7/u14_b2.5_L{L}_b2.5.npz')
        if (b, L) == (2.5, 8):
            p = np.load('D:/Dev/math/research/lattice-cpb/pilot_b2.5.npz')
            lvp = [int(x) for x in p['smears']]
            res[(b, L)] = table(tag, b, L, {s: p[f'O_{s}'] for s in lvp}, lvp)
            continue
        try:
            d = np.load(f)
        except FileNotFoundError:
            print(f'[MISSING] {f}'); continue
        lvv = [int(x) for x in d['levels']]
        res[(b, L)] = table(tag, b, L, {s: d[f'S{s}'] for s in lvv}, lvv)

    print('\n############# U14 HEADLINE TABLE: w1 vs L (m fixed 1.50) #############')
    print(f'  {"beta":>5s} {"L":>2s} {"N":>5s} | {"s=0":>16s} {"s=2":>16s} | '
          f'{"w1hat(1) s=0":>12s} {"w1hat(2) s=0":>12s}')
    for (b, L) in sorted(res):
        r = res[(b, L)]
        N = np.nan
        for s in r:
            pass
        a = r.get(0)
        c = r.get(2)
        if a is None:
            continue
        Ns = {4: {2.2: 2000, 2.5: 1500}, 6: {2.2: 4000}, 8: {2.2: 250, 2.5: 40}}
        N = Ns.get(L, {}).get(b, np.nan)
        sa = f'{a[0]:.3f}+-{a[1]:.3f}'
        sc = f'{c[0]:.3f}+-{c[1]:.3f}' if c else '   --   '
        print(f'  {b:5.1f} {L:2d} {N:5.0f} | {sa:>16s} {sc:>16s} | '
              f'{a[2].get(1, np.nan):12.3f} {a[2].get(2, np.nan):12.3f}')


if __name__ == '__main__':
    main()
