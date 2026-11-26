"""U14 — does w1 = f1^2/M0 depend on L?  (the one question of route 14)

Data, in order of trust:
  1. beta=2.2, L=6, N=4000 archive  research/unconditional/f1sq_b2.2_L6.npz
     (levels 0,2,4,8; rebuilds the p8 numbers 0.4838 / 0.8382 / 0.8569 / 0.7929)
  2. beta=2.2, L=4 (fresh, N as given) and L=8 (fresh, N as given): SAME operator,
     same sampler (w4_blocked.sample2), same APE levels -> the L-triple.
  3. beta=2.5, L=8: the existing 40-config pilot research/lattice-cpb/pilot_b2.5.npz
     (O_0, O_8 are already per-config zero-momentum time sums) -> w1 with a huge error.
  4. beta=2.5, L=4 fresh (cheap) -- but T=4 at beta=2.5 is the DECONFINED side
     (Nt=4 critical beta_c ~ 2.2986), flagged, not mixed into the confined trend.

Estimator (identical for every point):
  C(t) = <S(t)S(0)>_c ,  M0 = C(0) = Var(S) ,  A = amplitude of A*cosh(m(T/2-t)),
  w1 = A/M0 .  Also the finite-t estimator  w1hat(t) = e^{m t} C(t)/C(0) >= w1 (T2 sec 2.3),
  which needs no fit and is what a "no plateau" run really gives.
Errors: block bootstrap over configurations.

Usage: py -3.12 research/unconditional/w7/u14_w1L.py
"""
import os
import sys
import time

import numpy as np

CB = 'D:/Dev/math/research/lattice-cpb'
UN = 'D:/Dev/math/research/unconditional'
sys.path.insert(0, CB)
sys.path.insert(0, UN)
import w4_blocked as W            # noqa: E402  (SU(2) heat bath, sample2)
import _f1sq_su2_run as F         # noqa: E402  (mag_density_su2, ape_spatial2)
import _f1sq_su2_analyze as A     # noqa: E402  (fit_residue, block_boot_indices)

OUT = 'D:/Dev/math/research/unconditional/w7'


def measure(cfgs, levels, alpha=0.5):
    """per-config zero-momentum time sums S{n}[k,t] for the APE ladder."""
    out = {n: [] for n in levels}
    for U in cfgs:
        prev, Us = 0, U
        for n in levels:
            if n > prev:
                Us = F.ape_spatial2(Us, alpha, n - prev)
                prev = n
            out[n].append(F.mag_density_su2(Us).sum(axis=(1, 2, 3)))
    return {n: np.array(v) for n, v in out.items()}


def conn(S):
    """connected zero-momentum correlator C(t) from S[k,t]."""
    T = S.shape[1]
    mu = S.mean(axis=1)
    P = np.array([np.mean(S * np.roll(S, -t, axis=1), axis=1) for t in range(T)]).T
    m = mu.mean()
    return P.mean(axis=0) - m * m


def w1_stats(S, nboot=400, seed=7):
    """w1 = A/M0 with a cosh fit, plus w1hat(t).  Block bootstrap."""
    N, T = S.shape
    C = conn(S)
    M0 = C[0]
    tmin = 1 if T <= 4 else 2
    tmax = T // 2
    A_, m_ = A.fit_residue(C, T, tmin, tmax)
    w1 = A_ / M0 if M0 != 0 else np.nan
    # finite-t estimator with the fitted m (decreasing upper bound on w1)
    ts = list(range(1, tmax + 1))
    what = {t: float(np.exp(m_ * t) * C[t] / M0) for t in ts}

    rng = np.random.default_rng(seed)
    block = max(1, N // 20)
    bw, bm = [], []
    for idx in A.block_boot_indices(N, nboot, block, rng):
        Cb = conn(S[idx])
        if Cb[0] == 0:
            continue
        Ab, mb = A.fit_residue(Cb, T, tmin, tmax)
        bw.append(Ab / Cb[0]); bm.append(mb)
    return dict(T=T, N=N, M0=float(M0), C1C0=float(C[1] / M0),
                A=float(A_), m=float(m_), w1=float(w1),
                w1_se=float(np.nanstd(bw)), m_se=float(np.nanstd(bm)),
                what=what, tmin=tmin, tmax=tmax)


def fresh(L, beta, N, therm=600, gap=2, levels=(0, 2, 4, 8), seed=11, tag=''):
    rng = np.random.default_rng(seed)
    t0 = time.time()
    cfgs = W.sample2(L, beta, rng, therm, N, gap)
    print(f'  [{tag}] sampled L={L} beta={beta} N={N} in {time.time()-t0:.0f}s', flush=True)
    S = measure(cfgs, levels)
    np.savez(f'{OUT}/u14_{tag}_L{L}_b{beta}.npz', L=L, beta=beta, N=N,
             levels=np.array(levels), **{f'S{n}': v for n, v in S.items()})
    return S


def report(beta, L, S, tag):
    print(f'\n=== beta={beta}  L={L}  ({tag}) ===')
    print(f'  {"lvl":>3s} {"N":>5s} {"T":>2s} {"M0=C(0)":>10s} {"C1/C0":>7s} '
          f'{"A=f1^2":>9s} {"m":>5s} {"w1=A/M0":>14s} {"w1hat(t=1..)":>22s}')
    rows = {}
    for n in sorted(S):
        r = w1_stats(S[n])
        rows[n] = r
        wh = ' '.join('%.2f' % r['what'][t] for t in sorted(r['what']))
        print(f'  {n:3d} {r["N"]:5d} {r["T"]:2d} {r["M0"]:10.4g} {r["C1C0"]:7.4f} '
              f'{r["A"]:9.4g} {r["m"]:5.2f} {r["w1"]:8.3f}+-{r["w1_se"]:.3f} {wh:>22s}',
              flush=True)
    return rows


def main():
    t00 = time.time()
    res = {}

    # (1) archived beta=2.2 L=6, N=4000 -- reproduce p8 first
    d = np.load(f'{UN}/f1sq_b2.2_L6.npz')
    lv = [int(x) for x in d['levels']]
    S6 = {n: d[f'S{n}'] for n in lv}
    print(f'[archive] beta=2.2 L=6 N={d["N"]} levels={lv}')
    res[(2.2, 6)] = report(2.2, 6, S6, 'archive N=4000')

    # (2) fresh L=4 and L=8 at beta=2.2, same operator/ladder
    S4 = fresh(4, 2.2, 2000, levels=(0, 2, 4, 8), seed=41, tag='b2.2')
    res[(2.2, 4)] = report(2.2, 4, S4, 'fresh N=2000')
    S8 = fresh(8, 2.2, 250, levels=(0, 2, 8), seed=81, tag='b2.2')
    res[(2.2, 8)] = report(2.2, 8, S8, 'fresh N=250')

    # (3) the 40-config beta=2.5 L=8 pilot (already measured, O_0 / O_8)
    p = np.load(f'{CB}/pilot_b2.5.npz')
    Sp = {int(s): p[f'O_{int(s)}'] for s in p['smears']}
    print(f"\n[pilot] beta=2.5 L={int(p['L'])} N={int(p['ncfg'])}")
    res[(2.5, 8)] = report(2.5, 8, Sp, 'pilot N=40')

    # (4) fresh L=4 at beta=2.5 (cheap) -- deconfined side, flagged
    S4b = fresh(4, 2.5, 1500, levels=(0, 2), seed=45, tag='b2.5')
    res[(2.5, 4)] = report(2.5, 4, S4b, 'fresh N=1500 (deconfined side)')

    # ---- summary: the L-trend of w1 at fixed beta, level s=0 ----
    print('\n================ U14 SUMMARY: w1(L) at level s=0 ================')
    print(f'  {"beta":>5s} {"L":>2s} {"N":>5s} {"M0":>9s} {"w1":>7s} {"+-":>6s} '
          f'{"w1hat(1)":>8s} {"w1hat(2)":>8s}')
    for (b, L) in sorted(res):
        r = res[(b, L)].get(0)
        if r is None:
            continue
        w1h = [r['what'].get(t, np.nan) for t in (1, 2)]
        print(f'  {b:5.1f} {L:2d} {r["N"]:5d} {r["M0"]:9.4g} {r["w1"]:7.3f} '
              f'{r["w1_se"]:6.3f} {w1h[0]:8.3f} {w1h[1]:8.3f}')

    print('\n--- M0 (total source variance) L-scaling at beta=2.2, s=0 ---')
    for L in (4, 6, 8):
        r = res[(2.2, L)].get(0)
        if r:
            print(f'  L={L}: M0={r["M0"]:.4g}  f1^2={r["A"]:.4g}  m_G={r["m"]:.3f}')
    print(f'\ntotal {time.time()-t00:.0f}s')


if __name__ == '__main__':
    main()
