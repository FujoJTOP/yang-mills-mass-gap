#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gevp_rho2.py -- rho^2(L) = (e2/e1)^2 from the zero-momentum 0++ GEVP of 4D SU(2) Wilson
gauge theory at beta=2.5, L=6 and L=8 (ONE heavy job at a time: run L=6, then L=8).

Correspondence [T] (pure algebra, same as gevp_4d.py:7 and table_for:195-197):
    C(t) v_n = lam_n(t,t0) C(t0) v_n           (symmetric whitening of C(t0))
    e_n(t,t0) = -ln lam_n / (t-t0)             <- the two lowest GEVP "energies"
    rho      = e2/e1                           <- this IS e2/e1, NOT its square
    rho^2    = (e2/e1)^2                       <- the KT object; lower bound rho^2 >= 1+c
So "e2/e1" is rho (a mass ratio, expected 1.56 for the 0++*/0++ pair), and rho^2 is that
number squared (expected 2.43).  Reporting e2/e1 and printing both.

Sampler: numba heat-bath kernels of m1_fast.py (imported, not rewritten), one independent
chain per worker process (nchain processes, NUMBA/OMP threads pinned to 1 each).
Operators: gevp_4d.observe (o0 = plaquette, o2 = APE(0.5,2), o4 = APE(0.5,4)) + R = the
spatial 2x1 rectangle (geometry change: breaks the o0/o2/o4 collinearity that made the
3-operator C(t0) rank deficient in gevp_4d_L6.json).
Errors: hierarchical bootstrap (chains resampled, then configs inside a chain).

Usage:  py -3.12 -I w7/gevp_rho2.py measure 6 8000 [nchain] [therm] [gap]
        py -3.12 -I w7/gevp_rho2.py analyze 6
Labels [T]/[C]; no QED.
"""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import sys, time, json                              # noqa: E402
import numpy as np                                  # noqa: E402

sys.path.insert(0, r'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, r'D:/Dev/math/research/v3/w7')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import m1_fast as F                                 # noqa: E402  njit heat bath
import gevp_4d as G                                 # noqa: E402  observe/C_from_w/gevp
from multiprocessing import Pool                    # noqa: E402

BETA = 2.5
TAGS = ['o0_plaq', 'o2_ape2', 'o4_ape4', 'R_2x1']


def rect_t(U):
    """O(t) = sum_{i<j spatial} sum_x (1 - 1/2 Re tr [2x1 retrace in the (i,j) plane])."""
    L = U.shape[1]
    adj = G.W.adj
    op = np.zeros(L)
    for i in (1, 2, 3):
        for j in range(i + 1, 4):
            P = (U[i] @ np.roll(U[i], -1, axis=i) @ np.roll(U[j], -2, axis=i)
                 @ np.roll(np.roll(adj(U[i]), -1, axis=j), -1, axis=i)
                 @ np.roll(adj(U[i]), -1, axis=j) @ adj(U[j]))
            op += (1.0 - np.real(np.trace(P, axis1=-2, axis2=-1)) / 2.0).sum(axis=(1, 2, 3))
    return op


def observe(U):
    return np.vstack([G.observe(U), rect_t(U)[None, :]])     # (4, L)


def chain(arg):
    L, beta, nconf, therm, gap, seed = arg
    np.random.seed(seed)
    U = F.random_start(L, seed)
    for _ in range(therm):
        F.sweep(U, beta, L, 0.0)
    ops = np.zeros((nconf, len(TAGS), L))
    for i in range(nconf):
        for _ in range(gap):
            F.sweep(U, beta, L, 0.0)
        ops[i] = observe(U)
    return ops


def measure(L, N, nchain=8, therm=400, gap=10):
    base, rem = divmod(N, nchain)
    args = [(L, BETA, base + (1 if k < rem else 0), therm, gap, 7000 + 17 * k + 101 * L)
            for k in range(nchain)]
    t0 = time.time()
    with Pool(nchain) as p:
        parts = p.map(chain, args)
    ops = np.concatenate(parts, axis=0)
    wall = time.time() - t0
    np.save(os.path.join(HERE, f'gevp_L{L}_ops.npy'), ops)
    json.dump(dict(L=L, N=int(ops.shape[0]), bounds=[int(p.shape[0]) for p in parts],
                   beta=BETA, tags=TAGS, therm=therm, gap=gap, nchain=nchain,
                   wall_s=wall, engine='m1_fast njit heat bath, chain-parallel'),
              open(os.path.join(HERE, f'gevp_L{L}_meta.json'), 'w'), indent=1)
    print(f'# measured L={L} ops={ops.shape} in {wall:.1f}s '
          f'({ops.shape[0]/wall:.1f} conf/s, {1000*wall/ops.shape[0]:.2f} ms/conf, '
          f'{nchain} chains)', flush=True)
    return ops


SUBSETS = {'o0o2o4': (0, 1, 2), 'o0o2R': (0, 1, 3), 'o2o4R': (1, 2, 3),
           'all4': (0, 1, 2, 3), 'o0R': (0, 3), 'o2R': (1, 3)}


def readout(Cm, L, idx, t0, t):
    sub = Cm[np.ix_(idx, idx, np.arange(L))]
    lam, rank, cond = G.gevp(sub[:, :, t], sub[:, :, t0])
    if lam is None or len(lam) < 2 or lam[0] <= 0 or lam[1] <= 0:
        return None
    e1 = -np.log(lam[0]) / (t - t0)
    e2 = -np.log(lam[1]) / (t - t0)
    if not np.isfinite(e1 + e2) or e1 <= 0:
        return None
    rho = e2 / e1
    return dict(e1=float(e1), e2=float(e2), rho=float(rho), rho2=float(rho * rho),
                rank=int(rank), cond=float(cond))


def boot_counts(bounds, rng):
    st = np.concatenate([[0], np.cumsum(bounds)])
    c = np.zeros(int(st[-1]))
    for _ in range(len(bounds)):
        k = int(rng.integers(len(bounds)))
        a, b = int(st[k]), int(st[k + 1])
        c[rng.integers(a, b, b - a)] += 1
    return c.astype(float)


def analyze(L, boots=300):
    ops = np.load(os.path.join(HERE, f'gevp_L{L}_ops.npy'))
    meta = json.load(open(os.path.join(HERE, f'gevp_L{L}_meta.json')))
    bounds = meta['bounds']
    N, no, Lz = ops.shape
    assert Lz == L and no == len(TAGS)
    S, abar = G.per_config_outer(ops)
    C0 = G.C_from_w(np.ones(N), S, abar, L)

    acc = {}
    rng = np.random.default_rng(20260101)
    for _ in range(boots):
        Cb = G.C_from_w(boot_counts(bounds, rng), S, abar, L)
        for lab, idx in SUBSETS.items():
            for t0 in (1, 2):
                for t in range(t0 + 1, L // 2 + 1):
                    r = readout(Cb, L, idx, t0, t)
                    if r is not None:
                        acc.setdefault((lab, t0, t), []).append(r)

    out = dict(L=L, N=int(N), beta=BETA, tags=TAGS, bounds=bounds,
               meta=meta, boot=dict(kind='hierarchical (chains, then configs)', B=boots,
                                    seed=20260101),
               readout='rho=e2/e1=-ln(lam2)/-ln(lam1) per DOF; rho^2=(e2/e1)^2',
               rows=[], diag_meff=[])
    for i in range(no):
        out['diag_meff'].append([float(np.log(C0[i, i, t] / C0[i, i, t + 1]))
                                 if C0[i, i, t] > 0 and C0[i, i, t + 1] > 0 else None
                                 for t in range(1, L // 2 + 1)])
    print(f'# L={L} N={N} beta={BETA}  diag m_eff (log) per operator, t=1..{L//2}')
    for i, tg in enumerate(TAGS):
        print(f'  {tg:9s} ' + ' '.join(
            'nan' if v is None else f'{v:6.3f}' for v in out['diag_meff'][i]))
    print(f"  {'subset':8s} {'t0':>2} {'t':>2} {'rk':>2} {'e1':>7} {'e2':>7} "
          f"{'e2/e1':>7} {'+-':>7} {'rho2':>7} {'+-':>7}")
    for lab, idx in SUBSETS.items():
        for t0 in (1, 2):
            for t in range(t0 + 1, L // 2 + 1):
                r = readout(C0, L, idx, t0, t)
                row = dict(subset=lab, idx=list(idx), t0=t0, t=t, ok=r is not None)
                if r is None:
                    out['rows'].append(row)
                    print(f'  {lab:8s} {t0:2d} {t:2d}   -- singular / non-positive --')
                    continue
                row.update(r)
                arr = acc.get((lab, t0, t), [])
                if len(arr) >= 20:
                    r2 = np.array([a['rho2'] for a in arr], float)
                    rr = np.array([a['rho'] for a in arr], float)
                    r2 = r2[np.isfinite(r2)]
                    rr = rr[np.isfinite(rr)]
                    row['rho2_se'] = float(r2.std(ddof=1))
                    row['rho_se'] = float(rr.std(ddof=1))
                    row['rho2_q16'], row['rho2_q84'] = (float(np.percentile(r2, 16)),
                                                        float(np.percentile(r2, 84)))
                    row['n_boot'] = int(r2.size)
                out['rows'].append(row)
                print(f"  {lab:8s} {t0:2d} {t:2d} {r['rank']:2d} {r['e1']:7.4f} "
                      f"{r['e2']:7.4f} {r['rho']:7.4f} {row.get('rho_se', float('nan')):7.4f} "
                      f"{r['rho2']:7.3f} {row.get('rho2_se', float('nan')):7.3f}")
    json.dump(out, open(os.path.join(HERE, f'gevp_L{L}.json'), 'w'), indent=1)
    print(f'# -> w7/gevp_L{L}.json', flush=True)


if __name__ == '__main__':
    a = sys.argv[1:]
    cmd = a[0] if a else 'help'
    if cmd == 'measure':
        L = int(a[1]); N = int(a[2])
        measure(L, N, nchain=int(a[3]) if len(a) > 3 else 8,
                therm=int(a[4]) if len(a) > 4 else 400,
                gap=int(a[5]) if len(a) > 5 else 10)
    elif cmd == 'analyze':
        analyze(int(a[1]), boots=int(a[2]) if len(a) > 2 else 300)
    else:
        print(__doc__)
