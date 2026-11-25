#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gevp_4d_fast.py -- PARALLEL-BV-CHAIN version of gevp_4d.py (lead constraint: 10 min wall).

* Sampler: the **Numba heat-bath kernels of `w7/m1_fast.py`** (`sweep`, `heatbath_link`,
  `staple`, `random_start`), imported, NOT rewritten.
* Parallelism: BY CHAIN.  `Pool(12)`, one independent heat-bath chain per worker, its own
  seed, its own `np.random.seed`; `n = N/12` configs per chain.  Each worker runs with
  NUMBA_NUM_THREADS=1 / OMP_NUM_THREADS=1 (set before numba is imported) so 12 processes do
  not fight over cores.
* Operators: UNCHANGED definitions, imported from `gevp_4d.py`
  (`M.plaq_operator_t`, `W.ape_all2(alpha=0.5)`), i.e. o0 = space-summed plaquette,
  o2 = same after APE(0.5,2), o4 = same after APE(0.5,4).  Zero momentum, A1++.
* Storage: per-config operator array (n, nops, L) float64 -> `gevp_4d_L{L}_ops.npy`
  (the correlator matrix and the GEVP are then built by `gevp_4d.analyze`, unchanged).

Usage: py -3.12 -I gevp_4d_fast.py L N [therm] [gap] [nchain] [mode]
       mode = ym (default) | haar (beta=0) | free (Gaussian free massless field)
"""
import os
os.environ['NUMBA_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import sys, time, json
import numpy as np

sys.path.insert(0, r'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, r'D:/Dev/math/research/unconditional/w7')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import m1_fast as F            # noqa: E402  numba heat-bath kernels (unchanged)
import gevp_4d as G            # noqa: E402  operators + GEVP analysis (unchanged)
from multiprocessing import Pool  # noqa: E402

NOPS = 3


def observe(U):
    """The A1++ zero-momentum operators of one configuration, (nops, L).  Same as gevp_4d."""
    out = [G.M.plaq_operator_t(U)]
    for s in G.STEPS:
        out.append(G.M.plaq_operator_t(G.W.ape_all2(U, G.ALPHA, s)))
    return np.array(out)


def chain(arg):
    """One independent heat-bath chain -> its share of the operator array."""
    L, beta, nconf, therm, gap, seed, mode = arg
    np.random.seed(seed)
    if mode == 'free':
        rng = np.random.default_rng(seed)
        ops = np.zeros((nconf, NOPS, L))
        for i in range(nconf):
            ops[i] = observe(G.free_su2_config(L, rng, g=0.1))
        return ops
    U = F.random_start(L, seed) if F.random_start is not None else F.cold_start(L)
    for _ in range(therm):
        F.sweep(U, beta, L, 0.0)
    ops = np.zeros((nconf, NOPS, L))
    for i in range(nconf):
        for _ in range(gap):
            F.sweep(U, beta, L, 0.0)
        ops[i] = observe(U)
    return ops


def main():
    a = sys.argv[1:]
    L = int(a[0]); N = int(a[1])
    therm = int(a[2]) if len(a) > 2 else 200
    gap = int(a[3]) if len(a) > 3 else 5
    nch = int(a[4]) if len(a) > 4 else 12
    mode = a[5] if len(a) > 5 else 'ym'
    beta = 0.0 if mode == 'haar' else 2.5
    nch = max(1, min(nch, N))
    base, rem = divmod(N, nch)
    args = []
    for k in range(nch):
        n = base + (1 if k < rem else 0)
        if n:
            args.append((L, beta, n, therm, gap, 1000 + 17 * k + 101 * L, mode))
    print(f'# gevp_4d_fast  L={L} N={N} mode={mode} beta={beta} therm={therm} gap={gap} '
          f'chains={len(args)}', flush=True)
    t0 = time.time()
    with Pool(len(args)) as p:
        parts = p.map(chain, args)
    ops = np.concatenate(parts, axis=0)
    dt = time.time() - t0
    np.save(os.path.join(HERE, f'gevp_4d_{"free_" if mode=="free" else ""}'
                               f'{"haar_" if mode=="haar" else ""}L{L}_ops.npy'), ops)
    print(f'# measured {ops.shape} in {dt:.1f}s  ({ops.shape[0]/dt:.1f} conf/s, '
          f'{dt/max(1,ops.shape[0])*1000:.2f} ms/conf)', flush=True)
    plq = np.array([1.0 - 0.0])   # plaquette not stored in the fast path
    res = G.analyze(L, ops=ops, plq=None,
                    meta=dict(therm=therm, gap=gap, nchain=len(args), mode=mode,
                              engine='gevp_4d_fast.py + m1_fast.py njit heat bath',
                              wall_s=dt, cps=ops.shape[0] / dt),
                    tag=f'{mode}')
    G.show(res)
    name = ('gevp_4d_free_L%d.json' % L) if mode == 'free' else \
           ('gevp_4d_haar_L%d.json' % L) if mode == 'haar' else ('gevp_4d_L%d.json' % L)
    print('# ->', G.save_json(res, name))


if __name__ == '__main__':
    main()
