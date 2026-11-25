#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wg_l32_parallel.py -- L=32 gauge-invariant WG estimator in <=20 min via independent chains.

Why parallel chains: one L=32 heat-bath sweep is ~6 s; a single chain in 20 min gives only ~200
sweeps (N~40), too noisy to separate chi_inf~0.50 from a power law ~0.65.  Six independent chains
of N=50, therm=50, gap=2 (150 sweeps each) run concurrently and combine to N=300 with SE~0.04.

Estimator (same as WG): APE-smeared tr F^2 density O_s; Obar=(1/V)sum O_s;
chi_O = V * Var_c(Obar).  Gauge-invariant, total-normalized, non-self-averaging.

Run: py -3.12 wg_l32_parallel.py [nchains] [N] [therm] [gap]
Defaults: 6 50 50 2.  ONE process per chain; 6 processes total (memory ~2 GB).
Labels [V]/[O]; no QED.  Clay 0%.
"""
import sys, os, json, time
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
CPB = 'D:/Dev/math/research/lattice-cpb'
W7U = 'D:/Dev/math/research/unconditional/w7'
for p in (HERE, CPB, W7U):
    if p not in sys.path:
        sys.path.insert(0, p)

import w4_blocked as W
import m1_fast as F
import wg_invariant as G

L = 32
BETA = 2.5
NS_APE = 2


def chain(args):
    seed, N, therm, gap = args
    np.random.seed(seed)
    U = F.random_start(L, seed)
    F.sweep(U, BETA, L)  # JIT warm-up
    for _ in range(therm):
        F.sweep(U, BETA, L)
    obar = np.empty(N)
    q2 = np.empty(N)
    t0 = time.time()
    for i in range(N):
        for _ in range(gap):
            F.sweep(U, BETA, L)
        Us = W.ape_all2(U, alpha=0.5, nsteps=NS_APE)
        Os = G.scalar_density(Us)
        obar[i] = Os.mean()
        q2[i] = (Os * Os).mean()
    return dict(seed=seed, obar=obar.tolist(), q2=q2.tolist(),
                wall_s=time.time() - t0, N=N, therm=therm, gap=gap)


def main():
    nch = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    therm = int(sys.argv[3]) if len(sys.argv) > 3 else 50
    gap = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    t0 = time.time()
    args = [(3200 + 17 * i, N, therm, gap) for i in range(nch)]
    with Pool(nch) as pool:
        res = pool.map(chain, args)
    V = L ** 4
    obar = np.concatenate([np.array(r['obar']) for r in res])
    q2 = np.concatenate([np.array(r['q2']) for r in res])
    chiO = V * obar.var(ddof=1)
    se = V * obar.std(ddof=1) / np.sqrt(len(obar))
    C0 = float(q2.mean() - obar.mean() ** 2)
    out = dict(L=L, beta=BETA, ns_ape=NS_APE, nchains=nch, N_per_chain=N,
               therm=therm, gap=gap, N_total=len(obar),
               chi_O=float(chiO), chi_O_se=float(se), contact_C0=C0,
               chi_phys=float(chiO - C0),
               obar_mean=float(obar.mean()), obar_sd=float(obar.std(ddof=1)),
               chain_walls=[r['wall_s'] for r in res],
               wall_s=time.time() - t0,
               obar=obar.tolist())
    print(json.dumps({k: v for k, v in out.items() if k != 'obar'}, indent=1))
    with open(os.path.join(HERE, 'wg_l32_parallel_results.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
