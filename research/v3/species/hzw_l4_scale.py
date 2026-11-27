#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hzw_l4_scale.py -- HZW-L2 companion: the N_B trend of varpi at L=4, SU(2), beta=2.5.
Same action / same heat bath as hzw_l2_za.py. Blocks of b=1,2,4 lattice spacings
(N_B = 1, 16, 256). This is the decisive point for "O(1) vs ∝ N_B":
at fixed a, varpi(b) = N_B * Var(block average of q), q = 1 - (1/2) Re tr U_plaq.
"""
import numpy as np
import json, time, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hzw_l2_za as M


def run(L, beta, ntherm, nkeep, nstreams, seed0=9000):
    V = L ** 4
    c = np.array(np.unravel_index(np.arange(V), (L,) * 4)).T
    # block index for b=2 blocks: two sites per direction
    blk2 = ((c[:, 0] // 2) * 2 + (c[:, 1] // 2)) * 4 + (c[:, 2] // 2) * 2 + (c[:, 3] // 2)
    A = []
    out = []
    t0 = time.time()
    for s in range(nstreams):
        rng = np.random.default_rng(seed0 + s)
        U = M.random_su2(rng, (4, L, L, L, L))
        for _ in range(ntherm):
            M.sweep(U, L, beta, rng)
        rows = []
        for it in range(nkeep):
            M.sweep(U, L, beta, rng)
            q = M.q_planes(M.plaquettes(U))
            q1 = q[0].reshape(-1)
            q6 = q.sum(axis=0).reshape(-1)
            b2 = np.bincount(blk2, weights=q1, minlength=16) / (V // 16)
            b2_6 = np.bincount(blk2, weights=q6, minlength=16) / (V // 16)
            rows.append([q1.mean(), q6.mean(), q1.var(), q6.var(), *b2, *b2_6,
                         (1 - q).mean()])
        rows = np.array(rows)
        out.append(rows)
    return out, time.time() - t0


if __name__ == '__main__':
    L, beta = 4, 2.5
    ntherm, nkeep, nstreams = 500, 4000, 4
    out, dt = run(L, beta, ntherm, nkeep, nstreams)
    R = np.concatenate(out, axis=0)     # (ncfg, 2 + 2 + 16 + 16 + 1)
    n = len(R)
    qbar1, qbar6 = R[:, 0], R[:, 1]
    vq1, vq6 = R[:, 2], R[:, 3]
    B2_1 = R[:, 4:20]
    B2_6 = R[:, 20:36]
    Pmean = R[:, 36]
    V = L ** 4
    v1 = qbar1.var(ddof=1)
    v6 = qbar6.var(ddof=1)
    v2_1 = B2_1.var(axis=0, ddof=1).mean()
    v2_6 = B2_6.var(axis=0, ddof=1).mean()
    # errors: 8 blocks over configurations
    def be(x):
        nb = 8
        b = x[:len(x) // nb * nb].reshape(nb, -1)
        est = [b[i].var(ddof=1) for i in range(nb)]
        return np.std(est, ddof=1) / np.sqrt(nb)
    print("HZW L=4 beta=2.5 : n=%d configs, streams=%d x %d sweeps, wall %.1fs" % (n, nstreams, nkeep, dt))
    print("<P>            = %.5f +- %.5f      (repo single-config 0.6526)" % (Pmean.mean(), Pmean.std(ddof=1) / np.sqrt(n)))
    print("Var(q)  (b=1)  = %.6f +- %.6f" % (vq1.mean(), vq1.std(ddof=1) / np.sqrt(n)))
    print("Var(q6) (b=1)  = %.6f +- %.6f" % (vq6.mean(), vq6.std(ddof=1) / np.sqrt(n)))
    print("varpi(b=1)     = %.6f          (N_B=1)" % vq1.mean())
    print("varpi(b=2a)    = %.6f +- %.6f  (N_B=16)  = 16*Var(qbar_16)" % (16 * v2_1, 16 * be(B2_1.mean(axis=1)),))
    print("varpi(b=4a)    = %.6f +- %.6f  (N_B=256) = 256*Var(qbar_256)" % (V * v1, V * be(qbar1)))
    print("ratio varpi(4a)/varpi(2a) = %.3f   (N_B ratio = 16)" % (V * v1 / (16 * v2_1)))
    print("six-plane: varpi(2a)=%.5f  varpi(4a)=%.5f" % (16 * v2_6, V * v6))
    js = dict(L=L, beta=beta, ncfg=int(n), P=float(Pmean.mean()),
              varq1=float(vq1.mean()), varq6=float(vq6.mean()),
              varpi2a=float(16 * v2_1), varpi4a=float(V * v1),
              varpi2a_6=float(16 * v2_6), varpi4a_6=float(V * v6),
              e_varpi4a=float(V * be(qbar1)))
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hzw_l4_scale.json'), 'w') as f:
        json.dump(js, f, indent=1)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hzw_l4_scale.log'), 'w') as f:
        f.write(repr(js) + "\n")
