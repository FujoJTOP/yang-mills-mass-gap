"""MEMORY-SAFE streaming SU(2) 0++ tr F^2 correlator for the Clay residue f1^2.

Single-process, no multiprocessing.  Samples configs one at a time and stores
ONLY zero-momentum time-slice sums S_n[k,t] (shape N x L per smearing level),
never the 4*L^4 gauge fields -> RSS stays at tens of MB, independent of N.

O_s(t) = sum_x (1/3) sum_{1<=i<j<=3} (1/2) Re tr U_ij(x,t),  APE-smeared
spatial links at level s in LEVELS.  f1^2 = amplitude A of the lightest cosh:
    C(t) = A ( e^{-m t} + e^{-m (T-t)} ).
Reuses lattice-cpb/w4_blocked.py (SU(2) heat bath, staple, random_su2) and
su3.project_su3 (polar projection works for 2x2).
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
import w4_blocked as W   # noqa: E402
import su3 as S          # noqa: E402

adj = W.adj


def mag_density_su2(U):
    """s(t,x,y,z) = (1/3) sum_{spatial i<j} (1/2) Re tr U_ij.  axis0 = time."""
    L = U.shape[1]
    s = np.zeros((L, L, L, L))
    for i in (1, 2, 3):
        for j in range(i + 1, 4):
            Up = (U[i] @ np.roll(U[j], -1, axis=i)
                  @ adj(np.roll(U[i], -1, axis=j)) @ adj(U[j]))
            s += 0.5 * np.real(np.trace(Up, axis1=-2, axis2=-1))
    return s / 3.0


def ape_spatial2(U, alpha=0.5, n=1):
    """APE-smear the three spatial links (mu=1,2,3); temporal fixed."""
    U = U.copy()
    for _ in range(n):
        old = U.copy()
        for i in (1, 2, 3):
            Ss = np.zeros_like(U[i])
            for j in (1, 2, 3):
                if j == i:
                    continue
                Ss += old[j] @ np.roll(old[i], -1, axis=j) @ adj(np.roll(old[j], -1, axis=i))
                Ss += (adj(np.roll(old[j], 1, axis=j)) @ np.roll(old[i], 1, axis=j)
                       @ np.roll(np.roll(old[j], 1, axis=j), -1, axis=i))
            U[i] = S.project_su3((1.0 - alpha) * old[i] + alpha / 6.0 * Ss)
    return U


def run(L, beta, N, therm, gap, seed, alpha, levels, out, pilot=False):
    rng = np.random.default_rng(seed)
    kappa = beta / 2.0
    # random start (better decorrelation at strong coupling than cold)
    U = W.random_su2(rng, (4, L, L, L, L))

    t0 = time.time()
    for _ in range(therm):
        Uf = U.reshape(4, -1, 2, 2)
        L2 = L
        ns = L ** 4
        x0, x1, x2, x3 = np.meshgrid(*[np.arange(L)] * 4, indexing='ij')
        par = (x0 + x1 + x2 + x3) % 2
        idxp = [np.nonzero(par.reshape(-1) == c)[0] for c in (0, 1)]
        for c in (0, 1):
            ids = idxp[c]
            for mu in range(4):
                St = W.staple_sum2(U, mu).reshape(ns, 2, 2)[ids]
                A = Uf[mu, ids]
                M = A @ St
                v = S.su2_heatbath(kappa, M, 0, 1, rng)
                Uf[mu, ids] = v @ A
    t1 = time.time()
    print(f'[L={L} b={beta}] therm {therm} in {t1-t0:.1f}s '
          f'({(t1-t0)/max(therm,1):.3f}s/sweep)', flush=True)

    out_d = {n: [] for n in levels}
    plaq = []
    t1 = time.time()
    for k in range(N):
        for _ in range(gap):
            Uf = U.reshape(4, -1, 2, 2)
            ns = L ** 4
            x0, x1, x2, x3 = np.meshgrid(*[np.arange(L)] * 4, indexing='ij')
            par = (x0 + x1 + x2 + x3) % 2
            idxp = [np.nonzero(par.reshape(-1) == c)[0] for c in (0, 1)]
            for c in (0, 1):
                ids = idxp[c]
                for mu in range(4):
                    St = W.staple_sum2(U, mu).reshape(ns, 2, 2)[ids]
                    A = Uf[mu, ids]
                    M = A @ St
                    v = S.su2_heatbath(kappa, M, 0, 1, rng)
                    Uf[mu, ids] = v @ A
        plaq.append(float(np.mean(mag_density_su2(U))))
        prev = 0
        Us = U
        for n in levels:
            if n > prev:
                Us = ape_spatial2(Us, alpha, n - prev)
                prev = n
            out_d[n].append(mag_density_su2(Us).sum(axis=(1, 2, 3)))
        if (k + 1) % max(1, N // 10) == 0:
            print(f'  [{k+1}/{N}] {time.time()-t1:.0f}s', flush=True)
        if pilot and k + 1 >= 4:
            break

    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    np.savez(out, L=L, beta=beta, alpha=alpha, levels=np.array(levels),
             plaq=np.array(plaq), N=len(plaq), gap=gap, therm=therm, seed=seed,
             **{f'S{n}': np.array(v) for n, v in out_d.items()})
    print('saved', out, f'total {time.time()-t0:.1f}s  N={len(plaq)}', flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=6)
    ap.add_argument('--beta', type=float, required=True)
    ap.add_argument('--N', type=int, default=2000)
    ap.add_argument('--therm', type=int, default=1000)
    ap.add_argument('--gap', type=int, default=2)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--alpha', type=float, default=0.5)
    ap.add_argument('--levels', type=str, default='0,2,4,8')
    ap.add_argument('--out', type=str, required=True)
    ap.add_argument('--pilot', action='store_true')
    args = ap.parse_args()
    levels = tuple(int(x) for x in args.levels.split(','))
    run(args.L, args.beta, args.N, args.therm, args.gap, args.seed,
        args.alpha, levels, args.out, pilot=args.pilot)


if __name__ == '__main__':
    main()
