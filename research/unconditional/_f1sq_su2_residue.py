"""Reproduce the COMPUTE-F1SQ-SU2 residue table from saved SU(2) correlators.

For each npz produced by _f1sq_su2_mem.py:
  * per-level effective mass  m_eff(t)=arccosh[(C(t-1)+C(t+1))/2C(t)]
  * plateau mass              m_G = mean over smeared levels (s>=2) of m_eff(t=2)
  * per-level connected C(1)/C(0) with 40-block jackknife error + S/N
  * residue f1^2 = A from the single-cosh fit over t=1,2 with m fixed to m_G,
    jackknifed; degeneracy: A = (w.C)/(w.w), w = e^{-m t}+e^{-m(T-t)}.

Usage:  python _f1sq_su2_residue.py f1sq_b2.2_L6.npz f1sq_b1.0_L6.npz
"""
import sys
import numpy as np
import _f1sq_su2_analyze as FA


def analyse(path):
    d = np.load(path)
    levels = [int(x) for x in d['levels']]
    L = int(d['L']); T = L; beta = float(d['beta'])
    mu, P, T = FA.products(d, levels)
    N = mu.shape[0]; allidx = np.arange(N)
    C_of = lambda idx: FA.connected(idx, mu, P, T)
    C = C_of(allidx)

    def mg_of(Cm):
        ms = []
        for a, s in enumerate(levels):
            if s < 2:
                continue
            c = Cm[a, a]; r = (c[1] + c[3]) / (2 * c[2])
            if r >= 1 and np.isfinite(r):
                ms.append(np.arccosh(r))
        return float(np.mean(ms)) if ms else np.nan

    mG = mg_of(C)
    nb = 40; bs = [np.arange(i * (N // nb), (i + 1) * (N // nb)) for i in range(nb)]
    idxs = [np.concatenate([bs[j] for j in range(nb) if j != i]) for i in range(nb)]
    mg = np.array([mg_of(C_of(idx)) for idx in idxs])
    emG = float(np.sqrt((nb - 1) / nb * np.sum((mg - mg.mean()) ** 2)))

    ts = np.array([1, 2])
    w_of = lambda m: np.exp(-m * ts) + np.exp(-m * (T - ts))
    w = w_of(mG)

    print(f'# {path}  beta={beta} L={L} N={N}  m_G={mG:.3f}+-{emG:.3f}')
    print(f'  {"n":>2s} {"m_eff(1,2,3)":>18s} {"C1/C0":>9s} {"S/N":>5s} '
          f'{"f1^2":>10s} {"+-err":>8s} {"S/N":>5s}')
    for a, s in enumerate(levels):
        c = C[a, a]
        me = []
        for t in (1, 2, 3):
            if t + 1 < T:
                r = (c[t - 1] + c[t + 1]) / (2 * c[t])
                me.append(np.arccosh(r) if (np.isfinite(r) and r >= 1) else np.nan)
            else:
                me.append(np.nan)
        # C1/C0 jackknife
        r0 = c[1] / c[0]
        rv = np.array([C_of(idx)[a, a][1] / C_of(idx)[a, a][0] for idx in idxs])
        e_r = float(np.sqrt((nb - 1) / nb * np.sum((rv - rv.mean()) ** 2)))
        # residue with m fixed to mG
        A = float((w @ c[ts]) / (w @ w))
        Av = np.array([float((w @ C_of(idx)[a, a][ts]) / (w @ w)) for idx in idxs])
        eA = float(np.sqrt((nb - 1) / nb * np.sum((Av - Av.mean()) ** 2)))
        mesa = ','.join('%.2f' % x if np.isfinite(x) else '-' for x in me)
        print(f'  {s:2d} {mesa:>18s} {r0:+9.4f} {abs(r0)/max(e_r,1e-30):5.1f} '
              f'{A:+10.4f} {eA:8.4f} {A/eA if eA > 0 else float("nan"):5.1f}')
    print()


if __name__ == '__main__':
    for p in sys.argv[1:]:
        analyse(p)
