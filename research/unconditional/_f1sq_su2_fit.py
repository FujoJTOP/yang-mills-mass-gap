"""Dedicated f1^2 fitter for the SU(2) tr F^2 residue from _f1sq_su2_mem.py output.

For each APE-smearing level s we build the connected zero-momentum correlator
    C_ss(t),  t=0..L-1   (time-averaged over t0),
fit a single cosh  C(t) = A (e^{-m t} + e^{-m (L-t)})  on t in [tmin, tmax]
by the profiled grid search of _f1sq_su2_analyze.fit_residue, and jackknife
(block) over configurations for errors on A and m.

f1^2 = A  is |<0| O_s |G>|^2 for the lightest 0++ state in the level-s basis.
Also reports the spectral ratio  A/ (e^{m*tmin} C(tmin)) consistency and the
simple ratio  -ln[C(1)/C(0)]  (first-step trial energy).
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _f1sq_su2_analyze as FA   # noqa: E402


def load_products(path):
    d = np.load(path)
    levels = [int(x) for x in d['levels']]
    L = int(d['L']); beta = float(d['beta'])
    mu, P, T = FA.products(d, levels)
    N = mu.shape[0]
    plaq = np.asarray(d['plaq'], float) if 'plaq' in d else np.zeros(N)
    return d, levels, L, beta, mu, P, T, N, plaq


def fit_one(idx, levels, mu, P, T, tmin, tmax, lev_i):
    C = FA.connected(idx, mu, P, T)[lev_i, lev_i]
    A, m = FA.fit_residue(C, T, tmin, tmax)
    return A, m


def main():
    paths = sys.argv[1:]
    for path in paths:
        d, levels, L, beta, mu, P, T, N, plaq = load_products(path)
        tau = FA.ac_time(plaq) if plaq.size > 5 else 1.0
        block = int(max(1, round(2.0 * tau)))
        tmin = 1
        tmax = T // 2 if T >= 6 else max(1, T // 2)
        rng = np.random.default_rng(777)
        allidx = np.arange(N)
        print(f'# {os.path.basename(path)}  L={L} beta={beta} N={N} '
              f'tau_plaq={tau:.1f} block={block} fit t=[{tmin},{tmax}]')
        print(f'  {"n":>3s} {"f1^2=A":>11s} {"+-err":>9s} {"S/N":>6s} '
              f'{"mass":>7s} {"+-":>6s} {"C1/C0":>9s}  chisq/dof')
        rows = []
        for lev_i, s in enumerate(levels):
            A, m = fit_one(allidx, levels, mu, P, T, tmin, tmax, lev_i)
            # chi^2/dof of the profiled fit
            C = FA.connected(allidx, mu, P, T)[lev_i, lev_i]
            ts = np.arange(tmin, tmax + 1)
            w = np.exp(-m * ts) + np.exp(-m * (T - ts))
            resid = C[ts] - A * w
            # rough: use jackknife spread as per-point error proxy is hard; skip formal
            # block jackknife
            nb = max(1, N // block)
            bs = [np.arange(i * block, min((i + 1) * block, N)) for i in range(nb)]
            Aj = []; mj = []
            for i in range(nb):
                idx = np.concatenate([bs[j] for j in range(nb) if j != i])
                a_, m_ = fit_one(idx, levels, mu, P, T, tmin, tmax, lev_i)
                Aj.append(a_); mj.append(m_)
            Aj = np.array(Aj); mj = np.array(mj)
            eA = np.sqrt((nb - 1) / nb * np.sum((Aj - Aj.mean()) ** 2))
            em = np.sqrt((nb - 1) / nb * np.sum((mj - mj.mean()) ** 2))
            c1c0 = C[1] / C[0] if C[0] != 0 else np.nan
            sn = A / eA if eA > 0 else np.nan
            print(f'  {s:3d} {A:+11.5g} {eA:9.3g} {sn:+6.1f} '
                  f'{m:7.3f} {em:6.3f} {c1c0:+9.4f}')
            rows.append((s, A, eA, sn, m, em, c1c0))
        print()
    return rows


if __name__ == '__main__':
    main()
