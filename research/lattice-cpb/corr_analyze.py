"""Analyse the APE-smeared connected plaquette correlators (corr/)."""
import glob
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def conn(Sm):
    """connected correlator C(t) from time-slice sums, shape (n,L)."""
    n, L = Sm.shape
    X = Sm - Sm.mean()
    return np.array([np.mean(X * np.roll(X, -t, axis=1)) for t in range(L)])


def fit_xi(C, L, tmin=1, tmax=None):
    if tmax is None:
        tmax = L // 2
    ts = np.arange(tmin, tmax + 1)
    y = C[ts]
    ms = np.linspace(0.05, 4.0, 4000)
    Z = np.exp(-ms[:, None] * ts[None, :]) + np.exp(-ms[:, None] * (L - ts)[None, :])
    A = (Z @ y) / np.sum(Z * Z, axis=1)
    r = np.sum((y[None, :] - A[:, None] * Z) ** 2, axis=1)
    k = int(np.argmin(r))
    if k == 0 or k == len(ms) - 1:
        return float('nan'), float('nan')   # fit not constrained
    return float(ms[k]), float(A[k])


def meff(C, L):
    out = []
    for t in range(1, L // 2 + 1):
        cm, c0, cp = C[(t - 1) % L], C[t % L], C[(t + 1) % L]
        if min(cm, c0, cp) > 1e-14:
            out.append((t, float(np.arccosh((cm + cp) / (2 * c0)))))
        else:
            out.append((t, np.nan))
    return out


def tau_int(x):
    x = np.asarray(x, float); n = len(x); x = x - x.mean()
    v = np.dot(x, x) / n
    if v <= 0:
        return 0.5
    s = 0.0
    for t in range(1, min(n // 10, 500)):
        s += np.dot(x[:-t], x[t:]) / (n - t) / v
        if t >= 6 * (0.5 + s):
            break
    return max(0.5 + s, 0.5)


def jk(fn, n, b):
    b = max(b, int(np.ceil(n / 400.0)))
    nb = n // b
    bs = [np.arange(i * b, (i + 1) * b) for i in range(nb)]
    full = fn(np.arange(n))
    vals = [fn(np.concatenate([bs[j] for j in range(nb) if j != i])) for i in range(nb)]
    vals = np.array(vals)
    err = np.sqrt((nb - 1) / nb * np.sum((vals - vals.mean(axis=0)) ** 2, axis=0))
    return full, err


def main():
    lines = []
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for f in sorted(glob.glob(os.path.join(HERE, 'corr', '*.npz'))):
        d = dict(np.load(f, allow_pickle=True))
        L = int(d['L']); beta = float(d['beta'])
        b = max(1, int(np.ceil(2 * tau_int(d['plaq']))))
        for lev in [int(x) for x in d['levels']]:
            S = d[f'S{lev}']
            n = len(S)
            C, eC = jk(lambda idx: conn(S[idx]), n, b)
            mfit, _ = fit_xi(C, L)
            mfitj, emfit = jk(lambda idx: np.array([fit_xi(conn(S[idx]), L)[0]]), n, b)
            emfit = float(emfit[0])
            me = meff(C, L)
            m1 = me[0][1] if me else np.nan
            m1j, em1 = jk(lambda idx: np.array([meff(conn(S[idx]), L)[0][1]]), n, b)
            em1 = float(em1[0])
            lines.append(f"L={L} beta={beta:.3f} n_smear={lev}: N={n} bin={b}")
            lines.append(f"    C/C0 = {np.round(C / C[0], 5)}  +- {np.round(eC / C[0], 5)}")
            lines.append(f"    fit m = {mfit:.3f} +- {emfit:.3f} (window t=1..{L // 2}), "
                         f"xi/a = {1 / mfit:.3f} +- {emfit / mfit ** 2:.3f}")
            lines.append(f"    am_eff(1) = {m1:.3f} +- {em1:.3f}; "
                         f"m_eff(t) = {[(t, round(v, 3)) for t, v in me]}")
            axs[0].plot(np.arange(L), C / C[0], 'o-', label=f'L{L} b{beta} n{lev}')
            axs[1].plot([t for t, v in me], [v for t, v in me], 'o-',
                        label=f'L{L} b{beta} n{lev}')
    axs[0].set_yscale('log'); axs[0].set_xlabel('t')
    axs[0].set_ylabel('C(t)/C(0)'); axs[0].legend(fontsize=6); axs[0].grid(alpha=.3)
    axs[1].set_xlabel('t'); axs[1].set_ylabel('a m_eff'); axs[1].legend(fontsize=6)
    axs[1].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'plots', 'smeared_correlators.png'), dpi=140)
    txt = '\n'.join(lines)
    print(txt)
    open(os.path.join(HERE, 'corr_summary.txt'), 'w').write(txt + '\n')

    # append the smeared-correlator results to results.csv
    with open(os.path.join(HERE, 'results.csv'), 'a') as fh:
        for f in sorted(glob.glob(os.path.join(HERE, 'corr', '*.npz'))):
            d = dict(np.load(f, allow_pickle=True))
            L = int(d['L'])
            beta = float(d['beta'])
            n = len(d['plaq'])
            b = max(1, int(np.ceil(2 * tau_int(d['plaq']))))
            for lev in [int(x) for x in d['levels']]:
                S = d[f'S{lev}']
                mf, emf = jk(lambda idx: np.array([fit_xi(conn(S[idx]), L)[0]]), n, b)
                fh.write(f'xi_over_a_smeared{lev},{L},{beta:.4f},'
                         f'{1.0 / float(mf[0]):.8g},nan,{n},nan,{b}\n')
                fh.write(f'mass_a_smeared{lev},{L},{beta:.4f},'
                         f'{float(mf[0]):.8g},{float(emf[0]):.4g},{n},nan,{b}\n')


if __name__ == '__main__':
    main()
