"""Analyse trial-state energies R(psi) from _trial_state_run.py output.

For each APE-smearing level s (a trial operator) and for the GEVP-optimal
variational combination of levels, report:
  * connected zero-momentum correlator C(t) (jackknife error)
  * trial energy at the first step  R1 = -ln[C(1)/C(0)]
  * cosh effective mass  m_eff(t) = arccosh[(C(t-1)+C(t+1)) / 2C(t)]
  * GEVP: largest generalized eigenvalue of C(t=1) vs C(t=0); R_GEVP = -ln lambda1
"""
import glob
import os
import sys

import numpy as np
import scipy.linalg as sla

HERE = os.path.dirname(os.path.abspath(__file__))


def conn(Sm):
    """connected C(t) from time-slice sums Sm (nconf, L)."""
    n, L = Sm.shape
    X = Sm - Sm.mean()
    return np.array([np.mean(X * np.roll(X, -t, axis=1)) for t in range(L)])


def cross(Am, Bm):
    """C_AB(t) = <dA(0) dB(t)> from time-slice sums (nconf, L)."""
    n, L = Am.shape
    X = Am - Am.mean()
    Y = Bm - Bm.mean()
    return np.array([np.mean(X * np.roll(Y, -t, axis=1)) for t in range(L)])


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
    vals = np.array([fn(np.concatenate([bs[j] for j in range(nb) if j != i]))
                     for i in range(nb)])
    err = np.sqrt((nb - 1) / nb * np.sum((vals - vals.mean(axis=0)) ** 2, axis=0))
    return full, err


def meff_cosh(C, L):
    out = []
    for t in range(1, L // 2 + 1):
        cm, c0, cp = C[(t - 1) % L], C[t % L], C[(t + 1) % L]
        if min(cm, c0, cp) > 1e-16:
            arg = (cm + cp) / (2 * c0)
            out.append((t, float(np.arccosh(max(arg, 1.0)))))
        else:
            out.append((t, np.nan))
    return out


def gevp_energy(levels, data, idx):
    """largest GEVP eigenvalue at t=1 -> R = -ln lambda1, plus C(0),C(1)."""
    L = data['S%d' % levels[0]].shape[1]
    Xs = []
    for n in levels:
        v = data['S%d' % n][idx]
        Xs.append(v)
    m = len(levels)
    C0 = np.zeros((m, m))
    C1 = np.zeros((m, m))
    for a in range(m):
        for b in range(m):
            cab = cross(Xs[a], Xs[b])
            C0[a, b] = cab[0]
            C1[a, b] = cab[1 % L]
    # symmetrize
    C0 = 0.5 * (C0 + C0.T); C1 = 0.5 * (C1 + C1.T)
    try:
        w = sla.eigh(C1, C0, eigvals_only=True)
        lam = float(w[-1])
    except Exception:
        return np.nan, np.nan, np.nan
    return -np.log(max(lam, 1e-300)), C0[0, 0], C1[0, 0]


def main():
    lines = []
    for f in sorted(glob.glob(os.path.join(HERE, 'trial_runs', '*.npz'))):
        d = dict(np.load(f, allow_pickle=True))
        L = int(d['L']); beta = float(d['beta'])
        levels = [int(x) for x in d['levels']]
        n = len(d['plaq'])
        b = max(1, int(np.ceil(2 * tau_int(d['plaq']))))
        lines.append(f"### L={L} beta={beta:.3f}  N={n} bin={b} levels={levels}")
        # per-level trial energies
        for lev in levels:
            Sm = d['S%d' % lev]
            C, eC = jk(lambda i: conn(Sm[i]), n, b)
            C0 = C[0]
            R1 = -np.log(max(C[1] / C0, 1e-300)) if C[1] > 0 else np.nan
            R1j, eR1 = jk(lambda i: np.array([-np.log(max(conn(Sm[i])[1] / conn(Sm[i])[0], 1e-300))
                                              if conn(Sm[i])[1] > 0 else np.nan]), n, b)
            me = meff_cosh(C, L)
            ratio = C[1] / C0 if C0 != 0 else np.nan
            eratio = eC[1] / abs(C0) if C0 != 0 else np.nan
            lines.append(
                f"  n={lev:2d}: C1/C0={ratio:+.4f}+-{eratio:.4f} (sig {abs(ratio)/max(eratio,1e-30):.1f}s) "
                f"R1={R1:.3f}+-{float(eR1[0]):.3f}  meff={[(t, round(v,2)) for t,v in me]}")
        # GEVP over well-separated levels
        try:
            g = [lev for lev in (0, 4, 8) if lev in levels]
            if len(g) >= 2:
                Rg, C0g, C1g = gevp_energy(g, d, np.arange(n))
                Rgj, eRg = jk(lambda i: np.array([gevp_energy(g, d, i)[0]]), n, b)
                lines.append(f"  GEVP{g}: C1/C0={C1g/C0g:+.4f}  R_GEVP={Rg:.3f}+-{float(eRg[0]):.3f}")
        except Exception as ex:
            lines.append(f"  GEVP failed: {ex}")
        lines.append("")
    txt = '\n'.join(lines)
    print(txt)
    open(os.path.join(HERE, 'trial_runs', 'summary.txt'), 'w').write(txt + '\n')


if __name__ == '__main__':
    main()
