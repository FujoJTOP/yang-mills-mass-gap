"""
Combine independent chains per (beta,L) and extract the SU(2) 0++ gap + FSS.

For each (beta,L) group:
  * concatenate the O[n] config arrays of all matching-chain files,
  * use a FIXED smearing level n_op (default 8) so the operator is identical
    across volumes -> clean finite-size comparison; also report the
    highest-smearing level as a systematic.
  * connected correlator, cosh m_eff(t) for t=2..L//2, 1-cosh fit, jackknife.
Then build the FSS table Delta(beta,L)+-err and fit Delta_L = Delta_inf + c/L.

usage: python tm_gap_combine.py
"""
import glob
import os
import numpy as np
import tm_gap_gevp as G


def load_concat(files, level):
    """concatenate O[level] arrays across files; return (O, n_total, L, beta)."""
    Os = []; n = 0; L = None; beta = None
    for f in files:
        d = dict(np.load(f, allow_pickle=True))
        if f'O{level}' not in d:
            continue
        Os.append(d[f'O{level}'])
        n += d[f'O{level}'].shape[0]
        L = int(d['L']); beta = float(d['beta'])
    if not Os:
        return None
    return np.concatenate(Os, axis=0), n, L, beta


def conn(O):
    n, L = O.shape
    X = O - O.mean()
    return np.array([np.mean(X * np.roll(X, -t, axis=1)) for t in range(L)])


def meff(C, L):
    out = {}
    for t in range(1, L // 2 + 1):
        cm, c0, cp = C[(t - 1) % L], C[t % L], C[(t + 1) % L]
        if min(c0, cm + cp) > 1e-14 and (cm + cp) / (2 * c0) >= 1.0:
            out[t] = float(np.arccosh((cm + cp) / (2 * c0)))
        else:
            out[t] = np.nan
    return out


def fit1cosh(C, L, tmin=2, tmax=None):
    if tmax is None:
        tmax = L // 2
    ts = np.arange(tmin, tmax + 1)
    if len(ts) < 2:
        return np.nan
    y = C[ts]
    ms = np.linspace(0.05, 6.0, 6000)
    Z = np.exp(-ms[:, None] * ts[None, :]) + np.exp(-ms[:, None] * (L - ts)[None, :])
    A = (Z @ y) / np.sum(Z * Z, axis=1)
    r = np.sum((y[None, :] - A[:, None] * Z) ** 2, axis=1)
    k = int(np.argmin(r))
    return float(ms[k]) if 0 < k < len(ms) - 1 else np.nan


def extract(files, level, tau_per_bin=2):
    r = load_concat(files, level)
    if r is None:
        return None
    O, n, L, beta = r
    # bin: use a modest bin (2*tau_int ~ 3-4); configs from distinct chains
    b = 4
    ts = list(range(2, L // 2 + 1))

    def fn(idx):
        C = conn(O[idx])
        me = meff(C, L)
        vals = [me[t] for t in ts] + [fit1cosh(C, L, tmin=2)]
        return np.array(vals)

    full, err = G.jack(fn, n, b)
    return dict(L=L, beta=beta, n=n, level=level, ts=ts,
                me=full[:-1], errs=err[:-1], fit=full[-1], fiterr=err[-1], nfiles=len(files))


def pick(rec, relmax=0.5):
    cand = [(t, rec['me'][i], rec['errs'][i]) for i, t in enumerate(rec['ts'])
            if np.isfinite(rec['me'][i]) and np.isfinite(rec['errs'][i]) and rec['errs'][i] > 0]
    good = [c for c in cand if c[2] / max(abs(c[1]), 1e-9) < relmax]
    pick_ = good[-1] if good else (min(cand, key=lambda x: x[2]) if cand else None)
    if pick_ is None:
        return None
    t, v, e = pick_
    return t, v, e, (e / max(abs(v), 1e-9) > 0.35)


def fss_fit(points):
    pts = [(L, v, e) for (L, v, e) in points if np.isfinite(v) and np.isfinite(e) and e > 0]
    if len(pts) < 2:
        return None
    Ls = np.array([p[0] for p in pts], float); vs = np.array([p[1] for p in pts], float)
    es = np.array([p[2] for p in pts], float)
    X = np.vstack([np.ones_like(Ls), 1.0 / Ls]).T; W = np.diag(1.0 / es ** 2)
    cov = np.linalg.inv(X.T @ W @ X); bb = cov @ X.T @ W @ vs
    r = vs - X @ bb; chi2 = float(r.T @ W @ r) / max(len(pts) - 2, 1)
    return bb[0], np.sqrt(cov[0, 0]), chi2


def groups_auto():
    """Build (beta,L)->file list from base + chains."""
    G_ = {}
    base = glob.glob('tm_gap_data/L*_b*.npz')
    for f in base:
        d = dict(np.load(f, allow_pickle=True))
        L = int(d['L']); beta = float(d['beta'])
        G_.setdefault((round(beta, 2), L), []).append(f)
    # decorated: L6_b2.2_deep belongs to (2.2,6)
    for f in glob.glob('tm_gap_data/chains/*.npz'):
        d = dict(np.load(f, allow_pickle=True))
        L = int(d['L']); beta = float(d['beta'])
        G_.setdefault((round(beta, 2), L), []).append(f)
    return G_


def main():
    G_ = groups_auto()
    level = 8
    print(f"==== Per-(beta,L) extraction at FIXED smearing level n={level} ====\n")
    table = {}
    for (beta, L) in sorted(G_):
        rec = extract(G_[(beta, L)], level)
        if rec is None:
            continue
        pl = pick(rec)
        tagf = f"b={beta} L={L}"
        me_s = ", ".join(f"t{t}:{rec['me'][i]:.3f}({rec['errs'][i]:.2f})"
                         for i, t in enumerate(rec['ts']) if np.isfinite(rec['me'][i]))
        print(f"[{tagf}] N={rec['n']} nfiles={rec['nfiles']}  m_eff: {me_s}")
        print(f"            1cosh fit: {rec['fit']:.3f} +- {rec['fiterr']:.3f}")
        if pl:
            t, v, e, nz = pl
            print(f"            Delta = {v:.3f} +- {e:.3f} @ t={t}" + ("  [NOISY]" if nz else ""))
            table.setdefault(beta, {})[L] = (v, e, rec['fit'], rec['fiterr'], t, nz)
        print()

    allL = sorted({L for bb in table for L in table[bb]})
    print("########## FSS TABLE  Delta(beta,L)  [op level n=%d] ##########" % level)
    print("beta  | " + " | ".join(f" L={L}          " for L in allL))
    for beta in sorted(table):
        row = table[beta]; cells = []
        for L in allL:
            if L in row:
                v, e, fv, fe, t, nz = row[L]
                cells.append(f"{v:6.3f}+-{e:5.3f}" + ("*" if nz else " "))
            else:
                cells.append("      ---       ")
        print(f"{beta:>5} | " + " | ".join(cells))
    print("\n(* noise-dominated)\n")
    print("---- FSS fit Delta_L = Delta_inf + c/L (weighted) ----")
    for beta in sorted(table):
        pts = [(L, v, e) for L, (v, e, fv, fe, t, nz) in table[beta].items()]
        fit = fss_fit(pts)
        if fit:
            print(f"  beta={beta}: Delta_inf = {fit[0]:.3f} +- {fit[1]:.3f}  (chi2/dof={fit[2]:.2f}, n={len(pts)})")
        else:
            print(f"  beta={beta}: <2 points")


if __name__ == '__main__':
    main()
