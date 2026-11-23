"""Single-level inverse Monte Carlo (Swendsen) for the b=2 blocked SU(3) theory.

Model:  P_lam(V) = exp( sum_i lam_i s_i(V) ) / Z,
        s_i = mean over the loop family of (1/3) Re tr  (per-loop density),
        families in `ops` subset of {P,R,S,R3}.

Newton/Swendsen iteration:
        d<s>/dlam = +Cov_model(s)        (sign: S = sum lam (1-T))
        lam <- lam + damp * C_model^{-1} (s_true - s_model)

The true blocked ensemble is the one produced by gen_blocked.py.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              '..', 'lattice-cpb'))

import su3 as S
import model as M
import blk


# ---------------------------------------------------------------- helpers
def jackknife(fn, n, b):
    b = max(b, 1)
    nb = n // b
    if nb < 2:
        return fn(np.arange(n)), 0.0
    bs = [np.arange(i * b, (i + 1) * b) for i in range(nb)]
    full = fn(np.arange(n))
    vals = np.array([fn(np.concatenate([bs[j] for j in range(nb) if j != i]))
                     for i in range(nb)])
    err = np.sqrt((nb - 1) / nb * np.sum((vals - vals.mean(axis=0)) ** 2, axis=0))
    return full, err


def tau_int(x):
    x = np.asarray(x, float)
    n = len(x)
    x = x - x.mean()
    v = np.dot(x, x) / n
    if v <= 0:
        return 0.5
    s = 0.0
    for t in range(1, min(n // 10, 500)):
        s += np.dot(x[:-t], x[t:]) / (n - t) / v
        if t >= 6 * (0.5 + s):
            break
    return max(0.5 + s, 0.5)


FAM_PER_SITE = {'P': 6, 'R': 12, 'S': 6, 'R3': 12,
                'C6': 6, 'C8': 6, 'C10': 6}
DATA_KEY = {'P': 'op_P', 'R': 'op_R', 'S': 'op_S', 'R3': 'op_R3',
            'C6': 'c6', 'C8': 'c8', 'C10': 'c10'}


def family_counts(L, ops):
    """Number of loops of each family on the blocked lattice (L^4 sites)."""
    return np.array([FAM_PER_SITE[k] * L ** 4 for k in ops], float)


def op_matrix(d, ops):
    """Per-configuration operator sums S_i = sum over loops of (1/3) Re tr.

    Using the *sums* (not per-loop means) makes d<S_i>/d lam_j = Cov(S_i,S_j)
    exactly, which is the standard Swendsen form.
    """
    L = int(d['L']) // 2
    nf = family_counts(L, ops)
    X = np.array([d[DATA_KEY[k]] for k in ops]) * nf[:, None]
    return X


def model_run(lam, kinds, L, nsweep, therm, meas, seed, Uf0=None,
              verbose=False, full=False, width=1.0):
    """Sample the model at couplings lam (len == len(kinds)); return per-config
    observables and statistics."""
    tab = M.build_table(L, kinds)
    w = np.asarray(lam, float)
    args = (tab['t_first'], tab['t_nt'], tab['t_pos'], tab['t_kind'],
            tab['t_code'], tab['slot_link'], tab['slot_conj'], tab['link_off'],
            tab['link_term'], tab['order_idx'], tab['order_cnt'])
    if Uf0 is None:
        rng = np.random.default_rng(seed)
        U = S.random_su3(rng, (4, L, L, L, L))
        Uf = np.transpose(U, (1, 2, 3, 4, 0, 5, 6)).reshape(4 * L ** 4, 3, 3).copy()
    else:
        Uf = Uf0.copy()
    np.random.seed(seed)
    nacc = np.zeros(2, np.int64)
    for _ in range(therm):
        M.sweep(Uf, w, *args, nacc, width)
    orbits = blk.orbits_for(L)
    recs = []
    karr = np.arange(len(kinds))
    for k in range(nsweep):
        M.sweep(Uf, w, *args, nacc, width)
        if (k + 1) % meas:
            continue
        if full:
            U2 = np.transpose(Uf.reshape(L, L, L, L, 4, 3, 3),
                              (4, 0, 1, 2, 3, 5, 6))
            rec = blk.measure_all(U2, orbits)
        else:
            dens = M.loop_densities(Uf, w, tab['t_first'], tab['t_nt'],
                                    tab['t_pos'], tab['t_kind'], tab['t_code'],
                                    tab['slot_link'], tab['slot_conj'],
                                    tab['link_off'], tab['link_term'], karr)
            rec = {DATA_KEY[kd]: float(dens[i]) for i, kd in enumerate(kinds)}
        recs.append(rec)
    if verbose:
        print(f"    model lam={np.round(w,4)} acc={nacc[0]/max(1,nacc[1]):.2f}")
    return recs, acc(Uf, w, tab, kinds, nacc), Uf


def acc(Uf, w, tab, kinds, nacc):
    return nacc[0] / max(1, nacc[1])


# ---------------------------------------------------------------- fit
def fit(d, ops, kinds=None, iters=10, sweeps=2500, therm=400, meas=2,
        seed=1, damp=0.8, verbose=True, width=1.0):
    Lfine = int(d['L']); L = Lfine // 2; beta = float(d['beta'])
    kinds = list(ops) if kinds is None else list(kinds)
    X = op_matrix(d, ops)
    ncfg = X.shape[1]
    b = max(1, int(np.ceil(2 * tau_int(d['op_P']))))
    s_true = X.mean(axis=1)
    _, s_err = jackknife(lambda idx: X[:, idx].mean(axis=1), ncfg, b)
    Xb = X[:, :ncfg // b * b].reshape(len(ops), -1, b).mean(axis=2)
    C_true = np.atleast_2d(np.cov(Xb)) if Xb.shape[1] > 1 else np.eye(len(ops))
    ev = np.linalg.eigvalsh(C_true)
    hist = []
    lam = np.zeros(len(ops))
    Uf = None
    for it in range(iters):
        t0 = time.time()
        recs, a, Uf = model_run(lam, kinds, L, sweeps, therm, meas,
                                seed + 101 * it, Uf0=Uf, verbose=verbose,
                                width=width)
        nf = family_counts(L, ops)
        Ym = np.array([[r[DATA_KEY[k]] for r in recs] for k in ops]) * nf[:, None]
        s_mod = Ym.mean(axis=1)
        nmod = Ym.shape[1]
        bm = max(1, int(np.ceil(2 * tau_int(Ym[0]))))
        Ym_b = Ym[:, :nmod // bm * bm].reshape(len(ops), -1, bm).mean(axis=2)
        C_mod = np.atleast_2d(np.cov(Ym_b)) if Ym_b.shape[1] > 1 else np.eye(len(ops))
        r = s_true - s_mod
        step = np.linalg.solve(C_mod, r)
        nrm = np.linalg.norm(step)
        if nrm > 0.7:
            step *= 0.7 / nrm
        lam = lam + damp * step
        hist.append(dict(it=it, lam=lam.tolist(), s_mod=s_mod.tolist(),
                         residual=r.tolist(), step=step.tolist(), acc=a,
                         cond_mod=float(np.linalg.cond(C_mod)),
                         time=round(time.time() - t0, 1)))
        if verbose:
            msg = (f"  it{it} lam={np.round(lam,4)} r={np.round(r,4)} "
                   f"cond={np.linalg.cond(C_mod):.1f} {time.time()-t0:.0f}s")
            print(msg, flush=True)
            with open('progress.txt', 'a') as fh:
                fh.write(f"L{L} b{beta} {'/'.join(ops)} " + msg + "\n")
        if np.all(np.abs(r) < np.maximum(3 * s_err, 2e-4)):
            break
    out = dict(L=L, beta=beta, ops=list(ops), s_true=s_true.tolist(),
               s_err=s_err.tolist(), ncfg=ncfg, bin=b,
               C_true=C_true.tolist(), eig_true=ev.tolist(),
               cond_true=float(ev[-1] / ev[0]) if ev[0] > 0 else np.inf,
               lam=lam.tolist(), hist=hist)
    return out, Uf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--ops', default='P,R,S')
    ap.add_argument('--iters', type=int, default=10)
    ap.add_argument('--sweeps', type=int, default=2500)
    ap.add_argument('--meas', type=int, default=2)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    d = dict(np.load(args.data))
    ops = args.ops.split(',')
    out, _ = fit(d, ops, iters=args.iters, sweeps=args.sweeps,
                 meas=args.meas, seed=args.seed)
    json.dump(out, open(args.out, 'w'), indent=1)
    print('saved', args.out)


if __name__ == '__main__':
    main()
