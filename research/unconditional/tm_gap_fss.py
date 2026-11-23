"""
Finite-size-scaling consolidation for the SU(2) 0++ transfer-matrix gap.

Primary extractor (stable at low statistics):
  * the highest-smeared operator O_n (best ground-state overlap),
  * connected correlator C(t) = <O(t')O(t'+t)>_c, periodic time,
  * cosh effective mass m_eff(t) = arccosh[(C(t-1)+C(t+1))/2C(t)], t>=2
    (t=1 is excluded: C(0) carries the contact term and biases m_eff(1)),
  * 1-state cosh fit C(t)=A[e^{-mt}+e^{-m(L-t)}] over t=2..L//2.

Cross-check: GEVP lambda(t) with reference t0=2 over ALL smearing levels,
same cosh effective mass for t>=3.

Jackknife errors, bin = ceil(2*tau_int(plaquette)).

FSS: Delta_L = Delta_inf + c/L  (errors-in-weights) per beta.

Run:  python tm_gap_fss.py [files...]
"""
import sys
import numpy as np
import tm_gap_gevp as G


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
    if len(ts) < 1:
        return np.nan
    y = C[ts]
    ms = np.linspace(0.05, 6.0, 6000)
    Z = np.exp(-ms[:, None] * ts[None, :]) + np.exp(-ms[:, None] * (L - ts)[None, :])
    A = (Z @ y) / np.sum(Z * Z, axis=1)
    r = np.sum((y[None, :] - A[:, None] * Z) ** 2, axis=1)
    k = int(np.argmin(r))
    return float(ms[k]) if 0 < k < len(ms) - 1 else np.nan


def analyse(f, t0_gevp=2):
    d = dict(np.load(f, allow_pickle=True))
    L = int(d['L']); beta = float(d['beta']); levels = sorted(int(x) for x in d['levels'])
    n = d['plaq'].shape[0]
    b = max(1, int(np.ceil(2 * G.tau_int(d['plaq']))))
    hi = levels[-1]

    # pointwise m_eff (t>=2) and cosh fit for the highest-smeared operator
    def fn_hi(idx):
        C = conn(d[f'O{hi}'][idx])
        me = meff(C, L)
        ts = list(range(2, L // 2 + 1))
        vals = [me[t] for t in ts] + [fit1cosh(C, L, tmin=2)]
        return np.array(vals)
    full, err = G.jack(fn_hi, n, b)
    ts = list(range(2, L // 2 + 1))
    me_vals = full[:-1]; me_errs = err[:-1]
    fit_v, fit_e = full[-1], err[-1]

    # GEVP cross-check, reference t0=2, evaluate t>=3
    Olist = [d[f'O{lv}'] for lv in levels]
    def fn_gev(idx):
        lam = G.gevp_lambda([O[idx] for O in Olist], L, t0_gevp)
        de = G.delta_eff(lam, L)
        ts2 = list(range(3, L // 2 + 1))
        return np.array([de[t] for t in ts2])
    gts = list(range(3, L // 2 + 1))
    if gts and len(levels) >= 2:
        gfull, gerr = G.jack(fn_gev, n, b)
    else:
        gfull, gerr = np.array([]), np.array([])

    return dict(f=f, L=L, beta=beta, n=n, bin=b, hi=hi, levels=levels,
                ts=ts, me_vals=me_vals, me_errs=me_errs, fit_v=fit_v, fit_e=fit_e,
                gts=gts, gvals=gfull, gerrs=gerr)


def plateau(rec, relmax=0.5):
    """Delta at the largest t (>=2) with relative error < relmax; else best."""
    cand = [(t, rec['me_vals'][i], rec['me_errs'][i])
            for i, t in enumerate(rec['ts'])
            if np.isfinite(rec['me_vals'][i]) and np.isfinite(rec['me_errs'][i])
            and rec['me_errs'][i] > 0]
    good = [c for c in cand if c[2] / max(abs(c[1]), 1e-9) < relmax]
    pick = (good[-1] if good else (min(cand, key=lambda x: x[2]) if cand else None))
    if pick is None:
        return None
    t, v, e = pick
    noisy = (e / max(abs(v), 1e-9) > 0.35)
    return (t, v, e, noisy)


def fss_fit(points):
    pts = [(L, v, e) for (L, v, e) in points
           if np.isfinite(v) and np.isfinite(e) and e > 0]
    if len(pts) < 2:
        return None
    Ls = np.array([p[0] for p in pts], float)
    vs = np.array([p[1] for p in pts], float)
    es = np.array([p[2] for p in pts], float)
    X = np.vstack([np.ones_like(Ls), 1.0 / Ls]).T
    W = np.diag(1.0 / es ** 2)
    cov = np.linalg.inv(X.T @ W @ X)
    bb = cov @ X.T @ W @ vs
    r = vs - X @ bb
    chi2 = float(r.T @ W @ r) / max(len(pts) - 2, 1)
    return bb[0], np.sqrt(cov[0, 0]), chi2


def main():
    files = sys.argv[1:]
    if not files:
        import glob
        files = sorted(glob.glob('tm_gap_data/L*_b*.npz'))
    files = [f for f in files if '_smoke' not in f]
    table = {}
    for f in files:
        rec = analyse(f)
        tag = f.split('\\')[-1].split('/')[-1].replace('.npz', '')
        pl = plateau(rec)
        print(f"\n=== {tag}  L={rec['L']} beta={rec['beta']} N={rec['n']} "
              f"bin={rec['bin']} op n={rec['hi']}  levels={rec['levels']}")
        for i, t in enumerate(rec['ts']):
            v, e = rec['me_vals'][i], rec['me_errs'][i]
            if np.isfinite(v):
                fl = "  (noise)" if (e / max(abs(v), 1e-9) > 0.5) else ""
                print(f"   m_eff(t={t}) = {v:7.3f} +- {e:6.3f}{fl}")
        print(f"   1cosh fit (t>=2..L/2): {rec['fit_v']:.3f} +- {rec['fit_e']:.3f}")
        if len(rec['gts']):
            gs = ", ".join(f"t={t}:{rec['gvals'][i]:.3f}({rec['gerrs'][i]:.2f})"
                           for i, t in enumerate(rec['gts']) if np.isfinite(rec['gvals'][i]))
            print(f"   GEVP t0=2 cross-check: {gs}")
        if pl:
            t, v, e, noisy = pl
            print(f"   --> Delta(L={rec['L']},b={rec['beta']}) = {v:.3f} +- {e:.3f} at t={t}"
                  + ("  [NOISY]" if noisy else ""))
            table.setdefault(rec['beta'], {})[rec['L']] = (v, e, rec['fit_v'], rec['fit_e'], t, noisy)

    print("\n\n########## FSS SUMMARY (primary: highest-smeared cosh m_eff) ##########")
    allL = sorted({L for bb in table for L in table[bb]})
    hdr = "beta  | " + " | ".join(f"  L={L}        " for L in allL)
    print(hdr)
    for beta in sorted(table):
        row = table[beta]; cells = []
        for L in allL:
            if L in row:
                v, e, fv, fe, t, nz = row[L]
                cells.append(f"{v:6.3f}+-{e:5.3f}" + ("*" if nz else " "))
            else:
                cells.append("      ---       ")
        print(f"{beta:>5} | " + " | ".join(cells))
    print("\n(* = noise-dominated; Delta quoted at largest usable t>=2)\n")
    print("---- FSS: Delta_L = Delta_inf + c/L  (weighted) ----")
    for beta in sorted(table):
        pts = [(L, v, e) for L, (v, e, fv, fe, t, nz) in table[beta].items()]
        fit = fss_fit(pts)
        if fit:
            A, eA, chi = fit
            print(f"  beta={beta}:  Delta_inf = {A:.3f} +- {eA:.3f}   (chi2/dof={chi:.2f}, n={len(pts)})")
        else:
            print(f"  beta={beta}:  <2 finite points, no fit")


if __name__ == '__main__':
    main()
