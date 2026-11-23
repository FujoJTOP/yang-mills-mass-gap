"""
GEVP extraction of the 0++ transfer-matrix gap Delta from tm_gap_run.py data.

Operators O_n, n in {0,4,8} (APE smearing levels).  Build the connected
cross-correlator matrix  C_ab(t) = <O_a(t') O_b(t'+t)>_c  (periodic time),
solve  C(t) w = lambda(t) C(t0) w  for the leading eigenvalue lambda_0(t).

Because the time direction is PERIODIC (length L), every correlator carries
both e^{-Delta t} and e^{-Delta (L-t)} images, so a single-exponential
ln-ratio is wrong.  We therefore apply the standard cosh effective mass to
the leading GEVP eigenvalue (which has the same cosh shape as the scalar
correlator):
     Delta_eff(t) = arccosh( [lambda(t-1)+lambda(t+1)] / [2 lambda(t)] ),
evaluated at t=L//2 (the point where the periodic image is handled exactly).
The basis=[8] single-operator case reproduces the plain cosh m_eff; the
multi-operator GEVP is the variationally best 0++ interpolator.

Jackknife over configs (bins capped at 200).
"""
import glob
import sys

import numpy as np
import scipy.linalg as sla


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


def Cmat(Olist, t):
    """connected cross-correlator matrix C_ab(t), shape (n,n)."""
    n = len(Olist)
    means = [O.mean() for O in Olist]
    M = np.zeros((n, n))
    for a in range(n):
        Xa = Olist[a] - means[a]
        for b in range(n):
            Xb = Olist[b] - means[b]
            M[a, b] = np.mean(Xa * np.roll(Xb, -t, axis=1))
    return 0.5 * (M + M.T)


def gevp_lambda(Olist, L, t0=1):
    """leading GEVP eigenvalue lambda(t) vs reference C(t0), t=0..L-1."""
    C0 = Cmat(Olist, t0)
    lam = {}
    for t in range(L):
        Ct = Cmat(Olist, t)
        try:
            w, _ = sla.eigh(Ct, C0)
            lam[t] = float(np.max(w.real))
        except Exception:
            lam[t] = np.nan
    return lam


def delta_eff(lam, L):
    """cosh effective mass on lambda at t=1..L//2."""
    out = {}
    for t in range(1, L // 2 + 1):
        a = lam.get(t - 1, np.nan)
        b = lam.get(t, np.nan)
        c = lam.get((t + 1) % L, np.nan)
        if np.isfinite(a) and np.isfinite(b) and np.isfinite(c) and b != 0:
            arg = (a + c) / (2.0 * b)
            out[t] = float(np.arccosh(arg)) if arg >= 1.0 else np.nan
        else:
            out[t] = np.nan
    return out


def jack(fn, n, b, nbin_max=200):
    b = max(b, 1, n // nbin_max)
    nb = n // b
    if nb < 2:
        v = np.atleast_1d(fn(np.arange(n)))
        return v, np.full_like(v, np.nan)
    bs = [np.arange(i * b, min((i + 1) * b, n)) for i in range(nb)]
    full = np.atleast_1d(fn(np.arange(n)))
    vals = np.array([np.atleast_1d(fn(np.concatenate([bs[j] for j in range(nb) if j != i])))
                     for i in range(nb)])
    err = np.sqrt((nb - 1) / nb * np.sum((vals - vals.mean(axis=0)) ** 2, axis=0))
    return full, err


def analyse(f, t0=1, basis=None, verbose=True):
    d = dict(np.load(f, allow_pickle=True))
    L = int(d['L']); beta = float(d['beta']); levels = [int(x) for x in d['levels']]
    if basis is not None:
        levels = [lv for lv in levels if lv in basis]
    Olist = [d[f'O{lv}'] for lv in levels]
    n = Olist[0].shape[0]
    b = max(1, int(np.ceil(2 * tau_int(d['plaq']))))
    ts = list(range(1, L // 2 + 1))

    def fn(idx):
        lam = gevp_lambda([O[idx] for O in Olist], L, t0)
        de = delta_eff(lam, L)
        return np.array([de[t] for t in ts])

    full, err = jack(fn, n, b)
    tag = f.split('\\')[-1].split('/')[-1]
    if verbose:
        print(f"\n{tag}  L={L} beta={beta} N={n} bin={b} basis={levels} t0={t0}")
        for i, t in enumerate(ts):
            print(f"   GEVP Delta_eff(t={t}) = {full[i]:.3f} +- {err[i]:.3f}")
    # cleanest point = t=L//2 (last)
    return dict(L=L, beta=beta, t=ts[-1], val=float(full[-1]), err=float(err[-1]),
                ts=ts, full=full, errs=err, levels=levels)


if __name__ == '__main__':
    files = sys.argv[1:] or sorted(glob.glob('tm_gap_data/L*.npz'))
    for f in files:
        for basis in ([0, 4, 8], [4, 8], [8]):
            analyse(f, t0=1, basis=basis)
