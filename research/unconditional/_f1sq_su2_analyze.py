"""Analyse SU(2) connected tr F^2 correlators -> residue f1^2.

For each APE-smearing level s and each config k we formed S_s[k,t] (zero-
momentum time-slice sum).  The connected correlator is
    C_{ss'}(t) = <S_s(t) S_{s'}(0)> - <S_s><S_{s'}>
                = mean_k P^{ss'}_k(t) - mean_k mu^s_k * mean_k mu^{s'}_k
with P^{ss'}_k(t) = (1/T) sum_{t0} S_s[k,t0+t] S_{s'}[k,t0], mu^s_k = mean_t S_s[k].

f1^2 = |<0|O_s|G>|^2 is the amplitude A of the lightest cosh:
    C(t) = A ( e^{-m t} + e^{-m (T-t)} )   for  t in [tmin, T/2].
Errors: block bootstrap over configs (block size from plaquette autocorr).
"""
import sys
import numpy as np
from scipy.linalg import eigh as sp_eigh


_MGRID = np.linspace(0.02, 6.0, 6000)


def fit_residue(C1d, T, tmin, tmax, mgrid=_MGRID):
    """Robust 1-state cosh fit, A profiled out: for each trial m,
       A(m) = (w.y)/(w.w), w=exp(-m t)+exp(-m(T-t)); pick m minimising r(m).
       Returns (A, m).  Stability: grid search, no optimizer."""
    ts = np.arange(tmin, tmax + 1)
    y = np.asarray(C1d, float)[ts]
    W = np.exp(-np.outer(mgrid, ts)) + np.exp(-np.outer(mgrid, T - ts))
    denom = np.sum(W * W, axis=1)
    A = (W @ y) / denom
    r = np.sum((y[None, :] - A[:, None] * W) ** 2, axis=1)
    k = int(np.argmin(r))
    return float(A[k]), float(mgrid[k])


def products(d, levels):
    """mu (N, nlev) and P[(i,j)] = (N,T) config products."""
    N = d[f'S{levels[0]}'].shape[0]
    T = d[f'S{levels[0]}'].shape[1]
    mu = np.stack([d[f'S{s}'].mean(axis=1) for s in levels], axis=1)   # (N,nlev)
    P = {}
    for a, sa in enumerate(levels):
        Sa = d[f'S{sa}']
        for b, sb in enumerate(levels):
            Sb = d[f'S{sb}']
            P[(a, b)] = np.array([np.mean(Sa * np.roll(Sb, -t, axis=1), axis=1)
                                  for t in range(T)]).T                 # (N,T)
    return mu, P, T


def connected(idx, mu, P, T):
    """C matrix (nlev,nlev,T) for config subset idx."""
    nlev = mu.shape[1]
    m = mu[idx].mean(0)
    C = np.zeros((nlev, nlev, T))
    for (a, b), Pab in P.items():
        C[a, b] = Pab[idx].mean(0) - m[a] * m[b]
    return C


def ac_time(x):
    x = x - x.mean()
    v = np.dot(x, x) / x.size
    if v <= 0:
        return 1.0
    s = 1.0
    for lag in range(1, min(200, x.size // 4)):
        c = np.dot(x[:-lag], x[lag:]) / (x.size - lag) / v
        if c <= 0:
            break
        s += 2.0 * c
    return max(s, 1.0)


def block_boot_indices(N, nboot, block, rng):
    nb = int(np.ceil(N / block))
    for _ in range(nboot):
        starts = rng.integers(0, N, nb)
        idx = np.concatenate([(np.arange(s, s + block) % N) for s in starts])[:N]
        yield idx


def eff_mass(C1d, T):
    out = []
    for t in range(1, T - 1):
        a, b, c = C1d[t - 1], C1d[t], C1d[t + 1]
        if min(a, b, c) > 0:
            r = (a + c) / (2.0 * b)
            out.append(float(np.arccosh(r)) if r >= 1.0 else np.nan)
        else:
            out.append(np.nan)
    return np.array(out)


def gevp_mass(Cmat, T, t0=1, tmin=2, use=None):
    """GEVP C(t) v = lam C(t0) v  (generalized symmetric).  Ground (largest)
       lam(t) ~ e^{-m(t-t0)}; m = -ln(lam(t)/lam(t-1)).  use = level indices."""
    if use is None:
        use = list(range(Cmat.shape[0]))
    use = list(use)
    sub = Cmat[np.ix_(use, use)]            # (k,k,T)
    C0 = sub[:, :, t0]
    out = []
    lam_prev = None
    for t in range(t0, T // 2 + 1):
        Ct = sub[:, :, t]
        try:
            w = sp_eigh(Ct, C0, eigvals_only=True)
        except Exception:
            out.append(np.nan); lam_prev = None; continue
        lam = float(np.max(np.real(w)))
        if lam_prev is not None and lam > 0 and lam_prev > 0:
            out.append(float(np.log(lam_prev / lam)))
        else:
            out.append(np.nan)
        lam_prev = lam
    return np.array(out)


def analyse(path, nboot=500, verbose=True):
    d = np.load(path)
    levels = [int(x) for x in d['levels']]
    L = int(d['L']); beta = float(d['beta'])
    mu, P, T = products(d, levels)
    N = mu.shape[0]
    plaq = d['plaq'] if 'plaq' in d else np.array([0.0])
    tau = ac_time(np.asarray(plaq, float)) if plaq.size > 5 else 1.0
    block = int(max(1, round(2.0 * tau)))
    rng = np.random.default_rng(12345)

    allidx = np.arange(N)
    C = connected(allidx, mu, P, T)

    tmin = 1
    tmax = T // 2

    # central per-level fits + bootstrap errors
    res = {}
    for a, s in enumerate(levels):
        A0, m0 = fit_residue(C[a, a], T, tmin, tmax)
        me = eff_mass(C[a, a], T)
        # spectral-ratio fraction C(1)/C(0)
        c1c0 = C[a, a][1] / C[a, a][0] if C[a, a][0] != 0 else np.nan
        Ab = []; mb = []; c1b = []
        for idx in block_boot_indices(N, nboot, block, rng):
            Cb = connected(idx, mu, P, T)
            A_, m_ = fit_residue(Cb[a, a], T, tmin, tmax)
            Ab.append(A_); mb.append(m_)
            c1b.append(Cb[a, a][1] / Cb[a, a][0] if Cb[a, a][0] != 0 else np.nan)
        Ab = np.array(Ab); mb = np.array(mb); c1b = np.array(c1b)
        res[s] = dict(A=A0, A_se=float(np.nanstd(Ab)), m=m0, m_se=float(np.nanstd(mb)),
                      me=me, c1c0=c1c0, c1c0_se=float(np.nanstd(c1b)),
                      Acv=float(np.nanstd(Ab) / max(abs(A0), 1e-300)) if A0 else np.inf)

    # GEVP mass over smeared levels only (level>=2); unsmeared is contact noise
    guse = [a for a, s in enumerate(levels) if s >= 2]
    if len(guse) >= 2:
        gv = gevp_mass(C, T, use=guse)
    else:
        gv = np.array([])

    if verbose:
        print(f'# {path}')
        print(f'# SU(2) L={L} beta={beta} N={N} T={T} block={block} '
              f'(tau_plaq={tau:.1f}) tfit={tmin}..{tmax}')
        print(f'  {"n":>3s} {"A=f1^2":>12s} {"+-se":>10s} {"S/N":>7s} {"m":>7s} {"+-":>6s} '
              f'{"m_eff(1..)":>22s} {"C1/C0":>9s} {"+-":>7s}')
        for s in levels:
            r = res[s]
            me_s = ' '.join('%.2f' % x if np.isfinite(x) else '  - ' for x in r['me'][:4])
            sn = r['A'] / r['A_se'] if r['A_se'] > 0 else np.nan
            print(f'  {s:3d} {r["A"]:+12.5g} {r["A_se"]:10.3g} {sn:+7.1f} {r["m"]:7.3f} '
                  f'{r["m_se"]:6.3f} {me_s:>22s} {r["c1c0"]:+9.4f} {r["c1c0_se"]:7.4f}')
        print(f'  GEVP m_eff (t=2..T/2): ' +
              ' '.join('%.3f' % x if np.isfinite(x) else '-' for x in gv))
    return dict(levels=levels, L=L, beta=beta, N=N, T=T, block=block, res=res, gevp=gv)


if __name__ == '__main__':
    for fn in sys.argv[1:]:
        analyse(fn)
        print()
