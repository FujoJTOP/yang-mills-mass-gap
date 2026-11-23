"""
SU(2) pure-gauge string tension from Creutz ratios + 0++ glueball mass.

Convention (matches repo): S = beta * sum_p (1 - 1/2 Re tr U_p),  kappa = beta/2.
Heat bath: su2_heatbath (su3.py), sweep pattern from w4_blocked.sample2 (validated).

Wilson loops: time-space planes (mu=0, nu in {1,2,3}); R = spatial extent, T = temporal.
Creutz ratio  chi(R,T) = -ln[ W(R,T) W(R-1,T-1) / ( W(R,T-1) W(R-1,T) ) ].
Jackknife errors over configs (binned).

Glueball: 0++ operator O(t) = sum_{x-space} sum_{1<=i<j<=3} Re tr U_ij / 2  on
spatial-APE-smeared configs. Connected correlator, cosh fit.

usage:
  python su2_string_tension.py --mode timing --L 8 --beta 2.5
  python su2_string_tension.py --mode production --L 8 --beta 2.2 --ncfg 300 ...
"""
import sys, os, time, json, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w4_blocked as W
import su3 as S3
from su3 import su2_heatbath, ID2

adj = W.adj


# ---------------------------------------------------------------- measurement
def rect_loop_trace2(U, mu, nu, R, T):
    """Re tr(loop)/2 per site on the (mu,nu) plane, base-indexed array (L^4)."""
    L = U.shape[1]
    A = np.broadcast_to(ID2, (L, L, L, L, 2, 2))
    for k in range(R):
        A = A @ np.roll(U[mu], -k, axis=mu)
    for l in range(T):
        A = A @ np.roll(np.roll(U[nu], -R, axis=mu), -l, axis=nu)
    for k in range(R - 1, -1, -1):
        A = A @ adj(np.roll(np.roll(U[mu], -k, axis=mu), -T, axis=nu))
    for l in range(T - 1, -1, -1):
        A = A @ adj(np.roll(U[nu], -l, axis=nu))
    return np.real(np.trace(A, axis1=-2, axis2=-1)) / 2.0


def wilson_loop2(U, R, T):
    """time-space average (3 spatial directions, all base sites)."""
    vals = []
    for nu in (1, 2, 3):
        vals.append(rect_loop_trace2(U, 0, nu, R, T).mean())
    return float(np.mean(vals))


def mag_op_time(U):
    """0++ idensity per time slice: sum over each t of spatial plaquettes.
    O(t) = sum_{x1,x2,x3} sum_{1<=i<j<=3} Re tr U_ij / 2 ; returns array (L,)."""
    L = U.shape[1]
    s = np.zeros((L, L, L, L))
    for i in (1, 2, 3):
        for j in range(i + 1, 4):
            Up = (U[i] @ np.roll(U[j], -1, axis=i)
                  @ adj(np.roll(U[i], -1, axis=j)) @ adj(U[j]))
            s += np.real(np.trace(Up, axis1=-2, axis2=-1)) / 2.0
    # s has shape (L0,L1,L2,L3); sum over spatial x1,x2,x3 (axes 1,2,3)
    return s.sum(axis=(1, 2, 3))


# ---------------------------------------------------------------- sampler
def make_idxp(L):
    ns = L ** 4
    x0, x1, x2, x3 = np.meshgrid(*[np.arange(L)] * 4, indexing='ij')
    par = (x0 + x1 + x2 + x3) % 2
    return ns, [np.nonzero(par.reshape(-1) == c)[0] for c in (0, 1)]


def sweep2(U, kappa, rng, idxp, ns):
    Uf = U.reshape(4, ns, 2, 2)
    for c in (0, 1):
        ids = idxp[c]
        for mu in range(4):
            S = W.staple_sum2(U, mu).reshape(ns, 2, 2)[ids]
            A = Uf[mu, ids]
            M = A @ S
            v = su2_heatbath(kappa, M, 0, 1, rng)
            Uf[mu, ids] = v @ A


def measure_cfg(U, Rmax, Tmax, smears):
    """Return dict: ('W',R,T)->float, ('O',n)->array(L)."""
    out = {}
    for R in range(1, Rmax + 2):
        for T in range(1, Tmax + 2):
            out[('W', R, T)] = wilson_loop2(U, R, T)
    prev = 0
    for n in smears:
        Us = U if n == 0 else S3.ape_smear_spatial(U, 0.5, n) if n == 1 else None
        # avoid repeated smearing: handled by caller below
        out[('O', n)] = mag_op_time(Us)
    return out


# ---------------------------------------------------------------- driver
def run(beta, L, ncfg, ntherm, gap, seed, Rmax, Tmax, smears, outpath,
        timing=False):
    rng = np.random.default_rng(seed)
    kappa = beta / 2.0
    ns, idxp = make_idxp(L)
    U = W.random_su2(rng, (4, L, L, L, L))
    t0 = time.time()
    for _ in range(ntherm):
        sweep2(U, kappa, rng, idxp, ns)
    t_therm = time.time() - t0
    print(f"[beta={beta} L={L}] therm {ntherm} sweeps in {t_therm:.1f}s "
          f"({t_therm/ntherm*1e3:.1f} ms/sweep)", flush=True)

    Wvals = {}
    Ovals = {n: [] for n in smears}
    plaq = []
    tmeas = 0.0
    for k in range(ncfg):
        for _ in range(gap):
            sweep2(U, kappa, rng, idxp, ns)
        tm = time.time()
        for R in range(1, Rmax + 2):
            for T in range(1, Tmax + 2):
                Wvals.setdefault((R, T), []).append(wilson_loop2(U, R, T))
        # smearing chain
        Us = U
        prev = 0
        for n in sorted(smears):
            if n > prev:
                Us = S3.ape_smear_spatial(U, 0.5, n) if prev == 0 else \
                     S3.ape_smear_spatial(Us, 0.5, n - prev)
                prev = n
            Ovals[n].append(mag_op_time(Us))
        plaq.append(S3.mean_plaquette(U))
        tmeas += time.time() - tm
        if (k + 1) % max(1, ncfg // 10) == 0:
            print(f"  [{k+1}/{ncfg}] plaq={plaq[-1]:.5f} "
                  f"(meas {tmeas/(k+1):.2f}s/cfg)", flush=True)
        if timing and k >= 2:
            break
    print(f"[beta={beta}] total {time.time()-t0:.1f}s", flush=True)

    nsave = len(plaq)
    data = dict(beta=beta, L=L, ncfg=nsave, plaq=np.array(plaq),
                W=np.array([[Wvals[(R, T)] for T in range(1, Tmax + 2)]
                            for R in range(1, Rmax + 2)]),
                smears=np.array(sorted(smears)))
    for n in smears:
        data[f'O_{n}'] = np.array(Ovals[n])
    if outpath:
        np.savez(outpath, **data)
        print(f"  saved -> {outpath}", flush=True)
    return data


# ---------------------------------------------------------------- analysis
def creutz(w, R, T):
    """chi(R,T) from scalar W(R,T) values; returns nan if any <=0."""
    a = w.get((R, T)); b = w.get((R-1, T-1)); c = w.get((R, T-1)); d = w.get((R-1, T))
    if min(a, b, c, d) <= 0:
        return np.nan
    return -np.log(a * b / (c * d))


def jackknife_creutz(Warr, R, T, nbin):
    """Warr shape (Rmax,Tmax,ncfg). Binned jackknife over configs."""
    n = Warr.shape[-1]
    w = {}
    for RR in range(1, Warr.shape[0] + 1):
        for TT in range(1, Warr.shape[1] + 1):
            w[(RR, TT)] = Warr[RR-1, TT-1]
    full = creutz({k: v.mean() for k, v in w.items()}, R, T)
    nb = max(1, n // nbin)
    # split into nbin blocks as equal as possible
    edges = np.linspace(0, n, nbin + 1).astype(int)
    vals = []
    for b in range(nbin):
        idx = np.r_[0:edges[b], edges[b+1]:n]
        if len(idx) == 0:
            continue
        wb = {k: v[idx].mean() for k, v in w.items()}
        vals.append(creutz(wb, R, T))
    vals = np.array(vals)
    vals = vals[np.isfinite(vals)]
    if len(vals) < 2:
        return full, np.nan
    m = vals.mean()
    err = np.sqrt((len(vals) - 1) / len(vals) * np.sum((vals - m) ** 2))
    return full, err


def fit_sigma(Wmean, Rmax, Tmax, rmin=2):
    """Fit ln W(R,T) = c - sigma*R*T - p*(R+T) over plateau (R,T>=rmin, W>0)."""
    xs, area, perim, y = [], [], [], []
    for R in range(1, Rmax + 1):
        for T in range(1, Tmax + 1):
            v = Wmean[R-1, T-1]
            if v > 0 and R >= rmin and T >= rmin:
                xs.append((R, T)); area.append(R * T); perim.append(R + T)
                y.append(np.log(v))
    A = np.vstack([np.ones(len(y)), np.array(perim), np.array(area)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return -coef[2], coef  # sigma


def cosh_fit(C, Lt):
    """fit C(t) = A[exp(-m t)+exp(-m(Lt-t))] + b for t=0..Lt//2; return m."""
    ts = np.arange(Lt // 2 + 1)
    y = C[:Lt // 2 + 1]
    mask = np.isfinite(y) & (np.abs(y) > 0)
    ts, y = ts[mask], y[mask]
    if len(ts) < 3:
        return np.nan
    best = (np.inf, np.nan)
    for m in np.linspace(0.05, 3.0, 300):
        f = np.exp(-m * ts) + np.exp(-m * (Lt - ts))
        A = np.vstack([f, np.ones_like(f)]).T
        coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
        r = np.sum((A @ coef - y) ** 2)
        if r < best[0]:
            best = (r, m)
    return best[1]


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', default='timing')
    ap.add_argument('--L', type=int, default=8)
    ap.add_argument('--beta', type=float, default=2.5)
    ap.add_argument('--ncfg', type=int, default=100)
    ap.add_argument('--ntherm', type=int, default=400)
    ap.add_argument('--gap', type=int, default=10)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--Rmax', type=int, default=5)
    ap.add_argument('--Tmax', type=int, default=5)
    ap.add_argument('--smears', type=str, default='0,4,8')
    ap.add_argument('--out', type=str, default='')
    args = ap.parse_args()
    smears = [int(x) for x in args.smears.split(',')]
    out = args.out or f'D:/Dev/math/research/lattice-cpb/sigma_L{args.L}_b{args.beta}.npz'
    run(args.beta, args.L, args.ncfg, args.ntherm, args.gap, args.seed,
        args.Rmax, args.Tmax, smears, out, timing=(args.mode == 'timing'))
