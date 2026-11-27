"""GH3 -- operator-basis probe at beta=2.5 (small-alpha APE, p!=0 projection, area>1 loops).

One MC job; on the SAME configurations measure 6 operator levels x 5 observables.

Levels (all: only spatial links are smeared, time index never changes, SU(2) det-projection
        exact for real weights; alpha/n passed straight to run3_apebase._ape3d):
   0  raw                (spatial 1x1 only)
   1  alpha=0.10, n=1
   2  alpha=0.20, n=1
   3  alpha=0.30, n=1
   4  alpha=0.50, n=1   (== archived TAPE1 / APE1 control)
   5  alpha=0.10, n=3

Observables per (config, level), stored as row[t] of a (5, Lt) array:
   0  W11(x) = (1/3) sum_{i<j<=3} 1/2 Re tr U_p                 , p=0
   1  as 0 but with phase (-1)^x                                , p=(pi,0,0)
   2  as 0 but with phase (-1)^(x+y+z)                          , p=(pi,pi,pi)
   3  W12(x) = (1/6) sum_{ordered (mu,nu)} 1/2 Re tr R_{mu,nu}  , p=0
   4  W22(x) = (1/3) sum_{i<j<=3} 1/2 Re tr P_{i,j} (2x2)       , p=0
with S_obs(t) = sum_{x in slice t} (weighted) W(x,t); C(t) = (1/Lt) sum_t0 <S(t0+t)S(t0)>_c.

Self-checks before production:
  * obs 0 must equal run3_apebase._plaq_site summed over the 3 planes, bit-for-bit;
  * every loop trace must be in [-1,1] (path closure test);
  * on a HOT configuration all five means must be ~0 (control);
  * plaq (all 6 planes) must reproduce 0.6523-0.6527 at beta=2.5 (archive agreement).

Run:  py -3.12 gh3_probe.py --Ls 6 --Lt 12 --beta 2.5 --n 8000 --gap 4 --therm 400 --chains 8
Writes gh3_probe.{log,json,npz} next to this file.
"""
import os
import sys
import time
import json
import argparse

import numpy as np
from numba import njit, prange

from run3_apebase import (_sweep, _mul, _dag, _plaq_site, _ape3d, _plaq_all,
                          index_tables, hot, correlators, jack)

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []
NOBS = 5
NAMES = ['raw', 'a10n1', 'a20n1', 'a30n1', 'a50n1', 'a10n3']
LEVCFG = [(0.0, 0), (0.10, 1), (0.20, 1), (0.30, 1), (0.50, 1), (0.10, 3)]
OBSN = ['W11_p0', 'W11_px', 'W11_pxxx', 'W12_p0', 'W22_p0']


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


# ------------------------------------------------------------------ loop kernels
@njit(cache=True, fastmath=True, inline='always')
def _rect2(U, plus, i, mu, nu):
    """1x2 rectangle: U_mu(x) U_mu(x+mu) U_nu(x+2mu) U_mu(x+mu+nu)^dag U_mu(x+nu)^dag U_nu(x)^dag."""
    x1 = plus[mu, i]
    x2 = plus[mu, x1]
    a00, a01, a10, a11 = _mul(U[mu, i, 0, 0], U[mu, i, 0, 1], U[mu, i, 1, 0], U[mu, i, 1, 1],
                              U[mu, x1, 0, 0], U[mu, x1, 0, 1], U[mu, x1, 1, 0], U[mu, x1, 1, 1])
    a00, a01, a10, a11 = _mul(a00, a01, a10, a11,
                              U[nu, x2, 0, 0], U[nu, x2, 0, 1], U[nu, x2, 1, 0], U[nu, x2, 1, 1])
    xmn = plus[nu, x1]                                  # x+mu+nu
    b00, b01, b10, b11 = _dag(U[mu, xmn, 0, 0], U[mu, xmn, 0, 1], U[mu, xmn, 1, 0], U[mu, xmn, 1, 1])
    a00, a01, a10, a11 = _mul(a00, a01, a10, a11, b00, b01, b10, b11)
    xn = plus[nu, i]                                    # x+nu
    b00, b01, b10, b11 = _dag(U[mu, xn, 0, 0], U[mu, xn, 0, 1], U[mu, xn, 1, 0], U[mu, xn, 1, 1])
    a00, a01, a10, a11 = _mul(a00, a01, a10, a11, b00, b01, b10, b11)
    b00, b01, b10, b11 = _dag(U[nu, i, 0, 0], U[nu, i, 0, 1], U[nu, i, 1, 0], U[nu, i, 1, 1])
    a00, a01, a10, a11 = _mul(a00, a01, a10, a11, b00, b01, b10, b11)
    return 0.5 * (a00.real + a11.real)


@njit(cache=True, fastmath=True, inline='always')
def _sq22(U, plus, i, mu, nu):
    """2x2 square (8 links), corner at x, in the (mu,nu) plane."""
    x1 = plus[mu, i]
    x2 = plus[mu, x1]
    x3 = plus[nu, x2]          # x+2mu+nu
    x2n = plus[nu, plus[nu, i]]  # x+2nu
    xmn = plus[mu, x2n]        # x+mu+2nu
    xn = plus[nu, i]           # x+nu
    a00, a01, a10, a11 = _mul(U[mu, i, 0, 0], U[mu, i, 0, 1], U[mu, i, 1, 0], U[mu, i, 1, 1],
                              U[mu, x1, 0, 0], U[mu, x1, 0, 1], U[mu, x1, 1, 0], U[mu, x1, 1, 1])
    a00, a01, a10, a11 = _mul(a00, a01, a10, a11,
                              U[nu, x2, 0, 0], U[nu, x2, 0, 1], U[nu, x2, 1, 0], U[nu, x2, 1, 1])
    a00, a01, a10, a11 = _mul(a00, a01, a10, a11,
                              U[nu, x3, 0, 0], U[nu, x3, 0, 1], U[nu, x3, 1, 0], U[nu, x3, 1, 1])
    b00, b01, b10, b11 = _dag(U[mu, xmn, 0, 0], U[mu, xmn, 0, 1], U[mu, xmn, 1, 0], U[mu, xmn, 1, 1])
    a00, a01, a10, a11 = _mul(a00, a01, a10, a11, b00, b01, b10, b11)
    b00, b01, b10, b11 = _dag(U[mu, x2n, 0, 0], U[mu, x2n, 0, 1], U[mu, x2n, 1, 0], U[mu, x2n, 1, 1])
    a00, a01, a10, a11 = _mul(a00, a01, a10, a11, b00, b01, b10, b11)
    b00, b01, b10, b11 = _dag(U[nu, xn, 0, 0], U[nu, xn, 0, 1], U[nu, xn, 1, 0], U[nu, xn, 1, 1])
    a00, a01, a10, a11 = _mul(a00, a01, a10, a11, b00, b01, b10, b11)
    b00, b01, b10, b11 = _dag(U[nu, i, 0, 0], U[nu, i, 0, 1], U[nu, i, 1, 0], U[nu, i, 1, 1])
    a00, a01, a10, a11 = _mul(a00, a01, a10, a11, b00, b01, b10, b11)
    return 0.5 * (a00.real + a11.real)


@njit(cache=True, fastmath=True)
def _obs(U, plus, sx, sxyz, Lt, out):
    """out[obs, t] = slice sums for the 5 observables (see module docstring)."""
    ns = U.shape[1]
    for k in range(out.shape[0]):
        for t in range(Lt):
            out[k, t] = 0.0
    for i in range(ns):
        t = i % Lt
        p = 0.0
        for mu in range(3):
            for nu in range(mu + 1, 3):
                p += _plaq_site(U, plus, i, mu, nu)
        p /= 3.0
        out[0, t] += p
        out[1, t] += p * sx[i]
        out[2, t] += p * sxyz[i]
        r = 0.0
        for mu in range(3):
            for nu in range(3):
                if nu != mu:
                    r += _rect2(U, plus, i, mu, nu)
        out[3, t] += r / 6.0
        s = 0.0
        for mu in range(3):
            for nu in range(mu + 1, 3):
                s += _sq22(U, plus, i, mu, nu)
        out[4, t] += s / 3.0


@njit(cache=True, fastmath=True, parallel=True)
def _mchains(U, rng, plus, minus, kappa, per, gap, therm, Lt, sx, sxyz, alphas, nsteps, S, plaq):
    """prange over chains; S[j, lev, obs, t]. Never mutates U[c]."""
    nch = U.shape[0]
    nlev = S.shape[1]
    for c in prange(nch):
        _sweep(U[c], plus, minus, kappa, therm, rng, c)
        buf = np.empty_like(U[c])
        for k in range(per):
            j = c * per + k
            _sweep(U[c], plus, minus, kappa, gap, rng, c)
            plaq[j] = _plaq_all(U[c], plus, Lt)
            _obs(U[c], plus, sx, sxyz, Lt, S[j, 0])          # lev 0 = raw
            for lev in range(1, nlev):
                buf[:, :, :, :] = U[c]                       # copy: chain stays clean
                _ape3d(buf, plus, minus, alphas[lev], nsteps[lev])
                _obs(buf, plus, sx, sxyz, Lt, S[j, lev])


# ------------------------------------------------------------------ helpers
def phase_tables(Ls, Lt):
    ns = Ls * Ls * Ls * Lt
    idx = np.arange(ns)
    x = (idx // (Ls * Ls * Lt)) % Ls
    y = (idx // (Ls * Lt)) % Ls
    z = (idx // Lt) % Ls
    sx = np.where(x % 2 == 0, 1.0, -1.0)
    sxyz = np.where((x + y + z) % 2 == 0, 1.0, -1.0)
    return sx, sxyz


def corr_pair(Si, Sj, Lt):
    """C_ij(t) = (1/Lt) sum_t0 <S_i(t0+t) S_j(t0)>_c, zero mode removed from both."""
    a = Si - Si.mean(axis=0)
    b = Sj - Sj.mean(axis=0)
    C = np.zeros(Lt)
    for t in range(Lt):
        C[t] = np.mean(np.sum(a * np.roll(b, -t, axis=1), axis=1)) / Lt
    return C


def gewp(C, Lt, t0):
    """Generalised eigenvalues of C(t) w.r.t. C(t0); returns dict t -> array of lam."""
    G = C[t0]
    w, V = np.linalg.eigh(G)
    keep = w > 1e-2 * w.max()
    Wm = V[:, keep] / np.sqrt(w[keep])[None, :]
    out = {}
    for t in range(t0 + 1, Lt):
        M = Wm.T @ C[t] @ Wm
        M = 0.5 * (M + M.T)
        lam = np.linalg.eigvalsh(M)
        out[t] = np.sort(lam)[::-1]
    return out, keep.sum(), w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--Ls', type=int, default=6)
    ap.add_argument('--Lt', type=int, default=12)
    ap.add_argument('--beta', type=float, default=2.5)
    ap.add_argument('--n', type=int, default=8000)
    ap.add_argument('--gap', type=int, default=4)
    ap.add_argument('--therm', type=int, default=400)
    ap.add_argument('--chains', type=int, default=8)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--seed', type=int, default=20261126)
    ap.add_argument('--tag', type=str, default='probe')
    ap.add_argument('--budget', type=float, default=280.0)
    a = ap.parse_args()

    import numba
    numba.set_num_threads(a.threads)
    Ls, Lt, beta = a.Ls, a.Lt, a.beta
    nch = a.chains
    per = a.n // nch
    N = per * nch
    nlev = len(NAMES)

    P("=" * 96)
    P(f"GH3  SU(2) Wilson  beta={beta}  Ls={Ls}^3 x Lt={Lt}  N={N} ({nch} chains x {per})  "
      f"gap={a.gap} therm={a.therm}  budget={a.budget}s")
    P(f"levels: {', '.join(f'{NAMES[i]}(alpha={LEVCFG[i][0]},n={LEVCFG[i][1]})' for i in range(nlev))}")
    P(f"obs   : {', '.join(OBSN)}   (normalisation: W11,W22 -> /3 planes ; W12 -> /6 oriented loops)")
    P(f"numba={numba.__version__} threads={numba.get_num_threads()}")
    P("=" * 96)

    plus, minus = index_tables(Ls, Lt)
    sx, sxyz = phase_tables(Ls, Lt)

    # ---------------- self-checks on a hot config ----------------
    Uc = hot(1, Ls, Lt, 12345)
    o = np.zeros((NOBS, Lt))
    _obs(Uc[0], plus, sx, sxyz, Lt, o)
    P(f"[T] hot-config means (control, must be ~0): " +
      " ".join(f"{OBSN[k]}={o[k].mean()/ (Ls**3):+.4f}" for k in range(NOBS)))
    # exact agreement with run3_apebase._plaq_site on the same config
    ref = 0.0
    for i in range(Ls ** 3 * Lt):
        s = 0.0
        for mu in range(3):
            for nu in range(mu + 1, 3):
                s += _plaq_site(Uc[0], plus, i, mu, nu)
        ref += s / 3.0
    P(f"[T] obs0 total={o[0].sum():.10f}  run3 ref={ref:.10f}  diff={o[0].sum()-ref:+.2e}")
    P(f"[T] hot plaq(all 6 planes)={_plaq_all(Uc[0], plus, Lt):+.5f} (must be ~0)")
    # loop traces bounded: use a thermalised-ish config instead (hot is too symmetric), so
    # check the bound on several random configs
    worst = 0.0
    for s0 in range(3):
        U2 = hot(1, Ls, Lt, 777 + s0)
        o2 = np.zeros((NOBS, Lt))
        _obs(U2[0], plus, sx, sxyz, Lt, o2)
        for k in range(NOBS):
            v = np.abs(o2[k] / (Ls ** 3)).max()
            worst = max(worst, v)
    P(f"[T] |loop trace| bound on random configs: max |mean| = {worst:.4f} (must be <= 1)")

    # ---------------- production ----------------
    U = hot(nch, Ls, Lt, a.seed)
    rng = np.empty((nch, 2), np.uint64)
    ss = np.uint64(a.seed)
    for c in range(nch):
        ss = ss * np.uint64(6364136223846793005) + np.uint64(1442695040888963407)
        rng[c, 0] = ss | np.uint64(1)
        s1 = ss ^ np.uint64(0x9E3779B97F4A7C15)
        rng[c, 1] = s1 if s1 != np.uint64(0) else np.uint64(0x2545F4914F6CDD1D)
    S = np.zeros((N, nlev, NOBS, Lt))
    plaq = np.zeros(N)
    P(f"[T] sites={Ls**3*Lt} links={4*Ls**3*Lt} U={U.nbytes/1e6:.2f} MB/ch "
      f"S={S.nbytes/1e6:.1f} MB")

    lev_alphas = np.array([x[0] for x in LEVCFG], dtype=np.float64)
    lev_nsteps = np.array([x[1] for x in LEVCFG], dtype=np.int64)
    t = time.time()
    _mchains(U[:1], rng[:1], plus, minus, beta / 2.0, 2, a.gap, 2, Lt, sx, sxyz,
             lev_alphas, lev_nsteps, np.zeros((2, nlev, NOBS, Lt)), np.zeros(2))
    P(f"[T] compile+smoke: {time.time()-t:.1f}s")

    t = time.time()
    _mchains(U, rng, plus, minus, beta / 2.0, per, a.gap, a.therm, Lt, sx, sxyz,
             lev_alphas, lev_nsteps, S, plaq)
    dtw = time.time() - t
    P(f"[E] production: {N} cfg x (gap {a.gap} + {nlev} levels x {NOBS} obs) in {dtw:.1f}s "
      f"({N/dtw:.1f} cfg/s)  wall={time.time()-T0:.0f}s")
    P(f"[E] plaq={plaq.mean():.5f} +- {plaq.std(ddof=1)/np.sqrt(N):.5f}  "
      f"autocorr(lag1)={np.corrcoef(plaq[:-1], plaq[1:])[0,1]:.3f}")
    if time.time() - T0 > a.budget:
        P(f"[W] wall {time.time()-T0:.0f}s exceeded budget {a.budget}s")

    np.savez(os.path.join(HERE, f'gh3_{a.tag}.npz'), S=S, plaq=plaq, beta=beta, Ls=Ls, Lt=Lt,
             n=N, gap=a.gap, therm=a.therm, chains=nch, seed=a.seed,
             lev_alpha=np.array([x[0] for x in LEVCFG]), lev_n=np.array([x[1] for x in LEVCFG]))
    P(f"[E] wrote gh3_{a.tag}.npz ({S.nbytes/1e6:.1f} MB) -- raw data safe before analysis")

    # ---------------- analysis (tables only; GEVP lives in gh3_analyze.py) ----------------
    res = dict(n=N, Ls=Ls, Lt=Lt, beta=beta, gap=a.gap, therm=a.therm, chains=nch,
               levels=[dict(name=NAMES[i], alpha=LEVCFG[i][0], n=LEVCFG[i][1]) for i in range(nlev)],
               obs=OBSN, plaq=float(plaq.mean()), wall=float(time.time() - T0), lev={})
    for lev in range(nlev):
        res['lev'][NAMES[lev]] = {}
        P("")
        P(f"### {NAMES[lev]}  (alpha={LEVCFG[lev][0]}, n={LEVCFG[lev][1]})")
        P(f"{'obs':>9} {'C(0)':>9} " + " ".join(f"{'Phi'+str(k):>10}" for k in range(1, 5)) +
          f" {'S/N(1)':>7} {'m_eff(1)':>9} {'m_eff(2)':>9}")
        for k in range(NOBS):
            Sm = S[:, lev, k, :]
            C = correlators(Sm, Lt)
            if C[0] <= 0:
                P(f"{OBSN[k]:>9} C(0)={C[0]:.3e} <=0")
                continue
            err = np.zeros(Lt)
            for tt in range(Lt):
                _, e, _ = jack(Sm, Lt, 20, lambda Cc, tt=tt: Cc[tt])
                err[tt] = e
            phis = [C[tt] / C[0] for tt in range(1, 5)]
            snr = C[1] / err[1] if err[1] > 0 else float('nan')
            m1 = np.log(C[1] / C[2]) if C[2] > 0 else float('nan')
            m2 = np.log(C[2] / C[3]) if C[3] > 0 and C[2] > 0 else float('nan')
            P(f"{OBSN[k]:>9} {C[0]:9.4f} " + " ".join(f"{v:10.5f}" for v in phis) +
              f" {snr:7.2f} {m1:9.4f} {m2:9.4f}")
            P(f"{'':>9} {'+-':>9} " +
              " ".join(f"{err[tt]/C[0]:10.5f}" for tt in range(1, 5)))
            res['lev'][NAMES[lev]][OBSN[k]] = dict(
                C0=float(C[0]), phi=[float(v) for v in phis],
                err=[float(err[tt] / C[0]) for tt in range(1, 5)],
                snr1=float(snr), m_eff1=float(m1), m_eff2=float(m2),
                Cfull=C.tolist())

    # ---------------- GEVP is in gh3_analyze.py (runs off the npz, no MC rerun) ----------
    json.dump(res, open(os.path.join(HERE, f'gh3_{a.tag}.json'), 'w'), indent=1)
    P("")
    P(f"[wall] total {time.time()-T0:.1f}s ; wrote gh3_{a.tag}.npz / .json / .log")
    open(os.path.join(HERE, f'gh3_{a.tag}.log'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
