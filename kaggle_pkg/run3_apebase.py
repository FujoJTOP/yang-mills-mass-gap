"""RUN3 — fix the operator basis (keep time correlation) + N = 1.2e4 at L=8, prange over chains.

Two run1 limitations addressed:
  (1) S/N(2) at L=8 was 5.9 with N=4152 -> need N ~ 1.2e4.
  (2) 3D-APE n>=1 killed the time correlation (C(1)/C(0) ~ 1e-3, C(2)<0 for some n).
      New basis here: *spatial* APE, time direction untouched (operator at slice t is a
      function of the spatial links at slice t only), one/two steps, staple sum normalised
      by the actual staple count (4 in 3D, so alpha/4, NOT alpha/6 as in run1), and a
      guarded SU(2) projection.

Bases measured in the same Markov chain (no extra dynamics cost):
      lev 0: raw spatial plaquette (run1 baseline, same estimator)
      lev 1: T-APE n=1   (spatial only, alpha=0.5, /4 normalisation)
      lev 2: T-APE n=2   (same, iterated once more)

Estimators identical to run1:  S(t)=sum_{space} P(x,t),  C(t)=(1/Lt) sum_{t0}<S(t0+t)S(t0)>_c,
w1 = A/C(0) from C(t)=A(e^{-mt}+e^{-m(Lt-t)}) on t=1,2,3; errors = 20-block delete-one jackknife;
S/N(2) = C(2)/sigma_JK[C(2)].

Run:  py -3.12 run3_apebase.py --Ls 8 --Lt 12 --beta 2.5 --n 12000 --gap 3 --therm 400 --chains 8
Writes run3_apebase.{log,npz,json} next to this file.
"""
import os
import sys
import time
import json
import argparse

import numpy as np
from numba import njit, prange

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


# ------------------------------------------------------------------ RNG (per chain)
@njit(cache=True, inline='always')
def _u01(rng, c):
    s0 = rng[c, 0]
    s1 = rng[c, 1]
    x = s0
    y = s1
    rng[c, 0] = y
    x = x ^ (x << 23)
    rng[c, 1] = x ^ y ^ (x >> 17) ^ (y >> 26)
    return np.float64((rng[c, 1] + y) >> 11) * 1.1102230246251565e-16


@njit(cache=True, inline='always')
def _gauss(rng, c):
    u1 = _u01(rng, c)
    if u1 < 1e-300:
        u1 = 1e-300
    u2 = _u01(rng, c)
    return np.sqrt(-2.0 * np.log(u1)) * np.cos(6.283185307179586 * u2)


# ------------------------------------------------------------------ njit core
@njit(cache=True, inline='always')
def _mul(a00, a01, a10, a11, b00, b01, b10, b11):
    return (a00 * b00 + a01 * b10, a00 * b01 + a01 * b11,
            a10 * b00 + a11 * b10, a10 * b01 + a11 * b11)


@njit(cache=True, inline='always')
def _dag(a00, a01, a10, a11):
    return (np.conj(a00), np.conj(a10), np.conj(a01), np.conj(a11))


@njit(cache=True, fastmath=True, inline='always')
def _x0(lam, rng, c):
    """p(x) ~ sqrt(1-x^2) exp(lam x) on [-1,1] (Kennedy-Pendleton)."""
    if lam < 1e-6:
        while True:
            x = 2.0 * _u01(rng, c) - 1.0
            if _u01(rng, c) < np.sqrt(max(0.0, 1.0 - x * x)):
                return x
    e = np.exp(-2.0 * lam)
    while True:
        u = _u01(rng, c)
        x = 1.0 + np.log(u + (1.0 - u) * e) / lam
        r = _u01(rng, c)
        if r * r < 1.0 - x * x:
            return x


@njit(cache=True, fastmath=True)
def _sweep(U, plus, minus, kappa, nstep, rng, c):
    """U: (4, ns, 2, 2) complex128; time axis fastest. Sequential single-link heat bath."""
    ns = U.shape[1]
    for _ in range(nstep):
        for mu in range(4):
            for i in range(ns):
                s00 = 0j
                s01 = 0j
                s10 = 0j
                s11 = 0j
                for nu in range(4):
                    if nu == mu:
                        continue
                    ipm = plus[mu, i]
                    ipn = plus[nu, i]
                    imn = minus[nu, i]
                    imp = minus[nu, ipm]
                    x00, x01, x10, x11 = _mul(U[nu, ipm, 0, 0], U[nu, ipm, 0, 1],
                                              U[nu, ipm, 1, 0], U[nu, ipm, 1, 1],
                                              np.conj(U[mu, ipn, 0, 0]), np.conj(U[mu, ipn, 1, 0]),
                                              np.conj(U[mu, ipn, 0, 1]), np.conj(U[mu, ipn, 1, 1]))
                    y00, y01, y10, y11 = _mul(x00, x01, x10, x11,
                                              np.conj(U[nu, i, 0, 0]), np.conj(U[nu, i, 1, 0]),
                                              np.conj(U[nu, i, 0, 1]), np.conj(U[nu, i, 1, 1]))
                    s00 += y00
                    s01 += y01
                    s10 += y10
                    s11 += y11
                    x00, x01, x10, x11 = _mul(np.conj(U[nu, imp, 0, 0]), np.conj(U[nu, imp, 1, 0]),
                                              np.conj(U[nu, imp, 0, 1]), np.conj(U[nu, imp, 1, 1]),
                                              np.conj(U[mu, imn, 0, 0]), np.conj(U[mu, imn, 1, 0]),
                                              np.conj(U[mu, imn, 0, 1]), np.conj(U[mu, imn, 1, 1]))
                    y00, y01, y10, y11 = _mul(x00, x01, x10, x11,
                                              U[nu, imn, 0, 0], U[nu, imn, 0, 1],
                                              U[nu, imn, 1, 0], U[nu, imn, 1, 1])
                    s00 += y00
                    s01 += y01
                    s10 += y10
                    s11 += y11
                a00, a01, a10, a11 = _mul(U[mu, i, 0, 0], U[mu, i, 0, 1],
                                          U[mu, i, 1, 0], U[mu, i, 1, 1],
                                          s00, s01, s10, s11)
                b0 = 0.5 * (a00 + a11)
                b1 = 0.5 * (a01 + a10)
                b2 = 0.5j * (a01 - a10)
                b3 = 0.5 * (a00 - a11)
                q0 = b0.real
                q1 = -b1.imag
                q2 = -b2.imag
                q3 = -b3.imag
                an = np.sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3)
                x0 = _x0(2.0 * kappa * an, rng, c)
                g1 = _gauss(rng, c)
                g2 = _gauss(rng, c)
                g3 = _gauss(rng, c)
                gn = np.sqrt(g1 * g1 + g2 * g2 + g3 * g3)
                if gn < 1e-12:
                    g1, g2, g3, gn = 1.0, 0.0, 0.0, 1.0
                r = np.sqrt(max(0.0, 1.0 - x0 * x0))
                w1 = r * g1 / gn
                w2 = r * g2 / gn
                w3 = r * g3 / gn
                if an < 1e-14:
                    v00 = x0 + 1j * w3
                    v01 = 1j * w1 + w2
                    v10 = 1j * w1 - w2
                    v11 = x0 - 1j * w3
                else:
                    h0 = q0 / an
                    h1 = q1 / an
                    h2 = q2 / an
                    h3 = q3 / an
                    p00 = h0 + 1j * h3
                    p01 = 1j * h1 + h2
                    p10 = 1j * h1 - h2
                    p11 = h0 - 1j * h3
                    v00, v01, v10, v11 = _mul(p00, p01, p10, p11,
                                              x0 + 1j * w3, 1j * w1 + w2,
                                              1j * w1 - w2, x0 - 1j * w3)
                U[mu, i, 0, 0], U[mu, i, 0, 1], U[mu, i, 1, 0], U[mu, i, 1, 1] = _mul(
                    v00, v01, v10, v11,
                    U[mu, i, 0, 0], U[mu, i, 0, 1], U[mu, i, 1, 0], U[mu, i, 1, 1])
    return U


# ------------------------------------------------------------------ operators
@njit(cache=True, fastmath=True, inline='always')
def _plaq_site(U, plus, i, mu, nu):
    a00, a01, a10, a11 = _mul(U[mu, i, 0, 0], U[mu, i, 0, 1], U[mu, i, 1, 0], U[mu, i, 1, 1],
                              U[nu, plus[mu, i], 0, 0], U[nu, plus[mu, i], 0, 1],
                              U[nu, plus[mu, i], 1, 0], U[nu, plus[mu, i], 1, 1])
    q = plus[nu, i]
    b00, b01, b10, b11 = _dag(U[mu, q, 0, 0], U[mu, q, 0, 1], U[mu, q, 1, 0], U[mu, q, 1, 1])
    c00, c01, c10, c11 = _mul(a00, a01, a10, a11, b00, b01, b10, b11)
    d00, d01, d10, d11 = _dag(U[nu, i, 0, 0], U[nu, i, 0, 1], U[nu, i, 1, 0], U[nu, i, 1, 1])
    e00, e01, e10, e11 = _mul(c00, c01, c10, c11, d00, d01, d10, d11)
    return 0.5 * (e00.real + e11.real)


@njit(cache=True, fastmath=True)
def _slices(U, plus, Lt, row):
    """row(t) = (1/3) sum_{space} sum_{3 spatial planes} 1/2 Re tr U_p  (zero momentum)."""
    ns = U.shape[1]
    for t in range(Lt):
        row[t] = 0.0
    for i in range(ns):
        s = 0.0
        for mu in range(3):
            for nu in range(mu + 1, 3):
                s += _plaq_site(U, plus, i, mu, nu)
        row[i % Lt] += s / 3.0
    return row


@njit(cache=True, fastmath=True)
def _plaq_all(U, plus, Lt):
    """<1/2 Re tr U_p> over all 6 planes (sanity vs literature)."""
    ns = U.shape[1]
    tot = 0.0
    for i in range(ns):
        for mu in range(4):
            for nu in range(mu + 1, 4):
                tot += _plaq_site(U, plus, i, mu, nu)
    return tot / (6.0 * ns)


@njit(cache=True, fastmath=True)
def _ape3d(U, plus, minus, alpha, n):
    """Spatial-only APE (time direction untouched). Staple sum over the 2 transverse
    spatial directions (4 staples) normalised by alpha/4. Guarded SU(2) projection."""
    ns = U.shape[1]
    for _ in range(n):
        for mu in range(3):
            for i in range(ns):
                s00 = 0j
                s01 = 0j
                s10 = 0j
                s11 = 0j
                for nu in range(3):
                    if nu == mu:
                        continue
                    ipm = plus[mu, i]
                    ipn = plus[nu, i]
                    imn = minus[nu, i]
                    imp = minus[nu, ipm]
                    x00, x01, x10, x11 = _mul(U[nu, ipm, 0, 0], U[nu, ipm, 0, 1],
                                              U[nu, ipm, 1, 0], U[nu, ipm, 1, 1],
                                              np.conj(U[mu, ipn, 0, 0]), np.conj(U[mu, ipn, 1, 0]),
                                              np.conj(U[mu, ipn, 0, 1]), np.conj(U[mu, ipn, 1, 1]))
                    y00, y01, y10, y11 = _mul(x00, x01, x10, x11,
                                              np.conj(U[nu, i, 0, 0]), np.conj(U[nu, i, 1, 0]),
                                              np.conj(U[nu, i, 0, 1]), np.conj(U[nu, i, 1, 1]))
                    s00 += y00
                    s01 += y01
                    s10 += y10
                    s11 += y11
                    x00, x01, x10, x11 = _mul(np.conj(U[nu, imp, 0, 0]), np.conj(U[nu, imp, 1, 0]),
                                              np.conj(U[nu, imp, 0, 1]), np.conj(U[nu, imp, 1, 1]),
                                              np.conj(U[mu, imn, 0, 0]), np.conj(U[mu, imn, 1, 0]),
                                              np.conj(U[mu, imn, 0, 1]), np.conj(U[mu, imn, 1, 1]))
                    y00, y01, y10, y11 = _mul(x00, x01, x10, x11,
                                              U[nu, imn, 0, 0], U[nu, imn, 0, 1],
                                              U[nu, imn, 1, 0], U[nu, imn, 1, 1])
                    s00 += y00
                    s01 += y01
                    s10 += y10
                    s11 += y11
                o00 = (1.0 - alpha) * U[mu, i, 0, 0] + (alpha / 4.0) * s00
                o01 = (1.0 - alpha) * U[mu, i, 0, 1] + (alpha / 4.0) * s01
                o10 = (1.0 - alpha) * U[mu, i, 1, 0] + (alpha / 4.0) * s10
                o11 = (1.0 - alpha) * U[mu, i, 1, 1] + (alpha / 4.0) * s11
                det = o00 * o11 - o01 * o10
                ad = det.real * det.real + det.imag * det.imag
                if ad > 1e-12 and ad > 0.0:
                    f = 1.0 / np.sqrt(det)
                    U[mu, i, 0, 0] = o00 * f
                    U[mu, i, 0, 1] = o01 * f
                    U[mu, i, 1, 0] = o10 * f
                    U[mu, i, 1, 1] = o11 * f
    return U


# ------------------------------------------------------------------ geometry
def index_tables(Ls, Lt):
    ns = Ls * Ls * Ls * Lt
    idx = np.arange(ns)
    dims = (Ls, Ls, Ls, Lt)
    strd = [Ls * Ls * Lt, Ls * Lt, Lt, 1]
    plus = np.empty((4, ns), np.int64)
    minus = np.empty((4, ns), np.int64)
    for d in range(4):
        x = (idx // strd[d]) % dims[d]
        s = strd[d]
        plus[d] = np.where(x < dims[d] - 1, idx + s, idx - (dims[d] - 1) * s)
        minus[d] = np.where(x > 0, idx - s, idx + (dims[d] - 1) * s)
    return plus, minus


def hot(nch, Ls, Lt, seed):
    rng = np.random.default_rng(seed)
    ns = Ls * Ls * Ls * Lt
    z = (rng.standard_normal((nch, 4, ns, 2, 2)) + 1j * rng.standard_normal((nch, 4, ns, 2, 2)))
    q, r = np.linalg.qr(z)
    d = np.diagonal(r, axis1=-2, axis2=-1)
    q = q * (d / np.abs(d))[..., None, :]
    det = np.linalg.det(q)
    q[..., 0, :] *= np.conj(det)[..., None]
    return np.ascontiguousarray(q, dtype=np.complex128)


@njit(cache=True, fastmath=True, parallel=True)
def _chains(U, rng, plus, minus, kappa, per, gap, therm, alpha, nlev, Lt, S, plaq):
    nch = U.shape[0]
    for c in prange(nch):
        _sweep(U[c], plus, minus, kappa, therm, rng, c)
        buf = np.empty_like(U[c])
        for k in range(per):
            j = c * per + k
            _sweep(U[c], plus, minus, kappa, gap, rng, c)
            plaq[j] = _plaq_all(U[c], plus, Lt)
            _slices(U[c], plus, Lt, S[j, 0])
            buf[:, :, :, :] = U[c]                      # COPY: never mutate the chain
            for lev in range(1, nlev):
                _ape3d(buf, plus, minus, alpha, 1)
                _slices(buf, plus, Lt, S[j, lev])


# ------------------------------------------------------------------ analysis
def correlators(Sm, Lt):
    N = Sm.shape[0]
    x = Sm - Sm.mean(axis=0)
    C = np.zeros(Lt)
    for t in range(Lt):
        C[t] = np.mean(np.sum(x * np.roll(x, -t, axis=1), axis=1)) / Lt
    return C


def fit_A(C, Lt, ts=None):
    if ts is None:
        ts = np.arange(1, Lt // 2 + 1)
    y = C[ts]
    ms = np.linspace(0.05, 4.0, 4000)
    g = np.exp(-ms[:, None] * ts[None, :]) + np.exp(-ms[:, None] * (Lt - ts)[None, :])
    A = (g @ y) / np.sum(g * g, axis=1)
    r = np.sum((y[None, :] - A[:, None] * g) ** 2, axis=1)
    i = int(np.argmin(r))
    return float(ms[i]), float(A[i]), (i == 0 or i == len(ms) - 1)


def w1_of(C, Lt):
    m, A, _ = fit_A(C, Lt)
    return A / C[0], m, A


def A_at(C, Lt, m):
    ts = np.arange(1, 4)
    g = np.exp(-m * ts) + np.exp(-m * (Lt - ts))
    return float(np.sum(g * C[ts]) / np.sum(g * g))


def jack(Sm, Lt, nb=20, fn=None):
    N = Sm.shape[0]
    idx = np.array_split(np.arange(N), nb)
    vals = []
    for b in range(nb):
        keep = np.concatenate([idx[j] for j in range(nb) if j != b])
        try:
            vals.append(fn(correlators(Sm[keep], Lt)))
        except Exception:
            vals.append(np.nan)
    v = np.array(vals, float)
    m = np.nanmean(v)
    e = np.sqrt((nb - 1) / nb * np.nansum((v - m) ** 2))
    return m, e, int(np.sum(np.isfinite(v)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--Ls', type=int, default=8)
    ap.add_argument('--Lt', type=int, default=12)
    ap.add_argument('--beta', type=float, default=2.5)
    ap.add_argument('--n', type=int, default=12000)
    ap.add_argument('--gap', type=int, default=3)
    ap.add_argument('--therm', type=int, default=400)
    ap.add_argument('--chains', type=int, default=8)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--alpha', type=float, default=0.5)
    ap.add_argument('--nlev', type=int, default=3)
    ap.add_argument('--seed', type=int, default=20261125)
    ap.add_argument('--tag', type=str, default='apebase')
    ap.add_argument('--budget', type=float, default=900.0)
    a = ap.parse_args()

    import numba
    numba.set_num_threads(a.threads)
    Ls, Lt, beta = a.Ls, a.Lt, a.beta
    nch = a.chains
    per = a.n // nch
    N = per * nch
    tag = a.tag
    P("=" * 84)
    P(f"RUN3  SU(2) Wilson  beta={beta}  Ls={Ls}^3 x Lt={Lt}  N={N} ({nch} chains x {per}) "
      f"gap={a.gap} therm={a.therm}")
    P(f"basis: lev0=raw | lev1..{a.nlev-1}=T-APE(spatial only, alpha={a.alpha}, /4 staples, time untouched)")
    P(f"numba={numba.__version__} threads={numba.get_num_threads()}  w1=A/C(0)  S/N(2)=C(2)/sigma_JK")
    P("=" * 84)

    plus, minus = index_tables(Ls, Lt)
    U = hot(nch, Ls, Lt, a.seed)
    rng = np.empty((nch, 2), np.uint64)
    ss = np.uint64(a.seed)
    for c in range(nch):
        ss = ss * np.uint64(6364136223846793005) + np.uint64(1442695040888963407)
        rng[c, 0] = ss | np.uint64(1)
        s1 = ss ^ np.uint64(0x9E3779B97F4A7C15)
        rng[c, 1] = s1 if s1 != np.uint64(0) else np.uint64(0x2545F4914F6CDD1D)
    S = np.zeros((N, a.nlev, Lt))
    plaq = np.zeros(N)
    P(f"[T] Ls={Ls} Lt={Lt} sites={Ls**3*Lt} links={4*Ls**3*Lt} U bytes={U.nbytes/1e6:.2f} MB/ch")

    t = time.time()
    _chains(U[:1], rng[:1], plus, minus, beta / 2.0, 2, a.gap, 2, a.alpha, a.nlev, Lt,
            np.zeros((2, a.nlev, Lt)), np.zeros(2))          # compile (1 chain)
    P(f"[T] compile+smoke: {time.time()-t:.1f}s")

    t = time.time()
    _chains(U, rng, plus, minus, beta / 2.0, per, a.gap, a.therm, a.alpha, a.nlev, Lt, S, plaq)
    dtw = time.time() - t
    P(f"[E] production: {N} cfg x gap {a.gap} + therm {a.therm} in {dtw:.1f}s "
      f"({N/dtw:.1f} cfg/s, {nch} chains)  wall={time.time()-T0:.0f}s")
    P(f"[E] plaq={plaq.mean():.5f} +- {plaq.std(ddof=1)/np.sqrt(N):.5f}  "
      f"autocorr(lag1)={np.corrcoef(plaq[:-1], plaq[1:])[0,1]:.3f}")

    names = {0: 'raw'}
    for lev in range(1, a.nlev):
        names[lev] = f'TAPE{lev}'
    res = {'n': N, 'chains': nch, 'Ls': Ls, 'Lt': Lt, 'beta': beta, 'gap': a.gap,
           'therm': a.therm, 'alpha': a.alpha, 'levels': {}, 'plaq': float(plaq.mean())}
    P("")
    P(f"{'lvl':>6} {'C(0)':>10} {'C1/C0':>9} {'C2/C0':>9} {'m_fit':>7} {'w1=A/C0':>9} "
      f"{'+-':>8} {'w1@1.546':>9} {'S/N(2)':>7} {'Rbar':>7}")
    for lev in range(a.nlev):
        C = correlators(S[:, lev, :], Lt)
        if C[0] <= 0:
            P(f"{names[lev]:>6} C(0)={C[0]:.3e} not positive -> w1 undefined")
            res['levels'][names[lev]] = dict(C0=float(C[0]), defined=False)
            continue
        w, m, A = w1_of(C, Lt)
        wm, we, nok = jack(S[:, lev, :], Lt, 20, lambda Cc: w1_of(Cc, Lt)[0])
        _, s2, _ = jack(S[:, lev, :], Lt, 20, lambda Cc: Cc[2])
        _, s1, _ = jack(S[:, lev, :], Lt, 20, lambda Cc: Cc[1])
        snr = C[2] / s2 if s2 > 0 else float('inf')
        A154 = A_at(C, Lt, 1.546)
        r_ln = -np.log(max(C[1] / C[0], 1e-300))
        P(f"{names[lev]:>6} {C[0]:10.4f} {C[1]/C[0]:9.5f} {C[2]/C[0]:9.5f} {m:7.3f} "
          f"{w:9.4f} {we:8.4f} {A154/C[0]:9.4f} {snr:7.1f} {r_ln:7.4f}")
        P(f"       [E] sigma(C1)={s1:.2e} -> C1/C0 = {C[1]/C[0]:.5f} +- {s1/C[0]:.5f} ; "
          f"C2/C0 = {C[2]/C[0]:.5f} +- {abs(s2)/C[0]:.5f}")
        res['levels'][names[lev]] = dict(
            C0=float(C[0]), C1=float(C[1] / C[0]), C2=float(C[2] / C[0]),
            C3=float(C[3] / C[0]), C1_err=float(s1 / C[0]), C2_err=float(abs(s2) / C[0]),
            m=m, A=A, w1=w, w1_err=we, w1_jk=float(wm),
            w1_m1546=float(A154 / C[0]), snr2=float(snr), rbar=float(r_ln),
            defined=True, jack_ok=nok, Cfull=C.tolist(), sigma2=float(s2))

    # same-chain comparison  raw  vs  TAPE (correlated data, so also the pair difference)
    if 'raw' in res['levels'] and 'TAPE1' in res['levels']:
        C0 = correlators(S[:, 0, :], Lt)
        C1 = correlators(S[:, 1, :], Lt)
        w_r = A_at(C0, Lt, 1.546) / C0[0]
        w_t = A_at(C1, Lt, 1.546) / C1[0]
        res['w1_pair'] = dict(raw=float(w_r), tape1=float(w_t))
        P("")
        P(f"[C] same-chain w1@1.546: raw={w_r:.4f}  TAPE1={w_t:.4f}  diff={w_t-w_r:+.4f}")

    np.savez(os.path.join(HERE, f'run3_{tag}.npz'), S=S, plaq=plaq, beta=beta, Ls=Ls, Lt=Lt,
             n=N, gap=a.gap, therm=a.therm, alpha=a.alpha, chains=nch, seed=a.seed)
    json.dump(res, open(os.path.join(HERE, f'run3_{tag}.json'), 'w'), indent=1)
    P("")
    P(f"[wall] total {time.time()-T0:.1f}s ; wrote run3_{tag}.npz / .json / .log")
    open(os.path.join(HERE, f'run3_{tag}.log'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
