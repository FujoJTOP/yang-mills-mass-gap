"""RUN1 — w1 = A/C(0) and Rbar = M1/M0 on the minimal confining box Ls=6 x Lt=12, beta=2.5.

SU(2) Wilson action S = beta * sum_p (1 - 1/2 Re tr U_p), kappa = beta/2.
njit sequential single-link heat bath (Kennedy-Pendleton for the x0 density),
one job at a time (CLAUDE.md sec.4: never parallelise heavy lattice jobs).

Operator ladder: APE-3D (spatial links only, alpha=0.5) n = 1,2,3,4  (+ n=0 raw).
  P(x) = (1/3) sum_{spatial planes i<j} 1/2 Re tr U_p(x)
  S_n(t) = sum_{x in time slice t} P_n(x,t)          (zero momentum)
  C_n(t) = (1/T) sum_{t0} <S(t0+t) S(t0)>_c          (time-averaged, periodic)
  M0 = C(0),  A from C(t) = A(e^{-mt} + e^{-m(T-t)}),  w1 = A/C(0)
  Rbar = -C'(0)/C(0)  (first spectral cumulant, energy units)
       = sum_t t C(t) / sum_t C(t)  (mean Euclidean time, M1/M0 with M0 = sum C)

Run:  py -3.12 run1_b25_l6.py --n 3000 --gap 4 --therm 300 --budget 900
Writes run1_b25_l6.npz, run1_b25_l6.log (same directory).
"""
import os
import sys
import time
import json
import argparse

import numpy as np
from numba import njit

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


# ---------------------------------------------------------------- njit core
@njit(cache=True, fastmath=True, inline='always')
def _mul(a00, a01, a10, a11, b00, b01, b10, b11):
    return (a00 * b00 + a01 * b10, a00 * b01 + a01 * b11,
            a10 * b00 + a11 * b10, a10 * b01 + a11 * b11)


@njit(cache=True, fastmath=True)
def _x0(lam):
    """p(x) ~ sqrt(1-x^2) exp(lam x) on [-1,1] (Kennedy-Pendleton)."""
    if lam < 1e-6:
        while True:
            x = 2.0 * np.random.random() - 1.0
            if np.random.random() < np.sqrt(max(0.0, 1.0 - x * x)):
                return x
    e = np.exp(-2.0 * lam)
    while True:
        u = np.random.random()
        x = 1.0 + np.log(u + (1.0 - u) * e) / lam
        r = np.random.random()
        if r * r < 1.0 - x * x:
            return x


@njit(cache=True, fastmath=True)
def sweep(U, plus, minus, kappa, nstep):
    """U: (4, ns, 2, 2) complex128, ns flattened with the time axis fastest."""
    nlink, ns = U.shape[0], U.shape[1]
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
                    # A = U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag
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
                    # B = U_nu(x+mu-nu)^dag U_mu(x-nu)^dag U_nu(x-nu)
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
                # B = A_link @ S
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
                x0 = _x0(2.0 * kappa * an)
                g1 = np.random.standard_normal()
                g2 = np.random.standard_normal()
                g3 = np.random.standard_normal()
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
                    # Ahat = h0 I + i(h1 s1 + h2 s2 + h3 s3)
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


# ---------------------------------------------------------------- geometry
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


def hot(Ls, Lt, rng):
    z = (rng.standard_normal((4, Ls * Ls * Ls * Lt, 2, 2))
         + 1j * rng.standard_normal((4, Ls * Ls * Ls * Lt, 2, 2)))
    q, r = np.linalg.qr(z)
    d = np.diagonal(r, axis1=-2, axis2=-1)
    q = q * (d / np.abs(d))[..., None, :]
    det = np.linalg.det(q)
    q[..., 0, :] *= np.conj(det)[..., None]
    return np.ascontiguousarray(q, dtype=np.complex128)


def adj(A):
    return np.conj(np.swapaxes(A, -1, -2))


def plaq_avg(U, Ls, Lt):
    """<1/2 Re tr U_p> over the 6 planes (sanity vs literature)."""
    U = U.reshape(4, Ls, Ls, Ls, Lt, 2, 2)
    tot, n = 0.0, 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Up = (U[mu] @ np.roll(U[nu], -1, axis=mu)
                  @ adj(np.roll(U[mu], -1, axis=nu)) @ adj(U[nu]))
            tot += 0.5 * np.real(np.trace(Up, axis1=-2, axis2=-1)).sum()
            n += U[mu].size // 4
    return tot / n


def ape3d(U, Ls, Lt, alpha=0.5):
    """One APE step on the 3 spatial link directions, 3D (spatial) staples."""
    Um = U.reshape(4, Ls, Ls, Ls, Lt, 2, 2).copy()
    for mu in range(3):
        S = np.zeros_like(Um[mu])
        for nu in range(3):
            if nu == mu:
                continue
            Un = Um[nu]
            S += np.roll(Un, -1, axis=mu) @ adj(np.roll(Um[mu], -1, axis=nu)) @ adj(Un)
            S += (adj(np.roll(np.roll(Un, 1, axis=nu), -1, axis=mu))
                  @ adj(np.roll(Um[mu], 1, axis=nu)) @ np.roll(Un, 1, axis=nu))
        M = (1.0 - alpha) * Um[mu] + (alpha / 6.0) * S
        det = M[..., 0, 0] * M[..., 1, 1] - M[..., 0, 1] * M[..., 1, 0]
        Um[mu] = M / np.sqrt(det)[..., None, None]
    return Um.reshape(4, -1, 2, 2)


def slices(U, Ls, Lt):
    """S(t) = sum over the time slice of the spatial-plaquette density."""
    U4 = U.reshape(4, Ls, Ls, Ls, Lt, 2, 2)
    P = np.zeros((Ls, Ls, Ls, Lt))
    for mu in range(3):
        for nu in range(mu + 1, 3):
            Up = (U4[mu] @ np.roll(U4[nu], -1, axis=mu)
                  @ adj(np.roll(U4[mu], -1, axis=nu)) @ adj(U4[nu]))
            P += 0.5 * np.real(np.trace(Up, axis1=-2, axis2=-1))
    return P.sum(axis=(0, 1, 2)) / 3.0


# ---------------------------------------------------------------- Markov
def markov(Ls, Lt, beta, n, gap, therm, budget, seed, levels=(1, 2, 3, 4)):
    rng = np.random.default_rng(seed)
    np.random.seed(seed)
    plus, minus = index_tables(Ls, Lt)
    U = hot(Ls, Lt, rng)
    kappa = beta / 2.0
    nl = len(levels) + 1                       # level 0 = raw
    S = np.zeros((n, nl, Lt))
    plaq = np.zeros(n)
    t = time.time()
    sweep(U, plus, minus, kappa, therm)        # jit compile + thermalise
    P(f"[T] thermalised: {therm} sweeps in {time.time() - t:.1f}s "
      f"(Ls={Ls} Lt={Lt} ns={Ls**3*Lt} links={4*Ls**3*Lt})")
    t = time.time()
    k = 0
    for k in range(n):
        sweep(U, plus, minus, kappa, gap)
        plaq[k] = plaq_avg(U, Ls, Lt)
        S[k, 0] = slices(U, Ls, Lt)
        Us = U
        for j, lev in enumerate(levels):
            Us = ape3d(Us, Ls, Lt)
            S[k, j + 1] = slices(Us, Ls, Lt)
        if k % 100 == 0:
            dt = time.time() - t
            rate = (k + 1) / max(dt, 1e-9)
            P(f"  [cfg {k + 1}/{n}] plaq={plaq[k]:.4f}  {rate:.2f} cfg/s  "
              f"eta {dt / (k + 1) * (n - k - 1) / 60:.1f} min  wall={time.time() - T0:.0f}s")
        if time.time() - T0 > budget:
            k += 1
            P(f"[C] wall budget {budget}s hit at cfg {k}/{n}")
            break
    else:
        k = n
    return S[:k], plaq[:k], k


# ---------------------------------------------------------------- analysis
def correlators(Sm, Lt):
    """Sm: (N, Lt) time-slice sums -> C(t), time-averaged."""
    N = Sm.shape[0]
    x = Sm - Sm.mean(axis=0)
    C = np.zeros(Lt)
    for t in range(Lt):
        C[t] = np.mean(np.sum(x * np.roll(x, -t, axis=1), axis=1)) / Lt
    return C


def fit_A(C, Lt, ts=None):
    """C(t)=A(e^{-mt}+e^{-m(T-t)}): grid m, linear A."""
    if ts is None:
        ts = np.arange(1, Lt // 2 + 1)
    y = C[ts]
    ms = np.linspace(0.05, 4.0, 4000)
    g = np.exp(-ms[:, None] * ts[None, :]) + np.exp(-ms[:, None] * (Lt - ts)[None, :])
    A = (g @ y) / np.sum(g * g, axis=1)
    r = np.sum((y[None, :] - A[:, None] * g) ** 2, axis=1)
    i = int(np.argmin(r))
    edge = i == 0 or i == len(ms) - 1
    return float(ms[i]), float(A[i]), edge


def rbar(C):
    """Rbar = -C'(0)/C(0): two discretisations + mean-time cross-check."""
    r_ln = -np.log(max(C[1] / C[0], 1e-300))
    r_fd = (3.0 * C[0] - 4.0 * C[1] + C[2]) / (2.0 * C[0])
    tt = np.arange(len(C))
    r_mt = float((tt * C).sum() / C.sum())
    return r_ln, r_fd, r_mt


def jack(Sm, Lt, nb=20, fn=None):
    N = Sm.shape[0]
    idx = np.array_split(np.arange(N), nb)
    vals = []
    for b in range(nb):
        keep = np.concatenate([idx[j] for j in range(nb) if j != b])
        try:
            v = fn(correlators(Sm[keep], Lt))
        except Exception:
            v = np.nan
        vals.append(v)
    v = np.array(vals, float)
    m = np.nanmean(v)
    e = np.sqrt((nb - 1) / nb * np.nansum((v - m) ** 2))
    return m, e, int(np.sum(np.isfinite(v)))


def w1_of(C, Lt):
    m, A, _ = fit_A(C, Lt)
    return A / C[0], m, A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--Ls', type=int, default=6)
    ap.add_argument('--Lt', type=int, default=12)
    ap.add_argument('--beta', type=float, default=2.5)
    ap.add_argument('--n', type=int, default=3000)
    ap.add_argument('--gap', type=int, default=4)
    ap.add_argument('--therm', type=int, default=300)
    ap.add_argument('--budget', type=float, default=900.0)
    ap.add_argument('--seed', type=int, default=20261125)
    ap.add_argument('--levels', type=str, default='1,2,3,4')
    ap.add_argument('--tag', type=str, default='')
    ap.add_argument('--bench', type=int, default=0, help='benchmark sweeps only, then exit')
    a = ap.parse_args()
    Ls, Lt, beta = a.Ls, a.Lt, a.beta
    levels = tuple(int(x) for x in a.levels.split(','))
    tag = a.tag or f"b{beta}_Ls{Ls}Lt{Lt}"
    P("=" * 78)
    P(f"RUN1  SU(2) Wilson  beta={beta}  Ls={Ls}^3 x Lt={Lt}  target n={a.n} gap={a.gap} "
      f"therm={a.therm} budget={a.budget}s levels={levels} tag={tag}")
    P(f"w1 = A/C(0)  (atom fraction),  Rbar = -C'(0)/C(0) = M1/M0  (first cumulant)")
    P("=" * 78)
    if a.bench:
        plus, minus = index_tables(Ls, Lt)
        U = hot(Ls, Lt, np.random.default_rng(1))
        sweep(U, plus, minus, beta / 2.0, 5)                       # compile
        t = time.time()
        sweep(U, plus, minus, beta / 2.0, a.bench)
        dt = time.time() - t
        P(f"[BENCH] Ls={Ls} Lt={Lt} links={4*Ls**3*Lt} {a.bench} sweeps in {dt:.2f}s "
          f"= {1e3*dt/a.bench:.3f} ms/sweep")
        open(os.path.join(HERE, f'run1_bench_{tag}.log'), 'w').write("\n".join(OUT) + "\n")
        return
    S, plaq, k = markov(Ls, Lt, beta, a.n, a.gap, a.therm, a.budget, a.seed, levels)
    P(f"[E] configs kept: {k}  wall={time.time() - T0:.0f}s  "
      f"plaq={plaq.mean():.5f} +- {plaq.std(ddof=1) / np.sqrt(max(k, 1)):.5f}")
    P(f"[T] t0 autocorr of plaq series: "
      + ", ".join(f"lag{lg}:{np.corrcoef(plaq[:-lg], plaq[lg:])[0, 1]:.3f}"
                  for lg in (1, 2, 4, 8) if k > 4 * lg))
    levels = [0] + list(levels)
    names = {lv: ('raw' if lv == 0 else f'APE{lv}') for lv in levels}
    res = {}
    P("")
    P(f"{'lvl':>5} {'C(0)':>10} {'C1/C0':>8} {'C2/C0':>8} {'m_fit':>7} {'A':>8} "
      f"{'w1=A/C0':>9} {'w1_err':>8} {'w1@1.546':>8} {'SNR(t=2)':>9} {'Rbar':>7} {'<t>':>6}")
    for lv in levels:
        C = correlators(S[:, lv, :], Lt)
        ok = C[0] > 0
        if not ok:
            P(f"{names[lv]:>5} C(0)={C[0]:.3e} NOT positive -> w1 undefined")
            res[names[lv]] = dict(C0=float(C[0]), defined=False)
            continue
        w, m, A = w1_of(C, Lt)
        wm, we, nok = jack(S[:, lv, :], Lt, 20, lambda Cc: w1_of(Cc, Lt)[0])
        # SNR at t=2 (jackknife sigma of C(2))
        _, s2, _ = jack(S[:, lv, :], Lt, 20, lambda Cc: Cc[2])
        snr = C[2] / s2 if s2 > 0 else float('inf')
        r_ln, r_fd, r_mt = rbar(C)
        A154 = float(np.sum((np.exp(-1.546 * np.arange(1, 4)) + np.exp(-1.546 * (Lt - np.arange(1, 4))))
                            * C[1:4]) / np.sum((np.exp(-1.546 * np.arange(1, 4))
                                                + np.exp(-1.546 * (Lt - np.arange(1, 4)))) ** 2))
        P(f"{names[lv]:>5} {C[0]:10.4f} {C[1]/C[0]:8.5f} {C[2]/C[0]:8.5f} {m:7.3f} {A:8.4f} "
          f"{w:9.4f} {we:8.4f} {A154/C[0]:8.4f} {snr:9.1f} {r_ln:7.4f} {r_mt:6.3f}")
        res[names[lv]] = dict(C0=float(C[0]), C=float(C[1] / C[0]),
                              C2=float(C[2] / C[0]), m=m, A=A, w1=w, w1_err=we,
                              w1_m1546=float(A154 / C[0]),
                              snr2=float(snr), rbar=r_ln, rbar_fd=r_fd, tmean=r_mt,
                              defined=True, jack_ok=nok, Cfull=C.tolist())
    # cross-level near-orthogonality at t=2 (operator overlap matrix)
    P("")
    P("[C] overlap |rho| of the level operators on the t=2 slice (off-diag ~ orthogonality):")
    Sm2 = np.stack([S[:, lv, 2] for lv in levels], axis=1)
    R = np.corrcoef(Sm2.T)
    for i, lv in enumerate(levels):
        P("   " + names[lv].rjust(5) + " " + " ".join(f"{R[i, j]:6.3f}" for j in range(len(levels))))
    res['overlap_t2'] = R.tolist()
    np.savez(os.path.join(HERE, f'run1_{tag}.npz'), S=S, plaq=plaq, levels=levels,
             beta=beta, Ls=Ls, Lt=Lt, n=k, gap=a.gap, therm=a.therm, seed=a.seed)
    json.dump(res, open(os.path.join(HERE, f'run1_{tag}.json'), 'w'), indent=1)
    P("")
    P(f"[wall] total {time.time() - T0:.1f}s ; wrote run1_{tag}.npz / .json")
    open(os.path.join(HERE, f'run1_{tag}.log'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
