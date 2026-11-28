#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ww1_peak.py -- WW1: SU(2) 4D Wilson MC at L=8, beta in {2.40, 2.45, 2.50}.

Measured object (TT4-WILSON / VV2-MONO conventions):
    qbar_p = 1 - (1/2) Re tr U_p,      <qbar> = 1 - plaq,   plaq = <(1/2) Re tr U_p>
    chi_q  = -(1/6) d<qbar>/d beta     [definition, finite differences]
           = L^4 Var(qbar) = Var(Q)/L^4   [exact exponential-family identity, Q = sum_x q_x]
    Monotone decrease of chi_q on [2.4,2.5]  <=>  <qbar> convex  <=>  D2 > 0,
    D2 = <qbar>(2.40) - 2<qbar>(2.45) + <qbar>(2.50).

Engine: exact SU(2) heat bath, unit-quaternion links, one sweep = 1 measurement.
The 1-d samplers sample_w0 (density ~ sqrt(1-x^2) exp(lam x)) and unit3() are
imported from the existing validated module m1_fast.py in this directory.

Usage (one job at a time; ~3 min per point at N=24k):
    py -3.12 -I ww1_peak.py point 2.40 24000 600 7001
    py -3.12 -I ww1_peak.py point 2.45 24000 600 7002
    py -3.12 -I ww1_peak.py point 2.50 24000 600 7003
    py -3.12 -I ww1_peak.py combine

Outputs: ww1_peak_q<beta>.npz (per-point series), ww1_peak.json (combine report).
Console output is ASCII only (Windows codepage GBK).  Labels describe the data only.
"""
import os
import sys
import json
import time
import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from m1_fast import sample_w0, unit3          # noqa: E402

L = 8
BETAS = (2.40, 2.45, 2.50)
NBLK = 40
DELTA = 1.0 / 24.0                             # TT4-WILSON window
C_BOUND = 1.0 / (6.0 * DELTA)                  # = 4.0000


# ------------------------------------------------------------------ quaternion ops
@njit(cache=True, inline='always')
def qmul(a0, a1, a2, a3, b0, b1, b2, b3):
    return (a0 * b0 - (a1 * b1 + a2 * b2 + a3 * b3),
            a0 * b1 + b0 * a1 - (a2 * b3 - a3 * b2),
            a0 * b2 + b0 * a2 - (a3 * b1 - a1 * b3),
            a0 * b3 + b0 * a3 - (a1 * b2 - a2 * b1))


@njit(cache=True, inline='always')
def getq(U, mu, x0, x1, x2, x3):
    return U[mu, x0, x1, x2, x3, 0], U[mu, x0, x1, x2, x3, 1], \
           U[mu, x0, x1, x2, x3, 2], U[mu, x0, x1, x2, x3, 3]


@njit(cache=True, inline='always')
def putq(U, mu, x0, x1, x2, x3, a0, a1, a2, a3):
    U[mu, x0, x1, x2, x3, 0] = a0
    U[mu, x0, x1, x2, x3, 1] = a1
    U[mu, x0, x1, x2, x3, 2] = a2
    U[mu, x0, x1, x2, x3, 3] = a3


@njit(cache=True, inline='always')
def sh4(x0, x1, x2, x3, ax, d, L):
    if ax == 0:
        x0 = (x0 + d) % L
    elif ax == 1:
        x1 = (x1 + d) % L
    elif ax == 2:
        x2 = (x2 + d) % L
    else:
        x3 = (x3 + d) % L
    return x0, x1, x2, x3


@njit(cache=True)
def staple(U, x0, x1, x2, x3, mu, L):
    s0 = 0.0
    s1 = 0.0
    s2 = 0.0
    s3 = 0.0
    for nu in range(4):
        if nu == mu:
            continue
        y0, y1, y2, y3 = sh4(x0, x1, x2, x3, mu, 1, L)
        A0, A1, A2, A3 = getq(U, nu, y0, y1, y2, y3)
        y0, y1, y2, y3 = sh4(x0, x1, x2, x3, nu, 1, L)
        B0, B1, B2, B3 = getq(U, mu, y0, y1, y2, y3)
        C0, C1, C2, C3 = getq(U, nu, x0, x1, x2, x3)
        u0, u1, u2, u3 = qmul(A0, A1, A2, A3, B0, -B1, -B2, -B3)
        t0, t1, t2, t3 = qmul(u0, u1, u2, u3, C0, -C1, -C2, -C3)
        s0 += t0
        s1 += t1
        s2 += t2
        s3 += t3
        z0, z1, z2, z3 = sh4(x0, x1, x2, x3, nu, -1, L)
        y0, y1, y2, y3 = sh4(z0, z1, z2, z3, mu, 1, L)
        A0, A1, A2, A3 = getq(U, nu, y0, y1, y2, y3)
        B0, B1, B2, B3 = getq(U, mu, z0, z1, z2, z3)
        C0, C1, C2, C3 = getq(U, nu, z0, z1, z2, z3)
        u0, u1, u2, u3 = qmul(-A0, A1, A2, A3, -B0, B1, B2, B3)
        t0, t1, t2, t3 = qmul(u0, u1, u2, u3, C0, C1, C2, C3)
        s0 += t0
        s1 += t1
        s2 += t2
        s3 += t3
    return s0, s1, s2, s3


@njit(cache=True)
def heatbath_link(U, x0, x1, x2, x3, mu, beta, L):
    """Exact: weight exp( beta * (U_mu(x) . S) ),  (a . b) = a0 b0 - a1 b1 - a2 b2 - a3 b3."""
    s0, s1, s2, s3 = staple(U, x0, x1, x2, x3, mu, L)
    a0, a1, a2, a3 = getq(U, mu, x0, x1, x2, x3)
    m0, m1, m2, m3 = qmul(a0, a1, a2, a3, s0, s1, s2, s3)
    mn = np.sqrt(m0 * m0 + m1 * m1 + m2 * m2 + m3 * m3)
    if mn < 1e-13:
        g0 = np.random.standard_normal()
        g1 = np.random.standard_normal()
        g2 = np.random.standard_normal()
        g3 = np.random.standard_normal()
        gn = np.sqrt(g0 * g0 + g1 * g1 + g2 * g2 + g3 * g3)
        putq(U, mu, x0, x1, x2, x3, g0 / gn, g1 / gn, g2 / gn, g3 / gn)
        return
    lam = beta * mn
    w0 = sample_w0(lam)
    g0, g1, g2 = unit3()
    r = np.sqrt(max(1.0 - w0 * w0, 0.0))
    w1 = r * g0
    w2 = r * g1
    w3 = r * g2
    h0 = m0 / mn
    h1 = -m1 / mn
    h2 = -m2 / mn
    h3 = -m3 / mn
    v0, v1, v2, v3 = qmul(h0, h1, h2, h3, w0, w1, w2, w3)
    n0, n1, n2, n3 = qmul(v0, v1, v2, v3, a0, a1, a2, a3)
    putq(U, mu, x0, x1, x2, x3, n0, n1, n2, n3)


@njit(cache=True)
def sweep(U, beta, L):
    """One full sweep: 4 direction passes; links of one direction never share a
    plaquette, so the heat bath is exact for the whole pass."""
    for mu in range(4):
        for x0 in range(L):
            for x1 in range(L):
                for x2 in range(L):
                    for x3 in range(L):
                        heatbath_link(U, x0, x1, x2, x3, mu, beta, L)


@njit(cache=True)
def plaq_total(U, L):
    """sum over 6 planes and all x of (1/2) Re tr U_p  ( = s-part of the quaternion)."""
    tot = 0.0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            for x0 in range(L):
                for x1 in range(L):
                    for x2 in range(L):
                        for x3 in range(L):
                            y0, y1, y2, y3 = sh4(x0, x1, x2, x3, mu, 1, L)
                            A0, A1, A2, A3 = getq(U, mu, x0, x1, x2, x3)
                            B0, B1, B2, B3 = getq(U, nu, y0, y1, y2, y3)
                            y0, y1, y2, y3 = sh4(x0, x1, x2, x3, nu, 1, L)
                            C0, C1, C2, C3 = getq(U, mu, y0, y1, y2, y3)
                            D0, D1, D2, D3 = getq(U, nu, x0, x1, x2, x3)
                            u0, u1, u2, u3 = qmul(A0, A1, A2, A3, B0, B1, B2, B3)
                            v0, v1, v2, v3 = qmul(u0, u1, u2, u3, C0, -C1, -C2, -C3)
                            w0, w1, w2, w3 = qmul(v0, v1, v2, v3, D0, -D1, -D2, -D3)
                            tot += w0
    return tot


@njit(cache=True)
def seed_numba(s):
    np.random.seed(s)


def random_start(L, seed):
    rng = np.random.default_rng(seed)
    q = rng.standard_normal((4, L, L, L, L, 4))
    q /= np.linalg.norm(q, axis=-1, keepdims=True)
    return np.ascontiguousarray(q)


# ------------------------------------------------------------------ statistics
def jk_stats(x, nblk=NBLK):
    """mean and delete-one-block jackknife error."""
    n = x.size
    b = n // nblk
    blk = x[:b * nblk].reshape(nblk, b).mean(axis=1)
    jk = (blk.sum() - blk) / (nblk - 1)
    err = np.sqrt((nblk - 1.0) / nblk * np.sum((jk - jk.mean()) ** 2))
    return float(x[:b * nblk].mean()), float(err)


def jk_var(x, nblk=NBLK):
    """sample variance and its delete-one-block jackknife error."""
    n = x.size
    b = n // nblk
    m = x[:b * nblk].reshape(nblk, b)
    s1 = m.sum(axis=1)
    s2 = (m ** 2).sum(axis=1)
    s1t = s1.sum()
    s2t = s2.sum()
    ntot = b * nblk

    def var_of(a1, a2, n_):
        return (a2 - a1 * a1 / n_) / (n_ - 1.0)

    jk = np.empty(nblk)
    for k in range(nblk):
        jk[k] = var_of(s1t - s1[k], s2t - s2[k], ntot - b)
    v = var_of(s1t, s2t, ntot)
    err = np.sqrt((nblk - 1.0) / nblk * np.sum((jk - jk.mean()) ** 2))
    return float(v), float(err)


# ------------------------------------------------------------------ one point
def run_point(beta, ncfg, therm, seed, tag='WW1'):
    seed_numba(seed)
    np.random.seed(seed)
    U = random_start(L, seed)
    npl = 6.0 * L ** 4
    t0 = time.time()
    for _ in range(therm):
        sweep(U, beta, L)
    q = np.empty(ncfg)
    for i in range(ncfg):
        sweep(U, beta, L)
        q[i] = 1.0 - plaq_total(U, L) / npl
        if (i + 1) % 4000 == 0:
            print("  %s beta=%.2f  %5d/%d  %6.1fs  qbar=%.5f"
                  % (tag, beta, i + 1, ncfg, time.time() - t0, q[:i + 1].mean()),
                  flush=True)
    return q, time.time() - t0


def cmd_point(args):
    beta = float(args[0])
    ncfg = int(args[1]) if len(args) > 1 else 24000
    therm = int(args[2]) if len(args) > 2 else 600
    seed = int(args[3]) if len(args) > 3 else 7001
    q, wall = run_point(beta, ncfg, therm, seed)
    m, e = jk_stats(q)
    v, ve = jk_var(q)
    npl = 6.0 * L ** 4
    chi = npl * v / 6.0
    chie = npl * ve / 6.0
    sd = np.sqrt(v)
    neff = (sd / e) ** 2
    print("POINT beta=%.2f  N=%d  therm=%d  seed=%d  wall=%.1fs" % (beta, ncfg, therm, seed, wall))
    print("  plaq   = %.6f +- %.6f   (qbar mean %.6f +- %.6f)" % (1.0 - m, e, m, e))
    print("  chi_var= L^4 Var(qbar) = %.6f +- %.6f" % (chi, chie))
    print("  Neff   = %.1f   tau_int = %.2f sweeps" % (neff, ncfg / (2.0 * neff)))
    np.savez(os.path.join(HERE, 'ww1_peak_q%.2f.npz' % beta), q=q, beta=beta,
             ncfg=ncfg, therm=therm, seed=seed, L=L, wall=wall,
             qbar=m, qbar_err=e, chi_var=chi, chi_var_err=chie, neff=neff)


# ------------------------------------------------------------------ combine
def cmd_combine():
    npl = 6.0 * L ** 4
    rows = {}
    series = {}
    for b in BETAS:
        f = os.path.join(HERE, 'ww1_peak_q%.2f.npz' % b)
        d = np.load(f)
        q = d['q']
        series[b] = q
        m, e = jk_stats(q)
        v, ve = jk_var(q)
        sd = np.sqrt(v)
        neff = (sd / e) ** 2
        rows[b] = dict(qbar=m, qbar_err=e, plaq=1.0 - m, chi_var=npl * v / 6.0,
                       chi_var_err=npl * ve / 6.0, neff=neff,
                       tau=q.size / (2.0 * neff), ncfg=int(q.size))
        print("beta=%.2f  <qbar>=%.7f +- %.7f  plaq=%.6f  chi_var=%.6f +- %.6f "
              "Neff=%.0f tau=%.1f" % (b, m, e, 1.0 - m, rows[b]['chi_var'],
                                      rows[b]['chi_var_err'], neff, rows[b]['tau']))

    # finite-difference chi_q = -(1/6) d<qbar>/d beta, h = 0.05
    h = 0.5 * (BETAS[2] - BETAS[0])
    m0, e0 = rows[BETAS[0]]['qbar'], rows[BETAS[0]]['qbar_err']
    m1, e1 = rows[BETAS[1]]['qbar'], rows[BETAS[1]]['qbar_err']
    m2, e2 = rows[BETAS[2]]['qbar'], rows[BETAS[2]]['qbar_err']
    fd = {}
    fd[BETAS[0]] = ((m0 - m1) / (6.0 * h), np.hypot(e0, e1) / (6.0 * h))
    fd[BETAS[1]] = ((m0 - m2) / (12.0 * h), np.hypot(e0, e2) / (12.0 * h))
    fd[BETAS[2]] = ((m1 - m2) / (6.0 * h), np.hypot(e1, e2) / (6.0 * h))
    # midpoints of the two intervals
    fd_mid = ((m0 - m1) / (6.0 * h), np.hypot(e0, e1) / (6.0 * h),
              (m1 - m2) / (6.0 * h), np.hypot(e1, e2) / (6.0 * h))

    D2 = m0 - 2.0 * m1 + m2
    D2e = np.sqrt(e0 ** 2 + 4.0 * e1 ** 2 + e2 ** 2)
    print("")
    print("finite differences (-(1/6) d<qbar>/dbeta, h=0.05):")
    for b in BETAS:
        print("  beta=%.2f  chi_fd=%.6f +- %.6f   (chi_var=%.6f +- %.6f)"
              % (b, fd[b][0], fd[b][1], rows[b]['chi_var'], rows[b]['chi_var_err']))
    print("  midpoint 2.425: %.6f +- %.6f" % fd_mid[:2])
    print("  midpoint 2.475: %.6f +- %.6f" % fd_mid[2:])
    drop = fd_mid[0] - fd_mid[2]
    drope = np.hypot(fd_mid[1], fd_mid[3])
    print("  drop 2.425 -> 2.475: %+.6f +- %.6f  = %+.2f sigma"
          % (drop, drope, drop / drope))
    print("convexity D2 = <qbar>(2.40) - 2<qbar>(2.45) + <qbar>(2.50) = %+.3e +- %.3e "
          "(%+.2f sigma; D2>0 <=> chi_q down)" % (D2, D2e, D2 / D2e))

    # quadratic vertex through the three chi_var values (curvature>0 => minimum)
    c0 = rows[BETAS[0]]['chi_var']
    c1 = rows[BETAS[1]]['chi_var']
    c2 = rows[BETAS[2]]['chi_var']
    den = c0 - 2.0 * c1 + c2
    curv = den / (h * h)
    vertex = BETAS[1] + h * (c0 - c2) / (2.0 * den) if abs(den) > 0 else float('nan')
    print("quadratic vertex of chi_var(beta): %.4f  (2a=%+.4e; >0 => minimum)"
          % (vertex, curv))
    argmax_beta = max(BETAS, key=lambda b: rows[b]['chi_var'])

    chi_lo = rows[BETAS[0]]['chi_var']
    print("")
    print("chi_q(2.40)=%.6f  chi_q(2.45)=%.6f  chi_q(2.50)=%.6f   (variance estimator)"
          % (c0, c1, c2))
    print("grid argmax of chi_q at beta=%.2f ; TT4 bound chi_q(2.50) <= %.4f (delta=1/24)"
          % (argmax_beta, C_BOUND))
    print("chi_q(2.50)/C = %.6f" % (c2 / C_BOUND))

    out = dict(L=L, nblk=NBLK, betas=list(BETAS),
               rows=dict((('%.2f' % b), rows[b]) for b in BETAS),
               fd=dict((('%.2f' % b), list(fd[b])) for b in BETAS),
               fd_mid=dict(b2425=list(fd_mid[:2]), b2475=list(fd_mid[2:])),
               drop=drop, drop_err=drope, drop_sigma=drop / drope,
               D2=D2, D2_err=D2e, D2_sigma=D2 / D2e,
               vertex=vertex, curv=curv, argmax_beta=argmax_beta,
               C_bound=C_BOUND, delta=DELTA)
    with open(os.path.join(HERE, 'ww1_peak.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print("wrote ww1_peak.json")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return
    if sys.argv[1] == 'point':
        cmd_point(sys.argv[2:])
    elif sys.argv[1] == 'combine':
        cmd_combine()
    else:
        print(__doc__)


if __name__ == '__main__':
    main()
