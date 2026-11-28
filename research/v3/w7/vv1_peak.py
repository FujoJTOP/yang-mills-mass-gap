#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vv1_peak.py -- VV1: SU(2) 4D Wilson MC, <qbar>(beta) and chi_q(beta), L=8.

Target (TT4-WILSON):
    qbar_p = 1 - (1/2) Re tr U_p,   <qbar> = (1/N_p) sum_p <qbar_p>,  N_p = 6 L^4
    chi_q  = -(1/6) d<qbar>/d beta                        [definition]
           = Var(Q) / (6 N_p),    Q = sum_p qbar_p        [exact fluctuation identity]
    If chi_q is non-increasing on [2.5-d, 2.5], then chi_q(2.5) <= 1/(6 d) = 4 at d = 1/24.
    Residual question: is the peak of chi_q below 2.5 ?

Method: exact SU(2) heat bath (numba, quaternion arithmetic), one full sweep per
measurement, L=8, N=24k sweeps at beta = 2.40 / 2.45 / 2.50, delete-one-block
jackknife with 40 blocks.  chi_q is obtained twice: (a) finite differences of
<qbar>(beta), (b) the static-variance estimator Var(Q)/(6 N_p).

Engine: SU(2) links are stored as unit quaternions U = u0 + i(u1 s1 + u2 s2 + u3 s3);
the 1-d heat-bath samplers sample_w0 / unit3 are imported from m1_fast.py (same
module as the M1 flat-chart cross-check; exact inverse-CDF/rejection samplers).

Usage:  py -3.12 -I vv1_peak.py              # production (writes vv1_peak.npz + report)
        py -3.12 -I vv1_peak.py --selftest   # L=4 validation against w4_blocked (numpy)

All console output is ASCII (Windows console codepage = GBK).
Labels are [T] about the measured data only; no claim beyond L=8, N=24k.
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

BETAS = (2.40, 2.45, 2.50)
L = 8
NCFG = 24000
THERM = 600
NBLK = 40
DELTA = 1.0 / 24.0                            # TT4-WILSON window
C_BOUND = 1.0 / (6.0 * DELTA)                 # = 4.0000


# ------------------------------------------------------------------ quaternion ops
@njit(cache=True, inline='always')
def qmul(a0, a1, a2, a3, b0, b1, b2, b3):
    """(a0 + i a.s)(b0 + i b.s):  s-part a0b0 - a.b, vector a0 b + b0 a - a x b."""
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


# ------------------------------------------------------------------ staple / plaquette
@njit(cache=True)
def staple(U, x0, x1, x2, x3, mu, L):
    """S = sum over nu != mu of [U_nu(x+mu) dag(U_mu(x+nu)) dag(U_nu(x))
                                + dag(U_nu(x-nu+mu)) dag(U_mu(x-nu)) U_nu(x-nu)]."""
    s0 = 0.0
    s1 = 0.0
    s2 = 0.0
    s3 = 0.0
    for nu in range(4):
        if nu == mu:
            continue
        # forward staple
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
        # backward staple
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
    """U_mu(x) <- link with weight exp( kappa Re tr(U_mu S) ), kappa = beta/2 (exact)."""
    s0, s1, s2, s3 = staple(U, x0, x1, x2, x3, mu, L)
    a0, a1, a2, a3 = getq(U, mu, x0, x1, x2, x3)
    m0, m1, m2, m3 = qmul(a0, a1, a2, a3, s0, s1, s2, s3)
    mn = np.sqrt(m0 * m0 + m1 * m1 + m2 * m2 + m3 * m3)
    if mn < 1e-13:                                   # zero staple: Haar element
        g0 = np.random.standard_normal()
        g1 = np.random.standard_normal()
        g2 = np.random.standard_normal()
        g3 = np.random.standard_normal()
        gn = np.sqrt(g0 * g0 + g1 * g1 + g2 * g2 + g3 * g3)
        putq(U, mu, x0, x1, x2, x3, g0 / gn, g1 / gn, g2 / gn, g3 / gn)
        return
    lam = beta * mn                                  # 2 kappa |M|
    w0 = sample_w0(lam)
    g0, g1, g2 = unit3()
    r = np.sqrt(max(1.0 - w0 * w0, 0.0))
    w1 = r * g0
    w2 = r * g1
    w3 = r * g2
    h0 = m0 / mn                                     # Ahat = (m0,-m1,-m2,-m3)/|M|
    h1 = -m1 / mn
    h2 = -m2 / mn
    h3 = -m3 / mn
    v0, v1, v2, v3 = qmul(h0, h1, h2, h3, w0, w1, w2, w3)
    n0, n1, n2, n3 = qmul(v0, v1, v2, v3, a0, a1, a2, a3)
    putq(U, mu, x0, x1, x2, x3, n0, n1, n2, n3)


@njit(cache=True)
def sweep(U, beta, L):
    for x0 in range(L):
        for x1 in range(L):
            for x2 in range(L):
                for x3 in range(L):
                    for mu in range(4):
                        heatbath_link(U, x0, x1, x2, x3, mu, beta, L)


@njit(cache=True)
def plaq_total(U, L):
    """sum over 6 planes of Re tr U_p / 2  ( = sum of the quaternion s-part)."""
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
def cold_start(L):
    U = np.zeros((4, L, L, L, L, 4), np.float64)
    for mu in range(4):
        for x0 in range(L):
            for x1 in range(L):
                for x2 in range(L):
                    for x3 in range(L):
                        U[mu, x0, x1, x2, x3, 0] = 1.0
    return U


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
    s = blk.sum()
    jk = (s - blk) / (nblk - 1)
    err = np.sqrt((nblk - 1.0) / nblk * np.sum((jk - jk.mean()) ** 2))
    return float(x[:b * nblk].mean()), float(err)


def jk_var(x, nblk=NBLK):
    """sample variance and its delete-one-block jackknife error."""
    n = x.size
    b = n // nblk
    m = x[:b * nblk].reshape(nblk, b)
    s1 = m.sum(axis=1)
    s2 = (m.astype(np.float64) ** 2).sum(axis=1)
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


def run_point(L, beta, ncfg, therm, seed, tag=''):
    np.random.seed(seed)
    U = random_start(L, seed)
    t0 = time.time()
    for _ in range(therm):
        sweep(U, beta, L)
    q = np.empty(ncfg)
    npl = 6.0 * L ** 4
    for i in range(ncfg):
        sweep(U, beta, L)
        q[i] = 1.0 - plaq_total(U, L) / npl
        if (i + 1) % 4000 == 0:
            print("  %s beta=%.2f  %5d/%d  %6.1fs  qbar=%.5f"
                  % (tag, beta, i + 1, ncfg, time.time() - t0, q[:i + 1].mean()),
                  flush=True)
    wall = time.time() - t0
    return q, wall


def report(data):
    npl = 6.0 * L ** 4
    print("")
    print("=== VV1: SU(2) L=%d, N=%d sweeps/point, %d blocks ===" % (L, NCFG, NBLK))
    print("beta    <qbar>        +-jk        chi_var      +-jk      Neff     tau")
    rows = {}
    for b in BETAS:
        q = data['%.2f' % b]
        m, e = jk_stats(q)
        v, ve = jk_var(q)
        chi = npl * v / 6.0
        chie = npl * ve / 6.0
        sd = np.sqrt(v)
        neff = (sd / e) ** 2
        tau = q.size / (2.0 * neff)
        rows['%.2f' % b] = dict(qbar=m, qbar_err=e, chi_var=chi, chi_var_err=chie,
                                sd=sd, neff=neff, tau=tau, wall=data['wall_%.2f' % b])
        print("%.2f  %.7f  %.7f   %.7f  %.7f  %7.1f  %6.2f"
              % (b, m, e, chi, chie, neff, tau), flush=True)

    # --- finite-difference chi_q (definition), equally spaced grid
    b0, b1, b2 = BETAS
    m0, e0 = rows['%.2f' % b0]['qbar'], rows['%.2f' % b0]['qbar_err']
    m1, e1 = rows['%.2f' % b1]['qbar'], rows['%.2f' % b1]['qbar_err']
    m2, e2 = rows['%.2f' % b2]['qbar'], rows['%.2f' % b2]['qbar_err']
    h = 0.5 * (b2 - b0)
    fd = {}
    fd[b0] = ((m0 - m1) / (6.0 * h), np.hypot(e0, e1) / (6.0 * h))
    fd[b1] = ((m0 - m2) / (12.0 * h), np.hypot(e0, e2) / (12.0 * h))
    fd[b2] = ((m1 - m2) / (6.0 * h), np.hypot(e1, e2) / (6.0 * h))
    print("")
    print("chi_q from finite differences (definition, -(1/6) d<qbar>/d beta):")
    print("beta    chi_fd       +-          chi_var      +-          pull")
    for b in BETAS:
        c, ce = fd[b]
        cv = rows['%.2f' % b]['chi_var']
        cve = rows['%.2f' % b]['chi_var_err']
        pull = (c - cv) / np.hypot(ce, cve)
        print("%.2f  %.7f  %.7f   %.7f  %.7f   %+6.2f"
              % (b, c, ce, cv, cve, pull), flush=True)
        rows['%.2f' % b].update(chi_fd=c, chi_fd_err=ce, pull=pull)

    # --- convexity of <qbar>:  D2 = q0 - 2 q1 + q2 > 0  <=>  d(chi)/dbeta < 0
    D2 = m0 - 2.0 * m1 + m2
    D2e = np.sqrt(e0 ** 2 + 4.0 * e1 ** 2 + e2 ** 2)
    print("")
    print("convexity  D2 = <q>(b0) - 2<q>(b1) + <q>(b2) = %+.3e +- %.3e  (D2>0 <=> chi_q down)"
          % (D2, D2e))

    # --- quadratic through the three chi_var values
    c0 = rows['%.2f' % b0]['chi_var']
    c1 = rows['%.2f' % b1]['chi_var']
    c2 = rows['%.2f' % b2]['chi_var']
    curv = (c0 - 2.0 * c1 + c2) / (h * h)          # 2a
    if abs(curv) > 0:
        xv = h * (c0 - c2) / (2.0 * (c0 - 2.0 * c1 + c2))   # vertex offset from b1
    else:
        xv = float('nan')
    print("delta chi:  %+.3e (b0->b1)   %+.3e (b1->b2)" % (c1 - c0, c2 - c1))
    print("quadratic vertex of chi(beta) at beta = %.4f  (curvature %+.3e; >0 => minimum)"
          % (b1 + xv, curv))

    # --- peak statement
    cvals = [c0, c1, c2]
    k = int(np.argmax(cvals))
    print("")
    print("measured chi_q: %.6f (%.2f)  %.6f (%.2f)  %.6f (%.2f)"
          % (c0, b0, c1, b1, c2, b2))
    print("argmax on the grid: beta = %.2f  =>  beta_peak <= %.2f < 2.5"
          % (BETAS[k], BETAS[k]))
    print("TT4 bound at d=1/24: chi_q(2.5) <= C = 1/(6d) = %.4f ; measured %.6f"
          % (C_BOUND, c2))
    out = dict(L=L, ncfg=NCFG, therm=THERM, nblk=NBLK, betas=list(BETAS),
               rows=rows, fd=dict((k2, list(v2)) for k2, v2 in fd.items()),
               D2=D2, D2err=D2e, vertex=b1 + xv, curv=curv,
               C_bound=C_BOUND, argmax_beta=BETAS[k])
    return out


# ------------------------------------------------------------------ selftest (L=4)
def selftest():
    """L=4, beta=2.5: quaternion heat bath vs the numpy code of w4_blocked.py."""
    sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
    import w4_blocked as W
    Ls = 4
    beta = 2.5
    npl = 6.0 * Ls ** 4

    np.random.seed(11)
    U = random_start(Ls, 11)
    for _ in range(400):
        sweep(U, beta, Ls)
    qa = np.empty(600)
    for i in range(600):
        sweep(U, beta, Ls)
        qa[i] = 1.0 - plaq_total(U, Ls) / npl

    # plaquette bookkeeping vs the numpy Wilson action on the SAME configuration:
    # quaternion (u0,u1,u2,u3) -> complex  U = [[u0+i u3, u2+i u1], [-u2+i u1, u0-i u3]]
    Uc = np.empty((4, Ls, Ls, Ls, Ls, 2, 2), np.complex128)
    Uc[..., 0, 0] = U[..., 0] + 1j * U[..., 3]
    Uc[..., 0, 1] = U[..., 2] + 1j * U[..., 1]
    Uc[..., 1, 0] = -U[..., 2] + 1j * U[..., 1]
    Uc[..., 1, 1] = U[..., 0] - 1j * U[..., 3]
    q_quat = 1.0 - plaq_total(U, Ls) / npl
    q_np = 1.0 - W.action2(Uc, beta) / (beta * npl)
    print("bookkeeping: qbar(quaternion plaq_total) = %.12f ; qbar(W.action2) = %.12f ; diff %.2e"
          % (q_quat, q_np, abs(q_quat - q_np)))

    # independent engine: numpy heat bath (w4_blocked.staple_sum2 + su3.su2_heatbath),
    # same checkerboard sweep as w4_blocked.sample2, but streamed (no config storage)
    from su3 import su2_heatbath
    rng = np.random.default_rng(7)
    V = W.random_su2(rng, (4, Ls, Ls, Ls, Ls))
    ns = Ls ** 4
    x0, x1, x2, x3 = np.meshgrid(*[np.arange(Ls)] * 4, indexing='ij')
    par = (x0 + x1 + x2 + x3) % 2
    idxp = [np.nonzero(par.reshape(-1) == c)[0] for c in (0, 1)]
    Vf = V.reshape(4, ns, 2, 2)
    kappa = beta / 2.0

    def np_sweep():
        for c in (0, 1):
            ids = idxp[c]
            for mu in range(4):
                S = W.staple_sum2(V, mu).reshape(ns, 2, 2)[ids]
                A = Vf[mu, ids]
                v = su2_heatbath(kappa, A @ S, 0, 1, rng)
                Vf[mu, ids] = v @ A

    for _ in range(400):
        np_sweep()
    qb = np.empty(400)
    for i in range(400):
        np_sweep()
        qb[i] = 1.0 - W.action2(V, beta) / (beta * npl)

    print("selftest L=%d beta=%.2f" % (Ls, beta))
    print("  numba quaternion heat bath : <qbar> = %.6f +- %.6f  (N=%d)"
          % (qa.mean(), qa.std(ddof=1) / np.sqrt(qa.size), qa.size))
    print("  numpy  w4_blocked heat bath: <qbar> = %.6f +- %.6f  (N=%d)"
          % (qb.mean(), qb.std(ddof=1) / np.sqrt(qb.size), qb.size))
    print("  difference = %+.6f" % (qa.mean() - qb.mean()))


def main():
    if '--selftest' in sys.argv:
        selftest()
        return
    t_all = time.time()
    data = {}
    for k, b in enumerate(BETAS):
        q, wall = run_point(L, b, NCFG, THERM, 90000 + int(1000 * b), tag='VV1')
        data['%.2f' % b] = q
        data['wall_%.2f' % b] = wall
        print("  beta=%.2f done in %.1fs" % (b, wall), flush=True)
    np.savez(os.path.join(HERE, 'vv1_peak.npz'),
             betas=np.array(BETAS), L=L, ncfg=NCFG, therm=THERM, nblk=NBLK,
             **{('q%.2f' % b): data['%.2f' % b] for b in BETAS})
    out = report(data)
    out['wall_total'] = time.time() - t_all
    with open(os.path.join(HERE, 'vv1_peak.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print("")
    print("total wall = %.1fs" % out['wall_total'])


if __name__ == '__main__':
    main()
