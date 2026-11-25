#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
m1_fast.py -- Numba-accelerated 4D SU(2) Metropolis for the k=0 flat-chart susceptibility.

Purpose: a fast cross-check of the M1 heat-bath plateau (chi_flat = 3.81, 3.76, 3.82 at
L=4,6,8, beta=2.5).  This is a DIFFERENT algorithm (single-link Metropolis with a random
SU(2) proposal) on the SAME measure; it is not a replacement for the heat-bath result.
Validated against L=4/6 before the L=12 run.

Run: py -3.12 m1_fast.py [L] [N] [therm] [gap]
Defaults: L=12 N=500 therm=200 gap=5.  One process; no external modules beyond numpy/numba.
Labels [V]; no QED.  Clay 0%.
"""
import sys, time, json
import numpy as np
from numba import njit

# ---------------------------------------------------------------- 2x2 SU(2) algebra
@njit(cache=True, fastmath=False)
def mmul(A, B):
    C = np.empty((2, 2), np.complex128)
    C[0, 0] = A[0, 0] * B[0, 0] + A[0, 1] * B[1, 0]
    C[0, 1] = A[0, 0] * B[0, 1] + A[0, 1] * B[1, 1]
    C[1, 0] = A[1, 0] * B[0, 0] + A[1, 1] * B[1, 0]
    C[1, 1] = A[1, 0] * B[0, 1] + A[1, 1] * B[1, 1]
    return C

@njit(cache=True, fastmath=False)
def dag(A):
    C = np.empty((2, 2), np.complex128)
    C[0, 0] = np.conj(A[0, 0]); C[0, 1] = np.conj(A[1, 0])
    C[1, 0] = np.conj(A[0, 1]); C[1, 1] = np.conj(A[1, 1])
    return C

@njit(cache=True, inline='always')
def sh(x0, x1, x2, x3, ax, d, L):
    if ax == 0:
        x0 = (x0 + d) % L
    elif ax == 1:
        x1 = (x1 + d) % L
    elif ax == 2:
        x2 = (x2 + d) % L
    else:
        x3 = (x3 + d) % L
    return x0, x1, x2, x3

@njit(cache=True, fastmath=False)
def staple(U, x0, x1, x2, x3, mu, L):
    S = np.zeros((2, 2), np.complex128)
    for nu in range(4):
        if nu == mu:
            continue
        # S1 = U_nu(x+mu) @ U_mu(x+nu)^dag @ U_nu(x)^dag
        a0, a1, a2, a3 = sh(x0, x1, x2, x3, mu, 1, L)
        b0, b1, b2, b3 = sh(x0, x1, x2, x3, nu, 1, L)
        A = U[nu, a0, a1, a2, a3]
        B = U[mu, b0, b1, b2, b3]
        C = U[nu, x0, x1, x2, x3]
        S = S + mmul(mmul(A, dag(B)), dag(C))
        # S2 = U_nu(x-nu+mu)^dag @ U_mu(x-nu)^dag @ U_nu(x-nu)
        c0, c1, c2, c3 = sh(x0, x1, x2, x3, nu, -1, L)
        d0, d1, d2, d3 = sh(c0, c1, c2, c3, mu, 1, L)
        A2 = U[nu, d0, d1, d2, d3]
        B2 = U[mu, c0, c1, c2, c3]
        C2 = U[nu, c0, c1, c2, c3]
        S = S + mmul(mmul(dag(A2), dag(B2)), C2)
    return S

@njit(cache=True, fastmath=False)
def sample_w0(lam):
    if lam <= 2.0:
        while True:
            u = np.random.random()
            p = 2.0 * np.sqrt(u * (1.0 - u)) * np.exp(-2.0 * lam * (1.0 - u))
            if np.random.random() < p:
                return 2.0 * u - 1.0
    while True:
        g0 = np.random.standard_normal()
        g1 = np.random.standard_normal()
        g2 = np.random.standard_normal()
        z = (g0 * g0 + g1 * g1 + g2 * g2) / (2.0 * lam)
        if z <= 2.0:
            p = np.sqrt((2.0 - z) / 2.0)
            if np.random.random() < p:
                return 1.0 - z

@njit(cache=True, fastmath=False)
def unit3():
    n0 = 2.0 * np.random.random() - 1.0
    n1 = 2.0 * np.random.random() - 1.0
    n2 = 2.0 * np.random.random() - 1.0
    nn = np.sqrt(n0 * n0 + n1 * n1 + n2 * n2)
    while nn < 1e-12:
        n0 = 2.0 * np.random.random() - 1.0
        n1 = 2.0 * np.random.random() - 1.0
        n2 = 2.0 * np.random.random() - 1.0
        nn = np.sqrt(n0 * n0 + n1 * n1 + n2 * n2)
    return n0 / nn, n1 / nn, n2 / nn

@njit(cache=True, fastmath=False)
def heatbath_link(U, x0, x1, x2, x3, mu, beta, L):
    S = staple(U, x0, x1, x2, x3, mu, L)
    A = U[mu, x0, x1, x2, x3]
    B = mmul(A, S)
    b00 = B[0, 0]; b01 = B[0, 1]; b10 = B[1, 0]; b11 = B[1, 1]
    a0 = 0.5 * (b00 + b11).real
    a1 = 0.5 * (b01 + b10)
    a2 = 0.5j * (b01 - b10)
    a3 = 0.5 * (b00 - b11)
    av0 = a0; av1 = -a1.imag; av2 = -a2.imag; av3 = -a3.imag
    an = np.sqrt(av0 * av0 + av1 * av1 + av2 * av2 + av3 * av3)
    if an < 1e-14:
        n0, n1, n2, n3 = np.random.standard_normal(), np.random.standard_normal(), np.random.standard_normal(), np.random.standard_normal()
        nn = np.sqrt(n0 * n0 + n1 * n1 + n2 * n2 + n3 * n3)
        w0 = n0 / nn; w1 = n1 / nn; w2 = n2 / nn; w3 = n3 / nn
    else:
        lam = beta * an  # 2*kappa*an, kappa=beta/2
        w0 = sample_w0(lam)
        g0, g1, g2 = unit3()
        r = np.sqrt(max(1.0 - w0 * w0, 0.0))
        w1 = r * g0; w2 = r * g1; w3 = r * g2
        ah0 = av0 / an; ah1 = av1 / an; ah2 = av2 / an; ah3 = av3 / an
        Ahat = np.empty((2, 2), np.complex128)
        Ahat[0, 0] = ah0 + 1j * ah3
        Ahat[0, 1] = 1j * ah1 + ah2
        Ahat[1, 0] = 1j * ah1 - ah2
        Ahat[1, 1] = ah0 - 1j * ah3
        Wm = np.empty((2, 2), np.complex128)
        Wm[0, 0] = w0 + 1j * w3
        Wm[0, 1] = 1j * w1 + w2
        Wm[1, 0] = 1j * w1 - w2
        Wm[1, 1] = w0 - 1j * w3
        v = mmul(Ahat, Wm)
        U[mu, x0, x1, x2, x3] = mmul(v, A)
        return
    Wm = np.empty((2, 2), np.complex128)
    Wm[0, 0] = w0 + 1j * w3
    Wm[0, 1] = 1j * w1 + w2
    Wm[1, 0] = 1j * w1 - w2
    Wm[1, 1] = w0 - 1j * w3
    U[mu, x0, x1, x2, x3] = mmul(Wm, A)

@njit(cache=True, fastmath=False)
def sweep(U, beta, L, alpha=0.0):
    for x0 in range(L):
        for x1 in range(L):
            for x2 in range(L):
                for x3 in range(L):
                    for mu in range(4):
                        heatbath_link(U, x0, x1, x2, x3, mu, beta, L)

@njit(cache=True, fastmath=False)
def chi_flat(U, L):
    Nl = 4 * L ** 4
    A0 = 0.0; A1 = 0.0; A2 = 0.0
    for mu in range(4):
        for x0 in range(L):
            for x1 in range(L):
                for x2 in range(L):
                    for x3 in range(L):
                        a = U[mu, x0, x1, x2, x3, 0, 0]
                        b = U[mu, x0, x1, x2, x3, 0, 1]
                        ar = a.real
                        if ar > 1.0: ar = 1.0
                        if ar < -1.0: ar = -1.0
                        th = np.arccos(ar)
                        sn = np.sin(th)
                        # U = cos(th) I + i sin(th) n.sigma with th = theta/2;
                        # the log in U=exp(i xi.sigma/2) is xi = 2 th n.
                        fac = 2.0 if sn < 1e-7 else 2.0 * th / sn
                        A0 += fac * b.imag
                        A1 += fac * b.real
                        A2 += fac * a.imag
    A0 /= np.sqrt(Nl); A1 /= np.sqrt(Nl); A2 /= np.sqrt(Nl)
    return (A0 * A0 + A1 * A1 + A2 * A2) / 3.0

@njit(cache=True, fastmath=False)
def polyakov_chi(U, L):
    # P(x) = (1/2) Re tr prod_t U_0(t,x); returns <P^2> and <P>^2
    s1 = 0.0; s2 = 0.0; n = L ** 3
    for x1 in range(L):
        for x2 in range(L):
            for x3 in range(L):
                A = np.empty((2, 2), np.complex128)
                A[0, 0] = 1.0 + 0j; A[0, 1] = 0j
                A[1, 0] = 0j; A[1, 1] = 1.0 + 0j
                for t in range(L):
                    A = mmul(A, U[0, t, x1, x2, x3])
                p = 0.5 * A[0, 0].real
                s1 += p * p
                s2 += p
    return s1 / n, (s2 / n) ** 2

@njit(cache=True, fastmath=False)
def cold_start(L):
    U = np.zeros((4, L, L, L, L, 2, 2), np.complex128)
    for mu in range(4):
        for x0 in range(L):
            for x1 in range(L):
                for x2 in range(L):
                    for x3 in range(L):
                        U[mu, x0, x1, x2, x3, 0, 0] = 1.0 + 0j
                        U[mu, x0, x1, x2, x3, 1, 1] = 1.0 + 0j
    return U

def random_start(L, seed):
    rng = np.random.default_rng(seed)
    q = rng.standard_normal((4, L, L, L, L, 4))
    q /= np.linalg.norm(q, axis=-1, keepdims=True)
    U = np.empty((4, L, L, L, L, 2, 2), np.complex128)
    a = q[..., 0]; b1 = q[..., 1]; b2 = q[..., 2]; b3 = q[..., 3]
    U[..., 0, 0] = a + 1j * b3
    U[..., 0, 1] = b2 + 1j * b1
    U[..., 1, 0] = -b2 + 1j * b1
    U[..., 1, 1] = a - 1j * b3
    return U

def run(L, N, therm, gap, alpha=0.35, seed=12345, start='random'):
    np.random.seed(seed)
    U = random_start(L, seed) if start == 'random' else cold_start(L)
    t0 = time.time()
    for _ in range(therm):
        sweep(U, 2.5, L, alpha)
    chis = []
    for i in range(N):
        for _ in range(gap):
            sweep(U, 2.5, L, alpha)
        chis.append(chi_flat(U, L))
    chis = np.array(chis)
    p2, p1sq = polyakov_chi(U, L)
    dt = time.time() - t0
    return dict(L=L, N=N, therm=therm, gap=gap, alpha=alpha,
                chi_mean=float(chis.mean()), chi_se=float(chis.std(ddof=1) / np.sqrt(len(chis))),
                chi_sd=float(chis.std(ddof=1)), P2=float(p2), P1sq=float(p1sq),
                wall_s=float(dt), chi_series=chis.tolist())

if __name__ == '__main__':
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 500
    therm = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    gap = int(sys.argv[4]) if len(sys.argv) > 4 else 5
    out = run(L, N, therm, gap)
    print(json.dumps({k: v for k, v in out.items() if k != 'chi_series'}, indent=1))
    with open('m1_fast_results.json', 'w') as f:
        json.dump(out, f, indent=1)
