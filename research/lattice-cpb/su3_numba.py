"""Numba version of the SU(3) Cabibbo-Marinari heat bath (same algorithm as
su3.py, per-site loop instead of batched NumPy).

Kernel:  sweep(U, L, kappa)  updates U in place, one full heat-bath sweep
(2 parity colours x 4 directions, 3 SU(2) subgroups per link).

The stochastic trajectory is NOT bit-identical to the NumPy version (different
RNG); the deterministic parts (staples, M = U S, the transition probabilities
for a given RNG draw) are the same.  Validation therefore compares
  (a) staples / action / plaquettes bitwise on identical configurations,
  (b) the equilibrium distribution statistically (same <P>, W(R,T)).
"""
import math

import numpy as np
from numba import njit, prange


@njit(cache=True, inline='always', fastmath=True)
def _mm(A, B, C):
    for i in range(3):
        for j in range(3):
            s = 0.0 + 0.0j
            for k in range(3):
                s += A[i, k] * B[k, j]
            C[i, j] = s


@njit(cache=True, inline='always', fastmath=True)
def _mul_adj(A, B, C):          # C = A @ B^dag
    for i in range(3):
        for j in range(3):
            s = 0.0 + 0.0j
            for k in range(3):
                s += A[i, k] * np.conj(B[j, k])
            C[i, j] = s


@njit(cache=True, inline='always', fastmath=True)
def _adj_mul(A, B, C):          # C = A^dag @ B
    for i in range(3):
        for j in range(3):
            s = 0.0 + 0.0j
            for k in range(3):
                s += np.conj(A[k, i]) * B[k, j]
            C[i, j] = s


@njit(cache=True, inline='always')
def _shift(x0, x1, x2, x3, d, s, L):
    a0 = np.int64(x0)
    a1 = np.int64(x1)
    a2 = np.int64(x2)
    a3 = np.int64(x3)
    b = np.int64(s)
    LL = np.int64(L)
    if d == 0:
        return (a0 + b) % LL, a1, a2, a3
    elif d == 1:
        return a0, (a1 + b) % LL, a2, a3
    elif d == 2:
        return a0, a1, (a2 + b) % LL, a3
    else:
        return a0, a1, a2, (a3 + b) % LL


@njit(cache=True)
def _sample_x(lam):
    """density propto sqrt(1-x^2) exp(lam x) on [-1,1]"""
    if lam <= 2.0:
        while True:
            u = np.random.random()
            p = 2.0 * math.sqrt(u * (1.0 - u)) * math.exp(-2.0 * lam * (1.0 - u))
            if np.random.random() < p:
                return 2.0 * u - 1.0
    else:
        while True:
            g1 = np.random.normal()
            g2 = np.random.normal()
            g3 = np.random.normal()
            z = (g1 * g1 + g2 * g2 + g3 * g3) / (2.0 * lam)
            if z <= 2.0:
                if np.random.random() < math.sqrt(0.5 * (2.0 - z)):
                    return 1.0 - z


@njit(cache=True, inline='always', fastmath=True)
def _su2_update(M, A, i, j, kappa):
    """One SU(2) subgroup heat-bath update of rows/cols (i,j) of M and A."""
    b0 = 0.5 * (M[i, i] + M[j, j])
    b1 = 0.5 * (M[i, j] + M[j, i])
    b2 = 0.5j * (M[i, j] - M[j, i])
    b3 = 0.5 * (M[i, i] - M[j, j])
    a0 = b0.real
    a1 = -b1.imag
    a2 = -b2.imag
    a3 = -b3.imag
    an = math.sqrt(a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3)
    lam = 2.0 * kappa * an
    w0 = _sample_x(lam)
    g1 = np.random.normal()
    g2 = np.random.normal()
    g3 = np.random.normal()
    gn = math.sqrt(g1 * g1 + g2 * g2 + g3 * g3)
    r = math.sqrt(max(1.0 - w0 * w0, 0.0))
    w1 = r * g1 / gn
    w2 = r * g2 / gn
    w3 = r * g3 / gn
    if an < 1e-14:
        ah0 = 1.0
        ah1 = ah2 = ah3 = 0.0
    else:
        ah0 = a0 / an
        ah1 = a1 / an
        ah2 = a2 / an
        ah3 = a3 / an
    # Ahat = ah0 I + i(ah1 s1 + ah2 s2 + ah3 s3)
    A00 = ah0 + 1j * ah3
    A01 = 1j * ah1 + ah2
    A10 = 1j * ah1 - ah2
    A11 = ah0 - 1j * ah3
    # W = w0 I + i(w1 s1 + w2 s2 + w3 s3)
    W00 = w0 + 1j * w3
    W01 = 1j * w1 + w2
    W10 = 1j * w1 - w2
    W11 = w0 - 1j * w3
    v00 = A00 * W00 + A01 * W10
    v01 = A00 * W01 + A01 * W11
    v10 = A10 * W00 + A11 * W10
    v11 = A10 * W01 + A11 * W11
    for k in range(3):
        mi = v00 * M[i, k] + v01 * M[j, k]
        mj = v10 * M[i, k] + v11 * M[j, k]
        M[i, k] = mi
        M[j, k] = mj
        ai = v00 * A[i, k] + v01 * A[j, k]
        aj = v10 * A[i, k] + v11 * A[j, k]
        A[i, k] = ai
        A[j, k] = aj


@njit(cache=True, fastmath=True)
def sweep(U, L, kappa):
    S = np.empty((3, 3), np.complex128)
    T1 = np.empty((3, 3), np.complex128)
    T2 = np.empty((3, 3), np.complex128)
    M = np.empty((3, 3), np.complex128)
    A = np.empty((3, 3), np.complex128)
    for c in range(2):
        for mu in range(4):
            for x0 in range(L):
                for x1 in range(L):
                    for x2 in range(L):
                        for x3 in range(L):
                            if ((x0 + x1 + x2 + x3) & 1) != c:
                                continue
                            # ---- staple sum S
                            for i in range(3):
                                for j in range(3):
                                    S[i, j] = 0.0 + 0.0j
                            for nu in range(4):
                                if nu == mu:
                                    continue
                                # U_nu(x+mu), U_mu(x+nu), U_nu(x)
                                y0, y1, y2, y3 = _shift(x0, x1, x2, x3, mu, 1, L)
                                z0, z1, z2, z3 = _shift(x0, x1, x2, x3, nu, 1, L)
                                _mul_adj(U[nu, y0, y1, y2, y3], U[mu, z0, z1, z2, z3], T1)
                                _mul_adj(T1, U[nu, x0, x1, x2, x3], T2)
                                for i in range(3):
                                    for j in range(3):
                                        S[i, j] += T2[i, j]
                                # S2 = U_nu^dag(x-nu+mu) U_mu^dag(x-nu) U_nu(x-nu)
                                w0_, w1_, w2_, w3_ = _shift(x0, x1, x2, x3, nu, -1, L)
                                v0_, v1_, v2_, v3_ = _shift(w0_, w1_, w2_, w3_, mu, 1, L)
                                _mm(U[mu, w0_, w1_, w2_, w3_],
                                    U[nu, v0_, v1_, v2_, v3_], T1)   # U_mu(x-nu) U_nu(x-nu+mu)
                                _adj_mul(T1, U[nu, w0_, w1_, w2_, w3_], T2)
                                for i in range(3):
                                    for j in range(3):
                                        S[i, j] += T2[i, j]
                            # ---- M = A S
                            for i in range(3):
                                for j in range(3):
                                    A[i, j] = U[mu, x0, x1, x2, x3, i, j]
                            _mm(A, S, M)
                            # ---- three SU(2) subgroups
                            for (i, j) in ((0, 1), (0, 2), (1, 2)):
                                _su2_update(M, A, i, j, kappa)
                            for i in range(3):
                                for j in range(3):
                                    U[mu, x0, x1, x2, x3, i, j] = A[i, j]


@njit(cache=True)
def action_numba(U, L, beta):
    tot = 0.0
    P = np.empty((3, 3), np.complex128)
    Bm = np.empty((3, 3), np.complex128)
    C = np.empty((3, 3), np.complex128)
    for x0 in range(L):
        for x1 in range(L):
            for x2 in range(L):
                for x3 in range(L):
                    for mu in range(4):
                        for nu in range(mu + 1, 4):
                            y0, y1, y2, y3 = _shift(x0, x1, x2, x3, mu, 1, L)
                            z0, z1, z2, z3 = _shift(x0, x1, x2, x3, nu, 1, L)
                            _mm(U[mu, x0, x1, x2, x3], U[nu, y0, y1, y2, y3], Bm)
                            _mul_adj(Bm, U[mu, z0, z1, z2, z3], P)
                            _mul_adj(P, U[nu, x0, x1, x2, x3], C)
                            tot += (C[0, 0] + C[1, 1] + C[2, 2]).real
    return beta * (6.0 * L ** 4 - tot / 3.0)


@njit(cache=True)
def plaquette_mean(U, L):
    tot = 0.0
    P = np.empty((3, 3), np.complex128)
    Bm = np.empty((3, 3), np.complex128)
    C = np.empty((3, 3), np.complex128)
    n = 0
    for x0 in range(L):
        for x1 in range(L):
            for x2 in range(L):
                for x3 in range(L):
                    for mu in range(4):
                        for nu in range(mu + 1, 4):
                            y0, y1, y2, y3 = _shift(x0, x1, x2, x3, mu, 1, L)
                            z0, z1, z2, z3 = _shift(x0, x1, x2, x3, nu, 1, L)
                            _mm(U[mu, x0, x1, x2, x3], U[nu, y0, y1, y2, y3], Bm)
                            _mul_adj(Bm, U[mu, z0, z1, z2, z3], P)
                            _mul_adj(P, U[nu, x0, x1, x2, x3], C)
                            tot += (C[0, 0] + C[1, 1] + C[2, 2]).real
                            n += 1
    return tot / (3.0 * n)


@njit(cache=True)
def staple_at(U, x0, x1, x2, x3, mu, L):
    """Staple sum S(x,mu) exactly as used by sweep(); for bitwise comparison."""
    S = np.empty((3, 3), np.complex128)
    T1 = np.empty((3, 3), np.complex128)
    T2 = np.empty((3, 3), np.complex128)
    for i in range(3):
        for j in range(3):
            S[i, j] = 0.0 + 0.0j
    for nu in range(4):
        if nu == mu:
            continue
        y0, y1, y2, y3 = _shift(x0, x1, x2, x3, mu, 1, L)
        z0, z1, z2, z3 = _shift(x0, x1, x2, x3, nu, 1, L)
        _mul_adj(U[nu, y0, y1, y2, y3], U[mu, z0, z1, z2, z3], T1)
        _mul_adj(T1, U[nu, x0, x1, x2, x3], T2)
        for i in range(3):
            for j in range(3):
                S[i, j] += T2[i, j]
        w0_, w1_, w2_, w3_ = _shift(x0, x1, x2, x3, nu, -1, L)
        v0_, v1_, v2_, v3_ = _shift(w0_, w1_, w2_, w3_, mu, 1, L)
        _mm(U[mu, w0_, w1_, w2_, w3_], U[nu, v0_, v1_, v2_, v3_], T1)
        _adj_mul(T1, U[nu, w0_, w1_, w2_, w3_], T2)
        for i in range(3):
            for j in range(3):
                S[i, j] += T2[i, j]
    return S


@njit(cache=True, fastmath=True, parallel=True)
def sweep_par(U, L, kappa):
    """Same heat-bath sweep, parallel over x0 (sites of one parity are
    independent, so no race).  Per-thread work arrays are allocated inside
    the prange body."""
    for c in range(2):
        for mu in range(4):
            for sidx in prange(L * L):
                x0 = sidx // L
                x1g = sidx % L
                S = np.empty((3, 3), np.complex128)
                T1 = np.empty((3, 3), np.complex128)
                T2 = np.empty((3, 3), np.complex128)
                M = np.empty((3, 3), np.complex128)
                A = np.empty((3, 3), np.complex128)
                for x1 in range(x1g, x1g + 1):
                    for x2 in range(L):
                        for x3 in range(L):
                            if ((x0 + x1 + x2 + x3) & 1) != c:
                                continue
                            for i in range(3):
                                for j in range(3):
                                    S[i, j] = 0.0 + 0.0j
                            for nu in range(4):
                                if nu == mu:
                                    continue
                                y0, y1, y2, y3 = _shift(x0, x1, x2, x3, mu, 1, L)
                                z0, z1, z2, z3 = _shift(x0, x1, x2, x3, nu, 1, L)
                                _mul_adj(U[nu, y0, y1, y2, y3], U[mu, z0, z1, z2, z3], T1)
                                _mul_adj(T1, U[nu, x0, x1, x2, x3], T2)
                                for i in range(3):
                                    for j in range(3):
                                        S[i, j] += T2[i, j]
                                w0_, w1_, w2_, w3_ = _shift(x0, x1, x2, x3, nu, -1, L)
                                v0_, v1_, v2_, v3_ = _shift(w0_, w1_, w2_, w3_, mu, 1, L)
                                _mm(U[mu, w0_, w1_, w2_, w3_], U[nu, v0_, v1_, v2_, v3_], T1)
                                _adj_mul(T1, U[nu, w0_, w1_, w2_, w3_], T2)
                                for i in range(3):
                                    for j in range(3):
                                        S[i, j] += T2[i, j]
                            for i in range(3):
                                for j in range(3):
                                    A[i, j] = U[mu, x0, x1, x2, x3, i, j]
                            _mm(A, S, M)
                            for (i, j) in ((0, 1), (0, 2), (1, 2)):
                                _su2_update(M, A, i, j, kappa)
                            for i in range(3):
                                for j in range(3):
                                    U[mu, x0, x1, x2, x3, i, j] = A[i, j]


@njit(cache=True, inline='always', fastmath=True)
def _proj_polar(W, P):
    """Unitary polar factor of W by three Newton iterations, det -> 1."""
    # scale to RMS singular value ~1 so the Newton iteration converges
    tr = 0.0
    for i in range(3):
        for j in range(3):
            tr += (W[i, j] * np.conj(W[i, j])).real
    sc = 1.0 / np.sqrt(tr / 3.0) if tr > 0 else 1.0
    for i in range(3):
        for j in range(3):
            P[i, j] = W[i, j] * sc
    Q = np.empty((3, 3), np.complex128)
    Y = np.empty((3, 3), np.complex128)
    Z = np.empty((3, 3), np.complex128)
    for it in range(6):
        # Q = P^dag P
        for i in range(3):
            for j in range(3):
                s = 0.0 + 0.0j
                for k in range(3):
                    s += np.conj(P[k, i]) * P[k, j]
                Q[i, j] = s
        # Y = 1.5 I - 0.5 Q
        for i in range(3):
            for j in range(3):
                Y[i, j] = (1.5 if i == j else 0.0) - 0.5 * Q[i, j]
        # Z = P Y
        _mm(P, Y, Z)
        for i in range(3):
            for j in range(3):
                P[i, j] = Z[i, j]
    # det
    det = (P[0, 0] * (P[1, 1] * P[2, 2] - P[1, 2] * P[2, 1])
           - P[0, 1] * (P[1, 0] * P[2, 2] - P[1, 2] * P[2, 0])
           + P[0, 2] * (P[1, 0] * P[2, 1] - P[1, 1] * P[2, 0]))
    ph = det / abs(det)
    c = np.conj(ph) ** (1.0 / 3.0)
    for i in range(3):
        for j in range(3):
            P[i, j] = P[i, j] * c


@njit(cache=True, fastmath=True)
def ape_smear_num(U, L, alpha, nsteps):
    """APE smearing of the three spatial directions; returns a new array."""
    V = U.copy()
    W = np.empty((3, 3), np.complex128)
    P = np.empty((3, 3), np.complex128)
    T = np.empty((3, 3), np.complex128)
    S1 = np.empty((3, 3), np.complex128)
    S2 = np.empty((3, 3), np.complex128)
    for it in range(nsteps):
        old = V.copy()
        for i in (1, 2, 3):
            for x0 in range(L):
                for x1 in range(L):
                    for x2 in range(L):
                        for x3 in range(L):
                            for a in range(3):
                                for b in range(3):
                                    W[a, b] = 0.0 + 0.0j
                            for j in (1, 2, 3):
                                if j == i:
                                    continue
                                y0, y1, y2, y3 = _shift(x0, x1, x2, x3, j, 1, L)
                                z0, z1, z2, z3 = _shift(x0, x1, x2, x3, i, 1, L)
                                _mm(old[j, x0, x1, x2, x3], old[i, y0, y1, y2, y3], T)
                                _mul_adj(T, old[j, z0, z1, z2, z3], S1)
                                w0_, w1_, w2_, w3_ = _shift(x0, x1, x2, x3, j, -1, L)
                                v0_, v1_, v2_, v3_ = _shift(w0_, w1_, w2_, w3_, i, 1, L)
                                _mm(old[i, w0_, w1_, w2_, w3_], old[j, v0_, v1_, v2_, v3_], T)
                                _adj_mul(old[j, w0_, w1_, w2_, w3_], T, S2)
                                for a in range(3):
                                    for b in range(3):
                                        W[a, b] += S1[a, b] + S2[a, b]
                            for a in range(3):
                                for b in range(3):
                                    W[a, b] = (1.0 - alpha) * old[i, x0, x1, x2, x3, a, b]                                               + (alpha / 6.0) * W[a, b]
                            _proj_polar(W, P)
                            for a in range(3):
                                for b in range(3):
                                    V[i, x0, x1, x2, x3, a, b] = P[a, b]
    return V


@njit(cache=True, fastmath=True, parallel=True)
def ape_smear_num_par(U, L, alpha, nsteps):
    """APE smearing of the three spatial directions; returns a new array."""
    V = U.copy()
    for it in range(nsteps):
        old = V.copy()
        for i in (1, 2, 3):
            for x0 in prange(L):
                W = np.empty((3, 3), np.complex128)
                P = np.empty((3, 3), np.complex128)
                T = np.empty((3, 3), np.complex128)
                S1 = np.empty((3, 3), np.complex128)
                S2 = np.empty((3, 3), np.complex128)
                for x1 in range(L):
                    for x2 in range(L):
                        for x3 in range(L):
                            for a in range(3):
                                for b in range(3):
                                    W[a, b] = 0.0 + 0.0j
                            for j in (1, 2, 3):
                                if j == i:
                                    continue
                                y0, y1, y2, y3 = _shift(x0, x1, x2, x3, j, 1, L)
                                z0, z1, z2, z3 = _shift(x0, x1, x2, x3, i, 1, L)
                                _mm(old[j, x0, x1, x2, x3], old[i, y0, y1, y2, y3], T)
                                _mul_adj(T, old[j, z0, z1, z2, z3], S1)
                                w0_, w1_, w2_, w3_ = _shift(x0, x1, x2, x3, j, -1, L)
                                v0_, v1_, v2_, v3_ = _shift(w0_, w1_, w2_, w3_, i, 1, L)
                                _mm(old[i, w0_, w1_, w2_, w3_], old[j, v0_, v1_, v2_, v3_], T)
                                _adj_mul(old[j, w0_, w1_, w2_, w3_], T, S2)
                                for a in range(3):
                                    for b in range(3):
                                        W[a, b] += S1[a, b] + S2[a, b]
                            for a in range(3):
                                for b in range(3):
                                    W[a, b] = (1.0 - alpha) * old[i, x0, x1, x2, x3, a, b]                                               + (alpha / 6.0) * W[a, b]
                            _proj_polar(W, P)
                            for a in range(3):
                                for b in range(3):
                                    V[i, x0, x1, x2, x3, a, b] = P[a, b]
    return V




@njit(cache=True, fastmath=True)
def mag_slice(V, L):
    """S(t) = sum over spatial sites of the magnetic plaquette density."""
    out = np.zeros(L, dtype=np.float64)
    Bm = np.empty((3, 3), np.complex128)
    P = np.empty((3, 3), np.complex128)
    C = np.empty((3, 3), np.complex128)
    for x0 in range(L):
        acc = 0.0
        for x1 in range(L):
            for x2 in range(L):
                for x3 in range(L):
                    for i in (1, 2, 3):
                        for j in range(i + 1, 4):
                            y0, y1, y2, y3 = _shift(x0, x1, x2, x3, i, 1, L)
                            z0, z1, z2, z3 = _shift(x0, x1, x2, x3, j, 1, L)
                            _mm(V[i, x0, x1, x2, x3], V[j, y0, y1, y2, y3], Bm)
                            _mul_adj(Bm, V[i, z0, z1, z2, z3], P)
                            _mul_adj(P, V[j, x0, x1, x2, x3], C)
                            acc += (C[0, 0] + C[1, 1] + C[2, 2]).real / 3.0
        out[x0] = acc / 3.0
    return out
