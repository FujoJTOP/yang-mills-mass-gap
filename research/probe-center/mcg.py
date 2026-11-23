"""Z(3) maximal center gauge fixing, center projection, vortex extraction.

Gauge functional (maximal center gauge, MCG)
--------------------------------------------
    R[g] = sum_{x,mu} |Tr_F U_mu(x)|^2 ,   U -> g U g^dag

This is *identical* to the adjoint Landau-gauge functional, because for
M in SU(3)

    Tr_adj(M) = |Tr_F(M)|^2 - 1 .

So MCG = adjoint Landau gauge (this is the "indirect maximal center gauge"
of the vortex literature, whose functional is the direct MCG functional).

Local maximisation
------------------
For the objective f(g) = sum_a |Tr(g C_a)|^2 (a runs over the 8 staples that
touch the site: 4 outgoing U_mu(x) g^dag(x+mu) and 4 incoming
U_mu(x-mu)^dag g^dag(x-mu)), the Wirtinger gradient is

    d f / d g^dag = sum_a Tr(g C_a) C_a^dag ,

and the unitary maximiser of the *linear* model Re Tr(g G) is the unitary
polar factor of G.  We iterate g <- polar(gradient) with monotone-acceptance
(accept only if f increases).  This is a Riemannian ascent step on U(3); it
is the standard Landau-gauge fixed-point iteration lifted to the quartic
center-gauge functional.  det is corrected to +1.

Projection
----------
Each link U_mu(x) -> z_mu(x) I with z = omega^k, omega = exp(2 pi i/3),
k = round(3 theta / 2 pi) mod 3, theta = arg Tr U_mu(x).  Plaquette flux
k_{mu nu}(x) = k_mu(x) + k_nu(x+mu) - k_mu(x+nu) - k_nu(x) (mod 3).  Nonzero
flux on a lattice plaquette == one vortex plaquette (a dual plaquette).

Vortex clusters
---------------
Two vortex plaquettes belong to the same cluster iff they share at least one
lattice link (i.e. they are edge-adjacent tiles of the same worldsheet).
"""
import numpy as np
from numba import njit


@njit(cache=True, inline='always', fastmath=True)
def _mm(A, B, C):                      # C = A @ B
    for i in range(3):
        for j in range(3):
            s = 0.0 + 0.0j
            for k in range(3):
                s += A[i, k] * B[k, j]
            C[i, j] = s


@njit(cache=True, inline='always', fastmath=True)
def _mul_adj(A, B, C):                 # C = A @ B^dag
    for i in range(3):
        for j in range(3):
            s = 0.0 + 0.0j
            for k in range(3):
                s += A[i, k] * np.conj(B[j, k])
            C[i, j] = s


@njit(cache=True, inline='always', fastmath=True)
def _adj_mul(A, B, C):                 # C = A^dag @ B
    for i in range(3):
        for j in range(3):
            s = 0.0 + 0.0j
            for k in range(3):
                s += np.conj(A[k, i]) * B[k, j]
            C[i, j] = s


@njit(cache=True, inline='always')
def _shift(x0, x1, x2, x3, d, s, L):
    a0 = np.int64(x0); a1 = np.int64(x1)
    a2 = np.int64(x2); a3 = np.int64(x3)
    b = np.int64(s); LL = np.int64(L)
    if d == 0:
        return (a0 + b) % LL, a1, a2, a3
    elif d == 1:
        return a0, (a1 + b) % LL, a2, a3
    elif d == 2:
        return a0, a1, (a2 + b) % LL, a3
    else:
        return a0, a1, a2, (a3 + b) % LL


@njit(cache=True, inline='always', fastmath=True)
def _polar(W, P):
    """Unitary polar factor of W, det -> +1 (same scheme as su3_numba)."""
    tr = 0.0
    for i in range(3):
        for j in range(3):
            tr += (W[i, j] * np.conj(W[i, j])).real
    if tr < 1e-30:
        for i in range(3):
            for j in range(3):
                P[i, j] = 1.0 if i == j else 0.0
        return
    sc = 1.0 / np.sqrt(tr / 3.0)
    Q = np.empty((3, 3), np.complex128)
    Y = np.empty((3, 3), np.complex128)
    Z = np.empty((3, 3), np.complex128)
    for i in range(3):
        for j in range(3):
            P[i, j] = W[i, j] * sc
    for it in range(6):
        for i in range(3):
            for j in range(3):
                s = 0.0 + 0.0j
                for k in range(3):
                    s += np.conj(P[k, i]) * P[k, j]
                Q[i, j] = s
        for i in range(3):
            for j in range(3):
                Y[i, j] = (1.5 if i == j else 0.0) - 0.5 * Q[i, j]
        _mm(P, Y, Z)
        for i in range(3):
            for j in range(3):
                P[i, j] = Z[i, j]
    det = (P[0, 0] * (P[1, 1] * P[2, 2] - P[1, 2] * P[2, 1])
           - P[0, 1] * (P[1, 0] * P[2, 2] - P[1, 2] * P[2, 0])
           + P[0, 2] * (P[1, 0] * P[2, 1] - P[1, 1] * P[2, 0]))
    ph = det / abs(det)
    c = np.conj(ph) ** (1.0 / 3.0)
    for i in range(3):
        for j in range(3):
            P[i, j] = P[i, j] * c


@njit(cache=True)
def _site_f(g, C):
    """f = sum_{a=0..7} |Tr(g C_a)|^2"""
    f = 0.0
    for a in range(8):
        tr = 0.0 + 0.0j
        for i in range(3):
            for j in range(3):
                tr += g[i, j] * C[a][j, i]
        f += (tr.real * tr.real + tr.imag * tr.imag)
    return f


@njit(cache=True)
def mcg_sweep(U, g, L, maxit):
    """One Gauss-Seidel sweep of MCG; updates g in place."""
    C = np.empty((8, 3, 3), np.complex128)
    G = np.empty((3, 3), np.complex128)
    P = np.empty((3, 3), np.complex128)
    gl = np.empty((3, 3), np.complex128)
    for x0 in range(L):
        for x1 in range(L):
            for x2 in range(L):
                for x3 in range(L):
                    for mu in range(4):
                        xp = _shift(x0, x1, x2, x3, mu, 1, L)
                        xm = _shift(x0, x1, x2, x3, mu, -1, L)
                        _mul_adj(U[mu, x0, x1, x2, x3],
                                 g[xp[0], xp[1], xp[2], xp[3]], C[mu])
                        _adj_mul(U[mu, xm[0], xm[1], xm[2], xm[3]],
                                 g[xm[0], xm[1], xm[2], xm[3]], C[4 + mu])
                    for i in range(3):
                        for j in range(3):
                            gl[i, j] = g[x0, x1, x2, x3, i, j]
                    for it in range(maxit):
                        f_old = _site_f(gl, C)
                        for i in range(3):
                            for j in range(3):
                                G[i, j] = 0.0 + 0.0j
                        for a in range(8):
                            tr = 0.0 + 0.0j
                            for i in range(3):
                                for j in range(3):
                                    tr += gl[i, j] * C[a][j, i]
                            for i in range(3):
                                for j in range(3):
                                    G[i, j] += tr * np.conj(C[a][j, i])
                        _polar(G, P)
                        f_new = _site_f(P, C)
                        if f_new > f_old * (1.0 + 1e-12):
                            for i in range(3):
                                for j in range(3):
                                    gl[i, j] = P[i, j]
                        else:
                            break
                    for i in range(3):
                        for j in range(3):
                            g[x0, x1, x2, x3, i, j] = gl[i, j]


@njit(cache=True)
def apply_gauge(U, g, L):
    """Return U' with U'_mu(x) = g(x) U_mu(x) g^dag(x+mu)."""
    V = np.empty_like(U)
    tmp = np.empty((3, 3), np.complex128)
    for mu in range(4):
        for x0 in range(L):
            for x1 in range(L):
                for x2 in range(L):
                    for x3 in range(L):
                        xp = _shift(x0, x1, x2, x3, mu, 1, L)
                        _mm(g[x0, x1, x2, x3], U[mu, x0, x1, x2, x3], tmp)
                        _mul_adj(tmp, g[xp[0], xp[1], xp[2], xp[3]],
                                 V[mu, x0, x1, x2, x3])
    return V


@njit(cache=True)
def center_functional(U, L):
    """R = sum_{x,mu} |Tr U_mu(x)|^2 / (4 L^4);  max possible = 9."""
    tot = 0.0
    for mu in range(4):
        for x0 in range(L):
            for x1 in range(L):
                for x2 in range(L):
                    for x3 in range(L):
                        tr = (U[mu, x0, x1, x2, x3, 0, 0]
                              + U[mu, x0, x1, x2, x3, 1, 1]
                              + U[mu, x0, x1, x2, x3, 2, 2])
                        tot += (tr.real * tr.real + tr.imag * tr.imag)
    return tot / (4.0 * L ** 4)


@njit(cache=True)
def center_project(U, L):
    """Return integer link field n_mu(x) in {0,1,2} (Z(3) log) and the
    fraction of links whose phase is ill-defined (|Tr U| very small)."""
    n = np.zeros((4, L, L, L, L), np.int64)
    amb = 0
    for mu in range(4):
        for x0 in range(L):
            for x1 in range(L):
                for x2 in range(L):
                    for x3 in range(L):
                        re = (U[mu, x0, x1, x2, x3, 0, 0]
                              + U[mu, x0, x1, x2, x3, 1, 1]
                              + U[mu, x0, x1, x2, x3, 2, 2]).real
                        im = (U[mu, x0, x1, x2, x3, 0, 0]
                              + U[mu, x0, x1, x2, x3, 1, 1]
                              + U[mu, x0, x1, x2, x3, 2, 2]).imag
                        mod = np.sqrt(re * re + im * im)
                        if mod < 0.05:
                            amb += 1
                            n[mu, x0, x1, x2, x3] = 0
                            continue
                        th = np.arctan2(im, re)
                        k = np.int64(np.round(3.0 * th / (2.0 * np.pi)))
                        n[mu, x0, x1, x2, x3] = (k % 3)
    return n, amb / (4.0 * L ** 4)


@njit(cache=True)
def plaquette_flux(n, L):
    """k_pmu_nu(x) = n_mu(x)+n_nu(x+mu)-n_mu(x+nu)-n_nu(x) mod 3.
    Returns array (6, L,L,L,L) ordered (mu<nu)."""
    flux = np.zeros((6, L, L, L, L), np.int64)
    idx = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            for x0 in range(L):
                for x1 in range(L):
                    for x2 in range(L):
                        for x3 in range(L):
                            xm = _shift(x0, x1, x2, x3, mu, 1, L)
                            xn = _shift(x0, x1, x2, x3, nu, 1, L)
                            k = (n[mu, x0, x1, x2, x3]
                                 + n[nu, xm[0], xm[1], xm[2], xm[3]]
                                 - n[mu, xn[0], xn[1], xn[2], xn[3]]
                                 - n[nu, x0, x1, x2, x3]) % 3
                            flux[idx, x0, x1, x2, x3] = k
            idx += 1
    return flux


MU_NU = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def _plaq_links(x, mu, nu, L):
    """The four undirected lattice links (as keys) of plaquette (x;mu,nu)."""
    xm = list(x); xn = list(x)
    xm[mu] = (xm[mu] + 1) % L          # x+mu
    xn[nu] = (xn[nu] + 1) % L          # x+nu
    # links: (x,mu), (x+mu,nu), (x+nu,mu), (x,nu); canonical: (y,d) with y_d
    # as the lower endpoint -> use (y,d) directly, undirected so sort endpoints
    def edge(y, d):
        a = tuple(y)
        b = list(y); b[d] = (b[d] + 1) % L
        return (a, tuple(b))
    return [edge(tuple(x), mu), edge(tuple(xm), nu),
            edge(tuple(xn), mu), edge(tuple(x), nu)]


def vortex_clusters(flux, L):
    """Connected components of vortex plaquettes (nonzero flux), adjacency =
    share a lattice link.  Returns (labels list per plaquette, sizes array,
    plaquette list with flux value)."""
    plaqs = []
    f6 = flux
    for idx, (mu, nu) in enumerate(MU_NU):
        nz = np.argwhere(f6[idx] != 0)
        for x in nz:
            plaqs.append((tuple(int(v) for v in x), mu, nu, int(f6[idx][tuple(x)])))
    parent = list(range(len(plaqs)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    by_edge = {}
    for p, (x, mu, nu, k) in enumerate(plaqs):
        for e in _plaq_links(x, mu, nu, L):
            by_edge.setdefault(e, []).append(p)
    for e, lst in by_edge.items():
        for i in range(1, len(lst)):
            union(lst[0], lst[i])
    roots = {}
    labels = np.empty(len(plaqs), np.int64)
    sizes = []
    for p in range(len(plaqs)):
        r = find(p)
        if r not in roots:
            roots[r] = len(roots)
            sizes.append(0)
        labels[p] = roots[r]
        sizes[roots[r]] += 1
    return plaqs, labels, np.array(sizes, np.int64)


def cluster_charge(plaqs, labels, nclust):
    """Per cluster: signed orientation sum  sum (k=1:+1, k=2:-1)."""
    q = np.zeros(nclust, np.int64)
    for p, (x, mu, nu, k) in enumerate(plaqs):
        q[labels[p]] += (1 if k == 1 else -1)
    return q


# ---------------------------------------------------------------------------
# observables
# ---------------------------------------------------------------------------
@njit(cache=True, inline='always')
def _loop_trace(U, x0, x1, x2, x3, mu, nu, R, T, L):
    """Tr(prod loop)/3 for the R (in mu) x T (in nu) rectangle at base x."""
    A = np.empty((3, 3), np.complex128)
    B = np.empty((3, 3), np.complex128)
    for i in range(3):
        for j in range(3):
            A[i, j] = 1.0 if i == j else 0.0
    p0, p1, p2, p3 = x0, x1, x2, x3
    for s in range(R):                                   # +mu
        _mm(A, U[mu, p0, p1, p2, p3], B)
        for i in range(3):
            for j in range(3):
                A[i, j] = B[i, j]
        p0, p1, p2, p3 = _shift(p0, p1, p2, p3, mu, 1, L)
    for s in range(T):                                   # +nu
        _mm(A, U[nu, p0, p1, p2, p3], B)
        for i in range(3):
            for j in range(3):
                A[i, j] = B[i, j]
        p0, p1, p2, p3 = _shift(p0, p1, p2, p3, nu, 1, L)
    for s in range(R):                                   # -mu  (U^dag)
        p0, p1, p2, p3 = _shift(p0, p1, p2, p3, mu, -1, L)
        _mul_adj(A, U[mu, p0, p1, p2, p3], B)
        for i in range(3):
            for j in range(3):
                A[i, j] = B[i, j]
    for s in range(T):                                   # -nu  (U^dag)
        p0, p1, p2, p3 = _shift(p0, p1, p2, p3, nu, -1, L)
        _mul_adj(A, U[nu, p0, p1, p2, p3], B)
        for i in range(3):
            for j in range(3):
                A[i, j] = B[i, j]
    return (A[0, 0] + A[1, 1] + A[2, 2]) / 3.0


@njit(cache=True)
def wilson_loop_tensor(U, L, rmax):
    """Fundamental Wilson loops W[R,T] = <Re Tr/3>, averaged over the 6
    mu<nu planes and all base sites.  R,T = 1..rmax."""
    W = np.zeros((rmax + 1, rmax + 1))
    npl = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            npl += 1
            for R in range(1, rmax + 1):
                for T in range(1, rmax + 1):
                    acc = 0.0
                    for x0 in range(L):
                        for x1 in range(L):
                            for x2 in range(L):
                                for x3 in range(L):
                                    acc += _loop_trace(U, x0, x1, x2, x3,
                                                       mu, nu, R, T, L).real
                    W[R, T] += acc / (L ** 4)
    return W / npl


@njit(cache=True)
def wilson_loop_stats(U, L, rmax):
    """Returns (Tmean, Tabs2): Tmean[R,T] = <Tr/3> (complex),
    Tabs2[R,T] = <|Tr/3|^2> (real), over 6 planes and all base sites.
    Adjoint loop (normalised) = (9*Tabs2 - 1)/8."""
    Tm = np.zeros((rmax + 1, rmax + 1), np.complex128)
    T2 = np.zeros((rmax + 1, rmax + 1))
    npl = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            npl += 1
            for R in range(1, rmax + 1):
                for T in range(1, rmax + 1):
                    acc = 0.0 + 0.0j
                    acc2 = 0.0
                    for x0 in range(L):
                        for x1 in range(L):
                            for x2 in range(L):
                                for x3 in range(L):
                                    t = _loop_trace(U, x0, x1, x2, x3,
                                                    mu, nu, R, T, L)
                                    acc += t
                                    acc2 += (t.real * t.real
                                             + t.imag * t.imag)
                    Tm[R, T] += acc / (L ** 4)
                    T2[R, T] += acc2 / (L ** 4)
    return Tm / npl, T2 / npl


@njit(cache=True)
def expand_center_links(n, L):
    """Turn integer Z(3) log links into SU(3) matrices z I, z=exp(2 pi i k/3)."""
    U = np.zeros((4, L, L, L, L, 3, 3), np.complex128)
    for mu in range(4):
        for x0 in range(L):
            for x1 in range(L):
                for x2 in range(L):
                    for x3 in range(L):
                        k = n[mu, x0, x1, x2, x3] % 3
                        th = 2.0 * np.pi * k / 3.0
                        z = np.cos(th) + 1j * np.sin(th)
                        for i in range(3):
                            U[mu, x0, x1, x2, x3, i, i] = z
    return U


@njit(cache=True)
def polyakov_traces(U, L):
    """Tr(prod_t U_0(x,t)) as complex field (L,L,L) (temporal loop at x0=0)."""
    out = np.zeros((L, L, L), np.complex128)
    A = np.empty((3, 3), np.complex128)
    B = np.empty((3, 3), np.complex128)
    for x1 in range(L):
        for x2 in range(L):
            for x3 in range(L):
                for i in range(3):
                    for j in range(3):
                        A[i, j] = 1.0 if i == j else 0.0
                for t in range(L):
                    _mm(A, U[0, t, x1, x2, x3], B)
                    for i in range(3):
                        for j in range(3):
                            A[i, j] = B[i, j]
                out[x1, x2, x3] = A[0, 0] + A[1, 1] + A[2, 2]
    return out


def polyakov_corr(pf, L, rmax):
    """C(r) averaged over base sites and 3 spatial directions (no subtraction).
    pf is the complex field Tr P(x) (|pf|<=3).  Fundamental normalisation is
    C(r)/9."""
    C = np.zeros(rmax + 1, np.complex128)
    for r in range(rmax + 1):
        acc = 0j
        nn = 0
        for d in range(3):
            acc = acc + (pf * np.conj(np.roll(pf, -r, axis=d))).sum()
            nn += pf.size
        C[r] = acc / nn
    return C


def adjoint_polyakov_corr(pf, L, rmax):
    """Connected adjoint Polyakov-loop correlator.  a(x) = (|pf|^2 - 1)/8,
    C_A(r) = <a(x) a(x+r)> - <a>^2."""
    a = (np.abs(pf) ** 2 - 1.0) / 8.0
    am = a.mean()
    C = np.zeros(rmax + 1)
    for r in range(rmax + 1):
        acc = 0.0
        nn = 0
        for d in range(3):
            acc += (a * np.roll(a, -r, axis=d)).sum()
            nn += a.size
        C[r] = acc / nn - am * am
    return C
