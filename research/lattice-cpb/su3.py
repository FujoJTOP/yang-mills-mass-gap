"""
4D pure SU(3) lattice gauge theory on L^4 with periodic boundary conditions.

Wilson action   S = beta * sum_p (1 - Re tr U_p / 3),   beta = 6/g^2
Update          Cabibbo-Marinari heat bath:
                exact SU(2) heat bath on the three SU(2) subgroups
                (0,1), (0,2), (1,2), sequential scan, checkerboard in parity.

Conventions
-----------
Link variables U[mu][x0,x1,x2,x3] in SU(3).   Plaquette (mu<nu) at base x:
    U_p = U_mu(x) U_nu(x+mu) U_mu^dag(x+nu) U_nu^dag(x)
The single-link weight is exp( (beta/3) Re tr(U S) ), with S the sum of the
6 staples
    S1(nu) = U_nu(x+mu) U_mu^dag(x+nu) U_nu^dag(x)
    S2(nu) = U_nu^dag(x-nu+mu) U_mu^dag(x-nu) U_nu(x-nu)
(verified numerically in validate.py).
"""
import numpy as np

ID2 = np.eye(2, dtype=np.complex128)
ID3 = np.eye(3, dtype=np.complex128)
SIG = np.array([[[0, 1], [1, 0]],
                [[0, -1j], [1j, 0]],
                [[1, 0], [0, -1]]], dtype=np.complex128)


def adj(A):
    return np.conj(np.swapaxes(A, -1, -2))


def gell_mann():
    T = np.zeros((8, 3, 3), dtype=np.complex128)
    T[0] = [[0, 1, 0], [1, 0, 0], [0, 0, 0]]
    T[1] = [[0, -1j, 0], [1j, 0, 0], [0, 0, 0]]
    T[2] = [[1, 0, 0], [0, -1, 0], [0, 0, 0]]
    T[3] = [[0, 0, 1], [0, 0, 0], [1, 0, 0]]
    T[4] = [[0, 0, -1j], [0, 0, 0], [1j, 0, 0]]
    T[5] = [[0, 0, 0], [0, 0, 1], [0, 1, 0]]
    T[6] = [[0, 0, 0], [0, 0, -1j], [0, 1j, 0]]
    T[7] = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) / np.sqrt(3.0)
    return T / 2.0


T_GEN = gell_mann()


# ----------------------------------------------------------------------
# SU(3) helpers
# ----------------------------------------------------------------------
def random_su3(rng, shape=()):
    z = (rng.standard_normal(shape + (3, 3)) +
         1j * rng.standard_normal(shape + (3, 3))) / np.sqrt(2.0)
    q, r = np.linalg.qr(z)
    d = np.diagonal(r, axis1=-2, axis2=-1)
    ph = d / np.abs(d)
    q = q * ph[..., None, :]
    det = np.linalg.det(q)
    q[..., 0, :] *= np.conj(det)[..., None]
    return q


def unitaries_from_xi(xi, U0):
    """U = exp(i xi_a T^a) U0, xi has last axis 8."""
    A = np.einsum('...a,aij->...ij', xi, T_GEN)
    w, V = np.linalg.eigh(A)
    ph = np.exp(1j * w)
    E = np.einsum('...ik,...k,...jk->...ij', V, ph, np.conj(V))
    return E @ U0


def gauge_transform(U, g):
    """U(x,mu) -> g(x) U(x,mu) g(x+mu)^dag ; g (L,L,L,L,3,3)."""
    out = np.empty_like(U)
    for mu in range(4):
        out[mu] = g @ U[mu] @ adj(np.roll(g, -1, axis=mu))
    return out


# ----------------------------------------------------------------------
# staples and action
# ----------------------------------------------------------------------
def staple_sum(U, mu):
    L = U.shape[1]
    Um = U[mu]
    S = np.zeros((L, L, L, L, 3, 3), dtype=np.complex128)
    for nu in range(4):
        if nu == mu:
            continue
        Un = U[nu]
        # S1 = U_nu(x+mu) U_mu^dag(x+nu) U_nu^dag(x)
        S1 = np.roll(Un, -1, axis=mu) @ adj(np.roll(Um, -1, axis=nu)) @ adj(Un)
        # S2 = U_nu^dag(x-nu+mu) U_mu^dag(x-nu) U_nu(x-nu)
        D = np.roll(np.roll(Un, 1, axis=nu), -1, axis=mu)
        E = np.roll(Um, 1, axis=nu)
        F = np.roll(Un, 1, axis=nu)
        S += adj(D) @ adj(E) @ F
        S += S1
    return S


def plaquette_traces(U):
    """dict[(mu,nu)] -> (L,L,L,L) real array  Re tr U_p / 3 for each plaquette."""
    out = {}
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Up = (U[mu] @ np.roll(U[nu], -1, axis=mu)
                  @ adj(np.roll(U[mu], -1, axis=nu)) @ adj(U[nu]))
            out[(mu, nu)] = np.real(np.trace(Up, axis1=-2, axis2=-1)) / 3.0
    return out


def mean_plaquette(U):
    tr = plaquette_traces(U)
    return float(np.mean([a.mean() for a in tr.values()]))


def wilson_action(U, beta):
    tr = plaquette_traces(U)
    return beta * (6.0 * U.shape[1] ** 4 - sum(float(a.sum()) for a in tr.values()))


def magnetic_density(U):
    """s(x) = (1/3) sum_{1<=i<j<=3} Re tr U_ij(x)/3  (spatial plaquettes)."""
    s = np.zeros(U.shape[1:5], dtype=np.float64)
    for i in (1, 2, 3):
        for j in range(i + 1, 4):
            Up = (U[i] @ np.roll(U[j], -1, axis=i)
                  @ adj(np.roll(U[i], -1, axis=j)) @ adj(U[j]))
            s += np.real(np.trace(Up, axis1=-2, axis2=-1)) / 3.0
    return s / 3.0


def rect_loop_trace(U, mu, nu, R, T):
    L = U.shape[1]
    A = np.broadcast_to(ID3, (L, L, L, L, 3, 3))
    for k in range(R):
        A = A @ np.roll(U[mu], -k, axis=mu)
    for l in range(T):
        A = A @ np.roll(np.roll(U[nu], -R, axis=mu), -l, axis=nu)
    for k in range(R - 1, -1, -1):
        A = A @ adj(np.roll(np.roll(U[mu], -k, axis=mu), -T, axis=nu))
    for l in range(T - 1, -1, -1):
        A = A @ adj(np.roll(U[nu], -l, axis=nu))
    return np.real(np.trace(A, axis1=-2, axis2=-1)) / 3.0


def wilson_loop(U, R, T):
    """Average over the 6 orientations (mu<nu) and all base sites."""
    vals = []
    for mu in range(4):
        for nu in range(mu + 1, 4):
            vals.append(rect_loop_trace(U, mu, nu, R, T).mean())
    return float(np.mean(vals))


def polyakov_field(U):
    """Polyakov loop in direction 0, spatial field (L,L,L) complex/3."""
    L = U.shape[1]
    A = np.broadcast_to(ID3, (L, L, L, L, 3, 3))
    for t in range(L):
        A = A @ np.roll(U[0], -t, axis=0)
    return np.trace(A[0], axis1=-2, axis2=-1) / 3.0


# ----------------------------------------------------------------------
# exact SU(2) sampler  (density propto sqrt(1-x^2) exp(lam x) on [-1,1])
# ----------------------------------------------------------------------
def _sample_w0(lam, rng, counters=None):
    lam = np.asarray(lam, dtype=np.float64)
    N = lam.shape[0]
    x = np.empty(N)
    n_prop = 0

    lo = np.nonzero(lam <= 2.0)[0]
    if lo.size:
        act = lo
        while act.size:
            u = rng.random(act.size)
            n_prop += act.size
            p = 2.0 * np.sqrt(u * (1.0 - u)) * np.exp(-2.0 * lam[act] * (1.0 - u))
            acc = rng.random(act.size) < p
            x[act[acc]] = 2.0 * u[acc] - 1.0
            act = act[~acc]

    hi = np.nonzero(lam > 2.0)[0]
    if hi.size:
        act = hi
        while act.size:
            l = lam[act]
            g = rng.standard_normal((act.size, 3))
            n_prop += act.size
            z = np.sum(g * g, axis=1) / (2.0 * l)
            p = np.zeros(act.size)
            ok = z <= 2.0
            p[ok] = np.sqrt((2.0 - z[ok]) / 2.0)
            acc = rng.random(act.size) < p
            x[act[acc]] = 1.0 - z[acc]
            act = act[~acc]

    if counters is not None:
        counters['proposals'] = counters.get('proposals', 0) + n_prop
        counters['samples'] = counters.get('samples', 0) + N
    return x


def su2_heatbath(kappa, M, i, j, rng, counters=None):
    """Exact heat bath for the SU(2) subgroup acting on rows/cols (i,j) of M.

    Returns v (N,2,2) in SU(2) sampled from  exp( kappa * Re tr(v B) ),
    B = M[:, [i,j]][:, :, [i,j]].
    """
    B00 = M[:, i, i]
    B01 = M[:, i, j]
    B10 = M[:, j, i]
    B11 = M[:, j, j]
    b0 = 0.5 * (B00 + B11)
    b1 = 0.5 * (B01 + B10)
    b2 = 0.5j * (B01 - B10)
    b3 = 0.5 * (B00 - B11)
    a = np.stack([b0.real, -b1.imag, -b2.imag, -b3.imag], axis=1)
    an = np.sqrt(np.sum(a * a, axis=1))
    lam = 2.0 * kappa * an
    x0 = _sample_w0(lam, rng, counters)

    g = rng.standard_normal((len(x0), 3))
    g /= np.linalg.norm(g, axis=1, keepdims=True)
    r = np.sqrt(np.maximum(1.0 - x0 * x0, 0.0))
    w = np.stack([x0, r * g[:, 0], r * g[:, 1], r * g[:, 2]], axis=1)

    small = an < 1e-14
    an_safe = np.where(small, 1.0, an)
    ah = a / an_safe[:, None]
    Ahat = (ah[:, 0, None, None] * ID2
            + 1j * (ah[:, 1, None, None] * SIG[0]
                    + ah[:, 2, None, None] * SIG[1]
                    + ah[:, 3, None, None] * SIG[2]))
    W = (w[:, 0, None, None] * ID2
         + 1j * (w[:, 1, None, None] * SIG[0]
                 + w[:, 2, None, None] * SIG[1]
                 + w[:, 3, None, None] * SIG[2]))
    v = Ahat @ W
    if np.any(small):
        v[small] = W[small]
    return v


# ----------------------------------------------------------------------
# lattice
# ----------------------------------------------------------------------
class Lattice:
    def __init__(self, L, beta, rng, cold=True):
        self.L = L
        self.beta = beta
        self.kappa = beta / 3.0
        self.rng = rng
        self.U = self._init(cold)
        # flat index / parity tables
        idx = np.arange(L ** 4).reshape(L, L, L, L)
        x0, x1, x2, x3 = np.meshgrid(*[np.arange(L)] * 4, indexing='ij')
        par = (x0 + x1 + x2 + x3) % 2
        self.idx_par = [np.nonzero(par.reshape(-1) == c)[0] for c in (0, 1)]
        self.sampler_stats = {'proposals': 0, 'samples': 0}

    def _init(self, cold):
        L = self.L
        if cold:
            return np.broadcast_to(ID3, (4, L, L, L, L, 3, 3)).copy()
        U = random_su3(self.rng, (4, L, L, L, L))
        return U

    def sweep(self):
        U = self.U
        Uf = U.reshape(4, -1, 3, 3)
        kappa = self.kappa
        rng = self.rng
        for c in (0, 1):
            idx = self.idx_par[c]
            for mu in range(4):
                S = staple_sum(U, mu).reshape(-1, 3, 3)[idx]
                A = Uf[mu][idx]
                M = A @ S
                for (i, j) in ((0, 1), (0, 2), (1, 2)):
                    v = su2_heatbath(kappa, M, i, j, rng, self.sampler_stats)
                    Mi = M[:, i, :].copy()
                    Mj = M[:, j, :].copy()
                    M[:, i, :] = v[:, 0, 0, None] * Mi + v[:, 0, 1, None] * Mj
                    M[:, j, :] = v[:, 1, 0, None] * Mi + v[:, 1, 1, None] * Mj
                    Ai = A[:, i, :].copy()
                    Aj = A[:, j, :].copy()
                    A[:, i, :] = v[:, 0, 0, None] * Ai + v[:, 0, 1, None] * Aj
                    A[:, j, :] = v[:, 1, 0, None] * Ai + v[:, 1, 1, None] * Aj
                Uf[mu][idx] = A

    def sampler_efficiency(self):
        s = self.sampler_stats
        return s['samples'] / max(s['proposals'], 1)


# ----------------------------------------------------------------------
# cheap observables
# ----------------------------------------------------------------------
def link_hessian_min(U, beta, every=1, rng=None):
    """Per-link 8x8 action Hessian (left algebra basis); returns array of
    lowest eigenvalues (one per link, shape (4,L,L,L,L))."""
    L = U.shape[1]
    out = np.empty((4, L, L, L, L))
    for mu in range(4):
        S = staple_sum(U, mu)
        X = U[mu] @ S  # (L,L,L,L,3,3)
        H = np.einsum('aij,bjk,...ki->...ab', T_GEN, T_GEN, X).real
        H = 0.5 * (H + np.swapaxes(H, -1, -2))
        out[mu] = (beta / 3.0) * np.linalg.eigvalsh(H)[..., 0]
    return out


def action_gradient(U, beta):
    """dS_W/dxi_a for U = exp(i xi_a T^a) U0 (left variation)."""
    G = np.empty((4,) + U.shape[1:5] + (8,))
    for mu in range(4):
        X = U[mu] @ staple_sum(U, mu)
        G[mu] = (beta / 3.0) * np.einsum('aij,...ji->...a', T_GEN, X).imag
    return G


def project_su3(W):
    """Polar projection of a batch of 3x3 complex matrices onto SU(3)."""
    A, s, Bh = np.linalg.svd(W)
    U = A @ Bh
    det = np.linalg.det(U)
    U[..., 0, :] *= np.conj(det)[..., None]
    return U


def ape_smear_spatial(U, alpha=0.5, n=1):
    """APE-smear the three spatial link directions (mu=1,2,3); temporal fixed."""
    U = U.copy()
    for _ in range(n):
        old = U.copy()
        for i in (1, 2, 3):
            S = np.zeros_like(U[i])
            for j in (1, 2, 3):
                if j == i:
                    continue
                S += (old[j] @ np.roll(old[i], -1, axis=j) @ adj(np.roll(old[j], -1, axis=i))
                      + adj(np.roll(old[j], 1, axis=j)) @ np.roll(old[i], 1, axis=j)
                      @ np.roll(old[j], 1, axis=j)) if False else 0
                # forward staple: U_j(x) U_i(x+j) U_j^dag(x+i)
                S += old[j] @ np.roll(old[i], -1, axis=j) @ adj(np.roll(old[j], -1, axis=i))
                # backward staple: U_j^dag(x-j) U_i(x-j) U_j(x-j+i)
                S += (adj(np.roll(old[j], 1, axis=j)) @ np.roll(old[i], 1, axis=j)
                      @ np.roll(np.roll(old[j], 1, axis=j), -1, axis=i))
            U[i] = project_su3((1.0 - alpha) * old[i] + alpha / 6.0 * S)
    return U


def project_su3_cuberoot(W):
    """SU(3) projection by the uniform cube-root phase (matches su3_numba)."""
    A, s, Bh = np.linalg.svd(W)
    U = A @ Bh
    det = np.linalg.det(U)
    ph = det / np.abs(det)
    return U * (np.conj(ph) ** (1.0 / 3.0))[..., None, None]


def ape_smear_spatial_cr(U, alpha=0.5, n=1):
    """APE smearing with the cube-root SU(3) projection (reference for Numba)."""
    U = U.copy()
    for _ in range(n):
        old = U.copy()
        for i in (1, 2, 3):
            S = np.zeros_like(U[i])
            for j in (1, 2, 3):
                if j == i:
                    continue
                S += old[j] @ np.roll(old[i], -1, axis=j) @ adj(np.roll(old[j], -1, axis=i))
                S += (adj(np.roll(old[j], 1, axis=j)) @ np.roll(old[i], 1, axis=j)
                      @ np.roll(np.roll(old[j], 1, axis=j), -1, axis=i))
            U[i] = project_su3_cuberoot((1.0 - alpha) * old[i] + alpha / 6.0 * S)
    return U
