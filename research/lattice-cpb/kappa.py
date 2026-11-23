"""
kappa_eff: lowest eigenvalue of the Wilson-action Hessian restricted to the
gauge-invariant subspace.

Construction
------------
* Link-space coordinates: xi(x,mu,a),  U(x,mu) = exp(i xi_a T^a) U0(x,mu).
* H = d^2 S_W / d xi^2  evaluated by central finite differences (eps=2e-3).
* The gauge orbit is the range of the covariant derivative
      (D w)(x,mu) = w(x) - Ad_{U(x,mu)} w(x+mu),
  which H annihilates exactly (gauge invariance).  P = D (D^dag D + mu0)^-1 D^dag
  projects onto it (mu0 regularises the pseudoinverse; no gauge fixing needed).
* K = (1-P) H (1-P) + Lam * P .  Its lowest eigenvalue is the lowest physical
  curvature (orbit directions are lifted to Lam).
* Lanczos with full reorthogonalisation.

Limitations: eps finite => O(eps^2) truncation; the physical subspace is the
linearised gauge orbit only; the wavelength of a physical mode must fit in L^4.
"""
import numpy as np

import su3 as S

# structure constants [T^a,T^b] = i f^{abc} T^c
_tr3 = np.einsum('aij,bjk,cki->abc', S.T_GEN, S.T_GEN, S.T_GEN)
_F = -2j * (_tr3 - _tr3.transpose(1, 0, 2))
assert np.abs(_F.imag).max() < 1e-12
F_ABC = _F.real.copy()


def adjoint_tables(U):
    T = S.T_GEN
    L = U.shape[1]
    Ad = np.empty((4, L, L, L, L, 8, 8))
    for mu in range(4):
        X = np.einsum('aij,...jk->...aik', T, U[mu])
        Ad[mu] = 2.0 * np.einsum('...aik,bkl,...li->...ab', X, T, S.adj(U[mu])).real
    return Ad


def D_op(w, Ad):
    out = np.empty((4,) + w.shape)
    for mu in range(4):
        out[mu] = w - np.einsum('...ab,...b->...a', Ad[mu],
                                np.roll(w, -1, axis=mu))
    return out


def Dt_op(v, Ad):
    out = np.zeros(v.shape[1:])
    for mu in range(4):
        out += v[mu] - np.einsum('...ba,...b->...a',
                                 np.roll(Ad[mu], 1, axis=mu),
                                 np.roll(v[mu], 1, axis=mu))
    return out


def DtD_op(w, Ad):
    out = np.zeros_like(w)
    for mu in range(4):
        out += 2.0 * w
        out -= np.einsum('...ab,...b->...a', Ad[mu], np.roll(w, -1, axis=mu))
        out -= np.einsum('...ba,...b->...a', np.roll(Ad[mu], 1, axis=mu),
                         np.roll(w, 1, axis=mu))
    return out


def cg(matvec, b, tol=1e-8, maxit=1000, x0=None):
    x = np.zeros_like(b) if x0 is None else x0.copy()
    r = b - matvec(x)
    p = r.copy()
    rs = float(np.sum(r * r))
    bnorm = np.sqrt(float(np.sum(b * b))) + 1e-300
    it = 0
    for it in range(1, maxit + 1):
        Ap = matvec(p)
        denom = float(np.sum(p * Ap))
        if denom <= 0:
            break
        a = rs / denom
        x += a * p
        r -= a * Ap
        rs_new = float(np.sum(r * r))
        if np.sqrt(rs_new) <= tol * bnorm:
            break
        p = r + (rs_new / rs) * p
        rs = rs_new
    return x, np.sqrt(rs_new) / bnorm, it


class HessianFD:
    """Hessian-vector product from central differences of the analytic gradient.

    grad_a(x,mu) = dS_W/dxi_a = (beta/3) Im tr(T^a U(x,mu) S(x,mu)).
    """

    def __init__(self, U0, beta, eps=2e-3):
        self.U0 = U0
        self.beta = beta
        self.eps = eps
        self.shape = (4,) + U0.shape[1:5] + (8,)
        # antisymmetric BCH/connection part of the coordinate Jacobian
        self.sym = -(0.5) * np.einsum('abc,...c->...ab', F_ABC,
                                      S.action_gradient(U0, beta))

    def grad(self, xi):
        U = S.unitaries_from_xi(xi, self.U0)
        return S.action_gradient(U, self.beta)

    def __call__(self, v):
        v = np.asarray(v).reshape(self.shape)
        e = self.eps
        gp = self.grad(e * v)
        gm = self.grad(-e * v)
        # coordinate Jacobian  +  BCH correction => symmetric bi-invariant Hessian
        return ((gp - gm) / (2.0 * e)
                - np.einsum('...ab,...b->...a', self.sym, v)).ravel()


class ProjectedHessian:
    def __init__(self, U0, beta, eps=2e-3, mu0=1e-3, Lam=1e3):
        self.H = HessianFD(U0, beta, eps)
        self.Ad = adjoint_tables(U0)
        self.mu0 = mu0
        self.Lam = Lam
        self.shape = (4,) + U0.shape[1:5] + (8,)
        self.nproj = 0
        self.cg_res = []

    def project(self, v):
        v = np.asarray(v).reshape(self.shape)
        b = Dt_op(v, self.Ad)
        x, res, it = cg(lambda w: DtD_op(w, self.Ad) + self.mu0 * w, b)
        self.nproj += 1
        self.cg_res.append(res)
        return D_op(x, self.Ad)

    def matvec(self, v):
        v = v.reshape(self.shape)
        Pv = self.project(v)
        q = v - Pv
        Hq = self.H(q).reshape(self.shape)
        PHq = self.project(Hq)
        return ((Hq - PHq) + self.Lam * Pv).ravel()


def lanczos_min(matvec, shape, m=60, rng=None, v0=None):
    n = int(np.prod(shape))
    if v0 is None:
        v0 = rng.standard_normal(n) if rng is not None else np.ones(n)
    v = v0 / np.linalg.norm(v0)
    V = np.zeros((m + 1, n))
    V[0] = v
    alpha = np.zeros(m)
    beta = np.zeros(m)
    for j in range(m):
        w = matvec(V[j])
        a = float(V[j] @ w)
        w -= a * V[j]
        if j > 0:
            w -= beta[j - 1] * V[j - 1]
        # full reorthogonalisation (two passes)
        for _ in range(2):
            w -= V[:j + 1].T @ (V[:j + 1] @ w)
        alpha[j] = a
        bn = np.linalg.norm(w)
        if bn < 1e-14:
            m = j + 1
            break
        beta[j] = bn
        V[j + 1] = w / bn
    T = np.diag(alpha[:m])
    if m > 1:
        T += np.diag(beta[:m - 1], 1) + np.diag(beta[:m - 1], -1)
    th, s = np.linalg.eigh(T)
    return th, V[:m], s, m


def lowest_eigenvalue(U0, beta, m=60, seed=0, eps=2e-3, Lam=1e3, mu0=1e-3, return_all=False):
    op = ProjectedHessian(U0, beta, eps=eps, Lam=Lam, mu0=mu0)
    rng = np.random.default_rng(seed)
    th, V, s, used = lanczos_min(op.matvec, op.shape, m=m, rng=rng)
    theta = th[0]
    u = V.T @ s[:, 0]
    r = op.matvec(u) - theta * u
    res = np.linalg.norm(r) / max(abs(theta), 1e-12)
    out = dict(kappa=float(theta), ritz=th[:6], lanczos_res=res, m=used,
               cg_res_max=float(np.max(op.cg_res)) if op.cg_res else 0.0,
               lam=Lam, eps=eps)
    if return_all:
        out['op'] = op
    return out
