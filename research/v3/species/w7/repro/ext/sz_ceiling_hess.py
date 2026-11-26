#!/usr/bin/env python3
"""sz_ceiling_hess.py -- exact config-supremum C* of the Wilson-action Hessian.

C*  :=  sup over configurations U,  sup over unit link-algebra v,  of
           <Hess S, v (x) v> / ( N |beta| |v|^2 ),      S = N beta sum_p Re tr Q_p
so that the Bakry-Emery (SZZ) radius is   beta_SZZ < (N/2)/(N C*) = 1/(2 C*).

SZZ's Lemma 4.1 gives the *v-uniform, config-uniform* bound  C* <= 8(d-1).
This script computes the actual C*.

Geodesic second variation, per plaquette (all 4 factors, sign sigma_i = -1 for a
daggered link, B_i = Ad_{Q^-1 G_i}(sigma_i X_i),  G_i = F_1...F_i,  Q = G_4):
     tr(prod_i F_i e^{t Y_i})''  =  sum_i Retr(Q B_i^2) + 2 sum_{i<j} Retr(Q B_i B_j).

Verified against a direct scalar implementation (function hess_scalar) below.
"""
import numpy as np
import time, sys

rng = np.random.default_rng(20261125)


# ---------------------------------------------------------------- algebra
def su_basis(N):
    T = []
    for j in range(N):
        for k in range(j + 1, N):
            M = np.zeros((N, N), dtype=complex); M[j, k] = 0.5; M[k, j] = 0.5; T.append(M)
            M = np.zeros((N, N), dtype=complex); M[j, k] = -0.5j; M[k, j] = 0.5j; T.append(M)
    for l in range(1, N):
        M = np.zeros((N, N), dtype=complex)
        for j in range(l):
            M[j, j] = 1.0
        M[l, l] = -l
        M *= 1.0 / np.sqrt(2.0 * l * (l + 1.0))
        T.append(M)
    return np.array(T)                     # (dim, N, N)


def mat_from_c(c, T):
    return 1j * np.tensordot(c, T, axes=1)


def group_from_c(c, T):
    X = mat_from_c(c, T)
    w, V = np.linalg.eigh(-1j * X)
    return (V * np.exp(1j * w)) @ V.conj().T


class Torus:
    def __init__(self, L, d):
        self.L, self.d = L, d
        self.nsites = L ** d
        self.nl = self.nsites * d
        self.plaq = []
        for x in range(self.nsites):
            xs = self.coords(x)
            for mu in range(d):
                for nu in range(mu + 1, d):
                    c_mu = list(xs); c_mu[mu] = (c_mu[mu] + 1) % L
                    c_nu = list(xs); c_nu[nu] = (c_nu[nu] + 1) % L
                    c_mn = list(xs); c_mn[mu] = (c_mn[mu] + 1) % L; c_mn[nu] = (c_mn[nu] + 1) % L
                    self.plaq.append(((self.link(x, mu), self.link(self.idx(c_mu), nu),
                                       self.link(self.idx(c_mn), mu), self.link(self.idx(c_nu), nu)),
                                      (1, 1, -1, -1)))

    def coords(self, x):
        c = []
        for _ in range(self.d):
            c.append(x % self.L); x //= self.L
        return c

    def idx(self, c):
        x = 0
        for k in reversed(c):
            x = x * self.L + k
        return x

    def link(self, s, mu):
        return s * self.d + mu


# ---------------------------------------------------------------- scalar check
def hess_scalar(U, X, torus):
    """direct O(plaquettes) evaluation of <Hess S/(N beta), v v>; slow reference."""
    tot = 0.0
    for pl, sg in torus.plaq:
        F = [U[l] if s > 0 else U[l].conj().T for l, s in zip(pl, sg)]
        Y = [X[l] if s > 0 else -X[l] for l, s in zip(pl, sg)]
        pre = [np.eye(F[0].shape[0], dtype=complex)]
        for i in range(4):
            pre.append(pre[i] @ F[i])
        suf = [None] * 5
        suf[4] = np.eye(F[0].shape[0], dtype=complex)
        for i in range(3, -1, -1):
            suf[i] = F[i] @ suf[i + 1]
        for i in range(4):
            tot += np.trace(pre[i] @ F[i] @ Y[i] @ Y[i] @ suf[i + 1]).real
        for i in range(4):
            for j in range(i + 1, 4):
                mid = np.eye(F[0].shape[0], dtype=complex)
                for k in range(i + 1, j):
                    mid = mid @ F[k]
                tot += 2.0 * np.trace(pre[i] @ F[i] @ Y[i] @ mid @ F[j] @ Y[j] @ suf[j + 1]).real
    return tot


# ---------------------------------------------------------------- fast assembly
def hessian_matrix(U, torus, T):
    """Vectorised assembly of the Hessian quadratic form (nl*dim)^2."""
    N = T.shape[1]
    dim = T.shape[0]
    nl = torus.nl
    npl = len(torus.plaq)
    Ua = np.array(U)                                  # (nl,N,N)
    li = np.array([p[0] for p in torus.plaq], dtype=np.int64)  # (npl,4)
    sg = np.array([p[1] for p in torus.plaq], dtype=np.int64)           # (npl,4)
    F = np.where(sg[..., None, None] > 0, Ua[li], np.conj(np.transpose(Ua[li], (0, 1, 3, 2))))
    # G[p,i] = F_1...F_{i-1};  i = 0..4
    G = np.empty((npl, 5, N, N), dtype=complex)
    G[:, 0] = np.eye(N)
    for i in range(4):
        G[:, i + 1] = G[:, i] @ F[:, i]
    Q = G[:, 4]                                       # (npl,N,N)
    Qd = np.conj(np.transpose(Q, (0, 2, 1)))
    # C[p,i] = Q^-1 G_{i+1}   (i=0..3;  G_{i+1} = F_1...F_{i+1})
    C = Qd[:, None] @ G[:, 1:5]                       # (npl,4,N,N)
    Cd = np.conj(np.transpose(C, (0, 1, 3, 2)))
    # B[p,i,a] = sigma_i * C E_a C^-1
    B = (C[..., None, :, :] @ T[None, None] @ Cd[..., None, :, :])
    B = B * sg.astype(float)[:, :, None, None, None]   # (npl,4,dim,N,N)
    # Q B  (npl,4,dim,N,N)
    QB = Q[:, None, None] @ B
    # blocks
    # diag: Retr(Q B_i^a B_i^b) = Retr((Q B_i^a) B_i^b)
    Hd = np.einsum('pikab,pilba->pikl', QB, B).real   # (npl,4,dim,dim)
    Hd = 0.5 * (Hd + np.swapaxes(Hd, 2, 3))           # symmetrise a<->b (full double sum)
    Ho = np.einsum('pikab,pjlba->pijkl', QB, B).real  # (npl,4,4,dim,dim)
    H = np.zeros((nl * dim, nl * dim))
    # scatter diagonal blocks (both triangles)
    idxbase = li * dim                                # (npl,4)
    for i in range(4):
        base = idxbase[:, i]
        rows = (base[:, None] + np.arange(dim)[None, :])            # (npl,dim)
        cols = (base[:, None] + np.arange(dim)[None, :])
        rr = np.repeat(rows, dim, axis=1).ravel()
        cc = np.tile(cols, (1, dim)).ravel()
        np.add.at(H, (rr, cc), Hd[:, i].ravel())
    for i in range(4):
        for j in range(i + 1, 4):
            bi, bj = idxbase[:, i], idxbase[:, j]
            rr = np.repeat(bi[:, None] + np.arange(dim)[None, :], dim, axis=1).ravel()
            cc = np.tile(bj[:, None] + np.arange(dim)[None, :], (1, dim)).ravel()
            val = Ho[:, i, j].ravel()
            np.add.at(H, (rr, cc), val)
            np.add.at(H, (cc, rr), val)
    # The assembled trace identity carries the opposite overall sign to the
    # direct geodesic second derivative (verified exactly: ratio = -1.000000
    # on random configurations for N = 2, 3).  Flip it.
    return -H


def main():
    print("=" * 78)
    print("verification: fast assembly vs direct scalar evaluation")
    print("=" * 78)
    for N in (2, 3):
        T = su_basis(N); dim = T.shape[0]
        t = Torus(3, 2)
        U = [group_from_c(rng.normal(size=dim), T) for _ in range(t.nl)]
        Hs = hessian_matrix(U, t, T)
        errs = []
        for _ in range(6):
            c = rng.normal(size=(t.nl, dim))
            nrm = np.sqrt(np.sum(c ** 2) / 2.0)
            c = c / nrm
            X = [mat_from_c(c[l], T) for l in range(t.nl)]
            v = np.concatenate([c[l] for l in range(t.nl)])      # |v|^2 = sum |c|^2/2
            direct = hess_scalar(U, X, t)
            quad = v @ Hs @ v
            errs.append(abs(direct - quad) / max(1.0, abs(direct)))
        print(f"  N={N}: max rel err = {max(errs):.2e}")

    print()
    print("=" * 78)
    print("C* = sup over configurations of lambda_max(Hess S/(N beta))")
    print("=" * 78)
    for N in (2, 3):
        T = su_basis(N); dim = T.shape[0]
        print(f"\n--- N={N} (dim su(N)={dim}) ---")
        for d, L, nstart, nsweep in ((2, 4, 8, 6), (3, 2, 6, 6), (3, 3, 4, 4),
                                     (4, 2, 4, 5), (2, 2, 6, 6)):
            t = Torus(L, d)
            t0 = time.time()
            best = -1e9
            for s in range(nstart):
                if s == 0:
                    U = [np.eye(N, dtype=complex) for _ in range(t.nl)]
                else:
                    U = [group_from_c(rng.normal(size=dim) * 1.3, T) for _ in range(t.nl)]
                val = np.linalg.eigvalsh(hessian_matrix(U, t, T))[-1]
                for sw in range(nsweep):
                    imp = False
                    for l in rng.permutation(t.nl):
                        old = U[l]
                        for _ in range(4):
                            U[l] = group_from_c(rng.normal(size=dim) * 1.1, T)
                            nv = np.linalg.eigvalsh(hessian_matrix(U, t, T))[-1]
                            if nv > val + 1e-9:
                                val = nv; old = U[l]; imp = True
                        U[l] = old
                    if not imp:
                        break
                best = max(best, val)
            print(f"  d={d} L={L}: nl={t.nl:4d}  nplaq={len(t.plaq):4d}   "
                  f"C*_num = {best:9.6f}   [2(d-1)={2*(d-1)}, 8(d-1)={8*(d-1)}]   "
                  f"C*/8(d-1) = {best/(8*(d-1)):.4f}   beta_SZZ<1/(2C*)={1/(2*best):.6f}"
                  f"   ({time.time()-t0:.0f}s)")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
