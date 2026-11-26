#!/usr/bin/env python3
"""sz_ceiling_const.py -- SZZ radius-ceiling constant audit (B1 route).

Question: is  |beta| < 1/(16(d-1))  (=1/48 in d=4)  the *ceiling* of the
Shen-Zhu-Zhu (arXiv:2204.12737) Bakry-Emery method, or is the constant loose?

Structure of the SZZ constant (see research/v3/C1-SZZ-NSCALING.md, verified
against the e-print):
    K_S = Ric_coeff - Hess_coeff = N/2 - 8(d-1) N |beta| > 0
        => |beta| < (N/2) / (8(d-1)N) = 1/(16(d-1)).
Gain side  : Ric(v,v) = ((N+2)/2 - 1)|v|^2 = (N/2)|v|^2   -- EXACT (equality).
Cost side  : |Hess_S(v,v)| <= 8(d-1) N |beta| |v|^2      -- Lemma 4.1.

This script decides whether the cost-side constant 8(d-1) can be lowered
*as a bound valid uniformly in v and in the configuration* -- which is what
Bakry-Emery needs.  The answer is NO, and the extremal objects are:

  * extremal v      = the CONSTANT (k=0) mode,
  * extremal config = the maximal center-vortex configuration
                      (every plaquette holonomy = -1 in the centre Z_N).

Parts
  P0  arithmetic: radius, normalisation map beta_W = N^2 beta_SZZ, Ric, lambda_1
  P1  plaquette-pair counting identity  sum_p (sum_{i in p} a_i)^2
  P2  lambda_max of that quadratic form  (= the v-sup of the trace bound)
  P3  the actual Hessian quadratic form on the lattice: exact evaluation
  P4  centre-vortex saturation test (the crux)
  P5  random / optimised search for anything LARGER (upper-bound check)
"""
import numpy as np
import itertools, sys

rng = np.random.default_rng(20261125)
np.set_printoptions(precision=6, suppress=False)

# ----------------------------------------------------------------------
# Lie algebra helpers.  Metric: <X,Y> = -tr(XY);  |v|^2 = -tr(v^2).
# Basis: T_a hermitian, tr(T_a T_b) = delta_ab/2,  X = i c_a T_a,
#        |X|^2 = -tr(X^2) = |c|^2/2.
# ----------------------------------------------------------------------
def su_basis(N):
    T = []
    for j in range(N):
        for k in range(j + 1, N):
            M = np.zeros((N, N), dtype=complex)
            M[j, k] = 0.5; M[k, j] = 0.5
            T.append(M)
            M = np.zeros((N, N), dtype=complex)
            M[j, k] = -0.5j; M[k, j] = 0.5j
            T.append(M)
    for l in range(1, N):
        M = np.zeros((N, N), dtype=complex)
        for j in range(l):
            M[j, j] = 1.0
        M[l, l] = -l
        M *= 1.0 / np.sqrt(2.0 * l * (l + 1.0))
        T.append(M)
    return T


def alg_to_mat(c, T):
    """c (real vector) -> X = i c_a T_a  (anti-hermitian, traceless)."""
    X = np.zeros(T[0].shape, dtype=complex)
    for a, ca in enumerate(c):
        if ca:
            X += ca * T[a]
    return 1j * X


def group_from_alg(c, T):
    """exp(X) with X = i c_a T_a  (unitary, det 1)."""
    X = alg_to_mat(c, T)
    w, V = np.linalg.eigh(-1j * X)          # -iX hermitian  => eigenvalues real
    return (V * np.exp(1j * w)) @ V.conj().T


def unit_direction(dim, rng):
    v = rng.normal(size=dim)
    return v / np.linalg.norm(v)


def norm2(X):
    return -np.trace(X @ X).real


# ----------------------------------------------------------------------
# Lattice
# ----------------------------------------------------------------------
class Torus:
    """d-dimensional periodic lattice, L^d sites."""

    def __init__(self, L, d):
        self.L, self.d = L, d
        self.nsites = L ** d
        self.nl = self.nsites * d          # links
        self.plaq = []                     # list of (4 link indices, 4 signs)
        inc = [0] * d
        for x in range(self.nsites):
            xs = self.coords(x)
            for mu in range(d):
                for nu in range(mu + 1, d):
                    # holonomy = U(x,mu) U(x+mu,nu) U(x+mu+nu,mu)^dag U(x,nu)^dag
                    xs_mu = list(xs); xs_mu[mu] = (xs_mu[mu] + 1) % L
                    xs_nu = list(xs); xs_nu[nu] = (xs_nu[nu] + 1) % L
                    xs_mn = list(xs); xs_mn[mu] = (xs_mn[mu] + 1) % L
                    xs_mn[nu] = (xs_mn[nu] + 1) % L
                    l1 = self.link(x, mu)
                    l2 = self.link(self.idx(xs_mu), nu)
                    l3 = self.link(self.idx(xs_mn), mu)
                    l4 = self.link(self.idx(xs_nu), nu)
                    self.plaq.append(((l1, l2, l3, l4), (1, 1, -1, -1)))
        self.np_ = len(self.plaq)

    def coords(self, x):
        c, L, d = [], self.L, self.d
        for _ in range(d):
            c.append(x % L); x //= L
        return c

    def idx(self, c):
        x = 0
        for k in reversed(c):
            x = x * self.L + k
        return x

    def link(self, site, mu):
        return site * self.d + mu

    def site_of(self, l):
        return l // self.d

    def mu_of(self, l):
        return l % self.d


def plaquette_factors(U, pl, signs):
    """Return F_i and the direction sign sigma_i for the 4 factors."""
    Fs = []
    for l, s in zip(pl, signs):
        Fs.append(U[l] if s > 0 else U[l].conj().T)
    return Fs, signs


def hess_quadratic(U, X, torus):
    """<Hess S, v (x) v> / (N beta)  for S = N beta sum_p Re tr Q_p.

    U : list of N x N unitaries (per link)
    X : list of N x N anti-hermitian directions (per link), |v|^2 = sum -tr X^2
    Returns  sum_p  Re tr( Q_p(t) )'' |_{t=0}.
    Deformations: factor F_i(t) = F_i exp(t Y_i),  Y_i = sigma_i X_link_i
    (for a daggered link, U -> exp(tX)U gives F=U^dag -> F exp(-tX)).
    """
    tot = 0.0
    for pl, sg in torus.plaq:
        F = []
        Y = []
        for l, s in zip(pl, sg):
            F.append(U[l] if s > 0 else U[l].conj().T)
            Y.append(X[l] if s > 0 else -X[l])
        # prefix / suffix products
        pre = [np.eye(F[0].shape[0], dtype=complex)]
        for i in range(4):
            pre.append(pre[i] @ F[i])
        suf = [None] * 5
        suf[4] = np.eye(F[0].shape[0], dtype=complex)
        for i in range(3, -1, -1):
            suf[i] = F[i] @ suf[i + 1]
        val = 0.0
        # diagonal terms: (pre[i] F_i) Y_i^2 (suf[i+1])
        for i in range(4):
            A = pre[i] @ F[i]
            val += np.trace(A @ Y[i] @ Y[i] @ suf[i + 1]).real
        # cross terms 2 Re tr(pre[i]... )
        for i in range(4):
            for j in range(i + 1, 4):
                A = pre[i] @ F[i] @ Y[i]
                mid = np.eye(F[0].shape[0], dtype=complex)
                for k in range(i + 1, j):
                    mid = mid @ F[k]
                A = A @ mid @ F[j] @ Y[j] @ suf[j + 1]
                val += 2.0 * np.trace(A).real
        tot += val
    return tot


def vnorm2(X):
    return sum(norm2(x) for x in X)


# ----------------------------------------------------------------------
# P2:  lambda_max of the plaquette-pair matrix N_ij
#      (upper bound for the v-sup of the per-plaquette trace inequality)
# ----------------------------------------------------------------------
def pair_matrix(torus):
    n = torus.nl
    N = np.zeros((n, n))
    for pl, sg in torus.plaq:
        for a in range(4):
            for b in range(a + 1, 4):
                N[pl[a], pl[b]] += 1.0
                N[pl[b], pl[a]] += 1.0
    for l in range(n):
        N[l, l] = 2 * (torus.d - 1)
    return N


def main():
    print("=" * 78)
    print("P0  arithmetic of the SZZ constant")
    print("=" * 78)
    for d in (2, 3, 4):
        for N in (2, 3):
            Ric = (N + 2) / 2 - 1
            print(f"  d={d} N={N}: Ric_coeff=(N+2)/2-1={Ric:g} (=N/2={N/2:g}?)"
                  f"  8(d-1)={8*(d-1):g}   radius beta_SZZ<1/(16(d-1))={1/(16*(d-1)):.6f}"
                  f"   beta_W=N^2*beta={N*N/(16*(d-1)):.6f}")
    print("  N-independence check:  (N/2)/(8(d-1)N) = 1/(16(d-1))  [N cancels exactly]")
    print("  beta_W = N^2 beta_SZZ  (Wilson normalisation S_W=-(beta_W/N) sum Re tr):")
    print(f"     N=3, d=4: SZZ radius beta_W = {9/48:.6f}  vs head beta_W=2.5  (x{2.5/(9/48):.2f})")
    print(f"     sharp Haar constant lambda_1=(N^2-1)/N : N=2 -> {3/2:.4f}, N=3 -> {8/3:.4f},"
          f"  vs Ric N/2: N=2 -> {1.0:.4f}, N=3 -> {1.5:.4f}")

    print()
    print("=" * 78)
    print("P1/P2  plaquette-pair count: the identity behind 8(d-1)")
    print("=" * 78)
    for d in (2, 3, 4):
        L = 4 if d <= 3 else 2
        t = Torus(L, d)
        Nmat = pair_matrix(t)
        rows = Nmat.sum(axis=1)
        w = np.linalg.eigvalsh(Nmat)
        print(f"  d={d} L={L}: links={t.nl} plaq={t.np_}  "
              f"diag={Nmat[0,0]:.0f}=2(d-1)  rowsum min/max={rows.min():.0f}/{rows.max():.0f}"
              f"=8(d-1)={8*(d-1)}   lam_max(N)={w[-1]:.6f}   lam_min={w[0]:.3e}")
        print(f"        constant vector attains rowsum: ||N.1/1|| check "
              f"{np.allclose(Nmat.sum(axis=1), 8*(d-1))}")
    print("  => 8(d-1) = 4 (links/plaquette) x 2(d-1) (plaquettes/link) is NOT slack")
    print("     in the v-sup: the constant (k=0) mode saturates it exactly.")

    print()
    print("=" * 78)
    print("P3/P4  actual Hessian quadratic form; centre-vortex saturation")
    print("=" * 78)
    for N in (2, 3):
        T = su_basis(N)
        dim = len(T)
        print(f"\n  --- N={N} (dim su(N)={dim}) ---")
        for d in (2, 3, 4):
            L = 4 if d <= 3 else 2
            t = Torus(L, d)
            # --- centre-vortex configuration: every plaquette holonomy = z*1
            #     SU(2): z=-1 ; SU(3): z=exp(2 pi i/3); solve phase f_e in Z_N
            #     with (coboundary f)_p = 1 (mod modulus)
            z = -1.0 if N == 2 else np.exp(2j * np.pi / 3)
            mod = 2 if N == 2 else 3
            # solve linear system over Z_mod: A f = 1
            A = np.zeros((t.np_, t.nl))
            for pi, (pl, sg) in enumerate(t.plaq):
                for l in pl:
                    A[pi, l] += 1.0
            # gaussian elimination mod `mod`
            Ai = A.astype(int) % mod
            b = np.ones(t.np_, dtype=int) % mod
            rows, cols = Ai.shape
            r = 0
            piv = []
            for c in range(cols):
                p = None
                for rr in range(r, rows):
                    if Ai[rr, c] % mod:
                        p = rr; break
                if p is None:
                    continue
                Ai[[r, p]] = Ai[[p, r]]; b[[r, p]] = b[[p, r]]
                inv = pow(int(Ai[r, c]), -1, mod)
                Ai[r] = (Ai[r] * inv) % mod; b[r] = (b[r] * inv) % mod
                for rr in range(rows):
                    if rr != r and Ai[rr, c] % mod:
                        f = Ai[rr, c] % mod
                        Ai[rr] = (Ai[rr] - f * Ai[r]) % mod
                        b[rr] = (b[rr] - f * b[r]) % mod
                piv.append(c); r += 1
                if r == rows:
                    break
            inconsistent = False
            for rr in range(r, rows):
                if np.all(Ai[rr] % mod == 0) and b[rr] % mod:
                    inconsistent = True
            f = np.zeros(cols, dtype=int)
            for i, c in enumerate(piv):
                f[c] = b[i] % mod
            U = [np.linalg.matrix_power(z * np.eye(N), int(fe)) for fe in f] \
                if N == 2 else [(np.exp(2j * np.pi * fe / 3)) * np.eye(N) for fe in f]
            # verify
            hol = []
            for pl, sg in t.plaq:
                M = np.eye(N, dtype=complex)
                for l, s in zip(pl, sg):
                    M = M @ (U[l] if s > 0 else U[l].conj().T)
                hol.append(M)
            phase_ok = all(abs(h[0, 0] - z) < 1e-9 for h in hol)
            # --- uniform direction v: X_l = X0 for all links, |v|^2 = 1
            c0 = unit_direction(dim, rng)
            X0 = alg_to_mat(c0, T)
            # rescale so that |v|^2 = nl * |X0|^2 = 1
            X0 = X0 / np.sqrt(norm2(X0) * t.nl)
            X = [X0 for _ in range(t.nl)]
            val_u = hess_quadratic(U, X, t)
            pred = 8 * (d - 1) * vnorm2(X)
            # --- single-link direction v
            Xs = [np.zeros((N, N), dtype=complex) for _ in range(t.nl)]
            Xs[0] = X0 * np.sqrt(t.nl)          # |v|^2 = 1
            val_s = hess_quadratic(U, Xs, t)
            pred_s = 2 * (d - 1) * vnorm2(Xs)
            print(f"   d={d} L={L}: centre config consistent={phase_ok} "
                  f"(inconsistent={inconsistent})")
            print(f"      uniform v : <Hess>/(N beta) = {val_u: .6f}   "
                  f"8(d-1)|v|^2 = {pred: .6f}   ratio = {val_u/pred:.6f}")
            print(f"      single-lnk: <Hess>/(N beta) = {val_s: .6f}   "
                  f"2(d-1)|v|^2 = {pred_s: .6f}   ratio = {val_s/pred_s:.6f}")

    print()
    print("=" * 78)
    print("P5  search for a LARGER value (upper-bound sanity check)")
    print("=" * 78)
    for N in (2, 3):
        T = su_basis(N); dim = len(T)
        for d in (2, 3, 4):
            L = 2
            t = Torus(L, d)
            best = -1e9; bestcfg = None
            nl = t.nl
            for trial in range(200):
                U = [group_from_alg(rng.normal(size=dim) * 1.4, T) for _ in range(nl)]
                # v: uniform-ish + noise
                X = [alg_to_mat(0.7 * unit_direction(dim, rng)
                                + 0.5 * rng.normal(size=dim), T) for _ in range(nl)]
                nv = vnorm2(X)
                if nv < 1e-12:
                    continue
                X = [x / np.sqrt(nv) for x in X]
                val = hess_quadratic(U, X, t)
                if val > best:
                    best = val; bestcfg = (trial, U, X)
            # local refinement by coordinate ascent on the config (a few sweeps)
            trial, U, X = bestcfg
            val = best
            for sweep in range(30):
                improved = False
                for l in range(nl):
                    for _ in range(3):
                        cnew = rng.normal(size=dim) * 1.0
                        Unew = group_from_alg(cnew, T)
                        old = U[l]; U[l] = Unew
                        v2 = hess_quadratic(U, X, t)
                        if v2 > val + 1e-12:
                            val = v2; improved = True
                        else:
                            U[l] = old
                if not improved:
                    break
            print(f"  N={N} d={d} L=2: max over 200 random + ascent "
                  f"<Hess>/(N beta) = {val:.6f}   vs 8(d-1)={8*(d-1)}   "
                  f"ratio={val/(8*(d-1)):.4f}")

    print()
    print("=" * 78)
    print("P6  consequence for the radius")
    print("=" * 78)
    for N in (2, 3):
        for d in (2, 3, 4):
            for tag, C in (("SZZ 8(d-1)", 8 * (d - 1)), ("single-link 2(d-1)", 2 * (d - 1))):
                bA = 1.0 / (2 * C)
                print(f"  N={N} d={d}  cost coeff C={C:5g} ({tag:20s}) -> "
                      f"beta_SZZ < {bA:.6f}   beta_W = N^2 beta < {N*N*bA:.6f}")
            print()


if __name__ == "__main__":
    main()
