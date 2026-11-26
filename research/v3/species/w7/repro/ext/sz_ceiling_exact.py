#!/usr/bin/env python3
"""sz_ceiling_exact.py -- exact evaluation of C* on the extremal family, and the
proof that SZZ's 8(d-1) is OPTIMAL in d=2.

Upper bound (SZZ Lemma 4.1, and it is the exact v-supremum of the trace bound):
    C* <= 8(d-1) = 4 (links/plaquette) x 2(d-1) (plaquettes/link),
attained as a bound on the constant vector a_e = const  (Perron-Frobenius:
the plaquette-pair matrix N_ij has row sum 8(d-1) and N.1 = 8(d-1).1).

Extremal family used below:
  configuration : maximal centre vortex, every plaquette holonomy = z*1, z in Z_N
                  (solvable over Z_mod for L even; consistency needs L^d even per
                  2-plane, i.e. L even)
  direction     : plane wave X_{x,mu} = b * u^{x_mu} * T,  u = exp(i k_mu),
                  k_mu in {0,pi}, T a fixed algebra direction
At this config the plateau value is  q_p = Re[ z * (sum_i Y_i)^2 ] and, for a
Z_N - aligned Q, the whole form reduces to
    Q(v) = (16/d) * m (d-m)                     [m = #{mu : k_mu = 0}]
for a parallel direction T, normalized |v| = 1.  Checked below.
"""
import numpy as np
import itertools, sys

rng = np.random.default_rng(101)


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
    return np.array(T)


def hess_scalar(U, X, plaq):
    tot = 0.0
    for pl, sg in plaq:
        F = [U[l] if s > 0 else U[l].conj().T for l, s in zip(pl, sg)]
        Y = [X[l] if s > 0 else -X[l] for l, s in zip(pl, sg)]
        n = F[0].shape[0]
        pre = [np.eye(n, dtype=complex)]
        for i in range(4):
            pre.append(pre[i] @ F[i])
        suf = [None] * 5
        suf[4] = np.eye(n, dtype=complex)
        for i in range(3, -1, -1):
            suf[i] = F[i] @ suf[i + 1]
        for i in range(4):
            tot += np.trace(pre[i] @ F[i] @ Y[i] @ Y[i] @ suf[i + 1]).real
        for i in range(4):
            for j in range(i + 1, 4):
                mid = np.eye(n, dtype=complex)
                for k in range(i + 1, j):
                    mid = mid @ F[k]
                tot += 2.0 * np.trace(pre[i] @ F[i] @ Y[i] @ mid @ F[j] @ Y[j] @ suf[j + 1]).real
    return tot


def build(L, d):
    nsites = L ** d
    nl = nsites * d
    plaq = []

    def idx(c):
        x = 0
        for k in reversed(c):
            x = x * L + k
        return x

    def coords(x):
        c = []
        for _ in range(d):
            c.append(x % L); x //= L
        return c

    def link(s, mu):
        return s * d + mu

    for x in range(nsites):
        xs = coords(x)
        for mu in range(d):
            for nu in range(mu + 1, d):
                a = list(xs); a[mu] = (a[mu] + 1) % L
                b = list(xs); b[nu] = (b[nu] + 1) % L
                c = list(xs); c[mu] = (c[mu] + 1) % L; c[nu] = (c[nu] + 1) % L
                plaq.append(((link(x, mu), link(idx(a), nu), link(idx(c), mu), link(idx(b), nu)),
                             (1, 1, -1, -1)))
    return nsites, nl, plaq


def centre_config(L, d, N):
    """solve (coboundary f)_p = 1 over Z_mod, U_e = z^{f_e}."""
    nsites, nl, plaq = build(L, d)
    mod = 2 if N == 2 else 3
    A = np.zeros((len(plaq), nl), dtype=int)
    for pi, (pl, sg) in enumerate(plaq):
        for l in pl:
            A[pi, l] = (A[pi, l] + 1) % mod
    b = np.ones(len(plaq), dtype=int) % mod
    Ai = A.copy(); r = 0; piv = []
    for c in range(nl):
        p = next((rr for rr in range(r, len(plaq)) if Ai[rr, c] % mod), None)
        if p is None:
            continue
        Ai[[r, p]] = Ai[[p, r]]; b[[r, p]] = b[[p, r]]
        inv = pow(int(Ai[r, c]), -1, mod)
        Ai[r] = (Ai[r] * inv) % mod; b[r] = (b[r] * inv) % mod
        for rr in range(len(plaq)):
            if rr != r and Ai[rr, c] % mod:
                f = Ai[rr, c] % mod
                Ai[rr] = (Ai[rr] - f * Ai[r]) % mod
                b[rr] = (b[rr] - f * b[r]) % mod
        piv.append(c); r += 1
    f = np.zeros(nl, dtype=int)
    for i, c in enumerate(piv):
        f[c] = b[i] % mod
    ok = np.all((A @ f) % mod == 1 % mod)
    z = -1.0 if N == 2 else np.exp(2j * np.pi / 3)
    U = [np.linalg.matrix_power(z * np.eye(N), int(fe)) for fe in f]
    return U, ok, plaq, nl


def main():
    print("=" * 78)
    print("Exact C* on the extremal family  (centre vortex + {0,pi} plane wave)")
    print("=" * 78)
    print("prediction:  Q = (16/d) m(d-m)   for m = #{components k_mu = 0}")
    print("upper bound 8(d-1);  single-link family 2(d-1)\n")
    for N in (2, 3):
        T = su_basis(N)
        dim = T.shape[0]
        print(f"--- N={N} ---")
        for d in (2, 3, 4):
            L = 4 if d <= 2 else 2
            U, ok, plaq, nl = centre_config(L, d, N)
            nsites = L ** d
            # direction T0: eigenvalues (a,-a,0,...) -> allows the full trace bound
            lam = np.zeros(N)
            lam[0] = 1.0; lam[1] = -1.0
            T0 = 1j * np.diag(lam) / np.sqrt(np.sum(lam ** 2) / 2.0)
            best = (None, -1e9)
            for m in range(1, d):
                for zeros in itertools.combinations(range(d), m):
                    k = np.full(d, np.pi)
                    for z_ in zeros:
                        k[z_] = 0.0
                    # v: X_{x,mu} = b * exp(i k.x) T0 ; normalize |v|^2 = 1
                    X = []
                    for x in range(nsites):
                        xs = []
                        t = x
                        for _ in range(d):
                            xs.append(t % L); t //= L
                        ph = np.exp(1j * np.dot(k, xs))
                        for mu in range(d):
                            X.append(ph * T0)
                    nv = np.sqrt(sum(-np.trace(x @ x).real for x in X))
                    X = [x / nv for x in X]
                    val = hess_scalar(U, X, plaq)
                    pred = 16.0 * m * (d - m) / d
                    if val > best[1]:
                        best = ((tuple(np.round(k, 3)), m, pred), val)
            print(f"  d={d} L={L} centre-config solvable={ok}: "
                  f"best Q = {best[1]:.6f}   (predicted {16*best[0][1]*(d-best[0][1])/d:.4f} "
                  f"at m={best[0][1]}, k={best[0][0]})   "
                  f"C*/8(d-1) = {best[1]/(8*(d-1)):.4f}")

    print()
    print("=" * 78)
    print("Reference: the two provable bounds and the reading of the ceiling")
    print("=" * 78)
    print("   single-link v, any config        : 2(d-1)")
    print("   constant a_e, any config         : <= 8(d-1)  (Perron-Frobenius, sharp)")
    print("   =>  C* in [2(d-1), 8(d-1)];   ceiling beta_SZZ < 1/(2C*)")
    for d in (2, 3, 4):
        Cs = {"2": 8 * (d - 1)}
        print(f"   d={d}: 1/(2*8(d-1)) = {1/(2*8*(d-1)):.6f}   1/(2*2(d-1)) = {1/(2*2*(d-1)):.6f}")


if __name__ == "__main__":
    main()
