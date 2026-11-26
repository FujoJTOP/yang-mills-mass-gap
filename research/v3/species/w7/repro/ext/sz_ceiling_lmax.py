#!/usr/bin/env python3
"""sz_ceiling_lmax.py -- exact lambda_max of the Wilson Hessian at the extremal
configurations (centre vortex for N=2, reflection for N=3), and at controls.

This gives the *exact* supremum over v at a fixed configuration, hence a lower
bound C* >= lambda_max ; the upper bound C* <= 8(d-1) is SZZ Lemma 4.1.
"""
import numpy as np
import itertools, sys, time
import importlib.util

spec = importlib.util.spec_from_file_location("h", "research/v3/w7/sz_ceiling_hess.py")
h = importlib.util.module_from_spec(spec); spec.loader.exec_module(h)


def gf2_solve(A, b, mod):
    rows, cols = A.shape
    Ai = A.copy() % mod; bb = b.copy() % mod
    r = 0; piv = []
    for c in range(cols):
        p = next((rr for rr in range(r, rows) if Ai[rr, c] % mod), None)
        if p is None:
            continue
        Ai[[r, p]] = Ai[[p, r]]; bb[[r, p]] = bb[[p, r]]
        inv = pow(int(Ai[r, c]), -1, mod)
        Ai[r] = (Ai[r] * inv) % mod; bb[r] = (bb[r] * inv) % mod
        for rr in range(rows):
            if rr != r and Ai[rr, c] % mod:
                f = Ai[rr, c] % mod
                Ai[rr] = (Ai[rr] - f * Ai[r]) % mod
                bb[rr] = (bb[rr] - f * bb[r]) % mod
        piv.append(c); r += 1
    for rr in range(r, rows):
        if np.all(Ai[rr] % mod == 0) and bb[rr] % mod:
            return None
    f = np.zeros(cols, dtype=int)
    for i, c in enumerate(piv):
        f[c] = bb[i] % mod
    return f


def vortex_config(t, mod, elt, N):
    A = np.zeros((len(t.plaq), t.nl), dtype=int)
    for pi, (pl, sg) in enumerate(t.plaq):
        for l in pl:
            A[pi, l] = (A[pi, l] + 1) % mod
    f = gf2_solve(A, np.ones(len(t.plaq), dtype=int), mod)
    if f is None:
        return None
    U = [np.linalg.matrix_power(elt, int(fe)) for fe in f]
    ok = True
    for pl, sg in t.plaq:
        M = np.eye(N, dtype=complex)
        for l, s in zip(pl, sg):
            M = M @ (U[l] if s > 0 else U[l].conj().T)
        if not np.allclose(M, elt, atol=1e-9):
            ok = False; break
    return U, ok


def main():
    print("=" * 78)
    print("lambda_max of Hess S/(N beta) at extremal configurations")
    print("=" * 78)
    print("upper bound (SZZ) 8(d-1);  single-link 2(d-1)\n")
    for N in (2, 3):
        T = h.su_basis(N); dim = T.shape[0]
        print(f"--- N={N} ---")
        if N == 2:
            elt = -np.eye(2, dtype=complex); mod = 2
        else:
            elt = np.diag([-1.0, -1.0, 1.0]).astype(complex); mod = 2
        for d in (2, 3, 4):
            L = 4 if d <= 2 else 2
            t = h.Torus(L, d)
            res = vortex_config(t, mod, elt, N)
            if res is None:
                print(f"  d={d}: vortex config not solvable")
                continue
            Uv, ok = res
            lv = np.linalg.eigvalsh(h.hessian_matrix(Uv, t, T))[-1]
            Uf = [np.eye(N, dtype=complex) for _ in range(t.nl)]
            lf = np.linalg.eigvalsh(h.hessian_matrix(Uf, t, T))[-1]
            rng = np.random.default_rng(5)
            lrand = max(np.linalg.eigvalsh(h.hessian_matrix(
                [h.group_from_c(rng.normal(size=dim) * 1.5, T) for _ in range(t.nl)], t, T))[-1]
                for _ in range(4))
            print(f"  d={d} L={L} vortex-ok={ok}:  lambda_max(vortex) = {lv:.6f}   "
                  f"lambda_max(flat) = {lf:.6f}   lambda_max(random) = {lrand:.6f}")
            print(f"        C*/8(d-1) = {lv/(8*(d-1)):.4f}   ceiling 1/(2C*) = {1/(2*lv):.6f}"
                  f"   [SZZ 1/(16(d-1)) = {1/(16*(d-1)):.6f}]")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
