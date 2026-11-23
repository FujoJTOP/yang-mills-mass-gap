"""Restriction of the fine Hessian to the b=2 block-uniform subspace.

P_lift: N x Nc (Nc = 4*L'^4*8) sends a coarse field to the fine field that is
constant on each 2^4 block.  The "naive blocked Hessian" is P^T H P.  We report
its lowest eigenvalues.

This is NOT the marginal Hessian (which also has the -Var(dS) term), but it
quantifies how much negativity survives in the coarse sector.
"""
import itertools
import sys

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/cpb-attack/11-slc-blocked')
import su3 as S  # noqa: E402
from full_hessian import build_hessian, flat_idx  # noqa: E402


def lift_matrix(L, b=2):
    D = L // b
    rows, cols = [], []
    Nc = 4 * (D ** 4) * 8
    cidx = lambda mu, xp: (((mu * D + xp[0]) * D + xp[1]) * D + xp[2]) * D + xp[3]
    for mu in range(4):
        for xp in itertools.product(range(D), repeat=4):
            c0 = cidx(mu, xp) * 8
            for a in range(8):
                for d in itertools.product(range(b), repeat=4):
                    x = tuple(xp[k] * b + d[k] for k in range(4))
                    fi = flat_idx(mu, x, L) * 8 + a
                    rows.append(fi); cols.append(c0 + a)
    # value 1 (coarse value -> constant fine value); report P^T H P scaled by
    # 1/sqrt(b^4) so its eigenvalues are log-concavity-like per coarse coord
    vals = np.ones(len(rows)) / np.sqrt(b ** 4)
    return sparse.coo_matrix((vals, (rows, cols)), shape=(4 * L ** 4 * 8, Nc)).tocsr()


def main():
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    beta = 2.5
    rng = np.random.default_rng(13)
    lat = S.Lattice(L, beta, rng, cold=True)
    for _ in range(300):
        lat.sweep()
    P = lift_matrix(L, 2)
    for s in range(2):
        lat.sweep()
        H = build_hessian(lat.U, beta, L)
        Hc = (P.T @ (H @ P)).toarray()
        Hc = 0.5 * (Hc + Hc.T)
        ev = np.sort(eigsh(Hc, k=min(20, Hc.shape[0] - 2), which='SA',
                           return_eigenvectors=False))
        print(f"[cfg {s+1}] coarse Hessian (dim {Hc.shape[0]}) lowest 10: "
              f"{np.round(ev[:10], 4)}")
        print(f"          min={ev[0]:+.4f}  frac_neg(20)={np.mean(ev < 0):.3f}",
              flush=True)


if __name__ == '__main__':
    main()
