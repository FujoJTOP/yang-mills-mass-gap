"""Ensemble-average coarse-projected Hessian Hbar = E[P^T H P].

The marginal blocked Hessian is G''(A') = E[d2 S | A'] - Var(dS | A').  The
E-part's ensemble average is Hbar = E[P^T H P].  Along any FIXED coarse
direction v:
    E_{A'}[ v^T G''(A') v ] = v^T Hbar v - E_{A'}[ Var(v^T dS | A') ] <= v^T Hbar v,
since Var >= 0 only subtracts.  Hence if Hbar has a NEGATIVE eigenvalue, the
marginal Hessian is negative somewhere -> blocked-SLC is DEAD.

If min(Hbar) >= 0 this is inconclusive (the Var term could still push it
negative), so it is a one-sided decisive test (dead / inconclusive).
"""
import sys

import numpy as np
from scipy.sparse.linalg import eigsh  # noqa: F401

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/cpb-attack/11-slc-blocked')
import su3 as S  # noqa: E402
from full_hessian import build_hessian  # noqa: E402
from coarse_block_hessian import lift_matrix  # noqa: E402


def main():
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    ncfg = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    beta = 2.5
    rng = np.random.default_rng(21)
    lat = S.Lattice(L, beta, rng, cold=True)
    for _ in range(300):
        lat.sweep()
    P = lift_matrix(L, 2)

    Hbar = None
    permin = []
    for c in range(ncfg):
        lat.sweep()
        H = build_hessian(lat.U, beta, L)
        Hc = (P.T @ (H @ P)).toarray()
        Hc = 0.5 * (Hc + Hc.T)
        permin.append(float(np.linalg.eigvalsh(Hc)[0]))
        Hbar = Hc if Hbar is None else Hbar + Hc
        print(f"[cfg {c+1}/{ncfg}] per-config min(P^T H P) = {permin[-1]:+.4f}",
              flush=True)
    Hbar /= ncfg
    ev = np.linalg.eigvalsh(Hbar)
    print(f"\nHbar = E[P^T H P], dim={ev.size}")
    print(f"  min = {ev[0]:+.4f}   frac_neg = {np.sum(ev < 0)}/{ev.size}")
    print(f"  lowest 10 = {np.round(ev[:10], 4)}")
    print(f"  per-config mins = {np.round(permin, 4)}  "
          f"(negative in {np.sum(np.array(permin) < 0)}/{ncfg})")

    # sanity: cold vacuum U=1 must be positive definite
    Uc = np.broadcast_to(S.ID3, (4, L, L, L, L, 3, 3)).copy()
    H0 = build_hessian(Uc, beta, L)
    H0c = (P.T @ (H0 @ P)).toarray()
    ev0 = np.linalg.eigvalsh(0.5 * (H0c + H0c.T))
    print(f"  sanity cold U=1: min = {ev0[0]:+.4f}")

    verdict = "DEAD (marginal negative somewhere)" if ev[0] < 0 else "INCONCLUSIVE (min>=0)"
    print(f"\nVERDICT: {verdict}")


if __name__ == '__main__':
    main()
