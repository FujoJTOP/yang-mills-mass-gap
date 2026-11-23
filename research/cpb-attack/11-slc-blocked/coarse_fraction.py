"""Where do the negative Hessian modes live: coarse (block-uniform) or fluctuation?

Build the full fine Hessian (full_hessian.build_hessian), take its lowest
eigenvectors, and project them onto the b=2 block-uniform subspace
(lift = constant on each 2^4 block).

  coarse_frac = 16 * sum_blocks mean(v)^2 / sum v^2 .

A negative eigenvector with coarse_frac ~ 0 is a fluctuation (high-momentum)
mode -> the *blocked* action can still be convex.  coarse_frac ~ 1 means the
blocked action itself is non-convex.
"""
import sys

import numpy as np
from scipy.sparse.linalg import eigsh

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/cpb-attack/11-slc-blocked')
import su3 as S  # noqa: E402
from full_hessian import build_hessian  # noqa: E402


def coarse_frac(v, L, b=2):
    D = L // b
    V = v.reshape(4, L, L, L, L, 8)
    Y = V.reshape(4, D, b, D, b, D, b, D, b, 8)
    mean = Y.mean(axis=(2, 4, 6, 8))          # (4,D,D,D,D,8)
    num = (b ** 4) * np.sum(mean ** 2)
    return float(num / np.sum(v ** 2))


def main():
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    beta = 2.5
    ncfg = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    rng = np.random.default_rng(11)
    lat = S.Lattice(L, beta, rng, cold=True)
    for _ in range(300):
        lat.sweep()

    for s in range(ncfg):
        lat.sweep()
        H = build_hessian(lat.U, beta, L)
        ev, vec = eigsh(H, k=10, which='SA')
        order = np.argsort(ev)
        ev, vec = ev[order], vec[:, order]
        fracs = [coarse_frac(vec[:, i], L) for i in range(ev.shape[0])]
        print(f"[cfg {s+1}] lowest 10 eig:")
        for i in range(ev.shape[0]):
            print(f"   eig[{i}] = {ev[i]:+.4f}   coarse_frac = {fracs[i]:.4f}")
        print(flush=True)


if __name__ == '__main__':
    main()
