"""Gauge-mode removal from the coarse-projected Hessian P^T H P.

For each config:
  * build the coarse links U' (b=2 block average + polar projection);
  * build the coarse gauge image G_c (512 x 128): infinitesimal coarse gauge
    transformations of the coarse links, in LEFT algebra coordinates
        dxi^a_(x,mu) = omega^a(x) - (Ad_{U'_mu(x)})^a_b omega^b(x+mu);
  * sanity: Hc @ G_c ~ 0 (gauge invariance of the coarse Hessian);
  * project Hc onto the complement of span(G_c) and compute its spectrum there;
  * average the sorted physical spectra over configs.

If the averaged physical min stays positive -> the E-part is genuinely positive.
If it goes negative -> the positivity was a gauge artifact.
"""
import itertools
import sys

import numpy as np
from scipy.linalg import null_space

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/cpb-attack/11-slc-blocked')
import su3 as S  # noqa: E402
from full_hessian import build_hessian  # noqa: E402
from coarse_block_hessian import lift_matrix  # noqa: E402

T = S.T_GEN


def coarse_links(U, L, b=2):
    D = L // b
    Uc = np.empty((4, D, D, D, D, 3, 3), dtype=complex)
    for mu in range(4):
        Y = U[mu].reshape(D, b, D, b, D, b, D, b, 3, 3)
        Uc[mu] = S.project_su3(Y.mean(axis=(1, 3, 5, 7)))
    return Uc


def gauge_image(Uc, D):
    Nc = 4 * D ** 4 * 8
    cidx = lambda mu, xp: (((mu * D + xp[0]) * D + xp[1]) * D + xp[2]) * D + xp[3]
    cols = []
    sites = list(itertools.product(range(D), repeat=4))
    for site in sites:
        for a0 in range(8):
            v = np.zeros(Nc)
            for mu in range(4):
                Unu = Uc[mu]
                for y in sites:
                    Ulink = Unu[y]
                    xp1 = list(y); xp1[mu] = (xp1[mu] + 1) % D; xp1 = tuple(xp1)
                    idx = cidx(mu, y) * 8
                    # LEFT variation: dxi^a = (Ad_{U^{-1}} omega(x))^a - omega^a(x+mu)
                    if y == site:
                        for a in range(8):
                            M = np.trace(T[a] @ S.adj(Ulink) @ T[a0] @ Ulink).real / 2.0
                            v[idx + a] += M
                    if xp1 == site:
                        v[idx + a0] -= 1.0
            cols.append(v)
    return np.array(cols).T  # (Nc, 128)


def main():
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    ncfg = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    beta = 2.5
    D = L // 2
    rng = np.random.default_rng(21)
    lat = S.Lattice(L, beta, rng, cold=True)
    for _ in range(300):
        lat.sweep()
    P = lift_matrix(L, 2)
    Nc = 4 * D ** 4 * 8

    spectra, gauge_res = [], []
    for c in range(ncfg):
        lat.sweep()
        U = lat.U
        H = build_hessian(U, beta, L)
        Hc = (P.T @ (H @ P)).toarray()
        Hc = 0.5 * (Hc + Hc.T)
        Uc = coarse_links(U, L, 2)
        Gc = gauge_image(Uc, D)
        # orthonormal basis of the gauge span
        Umat, sv, _ = np.linalg.svd(Gc, full_matrices=False)
        r = int(np.sum(sv > 1e-8 * sv[0]))
        V = Umat[:, :r]
        gauge_res.append(float(np.linalg.norm(Hc @ V) / max(np.linalg.norm(Hc), 1e-12)))
        Q = null_space(V.T)                      # (Nc, Nc-r) complement basis
        Hcomp = Q.T @ Hc @ Q
        Hcomp = 0.5 * (Hcomp + Hcomp.T)
        ev = np.sort(np.linalg.eigvalsh(Hcomp))
        spectra.append(ev)
        print(f"[cfg {c+1}/{ncfg}] rank(gauge)={r}  complement dim={ev.size}  "
              f"physical min={ev[0]:+.4f}  gauge_res={gauge_res[-1]:.2e}", flush=True)

    m = min(e.size for e in spectra)
    A = np.mean([e[:m] for e in spectra], axis=0)   # averaged sorted physical spectrum
    print(f"\nAveraged physical spectrum (complement, {m} modes):")
    print(f"  min = {A[0]:+.4f}   lowest 10 = {np.round(A[:10], 4)}")
    print(f"  per-config physical mins = {np.round([e[0] for e in spectra], 4)}")
    print(f"  gauge invariance residual (||Hc V||/||Hc||) = "
          f"{np.round(gauge_res, 2e0 if False else 5)}")
    print(f"\nVERDICT: {'physical min POSITIVE (E-part genuinely positive)' if A[0] > 0 else 'physical min NEGATIVE (positivity was gauge artifact / E-part negative)'}")


if __name__ == '__main__':
    main()
