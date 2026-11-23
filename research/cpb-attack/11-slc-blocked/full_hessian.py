"""Full fine Wilson Hessian spectrum at the crossover (L=2 test, L=4 real).

The one-loop (Gaussian-fluctuation) marginal blocked Hessian is the Schur
complement of the fine Hessian.  A Schur complement of a positive-definite
matrix is positive-definite, so the decisive question is:

    Is the FULL fine Hessian positive definite at beta=2.5?

If yes -> the one-loop blocked Hessian is PD (SLC alive at the Gaussian level).
If no  -> the Gaussian marginal is ill-defined (the fluctuation integral
          diverges); a non-Gaussian treatment is required.

We build the exact fine Hessian by per-plaquette finite differences of
-(beta/3) Re tr U_p in the left variations U_e -> exp(i xi^a T^a) U_e.
"""
import itertools
import sys

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
import su3 as S  # noqa: E402

ADJ = S.adj
EPS = 1e-3


def flat_idx(mu, x, L):
    return ((mu * L + x[0]) * L + x[1]) * L * L + x[2] * L + x[3]


def build_hessian(U, beta, L):
    N = 4 * L ** 4 * 8
    rows, cols, vals = [], [], []
    xs = list(itertools.product(range(L), repeat=4))
    d = 32
    # coordinate -> (link_slot 0..3, algebra a)
    coords = [(e, a) for e in range(4) for a in range(8)]
    # precompute the finite-difference point list
    points = []          # each: dict of (coord_index) -> sign list, plus tag
    single = {}
    for i, (e, a) in enumerate(coords):
        single[(i, +1)] = len(points); points.append([(i, +1)])
        single[(i, -1)] = len(points); points.append([(i, -1)])
    mixed = {}
    for i in range(d):
        for j in range(i + 1, d):
            for si in (+1, -1):
                for sj in (+1, -1):
                    mixed[(i, si, j, sj)] = len(points)
                    points.append([(i, si), (j, sj)])
    P = len(points)

    for mu in range(4):
        for nu in range(mu + 1, 4):
            for x in xs:
                xm = list(x); xm[mu] = (xm[mu] + 1) % L
                xn = list(x); xn[nu] = (xn[nu] + 1) % L
                U1 = U[mu][x[0], x[1], x[2], x[3]]
                U2 = U[nu][xm[0], xm[1], xm[2], xm[3]]
                U3 = U[mu][xn[0], xn[1], xn[2], xn[3]]
                U4 = U[nu][x[0], x[1], x[2], x[3]]
                base = [U1, U2, U3, U4]
                linkflat = [flat_idx(mu, x, L), flat_idx(nu, xm, L),
                            flat_idx(mu, xn, L), flat_idx(nu, x, L)]

                X = np.zeros((P, 4, 8))
                for k, pts in enumerate(points):
                    for (ci, s) in pts:
                        e, a = coords[ci]
                        X[k, e, a] += s * EPS
                Up_list = [S.unitaries_from_xi(X[:, e, :], base[e]) for e in range(4)]
                Up = Up_list[0] @ Up_list[1] @ ADJ(Up_list[2]) @ ADJ(Up_list[3])
                f = np.real(np.trace(Up, axis1=-2, axis2=-1))  # (P,)
                coeff = -beta / 3.0

                Hp = np.zeros((d, d))
                f0 = np.real(np.trace(U1 @ U2 @ ADJ(U3) @ ADJ(U4)))
                for i in range(d):
                    fpi = f[single[(i, +1)]]; fmi = f[single[(i, -1)]]
                    Hp[i, i] = coeff * (fpi - 2 * f0 + fmi) / EPS ** 2
                    for j in range(i + 1, d):
                        a = f[mixed[(i, +1, j, +1)]]; b = f[mixed[(i, +1, j, -1)]]
                        c = f[mixed[(i, -1, j, +1)]]; dd = f[mixed[(i, -1, j, -1)]]
                        Hp[i, j] = Hp[j, i] = coeff * (a - b - c + dd) / (4 * EPS ** 2)

                for i in range(d):
                    ei, ai = coords[i]
                    gi = linkflat[ei] * 8 + ai
                    for j in range(d):
                        ej, aj = coords[j]
                        gj = linkflat[ej] * 8 + aj
                        if abs(Hp[i, j]) > 1e-12:
                            rows.append(gi); cols.append(gj); vals.append(Hp[i, j])

    H = sparse.coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsr()
    H = 0.5 * (H + H.T)
    return H


def main():
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    beta = 2.5
    sweeps = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    rng = np.random.default_rng(7)
    lat = S.Lattice(L, beta, rng, cold=True)
    for _ in range(300):
        lat.sweep()

    mins = []
    for s in range(sweeps):
        lat.sweep()
        H = build_hessian(lat.U, beta, L)
        N = H.shape[0]
        k = min(20, N - 2)
        ev = eigsh(H, k=k, which='SA', return_eigenvectors=False)
        ev = np.sort(ev)
        mins.append(float(ev[0]))
        if s < 5 or (s + 1) % 10 == 0:
            print(f"[{s+1}/{sweeps}] L={L} N={N} min={ev[0]:+.5f} "
                  f"2nd={ev[1]:+.5f} neg={np.sum(ev < 0)}/{k}", flush=True)
    mins = np.array(mins)
    print(f"\nL={L} beta={beta}: min-eig over {sweeps} configs: "
          f"mean={mins.mean():+.5f} min={mins.min():+.5f} max={mins.max():+.5f} "
          f"frac_configs_negative={np.mean(mins < 0):.3f}")


if __name__ == '__main__':
    main()
