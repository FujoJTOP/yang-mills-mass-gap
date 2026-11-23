"""W6 numeric check: does collective plane-wave positivity miss the negative modes?

On ONE thermalized SU(2) config at beta=2.5 (L=3):
  (1) full fine gauge Hessian H = d^2 S  -> min eigenvalue (expect < 0, indefinite)
  (2) plane-wave collective subspace (the 11-slc-blocked probe's object) -> min (expect > 0)
  (3) overlap of H's most-negative eigenvector with the plane-wave subspace (expect small)
  (4) Ric + Hess S along the most-negative direction, Ric = +3/2 I (generous upper bound
      for SU(2) bi-invariant metric in this normalisation) -> expect still negative.

Reproduces the probe's positivity AND shows it is a subspace statement, not a
pointwise Bakry-Emery bound.
"""
import sys
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
import w4_blocked as W


def build(U, beta, L):
    ns = L ** 4
    dim = 4 * ns * 3
    H = W.hessian2(U, beta)                     # dense (dim,dim), SU(2)
    coords = np.stack(np.meshgrid(*[np.arange(L)] * 4, indexing='ij'),
                      axis=-1).reshape(ns, 4)
    return H, coords, ns, dim


def planewave_basis(coords, ns, dim, L, nmax=1):
    """Collective plane-wave subspace: per direction nu and colour a, perturb ALL
    links of direction nu by cos(p.x) t^a and sin(p.x) t^a, p along one axis.
    Spans the 'collective mode' object of the 11-slc-blocked probe."""
    cols = []
    for nu in range(4):
        for axis in range(4):
            for n in range(1, nmax + 1):
                k = 2 * np.pi * n / L
                phase = np.exp(1j * (coords @ (k * np.eye(4)[axis])))
                for a in range(3):
                    for part in (phase.real, phase.imag):
                        v = np.zeros((4, ns, 3))
                        v[nu, :, a] = part / np.sqrt(ns)
                        cols.append(v.reshape(-1))
    M = np.array(cols).T
    Q, _ = np.linalg.qr(M)                       # orthonormalise (drop nulls)
    return Q


def main():
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    beta = 2.5
    rng = np.random.default_rng(11)
    U = W.sample2(L, beta, rng, ntherm=300, nkeep=1, gapevery=1)[0]
    H, coords, ns, dim = build(U, beta, L)

    ev, V = np.linalg.eigh(H)
    print(f"L={L} beta={beta}  dim={dim}")
    print(f"[full Hessian]   min={ev[0]:+.4f}  min/beta={ev[0]/beta:+.4f}  "
          f"frac_neg={np.mean(ev < 0):.4f}  n_neg={np.sum(ev < 0)}")

    pl = W.perlink_min2(U, beta)
    print(f"[per-link]       frac_neg={np.mean(pl < 0):.3f}  min/beta={pl.min()/beta:+.4f}")

    Qp = planewave_basis(coords, ns, dim, L, nmax=1)
    B = Qp.T @ H @ Qp
    B = 0.5 * (B + B.T)
    eb = np.linalg.eigvalsh(B)
    print(f"[plane-wave sub] dim={Qp.shape[1]:3d}  min={eb.min():+.4f}  "
          f"min/beta={eb.min()/beta:+.4f}  frac_neg={np.mean(eb < 0):.3f}")

    # overlap of the most-negative eigenvector with the plane-wave subspace
    w0 = V[:, 0]
    ow = np.linalg.norm(Qp.T @ w0) ** 2
    print(f"[overlap] most-negative eigvec weight in plane-wave subspace = {ow:.4f}")

    # BE curvature along the most-negative direction: Ric + Hess S, Ric=+1.5 I (generous)
    ric = 1.5
    print(f"[BE curvature] min eig of (Ric + Hess S), Ric=+{ric}: "
          f"{ev[0] + ric:+.4f}  (negfrac={np.mean(ev + ric < 0):.4f})")


if __name__ == '__main__':
    main()
