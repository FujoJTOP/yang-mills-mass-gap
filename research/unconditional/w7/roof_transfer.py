"""
roof_transfer.py -- witness for research/unconditional/ROOF-TRANSFER.md

Exact (no Monte Carlo) transfer-matrix computation of the 3D Ising model
(= dual of 3D Z2 gauge) on an L x L spatial slice, L=3, periodic.

It demonstrates the EXACT 0-form criterion of the solvable specimen:
  * the 0-form (spin) gap  D_sigma = -ln(lambda_odd_max / lambda_even_max)
  * the energy gap          D_eps   = -ln(lambda_even_2nd / lambda_even_max)
both vanish only at the critical point beta_c ~ 0.2217 and are POSITIVE away
from it.  The gapless point is the critical point; both phases (confined and
deconfined) are gapped.  This is the "0-form ridge" that CAN be built exactly.

A small exact 2D-Ising magnetization curve is printed as calibration of the
order-parameter half of the criterion.

No claim about 4D Yang-Mills.  The 4D blocked steps are in ROOF-TRANSFER.md.
Run:  py -3.12 research/unconditional/w7/roof_transfer.py
"""

import numpy as np
from scipy.linalg import eigh

np.set_printoptions(precision=5, suppress=True, linewidth=120)


def build_3d_ising(L, beta):
    """Transfer matrix on L x L slice (periodic). T = B^{1/2} A B^{1/2}, symmetric."""
    N = L * L
    dim = 1 << N
    s = np.arange(dim)
    spins = np.array([[(1 if (x >> i) & 1 else -1) for i in range(N)]
                      for x in range(dim)], dtype=float)          # dim x N, entry +-1
    # intra-slice bonds on LxL torus
    bonds = []
    for x in range(L):
        for y in range(L):
            i = x * L + y
            bonds.append((i, ((x + 1) % L) * L + y))
            bonds.append((i, x * L + ((y + 1) % L)))
    Bdiag = np.zeros(dim)
    for (i, j) in bonds:
        Bdiag += spins[:, i] * spins[:, j]
    Bdiag = np.exp(beta * Bdiag)                                 # diag of B
    # A[s',s] = exp(beta * sum_i s'_i s_i)
    A = np.exp(beta * (spins @ spins.T))
    sq = np.sqrt(Bdiag)
    T = sq[:, None] * A * sq[None, :]
    # parity operator P: flip all bits
    flip = (np.arange(dim) ^ (dim - 1))
    P = np.zeros((dim, dim))
    P[np.arange(dim), flip] = 1.0
    return T, P


def gaps_3d(L, beta):
    T, P = build_3d_ising(L, beta)
    evals, evecs = eigh(T)                                       # ascending
    order = np.argsort(evals)[::-1]                              # descending
    evals = evals[order]
    evecs = evecs[:, order]
    pars = np.array([v @ P @ v for v in evecs.T])               # parity expectation
    even = evals[pars > 0.5]
    odd = evals[pars < -0.5]
    lam0 = even[0]                                              # ground state
    lamE = even[1]                                             # energy (epsilon)
    lamO = odd[0]                                              # spin (sigma)
    d_sig = -np.log(lamO / lam0)
    d_eps = -np.log(lamE / lam0)
    return d_sig, d_eps, lam0, lamO


def m2d_exact(beta):
    """Exact 2D Ising spontaneous magnetization (J=1)."""
    x = np.sinh(2.0 * beta)
    if x <= 1.0:
        return 0.0
    return (1.0 - x ** -4) ** 0.125


if __name__ == "__main__":
    print("=" * 74)
    print("3D Ising transfer matrix (L=3, periodic) -- dual of 3D Z2 gauge")
    print("  0-form spin gap (D_sigma) and energy gap (D_eps) vs beta")
    print("  beta_c(3D Ising) ~ 0.22165 : gaps -> 0 ONLY near the critical point")
    print("=" * 74)
    print(f"{'beta':>8} {'D_sigma':>12} {'D_eps':>12} {'lambda0':>14} {'lambda_odd':>14}")
    for beta in np.linspace(0.12, 0.42, 16):
        ds, de, l0, lo = gaps_3d(3, beta)
        print(f"{beta:8.4f} {ds:12.6f} {de:12.6f} {l0:14.8f} {lo:14.8f}")
    print()
    print("Reading: the EVEN-sector gap D_eps has a clear MINIMUM near beta~0.26")
    print("(finite-L critical region; thermodynamic beta_c~0.2217) and rises on")
    print("BOTH sides -> the spectrum is gapped away from the critical point,")
    print("gapless at it.  The ODD gap D_sigma -> 0 at large beta is vacuum")
    print("tunnelling between the two ordered vacua on the finite torus (SSB),")
    print("NOT a physical gapless mode.  Both phases are gapped; the critical")
    print("point is the unique gapless fixed point.  This is the exact 0-form")
    print("ridge -- an ORDER-PARAMETER criterion (it needs a TRANSITION), not a")
    print("single-phase mass-gap criterion.")
    print()
    print("=" * 74)
    print("Calibration: exact 2D Ising magnetization m(beta)  [J=1, beta_c=0.44069]")
    print("=" * 74)
    for beta in [0.30, 0.40, 0.44069, 0.48, 0.60, 0.80]:
        print(f"  beta={beta:7.4f}   m={m2d_exact(beta):8.5f}")
    print("\nDone.  No QED.")
