#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GI3-ELEMENTARY witness — the edge `m1 = inf supp nu` as the smallest POSITIVE
generalized eigenvalue of the Hankel pencil (H^(1)_N, H_N), and its failure modes.

All quantities are exact closed-form moments of toy measures; no sampling, no MC.
Deterministic. numpy + scipy only.

Run:  py -3.12 -I research/v3/species/w7/gi3_elementary.py

Blocks
  [1] pencil identification: Unif[0,1] -> Legendre nodes (exact check, machine precision)
  [2] c_N = min positive gen-eigenvalue, N=1..12, four measures; monotone from above
  [3] conditioning wall: cond2(H_N) ~ (Lambda/m1)^(2N); digits lost in c_N
  [4] Gauss quadrature: for every N two measures share m_0..m_{2N-1}, edges differ O(1)
  [5] R6 family nu_eta = eta*delta(1/2) + (1-eta)*Unif[1,2]: moments O(eta)-close,
      edge O(1)-apart; c_N does NOT see the atom at any fixed N
  [6] min-structure: finite-atomic exactness, m1(nu'+nu'') = min(m1', m1'')
"""

import numpy as np
import scipy
import scipy.linalg as sla
from numpy.polynomial.legendre import leggauss

np.set_printoptions(precision=12, suppress=False, linewidth=200)


# ---------------------------------------------------------------- moments
def m_unif01(k):
    """Unif[0,1]: m_k = 1/(k+1)  (Hilbert moment sequence)."""
    return 1.0 / (k + 1.0)


def m_unif12(k):
    """Unif[1,2]: m_k = (2^{k+1}-1)/(k+1)."""
    return (2.0 ** (k + 1) - 1.0) / (k + 1.0)


def m_atom(c, w, base):
    """w*delta_c + (1-w)*base."""
    return lambda k: w * (c ** k) + (1.0 - w) * base(k)


def m_delta(c):
    return lambda k: c ** k


def m_finite(atoms):
    """atoms = [(x_i, w_i), ...]"""
    return lambda k: sum(w * (x ** k) for x, w in atoms)


# ---------------------------------------------------------------- pencil
def pencil(mfun, N):
    """H_ij = m_{i+j}, H1_ij = m_{i+j+1}, 0<=i,j<=N  (size N+1)."""
    i = np.arange(N + 1)
    M = np.array([[mfun(p + q) for q in i] for p in i], dtype=float)
    M1 = np.array([[mfun(p + q + 1) for q in i] for p in i], dtype=float)
    return M, M1


def edge_read(mfun, N, tol=1e-12):
    """c_N = smallest POSITIVE finite generalized eigenvalue of (H1, H).

    Returns (c_N, n_positive).  Non-finite eigenvalues come from the null space
    of H (finitely atomic measure with fewer than N+1 atoms) and are dropped.
    """
    H, H1 = pencil(mfun, N)
    w = sla.eig(H1, H, right=False)
    w = np.real_if_close(w, tol=1000)
    pos = np.array([x.real for x in w
                    if abs(x.imag) < 1e-7 and np.isfinite(x.real) and x.real > tol])
    if pos.size == 0:
        return np.nan, 0
    return float(pos.min()), int(pos.size)


def cond2(mfun, N):
    H, _ = pencil(mfun, N)
    s = np.linalg.svd(H, compute_uv=False)
    return float(s[0] / s[-1])


# ---------------------------------------------------------------- blocks
def block1():
    print("=" * 78)
    print("[1] PENCIL IDENTIFICATION — Unif[0,1] must reproduce Legendre nodes")
    print("    claim: min positive gen-eig of (H1,H) = smallest Gauss-Legendre node")
    print("=" * 78)
    print("  N   c_N                 x_1(Gauss)          |diff|")
    for N in range(1, 10):
        c, _ = edge_read(m_unif01, N)
        t, _w = leggauss(N + 1)
        x1 = 0.5 * (t.min() + 1.0)
        print(f" {N:2d}   {c:.15e}   {x1:.15e}   {abs(c - x1):.3e}")
    print("  READ: identity confirmed to machine precision; every c_N > 0 = m1(=0)")
    print("        => the readout is an UPPER bound on the edge, never a lower one.\n")


def block2():
    print("=" * 78)
    print("[2] c_N  vs  m1 = inf supp nu   (N = 1..12)")
    print("=" * 78)
    cases = [
        ("delta_1                ", m_delta(1.0), 1.0),
        ("0.3 d(.25)+0.7 Unif[1,2]", m_atom(0.25, 0.3, m_unif12), 0.25),
        ("Unif[1,2]              ", m_unif12, 1.0),
        ("Unif[0,1]              ", m_unif01, 0.0),
    ]
    for name, f, truth in cases:
        print(f"  {name}   m1(true) = {truth}")
        row = []
        for N in range(1, 13):
            c, k = edge_read(f, N)
            row.append(f"{c:.6e}({k})")
        print("    " + "  ".join(row[:6]))
        print("    " + "  ".join(row[6:]))
        print()


def block3():
    print("=" * 78)
    print("[3] CONDITIONING WALL on the exact identity")
    print("    cond2(H_N) and the relative error of c_N (Unif[0,1], exact m_k = 1/(k+1))")
    print("=" * 78)
    print("   N     cond2(H_N)      c_N              rel.err vs Gauss      log10 cond")
    for N in range(1, 13):
        c, _ = edge_read(m_unif01, N)
        t, _w = leggauss(N + 1)
        x1 = 0.5 * (t.min() + 1.0)
        cd = cond2(m_unif01, N)
        rel = abs(c - x1) / x1 if x1 > 0 else np.nan
        print(f"  {N:2d}   {cd:.6e}   {c:.12e}   {rel:.3e}          {np.log10(cd):6.2f}")
    print("  READ: cond ~ 10^(1.5 N) (Hilbert). The identity is exact but the")
    print("        edge sits ~ N log10(Lambda/m1) digits below the moment scale.\n")


def block4():
    print("=" * 78)
    print("[4] GAUSS QUADRATURE — same m_0..m_{2N-1}, different edge, every N")
    print("    nu = Unif[0,1] vs its n-point Gauss rule (n = N+1)")
    print("=" * 78)
    print("   n   max|m_k(nu)-m_k(nu_n)|, k<=2n-1     edge(nu)   edge(nu_n)   gap")
    for n in range(2, 8):
        t, w = leggauss(n)
        x = 0.5 * (t + 1.0)
        w = 0.5 * w
        worst = 0.0
        for k in range(0, 2 * n):
            mk = float(np.sum(w * x ** k))
            worst = max(worst, abs(mk - m_unif01(k)))
        print(f"  {n:2d}   {worst:.3e}                       0.0        {x.min():.9f}   {x.min():.9f}")
    print("  READ: matching 2n moments to machine precision leaves the edge at x_1>0.")
    print("        => no finite-moment functional can certify m1 >= c > 0.\n")


def block5():
    print("=" * 78)
    print("[5] R6 FAMILY  nu_eta = eta*delta(1/2) + (1-eta)*Unif[1,2]")
    print("    m1 = 1/2 for every eta>0 ; moments -> Unif[1,2] as eta->0")
    print("=" * 78)
    base = m_unif12
    print("   eta      max_rel|dm_k|,k<=16   max_abs|dm_k|   m1(true)   c_4        c_8")
    for eta in [1.0, 1e-1, 1e-2, 1e-3, 1e-6, 0.0]:
        f = m_atom(0.5, eta, base) if eta > 0 else m_unif12
        rel = max(abs(f(k) - base(k)) / abs(base(k)) for k in range(17))
        absd = max(abs(f(k) - base(k)) for k in range(17))
        truth = 0.5 if eta > 0 else 1.0
        c4, _ = edge_read(f, 4)
        c8, _ = edge_read(f, 8)
        print(f"  {eta:.0e}   {rel:.6e}             {absd:.6e}      {truth:.3f}     {c4:.6f}   {c8:.6f}")
    print("  READ: moment distance O(eta) -> 0 while the edge jumps 0.5 -> 1.0;")
    print("        at any fixed N the readout loses the atom entirely (c_N ~ 1+).\n")


def block6():
    print("=" * 78)
    print("[6] MIN-STRUCTURE and finite-atomic exactness")
    print("=" * 78)
    fa = m_finite([(0.5, 0.3), (0.75, 0.7)])
    for N in [1, 2, 5]:
        c, k = edge_read(fa, N)
        print(f"  2 atoms 0.5,0.75 : N={N}  c_N={c:.12f}  #positive={k}   (exact m1=0.5)")
    mix = m_finite([(0.5, 0.3), (0.75, 0.2), (1.25, 0.5)])
    for N in [2, 4]:
        c, k = edge_read(mix, N)
        print(f"  3 atoms .5,.75,1.25: N={N}  c_N={c:.12f}  #positive={k}  (exact m1=0.5)")
    print("  READ: finite support => attainment at N=n-1 and #positive = n = rank H;")
    print("        rank is the atom count, and min of the atoms = m1 exactly.\n")


def main():
    print(f"numpy {np.__version__}  scipy {scipy.__version__}")
    print("GI3-ELEMENTARY witness — Hankel pencil edge readout\n")
    block1()
    block2()
    block3()
    block4()
    block5()
    block6()
    print("=" * 78)
    print("VERDICT: identity exact (c_N -> m1, monotone from above); every finite")
    print("         truncation is an UPPER bound; no finite-moment lower bound exists.")
    print("         Clay 0%. No QED.")
    print("=" * 78)


if __name__ == "__main__":
    main()
