#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CCM (Constraint-Cycle Memory) kill tests.

Object under test:
  E_i = conditional expectations onto overlapping local algebras A_i
        (for a probability/Gibbs measure these are ORTHOGONAL PROJECTIONS in L^2).
  H[] = E4 E3 E2 E1 - E1 E2 E3 E4      (a 4-cycle of overlaps).

For self-adjoint projections P_i, H[] = A - A*,  A = P4 P3 P2 P1.
=> H[] is skew-adjoint, tr H[] = 0 identically, H[] != 0  iff  A is NOT self-adjoint
   iff the P_i fail to commute in the cyclic order.

Exact Gaussian/probability/2D computations, numpy only, one job, -I safe, no downloads.
"""

import numpy as np

np.set_printoptions(precision=6, suppress=True, linewidth=160)

# --------------------------------------------------------------------------
# Core: orthogonal projection onto span of coordinate rows, in covariance metric.
# For a zero-mean Gaussian with covariance C on R^n, the conditional expectation
# E[. | sigma(x_j : j in S)] restricted to the first (linear) Wiener chaos is the
# C-orthogonal projection onto span{e_j : j in S}.  If instead the local algebra
# is span{b_j : j in S} (b_j = image of site j in a reduced/physical frame), the
# projection is A_S (A_S^T C A_S)^{-1} A_S^T C with A_S = [b_j].
# --------------------------------------------------------------------------
def proj_S(A_S, C):
    """C-orthogonal projector onto span(columns of A_S)."""
    M = A_S.T @ C @ A_S
    return A_S @ np.linalg.pinv(M) @ A_S.T @ C


def cycle_gap(P):
    """Return A - A*, the skew-adjoint H[] for the cyclic order P1,P2,P3,P4."""
    P1, P2, P3, P4 = P
    A = P4 @ P3 @ P2 @ P1
    B = P1 @ P2 @ P3 @ P4
    return A - B


def cnorm(D, C):
    """Operator norm of D in the C (covariance) inner product, via an isometry."""
    # C = Lg Lg^T ; z = Lg^T x maps (C-inner product) -> Euclidean.
    Lg = np.linalg.cholesky((C + C.T) / 2)
    D_z = Lg.T @ D @ np.linalg.inv(Lg.T)
    return np.linalg.norm(D_z, 2), np.linalg.norm(D_z, 'fro')


def check_projector(P, C):
    """Verify P^2=P and C-self-adjoint (P^T C = C P)."""
    idem = np.linalg.norm(P @ P - P)
    sadj = np.linalg.norm(P.T @ C - C @ P)
    return idem, sadj


# ==========================================================================
# TEST 1a : FREE MASSLESS GAUSSIAN (SCALAR) FIELD on a cycle, zero mean.
#           massless == gap 0 (thermodynamic limit); check H[] explicitly.
# ==========================================================================
def massless_gff_cycle(n, regions):
    # periodic graph Laplacian on C_n
    L = np.zeros((n, n))
    for i in range(n):
        L[i, i] = 2.0
        L[i, (i - 1) % n] -= 1.0
        L[i, (i + 1) % n] -= 1.0
    C = np.linalg.pinv(L)                      # massless Green function (= L^+)
    # orthonormal basis of 1^perp  (= im C)
    P0 = np.eye(n) - np.ones((n, n)) / n
    _, _, Vt = np.linalg.svd(P0)
    B = Vt[:n - 1, :]                          # (n-1) x n, rows orthonormal, span 1^perp
    Cpsi = B @ C @ B.T                         # (n-1)x(n-1) SPD
    # site j in reduced frame: phi_j = B[:,j] . psi
    b = [B[:, j] for j in range(n)]
    P = [proj_S(np.column_stack([b[j] for j in S]), Cpsi) for S in regions]
    return P, Cpsi, C


# ==========================================================================
# TEST 1b / 3 : FREE MASSIVE Gaussian field (gapped) and massless scan (gapless).
# ==========================================================================
def gff_cycle_massive(n, m2, regions):
    L = np.zeros((n, n))
    for i in range(n):
        L[i, i] = 2.0
        L[i, (i - 1) % n] -= 1.0
        L[i, (i + 1) % n] -= 1.0
    C = np.linalg.inv(L + m2 * np.eye(n))
    X = [np.eye(n)[:, S] for S in regions]
    P = [proj_S(Xi, C) for Xi in X]
    return P, C


# ==========================================================================
# TEST 4 : GAPPED PRODUCT MODEL (independent sites).  C diagonal.
# ==========================================================================
def product_model(n, regions, var=1.0):
    C = var * np.eye(n)
    X = [np.eye(n)[:, S] for S in regions]
    P = [proj_S(Xi, C) for Xi in X]
    return P, C


# ==========================================================================
# 4D FREE MASSLESS U(1) GAUGE FIELD (Maxwell), L=2 torus.
#   F = d1 A ; F is white Gaussian on im(d1) (covariance = projector onto im d1,
#   i.e. Euclidean identity on that subspace after an orthonormal change of frame).
#   Local algebra of a region = plaquette components F_p, p in region.
# ==========================================================================
def free_maxwell_4d_regions(L, region_sites):
    dims = [L, L, L, L]
    nsite = L ** 4
    def sid(x):  # site -> linear index
        return ((x[0] * L + x[1]) * L + x[2]) * L + x[3]
    # links: (site, mu)
    links = [(s, mu) for s in range(nsite) for mu in range(4)]
    lidx = {lk: i for i, lk in enumerate(links)}
    # plaquettes: (site, mu, nu), mu<nu
    plaq = []
    for s in range(nsite):
        x = tuple((s // (L ** (3 - k))) % L for k in range(4))
        for mu in range(4):
            for nu in range(mu + 1, 4):
                plaq.append((x, mu, nu))
    pidx = {p: i for i, p in enumerate(plaq)}
    # boundary operator d1: each plaquette = sum of 4 links
    d1 = np.zeros((len(plaq), len(links)))
    for (x, mu, nu), pi in pidx.items():
        # corners: x, x+e_mu, x+e_nu, x+e_mu+e_nu ; oriented boundary
        def shift(x, mu):
            y = list(x); y[mu] = (y[mu] + 1) % L; return tuple(y)
        x0 = x
        xmu = shift(x, mu)
        xnu = shift(x, nu)
        xmunu = shift(shift(x, mu), nu)
        d1[pi, lidx[(sid(x0), mu)]] += 1.0
        d1[pi, lidx[(sid(xmu), nu)]] += 1.0
        d1[pi, lidx[(sid(x0), nu)]] -= 1.0
        d1[pi, lidx[(sid(xnu), mu)]] -= 1.0
    # orthonormal basis of im(d1)
    U, sv, _ = np.linalg.svd(d1, full_matrices=False)
    rank = int(np.sum(sv > 1e-9))
    U = U[:, :rank]                            # p x rank, iid coordinates psi
    # map site linear index -> its 6 plaquettes
    site_plaq = {s: [] for s in range(nsite)}
    for (x, mu, nu), pi in pidx.items():
        site_plaq[sid(x)].append(pi)
    regions_p = [sorted({pi for s in reg for pi in site_plaq[s]}) for reg in region_sites]
    # in the iid frame: coordinate F_p = (U[p,:]) . psi  => b_p = U[p,:]
    b = [U[p, :] for p in range(len(plaq))]
    C = np.eye(rank)
    P = [proj_S(np.column_stack([b[p] for p in S]), C) for S in regions_p]
    return P, C, rank, len(plaq)


# ==========================================================================
# 2D FREE U(1) GAUGE FIELD (Maxwell) on LxL torus: exact CCM of the field
# strength.  F = d1 A ; F is white Gaussian on im(d1) = {Sigma_p F_p = 0}
# (one global constraint on the torus), so couplings are a single global mode.
# ==========================================================================
def free_maxwell_2d(L, region_sites):
    nsite = L * L
    def sid(x): return x[0] * L + x[1]
    links = [(s, mu) for s in range(nsite) for mu in range(2)]
    lidx = {lk: i for i, lk in enumerate(links)}
    plaq = [(s, 0, 1) for s in range(nsite)]
    pidx = {p: i for i, p in enumerate(plaq)}
    d1 = np.zeros((len(plaq), len(links)))
    for (s, mu, nu), pi in pidx.items():
        x = [s // L, s % L]
        def shift(x, mu):
            y = list(x); y[mu] = (y[mu] + 1) % L; return tuple(y)
        x0 = tuple(x); xmu = shift(x, mu); xnu = shift(x, nu)
        d1[pi, lidx[(sid(x0), mu)]] += 1.0
        d1[pi, lidx[(sid(xmu), nu)]] += 1.0
        d1[pi, lidx[(sid(x0), nu)]] -= 1.0
        d1[pi, lidx[(sid(xnu), mu)]] -= 1.0
    U, sv, _ = np.linalg.svd(d1, full_matrices=False)
    rank = int(np.sum(sv > 1e-9))
    U = U[:, :rank]
    b = [U[p, :] for p in range(len(plaq))]
    C = np.eye(rank)
    regions_p = [sorted(set(reg)) for reg in region_sites]  # one plaquette per site
    P = [proj_S(np.column_stack([b[p] for p in S]), C) for S in regions_p]
    return P, C, rank, len(plaq)


# ==========================================================================
# RUN
# ==========================================================================
def banner(t):
    print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)


results = {}

n = 12
regions = [
    [0, 1, 2, 3],
    [3, 4, 5, 6],
    [6, 7, 8, 9],
    [9, 10, 11, 0],
]
# a control arrangement that COMMUTES (nested chain P1<=P2<=P3<=P4) -> H[] must be 0
regions_nested = [[0, 1, 2], [0, 1, 2, 3], [0, 1, 2, 3, 4], [0, 1, 2, 3, 4, 5]]
# disjoint-but-correlated coordinate blocks (do NOT commute for a correlated field)
regions_disjoint = [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]

banner("TEST 1a  FREE MASSLESS scalar GFF on C_%d (gap=0 in thermo limit)" % n)
P, Cpsi, Cfull = massless_gff_cycle(n, regions)
idem = max(check_projector(Pi, Cpsi)[0] for Pi in P)
sadj = max(check_projector(Pi, Cpsi)[1] for Pi in P)
D = cycle_gap(P)
opn, frn = cnorm(D, Cpsi)
print("projector checks: max|P^2-P|=%.2e  max|P^T C - C P|=%.2e" % (idem, sadj))
print("COMMUTATOR  P4P3P2P1 - P1P2P3P4 :  ||.||_C=%.6e   ||.||_F=%.6e" % (opn, frn))
print("tr(H[]) = %.3e (identically 0, skew-adjoint)" % np.trace(D))
print("H[] != 0 ? ->", "NONZERO" if opn > 1e-9 else "ZERO")
results['massless_scalar'] = opn

banner("control A: NESTED local algebras on the SAME massless field (must commute)")
Pc, _, _ = massless_gff_cycle(n, regions_nested)
Dc = cycle_gap(Pc)
opc, frc = cnorm(Dc, Cpsi)
print("COMMUTATOR (nested) ||.||_C=%.3e  ||.||_F=%.3e -> %s"
      % (opc, frc, "ZERO" if opc < 1e-9 else "NONZERO"))
results['massless_nested'] = opc

banner("control B: DISJOINT-but-CORRELATED blocks (do NOT commute: not about disjointness)")
Pd, _, _ = massless_gff_cycle(n, regions_disjoint)
Dd = cycle_gap(Pd)
opd, frd = cnorm(Dd, Cpsi)
print("COMMUTATOR (disjoint) ||.||_C=%.3e  ||.||_F=%.3e -> %s"
      % (opd, frd, "ZERO" if opd < 1e-9 else "NONZERO"))
results['massless_disjoint'] = opd

banner("TEST 1b/3  FREE MASSIVE scalar GFF (gapped) mass scan")
for m2 in [0.0, 1e-4, 1e-2, 0.1, 1.0, 4.0]:
    if m2 == 0.0:
        Pm = P; Cm = Cpsi
    else:
        Pm, Cm = gff_cycle_massive(n, m2, regions)
    Dm = cycle_gap(Pm)
    omn, _ = cnorm(Dm, Cm)
    print("m^2=%-7g  ||H[]||_C = %.6e   (%s)" %
          (m2, omn, "GAPPED" if m2 > 0 else "GAPLESS"))
    results['massive_m2_%g' % m2] = omn

banner("TEST 4  GAPPED PRODUCT model (independent sites, var=1)")
Pp, Cp = product_model(n, regions)
Dp = cycle_gap(Pp)
opp, fpp = cnorm(Dp, Cp)
print("COMMUTATOR ||.||_C=%.3e  ||.||_F=%.3e -> %s"
      % (opp, fpp, "ZERO" if opp < 1e-12 else "NONZERO"))
results['product_gapped'] = opp

banner("TEST 1(c)  4D FREE MASSLESS U(1) GAUGE FIELD (Maxwell), L=2 torus")
# 16 sites in a 4x4-cycle of overlapping site-arcs
Lg = 2
site_regions = [list(range(0, 5)), list(range(4, 9)), list(range(8, 13)), list(range(12, 16)) + [0]]
Pg, Cg, rank, npl = free_maxwell_4d_regions(Lg, site_regions)
idemg = max(check_projector(Pi, Cg)[0] for Pi in Pg)
Dg = cycle_gap(Pg)
opg, frg = cnorm(Dg, Cg)
print("hypercube: sites=16 links=64 plaquettes=%d  dim(im d1)=%d"
      % (npl, rank))
print("  (F=d1A lives in the %d-dim exact-2-form space of the %d plaquette coords)" % (rank, npl))
print("projector check max|P^2-P|=%.2e" % idemg)
print("COMMUTATOR ||.||_C=%.6e  ||.||_F=%.6e  -> %s"
      % (opg, frg, "NONZERO" if opg > 1e-9 else "ZERO"))
results['free_maxwell_4d'] = opg

banner("TEST 2  2D free U(1) gauge field: exact CCM of field strength F=d1A")
Lg2 = 8
regs2 = [list(range(0, 20)), list(range(16, 36)), list(range(32, 52)), list(range(48, 64)) + list(range(0, 4))]
P2, C2, rk2, np2 = free_maxwell_2d(Lg2, regs2)
D2 = cycle_gap(P2)
op2, fr2 = cnorm(D2, C2)
print("L=%d: plaquettes=%d  dim(im d1)=%d  (#global constraints = %d, the torus flux)"
      % (Lg2, np2, rk2, np2 - rk2))
print("COMMUTATOR ||.||_C=%.6e  ||.||_F=%.6e  -> %s"
      % (op2, fr2, "NONZERO" if op2 > 1e-9 else "ZERO"))
results['free_maxwell_2d'] = op2

banner("SUMMARY")
print("massless scalar (gapless) : H[] = %.4e   -> %s"
      % (results['massless_scalar'], "NONZERO" if results['massless_scalar'] > 1e-9 else "ZERO"))
print("gapped product model      : H[] = %.4e   -> %s"
      % (results['product_gapped'], "ZERO" if results['product_gapped'] < 1e-12 else "NONZERO"))
print("free massless 4D photon   : H[] = %.4e   -> %s"
      % (results['free_maxwell_4d'], "NONZERO" if results['free_maxwell_4d'] > 1e-9 else "ZERO"))
print("free massless 2D photon   : H[] = %.4e   -> %s"
      % (results['free_maxwell_2d'], "NONZERO" if results['free_maxwell_2d'] > 1e-9 else "ZERO"))
print()
print("=> CCM nonzero WITHOUT a gap (free massless): NOT NECESSARY (fires on free field).")
print("=> CCM zero WITH a gap (gapped product / 2D YM): NOT SUFFICIENT.")
print("No QED. Clay 0%.")
