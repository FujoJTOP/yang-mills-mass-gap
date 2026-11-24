#!/usr/bin/env python3
"""
s4_decimation.py -- witness for S4-DECIMATION.md

Illustrates the *logical relation* between (star) and (star dec) in the decimation step
of the UV RG map.  A witness, not a proof; everything here is elementary linear algebra /
Gaussian model, exactly checkable.

Sections
  S4.1  Gershgorin interface: (star) alone (diagonal = delta_W) does NOT make G'' positive
        definite; with an exponential score covariance (star dec), positivity needs
        kappa > kappa*(delta_W,C).  Row-sum is volume-free for kappa>0 (tempered).
  S4.2  (star) NOT=> (star dec): bounded diagonal variance (Var=A) + power-law covariance
        => no exponential decay for ANY kappa>0; row-sum diverges with volume.
  S4.3  (star dec) NOT=> (star): exponential decay (kappa>0) but Savvidy-negative mean
        Hessian => diagonal D_B < delta_W, so (star) fails while (star dec) holds.
  S4.4  Free-field acid test: free massless field has G''(0)=0 (star FAILS) while the
        score covariance is a finite-range (local) kernel (star dec HOLDS).

Deterministic; py -3.12, numpy.
"""

import numpy as np

np.set_printoptions(precision=4, suppress=True)


def torus_rowsum(kernel, M):
    """Max row-sum over an M x M torus of a radial kernel kernel(r) (r = minimum-image
    Euclidean distance), excluding the origin.  Translation invariant => one row enough.
    Enumerates separations a,b in [0, M//2] with the correct multiplicities."""
    half = M // 2
    a = np.arange(half + 1)
    # multiplicity of min-image separation a on a ring of length M
    w = np.where((a == 0) | (2 * a == M), 1.0, 2.0)
    A, B = np.meshgrid(a, a, indexing="ij")
    Wt = np.outer(w, w)
    R = np.sqrt(A ** 2 + B ** 2)
    val = kernel(R) * Wt
    # remove the origin site (a=b=0)
    val[0, 0] -= kernel(np.array([0.0]))[0]
    return val.sum()


def torus_dist_matrix(M):
    """Full minimum-image Euclidean distance matrix on an M x M torus (small M only)."""
    idx = np.arange(M)
    d1 = np.minimum(idx, M - idx).astype(float)
    # 2D coords
    X, Y = np.meshgrid(idx, idx, indexing="ij")
    coords = np.stack([X.ravel(), Y.ravel()], 1).astype(float)
    diff = coords[:, None, :] - coords[None, :, :]
    diff = np.abs(diff)
    diff = np.minimum(diff, M - diff)          # minimum image per component
    return np.sqrt((diff ** 2).sum(-1))


def section_41():
    print("[S4.1] Gershgorin interface  (d=2 torus, C=1, delta_W=0.5)")
    delta = 0.5
    C = 1.0
    # row-sum on increasing tori: exponential kernel  C exp(-kappa r)
    print("  volume-freeness of the row-sum R(M) = sum_{r != 0} C exp(-kappa |r|):")
    print("    kappa    M=16        M=32        M=64        M=96      (R -> R_inf?)")
    for kappa in [0.30, 0.50, 0.80, 1.20]:
        row = [torus_rowsum(lambda r: C * np.exp(-kappa * r), M) for M in [16, 32, 64, 96]]
        print("    %.2f   %9.4f  %9.4f  %9.4f  %9.4f" % (kappa, *row))

    # actual G'' = delta I - (Cov - diag) on an M x M torus.  The exponential kernel is
    # circulant, so eigenvalues of the off-diagonal part are the 2D DFT of its torus array.
    def off_rowsum_and_lammax(kappa, M):
        a = np.arange(M)
        da = np.minimum(a, M - a).astype(float)
        h = C * np.exp(-kappa * np.sqrt(da[:, None] ** 2 + da[None, :] ** 2))
        h[0, 0] = 0.0                         # off-diagonal part (diagonal handled by delta)
        R = h.sum()                           # total = each site's row-sum (translation invariant)
        lmax = np.fft.fft2(h).real.max()      # max eigenvalue of the off-diagonal block (=R at k=0)
        return R, lmax

    M = 24
    print("  spectrum on M=24:  G'' = delta_W I - (Cov - diag);  lambda_min = delta_W - lambda_max(off)")
    print("    kappa    rowsum R     lambda_min     Gershgorin lb   (star dec)->pos?")
    for kappa in [0.20, 0.40, 0.60, 0.70, 0.80, 1.00, 1.20, 1.60]:
        R, lmax = off_rowsum_and_lammax(kappa, M)
        lm = delta - lmax
        print("    %.2f   %9.4f   %9.4f     %9.4f       %s"
              % (kappa, R, lm, delta - R, "YES" if lm > 0 else "NO"))
    # empirical kappa* (sign flip of lambda_min) by bisection; use a big torus so R ~ R_inf
    Mbig = 240
    lo, hi = 0.2, 6.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        R, lmax = off_rowsum_and_lammax(mid, Mbig)
        if delta - lmax > 0:
            hi = mid
        else:
            lo = mid
    kstar_emp = 0.5 * (lo + hi)
    # asymptotic estimate  R_inf ~ c_d * C / kappa^d,  with the polar-integral constant
    # c_d = |S^{d-1}| * Gamma(d) = 2 pi^{d/2} Gamma(d) / Gamma(d/2).  For d=2 this is 2*pi
    # (int_{R^2} exp(-kappa r) d^2r = 2 pi / kappa^2).
    from math import pi, gamma
    d = 2
    c_d = 2.0 * pi ** (d / 2.0) * gamma(d) / gamma(d / 2.0)
    c_2 = c_d
    kstar_asym = np.sqrt(c_2 * C / delta)
    print("  empirical kappa* (lambda_min crosses 0, M=%d) = %.3f" % (Mbig, kstar_emp))
    print("  small-kappa continuum estimate sqrt(c_d C/delta_W) = %.3f  (c_2=2*pi; an UPPER heuristic)" % kstar_asym)
    print("  => (star) alone (diagonal = delta_W) does NOT give positivity; sign flip at kappa*.")
    print()


def section_42():
    print("[S4.2] (star) NOT=> (star dec): bounded diagonal, power-law covariance")
    A = 1.0
    alpha = 0.5  # d=2 massless-ish Green function exponent (alpha = d-2+eps -> eps=0.5)
    print("  same diagonal A = Var(s_B) = 1 (so (star) can hold), kernel A/(1+r)^alpha, alpha=%.2f" % alpha)
    print("    r       Cov       r       Cov      r       Cov")
    for r in [2, 4, 8, 16, 32]:
        print("    %3d   %8.4f" % (r, A / (1.0 + r) ** alpha))
    # log-log slope
    rr = np.array([2, 4, 8, 16, 32, 64], float)
    cc = A / (1.0 + rr) ** alpha
    slope = np.polyfit(np.log(rr), np.log(cc), 1)[0]
    print("  log-log slope = %.3f  (power law: NO exponential with any kappa>0)" % slope)
    # row-sum on increasing tori -> diverges
    print("  row-sum R(M) (grows with volume => (star dec) fails, no volume-free bound):")
    print("    M=16      M=32      M=64      M=96      M=128")
    vals = [torus_rowsum(lambda r: A / (1.0 + r) ** alpha, M) for M in [16, 32, 64, 96, 128]]
    print("   " + "  ".join("%8.3f" % v for v in vals))
    print("  => diagonal bounded (Var=A=1) yet (star dec) fails for every kappa>0.")
    print()


def section_43():
    print("[S4.3] (star dec) NOT=> (star): decay holds, but curvature comparison fails")
    kappa = 0.8
    C = 1.0
    beta = 2.5
    EH_plus_W = -1.1 * beta   # Savvidy floor (PROVE-Gpp-positive sec.3)
    var = C
    delta = 0.5
    Ddiag = EH_plus_W - var
    print("  exponential kernel C=1, kappa=%.1f (=> (star dec) HOLDS)" % kappa)
    print("  mean Hessian E[H]+W = -1.1*beta = %.2f ,  Var = %.2f ,  delta_W = %.2f"
          % (EH_plus_W, var, delta))
    print("  diagonal  D_B = E[H]+W - Var = %.2f  <  delta_W = %.2f   =>  (star) FAILS"
          % (Ddiag, delta))
    print("  => a uniform variance bound (from (star dec) at B=B') is NOT curvature.")
    print()


def section_44():
    print("[S4.4] free-field acid test (1D ring, factor-2 block averaging)")
    N = 128            # fine sites (even)
    ncoarse = N // 2
    # free Laplacian on the ring (massless):  (Q phi)_x = 2 phi_x - phi_{x-1} - phi_{x+1}
    Q = np.zeros((N, N))
    for x in range(N):
        Q[x, x] = 2.0
        Q[x, (x - 1) % N] = -1.0
        Q[x, (x + 1) % N] = -1.0
    # block-average factor 2, orthonormal:  (P phi)_I = (phi_{2I} + phi_{2I+1}) / sqrt(2)
    P = np.zeros((ncoarse, N))
    for I in range(ncoarse):
        P[I, 2 * I] = 1.0 / np.sqrt(2.0)
        P[I, 2 * I + 1] = 1.0 / np.sqrt(2.0)

    def gap_and_local(mass):
        Qm = Q + mass ** 2 * np.eye(N)        # SAME regularised action and covariance
        K = np.linalg.inv(Qm)
        Kc = P @ K @ P.T                      # coarse field covariance
        Gpp = np.linalg.inv(Kc)               # G'' = (P K P^T)^{-1}  (coarse precision)
        # free-Gaussian exact identity:  G'' = E[H] - Cov(s,s) + W,  E[H]=P Qm P^T, W=0
        #   =>  Cov(s,s) = P Qm P^T - G''
        Cov_s = P @ Qm @ P.T - Gpp
        idx = np.arange(ncoarse)
        sep = np.minimum(np.abs(idx[:, None] - idx[None, :]),
                         ncoarse - np.abs(idx[:, None] - idx[None, :]))
        return Gpp, Cov_s, sep

    print("  mass m     min eig G'' (=G''(k=0)=m^2)   max|Cov_s|   max|Cov_s(r>4 blocks)|")
    for mass in [0.5, 0.2, 0.1, 0.05, 0.02]:
        Gpp, Cov_s, sep = gap_and_local(mass)
        lam0 = np.linalg.eigvalsh(Gpp).min()
        maxall = np.abs(Cov_s).max()
        maxfar = np.abs(Cov_s[sep > 4]).max()
        print("    %.2f        %14.6e            %10.3e     %10.3e" % (mass, lam0, maxall, maxfar))
    # a decay profile for one mass, to show locality of (star dec)'s kernel
    Gpp, Cov_s, sep = gap_and_local(0.10)
    print("  Cov_s kernel decay profile (m=0.10): separation s -> max|Cov_s| at that s")
    prof = []
    for s in range(0, 12):
        m_ = sep == s
        prof.append(np.abs(Cov_s[m_]).max() if m_.any() else 0.0)
    print("    s:  " + "  ".join("%d" % s for s in range(12)))
    print("    |.: " + "  ".join("%.1e" % v for v in prof))
    print("  => as m -> 0: G''(k=0) = m^2 -> 0  (star FAILS on the free massless field),")
    print("     while Cov(s,s) stays short-range / negligible past a few blocks (star dec HOLDS).")
    print()


if __name__ == "__main__":
    section_41()
    section_42()
    section_43()
    section_44()
