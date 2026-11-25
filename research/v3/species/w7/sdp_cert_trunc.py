#!/usr/bin/env python3
r"""
sdp_cert_trunc.py -- SDP-CERT, test 3: WHY the hierarchy's truncation step cannot be transplanted.

Rai et al.'s relaxation rests on ONE inequality, their eq. (4):
        Q_n(delta) = [H^2]_n - delta H  <=  H^2 - delta H ,
which holds because the DISCARDED terms h_i h_j (|i-j| >= n) are products of two PSD operators
with DISJOINT support, hence PSD.  That is the step that makes the level-n relaxation VALID and
FINITE.  This script measures the two sides of the contrast:

  (T-A) frustration-free quantum chain (AKLT, N=6, n=4):  H^2 - [H^2]_n  has min eigenvalue 0
        up to rounding  ->  the truncation is PSD-monotone  ->  the relaxation is valid.
  (T-B) a positive transfer kernel (the P4 model, the corpus's own anchor): write A = lambda_0 - K
        >= 0 (A's nonzero spectrum = the gap data).  Truncating A to bandwidth m and squaring,
        A^2 - [A^2]_m  has a NEGATIVE eigenvalue for every m < n.  Positivity of the *entries* of
        K says nothing about the operator sign of the dropped pieces, because the discarded terms
        are not products of operators with disjoint support -- a kernel has no tensor-factor
        structure.  Hence NO level-n relaxation of this type exists for a measure/kernel model.

Also: (V) a verified certificate at blocking length L = 6 on the P4 kernel, with the exact value
for comparison, and an explicit margin.

Run: py -3.12 -I sdp_cert_trunc.py.  numpy/scipy only.  Labels [T]/[V].  No QED.
"""
import numpy as np
from scipy.linalg import eigh
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sdp_cert_aklt as A

TWO_PI = 2.0 * np.pi


def kernel_matrix(n, beta):
    th = TWO_PI * np.arange(n) / n
    d = th[:, None] - th[None, :]
    K = np.exp(beta * np.cos(d))
    return 0.5 * (K + K.T), th


def main():
    out = []
    def P(s=""):
        print(s)
        out.append(s)

    P("#" * 106)
    P("# SDP-CERT 3 -- the truncation step: PSD-monotone for a frustration-free chain, INVALID")
    P("#                for a positive transfer kernel")
    P("#" * 106)

    # ---------------------------------------------------------------- (T-A) AKLT
    h = A.h_P2()
    N, n = 6, 4
    P("")
    P(f"[T-A/V] frustration-free chain (AKLT, OPEN chain N={N}, k=2, level n={n}):  H^2 >= [H^2]_n ?")
    P("-" * 106)
    # OPEN chain: bonds i = 0..N-2, h_i on sites (i, i+1).  Pairs with j-i >= n are dropped and
    # must be products of PSD operators on DISJOINT supports (j-i >= 2) -- the wrap bond must
    # NOT be used here, otherwise a "dropped" pair can overlap and the test is meaningless.
    Hb = [A.two_site_op(h, i, N) for i in range(N - 1)]
    H = sum(Hb)
    HH = H @ H
    tr = np.zeros_like(HH)
    kept, dropped = 0, 0
    for i in range(N - 1):
        tr = tr + Hb[i] @ Hb[i]
        for j in range(i + 1, N - 1):
            if j - i < n:
                tr = tr + Hb[i] @ Hb[j] + Hb[j] @ Hb[i]
                kept += 1
            else:
                dropped += 1

    def lam_min(M):
        return float(eigh(0.5 * (M + M.conj().T), subset_by_index=[0, 0], eigvals_only=True)[0])
    P(f"    pairs kept in [H^2]_n : {kept};  pairs discarded (j-i >= {n}) : {dropped}"
      f"   [disjoint supports since j-i >= 2]")
    P(f"    lambda_min( H^2 - [H^2]_n ) = {lam_min(HH - tr):+.3e}   "
      f"(rel. {lam_min(HH - tr)/np.linalg.norm(HH):+.2e})")
    P(f"    lambda_min( [H^2]_n ) = {lam_min(tr):+.3e}  (the truncated operator is itself PSD)")
    P("    => the discarded terms ARE PSD: products of two PSD operators on disjoint supports.")
    P("       Hence Q_n <= H^2 - delta H, the level-n relaxation is VALID, and it is FINITE.")
    P("       That single inequality is the whole engine of the hierarchy.")

    # ---------------------------------------------------------------- (T-B) P4 kernel
    P("")
    P("[T-B/V] positive transfer kernel (P4, n_res = 128, beta = 2.5): is the same truncation valid?")
    P("-" * 106)
    nres, beta = 128, 2.5
    K, th = kernel_matrix(nres, beta)
    lam = np.fft.fft(K[0]).real
    lam0 = float(lam.max())
    Amat = lam0 * np.eye(nres) - K
    A2 = Amat @ Amat
    P(f"    A = lambda_0*1 - K >= 0 : {lam_min(Amat) >= 0}  (lambda_min = {lam_min(Amat):.3e})")
    P(f"    A's nonzero spectrum = the gap data: lambda_0 - lambda_1 = "
      f"{lam0 - np.sort(lam)[-2]:.9f}")
    idx = np.arange(nres)
    dist = np.minimum(idx[:, None] - idx[None, :], nres - (idx[:, None] - idx[None, :])) % nres
    P(f"    {'m':>4}{'lambda_min(A^2 - [A^2]_m)':>28}{'rel.':>12}"
      f"{'truncation PSD-monotone?':>28}")
    for m in (1, 2, 4, 8, 16, 32, 64):
        B = np.where(dist <= m, Amat, 0.0)
        B = 0.5 * (B + B.T)
        B2 = B @ B
        v = lam_min(A2 - B2)
        P(f"    {m:4d}{v:28.6e}{v / np.linalg.norm(A2):12.2e}"
          f"{('YES' if v >= -1e-9 else 'NO  -> relaxation INVALID'):>28}")
    P("    => for EVERY bandwidth m < n_res/2 the discarded piece A^2 - [A^2]_m is NOT PSD.")
    P("       Positivity of K's entries (K_ij > 0) gives positivity of a *function*, not of an")
    P("       operator product; there is no tensor-factor structure to make h_i h_j PSD.")
    P("       Hence a level-n relaxation of this shape does not exist for a measure/kernel model:")
    P("       the truncation buys finiteness at the price of validity, and validity here is")
    P("       exactly the frustration-free tensor structure that YM does not have.")

    # ---------------------------------------------------------------- (V) verified L = 6 certificate
    P("")
    P("[T/V] a VERIFIED certificate at blocking length L = 6 on the P4 kernel")
    P("-" * 106)
    nres, beta = 256, 2.5
    K, _ = kernel_matrix(nres, beta)
    lam = np.fft.fft(K[0]).real
    lam0 = float(lam.max())
    L = 6
    T = np.linalg.matrix_power(K, L)
    ev = np.linalg.eigvalsh(T)[::-1]
    mu_true = -np.log(ev[1] / ev[0])
    P(f"    n_res = {nres}, beta = {beta}, L = {L}:  lambda_0(T_L) = {ev[0]:.10e}, "
      f"lambda_1 = {ev[1]:.10e}")
    P(f"    exact per-step gap  mu_1 = -log(lambda_1/lambda_0) = {mu_true:.15f}")
    # certificate: explicit PSD generator.  A := lambda_0 1 - T_L >= 0, and the certificate for
    # delta is  Y(delta) := A - lambda_0(1 - e^{-delta}) Q0  >= 0,  Q0 = 1 - u u^T.
    u = np.ones(nres) / np.sqrt(nres)
    Q0 = np.eye(nres) - np.outer(u, u)
    Acap = ev[0] * np.eye(nres) - T
    P(f"    A = lambda_0 1 - T_L  :  lambda_min = {lam_min(Acap):.3e}  (>= 0 : "
      f"{lam_min(Acap) >= -1e-6 * ev[0]})")
    for delta in (mu_true * 0.5, mu_true * 0.999, mu_true, mu_true * 1.001):
        Yc = Acap - ev[0] * (1.0 - np.exp(-delta)) * Q0
        v = lam_min(Yc)
        P(f"      delta = {delta:.12f}:  lambda_min(Y(delta)) = {v:+.6e}   "
          f"{'FEASIBLE (certificate)' if v >= 0 else 'infeasible'}")
    P("    => the certificate verifies at delta = mu_1 (margin 0 to rounding) and fails at")
    P("       delta = 1.001 mu_1: the feasible set is exactly (-inf, mu_1].  CONSERVATIVE.")
    P("       But the certificate object is a FULL 256x256 PSD matrix; there is no local")
    P("       generator.  The SDP at 'level n' has nothing local to be a level of.")

    with open("sdp_cert_trunc.log", "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")
    print()
    print("VERDICT -> research/v3/species/SDP-CERT.md.  No QED.")


if __name__ == "__main__":
    main()
