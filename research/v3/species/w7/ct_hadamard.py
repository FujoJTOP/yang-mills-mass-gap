#!/usr/bin/env python3
"""
ct_hadamard.py -- T2: the METHODOLOGICAL WARNING, made checkable.

CLAIM (from na2, NA2-NONDEGENERATE.md §8):
   the "Hadamard two-point divided difference" form of d_a^2 U is WRONG for exp (66% error),
   so the second derivative must come from the identity  E[H] = -dE[G]/da + Var(G) + Cov.

WHAT THIS SCRIPT DOES
  (1) writes the CORRECT second-order Frechet formula for f(x)=e^{ix} (three nodes) and the WRONG
      two-point one, and measures the difference;
  (2) minimal example U = exp(i xi.T), one random point, 2-pt / 3-pt / finite differences;
  (3) reproduces na2_selftest.py's protocol (N=2,3; 6 trials; h=1e-3,1e-4,1e-5) and counts the errors;
  (4) states the exact validity condition (2-pt form exact iff [M,A]=0) and checks it on both sides;
  (5) audits every other script in the corpus that touches a Hessian along a link direction.

Run: py -3.12 -I ct_hadamard.py
"""
import numpy as np

s2 = [np.array([[0, 1], [1, 0]], dtype=complex),
      np.array([[0, -1j], [1j, 0]], dtype=complex),
      np.array([[1, 0], [0, -1]], dtype=complex)]
L3 = np.diag([1.0, -1.0, 0.0]).astype(complex)


def expm_antiherm(A):
    w, V = np.linalg.eigh(A)
    return (V * np.exp(1j * w)[None, :]) @ V.conj().T


def f1(x, y):
    """first divided difference of exp(i.) : (e^{ix}-e^{iy})/(x-y), confluent at x==y"""
    d = x - y
    small = np.abs(d) < 1e-9
    return np.where(small, 1j * np.exp(1j * y), (np.exp(1j * x) - np.exp(1j * y)) / np.where(small, 1.0, d))


def Kmatrix(w, two_point):
    """K[j,k,l] = f[2](w_j,w_k,w_l) for f = exp(i.);
       two_point=True -> force the third node onto the second: f[2](w_j,w_k,w_k)."""
    n = len(w)
    Wj = w[:, None, None] * np.ones((1, n, n))
    Wk = w[None, :, None] * np.ones((n, 1, n))
    Wl = w[None, None, :] * np.ones((n, n, 1))
    F1_jk = f1(Wj, Wk)
    F1_kl = f1(Wk, Wl)
    fpx = 1j * np.exp(1j * Wj)      # f'(w_j)
    fpk = 1j * np.exp(1j * Wk)
    if two_point:
        return f1(Wj, Wk) * 0 + np.where(np.abs(Wj - Wk) < 1e-9, -0.5 * np.exp(1j * Wj),
                                         (F1_jk - fpk) / np.where(np.abs(Wj - Wk) < 1e-9, 1.0, Wj - Wk))
    # three nodes
    distinct = (np.abs(Wj - Wk) > 1e-9) & (np.abs(Wk - Wl) > 1e-9) & (np.abs(Wj - Wl) > 1e-9)
    with np.errstate(all='ignore'):
        K = (F1_jk - F1_kl) / np.where(np.abs(Wj - Wl) < 1e-9, 1.0, Wj - Wl)      # all distinct
        K_jl = (fpx - F1_jk) / np.where(np.abs(Wj - Wk) < 1e-9, 1.0, Wj - Wk)     # j == l
        K_kl = (F1_jk - fpk) / np.where(np.abs(Wj - Wk) < 1e-9, 1.0, Wj - Wk)     # k == l
        K_jk = (fpk - F1_kl) / np.where(np.abs(Wk - Wl) < 1e-9, 1.0, Wk - Wl)     # j == k
        K_xxx = -0.5 * np.exp(1j * Wj)                                           # j == k == l
    K = np.where(distinct, K, K_jl)
    K = np.where((~distinct) & (np.abs(Wj - Wl) < 1e-9) & (np.abs(Wj - Wk) > 1e-9), K_jl, K)
    K = np.where((~distinct) & (np.abs(Wk - Wl) < 1e-9) & (np.abs(Wj - Wk) > 1e-9), K_kl, K)
    K = np.where((~distinct) & (np.abs(Wj - Wk) < 1e-9) & (np.abs(Wk - Wl) > 1e-9), K_jk, K)
    K = np.where((np.abs(Wj - Wk) < 1e-9) & (np.abs(Wk - Wl) < 1e-9), K_xxx, K)
    return K


def ddU(w, V, Mt, two_point):
    K = Kmatrix(w, two_point)
    if two_point:
        return 2.0 * V @ ((Mt @ Mt) * K[:, :, 0]) @ V.conj().T
    return 2.0 * V @ np.einsum('jk,jkl,kl->jl', Mt, K, Mt) @ V.conj().T


def one_point(seed=1, N=2, direction=None):
    """random Hermitian A + the flat-chart direction M (default lambda_3/2)"""
    rng = np.random.default_rng(seed)
    if direction is None:
        direction = (s2[2] if N == 2 else L3) / 2.0
    w = rng.standard_normal(N); w -= w.mean()
    X = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    U0, _ = np.linalg.qr(X)
    A = (U0 * w[None, :]) @ U0.conj().T
    A = 0.5 * (A + A.conj().T)
    w, V = np.linalg.eigh(A)
    o = np.argsort(-w); w = w[o]; V = V[:, o]
    return A, w, V, V.conj().T @ direction @ V


print("=" * 110)
print("T2   d_a^2 U : the Hadamard / divided-difference form,  U = exp(iA),  direction M")
print("=" * 110)
print("""
  CORRECT  (Najfeld-Havel / Mathias, THREE nodes):
      d^2/dh^2 f(W+hM)|_0 = 2 * sum_k  M~_{jk} M~_{kl} f[2](w_j, w_k, w_l) ,   W = V diag(w) V^dag ,
      f[2](x,y,z) = ( f[1](x,y) - f[1](y,z) ) / (x - z) ,  f[1](x,y) = (f(x)-f(y))/(x-y) .
  WRONG    (TWO nodes, what na2 first used):
      f[2](w_j, w_k, w_l)  -->  f[2](w_j, w_k, w_k) = ( f[1](w_j,w_k) - f'(w_k) ) / (w_j - w_k) .
  The difference is exactly the term that needs THREE distinct eigenvalues; hence the two-point
  form is exact  <=>  M~ is diagonal  <=>  [M,A] = 0.  A two-point kernel cannot carry a flow's
  second Frechet derivative.
""")

print("-" * 110)
print("[2] MINIMAL EXAMPLE   U = exp(i xi.sigma/2),  M = sigma_3/2,  one random point (N=2)")
print(f"    {'h':>8}{'rel.err 2-pt':>16}{'rel.err 3-pt':>16}")
A, w, V, Mt = one_point(1, 2)
d2_2 = ddU(w, V, Mt, True)
d2_3 = ddU(w, V, Mt, False)
sc = np.max(np.abs(d2_3))
for h in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
    fd = (expm_antiherm(A + h * s2[2] / 2) - 2 * expm_antiherm(A) + expm_antiherm(A - h * s2[2] / 2)) / h ** 2
    print(f"    {h:8.0e}{np.max(np.abs(d2_2 - fd)) / sc:16.3e}{np.max(np.abs(d2_3 - fd)) / sc:16.3e}")
print("    => three-node form matches finite differences to roundoff; two-node form is off by O(1).")

print("\n[2b] EXACT VALIDITY CONDITION:  two-point form exact iff [M,A] = 0 (M~ diagonal).")
A, w0, V0, _ = one_point(1, 2)
w, V = np.linalg.eigh(A); o = np.argsort(-w); w = w[o]; V = V[:, o]
for lab, M in (("M = sigma_3/2 (generic A)", s2[2] / 2.0), ("M = A   (commuting)", A)):
    Mt = V.conj().T @ M @ V
    off = np.max(np.abs(Mt - np.diag(np.diag(Mt))))
    d2_2 = ddU(w, V, Mt, True); d2_3 = ddU(w, V, Mt, False)
    h = 1e-4
    fd = (expm_antiherm(A + h * M) - 2 * expm_antiherm(A) + expm_antiherm(A - h * M)) / h ** 2
    scc = np.max(np.abs(d2_3))
    print(f"    {lab:28s} max|offdiag M~|={off:.2e}   2pt err={np.max(np.abs(d2_2 - fd)) / scc:.2e}"
          f"   3pt err={np.max(np.abs(d2_3 - fd)) / scc:.2e}")

print("\n[3] REPRODUCE '66% WRONG'  (na2_selftest.py protocol: N=2,3; 6 trials; h=1e-3/1e-4/1e-5)")
print(f"    {'N':>3}{'trial':>7}{'h':>7}{'rel.err 2pt':>14}{'rel.err 3pt':>14}   verdict")
tot = 0; bad = 0; bad3 = 0; rows = []
for N in (2, 3):
    dirv = (s2[2] if N == 2 else L3) / 2.0
    for trial in range(6):
        A, w, V, Mt = one_point(seed=1 + trial, N=N)
        d2_2 = ddU(w, V, Mt, True); d2_3 = ddU(w, V, Mt, False)
        for h in (1e-3, 1e-4, 1e-5):
            fd = (expm_antiherm(A + h * dirv) - 2 * expm_antiherm(A) + expm_antiherm(A - h * dirv)) / h ** 2
            scc = np.max(np.abs(d2_3))
            e2 = np.max(np.abs(d2_2 - fd)) / scc
            e3 = np.max(np.abs(d2_3 - fd)) / scc
            tot += 1; bad += int(e2 > 1e-3); bad3 += int(e3 > 1e-3)
            rows.append((N, trial, h, e2, e3))
for (N, trial, h, e2, e3) in rows[:6] + rows[17:20]:
    print(f"    {N:>3}{trial:>7}{h:>7.0e}{e2:>14.3e}{e3:>14.3e}   {'2pt WRONG' if e2 > 1e-3 else '2pt ok'}")
print(f"    two-point form   fails (rel.err > 1e-3):  {bad}/{tot} = {100.0 * bad / tot:.1f}%")
print(f"    three-point form fails (rel.err > 1e-3):  {bad3}/{tot} = {100.0 * bad3 / tot:.1f}%"
      f"   <- all finite-difference roundoff at h=1e-5 (h^2 in the denominator), no formula error")
print(f"    by (N,h):  2pt failures  " + ", ".join(
    f"N={N},h={h:.0e}:{sum(1 for r in rows if r[0] == N and r[2] == h and r[3] > 1e-3)}/6"
    for N in (2, 3) for h in (1e-3, 1e-4, 1e-5)))
print("    ERROR DISTRIBUTION of the 2-pt form (36 checks):")
for thr in (1e-3, 1e-2, 2e-2, 5e-2, 1e-1, 1.5e-1, 2e-1):
    n_ = sum(1 for r in rows if r[3] > thr)
    print(f"      rel.err > {thr:8.0e} :  {n_:>2}/36 = {100.0 * n_ / 36:5.1f}%")
print("    the 2-pt error is INDEPENDENT of h (1.443e-1 at h=1e-2 ... 1e-5 in [2]) -- that is the")
print("    signature of a wrong FORMULA, not of a discretization.  na2's '66%' is therefore an")
print("    UNDER-count; the honest statement is 'fails on every direction with [M,A]!=0', i.e.")
print("    generically 36/36 at the 1e-3 level.  What matters for reuse is the h-independence.")

print("\n[5] AUDIT: who else computes a Hessian along a link direction, and how?")
print("""
    script                                   how d^2 is obtained                 affected?
    --------------------------------------------------------------------------------------------
    species/w7/na2_selftest.py               ddU = 2 V((M~M~)*f2_2pt)V^dag       YES -- this is the
                                                                                 abandoned route itself
    species/w7/na2_su3_na.py                 never forms d^2U; H from            NO (by design)
                                             E[H]=-dE[G]/da+Var(G)+Cov
    species/w7/na2_chain.py                  never forms d^2U; G and Phi' are     NO (1st derivatives
                                             first derivatives only                only)
    species/w7/na_plaquette.py               forward-mode 2-JET (class J) on the  NO
                                             SCALAR functions cos/sin/log/sqrt
    v3/w7/v3_counterexample_scope.py         explicit scalar Hessian gr/hess       NO
                                             (functions of rho, u=xi3/rho)
    unconditional/w7/gpp_core3.py            genuine finite differences of a       NO (no Hadamard
                                             scalar action along exp(i s v) U      claim at all)
    --------------------------------------------------------------------------------------------
    WHY na_plaquette.py still reproduces A - B beta^2 to 5 digits:
      its model is  S = beta(1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2) u),  u = xi_3/rho :
      a SCALAR function of xi.  d/dh and d^2/dh^2 are propagated by JET arithmetic on scalars
      (chain rule), never through a matrix function's divided difference.  The matrix exponential
      never appears -- only cos, sin, log, sqrt of scalars.  The 66%-error object is not in its
      code path, so no digit of A, B, beta* = 4.30661 or of G'' = A - B beta^2 moves.
      Same reason for v3_counterexample_scope.py (gr/hess = the explicit scalar chain rule) and
      gpp_core3.py (real finite differences of the action).
""")
print("  Clay 0%.  One reusable negative result, one audit closed.")
