#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gy3_free_acid.py -- GY3 task (a): the free-field counterexample family, closed form.

The no-light-mode theorem ("no family with m1^2(L) <~ L^-2 and f1^2 >= r0 > 0",
HN1-H1.md (c-3)) must FAIL on the free massless field (acid F1, Q6-H1-STATUS.md:298-300).
This script computes the exact free-field KL atom tables that realize the failure:

  [1] box-gap acid family (single-particle channel, the corpus's free-massless row,
      M1-4D-RUN.md:89):  nu_L = f1^2 * delta_{m1^2},  m1^2 = lambda_min = 4 sin^2(pi/L)
      (a=1 lattice), f1^2 = 1 (unit carrier), chi_L = 1/m1^2 ~ L^2/(4 pi^2).
      Exact numbers + ratios; cross-check vs M1:89 (0.500 / 1.707 / 6.569 at L=4/8/16).
  [2] exact 4D lattice Gaussian verification (L=4, 256-dim, zero mode removed):
      Var(phi~_{k_min}) = 1/lambda_min; massive Var(a) = 1/m^2 independent of L.
  [3] HN1 two-atom family nu = (1/2)delta_{1/L^2} + (1/2)delta_4: chi_L = L^2/2 + 1/8
      (hn1h1_check.log [4]) -- exact arithmetic re-verified.
  [4] the 0^{++} F^2-analog composite channel (free scalar :phi^2:, cylinder T^3_L x R):
      space-summed carrier: atoms at (2 omega_k)^2, f1^2 = 3/omega_1^2 ~ 3L^2/(4pi^2) -> inf,
      chi = (1/8) sum 1/omega_k^4 ~ L^4 (IR);  point carrier: f1^2 = 18/(L^6 omega_1^2)
      ~ 4.5/(pi^2 L^4) -> 0.  Normalization hazard made explicit.

py -3.12 -I  deterministic, no MC.  Log: w7/gy3_free_acid.log
"""
import numpy as np

out = []
def P(*a):
    s = " ".join(str(x) for x in a)
    out.append(s); print(s, flush=True)

# ---------------------------------------------------------------- [1] box gap
P("=" * 72); P("[1] box-gap acid family (single-particle channel, a=1 lattice)")
P("    m1^2 = lambda_min = 4 sin^2(pi/L);  nu_L = 1 * delta_{m1^2};  chi_L = 1/m1^2")
P("=" * 72)
P(f"{'L':>4} | {'m1^2=4sin^2(pi/L)':>18} | {'m1^2*L^2/(4pi^2)':>16} | {'chi_L=1/m1^2':>12} | {'chi_L/L^2':>10}")
rows = {}
for L in (4, 6, 8, 12, 16, 32, 64):
    m1sq = 4 * np.sin(np.pi / L) ** 2
    chi = 1.0 / m1sq
    rows[L] = (m1sq, chi)
    P(f"{L:>4} | {m1sq:>18.10f} | {m1sq * L**2 / (4*np.pi**2):>16.10f} | {chi:>12.6f} | {chi / L**2:>10.7f}")
P("")
P("chi_L ratios vs the free L^2 law and vs the continuum limit (b/a)^2:")
for (a, b) in ((4, 6), (4, 8), (6, 8), (8, 12), (8, 16), (16, 32)):
    r = rows[b][1] / rows[a][1]
    cont = (b / a) ** 2
    P(f"  L={a}->{b}: chi ratio = {r:.6f}   continuum L^2 ratio = {cont:.6f}   "
      f"lattice/cont = {r/cont:.6f}")
P("")
P("M1-4D-RUN.md:89 free table check (massless Var at L=4/8/16 = chi_L):")
for L, ref in ((4, 0.500), (8, 1.707), (16, 6.569)):
    d = rows[L][1] - ref
    P(f"  L={L}: chi_L = {rows[L][1]:.6f}  vs M1 table {ref}  diff = {d:.2e}")
P("")
P("Continuum box gap: m1^2 = (2 pi/L)^2; lattice/continuum = (L sin(pi/L)/pi)^2 -> 1 - pi^2/(3L^2):")
for L in (4, 8, 16, 64):
    P(f"  L={L}: (L sin(pi/L)/pi)^2 = {(L*np.sin(np.pi/L)/np.pi)**2:.10f}   "
      f"1 - pi^2/(3L^2) = {1 - np.pi**2/(3*L**2):.10f}")
P("")
P("Transfer-matrix energy version (A1-FREE-TORUS.md:29, exact discrete time):")
P("  m1 = (2/a) arcsinh(sin(pi/L)) at a=1  =>  m1^2 = 4 arcsinh^2(sin(pi/L)):")
for L in (4, 8, 16):
    tm = 4 * np.arcsinh(np.sin(np.pi / L)) ** 2
    P(f"  L={L}: 4 arcsinh^2(sin(pi/L)) = {tm:.10f}   (Euclidean lambda_min = {rows[L][0]:.10f}, "
      f"ratio {tm/rows[L][0]:.6f} -> 1)")

# ------------------------------------------------- [2] exact lattice Gaussian
P(""); P("=" * 72)
P("[2] exact 4D lattice Gaussian checks (deterministic linear algebra)")
P("=" * 72)
def laplacian(L):
    """precision matrix of free scalar, action S = sum_x sum_mu (phi_{x+mu}-phi_x)^2
    (each undirected edge counted once in the sum, site degree 8)."""
    n = L ** 4
    A = np.zeros((n, n))
    idx = np.arange(n).reshape(L, L, L, L)
    lin = idx.ravel()
    for mu in range(4):
        nbr = np.roll(idx, -1, axis=mu).ravel()
        for i, j in zip(lin, nbr):
            A[i, j] += -1.0; A[j, i] += -1.0
    A[np.arange(n), np.arange(n)] += 8.0
    return A

for L in (4, 6):
    n = L ** 4
    A = laplacian(L)
    lam, V = np.linalg.eigh(A)
    lam0 = lam[0]
    assert abs(lam0) < 1e-8, lam0
    keep = lam > 1e-8
    Cov = ((V[:, keep] / lam[keep]) @ V[:, keep].T)
    lam_min = np.sort(lam)[1]
    P(f"  L={L}: spec(A) min nonzero = {lam_min:.12f}"
      f"   closed form 4 sin^2(pi/L) = {4*np.sin(np.pi/L)**2:.12f}"
      f"   diff = {abs(lam_min - 4*np.sin(np.pi/L)**2):.2e}")
    # Var(phi~_{k_min}): mode vector e = N^{-1/2} exp(i k.x); k_min = 2pi/L e_1
    kk = 2 * np.pi / L
    x = np.indices((L, L, L, L)).reshape(4, n)
    # complex mode phi~_k = N^{-1/2} sum e^{i k.x} phi_x ; Var = < |phi~_k|^2 > = 1/lambda_k
    e = np.exp(1j * kk * x[0]) / np.sqrt(n)
    v = float(np.real(e.conj() @ Cov @ e))
    P(f"  L={L}: Var(phi~_kmin) exact = {v:.12f}   1/lambda_min = {1/(4*np.sin(np.pi/L)**2):.12f}"
      f"   diff = {abs(v - 1/(4*np.sin(np.pi/L)**2)):.2e}")
    # massive: Var(a) = (1/N) 1^T (A + m^2)^{-1} 1 = 1/m^2  (L-independent)
    for m in (0.3, 0.5):
        Am = A + (m ** 2) * np.eye(n)
        one = np.ones(n)
        va = float(one @ np.linalg.solve(Am, one)) / n
        P(f"  L={L}, m={m}: Var(a) = {va:.12f}   1/m^2 = {1/m**2:.12f}   diff = {abs(va-1/m**2):.2e}")

# ---------------------------------------------------- [3] HN1 two-atom family
P(""); P("=" * 72)
P("[3] HN1-H1.md (c-1) family: nu_L = (1/2) delta_{1/L^2} + (1/2) delta_4")
P("    chi_L = int s^-1 dnu = L^2/2 + 1/8   (hn1h1_check.log [4])")
P("=" * 72)
ok = True
for L, ref in ((2, 17 / 8), (4, 65 / 8), (8, 257 / 8), (16, 1025 / 8), (64, 16385 / 8)):
    val = L ** 2 / 2 + 1 / 8
    ok &= abs(val - ref) < 1e-12
    P(f"  L={L:>3}: chi_L = L^2/2+1/8 = {val:.6f}   log ref = {ref:.6f}   match = {abs(val-ref)<1e-12}")
P(f"  all match: {ok}")
P("  m1^2 = 1/L^2 (inf supp), f1^2 = 1/2 >= r0, chi_L ~ L^2/2 unbounded")
P("  squeeze: f1^2/chi = (1/2)/(L^2/2+1/8) < 1/L^2 = m1^2  (strict, equality-free)  [E]")

# --------------------------------------- [4] F^2-analog composite (continuum)
P(""); P("=" * 72)
P("[4] free-scalar 0++ composite channel :phi^2: (F^2 analog), cylinder T^3_L x R")
P("    omega_n = sqrt((2 pi |n|/L)^2 + m^2), n in Z^3")
P("    space-summed carrier: C(t) = (1/2) sum_{k!=0} e^{-2 omega_k t} / omega_k^2")
P("        atoms at s_k = (2 omega_k)^2, weight w_k = 1/(2 omega_k^2) (per k in the sum)")
P("    point carrier: C(t) = (1/(2L^6)) sum_{k,k'} e^{-(omega_k+omega_k')t}/(omega_k omega_k')")
P("=" * 72)
def shell_mult(L, m2, r2max=8):
    """mult(n) for |n|^2 = r2 with omega = sqrt((2pi/L)^2 r2 + m2); r2 <= r2max."""
    shells = {}
    R = int(np.ceil(np.sqrt(r2max))) + 1
    for i in range(-R, R + 1):
        for j in range(-R, R + 1):
            for k in range(-R, R + 1):
                r2 = i * i + j * j + k * k
                if 0 < r2 <= r2max:
                    shells[r2] = shells.get(r2, 0) + 1
    return shells

sh = shell_mult(8, 0.0)
P(f"  |n|^2 shell multiplicities in Z^3 (r2<=8): " +
  ", ".join(f"|n|^2={r}:{c}" for r, c in sorted(sh.items())))

for m in (0.0, 0.1):
    P(f"  --- regulator m = {m} ---")
    P(f"  {'L':>4} | {'omega_1':>10} | {'m1^2=(2w1)^2':>13} | {'m1^2*L^2/(16pi^2)':>18} | "
      f"{'f1^2 (space-summed)':>20} | {'f1^2 (point)':>13}")
    for L in (4, 6, 8, 16, 32):
        w1 = np.sqrt((2 * np.pi / L) ** 2 + m ** 2)
        m1sq = (2 * w1) ** 2
        f_space = 3.0 / w1 ** 2                    # 6 * 1/(2 w1^2)
        f_point = 18.0 / (L ** 6 * w1 ** 2)        # 36 * 1/(2 L^6 w1^2)
        P(f"  {L:>4} | {w1:>10.6f} | {m1sq:>13.6f} | {m1sq*L**2/(16*np.pi**2):>18.8f} | "
          f"{f_space:>20.6e} | {f_point:>13.6e}")
P("")
P("  space-summed: f1^2 = 3/omega_1^2 ~ 3 L^2/(4 pi^2) -> inf  (ATOM, growing)")
P("  point:        f1^2 ~ 4.5/(pi^2 L^4) -> 0                (atom vanishes)")
P("  chi (space-summed) = (1/8) sum_{k!=0} 1/omega_k^4:  m=0 exact closed form:")
# sum_{n in Z^3, n!=0} 1/|n|^4  (Epstein zeta Z_3(2) = 16.53231596...)
R = 120
tot = 0.0
ij2s = []
for i in range(-R, R + 1):
    for j in range(-R, R + 1):
        ij2 = i * i + j * j
        if ij2 <= 4 * R * R:
            ks = np.arange(-R, R + 1) ** 2 + ij2
            msk = (ks > 0) & (ks <= 4 * R * R)
            tot += float(np.sum(1.0 / ks[msk] ** 2.0))
# tail: 4 pi int_R^inf r^2/r^4 dr = 4 pi / R
ep = tot + 4 * np.pi / R
P(f"    sum_(n!=0) 1/|n|^4  (cutoff {R}, tail-corrected) ~ {ep:.6f}   "
  f"(Epstein Z_3(2) reference 16.53231596, diff {abs(ep-16.53231596):.4f})")
for L in (4, 8, 16):
    chi0 = (1.0 / 8.0) * (L ** 4 / (16 * np.pi ** 4)) * ep
    P(f"    L={L}: chi = (1/8) L^4/(16 pi^4) * S4 = {chi0:.4f}   (chi/L^4 = {chi0/L**4:.6e})")
P("    with m>0 regulator, numeric sum (k-grid n in Z^3, |n|<=90):")
for m in (0.05, 0.1):
    for L in (4, 8):
        Rn = 90
        ns = np.arange(-Rn, Rn + 1)
        II, JJ, KK = np.meshgrid(ns, ns, ns, indexing="ij")
        r2 = II ** 2 + JJ ** 2 + KK ** 2
        w = np.sqrt((2 * np.pi / L) ** 2 * r2 + m ** 2)
        mask = r2 > 0
        chi_m = float(np.sum(1.0 / (8.0 * w[mask] ** 4)))
        P(f"    m={m}, L={L}: chi = {chi_m:.6f}")
P("")
P("VERDICT (a): the theorem-form must exclude (i) the single-particle box-mode carrier")
P("  [m1^2 = 4 sin^2(pi/L), f1^2 = 1, chi ~ L^2/4pi^2] and (ii) the HN1 two-atom family;")
P("  the free F^2-analog composite evades the theorem in the density normalization")
P("  (f1^2 -> 0) but violates it in the extensive (space-summed) normalization")
P("  (f1^2 ~ 3L^2/4pi^2 -> inf, chi ~ L^4).  Carrier normalization is load-bearing.")

with open("D:/Dev/math/research/v3/w7/gy3_free_acid.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(out) + "\n")
print("\n# log written: w7/gy3_free_acid.log", flush=True)
