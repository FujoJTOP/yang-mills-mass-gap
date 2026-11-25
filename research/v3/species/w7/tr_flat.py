#!/usr/bin/env python3
r"""
tr_flat.py -- TR species, test 3: the EXACTLY solvable flat-direction model.

This is the cleanest possible realisation of the P4 anchor "S_W is exactly flat on the constant
holonomy slice": a circulant kernel on the ring Z_n,

    K_ij = exp(beta cos(theta_i - theta_j)),      theta_k = 2 pi k / n ,

i.e. the Wilson weight of a pure *difference* (no frustration, unique... but n-fold degenerate
maximiser).  Everything is exact:
    spectrum  lambda_j = sum_m e^{beta cos(theta_m)} e^{2 pi i j m / n}   (= 2 pi I_j(beta) + O(1/n) )
    gap       Delta(beta) = -log(lambda_1/lambda_0)          [exact, n-independent to O(1/n)]
    tropical  W_ij = cos(theta_i-theta_j), rho_inf = 1, critical graph = ALL n self-loops
              => tropical Delta = 0 EXACTLY, for every beta, every n.
    squeeze   width = (1/beta) log sum_j e^{beta(cos theta_j - 1)} ~ (log n)/beta + O(1)

So: the tropical data predict Delta = 0 while the exact gap is O(1); the gap is entirely the
Maslov / Van Vleck (prefactor) term; and the squeeze width grows like log(n)/beta, i.e. it is NOT
uniform in the state-space resolution (equivalently: not uniform in (L,a)).

Also reports the large-beta asymptotics Delta ~ 1/(2 beta) and how far it is from the truth at
beta = 2.5 (the expansion parameter 1/beta = 0.4 is not small).

Run: py -3.12 -I tr_flat.py.  numpy/scipy only, deterministic.  Labels [T]/[V].  No QED.
"""
import numpy as np
from scipy.special import iv

TWO_PI = 2.0 * np.pi


def circulant_row(n, beta):
    m = np.arange(n)
    return np.exp(beta * np.cos(TWO_PI * m / n))


def spectrum(n, beta):
    row = circulant_row(n, beta)
    lam = np.fft.fft(row).real          # eigenvalues of the circulant (real, symmetric kernel)
    lam = np.sort(lam)[::-1]
    return lam


def main():
    print("#" * 100)
    print("# TR test 3 -- EXACT flat-direction (circulant) model: tropical Delta = 0 vs exact Delta > 0")
    print("#" * 100)

    # ---- validation: circulant spectrum vs 2 pi I_j(beta)
    n, beta = 512, 2.5
    lam = spectrum(n, beta)
    print(f"[check] n={n}, beta={beta}: lambda_j/n vs 2pi I_j(beta)")
    for j in (0, 1, 2, 3):
        print(f"    j={j}: lambda_j/n = {lam[j]/n:.10f}   2pi I_j = {TWO_PI*iv(j, beta):.10f}"
              f"   rel.err = {abs(lam[j]/n/(TWO_PI*iv(j, beta))-1):.2e}")

    # ---- exact gap vs beta (n large enough that it is n-independent)
    print("\n" + "=" * 100)
    print("[T/V] EXACT gap Delta(beta) = -log(l1/l0), the tropical prediction (0), the squeeze width,")
    print("      and the large-beta asymptote 1/(2 beta):")
    print("-" * 100)
    print(f"    {'beta':>6}{'l1/l0':>12}{'Delta_exact':>13}{'Delta_trop':>12}{'width':>10}"
          f"{'width/Delta':>13}{'1/(2beta)':>11}{'asym/exact':>12}")
    n = 4096
    for beta in (0.5, 1.0, 2.5, 5.0, 10.0, 20.0):
        lam = spectrum(n, beta)
        ratio = lam[1] / lam[0]
        D = -np.log(ratio)
        # exact tropical width for this beta: (1/beta) log sum_j e^{beta(cos theta_j - 1)}
        th = TWO_PI * np.arange(n) / n
        wid = np.log(np.exp(beta * (np.cos(th) - 1.0)).sum()) / beta
        print(f"    {beta:6.2f}{ratio:12.6f}{D:13.6f}{0.0:12.6f}{wid:10.4f}{wid/D:13.3f}"
              f"{1/(2*beta):11.4f}{(1/(2*beta))/D:12.3f}")

    # ---- n-dependence of the squeeze width at beta=2.5 (must be uniform in the resolution)
    print("\n" + "=" * 100)
    print("[V] SQUEEZE WIDTH vs state-space resolution n at beta=2.5 (the exact Delta is n-independent):")
    print("-" * 100)
    beta = 2.5
    print(f"    {'n':>6}{'Delta_exact':>13}{'width':>10}{'width/Delta':>13}{'log(n)/beta':>13}")
    for n in (16, 32, 64, 128, 256, 512, 1024, 4096):
        lam = spectrum(n, beta)
        D = -np.log(lam[1] / lam[0])
        th = TWO_PI * np.arange(n) / n
        wid = np.log(np.exp(beta * (np.cos(th) - 1.0)).sum()) / beta
        print(f"    {n:6d}{D:13.6f}{wid:10.4f}{wid/D:13.3f}{np.log(n)/beta:13.4f}")

    # ---- why: the critical graph is the whole diagonal (n-fold degeneracy), for every beta
    print("\n" + "=" * 100)
    print("[T] STRUCTURE: W = cos(theta_i-theta_j) has W_ii = 1 for ALL i => critical graph = entire")
    print("    diagonal (n-fold degenerate) => tropical spectral gap = 0 => tropical Delta = 0.")
    print("    The exact gap is carried entirely by the PREFACTOR: the Gaussian width 1/sqrt(beta)")
    print("    (Van Vleck) times the compact measure.  In the non-compact (free) case the same")
    print("    computation gives a CONTINUOUS spectrum with no gap at all -- the discretisation/")
    print("    compactness is what creates the gap, and compactness is a measure, not an action,")
    print("    property.  Neither is a tropical datum.")
    print()
    print("VERDICT -> TR-TROPICAL.md.  No QED.")


if __name__ == "__main__":
    main()
