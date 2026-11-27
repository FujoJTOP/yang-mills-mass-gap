#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gp9rb_form.py  --  GP9-RB : what is epsilon_a^form(H_a - H) in the Gaussian /
lattice-Laplacian sector?

Deterministic, pure numpy, no Monte Carlo, no random seeds.
Writes gp9rb_form.log (caller redirects stdout).

Sections
 [1] exact closed form  eps_free = 1 - 4/pi^2  (d-independent), numeric cross-check
 [2] attainment / weighted-average identity
 [3] massive case: eps_a(m) -> 1-4/pi^2 with O(a^2 m^2);  is it O(a^k), k>0 ?  NO
 [4] non-invariance under the choice of lattice kinetic symbol h
 [5] Shale-Stinespring / quasi-equivalence diagnostics
 [6] GO2 family counting function vs a physical (lattice) family
"""

import numpy as np

SEP = "-" * 78


def h_wilson(k):
    """Nearest-neighbour (Wilson) lattice Laplacian symbol, in units a=1.
    h_W(k) = 4 sum_i sin^2(k_i/2) = 2 sum_i (1 - cos k_i)."""
    return 2.0 * np.sum(1.0 - np.cos(k), axis=-1)


def ksq(k):
    return np.sum(k * k, axis=-1)


def sup_ratio_grid(hfun, d, n, extra=1.0, m2=0.0):
    """sup over the BZ of |h(k)-|k|^2| / (|k|^2 + m2), on an n^d grid
    (the grid always contains the corners |k_i| = pi, where the sup sits).
    `extra` multiplies hfun (for coupling-coefficient bookkeeping)."""
    ax = np.linspace(-np.pi, np.pi, n)
    grids = np.meshgrid(*([ax] * d), indexing="ij")
    k = np.stack([g.ravel() for g in grids], axis=-1)
    num = np.abs(extra * hfun(k) - ksq(k))
    den = ksq(k) + m2
    with np.errstate(divide="ignore", invalid="ignore"):
        # guard: |k|^2 must be resolvable, otherwise the ratio is 0/0 noise
        r = np.where(den > 1e-10, num / den, 0.0)
    return float(np.max(r))


def main():
    print("GP9-RB : epsilon_a^form(H_a - H), Gaussian sector")
    print("symbols from GN2-VARIATIONAL.md:52 : omega_a(p) = (2/a^2) sum_mu (1-cos a p_mu), |p|^2")
    print(SEP)

    # ------------------------------------------------------------------ [1]
    print("[1] exact closed form (d-independent)")
    print("    r(k) = 1 - 4 sin^2(k_i/2)/k_i^2 per axis; ratio is the w_i-weighted")
    print("    average of the per-axis r_i with w_i = k_i^2/sum_j k_j^2, hence")
    print("    ratio <= max_i r_i <= 1 - 4/pi^2, attained iff every k_i in {0,+-pi}.")
    c_exact = 1.0 - 4.0 / np.pi ** 2
    print("    1 - 4/pi^2 = %.16f" % c_exact)
    for d in (1, 2, 3, 4):
        n = 20001 if d == 1 else (601 if d == 2 else (121 if d == 3 else 41))
        s = sup_ratio_grid(h_wilson, d, n)
        print("    d=%d  numeric sup (n=%d) = %.12f   |sup-(1-4/pi^2)| = %.2e"
              % (d, n, s, abs(s - c_exact)))
    # corner value, exact by hand
    print("    corner k=(pi,...,pi): |4d - d pi^2|/(d pi^2) = %.16f  (same for every d)"
          % (abs(4.0 - np.pi ** 2) / np.pi ** 2))
    print("    => eps_free = 1 - 4/pi^2 = %.10f  in (0,1), a-INDEPENDENT (O(1), not O(a^k))"
          % c_exact)

    # ------------------------------------------------------------------ [2]
    print(SEP)
    print("[2] attainment and the weighted-average identity (numeric)")
    ax = np.linspace(-np.pi, np.pi, 201)
    K = np.stack(np.meshgrid(ax, ax, ax, indexing="ij"), axis=-1).reshape(-1, 3)
    num = np.abs(h_wilson(K) - ksq(K))
    r = np.where(ksq(K) > 0, num / ksq(K), 0.0)
    print("    max over 201^3 grid           : %.12f" % r.max())
    print("    argmax representative        : %s" % np.round(K[np.argmax(r)], 6))
    # weighted-average check on a few random-but-fixed points
    for kk in ([2.9, 0.4, 0.1], [1.5, 1.5, 1.5], [3.0, 3.0, 3.0], [0.05, 0.05, 0.05]):
        kk = np.array(kk)
        w = kk ** 2 / np.sum(kk ** 2)
        ri = 1.0 - 4.0 * np.sin(kk / 2.0) ** 2 / kk ** 2
        print("    k=%-18s ratio=%.12f  sum_i w_i r_i=%.12f"
              % (np.round(kk, 3), np.abs(h_wilson(kk) - ksq(kk)) / ksq(kk),
                 float(np.sum(w * ri))))

    # ------------------------------------------------------------------ [3]
    print(SEP)
    print("[3] massive case: eps_a(m) = sup_k |h(k)-k^2| / (k^2 + a^2 m^2)   (k = a p)")
    print("    a          eps_a(m=1)        eps_a(m=0.1)      1-4/pi^2 - eps   rel.corr")
    for a in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5):
        row = []
        for m in (1.0, 0.1):
            m2 = (a * m) ** 2
            row.append(sup_ratio_grid(h_wilson, 3, 201, m2=m2))
        print("    %-10.0e %.12f   %.12f   %.3e     %.2e"
              % (a, row[0], row[1], c_exact - row[0],
                 (c_exact - row[0]) / c_exact))
    print("    => eps_a(m) = c_exact - O(a^2 m^2):  NOT O(a), NOT O(a^k) for any k>0.")
    print("       (corner formula: d(pi^2-4)/(d pi^2 + a^2 m^2), d=3.)")
    for a in (1e-1, 1e-2, 1e-4):
        print("       a=%-8.0e corner = %.12f" % (a, 3.0 * (np.pi ** 2 - 4.0) / (3.0 * np.pi ** 2 + a * a)))
    print("    => the RATE version (eps_a + b_a -> 0) fails: liminf eps_a = 1-4/pi^2 > 0.")

    # ------------------------------------------------------------------ [4]
    print(SEP)
    print("[4] non-invariance: same continuum limit, different eps^form")
    print("    h_t = h_W + t*phi,  phi(k) = (sum_i (1-cos k_i))^2 = O(k^4) near 0")
    print("      (irrelevant, positive, dim>=6 addition: IR dispersion unchanged)")

    def h_t(k, t):
        return h_wilson(k) + t * np.sum(1.0 - np.cos(k), axis=-1) ** 2

    print("    t        sup |h_t-k^2|/k^2   > 1 ?")
    for t in (0.0, 0.05, 0.2, 0.5, 1.0, 5.0, 50.0):
        s = sup_ratio_grid(lambda k, t=t: h_t(k, t), 3, 201)
        print("    %-8.2f %.12f    %s" % (t, s, "YES" if s > 1.0 else "no"))
    print("    => eps^form is NOT an invariant of the continuum theory; it can be made")
    print("       arbitrarily large by irrelevant (continuum-invisible) additions.")
    print("    Symanzik-improved Laplacian  h_S = h_W - h_W^2/12 (d=3):")
    sB = sup_ratio_grid(lambda k: h_wilson(k) - h_wilson(k) ** 2 / 12.0, 3, 201)
    print("       sup = %.12f   (exactly 1: h_S(pi,pi,pi)=0 => corner mode becomes massless;" % sB)
    print("        hence also NOT an admissible 'same-continuum-limit' regulator)")


# ---------------------------------------------------------------------- [5]
def section5():
    print(SEP)
    print("[5] Shale-Stinespring / quasi-equivalence of lattice and continuum Gaussians")
    ax = np.linspace(-np.pi, np.pi, 401)
    K = np.stack(np.meshgrid(ax, ax, ax, indexing="ij"), axis=-1).reshape(-1, 3)
    k2 = ksq(K)
    r = np.where(k2 > 0, h_wilson(K) / k2, 1.0)
    d3k = (2 * np.pi / 400) ** 3
    m_mask = np.abs(r - 1.0) >= 0.1
    print("    fraction of BZ with |h/k^2 - 1| >= 0.1 : %.6f" % (m_mask.mean(),))
    I2 = float(np.sum((r - 1.0) ** 2) * d3k)
    print("    I2 = (1/(2pi)^3) int_BZ (h/k^2-1)^2 dk = %.6f" % I2)
    print("    #modes = (L/a)^3 * I2   --> infinity as a->0, for every L>0")
    for a in (1e-1, 1e-2, 1e-3):
        print("      a=%-8.0e  L=1:  #modes with |r-1|>=0.1 ~ %.3e"
              % (a, m_mask.mean() * (1.0 / a) ** 3))
    print("    => C_a^{-1/2} C C_a^{-1/2} - 1 is NOT Hilbert-Schmidt: the lattice and")
    print("       continuum Gaussian vacua are MUTUALLY SINGULAR (any a>0).")
    print("       Consequence: the identification iota_a cannot be canonical.")


# ---------------------------------------------------------------------- [6]
def section6():
    print(SEP)
    print("[6] GO2 family vs a physical (lattice) family: who carries the spectrum?")
    print("    RETRACTION (draft v1 of this script claimed N_a(E) -> oo at fixed E: FALSE).")
    print("    At fixed PHYSICAL energy E, any local QFT at fixed volume has an a-independent")
    print("    density of states; free field: N(E) = (L/2pi)^3 (4pi/3) E^3, no a.  The (L/a)^3")
    print("    growth sits only at the cutoff scale E ~ 1/a.  So counting does NOT separate")
    print("    the two families.  What separates them is WHO carries the a-dependence:")
    print("      GO2 (GO2-SOFTMODE.md:66): sigma(H_a) = N_0 u {eps_a}")
    print("         -> exactly ONE moving branch; the rest of the spectrum is a-independent.")
    print("      physical KS: g^2(a)/2 -> 0 and 1/(2 g^2(a)) -> oo (GO2:162) and the low")
    print("         spectrum is a collective property of 9(L/a)^3 link dof (SU(2)).")
    print("    link dof count for SU(2), 3 colours x (N^2-1) x (L/a)^3 = 9 (L/a)^3, L=1:")
    for a in (1e-1, 1e-2, 1e-3, 1e-4):
        print("      a=%-8.0e  links = %.3e      GO2 moving branches = 1"
              % (a, 9.0 * (1.0 / a) ** 3))
    print("    => the conversion GO2-family -> physical family fails at S2 (rank-one defect")
    print("       vs local plaquette sum), NOT at S1 (T-matrix identity, GO2:104 passes).")


if __name__ == "__main__":
    main()
    section5()
    section6()
    print(SEP)
    print("done. deterministic; no MC; no random numbers.")
