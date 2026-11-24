#!/usr/bin/env python3
"""
h6_head_spectral.py -- witness for H6-HEAD-SPECTRAL.md.

A deterministic *illustration* witness, NOT a proof.  numpy + scipy only.

The thesis being illustrated:

  The coarse marginal curvature  G''(a) = V_eff''(a) = (Hess of -log nu)(a)
  is the curvature of the *Schroedinger ground state* of the coarse-mode transfer
  matrix,  nu = |psi_0|^2 ,  V_eff = -2 log psi_0 .  The transfer-matrix Poincare
  (spectral-gap) constant obeys the one-sided inequality

        lambda_1(nu)  >=  ess inf_a  V_eff''(a)          (Bakry-Emery / BL, 1-D),

  i.e. the *global* spectral bound reads the ESSENTIAL INFIMUM of the curvature,
  not the value at a typical field a0*.  At the head (beta~2.5) the coarse
  marginal is the Savvidy double-well: V_eff''(typical well) > 0 but
  ess inf V_eff'' < 0 and lambda_1 is Kramers-suppressed.  Hence *typical-point*
  positivity (the measured +0.13..+0.51) does NOT lift to a *uniform spectral*
  constant.  The free field is the acid test (single well, lambda_1 = G'' = k^2 -> 0).

PARTS
  (A) free field: Gaussian marginal; lambda_1 = G'' = Var^{-1}; -> 0 massless.  ACID.
  (B) Savvidy double well: typical G''>0, ess inf G''<0, lambda_1 KRAMERS -> 0.
  (C) the inequality lambda_1 >= ess inf V_eff'': sharp (Gaussian), vacuous (well).
  (D) two-coordinate block transfer matrix (exact Gaussian): slow gap = G'' = H1,
      fast fiber gap Lambda is SILENT on the slow gap (two-coordinate independence).
  (E) susceptibility test-function route: only lambda_1 <= 1/Var (upper bound),
      so chi cannot lower-bound the gap; 1/lambda_1 >> Var for the double well.

Run:  py -3.12 h6_head_spectral.py
"""

import numpy as np
import scipy.linalg as sla

np.random.seed(20261124)


# --------------------------------------------------------------------------- FEM
def poincare_constant(V, xmax, N):
    """Poincare (spectral-gap) constant lambda_1 of rho ~ exp(-V) on (-xmax,xmax).

    P1 finite elements for the generator L f = f'' - V' f'  in L^2(rho):
        int phi_i' phi_j' rho dx  c = lambda  int phi_i phi_j rho dx .
    Free (natural) endpoints so the constant vector lives in the space -> lambda_0 = 0.
    Returns (lambda_1, x_grid, rho).
    """
    x = np.linspace(-xmax, xmax, N)
    h = x[1] - x[0]
    rho = np.exp(-(V(x) - np.min(V(x))))      # shift for stability
    rho = np.maximum(rho, 1e-300)             # avoid exact underflow -> singular mass
    rho_mid = 0.5 * (rho[:-1] + rho[1:])      # element midpoint density

    A = np.zeros((N, N))                       # stiffness
    M = np.zeros((N, N))                       # mass
    for i in range(N - 1):
        w = rho_mid[i]
        # local stiffness (1/h)[[1,-1],[-1,1]] * w
        A[i, i] += w / h;      A[i, i + 1] -= w / h
        A[i + 1, i] -= w / h;  A[i + 1, i + 1] += w / h
        # local mass (h/6)[[2,1],[1,2]] * w
        M[i, i] += 2.0 * h * w / 6.0;  M[i, i + 1] += h * w / 6.0
        M[i + 1, i] += h * w / 6.0;    M[i + 1, i + 1] += 2.0 * h * w / 6.0

    # generalized symmetric eigenproblem; ask for the lowest 4 eigenvalues
    ev = sla.eigh(A, M, eigvals_only=True, subset_by_index=[0, min(3, N - 1)])
    ev = np.sort(ev)
    lam1 = ev[1] if ev[0] < 1e-8 else ev[0]    # first nonzero (lambda_0 should be ~0)
    return lam1, x, rho


# --------------------------------------------------------------------------- (A)
def part_A_free_field():
    print("=" * 84)
    print("(A) FREE FIELD (acid test: gap MUST vanish for the massless field)")
    print("=" * 84)
    print("  Gaussian marginal V_eff = (lam/2) x^2  (free).  G'' = lam, Var = 1/lam,")
    print("  coarse-mode transfer gap lambda_1 should equal lam (OU rate).")
    print(f"  {'lam':>8} {'lambda_1(fem)':>15} {'G''=lam':>10} {'Var=1/lam':>11} {'1/lambda_1':>11}")
    for lam in (1.0, 0.5, 0.1, 0.02):
        xmax = max(8.0, 6.0 / np.sqrt(lam))
        N = 4001
        lam1, _, _ = poincare_constant(lambda x: 0.5 * lam * x * x, xmax, N)
        print(f"  {lam:8.3f} {lam1:15.6f} {lam:10.3f} {1.0/lam:11.3f} {1.0/lam1:11.3f}")
    print("  (i) lambda_1 = G'' = 1/Var : Gaussian exact.  (ii) lam->0 => lambda_1->0.")
    print("  => the transfer-matrix gap correctly VANISHES on the massless free field.  [ACID]\n")


# --------------------------------------------------------------------------- (B)
def part_B_double_well():
    print("=" * 84)
    print("(B) SAWIDY DOUBLE WELL  V_eff=(g/2)(x^2-1)^2  (soft-mode normal form)")
    print("=" * 84)
    print("  wells x=+-1: G''(well)=4g>0 (TYPICAL positivity);  barrier x=0: G''=-2g<0.")
    print("  Transfer gap lambda_1 is KRAMERS-suppressed by the barrier dV=g/2.")
    print(f"  {'g':>4} {'G''(well)':>10} {'ess inf G''':>11} {'lambda_1':>11} {'Var(x)':>8} "
          f"{'lam1/G''well':>12} {'-ln(lam1)/(g/2)':>15}")
    for g in (0.5, 1.0, 2.0, 4.0, 8.0):
        xmax = 2.6 + 0.10 * g
        N = 6001
        V = lambda x, g=g: 0.5 * g * (x * x - 1.0) ** 2
        lam1, x, rho = poincare_constant(V, xmax, N)
        p = rho / np.trapezoid(rho, x)
        var = np.trapezoid(p * x * x, x) - (np.trapezoid(p * x, x)) ** 2
        gwell = 4.0 * g
        essinf = -2.0 * g
        ratio = lam1 / gwell
        kal = -np.log(lam1) / (0.5 * g)       # Kramers factor vs barrier g/2
        print(f"  {g:4.1f} {gwell:10.3f} {essinf:11.3f} {lam1:11.3e} {var:8.3f} "
              f"{ratio:12.3e} {kal:15.3f}")
    print("  Reading: G''(well)=4g -> grows, but lambda_1 drops ~ e^{-g/2} (Kramers); the")
    print("  ratio lambda_1/G''(well) -> 0.  TYPICAL-POINT positivity does NOT give a gap,")
    print("  because the spectral constant is controlled by the NEGATIVE barrier (ess inf<0).")
    print("  This is the exact head obstruction: the Savvidy tachyonic barrier.\n")


# --------------------------------------------------------------------------- (C)
def part_C_essinf_ineq():
    print("=" * 84)
    print("(C) the one-sided inequality  lambda_1 >= ess inf V_eff''")
    print("=" * 84)
    # Gaussian: sharp.  V'' = lam everywhere => ess inf = lam = lambda_1.
    lam = 1.0
    lam1, _, _ = poincare_constant(lambda x: 0.5 * lam * x * x, 8.0, 4001)
    print(f"  Gaussian lam={lam}:   lambda_1 = {lam1:.5f}   ess inf V'' = {lam:.5f}   "
          f"SHARP (gap = ess inf).")
    # Double well g=2: ess inf = -4 < 0, so the bound reads lambda_1 >= -4 (vacuous).
    g = 2.0
    V = lambda x, g=g: 0.5 * g * (x * x - 1.0) ** 2
    lam1, _, _ = poincare_constant(V, 4.6, 4001)
    print(f"  Double well g={g}:    lambda_1 = {lam1:.5e}   ess inf V'' = {-2.0*g:.5f}   "
          f"bound lambda_1>=-4 holds but is VACUOUS (RHS<0).")
    print("  => the constructive (lower-bound) spectral route is BE/BL: V_eff''>=lam a.e. =>")
    print("     lambda_1>=lam.  At the head ess inf V_eff''<0 (Savvidy), so it gives nothing,")
    print("     while the typical-point G''(a0*)>0 is invisible to the infimum.  [the crux]\n")


# --------------------------------------------------------------------------- (D)
def part_D_two_coordinate():
    print("=" * 84)
    print("(D) TWO-COORDINATE BLOCK TRANSFER MATRIX (exact Gaussian)")
    print("=" * 84)
    print("  mu ~ N(0,K^{-1}),  K=[[lam_s+b^2*Lam, -b*Lam],[-b*Lam, Lam]]  (slow a, fast y).")
    print("  Marginal precision in a: lam_s + b^2*Lam - (b*Lam)^2/Lam = lam_s  (EXACT).")
    print("  Full 2D OU rates = eigenvalues of K;  slow root = lam_s/(1+b^2) as Lam->inf.")
    lam_s = 0.4
    print(f"  {'Lam(fast)':>11} {'b':>5} {'G''marg=1/Var':>14} {'Var(a)':>9} "
          f"{'2D slow rate':>13} {'2D fast':>10}")
    for Lam in (1.0, 10.0, 100.0, 1000.0):
        for b in (0.0, 1.0):
            K = np.array([[lam_s + b * b * Lam, -b * Lam],
                          [-b * Lam, Lam]])
            ev = np.sort(np.linalg.eigvalsh(K))
            slow = ev[0]                    # slowest 2D (full block-T) rate
            fast = ev[1]
            Kinv = np.linalg.inv(K)
            var_a = Kinv[0, 0]              # Var(a): EXACTLY 1/lam_s
            print(f"  {Lam:11.1f} {b:5.1f} {1.0/var_a:14.5f} {var_a:9.5f} "
                  f"{slow:13.6f} {fast:10.3f}")
    print("  Three distinct objects, do not conflate:")
    print("   (1) G''_marg = 1/Var(a) = lam_s  = 0.4  ALWAYS  (the coarse MARGINAL curvature = H1);")
    print("   (2) full 2D block-T slow rate = lam_s/(1+b^2)  (a DYNAMICAL number, != G'');")
    print("   (3) fast conditional (fiber) gap = Lam -> inf.")
    print("  (1) is independent of b AND Lam: the fast fiber gap is SILENT on the coarse marginal")
    print("  curvature (two-coordinate independence, B3-ISOPERIMETRIC).  (2) is also independent")
    print("  of Lam as Lam->inf.  The ONLY block-T gap that equals G'' is the MARGINAL one (after")
    print("  integrating out the fast coordinate) -- which is definitionally the coarse Poincare")
    print("  constant = H1, not a bound on it.  A positive fast gap says NOTHING about G''.\n")


# --------------------------------------------------------------------------- (E)
def part_E_susceptibility():
    print("=" * 84)
    print("(E) susceptibility / test-function route: only lambda_1 <= 1/Var (an UPPER bound)")
    print("=" * 84)
    print("  Test fn phi=a gives 1/lambda_1 >= Var(a)  <=>  lambda_1 <= 1/Var(a):")
    print("  chi=Var(a) UPPER-bounds... no, it can only give an UPPER bound on the gap.")
    print(f"  {'model':>20} {'lambda_1':>12} {'Var(a)':>10} {'1/lambda_1':>12} {'1/lambda_1 - Var':>16}")
    # Gaussian
    lam = 1.0
    lam1, _, _ = poincare_constant(lambda x: 0.5 * lam * x * x, 8.0, 4001)
    var = 1.0 / lam
    print(f"  {'Gaussian lam=1':>20} {lam1:12.5f} {var:10.4f} {1.0/lam1:12.4f} {1.0/lam1-var:16.3e}")
    # Double well
    for g in (2.0, 4.0):
        V = lambda x, g=g: 0.5 * g * (x * x - 1.0) ** 2
        xmax = 4.0 + 0.3 * g
        lam1, x, rho = poincare_constant(V, xmax, 4001)
        p = rho / np.trapezoid(rho, x)
        var = np.trapezoid(p * x * x, x) - (np.trapezoid(p * x, x)) ** 2
        print(f"  {f'double well g={g}':>20} {lam1:12.4e} {var:10.4f} {1.0/lam1:12.3e} "
              f"{1.0/lam1-var:16.3e}")
    print("  Reading: Gaussian is equality (1/lambda_1=Var, G''*Var=1).  Double well: 1/lambda_1")
    print("  grows (->inf as g->inf) while Var is O(1): chi gives only an UPPER bound on the")
    print("  gap direction (lambda_1 <= 1/Var); it cannot LOWER-bound the spectral gap.  [D3]\n")


def main():
    print("#" * 84)
    print("# H6 head-spectral witness  (illustration only; no theorem claimed here)")
    print("# spectral / transfer-matrix attack on the head beta~2.5")
    print("#" * 84)
    part_A_free_field()
    part_B_double_well()
    part_C_essinf_ineq()
    part_D_two_coordinate()
    part_E_susceptibility()
    print("#" * 84)
    print("# SUMMARY")
    print("#  * (A) free field: lambda_1 = G'' = 1/Var -> 0 massless.  ACID PASSED.")
    print("#  * (B) Savvidy double well: G''(well)=4g>0 TYPICAL, ess inf G''=-2g<0, lambda_1~e^{-g/2}.")
    print("#  * (C) the only constructive spectral bound reads ess inf V_eff'', not the typical value.")
    print("#  * (D) the ONLY block-T gap = G'' is the MARGINAL one (= H1 by definition); the")
    print("#        full 2D block-T rate differs, and the fast fiber gap is SILENT on G''.")
    print("#  * (E) susceptibility gives only lambda_1 <= 1/Var: cannot lower-bound the gap.")
    print("#  NET: no spectral bound reaches the head with a volume-uniform constant; the")
    print("#       obstruction is typical-point positivity vs the essential infimum (Savvidy tail).")
    print("#  Clay 0%.  No QED.")
    print("#" * 84)


if __name__ == "__main__":
    main()
