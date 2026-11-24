#!/usr/bin/env python3
"""
s6_dobrushin.py  --  witness for research/v3/S6-DOBRUSHIN.md

Deterministic toy checks for the Dobrushin / Brascamp-Lieb attack on (star):

  (star)  Var_{mu(.|U')}( <P^T grad S_j, v> )  <=  E< P^T H_j P v, v > + W_j(v) - d_W
  equivalently  G''_j(v) = E[H_j] - Var_j + W_j >= d_W > 0 .

Parts
  A. Dobrushin coefficient for the SU(2) 4D block fibre -> threshold in beta,
     compared with the SZZ 1/48 and the physical 2.5.
  B. Free Gaussian (flat) block, exact linear algebra:  for a Gaussian with
     H = Sigma^{-1} one has Var(<H x, v>) = <H v, v> EXACTLY (equipartition).
     Hence the flat defect D := Var - E[H] = 0 and W = 0  =>  G'' = 0  =>
     (star) FAILS.  The boundary shifts the MEAN of x, never the VARIANCE.
     This is the required free-field acid test.
  C. One-link SU(2) PURE HAAR in the flat chart (3D quadrature): the Haar
     budget W_geo = E[Phi''_VV] - Var(Phi'_V) > 0.  Since B shows the flat
     part contributes exactly 0, W_geo is the ONLY source of the strict
     margin d_W -- it must be positive, and it is.

Run:  py -3 s6_dobrushin.py
"""
import numpy as np

np.set_printoptions(precision=6, suppress=True)


# --------------------------------------------------------------------------
# A.  Dobrushin coefficient for the SU(2) fibre on a 4D block
# --------------------------------------------------------------------------
def part_A():
    d = 4
    n_plaq_through_link = 2 * (d - 1)          # = 6
    n_links_sharing      = 6 * (d - 1)         # <= 18 (upper bound, overlaps ignored)
    c_bond               = 4.0                 # single-bond coefficient, O(beta)

    def alpha(beta):
        return n_links_sharing * c_bond * beta

    beta_dobr = 1.0 / (n_links_sharing * c_bond)

    print("=" * 70)
    print("A.  Dobrushin coefficient, SU(2) 4D block fibre")
    print("=" * 70)
    print(f"  plaquettes through one link           : {n_plaq_through_link}")
    print(f"  links sharing a plaquette with l      : {n_links_sharing}")
    print(f"  single-bond coefficient c_ll'         : {c_bond:.3g} * beta")
    print(f"  alpha(beta) = {n_links_sharing*c_bond:.0f} * beta")
    print(f"  beta_Dobr  (alpha = 1)                : {beta_dobr:.5f}  = 1/{1/beta_dobr:.1f}")
    print(f"  SZZ sharp threshold (paper beta<1/48) : Hooft lambda>96  <=>  beta_W ~ 0.19")
    print(f"  physical head (standard Wilson beta)  : 2.5")
    print(f"  honest distance head/threshold        : ~13x in lambda (beta_W 0.19 vs 2.5)")
    print(f"  crude Dobrushin is a cruder member of the SAME strong-coupling family")
    for b in (beta_dobr, 0.02, 0.19, 2.5):
        tag = "OK  (alpha<1)" if alpha(b) < 1 else "FAILS (alpha>=1)"
        print(f"    alpha({b:9.5f}) = {alpha(b):8.3f}   {tag}")
    print("  -> alpha(2.5) = 180 >> 1: the whole class of local-influence")
    print("     (Dobrushin / disagreement-percolation) methods dies at beta~2.5.")
    print()
    return beta_dobr


# --------------------------------------------------------------------------
# B.  free Gaussian flat block: exact equipartition, defect = 0
# --------------------------------------------------------------------------
def dirichlet_laplacian(n):
    A = np.zeros((n, n))
    for i in range(n):
        A[i, i] = 2.0
        if i > 0:     A[i, i - 1] = -1.0
        if i < n - 1: A[i, i + 1] = -1.0
    return A


def part_B():
    n = 6
    H = dirichlet_laplacian(n)                     # Wilson / free Hessian
    Sigma = np.linalg.inv(H)                        # Gaussian covariance, H = Sigma^{-1}
    v = np.ones(n) / np.sqrt(n)
    E_H   = v @ H @ v                               # E< H v, v >
    Var_G = v @ H @ Sigma @ H @ v                   # Var(<H x, v>)
    print("=" * 70)
    print("B.  Free Gaussian flat block  (H = Sigma^{-1})")
    print("=" * 70)
    print(f"  E[H_vv] = <H v, v>                    : {E_H:.12f}")
    print(f"  Var(G_v) = <H v, Sigma H v>           : {Var_G:.12f}")
    print(f"  defect D = Var - E[H]                 : {Var_G-E_H:.2e}")
    print(f"  W (=0, no Haar)                       : {0.0:.2e}")
    print(f"  G'' = E[H] - Var + W                  : {E_H-Var_G:.2e}")
    # boundary (mean) dependence: Var is mean-blind, so boundary cannot help.
    b = 3.0
    e1 = np.zeros(n); e1[0] = 1.0
    mean_x = b * np.linalg.solve(H, e1)
    print(f"  boundary value b = {b}")
    print(f"  interior mean (first, last)           : {mean_x[0]:.4f}, {mean_x[-1]:.4f}")
    print(f"  Var(G_v) depends on b ?               : no  (variance is mean-blind)")
    print(f"  => (star): G'' = {E_H-Var_G:.2e} >= d_W > 0 ?  NO  (fails, G''=0)")
    print("  The flat part of (star) is an EQUALITY with zero margin: the strict")
    print("  sign can only come from the Haar budget W (part C).")
    print()
    return E_H, Var_G


# --------------------------------------------------------------------------
# C.  one-link SU(2) pure Haar, flat chart, 3D quadrature
#     Phi(theta) = -2 log( sin(theta/2) / (theta/2) ),  theta = |X| in (0,2pi)
#     flat chart, X in su(2) ~ R^3,  dU = e^{-Phi(X)} dX ,  e^{-Phi} = (sin(theta/2)/(theta/2))^2.
#     W_geo is the PURE-HAAR MARGINAL curvature at a=0 (campaign: beta=0 => G''=W_geo):
#       f(a) = marginal density of a=<v,X> ;  V_eff(a) = -log f(a) ;  W_geo=V_eff''(0).
#     By cylindrical symmetry  f(a) = 2*pi * int_{|a|}^{2pi} theta e^{-Phi(theta)} dtheta,
#     hence  V_eff''(0) = -f''(0)/f(0) = 1 / int_0^{2pi} theta e^{-Phi(theta)} dtheta.
# --------------------------------------------------------------------------
def part_C(n_grid=2000001):
    theta = np.linspace(1e-7, 2*np.pi - 1e-7, n_grid)
    # e^{-Phi}(theta) = (sin(theta/2)/(theta/2))^2 = 4 sin^2(theta/2)/theta^2
    emPhi = (np.sin(theta / 2.0) / (theta / 2.0))**2
    norm  = np.trapezoid(theta * emPhi, theta)         # = f(0)/(2*pi)
    W_geo = 1.0 / norm                                 # = V_eff''(0)
    print("=" * 70)
    print("C.  one-link SU(2) pure Haar: the budget W_geo = V_eff''(0) > 0")
    print("=" * 70)
    print(f"  int_0^2pi theta e^-Phi dtheta         : {norm:.6f}")
    print(f"  W_geo = V_eff''(0) = 1/norm           : {W_geo:.6f}")
    print(f"  (campaign measured V_eff''=0.248 at beta=0, E[Phi'']|0 >= 1/6)")
    print(f"  sign of W_geo                         : {'POSITIVE' if W_geo>0 else 'NON-POSITIVE'}")
    print("  => the Haar term is a genuine positive budget; with the flat part")
    print("     identically 0 (part B), it is the only source of the margin d_W.")
    print("     On the COMPACT group this curvature can never go negative (Prekopa: e^-Phi")
    print("     log-concave => marginal log-concave => V_eff''>=0); the free field has NO")
    print("     such term, which is exactly why (star) fails there.")
    print()
    return W_geo


if __name__ == "__main__":
    part_A()
    part_B()
    part_C()
    print("done.  (toy; deterministic; a witness, not a proof.)")
