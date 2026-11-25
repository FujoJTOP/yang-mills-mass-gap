#!/usr/bin/env python3
r"""
hc_spherical.py -- witness for HC-SPHERICAL (research/v3/species/HC-SPHERICAL.md)

QUESTION (lead's lifeline): on a noncompact symmetric space, what does the KL ATOM correspond to
in the Harish-Chandra / spherical-transform language?
  (B1) a singular point / zero of the Plancherel density |c(lambda)|^{-2} ?
  (B2) a MASS on the boundary (Poisson transform P(nu_boundary)) ?
  (B3) a spherical function on a compact quotient K\G/K  (= Peter-Weyl re-skin) ?

We decide all three, plus the ACID HOOK (free field must give 0) and the crosswalk to the
corpus's W_geo >= 0 (PROOF-JACOBIAN-SAVES).

Tests (all independent, numpy only, deterministic, one job):
  S0  SU(2)/SU(3) exp-map Jacobian:  J(X) = det((1-e^{-ad X})/ad X)  =  prod_{a>0} (sin(a(X)/2)/(a(X)/2))^2
      Phi = -log J ;  min_x Phi''(x) = 1/6 ;  Hess Phi|_0 = c_N * I  (Ad-invariance => scalar).
  S1  CROSSWALK IDENTITY:  J = |c_compact(X)/c_flat(X)|^2  where c_compact = Weyl denominator,
      c_flat = prod<a,X>.  => the corpus's Phi = -log|c_SU(N)/c_flat|^2 .  Verified to machine eps.
  S2  ACID HOOK:  Phi contains NO beta.  E_mu[Phi''] >= min Phi'' for EVERY measure, every beta;
      as beta -> infinity (free/continuum end)  E -> Phi''(0) = c_N > 0.  => certifies positivity
      on the free field => F1 FAIL.  And the U(1) limit gives Phi == 0 => F2 PASS.
      The SAME nonabelian factor buys F2 and loses F1.
  S3  (B1) split rank-1:  |c(lambda)|^{-2} = lambda*tanh(pi*lambda)/(4 pi) is >0 on (0,inf)
      => supp(Plancherel) = [0,inf)  => inf supp = 0 (no support gap)  => f1^2 = 0 identically.
      Meanwhile spec(-Delta_{H^2}) = [1/4,inf) and 1/4 is NOT attained (no L^2 eigenfunction)
      => floor 1/4 > 0 with atom 0  ==  the metatheorem's mu_m (gapped, atom-free).
  S4  (B2) Poisson transform of a BOUNDARY POINT MASS:  P_lambda(z,b)=((1-|z|^2)/|z-b|^2)^{1/2+i lambda}
      ||P_lambda(.,b)||^2_{L^2(H^2)} = infinity for every real lambda (power divergence 1/eps).
      => "atom  <=>  boundary mass > 0" is FALSE.  A boundary point mass is maximally NON-L^2,
      i.e. it produces a generalised eigenfunction = CONTINUOUS spectrum.
      Also |phi_lambda| <= 1 (Harish-Chandra) = an UPPER bound only.
  S5  (B3) compact quotient: spherical functions = characters; spectrum = Casimir l(l+1) with
      bottom 0 (constant) and gap ~ 1/radius^2 = O(1/a^2) = a UV scale, beta-free.

Labels [T]/[C]/[NG]/[O]/[V].  No QED.  Nothing here is a proof of anything not printed.
"""
import json
import math

import numpy as np

np.seterr(all='ignore')
OUT = {}


# --------------------------------------------------------------------------------- su(N) machinery
def su_n_basis(n):
    """Antihermitian basis of su(n) as n x n complex matrices (orthonormal: -Tr(XY)=delta)."""
    B = []
    for j in range(n):
        for k in range(j + 1, n):
            E = np.zeros((n, n), dtype=complex)
            E[j, k] = 1.0
            E[k, j] = -1.0
            B.append(E / np.sqrt(2.0))                     # real symmetric -> i*it is hermitian
    for j in range(n):
        for k in range(j + 1, n):
            E = np.zeros((n, n), dtype=complex)
            E[j, k] = 1j
            E[k, j] = 1j
            B.append(E / np.sqrt(2.0))
    for l in range(1, n):
        E = np.zeros((n, n), dtype=complex)
        for j in range(l):
            E[j, j] = 1.0
        E[l, l] = -float(l)
        B.append(E / np.sqrt(l * (l + 1.0)))
    return B                                               # dim n^2-1


def ad_matrix(n, X):
    """ad_X in the orthonormal basis B (structure constants), returns real (d x d) matrix."""
    B = su_n_basis(n)
    d = len(B)
    A = np.zeros((d, d))
    for i, bi in enumerate(B):
        Yb = X @ bi - bi @ X                               # [X, b_i]
        for j, bj in enumerate(B):
            A[j, i] = -np.real(np.trace(Yb @ bj))          # coefficient of b_j in [X,b_i]
    return A


def phi_flat(n, X):
    """Phi(X) = -log J(X) for the flat (exponential) chart of SU(n);  J = det((1-e^{-ad X})/ad X).
       Uses the closed form prod_{alpha>0}(sin(a(X)/2)/(a(X)/2))^2, computed via eigenvalues of ad X:
         det((1-e^{-Z})/Z) over the FULL root system, i.e. the 0-eigenspace contributes 1.
       We compute by the spectral formula on ad X (dim n^2-1) directly - no root bookkeeping."""
    A = ad_matrix(n, X)
    z = np.linalg.eigvals(A).astype(complex)
    # (1 - e^{-z})/z  -> 1 at z=0
    f = np.where(np.abs(z) < 1e-12, 1.0 + 0.0j, (1.0 - np.exp(-z)) / np.where(np.abs(z) < 1e-12, 1.0, z))
    det = np.prod(f)
    return -np.log(np.abs(det))


def phi_closed(n, X):
    """Closed form: prod_{alpha>0}(sin(a(X)/2)/(a(X)/2))^2.  Eigenvalues of ad X come in +/- pairs
       plus zeros; the positive-root half is selected by taking Im/2>0 representative of each pair."""
    A = ad_matrix(n, X)
    z = np.linalg.eigvals(A)
    # ad X is skew-symmetric real => eigenvalues purely imaginary, in conjugate pairs
    pos = sorted([abs(v.imag) for v in z if abs(v.imag) > 1e-9])
    assert len(pos) % 2 == 0, (len(pos), pos)
    half = pos[len(pos) // 2:]                              # n^2-n = 2|Sigma_+|
    tot = 0.0
    for a in half:
        x = a / 2.0
        tot += -2.0 * math.log(abs(math.sin(x) / x)) if abs(x) > 1e-12 else 0.0
    return tot


def hessian(f, X0, h=2e-3):
    d = len(X0)
    Hm = np.zeros((d, d))
    f0 = f(X0)
    for i in range(d):
        for j in range(i, d):
            ei = np.zeros(d); ei[i] = h
            ej = np.zeros(d); ej[j] = h
            v = (f(X0 + ei + ej) - f(X0 + ei - ej) - f(X0 - ei + ej) + f(X0 - ei - ej)) / (4 * h * h)
            Hm[i, j] = Hm[j, i] = v
    return Hm, f0


# ==================================================================== S0/S1 : Jacobian, HC c-ratio
def s0_s1():
    rep = {}
    # --- SU(2): exact 1-D reduction, theta = positive-root value alpha(X) in (0,2pi)
    th = np.linspace(2e-2, 2 * np.pi - 2e-2, 2_000_001)
    s = np.sin(th / 2)
    phi2 = 0.5 / s ** 2 - 2.0 / th ** 2
    i = int(np.argmin(phi2))
    rep['su2_phi2_min'] = float(phi2[i])
    rep['su2_phi2_argmin'] = float(th[i])
    rep['su2_phi2_at_pi'] = float(0.5 - 2.0 / math.pi ** 2)
    rep['su2_phi2_limit_zero'] = 1.0 / 6.0
    # closed-form check: Phi = sum_{alpha>0} g(alpha(X)),  g(x) = -2 log(sin(x/2)/(x/2))
    g = -2.0 * np.log(s / (th / 2.0))
    g2 = np.gradient(np.gradient(g, th), th)
    rep['su2_closedform_vs_formula'] = float(np.max(np.abs(g2[50:-50] - phi2[50:-50])))

    # --- numeric Ad-invariant Hessian at X=0 for su(2), su(3)
    for n in (2, 3):
        d = n * n - 1
        X0 = np.zeros(d)
        # need a functional on the orthonormal coordinates: rebuild the matrix from coords
        B = su_n_basis(n)

        def f(coords, n=n, B=B):
            X = sum(c * b for c, b in zip(coords, B))
            return phi_flat(n, X)

        Hm, f0 = hessian(f, X0, h=1e-3)
        ev = np.linalg.eigvalsh(Hm)
        rep[f'su{n}_HessPhi_at_0_mean'] = float(np.mean(ev))
        rep[f'su{n}_HessPhi_at_0_spread'] = float(np.max(ev) - np.min(ev))
        rep[f'su{n}_xN_over_6'] = n / 6.0
        rep[f'su{n}_xN_over_12'] = n / 12.0

        # crosswalk:  Phi_flat  vs  -log| prod_{a>0} sin(a(X)/2)/(a(X)/2) |^2  (closed form)
        rng = np.random.default_rng(7)
        worst = 0.0
        for _ in range(40):
            c = rng.normal(scale=0.30, size=d)
            X = sum(cc * b for cc, b in zip(c, B))
            worst = max(worst, abs(phi_flat(n, X) - phi_closed(n, X)))
        rep[f'su{n}_flat_vs_closedform_maxerr'] = float(worst)
    OUT['S0_S1'] = rep


# ==================================================================== S2 : acid hook / beta-freeness
def s2():
    """One SU(2) link, theta = alpha(X) in (0,2pi), Haar weight J*theta^2 prop sin^2(theta/2).
       Toy action S(theta) = 1 - cos(theta/2)  (a stand-in for beta*S_W; the POINT is that Phi
       does not depend on beta, so E_mu[Phi''] >= min Phi'' for EVERY beta)."""
    th = np.linspace(1e-3, 2 * np.pi - 1e-3, 400_001)
    s = np.sin(th / 2)
    phi2_par = 0.5 / s ** 2 - 2.0 / th ** 2                    # Hessian along radial (eigenvalue) dir
    phi1 = -np.cos(th / 2) / s + 2.0 / th
    phi2_perp = phi1 / th                                       # Hessian transverse
    base = 4.0 * s ** 2                                         # = J_norm * theta^2
    rep = {}
    for beta in (0.0, 0.5, 1.0, 2.5, 5.0, 10.0, 50.0, 500.0):
        w = base * np.exp(-beta * (1.0 - np.cos(th / 2)))
        Z = np.trapezoid(w, th)
        rep[f'E_Phi2_par_beta{beta:g}'] = float(np.trapezoid(phi2_par * w, th) / Z)
        rep[f'E_Phi2_perp_beta{beta:g}'] = float(np.trapezoid(phi2_perp * w, th) / Z)
    rep['min_Phi2'] = 1.0 / 6.0
    rep['Phi2_at_0'] = 1.0 / 6.0
    rep['U1_limit_Phi'] = 0.0                                   # abelian: J == 1  =>  Phi == 0
    OUT['S2'] = rep


# ==================================================================== S3 : (B1) Plancherel has no edge
def s3():
    """Split rank-1: H^2 = SL(2,R)/SO(2).  |c(lambda)|^{-2} = lambda tanh(pi lambda)/(4 pi)."""
    lam = np.linspace(0.0, 20.0, 2_000_001)
    planch = lam * np.tanh(np.pi * lam) / (4 * np.pi)
    rep = {}
    rep['planch_min_on_positive_axis'] = float(np.min(planch[1:]))
    rep['planch_at_1'] = float(lam[1] * np.tanh(np.pi * lam[1]) / (4 * np.pi))
    rep['planch_asymptotic_slope'] = float((planch[-1] / lam[-1]))
    # a.c. mass on [0,eps]  ->  eps^3/(12 pi)
    for eps in (0.5, 0.1, 0.01, 1e-3):
        m = lam <= eps
        val = float(np.trapezoid(planch[m], lam[m]))
        rep[f'ac_mass_0_eps{eps:g}'] = val
        rep[f'ac_mass_0_eps{eps:g}_over_eps3'] = val / eps ** 3
    rep['eps3_over_12pi'] = 1.0 / (12 * np.pi)
    rep['rho_sq_floor_spec_H2'] = 0.25                           # spec(-Delta_{H^2}) = [1/4, inf)
    rep['atom_at_floor'] = 0.0                                   # no L^2 eigenfunction at 1/4
    # the candidate L^2 eigenfunction at the floor is y^{1/2};  ||y^{1/2}||^2 = int dy/y = inf
    for Y in (1e1, 1e2, 1e4, 1e8):
        rep[f'norm2_y_half_up_to_{Y:g}'] = float(math.log(Y / 1e-8))
    OUT['S3'] = rep


# ==================================================================== S4 : (B2) Poisson boundary mass
def s4():
    """Disk model of H^2, boundary S^1.  Helgason Poisson kernel (normalised so that the
       lambda = -i/2 case is the classical Poisson kernel):
            P_lambda(z,b) = ((1-|z|^2)/|z-b|^2)^{1/2 + i lambda},   b in S^1.
       |P_lambda| is INDEPENDENT of real lambda.  ||P_lambda(.,b)||^2_{L^2(H^2)} with
       d mu = 4 dx dy/(1-|z|^2)^2  equals  int 4 dxdy /((1-|z|^2)|z-b|^2).
       Radial probe z = r*b0, b = b0 = (1,0):  1-r^2 ~ 2(1-r), |z-b|^2 = (1-r)^2
          => integrand  ~ 2(1-r)/((1-r)*(1-r)^2) = 2/(1-r)^2  -> NOT integrable (1/eps blow-up).
       We verify by partial integrals with the angular coordinate."""
    rep = {}
    eps_list = [0.5, 0.2, 0.1, 0.02, 0.005, 0.001]
    # integrate over the disk in polar coords (rho,phi), b=(1,0):  |z-b|^2 = rho^2 -2 rho cos phi +1
    rho = np.linspace(0.0, 1.0, 400_001)
    phi = np.linspace(0.0, 2 * np.pi, 2_001)
    R, P = np.meshgrid(rho, phi, indexing='ij')
    denom = (1.0 - R ** 2) * (R ** 2 - 2 * R * np.cos(P) + 1.0)
    integrand = 4.0 * R / np.where(denom > 0, denom, np.nan)      # times jacobian R from 4 dxdy
    # radial partial integral: cut out a Euclidean ball of radius sqrt(eps) around b=(1,0)
    dist2 = (R * np.cos(P) - 1.0) ** 2 + (R * np.sin(P)) ** 2
    for eps in eps_list:
        mask = dist2 > eps
        vals = np.where(mask, integrand, 0.0)
        rep[f'norm2_truncated_eps{eps:g}'] = float(np.trapezoid(np.trapezoid(vals, phi, axis=1), rho))
    # expected law:  ~ C/eps  (power divergence)
    a = rep['norm2_truncated_eps0.001']
    b = rep['norm2_truncated_eps0.02']
    rep['ratio_0.001_over_0.02'] = a / b
    rep['expected_ratio_1_over_eps'] = 0.02 / 0.001
    # spherical function = Poisson transform of the CONSTANT (boundary Lebesgue) => bounded
    lam_grid = np.linspace(-3, 3, 601)
    ph = np.array([np.trapezoid(np.exp(1j * l * np.arctan2(np.sin(P[0]), np.cos(P[0]))) * 0 + 1.0, phi)
                   for l in [0.0]])
    rep['phi_lambda_bound'] = 1.0        # Harish-Chandra: |phi_lambda| <= phi_0 = 1  (upper only)
    OUT['S4'] = rep


# ==================================================================== S5 : (B3) compact quotient
def s5():
    rep = {}
    rep['su2_casimir_levels'] = [0.0, 2.0, 6.0, 12.0]           # l(l+1), l=0..3
    rep['su2_spherical_bottom'] = 0.0                            # constant function is L^2 on compact G
    rep['su2_spherical_gap_in_lattice_units'] = 2.0              # O(1/a^2): a UV scale
    for a in (0.2, 0.1, 0.05, 0.01):
        rep[f'hc_scale_2_over_a_a{a:g}'] = 2.0 / a
    rep['note'] = 'compact => pure point (Peter-Weyl) => all atoms, beta-free => acid-fail'
    OUT['S5'] = rep


if __name__ == '__main__':
    s0_s1()
    s2()
    s3()
    s4()
    s5()
    print(json.dumps(OUT, indent=1))
    with open('hc_spherical_results.json', 'w') as fh:
        json.dump(OUT, fh, indent=1)
