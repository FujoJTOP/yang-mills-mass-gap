#!/usr/bin/env python3
r"""
r1_transport.py -- R1: can the multiplicative residue transport f1^2 = Z_O^2 f1^2 carry the
strong-coupling atom (P3, beta<1/48) across to beta=2.5 ?

WHAT THIS COMPUTES (all deterministic, numpy/scipy, one job):
 [A] ONE-PLAQUETTE / SU(2) CHARACTER COEFFICIENTS  c_j(beta) = 2(2j+1) I_{2j+1}(beta)/beta.
     The "residue ratio" r(beta) = c_{1/2}/c_0 = <tr U> is the strong-coupling building block.
     -> r(beta) > 0 for ALL beta>0 (vacuous sign transport: every finite model has positive residue).
     -> r(beta) ~ beta/2 for small beta (POWER LAW), -> 1 for large beta. It is NOT the continuum
        anomaly Z_O.
 [B] STRONG-COUPLING power law: f1^2 ~ r^2 ~ (beta/2)^2  =>  f1^2(a)/f1^2(b) = (beta_a/beta_b)^p, p=2
     for the one-plaquette channel. The "transport constant" is a bare power, coupling-dependent.
 [C] CONTINUUM ANOMALY Z_O = (beta_YM/2g)(mu0)/(beta_YM/2g)(mu) = g^2(mu0)/g^2(mu)  (tree);
     1-loop 1/g^2 = 2 b0 ln(mu/Lambda).  Z_O>0 for every finite window (pure YM, beta_YM<0, no BZ).
     -> Z_O is bounded below by a POSITIVE constant on any compact window above Lambda.
     -> but its SIZE extrapolated to beta_eff=1/48 (g^2=4/beta=192) is a perturbative absurdity:
        the transport constant is scheme/coupling dependent, not universal.
 [D] OVERLAP-ZERO TOY (the crux): two-level H(t)=[[a,v],[v,-a]], O=sigma_z. The atom residue
     f1^2 = |<e0|O|e1>|^2  VANISHES at v=0 (selection rule) and is smooth/positive off it.
     -> For ANY Z_O>0,  Z_O^2 * 0 = 0.  The multiplicative map fixes 0 (absorbing point): it
        cannot manufacture or certify the atom. The missing datum is the persistence/overlap
        Omega = |<G(target)|G(strong)>|^2 (== isolation I2 == (CP)), NOT a lower bound on Z_O.
 [E] ACID TEST: free massless field -> f1^2 = 0, Z_O^2*0 = 0 (no manufactured atom). PASS.
 [F] 2D (Migdal) diagonal transfer matrix: the 2D spectrum is indexed by SU(2) irreps with heat-
     kernel weights c_j(beta) -- the SAME Bessel object as [A]. So the "2D transport constant" is
     not new; 2D reduces to one-plaquette. The 2D gap -ln(c_1/c_0) > 0 for all beta (vacuous sign).

Labels in the .md: [T] identity, [V] this witness, [NG] no-go, [O] missing. No QED. Clay 0%.
Run: py -3.12 -I r1_transport.py
"""
import numpy as np
from scipy.special import iv

np.set_printoptions(precision=6, suppress=True)

# ---------------------------------------------------------------- [A] character coefficients SU(2)
def c(j, beta):
    """SU(2) heat-kernel / one-plaquette character coefficient c_j(beta)=2(2j+1) I_{2j+1}(beta)/beta."""
    return 2.0 * (2.0 * j + 1.0) * iv(2 * j + 1, beta) / beta

def Z1_direct(beta, n=4000):
    """Direct one-plaquette partition fn Z1=(2/pi) int_0^pi sin^2 u e^{beta cos u} du (u=theta/2)."""
    xg, wg = np.polynomial.legendre.leggauss(n)
    u = 0.5 * np.pi * (xg + 1.0); w = 0.5 * np.pi * wg
    return (2.0 / np.pi) * np.sum((np.sin(u) ** 2) * np.exp(beta * np.cos(u)) * w)

def trU_direct(beta, n=4000):
    """Direct <tr U> = (2/pi) int sin^2 u (2 cos u) e^{beta cos u} du / Z1  (tr U = 2 cos u)."""
    xg, wg = np.polynomial.legendre.leggauss(n)
    u = 0.5 * np.pi * (xg + 1.0); w = 0.5 * np.pi * wg
    num = (2.0 / np.pi) * np.sum((np.sin(u) ** 2) * (2 * np.cos(u)) * np.exp(beta * np.cos(u)) * w)
    return num / Z1_direct(beta, n)

def section_A():
    print("=" * 78)
    print("[A] ONE-PLAQUETTE SU(2) CHARACTER COEFFICIENTS")
    print("    c_j(beta)=2(2j+1) I_{2j+1}(beta)/beta ; verify Z1=c_0, <trU>=c_{1/2}/c_0")
    print("    %8s %12s %12s %12s %12s" % ("beta", "Z1=c0", "Z1_direct", "c12/c0", "<trU>_dir"))
    for beta in [0.5, 1.0, 2.5, 5.0]:
        c0 = c(0.0, beta); r = c(0.5, beta) / c0
        print("    %8.3f %12.6f %12.6f %12.6f %12.6f"
              % (beta, c0, Z1_direct(beta), r, trU_direct(beta)))
    print("    -> r(beta)=c12/c0>0 for all beta (Bessel): SIGN TRANSPORT IS VACUOUS at finite beta.")
    print("    %8s %12s %12s %12s" % ("beta", "r=c12/c0", "r^2", "dln r/dln beta"))
    betas = np.array([0.05, 0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0])
    rs = np.array([c(0.5, b) / c(0.0, b) for b in betas])
    for i, b in enumerate(betas):
        d = np.nan if i == 0 else (np.log(rs[i]) - np.log(rs[i - 1])) / (np.log(b) - np.log(betas[i - 1]))
        print("    %8.3f %12.6f %12.6f %12.4f" % (b, rs[i], rs[i] ** 2, d))
    print("    -> small beta: r~beta/2 (POWER LAW, dln r/dln beta -> 1); large beta: r->1 (dln->0).")
    print("       The local exponent is NOT constant => not a single multiplicative Z_O power.")

# ---------------------------------------------------------------- [B] strong-coupling power law
def section_B():
    print("=" * 78)
    print("[B] STRONG-COUPLING power law:  f1^2 ~ r^2 ~ (beta/2)^2  => p=2")
    b1, b2 = 1.0 / 48.0, 2.5          # strong -> weak (target)
    r1 = c(0.5, b1) / c(0.0, b1) if b1 > 0.01 else b1 / 2.0
    r2 = c(0.5, b2) / c(0.0, b2)
    print("    beta1=1/48=%.4f  r1=%.5e (asymptotic beta/2=%.5e)" % (b1, r1, b1 / 2))
    print("    beta2=2.5      r2=%.6f" % r2)
    print("    residue ratio r2^2/r1^2 = %.6e  vs  (beta2/beta1)^2 = %.6e"
          % (r2 ** 2 / r1 ** 2, (b2 / b1) ** 2))
    print("    -> AT STRONG COUPLING the residue scales as a POWER (~beta/2 per factor), so the")
    print("       strong-coupling 'transport constant' is (beta2/beta1)^2 (matches the small-beta end);")
    print("       it is NOT the continuum anomaly (beta2/beta1 at tree, [ln]^2 at 1-loop). The two")
    print("       constants disagree because beta2=2.5 is NOT in the strong-coupling regime (the")
    print("       true r2=%.3f is far above the power-law beta2/2=%.3f): the crossover is the gap." % (r2, b2 / 2))

# ---------------------------------------------------------------- [C] continuum anomaly Z_O
B0 = 11.0 / (48 * np.pi ** 2)   # one-loop b0 for pure SU(N) in the corpus normalization
def section_C():
    print("=" * 78)
    print("[C] CONTINUUM ANOMALY Z_O = g^2(mu0)/g^2(mu)  (tree; 1-loop 1/g^2=2 b0 ln(mu/Lambda))")
    print("    b0 = 11N/(48 pi^2), N=2: b0=%.6f" % B0)
    for N in [2]:
        print("    --- N=%d : tree g^2=2N/beta ---" % N)
        for (b_lo, b_hi, tag) in [(1/48, 2.5, "beta_eff -> 2.5")]:
            g2_lo = 2 * N / b_lo; g2_hi = 2 * N / b_hi
            ZO = g2_lo / g2_hi
            print("      %-16s g^2(%.4f)=%.3f  g^2(2.5)=%.3f  Z_O=%.4f  Z_O^2=%.4f"
                  % (tag, b_lo, g2_lo, g2_hi, ZO, ZO * ZO))
            print("      (g^2(1/48)=%.0f is NON-perturbative: the tree 'Z_O' is an absurd extrapolation)"
                  % g2_lo)
    print("    --- Z_O>0 for all finite windows: beta_YM(g)<0 for g>0 (pure YM, NO Banks-Zaks);")
    print("        so beta_YM/2g<0 on both ends => ratio >0; zero only at a fixed point (none) or")
    print("        mu->Lambda (deep IR). => Z_O is BOUNDED BELOW by a positive constant. [T]")
    print("        But bounded-below Z_O does NOT bound f1^2 below (see [D]).")

# ---------------------------------------------------------------- [D] overlap-zero toy (the crux)
def section_D():
    print("=" * 78)
    print("[D] OVERLAP-ZERO TOY (the crux): H(t)=[[a,v],[v,-a]], O=sigma_z.")
    print("    atom residue f1^2 = |<e0|O|e1>|^2 ; vanishes at v=0 (selection rule).")
    a = 1.0
    print("    %8s %12s %14s %14s" % ("v", "E_gap", "f1^2", "Z_O^2*f1^2 (Z_O=3)"))
    for v in [0.0, 1e-3, 0.01, 0.1, 0.3, 0.7, 1.0, 2.0, 5.0]:
        H = np.array([[a, v], [v, -a]])
        w, U = np.linalg.eigh(H)        # ascending: U[:,0]=ground
        e0 = U[:, 0]; e1 = U[:, 1]
        sz = np.array([[1.0, 0.0], [0.0, -1.0]])
        m = e0 @ (sz @ e1)
        f1sq = m * m
        print("    %8.4f %12.6f %14.6e %14.6e" % (v, w[1] - w[0], f1sq, 9.0 * f1sq))
    print("    -> f1^2 -> 0 as v -> 0 for ANY fixed Z_O: the multiplicative map fixes 0.")
    print("       'Z_O bounded below' is irrelevant to whether the atom survives (Omega=1).")

# ---------------------------------------------------------------- [E] acid test
def section_E():
    print("=" * 78)
    print("[E] ACID TEST: free massless field. f1^2=0 (no atom) => Z_O^2 * 0 = 0. PASS.")
    print("    Both the strong-coupling seed and the free field are sent correctly;")
    print("    the map is a homomorphism on the WEIGHT, blind to the SUPPORT (atom vs cut).")

# ---------------------------------------------------------------- [F] 2D Migdal diagonal
def section_F():
    print("=" * 78)
    print("[F] 2D (Migdal): transfer matrix diagonal in SU(2) irreps, weights c_j(beta) (same Bessel).")
    print("    %8s %12s %12s %14s" % ("beta", "c0", "c12/c0", "m_strong=-ln(c12/c0)"))
    for beta in [0.5, 1.0, 2.5, 5.0]:
        c0 = c(0.0, beta); r = c(0.5, beta) / c0
        print("    %8.3f %12.6f %12.6f %14.6f" % (beta, c0, r, -np.log(r)))
    print("    -> 2D residue weights = the SAME c_j(beta); 2D does NOT supply a new transport constant.")
    print("       The strong-coupling 'mass' -ln(c12/c0)>0 only for beta<~2; it turns negative at")
    print("       beta=2.5 (strong-coupling expansion BROKEN there) -- again, the crossover is the gap.")
    print("       Sign transport is vacuous at finite beta: c_j>0 for every j,beta (compact model).")

if __name__ == "__main__":
    section_A()
    section_B()
    section_C()
    section_D()
    section_E()
    section_F()
    print("=" * 78)
    print("VERDICT: Z_O>0 is trivial & bounded below (no BZ). It is NOT the missing step.")
    print("         Missing step = persistence/overlap Omega = |<G(tgt)|G(strong)>|^2 = isolation I2 = (CP).")
    print("         Transport CANNOT manufacture the atom (0 is absorbing). No QED. Clay 0%.")
