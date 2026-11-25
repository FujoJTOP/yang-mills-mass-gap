#!/usr/bin/env python3
"""
T2-TUBE-OVERLAP witness.  Deterministic, numpy only.  Run:  py -3.12 -I t2_tube_overlap.py

Object under study: the *normalised flux-tube overlap*
    w1 := f1^2 / M0 ,  f1^2 = |<G|O|Omega>|^2 ,  M0 = Var(O) = ||phi||^2 ,  phi = Q O|Omega> .
w1 in [0,1] is the fraction of the source's spectral weight sitting in the isolated
lightest tube.  It is dimensionless and invariant under the multiplicative transport
O -> Z_O O (numerator and M0 both scale by Z_O^2).  Honest, scheme-free form of the floor.

What this witness shows (NOT a proof of anything about YM; a check of the inequality
structure and of the numerical facts quoted in T2-TUBE-OVERLAP.md):
 (A) w1 is a normalised projection fraction  (exact identity, KL/RP measure dnu>=0).
 (B) finite-t estimator  w1hat(t) = e^{E1 t} Cbar(t)  is a DECREASING UPPER bound on w1;
     it collapses *toward* w1 with tube length, never below -> finite t cannot certify a floor.
 (C) reflection positivity / log-convexity of C(t) holds ALSO when w1 = 0
     (pure cut) -> RP gives no lower bound on the overlap (acid-test-blind).
 (D) threshold/dissolution: atom fraction w1 -> 0 as the bound state reaches threshold;
     free field (composite cut): w1 = 0 exactly (acid test).
 (E) residue-level-gap lower bound requires the second level Lambda = E2; drop it -> 0.
 (F) strong-coupling overlap floor a0 >= (2N)^-2 and the KP radius uK at beta=2.5.
"""
import numpy as np

np.set_printoptions(precision=6, suppress=True, linewidth=140)


def bar_c(t, E, w):
    """normalised connected correlator Cbar(t)=sum w_n e^{-E_n t}/sum w_n, w_n>=0 (KL)."""
    E = np.asarray(E, float); w = np.asarray(w, float)
    return np.sum(w[None, :] * np.exp(-E[None, :] * np.atleast_1d(t)[:, None]), axis=1) / w.sum()


def eff_mass(t, E, w):
    """m_eff(t) = -(d/dt) log Cbar(t)  (RP: decreasing from above)."""
    C = bar_c(t, E, w)
    dC = -np.sum(np.asarray(w)[None, :] * np.asarray(E)[None, :] *
                 np.exp(-np.asarray(E)[None, :] * np.atleast_1d(t)[:, None]), axis=1) / np.sum(w)
    return -dC / C


def second_diff(t, f):
    """discrete second difference of f on grid t (log-convexity test: >=0)."""
    d2 = f[:-2] - 2 * f[1:-1] + f[2:]
    # scale to a curvature per unit t^2 (grid uniform)
    h = t[1] - t[0]
    return d2 / h**2


print("=" * 78)
print("T2-TUBE-OVERLAP witness  (numpy %s)" % np.__version__)
print("=" * 78)

# ----------------------------------------------------------------------------
# (A) w1 is a normalised projection fraction, in [0,1].  Any KL measure.
# ----------------------------------------------------------------------------
E = np.array([1.00, 2.00, 2.10, 2.25, 2.45, 2.70, 3.00])
w = np.array([0.30, 0.20, 0.15, 0.13, 0.10, 0.07, 0.05])
w1_true = w[0] / w.sum()
print("\n(A) exact identity: w1 = weight_lightest / total weight  in [0,1]")
print("    spectrum atoms (E_n, w_n):", list(zip(E.tolist(), w.tolist())))
print("    w1_true = w_1 / sum_n w_n = %.6f   (0<=w1<=1 exact; sum_n w_n = M0)"
      % w1_true)

# ----------------------------------------------------------------------------
# (B) finite-t estimator is a DECREASING UPPER bound on w1.
# ----------------------------------------------------------------------------
t = np.arange(0.0, 60.0 + 1e-9, 2.0)
w1hat = np.exp(E[0] * t) * bar_c(t, E, w)
mono = np.all(np.diff(w1hat) <= 1e-12)
print("\n(B) finite-t overlap estimator  w1hat(t) = e^{E1 t} Cbar(t)  (E1=%.2f):" % E[0])
for i in range(0, len(t), 4):
    print("    t=%5.1f   w1hat=%.6f   (>= w1_true=%.6f ? %s)"
          % (t[i], w1hat[i], w1_true, w1hat[i] >= w1_true - 1e-9))
print("    monotone non-increasing in t : %s" % mono)
print("    limit w1hat(inf)             = %.6f  == w1_true : %s"
      % (w1hat[-1], abs(w1hat[-1] - w1_true) < 1e-3))
print("    READ: every finite tube length OVER-estimates the overlap; the true value")
print("          is the *limit*, never certified from finite t (upper bound only).")

# ----------------------------------------------------------------------------
# (C) reflection positivity / log-convexity holds also when w1 = 0.
# ----------------------------------------------------------------------------
print("\n(C) RP / log-convexity of Cbar(t)  (must hold, dnu>=0), for w1>0 AND w1=0:")
cases = [("w1>0 (atom + cut)", E, w),
         ("w1=0 (pure cut)   ", np.array([2.00, 2.10, 2.25, 2.45, 2.70, 3.00]),
                                np.array([0.20, 0.18, 0.20, 0.20, 0.12, 0.10]))]
tt = np.arange(0.0, 40.0 + 1e-9, 2.0)
for name, Ec, wc in cases:
    lc = np.log(bar_c(tt, Ec, wc))
    curv = second_diff(tt, lc)
    meff = eff_mass(tt, Ec, wc)
    print("    %s :  log Cbar convex? min 2nd-diff=%+.3e (>=0 => convex); "
          "m_eff decreasing? %s ; m_eff(0)=%.3f -> m_eff(end)=%.3f"
          % (name, curv.min(), np.all(np.diff(meff) <= 1e-9), meff[0], meff[-1]))
print("    READ: both measures are legal (RP/KL). log Cbar is convex in BOTH; RP pins")
print("          m_eff decreasing from above but is blind to whether w1 is 0. acid test.")

# ----------------------------------------------------------------------------
# (D) threshold / dissolution: atom fraction -> 0 as the level reaches threshold.
#     Threshold behaviour of a Feshbach residue: Z = kappa*sqrt(delta), delta=gap.
# ----------------------------------------------------------------------------
print("\n(D) dissolution at threshold: dnu = Z*delta_{E1} + A*sqrt(s-theta) ds, "
      "theta=E1+delta,")
print("    Z = kappa*sqrt(delta) (threshold residue).  A=1, W=2 (cut width), kappa=1, E1=1:")
theta0, W, A, kappa = 1.0, 2.0, 1.0, 1.0
def cut_mass(theta, W, A):
    # integral_theta^{theta+W} A sqrt(s-theta) ds = A * (2/3) W^{3/2}
    return A * (2.0 / 3.0) * W ** 1.5
print("      delta        Z        cut_mass     w1=Z/(Z+cut)")
for delta in [1.0, 0.5, 0.1, 1e-2, 1e-3, 1e-4, 1e-6]:
    Z = kappa * np.sqrt(delta)
    cm = cut_mass(theta0 + delta, W, A)
    print("    %9.1e  %9.4e  %9.4f    %10.6f" % (delta, Z, cm, Z / (Z + cm)))
print("    READ: w1 -> 0 continuously as binding delta -> 0; the measure stays RP.")
print("          NO floor survives the dissolution limit. (R6: f1^2 upper-semicontinuous.)")

# free field: composite tr F^2 is a two-particle CUT, no point mass.
print("\n    free field acid test: composite tr F^2 has dnu = rho(s) ds on [4m^2,inf),")
print("    no atom => w1 = nu({m^2}) = 0 exactly (independent of m).")

# ----------------------------------------------------------------------------
# (E) residue-level-gap lower bound requires Lambda=E2.
# ----------------------------------------------------------------------------
print("\n(E) residue-level-gap:  w1 >= (Lambda*M0 - M1)/(M0*(Lambda-m1^2)),  needs Lambda.")
E2 = np.array([1.0, 2.0, 2.5])
w2 = np.array([0.4, 0.35, 0.25])
M0 = w2.sum(); M1 = np.sum(w2 * E2); m1sq = E2[0] ** 2; w1t = w2[0] / M0
# interpret atoms as masses^2: use s_n = E_n (as if these are m^2 values)
M0 = w2.sum(); M1 = np.sum(w2 * E2); m1sq = E2[0]; w1t = w2[0] / M0
print("    true  w1 = %.6f" % w1t)
for Lam, tag in [(E2[1], "Lambda=E2 (sound)"), (E2[0] + 1e-6, "Lambda=m1^2+ (no 2nd level)"),
                 (3.0, "Lambda=3.0 (>E2, unsound)")]:
    bound = (Lam * M0 - M1) / (M0 * (Lam - m1sq))
    print("      %-26s -> bound w1 >= %+.6f %s"
          % (tag, bound, "" if tag.startswith("Lambda=E2") else
             "(vacuous/loose)" if bound <= w1t + 1e-9 else "(OVERCLAIMS!)"))
print("    READ: the only true lower bound CONSUMES Lambda=E2; without a second level")
print("          it collapses to 0. Same circularity, in overlap language.")

# ----------------------------------------------------------------------------
# (F) strong-coupling overlap floor a0>0 and KP radius at beta=2.5.
# ----------------------------------------------------------------------------
print("\n(F) strong coupling (SU(2) example): a0 >= (2N)^-2 ; w1 = f1^2/M0 >= w0 > 0.")
N = 2
a0 = (2 * N) ** (-2)
print("    N=2: a0 >= (2N)^-2 = %.6f > 0   (positive recoupling, t- and L-independent)" % a0)
print("    KP convergence needs  u*K < 1 ,  K <= 6.1e2 * N = %.0f ;  u=I2/I1(measured)" %
      (6.1e2 * N))
data = [(0.01, 5.0e-3), (1.0 / 48, 1.04e-2), (1.0, 4.33e-1), (2.5, 7.19e-1)]
K = 6.1e2 * N
print("      beta        u=I2/I1       u*K          inside KP radius u*K<1 ?")
for beta, u in data:
    print("    %7.4f   %10.4e   %10.1f      %s" % (beta, u, u * K, u * K < 1.0))
# beta* where crude u*K=1 (using SU(2) u~beta/4 leading order)
u_at = 1.0 / K
beta_star_crude = 4.0 * u_at
print("    crude radius (SU(2) u~beta/4): u=1/K=%.2e -> beta* ~ %.2e  "
      "(rigorous P3: beta1=1/48=%.4f)" % (u_at, beta_star_crude, 1.0 / 48))
print("    READ: at beta=2.5, u*K ~ %.0f >> 1  ==> the cluster expansion is OUTSIDE its" % (0.719 * K))
print("          convergence radius; the strong-coupling floor is not even defined there.")

print("\n" + "=" * 78)
print("Verdict of the witness (structural, not a YM proof):")
print("  * w1 = f1^2/M0 in [0,1] is a clean, Z_O-invariant overlap fraction.")
print("  * RP/log-convexity and PF give NO lower bound on w1 (acid-test-blind).")
print("  * finite-t (finite tube length) estimators are UPPER bounds: the residue is")
print("    t-independent, but its CERTIFICATION is not -- no finite t or L certifies w1>0.")
print("  * at dissolution (binding -> 0) w1 -> 0 continuously; free field w1 = 0 exactly.")
print("  * the only proven floor is strong coupling; at beta=2.5 u*K >> 1, outside radius.")
print("=" * 78)
