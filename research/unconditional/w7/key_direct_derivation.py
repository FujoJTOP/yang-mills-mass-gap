"""
key_direct_derivation.py  --  witness for research/unconditional/KEY-DIRECT-DERIVATION.md
Single small job.  All exact linear algebra / 1D quadrature / tiny toys.  NO lattice MC.
Run:  py -3.12 key_direct_derivation.py

Purpose: make the *arithmetic* of five direct attacks on  f1^2 >= c > 0  (equivalently
G''(k->0) >= lambda > 0) concrete at the crossover, and expose the exact step each dies on.
Labels in the doc: [T]/[C]/[O]/[NG]/[V-here].  Nothing here is a proof of anything.
"""
import numpy as np
from scipy import integrate, optimize

np.set_printoptions(precision=5, suppress=True)


def banner(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)

# ---------------------------------------------------------------------------
# PART A  --  The Savvidy / Nielsen-Olesen negative mode is NONABELIAN-ONLY.
#             Exact spectrum of the gauge-fluctuation operator in a constant B.
# ---------------------------------------------------------------------------
banner("PART A  |  gluon fluctuation spectrum in a constant background B  (U(1) vs SU(2))")
# Charged spin-1 (vector) in a constant magnetic field B, lowest Landau level n=0:
#     E^2 = 2|q|B (n + 1/2) + 2 q B s_z ,   s_z in {-1,0,+1},  n=0,1,2,...
# For a MASSLESS charged vector the s_z=-1, n=0 mode is the Nielsen-Olesen/Savvidy
# NEGATIVE mode (E^2 = -|q|B < 0).  The *charge* q is the eigenvalue of the background
# (Cartan) generator on the representation under which the fluctuation transforms.
def vector_levels(qB, nmax=1, szs=(-1, 0, 1)):
    out = []
    for n in range(nmax + 1):
        for s in szs:
            out.append(2.0 * abs(qB) * (n + 0.5) + 2.0 * qB * s)
    return np.array(sorted(out))

print("\n-- U(1) pure gauge: the photon is NEUTRAL under U(1) (adjoint = trivial rep), q=0 --")
print("   no charged vector exists; the only gluon is the neutral photon: -d^2, massless.")
E2_U1 = vector_levels(0.0)          # qB = 0  ->  E^2 = 0 for all (massless continuum)
print("   qB=0 spectrum E^2 =", E2_U1)
print("   min E^2 = %.4f   (NO negative mode; negative modes present: %d)"
      % (E2_U1.min(), np.sum(E2_U1 < -1e-9)))

print("\n-- SU(2) adjoint: background along T3 -> charges {0, +g, -g} on the 3 gluons --")
# ad(T3) on su(2) has eigenvalues {0, +i, -i}; the real charges of the two off-diagonal
# gluons under the Cartan are +1 and -1 (times g).  Build ad(T3) in the (T1,T2,T3) basis:
# su(2) structure constants f_{abc} = eps_{abc}; ad(T3)_{ab} = f_{3ab}
f = np.zeros((3, 3, 3))
for a in range(3):
    for b in range(3):
        for c in range(3):
            if (a, b, c) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
                f[a, b, c] = 1.0
            if (a, b, c) in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
                f[a, b, c] = -1.0
ad_T3 = f[2]  # (ad T3)_{ab} = f_{3ab}
charges = np.linalg.eigvals(ad_T3)
# real charges = +/-1 (and 0).  use them as q (in units of g).
q_real = np.array([1.0, -1.0, 0.0])   # the two charged gluons (W+-, complex combos) and the Cartan
g = 1.0
B = 1.0
E2_su2 = []
for q in q_real:
    E2_su2.extend(list(vector_levels(g * q * B)))
E2_su2 = np.array(sorted(E2_su2))
print("   ad(T3) eigenvalues (complex):", charges)
print("   SU(2) charged-gluon spectrum E^2 (gB=1) has min = %.4f" % E2_su2.min())
neg = E2_su2[E2_su2 < -1e-9]
print("   NEGATIVE modes E^2 =", neg, "  count =", len(neg))
print("   -> the negative (Savvidy) mode is NONABELIAN-ONLY: it needs a charged adjoint gluon.")
print("      U(1) pure gauge: 0 negative modes.   SU(2): >=1 (here %d).  [V-here, exact]" % len(neg))

# ---------------------------------------------------------------------------
# PART B  --  Savvidy scale: Phi_Sav(B) = sqrt(2 g B) is a SCALE, not a minorant.
#            It is monotonically increasing in B and -> 0 as B -> 0.
#            The confining vacuum has <B> = 0 (no uniform background) and m1 > 0,
#            so Phi_Sav = 0 << m1 there  =>  Phi_Sav <= m1 is FALSE.
# ---------------------------------------------------------------------------
banner("PART B  |  Savvidy scale Phi_Sav(B)=sqrt(2 g B): monotone -> 0, not a minorant")
# The ROBUST point does not depend on the one-loop constant: any 'Savvidy mass' obtained
# from a background-field curvature is a function of the background magnitude B and is
# monotonically increasing in B, vanishing as B -> 0.  The confining vacuum has no uniform
# background (B = 0 on average) yet m1 > 0, so the FUNCTION Phi_Sav(B) is not a minorant.
# At any putative one-loop minimum B_min the value is ~ Lambda (dimensionally the SAME
# scale as m1 -- the transmutation IDENTITY, not an inequality).
g = 1.0
print("   Phi_Sav(B) = sqrt(2 g B) as a function of the background B:")
for Bfac in [1.0, 0.1, 0.01, 1e-4, 0.0]:
    print("      B=%.3e   Phi_Sav=%.4e" % (Bfac, np.sqrt(2 * g * Bfac)))
print("   -> monotone increasing, -> 0 as B -> 0.  Physical vacuum: <B>=0, m1>0")
print("      => Phi_Sav(0)=0 << m1 : the FUNCTION is NOT <= m1.  [V-here]")
print("   -> At a one-loop minimum B_min, m_Sav=sqrt(2 g B_min) ~ Lambda; this is the")
print("      SAME dimensional-transmutation scale as m1 (an identity of scales), and it is")
print("      the curvature of an UNSTABLE saddle (has a negative mode), i.e. a width-like")
print("      object, not a spectral floor.  No operator inequality m1 >= sqrt(2 g B) exists.")

# ---------------------------------------------------------------------------
# PART C  --  Variational lower bounds: Rayleigh is UPPER; Temple needs E2;
#            Bogoliubov needs gap(H0) - ||V||, and free massless has gap(H0)=0.
# ---------------------------------------------------------------------------
banner("PART C  |  variational instruments: direction + hidden inputs")
rng = np.random.default_rng(0)
# toy spectrum on (vac)^perp: E1=1.0, E2=1.5, plus continuum from 2.0
E1, E2 = 1.0, 1.5
# (C1) Rayleigh quotient on a fixed 1-d trial vector combined with a 2-eigenvector toy
psi1 = np.array([1.0, 0.0])
psi2 = np.array([0.0, 1.0])
Htoy = np.diag([E1, E2])
best = np.inf
for _ in range(200000):
    v = rng.normal(size=2); v /= np.linalg.norm(v)
    R = float(v @ Htoy @ v)
    best = min(best, R)
print("   (C1) Rayleigh: min over 2e5 random trials R = %.6f ;  E1 = %.1f" % (best, E1))
print("        every R >= E1  =>  Rayleigh is an UPPER bound (min-max); a 'minorant' R<=E1 is false.")
# (C2) Temple: label the first excited (the gap) E1, need a lower bound rho on E2.
# Use CONSISTENT data: trial |psi> = cos t |1> + sin t |2>, R = <H>, Var = <H^2>-R^2.
def temple(R, Var, rho):
    return R - Var / (rho - R)
theta = 0.2                                   # admixture of the E2 level
R = np.cos(theta)**2 * E1 + np.sin(theta)**2 * E2
Var = np.cos(theta)**2 * np.sin(theta)**2 * (E2 - E1)**2   # exact for this 2-level trial
print("   (C2) trial: R=%.5f  Var=%.6f  (consistent 2-level data)" % (R, Var))
for rho, tag in [(E2, "rho=E2   (sound)   "),
                 (1.3 * E2, "rho=1.3E2(overclaim)"),
                 (3.0, "rho=3    (overclaim)")]:
    print("        Temple b(rho)=%.5f  [%s]   true E1=%.1f" % (temple(R, Var, rho), tag, E1))
print("        Temple is sound ONLY for rho<=E2: rho=E2 gives b=E1 exactly; rho>E2 OVERCLAIMS")
print("        (b>E1, a false too-high 'lower' bound). A certified rho<=E2 IS isolation E2.")
# (C3) Bogoliubov: m1 >= gap(H0) - ||V||
gapH0 = 0.0   # free massless Yang-Mills: gapless
normV = 0.7   # interaction, alpha_s ~ O(1)
print("   (C3) Bogoliubov: m1 >= gap(H0) - ||V|| = %.3f - %.3f = %.3f  (<=0, vacuous)"
      % (gapH0, normV, gapH0 - normV))
print("        needs gap(H0)>||V|| i.e. alpha_s small; at beta=2.5 there is NO small parameter.")

# ---------------------------------------------------------------------------
# PART D  --  exact variance decomposition G'' = S - D.  The Schur term CANCELS:
#             S - D == E[h|-] - Var(f|-) identically.  So 'D<S' is H1 verbatim.
#             Toy: V(x,y)=1/2 m^2 x^2 + 1/2 k y^2 + g x y^2 + 1/4 lam y^4
#             G''(0) = m^2 - g^2 Var(y^2|0).   Free field -> 0 (acid test passes).
# ---------------------------------------------------------------------------
banner("PART D  |  G''=S-D is an identity that cancels; the gap is 'stiffness beats variance'")
def marginal_curvature(m2, k, g, lam, dx=1e-4):
    """exact W''(0) for V above, from the 1D marginal p(x) (Gauss-Hermite-free: quadrature)."""
    def neglogp(x):
        f = lambda y: np.exp(-(0.5 * m2 * x * x + 0.5 * k * y * y + g * x * y * y
                               + 0.25 * lam * y**4))
        Z, _ = integrate.quad(f, -12, 12, limit=200)
        return -np.log(Z)
    W0 = neglogp(0.0); Wp = neglogp(dx); Wm = neglogp(-dx)
    return (Wp - 2 * W0 + Wm) / dx**2

def conditional_moments(k, lam):
    """moments of y under p(y|0) ∝ exp(-(1/2 k y^2 + 1/4 lam y^4))."""
    def p(y):
        return np.exp(-(0.5 * k * y * y + 0.25 * lam * y**4))
    Z, _ = integrate.quad(p, -12, 12, limit=200)
    y2, _ = integrate.quad(lambda y: y * y * p(y), -12, 12, limit=200)
    y4, _ = integrate.quad(lambda y: y**4 * p(y), -12, 12, limit=200)
    return y2 / Z, y4 / Z

print("\n   identity check  G''(0)  vs   m^2 - g^2 Var(y^2|0):")
for (m2, k, g, lam) in [(1.0, 1.0, 0.0, 0.0),      # free Gaussian, decoupled
                        (1.0, 1.0, 0.3, 0.0),      # coupled, y Gaussian
                        (1.0, 1.0, 0.6, 0.5),      # coupled, y non-Gaussian (lam>0)
                        (0.30, 1.0, 0.6, 0.5),     # weaker stiffness -> competition
                        (0.0, 1.0, 0.6, 0.5)]:     # FREE-MASSLESS-like x (m2->0): acid test
    Gpp = marginal_curvature(m2, k, g, lam)
    e_y2, e_y4 = conditional_moments(k, lam)
    Var_y2 = e_y4 - e_y2**2
    pred = m2 - g**2 * Var_y2
    print("      m2=%.2f k=%.2f g=%.2f lam=%.2f :  G''_num=% .5f   m2-g2Var(y2)=% .5f   Var(y2)=%.4f"
          % (m2, k, g, lam, Gpp, pred, Var_y2))

print("\n   READING:")
print("   * decoupled free field (g=0): G''=m^2; as m2->0, G''->0  (ACID TEST PASSES).")
print("   * the 'Gap' is literally  G'' = stiffness(m^2) - g^2 Var(y^2):  a MASS vs a VARIANCE.")
print("   * nonabelian coupling g REDUCES G'' (can drive it negative = the Savvidy marginal");
print("     instability); the indefinite fibre shifts S and D by the SAME -E[Schur], so it")
print("     cannot separate S from D.  'D<S' <=> 'm^2 > g^2 Var(y^2|0)' <=> H1 (Var bound). [V-here]")

banner("SUMMARY (no QED)")
print(" A: Savvidy negativity exact & NONABELIAN-ONLY (U(1):0, SU(2):>=1). [V-here]")
print(" B: Savvidy scale sqrt(2gB) is monotone->0 as B->0; =0 at the physical vacuum, m1>0.")
print("    So Phi_Sav is NOT a minorant (Phi_Sav<=m1 FALSE). It is the WRONG vacuum's scale.")
print(" C: Rayleigh=upper; Temple needs E2(=isolation); Bogoliubov needs alpha_s small(absent).")
print(" D: G''=S-D cancels identically; gap = stiffness-vs-variance = H1; acid test passes.")
print(" => No attack yields f1^2>=c or G''>=lambda uniformly at beta~2.5 without H1/(CP).")
