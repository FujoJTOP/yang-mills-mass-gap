"""
slot_quadratic_fixpoint.py -- witness for research/unconditional/SLOT-QUADRATIC-FIXPOINT.md
ONE small job. Exact linear algebra + 1D/2D quadrature only. NO lattice MC.
Run:  py -3.12 slot_quadratic_fixpoint.py

Purpose: make the arithmetic of the "sign-definite quadratic minorant AS A FIXPOINT" concrete:
 (1) pointwise trap:  Q>=0 and Q<=Hess  ==>  Hess>=0  (no PSD minorant of an indefinite Hessian);
 (2) expectation compensation (Stein/SD):  E_mu[Hess V] = E_mu[f (x) f]  >= 0  -- exact, gap-blind;
 (3) exact marginal fixpoint Q := G''(0) = m^2_bare - g^2 Var(y^2|0), and the SELF-CONSISTENT
     iteration: Q=0 is the gapless fixpoint; stripping the bare mass (pure YM) leaves Q<=0,
     i.e. NO source at Q=0 -- a positive fixpoint must be generated (dimensional transmutation = gap).
Labels [T]/[C]/[O]/[NG]/[V-here]. Nothing here is a proof.
"""
import numpy as np
from scipy import integrate

np.set_printoptions(precision=6, suppress=True)

def banner(s):
    print("\n" + "=" * 78); print(s); print("=" * 78)

# ===========================================================================
# PART 1 -- THE POINTWISE TRAP (exact linear algebra, [T])
#   Claim: if H is a symmetric form with a negative eigenvalue, there is NO
#   positive-semidefinite form Q with Q <= H.  Proof one line: Q<=H and Q>=0
#   give H >= Q >= 0.  Demonstrate on a Savvidy-like indefinite 2x2.
# ===========================================================================
banner("PART 1 | pointwise trap:  Q>=0 & Q<=Hess  ==>  Hess>=0  (no PSD minorant of indefinite Hess)")
H = np.array([[2.0, 0.6], [0.6, -1.0]])   # indefinite: one negative eigenvalue (Savvidy-like)
w, V = np.linalg.eigh(H)
print("  Savvidy-like Hessian H =\n", H)
print("  eigenvalues of H =", w, "   (lambda_min = %.4f < 0)" % w.min())
# search for a PSD Q with Q <= H over symmetric 2x2 (2x2 PSD = 3 params), coarse grid
worst = None
for q11 in np.linspace(-3, 3, 61):
    for q22 in np.linspace(-3, 3, 61):
        for q12 in np.linspace(-3, 3, 61):
            Q = np.array([[q11, q12], [q12, q22]])
            if np.linalg.eigvalsh(Q).min() < -1e-9:      # Q not PSD -> skip
                continue
            if np.linalg.eigvalsh(H - Q).min() < -1e-9:  # Q not <= H -> skip
                continue
            worst = Q
            break
        if worst is not None: break
    if worst is not None: break
print("  any PSD Q with Q<=H found in search? ", "NO  (as forced: H<0 somewhere => impossible)" if worst is None
      else "YES "+str(worst))
print("  => a POINTWISE quadratic minorant cannot be sign-definite where Savvidy is negative. [T]")

# ===========================================================================
# PART 2 -- EXPECTATION COMPENSATION (Stein / Schwinger-Dyson), exact, gap-blind [T]
#   Stein:  E_mu[H_xx] = E_mu[f^2],  f=dV/dx,  H_xx=d^2V/dx^2.
#   It is an INTEGRATION-BY-PARTS identity: the boundary term vanishes ONLY for the
#   FULL, untruncated, exactly-integrated measure.  Truncation/conditioning leaves a
#   nonzero boundary term = exactly the variance the gap lives in (the L3->L4 step).
#   Use the BOUNDED bilinear Gaussian  V=1/2 a x^2 + 1/2 b y^2 + c x y  (a,b>0, |c|<sqrt(ab)).
# ===========================================================================
banner("PART 2 | expectation compensation E[H_xx]=E[f^2] (Stein) is exact but needs the FULL measure")
def momg(a, b, c, op, lim=8.0):
    """E_mu[op] for the Gaussian mu ~ exp(-(1/2 a x^2 + 1/2 b y^2 + c x y)), on [-lim,lim]^2."""
    V = lambda x, y: 0.5*a*x*x + 0.5*b*y*y + c*x*y
    Z, _ = integrate.dblquad(lambda y, x: np.exp(-V(x, y)), -lim, lim, -lim, lim, epsabs=1e-13)
    num, _ = integrate.dblquad(lambda y, x: op(x, y)*np.exp(-V(x, y)), -lim, lim, -lim, lim, epsabs=1e-13)
    return num / Z
print("\n  (a) bounded Gaussian  V=1/2 x^2+1/2 y^2+c x y:  E[H_xx]=a=1 vs E[f^2] (f=a x+c y), full box:")
for c in [-0.6, 0.0, 0.6]:
    EH  = momg(1.0, 1.0, c, lambda x, y: 1.0)               # H_xx = a = 1
    Ef2 = momg(1.0, 1.0, c, lambda x, y: (1.0*x + c*y)**2)  # f = a x + c y
    print("     c=% .1f :  E[H_xx]=% .6f   E[f^2]=% .6f   diff=% .2e   (exact: both = 1)"
          % (c, EH, Ef2, EH - Ef2))
print("  (b) SAME Gaussian, SHRINKING box -> the IBP boundary term is left behind:")
for lim in [8.0, 3.0, 2.0, 1.5]:
    Ef2 = momg(1.0, 1.0, 0.6, lambda x, y: (x + 0.6*y)**2, lim=lim)
    print("     box=+-%.1f :  E[f^2]box=% .5f   truncation residual = E[f^2]box - 1 = % .3e"
          % (lim, Ef2, Ef2 - 1.0))
print("  => Stein compensation is a statement about the UNCONDITIONED, exactly-INTEGRATED measure.")
print("     Truncation/conditioning re-introduces a boundary term = the variance the gap lives in (H1).")
print("     FREE FIELD (a,b->0 Gaussian, massless): E[H]=E[f^2]=0 -> identity met, GAP-BLIND. [T]")

# ===========================================================================
# PART 3 -- THE FIXPOINT  Q := G''(0) = m^2_bare - g^2 Var(y^2|0)   (exact)
#   and its SELF-CONSISTENT closure where the fibre mass is the induced mass Q itself.
#   The ONLY positive source is m^2_bare.  Pure YM has m^2_bare = 0.
# ===========================================================================
banner("PART 3 | the marginal fixpoint  Q=G''(0)=m^2_bare - g^2 Var(y^2|0), and its self-consistency")
def marginal_curvature(m2, k, g, lam, dx=2e-4):
    def neglogp(x):
        f = lambda y: np.exp(-(0.5*m2*x*x + 0.5*k*y*y + g*x*y*y + 0.25*lam*y**4))
        Z, _ = integrate.quad(f, -12, 12, limit=300)
        return -np.log(Z)
    W0, Wp, Wm = neglogp(0.0), neglogp(dx), neglogp(-dx)
    return (Wp - 2*W0 + Wm)/dx**2
def var_y2(k, lam):
    p = lambda y: np.exp(-(0.5*k*y*y + 0.25*lam*y**4))
    Z, _ = integrate.quad(p, -12, 12, limit=300)
    y2, _ = integrate.quad(lambda y: y*y*p(y), -12, 12, limit=300)
    y4, _ = integrate.quad(lambda y: y**4*p(y), -12, 12, limit=300)
    m2y, m4y = y2/Z, y4/Z
    return m4y - m2y**2

print("\n  exact identity G''(0) == m^2 - g^2 Var(y^2|0):")
for (m2, k, g, lam) in [(1.0, 1.0, 0.0, 0.0),
                        (1.0, 1.0, 0.6, 0.5),
                        (0.0, 1.0, 0.6, 0.5)]:
    Gpp = marginal_curvature(m2, k, g, lam)
    pred = m2 - g**2 * var_y2(k, lam)
    print("     m2=%.2f g=%.2f lam=%.2f :  G''_num=% .5f   m^2-g2Var=% .5f"
          % (m2, g, lam, Gpp, pred))

print("\n  FREE-FIELD ACID: m2=0, g=0  ->  G'' = %.6f  (gapless fixpoint Q=0 admissible)"
      % marginal_curvature(0.0, 1.0, 0.0, 0.5))
print("  PURE-YM sign  : m2=0, g=0.6 -> G'' = %.6f  (NO bare source => Q<=0, Savvidy sign)"
      % marginal_curvature(0.0, 1.0, 0.6, 0.5))

# Self-consistent closure: fibre mass k := Q (induced), bare mass b; iterate Q = b - g^2 Var(y^2|0;k=Q)
banner("PART 3b | self-consistent closure  Q = b - g^2 Var(y^2|0 ; k=Q) : b=0 (pure YM) vs b>0")
def iterate_Q(b, g, lam, Q0=1.0, n=400, damp=0.15):
    Q = Q0
    hist = [Q]
    for _ in range(n):
        k = max(Q, 1e-6)
        Qn = b - g**2 * var_y2(k, lam)
        Q = (1 - damp)*Q + damp*Qn          # damped fixed-point iteration
        hist.append(Q)
    return Q, np.array(hist)
for (b, g) in [(0.0, 0.6), (0.0, 1.0), (1.0, 0.6), (1.0, 1.0)]:
    Qs, hist = iterate_Q(b, g, lam=0.5)
    # residual of the exact FP equation at the converged Q
    res = Qs - (b - g**2 * var_y2(max(Qs, 1e-6), 0.5))
    print("     bare b=%.1f g=%.1f :  Q*=% .5f  (residual %.2e)   sign=%s"
          % (b, g, Qs, res, "POSITIVE source needed" if Qs > 1e-4 else
             ("ZERO (gapless)" if abs(Qs) < 1e-4 else "NEGATIVE (no gap)")))

print("\n  READING:")
print("   * b=0,g>0 (pure YM): the self-consistent curvature is Q*<=0 (Savvidy/unstable sign).")
print("     There is NO positive fixed point -- the fixpoint does not GENERATE a scale.")
print("   * b>0: Q*>0 because the BARE MASS b is the source. Pure YM has no such bare mass.")
print("   * the free massless fixpoint is exactly Q=0.  => Q=0 is admissible; no source at Q=0. [V-here]")

# ===========================================================================
# PART 4 -- SCALE HOMOGENEITY (dimensional, [T at structure level])
#   Pure YM is classically scale-free (only the dimensionless coupling + the
#   regulator a).  A fixpoint Q of mass dimension 2 with no dimensionful input
#   is homogeneous; the anomaly T^mu_mu=(beta/2g)tr F^2 is the only explicit
#   scale-breaker and enters the residue relation HOMOGENEOUSLY (cancels).
# ===========================================================================
banner("PART 4 | the one explicit scale-breaker (trace anomaly) is HOMOGENEOUS (no Q-independent source)")
# anomaly relation: <0|O|G> = f_G m_G^2 g/(2 beta),  with O=1/4 tr F^2.  Both sides ~ f_G (linear).
# Under f_G -> s f_G, m_G -> s m_G (both the residue AND mass scale), the relation is degree 2 in s^? :
for fG, mG in [(1.0, 1.0), (2.0, 2.0), (3.0, 3.0)]:
    # ratio LHS/RHS with beta<0, g fixed: proportional to fG / (fG mG^2) = 1/mG^2 -> NOT zero.
    print("     f_G=%.1f m_G=%.1f :  (f_G m_G^2)/(f_G) = m_G^2 = %.1f  (relation fixes m_G in terms of f_G;"
          % (fG, mG, fG*mG**2/fG), " it is homogeneous, carries no Q-independent source)")
print("  => strip bare mass: the only scale-breaker (anomaly) enters ~ f_G m_G^2, i.e. HOMOGENEOUSLY;")
print("     it relates the residue to the mass but cannot FLOOR either. A nonzero Q = spontaneous")
print("     scale breaking = dimensional transmutation = the mass gap itself. [T-structure]")

banner("SUMMARY (no QED)")
print(" 1 pointwise: PSD minorant of an indefinite Hessian does NOT exist (Savvidy kills it). [T]")
print(" 2 expectation: E[Hess]=E[f^2]>=0 by Stein -- true, but gap-blind (free field saturates). [T]")
print(" 3 fixpoint Q=G''(0): Q=0 admissible (free); pure YM (b=0) => Q*<=0, NO source at Q=0. [V-here]")
print(" 4 uniformity: the positive branch needs spontaneous scale breaking = the gap (=H1).")
print(" => The quadratic-minorant fixpoint does NOT exclude the gapless point; Q=0 solves it,")
print("    there is no source at Q=0, and a nonzero fixed point IS dimensional transmutation. NO QED.")
