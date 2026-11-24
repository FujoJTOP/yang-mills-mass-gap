"""
invent_curvature.py  --  witness for INVENT-CURVATURE (2026-11-24)

Question: is a thermodynamic-geometry (Fisher-Rao / Ruppeiner) curvature functional
on the family {mu_beta} a genuinely NEW strut for the 0++ gap, or a relabelling of chi?

Pieces (all cheap, numpy only):
 (A) 1-parameter family -> scalar curvature is IDENTICALLY ZERO (1D flat). Pure YM at
     fixed volume has one coupling => the naive curvature of {mu_beta} is vacuous.
 (B) Gaussian family (free field, massless OR massive): Fisher curvature is CONSTANT
     -1/2 (scale-independent) => cannot distinguish massless from massive => acid-vacuous.
 (C) Fisher = susceptibility identity: for family mu_h ~ e^{-S + h O}, g_hh = d^2 logZ/dh^2
     = Var(O) = chi.  So the metric literally contains the susceptibility; any 2-param
     curvature R is a function of chi and its slopes => homogeneous, the same wall.
 (D) Non-Gaussian family (quartic): R deviates from -1/2, SIGN-INDEFINITE, and the
     deviation = non-Gaussianity (READ not FORCE).

Run: py -3.12 -I invent_curvature.py
"""
import numpy as np

# ---------------------------------------------------------------------------
# general 2D scalar curvature from a metric  g(u,v) = [[E,F],[F,G]] by finite diff
# ---------------------------------------------------------------------------
def curvature(gfun, u, v, h=2e-3):
    def g(uu, vv): return np.array(gfun(uu, vv), float)
    def dg(uu, vv, ax):
        du = h if ax == 0 else 0.0; dv = h if ax == 1 else 0.0
        return (g(uu+du, vv+dv) - g(uu-du, vv-dv)) / (2*h)
    def Gam(uu, vv):
        G = np.linalg.inv(g(uu, vv)); D = [dg(uu, vv, 0), dg(uu, vv, 1)]
        Gm = np.zeros((2, 2, 2))
        for k in range(2):
            for i in range(2):
                for j in range(2):
                    Gm[k, i, j] = 0.5*sum(G[k, l]*(D[i][j, l]+D[j][i, l]-D[l][i, j]) for l in range(2))
        return Gm
    G0 = Gam(u, v)
    dG = [(Gam(u+h, v)-Gam(u-h, v))/(2*h), (Gam(u, v+h)-Gam(u, v-h))/(2*h)]
    Riem = np.zeros((2, 2, 2, 2))
    for l in range(2):
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    Riem[l, i, j, k] = dG[i][l, j, k]-dG[j][l, i, k] + sum(G0[l, i, m]*G0[m, j, k]-G0[l, j, m]*G0[m, i, k] for m in range(2))
    Ricci = np.array([[sum(Riem[k, i, k, j] for k in range(2)) for j in range(2)] for i in range(2)])
    Gup = np.linalg.inv(g(u, v))
    # convention: Poincare half-plane (K=-1) must give scalar R=-2 => negate to match
    return -float(sum(Gup[i, j]*Ricci[i, j] for i in range(2) for j in range(2)))


def banner(s):
    print("\n" + "="*76); print(s); print("="*76)


# ---------------------------------------------------------------------------
banner("(A) 1-parameter family: scalar curvature == 0 (1D flat)")
print("  A 1x1 metric has no curvature tensor.  Pure YM (one coupling beta,")
print("  fixed volume) => the naive Fisher-Rao/Ruppeiner curvature of {mu_beta} is 0.")
print("  => VACUOUS as a strut (nothing to evaluate).   [T]")

# ---------------------------------------------------------------------------
banner("(B) Gaussian family: scalar Fisher curvature R == -1 CONSTANT (Gauss K == -1/2)")
# g = diag(1/v^2, 2/v^2),  u = mean mu, v = sigma.  Poincare half-plane form.
def gauss_metric(u, v):
    return [[1.0/v**2, 0.0], [0.0, 2.0/v**2]]
for (v, tag) in [(1.0, "sigma=1  (massive-like)"), (0.1, "sigma=0.1"), (10.0, "sigma=10 (approach massless: chi->inf)")]:
    R = curvature(gauss_metric, 0.0, v)
    print(f"  Gaussian {tag:34s}: scalar R = {R:+.5f}")
print("  => scalar curvature R == -1 (Gauss K == -1/2) for EVERY Gaussian")
print("     (massless AND massive).  It does NOT distinguish gapless from gapped")
print("     free field => ACID-VACUOUS.   [T]")
print("     (Analytic: Bonnet formula gives K = -1/2 exactly; scalar R = 2K = -1.)")

# ---------------------------------------------------------------------------
banner("(C) Fisher == susceptibility:  g_hh = d^2 logZ/dh^2 = Var(O) = chi")
# measure  p(x) ~ exp(-a x^2/2 - b x^4/4 + h x^2);  O = x^2;  chi = Var(O)
def logZ(a, b, h, L=12.0, n=600001):
    x = np.linspace(-L, L, n)
    lf = -a*x**2/2 - b*x**4/4 + h*x**2
    m = lf.max(); Z = np.exp(lf-m)
    return m + np.log(np.trapezoid(Z, x))
def VarO(a, b, h, L=12.0, n=600001):
    x = np.linspace(-L, L, n)
    lf = -a*x**2/2 - b*x**4/4 + h*x**2
    m = lf.max(); w = np.exp(lf-m); w /= np.trapezoid(w, x)
    O = x**2; mO = np.trapezoid(w*O, x)
    return np.trapezoid(w*(O-mO)**2, x)
for (a, b) in [(1.0, 0.0), (1.0, 0.08)]:
    hh = 1e-4
    g_hh = (logZ(a, b, hh) - 2*logZ(a, b, 0.0) + logZ(a, b, -hh))/hh**2
    chi = VarO(a, b, 0.0)
    print(f"  a={a}, b={b}:  g_hh = {g_hh:.6f}   Var(O)=chi = {chi:.6f}   diff={abs(g_hh-chi):.2e}")
print("  => the metric component in the source direction IS chi exactly.     [T]")
print("     A 2-parameter curvature R is built from this metric (incl. chi) +")
print("     its slopes => bounding |R| requires controlling chi == H1. HOMOGENEOUS.")

# ---------------------------------------------------------------------------
banner("(D) R depends on the CHOSEN 2-parameter family, not on mu alone (deeper issue)")
# exponential family p ~ exp(-a x^2/2 - b x^4/4); natural params (a,b); Fisher = Hess psi
def psi(a, b, L=12.0, n=400001):
    x = np.linspace(-L, L, n)
    lf = -a*x**2/2 - b*x**4/4
    m = lf.max(); Z = np.exp(lf-m)
    return m + np.log(np.trapezoid(Z, x))
def quartic_metric(a, b, d=1e-3):
    # Hessian of psi wrt (a,b)
    pa = lambda A, B: psi(A, B)
    g_aa = (pa(a+d, b) - 2*pa(a, b) + pa(a-d, b))/d**2
    g_bb = (pa(a, b+d) - 2*pa(a, b) + pa(a, b-d))/d**2
    g_ab = (pa(a+d, b+d) - pa(a+d, b-d) - pa(a-d, b+d) + pa(a-d, b-d))/(4*d**2)
    return [[g_aa, g_ab], [g_ab, g_bb]]
print("  scalar R at fixed a=1, varying non-Gaussianity b:")
print("  (b=0 slice IS a Gaussian in x, but a DIFFERENT 2D family than (mu,sigma):")
print("   its Gaussian-point curvature is NOT -1 -> R depends on the family, not mu.)")
R0 = None
for b in [0.0, 0.02, 0.05, 0.1, 0.2, 0.5]:
    R = curvature(lambda A, B: quartic_metric(A, B), 1.0, b)
    if R0 is None: R0 = R
    print(f"    b={b:4.2f}:  scalar R = {R:+.5f}    (R-R(b=0) = {R-R0:+.5f})")
print("  => (i) the Gaussian point of the (a,b) family gives R~-6.9, while the")
print("     (mu,sigma) Gaussian family gives R=-1: TWO embeddings of a Gaussian")
print("     measure give DIFFERENT curvature => R is NOT a functional of mu alone")
print("     (it needs a choice of 2 thermodynamic variables).")
print("     (ii) the drift from R(b=0) tracks NON-GAUSSIANITY (READ, not FORCE).")
print("     => not a well-defined mu-functional (P1) and not sign-definite. [NG]")

# ---------------------------------------------------------------------------
banner("VERDICT")
print("  (A) 1-param  -> curvature == 0 (flat): vacuous.")
print("  (B) Gaussian -> scalar R == -1 constant: acid-vacuous (massless == massive).")
print("  (C) g_hh == chi exactly: 2-param curvature is a function of chi => HOMOGENEOUS.")
print("  (D) R depends on the CHOSEN family, not mu alone (two Gaussian embeddings give")
print("      R=-1 vs R~-6.9); its drift is read non-Gaussianity, not sign-definite.")
print("  Thermodynamic curvature is NOT a kappa-independent sign-definite mu-functional;")
print("  it is (flat 1-param) or (homogeneous in chi) or (embedding-dependent).")
print("  Exact failing step: (C)+(D) -- R inherits chi, and needs a variable choice;")
print("  bounding R <=> bounding chi = H1.   Clay 0%.  No QED.")
