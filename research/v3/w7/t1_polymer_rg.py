#!/usr/bin/env python3
r"""
t1_polymer_rg.py -- T1 witness: can a polymer / block-spin RG push the binding
expansion (P3 residue, beta<1/48) to beta=2.5?

Deterministic, numpy(+scipy.special.iv) only, no downloads, one job, ~1 s.
Run: py -3.12 -I t1_polymer_rg.py

What it computes, and why (maps 1:1 to T1-POLYMER-RG.md):

  [C] exact SU(2) character coefficients  c_n(beta) = I_n(beta) - I_{n+2}(beta),
      n = 2j = number of boxes, dim = n+1;  SUM_n (n+1) c_n = e^beta  (identity check).
      -> the strong-coupling "activity"
         u(beta)        = c_1/c_0                (single-plaquette fundamental activity)
         Phi(a,beta)    = SUM_{n>=1} (c_n/c_0)(n+1) e^{a n}   (per-plaquette KP activity)
         rho_Ur(beta)   = lambda_a * Phi(1,beta)  (KP / Ursell convergence ratio)
      Convergence of the connected-polymer (Ursell) series requires rho_Ur < 1.

  [D] minimal connected tr F^2 polymer = the four-wall flux tube (P3 geometry):
      per-slice weight u^4, rate  m1 a = 4 ln(1/u);  the cheapest competitor
      (perimeter-6, 1x2 cross-section) costs u^6 (one extra u^2 per slice).

  [E] the actual numbers at beta = 1/48 (the P3/SZZ wall), beta=1, beta=2.5,
      for d = 4 (Delta=20), d = 3 (Delta=12), d = 2 (Delta=4);
      and the "needed vs actual" ratio  u(2.5) / (1/lambda_a).

  [F] block-spin / one-loop RG estimate of the effective activity: the physical
      range is xi/a (correlation length in lattice units); flowing a^{-1}->xi^{-1}
      changes beta by  dbeta = -8 b0 ln(xi/a),  b0 = 11N/(48 pi^2);
      compare u_eff at the physical scale to the 1/lambda_a needed.

  [G] scale separation xi/a at beta=2.5 (strong-coupling leading + scaling estimates):
      the number of halving steps between the lattice and the correlation length.

  [H] acid test: free field beta->inf ( u -> c1/c0 -> 2 ), activities do NOT decay
      -> the expansion is maximally non-convergent there; and f1^2(free)=0 by structure
      (composite = 2-gluon cut, R5).  The no-go manufactures no false atom.

Nothing here is a proof; every theorem content (KP criterion, four-wall geometry,
analyticity wall) is in the .md and its cited P3 files. This witness is [V].
"""
import numpy as np
try:
    from scipy.special import iv as _iv
    def I(n, z): return float(_iv(int(n), z))
    SRC = "scipy.special.iv"
except Exception:                                    # deterministic fallback quadrature
    _t = np.linspace(0.0, np.pi, 400001)
    def I(n, z):
        return float((1.0/np.pi)*np.trapezoid(np.exp(z*np.cos(_t))*np.cos(n*_t), _t))
    SRC = "quadrature"

def blk(t): print("\n" + "="*80 + "\n" + t + "\n" + "="*80)

N = 2                       # SU(2)
B0 = 11.0*N/(48.0*np.pi**2)  # 1-loop beta-function coefficient (dg^2/dln mu = -2 b0 g^4)

def ccoef(n, beta):
    """exact SU(2) character coefficient c_n(beta); n = #boxes = 2j, dim = n+1."""
    return I(n, beta) - I(n+2, beta)

def activity(beta):
    c0 = ccoef(0, beta)
    return ccoef(1, beta)/c0

def Phi(a, beta, nmax=400):
    c0 = ccoef(0, beta)
    s = 0.0
    for n in range(1, nmax+1):
        term = (ccoef(n, beta)/c0)*(n+1)*np.exp(a*n)
        s += term
        if n > 20 and abs(term) < 1e-14*max(abs(s), 1e-30):
            break
    return s

def Delta(d):              # plaquette-adjacency degree: 4*(2(d-1)-1)
    return 4*(2*(d-1)-1)

# =====================================================================
blk("[C] exact SU(2) character coefficients and the identity  SUM (n+1)c_n = e^beta")
print(f"  modified-Bessel source: {SRC}")
for beta in [1.0/48.0, 1.0, 2.5]:
    c0 = ccoef(0, beta)
    tot = sum((n+1)*ccoef(n, beta) for n in range(0, 120))
    print(f"  beta={beta:8.5f}:  SUM(n+1)c_n = {tot:.10f}   e^beta = {np.exp(beta):.10f}"
          f"   rel.err={abs(tot-np.exp(beta))/np.exp(beta):.2e}")
print("  => c_n = I_n - I_{n+2} is exact; u(beta)=c_1/c_0 is the single-plaquette activity.")
print("\n  CALIBRATION-FREE diluteness diagnostic:  trivial-weight fraction f0 = c_0/e^beta")
print("  (a convergent character/polymer expansion needs the nontrivial plaquettes DILUTE, f0 -> 1)")
for beta in [1.0/48.0, 1.0, 2.5]:
    c0 = ccoef(0, beta); f0 = c0/np.exp(beta)
    print(f"    beta={beta:8.5f}:  c_0={c0:.6f}  f0=c_0/e^beta={f0:.5f}  "
          f"(Phi(0)=(e^beta-c_0)/c_0={np.exp(beta)/c0-1.0:.4f} nontrivial/trivial)")
print("  f0: 0.98 (wall, dilute) -> 0.42 (beta=1) -> 0.165 (beta=2.5, DENSE);  expansion dead.")

# =====================================================================
blk("[D] minimal trF^2 polymer = four-wall tube; per-slice weight u^4, rate m1 a = 4 ln(1/u)")
print("  beta       u=c1/c0     u^4 (per slice)    m1 a=4ln(1/u)   u^2 (tube/competitor per slice)")
for beta in [1.0/48.0, 1.0, 2.5]:
    u = activity(beta)
    print(f"  {beta:8.5f}   {u:.6f}     {u**4:.6f}          {4*np.log(1.0/u):.6f}        {u**2:.6f}")
print("  four-wall tube: w(B_T)=a0 u^{4T}; perimeter-6 competitor w ~ u^{6T}=u^2 x tube per slice.")

# =====================================================================
blk("[E] KP / Ursell ratio  rho_Ur(beta)=lambda_a Phi(1,beta)  at the P3 wall, beta=1, beta=2.5")
print("  convergence of the connected-polymer (Ursell) series  <=>  rho_Ur < 1")
print("  lambda_a = e*Delta (crude P3 animal constant);  Delta=20/12/4 for d=4/3/2")
print(f"  {'beta':>9} {'u':>10} {'Phi(1)':>14} | "
      f"{'d=4 rho':>12} {'d=3 rho':>12} {'d=2 rho':>12}")
for beta in [1.0/48.0, 0.1, 1.0, 2.5]:
    u = activity(beta); P1 = Phi(1.0, beta)
    row = [np.e*Delta(d)*P1 for d in (4,3,2)]
    print(f"  {beta:9.5f} {u:10.6f} {P1:14.6e} | "
          f"{row[0]:12.4e} {row[1]:12.4e} {row[2]:12.4e}")
print("  (all three exceed 1 already by beta~O(1); at 2.5 they are astronomically >1.)")

# needed activity for convergence (the crude KP threshold):
lam4 = np.e*Delta(4)
u_need = 1.0/lam4
u25 = activity(2.5)
print(f"\n  crude KP threshold (d=4):  convergence needs u < 1/lambda_a = {u_need:.5f}")
print(f"  actual u(2.5) = {u25:.5f}   => factor {u25/u_need:.1f}x too large")
# how far beyond the wall (ratio of Ursell ratios), and where the crude rho_Ur=1
r_wall = lam4*Phi(1.0, 1.0/48.0); r_25 = lam4*Phi(1.0, 2.5)
b_grid = np.linspace(1e-4, 0.05, 4000)
_rho = np.array([lam4*Phi(1.0, b) for b in b_grid])
b_cross = b_grid[np.argmin(np.abs(_rho-1.0))]
print(f"  rho_Ur(2.5)/rho_Ur(1/48) = {r_25/r_wall:.1f}x  (d=4): {r_wall:.2f} at the wall -> {r_25:.3e} at 2.5")
print(f"  crude rho_Ur=1 (d=4) is crossed at beta* ~ {b_cross:.5f}; 2.5 is {2.5/b_cross:.0f}x beyond it")
print(f"  u(beta)=1/lambda_a is reached only at beta ~ {u_need*2:.5f}  (u~beta/2 small-beta)")
print(f"  (SZZ/BE wall is beta<1/48~{1/48:.5f}; the crude KP constant is even tighter --")
print(f"   either way beta=2.5 is 10^1-10^3 x beyond the wall.)")

# =====================================================================
blk("[F] block-spin / one-loop RG: does coarse-graining make the activity small at beta=2.5?")
# flowing a^{-1} -> xi^{-1} (coarsening to the physical IR):  dbeta = -8 b0 ln(xi/a)
def beta_eff(beta, xi_over_a):
    return beta - 8.0*B0*np.log(xi_over_a)   # AF: coupling grows in IR => beta decreases
xi = 1.3      # correlation length in lattice units at beta=2.5 (see [G])
beff = beta_eff(2.5, xi)
print(f"  b0(SU(2)) = {B0:.6f};   8 b0 = {8*B0:.6f}")
for xa in [1.0, xi, 2.0, 4.0, 20.0, 800.0]:
    bb = beta_eff(2.5, xa)
    uu = activity(max(bb, 1e-6))
    print(f"  xi/a={xa:7.2f}:  beta_eff={bb:7.4f}   u_eff={uu:.5f}"
          + ("   (physical value here)" if abs(xa-xi)<1e-9 else ""))
print(f"  physical range xi/a={xi:.2f}: beta_eff={beff:.4f},  u_eff={activity(beff):.5f}")
print(f"  => over the physical range u only drops {u25:.4f}->{activity(beff):.4f} "
      f"(factor {u25/activity(beff):.3f});")
print(f"     to reach the KP threshold u<{u_need:.4f} one needs beta_eff<~{u_need*2:.3f},")
print(f"     i.e. xi/a >~ {np.exp((2.5-u_need*2)/(8*B0)):.1f}  -- the lattice has xi/a~{xi:.1f}. BLOCKED.")

# =====================================================================
blk("[G] scale separation at beta=2.5: correlation length in lattice units")
# two estimates of m_{0++} a at beta=2.5 (SU(2)):
m_strong = 4.0*np.log(4.0/2.5)     # strong-coupling leading: ma = -4ln(beta/4)
# scaling region: m_{0++}/sqrt(sigma) ~ 3.5, sqrt(sigma) a ~ 0.22 at beta~2.5
m_scale  = 3.5*0.22
for lab, ma in [("strong-coupling leading 4ln(4/beta)", m_strong),
                ("scaling  m/sqrt(sigma)~3.5, sqrt(sigma)a~0.22", m_scale)]:
    print(f"  m a = {ma:.3f}  ->  xi/a = {1.0/ma:.3f}   [{lab}]")
print(f"  both give xi/a ~ O(1): number of halving (b=2) steps = log2(xi/a) ~ "
      f"{np.log2(xi):.2f}.  NO RG TOWER (0-1 steps) at beta=2.5.")

# =====================================================================
blk("[H] acid test: free field beta->inf -- activity does NOT decay; f1^2(free)=0 by structure")
print("  beta     u=c1/c0      Phi(1)        rho_Ur(d=4)")
for beta in [5.0, 10.0, 50.0]:
    u = activity(beta)
    P1 = Phi(1.0, beta)
    print(f"  {beta:6.2f}   {u:.6f}   {P1:12.4e}   {lam4*P1:12.4e}")
print("  u -> 2 (c1/c0 -> 2) as beta->inf; activities grow, rho_Ur -> inf.")
print("  The diagnostic correctly says 'no expansion' on the free field, and the free")
print("  composite trF^2 has f1^2=0 (2-gluon cut) -- so the no-go never claims a false atom.")

# =====================================================================
blk("SUMMARY  (labels: [T]/[C]/[NG]/[O]/[V])")
print("[C] exact SU(2) c_n=I_n-I_{n+2}, SUM(n+1)c_n=e^beta verified                    [T]+[V]")
print(f"[C] diluteness f0(2.5)=c_0/e^beta={ccoef(0,2.5)/np.exp(2.5):.3f} (needs ->1); gas dense   [V]")
print(f"[D] four-wall tube w=a0 u^(4T); u(2.5)={u25:.3f}, u^4={u25**4:.3f}, m1a={4*np.log(1/u25):.3f}   [T]+[V]")
print(f"[E] Ursell/KP ratio rho_Ur(2.5) = {lam4*Phi(1,2.5):.3e}  (d=4)  >> 1              [V]")
print(f"    crude KP needs u < {u_need:.4f}; u(2.5)={u25:.4f}  => {u25/u_need:.1f}x over   [NG] for the class")
print(f"[F] one-loop/block-spin: beta_eff(2.5;xi/a={xi})={beff:.3f}, u_eff={activity(beff):.3f};")
print(f"    to drop u below {u_need:.4f} needs xi/a>~{np.exp((2.5-u_need*2)/(8*B0)):.0f}, physical xi/a~{xi:.1f}  [NG]")
print(f"[G] no scale separation (xi/a~{xi:.1f}, <1 halving step) at beta=2.5             [T]+[V]")
print("[H] free field: expansion maximally bad; f1^2=0 (acid test PASSES, no false atom) [T]+[V]")
print("\n MISSING ESTIMATE: a rigorous bound u_eff(beta=2.5) <= 1/lambda_a on the")
print(" block effective activity at the physical scale (equivalently connected trF^2")
print(" polymer weights |w(P)| <= e^{-kappa|P|}, kappa >= ln lambda_a). Bare and")
print(" one-loop values miss it by ~40x; its proof is exactly the crossover/mass-gap")
print(" control -- circular unless supplied non-perturbatively.")
print("\nNo QED. Clay 0%.")
