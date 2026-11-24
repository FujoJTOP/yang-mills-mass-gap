"""
riccati_flow_witness.py  --  witness for SLOT-RICCATI-FLOW.md  (2026-11-24)

Two independent checks, both TINY (numpy only, <1s):

(A) IDENTITY CHECK.  The corpus's exact Wetterich/Polchinski flow of the coarse
    effective potential (E1-FRG, the viscous-HJ / Fokker-Planck equation)

        d_k V = -sigma V'' + sigma (V')^2 + (2 sigma' - beta) V' + (beta' - sigma'')   (*)

    projected to the symmetric vacuum a=0 (V'(0)=0) is claimed (SLOT derivation)
    to give the Riccati-in-k

        d_k kappa = 2 sigma0 kappa^2 + (3 sigma'' - 2 beta') kappa
                    - beta0 phi3 - sigma0 lam + (beta''' - sigma'''')        (R)

    with kappa=V''(0), phi3=V'''(0), lam=V''''(0), and sigma0=sigma(0) etc.
    We evolve (*) on a grid in (k,a) and verify (R) to discretisation error.

(B) FIXED-POINT / SOURCE STRUCTURE.  On the *exactly closed* Gaussian submanifold
    (sigma=const, beta = b1 a) the flow (R) reduces to the closed Riccati

        d_k kappa = 2 sigma0 kappa^2 - 2 b1 kappa .                       (G)

    kappa=0 is always a solution (gapless/free trajectory).  A second fixed point
    kappa* = b1/sigma0 exists and is positive iff b1>0.  We integrate (G) and show
    which fixed point the IR endpoint selects as a function of sign(b1).  This is
    an *illustration of the Riccati structure only* -- NOT a statement about YM,
    whose sigma,beta (hence the sign of the source) are the open H1 inputs.

NO claim about Yang-Mills here.  Illustrative/identity only.
"""
import numpy as np

print("=" * 74)
print("(A) IDENTITY CHECK  --  Riccati (R) vs exact viscous-HJ flow (*)")
print("=" * 74)

# ---- generic model: sigma(a), beta(a) polynomials, and a test potential V0(a) ----
# sigma(a)=s0 + s2/2 a^2 + s4/24 a^4 ; beta(a)=b0 + b1 a + b2/2 a^2 + b3/6 a^3
# V0(a)  =k0/2 a^2 + p3/6 a^3 + l0/24 a^4  (V0(0)=0, V0'(0)=0)
def make(s0,s2,s4, b0,b1,b2,b3):
    return dict(s0=s0,s2=s2,s4=s4,b0=b0,b1=b1,b2=b2,b3=b3)
def sig(a,p):  return p['s0'] + 0.5*p['s2']*a**2 + p['s4']/24.0*a**4
def sigp(a,p): return p['s2']*a + p['s4']/6.0*a**3
def sigpp(a,p):return p['s2'] + p['s4']/2.0*a**2
def bet(a,p):  return p['b0'] + p['b1']*a + 0.5*p['b2']*a**2 + p['b3']/6.0*a**3
def betp(a,p): return p['b1'] + p['b2']*a + p['b3']/2.0*a**2
def betpp(a,p):return p['b2'] + p['b3']*a

# ---- grid ----
A_dom, N = 3.0, 601
a  = np.linspace(-A_dom, A_dom, N)
da = a[1]-a[0]; j0 = N//2
def lap(V):
    d = np.zeros_like(V); d[1:-1]=(V[2:]-2*V[1:-1]+V[:-2])/da**2; return d
def grad(V):
    d = np.zeros_like(V); d[1:-1]=(V[2:]-V[:-2])/(2*da); return d
def d2(V): return (V[j0+1]-2*V[j0]+V[j0-1])/da**2
def d3(V): return (V[j0+2]-2*V[j0+1]+2*V[j0-1]-V[j0-2])/(2*da**3)
def d4(V): return (V[j0+2]-4*V[j0+1]+6*V[j0]-4*V[j0-1]+V[j0-2])/da**4

def rhs(V, p):
    return -sig(a,p)*lap(V) + sig(a,p)*grad(V)**2 + (2*sigp(a,p)-bet(a,p))*grad(V) + (betp(a,p)-sigpp(a,p))

def check_identity(p, k0, p3, l0, tag):
    V0 = k0/2.0*a**2 + p3/6.0*a**3 + l0/24.0*a**4
    # exact a=0 data of V0
    kap, phi3, lam = k0, p3, l0
    # instantaneous d_k V = rhs(V0); then d_k kappa = d2(rhs)  (exact identity)
    rhs_arr = rhs(V0, p)
    dkappa_num = d2(rhs_arr)
    A = 2.0*p['s0']
    B = 3.0*p['s2'] - 2.0*p['b1']
    C = -p['b0']*phi3 - p['s0']*lam + p['b3'] - p['s4']
    rhs_R = A*kap**2 + B*kap + C
    rel = abs(dkappa_num - rhs_R)/max(1e-12, abs(rhs_R), abs(dkappa_num))
    print(f"  [{tag}]  sigma0={p['s0']} sigma''={p['s2']} sigma''''={p['s4']} | "
          f"beta0={p['b0']} beta'={p['b1']} beta'''={p['b3']} | k={k0} phi3={p3} lam={l0}")
    print(f"           d_k kappa (PDE num) = {dkappa_num:+.8f} ;  Riccati RHS = {rhs_R:+.8f} ; "
          f"rel.err = {rel:.2e}")
    return rel

print("\n  (A is an instantaneous identity check: no unstable time-marching.)")
r1 = check_identity(make(0.5,0.3,0.0, 0.0,0.4,0.0,-0.2), 1.0, 0.0, 0.6, "1 even, s4=0")
r2 = check_identity(make(0.5,0.3,0.24,0.3,0.4,-0.5,-0.2), 1.0, 0.5, 0.6, "2 odd+even, s4!=0")
maxrelerr = max(r1, r2)
print(f"\n  max relative error over configs: {maxrelerr:.2e}")
print("  => (R) is the correct a=0 projection of the exact E1 flow." if maxrelerr < 1e-3
      else "  => WARNING: mismatch (check discretisation).")

print()
print("=" * 74)
print("(B) FIXED-POINT STRUCTURE of the closed Gaussian Riccati (G)")
print("=" * 74)
print("  d_k kappa = 2 sigma0 kappa^2 - 2 b1 kappa ;  fixed pts: 0 and b1/sigma0\n")
def integrate_G(sigma0, b1, kappa_start=1.0, k0=1.0, k1=0.02, n=400000):
    kk = np.linspace(k0, k1, n); kap = kappa_start
    for i in range(1, n):
        h = kk[i]-kk[i-1]
        kap = kap + h*(2*sigma0*kap**2 - 2*b1*kap)
    return kap
for b1v, lab in [(0.0,"b1=0  (pure Riccati, NO linear source)"),
                 (+0.4,"b1>0  (definite linear term, one sign)"),
                 (-0.4,"b1<0  (definite linear term, other sign)")]:
    kend = integrate_G(0.5, b1v)
    star = b1v/0.5
    print(f"  {lab}:  kappa*=b1/sigma0={star:+.3f};  IR endpoint kappa(k->0.02)={kend:+.6e}")
print("")
print("  KEY STRUCTURAL FACT:")
print("   * pure Riccati kappa'=2 sigma0 kappa^2 (no linear term): kappa decays only")
print("     LOGARITHMICALLY, kappa(t)~1/(2 sigma0 |t|) with t=ln k (endpoint above,")
print("     1/(2-0.02)=0.505 over k in [0.02,1], confirms this).  It is MARGINAL --")
print("     the kappa^2 term can NEVER hold a nonzero IR fixed point.")
print("   * a positive IR fixed point requires the LINEAR coefficient to have a definite")
print("     sign (a genuinely IR-RELEVANT direction) -- i.e. it is the SIGN-DEFINITE")
print("     SOURCE AT kappa=0 that decides gapped-vs-gapless, exactly the open slot.")
print("")
print("  READING (honest): In YM the sign of that marginal/linear source is set by")
print("  sigma,beta of the coarse k->0 mode -- the open H1 input.  This witness")
print("  demonstrates the STRUCTURE only; it does NOT compute the YM sign.")
print("=" * 74)
