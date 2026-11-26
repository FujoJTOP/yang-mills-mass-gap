# D07-G : is (G)  ∫_{beta1}^{2.5} A/Delta dbeta < inf  a LIVE route?
# ------------------------------------------------------------------
# (G)'s numerator uses the SHARP form of D04-ODE (1.4):
#     sum_{n>=1} |<n|S_W|G>|^2 = Var(S_W)                [T]  (only n>=1: mean drops)
#     A^2 = Var_O(Omega) * Var(S_W)                      [T]  (Cauchy-Schwarz, sharp)
# and the EXACT beta-mass identity (FTC + Var(S_W) = -d/dbeta <S_W>):
#     int_a^b Var(S_W) dbeta = #plaq * [ <s_p>(a) - <s_p>(b) ]      [T]
# ------------------------------------------------------------------
# Four models: free Gaussian / 2D YM SU(2) / M6 (no beta!) / single plaquette (no transfer)
import numpy as np
from scipy.integrate import quad

L2 = np.log(2.0)
def out(s): print(s, flush=True)

# ---------- exact SU(2) single-plaquette quadrature: dmu = sin^2 t e^{b cos t}, t in [0,pi]
def _Z(b):  return quad(lambda t: np.sin(t)**2*np.exp(b*np.cos(t)), 0, np.pi, epsabs=1e-14, limit=400)[0]
def mom(b, k):
    Z = _Z(b)
    return quad(lambda t: np.cos(t)**k*np.sin(t)**2*np.exp(b*np.cos(t)), 0, np.pi, epsabs=1e-14, limit=400)[0]/Z
def u_pl(b):    return mom(b, 1)              # <(1/2)Re tr U_p> = <cos t>
def var_s(b):   return mom(b, 2) - mom(b, 1)**2   # s=1-cos t  =>  Var(s)=Var(cos t)

out("="*78)
out("[0] exact SU(2) single plaquette (Haar+Gibbs), s = 1 - (1/2)Re tr U")
out(f"    u_p(2.5) = {u_pl(2.5):.9f}   (corpus anchor v06_radius.log:43 = 0.507195)")
out(f"    Var(s)(2.5) = {var_s(2.5):.9f}   Var(s)(0) = {var_s(0.0):.9f}  (=1/4 Haar, PW: <|tr U/2|^2>=1/4)")
Ivar = quad(lambda b: var_s(b), 0.0, 2.5, epsabs=1e-12, limit=400)[0]
out(f"    int_0^2.5 Var(s) dbeta = {Ivar:.9f}   vs  u_p(2.5)-u_p(0) = {u_pl(2.5)-u_pl(0.0):.9f}")
out(f"    => identity int Var = #plaq * [<s>(a)-<s>(b)]  VERIFIED (residual {abs(Ivar-(u_pl(2.5)-u_pl(0.0))):.2e})")
out("="*78)

B1, B2 = 0.0034, 2.5          # Z01: beta1 ~ 3.4e-3 ; target point 2.5
T = B2 - B1
# Var_O(Omega) = 1 placeholder (local O); A = sqrt(Var_O * Var S_W)
VO = 1.0
K_NEED = np.log(64.0/1e-6)    # D04 (4.1): needed exponent for L-uniform floor eps=1e-6  (r0=1/64)
out(f"[1] Grönwall budget: need int A/Delta <= log(r0/eps) = {K_NEED:.3f}   (eps=1e-6, r0=1/64)")
out(f"    (eps=1e-3 -> {np.log(64.0/1e-3):.3f} ; eps=1e-12 -> {np.log(64.0/1e-12):.3f})")

# ---------- M1 free Gaussian field: Var(S_W)=dim/2 exactly (dim=3V polarizations), Delta=2pi/L
def free_rows(L):
    V = L**4; var = 3.0*V/2.0
    A = np.sqrt(VO*var); Delta = 2*np.pi/L
    return A, Delta, T*A/Delta
# ---------- M2 2D YM SU(2): V=L^2 plaquettes, INDEPENDENT (Migdal) => Var(S_W)=V*Var(s) exact
#            exact gap m = 3g^2/8 with beta = 2N/g^2 = 4/g^2  =>  Delta(beta)=3/(2 beta)
def ym2d_G(L):
    V = L**2
    f = lambda b: np.sqrt(VO*V*var_s(b))/(3.0/(2.0*b))
    return quad(f, B1, B2, epsabs=1e-13, limit=400)[0]
# ---------- M4 4D SU(2): exact beta-mass int Var = 6V u_p(2.5) -> Jensen lower bound on int sqrt(Var)
def ym4d_rows(L):
    V = L**4
    mass = 6.0*V*(u_pl(2.5)-u_pl(B1))          # = int Var(S_W) dbeta   [T]
    int_sqrt_lb = np.sqrt(T*mass)               # Jensen (sqrt concave)  [T]
    A_lb = np.sqrt(VO)*int_sqrt_lb              # lower bound on int A dbeta
    triv = 12.0*V*T                             # D04 (1.4) trivial A=12V|O^dag Omega|
    return int_sqrt_lb, A_lb, triv

out("")
out("[2] FOUR-MODEL TABLE for  G = int A/Delta dbeta     (A = sqrt(Var_O * Var S_W))")
out(f"{'model':<26}{'L':>4}{'Delta':>12}{'Var(S_W)':>14}{'G':>14}   verdict")
rows=[]
for L in (4,6,8,16):
    A,D,G = free_rows(L); rows.append(("free Gaussian",L,D,3*L**4/2,G))
for L in (4,6,8,16):
    G = ym2d_G(L); rows.append(("2D YM SU(2)",L,3.0/(2*B2),L**2*var_s(B2),G))
for name,L,D,Var,G in rows:
    out(f"{name:<26}{L:>4}{D:>12.4g}{Var:>14.4g}{G:>14.5g}   "
        + ("DIVERGES" if G > 10 else "ok") + ("  <= free field: f1^2=0 (must fail)" if name.startswith('free') else "  <= target TRUE here [T]"))
out(f"{'M6 (osc ⊗ free)':<26}{'-':>4}{'undef':>12}{'undef':>14}{'undef':>14}   no beta at all => (G) unformable")
out(f"{'single plaquette':<26}{'1':>4}{'undef':>12}{var_s(B2):>14.4g}{'undef':>14}   no transfer direction => Delta undefined")
out("")
out("[3] 4D SU(2): what Delta would (G) DEMAND?   (int A dbeta from the exact identity + Jensen)")
out(f"{'L':>4}{'int sqrt(Var) lb':>18}{'int A lb':>12}{'Delta_req=int A/10':>20}{'vs 2D Delta=0.6':>17}")
for L in (4,6,8,16):
    s_lb, A_lb, triv = ym4d_rows(L)
    out(f"{L:>4}{s_lb:>18.5g}{A_lb:>12.5g}{A_lb/10:>20.5g}{A_lb/10/(3.0/(2*B2)):>17.4g}")
out("")
out("[4] trivial D04 bound A=12V|O^dag Omega| (no variance structure):")
for L in (4,6):
    s_lb, A_lb, triv = ym4d_rows(L)
    out(f"    L={L}: int A dbeta <= {triv:.4g} -> even with Delta=1 the exponent >= {triv:.4g}"
        f"  (budget {K_NEED:.2f}, ratio {triv/K_NEED:.3g}x)")
out("")
out("[5] required-gap scaling: (G) with K=10 needs sup Delta >= int A dbeta / 10  ~ 0.28 L^2")
out("    every model's Delta is L-BOUNDED (free 2pi/L ->0 ; 2D 3/(2beta) exact L-indep ; 4D physical O(1))")
out("    => (G) is NOT implied by a uniform gap: uniform gap => int A/Delta = int A/Delta ~ sqrt(V) -> inf")
