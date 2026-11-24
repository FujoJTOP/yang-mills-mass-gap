#!/usr/bin/env python3
r"""
q1_infrared_bound.py  --  Q1: can a reflection-positivity (RP) infrared bound reach beta=2.5?

Deliverable witness for research/v3/Q1-INFRARED-BOUND.md.

The H1 object is the zero-momentum (k=0) susceptibility
    chi = Var_mu(a) ,   a = (1/sqrt(N_links)) * sum_links xi_link ,  xi = log U (flat chart),
equivalently G''_{k=0} = 1/chi (Gaussian marginal).  H1 = "chi <= 1/lambda uniformly in (L,a)."

We test the four standard RP/compactness handles and show exactly where each dies.

 [A] FREE NONCOMPACT MASSLESS FIELD (the "free photon" acid object, d=4).
     The k=0 mode is the zero mode of -Delta: chi = 1/m^2 -> infinity as m->0, INDEPENDENT of L.
     This is the object FSS must NOT bound (acid test).  Also show the *other* object FSS does
     bound -- the average of the NON-constant modes -- is finite; the divergence is purely k=0.

 [B] COMPACT beta=0 SU(2): chi = Var_Haar(xi_3) = (4 pi^2 - 6)/9 = 3.7198, EXACT, L-independent.
     Compactness makes beta=0 finite (would be infinite for the noncompact photon).  [T]

 [C] WHY COMPACTNESS ALONE IS NOT ENOUGH (bounded-range vs bounded-variance).
     For bounded X_i in [-pi,pi], the *normalized* sum has |a| <= pi sqrt(N): boundedness gives
     only Var(a) <= pi^2 N (grows).  Independent Efron-Stein gives the finite sigma^2; a perfectly
     correlated model gives N sigma^2.  Bounded range <=> Var = O(V): the volume must be beaten by
     DECORRELATION, which is the gap, not compactness.

 [D] BOUNDED-DENSITY COMPARISON mu = e^{-S_W} Haar / Z, 0 <= S_W <= 2 beta N_p.
     Var_mu(a) <= ||d mu/d Haar||_inf Var_Haar(a) with inflation <= e^{2 beta N_p}.  For d=4,
     N_p = 6 L^4, so 2 beta N_p = 12 beta L^4 = Theta(L^4): volume-growing, not uniform.

 [E] DOBRUSHIN / INFLUENCE threshold (ties to S6).  alpha(beta) = 72 beta; alpha<1 iff
     beta < 1/72.  alpha(2.5)=180 >> 1.  The only volume-uniform variance route is strong coupling.

 [F] 1D COMPACT SU(2): the transfer operator is central (depends on U^dag U'), so by Peter-Weyl
     its eigenvalues are the SU(2) Fourier coefficients  fhat_j.  We compute fhat_0, fhat_{1/2}
     by exact quadrature, get r = fhat_{1/2}/fhat_0 < 1, gap Delta = -log r > 0, and
     chi_1D = (1+r)/(1-r) < infinity.  This is the compactness theorem that DOES work -- in 1D.
     It does NOT lift to 4D (the 3D-layer transfer operator is not compact; massless phases exist).

Run: py -3.12 q1_infrared_bound.py [quick|full].  numpy+scipy only, one job. Small sizes.
No third-party sources.  Labels [T]/[V]/[C].  No QED.  Clay 0%.
"""
import sys
import numpy as np
from scipy.special import iv

np.seterr(all='ignore')
MODE = sys.argv[1] if len(sys.argv) > 1 else 'quick'

def sec(t): print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)

# =====================================================================  [A] free noncompact
def part_A():
    sec("[A] FREE NONCOMPACT MASSLESS FIELD d=4: the k=0 susceptibility chi = 1/m^2 -> inf")
    print("  The H1 object chi = Var(a), a = (1/sqrt V) sum_x xi_x.  For a massive free field")
    print("  with covariance 1/(m^2+lambda_k), the k=0 (=constant) mode gives Var(a) = 1/m^2")
    print("  EXACTLY, for every L.  As m->0 it diverges.  (The free photon's constant gauge mode")
    print("  is the same flat direction: no restoring force.)  ACID OBJECT: must NOT be bounded.")
    print()
    print("     m        chi = 1/m^2")
    for m in [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]:
        print(f"   {m:6.3f}   {1.0/m**2:12.2f}")
    # FSS *does* bound the non-constant average: (1/V) sum_{k!=0} 1/lambda_k on the L^4 torus.
    # show it is finite and volume-bounded (this is the OTHER object, the field variance).
    print("\n  Contrast: the mean over NON-constant modes, (1/#k!=0) sum_{k!=0} 1/lambda_k,")
    print("  is the object FSS actually bounds.  On the L^4 torus it is finite and L-bounded:")
    for L in [4, 6, 8]:
        ks = 2*np.pi*np.arange(L)/L
        # lambda(k) = sum_mu 2(1-cos k_mu)
        lam = np.zeros((L, L, L, L))
        c = 2*(1 - np.cos(ks))
        for mu in range(4):
            shp = [1, 1, 1, 1]; shp[mu] = L
            lam = lam + c.reshape(shp)
        nz = lam > 1e-12
        avg = (1.0/np.sum(nz)) * np.sum(1.0/lam[nz])
        print(f"     L={L:2d}:  (1/#k+) sum 1/lambda_k = {avg:.4f}   (finite, volume-bounded)")
    print("  => FSS bounds the NON-constant average (finite); the k=0 mode it cannot touch.  [T]")

# =====================================================================  [B] compact beta=0
def part_B():
    sec("[B] COMPACT beta=0 SU(2): chi = Var_Haar(xi_3), EXACT and L-independent")
    # exact closed form for the exponential chart: Var_Haar(xi_3) = (4 pi^2 - 6)/9
    chi = (4*np.pi**2 - 6)/9.0
    print(f"  Var_Haar(xi_3) = (4 pi^2 - 6)/9 = {chi:.6f}   (L-independent, all beta=0 volumes)")
    print("  => G''_{k=0}(beta=0) -> 1/chi = %.6f  (the exact Haar floor, P1/D6)." % (1.0/chi))
    # exact quadrature cross-check: one-Cartesian-component Haar marginal of U=exp(i xi.sigma/2)
    # rho1(z) = pi ( Cin(2pi) - Cin(|z|) ),  |z|<=2pi  (from p1 witness)
    from scipy.special import sici
    EULER = 0.5772156649015328606
    def Cin(x):
        _, ci = sici(x)
        return EULER + np.log(np.where(x > 0, x, 1.0)) - ci
    z = np.linspace(-2*np.pi, 2*np.pi, 400001)
    w = np.full(z.size, 4*np.pi/(z.size-1)); w[[0, -1]] *= 0.5
    r1 = np.pi*(Cin(2*np.pi) - Cin(np.abs(z)))
    r1 = np.where(np.abs(z) <= 2*np.pi, r1, 0.0)
    dens = r1*w; Zc = dens.sum(); m1 = (z*dens).sum()/Zc
    v1 = ((z-m1)**2*dens).sum()/Zc
    print(f"  quadrature cross-check: Var = {v1:.6f}  (matches closed form to {abs(v1-chi):.1e})")
    print("  This finiteness is PURE COMPACTNESS: would be +inf for the noncompact photon ([A]).  [T]+[V]")

# =====================================================================  [C] bounded vs variance
def part_C():
    sec("[C] COMPACTNESS ALONE IS NOT ENOUGH: bounded range gives Var(a) = O(V), not O(1)")
    rng = np.random.default_rng(12345)
    print("  a = (1/sqrt N) sum X_i, X_i in [-pi,pi].  Bounded range => |a| <= pi sqrt(N).")
    print("  Independent Efron-Stein (product measure) gives the finite Var = sigma^2;")
    print("  perfect correlation gives N*sigma^2.  Same bounded range, opposite volume scaling:")
    print()
    print("     N      independent a       perfectly-correlated a     trivial bound pi^2*N")
    for N in [16, 64, 256, 1024]:
        X = rng.uniform(-np.pi, np.pi, size=(20000, N))
        a_ind = (X.sum(axis=1))/np.sqrt(N)
        # perfectly correlated: one variable repeated
        Y = rng.uniform(-np.pi, np.pi, size=(20000, 1))
        a_cor = (np.sqrt(N)*Y[:, 0])
        print(f"   {N:5d}   {np.var(a_ind):12.4f}      {np.var(a_cor):14.4f}         {np.pi**2*N:12.1f}")
    print("  => boundedness (compactness) caps Var at pi^2 N (grows).  Volume-uniform Var(a)=O(1)")
    print("     needs DECORRELATION of the sum -- i.e. the correlation length / gap.  [T]+[V]")

# =====================================================================  [D] bounded-density inflation
def part_D():
    sec("[D] BOUNDED-DENSITY COMPARISON: mu = e^{-S_W} Haar / Z  inflates Var by e^{2 beta N_p}")
    print("  0 <= S_W <= 2 beta N_p  (per plaquette 1 - (1/N)Re tr U_p in [0,2]); d=4, N_p = 6 L^4.")
    print("  Var_mu(a) <= ||d mu / d Haar||_inf * Var_Haar(a),  inflation <= e^{2 beta N_p}.")
    print("  At beta=2.5: exponent 2 beta N_p = 12 beta L^4 = 30 L^4.  Volume-growing:")
    print()
    print("     L     N_p=6L^4    2*beta*N_p(beta=2.5)     log10(inflation)")
    beta = 2.5
    for L in [2, 3, 4, 6, 8]:
        Np = 6*L**4
        expo = 2*beta*Np
        print(f"   {L:3d}   {Np:10d}      {expo:16.1f}        {expo/np.log(10):12.2f}")
    print("  => the Haar reference is finite, but the interaction reweighting is e^{Theta(L^4)}:")
    print("     NO volume-uniform transfer from the beta=0 bound.  [T]")

# =====================================================================  [E] Dobrushin
def part_E():
    sec("[E] DOBRUSHIN / INFLUENCE: the only volume-uniform variance route is strong coupling")
    kap0 = 4.0; deg = 6*(4-1)  # 6(d-1) neighbours in d=4
    coef = deg*kap0  # alpha(beta) = 18*... see S6: 6(d-1)*kap0 = 72
    coef = 6*(4-1)*kap0  # = 72
    print(f"  Dobrushin coefficient alpha(beta) = {coef:.0f} * beta ; alpha<1  <=>  beta < 1/{coef:.0f}")
    for b in [1/72, 0.19, 1.0, 2.5]:
        print(f"     beta = {b:7.4f}:  alpha = {coef*b:10.2f}   {'<1 (mixing)' if coef*b < 1 else '>>1 (no contraction)'}")
    print(f"  head beta=2.5: alpha = {coef*2.5:.0f}  (the influence/Dobrushin route is dead).  [T] (S6)")

# =====================================================================  [F] 1D compact SU(2)
def part_F():
    sec("[F] 1D COMPACT SU(2): central transfer operator => spectral gap => FINITE chi (works in 1D)")
    # transfer kernel T(U,U') = exp(beta * (1/2) tr(U^dag U')) ; central => Fourier coeffs fhat_j.
    # SU(2) class angle theta in [0,pi], Haar dU = (2/pi) sin^2 theta dtheta,
    # (1/2) tr U = cos theta,  chi_j(theta) = sin((2j+1)theta)/sin theta, dim = 2j+1.
    from numpy.polynomial.legendre import leggauss
    def fhat(beta, j):
        n = 600
        x, w = leggauss(n)               # on [-1,1]
        th = 0.5*np.pi*(x + 1.0)         # theta in [0,pi]
        w = w*0.5*np.pi
        sinth = np.sin(th)
        chi = np.where(sinth > 1e-12, np.sin((2*j+1)*th)/np.where(sinth>1e-12,sinth,1.0), 2*j+1)
        integrand = (2.0/np.pi)*sinth**2 * np.exp(beta*np.cos(th)) * chi
        return np.sum(integrand*w)/(2*j+1)
    print("  r = fhat_{1/2}/fhat_0 (fundamental decay ratio); gap Delta=-log r; chi=(1+r)/(1-r).")
    print()
    print("    beta     fhat_0       fhat_1/2      r         Delta=-log r    chi_1D")
    for beta in [0.5, 1.0, 2.5, 4.0]:
        f0 = fhat(beta, 0.0)
        fh = fhat(beta, 0.5)
        r = abs(fh/f0)
        Delta = -np.log(r) if r > 0 else np.inf
        chi = (1+r)/(1-r)
        print(f"   {beta:5.2f}   {f0:9.4f}   {fh:9.4f}   {r:8.5f}   {Delta:9.5f}   {chi:10.3f}")
    print("  => r<1 for every beta>0: the compact 1D transfer operator is gapped, chi_1D finite.  [T]+[V]")
    print("     The gap is a 1D COMPACTNESS statement.  It does NOT lift to 4D: the 3D-layer")
    print("     transfer operator is not Hilbert-Schmidt / compact, and massless phases (free photon)")
    print("     exist in d=4.  That lift is exactly H1.  [T] (structural)")

def main():
    print("Q1-INFRARED-BOUND witness  (mode=%s)" % MODE)
    part_A(); part_B(); part_C(); part_D(); part_E(); part_F()
    print("\n" + "=" * 78)
    print("VERDICT.  Compactness makes beta=0 finite (chi=3.7198, [B]) and the 1D chain gapped")
    print("([F]); RP/FSS bounds only the non-constant average ('free covariance' at k=0 is +inf, [A]);")
    print("the interaction reweighting is e^{Theta(L^4)} ([D]) and the influence route is")
    print("strong-coupling only ([E]).  NO RP bound reaches beta=2.5.  The missing lemma is a")
    print("volume-uniform finite-mass Gaussian-domination / non-Abelian FKG inequality at the")
    print("constant mode = H1.  No QED.  Clay 0%.")
    print("=" * 78)

if __name__ == "__main__":
    main()
