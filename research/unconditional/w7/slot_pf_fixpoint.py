#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
slot_pf_fixpoint.py  --  witness for  research/unconditional/SLOT-PF-FIXPOINT.md

Can the gauge-invariant Perron-Frobenius (strictly positive transfer kernel)
BE a k->0 sign-definite self-consistency that forces a gap rho = lambda1/lambda0 < 1,
and does it EXCLUDE the gapless fixed point rho = 1 ?

Single job, pure numpy/scipy (+ mpmath reused from finite_vol_certify). No torch/cvxpy.

Sections
  A. PF / Jentzsch on strictly-positive toy kernels:
       lambda0 simple, Omega>0, lambda1<lambda0  -- strict separation exists.
       Two-block kernel with cross coupling delta: rho->1 as delta->0 WHILE the
       kernel stays strictly positive (no uniform modulus from positivity alone).
  B. Sign-blindness of the positive cone (the central lemma, numerically):
       the deflated operator T_perp = T - lambda0 Omega Omega^T is sign-indefinite
       (min entry < 0); its spectral radius is exactly lambda1 but the lambda1
       eigenvector psi1 is sign-changing and the Collatz-Wielandt functional on
       T_perp FAILS to reach lambda1 (no sign-definite lower bracket).
  C. FREE MASSLESS FIELD acid test:
       T = exp(-a*M), M = -Delta_lattice on the L^3 torus -> manifestly strictly
       positive heat kernel (PF applies exactly), lambda0=1, lambda1=e^{-a mu1},
       rho = e^{-a mu1} -> 1 as L->inf.  Strictly positive kernel + PF +
       rho->1  =>  positivity does NOT exclude the gapless limit.
  D. Cheeger / Poincare sign-definite LOWER face:
       1-rho >= Phi^2/2 (Cheeger).  Phi_free bracketed by gap -> 0.  The only
       sign-definite lower mechanism is the conductance = Poincare constant = H1.
  E. SU(2) single-link at beta=2.5: Peter-Weyl ratio lambda_{1/2}/lambda0 =
       I2(beta_eff)/I1(beta_eff) (slow, close to 1).  The coarse-mode gap is not
       controlled by any local positivity datum (uniformity wall).

HONEST SCOPE.  Nothing here proves a uniform 4D gap; it pins the EXACT step where
the PF self-consistency fails (sign-blindness of the single cone + non-uniform
finite-V separation).  Labels in the .md: [T]/[C]/[O]/[NG]/[V-here].  No QED.
"""
import numpy as np
from scipy.linalg import expm
from scipy.special import iv

np.set_printoptions(precision=6, suppress=True)
SEP = "=" * 80
def sec(s):
    print("\n" + SEP + "\n" + s + "\n" + SEP)


# ======================================================================
# A.  PF / Jentzsch on strictly-positive toy kernels
# ======================================================================
sec("A. PF / Jentzsch: lambda0 simple, Omega>0, lambda1<lambda0  (strict separation)")

def pf_report(K, name):
    """Symmetric matrix with strictly positive entries: PF check."""
    assert np.all(K > 0), "kernel not strictly positive"
    lam, Q = np.linalg.eigh(K)               # ascending
    lam = lam[::-1]; Q = Q[:, ::-1]          # descending
    lam0, lam1 = lam[0], lam[1]
    Om = Q[:, 0] * np.sign(Q[0, 0])          # Perron vector, fix sign
    # multiplicity of lam0 (simple?)
    simple = (lam0 - lam1) > 1e-10 * abs(lam0)
    print(f"  {name}: n={K.shape[0]}, lambda0={lam0:.6f}, lambda1={lam1:.6f}, "
          f"rho=lambda1/lambda0={lam1/lam0:.6f}")
    print(f"      lambda0 simple = {simple};  Omega>0 everywhere = {np.all(Om > 0)}"
          f" (min Omega={Om.min():.3e});  min(K)={K.min():.3e}>0")
    return lam, Q, Om

# Toy 1: Gaussian (strongly gapped)
n = 8; x = np.arange(n)
K1 = np.exp(-(x[:, None] - x[None, :])**2 / (2 * 0.7**2))
pf_report(K1, "Gaussian kernel (bandwidth .7)")

# Toy 2: two dense blocks weakly coupled by delta>0  (strictly positive, rho->1)
def two_block(nb, delta):
    N = 2 * nb; K = np.full((N, N), delta)
    K[:nb, :nb] = 1.0; K[nb:, nb:] = 1.0
    return K
print("\n  Two-block kernel (cross coupling delta>0; strictly positive for all delta):")
print(f"    {'delta':>10} {'lambda0':>10} {'lambda1':>10} {'rho':>10} {'1-rho':>10} {'min(K)':>9}")
for delta in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6]:
    K = two_block(4, delta)
    lam = np.linalg.eigvalsh(K)[::-1]
    print(f"    {delta:>10.0e} {lam[0]:>10.5f} {lam[1]:>10.5f} "
          f"{lam[1]/lam[0]:>10.6f} {1-lam[1]/lam[0]:>10.2e} {K.min():>9.0e}")
print("  => strict positivity (min K=delta>0) gives lambda1<lambda0 at EVERY delta,")
print("     but 1-rho -> 0 as delta->0: positivity alone supplies NO uniform modulus.")

lam, Q, Om = pf_report(K1, "Gaussian (again, for B)")
lam0, lam1 = lam[0], lam[1]


# ======================================================================
# B.  Sign-blindness of the positive cone  (the central lemma)
# ======================================================================
sec("B. Deflation kills positivity: T_perp = T - lambda0 Omega Omega^T is sign-indefinite")
print("    and Collatz-Wielandt on T_perp cannot see lambda1 (no sign-definite lower bound).")

T = K1.copy()
Tperp = T - lam0 * np.outer(Om, Om)
print(f"  T strictly positive: min(T)={T.min():.3e}>0")
print(f"  T_perp = T - lambda0 Omega Omega^T :  min(T_perp)={Tperp.min():.3e}  "
      f"max(T_perp)={Tperp.max():.3e}")
print(f"  negative entries in T_perp: {np.sum(Tperp < 0)} / {Tperp.size}  "
      f"-> kernel is SIGN-INDEFINITE")

# eigenvalues of T_perp = {0} U {lambda1, lambda2, ...}
evp = np.linalg.eigvalsh(Tperp)[::-1]
print(f"  eig(T_perp) top three = {evp[:3]};  spectral radius = {evp[0]:.6f}  "
      f"(=|lambda1(T)|={lam1:.6f}: MATCH)")

# lambda1 eigenvector of T is sign-changing
psi1 = Q[:, 1]
print(f"  psi1 (lambda1-eigenvector of T): sign-changing = "
      f"{np.any(psi1 > 1e-9) and np.any(psi1 < -1e-9)}  "
      f"(min={psi1.min():.3f}, max={psi1.max():.3f});  <psi1,Omega>={psi1@Om:.2e}")

# Collatz-Wielandt functionals on a positive operator vs the sign-indefinite T_perp.
# For a NONNEGATIVE matrix A: rho(A) = max_{x>0} min_i (Ax)_i/x_i.
# Apply the SAME functional to T_perp: does it reach lambda1? (It must not.)
def cw_lower(A, rng, trials=4000, seed_vecs=()):
    """sup over positive x of  min_i (Ax)_i / x_i   (Collatz-Wielandt lower).
       seed_vecs: extra positive test vectors evaluated exactly (e.g. Perron)."""
    best = -np.inf
    for xv in seed_vecs:
        r = (A @ xv) / xv
        best = max(best, r.min())
    for _ in range(trials):
        xv = np.exp(rng.normal(0, 1.5, A.shape[0]))     # positive test vector
        r = (A @ xv) / xv
        val = r.min()
        if val > best:
            best = val
    return best

rng = np.random.default_rng(0)
# include the Perron vector (gives exactly lambda0 on T) as an exact trial
cw_T    = cw_lower(T,     rng, seed_vecs=(Om,))
cw_Tperp = cw_lower(Tperp, rng, seed_vecs=(Om,))
print(f"\n  Collatz-Wielandt lower functional  max_(x>0) min_i (A x)_i / x_i :")
print(f"    on T      (nonnegative): {cw_T:.6f}   ~ lambda0={lam0:.6f}  (works: CW sees the TOP)")
print(f"    on T_perp (sign-indef) : {cw_Tperp:.6f}   vs true spectral radius lambda1={lam1:.6f}")
print("  => CW on T_perp gives a value BELOW 0 (or far from lambda1): the sign-definite")
print("     lower bracket COLLAPSES. The single positive cone is blind to the 2nd level.")


# ======================================================================
# C.  Free massless field acid test: strictly-positive heat kernel, rho -> 1
# ======================================================================
sec("C. FREE MASSLESS FIELD: T=exp(-a M), M=-Delta_lat on L^3 torus (strictly positive)")

def neg_lap_position(L):
    """Laplacian M = D - A (integer, periodic, 7-point) on the L^3 torus, position space."""
    N = L ** 3
    def idx(i, j, k): return (i % L) * L * L + (j % L) * L + (k % L)
    M = np.zeros((N, N), dtype=float)
    for i in range(L):
        for j in range(L):
            for k in range(L):
                a = idx(i, j, k); M[a, a] += 6
                for di, dj, dk in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
                    M[a, idx(i+di, j+dj, k+dk)] -= 1
    return M

a_step = 1.0
print(f"  a (Euclidean step) = {a_step};  box gap (energy) = 2 sin(pi/L) [matches FINITE-VOL-CERTIFY]")
print(f"  {'L':>3} {'mu1':>10} {'2sin(pi/L)':>12} {'rho=e^{-a mu1}':>16} {'1-rho':>10} "
      f"{'min(T)>0?':>10}")
for L in [2, 3, 4, 6, 8]:
    M = neg_lap_position(L)
    mu = np.linalg.eigvalsh(M)                    # ascending, mu0=0
    mu1 = mu[1]
    Tf = expm(-a_step * M)                         # heat kernel: strictly positive
    minT = Tf.min()
    rho = np.exp(-a_step * mu1)
    print(f"  {L:>3} {mu1:>10.5f} {2*np.sin(np.pi/L):>12.5f} {rho:>16.6f} "
          f"{1-rho:>10.2e} {str(minT > 0):>10}")
    if L == 4:
        # explicit PF confirmation of the heat-kernel matrix at L=4.
        # The Perron vector of the heat kernel is the CONSTANT mode (M.1=0):
        # mathematically strictly positive; eigh returns it up to 1e-15 roundoff.
        lamf, Qf = np.linalg.eigh(Tf); lamf = lamf[::-1]
        Omf = Qf[:, -1] * np.sign(Qf[0, -1])   # eigenvector of the LARGEST eigenvalue
        const_ok = np.allclose(Omf, Omf[0]) and Omf[0] > 0
        print(f"      [L=4 explicit] lambda0={lamf[0]:.6f} simple, "
              f"lambda1={lamf[1]:.6f}, rho={lamf[1]/lamf[0]:.6f}, "
              f"Omega=const>0={const_ok} (min Omega={Omf.min():.2e}), min(T)={minT:.2e}")
print("  => the FREE MASSLESS transfer is a STRICTLY POSITIVE (PF) kernel with lambda1<lambda0")
print("     at every finite L, yet rho -> 1 as L -> inf: strict positivity does NOT exclude the")
print("     gapless limit. (Acid test: the PF mechanism manufactures no uniform mass.)")


# ======================================================================
# D.  Cheeger / Poincare sign-definite LOWER face:  1 - rho >= Phi^2/2
# ======================================================================
sec("D. Cheeger sign-definite lower bound: 1-rho >= Phi^2/2,  Phi_free bracketed -> 0")
print("  For a reversible chain gap g=1-lambda1, Cheeger's inequality gives")
print("     g <= 2*Phi  and  g >= Phi^2/2   =>   g/2 <= Phi <= sqrt(2 g).")
print(f"  {'L':>3} {'g=mu1':>10} {'Phi_lower=g/2':>14} {'Phi_upper=sqrt(2g)':>20}")
for L in [2, 3, 4, 6, 8, 16]:
    M = neg_lap_position(L); mu = np.linalg.eigvalsh(M)
    g = mu[1]                                      # continuous-time gap (=mu1, mu0=0)
    print(f"  {L:>3} {g:>10.5f} {g/2:>14.5f} {np.sqrt(2*g):>20.5f}")
print("  => Phi_free -> 0 (both brackets): the free field has NO uniform isoperimetry.")
print("     The ONLY sign-definite LOWER mechanism (Cheeger conductance = Poincare constant)")
print("     therefore yields nothing uniformly on the free field -- and closing it on YM")
print("     (Phi>=delta uniform at beta=2.5) is EXACTLY the coarse Poincare = H1")


# ======================================================================
# E.  SU(2) single-link at beta=2.5: local Peter-Weyl ratio (slow)
# ======================================================================
sec("E. SU(2) single-link at beta=2.5: Peter-Weyl ratio lambda_{1/2}/lambda0 (uniformity wall)")
beta = 2.5
# Single plaquette RP transfer (=2D SU(2) YM = heat kernel): leading ratio I2(beta)/I1(beta)
r_plaq = iv(2, beta) / iv(1, beta)
# Single link, aligned 6-staple background: exponent beta_eff * u0, beta_eff = 6*beta = 15
beta_eff = 6.0 * beta
r_link = iv(2, beta_eff) / iv(1, beta_eff)
print(f"  beta = {beta} (SU(2), Wilson).")
print(f"  single-plaquette (2D) ratio I2/I1         = {r_plaq:.6f}  (local gap 1-ratio = {1-r_plaq:.6f})")
print(f"  single-link (beta_eff=6beta={beta_eff:.0f}) I2/I1 = {r_link:.6f} (local gap = {1-r_link:.6f})")
print("  As beta_eff grows the local ratio -> 1 (slow): the local positivity datum is")
print("  NOT a modulus for the COARSE mode (k->0); block positivity is silent at k=0")
print("  (PROOF-COARSE-POINCARE).  No uniform deflated gap at beta=2.5 is obtainable from")
print("  local positivity; FINITE-VOL-CERTIFY: dense certification infeasible at L>=2.")


# ======================================================================
# F.  VERDICT
# ======================================================================
sec("F. VERDICT")
print("""  [T]  Strict positivity => lambda0 simple, Omega>0, lambda1<lambda0  (Jentzsch).
       => rho=1 EXCLUDED at each finite (V,a): a finite-volume gap EXISTS.
  [T]  Sign-blindness (central): every eigenvector psi (lambda<lambda0) is Omega-orthogonal
       with Omega>0 => psi SIGN-CHANGING => outside the positive cone.  Hence Collatz-
       Wielandt / Hilbert-metric / Birkhoff CANNOT see lambda1; the deflated kernel
       T_perp is sign-indefinite (Sec. B).  The single PF cone is structurally blind to the gap.
  [V]  Free massless field: strictly positive PF kernel, rho->1 as L->inf (Sec. C);
       Phi_free->0 (Sec. D).  => rho=1 IS admissible as a UNIFORM limit with strict positivity.
  [NG] Strict positivity forces separation only PER-(V,a), with NO uniform modulus
       (Sec. A two-block: 1-rho->0 while min K>0).
  [O]  The only sign-definite LOWER mechanism (Cheeger/Poincare conductance) == the
       coarse Poincare == H1 (Sec. D).  Uniform exclusion of rho=1 = H1.
  ------------------------------------------------------------------------------
  Is rho=1 a solution?   finite (V,a): NO (Jentzsch).  uniformly: YES (free field).
  Does strict positivity force separation?  YES per-(V,a) [T]; NO uniformly [NG].
  Verdict: the PF SELF-CONSISTENCY does not close -- exact failing step pinned:
           sign-blindness of the single cone + non-uniform finite-V separation;
           the uniform statement reduces to the Poincare/Cheeger constant = H1.
  No QED.  Clay open.
""")
