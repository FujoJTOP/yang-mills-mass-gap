"""
pnt_chi_certificate.py -- witness for research/unconditional/PNT-CHI-CERTIFICATE.md

Demonstrates four things (numpy only, <2 s):
  (A)  Lemma G (exact spectral identity): chi = sum_n w_n/(1-lam_n), w_n=| <f,phi_n>|^2 >= 0,
       1-lam_n >= gamma  =>  w_2/gamma <= chi <= Var(f)/gamma .   (so chi <= 1/lambda <=> gamma >= lambda*C)
  (B)  Floor-blindness: a tiny atom at small s changes chi = int dnu/s by O(1) while perturbing
       the polynomial moments M_j = int s^j dnu by O(epsilon).  => no UV/spectral datum bounds chi.
  (C)  Certificate equivalence: chi <= 1/lambda  <=>  gamma >= lambda*C  (C in [w_2, Var]).
  (D)  Non-Gaussianity defect (Ward/Stein, exact): kappa*chi = 1-eta ; Gaussian eta=0 => chi=1/kappa;
       interacting eta>0 => chi=(1-eta)/kappa  => certificate weakens to kappa >= lambda(1-eta).
NOT a Yang-Mills calculation -- reduced models that reproduce the exact logical structure. No QED.
Run:  py -3.12 -I research/unconditional/w7/pnt_chi_certificate.py
"""
import numpy as np

def banner(s): print("\n" + "=" * 76 + "\n" + s + "\n" + "=" * 76)

# ===========================================================================
# (A) Lemma G on a Gaussian transfer operator
#   eigenvalues lam_k = e^{-m k}, k=1,2,... ; KL weights w_k >= 0
#   Var = sum w_k ; chi = sum w_k/(1-lam_k) ; gamma = 1 - lam_2 = 1 - e^{-2m}? (index k=1 -> lam_1=e^{-m})
#   NOTE: gap = 1 - lam_2 in the ordering lam_1=1 (top), lam_2=e^{-m}; here we index subdominant k=1..:
#   take lam_k = e^{-m*k} for k=1,2,... with k=1 the SUBdominant (so 1-lam_1 is the smallest denom = gamma).
# ===========================================================================
def chi_and_bounds(m, w):
    k = np.arange(1, len(w) + 1)
    lam = np.exp(-m * k)               # lam_1 is the largest subdominant eigenvalue
    denom = 1.0 - lam                  # 1-lam_k, increasing in k; smallest = 1-lam_1 = gamma
    chi = np.sum(w / denom)
    Var = np.sum(w)
    gamma = denom[0]                   # 1 - e^{-m} = gap
    lo = w[0] / gamma                  # n=2 (subdominant) term  <= chi
    hi = Var / gamma                   # sum w_n/gamma           >= chi
    return chi, Var, gamma, lo, hi

w = np.array([1.0, 0.5, 0.25, 0.125, 0.0625])   # geometric KL weights (positive)

banner("(A) Lemma G:  w_2/gamma <= chi <= Var/gamma   (exact spectral identity)")
print(f"  KL weights w = {w},  Var = {w.sum():.4f}")
for m in [0.2, 0.5, 1.0, 2.0]:
    chi, Var, gamma, lo, hi = chi_and_bounds(m, w)
    ok = (lo - 1e-12 <= chi <= hi + 1e-12)
    print(f"  m={m:4.2f}  chi={chi:9.5f}  gamma={gamma:.5f}  Var={Var:.4f}  "
          f"[lo={lo:.5f}, hi={hi:.5f}]  in-range={ok}")

banner("(A') chi ~ Var/gamma :  0 < w_2/Var <= chi*gamma/Var <= 1   (chi never easier than the gap)")
for m in [0.2, 0.5, 1.0, 2.0]:
    chi, Var, gamma, lo, hi = chi_and_bounds(m, w)
    r = chi * gamma / Var
    print(f"  m={m:4.2f}  chi*gamma/Var = {r:.5f}   (<=1:{r<=1+1e-12},  >=w2/Var:{r>=w[0]/Var-1e-12})")

# ===========================================================================
# (A2) KL-atom form of Lemma G:   f1^2/m1^2 <= chi = int dnu/s <= M/m1^2
#      nu = atom f1^2 at m1^2  +  continuum on [m1^2, hi], total mass M.
# ===========================================================================
banner("(A2) KL-atom Lemma G:  f1^2/m1^2 <= chi=int dnu/s <= M/m1^2")
print("  chi >= f1^2/m1^2  (atom term)  ;  chi <= M/m1^2  (all weight >= m1^2)")
for m1sq, f1sq, M in [(1.0, 0.5, 1.0), (1.0, 0.3, 1.0), (4.0, 0.5, 1.0), (0.25, 0.5, 1.0)]:
    lo, hi = m1sq, 10.0 * m1sq
    sc = np.linspace(lo, hi, 200001)
    dens = (M - f1sq) / (hi - lo)            # uniform in s
    chi = f1sq / m1sq + np.trapezoid(dens / sc, sc)
    print(f"  m1^2={m1sq:5.2f} f1^2={f1sq:.2f} M={M:.2f}: chi={chi:7.5f}  "
          f"[f1^2/m1^2={f1sq/m1sq:.4f}, M/m1^2={M/m1sq:.4f}]  "
          f"in-range={f1sq/m1sq - 1e-9 <= chi <= M/m1sq + 1e-9}")

# ===========================================================================
# (B) Floor-blindness
# ===========================================================================
banner("(B) Floor-blindness: tiny atom at small s blows up chi=int dnu/s, moments barely move")
s0 = np.array([1.0, 4.0]); wt0 = np.array([0.5, 0.5])
js = [0, 1, 2, 3]
def moments(s, w): return np.array([np.sum(w * s**j) for j in js])
M_base = moments(s0, wt0); chi_base = np.sum(wt0 / s0)
print(f"  base: s={s0}, w={wt0}  chi={chi_base:.4f}  M={M_base}")
for eps in [1e-6, 1e-4, 1e-2]:
    sB = np.concatenate([s0, [eps]]); wtB = np.concatenate([wt0, [eps]])
    dM = np.max(np.abs(moments(sB, wtB) - M_base)); chi = np.sum(wtB / sB)
    print(f"  +atom(eps={eps:.0e}): |dM|max={dM:.2e}   chi={chi:.4f}  (jump {chi-chi_base:+.3f})")
print("  => chi is NOT a function of finitely many polynomial moments (the IR tail is free).")

# ===========================================================================
# (C) Certificate equivalence
# ===========================================================================
banner("(C) Certificate  chi <= 1/lambda  <=>  gamma >= lambda*C  (C in [w_2, Var])")
print("  geometry: chi >= w_2/gamma (=> necessary: gamma >= lam*w_2);")
print("            chi <= Var/gamma (=> sufficient: gamma >= lam*Var). C clamped in [w_2, Var].")
lam = 0.25
for m in [0.2, 0.5, 1.0, 2.0]:
    chi, Var, gamma, lo, hi = chi_and_bounds(m, w)
    print(f"  m={m:.2f}: chi={chi:.4f} (chi<=1/lam:{chi<=1/lam})  gamma={gamma:.4f}  "
          f"nec(g>=lam*w2):{gamma>=lam*w[0]}  suf(g>=lam*Var):{gamma>=lam*Var}")

# ===========================================================================
# (D) Non-Gaussianity defect: kappa*chi = 1-eta
# ===========================================================================
banner("(D) Ward/Stein defect:  kappa*chi = 1-eta.  Gaussian eta=0 -> chi=1/kappa;"
       " interacting eta>0 -> chi=(1-eta)/kappa")
def toy(k0, g):
    a = np.linspace(-12, 12, 400001)
    lw = -(0.5 * k0 * a**2 + 0.25 * g * a**4)
    lw -= lw.max(); p = np.exp(lw); p /= np.trapezoid(p, a)
    m1 = np.trapezoid(a * p, a); m2 = np.trapezoid(a**2 * p, a)
    Var = m2 - m1**2
    kappa = k0            # -p''/p|_0 = k0 for this quartic
    return kappa, Var, kappa * Var

for k0, g in [(1.0, 0.0), (1.0, 0.5), (1.0, 1.0), (1.0, 2.0)]:
    kappa, Var, prod = toy(k0, g)
    print(f"  k0={k0:.1f} g={g:.1f}: kappa={kappa:.3f} Var={Var:.4f} "
          f"kappa*Var={prod:.4f}  eta={1-prod:.4f}")

banner("SUMMARY")
print("  (A) chi = sum w_n/(1-lam_n):  w2/gamma <= chi <= Var/gamma   (exact)      [T]")
print("  (B) polynomial/UV moments do NOT determine chi=int dnu/s     (floor-blind) [T]")
print("  (C) chi <= 1/lambda <=> gamma >= lambda*C, C in [w2,Var]      (cert=gap)   [T]")
print("  (D) kappa*chi = 1-eta ; certificate = kappa >= lambda(1-eta)  (H1 mod eta) [T]")
print("  => certificate chi<=1/lambda is H1 verbatim (up to bounded C and uncontrolled eta).")
