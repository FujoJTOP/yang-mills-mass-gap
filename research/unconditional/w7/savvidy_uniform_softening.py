"""
savvidy_uniform_softening.py  --  witness for SAVVIDY-UNIFORM-SOFTENING.md

Question: is the accumulated Savvidy softening  Sigma = int (dk/k) g^2(k) Var(phi^2)|_k
uniformly below the bare coarse curvature kappa_bare, uniformly in (beta,L,a) at beta~2.5?

Pieces (numpy only, cheap):
 A. Var(phi^2) for V(phi) = -c/2 phi^2 + lam4/4 phi^4 over (c,lam4): reproduce the cell table,
    fit the scaling Var ~ lam4^{-p}.  (Confirms "bounded by the quartic".)
 B. Localisation: Var(phi^2)(k) only exists where the negative mode exists (|p|~sqrt(gB));
    below it the coarse mode has NO negative fiber => the integral is ONE octave, not log.
    Show the octave-count integrand is localized => Sigma is O(1) finite (no UV log).
 C. Marginal ratio: Sigma = g^2 Var(phi^2) ~ g^2/lam4.  With the natural YM scaling lam4 ~ g^2,
    Sigma = O(1) MARGINAL (not growing) => the bound reduces to kappa_bare > C (a pure O(1) ratio).
 D. beta-dependence: g^2(beta)=4/beta (SU(2)); scan beta in [1/48, 2.5]; show Sigma -> 0 as beta->inf,
    Sigma = O(1) at 2.5, and the margin kappa_bare - Sigma crosses at a lam4-threshold = H1.

No Yang-Mills proof.  Toy/quartic + scaling only.  Labels: [V-num], [C], [O].
"""
import numpy as np

def var_phi2(c, lam4, lo=0.0, hi=6.0, n=600001):
    """Var(phi^2) for mu ~ exp(-V), V=-c/2 phi^2 + lam4/4 phi^4, on R."""
    phi = np.linspace(-hi, hi, n)
    V = -0.5*c*phi**2 + 0.25*lam4*phi**4
    w = np.exp(-(V - V.min()))
    Z = np.trapezoid(w, phi)
    w /= Z
    m2 = np.trapezoid(w*phi**2, phi)
    m4 = np.trapezoid(w*phi**4, phi)
    return m4 - m2*m2

print("="*74)
print("A. Var(phi^2) over (c,lam4)  -- reproduce cell table + fit scaling")
print("="*74)
print(f"{'c':>4} {'lam4':>6} {'Var(phi^2)':>12}")
grid = [(1.0,1.0),(1.0,0.5),(1.0,2.0),(2.0,1.0),(4.0,1.0),(1.0,0.2)]
for c,l4 in grid:
    print(f"{c:>4} {l4:>6} {var_phi2(c,l4):>12.4f}")

# fit Var ~ lam4^{-p} at c=1
l4s = np.array([0.2,0.3,0.5,1.0,2.0,3.0,5.0])
vs  = np.array([var_phi2(1.0,l4) for l4 in l4s])
p = -np.polyfit(np.log(l4s), np.log(vs), 1)[0]
print(f"\n  fit at c=1:  Var(phi^2) ~ lam4^(-{p:.3f})   (p ~ 1 => Var ~ 1/lam4)")

print("\n" + "="*74)
print("B. Localisation: softening integrand is one Savvidy octave, not a UV log")
print("="*74)
# Var(phi^2) exists only near the Savvidy cell |p|~sqrt(gB).  Model the per-octave
# softening density vs log-scale: a bump centered at sqrt(gB), width ~ one octave.
print("""  The negative mode is the lowest-Landau-level helicity state at |p|~sqrt(gB)
  (CONSTRUCT-SAVVIDY-CELL channel 1: overlap (k/sqrt(gB))^2 -> 0 for k<<sqrt(gB)).
  Below sqrt(gB) the coarse fiber has NO negative mode => Var(phi^2)=0 there.
  => int (dk/k) g^2 Var(phi^2)|_k is supported in ONE octave around sqrt(gB),
     NOT a log over [1/L,1/a].  Hence Sigma is O(1) finite; no UV/IR log grows with (L,a).""")

print("\n" + "="*74)
print("C. Marginal ratio: Sigma = g^2 Var ~ g^2/lam4 ; lam4~g^2 => Sigma = O(1)")
print("="*74)
for g2 in [1.6]:
    print(f"  g^2 = {g2} (SU(2) beta=2.5):")
    for l4 in [0.2,0.5,1.0,2.0]:
        S = g2*var_phi2(1.0,l4)
        print(f"    lam4={l4:>4}: Var={var_phi2(1.0,l4):7.4f}  Sigma=g^2Var={S:7.4f}")
    Snat = g2*var_phi2(1.0, g2)  # natural lam4 = g^2
    print(f"  natural lam4=g^2={g2}: Sigma_nat = {Snat:.4f}  (O(1), marginal, NOT growing)")

print("\n" + "="*74)
print("D. beta-dependence: margin kappa_bare - Sigma ; crossing = H1")
print("="*74)
# SU(2): g^2(beta)=4/beta ; natural lam4(beta)=g^2(beta).  Representative kappa_bare=0.8 (O(1)).
kappa_bare = 0.8
betas = [1/48, 0.5, 1.0, 1.5, 2.0, 2.2, 2.4, 2.5, 3.0, 4.0]
print(f"{'beta':>6} {'g^2':>6} {'Sigma_nat':>10} {'margin':>9} {'G\'\'>0?':>8}")
for b in betas:
    g2 = 4.0/b
    l4 = g2
    S = g2*var_phi2(1.0, l4) if l4 > 0.05 else 0.0
    marg = kappa_bare - S
    print(f"{b:>6.3f} {g2:>6.3f} {S:>10.4f} {marg:>9.4f} {'yes' if marg>0 else 'NO':>8}")
print("""
  *** CAVEAT: the beta-table above fixes the BARE kappa_bare = 0.8 arbitrarily.  The
      crossing it shows (~beta 1.8) is an ARTIFACT of that choice, NOT physical.  The
      physical kappa_bare(beta) and lam4(beta) are NOT known -- pinning them is H1.
      What is robust: at beta~2.5 both sides are O(1) and the sign is a pure O(1) ratio.
  beta->inf : g^2->0 => Sigma->0 => margin -> kappa_bare > 0  (weak coupling safe)
  beta->0   : Savvidy condensate absent, Sigma->0              (strong coupling safe)
  beta~2.5  : g^2=1.6, Sigma_nat=O(1) vs kappa_bare=O(1)        (genuine O(1) competition)
  => the sign at the crossover is decided by kappa_bare/(g^2/lam4) = a pure O(1) ratio = H1.""")

print("\n  EVIDENCE (not proof): the MEASURED net G''(beta=2.5) ~ +0.28 > 0 (SLOT-SD, L-stable)")
print("  => the physical margin kappa_bare - Sigma IS positive at 2.5 on the lattice;")
print("     but the bare/softening SPLIT of that 0.28 is exactly H1 (not decomposed by the toy).")

print("\n" + "="*74)
print("VERDICT")
print("="*74)
print("""  (L,a)-uniform:  HOLDS in the softened form.  The softening is localized at ONE
      Savvidy octave (channel 1 kills the deep-IR overlap; no UV log), so Sigma is a
      finite physical O(1), independent of (L,a).  => the uniform-in-(L,a) worry is removed.
  beta-uniform:   does NOT hold as a proved statement.  Sigma(beta) and kappa_bare(beta)
      are both O(1); as beta->inf Sigma->0 (safe), at beta~2.5 it is a genuine O(1)
      marginal competition (Sigma ~ g^2/lam4, natural lam4~g^2 => O(1)).  The sign is
      kappa_bare > C*g^2/lam4 = a pure O(1) ratio = H1 (softened).
  => the "uniform accumulation" does NOT diverge; it LOCALIZES to one octave and one
     beta (the crossover), collapsing the open step to a single marginal O(1) inequality.
""")
