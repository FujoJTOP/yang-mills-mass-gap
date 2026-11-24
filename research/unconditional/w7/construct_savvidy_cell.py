"""
construct_savvidy_cell.py -- witness for CONSTRUCT-SAVVIDY-CELL.
Does the Savvidy negative mode (lambda_min = -g|B|) DESTROY the chi floor
    Var(a_{k->0}) <= 1/lambda,
or can its contribution be CONTROLLED (resummed into the floor)?

Two distinct channels by which the negative mode can reach Var(a_{k->0}):

  (1) DIRECT one-mode variance.  The negative mode is the lowest-Landau-level
      gluon helicity state; its fluctuation operator eigenvalue is
          (2n+1) gB + 2 gB * s ,  n=0, s=-1  =>  -gB ,
      so its wavefunction is a lowest-Landau-level mode localized on the
      magnetic length 1/sqrt(gB): its FOURIER support is |p| ~ sqrt(gB), NOT
      p=0.  The k->0 coarse mode has support |p| < k.  Their overlap is
      O((k/sqrt(gB))^2) -> 0 as k -> 0.  => the direct channel VANISHES.

  (2) RG-INDUCED SOFTENING.  Integrating out the negative mode through the
      cubic vertex g a phi^2 softens the effective coarse curvature:
          Gamma''(0) = kappa - g^2 Var(phi^2) ,
      where Var(phi^2) is bounded by the quartic (bimodal effective potential
          V(phi) = -1/2 c phi^2 + 1/4 lam4 phi^4,  c = |lambda_-| > 0 ).
      The floor (Var(a)=1/Gamma'' finite) survives iff kappa > g^2 Var(phi^2).
      The softening is BOUNDED (set by lam4, not by c) => NOT a divergence.

Verdict shape: controlled, not destructive.  Softened H1.  No QED.
Run:  py -3.12 -I construct_savvidy_cell.py
"""
import numpy as np

np.seterr(all='ignore')


def bimodal_var_phi2(c, lam4, N=800001):
    """Var(phi^2) under V = -1/2 c phi^2 + 1/4 lam4 phi^4 (c>0: double well)."""
    if c <= 0 or lam4 <= 0:
        return float('nan')
    phi0 = np.sqrt(c / lam4)               # minima at +/- phi0
    L = 8.0 * phi0 + 6.0
    phi = np.linspace(-L, L, N)
    V = -0.5 * c * phi**2 + 0.25 * lam4 * phi**4
    w = np.exp(-(V - V.min()))
    Z = np.trapezoid(w, phi)
    p2 = np.trapezoid(phi**2 * w, phi) / Z
    p4 = np.trapezoid(phi**4 * w, phi) / Z
    return p4 - p2**2


bar = "=" * 78
print(bar)
print("W2 | RG-INDUCED SOFTENING:  Gamma''(0) = kappa - g^2 Var(phi^2)")
print(bar)
print(f"{'c=|lam_-|':>10} {'lam4':>8} {'Var(phi^2)':>13} {'g^2 Var (g=1)':>14}")
rows = [(1.0, 1.0), (1.0, 0.5), (1.0, 2.0), (2.0, 1.0), (0.5, 1.0),
        (4.0, 1.0), (1.0, 0.2)]
for (c, lam4) in rows:
    v = bimodal_var_phi2(c, lam4)
    print(f"{c:>10.3f} {lam4:>8.3f} {v:>13.5f} {v:>14.5f}")
print()
print("  Reading: as c=|lam_-| grows the negative curvature grows, but Var(phi^2)")
print("  stays O(1/lam4)-bounded (NOT divergent).  The softening is set by the")
print("  quartic stabilizer lam4, not by the negative curvature c.")

print()
print(bar)
print("FLOOR TEST:  Var(a) = 1/(kappa - g^2 Var(phi^2));  soft-then-destroy")
print(bar)
c, lam4, g = 1.0, 1.0, 1.0
soft = g * g * bimodal_var_phi2(c, lam4)
print(f"  softening g^2 Var(phi^2) = {soft:.5f}   (c={c}, lam4={lam4}, g={g})")
for kappa in [3 * soft, 1.5 * soft, 1.05 * soft, 1.0 * soft, 0.9 * soft]:
    gap = kappa - soft
    if gap > 1e-9:
        var_a = 1.0 / gap
        tag = "finite floor"
    else:
        var_a = float('inf')
        tag = "FLOOR DESTROYED (Gamma''<=0)"
    print(f"    kappa={kappa:8.4f}   kappa-soft={gap:+10.5f}   "
          f"Var(a)={var_a:11.4f}   {tag}")
print()
print("  Reading: the negative mode SOFTENS by g^2 Var(phi^2) and destroys the")
print("  floor only when kappa drops BELOW the softening (a critical ratio),")
print("  i.e. it acts as a finite threshold shift -- not an automatic blow-up.")

print()
print(bar)
print("W1 | DIRECT one-mode overlap:  |<a_{k->0} | psi_->|^2 ~ (k/sqrt(gB))^2")
print(bar)
print("  Coarse mode: Fourier support |p|<k.  Negative mode: |p|=p_-=sqrt(gB).")
print("  Disjoint supports for k<p_-  =>  overlap identically 0 (Fourier).")
print("  Smearing hard cutoff at k: overlap^2 ~ min(1,(k/p_-)^2).")
p_neg = 1.0
for k in [0.05, 0.1, 0.25, 0.5, 1.0, 2.0]:
    ov = min(1.0, (k / p_neg) ** 2)
    print(f"    k={k:4.2f}   overlap^2 ~ (k/p_-)^2 = {ov:8.5f}")
print()
print("  Reading: the direct channel VANISHES as k->0.  The Savvidy negative")
print("  mode lives at the magnetic scale sqrt(gB); the k->0 variance is")
print("  orthogonal to it.  This is the analytic face of the W6-BE finding that")
print("  90-97% of the negative weight is high-momentum.")

print()
print(bar)
print("CONTINUUM OUTLOOK")
print(bar)
print("  IR (k << sqrt(gB)): direct channel protected (overlap->0); the induced")
print("  softening is bounded per scale => chi floor SURVIVES (softened H1).")
print("  Crossover (k ~ sqrt(gB)): overlap O(1), softening O(1) => a genuine O(1)")
print("  competition; this is exactly where the Savvidy cell can bite.")
print("  The floor is now a SOFTENED inequality:")
print("      kappa_bare  >  (accumulated Savvidy softening over scales)   =  H1.")
print("  The uniform bound on the accumulated softening is the open H1 step.")
print("  No QED; Clay 0%.")
