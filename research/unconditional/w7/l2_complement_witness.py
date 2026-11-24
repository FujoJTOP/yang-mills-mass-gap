"""
l2_complement_witness.py  --  [V-here] witnesses for E2-FLIP-L2-COMPLEMENT.

Three two-minute checks, each a *lemma illustration*, not a Clay computation:
  W1  thin-atom H1-blindness: the H1 *curvature* reading  G'' = 1/chi  stays
      bounded below while the KL edge m1 -> 0 (H1-curvature is NOT sufficient).
  W2  bounded/continuous NONLINEAR functionals of the spectral measure are
      edge-blind to a thin edge atom (value moves by O(w) while m1 -> 0),
      whereas only a SHARP-at-0 functional sees the edge -- and the uniform
      version of that sharp functional IS the anti-concentration (H1) datum.
  W3  2-pt / higher-pt orthogonality: variance (2-pt) and connected 4th
      cumulant (4-pt) are independent statistics -- fixing one does not fix
      the other, so a beyond-2pt functional is orthogonal to the gap's location.

Run:  py -3.12 l2_complement_witness.py
"""
import numpy as np

def chi(s, w):
    """zero-momentum susceptibility  chi = int s^{-1} dnu   (s = mass^2)."""
    s = np.asarray(s, float); w = np.asarray(w, float)
    return float(np.sum(w / s))

def edge(s, w):
    """m1^2 = inf supp nu."""
    s = np.asarray(s, float); w = np.asarray(w, float)
    return float(np.min(s[w > 0]))

print("=" * 70)
print("W1  thin-atom H1-blindness  (nu = w0*delta_eps + w1*delta_S, w0 = eps)")
print("=" * 70)
print(f"{'eps':>10} {'m1=sqrt(eps)':>14} {'chi':>12} {'G\"=1/chi':>12} {'H1(G\">=0.05)':>13}")
S = 4.0
for eps in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6]:
    w0 = eps; w1 = 1.0 - w0
    s = [eps, S]; w = [w0, w1]
    ch = chi(s, w); g = 1.0 / ch; m1 = np.sqrt(edge(s, w))
    print(f"{eps:>10.0e} {m1:>14.4e} {ch:>12.5f} {g:>12.5f} {str(g >= 0.05):>13}")
print("  -> chi -> 1 + 1/S = 1.25 (finite), G'' -> 0.8 (bounded), while m1 -> 0.")
print("     The H1 *curvature* bound is satisfied by a gapless (thin-atom) theory.")
print("     (Necessity is fine: gap => chi <= M0/lam. Sufficiency FAILS without the atom.)")

print()
print("=" * 70)
print("W2  nonlinear-but-bounded spectral functionals are edge-blind to thin atom")
print("=" * 70)
print(f"{'eps':>10} {'m1':>12} {'L_t=Int e^-s dnu':>18} {'particip.1/(sum w^2)':>22} {'entropy -sum wlnw':>18}")
S = 4.0; t = 0.7
for eps in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6]:
    w0 = eps; w1 = 1.0 - w0
    s = np.array([eps, S]); w = np.array([w0, w1])
    m1 = np.sqrt(eps)
    Lt = float(np.sum(w * np.exp(-t * s)))
    part = float((w.sum())**2 / np.sum(w**2))
    ent = float(-np.sum(w * np.log(w)))
    print(f"{eps:>10.0e} {m1:>12.3e} {Lt:>18.6f} {part:>22.6f} {ent:>18.6f}")
print("  -> L_t, participation, entropy all move by O(eps) or less while m1 -> 0.")
print("     Only a SHARP functional 'mass below delta' sees the edge; its uniform")
print("     lower-bounded version ('has support down to delta') IS anti-concentration=H1.")

print()
print("=" * 70)
print("W3  2-pt (variance) and connected-4th cumulant are independent statistics")
print("=" * 70)
rng = np.random.default_rng(0)
Nm = 4_000_000

def stats(x):
    m2 = np.mean(x**2)
    k4 = np.mean(x**4) - 3 * m2**2
    return m2, k4

# (a) standard Gaussian: var 1, kappa4 = 0
xg = rng.standard_normal(Nm)
# (b) two-point +-1 w.p. 1/2 (Rademacher): var 1, kappa4 = E[x^4]-3 = 1-3 = -2
xr = rng.choice([-1.0, 1.0], size=Nm)
# (c) scaled symmetric two-point with var 1 but different 4th moment:
#     x = +-a w.p. p/2 each and +-b otherwise, chosen so E[x^2]=1, E[x^4] != 3
def two_scale_moments(a, b, p):
    m2 = p * a**2 + (1 - p) * b**2
    m4 = p * a**4 + (1 - p) * b**4
    return m2, m4 - 3 * m2**2
# pick a=2, p=0.2 -> var = 0.2*4 + 0.8*b^2 = 1 -> b^2 = (1-0.8)/0.8 = 0.25, b=0.5
a, b, p = 2.0, 0.5, 0.2
xu = rng.choice([-a, a, -b, b], size=Nm, p=[p/2, p/2, (1-p)/2, (1-p)/2])

for name, x in [("Gaussian N(0,1)", xg), ("Rademacher +-1", xr), ("two-scale a=2,b=.5,p=.2", xu)]:
    m2, k4 = stats(x)
    print(f"  {name:>26}: var(2-pt) = {m2:8.4f}   kappa4(4-pt) = {k4:8.4f}")
print("  -> all three have var = 1 (identical 2-pt / identical gap location in a")
print("     single-mode toy) but kappa4 differs (0, -2, finite). A 'kappa4 >= c'")
print("     functional is therefore orthogonal to the 2-pt edge and cannot certify m1;")
print("     a gapless interacting point can have kappa4 != 0, so it also fails acid.")

print()
print("Factory-check summary:")
print("  W1: H1-curvature (1/chi) bounded below while m1->0  => H1-curvature insufficient.")
print("  W2: bounded nonlinear spectral functionals are edge-blind to a thin atom.")
print("  W3: 2-pt and 4-pt statistics are orthogonal => beyond-2pt is not a gap-location datum.")
