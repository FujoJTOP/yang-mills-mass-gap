"""
impossible_invertibility.py -- witness for IMPOSSIBLE-INVERTIBILITY.md

Base: gap = invertibility of (I - T_phys) on the EXCITATION space, via a
degree / index argument (operator space), NOT spectral positivity of a measure.

(A) degree = sign det_{Omega^perp}(I - T/lambda0) = prod_{i>=1} (1 - e^{-a(E_i-E0)})
    >= 0, == 0 iff gapless  =>  the degree is only 0 or +1 (Krein-Rutman) = the gap's indicator.
(B) Homotopy: the indicator drops to 0 at a level crossing => no nonzero integer index is
    conserved across beta (rigidity: the gap is a non-analytic jump, unprotected).
(C) Nonlinear/Leray-Schauder toy: F(k)=lambda*k/(1+k^2); F(0)=0 always (gapless); a
    nonzero fixed point is 'degree-forced' iff lambda>1 = the tachyon = the gap (circular).

numpy only.  Run:  py -3.12 -I research/unconditional/w7/impossible_invertibility.py
"""
import numpy as np

def banner(s):
    print("\n" + "=" * 72)
    print(s)
    print("=" * 72)

def analyze(spec, a=1.0):
    """det of (I - T/lambda0) on the EXCITATION space (vacuum/Perron mode removed)."""
    E = np.array(sorted(float(x) for x in spec))
    gaps = E[1:] - E[0]                       # excitation energies E_i - E0
    det = float(np.prod(1.0 - np.exp(-a * gaps)))
    return E[1] - E[0], det                    # gap, det

# ---------------------------------------------------------------- (A)
banner("(A) finite-dim: degree = sign det_{Omega^perp}(I - T/lambda0)  vs the gap")
print("    (vacuum/Perron mode removed; det = prod_{i>=1}(1 - e^{-a(E_i-E0)}))\n")
cases = [
    ("gapped        (delta=1)",    [0.0, 1.0, 4.0, 9.0]),
    ("near-gapless  (delta=1e-3)", [0.0, 1e-3, 4.0, 9.0]),
    ("gapless       (delta=0)",    [0.0, 0.0, 4.0, 9.0]),
]
for name, spec in cases:
    g, d = analyze(spec)
    print(f"  {name:28s}: delta={g:.4f}  det={d:+.4e}  degree={np.sign(d):+.0f}")
print("  -> det >= 0, ==0 exactly when delta=0.  degree in {0,+1}: it IS the gap's indicator.")
print("     (Krein-Rutman: a positive compact T has one simple Perron eigenvalue, all others <1.)")

# ---------------------------------------------------------------- (B)
banner("(B) homotopy: a level crossing; the indicator passes through 0 (no conserved index)")
def spectrum_at(s):
    return [0.0, s - 0.5, 1.5, 3.0]           # E1(s) crosses E0=0 at s=0.5
for s in np.linspace(0.0, 1.0, 11):
    g, d = analyze(spectrum_at(s))
    tag = "   <-- gap closes (crossing)" if abs(s - 0.5) < 1e-9 else ""
    print(f"  s={s:.1f}: E1-E0={g:+.3f}  degree={np.sign(d):+.0f}{tag}")
print("  -> the degree is +1 away from the crossing and 0 AT it: it does NOT stay nonzero.")
print("     No integer index is conserved across the homotopy (rigidity B3: the gap is a")
print("     non-analytic beta-jump; BASE-INDEX: indices are floor-blind).")

# ---------------------------------------------------------------- (C)
banner("(C) nonlinear / Leray-Schauder toy: F(k)=lambda*k/(1+k^2)")
def fixed_points(lam):
    out = [0.0]
    if lam > 1.0:
        out.append(np.sqrt(lam - 1.0))       # F(k)=k  <=>  k^2 = lambda-1
    return out
for lam in [0.5, 1.0, 2.0]:
    fps = fixed_points(lam)
    kind = "gapless only" if lam <= 1.0 else "NONZERO fixed point exists"
    print(f"  lambda={lam:.1f}: fixed points k*={['%.3f' % x for x in fps]}  ({kind})")
print("  -> a nonzero fixed point is 'degree-forced' iff lambda>1.")
print("     But lambda>1 is exactly the k->0 tachyon / the gap. So the nonlinear degree")
print("     argument forces the gap only after assuming the tachyon = the gap (circular).")
print("     And the physical F is the infinite DSE/ERG tower: NOT compact => Leray-Schauder")
print("     degree is undefined; any truncation admits the gapless fixed point k=0.")

# ---------------------------------------------------------------- acid
banner("acid: invertibility distinguishes free-massless from gapped -- but IS the gap")
for m in [0.0, 0.5]:
    g, d = analyze([0.0, m, 1.0, 4.0])       # m=0 -> E1=E0 gapless; m>0 -> gapped
    state = "singular (free/gapless)" if abs(d) < 1e-12 else "invertible (gapped)"
    print(f"  m={m}: delta={g:.3f}  det={d:+.4e}  -> {state}")
print("  -> acid-clean in FORM (det=0 iff gapless), but det = prod(1 - e^{-a*delta_i}) = the")
print("     gap. The degree carries no information beyond the sign of the gap: a reformulation,")
print("     not a handle.")
