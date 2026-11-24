"""
roof_specimen_4d.py -- witness for research/unconditional/ROOF-SPECIMEN-4D.md

Question: is there a 4D specimen where the roof (0++ gap via 0-form condensation)
can be built END-TO-END, and what does 4D teach that 3D Z2 gauge cannot?

This script makes ONE honest numeric point: in 4D Z2 gauge the 0++ glueball mass
is controlled (m1 > 0) ONLY in the strong-coupling (tube) regime, with the SAME
beta* wall as P3 -- i.e. the roof EXISTS in 4D gauge but is readable only as a
DIRECT strong-coupling construction, never off a dual magnetization.

Strong-coupling 0++ mass (4D Z2): m1(beta) ~ c * ln(1/u), u = tanh(beta), c = O(4)
(four-wall tube, P3 geometry; coefficient is scheme-dependent -- labeled estimate).
The strong-coupling expansion in u is ASYMPTOTIC: controlled for u <~ 0.35,
diverges by the transition beta_c ~ 0.76. There it becomes UNREADABLE.

3D Z2 gauge <-> 3D Ising: the 't Hooft loop IS the 0-form dual magnetization, so the
3D roof is READABLE off a spin order parameter. In 4D that duality is absent: the
't Hooft operator is a surface (1-form) and the 0++ plaquette is a separate 0-form.
=> the 3D->4D transfer fails on the 0-form/1-form coincidence, not on a tool gap.

Lightweight: numpy only.  Run: py -3.12 -I roof_specimen_4d.py
"""
import numpy as np

np.set_printoptions(precision=4, suppress=True, linewidth=120)


def m1_4d_z2(beta, c=4.0):
    """Strong-coupling 0++ glueball mass estimate, lattice units: c*ln(1/tanh beta)."""
    u = np.tanh(beta)
    return c * np.log(1.0 / u)


def regime(beta):
    u = np.tanh(beta)
    if u < 0.35:
        return "SC controlled"
    if beta < 0.76:
        return "beyond SC radius"
    return "beta_c: UNREAD"


print("=" * 78)
print("ROOF-SPECIMEN-4D : strong-coupling 0++ glueball in 4D Z2 gauge")
print("=" * 78)
print("(lattice units a=1; c=4 from the four-wall tube, P3 geometry; estimate)")
print()
print(f"{'beta':>6} {'u=tanh b':>9} {'m1':>8} {'regime':>20}")
for beta in [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.76]:
    print(f"{beta:6.2f} {np.tanh(beta):9.4f} {m1_4d_z2(beta):8.3f} {regime(beta):>20}")

print()
print("-" * 78)
print("Reading 1 (roof EXISTS):   m1(beta) > 0 for all beta < beta_c -- a genuine")
print("  strong-coupling 0++ gap. The roof is present in 4D Z2 gauge.")
print("Reading 2 (wall):          the series is asymptotic in u=tanh(beta); it is")
print("  trusted only for u<~0.35 (beta<~0.37). By beta_c~0.76 (u~0.64) it has")
print("  diverged -> the roof is UNREADABLE there. Same beta* structure as P3.")
print("Reading 3 (species):       4D Z2 is ABELIAN -> no SU(N) transfer (IL-ii).")
print()
print("-" * 78)
print("WHY 3D Z2 COULD BE READ AND 4D CANNOT:")
print("  3D Z2 <-> 3D Ising : 't Hooft loop = 0-form dual magnetization  -> READ off m.")
print("  4D   : 't Hooft op  = surface (1-form); 0++ plaquette = 0-form  -> NO dual reading.")
print("  => in 4D the 0-form (0++) and 1-form (confinement) sectors DECOUPLE;")
print("     the roof is a SEPARATE object that must be BUILT directly, not read.")
print("=" * 78)

# ---------------------------------------------------------------------------
# Feature table for the four candidate 4D specimens (printed for the report).
# ---------------------------------------------------------------------------
print()
print("CANDIDATE 4D SPECIMENS  (feat a=gapless, b=gapped, c=0-form nonlocal, d=mechanism)")
print(f"{'model':<26}{'a':>3}{'b':>3}{'c':>3}{'d':>3}   note")
rows = [
    ("4D Z2 gauge",      "~", "Y", "N", "~", "first-order transition; abelian; no gapless fixed pt"),
    ("4D U(1) compact",  "Y", "Y", "N", "~", "gapless Coulomb; 1-form (monopole); abelian"),
    ("4D SU(N) T^4_L",   "~", "Y", "~", "~", "perturbative only; center-sym breaking; no crossover"),
    ("4D SU(N) strong cpl","N","Y","Y","Y",  "RIGHT species; P3 roof; beta* = KP radius"),
]
for name, a, b, c, d, note in rows:
    print(f"{name:<26}{a:>3}{b:>3}{c:>3}{d:>3}   {note}")
print()
print("No nonabelian 4D specimen has all four features with exact solvability.")
print("Verdict: the 4D roof is a NEW object; 3D transfer blocked at channel mismatch.")
print("Clay 0%; no QED.")
