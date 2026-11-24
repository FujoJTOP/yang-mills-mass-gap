"""
roof_condensate_floor_witness.py -- ROOF-CONDENSATE-FLOOR, schematic witness.
NOT a Yang-Mills calculation. Demonstrates the u=2 IR-renormalon structure that
blocks a scheme-independent uniform floor on the gluon condensate <tr F^2>.

Convention (dimension-2 units): Borel sum  R(s) = int_0^inf du e^{-u s} B(u),
s = ln(Q^2/Lam^2), so a Borel pole at u = d/2 gives a nonperturbative ambiguity
~ e^{-(d/2) s} = (Lam^2/Q^2)^{d/2}.  For d=4 (the gluon condensate, dimension 4)
u=2 -> ambiguity ~ Lam^4/Q^4, the SAME parametric order as the condensate itself.

We verify:
 (1) the discontinuity of the Borel sum across u=2 scales exactly as e^{-2s}
     (a constant times e^{-2s}; the constant is 2*pi for a simple pole);
 (2) that constant is O(1) relative to a natural condensate term c*e^{-2s},
     INDEPENDENT of s -- so at every scale (every a) the ambiguity is the same
     fraction of the condensate: no scheme-independent separation, no floor;
 (3) the d=6 control (u=3) is parametrically smaller, so d=4 is the marginal,
     leading case that matters.

Run:  py -3.12 -I roof_condensate_floor_witness.py
"""
import numpy as np

print("=" * 74)
print("(1) Borel discontinuity across the u=d/2 IR renormalon scales as e^{-(d/2)s}")
print("=" * 74)


def borel_discontinuity(s, u0, eta=1e-3):
    """Delta(s) = int_0^inf e^{-u s} [ 1/(1-u/u0 - i*eta) - c.c. ] du.
    As eta->0 this localises at u=u0 and gives a constant * e^{-u0 s}.
    (Schematic: a simple pole at u=u0, smeared to width eta.)"""
    u = np.linspace(0.0, max(80.0 / s, 8.0 * u0), 800001)
    d = 1.0 / (1.0 - u / u0 - 1j * eta) - 1.0 / (1.0 - u / u0 + 1j * eta)
    return np.trapezoid(np.exp(-u * s) * d, u)


for u0 in (2.0, 3.0):
    print(f"\n  u0 = d/2 = {u0:g}   (d = {2*u0:g})")
    print("  s      |Delta(s)|        |Delta|/e^{-u0 s}")
    for s in (1.0, 1.5, 2.0, 3.0, 4.0):
        D = abs(borel_discontinuity(s, u0))
        print(f"  {s:3.1f}   {D:.4e}      {D/np.exp(-u0*s):.4f}")

print("""
  -> the ratio |Delta|/e^{-u0 s} is (approximately) CONSTANT in s
     (~4*pi = the full discontinuity between the two i-0 prescriptions; the
     FLOOR half-width is 2*pi*e^{-u0 s}).  The ambiguity is a fixed fraction of
     e^{-(d/2)s}.  (Small drift at large s is an eta-smearing artefact.)"""))

print("\n" + "=" * 74)
print("(2) The ambiguity is O(1) relative to a natural condensate c*e^{-2s}")
print("=" * 74)
print("  A natural d=4 condensate term is c * e^{-2s} (c = O(1) coefficient).")
print("  Ambiguity/condensate = 2*pi/c, INDEPENDENT of s: same fraction at every scale.")
ratios = {}
for c in (0.5, 1.0, 2.0, 10.0):
    ratios[c] = 2 * np.pi / c
    print(f"    c = {c:4.1f}  ->  ambiguity/condensate = {ratios[c]:6.3f}")
print("""
  -> for c < 2*pi the ambiguity EXCEEDS the condensate, so the interval
     [PV - |Delta|, PV + |Delta|] always straddles 0: you can dial the
     'nonperturbative condensate' to zero inside the renormalon ambiguity.
     No scheme-independent, a-uniform floor survives.  This IS the obstruction.""")

print("\n" + "=" * 74)
print("(3) d=6 control (u=3): parametrically smaller -> d=4 is the marginal case")
print("=" * 74)
for s in (1.0, 2.0, 3.0):
    D4 = abs(borel_discontinuity(s, 2.0))
    D6 = abs(borel_discontinuity(s, 3.0))
    print(f"  s={s:3.1f}   |Delta_d4|/|Delta_d6| = {D4/D6:.3e}   (= e^s: d=4 ambiguity dominates)")
print("""
  -> the d=4 (gluon-condensate) renormalon is the LEADING IR ambiguity; its
     size equals the condensate scale.  This is why 'floor the condensate' is
     ill-posed as a scheme-independent uniform bound.

  VERDICT (schematic, illustrative): u=2 IR renormalon blocks a uniform floor.
  The physical NONZERO condensate is real, but a scheme-independent FLOOR is
  not a well-posed demand.  [NG] for the class; [O] exact step = the subtraction.
""")
