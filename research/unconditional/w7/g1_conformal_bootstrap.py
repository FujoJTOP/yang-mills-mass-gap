"""
g1_conformal_bootstrap.py -- SCHEMATIC witness for G1-CONFORMAL-BOOTSTRAP.

This is NOT a real numerical bootstrap (a genuine 4D crossing/positivity bootstrap
is a semidefinite program, far beyond a cheap witness). It only prints the exact
Hofman-Maldacena collider band, the free-field central charges, and demonstrates that
crossing + unitarity admit a continuum of generalized-free-field (GFF) solutions --
so they neither single out nor exclude a hypothetical non-SUSY "pure-glue" CFT.

Run: py -3.12 -I research/unconditional/w7/g1_conformal_bootstrap.py
"""
from fractions import Fraction as F

print("=" * 72)
print("(1) Hofman-Maldacena collider bound (unitary 4D CFT, unique stress tensor):")
lo, hi = F(1, 3), F(31, 18)
print("    1/3 = %.5f  <=  a/c  <=  31/18 = %.5f   (width %.4f)" % (float(lo), float(hi), float(hi - lo)))
print("    -> a WIDE necessary band; a pure-glue CFT would sit inside it; excludes nothing.")

print()
print("(2) Free-field central charges (commonly-quoted; standard, illustrative):")
# free real scalar: a = c = 1/90  (well established)
a_s = c_s = F(1, 90)
print("    real scalar : a=%s, c=%s, a/c=%.3f" % (a_s, c_s, float(a_s / c_s)))
# free Dirac fermion: a = 11/360, c = 1/20 = 18/360  (standard)
a_f, c_f = F(11, 360), F(18, 360)
print("    Dirac ferm. : a=%s, c=%s, a/c=%.3f" % (a_f, c_f, float(a_f / c_f)))
print("    -> both inside the band; the bound does not discriminate.")

print()
print("(3) GFF crossing illustration -- crossing+unitarity admit a CONTINUUM:")
print("    a generalized free field of scalar dimension Delta is crossing-symmetric")
print("    for every Delta, with double-trace spectrum {2n.Delta}. Unitarity needs Delta>1.")
for D in [F(1, 2), 1.0, 1.2, 2.0, 3.0, 5.0]:
    Df = float(D)
    unit = Df > 1.0
    dt = [2 * n * Df for n in range(1, 4)]
    print("    Delta=%.2f : unitary=%s, double-trace dims=%s" % (Df, str(unit), ["%.2f" % x for x in dt]))
print()
print("=> crossing + positivity alone admit a continuum of GFF solutions and do NOT")
print("   exclude a multi-operator pure-glue CFT. A decisive SDP is beyond this witness.")
print("   Structurally, even an abstract-CFT verdict would not decide whether pure YM")
print("   FLOWS to it (= the dynamical gap). Thus (G1) is [O], not closed by the bootstrap.")
