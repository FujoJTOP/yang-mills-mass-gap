"""
blast_lattice_cert.py -- witness for BLAST-LATTICE-CERT.md
A rigorous, L-UNIFORM upper bound on Var(a_{k->0}) from compactness, vs the
measured value -- and the limit-localization (L is NOT the wall; a->0 is).

No heavy MC. Reads the existing k0 SU(2) results; computes the compactness ceiling.
"""
import json, math

# ---- (1) the rigorous compactness ceiling ----------------------------------
# a_{k->0} = <v, X>, X = principal log of a configuration in su(2) (branch-safe,
# center-twist corrected): every SU(2) angle theta in [0,pi], so the operator
# norm |X| <= pi, and |a| <= |X| * |v| <= pi (v a unit mode).
# Therefore  Var(a) <= (max a - min a)^2 <= (2*max|a|)^2 <= (2*pi)^2.
C = (2.0 * math.pi) ** 2
print("compactness ceiling  C = (2*pi)^2 =", round(C, 4), "  (rigorous, L-uniform, |a|<=pi)")

# ---- (2) the measured value (SU(2), beta=2.4, L=8, N=400) ------------------
d = json.load(open("research/unconditional/w7/k0_su2_b2.4_L8_n2_results.json"))
invs = [m["full"]["inv"] for m in d["per_mode"] if m["full"]["inv"] > 0]
var = [1.0 / x for x in invs]
print("n modes:", len(var))
print("measured Var(a):  min %.3f  max %.3f  mean %.3f" % (min(var), max(var), sum(var) / len(var)))
print("measured 1/Var(a): min %.3f  max %.3f" % (1/max(var), 1/min(var)))

# ---- (3) the widening gap matters only through the continuum --------------
print()
print("certified interval:   0 <= Var(a) <= %.3f   (contains truth, includes 0)" % C)
print("ceiling / measured-max = %.2fx   (loose by ~an order of magnitude)" % (C / max(var)))
# Gaussian-proxy floor on G'' = 1/Var (MODEL-dependent, NOT a certificate):
print("Gaussian-proxy floor  G'' >= 1/C = %.4f   (weak, model-dependent)" % (1.0 / C))

# ---- (4) why the a->0 step is the wall, not L ------------------------------
print()
print("L-uniformity:  the bound C is INDEPENDENT of L  -> L->infinity is NOT the obstruction.")
print("a->0:          H1 needs the PHYSICAL chi <= 1/lambda; the lattice ceiling C is in")
print("               lattice units, and the lattice->physical map Z(a) is uncontrolled by")
print("               compactness -> the a->0 (dimensional-transmutation) step is the wall.")
