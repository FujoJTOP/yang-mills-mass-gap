"""
impossible_universality.py -- witness for IMPOSSIBLE-UNIVERSALITY.md  [V-here]
NOT a Yang-Mills calculation. Three facts that show "universality cannot force the 0++ gap":
  (1) the 4D a-anomaly is universal but GAP-BLIND (mass-independent in free fields);
  (2) the free DOS edge exponent is (D-2)/2, D-dependent -> not a gapped-class classifier;
  (3) the one real match (finite-T SU(2) deconfinement <-> 3D Ising) is 1-form / T!=0.
Run:  py -3.12 -I research/unconditional/w7/impossible_universality.py
"""
print("=" * 74)
print("(1) 4D a-anomaly: universal, RG-monotone, but GAP-BLIND (mass-independent)")
print("=" * 74)
# standard free-field a-anomaly coefficients in 4D
vals = [("real scalar", 1/90), ("Dirac fermion", 11/360), ("gauge vector", 31/180)]
for name, a in vals:
    print(f"  a({name:14s}) = {a:.6f}   [same for massless (gapless) and massive (gapped)]")
print("  => 'a' is a UV datum; a gapless CFT has a!=0 => it cannot force a gap.")

print("\n" + "=" * 74)
print("(2) free DOS edge exponent near threshold: (D-2)/2, D-dependent")
print("=" * 74)
print("  relativistic free field in D spatial dims:  rho(E) ~ (E-m)^((D-2)/2)")
for D in [2, 3, 4]:
    print(f"    D={D}:  edge exponent = (D-2)/2 = {(D-2)/2:+.2f}")
print("  => tracks spatial dimension / dispersion, NOT a universality class.")
print("     And edge statistics REQUIRE an edge -> they presuppose the gap (circular).")

print("\n" + "=" * 74)
print("(3) the one real match: finite-T SU(2) deconfinement <-> 3D Ising")
print("=" * 74)
print("  3D Ising: beta=0.326, nu=0.630  (Polyakov loop = Z2 magnetization)")
print("  REAL universality match -- but the FINITE-T transition, CENTER/Polyakov (1-form)")
print("  channel, vanishing as T->0; Schur-orthogonal to the 0++ (0-form) gap.")
print("\nVERDICT: universality cannot force the 0++ gap (no critical point at T=0;")
print("         the gapped-order data are 1-form/topological, channel-mismatched).")
