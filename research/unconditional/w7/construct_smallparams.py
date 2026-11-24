"""
construct_smallparams.py -- witness for CONSTRUCT-CONSTRUCTIVE-QFT.md  [V-num]

Quantifies the "three small parameters" of constructive QFT at the PHYSICAL
crossover point (N=3, L=infinity, beta~2.5), and shows each is O(1) there.

  1. Strong coupling: SZZ rigorous gap radius lambda>96 (=g^2 N, d=4: 1/(16*3));
     physical crossover lambda_cross ~ 7.3 -> the expansion parameter is NOT small.
  2. 1/N: large-N master field needs N->inf; at N=3, 1/N = 0.33 -> O(1) correction.
  3. Small circle (Unsal-Yaffe / deformed): needs Lambda*L << 1; decompactifying
     to R^4 means Lambda*L >> 1 -> the small parameter is lost at the target.

No lattice MC. math only. Run: py -3.12 -I construct_smallparams.py
"""

import math

def banner(s):
    print("\n" + "=" * 72)
    print(s)
    print("=" * 72)

# ----------------------------------------------------------------------
# 1. Strong coupling: SZZ radius vs crossover
#    SZZ (arXiv:2204.12737, CMP 400 (2023) 805): mass gap for beta < 1/48,
#    i.e. lambda := g^2 N > 96   (d=4 constant 1/(16*(d-1)) = 1/48).
# ----------------------------------------------------------------------
banner("1 | strong-coupling expansion parameter at the crossover")
beta_szz = 1.0 / 48.0            # rigorous gap radius (d=4)
beta_cross = 2.5                 # physical crossover (bare beta)
lam_szz = 96.0                   # lambda > 96 equivalent (d=4)
lam_cross = 7.3                  # SU(3) Hooft coupling at the crossover
print(f"  rigorous mass-gap radius : beta < 1/48 = {beta_szz:.5f}")
print(f"  physical crossover       : beta ~ {beta_cross}   (lambda ~ {lam_cross})")
print(f"  ratio in beta            : {beta_cross / beta_szz:.1f}x")
print(f"  ratio in lambda          : {lam_szz / lam_cross:.1f}x")
print("  -> the cluster/LSI expansion parameter is FAR from its convergence")
print("     regime at the crossover (needs lambda LARGE; physical lambda is O(1)).")

# ----------------------------------------------------------------------
# 2. 1/N expansion at N=3
# ----------------------------------------------------------------------
banner("2 | large-N expansion parameter at N=3")
for N in (3, 4, 5, 8, 16, 1000):
    print(f"  N = {N:>4d} : 1/N = {1.0/N:.4f}   1/N^2 = {1.0/N**2:.4f}")
print("  -> at N=3, 1/N = 0.333 (O(1)); the leading large-N gap is not")
print("     quantitatively controlled at the physical gauge group.")

# ----------------------------------------------------------------------
# 3. Small circle: Lambda*L decompactification
#    A small-circle deformation has a small parameter only when Lambda*L << 1.
#    The physical target is R^4, i.e. Lambda*L -> infinity.
# ----------------------------------------------------------------------
banner("3 | small-circle deformation: the parameter is lost at the target")
print("  small circle controlled for  Lambda*L << 1  (KK scale >> Lambda)")
for LL in (0.1, 0.5, 1.0, 2.0, 5.0):
    print(f"  Lambda*L = {LL:>4.1f} : KK gap / Lambda ~ 1/(Lambda L) = {1.0/LL:.3f}"
          + ("   (controlled)" if LL < 0.5 else "   (NOT small)"))
print("  -> decompactification R^4 = Lambda*L -> infinity: no small parameter remains;")
print("     the limit back to R^4 is the transport wall (I2).")

# ----------------------------------------------------------------------
# summary
# ----------------------------------------------------------------------
banner("SUMMARY: all three small parameters are O(1) at the physical point")
print("  strong coupling : needs lambda>96, physical lambda~7.3     -> O(1) gap")
print("  1/N             : needs N->inf, physical N=3 (1/N=0.33)    -> O(1)")
print("  small circle    : needs Lambda*L<<1, physical R^4 (>>1)    -> O(1)")
print("  => constructive QFT has NO small parameter at the crossover.")
print("  Missing estimate = uniform (L,a) lower bound on m_0++ at beta=2.5 = H1.")
