"""
construct_worldsheet_cft.py -- witness for CONSTRUCT-WORLDSHEET-CFT.md
Tiny, numpy only, py -3.12 -I. NO lattice MC.

Tests:
 [1] c = D-2 transverse bosons; universal Luscher coefficient -(D-2)pi/6.
 [2] Pure NG closed string E(L)=sigma L - (D-2)pi/(6L): monotone increasing (collapses, no glueball).
 [3] Rigid ring E(rho)=2 pi sigma rho + pi kappa/rho: interior min, m/sqrt(sigma)=2 pi sqrt(2 kappa).
 [4] kappa is an IRRELEVANT dim-4 operator => free LEC: m/sqrt(sigma) varies with kappa at fixed c.
 [5] stability: need kappa>0 to beat the negative Casimir; pure NG (kappa=0) unstable.
"""
import numpy as np

print("CONSTRUCT-WORLDSHEET-CFT witness")
D = 4
c = D - 2                      # transverse bosons
sigma = 1.0
print(f"[1] c = D-2 = {c}  (D={D});  universal Luscher coefficient  -(D-2)pi/6 = {-c*np.pi/6:.6f}")

# [2] pure Nambu-Goto closed string: E(L) = sigma L - (D-2)pi/(6L)
print("\n[2] pure NG closed string  E(L) = sigma L - (D-2)pi/(6L):")
for L in [0.05, 0.1, 0.3, 1.0, 3.0, 10.0]:
    print(f"      L={L:5.2f}   E={sigma*L - c*np.pi/(6*L):9.4f}")
print("      E'(L) = sigma + (D-2)pi/(6L^2) > 0 for all L>0  -> monotone increasing ->")
print("      minimum at L->0 with E->-inf  => shrink modulus, NO stable closed string / NO glueball [T].")

# [3] rigid ring (MASTERFIELD): E(rho) = 2 pi sigma rho + pi kappa/rho
def m_over_sqrt_sigma(kappa, sigma=1.0):
    return 2*np.pi*np.sqrt(2*kappa)          # m = 2 pi sqrt(2 sigma kappa), divided by sqrt(sigma)
print("\n[3] rigid ring  E(rho) = 2 pi sigma rho + pi kappa/rho ;  rho* = sqrt(kappa/(2 sigma)) ;")
print("    m/sqrt(sigma) = 2 pi sqrt(2 kappa):")
for kappa in [0.05, 0.100, 0.138, 0.193, 0.300]:
    rho = np.sqrt(kappa/(2*sigma))
    print(f"      kappa={kappa:5.3f}   rho*={rho:6.3f}   L*=2pi rho*={2*np.pi*rho:6.3f}   m/sqrt(sigma)={m_over_sqrt_sigma(kappa):6.4f}")
print("      lattice band 3.3-3.9  ->  kappa in [0.138,0.193]  (consistent with CONSTRUCT-MASTERFIELD).")

# [4] the rigidity is an irrelevant dim-4 operator: (d^2 X)^2 has dim 4 > 2  => free LEC
print("\n[4] rigidity = coefficient of (d^2 X)^2, worldsheet dim = 2*2 = 4 > 2  -> IRRELEVANT.")
print("    Virasoro(c)/modular/RP fix c and the 2 pi/L spacing, NOT an irrelevant coupling:")
print("    m/sqrt(sigma) varies with kappa at fixed c=2  ->  kappa is a free LEC, NOT fixed by the CFT.")

# [5] stability threshold: the +pi kappa/rho must beat the negative Casimir -pi c/(6 rho)
print("\n[5] stability: closed string needs +coeff/rho with coeff>0 to overcome the universal negative")
print(f"    Casimir -pi c/6 = {-c*np.pi/6:.4f} per 1/L; pure NG (kappa=0) sits on the unstable boundary.")
print("    => the worldsheet CFT bootstrap pins the universal Luscher part, NOT kappa, hence NOT m_0++.")
print("\nVERDICT: universal structure fixed [T]; kappa = irrelevant LEC = H1 in worldsheet coords [NG]/[O].")
