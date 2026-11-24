"""
kappa_eft_match.py -- microscope-match of the rigid-string rigidity kappa.
Model (CONSTRUCT-MASTERFIELD): closed flux ring E(R)=2*pi*sigma*R + pi*kappa/R,
   R*=sqrt(kappa/(2 sigma)),  m_0++ = 2*pi*sqrt(2 sigma kappa)  =>  m/sqrt(sigma)=2*pi*sqrt(2 kappa).
We test three microscopic matching routes and the exact inversion.  Lightweight.
"""
import math

def kappa_from_m(m_over_rs):          # exact inversion of the model
    return (m_over_rs**2)/(8*math.pi**2)
def m_from_kappa(kappa):
    return 2*math.pi*math.sqrt(2*kappa)

print("="*72); print("[A] exact inversion kappa=(m/sqrt(sigma))^2/(8 pi^2)  (model identity)"); print("="*72)
for mr in [3.3,3.5,3.7,3.9]:
    print(f"  m/sqrt(sigma)={mr:4.1f}  ->  kappa={kappa_from_m(mr):.4f}")
print(f"  band 3.3-3.9 <-> kappa in [{kappa_from_m(3.3):.3f},{kappa_from_m(3.9):.3f}]  (matches [0.138,0.193])")

print("\n"+"="*72); print("[B] geometric bending-stiffness: kappa = 2 sigma R_core^2 = 2 c^2, c=R_core*sqrt(sigma)"); print("="*72)
for kappa in [0.138,0.193]:
    c = math.sqrt(kappa/2)
    print(f"  kappa={kappa:.3f}  ->  R_core = {c:.3f}/sqrt(sigma)   (c=O(1): core ~ confinement length)")

print("\n"+"="*72); print("[C] condensate route: tube width w0^2 and core radius carry ln(sigma/Lambda^2) (H1 input)"); print("="*72)
# dual/bag core: R_core^2 ~ (1/(2 pi sigma)) * ln( sigma / Lambda^2 )  (Lambda = IR/condensate scale)
sig = 0.44**2                          # GeV^2, sqrt(sig)=0.44 GeV
for Lam in [0.15,0.20,0.25]:           # GeV  (condensate/IR scale, itself = H1)
    R2 = (1.0/(2*math.pi*sig))*math.log(sig/(Lam**2))
    kap = 2*sig*R2
    print(f"  Lambda={Lam:.2f} GeV: ln(sigma/Lambda^2)={math.log(sig/Lam**2):+.2f}, R_core^2={R2:.3f}/sigma, kappa~{kap:.3f}")
print("  -> kappa drifts with Lambda (= condensate scale = H1); no lower-dimension pinning.")

print("\n"+"="*72); print("[D] pure-NG check: rigidity is redundant (Gauss-Codazzi) and string is unstable"); print("="*72)
print("  dE_0/dR = 2*pi*sigma - pi*kappa/R^2 ; for kappa->0, E(R)=2*pi*sigma*R monotone increasing")
print("  => NG closed string shrinks to R->0 (shrink modulus), NO stable glueball, NO physical kappa.")
print("  One-loop NG determinant fixes coefficients but at the (unstable) kappa=0 boundary.")

print("\nCONCLUSION: no lower-dimension coefficient fixes kappa;")
print("            kappa = (m/sqrt(sigma))^2/(8 pi^2) = the gap in EFT coordinates (H1).")
