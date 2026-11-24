"""blast_savvidy_close.py -- witness for BLAST-SAVVIDY-CLOSE.md
Assemble the Savvidy-cell H1 attack: Var(phi^2) bound, (L,a)-uniformity,
and the multi-cell IR-marginal accumulation. Decide close via kappa_bare > Sigma.
Run: py -3.12 -I blast_savvidy_close.py
"""
import numpy as np

C_RG = 24*np.pi**2/11.0   # SU(2) one-loop, = 21.53 (matches lean C_RG)
print(f"C_RG (SU(2)) = {C_RG:.4f}\n")

# ---- (1) Var(phi^2) bound: 2/lam4 + 2/(gB)^2 + ... -> 2/lam4, c-independent
def var_phi2(c_gB, lam4):
    # exact-ish saddle form from SAVVIDY-VAR-PHI2: 2/lam4 + 2/c^2 + O(c^-4)
    return 2.0/lam4 + 2.0/c_gB**2
print("(1) Var(phi^2) vs c=gB at lam4=1  (-> 2/lam4 = 2, c-INDEPENDENT):")
for c in [10,20,40,80,160]:
    print(f"    c={c:4d}  Var(phi^2)={var_phi2(c,1.0):.6f}")

# ---- (2) softening g^2 Var(phi^2) O(1), g-INDEPENDENT at fixed a = gB/sqrt(lam4)
# physical scaling lam4 = k*g^2, c=gB, a = c/sqrt(lam4) = B/sqrt(k)
print("\n(2) softening sig = g^2 Var(phi^2) at fixed a=1, varying g (k=1):")
for g in [1.0,0.7,0.5,0.4,0.3]:
    k=1.0; lam4=k*g*g; c=g*1.0/np.sqrt(k)   # a = c/sqrt(lam4) = B/sqrt(k) with B=1 -> keep a fixed
    # fix a: B such that a=1 => B=sqrt(k). c=g*B
    B=np.sqrt(k); c=g*B; lam4=k*g*g
    sig=g*g*var_phi2(c,lam4) if c>0 else np.nan
    print(f"    g={g:.1f}  lam4={lam4:.3f} c={c:.3f}  g^2 Var={sig:.5f}")

# ---- (3) multi-cell accumulation Sigma = (2 C_RG / lam4) * ln( ln(kUV/Lam)/ln(kIR/Lam) )
def Sigma_multicell(lam4, kUV_over_Lam, kIR_over_Lam):
    num=np.log(kUV_over_Lam); den=np.log(kIR_over_Lam)
    if den<=0: return np.inf
    return (2*C_RG/lam4)*np.log(num/den)
print("\n(3) multi-cell Sigma (lam4=1, kUV=100*Lam) vs kIR -> Landau pole (IR-MARGINAL):")
for kIR in [10,5,3,2,1.5,1.2,1.1,1.05,1.01]:
    S=Sigma_multicell(1.0,100.0,kIR)
    print(f"    kIR/Lam={kIR:5.2f}  Sigma={S:9.3f}")

# ---- (4) does it close: margin kappa_bare - Sigma for natural kappa_bare ~ O(1)
print("\n(4) margin kappa_bare - Sigma  (overshoot when < 0):")
for kb in [1.0,2.0,5.0]:
    row=[]
    for kIR in [2.0,1.5,1.2,1.1]:
        row.append(f"kIR={kIR}:{kb-Sigma_multicell(1.0,100.0,kIR):7.2f}")
    print(f"    kappa_bare={kb:.0f}:  "+"  ".join(row))

# ---- (5) ONE-OCTAVE alternative: if negative mode exists in ONE octave only,
# Sigma_oct = g^2 Var * ln(2)  (bounded) -- the competing picture
print("\n(5) ONE-OCTAVE alternative Sigma_oct = g^2 Var * ln2 (bounded):")
for sig in [0.5,0.95,2.02]:
    print(f"    g^2 Var={sig:.2f}  Sigma_oct={sig*np.log(2):.4f}")

# ---- crossover kIR where Sigma = kappa_bare (the exact step, for kappa_bare=1)
print("\n(6) crossover kIR* where Sigma(kIR*)=kappa_bare=1 (solve):")
from scipy.optimize import brentq
f=lambda kIR: Sigma_multicell(1.0,100.0,kIR)-1.0
kIR_star=brentq(f,5.0,99.99)
print(f"    kIR* = {kIR_star:.4f} * Lam   (Sigma=k_bare here; below it Sigma>k_bare -> OVERSHOOT)")
print(f"    => sub-criticality kappa_bare>Sigma holds ONLY for kIR>kIR*, i.e. NOT down to the Landau pole.")
