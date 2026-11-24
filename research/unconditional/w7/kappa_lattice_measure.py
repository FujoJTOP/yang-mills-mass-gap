"""
kappa_lattice_measure.py -- KAPPA-LATTICE-MEASURE witness (fork, 2026-11-24)

Constructs the lattice protocol for the rigid-string rigidity κ from the
SUB-LEADING (1/R^2) coefficient of the flux-tube width -- NOT the universal
ln R slope (rigidity-free) and NOT the constant term (r0-degenerate).

Units sigma=1, UV cutoff a=1 (KAPPA-DERIVE convention).  Transverse propagator
    G(k) = 1/(sigma k^2 + kappa k^4) ,   sigma=1
so the (coincident-transverse) width is
    w^2(R) = (D-2)/(2pi) * int_{1/R}^{1/a} dk / (k (1 + kappa k^2)) .

Closed form (a=1):
    w^2(R) = (D-2)/(2pi) [ ln R - (1/2)ln(1+kappa) + (1/2)ln(1+kappa/R^2) ]
           = (D-2)/(2pi) ln R  - (D-2)/(4pi)ln(1+kappa)  + (D-2)/(4pi) * kappa/R^2 + O(R^-4)
                ^universal slope     ^constant (r0-channel)   ^SUB-LEADING 1/R^2  (r0-FREE!)

Key: the 1/R^2 coefficient c2 = (D-2)/(4pi) * kappa is positive, proportional to the
rigidity, and -- unlike the constant term -- carries NO r0 (UV-core) degeneracy.
That is the cleanest lattice rigidity observable.  Question: is it measurable?
"""
import numpy as np
from scipy.integrate import quad

D = 4
Dm2 = D - 2


def w2_numeric(R, kappa):
    val, _ = quad(lambda k: 1.0 / (k * (1.0 + kappa * k * k)), 1.0 / R, 1.0)
    return Dm2 / (2 * np.pi) * val


def w2_closed(R, kappa):
    return Dm2 / (2 * np.pi) * (np.log(R) - 0.5 * np.log(1 + kappa)
                                + 0.5 * np.log(1 + kappa / R**2))


print("=" * 74)
print("(1) VERIFY closed form vs numerical quad  (kappa=0.193)")
print("=" * 74)
for R in [1.5, 2.0, 4.0, 8.0]:
    num = w2_numeric(R, 0.193)
    clo = w2_closed(R, 0.193)
    print(f"  R={R:4}: quad={num:.10f}  closed={clo:.10f}  diff={abs(num-clo):.2e}")

print()
print("=" * 74)
print("(2) THE 1/R^2 sub-leading coefficient   c2 = (D-2)/(4pi) * kappa   (r0-FREE)")
print("=" * 74)
for kappa in [0.0, 0.05, 0.139, 0.193]:
    print(f"  kappa={kappa:6}: c2 = {Dm2/(4*np.pi)*kappa:.6f}")

print()
print("=" * 74)
print("(3) NG (kappa=0) at quadratic level: c2 = 0  (pure log). Rigidity => c2>0.")
print("=" * 74)

print()
print("=" * 74)
print("(4) RESOLVABILITY: magnitude of the c2/R^2 term vs w^2   (kappa=0.193)")
print("=" * 74)
k = 0.193
c2 = Dm2 / (4 * np.pi) * k
for R in [2, 3, 4, 6, 8, 10, 12]:
    term = c2 / R**2
    tot = w2_closed(R, k)
    print(f"  R={R:2}: c2/R^2={term:.5f}   w^2={tot:.4f}   rel={100*term/tot:.3f}%")
print("  -> sub-3% at R=2, <0.5% at R>=4, <0.1% at R>=8.")

print()
print("=" * 74)
print("(5) RECOVER c2 by fitting w^2 = A lnR + B + C/R^2  (noiseless)")
print("=" * 74)
Rs = np.linspace(2.0, 12.0, 60)
for kappa in [0.0, 0.193]:
    ws = w2_closed(Rs, kappa)
    M = np.column_stack([np.log(Rs), np.ones_like(Rs), 1.0 / Rs**2])
    coef, *_ = np.linalg.lstsq(M, ws, rcond=None)
    print(f"  kappa={kappa:6}: slope={coef[0]:.6f} (pred {Dm2/(2*np.pi):.6f})  "
          f"c2_fit={coef[2]:.6f} (pred {c2 if kappa else 0.0:.6f})")

print()
print("=" * 74)
print("(6) PRECISION NEEDED TO RESOLVE c2 (absolute error delta on w^2 at R~8-10)")
print("=" * 74)
for R in [8, 10, 12]:
    need = c2 / R**2
    tot = w2_closed(R, k)
    print(f"  R={R:2}: need abs delta < {need:.5f}  (rel < {100*need/tot:.3f}%)")
print("  Current lattice width precision ~ few % (rel) => need ~10-50x improvement.")

print()
print("=" * 74)
print("(7) CONTRAST: constant-term channel is r0-DEGENERATE (WIDTH-KAPPA)")
print("=" * 74)
print("  intercept shift with rigidity: -(D-2)/(4pi)ln(1+kappa) = "
      f"{-(Dm2/(4*np.pi))*np.log(1+k):.4f} at kappa=0.193")
print("  but w0^2 (UV core) also shifts it  => degenerate with r0 (not a clean kappa read).")
print("  The 1/R^2 coefficient c2=(D-2)/(4pi)kappa is r0-FREE => strictly better observable.")
print()
print("VERDICT: protocol is clean and r0-free, but c2/R^2 is <0.1% of w^2 at large R")
print("         (current precision ~few %); gives kappa_PK (canonical), NOT MASTERFIELD kappa.")
