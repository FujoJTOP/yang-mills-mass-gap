"""
b1_nf_monotone.py -- Attack surface B: an N_f-monotonicity theorem for the SU(3) pure-YM gap.

Pure algebra + a toy no-go witness. Exact rational arithmetic (fractions), no lattice, no fitting.
Conventions (campaign / G1-ACID):  beta(g) = -b0 g^3 - b1 g^5,  SU(3):
    b0(Nf) = 11 - (2/3) Nf          (asymptotic freedom boundary b0=0 at Nf = 33/2 = 16.5)
    b1(Nf) = 102 - (38/3) Nf        (Banks-Zaks edge       b1=0 at Nf = 306/38 = 153/19)
The fixed point of beta = -b0 g^3 - b1 g^5 sits at g*^2_abs = -b0/b1 (absorbed coupling; scheme-dependent).
a := -b0/b1 = g*^2_MSbar/(16 pi^2) is the small parameter, sqrt(a) = g*_MSbar/(4 pi).
Leading IRFP mass anomalous dimension (MS-bar): gamma_m* = 6 C_F a + O(a^2).

Run:  py -3.12 research/v3/w7/b1_nf_monotone.py
"""
from fractions import Fraction as Q
import math

def b0(x):  return Q(11)  - Q(2, 3) * Q(x)
def b1(x):  return Q(102) - Q(38, 3) * Q(x)

def sign(v): return '+' if v > 0 else ('0' if v == 0 else '-')

def gstar2(x):
    """a = -b0/b1, the absorbed FP coupling^2 = g*^2_MSbar/(16 pi^2) (defined when b0>0 and b1<0)."""
    a, b = b0(x), b1(x)
    return (a > 0 and b < 0), (-a / b)

print("=== 1. two-loop SU(3) coefficients, exact rationals, Nf = 0..17 ===")
print(f"{'Nf':>4} {'b0':>9} {'b1':>10} {'sgn b0':>7} {'sgn b1':>7}  {'2-loop IRFP?':>12}")
for Nf in range(0, 18):
    has, g2 = gstar2(Nf)
    extra = f"a={float(g2):.4f} (sqrt a={math.sqrt(float(g2)):.3f})" if has else "none (beta<0 all g)"
    print(f"{Nf:>4} {float(b0(Nf)):>9.4f} {float(b1(Nf)):>10.4f} "
          f"{sign(b0(Nf)):>7} {sign(b1(Nf)):>7}  {extra:>12}")

print("\n=== 2. exact sign-change points and the Nf<=8 no-zero theorem ===")
print(f"  b0 = 0 exactly at Nf = 33/2 = {float(Q(33,2))}  (AF upper edge)")
print(f"  b1 = 0 exactly at Nf = 306/38 = 153/19 = {float(Q(153,19)):.10f}  (BZ lower edge)")
print("  For Nf <= 8 (=integer): b0(8)=%s=%.4f>0 and b1(8)=%s=%.4f>0 ;" %
      (b0(8), float(b0(8)), b1(8), float(b1(8))))
print("  b0,b1 both strictly decreasing in Nf => for all Nf<=8, b0>=17/3>0, b1>=2/3>0.")
print("  => beta(g) = -g^3 (b0 + b1 g^2) < 0 for every g>0 : NO positive zero (exact).")

print("\n=== 3. beta(g) has no positive zero for Nf<=8 (numerical confirmation) ===")
def beta(g, Nf): return -float(b0(Nf)) * g**3 - float(b1(Nf)) * g**5
grid = [i * 0.01 for i in range(1, 4001)]  # g in (0,40]
for Nf in (0, 2, 6, 8):
    allneg = all(beta(g, Nf) < 0 for g in grid)
    print(f"  Nf={Nf}: beta<0 on (0,40]? {allneg}   beta(1)={beta(1,Nf):.3e} beta(10)={beta(10,Nf):.3e}")

print("\n=== 4. EXACT monotonicity of a=-b0/b1 in the window: d a/dNf = -642/(38 Nf - 306)^2 < 0 ===")
print(f"{'Nf':>5} {'a':>12} {'sqrt a':>9} {'analytic slope':>18} {'finite-diff slope':>18}")
for Nf in [Q(82,10), Q(9), Q(10), Q(12), Q(14), Q(16), Q(164,10)]:  # 8.2,9,...,16.4
    has, g2 = gstar2(Nf)
    if not has: continue
    an = -Q(642) / (Q(38) * Nf - Q(306))**2                 # exact derivative
    # central finite difference with h=1/100 in exact rationals
    h = Q(1, 100)
    _, g2p = gstar2(Nf + h); _, g2m = gstar2(Nf - h)
    fd = (g2p - g2m) / (2 * h)
    print(f"{float(Nf):>5.2f} {float(g2):>12.6f} {math.sqrt(float(g2)):>9.4f} "
          f"{float(an):>18.6f} {float(fd):>18.6f}")
print("  a: +inf as Nf->153/19^+ ,  ->0 as Nf->33/2^- ; strictly DECREASING in Nf on the window.")

print("\n=== 5. universal IRFP datum gamma_m* (leading order) : increasing as Nf decreases ===")
CF = Q(4, 3)  # SU(3)
print(f"{'Nf':>5} {'a=-b0/b1':>14} {'gamma_m* ~ 6 CF a':>18}  note")
for Nf in [9, 10, 11, 12, 13, 14, 15, 16]:
    has, g2 = gstar2(Nf)
    a = g2  # this already is -b0/b1 = g*^2/(16pi^2) in the MS-bar convention
    gm = 6 * CF * a
    note = "pert. weak" if Nf >= 15 else ("near lower edge, O(a) unreliable" if Nf <= 10 else "")
    print(f"{Nf:>5} {float(a):>14.5f} {float(gm):>18.4f}  {note}")
print("  gamma_m*(Nf) vanishes at the Gaussian edge (Nf->16.5, g*->0) and GROWS as Nf falls.")
print("  gamma_m* is a SCHEME-INDEPENDENT CFT datum, but it exists ONLY where a fixed point does.")

print("\n=== 6. pointwise beta monotonicity in Nf (exact, 2-loop): beta_{Nf+2}(g) - beta_{Nf}(g) > 0 ===")
for g in (0.3, 1.0, 2.0):
    for Nf in (0, 6, 10):
        d = beta(g, Nf + 2) - beta(g, Nf)   # = (4/3) g^3 + (76/3) g^5 > 0
        print(f"  g={g} Nf={Nf}->{Nf+2}:  delta beta = {d:.5f}  (>0: {d > 0})")
print("  => at every fixed g>0 the 2-loop beta is strictly INCREASING in Nf (more screening).")
print("  This is the strongest EXACT Nf-monotonicity, but it is a statement about the truncated")
print("  polynomial, not about the nonperturbative zero structure (the G1-ACID lesson).")

print("\n=== 7. THE NO-GO TOY: an order parameter over a nonempty open window is NOT strictly monotone ===")
# schematic: gapped for Nf < Nfc ; conformal (gapless) on the open window [Nfc, 16.5].
# any order parameter O (gap, 0<O gapped, O=0 conformal) is 0 on the whole window:
Nfc = float(Q(153,19))          # schematic lower edge at the 2-loop value
def O(Nf):
    return max(0.0, (Nfc - Nf)) if Nf < Nfc else 0.0   # 0 on [Nfc, 16.5]
print(f"  schematic order parameter O(Nf) = max(0, Nfc-Nf), Nfc={Nfc:.4f}")
print(f"  O(0)={O(0):.4f}  O(8)={O(8):.4f}  O(12)={O(12):.4f}  O(15)={O(15):.4f}  O(16)={O(16):.4f}")
print("  O vanishes on the WHOLE window [Nfc,16.5] (an open interval) -> like the gap.")
# strict monotonicity would require O(Nf1) > O(Nf2) for ALL Nf1<Nf2:
viol = (O(12) > O(15))  # false, both 0
print(f"  strict-decreasing test O(12)>O(15): {viol}  (both 0 -> strict monotonicity VIOLATED)")
print("  The intended argument: O(Nf*)=0 (rigorous BZ) & O strictly decreasing => O(0)>0 (gapped).")
print("  But 'O strictly decreasing' is FALSE on any interval containing two window points.")
print("  Non-strict (O non-increasing) gives only O(0) >= 0  -- VACUOUS.")
print("  => STRICT monotonicity is REQUIRED for the conclusion, and is IMPOSSIBLE.")

print("\n=== 8. ACID TEST: Nf=0 (gapped) vs Nf=12 (conformal) at fixed N=3 ===")
for Nf in (0, 12):
    has, g2 = gstar2(Nf)
    phase = "GAPPED (evidence + SZZ strong coupling)" if Nf == 0 else "IR-CONFORMAL (lattice, contested)"
    if has:
        gm = 6 * CF * g2
        print(f"  Nf={Nf}: b0={float(b0(Nf)):+.2f} b1={float(b1(Nf)):+.2f} -> IRFP sqrt a={math.sqrt(float(g2)):.3f}, "
              f"gamma_m*={float(gm):.3f}; phase={phase}")
    else:
        print(f"  Nf={Nf}: b0={float(b0(Nf)):+.2f} b1={float(b1(Nf)):+.2f} -> NO IRFP (gamma_m* UNDEFINED); "
              f"phase={phase}")
print("  Any monotone must output: Nf=0 gapped, Nf=12 gapless. The ONLY scheme-independent")
print("  distinguisher is the PHASE QUALIFIER itself; every continuous observable that could carry")
print("  it is either an order parameter (vanishes on the window -> no-go above) or scheme-dependent.")
