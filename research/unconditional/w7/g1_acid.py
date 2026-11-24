"""
g1_acid.py -- ACID TEST of (G1): can a SINGLE-coupling AF theory have a beta-function zero?

Two-loop QCD beta (SU(N)):  beta(g) = -b0 g^3 - b1 g^5
  b0 = (11N - 2 Nf)/3
  b1 = 34 N^2/3 - (10N/3 + 2 Cf) Nf,   Cf = (N^2-1)/(2N)
A Banks-Zaks zero  g*^2 = -b0/b1 > 0  exists  iff  b0>0 and b1<0 (or both <0).

For FIXED Nf (one coupling) g* is a single number -> a zero is NOT forbidden by "single coupling".
Nf=0 (pure YM): b0,b1 both >0 -> NO zero.  Nf in the window: zero EXISTS.
=> the window is a real SINGLE-COUPLING AF gapless-conformal theory; "AF" / "single coupling"
   does NOT forbid a fixed point. (G1) for Nf=0 is a matter-content / dynamical fact, not kinematic.
Run:  py -3.12 research/unconditional/w7/g1_acid.py
"""
import math

def coeffs(N, Nf):
    b0 = (11 * N - 2 * Nf) / 3.0
    Cf = (N * N - 1) / (2.0 * N)
    b1 = 34 * N * N / 3.0 - (10 * N / 3.0 + 2 * Cf) * Nf
    return b0, b1

def bz_zero(N, Nf):
    b0, b1 = coeffs(N, Nf)
    if b1 == 0:
        return None, b0, b1
    g2 = -b0 / b1
    return (math.sqrt(g2) if g2 > 0 else None), b0, b1

print("SU(3) two-loop Banks-Zaks scan (one coupling per Nf):")
print(f"{'Nf':>6} {'b0':>9} {'b1':>10} {'g* (BZ zero)':>14}")
for Nf in [0, 2, 6, 8, 9, 10, 12, 14, 15, 16, 16.4]:
    z, b0, b1 = bz_zero(3, Nf)
    g = "-" if z is None else f"{z:.4f}"
    print(f"{Nf:>6} {b0:>9.4f} {b1:>10.4f} {g:>14}")

print()
print("CHECK acid claims:")
# (1) single coupling: Nf fixed -> one g, one beta -> a zero is just one condition (not forbidden)
z12, _, _ = bz_zero(3, 12)
z16, _, _ = bz_zero(3, 16)
print("  Nf=12 has a BZ zero:", z12 is not None, " g*=", None if z12 is None else round(z12, 4))
print("  Nf=16 has a BZ zero:", z16 is not None, " g*=", None if z16 is None else round(z16, 4))
# (2) Nf=0 has NO zero (2-loop): b0,b1 > 0
z0, b0_0, b1_0 = bz_zero(3, 0)
print("  Nf=0 has a zero:", z0 is not None, " (b0=%.2f, b1=%.2f)" % (b0_0, b1_0))
# (3) pure-glue/adjoint theories: same b0>0 for N_f=0 with matter in a rep R if (11N-2*T_R*Nf)>0
#     -> the exclusion of a zero at Nf=0 is a matter-content fact, not a single-coupling one.

# monotonic beta below g=0? show beta<0 for all g at Nf=0
def beta(g, N, Nf):
    b0, b1 = coeffs(N, Nf)
    return -b0 * g**3 - b1 * g**5

print()
print("Nf=0 beta(g) on g in (0,5]: sign =",
      "all negative (no zero)" if all(beta(g, 3, 0) < 0 for g in [i*0.05 for i in range(1, 101)]) else "has a zero")
print("Nf=12 beta(g) sign changes (zero present):",
      not all(beta(g, 3, 12) < 0 for g in [i*0.05 for i in range(1, 101)]))
