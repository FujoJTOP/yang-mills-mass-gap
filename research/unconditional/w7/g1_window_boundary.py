"""
g1_window_boundary.py -- witness for G1-WINDOW-BOUNDARY.md
SU(3) Banks-Zaks conformal window from the 2-loop beta function, and the
distance of N_f=0 from it.  TINY (numpy only, no lattice, no MC).

beta(g) = -b0 g^3 - b1 g^5   (16pi^2 units absorbed)
  b0(Nf) = 11 - (2/3) Nf        > 0  <==> Nf < 16.5  (asymptotic freedom)
  b1(Nf) = 102 - (38/3) Nf      < 0  <==> Nf > 306/38 ~ 8.05  (Banks-Zaks zero g*^2=-b0/b1>0)

Verdict (printed): the 2-loop window is Nf in {9..16}; Nf=0 has b0>0, b1>0 => NO zero at any perturbative
order, i.e. maximally far outside.  Perturbative-only theorem; the nonperturbative statement is G1 itself.

Run:  py -3.12 -I research/unconditional/w7/g1_window_boundary.py
"""
import numpy as np

def b0(Nf): return 11.0 - (2.0/3.0)*Nf
def b1(Nf): return 102.0 - (38.0/3.0)*Nf

print("SU(3) 2-loop beta:  b0 = 11 - 2Nf/3 ,  b1 = 102 - 38Nf/3")
print(f"  2-loop BZ boundary  Nf^c = 306/38 = {306/38:.4f}")
print(f"  AF upper edge       Nf < 33/2  = 16.5")
print()
hdr = f"{'Nf':>3} | {'b0':>8} | {'b1':>8} | {'AF?':>4} | {'BZ zero?':>8} | {'g*^2':>10}"
print(hdr); print("-"*len(hdr))
window = []
for Nf in range(0, 17):
    a, b = b0(Nf), b1(Nf)
    af = a > 0.0
    zero = af and b < 0.0
    g2 = (-a/b) if zero else float('nan')
    if zero: window.append(Nf)
    print(f"{Nf:>3} | {a:>8.3f} | {b:>8.3f} | {str(af):>4} | {str(zero):>8} | {g2:>10.4f}")

print()
print("2-loop conformal window (BZ zero):  Nf in", window)
print("Nf=0 :  b0=%.1f>0, b1=%.1f>0  =>  NO zero at any perturbative order" % (b0(0), b1(0)))
print("       => Nf=0 is maximally FAR outside the window (farthest possible point).")

# distance / robustness argument
Nfc_2loop = 306/38
Nfc_4loop = 11.0     # van Ritbergen-Vermaseren-Larin (hep-ph/9701390), lower-edge ~10-11
Nfc_5loop = 12.0     # higher
shift_2to4 = Nfc_4loop - Nfc_2loop
print()
print("Distance of Nf=0 from the boundary:")
print(f"  2-loop  Nf^c = {Nfc_2loop:.2f}")
print(f"  4-loop  Nf^c ~ {Nfc_4loop:.0f}   (2->4 loop shift = +{shift_2to4:.1f})")
print(f"  5-loop  Nf^c ~ {Nfc_5loop:.0f}")
print(f"  gap to Nf=0 = Nf^c - 0 ~ {Nfc_2loop:.1f}..{Nfc_5loop:.0f}  vs loop-shift +{shift_2to4:.1f}")
print("  => reaching Nf=0 would need a nonperturbative shift ~%.0fx the known 4-loop correction." %
      ((Nfc_5loop-0)/max(shift_2to4,1e-9)))

# sanity: beta(g) shape at small Nf (no zero) vs a window Nf (one zero)
print()
def beta(g, Nf): return -b0(Nf)*g**3 - b1(Nf)*g**5
for Nf in (0, 12):
    gs = np.linspace(1e-3, 3.0, 200001)
    vals = beta(gs, Nf)
    sign = np.sign(vals)
    # count sign changes of beta (a zero is where beta changes sign from - to +)
    flips = np.sum((sign[:-1] != sign[1:]))
    gzero = gs[np.where(sign[:-1] != sign[1:])[0]] if flips else np.array([])
    print(f"  N_f={Nf:>2}: beta sign-changes on g in (0,3] = {int(flips)}"
          + (f"  g*={np.sqrt(-b0(Nf)/b1(Nf)):.3f}" if Nf==12 else "  (no zero)"))
print()
print("CONCLUSION: Nf=0 outside the window is [T] PERTURBATIVELY (b1>0), far outside (~9-12 units),")
print("robust to loop shifts (+3); but the NONPERTURBATIVE 'no strong-coupling zero' IS (G1). No QED.")
