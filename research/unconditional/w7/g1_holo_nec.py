"""
g1_holo_nec.py -- witness for G1-HOLO-NEC.

Claim: the bulk NEC (holographic c-theorem, A'' <= 0) does NOT distinguish an
IR AdS throat (gapless CFT) from a gapped wall. Both satisfy NEC.

Domain-wall ansatz: ds^2 = e^{2A(r)} eta_{mn} dx^m dx^n + dr^2,
convention: UV at r=+inf, IR at r=-inf.  NEC <=> A''(r) <= 0 (A concave).
c(r) = 1/A'(r)^(d-1)  is monotone decreasing toward the IR (Freedman-Gubser-
Pilch-Warner 1999); a-theorem: a_UV >= a_IR.

Case 1 (IR AdS throat, gapless CFT): A' runs 1/Luv (UV) -> 1/Lir (IR) with
        Lir < Luv (a_IR < a_UV). A' increases toward IR => A' decreasing in r
        => A'' <= 0 (NEC). Both endpoints AdS (A' > 0).
Case 2 (hard wall, gapped): A' = 1/Luv constant, space caps at r_wall
        (no smooth IR throat). A'' = 0 <= 0 (NEC).

Both NEC-satisfied => NEC cannot decide IR-CFT vs wall.

Run: py -3.12 -I research/unconditional/w7/g1_holo_nec.py
"""
import numpy as np

d = 4  # boundary dimensions (bulk d+1 = 5)


def cfun(Ap):
    return 1.0 / np.abs(Ap) ** (d - 1)


# ---------------- Case 1: UV-AdS -> IR-AdS (gapless CFT IR) ----------------
r = np.linspace(-6, 6, 40001)          # UV at +6, IR at -6
Luv, Lir = 1.0, 0.5                    # Lir < Luv  =>  a_IR < a_UV (a-theorem)
s = 1.0 / (1.0 + np.exp(-3.0 * r))     # s: 0 at r=-inf (IR), 1 at r=+inf (UV)
# A'(IR at r=-inf) = 1/Lir = 2.0 ; A'(UV at r=+inf) = 1/Luv = 1.0
Ap = (1.0 / Lir) - (1.0 / Lir - 1.0 / Luv) * s
A2 = np.gradient(Ap, r)                # A'' = dA'/dr
iUV, iIR = -50, 50                     # index near r=+6 (UV), r=-6 (IR)
print("Case 1  UV-AdS -> IR-AdS  (gapless conformal IR):")
print(f"  A'(UV at r=+6) = {Ap[iUV]:.4f}   A'(IR at r=-6) = {Ap[iIR]:.4f}   "
      f"(both > 0 => both AdS: {Ap[iUV] > 0 and Ap[iIR] > 0})")
print(f"  max A'' = dA'/dr = {A2.max():+.3e}   (<=0 => NEC satisfied: {A2.max() <= 1e-9})")
print(f"  c(UV) = {cfun(Ap[iUV]):.4f}   c(IR) = {cfun(Ap[iIR]):.4f}   "
      f"(c decreases toward IR: {cfun(Ap[iUV]) >= cfun(Ap[iIR])})")

# ---------------- Case 2: hard wall (gapped) ----------------
# A' = 1/Luv constant; the space ends at r_wall (a cap), no smooth IR throat.
Ap_wall = np.full_like(r, 1.0 / Luv)
A2_wall = np.gradient(Ap_wall, r)
print()
print("Case 2  hard wall (gapped):")
print(f"  A' = 1/Luv = {Ap_wall[0]:.4f} constant; space caps at r_wall (no IR throat)")
print(f"  max A'' = dA'/dr = {A2_wall.max():+.3e}   (<=0 => NEC satisfied: {A2_wall.max() <= 1e-9})")

print()
print("Both cases satisfy NEC (A'' <= 0).")
print("=> The bulk NEC / c-theorem does NOT exclude an IR AdS throat (gapless CFT).")
print("   Forcing the wall needs an extra input (a dilaton potential / 'no IR CFT') = (G1) = H1.")
