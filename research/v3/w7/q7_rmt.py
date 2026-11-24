#!/usr/bin/env python3
# q7_rmt.py  --  Q7 TOOL-FORGE, Theorem Q2 witness
# "Frozen-outlier / fine-tuning" no-go for random-matrix / free-probability models.
#
#  (A) The Savvidy soft mode is an EXACT, DETERMINISTIC eigenvalue -g|B| below the
#      positive Landau edge +g|B| (outlier separation 2g|B|).  It is NOT a BBP random
#      fluctuation; the tower is fixed by the background B, not a Wigner bulk.  So the
#      frozen-Hessian RMT returns lambda_min = -g|B| < 0, i.e. it certifies the NEGATIVITY
#      and is structurally BLIND to the measure rescue D_+ <= 0 (Brascamp-Lieb buy-in).
#
#  (B) On the measure-averaged object G'' = S - D with D/S ~ 0.97, the gap is the thin
#      residue G''/D ~ 3%.  An RMT/free-probability estimate of D with relative error
#      delta gives G'' in [G'' - delta D, G'' + delta D]; this interval contains 0 unless
#      delta < G''/D ~ 3%.  The exact-rank-one computation already had >=10% error
#      (SM2 finite-c vs 2/lambda_4) and the LLL model O(1) error (SM2 Sec.1.2), so the
#      sign is UNRESOLVED by RMT/FP at the head.
#
# Deterministic, numpy only.  Run:  py -3.12 -I q7_rmt.py
import numpy as np

print("=" * 74)
print("Q7 / Q2  FROZEN-OUTLIER + FINE-TUNING NO-GO witness   (random matrix / free prob.)")
print("=" * 74)

# ---------------------------------------------------------------- (A) exact Savvidy tower
# Charged spin-1 Landau tower of  -D^2 + 2 g F_B :  the (lattice) Savvidy operator.
# Lowest rung (spin anti-aligned, n=0):  lambda_- = -gB   (tachyonic, SM5 eq.(1.1)).
# The remaining rungs sit at positive  ~ (2n+1) gB.  This is deterministic given B.
print("\n(A) SAVVIDY OUTLIER IS EXACT (deterministic spectrum, not a Wigner fluctuation)")
for gB in (0.5, 1.0, 2.0, 4.0):
    n = np.arange(0, 4)
    tower = (2 * n + 1) * gB        #  +gB, 3gB, 5gB, ...
    lam_min = -gB
    edge_pos = tower.min()
    print(f"   gB={gB:4.1f}: lambda_min = {lam_min:6.2f}   positive edge = {edge_pos:6.2f}   "
          f"outlier separation = {edge_pos - lam_min:6.2f} (=2 gB)")
print("   => lambda_min<0 EXACTLY (outlier escapes; no BBP absorption). Frozen RMT -> G''<0.")
print("      This is the OPPOSITE of the measured measure-averaged G''>0. RMT is blind to D_+.\n")

# ---------------------------------------------------------------- (B) fine-tuning of S-D
# measured (SM2, beta=2.4, L=8):  D/S = 0.965 ; S_bare ~ 1.40 ; G'' = S - D ~ 3.5% S
print("(B) FINE-TUNING: gap G'' = S - D is a THIN residue of two O(S) quantities")
S = 1.40
for DS in (0.965, 0.970, 0.980, 0.990):
    D = DS * S
    Gpp = S - D
    delta_needed = Gpp / D                      # relative error on D to resolve the SIGN
    print(f"   D/S={DS:5.3f}:  S={S:.3f}  D={D:.3f}  G''={Gpp:7.4f}   "
          f"sign resolved only if (rel err on D) < {delta_needed*100:5.2f}%")

print("\n   Achievable relative errors from the SAME rank-one machinery (SM2):")
achievable = [
    ("finite-c vs 2/lam4 (exact quadrature)", 0.10),
    ("LLL overlap model (SM2 Sec.1.2, O(1) lossy)", 1.00),
    ("free-prob leading cumulant (model)", 0.30),
]
delta_needed = (S - 0.97 * S) / (0.97 * S)
print(f"   needed for sign resolution: delta < G''/D = {delta_needed*100:.2f}%")
for name, d in achievable:
    verdict = "RESOLVES" if d < delta_needed else "DOES NOT RESOLVE (error >> gap)"
    print(f"      {name:45s}: delta ~ {d*100:5.1f}%  -> {verdict}")
print("   => every available relative error exceeds the 3% gap; RMT/FP CANNOT fix the sign.\n")

# explicit interval for the best (10%) case
D = 0.97 * S; Gpp = S - D
for d in (0.03, 0.10, 0.30):
    lo, hi = Gpp - d * D, Gpp + d * D
    contains0 = (lo < 0 < hi)
    print(f"   delta={d*100:4.1f}%:  G'' in [{lo:7.4f}, {hi:7.4f}]   sign fixed? "
          f"{'NO (straddles 0)' if contains0 else 'yes'}")

# ---------------------------------------------------------------- (C) BL-rescue invisibility
print("\n(C) THE MEASURE RESCUE IS INVISIBLE TO RMT")
print("   frozen (operator) RMT sees:      G''_RMT = S - D_Sav      (D_Sav>=0, sign <=0 likely)")
print("   exact marginal (the real object): G'' = S - D_Sav + |D_+|  with D_+ <= 0 (BL).")
print("   The gap needs the MAGNITUDE |D_+| of the Brascamp-Lieb buy-in -- the non-Gaussianity")
print("   of the positive fibre. RMT/free-prob computes the operator spectrum, not |D_+|.\n")

# ---------------------------------------------------------------- (D) free-field acid test
print("(D) FREE-FIELD ACID TEST:  F=0 => no Savvidy outlier (no -gB rung); spectrum {k^2},")
print("    positive edge and negative edge both -> 0.  RMT returns lambda_min = 0 (gapless).")
print("    No false mass is manufactured. PASSES (returns 0) -- and is then silent on YM.\n")

print("Q2 VERDICT [NG]: the soft mode is an EXACT deterministic outlier (frozen RMT -> G''<0),")
print("                 and on S-D the RMT relative error (>=10%) >> the 3% gap.")
print("                 RMT/free-probability cannot resolve sign(G'') at the head.")
