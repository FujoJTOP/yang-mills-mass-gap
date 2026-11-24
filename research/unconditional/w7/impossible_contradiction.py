"""
impossible_contradiction.py -- witness for IMPOSSIBLE-CONTRADICTION.md

Tests the EXACT inputs used to try to derive a contradiction from "pure YM is
gapless", WITHOUT any IR assumption:

  R1  beta-function: pure SU(N) YM at 2-loop has NO zero at g>0 (no Banks-Zaks
      window) -- only the UV fixed point g=0.  A gapless *CFT* IR would need
      beta(g*)=0 at g*>0.

  R2  a-theorem: a_UV = free Maxwell = 62(N^2-1)/360 > 0; along any flow
      a_IR < a_UV.  BOTH a_IR=0 (gapped) and 0<a_IR<a_UV (gapless CFT)
      satisfy it -> the theorem excludes neither.

  R3  in a scale-invariant vacuum <TrF^2>=0 (dimension-4 vev, no scale);
      the confined phase has <TrF^2> > 0.  No EXACT relation forbids either:
      "<TrF^2>=0 or >0" IS the gapped/gapless question.

Run:  py -3.12 research/unconditional/w7/impossible_contradiction.py
"""
import numpy as np

print("=" * 74)
print("R1  pure-SU(N) YM 2-loop beta-function: zeros at g>0 (BZ)?")
print("=" * 74)

def beta2(g, N):
    # dg/dln mu = -b0 g^3 - b1 g^5  (pure YM; b0,b1>0 for all N)
    b0 = 11.0 * N / (48 * np.pi**2)
    b1 = 34.0 * N**2 / (3 * (16 * np.pi**2)**2)
    return -b0 * g**3 - b1 * g**5

for N in (2, 3, 4):
    gs = np.linspace(1e-4, 8.0, 400001)
    b = beta2(gs, N)
    zero_at_gpos = np.any(np.diff(np.sign(b)) != 0)  # sign change => a zero
    print(f"  N={N}: b0={11*N/(48*np.pi**2):.5f} b1={34*N**2/(3*(16*np.pi**2)**2):.6f} "
          f"| beta<0 for all g>0? {np.all(b < 1e-14)} | zero at g>0? {zero_at_gpos}")
print("  => pure YM: NO Banks-Zaks zero at g>0 (2-loop).  A gapless CFT IR")
print("     would need beta(g*)=0, g*>0.  Contradiction ONLY IF 'gapless => CFT'.\n")

print("=" * 74)
print("R2  a-theorem: does it exclude the gapless branch?")
print("=" * 74)
for N in (2, 3):
    a_UV = 62.0 * (N**2 - 1) / 360.0
    print(f"  N={N}: a_UV={a_UV:.5f}")
    for name, a_IR in [("gapped  a_IR=0", 0.0), ("gapless a_IR=a_UV/2", a_UV/2)]:
        print(f"    {name:22s}: a_IR<a_UV (strict flow) AND a_IR>=0 ? "
              f"{0 <= a_IR < a_UV}")
print("  => BOTH branches satisfy the a-theorem.  It excludes NEITHER.\n")

print("=" * 74)
print("R3  <TrF^2>: scale-invariant (CFT) vs confined -- any EXACT obstruction?")
print("=" * 74)
print("  In a scale-invariant vacuum every dimension-4 vev vanishes:")
print("    CFT (gapless):  <TrF^2> = 0")
print("    confined:       <TrF^2> > 0  (gluon condensate, trace anomaly scale)")
print("  The trace anomaly Ward identity  <T^m_m> = (beta/2g)<TrF^2> + (a E4 - c W^2)")
print("  is an EXACT moment sum-rule; in flat space (E4=W^2=0) it reduces to")
print("  (beta/2g)<TrF^2> = 0.  In a CFT beta=0 makes this hold for ANY <TrF^2>,")
print("  including 0.  => no exact relation forbids <TrF^2>=0.")
print("  Choosing <TrF^2>=0 vs >0 (i.e. gapless vs gapped) is the open question.\n")

print("=" * 74)
print("VERDICT: no contradiction is derived. The sharpest route (R1) reduces to")
print("'gapless interacting 4D IR attractor => CFT (beta=0 at g*)'; that step is")
print("the exact gap. R2 (a-theorem) excludes neither branch. R3 is circular.")
print("=" * 74)
