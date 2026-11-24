"""
star_counterexample.py -- ACID TEST of (star):
 "(star) a gapless interacting 4D AF Wilsonian IR attractor of an asymptotically-free
  gauge theory is conformal."

Two parts:
 (A) Counterexample search: is there a gapless, interacting, 4D, AF gauge theory that is
     NOT a CFT?  (walking / limit cycle / gapless-non-conformal phase / dilaton / theta-CP)
 (B) THE HONEST HOLE: (star) is TRUE (~Polchinski + a-theorem) but NOT SUFFICIENT -- the
     contradiction route then needs "pure YM has NO nonperturbative IR fixed point", which
     is NOT the 2-loop beta<0 fact; it is the gap itself (circular).

Run: py -3.12 -I research/unconditional/w7/star_counterexample.py
"""
import numpy as np

print("=" * 74)
print("(A) COUNTEREXAMPLE SEARCH: a clean gapless interacting 4D AF non-CFT?")
print("=" * 74)
print("""  (i)   walking          -> flows TO an approximate fixed point (a CFT): gapless+CFT,
                            NOT a counterexample (IR limit IS the CFT).
  (ii)  limit cycle      -> a periodic flow; EXCLUDED in unitary 4D by the a-theorem
                            (a flows monotonically, cannot cycle).
  (iii) gapless non-conformal phase -> if it has a scale it is NOT an IR attractor
                            (the flow does not converge to it); if the limit converges,
                            it is a fixed point = scale invariant = (by Polchinski) CFT.
  (iv)  dilaton          -> massless dilaton => scale invariance => CFT (flat direction);
                            massive/unstable dilaton => not gapless. No counterexample.
  (v)   theta / CP       -> theta=0 breaks nothing; theta=pi is a PHASE, not a gapless
                            non-conformal interacting AF theory.
  => (star) SURVIVES: no clean gapless interacting 4D AF non-CFT exists.""")

print("\n" + "=" * 74)
print("(B) THE HONEST HOLE: 2-loop beta<0 does NOT exclude a nonperturbative IR fixed point")
print("=" * 74)
# SU(2) 2-loop beta: beta = -b0 g^3 - b1 g^5, b0,b1>0 => beta<0 for all g>0, NO zero.
b0 = 22.0 / (3.0 * 16.0 * np.pi**2)
b1 = 136.0 / (3.0 * (16.0 * np.pi**2)**2)
beta2 = lambda g: -b0 * g**3 - b1 * g**5
gs = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0])
print("2-loop beta(g) =", np.round(beta2(gs), 4), " (all < 0, no zero)")
print("  => contradiction route: 'no BZ zero => no IR fixed point => gapped'")

# A nonperturbative IR-STABLE zero is NOT excluded by the 2-loop truncation:
# beta_full(g) = -b0 g^3 (1 - g^2/g*^2): 2-loop part -b0 g^3 (<0, no zero),
# but the FULL beta has an IR-stable zero at g* > 0.
gstar = 5.0
beta_full = lambda g: -b0 * g**3 * (1.0 - g**2 / gstar**2)
print(f"\nfull beta(g) = -b0 g^3 (1 - g^2/g*^2), g* = {gstar}")
print("  g<g*: beta<0 (flows UP to g* in the IR);  g>g*: beta>0 (flows DOWN to g*)")
print("  => g* is an IR-STABLE fixed point = a strongly-coupled IR CFT")
print("  => 2-loop beta<0 does NOT exclude it; the contradiction route ASSUMES 'no such g*'")
print("     = the gap itself (CIRCULAR).")
print("\n  IR-flow check (dg/d(-ln mu) = -beta points toward g*):")
for g in [1.0, 3.0, 8.0]:
    print(f"   g={g:4.1f}: beta_full={beta_full(g):+.5f}  (sign of -beta sets IR direction)")

print("\nVERDICT: (star) PASSES the counterexample search (it is TRUE, ~Polchinski scale+unitarity")
print("         => conformal, + a-theorem excluding limit cycles).")
print("         BUT it does NOT yield the gap: the contradiction route needs 'pure YM has no")
print("         nonperturbative IR fixed point', which is NOT the 2-loop beta<0 fact -- it is")
print("         the gap itself. So (star) is NECESSARY but CIRCULAR as a gap proof.")
