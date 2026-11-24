#!/usr/bin/env python3
"""
e2_second_scale_witness.py  --  [V-here] witness for E2-FLIP-SECOND-SCALE.

Purpose: instantiate the *scale-anchor dichotomy* (Lemma SA) arithmetically.
It is NOT new evidence about QCD.  It is a transparent machine-check of three
logical facts about a single-scale theory:

  (F1) For any candidate anchor Q = q*Lam that is itself proportional to the
       single scale Lam, the anchor inequality  m1 >= alpha*Q  is EQUIVALENT to
       the pure-number inequality  alpha <= m1/Q .  So the best alpha one could
       ever prove is the (unknown) pure ratio m1/Q -- i.e. the gap in disguise.

  (F2) The whole statement is scale-invariant: replacing Lam -> lam*Lam leaves
       every ratio, and hence the truth of  m1 >= alpha*Q  for fixed alpha,
       unchanged.  A scale-covariant anchor therefore cannot *produce* a floor;
       it can only *restate* the ratio.

  (F3) A scale Q whose dependence on Lam is dressed by a pure factor N^a is no
       exception: m1/(N^a Lam) = (m1/Lam)/N^a is still a pure number, so the
       N-dressing buys nothing (the candidate-(ii) N^a*Lam proposal).

Converse (why the hunt targets "external" scales): if Q_ext is KNOWN and NOT
proportional to Lam (an input we set), then m1 >= alpha*Q_ext IS an absolute
floor.  Lemma SA says the only chance is therefore an external scale -- and the
six-candidate audit (doc sections) shows each external scale fails uniformity,
direction, or connection-to-pure-YM.

All numbers below are [R] (standard/as-reported lattice-hadron values, NOT
re-verified here); the arithmetic is the point, not the values.
"""
import math

# ---- [R] reference values (GeV), pure SU(3) Yang-Mills, Nf=0, zero T ----
# single dynamical scale: use sqrt(sigma) as the scheme-independent anchor.
sqrt_sigma = 0.440     # Morningstar-Peardon scaling; ~0.44 GeV
m_0pp      = 1.710     # lightest 0++ glueball
m_2pp      = 2.400     # 2++ glueball
m_1mp      = 2.560     # 1-+ glueball (exotic)
chi14      = 0.180     # topological susceptibility^(1/4)
Tc         = 0.270     # deconfinement temperature
Lam_MSbar  = 0.250     # Lambda_MSbar^{Nf=0}; SCHEME-DEPENDENT label of same scale

m1 = m_0pp             # the gap (lightest glueball in the 0++ channel)

print("="*72)
print("(F1) 'proportional to Lambda' anchors reduce to a pure-number ratio")
print("="*72)
anchors = {
    "sqrt(sigma)": sqrt_sigma,
    "m_2++":       m_2pp,
    "m_1-+":       m_1mp,
    "chi_t^(1/4)": chi14,
    "T_c":         Tc,
    "Lambda_MSbar":Lam_MSbar,
}
print(f"  gap m1 = m_0++ = {m1:.3f} GeV ; reference scale sqrt(sigma)={sqrt_sigma:.3f}")
print(f"  {'anchor Q':<14}{'Q':>8}{'Q/sqrt(sigma)':>16}{'m1/Q (= best alpha)':>22}")
for name, Q in anchors.items():
    print(f"  {name:<14}{Q:>8.3f}{Q/sqrt_sigma:>16.3f}{m1/Q:>22.3f}")
print("  -> every 'best alpha' is a pure number = the unknown ratio; proving")
print("     m1 >= alpha*Q is proving alpha <= m1/Q  = the gap.  No absolute floor.")
print("     (The ratio happens to be O(1)-nice only because it was MEASURED;")
print("      no theorem fixes it, and a proven lower bound on it IS the gap.)")

print()
print("="*72)
print("(F2) scale invariance: Lam -> lam*Lam with fixed alpha changes nothing")
print("="*72)
alpha = 0.10
for Q,name in [(sqrt_sigma,"sqrt(sigma)"),(m_2pp,"m_2++"),(chi14,"chi_t^(1/4)")]:
    print(f"  anchor Q={name}:")
    for lam in (0.1, 1.0, 10.0):
        m1_l, Q_l = m1*lam, Q*lam
        ok = (m1_l >= alpha*Q_l)
        # the *reason* ok holds is purely the ratio, independent of lam:
        print(f"     lam={lam:>4}: m1={m1_l:7.3f} Q={Q_l:7.3f}  "
              f"'m1>={alpha} Q' -> {ok}   (ratio m1/Q={m1_l/Q_l:.3f} unchanged)")

print()
print("="*72)
print("(F3) N-dressed scale N^a*Lam is still proportional to Lam  (candidate ii)")
print("="*72)
for a in (1, 2):
    for N in (3, 10):
        Q = (N**a) * sqrt_sigma
        print(f"  N={N:>2} a={a}: Q=N^a*sqrt(sigma)={Q:8.3f}  m1/Q={m1/Q:7.3f}"
              f"  (= (m1/sqrt(sigma))/N^a = {m1/sqrt_sigma/N**a:.3f})")
print("  -> N^a is a pure factor; the ratio is still unknown.  No anchor.")

print()
print("="*72)
print("CONVERSE: a KNOWN external scale (NOT proportional to Lam) does give a floor")
print("="*72)
Q_ext = 1.0   # imagine a known dimensionful input we fix by hand (units: GeV)
for alpha in (0.05, 0.20):
    print(f"  if Q_ext={Q_ext} is KNOWN & fixed (not proportional to Lam):"
          f"m1 >= {alpha}*Q_ext = {alpha*Q_ext:.3f} GeV  -> a genuine floor.")
print("  -> This is the ONLY winning form; the hunt's six candidates are exactly")
print("     the external scales, and each fails (uniformity / direction / connection).")

print()
print("="*72)
print("ACID TEST: free massless field has no gap (m1=0); every anchor is vacuous")
print("="*72)
m1_free = 0.0
for name,Q in anchors.items():
    bound = m1_free / Q if Q else float('nan')   # 'best alpha' -> 0
    print(f"  {name:<14}: best alpha = m1_free/Q = {bound:.3f}  "
          f"-> 'm1 >= 0' : correct number, certifies nothing (vacuous).")
print("  -> consistent with the charter acid test: a scale-covariant anchor is")
print("     silent (0) on the free field, i.e. it is not the YM-specific lever.")

print()
print("NOTE: these are witnesses of the dimensional lemma (F1-F3), NOT a proof")
print("      about QCD.  The structural failures of external scales are argued in")
print("      the document (volume: backward monotonicity; deformations: different")
print("      theory / degenerate limit; T>0: not Clay).  Clay stays open.")
