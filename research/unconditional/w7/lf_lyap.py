# -*- coding: utf-8 -*-
"""
lf_lyap.py -- witness for candidate (3) of the death-cause axis:
              non-autonomous / Floquet transfer operator with beta (or scale) as "time".

Question put to the numbers:
  (Q1) DIRECTION ACCOUNT FIRST: does the Lyapunov/dichotomy datum bound the gap from above
       or from below?
  (Q2) What does the beta-cocycle's exponent actually equal on the physical beta-path?
  (Q3) What replaces "eigenvalue" in the non-autonomous case, and what does a uniform
       lower bound then cost?
"""
import numpy as np

np.random.seed(20261125)
OUT = []


def say(s=""):
    print(s)
    OUT.append(str(s))


# ================================================================ W1  vacuum pins lambda_min
say("=" * 78)
say("W1  DIRECTION ACCOUNT.  T = e^{-aH}, spec(H) = {0} v [Delta, inf).  a=1.")
say("=" * 78)
Delta, M, a = 0.6444, 4.0, 1.0
spec = np.array([0.0, Delta, Delta + 0.7, 2.0, M])
T = np.exp(-a * spec)
say(f"  spec(H) = {np.round(spec,4)}")
say(f"  spec(T) = {np.round(T,6)}")
say(f"{'n':>6} {'||T^n||':>14} {'rate=-(1/(an))log||T^n||':>28}")
for n in (1, 5, 20, 100, 1000):
    nt = T.max() ** n              # norm = largest eigenvalue = 1 (vacuum)
    rate = -(1.0 / (a * n)) * np.log(nt)
    say(f"{n:6d} {nt:14.8f} {rate:28.3e}")
say("  -> the vacuum eigenvalue 1 pins the top of spec(T) at 1: the UNPROJECTED exponent is")
say("     identically 0.  Any 'lower bound on the gap' routed through it reads  Delta <= 0.")
say("     => without projecting out the vacuum, the Lyapunov datum is vacuous.")

# ================================================================ W2  projected: rate = gap
say()
say("=" * 78)
say("W2  after projecting the vacuum out: ||T^n P|| = e^{-a Delta n} EXACTLY")
say("    => the rate is not a bound on the gap, it IS the gap  (the D2 horn)")
say("=" * 78)
say(f"{'n':>6} {'||T^n P||':>16} {'rate':>12} {'Delta':>10}")
for n in (1, 5, 20, 100):
    nt = np.exp(-a * Delta * n)
    say(f"{n:6d} {nt:16.10f} {-(1/(a*n))*np.log(nt):12.8f} {Delta:10.6f}")
say("  DIRECTION: a *proven* decay bound  ||T^n P|| <= C e^{-a lam n} (all n)  implies")
say("             e^{-a Delta n} <= C e^{-a lam n}  =>  Delta >= lam - log C/(a n)  =>  Delta >= lam.")
say("  => the direction is CORRECT (a norm UPPER bound gives a gap LOWER bound).")
say("  ==> D1 (direction reversed) does NOT fire for this candidate.  This is the first")
say("      candidate since C3 (trace/Alon-Boppana) to pass the direction account.")
say("  But the hypothesis 'C e^{-a lam n} uniformly' is, by the displayed identity, equivalent")
say("  to  Delta >= lam  in the limit: proving it IS proving the gap  (D2).")

# ================================================================ W3  beta-cocycle
say()
say("=" * 78)
say("W3  beta as TIME: the exponent of the beta-cocycle on the physical path")
say("    asymptotic freedom in lattice units:  a*m(beta) = D0 exp(-c beta)  ->  0 as beta->inf")
say("=" * 78)
c = np.log(10.0) / 2.5          # gap drops by 10x between beta=0 and beta=2.5
D0 = 0.6444 * np.exp(c * 2.5)   # so that a*m(2.5) = 0.6444
say(f"  fitted c={c:.5f}, D0={D0:.5f}   (a*m(2.5)={D0*np.exp(-c*2.5):.4f})")
Dfun = lambda b: D0 * np.exp(-c * b)
b0 = 1.0 / 48.0
say(f"  a*m(1/48) = {Dfun(b0):.4f}   a*m(2.5) = {Dfun(2.5):.4f}   a*m(10) = {Dfun(10):.3e}")
say()
say(f"{'B (upper end of the beta-window)':>34} {'mean of a*m over [1/48,B]':>26}")
for B in (2.5, 3.0, 5.0, 10.0, 30.0, 100.0, 1e4):
    mean = (Dfun(b0) - Dfun(B)) / (c * (B - b0))
    say(f"{B:34.2f} {mean:26.6f}")
say("  -> the beta-average decays to 0 like 1/B: the Lyapunov/ergodic exponent of the")
say("     beta-cocycle is the beta-INFINITY value of the gap, and that value is 0 (free theory).")
say("     The physical path points AT the free end.  Wrong half of the axis.  (D3)")

# ================================================================ W4  window artefact
say()
say("=" * 78)
say("W4  the beta-average is a function of the WINDOW, not an invariant of the beta=2.5 theory")
say("    (numerically the window-average sweeps the whole admissible range)")
say("=" * 78)
target = Dfun(2.5)
say(f"  a*m(2.5) = {target:.6f}   a*m(1/48) = {Dfun(b0):.6f}")
say(f"{'B':>10} {'mean over [1/48,B]':>22} {'mean / a*m(2.5)':>18}")
for B in (1.0, 2.0, 2.5, 3.0, 5.0, 10.0, 30.0, 100.0, 1e4):
    mean = (Dfun(b0) - Dfun(B)) / (c * (B - b0))
    say(f"{B:10.2f} {mean:22.6f} {mean/target:18.4f}")
say("  -> by moving the window end the 'average datum' takes ANY value in (0, a*m(1/48)]:")
say("     it is a window artefact, not an invariant of the theory at beta=2.5.  The only")
say("     canonical window is the whole half-line, and there the average is 0 (W3).")
say("     => the average is a strictly different functional of the beta-family than the")
say("        point value  (D14); and the only canonical choice of it vanishes  (D3).")

# ================================================================ W5  non-uniform hyperbolicity
say()
say("=" * 78)
say("W5  Lyapunov existence  !=  uniform hyperbolicity (exponential dichotomy)")
say("    scalar cocycle Phi(n,0) = prod_{k<=n} e^{-Delta}(1+1/k) = e^{-Delta n}(n+1)")
say("=" * 78)
Dl = 0.7
say(f"{'n':>6} {'rate_n':>12} {'C_n=||Phi||e^{Delta n}':>24}")
for n in (1, 10, 100, 1000, 10000, 100000):
    rate = Dl - np.log(n + 1.0) / n          # = -(1/n) log( e^{-Dl n}(n+1) )
    Cn = float(n + 1)                        # exact: Phi*e^{Dl n} = n+1
    say(f"{n:6d} {rate:12.8f} {Cn:24.1f}")
say(f"  -> the exponent exists and equals -{Dl} (rate {Dl}: rate_n -> {Dl} as log(n+1)/n -> 0),")
say("     but C_n = n+1 -> inf: there is NO uniform C.  Lyapunov exponents are limits,")
say("     not uniform bounds.")
say("     A (L,a)-uniform lower bound needs the dichotomy CONDITION (uniform hyperbolicity),")
say("     which is strictly stronger -- and in the beta coordinate it is literally")
say("     'dist(spec_dichotomy, 0) >= a*lam', i.e. the gap restated in another coordinate.")

# ================================================================ W6  ergodic average
say()
say("=" * 78)
say("W6  and if beta is randomised, Birkhoff gives the ergodic MEAN, never a point value")
say("=" * 78)
for n in (2000, 20000):
    bet = np.random.uniform(b0, 2.5, n)
    d = Dfun(bet)
    say(f"  n={n:6d}:  mean={d.mean():.6f}   min={d.min():.6f}   max={d.max():.6f}"
        f"   (a*m(2.5)={Dfun(2.5):.6f})")
say("  mean >= min always, so a proven lower bound on the mean says NOTHING about the minimum,")
say("  and the point value at beta=2.5 is not an ergodic average of anything.")

say()
say("=" * 78)
say("SUMMARY (candidate 3 / non-autonomous, Floquet, Lyapunov)")
say("=" * 78)
say(" - DIRECTION: PASS.  A norm upper bound  ||T^n P|| <= C e^{-a lam n}  gives Delta >= lam.")
say("   D1 does not fire.  (Second candidate ever to pass the direction account, after C3.)")
say(" - unprojected the exponent is identically 0 (vacuum), i.e. vacuous (W1).")
say(" - projected, the rate IS the gap (W2): the required uniform decay bound is equivalent")
say("   to the conclusion -- the D2 horn.")
say(" - beta-as-time: the exponent is the beta-infinity value, which asymptotic freedom sends")
say("   to 0; the average decays like 1/B (W3) -- wrong half of the axis (D3).")
say(" - the beta-average is a WINDOW artefact, not an invariant: it sweeps (0, a*m(1/48)] as")
say("   the window end moves (W4); the only canonical window is the half-line, where it is 0.")
say("   Hence it is a strictly different functional of the beta-family than the point value (D14),")
say("   and its canonical value vanishes (D3).")
say(" - uniform hyperbolicity is strictly stronger than Lyapunov existence (W5) and is the")
say("   gap in the beta coordinate; ergodic means never control minima (W6).")
say(" - D12 (no scale separation) never becomes the main cause: it is not even reached,")
say("   because D3/D14/D2 fire first.  Prediction correction: D13 (carrier identically zero)")
say("   does NOT fire either -- the carriers here are non-degenerate.")
