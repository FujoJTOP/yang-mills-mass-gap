#!/usr/bin/env python3
r"""
v05_i2_toy.py -- I2 (Lambda* > m1^2) on EXACTLY SOLVABLE models: existence test.

DICTIONARY (R2-ISOLATION.md, T6-BINDING-ADVERSARIAL.md):
    C(t) = int e^{-s t} dnu(s),  dnu >= 0   (Kallen-Lehmann measure of the 0++ channel)
    m1^2 := inf supp nu ;  f1^2 := nu({m1^2}) ;  Lambda* := inf supp(nu - f1^2 delta_{m1^2})
    (I2)  Lambda* > m1^2      [s = the same coordinate as m^2; only RATIOS are convention-free]
NOTE ON UNITS: all models give Lambda*/m1^2 = a pure Casimir ratio, so the verdict is unit-free.
If the KL variable is instead E^2 (s = m^4), every ratio squares (8/3 -> 64/9) and no sign changes.

MODELS
  [M1] single plaquette (SU(2)): 0-dimensional matrix integral -- two readings:
       (a) as a model: no transfer direction => m1, Lambda* UNDEFINED (T03-TOY discipline)
       (b) on the only spectrum present (Laplace/Peter-Weyl tower of -Delta_G):
           m1^2 = C2(1/2) = 3/4,  Lambda* = C2(1) = 2  => I2 TRUE, margin 5/4 [units g^2/2]
  [M2] 2D YM (Migdal, exact): E_R = (g^2/2) C2(R); 0++ tower {C2(R)}.
       O = chi_{1/2} + chi_1 fuses from the vacuum into R=1/2 and R=1 => two atoms:
       m1^2 = 3g^2/8, Lambda* = g^2, margin 5g^2/8, rho^2 = 8/3 EXACT  => I2 TRUE
       (single-mode operator O = chi_{1/2}: nu is ONE atom => inf supp(empty) = +inf, vacuous)
  [M3] N=infinity single plaquette: C2(f)=(N^2-1)/(2N), C2(adj)=N => rho^2 = 2N^2/(N^2-1) -> 2
       margin (N^2+1)/(2N) -> infinity  => I2 TRUE, ratio -> 2 (NOT -> 1)
  [M4] M6 (oscillator (x) free, m=0): excitation spectrum {n w} (+) [0,inf) = [0,inf)
       => m1 = 0, Lambda* = 0, margin 0  => I2 FALSE.  Box regularization shows the margin
       is the box level spacing 6 pi^2/L^2 -> 0; a massive free module (mu) still fails (margin 0).
  [X]  abstract witness: nu = 0.6 delta_4 + 1_{(4,inf)} => f1^2 > 0 (I1 TRUE) but Lambda*=4=m1^2
       => having an atom at the bottom does NOT give I2.
Run: py -3.12 w7/v05_i2_toy.py   (writes w7/v05_i2_toy.log)
Labels [T]/[C].  No QED.  Clay 0%.
"""
import io
import sys
from fractions import Fraction as F
import numpy as np

LOG = io.StringIO()
TICK = {'T': 0, 'C': 0}


def say(*a):
    line = " ".join(str(x) for x in a)
    print(line, flush=True)
    LOG.write(line + "\n")


def frac(x):
    return f"{x} = {float(x):.6f}"


def row(name, m1sq, lam, i2, lab='T'):
    """one model row of the verdict table"""
    d = lam - m1sq
    say(f"  {name:<34} m1^2={float(m1sq):<10.6f} L*={float(lam):<10.6f} "
        f"L*-m1^2={float(d):<10.6f} I2={'TRUE ' if d > 0 else 'FALSE'}   [{lab}]")
    TICK[lab] = TICK.get(lab, 0) + 1
    return d


def main():
    say("=" * 100)
    say("[V05] I2  (Lambda* > m1^2)  on solvable models -- existence test")
    say("  nu = KL measure of 0++ channel; m1^2=inf supp nu; Lambda*=inf supp(nu - f1^2 delta)")
    say("=" * 100)

    # ------------------------------------------------------------------ [M1] single plaquette
    say("\n[M1] SINGLE PLAQUETTE (SU(2), 0-dimensional matrix integral)")
    say("  (a) as a model: 0 spacetime directions => no transfer => m1, Lambda* UNDEFINED  [T]")
    say("      (T03-TOY: writing 'undefined' as 0 or >0 fabricates a verdict in either direction)")
    say("  (b) on the only spectrum present -- Laplace (-Delta_G) / Peter-Weyl tower, 0++ = class fns:")
    c2 = lambda j: F(j) * (F(j) + 1)          # SU(2) Casimir, j = 0,1/2,1,...
    js = [F(1, 2), F(1), F(3, 2), F(2), F(5, 2)]
    say("      j      : " + "".join(f"{str(j):>10}" for j in js))
    say("      C2(j)  : " + "".join(f"{str(c2(j)):>10}" for j in js))
    gaps = [c2(js[i + 1]) - c2(js[i]) for i in range(len(js) - 1)]
    say("      gaps   : " + "".join(f"{str(g):>10}" for g in gaps) + "   (= j+1, strictly increasing)")
    m1 = c2(F(1, 2))
    lam = c2(F(1))
    say(f"      m1^2 = C2(1/2) = {frac(m1)}   Lambda* = C2(1) = {frac(lam)}   rho^2 = {lam/m1} = {float(lam/m1):.6f}")
    row("M1(b) plaquette, Laplace tower", m1, lam, True)
    say("      => I2 TRUE, but the axis is compactness (C2 strictly increasing, no accumulation),")
    say("         NOT dynamics: this IS 2D YM's leading order, not an independent model. [C]")
    TICK['C'] = TICK.get('C', 0) + 2       # two [C] caveats: M1(b) no-dynamics, M3 arc-not-spectrum
    say("  (c) 'one-mode' reading (a single level, so Lambda*=m1): that is the same statement as")
    say("      'there is no spectrum' => collapses into (a), not a counterexample. [T]")

    # ------------------------------------------------------------------ [M2] 2D YM
    say("\n[M2] 2D YM (Migdal exact): E_R = (g^2/2) C2(R), 0++ tower {C2(R)}")
    g2 = F(1)                                  # work in units g^2 = 1 (= 8/3 ratio is unit-free)
    E = lambda j: g2 / 2 * c2(j)
    e1, e2 = E(F(1, 2)), E(F(1))
    say(f"      E(1/2) = {frac(e1)}   E(1) = {frac(e2)}")
    say("      operator O = chi_{1/2}+chi_1 fuses from the vacuum to R=1/2,1 => nu has 2 atoms")
    say(f"      (both residues nonzero since <R|chi_R'|0> != 0 for R' = R) => f1^2 = |c_1|^2 > 0")
    row("M2 2D YM SU(2)", e1, e2, True)
    say(f"      rho^2 = Lambda*/m1^2 = C2(1)/C2(1/2) = 8/3 = {8/3:.6f} EXACT  (R2 W2, reconfirmed)")
    say("      single-mode operator O = chi_{1/2}: nu = |c|^2 delta_{E(1/2)} only, so")
    say("      supp(nu - f1^2 delta) = EMPTY => Lambda* = +inf by the inf-empty convention:")
    say("      I2 holds vacuously. VOID row, not evidence. [T]")
    c2s = lambda N: (N * F(N) - 1) / (2 * N)
    say(f"      SU(3) check: C2(f)={c2s(3)} C2(adj)=3 => rho^2 = {3/c2s(3)} = {float(3/c2s(3)):.6f}")

    # ------------------------------------------------------------------ [M3] N=infinity
    say("\n[M3] N=infinity SINGLE PLAQUETTE (Casimir tower at large N; 0++ = characters chi_R)")
    say(f"      {'N':>8}{'C2(f)':>14}{'C2(adj)':>12}{'rho^2=C2adj/C2f':>18}{'margin=C2adj-C2f':>19}")
    for N in (2, 3, 4, 10, 100, 1000):
        cf, ca = (F(N * N - 1, 2 * N)), F(N)
        say(f"      {N:>8}{float(cf):>14.6f}{float(ca):>12.6f}{float(ca / cf):>18.6f}{float(ca - cf):>19.6f}")
    say(f"      N=inf : rho^2 = 2N^2/(N^2-1) -> 2  (from above);  margin (N^2+1)/(2N) -> inf")
    row("M3 N=inf plaquette (limit)", F(1, 2), F(1), True)
    say("      => I2 TRUE at N=inf with ratio EXACTLY 2 (not 1): 1/N corrections do not kill it.")
    say("      Caveat: N=inf here is the Casimir tower, not the master field's eigenvalue arc")
    say("      (the arc is a continuum but it is an eigenvalue DENSITY, not the mass spectrum). [C]")

    # ------------------------------------------------------------------ [M4] M6
    say("\n[M4] M6 = oscillator (x) free module,  m = 0   (T03-TOY M6: p_a^2/2 + p_b^2/2 + w^2 a^2/2)")
    w = 2.0
    say(f"      excitation spectrum = {{n*w : n>=0}} (+) [0,inf) = [0,inf)   (w = {w})")
    say("      => m1 = 0 (bottom is the free continuum, no atom: f1^2 = 0)")
    say("      => Lambda* = inf supp(nu) = 0 = m1^2   margin 0      I2 FALSE")
    row("M4 M6 osc(x)free, m=0", F(0), F(0), False, lab='T')
    say("\n     box regularization (free direction in a ring of length L, k=2 pi n/L):")
    say(f"      {'L':>8}{'m1^2':>12}{'Lambda*':>12}{'margin':>12}{'6pi^2/L^2':>12}   I2")
    for L in (0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0, 300.0):
        lev = sorted({w * n + 2 * np.pi ** 2 * mm * mm / (L * L)
                      for n in range(0, 24) for mm in range(1, 60)})
        lev = [x for x in lev if x > 1e-12][:2]
        a, b = lev[0], lev[1]
        say(f"      {L:>8.2f}{a:>12.6f}{b:>12.6f}{b-a:>12.6f}{6*np.pi**2/L**2:>12.6f}"
            f"   {'TRUE' if b-a > 1e-9 else 'FALSE'}")
    say("      => margin = box level spacing ~ 6 pi^2/L^2 -> 0: the failure is the continuum at")
    say("         the bottom, not the size of m1. (Small box: I2 revives because the cut is gapped.)")
    say("      massive free module (mu>0): spectrum {n*w} (+) [mu^2,inf) => m1^2 = Lambda* = mu^2,")
    say("      margin 0 for EVERY mu: a gapped free module still fails I2. [T]")
    row("M4' massive free module", F(4), F(4), False, lab='T')

    # ------------------------------------------------------------------ [X] witness
    say("\n[X] abstract witness nu = 0.6 delta_4 + 1_{(4,inf)} e^{-(s-4)} ds  (R2 Lemma B / T6-4)")
    say("      f1^2 = 0.6 > 0 (I1 TRUE), m1^2 = 4, Lambda* = 4 => margin 0 => I2 FALSE")
    row("X  atom touching its cut", F(4), F(4), False, lab='T')

    # ------------------------------------------------------------------ verdict
    say("\n" + "=" * 100)
    say("[V] VERDICT TABLE")
    say("  model                          Lambda*      m1^2        Lambda*-m1^2   I2")
    say("  M1(a) plaquette as a model     UNDEFINED    UNDEFINED   -              (c) undefined")
    say("  M1(b) plaquette, Laplace tower 2            3/4         5/4            TRUE  (no dynamics)")
    say("  M2  2D YM SU(2)                g^2          3g^2/8      5g^2/8         TRUE  (exact)")
    say("  M3  N=inf plaquette            2m1^2        (N^2-1)/2N  (N^2+1)/2N     TRUE  (ratio->2)")
    say("  M4  M6 osc(x)free, m=0         0            0           0              FALSE")
    say("  M4' massive free module        mu^2         mu^2        0              FALSE")
    say("  X   atom on its cut            4            4           0              FALSE")
    say("  => (a) TRIGGERED: I2 is TRUE on 3 solvable models (2D YM and the Casimir towers),")
    say("     so I2 is not self-contradictory. (b) NOT triggered. (c) M1(a): undefined.")
    say("  => BUT every TRUE row is compactness of the gauge group (C2 strictly increasing,")
    say("     purely discrete tower, NO multi-particle cut in the 0++ channel). M6 and X show the")
    say("     margin is exactly 0, not small, when the bottom is a cut. So the toys do not test")
    say("     the 4D question, which is precisely whether a cut sits at m_g: existence only.")
    say("=" * 100)
    say(f"labels: [T]={TICK.get('T', 0)}  [C]={TICK.get('C', 0)}   Clay 0%.  No QED.")


if __name__ == '__main__':
    main()
    sys.stdout.flush()
    open(r"D:\Dev\math\w7\v05_i2_toy.log", "w", encoding="utf-8").write(LOG.getvalue())
