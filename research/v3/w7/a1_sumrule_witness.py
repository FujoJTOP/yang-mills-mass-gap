#!/usr/bin/env python3
"""
A1-SUMRULE-A numerical witness (numpy only, single job).

This is a *bookkeeping sanity check* of the algebra in research/v3/A1-SUMRULE-A.md,
NOT evidence for or against Yang-Mills. It confirms three elementary facts behind
the [NG] verdict on the monotone-to-residue bridge Delta_a vs f1^2:

  (1) a_UV^{(SU(N))} = (N^2-1)*31/180.  For N>=3, a_UV>1, so in the *project*
      normalization (f1^2 in [0, M0], M0=1) the inequality "Delta_a >= c f1^2"
      is WEAKER than the trivial cap f1^2 <= M0 for every c <= a_UV (useless).
      For SU(2), a_UV=0.517<1 so the cap can bite for c>a_UV, but it is still an
      UPPER bound and still permits f1^2=0 (no atom certified). Verdict unchanged.
  (2) KS sum rule  Delta_a = (f^4/pi) * Int_0^inf Im A(s)/s^3 ds  is HOMOGENEOUS
      in the dilaton decay constant f (which measures f1^2). With Delta_a FIXED and
      a positive toy density Im A = A0 s^3 exp(-s/M^2)  (so Int/s^3 = A0 M^2),
      the solved f = (pi Delta_a/(A0 M^2))^{1/4} ranges over orders of magnitude:
      Delta_a does not bound f  (<=> f1^2).  Only the SIGN (Im A>0) is content.
  (3) f1^2 = 0 is compatible with Delta_a = a_UV > 0 (the free field / gapped
      trivial endpoint). Hence NO constant c>0 in "Delta_a >= c f1^2" can certify
      a positive atom: the inequality returns 0 >= 0 at f1^2 = 0.
"""
import numpy as np

def aUV(N):
    """UV Euler anomaly for pure SU(N) in the standard normalization (per vector 31/180)."""
    return (N**2 - 1) * 31.0 / 180.0

def main():
    print("=" * 68)
    print("(1) a_UV(SU(N)) = (N^2-1)*31/180,  trivial cap M0 = 1")
    print("=" * 68)
    print(f"{'N':>3} {'a_UV':>10} {'a_UV/M0':>10}  note")
    for N in range(2, 7):
        a = aUV(N)
        if a > 1:
            note = "a_UV>1 => (S1) cap f1^2<=a_UV/c weaker than f1^2<=1 for c<=a_UV"
        else:
            note = "a_UV<1 => cap can bite for c>a_UV, but STILL only an upper bound; f1^2=0 allowed"
        print(f"{N:>3} {a:>10.5f} {a/1.0:>10.5f}  {note}")
    print("  => for SU(3): a_UV = %.5f. The (S1)-implied cap f1^2 <= a_UV/c exceeds the" % aUV(3))
    print("     trivial cap f1^2 <= M0=1 for every c <= %.3f; it never gives a FLOOR.\n" % aUV(3))

    print("=" * 68)
    print("(2) KS sum rule homogeneity in f  (Delta_a fixed, vary spectral norm A0,M)")
    print("=" * 68)
    Delta_a = aUV(3)                      # fix the anomaly drop (gapped pure YM value)
    A0s = np.array([1e-3, 1e-2, 1e-1, 1.0, 10.0])
    Ms  = np.array([0.5, 1.0, 2.0, 4.0])  # positive, physical toy scales
    fs = []
    for A0 in A0s:
        for M in Ms:
            integral = A0 * M**2          # Int_0^inf A0 s^3 e^{-s/M^2} / s^3 ds = A0 M^2
            f = (np.pi * Delta_a / integral) ** 0.25   # KS: Delta_a = (f^4/pi)*integral
            fs.append(f)
    fs = np.array(fs)
    print(f"  Delta_a fixed = {Delta_a:.5f}")
    print(f"  solved f ranges over  [{fs.min():.4g}, {fs.max():.4g}]  "
          f"(ratio {fs.max()/fs.min():.1e}) for the SAME Delta_a.")
    # perfect log anticorrelation between spectral norm and solved f => no constraint
    logs_int = np.log([A0 * M**2 for A0 in A0s for M in Ms])
    corr = np.corrcoef(logs_int, np.log(fs))[0, 1]
    print(f"  corr( log(A0 M^2), log f ) = {corr:+.4f}   (=> Delta_a sets f only via the")
    print("  spectral normalization; it does NOT bound f, i.e. does NOT bound f1^2).")
    print("  Only the SIGN in 'Im A >= 0' is physical content (Delta_a > 0).\n")

    print("=" * 68)
    print("(3) f1^2 = 0 is compatible with Delta_a = a_UV > 0  => no atom certified")
    print("=" * 68)
    print(f"  gapped trivial endpoint: f1^2 = 0, Delta_a = a_UV(SU3) = {aUV(3):.5f} > 0.")
    print("  'Delta_a >= c f1^2' at this point reads: %.5f >= 0. True for all c>0," % aUV(3))
    print("  but produces NO lower bound on f1^2 (0 is allowed). The inequality is")
    print("  one-sided (upper bound on f1^2), exactly the wrong direction for a residue FLOOR.")
    print("\nVERDICT consistent with A1-SUMRULE-A.md: bridge [NG], not [T]/[C].")

if __name__ == "__main__":
    main()
