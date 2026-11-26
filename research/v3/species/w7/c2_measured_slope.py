"""C2 — measured `d ln(m a)/d beta` from the corpus, vs the asymptotic 12*pi^2/(11 N^2).

Rules of this script:
  * every (beta, m a, err) is TRANSCRIBED from the corpus with `file:line`; nothing invented;
  * normalization declared per point: all betas are Wilson `beta_W = 2N/g^2`
    (SU(2): 4/g^2 = 2N/g^2; SU(3): 6/g^2 = 2N/g^2);
  * requirement (MAGNITUDE) 12*pi^2/(11 N^2): SU(2) 2.692, SU(3) 1.197.
    Sign: a ~ exp(-1/(2 b0 g^2)) => m a DECREASES with beta, so the signed
    requirement is -12*pi^2/(11 N^2). Ratios below are |measured| / requirement.

usage: py -3.12 c2_measured_slope.py  (writes nothing; stdout is tee'd to c2_measured_slope.log)
"""
import math

N_SERIES = {"SU(2)": 2, "SU(3)": 3}
B0 = lambda N: 11.0 * N / (48.0 * math.pi**2)          # 1/(2 b0 g^2) = 12 pi^2 beta/(11 N^2)
B1 = lambda N: 34.0 * N**2 / (3.0 * (16.0 * math.pi**2) ** 2)


def req(N):
    """asymptotic |d ln(ma)/d beta|"""
    return 12.0 * math.pi**2 / (11.0 * N * N)


# ---------------------------------------------------------------- corpus points [E]
# (series, label, beta, m a, err, source)
PTS = [
    ("SU(3)", "smear n=8, L=6, t=1-3", 6.00, 1.073, 0.088,
     "research/lattice-cpb/REPORT.md:357"),
    ("SU(3)", "smear n=8, L=24 pilot", 6.30, 0.886, 0.159,
     "research/lattice-cpb/REPORT.md:1060"),
    ("SU(2)/Ntau4", "GEVP t=2", 2.2, 1.1506, 0.0406,
     "research/unconditional/tm_gap_results.csv:2"),
    ("SU(2)/Ntau4", "GEVP t=2", 2.5, 1.2562, 0.1008,
     "research/unconditional/tm_gap_results.csv:3"),
    ("SU(2)/Ntau4", "GEVP t=2", 2.8, 1.3975, 0.0666,
     "research/unconditional/tm_gap_results.csv:4"),
    ("SU(2)/Ntau6", "GEVP t=3", 2.2, 1.0800, 0.6001,
     "research/unconditional/tm_gap_results.csv:5"),
    ("SU(2)/Ntau6", "GEVP t=3", 2.5, 0.8467, 0.1481,
     "research/unconditional/tm_gap_results.csv:6"),
    ("SU(2)/Ntau6", "GEVP t=3", 2.8, 0.7153, 0.2257,
     "research/unconditional/tm_gap_results.csv:7"),
]
GROUP_OF = {"SU(2)/Ntau4": "SU(2)", "SU(2)/Ntau6": "SU(2)", "SU(3)": "SU(3)"}


def slope(p1, p2):
    """ln((ma)2/(ma)1)/(beta2-beta1) with error propagation (errs independent)."""
    (s, l1, b1_, m1, e1, _), (_, l2, b2_, m2, e2, _) = p1, p2
    N = N_SERIES[GROUP_OF[s]]
    s_ = math.log(m2 / m1) / (b2_ - b1_)
    e_ = math.sqrt((e1 / m1) ** 2 + (e2 / m2) ** 2) / abs(b2_ - b1_)
    return s_, e_, req(N)


# ---------------------------------------------------------------- main table
print("=" * 100)
print("C2  measured  d ln(m a)/d beta  vs  asymptotic  12 pi^2/(11 N^2)   [beta = 2N/g^2 everywhere]")
print("=" * 100)
rows = []
for i in range(len(PTS) - 1):
    p1, p2 = PTS[i], PTS[i + 1]
    if p1[0] != p2[0]:
        continue
    s_, e_, r_ = slope(p1, p2)
    rows.append((GROUP_OF[p1[0]], p1[2], p2[2], s_, e_, r_, abs(s_) / r_, p1[0]))

hdr = "%-13s %-11s %10s %8s %8s %6s" % ("series", "beta1->b2", "slope", "err", "req", "ratio")
print(hdr); print("-" * len(hdr))
for g, b1_, b2_, s_, e_, r_, ra, ser in rows:
    lab = ser if "/" in ser else ser + " (L=6->L=24)"
    print("%-13s %4.2f->%-4.2f %10.4f %8.4f %8.4f %6.2f  [E]" % (lab, b1_, b2_, s_, e_, r_, ra))

# ---------------------------------------------------------------- systematic band (SU(3) only)
print("\n-- SU(3) operator/volume systematic on the SAME beta pair [E for inputs, C for the slope]")
alt = [("L=6 unsmeared n=0", 1.68, 0.41, "research/lattice-cpb/REPORT.md:356"),
       ("L=10 unsmeared xi/a=1.22 -> am", 1.0 / 1.22, 0.0, "research/lattice-cpb/REPORT.md:378")]
for lbl, m_alt, e_alt, src in alt:
    s_ = math.log(0.886 / m_alt) / 0.3
    print("  beta=6.0 endpoint taken as am=%.3f (%-34s %s) -> slope %+.4f, ratio %5.2f"
          % (m_alt, lbl, src.split(":")[-1], s_, abs(s_) / req(3)))

# ---------------------------------------------------------------- strong-coupling bound point
print("\n-- 2-sigma bound point beta<=3.3 (am > 4.4): only an inequality, no slope value [E]")
s_b = math.log(1.073 / 4.4) / (6.0 - 3.3)
print("  |slope| over [3.3,6.0] > %.4f  => ratio > %5.2f  (2 sigma)" % (abs(s_b), abs(s_b) / req(3)))

# ---------------------------------------------------------------- extrapolation to 90% [C]
print("\n-- extrapolation to 90% of the requirement, SU(3), N=3 [C -- crude model, formula below]")
b0, b1c = B0(3), B1(3)
B = b1c / (2 * b0**2)                    # slope drift of the 2-loop a^-1: S2(b)=1.197+B/b
S2 = lambda b: req(3) + B / b
b_bar, S_bar = 6.15, abs(rows[0][3])
A = b_bar * (1.0 - S_bar / S2(b_bar))    # ansatz S(b) = S2(b)*(1 - A/b)
print("  b0=%.6f b1=%.6f  B=b1/(2b0^2)=%.4f  S2(b)=1.1963+%.4f/b" % (b0, b1c, B, B))
print("  fit at beta_bar=%.2f: |S|=%.4f, S2=%.4f -> A=%.3f   (S(b)=S2(b)(1-A/b))"
      % (b_bar, S_bar, S2(b_bar), A))
target = 0.9 * req(3)
lo, hi = 6.3, 400.0
for _ in range(200):
    mid = 0.5 * (lo + hi)
    lo, hi = (mid, hi) if S2(mid) * (1 - A / mid) < target else (lo, mid)
b90 = 0.5 * (lo + hi)
print("  => S(beta90)=0.9*1.197=%.4f at beta90 = %.1f" % (target, b90))
c = b_bar * (1.0 - S_bar / req(3))       # simpler ansatz |S| = req*(1 - c/b)
print("  simplified ansatz |S|=1.197*(1-%.3f/b)  -> beta90 = %.1f" % (c, c / 0.1))
d = b_bar**2 * (1.0 - S_bar / req(3))    # alternative: deficit ~ 1/b^2 (O(g^4))
print("  alternative deficit ~1/b^2 (O(g^4))       -> beta90 = %.1f  [spread matters]" % math.sqrt(d / 0.1))
xi_slope = 1.0 / (12.0 * b0)             # d ln(xi/a)/d beta at 2 loops (== 1.1963)
print("  same point on the 2-loop trajectory: xi/a = 0.93*exp(%.4f*(%.1f-6)) = %.2e  [C]"
      % (xi_slope, b90, 0.93 * math.exp(xi_slope * (b90 - 6.0))))

print("\n-- sign check: requirement is NEGATIVE (a -> 0 as beta -> inf); SU(2)/Ntau4 rises =>")
print("   that series has the WRONG SIGN (finite-T thermal mass above deconfinement).")
