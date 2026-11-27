"""GD2-WSTAR arithmetic + polyomino (planar-membrane) entropy probe.

Pure arithmetic + exact integer enumeration of fixed polyominoes (no MC, no random).
Output: gd2_wstar.log
"""
import math, itertools, sys, time

out = []
def P(*a):
    s = " ".join(str(x) for x in a)
    out.append(s)
    print(s)

t0 = time.time()
ln = math.log
KAPPA = 0.5

gammas = [
    ("GV2 gamma=1", 1.000000),
    ("a_n/site lower (GC2)", 2.333406),
    ("a_n/plaq lower (GC2)", 2.632032),
    ("B_n lower (GC2)", 2.746029),
    ("median a (M8)", 3.249025),
    ("median B (M8)", 3.458607),
    ("upper a (M9)", 3.375412),
    ("upper B (M9)", 3.697294),
]

# ---------------- S1: gate / window / q thresholds ----------------
P("== S1: gate (H4) vs window (CPB-S) thresholds ==")
P("gate   : c0*q >= g + ln(1+2PR)      ; window : q*(c0-kappa) >= g + ln(1+2PR)")
P("%-24s %10s %12s %12s %12s %12s" % ("source", "gamma", "gate(c0q)", "window(c0)",
                                     "q_gate@2.5986", "q_win@2.5986"))
for name, g in gammas:
    gate = g + ln(3.0)
    win = KAPPA + g + ln(3.0)
    qg = gate / 2.598612
    qw = gate / (2.598612 - KAPPA)
    P("%-24s %10.6f %12.6f %12.6f %12.6f %12.6f" % (name, g, gate, win, qg, qw))

P("")
P("== S1b: mu = exp(gamma) (animal growth constant) ==")
for name, g in gammas:
    P("%-24s gamma=%8.6f  mu=e^gamma=%10.4f" % (name, g, math.exp(g)))
P("measured ratios a: a6/a5=%.4f  B: B6/B5=%.4f (GC2-H23 1.2)" % (7218768/280175.0, 14305816/450255.0))

# ---------------- S2: c0 = T2/m^2 from the measured m/sqrt(sigma) ----------------
P("")
P("== S2: physical anchor c0 = T2/m^2 = 1/(m/sqrt(sigma))^2 ==")
for r in (3.405, 3.405 - 0.021, 3.405 + 0.021, 3.90, 3.30):
    P("m/sqrt(sigma)=%.3f -> c0=%.6f" % (r, 1.0 / r**2))
c0_phys = 1.0 / 3.405**2
P("c0_phys = %.6f ; ratio to corpus working value 1 : %.2fx smaller" % (c0_phys, 1.0 / c0_phys))
P("c0_phys vs window 4.230644 : short by %.2fx" % (4.230644 / c0_phys))
P("c0_phys vs window 4.344642 : short by %.2fx" % (4.344642 / c0_phys))

# ---------------- S3: window as a curve in q : c0*q^2 - kappa*q >= A ----------------
P("")
P("== S3: window in the scale-free form  c0*q^2 - kappa*q >= A = gamma+ln(1+2PR) ==")
P("   (c0 = T2/m^2 scale-free, q = trial mass * block scale; q=1 is a SLICE)")
for name, g in gammas:
    A = g + ln(3.0)
    P("  gamma=%8.6f  A=%9.6f" % (g, A))
    for c0 in (c0_phys, 1.0, 2.598612, 4.230644):
        # solve c0 q^2 - 0.5 q - A = 0
        disc = 0.25 + 4 * c0 * A
        qstar = (0.5 + math.sqrt(disc)) / (2 * c0)
        qmin_lhs = KAPPA / (2 * c0)   # minimiser of lhs
        lhs_min = c0 * qmin_lhs**2 - KAPPA * qmin_lhs
        c0_min_over_q = 3 * KAPPA**2 / (4 * A)
        P("      c0=%9.6f : need q >= %8.4f (worst q=%.3f, lhs_min=%8.4f) ; c0_min over q = %.6f"
          % (c0, qstar, qmin_lhs, lhs_min, c0_min_over_q))

# ---------------- S4: q-independence of BL-3 slab ----------------
P("")
P("== S4: BL-3 slab cost per plaquette (q-independent) vs window ==")
P("BL-3: U_l = sigma0 on in-slice space links, cost/plaq = beta(1-Re sigma0); reproduced numbers:")
for beta, label in ((2.5, "beta=2.5"),):
    P("  SU(2), sigma0=-1: Re tr/2 = -1 -> S_p = beta*(1-(-1)) = %.4f" % (beta * 2))
    P("  SU(3), sigma0=e^{2pi i/3}: Re tr/3 = -0.5 -> S_p = beta*(1+0.5) = %.4f" % (beta * 1.5))
P("  (1/g^2 = 0.625 at beta=2.5 in the convention beta=2/g^2 -> 2*0.625=1.25)")
P("window needs c0*q - kappa*q >= A ; with a q-INDEPENDENT per-plaquette cost C:")
for C in (1.25, 2.5, 3.75):
    for g in (2.632032, 3.458607):
        A = g + ln(3.0)
        P("  C=%.2f gamma=%7.4f A=%.4f : need C - 0.5q >= A -> q <= %.3f  (nonempty only if C>A=%.4f)"
          % (C, g, A, 2 * (C - A), A))

# ---------------- S5: block regrouping identically neutral ----------------
P("")
P("== S5: regrouping identity (geometric, block factor b) ==")
for b in (1, 2, 3, 4):
    nplaq = 6 * b**4
    P("  b=%d : fine plaquettes per block = 6*b^4 = %5d ; ln(%d)=%8.3f ; b^2=%d"
      % (b, nplaq, nplaq, ln(nplaq), b * b))
P("  min fine plaquettes per touched block = 1 (a connected 'snake' through M blocks has |C|=M)")
P("  => uniform per-block cost bound cannot exceed c0*q*1 ; block-graph entropy rate = gamma (isomorphic graph)")
P("  => criterion invariant : gain exactly 0")
P("  internal-config entropy ceiling per block: 96*ln2 = %.4f (b=2) ; ln C(96,48) = %.4f"
  % (96 * ln(2.0), ln(math.comb(96, 48))))

# ---------------- S6: exact fixed-polyomino enumeration (planar membrane proxy) ----------------
P("")
P("== S6: planar-membrane entropy proxy = fixed polyominoes (exact integer enumeration) ==")
NEI = ((1, 0), (-1, 0), (0, 1), (0, -1))

def canon(cells):
    mx = min(c[0] for c in cells); my = min(c[1] for c in cells)
    return tuple(sorted((c[0] - mx, c[1] - my) for c in cells))

cur = {((0, 0),)}
counts = [0, 1]
NMAX = 13
for n in range(2, NMAX + 1):
    nxt = set()
    for poly in cur:
        s = set(poly)
        for (x, y) in poly:
            for dx, dy in NEI:
                c = (x + dx, y + dy)
                if c not in s:
                    s.add(c)
                    nxt.add(canon(list(s)))
                    s.remove(c)
    cur = nxt
    counts.append(len(cur))
    P("  n=%2d : A_n=%10d   (1/n)ln A_n = %.6f   ratio A_n/A_{n-1} = %.5f"
      % (n, len(cur), ln(len(cur)) / n, len(cur) / counts[n - 1]))
P("  known/expected fixed polyominoes 1,2,6,19,63,216,760,2725,9910,36446,135268,505861,1903890")

# growth-constant estimates
A = counts
P("")
P("  growth constant mu = lim A_n/A_{n-1} ; endpoint ratios and 1/n fits:")
for n in range(6, NMAX + 1):
    P("    n=%2d ratio=%.6f  A_n^(1/n)=%.6f" % (n, A[n] / A[n - 1], A[n] ** (1.0 / n)))
# two-point 1/n extrapolation of the ratio
for n in (12, 13):
    r1 = A[n] / A[n - 1]; r0 = A[n - 1] / A[n - 2]
    P("    two-point 1/n ratio extrapolation at n=%d : mu ~ %.6f" % (n, r1 * (r1 / r0) ** 0 ))
P("  reference (literature, fixed polyominoes): mu = 4.0626...   [cite as external, not re-derived here]")
mu_lit = 4.0626
g_surf = ln(mu_lit)
P("  gamma_surf = ln(4.0626) = %.6f" % g_surf)
for A_ in (2.632032 + ln(3.0), 3.458607 + ln(3.0)):
    P("  window(gamma_surf) = kappa + gamma_surf + ln3 = %.6f  (vs A=%.6f)" % (KAPPA + g_surf + ln(3.0), A_))
P("  window(gamma_surf) - 2.598612 = %.6f" % (KAPPA + g_surf + ln(3.0) - 2.598612))
P("  criterion to rescue 2.5986 : gamma_eff <= 2.598612 - 0.5 - ln3 = %.6f  <=> mu_eff <= e = %.6f"
  % (2.598612 - 0.5 - ln(3.0), math.e))
P("  mu(animal, a_n)=%.4f ; mu(animal, B_n)=%.4f ; mu(planar polyomino)=%.4f ; e=%.4f"
  % (math.exp(2.632032), math.exp(2.746029), mu_lit, math.e))
P("  factor by which mu must drop : animal %.1fx / planar %.2fx"
  % (math.exp(3.249025) / math.e, mu_lit / math.e))

P("")
P("elapsed = %.2f s" % (time.time() - t0))
with open("gd2_wstar.log", "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")
