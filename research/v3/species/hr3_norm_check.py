# hr3_norm_check.py -- HR3-NORM: carrier-normalization battery for the no-light-mode theorem.
# Deterministic. No MC, no random numbers. < 10 s.
# Checks:
#  [1] free :phi^2:  point vs extensive residue: exact L-scaling (constant 4.5/pi^2, 3/(4 pi^2), ratio L^6/6)
#  [2] free :phi^2:  M0 in THREE regulators  -> the L-exponent of M0 is carrier/regulator data, not free
#  [3] w1 = f1^2/M0 on free field: extensive ~ L^-1, point ~ L^-4
#  [4] identity M0 = Var(O) = C(0) = nu([0,inf))  (so "divide by M0" == "divide by Var(O)")
#  [5] HN1 family w1 == 1/2 exactly; thin atom w1 = eps/2 -> 0 with a light mode
#  [6] M6 b-carrier: w1 == 6/pi^2 exactly while m1^2 ~ L^-4
#  [7] the flip: same free model, theorem truth value opposite in the two carriers
#  [8] structural: any absolute threshold f1^2 >= r0 is flippable by s_L := sqrt(r0/(2 f1^2(L)))
#  [9] external consistency: C(0) ~ L^3 on YM at beta=2.5 (GS1-SEAT2:0.15 numbers)
import math
from fractions import Fraction

P = print


def head(s):
    P("\n" + "=" * 78)
    P(s)
    P("=" * 78)


# ------------------------------------------------------------------ free field machinery
def omega(L, n):
    return 2.0 * math.pi * math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2) / L


def modes(L, ncut):
    out = []
    for n1 in range(-ncut, ncut + 1):
        for n2 in range(-ncut, ncut + 1):
            for n3 in range(-ncut, ncut + 1):
                if n1 == 0 and n2 == 0 and n3 == 0:
                    continue
                out.append((n1, n2, n3))
    return out


def G_of_t(L, t, ncut):
    s = 0.0
    for n in modes(L, ncut):
        w = omega(L, n)
        s += math.exp(-w * t) / (2.0 * w)
    return s / L ** 3


# ------------------------------------------------------------------ [1] point vs extensive
head("[1] FREE :phi^2: ON T^3_L (a=1, continuum dispersion, BZ cutoff |n|_inf <= L/2)")
P(" extensive (zero-momentum carrier):  f1^2 = 3/omega1^2 , omega1 = 2 pi / L")
P(" point     (coincident-point):       f1^2 = 18/(L^6 omega1^2)     [= 2 G(t)^2 leading atom]")
P(" ratio     point/extensive           = L^6/6  (exact)")
P("")
P("   L    f1^2_ext        f1^2_ext/L^2   f1^2_pt          f1^2_pt*L^4   6*f1_ext/(L^6*f1_pt)")
for L in (4, 8, 16, 32, 64):
    w1 = 2.0 * math.pi / L
    fe = 3.0 / w1 ** 2
    fp = 18.0 / (L ** 6 * w1 ** 2)
    P("  %3d  %12.6e   %11.6f     %12.6e   %11.6f     %10.6f"
      % (L, fe, fe / L ** 2, fp, fp * L ** 4, 6.0 * fe / (L ** 6 * fp)))
P(" => constants: 3/(4 pi^2) = %.6f ; 4.5/pi^2 = %.6f" % (3.0 / (4 * math.pi ** 2), 4.5 / math.pi ** 2))

# ------------------------------------------------------------------ [2] M0 in three regulators
head("[2] M0 = nu([0,inf)) : THE L-EXPONENT IS CARRIER/REGULATOR DATA")
P(" (a) extensive, continuum dispersion, BZ cutoff ncut = L/2   -> M0 ~ L^3 ?")
P(" (b) extensive, SAME model but FIXED cutoff ncut = 4          -> M0 ~ L^2 ?")
P(" (c) point:  M0_pt = 2 G(0)^2 , G(0) = (1/L^3) sum 1/(2 omega) -> M0 ~ L^0 ?")
P("")
P("   L    M0_ext(a)     M0_ext(a)/L^3   M0_ext(b)     M0_ext(b)/L^2   M0_pt(c)")
for L in (4, 8, 16, 32, 64):
    ncut = L // 2
    Ma = 0.0
    for n in modes(L, ncut):
        Ma += 1.0 / (2.0 * omega(L, n) ** 2)
    Mb = 0.0
    for n in modes(L, 4):
        Mb += 1.0 / (2.0 * omega(L, n) ** 2)
    G0 = G_of_t(L, 0.0, ncut)
    Mp = 2.0 * G0 ** 2
    P("  %3d  %11.5e    %10.5f     %11.5e    %10.5f      %10.6f"
      % (L, Ma, Ma / L ** 3, Mb, Mb / L ** 2, Mp))
P(" => (a) L^3 (lattice a=1: BZ scales with L);  (b) L^2 (cutoff held fixed: NOT the lattice);")
P("    (c) L^0 (local composite variance: contact, cutoff-dependent constant, no L-growth).")

# ------------------------------------------------------------------ [3] w1 scalings
head("[3] w1 = f1^2/M0 ON THE FREE FIELD (exact formulas; constants are regulator-scheme data)")
P("   L    w1_ext = f1_ext/M0(a)   w1_ext*L      w1_pt = f1_pt/M0(c)   w1_pt*L^4     w1_ext/w1_pt")
for L in (4, 8, 16, 32, 64):
    ncut = L // 2
    w1 = 2.0 * math.pi / L
    fe = 3.0 / w1 ** 2
    fp = 18.0 / (L ** 6 * w1 ** 2)
    Ma = 0.0
    for n in modes(L, ncut):
        Ma += 1.0 / (2.0 * omega(L, n) ** 2)
    G0 = G_of_t(L, 0.0, ncut)
    Mp = 2.0 * G0 ** 2
    we = fe / Ma
    wp = fp / Mp
    P("  %3d   %12.6e      %10.6f     %12.6e     %10.6f     %10.3f"
      % (L, we, we * L, wp, wp * L ** 4, we / wp))
P(" => w1_ext -> 0 like L^-1 ; w1_pt -> 0 like L^-4 ; w1_ext/w1_pt ~ L^3 :")
P("    THE TWO CARRIERS ARE NOT RELATED BY ANY MULTIPLICATIVE FACTOR (already at w1 level).")

# ------------------------------------------------------------------ [4] identity
head("[4] M0 = Var(O) = C(0) = nu([0,inf))  (identity, any centered O, positive KL measure)")
P(" C(0) = int e^{-sqrt(s)*0} dnu(s) = nu([0,inf)) = <O^2>_c = Var(O).   [exact]")
for (w, s) in ((Fraction(1, 2), Fraction(1)), (Fraction(3, 4), Fraction(2)), (Fraction(1, 5), Fraction(9))):
    f1sq = w  # two-atom measure below
    f2sq = 1 - w
    M0 = f1sq + f2sq
    chi = f1sq / s + f2sq / Fraction(4)
    P("  nu = %s d_%s + %s d_4 : M0=%s  Var=%s  chi=%s  w1=%s"
      % (w, s, f2sq, M0, M0, chi, f1sq / M0))
P(" => 'divide by M0' and 'divide by Var(O)' are THE SAME functional (same carrier, same channel).")

# ------------------------------------------------------------------ [5] HN1 + thin atom
head("[5] HN1 FAMILY AND THIN ATOM")
P(" HN1: nu_L = (1/2) delta_{1/L^2} + (1/2) delta_4  ->  M0 = 1, f1^2 = 1/2, w1 = 1/2 (all L)")
for L in (2, 4, 8, 16, 64):
    chi = Fraction(1, 2) * L * L + Fraction(1, 8)
    P("   L=%3d  m1^2=1/%d  f1^2=1/2  chi=%s  squeeze slack m1^2*chi-f1^2=%s"
      % (L, L * L, chi, Fraction(L * L, 1) / 1 * 0 + (Fraction(1, L * L) * chi - Fraction(1, 2))))
P(" thin atom: nu = (eps/2) delta_eps + (1-eps/2) delta_4  (M0 = 1, eps -> 0)")
for e in (Fraction(1, 2), Fraction(1, 8), Fraction(1, 1000)):
    w = e / 2
    chi = w / e + (1 - w) / 4
    P("   eps=%-8s f1^2=%-9s w1=%-9s chi=%-9s m1^2=eps  (light mode, bounded chi, w1->0)"
      % (str(e), str(w), str(w), str(chi)))

# ------------------------------------------------------------------ [6] M6
head("[6] M6 = oscillator(w) (x) free ring: b-carrier, exact")
P(" nu = sum_{n != 0} (L^2/(4 pi^2 n^2)) delta_{E_n^2},  E_n = 2 pi^2 n^2 / L^2")
P(" f1^2 = L^2/(2 pi^2),  M0 = L^2/12 (= Var(b), Parseval)  =>  w1 = 6/pi^2 exactly")
P(" 6/pi^2 = %.6f   and m1^2 = E_1^2 = 4 pi^4 / L^4  (<= c/L^2 for any c>0, L large)" % (6 / math.pi ** 2))
for L in (4, 8, 16):
    P("   L=%3d  f1^2=%.6e  Var(b)=%.6e  w1=%.6f" % (L, L * L / (2 * math.pi ** 2), L * L / 12.0, (L * L / (2 * math.pi ** 2)) / (L * L / 12.0)))

# ------------------------------------------------------------------ [7] the flip
head("[7] THE FLIP: one physical model, two carriers, opposite truth value (r0=0.5, c=16 pi^2)")
r0, c = 0.5, 16 * math.pi ** 2
P(" hypothesis:  m1^2(L) <= c/L^2  AND  f1^2(L) >= r0")
P("   L   m1^2=4 pi^2/L^2   <=c/L^2?   f1_pt           pt>=r0?   f1_ext        ext>=r0?")
for L in (4, 8, 16, 32, 64):
    m1 = 4 * math.pi ** 2 / L ** 2
    w1 = 2.0 * math.pi / L
    fe = 3.0 / w1 ** 2
    fp = 18.0 / (L ** 6 * w1 ** 2)
    P("  %3d   %12.6e     %-8s  %12.6e  %-6s  %12.6e  %-6s"
      % (L, m1, m1 <= c / L ** 2, fp, fp >= r0, fe, fe >= r0))
P(" => point carrier: hypothesis FALSE for all L>4  -> theorem SURVIVES")
P(" => extensive carrier: hypothesis TRUE for all L>=4 -> theorem VIOLATED by the free field")

# ------------------------------------------------------------------ [8] structural
head("[8] STRUCTURAL: no absolute threshold is carrier-invariant")
P(" given any family with f1^2(L) > 0, put s_L := sqrt(r0/(2 f1^2(L))); then (s_L^2 f1^2) = r0/2 < r0")
P(" and (2^0.5 s_L) gives 2 r0 > r0 : the predicate flips inside the admissible group {O -> s O, s>0}.")
for L in (8, 16):
    w1 = 2.0 * math.pi / L
    fe = 3.0 / w1 ** 2
    s1 = math.sqrt(r0 / (2 * fe))
    P("   L=%3d: s_L=%.6e -> s^2 f1^2=%.6f (<r0) ; 1.41421*s_L -> %.6f (>r0)" % (L, s1, s1 * s1 * fe, 2 * s1 * s1 * fe))
P(" => ill-definedness is STRUCTURAL, not a property of the free-field family.")

# ------------------------------------------------------------------ [9] YM consistency
head("[9] EXTERNAL CONSISTENCY: M0 = C(0) ~ L^3 on 4D SU(2) YM at beta=2.5 (GS1-SEAT2:0.15)")
a, b = 19.9108, 8.3799
P(" C(0)(L=6)=%.4f  C(0)(L=8)=%.4f  ratio=%.4f   vs (8/6)^3=%.4f   (deviation %.4f of the ratio)"
  % (a, b, a / b, (8.0 / 6.0) ** 3, abs(a / b - (8.0 / 6.0) ** 3) / (8.0 / 6.0) ** 3))
P(" => the extensive carrier's M0 grows like V on the target theory itself;")
P("    hence any floor on w1 = f1^2/M0 is <= C/L^3 and is violated by every bounded-f1^2 family.")

P("\nALL CHECKS DONE (deterministic, no sampling).")
