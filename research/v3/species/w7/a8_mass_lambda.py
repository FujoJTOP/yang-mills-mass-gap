# -*- coding: utf-8 -*-
"""A8-MASS-LAMBDA: m/Lambda from the corpus' measured m a, two-loop running.

Question: is m/Lambda a pure number c (constant-type shape m >= c*Lambda) ?
Every normalization step is checked against the standard beta-function form.
"""
import numpy as np

PI2 = np.pi ** 2
PI4 = np.pi ** 4

print("=" * 78)
print("A8  m/Lambda from corpus m a ; two-loop asymptotic scaling")
print("=" * 78)

# ---------------------------------------------------------------- step 0
print("\n[0] beta-function normalization, standard form")
print("    beta(g) = dg/dln mu = -b0 g^3 - b1 g^5 + O(g^7)   [convention]")
print("    b0 = 11N/(3 (4pi)^2) = 11N/(48 pi^2)")
print("    b1 = 34N^2/(3 (4pi)^4) = 34N^2/(768 pi^4) = 34N^2/(3 (16 pi^2)^2)")

def b0_of(N):
    return 11.0 * N / (3.0 * (4 * np.pi) ** 2)

def b1_of(N):
    return 34.0 * N ** 2 / (3.0 * (4 * np.pi) ** 4)

print("\n    N |        b0 |        b1 |  b1/b0^2 | 102/121 | b0*48pi^2/11N | b1*768pi^4/34N^2")
for N in (2, 3):
    b0, b1 = b0_of(N), b1_of(N)
    print(f"    {N} | {b0:.6f} | {b1:.6f} | {b1/b0**2:.6f} | {102/121:.6f} |"
          f" {b0*48*PI2/(11*N):.12f} | {b1*768*PI4/(34*N*N):.12f}")

print("\n    cross-check vs textbook Nf=0 values (units of 1/(4pi)^2k):")
for N in (2, 3):
    print(f"    N={N}: 11N/3 = {11*N/3:.6f}  vs b0*(4pi)^2 = {b0_of(N)*16*PI2:.6f}"
          f"   ;  34N^2/3 = {34*N*N/3:.6f}  vs b1*(4pi)^4 = {b1_of(N)*256*PI4:.6f}")
print("    [T] b1/b0^2 = (34/3)(9/121) = 102/121 exactly, N-independent.")

# ---------------------------------------------------------------- step 1
def aLambda(g2, N, loops=2):
    b0, b1 = b0_of(N), b1_of(N)
    one = np.exp(-1.0 / (2 * b0 * g2))
    if loops == 1:
        return one
    return one * (2 * b0 * g2) ** (-b1 / (2 * b0 ** 2))

def aLambda_exact2loop(g2, N):
    """closed 2-loop solution: 1/(2b0g2) + (b1/2b0^2) ln[b0g2/(1+(b1/b0)g2)] = ln(mu/Lam)"""
    b0, b1 = b0_of(N), b1_of(N)
    c = b1 / b0
    L = b1 / (2 * b0 ** 2)
    return np.exp(-1.0 / (2 * b0 * g2)) * (b0 * g2 / (1 + c * g2)) ** (-L)

print("\n[1] a Lambda = exp(-1/(2b0g2)) (2 b0 g2)^(-b1/(2b0^2))")
print("    g2 convention: beta_W = 2N/g2  =>  g2 = 2N/beta_W  (SU2: 4/b, SU3: 6/b)")
print("\n    N=3, g2=1:  b0g2=%.6f  2b0g2=%.6f  1/(2b0g2)=%.5f  exp=%.6e"
      % (b0_of(3), 2 * b0_of(3), 1 / (2 * b0_of(3)), np.exp(-1 / (2 * b0_of(3)))))
print("    N=3, g2=1:  (2b0g2)^(-b1/2b0^2)=%.6e   -> aLam=%.6e"
      % ((2 * b0_of(3)) ** (-b1_of(3) / (2 * b0_of(3) ** 2)), aLambda(1.0, 3)))
print("    truncation check, exact 2-loop closed form vs task formula (N=3,g2=1):")
print("      exact=%.6e  task=%.6e  ratio=%.4f"
      % (aLambda_exact2loop(1.0, 3), aLambda(1.0, 3), aLambda(1.0, 3) / aLambda_exact2loop(1.0, 3)))

# ---------------------------------------------------------------- step 2
# corpus data  (A4-SCALE-SEP.md  table rows 1,2,3,6,7)
data = [
    # group, beta_W, N, ma, sigma_ma, source
    ("SU(2)", 2.50, 2, 1.2562, 0.1008, "tm_gap_results.csv:3 (L4, Ntau=4)"),
    ("SU(2)", 2.50, 2, 0.8467, 0.1481, "tm_gap_results.csv:6 (L6, Ntau=6)"),
    ("SU(2)", 2.50, 2, 0.770,  0.000,  "T1-POLYMER-RG.md:241-242 (ratio chain, no err)"),
    ("SU(3)", 6.00, 3, 1.073,  0.088,  "lattice-cpb/REPORT.md:357 (n=8, L=6)"),
    ("SU(3)", 6.30, 3, 0.886,  0.159,  "lattice-cpb/REPORT.md:1060 (n=8, L=24)"),
]

print("\n[2] corpus -> Lambda a -> m/Lambda")
hdr = ("group | beta_W |   g2   |   Lambda a  |  m a         |   m/Lambda    | err")
print("    " + hdr)
print("    " + "-" * len(hdr))
rows = []
for grp, beta, N, ma, sma, src in data:
    g2 = 2 * N / beta
    La = aLambda(g2, N)
    ratio = ma / La
    errm = 0 if sma == 0 else (sma / ma) * ratio
    rows.append((grp, beta, g2, La, ma, ratio, errm))
    print(f"    {grp} | {beta:6.2f} | {g2:6.4f} | {La:11.4e} | {ma:5.4f}      |"
          f" {ratio:11.1f}    | {'-' if errm==0 else f'{errm:.1f}'}")

# ---------------------------------------------------------------- step 3
print("\n[3] scaling test")
print("    SU(2): only one coupling beta=2.5 -> CANNOT TEST (three masses share one Lambda a)")
r_su2 = [r[5] for r in rows if r[0] == "SU(2)"]
print(f"      m/Lambda(SU2) = {', '.join(f'{x:.0f}' for x in r_su2)}  [same ^]")
print(f"      spread from mass estimates only: {max(r_su2)/min(r_su2):.2f}x")
s30 = [r for r in rows if r[0] == "SU(3)" and r[1] == 6.0][0]
s63 = [r for r in rows if r[0] == "SU(3)" and r[1] == 6.3][0]
print(f"    SU(3): m/Lambda(6.00) = {s30[5]:.1f} +- {s30[6]:.1f} ;"
      f"  m/Lambda(6.30) = {s63[5]:.1f} +- {s63[6]:.1f}")
sig = abs(s30[5] - s63[5]) / np.hypot(s30[6], s63[6])
print(f"      ratio 6.00/6.30 = {s30[5]/s63[5]:.2f}   (errors ~{s30[6]/s30[5]*100:.0f}%,"
      f" {s63[6]/s63[5]*100:.0f}% -> discrepancy {sig:.2f} sigma)")
print(f"      Lambda a ratio aLam(6.30)/aLam(6.00) = {s63[3]/s30[3]:.2f}"
      f"   vs m a ratio {s63[4]/s30[4]:.3f}  (m a is flat, Lambda a is NOT)")

# ---------------------------------------------------------------- step 4
print("\n[4] structure of aLambda(g2): steep branch vs plateau")
g2s = np.array([0.5, 0.7, 0.953, 1.0, 1.6, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 50.0])
print("     g2  |  aLambda(SU2) |  aLambda(SU3) | m/Lambda if m a = 1")
for g2 in g2s:
    print(f"    {g2:5.3f} | {aLambda(g2,2):11.3e}  | {aLambda(g2,3):11.3e}   |"
          f"   {1/aLambda(g2,3):10.2f} (SU3)")
gg = np.linspace(1e-3, 200, 400001)
am = aLambda(gg, 3)
imax = int(np.argmax(am))
print(f"    max aLambda over g2 in (0,200]: {am[imax]:.4f} at g2 = {gg[imax]:.3f}   [N=3]")
am2 = aLambda(gg, 2)
imax2 = int(np.argmax(am2))
print(f"    max aLambda over g2 in (0,200]: {am2[imax2]:.4f} at g2 = {gg[imax2]:.3f}   [N=2]")
print("    => m/Lambda >= (m a)/%.3f : with m a ~ 1, m/Lambda >= %.2f ALWAYS [N=3]"
      % (am[imax], 1 / am[imax]))

print("\n[5] what the corpus would need for m/Lambda = O(1)  (i.e. aLambda ~ 0.3-1)")
for target in (0.3, 0.5, 1.0):
    idx = np.arange(len(gg))[am > target * 0.999]
    if len(idx):
        print(f"    aLambda >= {target}: g2 >= {gg[idx[0]]:.2f}"
              f"  (beta_W <= {6/gg[idx[0]]:.2f} for SU3)")
    else:
        print(f"    aLambda >= {target}: IMPOSSIBLE for N=3 (max {am[imax]:.3f})")

print("\n[6] consistency: does ma ~ 1 fit the perturbative branch at all?")
print("    if the point were in the 2-loop scaling window with m/Lambda = 3 (say):")
for grp, beta, N, ma, sma, src in data:
    g2 = 2 * N / beta
    pred = 3.0 * aLambda(g2, N)
    print(f"      {grp} beta={beta:5.2f}: predicted m a = 3*aLam = {pred:.3e}"
          f"   vs measured {ma:.4f}   ratio = {ma/pred:9.0f}")

print("\n[7] one-loop cross-check (no b1): same order of magnitude?")
for grp, beta, N, ma, sma, src in data:
    g2 = 2 * N / beta
    print(f"      {grp} beta={beta:5.2f}: m/Lambda_1loop = {ma/aLambda(g2,N,loops=1):10.1f}"
          f"   (2-loop: {ma/aLambda(g2,N):10.1f})")

print("\n[8] scheme sensitivity: c = m/Lambda is NOT scheme independent")
print("    a scheme change is a constant shift d(1/g2);  d ln(m/Lambda)/d(1/g2) = 1/(2 b0)")
for N in (2, 3):
    print(f"      N={N}: 1/(2b0) = {1/(2*b0_of(N)):.4f}  =>"
          f"  d(1/g2)=0.5 moves m/Lambda by e^{0.5/(2*b0_of(N)):.3f} ="
          f" {np.exp(0.5/(2*b0_of(N))):.1f}x")
r_su3 = [r[5] for r in rows if r[0] == "SU(3)"]
print(f"    to pull m/Lambda = {s30[5]:.0f} down to 1 needs d(1/g2) ="
      f" ln({s30[5]:.0f})*2b0 = {np.log(s30[5])*2*b0_of(3):.3f}   (an ordinary O(1) matching constant)")
print("    corpus has NO scheme matching (bare lattice g2 -> MS-bar): unchecked piece.")

print("\n[9] internal a-ratio test (normalization-free)")
print("    two-loop: a(6.30)/a(6.00) = aLam(6.30)/aLam(6.00) = %.4f" % (s63[3] / s30[3]))
print("    measured: m a(6.30)/m a(6.00)     = %.4f   (same physical m assumed)"
      % (s63[4] / s30[4]))
print("    => agreement within the ~18%% error: the beta-dependence of the formula is OK.")
print("    physical anchor (probe-m/REPORT.md:150, Necco-Sommer r0/a=5.368,7.380,9.660):")
print("      a^-1(6.0)=2.91 GeV  =>  m = 1.073*2.91 = %.2f GeV ; m/Lambda_2loop with"
      " Lambda=%.1f MeV" % (1.073 * 2.91, 2.91 * s30[3] * 1000))

print("\n[10] verdict")
print("    (c) fires ON ITS FACE: m/Lambda = 3e2..7e2, i.e. 2.5-2.9 orders above 1.")
print("    (b) DOES NOT fire: SU(3) 6.00 vs 6.30 -> %.0f vs %.0f, diff %.0f +- %.0f = %.2f sigma."
      % (s30[5], s63[5], s30[5] - s63[5], np.hypot(s30[6], s63[6]),
         abs(s30[5] - s63[5]) / np.hypot(s30[6], s63[6])))
print("    (a) does NOT hold: the value is ~6e2, not O(0.1-3).")
print("    normalization audit: b0,b1 = 11N/3, 34N^2/3 in 1/(4pi)^2k units, exact;")
print("      b1/b0^2 = 102/121 exactly; b1/(2b0^2) = 0.4215 (a lost factor 2 here would")
print("      move m/Lambda by (2b0g2)^{0.4215} = %.2f -> checked); exact 2-loop closed form"
      % (2 * b0_of(3)) ** (0.42149))
print("      vs the task formula differs by %.1f%% at g2=1. No factor-2/pi^2 slip can"
      % (100 * (aLambda(1.0, 3) / aLambda_exact2loop(1.0, 3) - 1)))
print("      produce 6e2 from O(1).")
print("    => c = m/Lambda is a scheme-conditional number at g2 ~ 1; not a pure number.")

print("\n[11] SUMMARY")
print("    (a) NOT supported: m/Lambda = 2.9e2..7.1e2, not O(0.1-3).")
print("    (b) NOT fired: SU(3) two points agree (0.71 sigma); a-ratio 0.71 vs 0.83.")
print("    (c) FIRED: m/Lambda > 1 by 2.5-2.9 orders -> normalization/interpretation check,")
print("        done in [0],[1],[8]: arithmetic exact, the gap is the unsupplied scheme")
print("        constant d(1/g2) ~ 0.9 (Lambda_2loop = 5.1 MeV vs physical ~2e2 MeV).")
print("    => c = m/Lambda is NOT a pure number from these data; the corpus fixes its")
print("       shape (beta-constant) but not its normalization.")
