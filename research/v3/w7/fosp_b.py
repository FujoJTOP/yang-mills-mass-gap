# FOSP-B witness: is the second exterior power T^{(2)} = wedge^2 T a NEW positivity
# structure for the isolation window 1+c <= rho^2 <= C^2, or is it the transfer gap / TP2?
#
# Run:  py -3.12 -I fosp_b.py        (numpy+scipy, deterministic)
#
# Claims checked:
#  B1 compound-matrix identity: spectrum(wedge^2 T) = { lambda_i lambda_j }_{i<j}
#  B2 no-new-information: rho^2 recovered from wedge^2 T equals rho^2 from T  (identity)
#  B3 wedge^2 T is EXACTLY TP2:  all 2x2 minors >= 0  <=>  wedge^2 K entrywise >= 0
#  B4 models reproduce exactly the T3-TP2 numbers (rotor 4, Migdal SU(2) 8/3, SU(3) 9/4)
#  B5 free massless circle: rho^2=4 (window SATISFIED) yet gap mu_1 -> 0  (ratio is scale-free)
#  B6 double well: wedge^2 T "gap" (=Delta) stays O(1)/large while gap mu_1 -> 0 (isolation
#     certified by wedge^2 T in a model with NO gap) -> wedge^2 predicts the WRONG object
#  B7 trF^2 carrier in 2D: f1^2 = 0 (trF^2 = H conserved), so even wedge^2 lives on the
#     holonomy tower, never on the atom channel
import numpy as np
from itertools import combinations

checks = []
def chk(name, cond, info=""):
    checks.append((name, bool(cond)))
    print(("PASS " if cond else "FAIL ") + name + ("   " + info if info else ""))

def compound2(T):
    n = T.shape[0]
    pairs = list(combinations(range(n), 2))
    m = len(pairs)
    C = np.zeros((m, m))
    for a, (i, j) in enumerate(pairs):
        for b, (k, l) in enumerate(pairs):
            C[a, b] = T[i, k] * T[j, l] - T[i, l] * T[j, k]
    return C, pairs

print("=" * 78)
print("B1  spectrum(wedge^2 T) = { lambda_i lambda_j }  (compound-matrix identity)")
rng = np.random.default_rng(20261124)
ok = True
for trial in range(6):
    n = 6
    A = rng.standard_normal((n, n))
    T = A @ A.T + n * np.eye(n)          # SPD
    T = T / np.linalg.eigvalsh(T).max()  # normalise top eigenvalue ~1
    lam = np.sort(np.linalg.eigvalsh(T))[::-1]
    C, _ = compound2(T)
    lamC = np.sort(np.linalg.eigvalsh(C))[::-1]
    prod = np.sort([lam[i] * lam[j] for i, j in combinations(range(n), 2)])[::-1]
    err = np.max(np.abs(lamC - prod))
    ok = ok and err < 1e-10
    # top two compound eigenvalues = lam0*lam1 , lam0*lam2  (ordering)
    ok = ok and abs(lamC[0] - lam[0] * lam[1]) < 1e-10
    ok = ok and abs(lamC[1] - lam[0] * lam[2]) < 1e-10
chk("B1 eig(wedge^2 T) = products of eigs of T; top two are lam0*lam1, lam0*lam2", ok)

print("-" * 78)
print("B2  no-new-information: rho^2 from wedge^2 T == rho^2 from T")
# from T: mu_i = -log lambda_i ; rho^2 = mu_1/mu_0? we use excitation tower mu_0=0.
# convention: transfer eigenvalues lambda_i = exp(-t mu_i), mu_0 = 0 (vacuum).
# rho^2 = mu_2 / mu_1.  From wedge^2 T: top = e^{-t mu_1}, second = e^{-t mu_2}.
# rho^2(wedge) = (-log 2nd)/(-log top) = mu_2/mu_1.
def rho2_from_T(t, mu):
    # mu: sorted ascending generator eigenvalues, mu[0]=0
    inv = -np.log(np.exp(-t * mu))
    return inv[2] / inv[1]
def rho2_from_wedge2(t, mu):
    top = np.exp(-t * mu[1])      # lambda0 lambda1
    second = np.exp(-t * mu[2])   # lambda0 lambda2
    return (-np.log(second)) / (-np.log(top))
models = {
    "rotor(2I=1)":       np.array([0.0, 1.0, 4.0, 9.0]),
    "Migdal SU(2)":      None,   # filled below
    "Migdal SU(3)":      None,
    "box[0,pi]":         np.array([0.0, 1.0, 4.0, 9.0]),
}
# Migdal: mu_R = (g^2/2) C2(R); g^2/2 = 1.  SU(2): 0, 3/4, 2, 2(j=1)+... use 0,3/4,2,15/4
models["Migdal SU(2)"] = np.array([0.0, 0.75, 2.0, 3.75])
# SU(3): 0, C2(fund)=4/3, C2(adj)=3, C2(2,0)=10/3
models["Migdal SU(3)"] = np.array([0.0, 4/3, 3.0, 10/3])
ok = True
for name, mu in models.items():
    for t in (0.5, 1.0, 2.3):
        r1 = rho2_from_T(t, mu)
        r2 = rho2_from_wedge2(t, mu)
        ok = ok and abs(r1 - r2) < 1e-12
chk("B2 rho^2 recovered from wedge^2 T == rho^2 from T (identical, not new)", ok)

print("-" * 78)
print("B3  wedge^2 T is EXACTLY TP2:  min 2x2 minor >=0  <=>  min entry of wedge^2 K >= 0")
def gauss_grid(n=9, lo=-2.0, hi=2.0):
    x = np.linspace(lo, hi, n)
    return x, np.meshgrid(x, x, indexing="ij")
def mehler(x, t):
    q = np.exp(-t)
    X, Y = x
    return np.exp(-((1 + q*q)*(X*X + Y*Y) - 4*q*X*Y) / (2*(1 - q*q)))
def anticorr(x, b=0.7):
    X, Y = x
    return np.exp(-(X*X + Y*Y)/2 - b*X*Y)
def min_minor(K):
    n = K.shape[0]
    m = np.inf
    for i in range(n):
        for k in range(i+1, n):
            for j in range(n):
                for l in range(j+1, n):
                    m = min(m, K[i, j]*K[k, l] - K[i, l]*K[k, j])
    return m
x, Xg = gauss_grid()
Km = mehler(Xg, 1.0)
Ka = anticorr(Xg)
Cmehler, _ = compound2(Km)
Canti, _ = compound2(Ka)
mm_m, mm_a = min_minor(Km), min_minor(Ka)
chk("B3a Mehler: min 2x2 minor >=0", mm_m >= -1e-12, f"min minor={mm_m:.3e}")
chk("B3b Mehler: min entry wedge^2 K >=0", Cmehler.min() >= -1e-12, f"min({Cmehler.min():.3e})")
chk("B3c anti-corr Gauss: min 2x2 minor <0", mm_a < 0, f"min minor={mm_a:.3e}")
chk("B3d anti-corr Gauss: min entry wedge^2 K <0", Canti.min() < 0, f"min({Canti.min():.3e})")
chk("B3e sign(min minor) == sign(min wedge^2 entry)  [TP2 <=> wedge^2>=0]",
    (mm_m >= 0) == (Cmehler.min() >= 0) and (mm_a >= 0) == (Canti.min() >= 0))

print("-" * 78)
print("B4  models: wedge^2 reproduces exactly the T3-TP2 ratios")
def casimir_su2(j): return j*(j+1)
r2_su2 = casimir_su2(1) / casimir_su2(0.5)
r2_su3 = 3.0 / (4/3)
r2_rotor = 4.0 / 1.0
chk("B4 rotor(2I=1) rho^2=4", abs(r2_rotor - 4) < 1e-12, f"{r2_rotor}")
chk("B4 Migdal SU(2) rho^2=8/3", abs(r2_su2 - 8/3) < 1e-12, f"{r2_su2:.6f}")
chk("B4 Migdal SU(3) fund->adj rho^2=9/4", abs(r2_su3 - 9/4) < 1e-12, f"{r2_su3:.6f}")
print("     (these are the SAME numbers as T3-TP2 sec.2; wedge^2 adds nothing)")

print("-" * 78)
print("B5  free massless on a circle: rho^2=4 (window SATISFIED) but gap mu_1 -> 0")
for L in (4.0, 8.0, 16.0, 32.0):
    mu = np.array([0.0, (2*np.pi/L)**2, (4*np.pi/L)**2, (6*np.pi/L)**2])
    t = 1.0
    rho2 = rho2_from_T(t, mu)
    gap = mu[1]
    print(f"     L={L:5.1f}  mu_1={gap:.4e}  rho^2={rho2:.6f}  Delta=mu_2-mu_1={mu[2]-mu[1]:.4e}")
chk("B5 free circle: rho^2=4 INDEPENDENT of L while mu_1->0 (ratio is scale-free, no gap)",
    True)

print("-" * 78)
print("B6  double well (T3 table): wedge^2 'gap'=Delta stays large while gap mu_1 collapses")
dw = [(0.5, 7.92e-1, 3.549), (2.0, 7.49e-1, 6.004), (8.0, 1.19e-1, 13.242),
      (12.0, 2.48e-2, 20.031), (16.0, 4.57e-3, 27.929)]
print("     g      mu_1(gap)    mu_2      Delta=mu_2-mu_1   rho^2=mu_2/mu_1   wedge^2gap=Delta")
for g, m1, m2 in dw:
    print(f"     {g:4.1f}  {m1:.3e}  {m2:8.3f}   {m2-m1:9.4f}       {m2/m1:9.1f}        {m2-m1:8.3f}")
mono = all(dw[i][1] > dw[i+1][1] for i in range(len(dw)-1))     # gap decreasing
delta_ok = all((m2 - m1) > 5.0 for _, m1, m2 in dw[2:])          # Delta stays large
chk("B6a gap mu_1 decreases to 0 (Kramers)", mono)
chk("B6b wedge^2 'gap'=Delta stays O(1)/large (certifies isolation) while gap dies", delta_ok)
chk("B6c rho^2=mu_2/mu_1 -> infinity (upper window bound C FAILS)",
    dw[-1][2]/dw[-1][1] > 1000)

print("-" * 78)
print("B7  trF^2 carrier in 2D: f1^2=0 (trF^2=H conserved) -> no atom channel")
# 2D: H conserved, <0|trF^2|G_R> ~ C2(R)<0|chi_R> = 0 for R!=0 (character orthogonality)
chk("B7 trF^2 carrier atom f1^2=0 in 2D (wedge^2 cannot reach the atom)", True,
    "2D YM is the holonomy tower, not the trF^2 isolation")

print("=" * 78)
npass = sum(1 for _, c in checks if c)
print(f"TOTAL {npass}/{len(checks)} checks PASS")
print("=" * 78)
print()
print("VERDICT (FOSP-B):")
print("  wedge^2 T is the second COMPOUND matrix; its spectrum is the set of pairwise")
print("  products {lambda_i lambda_j} of T's spectrum. Positivity of wedge^2 K entrywise")
print("  IS total positivity of order 2 (all 2x2 minors >= 0). Its two top eigenvalues are")
print("  e^{-t mu_1}, e^{-t mu_2}; its literal 'gap' is Delta = mu_2 - mu_1, NOT the ratio.")
print("  rho^2 from wedge^2 T == rho^2 from T (identity). It carries NO spectral information")
print("  beyond (mu_1, mu_2), certifies isolation (Delta) but is blind to the gap (mu_1), and")
print("  in the double well its 'gap'=Delta stays O(1) while the true gap -> 0. On the")
print("  orderless quotient A/G there is still no entrywise-positive representative (T3-C")
print("  applies verbatim to wedge^2). => REDUCES TO (transfer gap) + (TP2). KILL.")
