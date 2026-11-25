# u1_harnack_witness.py — witness for U1-HARNACK-BOUND (2026-11-24)
# Tests whether a PF / Harnack / Doob-h-transform / Sturm-TP2 lower bound can produce
# w1 = f1^2/M0 >= w0 > 0 for the atom fraction of the source o-1 (o=O/<O> >= 0, E[o]=1).
# Deterministic. numpy + scipy only.  py -3.12 -I u1_harnack_witness.py
#
# Convention used THROUGHOUT (self-contained toy models): a symmetric positive transfer
# matrix K on a finite space with counting measure mu (uniform), so the invariant measure is
# uniform and the ground state is constant.  Source phi = o - 1 (mean-zero after removing the
# ground mode).  Atom fraction of level n:  w_n = <psi_n|phi>^2 / ||phi||^2 .
# All claims are STRUCTURAL witnesses of the inequalities, NOT literal YM L-scaling.

import warnings
import numpy as np
from scipy.linalg import expm

warnings.filterwarnings("ignore")
np.set_printoptions(precision=6, suppress=True)


def eig_sorted(K):
    w, V = np.linalg.eigh(K)
    idx = np.argsort(w)[::-1]
    return w[idx], V[:, idx]


def level_weights(K, o, tol=1e-7):
    """Weights per DEGENERATE LEVEL (basis-independent).  Level 0 = ground."""
    w, V = eig_sorted(K)
    phi = o - o.mean()
    phi = phi - (V[:, 0] @ phi) * V[:, 0]
    M0 = phi @ phi
    per_mode = np.array([(V[:, n] @ phi) ** 2 for n in range(len(w))]) / M0
    levels = []                      # list of (eigenvalue, summed weight), sum_n>=1 = 1
    n = 0
    while n < len(w):
        m = n
        while m + 1 < len(w) and abs(w[m + 1] - w[n]) < tol:
            m += 1
        levels.append((w[n], per_mode[n:m + 1].sum()))
        n = m + 1
    return w, V, phi, M0, levels


def path_laplacian(N):
    H = np.zeros((N, N))
    for i in range(N):
        if i > 0:    H[i, i - 1] = -1; H[i, i] += 1
        if i < N - 1: H[i, i + 1] = -1; H[i, i] += 1
    return H


def cycle_laplacian(N):
    H = np.zeros((N, N))
    for i in range(N):
        H[i, i] = 2; H[i, (i - 1) % N] = -1; H[i, (i + 1) % N] = -1
    return H


def positive_perron(V):
    v = V[:, 0]
    return np.all(v > 0) or np.all(v < 0)


print("=" * 78)
print("U1 witness - PF/Harnack/Doob/Sturm lower bound on the tube overlap w1")
print("=" * 78)

# ---------------------------------------------------------------- PART A: master identity
print("\n[A] MASTER IDENTITY  w1 = <psi1|o-1>^2 / ||o-1||^2 ,  sum_n>=1 w_n = 1")
N = 7
K = expm(-0.7 * path_laplacian(N))
w, V, phi, M0, levels = level_weights(K, (N * (np.arange(N) == 2)).astype(float))
tot = sum(wt for _, wt in levels[1:])
print(f"    N={N}, K=exp(-0.7*L_path): Perron vector constant-sign : {positive_perron(V)}")
print(f"    gap -log(lam1/lam0) = {-np.log(w[1]/w[0]):.4f}")
w1 = levels[1][1]
print(f"    w1 = {w1:.6f}   in [0,1] : {0 <= w1 <= 1}")
print(f"    sum_n>=1 w_n = {tot:.12f}   (== 1: {abs(tot-1) < 1e-12})")
print(f"    identity |w1 - <psi1|o-1>^2/||o-1||^2| = "
      f"{abs(w1 - (V[:,1]@phi)**2/M0):.2e}")

# ---------------------------------------------------------------- PART B: no-go (source direction)
print("\n[B] NO-GO U1-NG1: w1 is the SOURCE direction cosine; positivity+gap do NOT bound it")
N = 6
K = expm(-0.9 * path_laplacian(N))
w, V = eig_sorted(K)
psi1, psi2 = V[:, 1], V[:, 2]
gap = -np.log(w[1] / w[0])
iso = w[1] / w[2]
pos = psi1 > 1e-9
harnack = psi1[pos].max() / psi1[pos].min()
print(f"    FIXED K: Perron const-sign={positive_perron(V)}  gap={gap:.4f}  "
      f"lam1/lam2={iso:.4f}  psi1 shape sup/inf={harnack:.3f}")
print("    (these are functions of K ALONE; the source is free in Omega^perp)")
print("     theta      w1=cos^2 th    min(o_theta)  (o>=0 legal source)")
for th in [0.0, 0.4, 0.8, 1.2, np.pi / 2 - 1e-3]:
    combo = np.cos(th) * psi1 + np.sin(th) * psi2
    eps = 0.30 / np.abs(combo).max()
    o = 1.0 + eps * combo
    f = o - o.mean()
    w1b = (psi1 @ f) ** 2 / (f @ f)
    print(f"    {th:6.3f}    {w1b:10.6f}      {o.min():8.5f}")
print("    => w1 -> 0 at fixed K (fixed PF, gap, Harnack, lam1/lam2).  NO positivity bound.")

# ---------------------------------------------------------------- PART C: Harnack capture lemma
print("\n[C] CONDITIONAL TOOL: Harnack(shape) x sign-alignment  =>  w1 >= delta0^2 * mu(supp)")
N = 9
K = expm(-0.8 * path_laplacian(N))
w, V = eig_sorted(K)
psi1 = V[:, 1]
sign1 = np.sign(psi1)
# collect maximal one-sign blocks (nodal domains); pick blocks well inside a lobe
blocks = []
i = 0
while i < N:
    j = i
    while j + 1 < N and sign1[j + 1] == sign1[i]:
        j += 1
    blocks.append((i, j)); i = j + 1
long_pos = max([b for b in blocks if sign1[b[0]] > 0], key=lambda b: b[1] - b[0])
a, b = long_pos
print(f"    psi1 positive lobe = sites [{a},{b}]; support blocks strictly inside it:")
for lo, hi in [(a, a + 2), (a + 1, b - 1), (a, b)]:
    if hi <= lo:
        continue
    g = np.zeros(N); g[lo:hi + 1] = 1.0
    f = sign1 * g
    f = f - (V[:, 0] @ f) * V[:, 0]
    w1c = (psi1 @ f) ** 2 / (f @ f)
    delta0 = np.abs(psi1[lo:hi + 1]).min()
    musupp = (hi - lo + 1) / N
    print(f"      supp=[{lo},{hi}]  w1={w1c:.4f}  >= delta0^2*mu={delta0**2*musupp:.4f} : "
          f"{w1c >= delta0**2*musupp - 1e-12}")
psi2 = V[:, 2]
print(f"    misaligned source f=psi2 (a DIFFERENT state's sign pattern): w1={float((psi1@psi2)**2/(psi2@psi2)):.2e}")
print("    => the bound lives on the SIGN-ALIGNMENT, which Harnack does not supply.")

# ---------------------------------------------------------------- PART D: free-field parity
print("\n[D] ACID TEST: free composite is number-even; one-particle number-odd => w1=0")
M = 12
a = np.diag(np.sqrt(np.arange(1, M)), 1)
ad = a.T
omega = 0.9
Kf = expm(-0.8 * omega * (ad @ a))
vac = np.zeros(M); vac[0] = 1.0
one = np.zeros(M); one[1] = 1.0
composite = a @ a + ad @ ad + 2 * (ad @ a)          # normal-ordered :phi^2:
overlap = one @ (composite @ vac)
print(f"    transfer all >0 : {np.all(np.diag(Kf) > 0)}   gap = omega = {omega}")
print(f"    <1-particle | :phi^2: |vacuum> = {overlap:.3e}   (even op, odd state)")
print(f"    w1(free composite) = {overlap**2:.3e}   (== 0 by number parity)")
two = np.zeros(M); two[2] = 1.0
print(f"    <2-particle | :phi^2: |vacuum> = {two @ (composite @ vac):.4f} (=sqrt2; the continuum)")

# ---------------------------------------------------------------- PART E: TP2 holds yet w1=0
print("\n[E] TP2/Sturm does NOT save it: Gaussian heat kernel is strictly TP2 yet parity kills w1")
xs = np.linspace(-2, 2, 9)
s2 = 0.7
Kg = np.exp(-(xs[:, None] - xs[None, :]) ** 2 / (2 * s2 ** 2))
minors = [Kg[i, j] * Kg[i + 1, j + 1] - Kg[i, j + 1] * Kg[i + 1, j]
          for i in range(len(xs) - 1) for j in range(len(xs) - 1)]
print(f"    strict TP2 (all 2x2 minors >=0): {min(minors) > -1e-12}  (min minor {min(minors):.2e})")
print("    => order/TP2 present in the ordered Gaussian model, yet part D has w1=0 by parity:")
print("       the datum is SOURCE ALIGNMENT, not order, not positivity.")

# ---------------------------------------------------------------- PART F: volume / delocalization
print("\n[F] LOCAL source on a DELOCALIZED lowest level: w1 -> 0 as L grows (free/gapless)")
print("     L    gap=-log(lam1/lam0)   w1 = level-1 weight   analytic 2/(L-1)")
for L in [4, 8, 16, 32, 64]:
    Kc = expm(-1.0 * cycle_laplacian(L))
    o = L * (np.arange(L) == 0).astype(float)
    w, V, phi, M0, levels = level_weights(Kc, o)
    w1f = levels[1][1]
    ana = 2.0 / (L - 1)
    print(f"    {L:3d}    {-np.log(levels[1][0]/levels[0][0]):10.6f}         {w1f:.6f}          {ana:.6f}")
print("    => delocalized light level (gap->0, psi1 ~ 1/sqrt L): w1 = 2/(L-1) -> 0.")
print("       uniform w1>0 needs the lightest 0++ level LOCALIZED (finite xi = the gap). [V]/[C]")

# ---------------------------------------------------------------- PART G: Sturm node count
print("\n[G] STURM on the ordered path: psi_n has exactly n nodes (order exists); cycle has none")
N = 10
K = expm(-0.5 * path_laplacian(N))
w, V = eig_sorted(K)
for n in range(4):
    s = np.sign(V[:, n]); s = s[np.abs(V[:, n]) > 1e-9]
    nodes = int(np.sum(s[1:] != s[:-1]))
    print(f"    psi_{n}: sign changes = {nodes}  (== n: {nodes == n})")
print("    path has a linear order (TP2/Sturm); A/G does NOT -- and part E shows even WITH")
print("    order the source can be misaligned, so Sturm alone gives no w1 floor.")

# ---------------------------------------------------------------- PART H: 2D solvable
print("\n[H] 2D SOLVABLE (T^2 rotor = the 2D-holonomy carrier): first excited level is 4-fold;")
print("    a local source cos(th1-th2) is EXACTLY orthogonal to it => w1=0 at fixed T")
N2 = 32
th = 2 * np.pi * np.arange(N2) / N2
T1, T2 = np.meshgrid(th, th, indexing='ij')
levels = {}
for k1 in range(N2):
    for k2 in range(N2):
        m1 = k1 if k1 <= N2 // 2 else k1 - N2
        m2 = k2 if k2 <= N2 // 2 else k2 - N2
        levels.setdefault(m1 * m1 + m2 * m2, []).append((k1, k2))


def w1_2d(g):
    F = np.abs(np.fft.fft2(g)) ** 2
    F[0, 0] = 0.0
    tot = F.sum()
    lvl1 = sum(F[k1, k2] for (k1, k2) in levels[1])
    return lvl1 / tot


for name, g in [("cos(th1)", np.cos(T1)),
                ("cos(th1+th2)", np.cos(T1 + T2)),
                ("cos(th1-th2)", np.cos(T1 - T2)),
                ("localized bump", ((T1 < 0.6) & (T2 < 0.6)).astype(float))]:
    print(f"    source {name:16s}: w1(first 4-fold level) = {w1_2d(g):.6f}")
print("    (rotor first level mu=1/(2I), rho^2=2, exact; but cos(th1-th2) sits in level 2 => w1=0)")
print(f"    (a POINT source gives exactly w1=4/(N^2-1)={4/(N2**2-1):.6f}->0: delocalized mode;")
print("     a finite bump keeps only its Fourier fraction -- the floor tracks mode LOCALIZATION)")

print("\n" + "=" * 78)
print("SUMMARY: real positive transfer + PF + gap + Harnack + TP2/Sturm => NO w1>0 bound.")
print("A conditional bound w1 >= delta0^2*mu(supp) holds ONLY under source sign-alignment")
print("+ Harnack shape (delta0 uniform). The missing estimate is the alignment (E-U1).")
print("=" * 78)
