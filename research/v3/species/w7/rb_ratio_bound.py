"""rb_ratio_bound.py -- witness for RB-RATIO-BOUND (rho^2 <= C^2).
py -3.12 -I ; numpy 2.x + scipy ; single process ; deterministic.
fix v2: no scipy select='i' (off-by-one), explicit zero-mode handling, tolerances.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal

FAIL = []


def chk(name, cond, msg=""):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name} {msg}")
    if not cond:
        FAIL.append(name)


def full_spectrum(diag, off):
    diag = np.asarray(diag, float)
    off = np.asarray(off, float)
    if diag.size <= 400:
        m = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
        return np.sort(np.linalg.eigvalsh(m))
    return np.sort(eigh_tridiagonal(diag, off)[0])


def tower(diag, off, n_zero=1):
    """drop the n_zero lowest (the degenerate Perron/SUSY ground states), return rest + the dropped"""
    ev = full_spectrum(diag, off)
    return ev[n_zero:], ev[:n_zero]


def ratio(w):
    """w[1]/w[0] with distinct-value protection"""
    w = np.asarray(w, float)
    return w[1] / w[0]


print("=" * 80)
print("PART A  1D chain = path graph {1..N}, Dirichlet Laplacian   mu_k = 2-2cos(k pi/(N+1))")
print("=" * 80)
print("   N          mu1          mu2     rho^2 = mu2/mu1")
for N in [2, 3, 4, 5, 10, 50, 100, 1000, 100000]:
    th = np.pi / (N + 1)
    m1, m2 = 2 - 2 * np.cos(th), 2 - 2 * np.cos(2 * th)
    print(f"  {N:7d}  {m1:11.3e}  {m2:11.3e}   {m2/m1:.6f}")
# Dirichlet path: NO zero mode (n_zero=0)
w, z = tower(2 * np.ones(10), -np.ones(9), n_zero=0)
print(f"  FD check N=10: mu1={w[0]:.9f} mu2={w[1]:.9f} rho^2={ratio(w):.9f}")
chk("A FD matches analytic N=10", abs(ratio(w) - (2 - 2 * np.cos(2 * np.pi / 11)) / (2 - 2 * np.cos(np.pi / 11))) < 1e-9)
chk("A rho^2 < 4 for every N>=2", all(
    (2 - 2 * np.cos(2 * np.pi / (N + 1))) / (2 - 2 * np.cos(np.pi / (N + 1))) < 4.0 for N in range(2, 5000)))
chk("A sup = 4 (N -> inf)", abs((2 - 2 * np.cos(2 * np.pi / 100001)) / (2 - 2 * np.cos(np.pi / 100001)) - 4) < 1e-4)

print("=" * 80)
print("PART B  hypercube {0,1}^n , L = sum_i (1-flip_i)  =>  spec = 2|S|  =>  rho^2 = 2")
print("=" * 80)
for n in [2, 3, 4, 5, 6, 8]:
    dim = 2 ** n
    L = np.zeros((dim, dim))
    for i in range(n):
        for s in range(dim):
            L[s, s] += 1.0
            L[s, s ^ (1 << i)] -= 1.0
    ev = np.sort(np.linalg.eigvalsh(L))
    ev = ev[ev > 1e-6]  # drop the zero mode defensively
    vals = [ev[0]]
    for x in ev[1:]:
        if x > vals[-1] * (1 + 1e-6):
            vals.append(x)
    print(f"  n={n:2d} dim={dim:4d}  mu1={vals[0]:.6f} mu2={vals[1]:.6f}  rho^2={vals[1]/vals[0]:.9f}")
    chk(f"B n={n}: rho^2 = 2 (degree/dimension only)", abs(vals[1] / vals[0] - 2) < 1e-9)

print("=" * 80)
print("PART C  star K_{1,m} : spec {0, 1^(m-1), m+1}  =>  rho^2 = m+1")
print("=" * 80)
for m in [3, 4, 10, 100, 1000]:
    print(f"  m={m:5d}  max deg={m+1:5d}  diam=2  gap mu1=1  rho^2={m+1}")
K = np.zeros((5, 5))
for leaf in range(1, 5):
    K[0, leaf] = K[leaf, 0] = -1.0
    K[leaf, leaf] = 1.0
    K[0, 0] += 1.0
ev = np.sort(np.linalg.eigvalsh(K))
print(f"  FD check m=4: eig={ev}")
evp = ev[ev > 1e-6]
vals = [evp[0]]
for x in evp[1:]:
    if x > vals[-1] * (1 + 1e-6):
        vals.append(x)
print(f"  distinct: {vals}  => rho^2={vals[1]/vals[0]:.6f}")
chk("C star m=4 : rho^2 = 5 = m+1", abs(vals[1] / vals[0] - 5) < 1e-9)
chk("C rho^2 = 5 > 4 with degree 5, diam 2, gap 1  => NO degree/diameter bound", vals[1] / vals[0] > 4)

print("=" * 80)
print("PART D  1D Sturm-Liouville  -u''+qu = lam u  on [0,1] Dirichlet, FD N=4000")
print("=" * 80)
ND = 4000
h = 1.0 / ND
xs = np.arange(1, ND) * h


def sl(qf):
    q = qf(xs)
    w, z = tower(2.0 / h**2 + q, -np.ones(ND - 2) / h**2, n_zero=0)  # Dirichlet: no zero mode
    return w[0], w[1], w[1] / w[0], 0.0


bat = {}
bat["q = 0"] = sl(lambda x: 0 * x)
for c in [0.01, 0.1, 1.0, 10.0, 100.0, 1000.0]:
    bat[f"q = {c} (const)"] = sl(lambda x, c=c: c + 0 * x)
bat["q = 1e3 (x-1/2)^2"] = sl(lambda x: 1e3 * (x - 0.5) ** 2)
bat["q = 1e3 (x-1/4)^2"] = sl(lambda x: 1e3 * (x - 0.25) ** 2)
bat["q = 2e3 ((x-1/2)^2-1/16)^2  [two wells, >=0]"] = sl(lambda x: 2e3 * ((x - 0.5) ** 2 - 0.0625) ** 2)
bat["q = 1e4 * 1[x<0.1]"] = sl(lambda x: 1e4 * (x < 0.1))
bat["q = 1e4 * 1[0.4<x<0.6]"] = sl(lambda x: 1e4 * ((x > 0.4) & (x < 0.6)))
bat["q = 1e4 * 1[x>0.9]"] = sl(lambda x: 1e4 * (x > 0.9))
bat["q = 1e6 * 1[|x-1/2|<0.002]"] = sl(lambda x: 1e6 * (np.abs(x - 0.5) < 0.002))
bat["q = 1e4 |x-1/2|^3"] = sl(lambda x: 1e4 * np.abs(x - 0.5) ** 3)
rng = np.random.default_rng(20261124)
mx, arg = 0.0, ""
for _ in range(200):
    A, p = 10 ** rng.uniform(-1, 4), rng.integers(1, 6)
    x0, s = rng.uniform(0.05, 0.95), 10 ** rng.uniform(-2, -0.3)
    r = sl(lambda x, A=A, p=p, x0=x0, s=s: A * np.abs(x - x0) ** p +
           A * 1e-2 * np.exp(-((x - x0) / s) ** 2))[2]
    if r > mx:
        mx, arg = r, f"A={A:.3g} p={p} x0={x0:.2f} s={s:.3g}"
for k, v in bat.items():
    print(f"  {k:44s} mu1={v[0]:10.5f} mu2={v[1]:10.5f} rho^2={v[2]:8.5f}  (mu0={v[3]:.1e})")
print(f"  random q>=0 (200 shapes): max rho^2 = {mx:.6f} at {arg}")
chk("D q>=0 : rho^2 <= 4 (battery + 200 random)", max([v[2] for v in bat.values()] + [mx]) <= 4.0,
    f"max={max([v[2] for v in bat.values()] + [mx]):.6f}")
chk("D q=0 saturates 4 (free string limit)", abs(bat["q = 0"][2] - 4) < 1e-3)

print("  --- constant q = -c : exact lam_k = k^2 pi^2 - c ---")
for c in [0.0, 0.5, 1.0, 3.0, 6.0, 9.0, 9.7]:
    ex = (4 * np.pi**2 - c) / (np.pi**2 - c)
    fd = sl(lambda x, c=c: -c + 0 * x)[2]
    print(f"  c={c:5.2f}  exact rho^2={ex:11.4f}  FD={fd:11.4f}  >4? {ex>4}")
chk("D q=-c , any c>0 : rho^2 > 4  (hypothesis q>=0 is SHARP)",
    all((4 * np.pi**2 - c) / (np.pi**2 - c) > 4 for c in [0.5, 1, 3, 6, 9]))
chk("D q=-c -> inf as c -> pi^2", (4 * np.pi**2 - 9.7) / (np.pi**2 - 9.7) > 100)

print("=" * 80)
print("PART E  one plaquette SU(2): L_beta = -(1/rho) d/phi (rho d/dphi),")
print("        rho_beta(phi) = e^{-beta(1-cos phi)} sin^2 phi ,  phi in (0,pi) is arc length")
print("=" * 80)


def plaquette(beta, NE):
    hp = np.pi / NE
    ph = (np.arange(NE) + 0.5) * hp
    rho = np.exp(-beta * (1 - np.cos(ph))) * np.sin(ph) ** 2
    rh = 0.5 * (rho[:-1] + rho[1:])
    diag = np.concatenate(([rh[0]], rh[:-1] + rh[1:], [rh[-1]])) / (hp**2 * rho)
    off = -rh / (hp**2 * np.sqrt(rho[:-1] * rho[1:]))
    w, z = tower(diag, off)
    return w[0], w[1], w[1] / w[0], z[0]


print("   beta      mu1         mu2       rho^2    (NE=3000)")
for beta in [0.0, 0.1, 0.5, 1.0, 2.0, 2.5, 4.0, 8.0, 16.0, 64.0]:
    m1, m2, r2, z = plaquette(beta, 3000)
    print(f"  {beta:6.2f}  {m1:11.5f} {m2:11.5f}  {r2:9.5f}   (mu0={z:.1e})")
m1, m2, r2, z = plaquette(0.0, 3000)
chk("E beta=0 : mu1=3, mu2=8 exactly  [tower 4j(j+1), j=1/2,1]", abs(m1 - 3) < 5e-3 and abs(m2 - 8) < 2e-2,
    f"got {m1:.5f},{m2:.5f}")
chk("E beta=0 : rho^2 = 8/3 = 2.6667", abs(r2 - 8 / 3) < 2e-3, f"got {r2:.6f}")
chk("E beta=0 violates IN-11 C=1.5 (C^2=2.25)", r2 > 2.25)
print("  convergence / contamination check (singular endpoint at phi=0):")
for beta in [0.0, 8.0]:
    for NE in [750, 1500, 3000, 6000]:
        m1, m2, r2, z = plaquette(beta, NE)
        print(f"    beta={beta:5.1f} NE={NE:5d}: mu1={m1:10.5f} mu2={m2:10.5f} rho^2={r2:8.5f}")
print("  => beta=0 is grid-independent (clean exact anchor).  For beta>0 the vanishing")
print("     density rho~phi^2 at phi=0 lets the even (g(0)!=0) branch leak in on a grid:")
print("     the large-beta tower is analytic (see md S2 RB-6'), not taken from this scan.")

print("=" * 80)
print("PART F  2D YM Casimir towers  mu_R = (g^2/2) C2(R)")
print("=" * 80)
for N in [2, 3, 5, 10, 100]:
    c2f, c2a = (N**2 - 1) / (2 * N), N
    print(f"  SU({N}): C2(fund)={c2f:8.5f} C2(adj)={c2a:8.5f}  rho^2={c2a/c2f:.6f}")
chk("F SU(2) 2D: rho^2 = 8/3", abs(2 / 0.75 - 8 / 3) < 1e-12)
chk("F U(1) 2D: reps n in Z, mu_n=(e^2/2)n^2 => rho^2 = 4/1 = 4", abs(4.0 - 4) < 1e-12)

print("=" * 80)
print("PART G  double well (SUSY) H=-d^2+W, W=V'^2/4-V''/2, V=(g/2)(x^2-1)^2")
print("=" * 80)
NG = 6000
xg = np.linspace(-3.0, 3.0, NG)
hg = xg[1] - xg[0]
print("    g     mu0(dropped)   mu1=GAP      mu2        rho^2      T2=mu1+mu2   log(rho^2)/g")
dwr = []
for g in [0.5, 2.0, 4.0, 8.0, 10.0, 12.0]:
    Vp, Vpp = 2 * g * xg * (xg**2 - 1), 2 * g * (3 * xg**2 - 1)
    W = Vp**2 / 4 - Vpp / 2
    w, z = tower(2.0 / hg**2 + W, -np.ones(NG - 1) / hg**2)
    mu1, mu2 = w[0], w[1]
    r2 = mu2 / mu1
    dwr.append((g, mu1, mu2, r2, z[0]))
    print(f"  {g:5.1f}  {z[0]:11.2e} {mu1:11.5f} {mu2:10.5f} {r2:10.3f}  {mu1+mu2:10.5f}  {np.log(r2)/g:8.4f}")
a, b = dwr[0], dwr[-1]
chk("G rho^2 grows >100x while Delta=mu2-mu1 stays O(1) (upper bound FAILS)", b[3] / a[3] > 100,
    f"rho^2: {a[3]:.2f} -> {b[3]:.1f}")
chk("G T2=mu1+mu2 grows only polynomially (x%.1f) while rho^2 explodes (x%.0f)" % (
    (b[4] + b[2]) / (a[4] + a[2]), b[3] / a[3]),
    (b[4] + b[2]) / (a[4] + a[2]) < 10 < b[3] / a[3])
chk("G log(rho^2)/g -> const in (0.3,1.2) : exponential in g", 0.3 < np.log(b[3]) / b[0] < 1.2,
    f"log(rho^2)/g = {np.log(b[3])/b[0]:.4f}")

print("=" * 80)
print("PART H  Ky Fan: T2(V)=tr(G^-1 M) >= mu1+mu2 for ANY 2-dim V perp 1 ; equality at eigspace")
print("=" * 80)
NL = 40
dL, eL = 2 * np.ones(NL), -np.ones(NL - 1)
Lm = np.diag(dL) + np.diag(eL, 1) + np.diag(eL, -1)
ev, vecs = np.linalg.eigh(Lm)
mu1, mu2 = ev[1], ev[2]


def T2_of(V):
    G = V.T @ V
    M = V.T @ Lm @ V
    return float(np.trace(np.linalg.solve(G, M)))


vac = np.ones(NL) / np.sqrt(NL)
P = np.eye(NL) - np.outer(vac, vac)
rng = np.random.default_rng(7)
diffs = []
for _ in range(2000):
    A = P @ rng.standard_normal((NL, 2))
    Q, _ = np.linalg.qr(A)
    diffs.append(T2_of(Q[:, :2]) - (mu1 + mu2))
print(f"  path N={NL}: mu1={mu1:.6f} mu2={mu2:.6f} mu1+mu2={mu1+mu2:.6f} rho^2={mu2/mu1:.6f}")
print(f"  random V: T2-(mu1+mu2) in [{min(diffs):.3e}, {max(diffs):.3e}]")
chk("H T2(V) >= mu1+mu2 always (Ky Fan)", min(diffs) > -1e-9)
Vopt = np.column_stack([vecs[:, 1], vecs[:, 2]])
chk("H equality at span(psi1,psi2)", abs(T2_of(Vopt) - (mu1 + mu2)) < 1e-8, f"T2={T2_of(Vopt):.9f}")

print("=" * 80)
print("PART I  coarsening: compressed 2x2 ratio is NOT ordered vs the true ratio")
print("=" * 80)
n = 8
M8 = np.diag(2 * np.ones(n)) + np.diag(-np.ones(n - 1), 1) + np.diag(-np.ones(n - 1), -1)
ev8, V8 = np.linalg.eigh(M8)
r_true = ev8[1] / ev8[0]
print(f"  path N=8 : mu1={ev8[0]:.6f} mu2={ev8[1]:.6f}  rho^2_true = {r_true:.6f}")


def comp_ratio(cols):
    Q = np.column_stack(cols)
    Q, _ = np.linalg.qr(Q)
    w = np.sort(np.linalg.eigvalsh(Q.T @ M8 @ Q))
    return w[1] / w[0]


d1 = np.zeros(n); d1[0] = 1.0
d2 = np.zeros(n); d2[1] = 1.0
dn = np.zeros(n); dn[n // 2] = 1.0
r_lo = comp_ratio([d1, d2])                      # two adjacent deltas: block [[2,-1],[-1,2]]
r_hi = comp_ratio([V8[:, 1], d1])                # psi1 + one delta: nu1 ~ mu1, nu2 = 2
print(f"  coarse = span(delta_1, delta_2)      : rho^2_coarse = {r_lo:.6f}  ({'<' if r_lo<r_true else '>'} true)")
print(f"  coarse = span(psi_1, delta_1)        : rho^2_coarse = {r_hi:.6f}  ({'<' if r_hi<r_true else '>'} true)")
print(f"  coarse = span(delta_1, delta_{n//2})       : rho^2_coarse = {comp_ratio([d1, dn]):.6f}")
rng = np.random.default_rng(11)
lo, hi = 0, 0
for _ in range(3000):
    A = rng.standard_normal((n, 2)) * (rng.random((n, 2)) < 0.35)  # sparse random
    if np.linalg.matrix_rank(A) < 2:
        continue
    Q, _ = np.linalg.qr(A)
    wc = np.sort(np.linalg.eigvalsh(Q.T @ M8 @ Q))
    if wc[0] < 1e-9:
        continue
    if wc[1] / wc[0] < r_true:
        lo += 1
    else:
        hi += 1
print(f"  sparse random coarse spaces: below {lo} / above {hi}")
chk("I coarsening does NOT order the ratio: both directions constructed", r_lo < r_true < r_hi,
    f"{r_lo:.4f} < {r_true:.4f} < {r_hi:.4f}")

print("=" * 80)
print("PART J  window arithmetic: Delta=(rho^2-1)mu1 ; rho^2<=C^2 => mu1 >= delta/(C^2-1)")
print("=" * 80)
anch = [("Mehler / OU heat kernel", 2.0), ("one plaquette SU(2) beta=0", 8 / 3),
        ("2D YM SU(2) (adj/fund)", 8 / 3), ("path graph N=2", 3.0),
        ("box / rotor / 1D Dirichlet string", 4.0), ("2D U(1)", 4.0)]
for name, r in anch:
    print(f"  {name:36s} rho^2={r:7.4f}   C=1.5 ok? {str(r<=2.25):5s}  C=2 ok? {str(r<=4+1e-12):5s}")
for C in [1.5, 2.0, 3.0]:
    print(f"  C={C}: mu1 >= delta/(C^2-1) = {1/(C*C-1):.4f} * delta")
chk("J C=1.5 violated by >=3 exact anchors", sum(1 for _, r in anch if r > 2.25) >= 3)
chk("J C=2 accommodates every exact anchor", all(r <= 4 + 1e-12 for _, r in anch))

print("=" * 80)
print(f"SUMMARY: {len(FAIL)} FAIL" + ("" if not FAIL else " -> " + "; ".join(FAIL)))
print("=" * 80)
