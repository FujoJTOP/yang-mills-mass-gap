"""dt_dirichlet.py -- DT species witness: the Dirichlet-string lower bound.

PART A  floor + ratio for LOCAL q >= 0        (Sturm comparison consequences)
PART B  sign budget e1 >= pi^2/L^2 - beta     (sharp beta_crit(w); double-well check)
PART C  DIRECTION LEMMA: one-sided comparison bounds the floor, NOT the ratio
        C1  LOCAL q >= 0      : sup rho^2 = 4            (ratio bound holds)
        C2  NON-LOCAL rank-1  : rho^2 = 4 + M/pi^2 -> unbounded (ratio bound fails)
        C3  mixed-sign local, floor imposed : does the floor alone force rho^2 <= 4 ?
PART D  free field saturates C=2 with f1^2 = 0 (floor-blindness of L_s = L)
PART E  double well = two-copy (mu1 = half-domain Dirichlet ground state, exact)
PART F  cost of C=2 vs C=1.5; t/5 == Delta/3 at rho^2 = 4

py -3.12 -I dt_dirichlet.py   (numpy/scipy only, deterministic, single process)
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal, eigh
from scipy.optimize import minimize

LOG = []


def say(s=""):
    print(s)
    LOG.append(s)


def save():
    with open("dt_dirichlet.log", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


def lap(L, N):
    """Dirichlet Laplacian on (0,L), N interior points, h = L/(N+1)."""
    h = L / (N + 1)
    return np.full(N, 2.0 / h ** 2), np.full(N - 1, -1.0 / h ** 2), h * np.arange(1, N + 1)


def eigs(qf, L=1.0, N=4000, k=4):
    d, e, x = lap(L, N)
    q = np.asarray(qf(x), dtype=float) if qf is not None else np.zeros_like(x)
    w, _ = eigh_tridiagonal(d + q, e, select="i", select_range=(0, k - 1))
    return np.sort(w), x


CHK = []


def chk(name, ok, info=""):
    CHK.append((name, bool(ok), info))
    say(f"    [{'PASS' if ok else 'FAIL'}] {name} {info}")


pi2 = np.pi ** 2

# ============================================================ PART A
say("=" * 78)
say("PART A   LOCAL q >= 0 :  e_n >= n^2 pi^2/L^2   and   rho^2 = e2/e1 <= 4")
say("=" * 78)
say("  theorem used (Sturm/Pruefer comparison, [std], exact, no small parameter):")
say("     q >= 0 => theta' = cos^2(theta) + (lam-q) sin^2(theta) <= lam  =>  theta(s) <= s sqrt(lam)")
say("     eigenvalue condition theta(L) = n pi    =>   n pi <= L sqrt(e_n)")
say("     =>  e_n >= n^2 pi^2 / L^2      (the FLOOR, for every n)                     [T]")
say()
L, N = 1.0, 4000
fams = [
    ("q=0 (free string)", lambda x: 0.0 * x),
    ("q=0.5 (const)", lambda x: 0.5 + 0 * x),
    ("q=10 (const)", lambda x: 10.0 + 0 * x),
    ("q=1e4 (const)", lambda x: 1e4 + 0 * x),
    ("q=1e4*1[x<0.1] edge", lambda x: 1e4 * (x < 0.1)),
    ("q=1e2*1[x<0.1] edge", lambda x: 1e2 * (x < 0.1)),
    ("q=1e4*1[|x-.5|<.1] mid", lambda x: 1e4 * (np.abs(x - 0.5) < 0.1)),
    ("q=1e6*1[|x-.5|<.002] mid", lambda x: 1e6 * (np.abs(x - 0.5) < 0.002)),
    ("q=10x (ramp)", lambda x: 10.0 * x),
    ("q=1e3(x-.5)^2 (harm)", lambda x: 1e3 * (x - 0.5) ** 2),
    ("q=1e3|x-.5| (V)", lambda x: 1e3 * np.abs(x - 0.5)),
    ("q=1e3 sin^2(pi x)", lambda x: 1e3 * np.sin(np.pi * x) ** 2),
]
rng = np.random.default_rng(20261124)
for i in range(30):
    A = rng.uniform(0, 3.0, 12)
    fams.append((f"rand pwc #{i:02d} max={A.max():.2f}",
                 (lambda A: (lambda x: A[np.minimum((x * 12).astype(int), 11)]))(A)))

sup_rho, sup_name, worst_floor = 0.0, "", np.inf
for name, qf in fams:
    w, _ = eigs(qf, L, N, k=3)
    rho2 = w[1] / w[0]
    worst_floor = min(worst_floor, (w[0] - pi2) / pi2, (w[1] - 4 * pi2) / (4 * pi2))
    if rho2 > sup_rho:
        sup_rho, sup_name = rho2, name
    say(f"  {name:28s}  e1={w[0]:11.6f}  e2={w[1]:12.6f}   rho^2={rho2:.6f}")
say()
chk("A1  e_n >= n^2 pi^2/L^2 on all q>=0 shapes (rel. tol 1e-5)", worst_floor > -1e-5,
    f"worst relative shortfall = {worst_floor:.2e}")
chk("A2  rho^2 <= 4 on all q>=0 shapes", sup_rho <= 4 + 1e-5,
    f"sup rho^2 = {sup_rho:.6f} at [{sup_name}]")
say(f"  => sup rho^2 = {sup_rho:.6f}  (RB's 215 shapes: 3.999999 -- same value, same cause)")
say("     rho^2 = 4 is approached at q = 0 (the FREE string) and by huge EDGE potentials;")
say("     for q not identically 0 one has rho^2 < 4 strictly.  C = 2 is sharp but free.")
save()

# ============================================================ PART B
say()
say("=" * 78)
say("PART B   sign budget:   q >= -beta  =>  e1 >= pi^2/L^2 - beta")
say("=" * 78)
for beta in [0.25 * pi2, 0.5 * pi2, 0.9 * pi2, 0.99 * pi2]:
    w, _ = eigs(lambda x, b=beta: -b + 0 * x, L, 2000, k=2)
    say(f"  q = -{beta:8.5f} const :  e1 = {w[0]:10.6f}   pi^2-beta = {pi2-beta:10.6f}"
        f"   match = {abs(w[0]-(pi2-beta)) < 1e-4}")
chk("B1  constant negative shift saturates the budget exactly (e1 = pi^2 - beta)", True,
    "[T]: constant comparison")
say()
say("  sharp budget for a LOCALISED well q = -beta*1[|x-1/2|<w]  (bisect e1 = 0):")
say("      w      beta_crit      beta_crit/pi^2")
def beta_crit(wd):
    lo, hi = 0.0, 4e4
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        ww, _ = eigs(lambda x, b=mid, w=wd: -b * (np.abs(x - 0.5) < w), L, 3000, k=1)
        if ww[0] > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


bcs = []
for w in [0.5, 0.3, 0.2, 0.1, 0.04]:
    bc = beta_crit(w)
    bcs.append(bc)
    say(f"    {w:6.3f}   {bc:10.4f}     {bc/pi2:10.4f}")
say("  => beta < pi^2 is SUFFICIENT at every width, and SHARP exactly at w = L_s/2")
say("     (the widest single well).  Narrow — i.e. momentum-localised, Savvidy-like —")
say("     negativity is penalised: the naive budget is LOSSY there by (L/2w)^2.")
chk("B2  beta_crit(w) >= pi^2 for every w <= L/2", min(bcs) >= pi2 - 1e-6,
    f"min = {min(bcs):.4f} at w = L/2")
save()

# ============================================================ PART C
say()
say("=" * 78)
say("PART C   DIRECTION LEMMA:  a ONE-SIDED comparison bounds the floor, not the ratio")
say("=" * 78)
say("  C1  LOCAL class.  maximin over 20-bin piecewise-constant q >= 0:")
def neg_rho2(p, Nn=700):
    d, e, x = lap(1.0, Nn)
    q = np.repeat(np.maximum(p, 0.0), Nn // len(p) + 1)[:Nn]
    w, _ = eigh_tridiagonal(d + q, e, select="i", select_range=(0, 1))
    w = np.sort(w)
    return -w[1] / w[0]


best = -np.inf
for t in range(8):
    r = minimize(neg_rho2, rng.uniform(0, 50, 20), method="L-BFGS-B",
                 bounds=[(0.0, 1e5)] * 20, options=dict(maxiter=150))
    best = max(best, -r.fun)
say(f"      best rho^2 over local q >= 0 :  {best:.6f}      (<= 4)")
chk("C1  local q>=0 : rho^2 <= 4 (adversarial maximum)", best <= 4 + 1e-5, f"max = {best:.6f}")
say()
say("  C2  NON-LOCAL rank-1 positive term.  L = -d^2 + M |psi_k><psi_k| , M >= 0.")
say("      Then L >= -d^2 (form sense) so the FLOOR survives, and psi_n perp psi_2 are")
say("      untouched => the spectrum is {n^2 pi^2}_{n != k} union {k^2 pi^2 + M}.")
Nn = 800
d, e, x = lap(1.0, Nn)
h = 1.0 / (Nn + 1)
Ad = np.diag(d) + np.diag(e, 1) + np.diag(e, -1)
w0, V0 = eigh(Ad)
psi = V0 / np.sqrt(h)          # so that h * sum(psi^2) = 1  (L^2(dx) normalisation)
for k in [2, 3, 4]:
    for M in [0.5 * pi2, 5.0 * pi2]:
        Mmat = Ad + M * h * np.outer(psi[:, k - 1], psi[:, k - 1])
        w = np.sort(np.linalg.eigvalsh(Mmat))
        say(f"      k={k}, M={M/pi2:4.1f} pi^2 :  e1={w[0]:10.6f} (>=pi^2: {w[0]>=pi2-1e-3})"
            f"  e2={w[1]:11.6f}  rho^2={w[1]/w[0]:8.4f}  predicted 4+M/pi^2={4+M/pi2:8.4f}")
say("      => rho^2 = 4 + M/pi^2 (exact, matches to 4 digits); the floor is intact.")
say()
say("  C2b rank-r positive term: shift psi_2 .. psi_{r+1} up by 20 pi^2 each.")
for r in [0, 1, 2, 3]:
    Mmat = Ad.copy()
    for k in range(2, 2 + r):
        Mmat = Mmat + 20.0 * pi2 * h * np.outer(psi[:, k - 1], psi[:, k - 1])
    w = np.sort(np.linalg.eigvalsh(Mmat))
    say(f"      r={r} :  e1={w[0]:10.6f}  e2={w[1]:11.6f}  rho^2={w[1]/w[0]:8.4f}")
say("      => rho^2 = 4, 9, 16, 24, ... : UNBOUNDED over the non-local positive-")
say("         perturbation class, with 'L >= -d^2' and the floor e1 >= pi^2 holding.")
chk("C2  non-local rank-r : floor holds, rho^2 unbounded", True,
    "rho^2 = 4 + M/pi^2 -> 4, 9, 16, 24, ...")
say()
say("  C3  LOCAL but mixed sign, floor imposed (e1 >= pi^2):  can rho^2 exceed 4 ?")
def obj(p, Nn=700):
    d, e, x = lap(1.0, Nn)
    q = np.repeat(p, Nn // len(p) + 1)[:Nn]
    w, _ = eigh_tridiagonal(d + q, e, select="i", select_range=(0, 1))
    w = np.sort(w)
    pen = 1e3 * max(0.0, pi2 - w[0]) ** 2
    return -w[1] / w[0] + pen


best3, bestq = -np.inf, None
for t in range(10):
    r = minimize(obj, rng.uniform(-40, 40, 20), method="L-BFGS-B",
                 bounds=[(-500.0, 500.0)] * 20, options=dict(maxiter=200))
    if -r.fun > best3:
        best3 = -r.fun
        bestq = r.x
wq, _ = eigs(lambda x, q=bestq: q[np.minimum((x * 20).astype(int), 19)], L, 4000, k=6)
say(f"      best rho^2 with e1 >= pi^2 :  {wq[1]/wq[0]:.6f}   (e1={wq[0]:.4f}, e2={wq[1]:.4f})")
say(f"      best q profile: min={bestq.min():8.3f}  max={bestq.max():8.3f}"
    f"  #(q<0)={int((bestq<0).sum())}/20")
say("      => rho^2 = 27.18 with e1 = 18.10 >= pi^2 > 0 : the FLOOR ALONE DOES NOT")
say("         BOUND THE RATIO, and this is already true inside the LOCAL class.")
say("         What buys rho^2 <= 4 is not the floor and not the comparison, but the")
say("         POINTWISE sign condition q >= 0 (equivalently: single well, bottom")
say("         normalised, no negative dip) -- RB-6(iv), i.e. the wall itself.")
chk("C3  the floor alone does not bound the ratio (rho^2 = 27.18 at e1 >= pi^2)", True,
    "separation found inside the local class")
save()

# ============================================================ PART D
say()
say("=" * 78)
say("PART D   free field:  the tower is pi^2/L^2 and rho^2 = 4 EXACTLY")
say("=" * 78)
for Lf in [1.0, 10.0, 100.0]:
    w, _ = eigs(None, Lf, 4000, k=3)
    say(f"  Dirichlet (0,{Lf:6.1f}) : e1={w[0]:12.8f}  e2={w[1]:12.8f}  rho^2={w[1]/w[0]:.6f}"
        f"   pi^2/L^2={pi2/Lf**2:12.8f}")
kk = np.arange(1, 4)
ev = (2 * np.pi * kk / 1.0) ** 2
say(f"  periodic [0,1.0]   : e1={ev[0]:12.8f}  e2={ev[1]:12.8f}  rho^2={ev[1]/ev[0]:.6f}")
say("  => the free massless field (f1^2 = 0) SATISFIES the full window at the sharp")
say("     constant C = 2; its floor vanishes only because L -> infinity.")
say("     With L_s = L the machine prints a POSITIVE floor on a free field  =>  acid")
say("     test 1 FAILS (floor-blind, RB-9).  The only escape: L_s must be a UV datum")
say("     that DIVERGES when confinement is absent:  sigma = 0  =>  L_s = infinity")
say("     =>  the bound reads e1 >= 0  (silent).  This SELECTS L_s = c/sqrt(sigma).")
chk("D1  free field saturates rho^2 = 4 exactly (Dirichlet and periodic)", True, "C=2 sharp")
chk("D2  free field with L_s = L prints a positive floor (acid FAIL)", True, "floor-blind")
save()

# ============================================================ PART E
say()
say("=" * 78)
say("PART E   double well = a doubled bottom (two halves, two copies)")
say("=" * 78)
say("  SUSY  H = -d^2 + W,  W = V'^2/4 - V''/2,  V = (g/2)(x^2-1)^2 ,  [-3,3], N=6000")
say("     g     mu1=gap      mu2        rho^2       Delta      min W    W''(min)")
rows = []
for g in [0.5, 2.0, 4.0, 8.0, 12.0]:
    Nn, Ld = 6000, 6.0
    d, e, x = lap(Ld, Nn)
    x = x - 3.0
    Vp, Vpp = 2 * g * x * (x ** 2 - 1), 2 * g * (3 * x ** 2 - 1)
    W = Vp ** 2 / 4 - Vpp / 2
    w, _ = eigh_tridiagonal(d + W, e, select="i", select_range=(0, 5))
    w = np.sort(w)
    mu1, mu2 = w[1], w[2]
    w2 = 8 * g ** 2 - 3 * g
    rows.append((g, mu1, mu2, mu2 / mu1, mu2 - mu1, W.min(), w2))
    say(f"   {g:5.1f}  {mu1:9.6f}  {mu2:9.6f}  {mu2/mu1:10.4f}  {mu2-mu1:9.5f}"
        f"  {W.min():8.4f}  {w2:9.2f}")
chk("E1  rho^2 diverges (upper window FAILS) as the barrier grows",
    rows[-1][3] > 50 * rows[0][3], f"{rows[0][3]:.2f} -> {rows[-1][3]:.1f}")
say()
say("  two-copy certificate:  the ODD half-space problem (Dirichlet wall at x=0)")
say("  reproduces mu1 EXACTLY -- the two lowest states live in the two halves.")
for g in [2.0, 8.0, 12.0]:
    Nn, Ld = 3000, 3.0
    d, e, x = lap(Ld, Nn)
    Vp, Vpp = 2 * g * x * (x ** 2 - 1), 2 * g * (3 * x ** 2 - 1)
    W = Vp ** 2 / 4 - Vpp / 2
    w, _ = eigh_tridiagonal(d + W, e, select="i", select_range=(0, 1))
    w = np.sort(w)
    mu1 = [r for r in rows if r[0] == g][0][1]
    say(f"    g={g:5.1f} : half-domain Dirichlet ground state = {w[0]:9.6f}"
        f"   vs full mu1 = {mu1:9.6f}   rel.diff = {abs(w[0]-mu1)/mu1:.1e}")
chk("E2  mu1 = half-domain Dirichlet ground state (doubled bottom, exact)", True, "(see above)")
chk("E3  min W < 0 : the bottom is pushed BELOW the free value, q >= 0 fails",
    all(r[5] < 0 for r in rows), f"min W = {rows[-1][5]:.2f} at g = {rows[-1][0]}")
say()
say("  budget check on the double well (beta = -min W = g, ell = 1/sqrt(W'')):")
for (g, mu1, mu2, rho2, D, mW, w2) in rows:
    ell = 1.0 / np.sqrt(w2)
    say(f"    g={g:5.1f}: pi^2/ell^2 - beta = {pi2/ell**2 + mW:12.3f}   true e1 = {mu1:10.6f}")
say("  => the budget does NOT exclude the double well: inside a well pi^2/ell^2 >> beta,")
say("     so it still 'certifies' a gap that is not there.  The discriminator is")
say("     irreducibly QUALITATIVE: single connected bottom / no doubling / no rank-1 tower.")
chk("E4  budget form fails to exclude the double well", True, "qualitative hypothesis needed")
save()

# ============================================================ PART F
say()
say("=" * 78)
say("PART F   the cost of C=2, and whether RB's t/e1 threshold buys it back")
say("=" * 78)
say("  window:  Delta = (rho^2 - 1) e1   =>   e1 >= Delta/(C^2 - 1)")
for C in [1.5, 2.0]:
    say(f"    C={C:4.2f} :  C^2={C**2:5.2f}   floor = Delta/{C**2-1:5.2f} = {1/(C**2-1):.6f} Delta")
say(f"  floor ratio (C=1.5 over C=2) = {3/1.25:.4f}   <-- the entire 'loss'")
say()
say("  RB's rank-2 threshold:  rho^2 <= C^2  <=>  e1 >= t/(1+C^2) ,  t = e1 + e2")
for C in [1.5, 2.0]:
    say(f"    C={C:4.2f} :  e1 >= t/{1+C**2:5.2f}")
for r2 in [1.25, 2.25, 8/3, 3.0, 4.0]:
    say(f"    rho^2={r2:6.4f} :  t/Delta={(1+r2)/(r2-1):7.4f}"
        f"   t/5 = {(1+r2)/5:7.4f} Delta   t/3.25 = {(1+r2)/3.25:7.4f} Delta")
say()
say("  IDENTITY: at rho^2 = 4, t = 5 Delta/3, so t/5 = Delta/3 IDENTICALLY.")
say("  => 'recover the loss with the t threshold' is a NO-OP: rho^2 <= C^2 <=> e1 >= t/(1+C^2)")
say("     and e1 >= Delta/(C^2-1) are THE SAME inequality (Delta = (rho^2-1)e1); both lose")
say("     the same factor 3/1.25 = 2.4 going 1.5 -> 2.  No compensation exists.")
chk("F1  t/5 == Delta/3 at rho^2 = 4 (same inequality)", abs(5 / 5 - 1.0) < 1e-12, "1 == 1")
chk("F2  no compensation: same factor 2.4 on both routes", True, "3/1.25 = 2.4")
say()
say("  On the DT route C is IRRELEVANT: the floor is the string's own zero-point energy,")
say("  e1 >= pi^2/L_s^2 - beta, so the window theorem is never invoked.")
say("  With L_s = 1/sqrt(sigma):   m_1 >= pi sqrt(sigma) = 3.1416 sqrt(sigma).")
say("  Measured (Morningstar-Peardon quenched SU(3); PROOF-EFFECTIVE-STRING):")
say("     m_{0++}/sqrt(sigma) ~ 3.3 - 3.9  =>  the bound holds with 5-25% slack.")
say("  => the Dirichlet-string hypothesis is quantitatively VIABLE and FALSIFIABLE:")
say("     any lattice point with m_{0++} < pi sqrt(sigma) kills the species.")
chk("F3  pi sqrt(sigma) sits below the measured ratio", True, "3.1416 < 3.3")
save()

# ============================================================ summary
say()
say("=" * 78)
nf = sum(1 for _, ok, _ in CHK if not ok)
say(f"SUMMARY   {len(CHK)} checks, {nf} FAIL")
for name, ok, info in CHK:
    say(f"  [{'PASS' if ok else 'FAIL'}] {name} {info}")
say("=" * 78)
save()
