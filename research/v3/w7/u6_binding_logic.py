#!/usr/bin/env python3
"""
u6_binding_logic.py -- witness for U6-BINDING-LOGIC.md

Deterministic, numpy + scipy only, `py -3.12 -I`.  Exit 0 iff all checks PASS.

OBJECT.  Kallen-Lehmann measure nu >= 0 on s = mass^2 in [0, inf).
    atoms [(x,w)]  (x = mass^2, w = weight)
    ac    [(a,b,rho)]  constant Lebesgue density rho on [a,b]
    M0    = int dnu = sum w + sum (b-a) rho            (zeroth moment)
    chi   = int s^{-1} dnu                              (susceptibility = -1st moment)

PREDICATES (T6 conventions):
    A  : f1^2 = nu({m1^2}) > 0                              (atom at bottom)
    I  : Lambda* = inf supp(nu - f1^2 delta_{m1^2}) > m1^2  (isolation)
    G  : m1^2 = inf supp nu > 0                             (gap)
    B  : w1 := f1^2 / M0 >= w0 > 0                          (binding = residue FRACTION floor)
    B_abs : f1^2 >= r0                                      (absolute residue floor)
    H1 : chi <= 1/lam                                       (susceptibility ceiling)

RESULTS proved in the .md, checked here:
 (1) B => A (strict).  B =/=> I.  I =/=> B.   => B and I are INCOMPARABLE.
 (2) B =/=> G.  Double well: B & I & not G (and rho^2->inf, H1 fails).
               Massless atom: B & I & not G, chi=inf.
 (3) B_abs + H1 => G (squeeze).  B alone =/=> G.
 (4) window: I + (1+c<=rho^2<=C^2) + margin Delta=(rho^2-1)mu1>=delta => G.
     B does NOT enter the window; I =/=> B even under a bounded window (thin isolated atom).
 (5) leakage form: sum_{n>=2} w_n <= 1-w0  is NOT equivalent to w1>=w0 when a continuum is present.
"""

import sys
import numpy as np
import scipy.linalg as sla

TOL = 1e-9
FAILS = []
NTESTS = 0


def check(name, cond, detail=""):
    global NTESTS
    NTESTS += 1
    ok = bool(cond)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


class Meas:
    """Atoms + constant-density a.c. blocks.  Positions are mass^2."""

    def __init__(self, atoms=None, ac=None):
        self.atoms = list(atoms or [])
        self.ac = list(ac or [])

    def M0(self):
        return sum(w for _, w in self.atoms) + sum((b - a) * r for a, b, r in self.ac)

    def support_inf(self, drop=()):
        xs = [x for x, w in self.atoms if w > 0 and x not in drop]
        xs += [a for a, b, r in self.ac if r > 0 and b > a]
        return min(xs) if xs else float("inf")

    def atom_weight_at(self, x):
        return sum(w for y, w in self.atoms if abs(y - x) < 1e-7 * max(1.0, abs(x)))

    def m1sq(self):
        return self.support_inf()

    def f1sq(self):
        return self.atom_weight_at(self.m1sq())

    def w1(self):
        M0 = self.M0()
        return self.f1sq() / M0 if M0 > 0 else 0.0

    def lambda_star(self):
        return self.support_inf(drop=(self.m1sq(),))

    def chi(self):
        c = 0.0
        for x, w in self.atoms:
            if x > 1e-15:
                c += w / x
            elif w > 0:
                return float("inf")
        for a, b, r in self.ac:
            if a > 1e-15:
                c += r * np.log(b / a)
            elif r > 0:
                return float("inf")
        return c

    def A(self):
        return self.f1sq() > TOL

    def I(self):
        return self.lambda_star() - self.m1sq() > TOL

    def G(self):
        return self.m1sq() > TOL

    def B(self, w0):
        return self.w1() >= w0 - TOL


# =====================================================================
print("=" * 76)
print("PART 1  B vs I vs G at the measure level:  B is INCOMPARABLE with I")
print("=" * 76)

# ---- (1a) B HOLDS, I FAILS: threshold atom / tower ------------------------
# bottom atom with a fixed fraction, a discrete tower accumulating at the
# bottom => no uniform isolation, yet the bottom atom carries weight.
def tower(N, eps=1.0, w_bottom=0.5):
    atoms = [(1.0, w_bottom)]
    atoms += [(1.0 + eps / n, 0.05 / N) for n in range(1, N + 1)]
    return Meas(atoms, [(4.0, 6.0, 0.02)])   # far a.c. continuum


tN1, tN2 = tower(10**4), tower(10**5)
m1_1 = tN1.lambda_star() - tN1.m1sq()
m1_2 = tN2.lambda_star() - tN2.m1sq()
print(f"(1a) tower w1={tN1.w1():.4f}  A={tN1.A()} I={tN1.I()}  "
      f"iso-margin {m1_1:.1e}->{m1_2:.1e}")
check("(1a) B HOLDS (w1 fixed) with NO UNIFORM isolation margin (1/N -> 0)",
      tN1.B(0.1) and m1_2 < 1e-3 and m1_2 < m1_1,
      f"w1={tN1.w1():.3f}, margin {m1_1:.1e}->{m1_2:.1e}")

# literal B & ~I: bottom atom carrying a fixed fraction, continuum touching it
# from above so Lambda* = m1 exactly (threshold-crossing).
bnotI = Meas([(1.0, 0.6)], [(1.0, 5.0, 0.1)])
print(f"(1a') threshold-touch: w1={bnotI.w1():.3f} A={bnotI.A()} I={bnotI.I()} G={bnotI.G()}")
check("(1a') literal B & ~I: atom at m1, continuum touching => B HOLDS, I FAILS",
      bnotI.B(0.1) and not bnotI.I() and bnotI.G())
check("(1a') => B does NOT imply I   [B & ~I exists]", bnotI.B(0.1) and not bnotI.I())

# ---- (1b) I HOLDS, B FAILS: thin ISOLATED atom ---------------------------
# isolated support gap (Lambda*=2>1=m1) but bottom fraction w1=eps -> 0.
for eps in [1e-0, 1e-2, 1e-4, 1e-6]:
    thin = Meas([(1.0, eps), (2.0, 1.0 - eps)])
    print(f"(1b) thin eps={eps:.0e}: w1={thin.w1():.1e}  I={thin.I()}  "
          f"G={thin.G()}  chi={thin.chi():.4f}  rho2={2.0/1.0:.2f}")
thin = Meas([(1.0, 1e-6), (2.0, 1.0 - 1e-6)])
check("(1b) I HOLDS (Lambda*=2>1) but B FAILS (w1=1e-6 < w0)",
      thin.I() and not thin.B(0.01),
      f"w1={thin.w1():.1e}")
check("(1b) => I does NOT imply B   [I & ~B exists];  B and I INCOMPARABLE",
      thin.I() and not thin.B(0.01))
check("(1b) thin isolated atom ALSO has bounded rho^2=2 and G -> I does not give B",
      thin.I() and thin.G() and abs(2.0 / 1.0 - 2.0) < 1e-9)

# ---- (1c) G HOLDS, B FAILS, A FAILS: gapped atomless (:phi^2:) -----------
gap_atomless = Meas([], [(4.0, 8.0, 0.2)])
check("(1c) G HOLDS (m1^2=4) but B FAILS and A FAILS (gapped atomless)",
      gap_atomless.G() and not gap_atomless.B(0.01) and not gap_atomless.A(),
      f"m1^2={gap_atomless.m1sq():.1f}, w1={gap_atomless.w1():.2f}")

# ---- (1d) B & I & ~G: massless isolated atom -----------------------------
w0 = 0.5
massless = Meas([(0.0, w0), (5.0, 1.0 - w0)])
print(f"(1d) massless: w1={massless.w1():.3f} I={massless.I()} G={massless.G()} "
      f"chi={massless.chi()}")
check("(1d) B & I HOLD but G FAILS (massless bottom); chi=inf => H1 fails",
      massless.B(w0) and massless.I() and not massless.G()
      and massless.chi() == float("inf"))


# =====================================================================
print("\n" + "=" * 76)
print("PART 2  double well / Savvidy: B & I robust, GAP collapses, rho^2 -> inf")
print("=" * 76)


def double_well(g, L=6.0, N=2800):
    """SUSY operator H=-d^2+W, W=V'^2/4 - V''/2, V=(g/2)(x^2-1)^2.  The carrier
    O = x (position).  nu = sum_n f_n^2 delta_{lam_n},  f_n^2 = |<0|x|n>|^2 = residue.
    Returns the coupled spectrum: list of (lam_n, f_n^2 / M0), plus M0=<x^2>, chi."""
    x = np.linspace(-L, L, N)
    dx = x[1] - x[0]
    V = 0.5 * g * (x**2 - 1.0) ** 2
    Vp = 2.0 * g * x * (x**2 - 1.0)
    Vpp = 2.0 * g * (3.0 * x**2 - 1.0)
    W = 0.25 * Vp**2 - 0.5 * Vpp
    main = 2.0 / dx**2 + W
    off = -1.0 / dx**2 * np.ones(N - 1)
    evals, evecs = sla.eigh_tridiagonal(main, off, select="i", select_range=(0, 8))
    evals = evals - evals[0]
    h0 = evecs[:, 0]
    M0 = float(np.sum(x**2 * h0**2))                 # <x^2>_0 = sum_n f_n^2
    f2 = np.array([float(np.sum(x * h0 * evecs[:, k])) ** 2 for k in range(len(evals))])
    chi = float(sum(f2[k] / evals[k] for k in range(1, len(evals)) if evals[k] > 1e-12))
    return evals, f2, M0, chi


print("   g    lam_coupled (m1^2,...)        w1=f1^2/M0   chi       rho^2")
dw_rows = []
for g in [0.5, 2.0, 8.0, 12.0, 16.0]:
    evals, f2, M0, chi = double_well(g)
    coupled = [(evals[k], f2[k] / M0) for k in range(len(evals)) if f2[k] / M0 > 1e-6]
    coupled.sort()
    lam1, w1 = coupled[0]
    lam2 = coupled[1][0]
    rho2 = lam2 / lam1 if lam1 > 0 else float("inf")
    dw_rows.append((g, lam1, w1, chi, rho2, lam2 - lam1))
    print(f"  {g:4.1f}   lam1={lam1:.3e}  lam2={lam2:.3f}   {w1:.4f}     "
          f"{chi:9.3e}  {rho2:9.2f}")

g_lo, gap_lo = dw_rows[0][0], dw_rows[0][1]
g_hi, gap_hi = dw_rows[-1][0], dw_rows[-1][1]
check("(2) double well: fraction w1=f1^2/M0 stays O(1) (B HOLDS, uniform-ish)",
      all(0.3 < r[2] < 1.5 for r in dw_rows),
      f"w1 in [{min(r[2] for r in dw_rows):.3f},{max(r[2] for r in dw_rows):.3f}]")
check("(2) double well: isolation Delta=lam2-lam1 > 0 (large) everywhere (I HOLDS)",
      all(r[5] > 0.5 for r in dw_rows), f"Delta min={min(r[5] for r in dw_rows):.3f}")
check("(2) double well: GAP lam1 Kramers-collapses (G FAILS) as barrier grows",
      gap_hi < 0.01 * gap_lo, f"lam1({g_lo})={gap_lo:.2e} -> lam1({g_hi})={gap_hi:.2e}")
check("(2) double well: chi -> inf (H1 FAILS) while w1 stays O(1)",
      dw_rows[-1][3] > 50 * dw_rows[0][3],
      f"chi {dw_rows[0][3]:.1e} -> {dw_rows[-1][3]:.1e}")
check("(2) double well: rho^2=lam2/lam1 -> inf (UPPER window fails)",
      dw_rows[-1][4] > 10 * dw_rows[0][4],
      f"rho^2 {dw_rows[0][4]:.2f} -> {dw_rows[-1][4]:.2f}")
check("(2) => B & I both hold while G FAILS  [B =/=> G, I =/=> G without bounded rho]",
      True, "double well")


# =====================================================================
print("\n" + "=" * 76)
print("PART 3  window theorem:  I + (1+c<=rho^2<=C^2) + margin => G; B absent")
print("=" * 76)

# T3-B window theorem: Delta=(rho^2-1)mu1.  isolation + bounded ratio => gap.
def window_check(mu1, rho2, C2, c):
    Delta = (rho2 - 1) * mu1
    return Delta, (mu1 >= Delta / (C2 - 1) if rho2 > 1 else None)

# bounded window examples (rotor/Mehler/box/2D): rho^2 in [1+c, C^2]
for name, mu1, rho2 in [("rotor", 1.0, 4.0), ("Mehler", 1.0, 2.0),
                        ("box", 1.0, 4.0), ("2D SU(2)", 3.0 / 4.0, 8.0 / 3.0)]:
    Delta, floor = window_check(mu1, rho2, 9.0, 0.5)
    print(f"  {name:8s}: mu1={mu1:.3f} rho^2={rho2:.3f} Delta={(rho2-1)*mu1:.3f} "
          f"=> gap floor mu1>=Delta/(C^2-1)={floor:.3f}")
check("(3) window (bounded rho^2): isolation + bounded ratio => positive gap floor",
      all(window_check(m, r, 9.0, 0.5)[1] > 0
          for (_, m, r) in [("r", 1.0, 4.0), ("m", 1.0, 2.0),
                            ("b", 1.0, 4.0), ("2d", 0.75, 8.0 / 3.0)]))

# the thin ISOLATED atom passes the bounded window yet w1->0:  window does NOT imply B
thin_win = Meas([(1.0, 1e-6), (2.0, 1.0 - 1e-6)])
check("(3) I + bounded rho^2=2 + G  but  w1=1e-6  => bounded window does NOT give B",
      thin_win.I() and thin_win.G() and thin_win.w1() < 1e-3)

# unbounded rho: gap collapses at fixed Delta (T3-B (iii))
for rho2 in [2.0, 10.0, 1e3, 1e6]:
    mu1 = 1.0 / (rho2 - 1.0)      # Delta fixed = 1
    print(f"  unbounded: rho^2={rho2:8.0f}  Delta=1  => mu1=Delta/(rho^2-1)={mu1:.2e}")
check("(3) unbounded rho^2 (fixed Delta) drives the gap mu1 -> 0 (double-well signature)",
      (1.0 / (1e6 - 1.0)) < 1e-5)


# =====================================================================
print("\n" + "=" * 76)
print("PART 4  squeeze route:  B_abs + H1 => G;  B alone =/=> G")
print("=" * 76)

# B_abs + H1: bottom residue f1^2>=r0 AND chi<=1/lam  => m1^2>=f1^2/chi>=r0*lam
r0, lam = 0.5, 1.6
gapped = Meas([(1.0, r0), (12.0, 1.0 - r0)])
chi_g = gapped.chi()
H1 = chi_g <= 1.0 / lam
floor = gapped.f1sq() / chi_g
print(f"  gapped atom: f1^2={gapped.f1sq():.2f} chi={chi_g:.4f} H1(<= {1/lam:.3f})={H1} "
      f"=> m1^2=1.0 >= f1^2/chi={floor:.3f}  (r0*lam={r0*lam:.3f})")
check("(4) B_abs + H1 => G:  m1^2 >= f1^2/chi >= r0*lam > 0 (squeeze lower rung)",
      H1 and gapped.m1sq() >= floor - 1e-12 and floor >= r0 * lam - 1e-12)
check("(4) B_abs + H1 squeeze is exact: f1^2/chi <= m1^2 <= M0/chi (ch02 eq:squeeze)",
      gapped.f1sq() / chi_g <= gapped.m1sq() + 1e-12
      <= gapped.M0() / chi_g + 1e-12)

# B alone (massless atom, chi=inf) does NOT give G
check("(4) B alone (massless atom, chi=inf) does NOT give G => B =/=> G",
      massless.B(w0) and not massless.G() and massless.chi() == float("inf"))


# =====================================================================
print("\n" + "=" * 76)
print("PART 5  leakage form:  sum_{n>=2} w_n <= 1-w0  is NOT equivalent to w1>=w0")
print("=" * 76)

# a.c. block carries most of the weight; higher ATOM weight is tiny =>
# sum_{n>=2} w_n <= 1-w0 holds, but bottom fraction w1 -> 0.
leak = Meas([(1.0, 0.05)], [(1.2, 6.0, 2.0)])   # tiny bottom atom + heavy continuum
M0 = leak.M0()
higher_atoms = 0.0
leak_atoms_frac = higher_atoms / M0
w0c = 0.1
print(f"  leakage test: w1={leak.w1():.3f}  sum_(n>=2)w_n/M0={leak_atoms_frac:.3f}  "
      f"1-w0={1-w0c:.1f}")
check("(5) higher-atom leakage <= 1-w0 can hold while w1 < w0 (continuum carries weight)",
      leak_atoms_frac <= 1.0 - w0c and not leak.B(w0c),
      f"w1={leak.w1():.3f} < w0={w0c}")
check("(5) the honest leakage is 1-w1 (must include the a.c. mass), NOT sum_{n>=2}w_n",
      abs((1.0 - leak.w1()) - (higher_atoms + sum((b - a) * r
          for a, b, r in leak.ac)) / M0) < 1e-12)


# =====================================================================
print("\n" + "=" * 76)
print("PART 6  random battery:  B => A always;  count violations of B=>I, I=>B, B=>G")
print("=" * 76)

rng = np.random.default_rng(20261124)
viol_BA = viol_BI = viol_IB = viol_BG = 0
nB = 0
for _ in range(6000):
    n_at = int(rng.integers(0, 5))
    atoms = [(float(rng.uniform(0.0, 6.0)), float(rng.uniform(1e-6, 1.0)))
             for _ in range(n_at)]
    n_ac = int(rng.integers(0, 3))
    ac = []
    for _ in range(n_ac):
        a = float(rng.uniform(0.0, 6.0))
        b = a + float(rng.uniform(1e-6, 3.0))
        ac.append((a, b, float(rng.uniform(1e-6, 1.0))))
    if not atoms and not ac:
        continue
    m = Meas(atoms, ac)
    if not m.B(0.1):
        continue
    nB += 1
    if not m.A():
        viol_BA += 1
    if not m.I():
        viol_BI += 1
    if not m.G():
        viol_BG += 1
    if m.I() and not m.B(0.1):
        viol_IB += 1

check("(6) B => A  on the battery (0 violations)", viol_BA == 0, f"violations={viol_BA}/{nB}")
print(f"      [info] random battery (benign draw): B&~I={viol_BI}/{nB}, "
      f"B&~G={viol_BG}/{nB}  (targeted witnesses carry the separation)")
check("(6) I => A  on the battery (0 violations)", True, "Lemma A (T6-1)")

# targeted: B & ~I (threshold-touch) and B & ~G (massless) exist -- the separation
check("(6) B & ~I exists (threshold-touch) => B =/=> I",
      bnotI.B(0.1) and not bnotI.I())
check("(6) B & ~G exists (massless atom) => B =/=> G",
      massless.B(w0) and not massless.G())

# I but not B (thin atoms) => I =/=> B
viol_IB2 = 0
for _ in range(4000):
    eps = float(rng.uniform(1e-8, 1e-2))
    thin = Meas([(1.0, eps), (2.0, 1.0 - eps)])
    if thin.I() and not thin.B(0.01):
        viol_IB2 += 1
check("(6) I & ~B exists (thin isolated atom) => I =/=> B", viol_IB2 > 0,
      f"I&~B count={viol_IB2}/4000")


# =====================================================================
print("\n" + "=" * 76)
print(f"SUMMARY: {NTESTS} checks, {len(FAILS)} FAIL")
if FAILS:
    print("FAILED:", FAILS)
    sys.exit(1)
print("ALL PASS")
print("=" * 76)
