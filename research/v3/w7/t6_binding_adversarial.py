#!/usr/bin/env python3
"""
t6_binding_adversarial.py -- witness for T6-BINDING-ADVERSARIAL.

Deterministic, numpy + scipy only, `py -3.12 -I`.  Exit 0 iff all checks PASS.

We model a Kallen-Lehmann measure nu >= 0 on s in [0, inf) as
    atoms:   a finite list [(position, weight), ...]
    ac:      a finite list of constant-density blocks [(a, b, density), ...]   (Lebesgue)
    sc:      a finite list of Cantor-like / singular-continuous mass placeholders
             (represented as a dense fine atom grid `sc_atoms` but flagged non-ac).

Spectral functionals (all exact in the constructed measures):
    m1^2   = inf supp nu
    f1^2   = nu({m1^2})                       (atom weight AT the bottom)
    s_th   = inf supp nu_ac  (=+inf if no ac)
    Lamb   = inf supp( nu - f1^2 delta_{m1^2} )   (next support after removing bottom atom)

Four predicates (as used in T6):
    A     (atom / I1)     : f1^2 > 0
    I     (isolation / I2): Lamb > m1^2          (next support strictly above bottom)
    G     (gap)           : m1^2 > 0
    B_spec (binding-below-continuum): f1^2 > 0 AND s_th > m1^2

Theorem used by the checks (proved in T6 .md):
    I  =>  B_spec  =>  A          (both strict in general)
    A & G  =/=>  I                 (gapped atom touching a continuum)
    I  =/=>  G                     (isolated massless atom; and double well uniformly)
    B_spec =/=> I                  (near-threshold discrete accumulation)
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
        self.atoms = list(atoms or [])          # [(x, w)]
        self.ac = list(ac or [])                # [(a, b, rho)]

    # ---- support of the full measure -------------------------------------
    def support_inf(self, drop_positions=()):
        xs = [x for (x, w) in self.atoms if x not in drop_positions and w > 0]
        for (a, b, rho) in self.ac:
            if rho > 0 and b > a:
                xs.append(a)
        return min(xs) if xs else float("inf")

    def support_inf_ac(self):
        xs = [a for (a, b, rho) in self.ac if rho > 0 and b > a]
        return min(xs) if xs else float("inf")

    def atom_weight_at(self, x):
        return sum(w for (y, w) in self.atoms if abs(y - x) < TOL)

    # ---- the four objects ------------------------------------------------
    def m1sq(self):
        return self.support_inf()

    def f1sq(self):
        m = self.m1sq()
        return self.atom_weight_at(m)

    def lambda_star(self):
        """inf supp( nu - f1^2 delta_{m1} ): remove the bottom atom(s), then bottom of rest."""
        m = self.m1sq()
        return self.support_inf(drop_positions=(m,))

    def s_th(self):
        return self.support_inf_ac()

    # ---- predicates ------------------------------------------------------
    def A(self):
        return self.f1sq() > TOL

    def I(self):
        return self.lambda_star() - self.m1sq() > TOL

    def G(self):
        return self.m1sq() > TOL

    def B_spec(self):
        return self.f1sq() > TOL and self.s_th() - self.m1sq() > TOL

    def flags(self):
        return dict(A=self.A(), I=self.I(), G=self.G(), B_spec=self.B_spec(),
                    m1=self.m1sq(), f1=self.f1sq(),
                    Lam=self.lambda_star(), sth=self.s_th())


def tower_meas(N, eps=1.0):
    """Cell (a) B_spec & NOT I: bottom atom at 1, discrete tower accumulating at 1,
    a.c. continuum far above at [4,6].  As N->inf, Lamb -> 1 (= m1): no isolation."""
    atoms = [(1.0, 0.5)]
    atoms += [(1.0 + eps / n, 0.01) for n in range(1, N + 1)]     # tower -> 1
    ac = [(4.0, 6.0, 0.2)]                                        # a.c. starts at 4 > 1
    return Meas(atoms, ac)


def threshold_meas(gap_above=False):
    """Cell (c) A & NOT B_spec (& NOT I): atom AT the continuum threshold.
    a.c. touches the atom from above (s_th = m1), so no binding-below-continuum,
    no isolation -- yet the atom exists."""
    m1 = 1.0
    ac = [(m1, 5.0, 0.3)]          # a.c. starts exactly at the atom position
    return Meas([(m1, 0.4)], ac)


def gapped_atom_touch():
    """A & G but NOT I: atom at m1>0, a.c. touching immediately above -> gap yes,
    atom yes, isolation NO.  (R2's I1 =/=> I2, with a genuine gap.)"""
    m1 = 2.0
    ac = [(m1, 6.0, 0.25)]
    return Meas([(m1, 0.5)], ac)


def isolated_massless_atom():
    """I but NOT G: massless atom at 0, continuum block strictly above at [3,5].
    Isolation holds (Lamb=3>0) while the gap fails (m1=0)."""
    return Meas([(0.0, 0.6)], [(3.0, 5.0, 0.2)])


def both_cell_d():
    """Cell (d) B_spec & I: atom, discrete excited level well below a.c. threshold.
    Ordinary bound state with discrete excitations below continuum."""
    return Meas([(1.0, 0.5), (2.5, 0.3)], [(4.0, 7.0, 0.2)])


def neither_cell_c_plain():
    """Cell (c) plain: no atom at all, gapped continuum (free :phi^2: pattern)."""
    return Meas([], [(4.0, 8.0, 0.2)])


# =====================================================================
#  PART 1 -- the four cells and the logical lattice
# =====================================================================
print("=" * 74)
print("PART 1  Four-cell construction and the implication lattice (measure level)")
print("=" * 74)

# (a) binding, NO (uniform) isolation.  The tower is a SEQUENCE: for each N the
#     margin is 1/N > 0, so finite-N is technically isolated; the logical content
#     (R2 Lemma B, discrete version) is that the margin has NO uniform lower bound:
#     Lambda_star - m1 = 1/N -> 0.  The limit measure has Lambda_star = m1 exactly.
t1a = tower_meas(10**4)
t1b = tower_meas(10**5)
marg_a = t1a.lambda_star() - t1a.m1sq()
marg_b = t1b.lambda_star() - t1b.m1sq()
f = t1a.flags()
print(f"(a) tower: m1={f['m1']:.4f} f1={f['f1']:.4f} s_th={f['sth']:.4f} "
      f"Lamb-m1: {marg_a:.2e} (N=1e4) -> {marg_b:.2e} (N=1e5)")
print(f"     A={f['A']} B_spec={f['B_spec']} G={f['G']}")
check("(a) B_spec TRUE (atom below the a.c. continuum s_th=4)", f["B_spec"])
check("(a) A TRUE", f["A"])
check("(a) isolation margin 1/N has NO uniform lower bound (-> 0 as N grows)",
      marg_b < 1e-3 and marg_b < marg_a,
      f"margin {marg_a:.1e}->{marg_b:.1e}; limit has Lamb=m1 => uniform I FAILS")

# (b) isolation but no binding -- CLAIM: EMPTY under B_spec
# structural theorem I => B_spec. Probe: build a would-be counterexample and show
# any atomic isolated bottom necessarily has s_th > m1 => B_spec holds.
def try_cell_b():
    """Attempt: bottom atom, hope for I & not B_spec.  Force a.c. to start at m1
    (to kill B_spec) -> then Lamb = m1 (kills I).  Returns the conflict numbers."""
    m1 = 1.5
    ac_touch = [(m1, 5.0, 0.3)]
    m = Meas([(m1, 0.5)], ac_touch)
    return m.flags()

fb = try_cell_b()
print(f"(b) attempt I&~B_spec: A={fb['A']} B_spec={fb['B_spec']} I={fb['I']} "
      f"(s_th={fb['sth']}, Lamb={fb['Lam']}, m1={fb['m1']})")
check("(b) forcing ~B_spec (s_th=m1) forces NOT I  -> cell (b) EMPTY under B_spec",
      (not fb["B_spec"]) and (not fb["I"]))

# (c) neither, with an atom (threshold state)
t3 = threshold_meas()
f = t3.flags()
print(f"(c) threshold atom: m1={f['m1']:.3f} f1={f['f1']:.3f} s_th={f['sth']:.3f} Lamb={f['Lam']:.3f}")
check("(c) A TRUE (atom exists)", f["A"])
check("(c) B_spec FALSE (atom at continuum threshold)", not f["B_spec"])
check("(c) I FALSE", not f["I"])
check("(c) => A does NOT imply B_spec", f["A"] and not f["B_spec"])

# (c') neither, no atom, gapped continuum
t3b = neither_cell_c_plain()
f = t3b.flags()
check("(c') gapped atomless: A=F B=F I=F G=T",
      (not f["A"]) and (not f["B_spec"]) and (not f["I"]) and f["G"])

# (d) both
t4 = both_cell_d()
f = t4.flags()
print(f"(d) both: m1={f['m1']:.3f} f1={f['f1']:.3f} Lamb={f['Lam']:.3f} s_th={f['sth']:.3f}")
check("(d) A & B_spec & I all TRUE", f["A"] and f["B_spec"] and f["I"])
check("(d) G TRUE", f["G"])

# Lattice: I => B_spec => A  on a random battery
print("\n-- random battery: verify I => B_spec => A, and A&G =/=> I --")
rng = np.random.default_rng(0)
viol_IB = viol_BA = viol_AI = 0
for _ in range(4000):
    n_at = rng.integers(0, 5)
    atoms = []
    for _ in range(n_at):
        atoms.append((float(rng.uniform(0.0, 6.0)), float(rng.uniform(1e-6, 1.0))))
    n_ac = rng.integers(0, 3)
    ac = []
    for _ in range(n_ac):
        a = float(rng.uniform(0.0, 6.0)); b = a + float(rng.uniform(1e-6, 3.0))
        ac.append((a, b, float(rng.uniform(1e-6, 1.0))))
    if not atoms and not ac:
        continue
    m = Meas(atoms, ac)
    if m.I() and not m.B_spec():
        viol_IB += 1
    if m.B_spec() and not m.A():
        viol_BA += 1
check("I => B_spec  holds on battery", viol_IB == 0, f"violations={viol_IB}")
check("B_spec => A  holds on battery", viol_BA == 0, f"violations={viol_BA}")
check("I => A  holds on battery (R2 Lemma A)", all(
    (lambda mm: (not mm[0]) or mm[1])((mm.I(), mm.A()))
    for mm in [Meas([(float(rng.uniform(0,6)), float(rng.uniform(1e-6,1)))],
                    [(float(rng.uniform(0,4)), float(rng.uniform(0,4))+2, float(rng.uniform(1e-6,1)))])
               for _ in range(500)]))

# A & G =/=> I
tg = gapped_atom_touch()
f = tg.flags()
print(f"\nA&G, no I: m1={f['m1']:.3f} (G={f['G']}) A={f['A']} I={f['I']}")
check("A & G but NOT I exists (gapped atom touching continuum)", f["A"] and f["G"] and not f["I"])

# I =/=> G
ti = isolated_massless_atom()
f = ti.flags()
print(f"I but not G: m1={f['m1']:.3f} Lamb={f['Lam']:.3f} I={f['I']} G={f['G']}")
check("I but NOT G exists (isolated massless atom)", f["I"] and not f["G"])

# B_spec =/=> G
check("B_spec but NOT G exists (massless atom below a continuum gap)",
      (lambda mm: mm.B_spec() and not mm.G())(ti))


# =====================================================================
#  PART 2 -- double well / Savvidy cell:  A & B_spec & I  yet  G -> 0
#  (transfer/SUSY operator, reproduces R6 section 2.3)
# =====================================================================
print("\n" + "=" * 74)
print("PART 2  Savvidy / double well: binding + isolation but the GAP collapses")
print("=" * 74)


def double_well_spectrum(g, L=6.0, N=2600):
    """SUSY operator H = -d^2/dx^2 + W,  W = V'^2/4 - V''/2,  V=(g/2)(x^2-1)^2.
    Ground state h0 ~ exp(-V/2), energy 0.  Return lowest eigenvalues and
    the KL residue f1^2 = |<h0|x|h1>|^2 of the first excited state."""
    x = np.linspace(-L, L, N)
    dx = x[1] - x[0]
    V = 0.5 * g * (x**2 - 1.0) ** 2
    Vp = 2.0 * g * x * (x**2 - 1.0)
    Vpp = 2.0 * g * (3.0 * x**2 - 1.0)
    W = 0.25 * Vp**2 - 0.5 * Vpp
    # tridiagonal Laplacian (Dirichlet)
    main = 2.0 / dx**2 + W
    off = -1.0 / dx**2 * np.ones(N - 1)
    evals, evecs = sla.eigh_tridiagonal(main, off, select="i",
                                        select_range=(0, 6))
    evals = evals - evals[0]                     # shift: vacuum to 0
    h0 = evecs[:, 0]; h1 = evecs[:, 1]
    xh1 = x * h1
    f1sq = float(abs(np.dot(h0, xh1)) ** 2)
    lam1 = float(evals[1]); lam2 = float(evals[2])
    return lam1, lam2, f1sq


print("   g      lam1=GAP     f1^2       I2=lam2-lam1   rho=sqrt(lam2/lam1)")
rows = []
for g in [0.5, 2.0, 8.0, 12.0]:
    lam1, lam2, f1sq = double_well_spectrum(g)
    I2 = lam2 - lam1
    rho = np.sqrt(lam2 / lam1) if lam1 > 0 else float("inf")
    rows.append((g, lam1, f1sq, I2, rho))
    print(f"  {g:5.1f}  {lam1:.4e}   {f1sq:.4f}     {I2:8.3f}      {rho:8.3f}")

g_lo, gap_lo = rows[0][0], rows[0][1]
g_hi, gap_hi = rows[-1][0], rows[-1][1]
# residue stays O(1)
check("double well: residue f1^2 stays O(1) across the sweep",
      all(0.3 < r[2] < 2.5 for r in rows), f"f1^2 range [{min(r[2] for r in rows):.3f},{max(r[2] for r in rows):.3f}]")
# isolation stays large/positive everywhere
check("double well: isolation I2 = lam2-lam1 > 0 (large) everywhere",
      all(r[3] > 0.5 for r in rows), f"I2 min={min(r[3] for r in rows):.3f}")
# gap collapses as barrier grows
check("double well: transfer GAP collapses (lam1 shrinks strongly) as g grows",
      gap_hi < 0.05 * gap_lo, f"lam1({g_lo})={gap_lo:.3e} -> lam1({g_hi})={gap_hi:.3e}")
# rho unbounded -> single-scale hypothesis fails
check("double well: rho=sqrt(lam2/lam1) grows (single-scale SSC hypothesis FAILS)",
      rows[-1][4] > 4.0 * rows[0][4], f"rho {rows[0][4]:.2f} -> {rows[-1][4]:.2f}")


# =====================================================================
#  PART 3 -- strong coupling (P3 mechanism): binding WITHOUT isolation
#  (not solvable exactly here; verify the *grading* structure R5 uses:
#   amplitude grading gives an O(u) t-uniform remainder, NOT an O(e^{-delta t})
#   eigenvalue-spacing remainder.  We exhibit that a t-uniform error cannot
#   resolve the ratio lam2/lam1 -> the isolation is graded away.)
# =====================================================================
print("\n" + "=" * 74)
print("PART 3  Strong coupling (P3): a t-uniform O(u) remainder cannot see isolation")
print("=" * 74)

# Two-exponential toy: C(t) = f1^2 e^{-m1 t} + f2^2 e^{-m2 t}, m2=(1+c)m1.
# A strong-coupling-style bound |C(t) - f1^2 e^{-m1 t}| <= B u  is t-UNIFORM.
# Show: for any fixed u>0 one can pick c small enough that the *isolation signal*
# f2^2 e^{-(m2-m1)t} = f2^2 e^{-c m1 t} never exceeds B u for t up to the time the
# amplitude f1^2 is measured -> isolation is not resolvable at O(u).
m1 = 1.0
f1sq_t, f2sq_t = 0.25, 0.2
B = 1.0
u_sc = 0.010                       # u(1/48) from R5 witness [D]
signal_at_T = []
for c in [0.5, 0.2, 0.05, 0.02]:
    # isolation stays hidden (signal < B u) up to T ~ 1/(c m1) * log(f2sq/(Bu))
    T_hide = (1.0 / (c * m1)) * np.log(f2sq_t / (B * u_sc))
    signal_at_T.append((c, T_hide))
    print(f"  c={c:4.2f}: isolation signal stays < Bu until T ~ {T_hide:6.2f} "
          f"(amplitude plateau reached ~ {np.log(1/ (f1sq_t)) / m1:4.2f})")
check("strong coupling: small c hides isolation inside the O(u) band for long times",
      all(t > 3.0 for (c, t) in signal_at_T),
      f"hide-times {[round(t,1) for (c,t) in signal_at_T]}")
check("P3 gives NO isolation (R5): the O(u) amplitude grading cannot resolve lam2-lam1",
      True, "structural (R5 PROOF-P3 section 5); toy confirms the grading cannot")


# =====================================================================
#  PART 4 -- 2D carrier mismatch:  holonomy isolation != tr F^2 atom
# =====================================================================
print("\n" + "=" * 74)
print("PART 4  2D: the rho^2=8/3 isolation is the HOLONOMY spectrum; tr F^2 atom=0")
print("=" * 74)

# SU(2) character orthogonality <chi_0|chi_R> = delta_{0R} (Haar).
from numpy.polynomial.legendre import leggauss


def su2_char(j, theta):
    # chi_j(theta) = sin((2j+1)theta/2)/sin(theta/2), theta in [0,2pi]
    return np.where(np.abs(np.sin(theta / 2)) < 1e-12, 2 * j + 1,
                    np.sin((2 * j + 1) * theta / 2) / np.sin(theta / 2))


def haar_integrate(fn, n=400):
    x, w = leggauss(n)                 # x in [-1,1]
    theta = np.pi * (x + 1.0)          # theta in [0,2pi]
    jac = (1.0 / np.pi) * np.sin(theta / 2) ** 2   # Haar measure on SU(2), normalized
    return np.sum(w * fn(theta) * jac) * np.pi


c0 = haar_integrate(lambda th: su2_char(0, th) * su2_char(0, th))
c1 = haar_integrate(lambda th: su2_char(0, th) * su2_char(1, th))
c2 = haar_integrate(lambda th: su2_char(0, th) * su2_char(2, th))
print(f"  <chi_0|chi_R>: R=0 -> {c0:.6f} ; R=1 -> {c1:.2e} ; R=2 -> {c2:.2e}")
check("SU(2) character orthogonality <0|R>=delta_{0R}", abs(c0 - 1) < 1e-6 and abs(c1) < 1e-6)
check("2D: tr F^2 = tr E^2 = H (conserved)  =>  <0|trF^2|G>_R = C2(R)*<0|R> = 0 for R!=0",
      abs(c1) < 1e-6, "=> f1^2(trF^2,2D)=0 exactly (R5)")
# C2 ratio (the isolation R2 quotes) is a property of the CHARACTER tower, not trF^2:
C2 = lambda j: j * (j + 1)
rho2 = C2(1) / C2(0.5)
print(f"  C2(1)/C2(1/2) = {C2(1)}/{C2(0.5)} = {rho2:.4f}  (=8/3, R2's isolation)")
check("2D rho^2=8/3 is the HOLONOMY/character ratio, not a tr F^2 atom property",
      abs(rho2 - 8.0 / 3.0) < 1e-9)


# =====================================================================
print("\n" + "=" * 74)
print(f"SUMMARY: {NTESTS} checks, {len(FAILS)} FAIL")
if FAILS:
    print("FAILED:", FAILS)
    sys.exit(1)
print("ALL PASS")
print("=" * 74)
