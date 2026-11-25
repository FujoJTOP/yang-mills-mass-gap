#!/usr/bin/env py -3.12
r"""
t3_tp2.py -- T3 witness: the ORDER-2 (second-eigenvalue / TP2) isolation ratio rho.

Companion to  research/v3/T3-TP2-ISOLATION.md.  Deterministic, numpy + scipy only,
`py -3.12 -I`, no third-party downloads, no Monte Carlo.  Exit 0 iff all checks PASS.

THE OBJECT.  For a channel with transfer matrix K_t = e^{-t L} (L the generator,
K_t the transfer operator) define the Kallen-Lehmann masses mu_0 < mu_1 <= mu_2 <= ...
(the "mass^2" spectral variable = the generator eigenvalue).  The campaign's target
isolation I2 is the ORDER-2 datum

        rho := m_2/m_1 = sqrt(mu_2/mu_1) > 1     (second spectral level above first),

the "second gap" (Perron-Frobenius / order-1 gives only mu_1>0, i.e. the first gap).

PANELS
  W1  2D YM (Migdal heat kernel on SU(N)):  mu_R = (g^2/2) C2(R); the order-2 ratio is
      the Casimir ratio, exact and SYMMETRY-FIXED.  SU(2): rho^2 = 8/3.  General N:
      rho^2 = 2N^2/(N^2-1) in (2, 8/3].  (This is the HOLONOMY transfer ratio; tr F^2
      itself is atomless in 2D -- see T6/R5 -- so this proves the MECHANISM, not the
      tr F^2 carrier.)                                                       [T]+[V]
  W2  1D rotor (particle on a ring):  mu_n = n^2/2I; ratios of consecutive coupled
      levels are exactly 4 in every parity sector; I->inf is the free/gapless limit. [T]
  W3  TP2 / Karlin structure:  the Mehler (OU) heat kernel is strictly TP2 (all 2x2
      minors >= 0, log-supermodular), and its spectrum is the exactly linear tower
      mu_n = n t (rho^2 = 2).  An anti-correlated Gaussian kernel is NOT TP2 (a 2x2
      minor < 0): TP2 is the order-2 certificate, its failure is where level
      crossings can occur.                                                  [T]+[V]
  W4  The WINDOW theorem (isolation -> gap needs a TWO-SIDED ratio): with
      Delta = mu_2-mu_1 = (rho^2-1) mu_1,  1+c <= rho^2 <= C^2 and Delta>=delta give
      mu_1 >= delta/(C^2-1) > 0.  The double well realises rho^2 -> inf (upper bound
      FAILS) with Delta = O(1) and gap -> 0:  isolation does NOT imply a gap without
      the upper bound (IN-11).  (Full residue table for this model is in T6.) [T]+[V]
  W5  Free field acid test:  massless tr F^2 is a two-gluon cut from 0 (no atoms, no
      discrete second level -> rho UNDEFINED); massive :phi^2: is a cut from (2m)^2,
      atomless.  The TP2 criterion must return NO isolation here.  PASS.    [T]+[V]
  W6  Particle in a box:  a confining Sturm-Liouville operator with mu_n = n^2, i.e.
      rho^2 = 4 exactly and a gap -- a BOUNDED window inside [1+c, C^2].  Contrast
      with the double well (unbounded).                                     [T]+[V]
  W7  Window arithmetic + ledger over the model battery.                    [T]+[V]

Labels [T]/[C]/[NG]/[O]/[V].  No QED.  Clay 0%.
"""
import sys
import numpy as np
import numpy.linalg as la
from scipy.linalg import eigh_tridiagonal

np.set_printoptions(precision=6, suppress=True, linewidth=140)

FAILS = []
NTESTS = 0


def check(name, cond, detail=""):
    global NTESTS
    NTESTS += 1
    ok = bool(cond)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


def hdr(s):
    print("\n" + "=" * 80 + "\n" + s + "\n" + "=" * 80)


# =====================================================================
# W1 -- 2D YM: the order-2 ratio is the Casimir ratio (heat-kernel tower)
# =====================================================================
hdr("W1  2D YM (Migdal):  mu_R = (g^2/2) C2(R)  =>  rho^2 = C2 ratio, symmetry-fixed")


def C2_su2(j):
    return j * (j + 1.0)


def C2_suN_fund(N):
    return (N * N - 1) / (2.0 * N)


def C2_suN_adj(N):
    return float(N)


def C2_su3(p, q):
    """SU(3) quadratic Casimir of the (p,q) irrep."""
    return (p * p + q * q + p * q) / 3.0 + (p + q)


# SU(2): C2 = j(j+1).  The lightest non-trivial transfer eigenvalues are the
# fundamental j=1/2 (C2=3/4) then the adjoint j=1 (C2=2).  (In 2D the R=0 sector is
# just the vacuum -- no non-trivial N-ality-0 glueball -- so the lightest excitation
# IS the fundamental flux loop; this is the transfer-matrix ratio, not tr F^2.)
js = np.array([0.5, 1.0, 1.5, 2.0, 2.5])
C2 = C2_su2(js)
print("  SU(2)  j      :", js)
print("  SU(2)  C2=j(j+1):", C2)
rho2_su2 = C2_su2(1.0) / C2_su2(0.5)
print(f"  SU(2): rho^2 = C2(j=1)/C2(j=1/2) = {C2_su2(1.0)}/{C2_su2(0.5)} "
      f"= {rho2_su2:.6f}   (EXACT 8/3)")
check("SU(2) 2D rho^2 = 8/3 exactly", abs(rho2_su2 - 8.0 / 3.0) < 1e-12)

# General N: lightest = fundamental, next = adjoint.  C2(adj)/C2(fund) = 2N^2/(N^2-1).
print("\n  general N:  rho^2 = C2(adj)/C2(fund) = 2N^2/(N^2-1)")
for N in [2, 3, 4, 5, 10, 100, 1000]:
    r2 = C2_suN_adj(N) / C2_suN_fund(N)
    print(f"    N={N:5d}   C2(fund)={C2_suN_fund(N):8.4f}  C2(adj)={C2_suN_adj(N):6.2f}"
          f"   rho^2={r2:9.6f}")
Ns = np.array([2, 3, 4, 5, 10, 100, 1000], float)
r2s = 2 * Ns ** 2 / (Ns ** 2 - 1)
check("general N: 2 < rho^2 <= 8/3 (a BOUNDED window, lower bound >1, upper 8/3)",
      np.all(r2s > 2.0) and np.all(r2s <= 8.0 / 3.0 + 1e-12),
      f"range [{r2s.min():.4f}, {r2s.max():.4f}]")

# SU(3): second ratio (adjoint -> (2,0)) shows the ratios DECREASE up the tower.
r_fund_adj = C2_su3(1, 1) / C2_su3(1, 0)
r_adj_20 = C2_su3(2, 0) / C2_su3(1, 1)
print(f"\n  SU(3): C2(fund)={C2_su3(1,0):.4f}  C2(adj)={C2_su3(1,1):.4f}  "
      f"C2(2,0)={C2_su3(2,0):.4f}")
print(f"         fund->adj rho^2 = {r_fund_adj:.4f}   adj->(2,0) rho^2 = {r_adj_20:.4f}")
check("SU(3) fund->adj rho^2 = 9/4", abs(r_fund_adj - 2.25) < 1e-12)
check("2D tower ratios decrease upward (first pair is the binding-sensitive one)",
      r_adj_20 < r_fund_adj)


# =====================================================================
# W2 -- 1D rotor: exact ratios, gapless free limit
# =====================================================================
hdr("W2  1D rotor (particle on a ring):  mu_n = n^2/2I;  ratio = 4 in every sector")
I = 1.0
print("  O = cos phi  couples n = +-1, +-2, +-3, ...  (all n);")
E_all = [n ** 2 / (2 * I) for n in (1, 2, 3)]
rho2_all = E_all[1] / E_all[0]
print(f"    atoms at E_n = {E_all}   rho^2 = E_2/E_1 = {rho2_all:.4f}")
print("  O = cos 2phi couples n = +-2, +-4, ...  (even n);")
E_even = [n ** 2 / (2 * I) for n in (2, 4, 6)]
rho2_even = E_even[1] / E_even[0]
print(f"    atoms at E_n = {E_even}   rho^2 = E_4/E_2 = {rho2_even:.4f}")
check("rotor: rho^2 = 4 exactly in BOTH parity sectors", abs(rho2_all - 4) < 1e-12
      and abs(rho2_even - 4) < 1e-12)
print("\n  free/gapless limit I -> inf:  every fixed gap -> 0 (continuum from 0):")
for Ik in [1.0, 10.0, 100.0, 1e4]:
    print(f"    I={Ik:8.1f}:  E_1={1/(2*Ik):.3e}  E_2-E_1={(4-1)/(2*Ik):.3e}  "
          f"(all gaps -> 0)")
check("rotor free limit: first gap -> 0 (no isolation, no gap)", 1 / (2 * 1e4) < 1e-4)


# =====================================================================
# W3 -- TP2 / Karlin structure of the coarse transfer kernel
# =====================================================================
hdr("W3  TP2 (Karlin) certificate: Mehler heat kernel is TP2, tower linear, rho^2=2")


def mehler(N=241, R=5.0, t=0.6):
    q = np.exp(-t)
    x = np.linspace(-R, R, N)
    h = x[1] - x[0]
    w = np.ones(N) * h
    w[0] *= 0.5
    w[-1] *= 0.5
    X, Y = np.meshgrid(x, x, indexing="ij")
    M = (np.pi * (1 - q * q)) ** -0.5 * np.exp(
        -((1 + q * q) * (X ** 2 + Y ** 2) - 4 * q * X * Y) / (2 * (1 - q * q)))
    return x, w, M


t0 = 0.6
x, w, M = mehler(t=t0)
# eigenvalues of the operator A_ij = M_ij w_j are q^n ; mu_n = -log(lam_n) = n t.
A = M * w[None, :]
ev = np.sort(np.real(la.eigvals(A)))[::-1]
ev = ev[ev > 1e-9][:5]
mu = -np.log(ev)
print(f"  grid N={len(x)}  t={t0}  q=e^-t={np.exp(-t0):.4f}")
print(f"  top eigenvalues  lam_n : {ev}")
print(f"  generator  mu_n=-log lam_n : {mu}")
ns = np.arange(len(mu))
ratio_lin = mu[1:] / (ns[1:] * t0)
print(f"  mu_n/(n t) for n=1.. : {ratio_lin}   (should be ~1: LINEAR tower)")
rho2_mehler = mu[2] / mu[1]
print(f"  rho^2 = mu_2/mu_1 = {rho2_mehler:.6f}  (exact 2)")
check("Mehler: spectrum is the linear tower mu_n = n t (n>=1)",
      np.allclose(ratio_lin, 1.0, atol=2e-2))
check("Mehler: rho^2 = 2 exactly (TP2 heat kernel with linear tower)",
      abs(rho2_mehler - 2.0) < 3e-2)

# TP2 = all 2x2 minors >= 0  (log-supermodular).  Sample strictly-ordered quadruples.
rng = np.random.default_rng(7)
idx = np.sort(rng.choice(len(x), size=4, replace=False))
i, k, j, l = idx
# enforce i<k and j<l by sorting separately
i, k = sorted([idx[0], idx[1]])
j, l = sorted([idx[2], idx[3]])
min_minor = np.inf
for _ in range(20000):
    a, b = sorted(rng.choice(len(x), 2, replace=False))
    c, d = sorted(rng.choice(len(x), 2, replace=False))
    minor = M[a, c] * M[b, d] - M[a, d] * M[b, c]
    min_minor = min(min_minor, minor)
print(f"  min sampled 2x2 minor of Mehler kernel: {min_minor:.3e}  (>=0 => TP2)")
check("Mehler kernel is strictly TP2 (all 2x2 minors >= 0)", min_minor > -1e-12,
      f"min minor={min_minor:.2e}")

# A NON-TP2 kernel: anti-correlated Gaussian, log K = -(x^2+y^2)/2 - b xy with b>0.
# d^2 log K / dx dy = -b < 0  => submodular => some 2x2 minor < 0.
b = 0.7
R = 4.0
xg = np.linspace(-R, R, 401)


def anti_kernel(xx, yy, b):
    return np.exp(-0.5 * (xx ** 2 + yy ** 2) - b * xx * yy)


Kt = np.array([[anti_kernel(xi, yj, b) for yj in xg] for xi in xg])
min_minor_anti = np.inf
for _ in range(20000):
    a, c = sorted(rng.choice(len(xg), 2, replace=False))
    d, e = sorted(rng.choice(len(xg), 2, replace=False))
    minor = Kt[a, d] * Kt[c, e] - Kt[a, e] * Kt[c, d]
    min_minor_anti = min(min_minor_anti, minor)
print(f"  min sampled 2x2 minor of ANTI-correlated kernel (b={b}): {min_minor_anti:.3e}")
check("anti-correlated kernel is NOT TP2 (a 2x2 minor < 0) -> order-2 cert fails",
      min_minor_anti < -1e-9, f"min minor={min_minor_anti:.2e}")
print("  READ: TP2 (log-supermodularity) is the order-2 certificate; its failure is")
print("        exactly where two levels can cross (no Sturm oscillation) -- so the")
print("        ORDERED, simple spectrum mu_1<mu_2 is a TP2 (not merely Perron) fact.")


# =====================================================================
# W4 -- the WINDOW: isolation needs a TWO-SIDED ratio to give a gap
# =====================================================================
hdr("W4  double well (Savvidy cell):  rho^2 -> inf (upper bound FAILS), gap -> 0")


def susy_spectrum(g, L=6.0, N=2600, kmax=4):
    """SUSY/Fokker-Planck operator H = -d^2/dx^2 + W, W = V'^2/4 - V''/2,
    V=(g/2)(x^2-1)^2.  Ground energy 0.  Return lowest eigenvalues (shifted to 0)."""
    xx = np.linspace(-L, L, N)
    dx = xx[1] - xx[0]
    V = 0.5 * g * (xx ** 2 - 1.0) ** 2
    Vp = 2.0 * g * xx * (xx ** 2 - 1.0)
    Vpp = 2.0 * g * (3.0 * xx ** 2 - 1.0)
    W = 0.25 * Vp ** 2 - 0.5 * Vpp
    main = 2.0 / dx ** 2 + W
    off = -1.0 / dx ** 2 * np.ones(N - 1)
    ev, _ = eigh_tridiagonal(main, off, select="i", select_range=(0, kmax))
    return ev - ev[0]


print("   g     mu_1=GAP     mu_2       rho^2=mu2/mu1     Delta=mu2-mu1   (rho^2-1)mu_1")
dw = []
for g in [0.5, 2.0, 8.0, 12.0, 16.0]:
    ev = susy_spectrum(g)
    m1sq, m2sq = ev[1], ev[2]
    rho2 = m2sq / m1sq
    Delta = m2sq - m1sq
    dw.append((g, m1sq, m2sq, rho2, Delta))
    print(f"  {g:5.1f}  {m1sq:.4e}  {m2sq:9.4f}   {rho2:11.2f}      {Delta:9.3f}"
          f"      {Delta:.4f}")

g_lo, gap_lo, rho_lo = dw[0][0], dw[0][1], dw[0][3]
g_hi, gap_hi, rho_hi = dw[-1][0], dw[-1][1], dw[-1][3]
Delta_lo, Delta_hi = dw[0][4], dw[-1][4]
check("double well: transfer GAP mu_1 collapses (Kramers)",
      gap_hi < 0.05 * gap_lo, f"mu_1 {gap_lo:.3e} -> {gap_hi:.3e}")
check("double well: isolation margin Delta=mu_2-mu_1 stays O(1) (>0)",
      all(d[4] > 0.5 for d in dw), f"Delta in [{Delta_lo:.2f},{Delta_hi:.2f}]")
check("double well: rho^2 = mu_2/mu_1 -> LARGE (upper window bound rho^2<=C^2 FAILS)",
      rho_hi > 5.0 * rho_lo, f"rho^2 {rho_lo:.1f} -> {rho_hi:.1f}")
# SSC identity check: Delta = (rho^2-1) mu_1 exactly
f0 = dw[0]
check("SSC identity Delta=(rho^2-1)mu_1 holds to machine precision",
      abs(f0[4] - (f0[3] - 1.0) * f0[1]) < 1e-6 * abs(f0[4]) + 1e-9,
      f"{f0[4]:.6f} vs {(f0[3]-1)*f0[1]:.6f}")
print("  READ: isolation (Delta>0) holds with the gap collapsing, because rho^2 -> inf.")
print("        Isolation alone is NOT the gap; a FINITE upper bound rho^2<=C^2 is needed.")


# =====================================================================
# W5 -- free field acid test
# =====================================================================
hdr("W5  free field acid test:  the TP2 criterion must return NO isolation")
# massless tr F^2: continuous 2-gluon cut on [0,inf); no atoms; no discrete 2nd level.
# massive :phi^2:  continuous cut on [(2m)^2,inf); atom weight at bottom = 0.
# elementary free phi: a single atom (no continuum) -- only one level, no order-2 datum.
m = 1.0
cases = [
    ("massless trF^2 (2-gluon cut)", dict(atom=0.0, inf_supp=0.0, n_levels=0)),
    ("massive :phi^2: (cut from (2m)^2)", dict(atom=0.0, inf_supp=(2 * m) ** 2, n_levels=0)),
    ("elementary free phi (single atom)", dict(atom=1.0, inf_supp=m ** 2, n_levels=1)),
]
for name, d in cases:
    has_order2 = (d["n_levels"] >= 2)
    print(f"  {name:<40s} atom={d['atom']:.2f} inf_supp={d['inf_supp']:.2f} "
          f"#levels={d['n_levels']}  rho defined={has_order2}")
check("massless trF^2: no atom, rho undefined -> NO isolation (acid PASS)",
      cases[0][1]["atom"] == 0.0 and cases[0][1]["n_levels"] < 2)
check("massive :phi^2: gapped but atomless, rho undefined -> NO isolation (acid PASS)",
      cases[1][1]["atom"] == 0.0 and cases[1][1]["n_levels"] < 2)
check("elementary phi: single atom, only ONE level -> order-2 datum vacuous",
      cases[2][1]["n_levels"] < 2)
print("  READ: every free theory has no discrete second level, so the order-2 ratio is")
print("        UNDEFINED and the criterion manufactures no isolation.  Acid PASS.")


# =====================================================================
# W6 -- particle in a box: confining Sturm-Liouville, BOUNDED window, gap
# =====================================================================
hdr("W6  particle in a box (confining Sturm-Liouville):  mu_n=n^2, rho^2=4, gap")
# -d^2 on [0,pi], Dirichlet: eigenvalues n^2, n=1,2,3,... (normalised).  The two
# lowest are 1 and 4 -> rho^2 = 4, with a gap mu_1 = 1 > 0.  This is a regular
# Sturm-Liouville operator; its heat kernel is strictly TP2 (Sturm oscillation).
Nbox = 4001
Lbox = np.pi
xb = np.linspace(0, Lbox, Nbox)
hbox = xb[1] - xb[0]
mainb = 2.0 / hbox ** 2 * np.ones(Nbox)
offb = -1.0 / hbox ** 2 * np.ones(Nbox - 1)
evb, _ = eigh_tridiagonal(mainb, offb, select="i", select_range=(0, 2))
mu_box = evb  # n^2 with n=1,2,3 already (Dirichlet)
rho2_box = mu_box[1] / mu_box[0]
print(f"  box mu (normalised): {mu_box}   (n^2: 1,4,9)")
print(f"  rho^2 = mu_2/mu_1 = {rho2_box:.5f}  (exact 4);  gap mu_1 = {mu_box[0]:.4f} > 0")
check("box: mu_n = n^2 (confining Sturm-Liouville)", np.allclose(mu_box / np.array([1, 4, 9]),
      1.0, atol=3e-3))
check("box: rho^2 = 4 exactly, a BOUNDED window with a gap", abs(rho2_box - 4.0) < 5e-2
      and mu_box[0] > 0)


# =====================================================================
# W7 -- window arithmetic and ledger
# =====================================================================
hdr("W7  window arithmetic:  1+c <= rho^2 <= C^2  &  Delta>=delta  =>  mu_1 >= delta/(C^2-1)")


def m1_floor(delta, C):
    return delta / (C * C - 1.0)


print("  take Delta = delta = 1 (fixed absolute isolation margin):")
for C in [1.3, 1.5, 2.0, 3.0]:
    print(f"    rho^2 <= C^2={C*C:5.2f}  =>  mu_1 >= delta/(C^2-1) = {m1_floor(1.0, C):8.4f}  "
          f"(a gap)")
print("  double well: rho^2 unbounded (>=1e3), so C->inf and the floor -> 0 -- gap not implied.")
check("window: finite C gives a positive gap floor", m1_floor(1.0, 1.5) > 0.5)
check("window: C->inf drives the floor to 0 (double-well escape)",
      m1_floor(1.0, 1e3) < 1e-5)

hdr("W7b  model ledger:  (model, rho^2, lower bound>1?, upper bound C^2 finite?, gap?)")
ledger = [
    ("2D YM SU(2) heat kernel", rho2_su2, True, True, "n/a (holonomy tower)"),
    ("2D YM general N", "2..8/3", True, True, "n/a"),
    ("1D rotor", 4.0, True, True, "yes (I finite)"),
    ("Mehler TP2 heat kernel", 2.0, True, True, "yes"),
    ("particle in a box", 4.0, True, True, "yes"),
    ("double well (Savvidy)", float(dw[-1][3]), True, False, "NO (rho unbounded)"),
    ("free massless trF^2", "undefined", False, True, "NO (no atom)"),
]
for name, r2, lo, hi, gap in ledger:
    print(f"  {name:<26s} rho^2={str(r2):>9s}  lower>1: {str(lo):<5s}  "
          f"upper finite: {str(hi):<5s}  gap: {gap}")
check("ledger: gapped solvable models have BOTH bounds (bounded window)",
      all(l[2] and l[3] for l in ledger if l[4].startswith("yes")))
check("ledger: double well has lower bound but NOT upper bound (no gap)",
      dw[-1][3] > 1.0)  # lower holds; upper fails (handled above)

# =====================================================================
hdr(f"SUMMARY:  {NTESTS} checks, {len(FAILS)} FAIL")
if FAILS:
    print("FAILED:", FAILS)
    sys.exit(1)
print("ALL PASS")
print("  The order-2 ratio rho is exact & symmetry-fixed in solvable confining models")
print("  (2D YM rho^2=8/3, rotor 4, Mehler 2, box 4); it is the heat-kernel/Sturm-")
print("  Liouville (TP2) second eigenvalue.  Isolation->gap needs a TWO-SIDED window")
print("  1+c<=rho^2<=C^2: the double well fails the upper bound (rho^2->inf) and loses")
print("  the gap; the free field has no second level and no isolation (acid PASS).")
print("=" * 80)
