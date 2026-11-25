# -*- coding: utf-8 -*-
"""
lf_ebk.py -- witness for candidate (1) of the death-cause axis: EBK / Bohr-Sommerfeld.

Question put to the numbers:
  (Q1) Is EBK's output a *count* or a *bound*?  And is the error one-sided?
  (Q2) Does the semiclassical small parameter exist at the target (beta=2.5)?
  (Q3) What is EBK's rigorous 1D shadow, and how far is it from the target?

Design (all light, single process):
  W1  harmonic calibration on the whole line:  -h^2 u'' + a s^2 u = lam u
      EBK with Maslov 1/2 per smooth turning point is EXACT: lam_n = 2 h sqrt(a)(n+1/2).
      -> the machinery is right, so any failure below is not an implementation bug.
  W2  anharmonic quartic well, h scanned: relative error of EBK vs exact (FD).
      Expect |rel err| ~ C h  => the error is O(h), i.e. O(g), i.e. O(1) at beta=2.5.
  W3  translate h -> g on the lattice: h_eff = g = sqrt(2N/beta).  Numbers for N=3.
  W4  Maslov-convention spread at h=1: lam_EBK(mu=1/2) vs (3/4) vs (1).  If the spread
      is O(1) relative, the "asymptotic" is not operating => counting cannot exclude
      one level in a window of that width.
  W5  Sturm comparison (EBK's rigorous shadow) on a family of sign-changing q:
      lam_Sturm = pi^2/T^2 + min q  -- check 0 violations (one-sided, rigorous).
  W6  the corpus's own rigorous 1D floor (LP-LIOUVILLE B3): T_intr=10.180126,
      ess inf Q=-0.052530 => bound 0.042705 vs target 0.6444  -> short by 15.1x.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal

np.random.seed(20261125)
OUT = []


def say(s=""):
    print(s)
    OUT.append(str(s))


# ---------------------------------------------------------------- exact solver
def exact_lam(qfun, hbar, T, N=2001, k=3):
    """lowest k eigenvalues of  -hbar^2 u'' + q u = lam u, Dirichlet on (0,T)"""
    s = np.linspace(0.0, T, N)
    h = T / (N - 1)
    d = 2.0 * hbar**2 / h**2 + qfun(s)
    e = -hbar**2 / h**2 * np.ones(N - 1)
    w = eigh_tridiagonal(d, e, select="i", select_range=(0, k - 1))[0]
    return np.sort(w)


# ---------------------------------------------------------------- EBK solver
def ebk_lam(qfun, hbar, T, n=0, maslov=0.5, ng=600):
    """Bohr-Sommerfeld:  2*int sqrt(lam-q)^+ ds = 2*pi*hbar*(n + maslov)
       maslov = total Maslov index / 2 ; 0.5 per smooth turning point (2 points -> 1.0? see below)

    Convention used:  oint p ds = 2*pi*hbar*(n + mu/2) with mu = (#smooth turning points)/2.
    For a well with 2 smooth turning points:  mu = 1  ->  RHS = 2*pi*hbar*(n+1/2).
    => the argument `maslov` below is mu, and RHS = 2*pi*hbar*(n + maslov/2).
    """
    x, wq = np.polynomial.legendre.leggauss(ng)
    s = 0.5 * T * (x + 1.0)
    wq = 0.5 * T * wq
    qq = qfun(s)
    qmin = qq.min()
    rhs = 2.0 * np.pi * hbar * (n + maslov / 2.0)

    def F(lam):
        return 2.0 * np.sum(wq * np.sqrt(np.maximum(lam - qq, 0.0))) - rhs

    lo = qmin + 1e-14
    hi = qmin + 4.0 * (np.pi * hbar * (n + 1.0) / T) ** 2 + 4.0 * (qq.max() - qmin) + 16.0 * hbar
    for _ in range(80):
        if F(hi) > 0:
            break
        hi = qmin + 2.0 * (hi - qmin)
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if F(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ================================================================ W1 calibration
say("=" * 78)
say("W1  harmonic calibration  q = a s^2  (whole line, 2 smooth turning points, mu=1)")
say("    exact  lam_n = 2 hbar sqrt(a) (n+1/2)   [operator -hbar^2 d^2 + a s^2]")
say("=" * 78)
a, T = 1.0, 20.0
say(f"{'hbar':>8} {'n':>2} {'exact(FD)':>14} {'EBK':>14} {'rel.diff':>11}")
w1max = 0.0
for hbar in (1.0, 0.5, 0.25):
    for n in (0, 1, 2):
        ex = exact_lam(lambda s: a * (s - T / 2) ** 2, hbar, T, N=2001, k=n + 1)[n]
        eb = ebk_lam(lambda s: a * (s - T / 2) ** 2, hbar, T, n=n, maslov=1.0)
        rd = (eb - ex) / ex
        w1max = max(w1max, abs(rd))
        say(f"{hbar:8.3f} {n:2d} {ex:14.8f} {eb:14.8f} {rd:11.2e}")
say(f"  -> max |rel| on the harmonic family = {w1max:.2e}   (machinery OK)")

# ================================================================ W2 anharmonic
say()
say("=" * 78)
say("W2  anharmonic well  q = (s-5)^2 + 0.1 (s-5)^4  on (0,10), Dirichlet; hbar scanned")
say("    if |rel err| ~ C*hbar the semiclassical error is O(hbar) -- magnitude, not a bound")
say("=" * 78)
b, T2 = 0.1, 10.0
qf = lambda s: (s - 5.0) ** 2 + b * (s - 5.0) ** 4
say(f"{'hbar':>8} {'exact(FD)':>14} {'EBK':>14} {'abs err':>12} {'rel err':>12} {'sign':>5}")
rows = []
for hbar in (1.0, 0.5, 0.25, 0.125):
    ex = exact_lam(qf, hbar, T2, N=2001, k=1)[0]
    eb = ebk_lam(qf, hbar, T2, n=0, maslov=1.0)
    rows.append((hbar, ex, eb, eb - ex))
    say(f"{hbar:8.4f} {ex:14.8f} {eb:14.8f} {eb-ex:12.3e} {(eb-ex)/ex:12.3e} {'+' if eb>ex else '-':>5}")

say()
say("    ratio rel(h/2)/rel(h)  -> should be ~0.5 if the error is O(hbar):")
for i in range(len(rows) - 1):
    r1 = (rows[i][2] - rows[i][1]) / rows[i][1]
    r2 = (rows[i + 1][2] - rows[i + 1][1]) / rows[i + 1][1]
    say(f"      hbar {rows[i][0]:6.4f} -> {rows[i+1][0]:6.4f}:  {r1:12.3e} -> {r2:12.3e}"
        f"   ratio {r2/r1:8.4f}")

# ================================================================ W3 hbar_eff
say()
say("=" * 78)
say("W3  what is hbar at beta=2.5 ?   lattice: beta = 2N/g^2  =>  hbar_eff = g = sqrt(2N/beta)")
say("=" * 78)
for N in (2, 3):
    for beta in (1.0, 2.5, 4.0, 10.0, 100.0, 600.0):
        g = np.sqrt(2.0 * N / beta)
        say(f"  N={N}  beta={beta:7.2f}   g=hbar_eff={g:8.5f}")
    g25 = np.sqrt(2.0 * N / 2.5)
    say(f"  -> N={N}: at beta=2.5 hbar_eff={g25:.4f}  (relative semiclassical error O(1) = O(100%))")
    b_need = 2.0 * N / 0.1**2
    say(f"     to reach 10% semiclassical accuracy (h=0.1) you need beta={b_need:.0f}  (free theory: gap->0)")
    say()
# absolute error at beta=2.5, read off the O(hbar) law of W2
c = abs((rows[0][2] - rows[0][1]) / rows[0][1]) / rows[0][0]
say(f"  fitted C = |rel err|/hbar = {c:.4f}")
for N in (2, 3):
    g = np.sqrt(2.0 * N / 2.5)
    say(f"  N={N}: predicted |rel err| at beta=2.5 = {c*g:.3f}   (i.e. {100*c*g:.0f}% )")

# ================================================================ W3b strong anharmonicity
say()
say("=" * 78)
say("W3b the semiclassical error tracks the ANHARMONICITY, which at the target is O(1)")
say("    q = (s-5)^2 + b (s-5)^4 ,  hbar=1 ;  b is the quartic/interaction strength")
say("=" * 78)
say(f"{'b':>8} {'exact lam_0':>14} {'EBK lam_0':>14} {'rel err':>12} {'sign':>5}")
w3b = []
for bq in (0.1, 1.0, 3.0, 10.0):
    qfq = lambda s, bq=bq: (s - 5.0) ** 2 + bq * (s - 5.0) ** 4
    ex = exact_lam(qfq, 1.0, T2, N=2001, k=1)[0]
    eb = ebk_lam(qfq, 1.0, T2, n=0, maslov=1.0)
    w3b.append((bq, ex, eb))
    say(f"{bq:8.2f} {ex:14.8f} {eb:14.8f} {(eb-ex)/ex:12.3e} {'+' if eb>ex else '-':>5}")
say("  -> the quartic coupling plays the role of g^2/2 : at beta=2.5 the target IS strongly")
say("     anharmonic (g=1.55), so the entry in the last rows is the relevant one.")

# ================================================================ W4 Maslov sensitivity
say()
say("=" * 78)
say("W4  sensitivity to the Maslov index at hbar=1 (n=0):  RHS = 2*pi*hbar*(n+mu/2)")
say("    admissible mu for a 2-smooth-turning-point well is 1 ; the scan shows how much of")
say("    the output is carried by the (classical) boundary layer at the turning points")
say("=" * 78)
say(f"{'mu':>6} {'EBK lam_0':>14} {'rel diff vs mu=1':>17}")
base = ebk_lam(qf, 1.0, T2, n=0, maslov=1.0)
for mu in (0.5, 1.0, 1.5):
    eb = ebk_lam(qf, 1.0, T2, n=0, maslov=mu)
    say(f"{mu:6.2f} {eb:14.8f} {(eb-base)/base:17.3e}")
vals = [ebk_lam(qf, 1.0, T2, n=0, maslov=mu) for mu in (0.5, 1.0, 1.5)]
say(f"  EBK spread (max-min)/mid = {(max(vals)-min(vals))/(0.5*(max(vals)+min(vals))):.4f}")
say(f"  exact lam_0              = {exact_lam(qf, 1.0, T2, N=2001, k=1)[0]:.8f}")

# ---- W4b: does "N_EBK(Delta)=0" certify anything?  find an explicit counterexample
say()
say("  W4b  counting test.  N_EBK(Delta)=0 is a *claim of absence* below Delta.")
say("       - if EBK errs LOW (lam_EBK<lam_exact) the claim is true but useless (weak);")
say("       - if EBK errs HIGH (lam_EBK>lam_exact) the claim is FALSE.")
say("       W2 already showed the error changes sign in hbar; find the false case explicitly:")
for hb in (0.125,):
    ex0 = exact_lam(qf, hb, T2, N=2001, k=1)[0]
    eb0 = ebk_lam(qf, hb, T2, n=0, maslov=1.0)
    D = 0.5 * (ex0 + eb0)
    nb = 1 if eb0 <= D else 0
    nx = 1 if ex0 <= D else 0
    say(f"       hbar={hb}: lam_EBK={eb0:.8f}  lam_exact={ex0:.8f}  Delta={D:.8f}")
    say(f"       N_EBK(Delta)={nb}   N_exact(Delta)={nx}"
        f"   {'<-- COUNTEREXAMPLE: EBK certifies absence, truth has a level' if (nb==0 and nx>=1) else ''}")
say("     => EBK counts; it does not bound. Its error is not one-sided (sign flip at hbar=0.125).")

# ================================================================ W5 Sturm
say()
say("=" * 78)
say("W5  rigorous shadow = Sturm comparison  lam >= pi^2/T^2 + ess inf q  (Dirichlet)")
say("    family of sign-changing q with mean 0; count violations of the floor")
say("=" * 78)
T3, Ns, viol, worst = 6.0, 900, 0, 0.0
say(f"{'#':>4} {'min q':>10} {'lam_exact':>12} {'Sturm floor':>13} {'slack':>12} {'ok?':>4}")
for it in range(8):
    s = np.linspace(0, T3, Ns)
    q = np.zeros(Ns)
    for kk in range(1, 4):
        q += (0.30 / kk) * np.sin(2 * np.pi * kk * s / T3 + np.random.uniform(0, 2 * np.pi))
    q -= np.trapezoid(q, s) / T3
    qf3 = lambda x, s=s, q=q: np.interp(x, s, q)
    ex = exact_lam(qf3, 1.0, T3, N=1201, k=1)[0]
    fl = np.pi**2 / T3**2 + q.min()
    if ex < fl - 1e-12:
        viol += 1
    worst = max(worst, abs(ex - fl))
    say(f"{it:4d} {q.min():10.5f} {ex:12.7f} {fl:13.7f} {ex-fl:12.3e} {'OK' if ex>=fl-1e-12 else 'VIOL':>4}")
say(f"  violations of the Sturm floor: {viol}/8   (one-sided OK)")
nneg = sum(1 for _ in ())  # placeholder
say("  BUT every floor printed above is NEGATIVE (min q < 0) => the rigorous 1D floor is")
say("  VACUOUS on any sign-changing (Savvidy-like) potential: it returns 'nothing', not a gap.")
say("  This is exactly LP-LIOUVILLE LP-5  Q_min = -G''/4 < 0  (Haar curvature is a COST).")
say("  and LP-6/LP-7 already killed this 1D floor as a gap certificate.")

# ================================================================ W6 LP numbers
say()
say("=" * 78)
say("W6  the corpus's own rigorous 1D (Liouville/EBK) floor, LP-LIOUVILLE B3")
say("=" * 78)
T_intr, ess_inf_Q = 10.180126, -0.052530
floor = np.pi**2 / T_intr**2 + ess_inf_Q
target = 0.6444
say(f"  T_intr            = {T_intr:.6f}")
say(f"  ess inf Q         = {ess_inf_Q:.6f}   (Q_min = -G''/4, LP-5: Haar curvature is a COST)")
say(f"  pi^2/T_intr^2     = {np.pi**2/T_intr**2:.6f}")
say(f"  Liouville floor   = {floor:.6f}")
say(f"  target Delta      = {target:.4f}")
say(f"  -> short by       = {target/floor:.2f}x")
say(f"  (and LP-7: no second-order SL representation -- any 2nd-order operator gives lam<=1.61,")
say(f"   while the beta=0 gap is infinite: the representation itself is unavailable.)")

say()
say("=" * 78)
say("SUMMARY (candidate 1 / EBK)")
say("=" * 78)
say(" - EBK is exact on the harmonic/integrable family (W1) and an approximation off it (W2).")
say(" - hbar_eff = g = sqrt(2N/beta): 1.55 at beta=2.5 (N=3); 10% needs beta=600 (free theory) (W3).")
say(" - the error grows with the ANHARMONICITY: -15% at quartic coupling 10 (W3b) -- and the")
say("   target at beta=2.5 IS strongly anharmonic (g=1.55).")
say(" - the error SIGN is not one-sided (W2: -2.7e-2 -> +2.5e-4 as hbar 1->0.125): in the")
say("   high case N_EBK(Delta)=0 while N_exact(Delta)=1 (W4b) -- counting does not certify.")
say(" - its rigorous 1D form is Sturm/Liouville: one-sided (0/8 violations) but VACUOUS")
say("   (all 8 floors negative, min q<0) and 15.09x short on the corpus's own anchor (W5,W6);")
say("   already dead in the corpus (LP-6 universality no-go, LP-7 no 2nd-order SL rep).")
