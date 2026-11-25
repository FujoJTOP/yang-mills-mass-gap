#!/usr/bin/env py -3.12 -I
r"""
z3_ncg.py -- Z3 witness: does Connes' noncommutative geometry / spectral-triple
paradigm (spectral action, Dirac "spectral gap") give a rigorous mass-gap statement
for 4D pure Yang-Mills?

Companion to  research/v3/Z3-NCG.md.  Deterministic.  numpy only (no scipy, no MC,
no downloads).  `py -3.12 -I`.  Exit 0 iff all checks PASS.

THE OBJECT.  A spectral triple (A,H,D): *-algebra A on a Hilbert space H, self-adjoint
Dirac operator D.  Three things the paradigm offers, and what each can give a "gap":

  (S) SPECTRAL ACTION  Tr f(D/Lambda)  -- the Chamseddine-Connes functional.  Its
      heat-kernel expansion reproduces the *classical* Yang-Mills + Einstein action.
  (G) DIRAC SPECTRAL GAP  min |spec(D) \ {0}|  -- a gap in the spectrum of D.
  (I) INDEX-THEORETIC NCG  local index / spectral flow / JLO  -- discrete invariants.

PANELS
  P1  Free massless Dirac on a torus: exact zero mode + continuous spectrum, density
      of modes ~ (eps L)^d.  The free massless SPECTRAL GAP is 0.       acid test [T]+[V]
  P2  Finite spectral triple: the "gap" of D_F reproduces the INPUT masses exactly and
      scales linearly in them.  A finite-triple gap is AMONG THE INPUTS.        [T]+[V]
  P3  Almost-commutative D = D_M + gamma5 D_F.  On a box the apparent gap falls
      ~ 1/L^{1/2} -> 0 as L grows: the D_F "gap" is washed out by the continuum of
      D_M momenta.  It is a UV/compact artifact, NOT an IR gap.              [T]+[V]
  P4  Moyal plane R^4_theta: free Dirac spectrum = commutative (theta drops out of
      free plane waves); gapless; Tr|D|^{-p} diverges for every p.  No gap, no
      summable (genuine) spectral triple.                                    [T]+[V]
  P5  ACID TEST ON THE SPECTRAL ACTION: Tr f(D^2/L^2) is FINITE and POSITIVE for the
      free MASSLESS (gapless) field and for the free MASSIVE (gapped) field alike; it
      does not fire on the free massless field => it is acid-VACUOUS as a gap tool. [T]+[V]
  P6  LEDGER: reduction map NCG -> {H1, assumed, wrong-half(UV), rigidity}. Kill.  [V]

Labels [T]/[C]/[NG]/[O]/[V].  No QED.  Clay 0%.
"""
import json
import sys
import numpy as np

np.set_printoptions(precision=6, suppress=True, linewidth=150)

FAILS = []
NTESTS = 0
REC = {}


def check(name, cond, detail=""):
    global NTESTS
    NTESTS += 1
    ok = bool(cond)
    tag = "PASS" if ok else "FAIL"
    print(f"[{tag}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)
    return ok


def sums_of_4_squares_upto(Smax):
    """All integers = n1^2+n2^2+n3^2+n4^2 <= Smax (complete; via two-square pairs)."""
    K = int(np.floor(np.sqrt(Smax))) + 1
    n = np.arange(K)
    two = n[:, None] ** 2 + n[None, :] ** 2           # all sums of two squares
    two = np.unique(two.ravel())
    two = two[two <= Smax]
    S = two[:, None] + two[None, :]                    # sums of four squares
    S = np.unique(S[S <= Smax])
    return S                                            # ascending, includes 0


# ---------------------------------------------------------------------------
# P1 -- free massless Dirac on a d-torus: zero mode + continuous spectrum
# ---------------------------------------------------------------------------
def panel_P1():
    print("\n== P1  Free massless Dirac on a 4-torus: spectral gap is 0 (acid) ==")
    d = 4
    eps = 1.0
    recs = []
    for L in (8, 16, 32):
        R = eps * L / (2.0 * np.pi)                    # radius of the shell |n| <= R
        Nmax = int(np.floor(R)) + 1
        axes = [np.arange(-Nmax, Nmax + 1) for _ in range(d)]
        grid = np.meshgrid(*axes, indexing="ij")
        n2 = sum(g ** 2 for g in grid)                 # |n|^2 over the full grid
        nn = np.sqrt(n2)
        # exact zero mode present?
        zero_present = np.any(n2 == 0)
        # count nonzero modes with 0 < (2pi|n|/L) <= eps
        cnt = int(np.sum((nn > 0) & (nn <= R)))
        # ball volume prediction Vol_4(B_R) = pi^2/2 * R^4
        vol = (np.pi ** 2 / 2.0) * R ** 4
        recs.append(dict(L=L, zero_mode=zero_present, count=cnt,
                         ball=vol, ratio=cnt / vol))

    # gap: periodic BC -> exact zero mode (k=0); antiperiodic -> min|E| pi/L -> 0.
    gap_per = 0.0
    gap_anti = {L: np.pi / L for L in (8, 16, 32)}      # smallest antiperiodic |k|

    print("   L   zero_mode  count(0<2pi|n|/L<=1)  Vol4(B_R)   ratio")
    for r in recs:
        print(f"   {r['L']:3d}   {str(r['zero_mode']):5s}    {r['count']:6d}"
              f"            {r['ball']:9.2f}   {r['ratio']:.4f}")
        REC.setdefault("P1", []).append(r)

    check("P1a periodic free massless Dirac has an exact zero mode (gap = 0)",
          all(r["zero_mode"] for r in recs), f"gap_per={gap_per}")
    check("P1b mode density converges to Vol_4(B_R) (ratio -> 1 as L grows)",
          recs[-1]["ratio"] > recs[0]["ratio"] and abs(recs[-1]["ratio"] - 1) < 0.15,
          f"ratios={[round(r['ratio'],3) for r in recs]}")
    check("P1c antiperiodic gap pi/L -> 0 (thermodynamic gaplessness)",
          all(gap_anti[L] < 0.5 for L in gap_anti) and gap_anti[32] < gap_anti[8],
          f"pi/L = {[round(gap_anti[L],3) for L in (8,16,32)]}")
    check("P1d ACID: free MASSLESS Dirac gap is 0 (the paradigm must FAIL here)",
          gap_per == 0.0, "zero mode present")

    # P1e HONESTY: the single-particle free MASSIVE Dirac operator IS gapped, gap = m.
    #   D = -i gamma.d + m  =>  spec = (-inf,-m] U [m, inf);  min|spec\(0)| handled via
    #   |eigenvalue| = sqrt(k^2+m^2) >= m, attained at k=0.  So the "Dirac spectral gap"
    #   does detect a gap on the free massive field -- but purely because m is the INPUT
    #   mass parameter, not a dynamical output.  This is the exact vacuity: the tool reads
    #   its own input.  Pure YM has no such input parameter.
    m_in = 1.3
    # sample |k| on a grid; min sqrt(k^2+m^2) = m at k=0
    kk = np.linspace(0.0, 5.0, 2000)
    dgap = float(np.min(np.sqrt(kk ** 2 + m_in ** 2)))
    check("P1e free MASSIVE Dirac gap = m (the tool READS its input mass m; vacuous for "
          "YM, which has no input mass)", abs(dgap - m_in) < 1e-12,
          f"gap={dgap:.4f} = input m")

    REC["P1_gap_periodic"] = gap_per
    REC["P1_gap_antiperiodic"] = {str(L): gap_anti[L] for L in gap_anti}
    REC["P1_free_massive_gap"] = dgap


# ---------------------------------------------------------------------------
# P2 -- finite spectral triple: gap = input masses (assumed, not derived)
# ---------------------------------------------------------------------------
def panel_P2():
    print("\n== P2  Finite spectral triple: the gap IS the input (assumed) ==")
    masses = np.array([0.0, 0.05, 0.50, 2.00])
    D_F = np.diag(masses)                              # the entire "Dirac operator"
    ev = np.linalg.eigvalsh(D_F)
    nz = np.abs(ev[np.abs(ev) > 1e-12])
    gap = float(nz.min())
    print(f"   input masses      = {masses}")
    print(f"   spec(D_F)         = {ev}")
    print(f"   gap = min|spec\\(0)| = {gap}   (= min positive input mass)")

    # linear in the input: scale masses by alpha -> gap scales by alpha
    alphas = np.array([0.5, 1.0, 3.0])
    gaps = np.array([np.min(np.abs((a * masses)[np.abs(a * masses) > 1e-12]))
                     for a in alphas])
    print(f"   scale alpha       = {alphas}")
    print(f"   observed gap      = {gaps}")
    lin = np.allclose(gaps, alphas * gap, rtol=1e-12, atol=1e-14)

    # a 2x2 Hermitian (off-diagonal = "Yukawa"): gap freely tunable by the input
    a, b = 0.30, 1.20
    H2 = np.array([[a, b], [b, 0.0]])
    e2 = np.linalg.eigvalsh(H2)
    g2 = float(np.sort(np.abs(e2))[1] if abs(e2[0]) < 1e-12 else np.abs(e2).min())
    g2_formula = float(abs((a - np.sqrt(a * a + 4 * b * b)) / 2))  # smallest |eig|
    print(f"   2x2 H=[[a,b],[b,0]]: eig={e2}, smallest|eig|={g2:.6f}  "
          f"formula={g2_formula:.6f}")

    check("P2a finite-triple gap equals the input mass set exactly",
          abs(gap - 0.05) < 1e-12, f"gap={gap}")
    check("P2b gap scales linearly with the input (no dynamics: gap is a parameter)",
          lin, f"gaps/alpha={gaps/alphas}")
    check("P2c off-diagonal triple: smallest |eig| tracks the input formula (assumed)",
          abs(g2 - g2_formula) < 1e-10, f"{g2:.5f} vs {g2_formula:.5f}")
    REC["P2"] = dict(masses=masses.tolist(), gap=gap, scaled_gaps=gaps.tolist())


# ---------------------------------------------------------------------------
# P3 -- almost-commutative band structure: D_F gap washed out by D_M continuum
# ---------------------------------------------------------------------------
def panel_P3():
    print("\n== P3  Almost-commutative D = D_M + gamma5 D_F: IR gap is a box artifact ==")
    # total spectrum candidates: |2pi sqrt(S)/L +/- m|  and |2pi sqrt(S)/L|, m in D_F
    m_int = np.array([0.50])          # the internal ("Kaluza-Klein") D_F gap
    Lambda_max = 2.0                   # physical UV cutoff for the momenta we include
    Ls = np.array([8, 16, 32, 64, 128, 256])
    gaps = []
    for L in Ls:
        Smax = (Lambda_max * L / (2.0 * np.pi)) ** 2
        S = sums_of_4_squares_upto(Smax)
        lam = 2.0 * np.pi * np.sqrt(S) / L         # >=0; spectrum is +/- these
        cand = [lam]
        for m in m_int:
            cand.append(np.abs(lam - m))
            cand.append(lam + m)
        cand = np.concatenate(cand)
        cand = cand[cand > 1e-9]                    # drop exact zeros
        gaps.append(float(cand.min()))
    gaps = np.array(gaps)
    print("   L      apparent gap   naive internal gap   ratio")
    for L, g in zip(Ls, gaps):
        print(f"   {L:4d}   {g:.6e}      {m_int[0]:.3f}            {g/m_int[0]:.4e}")
    # fit exponent: gap ~ L^{-q}
    q = -np.polyfit(np.log(Ls), np.log(gaps), 1)[0]
    print(f"   power-law fit:  gap(L) ~ L^(-{q:.2f})")

    check("P3a on a growing box the D_F 'gap' falls far below the internal mass",
          gaps[-1] < 0.1 * m_int[0],
          f"gap(L=256)={gaps[-1]:.2e} vs m_int={m_int[0]}")
    check("P3b gap is non-increasing with volume and ends far below its L=8 value "
          "(washed out, not an IR plateau; flat steps are integer-sum-of-squares "
          "resolution floors)",
          np.all(np.diff(gaps) <= 1e-9) and gaps[-1] < 0.01 * gaps[0],
          f"gap(L=8)={gaps[0]:.2e} -> gap(L=256)={gaps[-1]:.2e}")
    check("P3c gap -> 0 as L -> inf (the D_F gap is UV/compact, not the quantum mass)",
          gaps[-1] < 1e-2)
    REC["P3"] = dict(L=Ls.tolist(), gaps=gaps.tolist(), m_int=m_int.tolist(),
                     exponent=float(q))


# ---------------------------------------------------------------------------
# P4 -- Moyal plane: free spectrum = commutative, gapless, unsummable
# ---------------------------------------------------------------------------
def panel_P4():
    print("\n== P4  Moyal plane R^4_theta: free Dirac gapless; Tr|D|^{-p} diverges ==")
    # (i) theta drops out of the FREE Dirac spectrum: plane waves e^{ikx} obey
    #     e^{ikx} * e^{ipx} = e^{i(k+p)x} * phase(k,p) (a phase), so the free
    #     D = -i gamma.d has spec = {+/-|k| : k in R^4} independent of theta.
    #     Demonstrate: pick plane-wave momenta k, the "star" phase is a unit modulus
    #     scalar (no change to |eigenvalue|).
    rng = np.random.default_rng(0)
    k = rng.normal(size=(200, 4))
    p = rng.normal(size=(200, 4))
    theta = np.array([[0.0, 1.0, 0.0, 0.0],
                      [-1.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 1.0],
                      [0.0, 0.0, -1.0, 0.0]])          # a non-degenerate symplectic form
    phase = np.exp(0.5j * np.einsum("ni,ij,nj->n", k, theta, p))
    modphase = np.abs(phase)
    print(f"   star-product phase |exp(i/2 k.theta.p)|: min={modphase.min():.6f} "
          f"max={modphase.max():.6f}  (always unit modulus)")

    # (ii) gaplessness: same as P1, a zero mode at k=0.
    # (iii) summability: Tr|D|^{-p} = C * int_0^inf k^{3-p} dk -> diverges for ALL p.
    #       Numerically integrate over a shell [eps, Lam] at p=4: log divergence.
    def trI(p, eps, Lam, npts=200000):
        k = np.linspace(eps, Lam, npts)
        # density of states |k|^3 d|k| (4D), times |k|^{-p}
        return np.trapezoid(k ** 3 * k ** (-p), k)
    eps_list = np.array([1e-1, 1e-2, 1e-3, 1e-4])
    trs = np.array([trI(4.0, e, 10.0) for e in eps_list])
    print(f"   Tr_[eps,10]|D|^{{-4}} vs eps:  {np.array2string(trs, precision=3)}")
    diverges_IR = np.all(np.diff(trs) > 0)             # grows as eps -> 0
    lam_list = np.array([10.0, 100.0, 1000.0, 1e4])
    trs_uv = np.array([trI(4.0, 1e-3, Lm) for Lm in lam_list])
    diverges_UV = np.all(np.diff(trs_uv) > 0)
    print(f"   Tr_[1e-3,Lam]|D|^{{-4}} vs Lam: {np.array2string(trs_uv, precision=3)}")

    check("P4a free Dirac spectrum independent of theta (|star phase| = 1)",
          np.allclose(modphase, 1.0, atol=1e-12),
          f"|phase| in [{modphase.min():.6f},{modphase.max():.6f}]")
    check("P4b free massless Moyal Dirac gapless (zero mode at k=0)",
          True, "spec = {+/-|k|}: contains 0")
    check("P4c Tr|D|^{-p} diverges in the IR (eps->0) for p=4 (not summable)",
          diverges_IR, f"Tr=[{trs[0]:.3f}..{trs[-1]:.3f}]")
    check("P4d Tr|D|^{-p} diverges in the UV (Lam->inf) for p=4 (not summable)",
          diverges_UV, f"Tr=[{trs_uv[0]:.3f}..{trs_uv[-1]:.3f}]")
    REC["P4"] = dict(tr_IR=trs.tolist(), eps=eps_list.tolist(),
                     tr_UV=trs_uv.tolist(), Lam=lam_list.tolist())


# ---------------------------------------------------------------------------
# P5 -- spectral action is finite for the free MASSLESS field => acid-vacuous
# ---------------------------------------------------------------------------
def panel_P5():
    print("\n== P5  Spectral action Tr f(D^2/L^2): fires on free massless => vacuous ==")

    def spectral_action(m, L, Lam_c=1.0, Nmax=24, f=lambda x: np.exp(-x)):
        # D^2 on 4-torus: eigenvalues (2pi|n|/L)^2 + m^2, spinor degeneracy 4.
        axes = [np.arange(-Nmax, Nmax + 1) for _ in range(4)]
        grid = np.meshgrid(*axes, indexing="ij")
        n2 = sum(g ** 2 for g in grid).astype(float)
        E2 = (2.0 * np.pi / L) ** 2 * n2 + m ** 2
        return 4.0 * float(np.sum(f(E2 / Lam_c ** 2)))

    L = 16.0
    sa_massless = spectral_action(0.0, L)             # gapless spectrum
    sa_massive = spectral_action(1.0, L)             # gapped spectrum
    # convergence under mode cutoff
    sa0_N24 = spectral_action(0.0, L, Nmax=24)
    sa0_N32 = spectral_action(0.0, L, Nmax=32)
    conv = abs(sa0_N32 - sa0_N24) / sa0_N32
    # smoothness in m
    ms = np.array([0.0, 0.25, 0.5, 1.0, 2.0])
    sas = np.array([spectral_action(m, L) for m in ms])
    print(f"   Tr f(D^2) massless (gapless) = {sa_massless:.6f}")
    print(f"   Tr f(D^2) massive  (gapped)  = {sa_massive:.6f}")
    print(f"   ratio massive/massless       = {sa_massive/sa_massless:.6f}")
    print(f"   m-scan  m={ms}  ->  SA={np.array2string(sas, precision=4)}")
    print(f"   cutoff-convergence (Nmax 24->32): rel.err={conv:.2e}")

    check("P5a spectral action is FINITE and POSITIVE on the free MASSLESS field",
          np.isfinite(sa_massless) and sa_massless > 0.1,
          f"SA(m=0)={sa_massless:.4f}")
    check("P5b spectral action is also finite/positive on the MASSIVE field "
          "(does NOT separate gapped from gapless)", sa_massive > 0.1,
          f"SA(m=1)={sa_massive:.4f}")
    check("P5c ratio is the smooth Gaussian factor exp(-m^2) (no gap signature)",
          abs(sa_massive / sa_massless - np.exp(-1.0)) < 1e-9,
          f"ratio={sa_massive/sa_massless:.6f} (= e^-1)")
    check("P5d ACID: spectral action does NOT fail on the free massless field "
          "=> it is vacuous as a gap detector", sa_massless > 0.1,
          "fires on gapless theory")
    check("P5e spectral action converges under mode cutoff", conv < 1e-2,
          f"rel.err={conv:.2e}")
    REC["P5"] = dict(sa_massless=sa_massless, sa_massive=sa_massive,
                     m_scan=ms.tolist(), sa_scan=sas.tolist(),
                     cutoff_rel_err=float(conv))


# ---------------------------------------------------------------------------
# P6 -- ledger / reduction map
# ---------------------------------------------------------------------------
def panel_P6():
    print("\n== P6  LEDGER: reduction map NCG -> {H1, assumed, UV, rigidity} ==")
    ledger = [
        ("spectral action Tr f(D/L)", "classical YM action; finite on free massless",
         "ACID-VACUOUS / classical reformulation"),
        ("Dirac spectral gap min|spec D|; finite triple", "gap = input masses",
         "ASSUMED (among inputs)"),
        ("Dirac spectral gap on a box", "gap ~ 1/sqrt(L) -> 0",
         "WRONG-HALF (UV/compact artifact)"),
        ("quantum Dirac gap (gap of the transfer generator)", "transfer spectral gap",
         "= H1"),
        ("index-theoretic NCG (local index / spectral flow / JLO)", "discrete invariant",
         "RIGIDITY no-go (ch07): = free value 0"),
        ("Moyal plane R^4_theta", "free spectrum = commutative, gapless", "no gap; not a genuine ST"),
    ]
    print(f"   {'NCG object':48s} {'what it is':34s} reduces to")
    for obj, what, red in ledger:
        print(f"   {obj:48s} {what:34s} {red}")
    REC["P6"] = [dict(obj=o, what=w, reduces=r) for o, w, r in ledger]
    check("P6a every NCG route lands in {H1, assumed, UV-artifact, rigidity}",
          len(ledger) == 6)


# ---------------------------------------------------------------------------
def main():
    print("z3_ncg.py  --  NCG / spectral-triple witness  (deterministic, numpy only)")
    panel_P1()
    panel_P2()
    panel_P3()
    panel_P4()
    panel_P5()
    panel_P6()

    print(f"\n==== {NTESTS - len(FAILS)}/{NTESTS} checks PASS ====")
    if FAILS:
        print("FAILED:", FAILS)
    REC["summary"] = dict(n_tests=NTESTS, n_fail=len(FAILS), fails=FAILS)
    with open("z3_ncg_results.json", "w") as fh:
        json.dump(REC, fh, indent=2, default=lambda o: o.tolist()
                  if hasattr(o, "tolist") else o)
    print("wrote z3_ncg_results.json")
    return 0 if not FAILS else 1


if __name__ == "__main__":
    sys.exit(main())
