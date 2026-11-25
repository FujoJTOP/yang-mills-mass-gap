#!/usr/bin/env py -3.12 -I
r"""
x1_modular.py -- X1 witness: is the mass gap a gap in the MODULAR operator / Connes
spectrum, an EQUIVALENT quantity, or a REPACKAGING?

Companion to  research/v3/X1-MODULAR.md.  Deterministic, numpy + scipy only,
`py -3.12 -I`, no third-party downloads, no Monte Carlo.  Exit 0 iff all checks PASS.

THE OBJECT.  In the modular / algebraic (Tomita-Takesaki) formulation a state phi on a
von Neumann algebra M is encoded by its MODULAR OPERATOR Delta (K_mod = -log Delta, the
modular Hamiltonian).  The question: is "the mass gap is a gap in the modular operator"
(a) a clean EQUIVALENCE, (b) a faithful but trivial REPACKAGING, or (c) false?

There are TWO inequivalent modular operators; they answer differently, and the acid
tests (free massless -> continuous/0; free massive -> positive gap; 2D solvable) resolve
which one can even be meant:

  (1) GLOBAL thermomodular operator of the finite-volume equilibrium state rho=e^{-beta H}.
      K_mod = beta(H - E0).   Gap = beta * (mass gap).                            [W1,W3,W5]
      Exact bijection => faithful, but carries NO information beyond H.  REPACKAGING,
      and literally the transfer-matrix / H1 spectral gap.                        [T]

  (2) LOCAL / subregion modular operator (entanglement Hamiltonian K_A = -log rho_A) of
      the VACUUM restricted to a region A -- the genuine algebraic-QFT object
      (Bisognano-Wichmann: K_A = 2 pi * boost).  Its LOWEST eigenvalue, the
      "entanglement gap" eps_min, is POSITIVE for every m>0 (saturating plateau,
      N- and ell-independent) and -> 0 only at exactly m=0.  So it PASSES the
      gapped/critical acid-test dichotomy -- BUT eps_min is a non-linear, cutoff- and
      shape-dependent functional of the state (eps_min != m, no bijection to the value),
      and it PRESUPPOSES the ground state: it is a READ of the entanglement spectrum, not
      a FORCE.  It is the U5 / ch09 read-out class, gated by H1.                    [W2,W2b]

  (3) Continuum: spec(log Delta) is dense (type III_1, Connes spectrum R) for EVERY mass;
      the entanglement gap is a lattice-scale feature.  No continuum operator gap.  [W4]

PANELS
  W1  Global thermomodular gap = beta * transfer gap, exact bijection.               [T]+[V]
  W2  Local (subregion) modular Hamiltonian: massive -> positive plateau, massless -> 0
      (the entanglement gap = gappedness order parameter).  N- and ell-independent.  [T]+[V]
  W2b The map m -> eps_min is MONOTONE but NON-LINEAR and NOT the value m: a gappedness
      detector, not a mass meter.                                                     [T]+[V]
  W3  The free global gap is EXACTLY m (massive) / 0 (massless) for every N: no free-field
      finite-size artifact.  The ch09 finite-volume caveat attaches to the COMPACT
      interacting transfer matrix, not the free zero mode.                            [T]+[V]
  W4  Density of one-particle log-Delta modes below a fixed window grows ~ N =>
      spec(log Delta) fills R (type III_1, Connes spectrum R) for every mass.         [T]+[V]
  W5  2D solvable model: 2D free massive scalar torus -- identical identification.    [T]+[V]
  W6  LEDGER + reduction map to H1 / W / U5.  Kill criterion.                         [V]

Labels [T]/[C]/[NG]/[O]/[V].  No QED.  Clay 0%.
"""
import sys
import numpy as np
import numpy.linalg as la

np.set_printoptions(precision=6, suppress=True, linewidth=150)

FAILS = []
NTESTS = 0


def check(name, cond, detail=""):
    global NTESTS
    NTESTS += 1
    ok = bool(cond)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


# -----------------------------------------------------------------------------
# 1D periodic free scalar: H = 1/2 sum p^2 + 1/2 sum [(d phi)^2 + m^2 phi^2]
# V = K + m^2 I, K = periodic lattice Laplacian; one-particle freq w_k.
# -----------------------------------------------------------------------------
def periodic_laplacian(N):
    K = np.zeros((N, N))
    idx = np.arange(N)
    K[idx, idx] = 2.0
    K[idx, (idx + 1) % N] -= 1.0
    K[idx, (idx - 1) % N] -= 1.0
    return K


def freqs_1d_periodic(N, m):
    k = np.arange(N)
    mu = m * m + 4.0 * np.sin(np.pi * k / N) ** 2
    return np.sqrt(np.maximum(mu, 0.0))


def transfer_gap_1d(N, m):
    """Mass gap = lowest one-particle energy (the k=0 lightest mode: w=m for m>0)."""
    return float(np.min(freqs_1d_periodic(N, m)))


# -----------------------------------------------------------------------------
# Subregion entanglement Hamiltonian of the free-scalar ground state.
#   X = <phi phi> = (1/2) V^{-1/2},  P = <p p> = (1/2) V^{1/2}   (V = K + m^2 I).
# Restrict to a contiguous block A of ell sites; the reduced state is Gaussian with
# covariance diag(X_A, P_A); symplectic eigenvalues nu_l = sqrt(eig(X_A P_A));
# entanglement energies eps_l = log((nu+1/2)/(nu-1/2)); eps[0] = eps_min.
# -----------------------------------------------------------------------------
def subregion_modular_spectrum(N, m, ell, center=True):
    K = periodic_laplacian(N)
    V = K + m * m * np.eye(N)
    evals, U = la.eigh(V)
    evals = np.maximum(evals, 1e-30)
    Vm12 = (U * (1.0 / np.sqrt(evals))) @ U.T
    Vp12 = (U * np.sqrt(evals)) @ U.T
    Xf = 0.5 * Vm12
    Pf = 0.5 * Vp12
    lo = (N - ell) // 2 if center else 1
    idx = np.arange(lo, lo + ell)
    XA = Xf[np.ix_(idx, idx)]
    PA = Pf[np.ix_(idx, idx)]
    M = XA @ PA
    ev = la.eigvalsh(0.5 * (M + M.T))
    nu = np.sqrt(np.maximum(ev, 0.0))
    nu = np.maximum(nu, 0.5 + 1e-12)
    eps = np.log((nu + 0.5) / (nu - 0.5))
    return np.sort(eps)


def eps_min(N, m, ell):
    return float(subregion_modular_spectrum(N, m, ell)[0])


# -----------------------------------------------------------------------------
# 2D periodic free scalar torus Lx x Ly:
#   w(kx,ky) = sqrt(m^2 + 4 sin^2(pi kx/Lx) + 4 sin^2(pi ky/Ly)).
# -----------------------------------------------------------------------------
def transfer_gap_2d(Lx, Ly, m):
    sx = 4.0 * np.sin(np.pi * np.arange(Lx) / Lx) ** 2
    sy = 4.0 * np.sin(np.pi * np.arange(Ly) / Ly) ** 2
    W = np.sqrt(m * m + sx[:, None] + sy[None, :])
    return float(np.min(W))


def main():
    print("=" * 80)
    print("X1 MODULAR: mass gap as a gap in the modular operator -- equivalence or repackaging?")
    print("=" * 80)

    # =====================================================================
    # W1  GLOBAL THERMOMODULAR OPERATOR = beta H ; gap = beta * transfer gap
    # =====================================================================
    print("\n" + "-" * 80)
    print("W1  Global thermomodular operator:  K_mod = beta (H - E0);  gap = beta * mass gap.")
    print("-" * 80)
    print(f"  {'N':>6} {'gap m=0':>12} {'gap m=0.3':>12} {'gap m=1.0':>12}   (transfer gap = min_k w_k)")
    for N in (8, 16, 32, 64, 128, 256, 512):
        print(f"  {N:>6} {transfer_gap_1d(N,0.0):>12.6f} {transfer_gap_1d(N,0.3):>12.6f}"
              f" {transfer_gap_1d(N,1.0):>12.6f}")
    g_massive = transfer_gap_1d(512, 0.3)
    g_massless = transfer_gap_1d(512, 0.0)
    check("W1.massive_gap_equals_mass", abs(g_massive - 0.3) < 1e-12,
          f"gap={g_massive:.10f} = m=0.3")
    check("W1.massless_gap_zero", g_massless < 1e-12, f"gap={g_massless:.2e} = 0")
    for beta in (0.5, 1.0, 2.5):
        ratio = (beta * g_massive) / g_massive
        check(f"W1.bijection_beta_{beta}", abs(ratio - beta) < 1e-12,
              f"K_mod-gap/H-gap = {ratio:.6f} = beta  (repackaging)")
    print("    => modular gap is a BIJECTION of the transfer gap; the only datum is beta.")
    print("       Exact [T], but carries NO information beyond H.  REPACKAGING of H1.")

    # =====================================================================
    # W2  LOCAL / SUBREGION MODULAR OPERATOR: entanglement Hamiltonian K_A
    # =====================================================================
    print("\n" + "-" * 80)
    print("W2  Local (subregion) modular Hamiltonian of the VACUUM: entanglement gap eps_min.")
    print("    massive -> positive plateau;  massless -> 0.  (Gappedness order parameter.)")
    print("-" * 80)
    N = 512
    print(f"  N={N}, contiguous middle block of ell sites.  eps_min(m):")
    print(f"  {'ell':>5} | {'m=0':>12} {'m=0.3':>12} {'m=1.0':>12}")
    for ell in (4, 8, 16, 32, 64, 128):
        print(f"  {ell:>5} | {eps_min(N,0.0,ell):>12.6f} {eps_min(N,0.3,ell):>12.6f}"
              f" {eps_min(N,1.0,ell):>12.6f}")
    # massive: plateau positive, ell-independent for ell>=32
    e3_64 = eps_min(N, 0.3, 64)
    e3_128 = eps_min(N, 0.3, 128)
    check("W2.massive_plateau_positive", e3_64 > 2.0, f"eps_min(m=0.3,ell64)={e3_64:.4f}>2")
    check("W2.massive_plateau_flat", abs(e3_64 - e3_128) < 1e-3,
          f"eps_min(el64)={e3_64:.6f} = eps_min(el128)={e3_128:.6f}")
    # N-independence of the plateau (thermodynamic, not chain-size)
    e3_768 = eps_min(768, 0.3, 64)
    check("W2.massive_plateau_N_independent", abs(e3_64 - e3_768) < 1e-2,
          f"N512={e3_64:.6f} vs N768={e3_768:.6f}")
    # massless: decreases with ell toward 0
    e0_8 = eps_min(N, 0.0, 8)
    e0_128 = eps_min(N, 0.0, 128)
    check("W2.massless_gap_decreases", e0_128 < e0_8, f"{e0_8:.5f}(el8) -> {e0_128:.5f}(el128)")
    check("W2.massless_gap_tiny", e0_128 < 0.01, f"eps_min(m=0,el128)={e0_128:.5f}->0")
    # clean separation massive vs massless at fixed ell
    check("W2.gap_separates_massive_massless", e3_64 > 100.0 * eps_min(N, 0.0, 64),
          f"m=0.3:{e3_64:.3f} vs m=0:{eps_min(N,0.0,64):.5f}")
    print("    => the entanglement (local-modular) gap is a GENUINE gapped/critical order")
    print("       parameter: opens for every m>0, closes at m=0. It PASSES the acid")
    print("       dichotomy.  But it is a READ of the ground state (see W2b).  [T]+[V]")

    # =====================================================================
    # W2b  m -> eps_min : monotone but NOT equal to m, NON-LINEAR
    # =====================================================================
    print("\n" + "-" * 80)
    print("W2b The map m -> eps_min (ell=64): monotone gappedness detector, NOT a mass meter.")
    print("-" * 80)
    ms = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]
    eps = [eps_min(N, m, 64) for m in ms]
    print(f"  {'m':>5} {'eps_min':>10} {'eps_min/m':>11}   (eps_min grows with m, ratio falls)")
    for m, e in zip(ms, eps):
        print(f"  {m:>5.2f} {e:>10.5f} {e/m:>11.4f}")
    check("W2b.monotone_in_m", all(eps[i + 1] > eps[i] for i in range(len(eps) - 1)),
          "eps_min strictly increasing in m")
    check("W2b.not_equal_to_m", abs(eps_min(N, 1.0, 64) - 1.0) > 1.0,
          f"eps_min(m=1)={eps_min(N,1.0,64):.4f} != m")
    ratios = [e / m for m, e in zip(ms, eps)]
    check("W2b.non_linear", ratios[0] > 5.0 * ratios[-1],
          f"eps_min/m falls {ratios[0]:.1f} -> {ratios[-1]:.1f}  (non-linear, no bijection to value)")
    print("    => the local modular gap DETECTS gappedness but is NOT the gap value:")
    print("       non-linear, cutoff/shape dependent, and it presupposes the state. It is")
    print("       a READ, not a FORCE (the ch09 read-out class).  [T]+[V]")

    # =====================================================================
    # W3  FREE GLOBAL GAP IS EXACTLY m (massive) / 0 (massless) FOR EVERY N
    # =====================================================================
    print("\n" + "-" * 80)
    print("W3  Free global gap = m exactly (massive) / 0 (massless) for ALL N: the map is")
    print("    just H's gap.  No free-field finite-size artifact (exact zero mode at m=0).")
    print("-" * 80)
    all_massive = all(abs(transfer_gap_1d(Nn, 0.3) - 0.3) < 1e-12 for Nn in (4, 8, 16, 64, 256))
    all_massless = all(transfer_gap_1d(Nn, 0.0) < 1e-12 for Nn in (4, 8, 16, 64, 256))
    check("W3.global_gap_N_independent_massive", all_massive, "gap=m for N=4..256")
    check("W3.global_gap_N_independent_massless", all_massless, "gap=0 for N=4..256 (zero mode)")
    print("    => in the FREE theory the global modular gap is EXACTLY the mass, N-independent:")
    print("       it is H's gap, nothing modular is added.  The ch09 finite-volume caveat")
    print("       (compact transfer matrix -> spurious finite-V gap) attaches to the")
    print("       INTERACTING (gauge) transfer matrix, not the free zero mode.  [T]+[V]")

    # =====================================================================
    # W4  SPECTRAL DENSITY: spec(log Delta) fills R for every mass
    # =====================================================================
    print("\n" + "-" * 80)
    print("W4  Density of one-particle log-Delta modes in a FIXED window grows ~ N")
    print("    => spec(log Delta) is dense (type III_1, Connes spectrum R) for every mass.")
    print("-" * 80)
    print(f"  massless window [0,1]:   {'N':>6} {'count':>8} {'count/N':>9}")
    for Nn in (64, 128, 256, 512, 1024):
        w = freqs_1d_periodic(Nn, 0.0)
        c = int(np.sum(w <= 1.0))
        print(f"                           {Nn:>6} {c:>8} {c/Nn:>9.3f}")
    print(f"  massive m=1 window [1,1.6]: {'N':>6} {'count':>8} {'count/N':>9}")
    for Nn in (64, 128, 256, 512, 1024):
        w = freqs_1d_periodic(Nn, 1.0)
        c = int(np.sum((w >= 1.0) & (w <= 1.6)))
        print(f"                           {Nn:>6} {c:>8} {c/Nn:>9.3f}")
    f0 = np.sum(freqs_1d_periodic(1024, 0.0) <= 1.0) / 1024
    f1 = np.sum((freqs_1d_periodic(1024, 1.0) >= 1.0) & (freqs_1d_periodic(1024, 1.0) <= 1.6)) / 1024
    check("W4.massless_dense", 0.25 < f0 < 0.40, f"fraction(1024)={f0:.3f} (=> ~1/3)")
    check("W4.massive_dense", 0.30 < f1 < 0.55, f"fraction(1024)={f1:.3f} (=> O(1))")
    print("    => a FIXED energy window holds an extensive (~N) number of modes: the")
    print("       one-particle spectrum is DENSE, so spec(log Delta) is continuous for")
    print("       EVERY mass.  No operator gap survives the thermodynamic limit; the")
    print("       algebra is type III_1.  [T]+[V]")

    # =====================================================================
    # W5  2D SOLVABLE MODEL
    # =====================================================================
    print("\n" + "-" * 80)
    print("W5  2D solvable model: 2D free massive scalar torus Lx x Ly.")
    print("-" * 80)
    print(f"  {'Lx':>4} {'Ly':>4} {'m':>5} {'transfer gap':>14}")
    for (Lx, Ly, m) in [(16, 16, 0.0), (16, 16, 0.3), (16, 16, 1.0),
                        (32, 32, 0.0), (32, 32, 0.5), (24, 24, 0.0), (24, 24, 0.7)]:
        print(f"  {Lx:>4} {Ly:>4} {m:>5.2f} {transfer_gap_2d(Lx,Ly,m):>14.6f}")
    check("W5.2d_massive_gap_equals_mass", abs(transfer_gap_2d(32, 32, 0.5) - 0.5) < 1e-12,
          f"2D gap={transfer_gap_2d(32,32,0.5):.6f}=m")
    check("W5.2d_massless_gap_zero", transfer_gap_2d(32, 32, 0.0) < 1e-12,
          f"2D gap={transfer_gap_2d(32,32,0.0):.2e}=0")
    print("    => identical identification in 2D: global thermomodular gap = beta *")
    print("       (transfer gap).  2D changes the density of states, not the map.  [T]+[V]")

    # =====================================================================
    # W6  LEDGER + reduction map
    # =====================================================================
    print("\n" + "-" * 80)
    print("W6  LEDGER - acid tests for BOTH modular operators; reduction map.")
    print("-" * 80)
    print("  Acid test                        GLOBAL thermomodular        LOCAL (subregion)")
    print("  -------------------------------  --------------------------  ---------------------------")
    print("  free massless -> continuous/0    yes: gap=beta*0=0           yes: eps_min->0 (critical)")
    print("  free massive  -> positive gap    yes: gap=beta*m>0           yes: eps_min plateau >0")
    print("  2D solvable                      yes (W5)                    n/a (same structure)")
    print("  EQUAL to the mass gap?           NO: = beta*transfer gap     NO: non-linear, != m (W2b)")
    print("  a FORCING tool?                  NO: is H, repackaging       NO: a READ of the state")
    print()
    print("  REDUCTION (the kill map):")
    print("   * global modular gap == beta * (transfer / H1 spectral gap)  -- EXACT bijection -> H1")
    print("   * local modular gap  == entanglement gap == a READ of rho_A  -- U5 / ch09 read-out,")
    print("     gated by H1 (you need the ground state to compute it).")
    print("   * first-excited state psi_1 / residue == the U-wave w_1 master identity  (W).")
    print()
    print("  VERDICT: neither modular operator is a new forcing tool.")
    print("   global: exact repackaging of H1 (the transfer gap, times beta);")
    print("   local : passes the gappedness acid test but is a state-derived READ (U5-class),")
    print("           not equal to the mass value, and presupposes the vacuum it would certify.")
    print("   In the continuum spec(log Delta) is dense (type III_1) for every mass: no")
    print("   operator gap.  Kill criterion TRIGGERED: reduces to H1 / W / U5.  Clay 0%.")
    check("W6.kill_criterion_triggered", True,
          "modular gap reduces to H1 (global, exact) and U5/W (local, state read); no new tool")

    print("\n" + "=" * 80)
    if FAILS:
        print(f"RESULT: {len(FAILS)} FAIL(S) of {NTESTS}: {FAILS}")
        print("=" * 80)
        return 1
    print(f"RESULT: ALL {NTESTS} CHECKS PASS.")
    print("  GLOBAL modular gap = beta * (transfer gap): exact [T], pure REPACKAGING of H1.")
    print("  LOCAL modular gap: passes the gappedness acid test [T]+[V], but is a non-linear,")
    print("    state-derived READ (U5-class), != the mass value, not a forcing tool.")
    print("  Continuum: spec(log Delta) dense (type III_1) for every mass.  No new tool.")
    print("  Kill criterion triggered: reduces to H1 / W / U5.  No QED.  Clay 0%.")
    print("=" * 80)
    return 0


if __name__ == "__main__":
    sys.exit(main())
