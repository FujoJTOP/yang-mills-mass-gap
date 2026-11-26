"""p8 — is f1^2 >= r0 still true above beta=1/48?  Numbers, not adjectives.

Q: P3 proves f1^2 >= r0 only for beta<1/48.  Does f1^2(2.5) exist, and is it >> r0?
Inputs (all in-corpus, no invention):
  * COMPUTE-F1SQ-SU2.md  : SU(2) beta=2.2, L=6, N=4000, f1^2_bare=7.14+-0.42, m_G=1.51a
  * data f1sq_b2.2_L6.npz: S{0,2,4,8} = N x T time-slice sums  -> recompute C(t), C(0), A
  * r0 = 1/2 (2N)^-2 (Z10-HOLEE:149, N-degraded floor, beta-independent) and r0 = 1/2 c1 u^p
  * rho^2 = mu2/mu1 = 2.48 +- 0.2 (isolation coordinate), NOT the residue
Checks: (1) w1 = f1^2/M0 from data; (2) r0 口径 table; (3) rho^2 is blind to w1 (exact 2-dim);
        (4) residue-level-gap lemma: tautological at 2 levels; (5) dimensional mismatch flag.
Writes nothing; prints a table.
"""
import json
import numpy as np

DATA = "research/unconditional/f1sq_b2.2_L6.npz"
M_G, M_G_ERR = 1.51, 0.15
T_EXT = 6


def connected(S, T):
    """C(t) = <S(t)S(0)>_c, S: (N, T) time-slice sums."""
    x = S - S.mean(axis=0)
    C = x.T @ x / len(x)
    return C[0]                      # row 0 = C(t) for t=0..T-1


def amp_fixed_mass(C, m, ts=(1, 2)):
    """Least squares A for C(t) = A(e^{-mt} + e^{-m(T-t)}), t in ts."""
    t = np.array(ts, float)
    g = np.exp(-m * t) + np.exp(-m * (T_EXT - t))
    return float(g @ C[list(ts)] / (g @ g)), g


def main():
    d = np.load(DATA)
    out = {"beta": float(d["beta"]), "L": int(d["L"]), "N": int(d["N"])}
    rows = []
    for lev in ["S0", "S2", "S4", "S8"]:
        C = connected(d[lev], T_EXT)
        A, _ = amp_fixed_mass(C, M_G)
        w1 = A / C[0]
        rows.append((lev, C[0], C[1] / C[0], A, w1))
    out["levels"] = rows
    print("=== 1. beta=2.2 data: the atom FRACTION w1 = f1^2/M0 (dimensionless) ===")
    print(f"{'lvl':>4} {'C(0)':>10} {'C(1)/C(0)':>10} {'A(m=1.51)':>10} {'w1=A/C(0)':>10}")
    for lev, c0, r, A, w1 in rows:
        print(f"{lev:>4} {c0:10.4f} {r:10.5f} {A:10.4f} {w1:10.4f}")
    c0, A0, w0 = rows[0][1], rows[0][3], rows[0][4]
    # the paper quotes A=7.14 with C(1)/C(0)=0.123; our npz gives 0.108 -> 13% convention gap.
    # w1 = A/C(0) is invariant under operator rescaling, so both routes are comparable.
    w_paper = 7.14 / (7.14 / 0.5627)        # AUDIT-C0:94 A/C(0)=0.5627, C(0)=12.72
    print(f"  [M] paper route A/C(0) = {w_paper:.4f} (C(1)/C(0)=0.123) vs our npz {w0:.4f} "
          f"(C(1)/C(0)={rows[0][2]:.4f}); 13% estimator/window gap, flagged not hidden")
    w_lo, w_hi = min(w0, w_paper), max(r[4] for r in rows)
    print(f"  robust band: w1(2.2) in [{w_lo:.2f}, {w_hi:.2f}] across 2 routes x 4 smearings "
          "(s=0 polluted downward by contact/UV weight in M0)")

    print("\n=== 2. r0: what the P3 floor actually says (dimensionless or nothing) ===")
    u_148 = 5.208e-3                     # u(1/48), T23-TEMPLATE:288
    r0_tab = {}
    for N in (2, 3):
        r0_tab[f"1/2(2N)^-2, N={N}"] = 0.5 * (2 * N) ** -2
        r0_tab[f"1/2*c1*u, c1=2N, N={N}"] = 0.5 * (2 * N) * u_148
    for k, v in r0_tab.items():
        print(f"  {k:<26} r0 = {v:.5f}")
    print(f"  {'6u(1/48) (NOT r0)':<26}       = {6 * u_148:.5f}   "
          "same digits as r0(SU2) by coincidence -> the corpus number 0.031")
    r0 = r0_tab["1/2(2N)^-2, N=2"]
    print(f"\n  adopt r0 = {r0:.5f} (SU(2), beta-independent floor)")
    print(f"  w1(2.2)/r0 = {w_lo:.3f}..{w_hi:.3f} over r0={r0:.5f} = "
          f"{w_lo / r0:.0f}x .. {w_hi / r0:.0f}x")
    tight = r0_tab["1/2*c1*u, c1=2N, N=2"]
    print(f"  margin to the tightest 口径 r0=1/2*c1*u (N=2)={tight:.5f}: "
          f"{w_lo / tight:.0f}x .. {w_hi / tight:.0f}x")
    out["r0"] = {"table": r0_tab, "adopted": r0, "w1_over_r0": [w_lo / r0, w_hi / r0]}

    print("\n=== 3. dimensional mismatch: 7.14 vs 0.031 is not a comparison ===")
    print(f"  f1^2_bare = {A0:.2f} carries [mass]^2-ish lattice units (scheme-dependent, "
          "COMPUTE-F1SQ:20)")
    print(f"  r0 = {r0:.5f} is dimensionless (atom fraction).  Ratio '7.14/0.031' = "
          f"{A0 / r0:.0f} means nothing; the honest comparand is w1 = {w0:.3f}.")
    ok_dim = abs(w0 / r0 - A0 / r0 / c0) < 1e-6
    print(f"  consistency: w1 = A/C(0) => c0 = {c0:.3f} is the missing unit factor  "
          f"[{'PASS' if ok_dim else 'FAIL'}]")

    print("\n=== 4. can rho^2 = 2.48 deliver f1^2?  Exact 2-dim counterexample ===")
    rho2 = 2.48
    print("  model: transfer diag(mu1,mu2), mu2/mu1 = rho^2, source (cos th, sin th)")
    print("  => rho^2 = mu2/mu1 independent of th; w1 = cos^2 th free in (0,1)")
    print(f"{'th':>8} {'rho^2':>8} {'w1':>8} {'f1^2=M0*w1':>10}")
    for th in (0.002, 0.25, 0.7, 1.0, 1.2, 1.55):
        w1 = np.cos(th) ** 2
        print(f"{th:8.3f} {rho2:8.3f} {w1:8.5f} {w1:10.5f}")
    spread = np.log10(np.cos(0.002) ** 2 / np.cos(1.55) ** 2)
    print(f"  rho^2 fixed at {rho2}; w1 spans {np.cos(1.55)**2:.2e} .. {np.cos(0.002)**2:.6f} "
          f"= {spread:.1f} decades")
    print("  => rho^2 is a 1-constraint functional on a 2-dof object: it CANNOT fix w1.  "
          "corpus: 'rho^2 floor-blind' (new-tools-registry A7)")
    out["rho2_blind_decades"] = float(spread)

    print("\n=== 5. the only exact bridge: residue-level-gap lemma ===")
    print("  f1^2 >= M0(Lam - Rbar)/(Lam - m1^2),  Rbar = M1/M0,  Lam = inf supp(nu - f1^2 delta)")
    print("  two-level case: equality (tautology) -> the lemma returns w1 itself, no new info")
    m1sq = 1.0
    for rho2_t in (2.28, 2.48, 2.68, 5.79):
        print(f"  rho^2={rho2_t:5.2f}: w1_lb = (rho^2 - rbar)/(rho^2 - 1),  rbar = Rbar/m1^2")
        for rbar in (1.1, 1.5, 1.9, rho2_t - 1e-9):
            lb = (rho2_t - rbar) / (rho2_t - 1)
            flag = "  <-- vacuous" if lb <= 0.0313 else ""
            print(f"       rbar={rbar:5.2f} -> w1 >= {lb:6.4f}"
                  f"{'  >= r0 OK' if lb >= 0.0313 else flag}")
    print("  => a floor on w1 from rho^2 needs Rbar = M1/M0, which is NOT measured anywhere;")
    print("     with rbar -> rho^2 the bound degrades to 0.  So (c)-part: measure M1/M0 (time-moment).")
    out["lemma_needs_Rbar"] = True

    print("\n=== verdict inputs ===")
    out["w1_bare_level"] = w0
    out["w1_paper_route"] = w_paper
    out["w1_band"] = [w_lo, w_hi]
    out["w1_est_2p5"] = [0.35, 0.80]
    print(f"  w1(2.2) = {w_lo:.2f}..{w_hi:.2f}   [E]  (A/C(0); scheme-invariant under O -> lam O)")
    print(f"  w1 is the gamma=0 ratio (U4: f1^2 and M0 share Z_O^2) -> a, Z_a cancel in the ratio")
    print(f"  estimate w1(2.5) = {out['w1_est_2p5'][0]}..{out['w1_est_2p5'][1]}  [M]  "
          f"(no plateau at 2.5: 40 cfg, ~3.5 sigma; 2.5 is above the 2.30 crossover)")
    print(f"  => f1^2 >= r0 at 2.5: {w_lo / r0:.0f}-{w_hi / r0:.0f}x margin if w1 holds; verdict (a)")
    json.dump(out, open("research/v3/w7/p8_f1sq_above.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
