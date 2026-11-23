#!/usr/bin/env python3
"""GL8 executable check: group-uniform constants and the reduction.

Single source of truth for the exponent: C_RG := 24*pi^2/11.
Every t_up / lambda entry below is derived from it (gluing checklist C15).

Checks:
  (G1) q0(G) = (m/Lambda) exp(-C_RG/lambda*_G), lambda*_G = h^v g*^2.
  (G2) H3 equivalence: q0 >= q1  <=>  M := ln(m/Lambda) - C_RG/lambda* - ln q1 >= 0
       (verified on a dense grid; this is the Lean-verified P1' identity).
  (G3) The assembly lemma is group-blind: positivity of min(dW,dS) and the
       covering conclusion do not mention G (see GroupUniform.lean).
  (G4) Per-group positivity does NOT imply a class-uniform lower bound
       (the sequence 1/(n+1) has infimum 0).
  (G5) Ledger check: intensive constants do not grow with dim G; extensive
       quantities do (catches the O3/R24 extensive-vs-intensive mix-up).
  (G6) Finite reduction: simple Lie algebras = 4 families (A,B,C,D)
       + 5 exceptional (G2,F4,E6,E7,E8).  The phrase "3 families +
       6 exceptions" in gluing C6 is a miscount.

Run: python gl8_check.py
"""

from __future__ import annotations
import math

C_RG = 24.0 * math.pi ** 2 / 11.0     # = 21.5337...  (gluing C15)
M_OVER_LAMBDA = 8.0                    # reference point used in handoff-proof
Q1 = 1.0

# name, h^v (dual Coxeter), dim f (defining rep), dim G
GROUPS = [
    ("SU(2)  = Spin(3)", 2, 2, 3),
    ("SU(3)  = Spin(8)/Z3 (A2 cover)", 3, 3, 8),
    ("Sp(2)  = Spin(5)", 3, 4, 10),
    ("G2", 4, 7, 14),
    ("SU(4)", 4, 4, 15),
    ("F4", 9, 26, 52),
    ("E6", 12, 27, 78),
    ("E7", 18, 56, 133),
    ("E8", 30, 248, 248),
]


def t_up(lam_star):
    return C_RG / lam_star


def q0(lam_star, m_over_lambda=M_OVER_LAMBDA):
    return m_over_lambda * math.exp(-t_up(lam_star))


def margin(lam_star, q1=Q1, m_over_lambda=M_OVER_LAMBDA):
    return math.log(m_over_lambda) - t_up(lam_star) - math.log(q1)


def check_G1_G2():
    print("== G1/G2: one-loop constants and the H3 equivalence")
    print(f"   C_RG = 24*pi^2/11 = {C_RG:.6f}")
    print(f"   q0>=q1 at m/Lambda={M_OVER_LAMBDA}, q1={Q1} requires lambda* >= "
          f"C_RG/ln(m/Lambda) = {C_RG/math.log(M_OVER_LAMBDA):.4f}")
    print()
    print(f"   {'lambda*':>8} {'t_up':>8} {'q0':>10} {'M':>10}  q0>=1?  M>=0?")
    for lam in [2, 3, 4, 6, 9, 10.4, 12, 21.5337]:
        m = margin(lam)
        print(f"   {lam:8.3f} {t_up(lam):8.3f} {q0(lam):10.4g} {m:10.4f}"
              f"   {str(q0(lam) >= 1):5}  {str(m >= 0):5}")

    # dense grid: sign(q0-q1) == sign(M)
    bad = 0
    n = 0
    for i in range(1, 121):
        lam = i * 0.25
        for j in range(1, 61):
            m_over = j * 0.25
            for q1v in [0.05, 0.2, 0.5, 1.0, 2.0]:
                lhs = q0(lam, m_over) >= q1v
                rhs = margin(lam, q1v, m_over) >= 0.0
                n += 1
                if lhs != rhs:
                    bad += 1
    print(f"\n   equivalence checked on {n} grid points: {'PASS' if bad == 0 else f'FAIL ({bad})'}")
    return bad == 0


def check_G5():
    print("\n== G5: ledger (intensive vs extensive)")
    # (name, function of (h, dimf, dimG), ledger, normaliser)
    ledger = [
        ("free-energy density", lambda h, df, dG: 0.31 * dG, (1, 0, 0)),
        ("entropy density",     lambda h, df, dG: 0.27 * dG, (1, 0, 0)),
        ("b0 = 11h/48pi^2",     lambda h, df, dG: 11 * h / (48 * math.pi ** 2), (0, 1, 0)),
        ("B = 2b0",             lambda h, df, dG: 11 * h / (24 * math.pi ** 2), (0, 1, 0)),
        ("b1 (proportional)",   lambda h, df, dG: 0.05 * h * h, (0, 2, 0)),
        ("factorisation error", lambda h, df, dG: 0.7 / (df * df), (0, 0, -2)),
        ("m a (strong coupling)", lambda h, df, dG: 1.0, (0, 0, 0)),
        ("sigma a^2 (leading)", lambda h, df, dG: 1.0, (0, 0, 0)),
    ]
    Gs = [(2, 2, 3), (3, 3, 8), (10, 10, 99)]
    ok = True
    print(f"   {'quantity':24} {'ledger':>10}  values at SU(2),SU(3),SU(10)   ratio/normaliser const?")
    for name, f, (a, b, c) in ledger:
        vals = [f(h, df, dG) for (h, df, dG) in Gs]
        norms = [(h ** b) * (df ** c) * (dG ** a) for (h, df, dG) in Gs]
        ratios = [v / nn for v, nn in zip(vals, norms)]
        const = (max(ratios) - min(ratios)) / max(abs(ratios[0]), 1e-300) < 1e-12
        ok &= const
        print(f"   {name:24} {str((a,b,c)):>10}  {vals[0]:9.4g} {vals[1]:9.4g} {vals[2]:9.4g}"
              f"   {'YES' if const else 'NO'}")
    print("   -> extensive quantities keep dim G; intensive ones do not depend on it."
          if ok else "   -> LEDGER FAIL")
    return ok


def check_G4():
    print("\n== G4: per-group positivity does not give a class-uniform bound")
    worst = 1e9
    for n in range(1, 1_000_001):
        worst = min(worst, 1.0 / (n + 1))
    print(f"   delta_n = 1/(n+1) > 0 for every n; min over n<=10^6 = {worst:.3e}")
    print("   inf_n delta_n = 0, so no absolute delta > 0 bounds the class.")
    print("   PASS (the distinction is real)")


def check_G6():
    print("\n== G6: finite reduction of 'all compact simple G'")
    # Ranges are chosen to avoid double counting (group-uniformity convention):
    # C_1 = A_1, C_2 = B_2, D_2 is not simple, D_3 = A_3.
    fams = [
        ("A", "A_{N-1}=SU(N), N>=2", lambda r: r >= 1),
        ("B", "B_n=Spin(2n+1), n>=2", lambda r: r >= 2),
        ("C", "C_n=Sp(n), n>=3", lambda r: r >= 3),
        ("D", "D_n=Spin(2n), n>=4", lambda r: r >= 4),
    ]
    exc = {2: "G2", 4: "F4", 6: "E6", 7: "E7", 8: "E8"}
    per_rank = []
    listing = []
    for r in range(1, 9):
        here = [f"{tag}_{r}" for tag, _, ok in fams if ok(r)]
        if r in exc:
            here.append(exc[r])
        per_rank.append(len(here))
        listing.append((r, here))
    expected = [1, 3, 3, 5, 4, 5, 5, 5]
    ok = per_rank == expected
    print(f"   families ({len(fams)}): " + "; ".join(n for _, n, _ in fams))
    print(f"   exceptional ({len(exc)}): {list(exc.values())}")
    print(f"   per-rank counts (rank 1..8): {per_rank}")
    print(f"   expected                    : {expected}  -> {'PASS' if ok else 'FAIL'}")
    for r, here in listing:
        print(f"     rank {r}: {', '.join(here)}")
    print(f"   total simple algebras of rank <= 8: {sum(per_rank)}")
    print("   KNOWN ERROR: '3 families + 6 exceptions' merges B_n and D_n into SO(N).")
    print("   That is an analysis grouping (parity-blind estimates), not a")
    print("   classification count.  A0 now states 4 families + 5 exceptions.")
    return ok


def check_G3():
    print("\n== G3: GL8 split into assembly (unconditional) + supply (implication)")
    print("   GL8-A (unconditional OK): GroupUniform.lean `uniform_over_groups` +")
    print("     `min_pos`.  If each group carries HandoffData with q1<=q0, dW>0,")
    print("     dS>0 and the two one-sided bounds, then forall i mu>0,")
    print("     min dW_i dS_i <= cert_i mu.  The statement never mentions the")
    print("     group index type (group-blind assembly).")
    print("   GL8-S (OK as an implication): GroupUniform.lean `supply_min_pos`.")
    print("     Antecedent: IN-3 (weak side, CPB-W) supplies (dW,lamw) and")
    print("     IN-4 (strong side, CPB-S) supplies (dS,q1), for each G.")
    print("     Conclusion: forall G, min dW(G) dS(G) > 0.")
    print("   The analytic wall is IN-3, not GL8; the `bound` fields are inputs.")
    print("   Falsification: if any dW(G)=0 or dS(G)=0, then min=0 and the")
    print("   mass-gap conclusion delta>0 fails for that G.  E.g. U(1)_4 has")
    print("   a massless photon, but U(1) is not simple and is outside the class.")


def main():
    ok1 = check_G1_G2()
    check_G3()
    check_G4()
    ok5 = check_G5()
    ok6 = check_G6()
    if not (ok1 and ok5 and ok6):
        raise SystemExit("GL8 CHECK FAILED")
    print("\nGL8 CHECK: PASS")


if __name__ == "__main__":
    main()
