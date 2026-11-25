# -*- coding: utf-8 -*-
"""
BC-BETHE / T2 : SU(2) on the Bethe lattice  vs  exact 2D SU(2) lattice gauge theory.

------------------------------------------------------------------------------
TREE SIDE (exact BP / recursive-distribution fixed point, no truncation):
  kernel on a link:   K(U,V) = exp(beta Re tr(U V^dag))   (class function of U V^dag)
  character expansion:  K = sum_j (2j+1) c_j chi_j ,
        c_j = (1/(2j+1)) int dU K(U) chi_j(U) = 2 I_{2j+1}(2 beta)/(2 beta)
  BP on a q-regular link tree (q-1 non-backtracking branches feeding a link):
        eta_{n+1}(U)  ~  [ (K eta_n)(U) ]^{q-1}
  trivial fixed point eta = const; write eta = c_0 + a chi_{1/2} + ... ; the
  linearised recursion multiplies a by (q-1) c_{1/2}/c_0, hence
        m_tree(beta) = -log( (q-1) c_{1/2}/c_0 )      (if the argument < 1)
                     = 0                              (ordered phase)
        beta_c^tree : (q-1) c_{1/2}/c_0 = 1 .
  q-1 = 7  <=> 8 neighbours per link, matched to Z^4 (2d = 8).
  q-1 = 3  <=> 4 neighbours per link, matched to Z^2 (2d = 4)  [control].
  (The 4-body plaquette tree is NOT this model; see the .md for why that
   distinction is itself part of the kill: "the tree gap" is not canonical.)

------------------------------------------------------------------------------
LATTICE SIDE (exact, 2D SU(2) Wilson LGT, temporal gauge U_t = 1):
  columns decouple; one column = chain of L_t links with kernel K.
  Transfer matrix eigenvalues on class functions = c_j (multiplicity (2j+1)^2),
  hence   m_lat(beta) = -log( c_{1/2}/c_0 ) = -log( I_2(2 beta)/I_1(2 beta) ).
  Exact for every beta; small-beta check:  m_lat = -log(beta/2) + O(beta^2).

No truncation, no fit, no Monte Carlo anywhere in this file.
"""
import json
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import iv


def cj(beta, j):
    """character coefficient of exp(beta Re tr U) on SU(2); exact (Bessel)."""
    z = 2.0 * beta
    return 2.0 * iv(int(2 * j + 1), z) / z


def cj_quad(beta, j):
    """independent quadrature check of cj (Haar: (2/pi) sin^2 theta dtheta)."""
    f = lambda th: (2.0 / math.pi) * math.sin(th) ** 2 * math.exp(2 * beta * math.cos(th)) \
        * math.sin((2 * j + 1) * th) / math.sin(th)
    val, _ = quad(f, 0.0, math.pi, limit=400)
    return val / (2 * j + 1)


def ratio(beta):
    return cj(beta, 0.5) / cj(beta, 0.0)          # = I_2(2b)/I_1(2b)


def m_lat(beta):
    return -math.log(ratio(beta))


def m_tree(beta, cav):
    x = cav * ratio(beta)
    return -math.log(x) if x < 1.0 else 0.0


def beta_c_tree(cav):
    f = lambda b: cav * ratio(b) - 1.0
    return brentq(f, 1e-6, 50.0, xtol=1e-14)


def main():
    out = {}

    print("=" * 78)
    print("BC-BETHE T2  --  SU(2): exact BP on the Bethe lattice vs exact 2D LGT")
    print("=" * 78)

    # ---- [0] verify the Bessel closed form against quadrature -----------------
    print("\n[0] check  c_j = 2 I_{2j+1}(2 beta)/(2 beta)  against Haar quadrature")
    err = 0.0
    for b in (0.3, 1.0, 2.5):
        for j in (0.0, 0.5, 1.0):
            e = abs(cj(b, j) - cj_quad(b, j)) / cj(b, j)
            err = max(err, e)
    print(f"    max relative error over beta in {{0.3,1,2.5}}, j in {{0,1/2,1}} : {err:.3e}")
    out["bessel_vs_quadrature_maxrelerr"] = err

    # ---- [1] critical couplings ---------------------------------------------
    print("\n[1] tree critical couplings (BP), and the 2D lattice (which has none: m_lat>0")
    print("    for every finite beta, m_lat -> 0 only as beta -> infinity)")
    bcs = {}
    for cav, tag in ((7.0, "q-1=7 (Z^4-matched)"), (3.0, "q-1=3 (Z^2-matched)")):
        b = beta_c_tree(cav)
        bcs[tag] = b
        print(f"    {tag:22s}  beta_c^tree = {b:.6f}   (small-beta law: beta_c ~ 2/(q-1) = {2.0/cav:.6f})")
    out["beta_c_tree"] = bcs

    # ---- [2] mass table ------------------------------------------------------
    print("\n[2] masses and ratio")
    print(f"{'beta':>7} {'m_lat(2D)':>11} {'m_tree(7)':>11} {'ratio7':>8} {'m_tree(3)':>11} {'ratio3':>8}")
    rows = []
    betas = [0.10, 0.20, 0.28, 0.30, 0.40, 0.60, 0.80, 1.00, 1.50, 2.00, 2.50, 3.00]
    for b in betas:
        ml = m_lat(b)
        mt7, mt3 = m_tree(b, 7.0), m_tree(b, 3.0)
        r7 = mt7 / ml
        r3 = mt3 / ml
        rows.append({"beta": b, "m_lat": ml, "m_tree_cav7": mt7, "ratio7": r7,
                     "m_tree_cav3": mt3, "ratio3": r3})
        print(f"{b:>7.2f} {ml:>11.6f} {mt7:>11.6f} {r7:>8.4f} {mt3:>11.6f} {r3:>8.4f}")
    out["mass_table"] = rows

    # ---- [3] the physical coupling beta = 2.5 --------------------------------
    print("\n[3] at the physical coupling beta = 2.5 (the point the charter cares about)")
    b = 2.5
    ml, mt7, mt3 = m_lat(b), m_tree(b, 7.0), m_tree(b, 3.0)
    print(f"    m_lat = {ml:.6f} ,  m_tree(cav=7) = {mt7:.6f}  (ordered: cav*ratio = {7*ratio(b):.4f} > 1)")
    print(f"    m_tree(cav=3) = {mt3:.6f}    =>  the tree bound m_lat >= m_tree is 0 : VACUOUS")
    print(f"    (the bound is non-trivial only for beta < beta_c^tree = {beta_c_tree(7.0):.4f}, i.e. deep")
    print("     strong coupling, where the tree value DIVERGES as beta -> 0 and is useless near 2.5)")
    out["beta_2p5"] = {"m_lat": ml, "m_tree_cav7": mt7, "m_tree_cav3": mt3,
                       "cav7_times_ratio": 7 * ratio(b)}

    # ---- [4] small-beta law --------------------------------------------------
    print("\n[4] small-beta check (strong coupling): m_lat + log(beta/2) -> 0 , m_tree - m_lat -> log(cav)")
    print(f"{'beta':>7} {'m_lat+log(b/2)':>16} {'m_tree7-m_lat':>15} {'log7':>8}")
    for b in (0.05, 0.10, 0.20):
        print(f"{b:>7.3f} {m_lat(b)+math.log(b/2):>16.6e} {m_tree(b,7.0)-m_lat(b):>15.6f} {math.log(7):>8.4f}")
    out["small_beta_shift"] = math.log(7.0)

    # ---- [5] the gap in the bound as a function of beta ----------------------
    print("\n[5] the tree bound  B(beta) := max(0, m_tree(beta))  vs beta")
    print("    B > 0 only for beta < beta_c^tree;  B = 0 for every beta above it,")
    print("    in particular for the whole scaling window beta in [2.0, 3.0].")
    out["bound_support"] = {"cav7": bcs["q-1=7 (Z^4-matched)"], "cav3": bcs["q-1=3 (Z^2-matched)"]}

    print("\n[6] reading")
    print("    (i)  in the window where the tree is disordered, m_tree < m_lat ALWAYS, with the")
    print("         beta-independent shift  m_lat - m_tree = log(cav)  (here log 7 = 1.9459);")
    print("    (ii) the shift is exactly the connectivity multiplicity, i.e. the tree is again")
    print("         the 'worst case' (fast counting, slow decay, small mass);")
    print("    (iii) at beta_c^tree the tree's own gap closes; above it the tree bound is 0,")
    print("         so at beta = 2.5 the comparison yields  m_lat >= 0  : vacuous.")

    with open(__file__.replace(".py", ".json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
