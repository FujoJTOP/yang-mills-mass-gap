# -*- coding: utf-8 -*-
"""
BC-BETHE / T0+T3 : exact Gaussian (free scalar) field on the Bethe lattice vs Z^d.

Model:  H = 1/2 sum_x phi_x^2 - beta sum_{<xy>} phi_x phi_y
  Z^d     : quadratic form 1 - 2 beta sum_{i=1..d} cos k_i ; min at k=0
            => m^2 = 1 - 4 d beta ,   beta_c = 1/(4d)
  tree(z) : adjacency spectral radius 2 sqrt(z-1) (z-regular tree, z = 2d neighbours)
            => m^2 = 1 - 4 beta sqrt(z-1) ,  beta_c = 1/(4 sqrt(z-1))
  m = sqrt(max(0,m^2)) : smallest excitation energy = mass gap.  m=0 <=> gapless.

Everything here is an exact formula (no truncation, no fit).  Purpose:
 (i)  fix the SIGN of the tree-vs-lattice mass comparison in a model that is
      exactly solvable on BOTH sides and is the free massless field at beta_c^lat;
 (ii) run the acid test: at beta = beta_c^lat the lattice is the free massless
      field (m_lat = 0) while the tree has m_tree > 0  =>  "m_tree > 0 => m_lat > 0"
      is FALSE, i.e. the tree can never certify a lattice gap.
"""
import json
import math

import numpy as np


def m_lat(d, beta):
    """Z^d free scalar mass gap (exact)."""
    return math.sqrt(max(0.0, 1.0 - 4.0 * d * beta))


def m_tree(d, beta):
    """(2d)-regular tree free scalar mass gap (exact)."""
    z = 2 * d
    return math.sqrt(max(0.0, 1.0 - 4.0 * beta * math.sqrt(z - 1.0)))


def bc_lat(d):
    return 1.0 / (4.0 * d)


def bc_tree(d):
    z = 2 * d
    return 1.0 / (4.0 * math.sqrt(z - 1.0))


def main():
    out = {"model": "H = 1/2 sum phi^2 - beta sum_edges phi phi  (free scalar / Gaussian)",
           "exact_formulae": {
               "Z^d": "m_lat = sqrt(max(0, 1 - 4 d beta)), beta_c = 1/(4d)",
               "tree(z=2d)": "m_tree = sqrt(max(0, 1 - 4 beta sqrt(z-1))), beta_c = 1/(4 sqrt(z-1))"}}

    print("=" * 78)
    print("BC-BETHE T0/T3  -- exact Gaussian/free-scalar comparison, tree(z=2d) vs Z^d")
    print("=" * 78)
    print("\n[1] critical couplings (exact)")
    print(f"{'d':>2} {'z':>3} {'beta_c^lat':>12} {'beta_c^tree':>12} {'ratio tr/lat':>13}")
    rates = {}
    for d in (2, 3, 4):
        bl, bt = bc_lat(d), bc_tree(d)
        rates[d] = {"beta_c_lat": bl, "beta_c_tree": bt, "ratio": bt / bl}
        print(f"{d:>2} {2*d:>3} {bl:>12.6f} {bt:>12.6f} {bt/bl:>13.4f}")
    out["critical_couplings"] = rates

    print("\n[2] mass ratio m_tree/m_lat at matched beta (tree ABOVE lattice everywhere")
    print("    in the disordered window: tree orders LAST in this model)")
    print(f"{'d':>2} {'beta':>8} {'m_lat':>10} {'m_tree':>10} {'tree/lat':>10}")
    rows = []
    for d in (2, 3, 4):
        bl = bc_lat(d)
        for frac in (0.25, 0.5, 0.75, 1.0, 1.25, 1.5):
            b = frac * bl
            ml, mt = m_lat(d, b), m_tree(d, b)
            r = (mt / ml) if ml > 0 else float("inf")
            rows.append({"d": d, "frac_of_beta_c_lat": frac, "beta": b,
                         "m_lat": ml, "m_tree": mt, "ratio_tree_over_lat": r})
            print(f"{d:>2} {b:>8.5f} {ml:>10.5f} {mt:>10.5f} {r:>10.4f}")
    out["mass_table"] = rows

    print("\n[3] ACID TEST at the free-massless-field point beta = beta_c^lat")
    print("    (Z^d is the free massless field: m_lat = 0 exactly)")
    print(f"{'d':>2} {'m_lat(beta_c^lat)':>18} {'m_tree(beta_c^lat)':>19} {'tree>0 & lat=0 ?':>18}")
    acid = []
    for d in (2, 3, 4):
        b = bc_lat(d)
        ml, mt = m_lat(d, b), m_tree(d, b)
        flag = (abs(ml) < 1e-15) and (mt > 0)
        acid.append({"d": d, "m_lat": ml, "m_tree": mt, "acid_counterexample": bool(flag)})
        print(f"{d:>2} {ml:>18.3e} {mt:>19.6f} {str(flag):>18}")
    out["acid_test"] = acid

    print("\n[4] reading")
    print("    * tree mass is strictly LARGER than the lattice mass at every matched beta;")
    print("    * at beta_c^lat the lattice is gapless (m=0) but the tree is massive")
    print("      => the implication 'm_tree>0 => m_lat>0' is FALSE in an exactly solvable model.")
    print("    * direction of the surviving strict inequality:  m_tree >= m_lat  (UPPER bound")
    print("      on the lattice mass out of a tree mass)  -> useless for a gap lower bound.")
    print("    * NOTE the contrast with the ferromagnetic Ising case (bc_ising.py), where the")
    print("      sign is OPPOSITE: the sign is not a geometric fact, it is FKG-dependent.")

    with open(__file__.replace(".py", ".json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
