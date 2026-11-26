#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
h09_h1_falsify.py  --  H09 / falsify "H1 => m > 0".

Target (H1, verbatim ch02-squeeze.tex:87-92 / Q6-H1-STATUS.md:76):
    H1 : chi := Var_mu(xi_{k->0}) = int_0^inf s^{-1} dnu(s) <= 1/lambda ,
         uniform in (L,a).                                    [nu = KL measure of 1/4 tr F^2]
Bare H1 (ch18-scalar-h1.tex:57-58): G''(0) > 0 ; family H1(d) <=> G''(0) >= d.
Physical member (app.tex:190-191, P03-EQUIV:253-258):
    H1 <=> Delta_dyn >= a^2 Delta_0^2 ,  delta_curv := a^2 Delta_0^2 + W_full ,
    with the split identity G'' = Delta_dyn + W_geo + X , W_full := W_geo + X   [ch18:28-31].
    Phi = 0 (Gaussian) ==> W_geo = X = 0 ==> W_full = 0 ==> Delta_dyn = G''(0).
Mass: m_1^2 := inf supp nu ,  m := sqrt(m_1^2) ;  f_1^2 := nu({m_1^2}) ;  M_0 := nu([0,inf)).
Gaussian marginal identity (Q6:83-86): G''(k=0) = 1/chi  (exactly, iff nu single atom; here we
use it as the *definition* of the one-coordinate curvature, i.e. Var = chi => G''=1/Var).

Everything below is exact rational arithmetic (fractions.Fraction) except the quoted
closed forms checked by quadrature.  No model is fitted.  Output is labelled [T]/[E]/[V].
"""
from fractions import Fraction as F
import math

BAR = "=" * 78
def hdr(t): print("\n" + BAR + "\n" + t + "\n" + BAR)

# ---------------------------------------------------------------- helpers
def chi_of(points):
    """exact chi = sum w_i / s_i  for an atomic measure given as [(w,s),...]"""
    return sum(w / s for w, s in points)

def weight_below(points, eps):
    return sum(w for w, s in points if s <= eps)

# ---------------------------------------------------------------- [E1] thin atom, Q02 form
def e1_thin_atom_Q02(lam=F(1), eps_list=(F(1,100), F(1,10**4), F(1,10**6), F(1,10**9))):
    hdr("[E1] thin atom, Q02/paper form:  w = eps/(2*lam),  s2 = 2*lam*(1-w)")
    print(f"{'eps':>12} {'w=f1^2':>14} {'s2':>18} {'chi':>10} {'1/lam':>8} {'H1':>6} "
          f"{'G\"=1/chi':>12} {'Delta_dyn':>12} {'m1^2':>10} {'f1^2/chi<=m1^2<=M0/chi':>24}")
    rows = []
    for eps in eps_list:
        w = eps / (2 * lam)
        s2 = 2 * lam * (1 - w)
        pts = [(w, eps), (1 - w, s2)]
        chi = chi_of(pts)
        G2 = 1 / chi                      # Phi=0 => W_full=0 => Delta_dyn = G''
        Ddyn = G2 - 0
        m1sq = eps
        h1 = chi <= 1 / lam
        lo = w / chi <= m1sq              # ch02 squeeze lower rung
        hi = m1sq <= 1 / chi              # M0 = 1 upper rung
        print(f"{float(eps):12.3e} {float(w):14.6e} {float(s2):18.12f} {float(chi):10.6f} "
              f"{float(1/lam):8.4f} {str(h1):>6} {float(G2):12.6f} {float(Ddyn):12.6f} "
              f"{float(m1sq):10.3e} {str(lo and hi):>24}")
        rows.append(dict(eps=eps, w=w, s2=s2, chi=chi, G2=G2, h1=h1, lo=lo, hi=hi))
    print("exact chi for eps=1e-9 :", rows[-1]['chi'], " == 1/lam ?", rows[-1]['chi'] == 1 / lam)
    print("=> [T] chi == 1/lam EXACTLY at every eps; H1 holds with EQUALITY (delta_max = G'' = lam).")
    print("=> [T] m1^2 = eps -> 0 while H1 true:  no implication (measure level), residue f1^2 = w ~ eps.")
    return rows

# ---------------------------------------------------------------- [E2] strict-margin version
def e2_thin_atom_strict(lam=F(1), eps_list=(F(1,100), F(1,10**4), F(1,10**6), F(1,10**9))):
    hdr("[E2] NEW: thin atom with a STRICT H1 margin:  w = eps/(4*lam),  s2 = 2*lam*(1-w)")
    print("(removes the 'boundary case / equality' objection to E1)")
    print(f"{'eps':>12} {'w=f1^2':>14} {'s2':>20} {'chi':>12} {'1/lam':>8} {'margin':>10} "
          f"{'G\"':>10} {'Delta_dyn':>12} {'m1^2':>10}")
    out = []
    for eps in eps_list:
        w = eps / (4 * lam)
        s2 = 2 * lam * (1 - w)
        pts = [(w, eps), (1 - w, s2)]
        chi = chi_of(pts)                      # = 1/(4lam) + 1/(2lam) = 3/(4lam)
        G2 = 1 / chi                           # = 4lam/3
        marg = 1 / lam - chi
        print(f"{float(eps):12.3e} {float(w):14.6e} {float(s2):20.12f} {float(chi):12.8f} "
              f"{float(1/lam):8.4f} {float(marg):10.6f} {float(G2):10.6f} {float(G2):12.6f} "
              f"{float(eps):10.3e}")
        assert chi == 3 * lam / 4, chi
        assert G2 == 4 * lam / 3
        assert chi < 1 / lam
        out.append(dict(eps=eps, w=w, s2=s2, chi=chi, G2=G2, marg=marg))
    print("=> [T] chi = 3/(4lam) < 1/lam  (H1 strict, margin lam/4),  G'' = 4lam/3 >= lam,")
    print("      Delta_dyn = G'' = 4lam/3 > 0  (Phi=0),  m1^2 = eps -> 0,  f1^2 = eps/(4lam).")
    print("=> so H1(lam) strictly true AND the physical member Delta_dyn >= a^2Delta_0^2 (W_full=0)")
    print("   holds strictly, yet m = 0.  This is a counterexample to BOTH forms of H1.")
    return out

# ---------------------------------------------------------------- [E3] near-degenerate pair
def e3_near_degenerate(lam=F(1), eps_list=(F(1,10), F(1,100), F(1,1000), F(1,10**4))):
    hdr("[E3] candidate (2): two NEARLY DEGENERATE atoms, O(1) weights -- does it satisfy H1?")
    print(" nu = p*delta_eps + (1-p)*delta_(eps+delta),  p = 1/2,  delta = eps  (mass ratio 1:2)")
    print(f"{'eps':>10} {'p':>6} {'chi_exact':>16} {'1/lam':>8} {'H1':>6} {'lb 1/(2eps)':>14}")
    bad = 0
    for eps in eps_list:
        pts = [(F(1, 2), eps), (F(1, 2), 2 * eps)]
        chi = chi_of(pts)                      # = 1/(2eps) + 1/(4eps) = 3/(4eps)
        h1 = chi <= 1 / lam
        lb = 1 / (2 * eps)
        print(f"{float(eps):10.3e} {0.5:6.2f} {str(chi):>16} {float(1/lam):8.4f} {str(h1):>6} "
              f"{float(lb):14.3e}")
        if not h1: bad += 1
    print(f"=> H1 false on {bad}/{len(eps_list)} points.")
    print("=> [T] general bound:  nu([0,2eps]) >= c  ==>  chi >= c/(2eps).  Near-degeneracy at the")
    print("   bottom with O(1) residue is FORBIDDEN by H1.  So candidate (2) does NOT kill the target:")
    print("   it dies on H1 itself, not merely on m>0.")

# ---------------------------------------------------------------- [E4] continuous edge
def e4_continuous_edge(lam=F(1)):
    hdr("[E4] candidate (3): CONTINUOUS spectrum, no atom at all:  dnu = s^a e^{-s}/Gamma(a+1) ds")
    print("closed form: M0 = 1,  chi = Gamma(a)/Gamma(a+1) = 1/a,  G'' = a,  m1^2 = inf supp = 0,")
    print("             f1^2 = 0 (no atom),  nu([0,eps]) = int_0^eps s^a e^{-s}/Gamma(a+1) ds.")
    try:
        from scipy.integrate import quad
        quad_ok = True
    except Exception as e:                                   # pragma: no cover
        quad_ok = False
        print("  [scipy unavailable:", e, "]  -> closed forms only")
    print(f"{'a':>6} {'chi=1/a':>10} {'quad chi':>12} {'1/lam':>8} {'H1(<=)':>8} {'strict':>8} "
          f"{'G\"=a':>8} {'m1^2':>6} {'f1^2':>6} {'nu([0,1e-3])':>14}")
    for a in (F(1, 4), F(1, 2), F(1), F(2), F(4)):
        af = float(a)
        chi = 1 / a
        if quad_ok:
            G = math.gamma(af + 1)
            q, _ = quad(lambda s: s ** (af - 1) * math.exp(-s) / G, 0, math.inf, limit=200)
        else:
            q = float('nan')
        mass_small, _ = (quad(lambda s: s ** af * math.exp(-s) / math.gamma(af + 1),
                              0, 1e-3, limit=200) if quad_ok else (float('nan'), 0))
        print(f"{af:6.2f} {float(chi):10.6f} {q:12.9f} {float(1/lam):8.4f} "
              f"{str(chi <= 1/lam):>8} {str(chi < 1/lam):>8} {af:8.2f} {0:6.1f} {0:6.1f} "
              f"{mass_small:14.3e}")
    print("=> [T] for every a >= lam:  H1 holds STRICTLY (chi = 1/a <= 1/lam), G'' = a >= lam,")
    print("   Delta_dyn = G'' > 0 (Phi=0), and m1^2 = 0 EXACTLY with NO atom (f1^2 = 0).")
    print("=> the mass in [0,eps] scales like eps^(a+1): the edge residue vanishes; the bottom is")
    print("   at 0 for every eps>0 -- H1 never sees it.  This is not a limit: m = 0 identically.")

# ---------------------------------------------------------------- [E5] for EVERY lambda
def e5_for_every_lambda():
    hdr("[E5] 'for every lambda':  pick a = 2*lam  ==>  chi = 1/(2lam) <= 1/lam (margin 1/2), m1^2 = 0")
    print(f"{'lam':>8} {'a=2lam':>8} {'chi=1/a':>14} {'1/lam':>10} {'chi/(1/lam)':>12} "
          f"{'G\"=a':>10} {'G\"/lam':>10} {'m1^2':>6}")
    for lam in (F(1, 10), F(1), F(10), F(1000)):
        a = 2 * lam
        chi = 1 / a
        print(f"{float(lam):8.1f} {float(a):8.1f} {float(chi):14.6e} {float(1/lam):10.6e} "
              f"{float(chi*lam):12.6f} {float(a):10.1f} {float(a/lam):10.2f} {0:6.1f}")
    print("=> [T] no choice of the floor lambda makes H1 imply m>0: the witness exists for each lam.")

# ---------------------------------------------------------------- [E6] the one-line mechanism
def e6_mechanism():
    hdr("[E6] the mechanism theorem (this is the ONE link that breaks)")
    print("For any positive measure nu and any eps>0:")
    print("    chi = int s^{-1} dnu >= (2eps)^{-1} nu([0,2eps])   ==>   nu([0,2eps]) <= 2 eps chi.")
    print("So H1 (chi <= 1/lam) forces     nu([0,2eps]) <= 2eps/lam   for every eps:  the edge")
    print("residue must vanish at rate O(eps).  H1 therefore constrains the RESIDUE DENSITY at the")
    print("bottom, never its LOCATION.  Any claim 'H1 => m>0' must convert a residue ceiling into a")
    print("position floor; those are different coordinates on the same cone of measures.  [T]")
    print("Numerical check (weight in [0,2eps] vs the H1 ceiling 2eps/lam, lam=1):")
    eps = F(1, 1000)
    lam = F(1)
    cases = {
        "E1 thin atom (w=eps/2lam)": [(eps / (2 * lam), eps), (1 - eps / (2 * lam), 2 * lam * (1 - eps / (2 * lam)))],
        "E2 thin atom (w=eps/4lam)": [(eps / (4 * lam), eps), (1 - eps / (4 * lam), 2 * lam * (1 - eps / (4 * lam)))],
        "E3 near-degenerate 1/2,1/2": [(F(1, 2), eps), (F(1, 2), 2 * eps)],
    }
    for name, pts in cases.items():
        chi = chi_of(pts)
        print(f"  {name:32s} chi={float(chi):10.4f}  nu([0,2eps])={float(weight_below(pts,2*eps)):9.6f}"
              f"  ceiling 2eps/lam={float(2*eps/lam):9.6f}  H1={chi <= 1/lam}")

# ---------------------------------------------------------------- main
if __name__ == "__main__":
    print("H09 -- falsify 'H1 => m>0' : exact rational witnesses")
    print("H1 : chi = int s^{-1} dnu <= 1/lambda (uniform in (L,a));  m1^2 = inf supp nu;  M0 = 1")
    print("Phi=0 (Gaussian) => W_full = 0 => Delta_dyn = G''(0) = 1/chi   [ch18:28-31 + Q6:83-86]")
    e1_thin_atom_Q02()
    e2_thin_atom_strict()
    e3_near_degenerate()
    e4_continuous_edge()
    e5_for_every_lambda()
    e6_mechanism()
    hdr("VERDICT")
    print("(a) counterexample FOUND.  Three witnesses, one mechanism:")
    print("    E2 (thin atom, strict margin) and E4 (pure continuum edge, no atom) both satisfy H1")
    print("    in the lambda form AND the physical member (Delta_dyn = G'' > 0 because W_full = 0),")
    print("    with m1^2 -> 0 / m1^2 = 0.  E3 (near-degenerate pair) FAILS H1 -- it is not a witness.")
    print("    => 'H1 => m>0' is FALSE as a predicate implication, for every lambda (E5).")
    print("    Scope: witnesses are KL-measure / Gaussian (generalized-free-field) models; that is the")
    print("    level at which H1 is stated (ch02 eq:h1) and machine-checked (H1Sharpening.lean).")
