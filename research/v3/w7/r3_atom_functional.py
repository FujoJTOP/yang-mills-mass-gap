#!/usr/bin/env python3
r"""
r3_atom_functional.py -- R3: can a weak-* DISCONTINUOUS functional READ the edge atom
        f1^2 = nu({m1^2})  (m1^2 = inf supp nu)  without already knowing the gap?

Build the discontinuous atom-reading functionals named in the task and TEST each on solvable /
adversary models, to decide whether any is (a) discontinuous (acid-capable), (b) computable from
the PHYSICAL measure, and (c) NOT circular.

Candidates:
  [C1] Stieltjes boundary inversion.  Correct form is the DOUBLE limit
         f1^2 = lim_{d->0} lim_{e->0} (1/pi) Int_{m-d}^{m+d} Im Pi(u+ie) du .
       (The single-window form with d=e, as literally printed in ch09 eq.(atomweight), returns
        w/2 -- a factor-2 slip; we record and fix it here.)   Needs: location m AND d,e->0.
  [C2] "sup-over-x" atom functional  A(nu)=sup_x nu({x})  (location-free, but NOT the EDGE atom).
  [C3] Finite-data blindness + Stone-Weierstrass resolution bound: atom vs atomless bump with
       identical low moments, f1^2 jumps O(1); sup|E-g| >= w/2 for any continuous (finite-data) g.
  [C4] Thin atom: chi=1/L for ALL e, yet edge atom -> 0  (no finite-data floor).
  [C5] Hankel / Christoffel function  lambda_n(x)=1/(v^T G_n^{-1} v), G_n the moment-Hankel matrix.

Solvable-model acid tests (tool MUST return f1^2=0):
  [S1] SU(2) one-plaquette Haar marginal rho1(z) (exact closed form) -- SMOOTH, atomless.
  [S2] free massless cut proxy rho(s)=const on [0,1]                -- atomless, f1^2=0.

numpy only, one job, -I safe.  Labels [T]/[C]/[NG]/[O]/[V].  No QED.  Clay 0%.
"""
import numpy as np
from scipy.special import sici

np.seterr(all='ignore')
EULER = 0.5772156649015328606
TWO_PI = 2.0 * np.pi

# =====================================================================  measure helpers
# measure = (atoms, dens):  atoms = list[(w, s)],  dens = (grid s[], rho[]) or None

def mass_in_window(atoms, dens, x, delta, eps=1e-6):
    """(1/pi) Int_{x-delta}^{x+delta} Im Pi(u+i eps) du   ->  nu([x-delta,x+delta]) as eps->0.
       Atom term exact: (1/pi)[arctan((x+d-s)/e) - arctan((x-d-s)/e)] -> 1_{|s-x|<d}.
       Density term: same kernel applied to rho(s) (stable; rho smooth)."""
    def A(s):
        return (1.0 / np.pi) * (np.arctan((x + delta - s) / eps) - np.arctan((x - delta - s) / eps))
    tot = 0.0
    if atoms:
        for w, s in atoms:
            tot += w * A(s)
    if dens is not None:
        g, rho = dens
        tot += np.trapezoid(rho * A(g), g)
    return tot

def moments_atoms_dens(atoms, dens, K):
    m = np.zeros(K)
    if atoms:
        for w, s in atoms:
            for k in range(K):
                m[k] += w * s ** k
    if dens is not None:
        g, rho = dens
        for k in range(K):
            m[k] += np.trapezoid(rho * g ** k, g)
    return m

def unif_moment(x, h, k):
    """k-th moment of Unif[x-h, x+h], mass 1:  (1/2h) Int_{-h}^{h}(x+t)^k dt."""
    from math import comb
    tot = 0.0
    for j in range(0, k + 1, 2):          # only even j survive the symmetric average
        tot += comb(k, j) * x ** (k - j) * h ** j / (j + 1)
    return tot

# =====================================================================  exact SU(2) Haar marginal
def Cin(x):
    x = np.asarray(x, float); _, ci = sici(x)
    out = EULER + np.log(np.where(x > 0, x, 1.0)) - ci
    return np.where(x > 0, out, 0.0)

def rho1_haar(z):
    z = np.asarray(z, float); az = np.abs(z)
    return np.where(az <= TWO_PI, np.pi * (Cin(TWO_PI) - Cin(np.minimum(az, TWO_PI))), 0.0)

def hdr(t):
    print("\n" + "=" * 100 + "\n" + t + "\n" + "=" * 100)

# =====================================================================  main
def main():
    print("#" * 100)
    print("# R3 -- discontinuous atom-reading functionals: build, test, decide circularity")
    print("#" * 100)

    # ------------------------------------------------------------------ [C1] Stieltjes boundary
    hdr("[C1] STIELTJES BOUNDARY INVERSION  f1^2 = lim_{d->0}lim_{e->0}(1/pi)Int_{m-d}^{m+d}ImPi(u+ie)du")
    atoms = [(0.30, 1.0), (0.70, 3.0)]
    print("  measure: 0.30*delta_1 + 0.70*delta_3   (edge atom at m=1, weight w=0.30)")
    print("  mass-in-window M(d) at edge x=1, vs a non-atom point x=2:")
    print("  %10s %18s %18s" % ("d", "M(d)@x=1(atom)", "M(d)@x=2(no atom)"))
    for d in (0.5, 0.1, 0.03, 0.01, 0.003):
        print("  %10.4g %18.6f %18.6f" % (d, mass_in_window(atoms, None, 1.0, d),
                                                mass_in_window(atoms, None, 2.0, d)))
    print("  -> M(d)@x=1 -> w=0.30 (atom);  M(d)@x=2 -> 0.  The atom read is EXACT. [T]")
    print("  CIRCULARITY: evaluation sets x = m1^2 = 1  ->  presupposes the gap. (R1)")
    print("  *** NOTE (correction to ch09 eq.atomweight): with d=e (window=kernel width) the INTEGRAL")
    print("      is w/2, not w -- the correct object is the double limit d->0 AFTER e->0. ***")

    print("\n  ATOMLESS BUMP: total mass 0.30, half-width h=0.04 centred at 1 (same mass & place):")
    h = 0.04; gb = np.linspace(0.5, 1.5, 200001); rb = np.where(np.abs(gb - 1.0) <= h, 0.30 / (2 * h), 0.0)
    print("  %10s %22s" % ("d", "M(d)@x=1 (bump)"))
    for d in (0.5, 0.1, 0.04, 0.02, 0.005):
        print("  %10.4g %22.6f" % (d, mass_in_window([], (gb, rb), 1.0, d)))
    print("  -> for d < h the window drops below 0.30; f1^2(bump)=0, f1^2(atom)=0.30 : the jump is REAL")
    print("     but only visible once resolution d < h.  A finite-resolution observer d_min > h sees")
    print("     the SAME M(d_min) for atom and bump -> cannot certify the atom. = A3/infinite. (R2)")

    # ------------------------------------------------------------------ [C2] sup-over-x atom functional
    hdr("[C2] 'SUP-OVER-x' ATOM FUNCTIONAL  A(nu)=sup_x nu({x})  (location-free, WRONG atom)")
    print("  thin atom: nu_e=(e/2)delta_e+(1-e/2)delta_2 ;  EDGE atom w_e=e/2 ,  UPPER atom ~1")
    print("  %10s %14s %16s %18s" % ("eps", "edge m=e", "edge w=e/2", "sup_x nu({x})"))
    for e in (0.5, 0.2, 0.05, 0.01, 0.001):
        print("  %10.4g %14.4g %16.6g %18.6f" % (e, e, e / 2, max(e / 2, 1 - e / 2)))
    print("  -> sup_x ~= 1 (upper atom); EDGE atom f1^2=e/2 -> 0. Targets the WRONG atom. [T]")
    print("     Also >0 on any FINITE-volume free field (discrete spectrum) -> FAILS the acid test. [NG]")

    # ------------------------------------------------------------------ [C3] finite-data + SW bound
    hdr("[C3] FINITE-DATA BLINDNESS + STONE-WEIERSTRASS RESOLUTION BOUND")
    A, x1, B, x2 = 0.45, 1.0, 0.55, 3.3
    mA = moments_atoms_dens([(A, x1), (B, x2)], None, 5)
    print("  nu_A=%.2f delta_%.1f + %.2f delta_%.1f (edge atom f1^2=%.2f); nu_B = same atoms -> Unif"
          % (A, x1, B, x2, A))
    print("  %10s %14s %14s %14s %14s %14s" % ("h", "|dM0|", "|dM1|", "|dM2|", "|dM3|", "|dM4|"))
    for h in (0.20, 0.10, 0.05, 0.02, 0.01, 0.001):
        mB = np.array([A * unif_moment(x1, h, k) + B * unif_moment(x2, h, k) for k in range(5)])
        d = np.abs(mA - mB)
        print("  %10.4g %14.3e %14.3e %14.3e %14.3e %14.3e" % (h, d[0], d[1], d[2], d[3], d[4]))
    print("  -> dM0=dM1=0 EXACTLY; dM2=(A+B)h^2/3, dM3=(Ax1+Bx2)h^2, ... all -> 0; yet f1^2: %.2f->0. [T]" % A)
    print("  STONE-WEITERSTRASS: any weak-* CONTINUOUS g (= any finite-data functional) has")
    print("     sup_nu |E(nu)-g(nu)| >= w/2 = %.3f  (a continuous instrument cannot resolve a jump" % (A / 2))
    print("     of size w better than w/2).  resolution >= signal  -> circular. [T]")

    # ------------------------------------------------------------------ [C4] thin atom
    hdr("[C4] THIN ATOM -- chi=1/L for ALL e, edge atom -> 0  (no finite-data floor)")
    Lam = 1.0
    print("  nu_e=(e/2L)delta_e+(1-e/2L)delta_{s2}, s2=2L(1-w), w=e/2L")
    print("  %10s %12s %16s %16s %14s" % ("eps", "w=e/2L", "chi=int s^-1 dnu", "edge mass w", "edge m"))
    for e in (0.5, 0.2, 0.05, 0.01, 0.001):
        w = e / (2 * Lam); s2 = 2 * Lam * (1 - w); chi = w / e + (1 - w) / s2
        print("  %10.4g %12.6f %16.6f %16.6f %14.4g" % (e, w, chi, w, e))
    print("  -> chi=1/L EXACTLY for every e [T]; edge atom w->0, edge m->0. Any functional of chi")
    print("     returns the SAME value for all e and cannot lower-bound w. [NG]")

    # ------------------------------------------------------------------ [C5] Hankel / Christoffel
    hdr("[C5] HANKEL / CHRISTOFFEL  lambda_n(x)=1/(v^T G_n^{-1} v),  G_n=(M_{j+k}), v=(1,x,..,x^n)")
    # support rescaled to [0,1] for conditioning:  t=s/3 ; atom at 1/3 (w .30), continuum [0.4,1] (w .70)
    def M_k(k):
        if k == 0: return 1.0
        return 0.30 * (1.0 / 3) ** k + 0.70 * (1 ** (k + 1) - 0.4 ** (k + 1)) / ((k + 1) * (1 - 0.4))
    def christoffel(t, n):
        G = np.array([[M_k(j + l) for l in range(n + 1)] for j in range(n + 1)])
        v = np.array([t ** j for j in range(n + 1)])
        return 1.0 / (v @ np.linalg.solve(G, v))
    ta, tg, tc = 1.0 / 3, 0.37, 0.6     # scaled atom (1/3), gap (0.37 ~ 1.1/3), continuum (0.6 ~ 1.8/3)
    print("  measure (t=s/3 in [0,1]): 0.30 delta_{1/3} + 0.70*Unif[0.4,1]  (edge atom, gap, continuum)")
    print("  %6s %18s %18s %18s" % ("n", "lambda_n(atom)", "lambda_n(gap)", "lambda_n(cont)"))
    for n in (2, 4, 8, 12):
        print("  %6d %18.6e %18.6e %18.6e" % (n, christoffel(ta, n), christoffel(tg, n), christoffel(tc, n)))
    print("  (n>12 needs high-precision arithmetic: the moment-Hankel G_n is exponentially ill-conditioned.)")
    print("  -> lambda_n(atom) -> w=0.30 (stays bounded below; rigorous liminf>=nu({x}));")
    print("     lambda_n(gap) -> 0 (isolation); lambda_n(cont) -> density.  Same atom read as [C1],")
    print("     but as n->inf (infinite data) AND needs x=m1^2 (R1+R2).  Hankel TOEPLITZ reading, [T]/[NG].")

    # ------------------------------------------------------------------ [S1] one-plaquette
    hdr("[S1] SOLVABLE ACID -- SU(2) one-plaquette Haar marginal rho1(z) (exact, SMOOTH -> f1^2=0)")
    g = np.linspace(-TWO_PI, TWO_PI, 400001); r = rho1_haar(g)
    print("  rho1(0)=%.6f finite; rho1'(0)=%.2e (smooth): CONTINUOUS density, no atom."
          % (rho1_haar(0.0), (rho1_haar(1e-4) - rho1_haar(-1e-4)) / 2e-4))
    print("  %10s %20s" % ("d", "M(d)@x=0 (should -> 0)"))
    for d in (0.5, 0.2, 0.05, 0.01):
        print("  %10.4g %20.6f" % (d, mass_in_window([], (g, r), 0.0, d)))
    print("  -> M(d) ~ 2 d rho1(0) -> 0 : a CORRECT tool returns 0 on the one-plaquette (atomless). [V]")

    # ------------------------------------------------------------------ [S2] free cut
    hdr("[S2] SOLVABLE ACID -- free massless cut proxy rho(s)=const on [0,1] (atomless -> f1^2=0)")
    gf = np.linspace(0.0, 1.0, 400001); rf = np.full_like(gf, 1.0)
    print("  %10s %20s" % ("d", "M(d)@x=0 (should -> 0)"))
    for d in (0.2, 0.05, 0.01, 0.003):
        print("  %10.4g %20.6f" % (d, mass_in_window([], (gf, rf), 0.0, d)))
    print("  -> M(d) ~ d -> 0 : ACID PASS in the continuum. [V]")
    print("     CAVEAT: at finite L the free field is finitely atomic (lowest mode ~ (2pi/L)^2),")
    print("     f1^2>0 at finite L, -> 0 only as L->inf.  The acid test is an L->inf statement. [T]")

    # ------------------------------------------------------------------ verdict
    hdr("VERDICT")
    print("  Discontinuous atom-reading functionals EXIST and are EXACT [T]: boundary jump [C1],")
    print("  sup-atom [C2], finite-moment blindness [C3], Christoffel/Hankel [C5].  They PASS the")
    print("  acid test (return 0 on atomless) -- but only through a LIMIT.  Failure modes:")
    print("    R1 location: [C1],[C5] must set x = m1^2 (the gap) before reading;")
    print("    R2 limit   : all need e->0 / d->0 / n->inf = infinite data (A3);")
    print("    R3 family  : [C2] targets a non-edge atom; on any finite-data class the atom is unreachable.")
    print("  R1=R2=R3 are ONE obstruction: atom = weak-* discontinuity; finite data continuous;")
    print("  atomless dense.  => NO discontinuous functional is computable from finite (L,a)-data")
    print("  without presupposing the gap/isolation.  [T]/[NG]")
    print("  EXACT MISSING INPUT: a uniform-in-(L,a) RESOLUTION bound d < I2 (isolation) at the edge,")
    print("  i.e. the atom floor f1^2>=r0 WITH isolation = H* itself, not a consequence. [O]")
    print("  Clay 0%.  No QED.")

if __name__ == "__main__":
    main()
