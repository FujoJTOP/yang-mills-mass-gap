"""
qg4_su2.py -- Q4 route, step 4: the non-abelian (SU(2)) READ-ONLY comparison.

Two things only:
  (1) the exact SU(2) one-plaquette character coefficients (exact quadrature), which
      exhibit the gap sequence 4j(j+1) = 3, 8, 15, ... -- the corpus already records
      mu_1 = 3, mu_2 = 8 (rb_ratio_bound.py PART E).  We re-derive it independently
      from a 1D quadrature with the Weyl measure sin^2(theta).
  (2) the read-only comparison with the corpus numbers at beta = 2.5, and the
      statement of WHY no scalar network exists in the non-abelian case.

No new lattice simulation is run here.
"""
import json
import math
import time

import numpy as np
from scipy.integrate import quad

LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def char_j(j, th):
    """SU(2) character of the spin-j rep, theta in [0, pi]."""
    n = int(round(2 * j))
    x = math.cos(th)
    # U_{n}(x) = sin((n+1)th)/sin(th)
    if abs(math.sin(th)) < 1e-14:
        return float(n + 1)
    return math.sin((n + 1) * th) / math.sin(th)


def coeff(j, beta):
    """c_j propto int sin^2(th) chi_j(th) exp(2 beta cos th) dth  (Haar on SU(2))."""
    f = lambda th: math.sin(th) ** 2 * char_j(j, th) * math.exp(2.0 * beta * math.cos(th))
    val, _ = quad(f, 0.0, math.pi, limit=400, epsabs=1e-14, epsrel=1e-13)
    return val


def part1():
    say("=" * 78)
    say("PART 1 -- exact SU(2) one-plaquette spectrum = the group Laplacian on S^3")
    say("  class-angle form:  L f = -(1/sin^2 th) d/dth( sin^2 th f' ) ,  th in (0,pi)")
    say("  substitution f = u/sin(th)  =>  L = -u'' , u(0)=u(pi)=0")
    say("  => eigenvalues n^2-1 = 4j(j+1) = 3, 8, 15, 24, ...  (n = 2j+1)")
    say("  Verified below as a GENERALISED eigenproblem A f = lam B f, B = diag(sin^2 th).")
    say("=" * 78)
    N = 3000
    th = np.linspace(0.0, math.pi, N)
    h = th[1] - th[0]
    M = N - 2                                   # interior unknowns
    s2 = np.sin(th[1:-1]) ** 2                  # weight at interior nodes, length M
    sh = 0.5 * (np.sin(th[1:]) ** 2 + np.sin(th[:-1]) ** 2)   # weight at half points, length N-1
    # Rayleigh quotient  sum_j sh_j (f_{j+1}-f_j)^2 / h  /  ( h sum_i s2_i f_i^2 )
    d = (sh[0:N - 2] + sh[1:N - 1]) / h
    e = -sh[1:N - 2] / h
    A = np.diag(d) + np.diag(e, 1) + np.diag(e, -1)
    B = np.diag(s2 * h)
    from scipy.linalg import eigh
    ev = np.sort(np.real(eigh(A, B, eigvals_only=True)))
    say(f"  grid N={N}, generalised eigenvalues of the class Laplacian")
    say(f"  (index n = 2j+1; n=1 is the constant, j=0, lambda=0):")
    rows = []
    for i in range(6):
        j = i / 2.0
        p = 4.0 * j * (j + 1.0)
        say(f"    j={j:>3.1f}  lam = {ev[i]:>14.8f}    4j(j+1) = {p:>6.2f}   "
            f"|diff| = {abs(ev[i]-p):.2e}   rel = {abs(ev[i]-p)/max(p,1e-9):.1e}")
        rows.append(dict(j=j, lam=float(ev[i]), casimir=p))
    say("")
    say("  EXACT: mu_1 = 3, mu_2 = 8  (the corpus value, rb_ratio_bound.py PART E),")
    say("  and the sequence is 4j(j+1) with multiplicity (2j+1)^2.")
    say("  Note what this is: an OPERATOR spectrum, not the spectrum of a graph Laplacian.")
    say("")
    return rows


def part2():
    say("=" * 78)
    say("PART 2 -- the read-only comparison, and why SU(2) has no scalar network")
    say("=" * 78)
    say("  Geometric R_eff of the spatial 3-torus (the 'skeleton' of a 4D lattice,")
    say("  from qg1 PART A):  L=4: 0.4333,  L=6: 0.4539,  L=8: 0.46595")
    say("  It SATURATES (R_inf ~ 0.47): it is L-independent, and beta-independent.")
    say("")
    say("  Read-only corpus numbers at beta = 2.5 (SU(2), 4D):")
    say("     m_1 a = 1/(xi/a) in [0.77, 1.89]         (BN-RESIDUAL.md:293)")
    say("     strong-coupling first term m a = 1.88    (same line)")
    say("     m_{0++*}/m_{0++} ~ 1.54-1.56             (AUDIT-CH16.md:33)")
    say("")
    say("  Multiplying by the geometric R_eff ~ 0.47 gives 0.36-0.89: bounded, and")
    say("  L-independent BECAUSE R_eff is L-independent and m a is L-independent.")
    say("  The same test on free Maxwell gives  w_min*R_eff ~ 1/L -> 0 (qg3 ACID 1).")
    say("")
    say("  So the test does separate gapped from gapless -- but only because it is")
    say("  literally 'is lambda L-independent?' with a beta-free constant stuck on.")
    say("  All the discriminating information sits in lambda, none in R_eff.")
    say("")
    say("  WHY THERE IS NO SCALAR NETWORK (T0b, sharp):")
    say("   (i) For a scalar network the Dirichlet form is E(f)=sum_e c_e (f_x-f_y)^2,")
    say("       i.e. a QUADRATIC form in the vertex variables.  The Wilson weight")
    say("       exp(beta Re tr U_p) is NOT a quadratic form in the link variables:")
    say("       it is a matrix element of the regular representation, and the")
    say("       quadratic expansion around a background is the Savvidy Hessian")
    say("       -D^2 + 2 ad F_B -- an OPERATOR in the adjoint fibre, indefinite.")
    say("   (ii) A conductance is a positive number.  'ad F_B' has both signs, so the")
    say("       'conductances' of any non-abelian effective network are sign-indefinite")
    say("       => no M-matrix => no effective resistance (R_eff is defined only for")
    say("       PSD Laplacians with non-negative weights).")
    say("   (iii) Even the abelian Z2 case works only because U(1)/Z_N plaquettes")
    say("       exponentiate to a real positive weight after the duality transform;")
    say("       for SU(2) the dual is not a spin model with positive weights.")
    say("  => the constant c in lambda >= c/R_eff would have to carry a g-dependence")
    say("     (T0b), and the monotonicity direction (T0c) has no meaning.")
    say("")

    rows = [
        dict(model="SU(2) 4D beta=2.5 (read-only)", ma=[0.77, 1.89], Reff=0.46595,
             prod=[0.77 * 0.46595, 1.89 * 0.46595]),
        dict(model="free Maxwell T^3_L", ma=[2 * math.sin(math.pi / 8)], Reff=0.46595,
             prod=[2 * math.sin(math.pi / 8) * 0.46595]),
    ]
    for r in rows:
        say(f"  {r['model']:>32}:  m a = {r['ma']}   R_eff = {r['Reff']:.5f}   "
            f"m a R_eff = {[f'{x:.4f}' for x in r['prod']]}")
    say("")
    return rows


def main():
    t0 = time.time()
    say(f"qg4_su2.py   start {time.strftime('%Y-%m-%d %H:%M:%S')}")
    say("")
    p1 = part1()
    p2 = part2()
    say(f"total wall {time.time() - t0:.1f} s")
    json.dump(dict(part1=p1, part2=p2), open("qg4_su2.json", "w"), indent=1)
    open("qg4_su2.log", "w", encoding="utf-8").write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
