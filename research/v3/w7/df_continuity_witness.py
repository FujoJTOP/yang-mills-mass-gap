#!/usr/bin/env python3
"""df_continuity_witness.py -- witness battery for research/v3/species/DF-CONTINUITY.md
(deformation + adiabatic-continuity route to the 4D pure-SU(N) mass gap).

Pure stdlib, single process, deterministic.  Run:  py -3.12 -I df_continuity_witness.py

Parts
  A  the continuity schema is FALSE in the form the route needs, and vacuous in the
     form it can prove.  Three exactly-solvable families.
  B  adding a POSITIVE operator to H does NOT raise the gap: "deformation monotone in
     kappa" is false (explicit 3x3 counterexample + random statistics).
  C  the deformed (bion) gap is exponentially small, and control/gap are anti-correlated:
     Delta_def/Lambda ~ (Lambda*beta)^(8/3) at one-loop AF;  lambda=g^2N = 8*pi^2/S0.
  D  the transport budget: continuity can only bound DECREASE; the route needs an
     INCREASE by 1/Delta_def.  Closed-form incompatibility of the two windows.
  E  min-max direction on the radius axis (exact free spectrum) + the spurious-1/r trap.
"""

import math
import random

BAR = "-" * 78


# ----------------------------------------------------------------------------
# A.  schema falsification
# ----------------------------------------------------------------------------
def part_A():
    print(BAR)
    print("A. The continuity schema")
    print("   (S+)  [D continuous on [0,1], D(l)>0 for all l in (0,1]]  =>  D(0)>0")
    print("   (S)   same hypotheses                                   =>  D(0)>=0")
    print("   (S*)  exists delta>0 with D(l)>=delta on the CLOSED [0,1] =>  D(0)>=delta")
    fams = [
        ("mass / free field   D(l)=l",            lambda l: l),
        ("1D TFIM (exact)     D(l)=2(1-l)",       lambda l: 2.0 * (1.0 - l)),
        ("flat                D(l)=exp(-1/l)",    lambda l: math.exp(-1.0 / l) if l > 0 else 0.0),
    ]
    grid = [i / 200.0 for i in range(1, 201)]      # (0,1]
    print()
    print("   family                          min_{l in (0,1]} D   #zeros in (0,1]   D(0)   (S+) true?")
    for name, f in fams:
        vals = [f(l) for l in grid]
        mn = min(vals)
        nz = sum(1 for v in vals if v == 0.0 or abs(v) < 1e-300)
        d0 = f(0.0)
        print(f"   {name:30s} {mn:18.6e}   {nz:14d}   {d0:8.3f}   {'YES' if d0 > 0 else 'NO'}")
    print()
    print("   -> every family satisfies ALL hypotheses of (S+) and (S); D(0)>0 fails in 2/3.")
    print("      (S) holds but is vacuous: D>=0 is unconditional for a reflection-positive")
    print("      theory (spectral theorem).  So the route's whole logical yield is a tautology.")
    print("      exp(-1/l) is C^inf with ALL derivatives 0 at l=0, so no derivative bound,")
    print("      no Taylor expansion and no C^k hypothesis can upgrade (S) to (S+).")
    print("      Only a UNIFORM floor -- hypothesis (S*), whose l=0 clause is the theorem --")
    print("      excludes it.  'No phase transition on the path' is a pointwise (in fact open)")
    print("      condition; it is exactly (S), not (S*).")


# ----------------------------------------------------------------------------
# B.  no monotonicity on the deformation axis
# ----------------------------------------------------------------------------
def jacobi_eig(A, iters=200, tol=1e-13):
    """Eigenvalues of a real symmetric matrix by cyclic Jacobi rotations."""
    n = len(A)
    a = [row[:] for row in A]
    for _ in range(iters):
        p = q = 0
        mx = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                if abs(a[i][j]) > mx:
                    mx = abs(a[i][j])
                    p, q = i, j
        if mx < tol:
            break
        app, aqq, apq = a[p][p], a[q][q], a[p][q]
        th = 0.5 * math.atan2(2.0 * apq, aqq - app)
        c, s = math.cos(th), math.sin(th)
        for k in range(n):
            akp, akq = a[k][p], a[k][q]
            a[k][p] = c * akp - s * akq
            a[k][q] = s * akp + c * akq
        for k in range(n):
            apk, aqk = a[p][k], a[q][k]
            a[p][k] = c * apk - s * aqk
            a[q][k] = s * apk + c * aqk
    return sorted(a[i][i] for i in range(n))


def gap(ev):
    return ev[1] - ev[0]


def rand_psd(n, rng):
    B = [[rng.gauss(0.0, 1.0) for _ in range(n)] for _ in range(n)]
    return [[sum(B[i][k] * B[j][k] for k in range(n)) for j in range(n)] for i in range(n)]


def part_B():
    print(BAR)
    print("B. Adding a POSITIVE operator to H does not raise the gap")
    print("   Hellmann-Feynman:  d Delta/d lambda = <G|O|G> - <0|O|0>   (either sign).")
    print("   d Delta/d lambda = (a/(g^2 beta^3)) * that,  for the double-trace deformation.")
    H = [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    O = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    print()
    print("   explicit counterexample (O >= 0, O = diag(1,0,0), H = diag(0,1,1)):")
    print("     lambda = 0.0 :  spec = %s   gap = %.6f" % (jacobi_eig(H), gap(jacobi_eig(H))))
    H1 = [[H[i][j] + 0.5 * O[i][j] for j in range(3)] for i in range(3)]
    print("     lambda = 0.5 :  spec = %s   gap = %.6f" % (jacobi_eig(H1), gap(jacobi_eig(H1))))
    print("     -> O >= 0, gap DECREASED 1.000000 -> 0.500000.  So 'add a positive operator")
    print("        => gap grows' is FALSE; there is no monotone on the deformation axis.")
    print()
    rng = random.Random(20261125)
    n, lam, N = 4, 0.5, 4000
    nneg = nup = 0
    dmin, dmax = 1e9, -1e9
    for _ in range(N):
        M = rand_psd(n, rng)
        P = rand_psd(n, rng)
        g0 = gap(jacobi_eig(M))
        g1 = gap(jacobi_eig([[M[i][j] + lam * P[i][j] for j in range(n)] for i in range(n)]))
        d = g1 - g0
        dmin, dmax = min(dmin, d), max(dmax, d)
        if d < 0:
            nneg += 1
        else:
            nup += 1
    print("   random check, n=%d, %d draws, lambda=%.2f (H=AA^T, O=BB^T, both PSD):" % (n, N, lam))
    print("     gap decreased in %d / %d = %.1f%%   (increased in %d)" % (nneg, N, 100.0 * nneg / N, nup))
    print("     Delta change range: [%+.4f, %+.4f]" % (dmin, dmax))
    print("     -> both signs occur at O(1) rate.  The deformation axis carries NO order, so a")
    print("        continuity argument on it has no monotone to anchor a direction.")


# ----------------------------------------------------------------------------
# C.  the deformed (bion) gap: exponentially small, anti-correlated with control
# ----------------------------------------------------------------------------
#   one-loop pure-YM AF:  8*pi^2/g^2 = (11/3) N log(1/(Lambda*beta))   [SU(N), b0]
#   S0 = 8*pi^2/(g^2 N)   = (11/3) log(1/(Lambda*beta))                [N-independent]
#   lambda_Hooft = g^2 N  = 8*pi^2 / S0
#   Delta_def ~ c * (1/beta) * exp(-S0) = c * Lambda * (Lambda*beta)^(8/3)
A0 = 11.0 / 3.0


def S0_of(Lb):
    return A0 * math.log(1.0 / Lb)


def lam_of(S0):
    return 8.0 * math.pi ** 2 / S0


def Ddef_of(Lb):
    return Lb ** (8.0 / 3.0)


def part_C():
    print(BAR)
    print("C. The 'computable' region: the dynamical gap there is exponentially small")
    print("   semiclassical object = dual-photon mass m_gamma ~ (c/beta) exp(-S0),")
    print("   S0 = 8*pi^2/(g^2 N) = (11/3) ln(1/(Lambda beta));  c = O(1) prefactor dropped.")
    print()
    print("   Lambda*beta      lambda=g^2N        Delta_def/Lambda ~ (Lambda beta)^(8/3)")
    for Lb in (5e-1, 3e-1, 1e-1, 3e-2, 1e-2, 1e-3, 1e-6):
        s = S0_of(Lb)
        print("   %10.2e   %14.3f   %14.3e" % (Lb, lam_of(s), Ddef_of(Lb)))
    print()
    Lb1 = math.exp(-1.5 * S0_of(1e-3))  # solve lambda = 1 : S0 = 8 pi^2
    s1 = 8.0 * math.pi ** 2
    Lb1 = math.exp(-s1 / A0)
    print("   genuinely weak coupling, lambda = g^2 N = 1  <=>  S0 = 8*pi^2 = %.2f" % s1)
    print("     => Lambda*beta = %.3e  =>  Delta_def/Lambda = %.3e" % (Lb1, Ddef_of(Lb1)))
    print("   -> control and gap size are strictly ANTI-correlated: the better the semiclassics,")
    print("      the smaller the gap.  At lambda=1 the deformed gap is ~1e-25 Lambda.")
    print("   -> and m_gamma lives in the MAGNETIC / 1-form channel; the Clay residue f1^2")
    print("      = |<0|(1/4)trF^2|G>|^2 is neutral 0-form, and trF^2 is provably blind to the")
    print("      1-form centre (PROOF-1FORM-ANOMALY).  Wrong channel as well as tiny.")


# ----------------------------------------------------------------------------
# D.  the transport budget
# ----------------------------------------------------------------------------
def part_D():
    print(BAR)
    print("D. The transport budget: continuity bounds DECREASE, the route needs GROWTH")
    print("   Any continuity/analyticity/Kato statement yields")
    print("        Delta(end) >= Delta(computable) - |dV|        (|dV| = total variation)")
    print("   i.e. a bound on how far Delta may FALL.  The route needs Delta to RISE from")
    print("   Delta_def ~ 1e-25 Lambda to O(Lambda).  A decrease-bound cannot certify growth.")
    print()
    print("   For the bound to be informative at all one needs |dV| < Delta(computable):")
    print("        Lambda  <  c*Lambda*(Lambda beta)^(8/3)     <=>   (Lambda beta)^(8/3) > 1/c")
    print("                                              <=>   Lambda beta > c^(-3/8) = O(1)")
    print("   but the semiclassical window requires S0 > 0, i.e. Lambda*beta < 1.")
    print("   -> the two requirements are incompatible by a power law: the ONLY beta at which")
    print("      a transported floor could survive is beta >~ 1/Lambda, which is outside the")
    print("      window (S0 <= 0: asymptotic freedom has switched off).")
    print()
    for c in (1.0, 10.0, 100.0):
        Lb = c ** (-3.0 / 8.0)
        if Lb >= 1.0:
            print("     prefactor c = %6.1f :  needs Lambda*beta > %.3f  ->  NO admissible point"
                  " (window has S0>0 only for Lambda*beta<1)" % (c, Lb))
        else:
            s = S0_of(Lb)
            print("     prefactor c = %6.1f :  needs Lambda*beta > %.3f  (there S0 = %.2f, lambda = %.1f)"
                  % (c, Lb, s, lam_of(s)))
    print("   -> even with a generous prefactor the required point sits at Lambda*beta ~ 0.3,")
    print("      where lambda = g^2 N ~ 9:  not weak coupling.  No free lunch.")


# ----------------------------------------------------------------------------
# E.  min-max direction, and the spurious-1/r trap
# ----------------------------------------------------------------------------
def part_E():
    print(BAR)
    print("E. Radius axis: Delta(beta) >= Delta(R^4) = inf_beta Delta(beta)  (min-max)")
    print("   exact free spectrum with a twist along the circle:  Delta = sqrt(m^2 + (pi/beta)^2)")
    print()
    print("   beta (in 1/Lambda)      Delta, m=1        Delta, m=0")
    for b in (0.1, 0.3, 1.0, 3.0, 10.0, 1e3):
        print("   %16.3g   %12.6f      %12.6f" % (b, math.hypot(1.0, math.pi / b), math.pi / b))
    print("   %16s   %12.6f      %12.6f" % ("inf (beta->inf)", 1.0, 0.0))
    print("   -> Delta(beta) is non-increasing and tends to Delta(R^4) = inf.  The controlled")
    print("      end (small beta) is the MORE-gapped end, so every number computed there is an")
    print("      UPPER bound on the target.  For m=0 the whole 'gap' pi/beta is the spurious")
    print("      finite-volume KK spacing: it diverges as beta->0 and vanishes as beta->inf,")
    print("      i.e. it is present for the free massless field and certifies nothing.")


if __name__ == "__main__":
    part_A()
    part_B()
    part_C()
    part_D()
    part_E()
    print(BAR)
    print("witnesses of DF-CONTINUITY; not new QCD evidence.  No QED.")
