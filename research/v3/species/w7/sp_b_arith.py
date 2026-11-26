# -*- coding: utf-8 -*-
"""sp_b_arith.py -- SP-B interface witness (3 narrow checks, deterministic, single process).

[C1] the three R-classes <-> growth of the beta-Taylor coefficients of log W  (explicit dictionary)
[C2] the corpus u = c_1/c_0 sequence: reproduce / locate the crossing u=1
[C3] the middle class f = R^{-p}(log R)^q against (star_1): exponent matching + coefficient matching,
     under BOTH bookkeepings of the R-power (file's c1 f'' = c2 g2N f^2/R^3, and the dimensionally
     consistent c1 f''/R = c2 g2N f^2/R^3).
"""
import numpy as np
from scipy.special import iv
import mpmath as mp

mp.mp.dps = 30
OUT = []
def say(s=""):
    print(s); OUT.append(s)

say("=" * 76)
say("[C2] corpus u = c_1/c_0 in the R5-AUDIT convention  c_n = I_n(b) - I_{n+2}(b)")
say("=" * 76)
def c(n, b):
    return iv(n, b) - iv(n + 2, b)
def u(b):
    return c(1, b) / c(0, b)
say("  b        u=c1/c0         c0            note")
for b in [0.02, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 2.63, 3.0]:
    say("  %-6.2f  %.6f      %.6f" % (b, u(b), c(0, b)))
say("  lead's list 0.0104/0.247/0.480/1.014 vs mine: u(0.02)=%.4f u(0.5)=%.4f u(1.0)=%.4f u(2.5)=%.4f"
    % (u(0.02), u(0.5), u(1.0), u(2.5)))
# crossing u = 1
from scipy.optimize import brentq
bs = brentq(lambda x: u(x) - 1.0, 1.0, 4.0)
say("  ==> u = 1 at  b* = %.6f ;  u(2.5) = %.6f  (reproduces TRG-EKR's 1.014390 to 6 digits)"
    % (bs, u(2.5)))
say("  ==> c_0(2.5) = %.6f ; c_0/e^2.5 = %.6f (reproduces TRG-EKR's 0.165268)"
    % (c(0, 2.5), c(0, 2.5) / np.exp(2.5)))
say("  n-growth at b = 2.5 :  c_n vs geometric/factorial reference")
for n in [0, 2, 4, 8, 12, 16, 24]:
    v = c(n, 2.5)
    fact = (1.25 ** n) / float(mp.factorial(n))
    say("    n=%2d  c_n=%.6e   (b/2)^n/n!=%.6e   c_n/((b/2)^n/n!)=%.4f"
        % (n, v, fact, v / fact))
say("  => c_n decays FACTORIALLY in n (ratio -> const * (b/2)^n/n!), not geometrically, not as a power.")
say("     In the dictionary of [C1] that is sigma_0 = infinity, i.e. the deepest possible class (E).")

say()
say("=" * 76)
say("[C3] middle class f = R^{-p}(log R)^q against (star_1):  full algebra")
say("=" * 76)
import sympy as sp
R, p, q, A, C = sp.symbols("R p q A C", positive=True)
L = sp.log(R)

def f_of(p, q):
    return A * R ** (-p) * L ** q

def test(fun, lhs_kind):
    fpp = sp.diff(fun, R, 2)
    if lhs_kind == "file":            # c1 f''      = C f^2 / R^3
        lhs, rhs = fpp, C * fun ** 2 / R ** 3
    else:                             # c1 f''/R    = C f^2 / R^3   (dimensionally consistent)
        lhs, rhs = fpp / R, C * fun ** 2 / R ** 3
    ratio = sp.simplify(lhs / rhs)
    return sp.simplify(ratio)

for kind in ["file", "consistent"]:
    say("  bookkeeping: %s" % kind)
    for (pp, qq) in [(1, 0), (0, 0), (-1, 0), (-1, 1), (-1, -1), (2, 0.5)]:
        r = test(f_of(sp.Rational(pp) if isinstance(pp, int) else pp,
                      sp.Rational(qq) if isinstance(qq, int) else qq), kind)
        say("    p=%-4s q=%-5s   LHS/RHS = %s" % (pp, qq, sp.simplify(r)))

say()
say("  --- exponent matching in general (log-power bookkeeping), file's (star_1) ---")
say("  f'' ~ R^{-p-2}[ p(p+1)L^q - q(2p+1)L^{q-1} + q(q-1)L^{q-2} ]")
say("  RHS ~ C R^{-2p-3} L^{2q}")
say("  (i)  R-powers equal : -p-2 = -2p-3  =>  p = -1")
say("  (ii) with p=-1 the L^q term has coefficient p(p+1)=0, so LHS log-degrees = {q-1, q-2},")
say("       RHS log-degree = 2q.  Need 2q=q-1 AND 2q=q-2  =>  contradiction for every real q.")
say("  (iii) coefficient check at p=-1: f=R, f''=0 => LHS=0 while RHS=C g2N R^2/R^3 = C g2N/R != 0.")
say("  => even p=-1 (the file's 'only exponent') is NOT a solution: LHS vanishes identically.")
say()
say("  --- same test, dimensionally consistent bookkeeping c1 f''/R = C f^2/R^3 ---")
say("  (i)  R-powers equal : -p-3 = -2p-3  =>  p = 0")
say("  (ii) then f = A L^q, and LHS log-degrees {q, q-1, q-2}, RHS 2q  => again contradiction")
say("       (and for p=0 exactly, the p(p+1) term vanishes identically: LHS = q R^{-3}(...) vs RHS R^{-3})")
say("  => under BOTH bookkeepings: no solution in the whole two-parameter power-log family.")
say()
say("  --- area law check (file's sec.5 item 2) ---")
sig = sp.symbols("sigma", positive=True)
area = sp.exp(-sig * R ** 2)
say("  f = e^{-sigma R^2}:  LHS/RHS (file) = %s" % sp.simplify(sp.diff(area, R, 2) / (C * area ** 2 / R ** 3)))
say("  -> ratio = e^{+sigma R^2}*(4 sigma^2 R^2-2 sigma)*R^3/C  -> +infinity: the area law fails (star_1) too.")
say("  NOTE: for the area law the bilinear factor is NOT f^2 -- the two arcs' minimal surfaces add,")
say("        W(C_xy)W(C_yx) = W(C) = f.  So 'f^2' is class-dependent (right for power, wrong for area).")

say()
say("=" * 76)
say("[C1] the dictionary, written out")
say("=" * 76)
say("  Let kappa = log(entropy of surfaces of unit area)  [a combinatorial constant of the 4D lattice,")
say("  independent of beta and NOT a functional of c_n].  Let sigma_0 = lim_n (-log c_n)/n and")
say("  alpha the sub-exponential exponent of c_n ~ e^{-kappa n} n^{alpha-1}.")
say("     (E) W ~ e^{-sigma R^2}                <=> sigma_0 > kappa          (a GAP in the ratio test)")
say("     (P) W ~ R^{-p}                        <=> sigma_0 = kappa exactly, and c_n ~ e^{-kappa n} n^{-p/2-1}")
say("     (I) W ~ R^{-p}(log R)^q               <=> sigma_0 = kappa exactly, and c_n ~ e^{-kappa n} n^{-p/2-1}(log n)^q")
say("  Derivation of the exponent map: W(R) = sum_{A <= A_max(R)} Omega_A prod_p c_{j_p}; with Omega_A ~ e^{kappa A}")
say("  and c ~ e^{-kappa A} A^{alpha-1}, the summand is A^{alpha-1}, and A_max ~ R^2, so W ~ R^{2 alpha}.")
say("  Matching R^{-p} = R^{2 alpha} gives alpha = -p/2, i.e. c_n ~ e^{-kappa n} n^{-p/2-1}.  q passes through as-is.")
say("  What fixes (P) vs (I) is alpha and q per se, NOT kappa; but the (E) vs (P),(I) border is sigma_0 vs kappa,")
say("  and kappa is NOT computable from c_n.")

open(__file__.replace(".py", ".log"), "w", encoding="utf-8").write("\n".join(OUT) + "\n")
