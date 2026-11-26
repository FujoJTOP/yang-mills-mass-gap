"""tn1 part 3: rigorous upper bound for the rooted tree-sum t_n, and the numbers.

Claim [T]:  t_n <= B^n with B = 19^19/18^18, exactly, where Delta = 20 is the
plaquette degree.  Reason: every tree animal (Gamma,T) with |Gamma| = n is an
embedded subtree of the plaquette graph G rooted at p0; T -> its lift in the
universal cover (the 20-regular tree) is injective (T is simply connected, the
lift is unique, and T is the image of the lift), so
   t_n <= #{subtrees of the 20-regular tree with n vertices containing a fixed vertex}
       = [x^n] F(x),  F = x(1+G)^20,  G = x(1+G)^19,
and F(x) = y(y-1) with x = (y-1)/y^19 has its (unique, dominant) singularity at
y = 19/18, x_s = (1/18)(18/19)^19  =>  growth constant 1/x_s = 19^19/18^18.
"""
from fractions import Fraction as F
import math
from decimal import Decimal, getcontext

getcontext().prec = 40


def main():
    B = F(19) ** 19 / F(18) ** 18
    print("B = 19^19/18^18 =", B, "=", float(B))
    print("B (40 digits) =", Decimal(B.numerator) / Decimal(B.denominator))
    Delta = 20
    eD = Decimal(1).exp() * Delta
    print("e*Delta =", eD)
    print("B/eD   =", (Decimal(B.numerator) / Decimal(B.denominator)) / eD)
    print("48/B   =", Decimal(48) / (Decimal(B.numerator) / Decimal(B.denominator)))
    print("B/48   =", (Decimal(B.numerator) / Decimal(B.denominator)) / 48)
    print()
    t = [1, 20, 570, 18776, 669215]
    a = [1, 20, 458, 11132, 280175]
    Be = float(B)
    print(" n |          t_n | t_n^(1/n) |    B^n bound | log t_n/log B |  t_n/t_{n-1}")
    for i, tn in enumerate(t):
        n = i + 1
        print("%2d | %12d | %9.6f | %12.4e | %13.6f | %s"
              % (n, tn, tn ** (1.0 / n),
                 Be ** n,
                 math.log(tn) / math.log(Be) if n > 1 else 0.0,
                 ("%.4f" % (tn / t[i - 1])) if i else "-"))
    print()
    print(" mean tau = t_n/a_n:", ["%d/%d=%.6f" % (t[i], a[i], t[i] / a[i]) for i in range(5)])
    print(" 1/n-corrected extrapolation of ratios r_n = t_n/t_{n-1}:")
    r = [t[i] / t[i - 1] for i in range(1, len(t))]
    for i in range(1, len(r)):
        # L - A/n through (n-1, r[i-1]), (n, r[i])  with n = i+1
        n1, n2 = i, i + 1
        A = (r[i] - r[i - 1]) / (1.0 / n1 - 1.0 / n2)
        L = r[i] + A / n2
        print("   from n=%d,%d : L = %.4f   (A = %.3f)" % (n1, n2, L, A))
    print()
    ra = [a[i] / a[i - 1] for i in range(1, len(a))]
    for i in range(1, len(ra)):
        n1, n2 = i, i + 1
        A = (ra[i] - ra[i - 1]) / (1.0 / n1 - 1.0 / n2)
        print("   animals: from n=%d,%d : L_a = %.4f" % (n1, n2, ra[i] + A / n2))
    print()
    print("threshold 2.4*Delta =", 2.4 * Delta)


main()
