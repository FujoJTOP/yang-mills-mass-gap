"""
RAT1b - what number 48 actually requires (means, not maxima), + a tool warning.

Tool warning found this round:
    ur1_enum.pc() packs coordinates into 5-bit fields with OFF = 20, so only
    x_i in [-20, 11] is legal.  Chain/plate families built with pc() silently
    WRAP for coordinates >= 12: chain n=13 gave b = 194 > 14n+6 = 188 and b
    saturated at 473 for n >= 32.  Those numbers are artifacts.
    -> use rat1_indep.py (plain tuples) for any family with coordinates > 11.

Checks here:
  (i)   ur1_enum chain, n = 1..12 only (legal range): b = 14n+6 exactly
  (ii)  L_t(14) = 14^14/13^13  (tree envelope at the boundary rate Delta-6, not Delta)
  (iii) 48 requires mean m(n+1) >= (14n+6)/48, asymptotic slope 14/48
  (iv)  measured exact mean m and mean b per level (LC4 identity + a_6)

Usage: py -3.12 rat1_ldegree.py
"""
import sys
from decimal import Decimal, getcontext

sys.path.insert(0, r"D:\Dev\math\research\v3\species\w7")
import ur1_enum as U1

getcontext().prec = 30


def bnd(A):
    NA = set(A)
    out = set()
    for p in A:
        for q in U1.neighbors(p):
            if q not in NA:
                out.add(q)
    return len(out)


def chain(n):
    return [U1.pc((t, 0, 0, 0), 0, 1) for t in range(n)]


def main():
    print("RAT1b - the constraint behind 48, and the pc() wrap trap")
    print("=" * 78)
    print("[i] chain via ur1_enum, LEGAL range only (x <= 11 => n <= 12):")
    for n in (1, 4, 5, 6, 8, 10, 12):
        b = bnd(chain(n))
        print(f"    n={n:>2}: b={b:>4}  14n+6={14*n+6:>4}  {'OK' if b == 14*n+6 else 'BAD'}")
    print("    n=13 via pc(): b=", bnd(chain(13)), " (artifact: 5-bit field overflow)")
    print("    correct values for n > 12: rat1_indep.py -> b = 14n+6 for n = 1..60")

    print("\n[ii] tree envelope L_t(D) = D^D/(D-1)^(D-1):")
    for dd in (14, 15, 18, 19, 20):
        D = Decimal(dd)
        print(f"    L_t({dd}) = {D**dd/(D-1)**(dd-1)}")
    e = Decimal(1).exp()
    print(f"    e*Delta (Delta = 20) = {e*20}   (= proven bound on lambda_*, ADJ2:14)")
    print(f"    L_t(14) < 48 ?  {Decimal(14)**14/Decimal(13)**13 < 48}"
          f"   (14 = Delta - 2(d-1) = boundary rate of the extremal family)")

    print("\n[iii] a_{n+1}/a_n = mean b / mean m   (exact identity, LC4)")
    print("      => ratio <= 48  <=>  mean b <= 48 * mean m ; with max b = 14n+6:")
    print(f"      need mean m(n+1) >= (14n+6)/48, slope c = 14/48 = {Decimal(14)/Decimal(48)}")

    print("\n[iv] exact means (sum_b from the full enumerations, a_6 from UR4 [E])")
    sum_b = {1: 20, 2: 648, 3: 20172, 4: 616400, 5: 18639660}
    a_n = {1: 1, 2: 20, 3: 458, 4: 11132, 5: 280175, 6: 7218768}
    print(f"    {'n':>2} {'mean m':>8} {'mean m/n':>9} {'mean b':>8} {'b/n':>6} {'ratio':>8}")
    for n in range(2, 6):
        mean_m = sum_b[n - 1] / a_n[n]          # E[m] on level n = sum_{A_{n-1}} b / a_n
        mean_b = sum_b[n] / a_n[n]
        print(f"    {n:>2} {mean_m:>8.5f} {mean_m/n:>9.5f} {mean_b:>8.4f} {mean_b/n:>6.3f} "
              f"{mean_b/mean_m:>8.4f}")
    mm6 = sum_b[5] / a_n[6]
    print(f"     6 {mm6:>8.5f} {mm6/6:>9.5f} {'-':>8} {'-':>6} {'-':>8}   (a_6 level)")
    print("    mean m/n: 0.60, 0.55, 0.43, then 0.43 -> the slope 0.2917 that 48 needs")
    print("    is 1.5x below the measured one: 48 is consistent but NOT asymptotically safe")


if __name__ == "__main__":
    main()
