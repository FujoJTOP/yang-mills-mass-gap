"""
LC4 — growth positions of plaquette animals (Delta = 20), target lambda_* <= 48

Reuses UR1's exact enumeration (same dir):  a_n = # connected plaquette sets
containing the fixed plaquette p0=(0,(0,1))  [4D hypercubic, edge = shared link].

New quantities (this script):
  b(A)  := |boundary of A| = #{q not in A : q adjacent to some p in A}   ("growth positions")
  m(B)  := #{p in B, p != p0 : B\{p} connected}                          ("parents")
  identity (exact):   sum_{A in A_n} b(A)  =  sum_{B in A_{n+1}} m(B)
  => a_{n+1}/a_n = mean b / mean m     (so NO uniform "growth position" bound exists
     unless sup_n max b <= C ; measure max b).

Arithmetic: e*Delta, L_t(Delta)=Delta^Delta/(Delta-1)^(Delta-1) for Delta=18..21,
Euler-number-free check of "lambda* <= 48  <=>  effective degree <= 18".

Usage: py -3.12 lc4_growth.py [maxn]
"""
import sys, time, json
from decimal import Decimal, getcontext
import ur1_enum as U1

getcontext().prec = 40


def boundary(A):
    """b(A): vertices outside A adjacent to A"""
    NA = set(A)
    out = set()
    for p in A:
        for q in U1.neighbors(p):
            if q not in NA:
                out.add(q)
    return len(out)


ROOT = U1.pc((0, 0, 0, 0), 0, 1)


def parents(B):
    r"""m(B): p in B, p != root, with B\{p} connected"""
    root = ROOT
    assert root in B
    S = set(B)
    cnt = 0
    for p in B:
        if p == root:
            continue
        rest = S - {p}
        seen = {root}
        stack = [root]
        while stack:
            v = stack.pop()
            for q in U1.neighbors(v):
                if q in rest and q not in seen:
                    seen.add(q)
                    stack.append(q)
        if len(seen) == len(rest):
            cnt += 1
    return cnt


def main():
    maxn = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    t0 = time.time()
    print(f"LC4 growth positions, 4D plaquette graph, Delta=20, maxn={maxn}")
    levels = U1.animals(maxn)                      # levels[n-1] = A_n
    a = [len(L) for L in levels]
    print("a_n =", a, f"  [{time.time()-t0:.1f}s]\n")

    rows = []
    for n in range(1, maxn + 1):
        t1 = time.time()
        L = levels[n - 1]
        bs = [boundary(A) for A in L]
        rows.append(dict(n=n, a_n=len(L), sum_b=sum(bs), mean_b=sum(bs) / len(bs),
                         max_b=max(bs), min_b=min(bs),
                         hist={str(k): bs.count(k) for k in sorted(set(bs))}))
        r = rows[-1]
        ratio = (a[n] / a[n - 1]) if n >= 1 and n < len(a) else float('nan')
        print(f"n={n}: a_n={r['a_n']:>8d}  |dA| mean={r['mean_b']:8.4f} "
              f"max={r['max_b']:>4d} min={r['min_b']:>3d}  sum={r['sum_b']:>10d}  "
              f"a_n+1/a_n={ratio:8.4f}  [{time.time()-t1:.1f}s]", flush=True)

    # identity check:  sum_{A_n} b = sum_{A_{n+1}} m   (only where level n+1 known)
    print("\n--- identity  sum_{A_n} b(A) = sum_{A_{n+1}} m(B) ---")
    for n in range(1, maxn):
        t1 = time.time()
        L2 = levels[n]                              # A_{n+1}
        ms = [parents(B) for B in L2]
        lhs, rhs = rows[n - 1]["sum_b"], sum(ms)
        print(f"n={n}: sum_b(A_{n})={lhs:>10d}   sum_m(A_{n+1})={rhs:>10d}   "
              f"mean m={sum(ms)/len(ms):.5f}  max m={max(ms)}  "
              f"{'OK' if lhs == rhs else 'MISMATCH'}  [{time.time()-t1:.1f}s]", flush=True)

    # mean m on level 6 free of charge:  sum_b(A_5)/a_6   (a_6 = UR1, ur1_enum_n6.log)
    if maxn >= 5:
        A6 = 7218768
        print(f"\nmean m(A_6) = sum_b(A_5)/a_6 = {rows[-1]['sum_b']}/{A6} = "
              f"{rows[-1]['sum_b']/A6:.5f}   [must be >= 1 : a check on a_6 = 7218768]")

    # ---- constants ----
    print("\n--- constants (40 dps) ---")
    D = Decimal(20)
    e = Decimal(1).exp()
    print("Delta                     =", D)
    print("e*Delta                   =", (e * D))
    for dd in (17, 18, 19, 20, 21):
        d = Decimal(dd)
        lt = d ** dd / (d - 1) ** (dd - 1)
        print(f"L_t({dd}) = {dd}^{dd}/{dd-1}^{dd-1}".ljust(28), "=", lt)
    print("(Delta-1)*e               =", (D - 1) * e)
    print("e*18                      =", e * 18)
    print("48                        =", Decimal(48))
    print("48/Delta                  =", Decimal(48) / D)
    print("2*Delta+8                 =", 2 * D + 8)
    print("ratio e*Delta/48          =", (e * D) / 48)
    print("ratio L_t(20)/48          =", (Decimal(20) ** 20 / Decimal(19) ** 19) / 48)
    print("ratio L_t(19)/48          =", (Decimal(19) ** 19 / Decimal(18) ** 18) / 48)
    print("\ntotal", f"{time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
