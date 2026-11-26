"""
RAT1c - INDEPENDENT re-implementation (own encoding, no ur1_enum bitfields).

Why: ur1_enum.pc() packs coordinates into 5-bit fields with OFF=20, i.e. only
x_i in [-20, 11] is legal; the chain family silently wraps for t >= 12
(measured b(n=13) = 194 > 14n+6 = 188, and b saturates at 473 for n >= 32).
This script uses plain tuples so large-n families are exact, and reproduces
the small-n enumerations as a cross-check of the neighbour rule.

Checks:
  (A) a_n = 1, 20, 458, 11132, 280175  (n = 1..5, root = plaquette (0,0,0,0;0,1))
  (B) max_A b(A) = 20, 34, 48, 62, 76  (n = 1..5)  and  b <= 14n+6 everywhere
  (C) chain n = 1..60: b = 14n+6 exactly      (sharpness of the theorem)
  (D) plate k x k: b  (k = 1..6)

Usage: py -3.12 rat1_indep.py [maxn]
"""
import sys
import time

DIRS = ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1))


def P(x, mu, nu):
    if mu > nu:
        mu, nu = nu, mu
    return (x[0], x[1], x[2], x[3], mu, nu)


def add(x, k):
    d = DIRS[k]
    return (x[0] + d[0], x[1] + d[1], x[2] + d[2], x[3] + d[3])


def sub(x, k):
    d = DIRS[k]
    return (x[0] - d[0], x[1] - d[1], x[2] - d[2], x[3] - d[3])


def neighbors(p):
    x = p[:4]
    mu, nu = p[4], p[5]
    links = ((x, mu), (add(x, mu), nu), (add(x, nu), mu), (x, nu))
    out = set()
    for (y, d) in links:
        for k in range(4):
            if k == d:
                continue
            out.add(P(y, d, k))          # plaquette at y in plane (d,k)
            out.add(P(sub(y, k), k, d))  # plaquette at y-k in plane (k,d)
    out.discard(p)
    return out


def boundary(A):
    SA = set(A)
    out = set()
    for p in A:
        out |= (neighbors(p) - SA)
    return out


def bnd(A):
    return len(boundary(A))


ROOT = P((0, 0, 0, 0), 0, 1)


def animals(maxn):
    levels = [[(ROOT,)]]
    for _ in range(1, maxn):
        nxt = set()
        for S in levels[-1]:
            base = set(S)
            for p in S:
                for q in neighbors(p):
                    if q not in base:
                        nxt.add(tuple(sorted(S + (q,))))
        levels.append([tuple(s) for s in nxt])
        print(f"   level n={len(levels)}: {len(nxt)} animals", flush=True)
    return levels


def chain(n):
    return [P((t, 0, 0, 0), 0, 1) for t in range(n)]


def plate(k):
    return [P((i, j, 0, 0), 0, 1) for i in range(k) for j in range(k)]


def main():
    maxn = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    t0 = time.time()
    print("RAT1c - independent encoding;  b(A) = boundary, Delta = 20")
    print("=" * 74)
    levels = animals(maxn)
    a = [len(L) for L in levels]
    corpus = [1, 20, 458, 11132, 280175]
    print(f"(A) a_n = {a}   corpus {corpus[:maxn]}  match={a == corpus[:maxn]}")
    print(f"(B) n, min b, mean b, max b, 14n+6, #viol(b>14n+6)")
    for n in range(1, maxn + 1):
        bs = [bnd(A) for A in levels[n - 1]]
        viol = sum(1 for x in bs if x > 14 * n + 6)
        print(f"    {n:>2}  {min(bs):>4} {sum(bs)/len(bs):>9.4f} {max(bs):>5} "
              f"{14*n+6:>6} {viol:>6}   #b<=48 = {sum(1 for x in bs if x <= 48)}", flush=True)

    print("\n(C) chain family, n = 1..60 (own encoding, no wrap):")
    ok = True
    for n in range(1, 61):
        b = bnd(chain(n))
        if b != 14 * n + 6:
            ok = False
            print(f"    n={n}: b={b}  14n+6={14*n+6}  MISMATCH")
    print(f"    b = 14n+6 exactly for all n = 1..60 : {ok}"
          f"   (b(60) = {bnd(chain(60))})")
    print("    sharpness =>  sup_n max_A b(A) = +infinity ;  b > 48 for every n >= 4")

    print("\n(D) plate k x k (coplanar):")
    for k in range(1, 7):
        n = k * k
        b = bnd(plate(k))
        print(f"    k={k}  n={n:>2}  b={b:>4}  b/n = {b/n:.3f}")
    print(f"\ntotal {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
