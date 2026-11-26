"""
RAT1 - is the local-ratio lemma  a_{n+1} <= 48 a_n  true / provable ?

Object: 4D plaquette graph (vertex = plaquette, edge = shared link), Delta = 20.
a_n = #{A connected, p0 in A, |A| = n}.   a_n = 1,20,458,11132,280175,7218768.

Quantities (exact, infinite-lattice boundary - no torus involved):
  b(A) = |{q not in A : q adjacent to A}|          ("growth positions" / boundary tiles)
  m(B) = #{p in B, p != p0 : B\{p} connected}      ("parents")
  identity  sum_{A in A_n} b(A) = sum_{B in A_{n+1}} m(B)
  =>  a_{n+1}/a_n = mean b / mean m   (so a uniform ratio bound needs a GLOBAL
      statement "mean b <= 48 * mean m", NOT a per-set bound "max b <= 48")

Output: worst-case boundary max_n b, the slope law, the provable envelope
min(20n, mean b), the 48 = 8*6 arithmetic, and the effective-degree crossing.

Usage: py -3.12 rat1_bounds.py [maxn]
"""
import sys, time
from decimal import Decimal, getcontext

sys.path.insert(0, r"D:\Dev\math\research\v3\species\w7")
import ur1_enum as U1

getcontext().prec = 30
ROOT = U1.pc((0, 0, 0, 0), 0, 1)
A_CORPUS = [1, 20, 458, 11132, 280175, 7218768]   # n = 1..6 (UR4:51, ur4_enum_n6.log:19-21)


def boundary(A):
    NA = set(A)
    out = set()
    for p in A:
        for q in U1.neighbors(p):
            if q not in NA:
                out.add(q)
    return out


def bnd(A):
    return len(boundary(A))


def parents(B):
    S = set(B)
    cnt = 0
    for p in B:
        if p == ROOT:
            continue
        rest = S - {p}
        seen = {ROOT}
        stack = [ROOT]
        while stack:
            v = stack.pop()
            for q in U1.neighbors(v):
                if q in rest and q not in seen:
                    seen.add(q)
                    stack.append(q)
        if len(seen) == len(rest):
            cnt += 1
    return cnt


def nedges(A):
    S = list(A)
    nb = {p: set(U1.neighbors(p)) for p in S}
    return sum(1 for i in range(len(S)) for j in range(i + 1, len(S)) if S[j] in nb[S[i]])


# ---------- explicit families on the infinite lattice (exact boundary) ----------
def fam_chain(n):
    """straight chain in plane (0,1), consecutive tiles share exactly one link"""
    return [U1.pc((t, 0, 0, 0), 0, 1) for t in range(n)]


def fam_plate(k):
    """k x k plate in plane (0,1) (all tiles coplanar)"""
    return [U1.pc((i, j, 0, 0), 0, 1) for i in range(k) for j in range(k)]


def fam_tube(k):
    """k of the 6 plaquettes through the common link {0, e_0} (max k = 6)"""
    x = (0, 0, 0, 0)
    all6 = []
    for kk in (1, 2, 3):
        all6.append(U1.pc(x, 0, kk))
        all6.append(U1.pc(U1.sub(x, U1.DIRS[kk]), kk, 0))
    return all6[:k]


def fam_slab(m):
    """2 x m slab in plane (0,1)  (n = 2m tiles)"""
    return [U1.pc((i, 0, j, 0), 0, 1) for i in range(m) for j in range(2)]


def main():
    maxn = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    t0 = time.time()
    print("RAT1 - local-ratio lemma a_{n+1} <= 48 a_n, 4D plaquette graph, Delta=20")
    print("=" * 78)
    levels = U1.animals(maxn)
    lists = [list(L) for L in levels]
    a = [len(L) for L in lists]
    print(f"a_n = {a}   (corpus {A_CORPUS[:maxn]}) match={a == A_CORPUS[:maxn]}"
          f"   [{time.time()-t0:.1f}s]")

    # ---- 1. boundary statistics ----
    print("\n[1] boundary b(A) = #growth positions, exact over ALL animals")
    print(f"{'n':>2} {'a_n':>8} {'min b':>6} {'mean b':>9} {'max b':>6} {'max/n':>6} "
          f"{'sum b':>10} {'#b<=48':>7} {'#b>48':>7} {'maxb-(14n+6)':>13}")
    rows = []
    for n in range(1, maxn + 1):
        L = lists[n - 1]
        bs = [bnd(A) for A in L]
        s = sum(bs)
        rows.append(dict(n=n, a_n=len(L), bs=bs, sum_b=s, mean_b=s / len(bs),
                         min_b=min(bs), max_b=max(bs),
                         le48=sum(1 for x in bs if x <= 48), gt48=sum(1 for x in bs if x > 48)))
        r = rows[-1]
        print(f"{n:>2} {r['a_n']:>8} {r['min_b']:>6} {r['mean_b']:>9.4f} {r['max_b']:>6} "
              f"{r['max_b']/n:>6.2f} {r['sum_b']:>10} {r['le48']:>7} {r['gt48']:>7} "
              f"{r['max_b']-(14*n+6):>13}", flush=True)

    print("\n  law check: max b(n) - (14n+6) -- (14 = Delta - 6 = 20 - 2(d-1), d=4)")
    print("  max b/n  =", [f"{rows[i]['max_b']/(i+1):.3f}" for i in range(len(rows))])
    print("  min b/n  =", [f"{rows[i]['min_b']/(i+1):.3f}" for i in range(len(rows))])

    # argmax structure
    print("\n  argmax sets (n, #argmax, #edges of an argmax set, b):")
    for r in rows:
        n, bs = r["n"], r["bs"]
        mx = r["max_b"]
        idx = [i for i, x in enumerate(bs) if x == mx]
        A = lists[n - 1][idx[0]]
        print(f"    n={n}: #argmax={len(idx):>6}  edges={nedges(A)}  "
              f"b={mx}  (tree would have {n-1} edges)")

    # ---- 2. the identity and the ratio ----
    print("\n[2] identity sum_{A_n} b = sum_{A_{n+1}} m   =>   a_{n+1}/a_n = mean b / mean m")
    print(f"{'n':>2} {'sum_b(A_n)':>11} {'sum_m(A_n+1)':>13} {'mean m':>8} "
          f"{'ratio a_n+1/a_n':>16} {'mean_b/mean_m':>14}")
    for n in range(1, maxn):
        L2 = lists[n]                                   # A_{n+1}
        ms = [parents(B) for B in L2]
        sm = sum(ms)
        lhs = rows[n - 1]["sum_b"]
        rat = a[n] / a[n - 1]
        print(f"{n:>2} {lhs:>11} {sm:>13} {sm/len(ms):>8.5f} {rat:>16.4f} "
              f"{rows[n-1]['mean_b']/(sm/len(ms)):>14.4f}  "
              f"{'OK' if lhs == sm else 'MISMATCH'}", flush=True)
    if maxn >= 5:
        mm6 = rows[4]["sum_b"] / A_CORPUS[5]
        print(f"    mean m(A_6) = sum_b(A_5)/a_6 = {rows[4]['sum_b']}/{A_CORPUS[5]} = {mm6:.5f}"
              f"  [a_6 from UR4, [E]]")
        print(f"    ratio a_6/a_5 = {A_CORPUS[5]/A_CORPUS[4]:.4f} = mean b(A_5)/mean m(A_6)"
              f" = {rows[4]['mean_b']:.4f}/{mm6:.4f}")

    # ---- 3. ratios vs 48 ----
    print("\n[3] the lemma  a_{n+1} <= 48 a_n  on data")
    print(f"{'n':>2} {'a_n+1/a_n':>10} {'<=48?':>7} {'slack 48/ratio':>14} {'48*m_mean - mean_b':>19}")
    for n in range(1, maxn):
        rat = a[n] / a[n - 1]
        mm = (rows[n]["sum_b"] / A_CORPUS[n + 1]) if n + 1 < len(A_CORPUS) else float('nan')
        print(f"{n:>2} {rat:>10.4f} {str(rat <= 48):>7} {48/rat:>14.4f} "
              f"{48*mm - rows[n-1]['mean_b']:>19.4f}")

    # ---- 4. the provable local envelope ----
    print("\n[4] provable 'one-deletion' envelope:  a_{n+1} <= B_eff(n) * a_n")
    print("    per-set: b(A) <= Delta*n = 20n ;  with the identity: mean b(A_n) = B_eff")
    print(f"{'n':>2} {'20n':>5} {'mean b':>8} {'min(20n,mean b)':>16} {'n with 20n<=48':>15}")
    for r in rows:
        n = r["n"]
        print(f"{n:>2} {20*n:>5} {r['mean_b']:>8.4f} {min(20*n, r['mean_b']):>16.4f}")
    print("    20n <= 48  <=>  n <= 2 ;   mean b <= 48  <=>  n <= 3.37 (measured)")
    print(f"    min(20n, mean b) at n=5 = {min(100, rows[4]['mean_b']):.4f}  (>> 48)")
    print("    worst case  max_A b(A) <= 48  <=>  n <= 3      (max b = 20,34,48,62,76)")

    # ---- 5. 48 = 8 x 6 arithmetic ----
    print("\n[5] 48 = 8 x 6 : what are 8 and 6 ?")
    d = 4
    print(f"    d = {d}")
    print(f"    8 = 2d            = {2*d}       (site link-degree: links at a site)")
    print(f"    6 = 2(d-1)        = {2*(d-1)}       (plaquettes through one link)")
    print(f"    Delta = 4*(6-1)   = {4*(2*(d-1)-1)}      (plaquettes adjacent to one plaquette)")
    print(f"    => 48 = (links/site) x (plaquettes/link) = {(2*d)*(2*(d-1))}"
          f" : a SITE-link-plaquette incidence count")
    print(f"    boundary increment of the extremal family: Delta - 6 = {20-6}"
          f"  (6 = plaquettes through the shared link, already boundary via the parent)")
    print(f"    max b(3) = 48 = (Delta-6)*3 + 6 = {(20-6)*3+6}  -> 48 is the n=3 extreme, not a constant")
    print(f"    other decompositions: 48 = 2*Delta+8 = {2*20+8} ; 48/Delta = {Decimal(48)/Decimal(20)}")

    # ---- 6. effective degree / tree envelope ----
    print("\n[6] tree envelope  L_t(D) = D^D/(D-1)^(D-1)  (proven: lambda_* <= L_t(Delta) = e*Delta)")
    for dd in (17, 18, 19, 20):
        Dd = Decimal(dd)
        print(f"    L_t({dd}) = {Dd**dd / (Dd-1)**(dd-1)}")
    lo, hi = Decimal(18), Decimal(19)
    for _ in range(80):
        mid = (lo + hi) / 2
        if mid**mid / (mid - 1)**(mid - 1) < 48:
            lo = mid
        else:
            hi = mid
    print(f"    D_eff with L_t(D_eff) = 48 :  D_eff = {lo}")
    e = Decimal(1).exp()
    print(f"    e*Delta = {e*20} ;  e*Delta/48 = {e*20/48}")

    # ---- 7. explicit families: boundary grows like n ----
    print("\n[7] explicit families on the infinite lattice (exact boundary, no torus)")
    print(f"{'n':>3} {'chain b':>8} {'b/n':>6} | {'plate kxk b':>11} {'b/n':>6} | "
          f"{'tube b':>7} | {'2 x n slab':>10}")
    plates = {k * k: bnd(fam_plate(k)) for k in (1, 2, 3, 4)}
    for n in range(1, 13):
        bc = bnd(fam_chain(n))
        bt = bnd(fam_tube(n)) if n <= 6 else 0
        bp = plates.get(n, None)
        bs2 = bnd(fam_slab(n // 2)) if n % 2 == 0 and n // 2 >= 1 else None
        print(f"{n:>3} {bc:>8} {bc/n:>6.2f} | "
              f"{(str(bp) if bp else '-'):>11} {(f'{bp/n:.2f}' if bp else '-'):>6} | "
              f"{bt:>7} | {(str(bs2) if bs2 else '-'):>10}")
    print("    plate 4x4 (n=16):", plates[16], " b/n =", plates[16] / 16)
    print("    => sup_n max_A b(A) = infinity  (both families are linear with slope ~14-16)")

    print(f"\ntotal {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
