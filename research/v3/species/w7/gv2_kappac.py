"""
GV2 -- kappa_c as a functional of the Ursell shell masses {m_n}.

Hypothesis under test (GU4-NECESSARY:53 / TT2:16-19):
    sign(phi^U(Gamma)) = (-1)^(|Gamma|-1) * |phi^U(Gamma)|      ("sign lives on n, not on shape")

If TRUE then, for a fixed boundary dC0,
    Z_signed(n) = | sum_m s_m A_m(n) | ,   Z_abs(n) = sum_m A_m(n) ,   A_m >= 0
i.e. kappa_c is the alternating partial sum of the shell masses => a functional of {m_n}
=> the SAME computation as lambda_Ur (gv1 / ur1_enum).

This script
  (1) re-derives the shell masses m_n on the 4D plaquette graph (Delta=20), n<=5, from scratch,
      and verifies sign(phi^U) == (-1)^(n-1) for BOTH bond-activity conventions;
  (2) imports the n=6 layer from ur1_enum_n6.json (exact integers) with an internal consistency check;
  (3) tabulates kappa_c(n) in two normalisations + a rigorous reverse-triangle lower bound;
  (4) u-scan: shows kappa_c(n) at small n is dominated by accidental low-order cancellation;
  (5) geometric-tail model  R = |1-rho|/(1+rho)  vs exact R_n.

Conventions (declared, not hidden):
  u   = beta/N^2   (repo convention: SU(2), beta=2.5 -> u=0.625; see TT1:60 / L1SIGN:54)
  phi^U(Gamma) = sum_{A subset E, (Gamma,A) connected} (-1)^|A| prod_{e in A} f_e
  convention A: f_e = 1    => phi^U = u^n * mu(Gamma),  mu integer (= c_Gamma of UR6/GU4 4.5)
  convention B: f_e = u    => phi^U = u^n * C_Gamma(u)
Output: gv2_kappac.json
"""
import json
import math
import time

DIRS = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
PLANES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
PMAP = {p: i for i, p in enumerate(PLANES)}
OFF = 20
NMAX = 5                      # re-enumerated here (cheap: 11 s in ur1_enum.py)
U_REPO = 0.625                # beta=2.5, N=2: u = beta/N^2
U_SU3 = 2.5 / 9.0             # beta=2.5, N=3: u = beta/N^2


def sc(x):
    return (((x[0] + OFF) << 15) | ((x[1] + OFF) << 10)
            | ((x[2] + OFF) << 5) | (x[3] + OFF))


def pc(x, mu, nu):
    if mu > nu:
        mu, nu = nu, mu
    return (sc(x) << 3) | PMAP[(mu, nu)]


def dec(c):
    s, pl = c >> 3, c & 7
    x = (((s >> 15) & 31) - OFF, ((s >> 10) & 31) - OFF,
         ((s >> 5) & 31) - OFF, (s & 31) - OFF)
    return x, PLANES[pl]


def add(x, d):
    return (x[0] + d[0], x[1] + d[1], x[2] + d[2], x[3] + d[3])


def sub(x, d):
    return (x[0] - d[0], x[1] - d[1], x[2] - d[2], x[3] - d[3])


_NB = {}


def neighbors(c):
    r = _NB.get(c)
    if r is not None:
        return r
    x, (mu, nu) = dec(c)
    links = [(x, mu), (add(x, DIRS[mu]), nu), (add(x, DIRS[nu]), mu), (x, nu)]
    out = set()
    for (y, d) in links:
        for k in range(4):
            if k == d:
                continue
            out.add(pc(y, d, k))
            out.add(pc(sub(y, DIRS[k]), k, d))
    out.discard(c)
    _NB[c] = tuple(sorted(out))
    return _NB[c]


def animals(maxn, root=None):
    """connected n-plaquette sets containing the root plaquette p0 (plane (0,1), x=0)."""
    if root is None:
        root = pc((0, 0, 0, 0), 0, 1)
    levels = [{(root,)}]
    for _ in range(1, maxn):
        prev, nxt = levels[-1], set()
        for S in prev:
            base = set(S)
            for p in S:
                for q in neighbors(p):
                    if q not in base:
                        nxt.add(tuple(sorted(S + (q,))))
        levels.append(nxt)
    return levels


def adjacency(S):
    n = len(S)
    nbrsets = [set(neighbors(p)) for p in S]
    adj = [0] * n
    mE = 0
    for i in range(n):
        for j in range(i + 1, n):
            if S[j] in nbrsets[i]:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
                mE += 1
    return adj, mE


def conn_sum(adj, n, x):
    """sum over connected spanning A subset E of (-x)^|A|  (x=1 -> exact int)."""
    m = 1 << n
    ec = [0] * m
    for T in range(1, m):
        low = T & -T
        i = low.bit_length() - 1
        R = T ^ low
        ec[T] = ec[R] + bin(adj[i] & R).count("1")
    g = [((1 - x) ** k) for k in ec]
    f = [0] * m
    for T in range(1, m):
        low = T & -T
        rest = T ^ low
        tot = g[T]
        Uu = rest
        while True:
            if Uu != rest:
                P = low | Uu
                tot -= f[P] * g[T ^ P]
            if Uu == 0:
                break
            Uu = (Uu - 1) & rest
        f[T] = tot
    return f[m - 1], ec[m - 1]


def main():
    t0 = time.time()
    print("=" * 78)
    print("GV2 step 1 -- shell masses and the SIGN of phi^U  (4D plaquette graph, Delta=20)")
    print("=" * 78)
    levels = animals(NMAX)
    rows = []
    for n in range(1, NMAX + 1):
        t1 = time.time()
        nB = 0                    # sum |C_Gamma(u)|   (convention B magnitudes)
        B = 0                     # sum |mu|           (convention A magnitudes)
        violA = violB = 0
        sgn = (-1) ** (n - 1)
        histA, histB = {}, {}
        for S in levels[n - 1]:
            adj, mE = adjacency(S)
            mu, _ = conn_sum(adj, n, 1)
            cu, _ = conn_sum(adj, n, U_REPO)
            B += abs(mu)
            nB += abs(cu)
            sa = 1 if mu > 0 else (-1 if mu < 0 else 0)
            if sa != sgn:
                violA += 1
            sb = 1 if cu > 0 else (-1 if cu < 0 else 0)
            if sb != sgn:
                violB += 1
            histA[mu] = histA.get(mu, 0) + 1
            histB[round(cu, 12)] = histB.get(round(cu, 12), 0) + 1
        rows.append(dict(n=n, count=len(levels[n - 1]), B=B, sum_abs_Cu=nB,
                         sign_viol_A=violA, sign_viol_B=violB, sgn=sgn,
                         histA={str(k): v for k, v in sorted(histA.items())},
                         histB={str(k): v for k, v in sorted(histB.items())}))
        print(f"n={n}: count={len(levels[n-1]):>7d}  B_n=sum|mu|={B:>9d}  "
              f"sum|C(u)|={nB:>12.4f}  sign-viol A={violA:>6d} B={violB:>7d}  "
              f"[{time.time()-t1:.1f}s]")

    # ---- n=6 layer: exact integers from the pre-existing run (ur1_enum_n6.json) ----
    n6 = json.load(open("ur1_enum_n6.json"))
    chk_cnt = sum(n6["hist"].values())
    chk_B = sum(abs(int(k)) * v for k, v in n6["hist"].items())
    allneg = all(int(k) < 0 for k in n6["hist"])
    print(f"\nn=6 (from ur1_enum_n6.json): count={n6['count']} (hist sum {chk_cnt}, "
          f"{'OK' if chk_cnt == n6['count'] else 'MISMATCH'})")
    print(f"     B_6 = {n6['B']} (hist |sum| {chk_B}, "
          f"{'OK' if chk_B == n6['B'] else 'MISMATCH'});  all mu<0 for n=6 even: {allneg}")
    rows.append(dict(n=6, count=n6["count"], B=n6["B"], sum_abs_Cu=None,
                     sign_viol_A=0 if allneg else -1, sign_viol_B=None, sgn=(-1) ** 5,
                     source="ur1_enum_n6.json"))

    B_seq = [r["B"] for r in rows]
    print(f"\nB_n (n=1..6) = {B_seq}")
    print(f"lam_n = B_n^(1/n) = " + "  ".join(f"{b ** (1.0 / i):.6f}"
                                              for i, b in enumerate(B_seq, start=1)))

    # ---------------------------------------------------------------- step 2
    print("\n" + "=" * 78)
    print("GV2 step 2 -- kappa_c(n) := (|sum_k<=n s_k m_k| / sum_k<=n m_k)^(1/n),  m_k = u^k B_k")
    print("=" * 78)

    def kappa_table(u, label):
        m = [u ** k * B_seq[k - 1] for k in range(1, 7)]
        S = 0.0
        Z = 0.0
        out = []
        for k, mk in enumerate(m, start=1):
            S += (-1) ** (k - 1) * mk
            Z += mk
            R = abs(S) / Z
            # rigorous reverse-triangle lower bound on R (needs m_k to dominate):
            Rlb = max(0.0, (2.0 * mk - Z) / Z)
            out.append(dict(n=k, m_n=mk, S=abs(S), Z=Z, R=R,
                            kappa=R ** (1.0 / k), kappa_lb=Rlb ** (1.0 / k),
                            req80=0.80 ** k, req81=0.81 ** k,
                            ratio_to_req=R / (0.80 ** k)))
        print(f"\n--- {label}: u={u:.10f} ---")
        print(" n |        m_n        |      Z_abs     |  R=Sz/Zabs | kappa_c(n) | LB   | 0.80^n req | R/req")
        for r in out:
            print(f" {r['n']} | {r['m_n']:18.6e} | {r['Z']:14.6e} | {r['R']:10.8f} | "
                  f"{r['kappa']:10.8f} | {r['kappa_lb']:.4f} | {r['req80']:10.3e} | "
                  f"{r['ratio_to_req']:9.3e}")
        return out

    tab_su3 = kappa_table(U_SU3, "SU(3), beta=2.5, u=beta/N^2")
    tab_su2 = kappa_table(U_REPO, "SU(2), beta=2.5, u=beta/N^2 (repo convention, TT1:60)")

    # ---------------------------------------------------------------- step 3
    print("\n" + "=" * 78)
    print("GV2 step 3 -- u-scan: kappa_c(n) at small n is an ACCIDENTAL-CANCELLATION artifact")
    print("=" * 78)
    us = [0.02 * j for j in range(1, 36)]
    print("   u     | kappa_c(2) | kappa_c(3) | kappa_c(4) | kappa_c(5) | kappa_c(6)")
    uscan = []
    for u in us:
        m = [u ** k * B_seq[k - 1] for k in range(1, 7)]
        S = 0.0
        Z = 0.0
        row = []
        for k, mk in enumerate(m, start=1):
            S += (-1) ** (k - 1) * mk
            Z += mk
            row.append((abs(S) / Z) ** (1.0 / k))
        uscan.append(dict(u=u, k=row))
        print(f" {u:5.2f}  |  {row[1]:8.5f}  |  {row[2]:8.5f}  |  {row[3]:8.5f}  |"
              f"  {row[4]:8.5f}  |  {row[5]:8.5f}")
    k6 = [r["k"][5] for r in uscan]
    k4 = [r["k"][3] for r in uscan]
    print(f"\nmin over u-grid of kappa_c(4) = {min(k4):.6f} (at u={us[k4.index(min(k4))]:.2f})")
    print(f"min over u-grid of kappa_c(6) = {min(k6):.6f} (at u={us[k6.index(min(k6))]:.2f})")
    print("u*(n=2 cancellation point) = 1/B_2 = %.10f" % (1.0 / B_seq[1]))

    # ---------------------------------------------------------------- step 4
    print("\n" + "=" * 78)
    print("GV4 step 4 -- geometric-tail model: R -> |1-rho|/(1+rho),  rho = m_n/m_(n-1)")
    print("=" * 78)
    print(" n | rho_n = m_n/m_(n-1) | exact R_n | |1-rho|/(1+rho) | (1-rho)/(1+rho) if rho<1")
    geo = []
    for u, lbl in ((U_SU3, "SU3"), (U_REPO, "SU2")):
        m = [u ** k * B_seq[k - 1] for k in range(1, 7)]
        S = 0.0
        Z = 0.0
        for k, mk in enumerate(m, start=1):
            S += (-1) ** (k - 1) * mk
            Z += mk
            rho = mk / m[k - 2] if k > 1 else float("nan")
            mod = abs(1 - rho) / (1 + rho) if k > 1 else float("nan")
            geo.append(dict(u=u, n=k, rho=rho, R=abs(S) / Z, model=mod))
            print(f" {k} |      {rho:14.8f} | {abs(S)/Z:9.6f} |     {mod:12.8f} |  "
                  f"{lbl}")
        print()

    # ---------------------------------------------------------------- step 5
    print("=" * 78)
    print("GV2 step 5 -- verdict inputs")
    print("=" * 78)
    print(f"kappa_c(6), u=beta/9 (SU3,beta=2.5) = {tab_su3[-1]['kappa']:.8f}"
          f"   [LB {tab_su3[-1]['kappa_lb']:.6f}]")
    print(f"kappa_c(6), u=0.625  (SU2,beta=2.5) = {tab_su2[-1]['kappa']:.8f}"
          f"   [LB {tab_su2[-1]['kappa_lb']:.6f}]")
    print("criterion: kappa_c >= 0.81 => CPB-W and CPB-S both CLOSE (BianchiLemma.md:24)")
    print("required for a live route: R_n <= 0.80^n, i.e. kappa_c(n) <= 0.80")

    json.dump(dict(rows=rows, B_seq=B_seq, tab_su3=tab_su3, tab_su2=tab_su2,
                   uscan=uscan, geo=geo,
                   thresholds=dict(close=0.81, live=0.80)),
              open("gv2_kappac.json", "w"), indent=1)
    print(f"\njson -> gv2_kappac.json     total {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
