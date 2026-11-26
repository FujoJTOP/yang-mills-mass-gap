"""C07 hand-check appendix: n=3,4 identity witnesses, term counts, unit-activity ratios.
Reproduces exactly the numbers quoted in C07-AUDIT.md §1.
"""
import itertools
import sympy as sp


def conn(n, ed):
    par = list(range(n))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for i, j in ed:
        a, b = f(i), f(j)
        if a != b:
            par[a] = b
    return len({f(i) for i in range(n)}) == 1


def tau(n, ed):
    L = sp.zeros(n, n)
    for i, j in ed:
        L[i, i] += 1
        L[j, j] += 1
        L[i, j] -= 1
        L[j, i] -= 1
    return int(L[:n - 1, :n - 1].det())


for n in (3, 4):
    Eall = list(itertools.combinations(range(n), 2))
    sym = sp.symbols(f"a0:{len(Eall)}")
    am = dict(zip(Eall, sym))
    LHS = sp.Integer(0)
    RHS = sp.Integer(0)
    RHS2 = sp.Integer(0)
    print("=" * 70)
    print(f"n={n}, edges indexed {list(enumerate(Eall))}")
    rows = []
    for m in range(1 << len(Eall)):
        ed = [Eall[i] for i in range(len(Eall)) if (m >> i) & 1]
        if not conn(n, ed):
            continue
        t = tau(n, ed)
        mono = sp.prod([am[e] for e in ed])
        LHS += t * mono
        rows.append((ed, t, len(ed)))
        if len(ed) == n - 1:
            RHS += mono
            RHS2 += mono * sp.prod([1 + am[e] for e in Eall if e not in ed])
    for ed, t, m in sorted(rows, key=lambda r: (r[2], r[0])):
        print(f"   edges={[Eall.index(e) for e in ed]}  |E|={m}  tau={t}")
    d = sp.expand(LHS - RHS)
    terms = sp.Add.make_args(d)
    print("   LHS          =", sp.expand(LHS))
    print("   RHS (trees)  =", sp.expand(RHS))
    print("   LHS - RHS    =", d)
    print(f"   #terms of difference = {len(terms)}, "
          f"degree range [{min(len(t.free_symbols) for t in terms)},"
          f"{max(len(t.free_symbols) for t in terms)}]")
    print("   LHS - Sum_T a^T prod(1+a_e) =", sp.expand(LHS - RHS2), " (0 => corrected identity)")
    # unit activities
    u = {e: 1 for e in Eall}
    L1 = sum(tau(n, ed) for ed in [e_ for e_ in ([Eall[i] for i in range(len(Eall)) if (m >> i) & 1] for m in range(1 << len(Eall)))] if ed and conn(n, ed))
    R1 = sum(1 for m in range(1 << len(Eall)) if bin(m).count("1") == n - 1 and conn(n, [Eall[i] for i in range(len(Eall)) if (m >> i) & 1]))
    print(f"   unit activity: LHS={L1}  RHS={R1}  ratio={L1//R1}  "
          f"2^C(n-1,2)={2**sp.binomial(n - 1, 2)}")
print("=" * 70)
print("hand-check done.")
