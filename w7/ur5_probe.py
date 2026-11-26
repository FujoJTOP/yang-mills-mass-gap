"""ur5_probe.py -- is the S_n class-function Weingarten ansatz legit for SU(2)?

Ansatz A(n):  Int prod_p U_{a_p b_p} prod_p Ubar_{c_p d_p}
            = sum_{s,t in S_n} Wg(s^-1 t) prod_p d_{a_p d_{s(p)}} prod_p d_{b_p c_{t(p)}}
with Wg a function of the class only (p(n) unknowns).

Test: build the exact linear system over many index configurations (indices in {0,1}, N=2,
exact via the S^3 hafnian) and compare rank(M) with rank([M|b]):
  rank(M)==rank([M|b])==p(n)  -> ansatz holds, Wg = unique solution
  rank(M) < p(n) but consistent -> ansatz under-determined by these probes
  rank([M|b]) > rank(M)        -> ansatz FALSE (extra invariant tensors: epsilon / Brauer sector)
"""
import itertools, random, sys
from fractions import Fraction as F
from sympy import Matrix, Rational
import ur5_wg as W


def build(n, ncfg, seed=3):
    perms = list(itertools.permutations(range(n)))
    ctlist, mult = W.classes_of(n)
    cidx = {ct: i for i, ct in enumerate(ctlist)}
    inv = {p: W.inverse(p) for p in perms}
    cache = {}

    def classrow(s):
        if s not in cache:
            si = inv[s]
            cache[s] = [cidx[W.cycle_type(W.compose(si, t))] for t in perms]
        return cache[s]

    rng = random.Random(seed)
    rows, rhs = [], []
    for _ in range(ncfg):
        ab = [rng.randint(0, 1) for _ in range(n)]
        cd = [rng.randint(0, 1) for _ in range(n)]
        ab2 = [rng.randint(0, 1) for _ in range(n)]
        cd2 = [rng.randint(0, 1) for _ in range(n)]
        vecs = [W.mvec(ab[p], cd[p]) for p in range(n)] + [W.mvec(ab2[p], cd2[p], True) for p in range(n)]
        lhs = W.integ_exact(vecs)
        row = [0] * len(ctlist)
        for s in perms:
            if any(ab[p] != ab2[s[p]] for p in range(n)):
                continue
            cr = classrow(s)
            for ti, t in enumerate(perms):
                if all(cd[p] == cd2[t[p]] for p in range(n)):
                    row[cr[ti]] += 1
        rows.append(row)
        rhs.append(lhs)
    return ctlist, rows, rhs


for n in (1, 2, 3, 4):
    ncfg = {1: 20, 2: 200, 3: 2000, 4: 4000}[n]
    ctlist, rows, rhs = build(n, ncfg)
    M = Matrix(rows)
    Maug = M.row_join(Matrix(rhs))
    r, raug = M.rank(), Maug.rank()
    print("n=%d  p(n)=%d  #eq=%d   rank(M)=%d  rank([M|b])=%d  -> %s" %
          (n, len(ctlist), len(rows), r, raug,
           "ANSATZ-OK (unique)" if (r == raug == len(ctlist)) else
           ("ANSATZ-UNDERDETERMINED" if r == raug else "ANSATZ-FALSE (Brauer/epsilon sector)")))
    if r == raug and r == len(ctlist):
        sol, params = M.gauss_jordan_solve(Matrix(rhs))
        print("    solution:", {W.cname(ct): sol[i] for i, ct in enumerate(ctlist)})
    elif r == raug:
        sol, params = M.gauss_jordan_solve(Matrix(rhs))
        print("    one solution:", {W.cname(ct): sol[i] for i, ct in enumerate(ctlist)})
    sys.stdout.flush()
