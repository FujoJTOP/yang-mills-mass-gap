"""
GL3-HIGHORDER numerics (v4).  Exact rational arithmetic, no MC.

loop   : closed step sequence on Z^4 based at 0, reduced (no immediate backtrack).
         step (mu,-1) at p is U_mu(p-e_mu)^dag  =>  field -A_mu(p-e_mu).
w_mu   : sum of (+-1) z^{field base},  z_mu = exp(i p_mu).
ident  : sum_mu w_mu (z_mu-1) = 0                       [checked: 0 violations]
orbit  : point group x cyclic rotation; summed over DISTINCT images (each loop and its
         reverse both occur => a uniform factor 2, cancelled by the H_1 normalization).
axis   : p=(iE,0,0,0), x=aE, H(x) = sum_g sum_{i=1,2,3} w_{gi}(z0) w_{gi}(1/z0).
         kappa_j := H_j/H_1,  H(x) = -x^2 (1 + sum_{j>=2} kappa_j x^{2j}).
action : normalized coeffs d_O (sum d = 1), conditions sum d_O kappa_j^O = 0.
free   : k = 2(j_first-1);  Delta = m - kappa_j a^{2(j-1)} m^{2(j-1)+1}/2 + ...
checks : plaquette kappa_2 = 1/12 ; rectangle kappa_2 = 5/24 ;
         LW (d_plaq=5/3, d_rect=-2/3) -> kappa_2 = 0, kappa_3 = -1/90
         => Delta = m + a^4 m^5/180   (= GH3/GK3 known value).
"""
import itertools, math
from fractions import Fraction
from collections import defaultdict

PERMS = list(itertools.permutations(range(4)))
SIGNS = [s for s in itertools.product((1, -1), repeat=4)]
JMAX = 8


def gen_loops(maxlen):
    res = []

    def rec(steps, pos):
        rem = maxlen - len(steps)
        if rem == 0:
            if pos == (0, 0, 0, 0):
                res.append(tuple(steps))
            return
        if sum(abs(c) for c in pos) > rem:
            return
        for mu in range(4):
            for s in (1, -1):
                if steps:
                    m0, s0 = steps[-1]
                    if mu == m0 and s == -s0:
                        continue
                p = list(pos); p[mu] += s
                steps.append((mu, s)); rec(steps, tuple(p)); steps.pop()

    rec([], (0, 0, 0, 0))
    return res


def cyc_canon(st):
    n = len(st)
    return min(st[k:] + st[:k] for k in range(n))


def cheap_key(steps):
    cnt = [[0, 0] for _ in range(4)]
    for mu, s in steps:
        cnt[mu][0 if s == 1 else 1] += 1
    best = None
    for perm in PERMS:
        for sg in SIGNS:
            new = [None] * 4
            for mu in range(4):
                a, b = cnt[mu]
                if sg[mu] == -1:
                    a, b = b, a
                new[perm[mu]] = (a, b)
            k = tuple(sorted(new))
            if best is None or k < best:
                best = k
    return best


def group_imgs(st):
    return [cyc_canon(tuple((perm[mu], s * sg[mu]) for (mu, s) in st)) for perm in PERMS for sg in SIGNS]


def links_of(steps):
    pos = [0, 0, 0, 0]; out = []
    for (mu, s) in steps:
        out.append((mu, s, tuple(pos))); pos[mu] += s
    assert pos == [0, 0, 0, 0]
    return out


def wvec(steps):
    w = [defaultdict(int) for _ in range(4)]
    for (mu, s, pos) in links_of(steps):
        b = list(pos)
        if s == -1:
            b[mu] -= 1
        w[mu][tuple(b)] += s
    return [dict(d) for d in w]


def transv_res(w):
    acc = defaultdict(int)
    for mu in range(4):
        for a, c in w[mu].items():
            acc[a] -= c
            a2 = list(a); a2[mu] += 1
            acc[tuple(a2)] += c
    return {k: v for k, v in acc.items() if v}


def orbit_profile(rep):
    imgs = sorted(set(group_imgs(rep)))
    d = defaultdict(int)
    for g in imgs:
        w = wvec(g)
        for i in (1, 2, 3):
            for a, ca in w[i].items():
                for b, cb in w[i].items():
                    if a[1] or a[2] or a[3] or b[1] or b[2] or b[3]:
                        continue
                    d[a[0] - b[0]] += ca * cb
    return {k: v for k, v in d.items() if v}, len(imgs)


def taylor(d, nmax=JMAX):
    out = {}
    for j in range(1, nmax + 1):
        s = Fraction(0)
        for m, c in d.items():
            if m:
                s += Fraction(2 * c) * Fraction(m ** (2 * j), math.factorial(2 * j))
        out[j] = s
    return out


def rank_q(rows):
    m = [list(r) for r in rows]
    nr = len(m); nc = len(m[0]) if nr else 0
    rk = 0
    for c in range(nc):
        piv = next((r for r in range(rk, nr) if m[r][c] != 0), None)
        if piv is None:
            continue
        m[rk], m[piv] = m[piv], m[rk]
        pv = m[rk][c]; m[rk] = [x / pv for x in m[rk]]
        for r in range(nr):
            if r != rk and m[r][c] != 0:
                f = m[r][c]; m[r] = [x - f * y for x, y in zip(m[r], m[rk])]
        rk += 1
    return rk


def solve_q(A, b):
    nr = len(A); nc = len(A[0])
    M = [list(A[i]) + [b[i]] for i in range(nr)]
    piv_cols, r = [], 0
    for c in range(nc):
        piv = next((rr for rr in range(r, nr) if M[rr][c] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]; M[r] = [x / pv for x in M[r]]
        for rr in range(nr):
            if rr != r and M[rr][c] != 0:
                f = M[rr][c]; M[rr] = [x - f * y for x, y in zip(M[rr], M[r])]
        piv_cols.append(c); r += 1
        if r == nr:
            break
    for rr in range(r, nr):
        if all(x == 0 for x in M[rr][:nc]) and M[rr][nc] != 0:
            return None, 0
    x = [Fraction(0)] * nc
    for i, c in enumerate(piv_cols):
        x[c] = M[i][nc]
    return x, nc - len(piv_cols)


# ============================================================ A. orbits
print("=" * 78)
print("A. loop orbits (reduced, based at 0), perimeter <= 8")
print("=" * 78)
orbits = []
for L in (4, 6, 8):
    raw = gen_loops(L)
    bucket = defaultdict(list)
    for st in raw:
        bucket[cheap_key(st)].append(st)
    reps = []
    for _k, lst in bucket.items():
        seen = set()
        for st in lst:
            c = cyc_canon(st)
            if c in seen:
                continue
            imgs = set(group_imgs(st))
            seen |= imgs
            reps.append(min(imgs))
    reps.sort()
    bad = sum(1 for rep in reps for g in set(group_imgs(rep)) if transv_res(wvec(g)))
    print(f"perimeter {L:2d}: raw = {len(raw):6d} ; orbits = {len(reps):3d} ; "
          f"transversality violations = {bad}")
    for rep in reps:
        d, nimg = orbit_profile(rep)
        t = taylor(d)
        if t[1] == 0:
            orbits.append(dict(L=L, rep=rep, d=d, kap=None, mext=max((abs(m) for m in d), default=0)))
            continue
        kap = {j: t[j] / t[1] for j in range(2, JMAX + 1)}
        orbits.append(dict(L=L, rep=rep, d=d, kap=kap, mext=max((abs(m) for m in d), default=0)))

print("\naxis profile per orbit  (kappa_j, +- = sign of H_1):")
for i, o in enumerate(orbits):
    tag = "  (axis-inactive)" if o["kap"] is None else ""
    print(f"[{i:2d}] P={o['L']} mext={o['mext']:d}{tag}")
    if o["kap"]:
        print("      kappa2..6 = " + " ".join(f"{str(o['kap'][j]):>11s}" for j in range(2, 7)))

# ---- checks
P = [o for o in orbits if o["L"] == 4][0]
REC = [o for o in orbits if o["L"] == 6 and o["mext"] == 2]
print("\nCHECK plaquette kappa_2 =", P["kap"][2], " (expect 1/12)")
for r in REC:
    print("CHECK mext=2 orbit kappa_2 =", r["kap"][2], " kappa_3 =", r["kap"][3], " (rect expect 5/24, 17/720)")
if REC:
    r = REC[0]
    d_p, d_r = Fraction(5, 3), Fraction(-2, 3)
    kv = {j: d_p * P["kap"][j] + d_r * r["kap"][j] for j in range(2, JMAX + 1)}
    print("CHECK LW (d=5/3,-2/3): kappa_2 =", kv[2], " kappa_3 =", kv[3],
          " => Delta/m - 1 = ", -kv[3] / 2, "* a^4 m^4  (expect +1/180)")

# ============================================================ B. conditions
print()
print("=" * 78)
print("B. axis conditions per perimeter budget:  sum d = 1, sum d kappa_j = 0")
print("=" * 78)
for Lmax in (4, 6, 8):
    idx = [i for i, o in enumerate(orbits) if o["L"] <= Lmax and o["kap"] is not None]
    N = len(idx)
    # dedupe identical profiles
    uniq, seenp = [], set()
    for i in idx:
        key = tuple(orbits[i]["kap"][j] for j in range(2, JMAX + 1))
        if key not in seenp:
            seenp.add(key); uniq.append(i)
    rows = [[Fraction(1)] * len(uniq)]
    rhs = [Fraction(1)]
    killed = 1
    for j in range(2, JMAX + 1):
        rows.append([orbits[i]["kap"][j] for i in uniq]); rhs.append(Fraction(0))
        sol, nd = solve_q(rows, rhs)
        if sol is None:
            rows.pop(); rhs.pop(); break
        killed = j
    sol, nd = solve_q(rows, rhs)
    kv = [sum(sol[k] * orbits[i]["kap"][j] for k, i in enumerate(uniq)) for j in range(2, JMAX + 1)]
    first = 2 + next(t for t, v in enumerate(kv) if v != 0)
    print(f"perimeter<={Lmax}: active orbits N={N:3d}, distinct profiles M={len(uniq):2d}; "
          f"conditions imposed = {killed}; nullspace = {nd}")
    print(f"   -> k = {2*(first-1)}   (first nonzero kappa_{first} = {kv[first-2]});  "
          f"Delta = m {'+' if -kv[first-2] < 0 else '-'} {abs(-kv[first-2] / 2)} a^{2*(first-1)} m^{2*(first-1)+1}")
    print(f"   d (on distinct profiles) = {[str(x) for x in sol]}")
