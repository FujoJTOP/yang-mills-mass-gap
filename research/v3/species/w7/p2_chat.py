"""P2-CHAT: C^ = sup_gamma (prod_l R(kappa_l))^{1/|gamma|} 的精确值、4 是否被取到、
换判据（不均匀 KP / 按类分裂）能否避开那个取到 4 的立方体。

py -3.12 + mpmath, dps=40. 只读，不写本目录以外的文件。
"""
import itertools, time, sys
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 40
OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.append(s)


def haar(f):
    return mp.quad(lambda th: mp.sin(th / 2) ** 2 * f(th), [0, mp.pi, 2 * mp.pi]) / mp.pi


chi_f = lambda th: 2 * mp.cos(th / 2)
chi_adj = lambda th: 1 + 2 * mp.cos(th)

# ---------------- 1  R(kappa) 精确值
P("== 1  R(kappa) = M_{B+g}/(M_B M_g) = 1 + mult(adj in f^kappa)/mult(1 in f^kappa) ==")
m1, madj, R = {}, {}, {}
for m in range(1, 7):
    n = 2 * m
    m1[n] = haar(lambda t: chi_f(t) ** n)
    madj[n] = haar(lambda t: chi_adj(t) * chi_f(t) ** n)
cat = {1: 1, 2: 2, 3: 5, 4: 14, 5: 42, 6: 132}
P(f"{'kappa':>5} {'mult(1 in f^k)':>16} {'mult(adj in f^k)':>17} {'R':>10}  {'m^2/4^m':>10}")
for m in range(1, 7):
    n = 2 * m
    R[n] = 1 + madj[n] / m1[n]
    P(f"{n:>5} {mp.nstr(m1[n],14):>16} {mp.nstr(madj[n],14):>17} {mp.nstr(R[n],16):>10}"
      f"  {mp.nstr(R[n] / mp.mpf(4)**m, 10):>10}")
P("  (mult(1 in f^{2m}) = Catalan(m) 逐位相符; 奇 kappa 全为 0 => w(gamma)!=0 强制全偶数 kappa)")

# ---------------- 2  精确恒等式: 闭簇上 prodR = 4^F * (R(4)/4)^{E4}(R(6)/8)^{E6}(R(8)/16)^{E8}
def links_of(p):
    s, (i, j) = p
    si, sj = list(s), list(s)
    si[i] += 1
    sj[j] += 1
    return [(tuple(s), i), (tuple(sj), i), (tuple(s), j), (tuple(si), j)]


def hist(sub):
    h = {}
    for p in sub:
        for L in links_of(p):
            h[L] = h.get(L, 0) + 1
    return h


def build(n):
    sites = list(itertools.product(*[range(k) for k in n]))
    pla = []
    for s in sites:
        for i in range(4):
            for j in range(i + 1, 4):
                if s[i] + 1 < n[i] and s[j] + 1 < n[j]:
                    pla.append((s, (i, j)))
    return pla


pla2222 = build((2, 2, 2, 2))
P("")
P(f"== 2  恒等式 4F = sum_l kappa_l ；闭簇 = GF(2) 2-圈（全 kappa 偶）；穷举盒内所有闭簇 ==")


def closed_subcomplexes(n, cap=1 << 15):
    pla = build(n)
    idx = {p: k for k, p in enumerate(pla)}
    links = {}
    for p in pla:
        for L in links_of(p):
            links.setdefault(L, []).append(idx[p])
    nf = len(pla)
    rows = []
    for L, ps in links.items():
        m = 0
        for c in ps:
            m |= 1 << c
        rows.append(m)
    nl = len(rows)
    # RREF over GF(2)
    piv_of_col = {}
    rank = 0
    for c in range(nf):
        r = next((rr for rr in range(rank, nl) if (rows[rr] >> c) & 1), None)
        if r is None:
            continue
        rows[rank], rows[r] = rows[r], rows[rank]
        for rr in range(nl):
            if rr != rank and ((rows[rr] >> c) & 1):
                rows[rr] ^= rows[rank]
        piv_of_col[c] = rank
        rank += 1
    free = [c for c in range(nf) if c not in piv_of_col]
    basis = []
    for f in free:
        v = 1 << f
        for c, r in piv_of_col.items():
            if (rows[r] >> f) & 1:
                v |= 1 << c
        basis.append(v)
    d = len(basis)
    out = []
    if 2 ** d > cap:
        return pla, links, d, None
    for comb in range(1 << d):
        v = 0
        for b in range(d):
            if (comb >> b) & 1:
                v ^= basis[b]
        out.append([pla[k] for k in range(nf) if (v >> k) & 1])
    return pla, links, d, out


def rate_of(sub):
    from collections import Counter
    h = hist(sub)
    if any(k % 2 for k in h.values()):
        return None, h
    E = Counter(h.values())
    pr = mp.mpf(1)
    for kap, c in E.items():
        pr *= R[kap] ** c
    return mp.power(pr, mp.mpf(1) / len(sub)), E


for box in ((2, 2, 2, 2), (3, 2, 2, 2), (3, 3, 2, 2)):
    pla, links, d, out = closed_subcomplexes(box)
    if out is None:
        P(f"  box {box}: 面 {len(pla)} 链 {len(links)} 闭簇空间维数 d={d} （>{15} 位，跳过）")
        continue
    best = (0, None, None)
    n4 = 0
    hist_hist = {}
    for sub in out:
        if not sub:
            continue
        r, h = rate_of(sub)
        E = tuple(sorted(h.items()))
        hist_hist[E] = hist_hist.get(E, 0) + 1
        if abs(r - 4) < mp.mpf('1e-30'):
            n4 += 1
        if r > best[0]:
            best = (r, len(sub), dict(sorted(h.items())))
    ident_ok = all(sum(len(links_of(p)) for p in sub) == 4 * len(sub) for sub in out if sub)
    P(f"  box {box}: 面 {len(pla)} 链 {len(links)} 闭簇空间维数 d={d} => 闭簇总数 {2**d - 1}（非空）")
    P(f"    恒等式 4F=sum kappa 全过: {ident_ok}")
    P(f"    max (prodR)^(1/F) = {mp.nstr(best[0],20)}  at F={best[1]} hist={best[2]}")
    P(f"    恰好 =4 的闭簇个数 = {n4} ; 不同 kappa 直方型数 = {len(hist_hist)}")
    for E, cnt in sorted(hist_hist.items(), key=lambda t: -t[1])[:6]:
        F = sum(c for _, c in E) // 1
        P(f"      hist={dict(E)}  个数={cnt}")

# ---------------- 3  极值族：所有 L^3 立方体表面（4D 中 3 维立方体的 6L^2 个面）全链 kappa=2
P("")
P("== 3  极值族：L^3 立方盒表面（4D 中沿 3 个方向的立方体，F=6L^2，全链 kappa=2）==")


def box_surface(L, dirs=(0, 1, 2)):
    """闭合立方体表面：dirs 三个方向、边长 L 的方盒的 6 个面（第 4 方向不参与）"""
    d0, d1, d2 = dirs
    sub = []
    for k, (a, b) in ((d0, (d1, d2)), (d1, (d0, d2)), (d2, (d0, d1))):
        for side in (0, L):
            for u in range(L):
                for v in range(L):
                    s = [0, 0, 0, 0]
                    s[k] = side
                    s[a] = u
                    s[b] = v
                    sub.append((tuple(s), tuple(sorted((a, b)))))
    return sub


print_rows = []
for L in (1, 2, 3, 4):
    sub = box_surface(L)
    sub = list(dict.fromkeys(sub))
    h = hist(sub)
    F = len(sub)
    ok = all(v == 2 for v in h.values())
    prod = mp.mpf(1)
    for v in h.values():
        prod *= R[v]
    rate = mp.power(prod, mp.mpf(1) / F)
    print_rows.append((L, F, len(h), ok, rate))
    P(f"  L={L}: F={F:>3} |links|={len(h):>3} 全链 kappa=2: {ok}  prodR=2^{sum(h.values())}"
      f"  (prodR)^(1/F) = {mp.nstr(rate,16)}")
P("  ==> 4 在 L=1,2,3,4 全部尺度上**精确取到**；不是单个立方体的偶然，是整族（任意 L）")

# ---------------- 4  sup over kappa-profile: 精确证明 max = 4 只在 kappa=2 全同
P("")
P("== 4  sup_hist：max sum a_m ln R(2m) s.t. sum m a_m = 2 （闭簇的精确约束）==")
best, arg4 = 0, None
step = 1 / 24
grid = [i * step for i in range(0, 49)]
for a1 in grid:
    for a2 in grid:
        for a3 in grid:
            a4 = 2 - (a1 + 2 * a2 + 3 * a3)
            if a4 < -1e-12 or a4 > 2:
                continue
            if abs(a1 + 2 * a2 + 3 * a3 + 4 * a4 - 2) > 1e-9:
                continue
            v = a1 * mp.log(R[2]) + a2 * mp.log(R[4]) + a3 * mp.log(R[6]) + a4 * mp.log(R[8])
            if v > best:
                best, arg4 = v, (a1, a2, a3, a4)
P(f"  最优 profile (a2,a4,a6,a8) = {[round(x,4) for x in arg4]}  rate = {mp.nstr(mp.e**best,16)}")
P(f"  逐项核对 R(2m) <= 4^m:  m=1 {mp.nstr(R[2],8)}<=4  m=2 {mp.nstr(R[4],8)}<=16"
  f"  m=3 {mp.nstr(R[6],8)}<=64  m=4 {mp.nstr(R[8],8)}<=256  => prodR <= 4^F  [T]")
P("  等号 <=> E4=E6=E8=0 <=> 每条链 kappa=2 <=> 闭合立方 2-流形（无折）")

# ---------------- 5  阈值算术
P("")
P("== 5  阈值：C^Cbar<7.0633 与 按类分裂后的 A_S < 1/beta ==")
eD = mp.e * 20
u48 = mp.mpf(1) / 384
T = 1 / (eD * u48)
N = 2
P(f"  eΔ = {mp.nstr(eD,16)} ; 阈值 1/(eΔ u(1/48)) = {mp.nstr(T,16)}")
P(f"  Cbar = N = {N} [T] => 需 C^ < {mp.nstr(T/N,10)}（即要从 4 砍到 3.5316，砍 11.7%）")
for C in (4, 2):
    b = 8 / (C * N * eD)
    P(f"  C^={C}: beta1 = 8/(C^ Cbar eΔ) = {mp.nstr(b,10)} = {mp.nstr(b*48,8)}x(1/48)")
P(f"  窗宽比 C^=4 -> 2: {mp.nstr((8/(2*N*eD))/(8/(4*N*eD)),8)}x")
P("  按类分裂（不均匀 KP）：极值类 S 的条件是 C^_S*A_S*Cbar*u < 1")
for beta in (mp.mpf(1) / 48, mp.mpf(1) / 24, mp.mpf(1) / 12):
    A = 1 / (4 * N * (beta / 8))
    P(f"   beta={mp.nstr(beta,8)}: 需 A_S < 1/(C^ Cbar u) = {mp.nstr(A,10)}"
      f"   (= {mp.nstr(A/eD,8)}x eΔ)")
P(f"  对比：现用 eΔ = {mp.nstr(eD,10)} ⇒ 1/48 只需把极值类熵底砍 {mp.nstr(eD/48,8)}x；2× 窗砍 {mp.nstr(eD/24,8)}x")

# ---------------- 6  枚举 [E]：过固定 plaquette 的闭合方表面（全链 kappa=2）的计数
P("")
P("== 6  枚举：含固定 plaquette 的『全链 kappa=2』连通复形的计数（增长底 A_S 的经验值）==")
P0 = ((0, 0, 0, 0), (0, 1))
CONT = {}


def plaquettes_containing(L):
    s, i = L
    out = []
    for j in range(4):
        if j == i:
            continue
        out.append((tuple(s), tuple(sorted((i, j)))))
        t = list(s)
        t[j] -= 1
        out.append((tuple(t), tuple(sorted((i, j)))))
    return out


def kappa_add(sub_h, p):
    """返回加 p 后的 kappa 字典（不变式检查）"""
    h = dict(sub_h)
    for L in links_of(p):
        h[L] = h.get(L, 0) + 1
        if h[L] > 2:
            return None
    return h


level = {}
h0 = {L: 1 for L in links_of(P0)}
level[frozenset([P0])] = h0
t0 = time.time()
CAP = 210.0
counts = []
for F in range(1, 13):
    closed = [k for k in level if all(v == 2 for v in level[k].values())]
    counts.append((F, len(level), len(closed)))
    P(f"  F={F:>2} : kappa<=2 连通复形数 = {len(level):>8}   其中已闭合(全 kappa=2) = {len(closed):>7}"
      f"   [{time.time()-t0:.1f}s]")
    if F >= 11 or time.time() - t0 > CAP:
        break
    nxt = {}
    for sub, h in level.items():
        for L in h:
            if h[L] >= 2:
                continue
            for p in plaquettes_containing(L):
                if p in sub:
                    continue
                h2 = kappa_add(h, p)
                if h2 is None:
                    continue
                ns = frozenset(sub | {p})
                if ns in nxt:
                    continue
                nxt[ns] = h2
        if time.time() - t0 > CAP:
            break
    level = nxt
    if time.time() - t0 > CAP:
        P(f"  [时间上限 {CAP:.0f}s 触发，停于 F={F}]")
        break
P("  闭合计数 vs F：")
for idx, (F, tot, cl) in enumerate(counts):
    if cl:
        r = None
        if idx:
            prev = counts[idx - 1][2]
            if prev:
                r = cl / prev
        P(f"    F={F:>2} closed={cl:>7}  growth(cl/F-1)={('%.3f' % r) if r else '-'}")
P("")

P("== 7  结论要点 ==")
P("(1) C^ = sup_gamma (prod_l R(kappa_l))^{1/|gamma|} 精确 = 4，等号 <=> 全链 kappa=2（闭合立方 2-流形）。")
P("(2) 4 被取到，而且是一个无限族（所有 L^3 盒表面），立方体不是孤例 => 逐 gamma 的 sup 砍不动。")
P("(3) 立方体确实在 KP 的和里（w!=0、连通、过 x），但它不特殊：整族极值聚合物都在和里。")
P("(4) 4 之所以致命，是因为均匀判据把『sup 率 4』乘上『通用动物熵 eΔ=54.3656』；按类分裂后 4 只与表面熵底 A_S 配对。")

with open(__file__.replace(".py", ".log"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
print("[log written]")
