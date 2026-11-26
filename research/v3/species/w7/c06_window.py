"""C06-WINDOW: KP 要的是逐 gamma 还是逐 plaquette 包络？精确 2 还是包络 4？
合并常数 K^ = sup_gamma (prod_l R_l * tau_gamma)^{1/|gamma|} 的可证界。
py -3.12 + mpmath 1.3.0, dps=40. 只读，不写仓库其他文件。
"""
import itertools, random
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 40
OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.append(s)


# ---------------- 1  Weyl: mult(1 in f^n), mult(adj in f^n), R(kappa)
def haar(f):
    """(1/pi) int_0^{2pi} sin^2(th/2) f(th) dth"""
    return mp.quad(lambda th: mp.sin(th / 2) ** 2 * f(th), [0, mp.pi, 2 * mp.pi]) / mp.pi


chi_f = lambda th: 2 * mp.cos(th / 2)
chi_adj = lambda th: 1 + 2 * mp.cos(th)

P("== 1  SU(2) Weyl (dps=40) ==")
P("check int chi_f^2 =", mp.nstr(haar(lambda t: chi_f(t) ** 2), 20), " (need 1)")
P("check int chi_adj^2 =", mp.nstr(haar(lambda t: chi_adj(t) ** 2), 20), " (need 3)")
m1, madj, R = {}, {}, {}
for n in range(1, 9):
    m1[n] = haar(lambda t: chi_f(t) ** n)
    madj[n] = haar(lambda t: chi_adj(t) * chi_f(t) ** n)
P(f"{'n':>2} {'mult(1 in f^n)':>18} {'mult(adj in f^n)':>18}  Catalan?"""
  )
cat = [1, 1, 2, 5, 14, 42, 132]
for n in range(1, 9):
    ok = ""
    if n % 2 == 0:
        ok = "Catalan(%d)=%d  match=%s" % (n // 2, cat[n // 2], abs(m1[n] - cat[n // 2]) < mp.mpf('1e-30'))
    P(f"{n:>2} {mp.nstr(m1[n], 12):>18} {mp.nstr(madj[n], 12):>18}  {ok}")
for n in (2, 4, 6):
    R[n] = 1 + madj[n] / m1[n]
    P(f"R(kappa={n}) = 1 + {mp.nstr(madj[n],10)}/{mp.nstr(m1[n],10)} = {mp.nstr(R[n], 20)}")
P("2^{k/2} 上界核对: R(2)=2 <= 2^1=2 ; R(4)=2.5 <= 2^2=4 ; R(6)=%.6f <= 2^3=8" % R[6])
P("")

# ---------------- 2  4D plaquette complexes: Euler/handshake identity
# box (n1,n2,n3,n4) sites; plaquettes = (site, plane 0<=i<j<4)
def build(n):
    sites = list(itertools.product(*[range(k) for k in n]))
    pla = []
    for s in sites:
        for i in range(4):
            for j in range(i + 1, 4):
                if s[i] + 1 < n[i] and s[j] + 1 < n[j]:
                    pla.append((s, (i, j)))
    return pla


def links_of(p):
    s, (i, j) = p
    si, sj = list(s), list(s)
    si[i] += 1
    sj[j] += 1
    return [(tuple(s), i), (tuple(sj), i), (tuple(s), j), (tuple(si), j)]


def kappa_hist(sub):
    h = {}
    for p in sub:
        for L in links_of(p):
            h[L] = h.get(L, 0) + 1
    return h


def stats(sub):
    h = kappa_hist(sub)
    F = len(sub)
    E = {}
    for v in h.values():
        E[v] = E.get(v, 0) + 1
    return F, E


def trace_parts(sub):
    """connected components by shared link"""
    lp = [set(links_of(p)) for p in sub]
    seen, comps = set(), []
    for a in range(len(sub)):
        if a in seen:
            continue
        st, comp = [a], []
        seen.add(a)
        while st:
            x = st.pop()
            comp.append(sub[x])
            for b in range(len(sub)):
                if b not in seen and lp[x] & lp[b]:
                    seen.add(b)
                    st.append(b)
        comps.append(comp)
    return comps


def khat(sub, F, E):
    """prod_l R(kappa) ; None if a kappa=1 link kills weight"""
    pr = mp.mpf(1)
    for k, cnt in E.items():
        if k == 1 or k == 3 or k == 5:
            return None
        pr *= R[k] ** cnt
    return pr


P("== 2  手数恒等式 4F = E1 + 2E2 + 4E4 + 6E6 与 E2/F <= 2 ==")
pla = build((2, 2, 2, 2))
P("box (2,2,2,2): plaquettes =", len(pla))
random.seed(7)
bad = 0
mx = 0.0
for _ in range(20000):
    k = random.randint(1, 12)
    sub = random.sample(pla, k)
    F, E = stats(sub)
    s = sum(j * E.get(j, 0) for j in E)
    if s != 4 * F:
        bad += 1
    mx = max(mx, E.get(2, 0) / F)
P(f"random 20000 subsets: identity violations = {bad} ; max E2/F = {mx:.4f}")
# 3-cube 的 6 个面（box (2,2,2,1) 全部 6 片）
ple = build((2, 2, 2, 1))
P("box (2,2,2,1) plaquettes =", len(ple), "(= 3-cube 的 6 面)")
F, E = stats(ple)
P("cube surface: F=%d  hist kappa=%s" % (F, dict(sorted(E.items()))))
P("  E2/F = %.6f (=2 ⟺ 全闭) ; prod R = 2^%d = %s ; 几何均值 = %s"
  % (E[2] / F, E[2], mp.nstr(2 ** E[2], 12), mp.nstr(mp.power(2, mp.mpf(E[2]) / F), 12)))
P("  ==> 逐 gamma 的 sup 已达 4（不是 2）")
# closed subcomplexes only: search for max E2/F
best = (0, None)
for _ in range(20000):
    sub = random.sample(pla, random.randint(2, 10))
    F, E = stats(sub)
    if E.get(1, 0) == 0 and E.get(3, 0) == 0 and E.get(5, 0) == 0:
        r = E.get(2, 0) / F
        if r > best[0]:
            best = (r, sub)
P("closed random search: max E2/F =", best[0], "(theoretical sup = 2)")
P("")

# ---------------- 3  tau (Kirchhoff, exact) + 合并常数 K^
def tau(exact_int, sub):
    """spanning-tree count of the plaquette adjacency multigraph; exact Fraction"""
    n = len(sub)
    if n == 1:
        return Fraction(1)
    lap = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    lp = [set(links_of(p)) for p in sub]
    for a in range(n):
        for b in range(a + 1, n):
            w = len(lp[a] & lp[b])
            if w:
                lap[a][a] += w
                lap[b][b] += w
                lap[a][b] -= w
                lap[b][a] -= w
    M = [row[:-1] for row in lap[:-1]]          # delete last row/col
    det = Fraction(1)
    m = len(M)
    for c in range(m):                          # fraction-free-ish Gaussian
        piv = next((r for r in range(c, m) if M[r][c] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            det = -det
        det *= M[c][c]
        for r in range(c + 1, m):
            f = M[r][c] / M[c][c]
            for cc in range(c, m):
                M[r][cc] -= f * M[c][cc]
    return det


P("== 3  cube 的 tau 与合并常数 ==")
F, E = stats(ple)
t = tau(None, ple)
P("cube adjacency graph = 八面体 K222 : tau =", t, " (known 384)")
kc = khat(ple, F, E)
P("cube: prod R =", mp.nstr(kc, 10), " ; tau =", t, " ; (prod R * tau)^(1/F) =",
  mp.nstr(mp.power(kc * t, mp.mpf(1) / F), 12))

# ---------------- 4  混合型 上界优化：(prod R)^{1/F} * (2*Egraph/F)
P("")
P("== 4  合并常数可证界：按 link 型混合优化 ==")
P("约束 sum_l kappa_l = 4F ；a_k = n_k/F ；G = exp(sum a_k ln R_k) * (2*sum a_k*C(k,2))")
bestg, bestmix = 0, None
for a2 in [i / 40 for i in range(0, 81)]:
    for a4 in [i / 40 for i in range(0, 81)]:
        for a6 in [i / 40 for i in range(0, 81)]:
            if abs(2 * a2 + 4 * a4 + 6 * a6 - 4) > 1e-9:
                continue
            eg = 2 * (a2 * 1 + a4 * 6 + a6 * 15)     # mean adjacency degree
            if eg <= 0:
                continue
            g = mp.e ** (a2 * mp.log(R[2]) + a4 * mp.log(R[4]) + a6 * mp.log(R[6])) * eg
            if g > bestg:
                bestg, bestmix = g, (a2, a4, a6, eg)
P("最优混合 a2,a4,a6 =", [round(x, 3) for x in bestmix[:3]], " 平均邻接度 =", round(bestmix[3], 4))
P("K^ 的可证上界 ≈", mp.nstr(bestg, 8), "   （对比 4*20=80、2*20=40）")
for mix in [(2, 0, 0), (0, 1, 0), (0, 0, 2 / 3)]:
    a2, a4, a6 = mix
    eg = 2 * (a2 + 6 * a4 + 15 * a6)
    g = mp.e ** (a2 * mp.log(R[2]) + a4 * mp.log(R[4]) + a6 * mp.log(R[6])) * eg
    P(f"  纯 kappa 混合 a=({a2:.2f},{a4:.2f},{a6:.2f}): 几何均值={mp.nstr(g/eg,6)} 平均度={eg:.1f} G={mp.nstr(g,8)}")

# ---------------- 5  阈值与窗口
P("")
P("== 5  阈值 / 窗口（40 位） ==")
eD = mp.e * 20
u48 = mp.mpf(1) / 384                       # u = beta/8, beta=1/48
T = 1 / (eD * u48)
P("eΔ =", mp.nstr(eD, 16))
P("阈值 Ĉ*C̄ < 1/(eΔ u(1/48)) =", mp.nstr(T, 12))
P("  Ĉ=4 ⇒ C̄ <", mp.nstr(T / 4, 10), " ;  Ĉ=2 ⇒ C̄ <", mp.nstr(T / 2, 10))
b48 = mp.mpf(1) / 48
rows = [("4 x C̄=20 (Z12,A01)", mp.mpf(80)), ("2 x C̄=20", mp.mpf(40)),
        ("合并 K^ 上界(本节)", bestg), ("合并 K^ 下界(cube 见证)", kc * t and mp.power(kc * t, mp.mpf(1) / F)),
        ("4 x C̄=4 (闭簇 Kirchhoff)", mp.mpf(16)), ("4 x C̄=1 (z02 假定)", mp.mpf(4)),
        ("1 x C̄=1", mp.mpf(1))]
P(f"{'打包':<28}{'K':>12}{'<7.0633?':>10}{'beta1=8/(K eΔ)':>18}{'/(1/48)':>10}")
for name, K in rows:
    b1 = 8 / (K * eD)
    P(f"{name:<28}{mp.nstr(K,8):>12}{('YES' if K<T else 'NO'):>10}{mp.nstr(b1,8):>18}{mp.nstr(b1/b48,6):>10}")
P("")
P("== 6  结论要点 ==")
P("(1) KP(2.1) 只需一个逐 gamma 的几何均值常数；其值(可证且达到) = 4，不是 2。")
P("(2) '2' 是每因子化点(每 kappa=2 link)的值；因子化点/面积 = E2/F = 2 (闭簇, Euler)，故 2^2=4。")
P("(3) 阈值 7.0633 在 '均匀逐 gamma 率' 打包里：可证上界 39.62、cube 见证 10.78 —— 都 > 7.0633" % ())
P("")

with open(__file__.replace(".py", ".log"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
print("[log written]")
