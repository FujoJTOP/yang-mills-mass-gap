# w7/uu2_sign_check.py — UU2 校验（零 MC：纯组合 + Bessel + Weyl 求积）
# 1) 同点 6 面 15 对的共用 link 计数  2) 超立方轨道
# 3) 2 片 ladder 的精确恒等式 Cov≡0     4) Weyl 求积核 int chi_j phi = delta_{j,1/2}
import itertools, math, sys

def add(x, a):
    return tuple(x[i] + (1 if i == a else 0) for i in range(4))

def plaq_links(x, mu, nu):
    # link 键 = (base, dir)，base 取该方向坐标较小的端点
    return {(x, mu), (x, nu), (add(x, mu), nu), (add(x, nu), mu)}

planes = [(m, n) for m in range(4) for n in range(m + 1, 4)]
x0 = (0, 0, 0, 0)
Lk = {p: plaq_links(x0, *p) for p in planes}
pairs = list(itertools.combinations(planes, 2))

print("[1] planes =", planes, " pairs = C(6,2) =", len(pairs))
hist = {}
for a, b in pairs:
    k = len(Lk[a] & Lk[b])
    hist[k] = hist.get(k, 0) + 1
print("[1] shared-link histogram =", dict(sorted(hist.items())))
print("[1] shared=0 pairs =", [(a, b) for a, b in pairs if not (Lk[a] & Lk[b])])
print("[1] shared=1 pairs =", [(a, b, sorted(Lk[a] & Lk[b])) for a, b in pairs if len(Lk[a] & Lk[b]) == 1])

# 超立方点群（带符号置换）在 6 面上的作用 = S4；数 15 对的轨道
seen, orbits = set(), []
for a, b in pairs:
    if (a, b) in seen:
        continue
    orb = set()
    for p in itertools.permutations(range(4)):
        aa = tuple(sorted((p[a[0]], p[a[1]])))
        bb = tuple(sorted((p[b[0]], p[b[1]])))
        if aa > bb:
            aa, bb = bb, aa
        orb.add((aa, bb))
    seen |= orb
    orbits.append(sorted(orb))
print("[2] orbits =", len(orbits), " sizes =", [len(o) for o in orbits])
for o in orbits:
    print("    orbit size", len(o), " rep:", o[0])

try:
    from scipy.special import iv
except Exception as e:
    print("NO-SCIPY:", e)
    sys.exit(0)

# e^{a tr X} = sum_j c_j chi_j(X),  c_j = (2j+1) I_{2j+1}(2a)/a,  a=beta/2
def cj(j, beta):
    a = beta / 2.0
    return (2 * j + 1) * iv(2 * j + 1, 2 * a) / a

print("[3] ladder: <phi> = c_(1/2)/c_0,  <phi1 phi2> = (c_(1/2)/c_0)^2  =>  Cov(q1,q2)=0")
for beta in (0.5, 2.5):
    c0, ch, c1 = cj(0, beta), cj(0.5, beta), cj(1.0, beta)
    r = ch / c0
    print(f"    beta={beta}: c0={c0:.12g} c(1/2)={ch:.12g} c1={c1:.12g}  "
          f"ratio={r:.12g}  Cov_ladder={0.25*(r*r-r*r):.1e}")

# Weyl 求积: int f dU = (2/pi) int_0^pi f(theta) sin^2 theta d theta
def chi(j, th):
    s = math.sin(th)
    return (2 * j + 1.0) if abs(s) < 1e-14 else math.sin((2 * j + 1) * th) / s

N = 200000
def quad(f):
    h = math.pi / N
    s = 0.0
    for i in range(N + 1):
        th = i * h
        w = 1 if i in (0, N) else (4 if i & 1 else 2)
        s += w * f(th) * math.sin(th) ** 2
    return s * h / 3 * (2 / math.pi)

for j in (0.0, 0.5, 1.0, 1.5, 2.0):
    print(f"[4] int chi_{j} * phi dU = {quad(lambda t, j=j: chi(j, t) * chi(0.5, t)): .10f}")
print("[4] int phi^2 =", round(quad(lambda t: chi(0.5, t) ** 2), 10))
print("[4] int phi^3 =", round(quad(lambda t: chi(0.5, t) ** 3), 10))
print("[4] int phi^4 =", round(quad(lambda t: chi(0.5, t) ** 4), 10))
print("[5] OK exit0")
