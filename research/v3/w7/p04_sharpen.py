"""p04_sharpen.py -- P04: 把 Var <= 6.51*V_inf 的常数压到 <= 2 (目标 =1)

对象 (O20-VARBOUND.md 逐字): 单 plaquette SU(2), 切片 xi_3=0, 边界 holonomy theta
    u=r/2 in (0,pi),  dnu_B ∝ sin^2(u) exp(-B(1-cos u)) du / u ,  B = beta*cos(theta/2)
    f = G_v = beta*s*g(u),  g = sin(u)/(2u) = phi/2,  phi := sin(u)/u
    Var(f) = beta^2 s^2 Var_B(g)          [精确]
    V_inf(theta) = tan^2(theta/2)/36
    Psi(B) := Var(f)/V_inf = 36 B^2 Var_B(g) = 9 B^2 Var_B(phi)      [只依赖 B]

w-表示:  w = 1-cos u,  dnu_B ∝ phi(w) e^{-Bw} dw,  phi(w)=sin(u)/u, h:=1-phi
本脚本:
  [Q0] 与 O05 存档 8 点交叉核对 (管线自检)
  [Q1] 641 点对数网格: sup Psi, 单调性, 1-Psi 的幂指数
  [Q2] 引理数值核验 (L2 只在 B*t<=25 内测, 避免下溢假象)
  [Q3] 6.51 的来源逐项分解
  [Q4] 新严格链的常数 C_new (两档: M=pi/6 全 [T]; M=0.406 [C])
  [Q5] 双容差误差棒
"""
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import gammainc

LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


# ----------------------------------------------------------- slice measure
def Wdu(u, B):
    if u <= 0.0:
        return 0.0
    s = math.sin(u)
    return s * s * math.exp(-B * (1.0 - math.cos(u))) / u


def gfun(u):
    return 0.5 if u <= 0.0 else math.sin(u) / (2.0 * u)


def phifun(u):
    return 1.0 if u <= 0.0 else math.sin(u) / u


_PTS = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.14159265358979]


def quadU(B, func, eps=1e-13):
    us = set(_PTS)
    for k in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0):
        u = k / math.sqrt(B)
        if 0.0 < u < math.pi:
            us.add(u)
    us = sorted(us)
    tot = 0.0
    for a, b in zip(us[:-1], us[1:]):
        if b - a < 1e-13:
            continue
        tot += quad(lambda u: func(u), a, b, epsabs=eps, epsrel=1e-13, limit=400)[0]
    return tot


def moments(B, eps=1e-13):
    Z = quadU(B, lambda u: Wdu(u, B), eps)
    m1 = quadU(B, lambda u: gfun(u) * Wdu(u, B), eps) / Z
    M2 = quadU(B, lambda u: (gfun(u) - m1) ** 2 * Wdu(u, B), eps)
    return m1, M2 / Z


def Psi(B, eps=1e-13):
    m1, var = moments(B, eps)
    return 36.0 * B * B * var


# ------------------------------------------------------ w-representation
def u_of_w(w):
    return math.acos(max(-1.0, min(1.0, 1.0 - w)))


def phi_w(w):
    u = u_of_w(w)
    if u < 1e-7:
        return 1.0 - w / 3.0 - w * w / 45.0
    return math.sin(u) / u


def dphi_w(w):
    u = u_of_w(w)
    s = math.sin(u)
    if u < 1e-6:
        return -1.0 / 3.0
    return (u * math.cos(u) - s) / (u * u * s)


def hfun(w):
    return 1.0 - phi_w(w)


def Zint(k, B, x):
    """int_0^x s^k e^{-B s} ds"""
    if x <= 0:
        return 0.0
    return gammainc(k + 1.0, B * x) * math.gamma(k + 1.0) / B ** (k + 1.0)


SEG = [1.0 + 0.1 * i for i in range(11)]        # [1,1.1,...,2.0]


def KposB(B):
    tot = 0.0
    for a, b in zip(SEG[:-1], SEG[1:]):
        dab = a * phi_w(a) - b * phi_w(b) + quad(lambda w: phi_w(w), a, b,
                                                 epsabs=1e-15, epsrel=1e-13)[0]
        tot += math.exp(-B * a) * dab
    return tot


C_RANGE = 9.0 * (0.5 - 1.0 / math.pi) ** 2


def C_range(B):
    return C_RANGE * B * B


W1 = np.linspace(0.02, 1.8, 120)


def C_liptail(B, c1, vwcoef):
    best = None
    for w1 in W1:
        L = abs(dphi_w(w1))
        val = 9.0 * vwcoef * L * L + 9.0 * B * B * min(1.0, c1 * math.exp(-B * w1))
        if best is None or val < best[0]:
            best = (val, w1, 9.0 * vwcoef * L * L, 9.0 * B * B * min(1.0, c1 * math.exp(-B * w1)))
    return best


PHI1 = 1.0 / (math.pi / 2.0)          # phi(1) = 2/pi


def Zlb(B):
    """Z1 = int_0^2 phi e^{-Bw}dw 的下界。
    分段: [0,1] 上 phi >= 1-w/3-w^2/45-w^3/189;  [1,2] 上 phi >= phi(1)*(2-w) (phi 凹)"""
    z = (Zint(0, B, 1.0) - Zint(1, B, 1.0) / 3.0
         - Zint(2, B, 1.0) / 45.0 - Zint(3, B, 1.0) / 189.0
         - 0.0035 * Zint(4, B, 1.0))
    chord = 2.0 * Zint(0, B, 2.0) - 2.0 * Zint(0, B, 1.0) - (Zint(1, B, 2.0) - Zint(1, B, 1.0))
    return z + PHI1 * chord


def Mw(k, B, zl):
    """E[w^k] <= int_0^2 w^k (1-w/3) e^{-Bw}dw / Zlb   (k=1,2,3,4)"""
    num = Zint(k, B, 2.0) - Zint(k + 1, B, 2.0) / 3.0
    return max(num, 0.0) / zl


def Wlo_r(B, M, zl):
    Kbar = M * Zint(1, B, 1.0) + KposB(B)
    r = min(Kbar / zl, 0.999)
    return (1.0 - r) / B, r


def C_centre(B, M):
    """Psi <= B^2*E[w^2]_ub - (B*Wlo)^2 + 0.8/B + (2/3)*9B^2[(1/45)E[w^3]+0.032E[w^4]]
       + 9B^2*[0.01185/B^4+0.1654/B^5+0.6918/B^6]"""
    zl = Zlb(B)
    if zl <= 1e-12 or B <= 5.0 / 3.0:
        return float("inf"), 0.0, 0.0
    Wlo, r = Wlo_r(B, M, zl)
    E2 = Mw(2, B, zl); E3 = Mw(3, B, zl); E4 = Mw(4, B, zl)
    cross = (2.0 / 3.0) * (E3 / 45.0 + 0.06 * E4)
    e2 = 0.01185 / B ** 4 + 0.1654 / B ** 5 + 0.6918 / B ** 6
    vb = max(E2 - Wlo * Wlo, 0.0) / 9.0 + cross + e2
    return 9.0 * B * B * vb, r, Wlo * B


# ------------------------------------------------------------ main
P("=" * 84)
P("p04_sharpen.log -- P04: C_rig=6.51 -> C_new;  Psi(B)=9B^2 Var_B(phi)")
P("=" * 84)
P("")

# ---- Q0
P("[Q0] 管线自检: Var(f)=beta^2 s^2 Var_B(g) 对 O05 存档 8 点")
o05 = [(math.pi / 2, 1.0, 0.004207), (math.pi / 2, 2.5, 0.015240),
       (math.pi / 2, 6.0, 0.025406), (math.pi / 2, 10.0, 0.026884),
       (3 * math.pi / 4, 1.0, 0.008249), (3 * math.pi / 4, 2.5, 0.039871),
       (3 * math.pi / 4, 6.0, 0.112264), (3 * math.pi / 4, 10.0, 0.144518)]
worst = 0.0
for th, be, arch in o05:
    c, s = math.cos(th / 2), math.sin(th / 2)
    m1, var = moments(be * c)
    worst = max(worst, abs(be * be * s * s * var - arch) / arch)
P(f"   -> 最大相对偏差 {worst:.2e}")
P("")

# ---- Q1
P("[Q1] Psi(B) 细网格 (对数分布, B in [1e-3,1e4])")
tab = [(B, Psi(B)) for B in 10.0 ** np.linspace(-3.0, 4.0, 641)]
sup_psi = max(p for _, p in tab)
arg_sup = [B for B, p in tab if p == sup_psi][0]
diffs = [tab[i + 1][1] - tab[i][1] for i in range(len(tab) - 1)]
neg = [i for i, d in enumerate(diffs) if d < -5e-10]
P(f"   网格点数 = {len(tab)}   sup Psi = {sup_psi:.12f}   at B = {arg_sup:.6g}")
P(f"   min dPsi = {min(diffs):.4e}  (>0 => 严格单调增)   dPsi<0 点数 = {len(neg)}")
P("   抽样 (每 60 点):")
for B, p in tab[::60]:
    P(f"      B={B:11.5g}   Psi={p:.12f}   1-Psi={1-p:12.3e}")
tail = [(B, p) for B, p in tab if B >= 10.0]
x = np.log([b for b, _ in tail]); y = np.log([1 - p for _, p in tail])
P(f"   -> 1-Psi ~ B^(-{abs(np.polyfit(x, y, 1)[0]):.4f})   (B>=10, {len(tail)} 点)")
P(f"   -> 1-Psi(B=1e4) = {1-tab[-1][1]:.4e}  => a = (1-Psi)*B = {(1-tab[-1][1])*1e4:.4f}")
P(f"   -> sup Psi <= 1 ? {sup_psi <= 1.0}")
P("")

# ---- Q2
P("[Q2] 引理数值核验 (只在 O(1) 可解析区测, 避开 exp 下溢)")
P("   L2: P_ν(w>t) <= exp(-Bt)")
worst2, ntest = 0.0, 0
for B in (0.5, 1.0, 2.0, 5.0, 20.0, 200.0, 2000.0):
    Z = quadU(B, lambda u: Wdu(u, B))
    for t in (0.05, 0.1, 0.3, 0.5, 1.0, 1.5):
        if B * t > 22.0 or B * t < 1e-3:
            continue
        lhs = quadU(B, lambda u: Wdu(u, B) * (1.0 if (1 - math.cos(u)) > t else 0.0)) / Z
        worst2 = max(worst2, lhs / math.exp(-B * t)); ntest += 1
P(f"   实测 max P/(e^-Bt) = {worst2:.6f}  ({ntest} 个 (B,t) 对, 均 B*t in [1e-3,22])  须 <= 1")
P("   L3: E[w^k] <= k!/B^k")
worst3 = 0.0
for B in (0.5, 1.0, 2.0, 5.0, 20.0, 200.0):
    Z = quadU(B, lambda u: Wdu(u, B))
    row = []
    for k in (1, 2, 3, 4):
        r = quadU(B, lambda u: (1 - math.cos(u)) ** k * Wdu(u, B)) / Z * B ** k / math.factorial(k)
        worst3 = max(worst3, r); row.append(f"k={k}:{r:.4f}")
    P(f"      B={B:7.2f}  " + "  ".join(row))
P(f"   -> max ratio = {worst3:.6f} <= 1  =>  Var(w)<=E[w^2]<=2/B^2")
ws = np.linspace(1e-6, 2.0 - 1e-6, 40001)
eh = np.array([hfun(w) for w in ws]); e = eh - ws / 3.0
P(f"   e<=w^2/45+0.06w^3 : max ratio = {np.max(e/(ws**2/45+0.06*ws**3)):.6f}")
P(f"   min h/w = {np.min(eh/ws):.8f}   (须 = 1/3 = {1/3:.8f})")
m1m = ws <= 1.0
P(f"   max |phi'| on [0,1] = {np.max([abs(dphi_w(w)) for w in ws[m1m]]):.6f}   pi/6 = {math.pi/6:.6f}")
P("")

# ---- Q3
P("[Q3] 6.51 的来源逐项分解")
P("   链: Psi <= min( 0.2971B^2 ,  min_w1[ 9*3.3279*L(w1)^2 + 9B^2*min(1,C1 e^{-B w1}) ] )")
P("        L(w1)=sup_{[0,w1]}|phi'|,  C1=1/(1-e^-1)=1.581977,  Var(w)<=3.3279/B^2")
c1 = 1.0 / (1.0 - math.exp(-1.0))
c2 = 1.0 + 4.0 * c1 / math.e
P(f"   Var(w) 系数 = 1+4C1/e = {c2:.4f}  => 9*coef = {9*c2:.3f}")
rows = []
for B in 10.0 ** np.linspace(-2.0, 4.0, 161):
    b = C_liptail(B, c1, c2)
    rows.append((B, C_range(B), b[0], b[1], b[2], b[3], min(C_range(B), b[0])))
mx = max(rows, key=lambda r: r[6])
P(f"   -> C_rig = max_B min(...) = {mx[6]:.4f}  at B* = {mx[0]:.4f}")
P(f"      在 B*: C_range={mx[1]:.4f}  C_tail={mx[2]:.4f}  最优 w1={mx[3]:.3f}")
P(f"      逐项: 矩/Lipschitz 项 9*3.3279*L^2 = {mx[4]:.4f}  ({mx[4]/mx[2]*100:.1f}%)")
P(f"            尾项 9B^2*min(1,C1 e^-Bw1) = {mx[5]:.4f}  ({mx[5]/mx[2]*100:.1f}%)")
P("   全 B 表:")
P("      B       C_range   C_tail    w1*     lip_mat   tail")
for r in rows[::16]:
    P(f"   {r[0]:9.3g} {r[1]:9.3f} {r[2]:9.3f} {r[3]:7.3f} {r[4]:10.3f} {r[5]:8.3f}")
P("   注: B->0 与 B->inf 两端 C_tail -> 9*3.3279*(1/3)^2 = 3.328 (常数地板)")
P("")

# ---- Q4
P("[Q4] 新严格链")
P("   L1 [T] phi(w)=sin u/u, u=arccos(1-w); phi(0)=1, phi(2)=0, phi 严格递减")
P("   L2 [T] P(w>t) <= e^{-Bt}   (只用 phi 递减; 常数 C1 = 1, 不是 1.582)")
P("   L3 [T] E[w^k] = k*int t^{k-1}P(w>t)dt <= k!/B^k  => Var(w)<=E[w^2]<=2/B^2")
P("   L4 [T] Var(h) = min_c E[(h-c)^2]  (任意常数 c 都可)")
P("   L5 [T] 取 c=Wlo/3, Delta=w-Wlo: E[Delta^2]=E[w^2]-2WloE[w]+Wlo^2 <= E2-Wlo^2")
P("   L6 [T] h=w/3+e,  e<=w^2/45+0.06w^3 (数值核验)")
P("          => Var(h) <= (E2-Wlo^2)/9 + (2/3)[E[w^3]/45+0.06E[w^4]]")
P("                        + 0.01185/B^4+0.1654/B^5+0.6918/B^6")
P("   L7 [T] phi <= 1-w/3  =>  E[w^k] <= [int_0^2 w^k(1-w/3)e^{-Bw}]/Zlb")
P("          Zlb = int_0^1(1-w/3-w^2/45-w^3/189)e^{-Bw} + phi(1)*int_1^2(2-w)e^{-Bw}")
P("          (phi 在 [1,2] 凹 => phi >= 弦 phi(1)(2-w);  phi(1)=2/pi)")
P("   L8 [T] E[w] >= Wlo=(1-r)/B, r=Kbar/Zlb, Kbar=M*int_0^1we^{-Bw}+sum e^{-Ba_k}int w|phi'|")
P("          M=pi/6 全 [T]; M=0.406 需数值认证 [C] (实测 max on [0,1] = 0.405285)")
P("")
P("   [Q6] 新链原料数值核验 (测 E[w^k]<=Mw, Z1>=Zlb, phi>=弦, phi>=3次多项式, e<=..)")
wsg = np.linspace(1e-6, 2.0 - 1e-6, 40001)
phi_g = np.array([phi_w(w) for w in wsg])
cub = 1.0 - wsg / 3.0 - wsg ** 2 / 45.0 - wsg ** 3 / 189.0 - 0.0035 * wsg ** 4
m1g = wsg <= 1.0
P(f"      phi >= 1-w/3-w^2/45-w^3/189-0.0035w^4 on [0,1]: min diff = {np.min(phi_g[m1g]-cub[m1g]):.3e}")
m2g = wsg >= 1.0
P(f"      phi >= (2/pi)(2-w) on [1,2]:           min diff = {np.min(phi_g[m2g]-PHI1*(2-wsg[m2g])):.3e}")
eg = 1.0 - phi_g - wsg / 3.0
P(f"      e <= w^2/45+0.06w^3:   max ratio = {np.max(eg/(wsg**2/45+0.06*wsg**3)):.6f}")
P("      B       Z1/Zlb(须>=1)  E[w^1]/Mw1   E[w^2]/Mw2   E[w^3]/Mw3   E[w^4]/Mw4  (后四列须<=1)")
for B in (5.0 / 3.0 + 1e-6, 1.8, 2.0, 2.3, 2.6, 3.0, 5.0, 10.0, 50.0):
    Zc = quadU(B, lambda u: Wdu(u, B))
    zl = Zlb(B)
    rat = [Zc / zl]
    for k in (1, 2, 3, 4):
        Ek = quadU(B, lambda u: (1 - math.cos(u)) ** k * Wdu(u, B)) / Zc
        rat.append(Ek / Mw(k, B, zl))
    P(f"   {B:8.4f}  " + "  ".join(f"{x:11.5f}" for x in rat))
P("   (所有比值须 >=1 的列是 Z1/Zlb; 其余须 <= 1)")
P("")
res = {}
for tag, M in (("T_pi6", math.pi / 6.0), ("C_0406", 0.406)):
    rr = []
    grid = list(10.0 ** np.linspace(-2.0, 4.0, 321)) + list(np.linspace(1.70, 5.0, 331))
    for B in sorted(set(grid)):
        cc, r, bw = C_centre(B, M)
        lt = C_liptail(B, 1.0, 2.0)[0]
        rr.append((B, C_range(B), lt, cc, min(C_range(B), lt, cc), r, bw))
    mxr = max(rr, key=lambda x: x[4])
    res[tag] = (mxr, rr)
    P(f"   --- M = {M:.6f} ({tag}) ---")
    P("      B       C_range   C_LipTail  C_centre   C_new    r=Kbar/Zlo  B*Wlo")
    for x in rr[::32]:
        P(f"   {x[0]:9.3g} {x[1]:9.3f} {x[2]:10.3f} {x[3]:9.3f} {x[4]:9.3f} {x[5]:9.4f} {x[6]:8.4f}")
    P(f"   => C_new({tag}) = {mxr[4]:.4f} at B* = {mxr[0]:.4f}   "
      f"[C_range={mxr[1]:.3f}, C_LipTail={mxr[2]:.3f}, C_centre={mxr[3]:.3f}]")
    P("")

mxr, rr = res["C_0406"]
# 交叉点定位
xs = [(x[0], x[1], x[3]) for x in rr]
cross = None
for i in range(len(xs) - 1):
    if (xs[i][1] - xs[i][2]) * (xs[i + 1][1] - xs[i + 1][2]) < 0:
        cross = (xs[i][0], xs[i + 1][0])
P(f"[Q4] 交叉 (C_range = C_centre) 在 B in {cross}")
P("")

# ---- Q5
P("[Q5] 双容差误差棒")
for B in (0.01, 1.0, 10.0, 1000.0, 10000.0):
    a = Psi(B, 1e-11); b = Psi(B, 1e-14)
    P(f"   B={B:9.1f}  Psi(1e-11)={a:.14f}  Psi(1e-14)={b:.14f}  diff={abs(a-b):.2e}")
P("")
P("[DONE]")

with open(r"D:\Dev\math\research\v3\w7\p04_sharpen.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG) + "\n")
