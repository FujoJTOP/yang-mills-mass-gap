"""o20_varbound.py -- O20: Var(G_v|0) 的 beta-一致上界

单 plaquette SU(2), 切片 xi_3=0, 边界 holonomy theta.
切片律 (H1A-EXACT / O05 逐字):  u = r/2 in (0,pi)
    rho_B(u) du ∝ sin^2(u) exp(-B(1-cos u)) du / u ,    B := beta*cos(theta/2)
    f := G_v = beta*sin(theta/2)*sin(r/2)/r = beta*s*sin(u)/(2u)
        => f = beta*s*g(u),   g(u) = sin(u)/(2u)
    => Var(f) = beta^2 s^2 Var_B(g)     (精确, 只依赖 B)

本脚本:
  [P1] 缩放恒等式 + 与 O05 存档 8 个数交叉核对
  [P2] Psi(B) := 36 B^2 Var_B(g) 的网格: 单调性, 上确界, ->1 的幂
  [P3] 闭式 V_inf = tan^2(theta/2)/36 的核对
  [P4] 严格链的常数: L(w1)=|phi'|, phi(w1), E[w^2]B^2, 尾界常数, 统一常数 C_rig
  [P5] 切片上 Hess V 的负模 (BL 无 beta-一致 lambda 的直接证据)
  [P6] 双容差误差棒
"""
import math
import numpy as np
from scipy.integrate import quad

LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


# ---------------------------------------------------------------- reduced moments
def _f(B, power, y):
    if y <= 0.0:
        return 0.0
    u = y / math.sqrt(B)
    s = math.sin(u)
    dens = s * s * math.exp(-B * (1.0 - math.cos(u))) / y
    if power == 0:
        return dens
    g = s / (2.0 * u)
    return dens * g if power == 1 else dens * g * g


def _quad(B, power, eps):
    ymax = math.pi * math.sqrt(B)
    pts = [0.0] + ([20.0] if ymax > 20.0 else []) + [ymax]
    tot = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        tot += quad(lambda y: _f(B, power, y), a, b, epsabs=eps, epsrel=1e-14,
                    limit=400)[0]
    return tot


def var_g(B, eps=1e-13):
    Z0 = _quad(B, 0, eps)
    Z1 = _quad(B, 1, eps)
    Z2 = _quad(B, 2, eps)
    m1, m2 = Z1 / Z0, Z2 / Z0
    return m2 - m1 * m1


def Psi(B, eps=1e-13):
    return 36.0 * B * B * var_g(B, eps)


# ------------------------------------------------------- w-representation (rigor)
def phi(w):
    """density shape in w = 1-cos u:  dnu ∝ phi(w) e^{-B w} dw,  g = phi/2."""
    return math.sqrt(w * (2.0 - w)) / math.acos(1.0 - w)


def dphi(w, h=1e-7):
    return (phi(min(2.0, w + h)) - phi(max(0.0, w - h))) / (min(2.0, w + h) - max(0.0, w - h))


# ---------------------------------------------------------------------- main
P("=" * 78)
P("o20_varbound.log -- O20: Var(G_v|0) 的 beta-一致上界 (单 plaquette SU(2), 切片 xi_3=0)")
P("=" * 78)
P("reduced: u=r/2 in (0,pi); rho_B ∝ sin^2 u exp(-B(1-cos u))/u; B=beta*cos(theta/2)")
P("         g(u)=sin u/(2u);  Var(f)=beta^2 sin^2(theta/2) Var_B(g)   [exact]")
P("")

# ---- [P1] scaling identity vs O05 archive
P("[P1] 缩放恒等式核对: Var(f) = beta^2 s^2 Var_B(g), B=beta*c")
o05 = [  # (theta, beta, Var_archived)
    (math.pi / 2, 1.0, 0.004207), (math.pi / 2, 2.5, 0.015240),
    (math.pi / 2, 6.0, 0.025406), (math.pi / 2, 10.0, 0.026884),
    (3 * math.pi / 4, 1.0, 0.008249), (3 * math.pi / 4, 2.5, 0.039871),
    (3 * math.pi / 4, 6.0, 0.112264), (3 * math.pi / 4, 10.0, 0.144518),
]
P("   theta      beta     B=beta*c    Var(calc)     Var(O05)     rel.diff")
worst = 0.0
for th, be, arch in o05:
    c, s = math.cos(th / 2), math.sin(th / 2)
    B = be * c
    v = be * be * s * s * var_g(B)
    rel = abs(v - arch) / arch
    worst = max(worst, rel)
    P(f"   {th/math.pi:6.3f}pi  {be:6.2f}  {B:8.4f}  {v:11.6f}  {arch:11.6f}  {rel:10.2e}")
P(f"   -> 最大相对偏差 {worst:.2e}  (O05 存档只有 4-6 位有效, 故 1e-3 量级为正常)")
P("")

# ---- [P2] Psi(B) grid
P("[P2] Psi(B) = 36 B^2 Var_B(g) 网格 (should be <=1, increasing, ->1)")
Bs = [10 ** x for x in np.linspace(-3, 4, 36)]
tab = []
for B in Bs:
    tab.append((B, Psi(B)))
prev = None
mono = True
for B, p in tab:
    if prev is not None and p < prev - 5e-9:
        mono = False
    prev = p
for B, p in tab[::4]:
    P(f"   B={B:10.4g}   Psi={p:.12f}   1-Psi={1-p:12.3e}")
P(f"   -> max Psi over grid = {max(p for _, p in tab):.12f}  (<=1 ? {max(p for _,p in tab) <= 1.0})")
P(f"   -> 单调不减(容差5e-9)? {mono}")
# fit exponent of 1-Psi vs B on the large-B tail
tail = [(B, p) for B, p in tab if B >= 20.0 and 1 - p > 1e-12]
if len(tail) >= 3:
    x = np.log([b for b, _ in tail]); y = np.log([1 - p for _, p in tail])
    q = np.polyfit(x, y, 1)
    P(f"   -> 1-Psi ~ B^(-{abs(q[0]):.3f})  (大 B 端, {len(tail)} 点)")
P("")

# ---- [P3] closed form
P("[P3] 闭式: V_inf(theta) = tan^2(theta/2)/36")
P("   theta      tan^2(th/2)/36    O05 'Var_inf'    rel.diff")
for th, arch in ((math.pi / 2, 0.0278), (3 * math.pi / 4, 0.1618)):
    cf = math.tan(th / 2) ** 2 / 36.0
    P(f"   {th/math.pi:6.3f}pi   {cf:.8f}        {arch:.4f}        {abs(cf-arch)/arch:.2e}")
P(f"   theta=pi/2 精确值 1/36 = {1/36:.12f}")
P(f"   theta=3pi/4 精确值 5.8284271.../36 = {math.tan(3*math.pi/8)**2/36:.12f}")
P("   注: O05 的 0.0278/0.1618 是有限 beta 外推值; 闭式是 beta->inf 的极限")
P("")

# ---- [P4] rigorous chain constants
P("[P4] 严格链常数 (w = 1-cos u, dnu ∝ phi(w)e^{-Bw}dw, g=phi/2)")
P("   w1      phi(w1)      |phi'(w1)|   30*L^2")
best = None
for w1 in (0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 1.7):
    L = abs(dphi(w1))
    P(f"   {w1:4.2f}   {phi(w1):9.6f}   {L:10.6f}   {30*L*L:8.4f}")
    if best is None:
        best = (w1, L)
P("")
P("   尾界: nu((eps,2]) <= C1 exp(-B eps)/(1-e^-1),  C1=1/(1-e^-1)=2/... ")
c1 = 1.0 / (1.0 - math.exp(-1.0))
P(f"   C1 = 1/(1-1/e) = {c1:.6f}")
P("   E[w^2] <= 2[1/(2B^2)] + 2C1*2/e/B^2 = (1+4C1/e)/B^2 :")
c2 = 1.0 + 4.0 * c1 / math.e
P(f"   coeff = {c2:.6f}  => Var(w) <= {c2:.4f}/B^2")
# verify E[w^2] numerically at sample B
P("   B      E[w^2]*B^2(实测)   (须 <= %.3f)" % c2)


def ew2(B):
    def f0(w):
        return phi(w) * math.exp(-B * w)
    def f2(w):
        return w * w * phi(w) * math.exp(-B * w)
    Z0 = quad(f0, 0, 2, epsabs=1e-15, epsrel=1e-14, limit=400)[0]
    Z2 = quad(f2, 0, 2, epsabs=1e-15, epsrel=1e-14, limit=400)[0]
    return Z2 / Z0 * B * B


for B in (1.0, 2.0, 5.0, 20.0, 100.0):
    P(f"   {B:5.1f}   {ew2(B):12.6f}")
P("")
P("   统一常数: C_rig(B) = min( 0.2971*B^2 [range bound],  min_w1[30 L^2 + 4*C1*B^2 e^{-B w1}] )")
rng_c = 9.0 * (0.5 - 1.0 / math.pi) ** 2
P(f"   range-bound 系数 9*(1/2-1/pi)^2 = {rng_c:.6f}  (Var(f) <= (fmax-fmin)^2/4)")
rows = []
for B in [10 ** x for x in np.linspace(-2, 4, 61)]:
    cb = rng_c * B * B
    bestc = min(29.95 * abs(dphi(w1)) ** 2
                + 9.0 * B * B * min(1.0, c1 * math.exp(-B * w1))
                for w1 in np.linspace(0.05, 1.8, 60))
    rows.append((B, cb, bestc, min(cb, bestc)))
P("   B      C_range      C_tail      C_used")
for B, cb, ct, cu in rows[::6]:
    P(f"   {B:8.3g}  {cb:10.3f}  {ct:10.3f}  {cu:9.3f}")
P(f"   -> C_rig(uniform) = max_B min(...) = {max(r[3] for r in rows):.3f}")
P("")

# ---- [P5] indefinite Hessian on the slice
P("[P5] 切片上 Hess V 的特征值 (V = beta*c*(1-cos u) - 2 log(sin u/u), 2D 径向)")
P("   列 = 径向/角向特征值的下确界 (含 Phi 部分);  末列 = 仅 S_W 部分 (对照)")
P("   theta      beta    min eig_rad    min eig_ang   (仅 S_W)")


def eig_tot(B):
    us = np.linspace(1e-5, math.pi - 1e-5, 400001)
    sw = B * np.cos(us)
    ang = B * np.sin(us) / us
    # Phi = -2 log(sin u/u);  psi'' = 2(csc^2 u - 1/u^2), psi'/u = 2(1/u^2 - cot u/u)
    ps2 = 2.0 * (1.0 / np.sin(us) ** 2 - 1.0 / us ** 2)
    ps1u = 2.0 * (1.0 / us ** 2 - np.cos(us) / (np.sin(us) * us))
    rad_tot = 0.25 * (sw + ps2)
    ang_tot = 0.25 * (ang + ps1u)
    return rad_tot.min(), ang_tot.min(), (0.25 * sw).min()


for th, be in ((math.pi / 2, 1.0), (math.pi / 2, 10.0), (math.pi / 2, 1000.0),
               (3 * math.pi / 4, 10.0)):
    B = be * math.cos(th / 2)
    e1, e2, e3 = eig_tot(B)
    P(f"   {th/math.pi:6.3f}pi  {be:7.1f}   {e1:12.6f}  {e2:12.6f}   (仅 S_W: {e3:.6f})")
P("   => 总 Hess V 的径向/角向下确界 >=0 只在小 B 成立; 大 B 时径向变负 (u>pi/2).")
P("      故切片上无 beta-一致 BL 常数: lambda_cond 只能给同 beta 的界 (与 O05 sec.2 一致).")
P("")

# ---- [P6] error bars
P("[P6] 误差棒 (双容差)")
for B in (1.0, 10.0, 1000.0):
    a = Psi(B, 1e-11)
    b = Psi(B, 1e-15)
    P(f"   B={B:8.1f}  Psi(1e-11)={a:.14f}  Psi(1e-15)={b:.14f}  diff={abs(a-b):.2e}")
P("")
P("[DONE]")

with open(r"D:\Dev\math\research\v3\w7\o20_varbound.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG) + "\n")
