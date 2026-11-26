"""p13_acid.py -- P13 酸测: O22-THRESHOLD 在自由无质量场 (S_W -> 0, beta -> 0) 上的行为

模型 (O05/O20/O21 逐字, 单 plaquette SU(2), 切片 xi_3=0, 边界 holonomy theta):
    u = r/2 in (0,pi),  rho_B(u) du  ∝  sin^2(u) exp(-B(1-cos u)) du / u ,   B := beta*cos(theta/2)
    g(u) = sin u/(2u),  f = G_v = beta*sin(theta/2)*g(u)
    Var(f) = beta^2 sin^2(theta/2) Var_B(g) = V_inf * Psi(B),  Psi := 36 B^2 Var_B(g)
    H_vv  = (B/2) g(u)          [O20 [P5] 的 ang 结构: 0.25*(B sin u/u + psi'/u), sin u/u = 2g]
    Phi'' = 1/(2u^2) - cot(u)/(2u)  >= 0   [O20 sec.1]
      => 自由端: H ∝ B -> 0;  Phi'' -> 1/6 有 beta-无关地板;  Var -> 0 (∝ B^2)
    H1(M): Var(G|0) <= E[H_vv|0] + W_full - delta,  切片上 Phi'==0 => W_full = E[Phi'']

判据 (项目宪章): 工具必须在自由无质量场上失效, 否则 vacuous.
本脚本量化: 自由场上 V_inf / c_inf / beta_0 / 前提 / 结论 各是什么.
"""
import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


# ------------------------------------------------------------------ moments
def Phi2(u):
    """Phi'' = 2/r^2 - cot(r/2)/r,  r=2u;  小 u 用级数 1/6 + u^2/90 + O(u^4)."""
    if u < 5e-2:
        return 1.0 / 6.0 + u * u / 90.0
    return 1.0 / (2.0 * u * u) - (math.cos(u) / math.sin(u)) / (2.0 * u)


def gfun(u):
    return math.sin(u) / (2.0 * u) if u > 0 else 0.5


def mom(B, h, eps=1e-14):
    """E_{nu_B}[h],  dnu_B ∝ sin^2 u exp(-B(1-cos u)) du / u."""
    def f(u):
        s = math.sin(u)
        return s * s * math.exp(-B * (1.0 - math.cos(u))) / u * h(u)
    pts = {0.0, math.pi / 2, math.pi}
    if B > 1.0:
        pts.add(min(math.pi / 2, 12.0 / math.sqrt(B)))
    pts = sorted(p for p in pts if 0.0 <= p <= math.pi)
    tot = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        if b - a > 1e-12:
            tot += quad(f, a, b, epsabs=eps, epsrel=1e-13, limit=500)[0]
    return tot


def stats(B):
    Z0 = mom(B, lambda u: 1.0)
    Z1 = mom(B, gfun)
    Z2 = mom(B, lambda u: gfun(u) ** 2)
    m1, m2 = Z1 / Z0, Z2 / Z0
    return m1, m2 - m1 * m1, mom(B, Phi2) / Z0   # E[g], Var(g), E[Phi'']


def Psi(B):
    return 36.0 * B * B * stats(B)[1]


THETAS = [math.pi / 2, 3 * math.pi / 4, 7 * math.pi / 8]
BETAS = [0.0, 1e-6, 1e-3, 1e-2, 0.05, 0.1571348, 0.5, 2.5, 6.0, 50.0, 1e3, 1e4]

P("=" * 100)
P("p13_acid.log -- P13 酸测: O22-THRESHOLD 在自由无质量场 (S_W->0, beta->0) 上")
P("=" * 100)
P("模型: u=r/2, rho_B ∝ sin^2 u e^{-B(1-cos u)}/u; B=beta*cos(theta/2); g=sin u/(2u)")
P("      f=beta*sin(theta/2)*g;  Var(f)=V_inf*Psi(B);  H_vv=(B/2)g;  Phi''>=0")
P("      H1(M): Var <= E[H_vv]+E[Phi''] (切片上 Phi'==0 故 W_full=E[Phi''])")
P("")

# ---- [P0] 锚: 与 O22 sec.5 / H1A 的存档值对账 (theta=pi <=> B=0 <=> Haar 律)
P("[P0] 锚点对账 (theta=pi 即 B=0, Haar 律; 存档: Var/beta^2=0.0110593, E[Phi'']=0.2051153)")
m1h, w0, ef0 = stats(1e-14)
P(f"   Var_Haar(g) = {w0:.10f}   (存档 V0 = 0.0110593)")
P(f"   E_Haar[Phi''] = {ef0:.10f}   (存档 W  = 0.2051153)")
P(f"   E_Haar[g]    = {m1h:.10f}   (< 1/2, 因 g(u)=sin u/(2u) <= 1/2, 等号仅 u=0)")
P("")

# ---- [P1] 自由端: Psi(B) -> 0
P("[P1] 自由端 Psi(B)=Var/V_inf (B->0 = 自由场);  B=beta*cos(theta/2)")
P("   B            Psi(B)        1-Psi        (Var/V_inf)")
for B in (1e-6, 1e-3, 1e-1, 1.0, 1e1, 1e2, 1e3, 1e4):
    p = Psi(B)
    P(f"   {B:<12.4g} {p:.10f}   {1-p:11.3e}")
P("   => B->0: Psi->0  (自由场方差比 V_inf 小到 0, 上界不紧但成立)")
P("")

# ---- [P2] 三常数在自由场上
P("[P2] 自由场上的 V_inf / c_inf / beta_0 (=V_inf/c_inf=(1-x^2)/(9x^3), x=cos(theta/2))")
P("   theta/pi   V_inf        c_inf       beta_0      beta_0<=2.5?")
for th in THETAS:
    c, s = math.cos(th / 2), math.sin(th / 2)
    V = s * s / (36 * c * c)
    ci = c / 4.0
    b0 = V / ci
    P(f"   {th/math.pi:8.4f}   {V:.8f}   {ci:.8f}   {b0:10.5f}   {b0<=2.5}")
P("   (V_inf 是 beta->inf 的极限, c_inf 是 c(beta) 的 beta-一致下确界=beta->inf 极限)")
P("   (两者都与 beta 无关地是 theta 的常数 => 自由场不改它们的值)")
P("")

# ---- [P3] 关键: E[H]/beta vs c_inf  (定理第二条腿, 字面版: 只含 H)
P("[P3] 定理第二条腿 (字面 'E[H]>=c_inf*beta', H 只来自 S_W):  E[H]/beta = (c/2)E_B[g]")
P("   theta/pi  beta        B          E[g]        E[H]/beta   c_inf      比值     前提成立?")
for th in THETAS[:2]:
    c, s = math.cos(th / 2), math.sin(th / 2)
    ci = c / 4.0
    for be in BETAS:
        B = be * c
        m1, w, ef = stats(B)
        r = (B / 2.0) * m1 / be if be > 0 else float("nan")
        P(f"   {th/math.pi:7.3f}  {be:<10.4g} {B:<10.4g} {m1:.8f}  {r:11.6f}  {ci:.6f}  "
          f"{r/ci:8.4f}   {r>=ci}")
P("   => 每个有限 beta 上 E[H]/beta < c_inf (E[g]<1/2), 只在 beta->inf 取等")
P("")

# ---- [P4] 总曲率版 (O21 实际测的对象 E[H+Phi'']) 与 H1
P("[P4] 总曲率版 E[H+Phi''] 与 H1(M): 自由场上是否满足前提与结论?")
P("   theta/pi  beta     Var        E[H]+E[Phi'']  H1-slack   c(beta)=tot/beta  c(beta)>=c_inf?")
for th in THETAS:
    c, s = math.cos(th / 2), math.sin(th / 2)
    V = s * s / (36 * c * c)
    ci = c / 4.0
    for be in BETAS:
        B = be * c
        m1, w, ef = stats(B)
        var = V * Psi(B)
        tot = (B / 2.0) * m1 + ef
        cb = tot / be if be > 0 else float("inf")
        P(f"   {th/math.pi:7.3f}  {be:<8.4g} {var:10.3e}  {tot:11.6f}  {tot-var:10.3e}  "
          f"{cb:13.6f}    {cb>=ci}")
P("   => 自由端 (beta->0): Var->0, E[H]+E[Phi'']->0.2051, H1 以大余量成立; c(beta)->inf >= c_inf")
P("   => 前提与结论在自由场上都成立 => 该定理分不出 YM 与自由场")
P("")

# ---- [P5] 自由场的 '0/inf' 类比比
P("[P5] 若强行在自由端读同样的比值: V_free:=lim_{b->0}Var=0, c_free:=lim_{b->0}c(b)=+inf")
P("   => beta_0^free = 0/inf = 0 (不是 0/0);  自由场整体落在 'beta>=beta_0' 内")
for th in THETAS[:1]:
    c, s = math.cos(th / 2), math.sin(th / 2)
    P(f"   theta=pi/2: Var(1e-6)={s*s/36*Psi(1e-6*c):.3e}, c(1e-6)={((1e-6*c/2)*stats(1e-6*c)[0]+stats(1e-6*c)[2])/1e-6:.3e}")
P("")

# ---- [P6] beta_0 = 2.5 / 6 的分界 (复核 O22)
P("[P6] 复核 beta_0 分界 (O22 sec.0): (1-x^2)/(9x^3)=2.5, 6")
for target in (2.5, 6.0):
    f = lambda x: (1 - x * x) / (9 * x ** 3) - target
    x = brentq(f, 1e-3, 0.999999)
    P(f"   beta_0={target}: x=cos(theta/2)={x:.6f}  theta={2*math.acos(x)/math.pi:.6f} pi")
P("")
P("[DONE]")

with open(r"D:\Dev\math\research\v3\w7\p13_acid.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG) + "\n")
