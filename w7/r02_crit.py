# -*- coding: utf-8 -*-
"""
r02_crit.py — q04 判据 2*beta*t^2*Var(w) <= c*<w> 的解析证明核对
w(rho)=sin(rho/2)/rho, t=sin(th/2), c=cos(th/2)
mu ∝ exp(-beta(1-c*cos(rho/2))) sin^2(rho/2)/rho drho, rho in [0,22pi]

解析结论（本脚本做数值核对，见 R02-CRIT.md）:
  (1) beta=0: Z=<分母>, A=Z<w>, B=Z<w^2> 全部是 Cin/Ci 的显式组合
  (2) |w| <= 1/2  =>  w^2 <= |w|/2  =>  2 E[w^2] <= E|w|
  (3) 驼峰配对: 每个 pi-块与后一块关于 (2j+1)pi 镜像，exp 因子与 sin^2 精确对称，
      1/rho^2 严格递降 => E|w| <= <w>
  => 2Var(w) <= <w> 无条件 [T]; 判据在 beta t^2 <= c 下 [T]
纯 numpy/scipy，单进程。
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import sici
from math import pi, log, sin, cos, sqrt

LAM = 22*pi
GAM = 0.5772156649015328606
CUTS = [2*pi*k for k in range(1, 12)]        # 11 个 sin^2 周期边界


def quad1(f, Lam=LAM):
    return quad(f, 0.0, Lam, epsabs=1e-15, epsrel=1e-13,
                limit=4000, points=CUTS)[0]


def mom(beta, c, Lam=LAM):
    """返回 (Z, A=Z<w>, B=Z<w^2>, C=Z<|w|>)"""
    def f(r):
        return np.exp(-beta*(1.0-c*np.cos(0.5*r)))*np.sin(0.5*r)**2/r
    def w(r):
        return np.sin(0.5*r)/r
    Z = quad1(f, Lam)
    A = quad1(lambda r: f(r)*w(r), Lam)
    B = quad1(lambda r: f(r)*w(r)**2, Lam)
    C = quad1(lambda r: f(r)*abs(w(r)), Lam)
    return Z, A, B, C


def ci(x):
    return sici(x)[1]


def cin(x):
    return GAM + log(x) - ci(x)


def exact_b0(Lam=LAM):
    """beta=0 的 Cin/Ci 闭式（rho 截断 Lam, X=Lam/2）"""
    X = 0.5*Lam
    Z0 = 0.5*cin(Lam)
    A0 = 0.5*(-sin(X)**3/X + 0.75*(log(3.0) + ci(X) - ci(3*X)))
    B0 = 0.25*(-sin(X)**4/(2*X*X) - 4*sin(X)**3*cos(X)/(2*X)
               + (log(2.0) + ci(2*X) - ci(4*X)))
    return Z0, A0, B0


def bessel_series(beta, c, K=60, Lam=LAM):
    """Z = e^{-beta} sum_k (beta c)^k/k! M_k,  M_k = ∫_0^{Lam/2} cos^k x sin^2x/x dx"""
    X = 0.5*Lam
    cuts = [pi*k for k in range(1, int(X/pi)+1)]
    M = np.empty(K+1)
    for k in range(K+1):
        M[k] = quad(lambda x: cos(x)**k*sin(x)**2/x, 0.0, X, epsabs=1e-15,
                    epsrel=1e-13, limit=4000, points=cuts)[0]
    from math import lgamma
    k = np.arange(K+1)
    coef = np.exp(k*np.log(beta*c) - np.array([lgamma(i+1.0) for i in k]))
    return np.exp(-beta)*float(np.sum(coef*M))


def w_dens(beta, c, Lam=LAM):
    def f(r):
        return np.exp(-beta*(1.0-c*np.cos(0.5*r)))*np.sin(0.5*r)**2/r
    Z = quad1(f, Lam)
    A = quad1(lambda r: f(r)*np.sin(0.5*r)/r, Lam)
    B = quad1(lambda r: f(r)*(np.sin(0.5*r)/r)**2, Lam)
    C = quad1(lambda r: f(r)*abs(np.sin(0.5*r)/r), Lam)
    mean = A/Z
    m2 = B/Z
    return mean, m2, m2-mean*mean, C/Z


def pair_nets(beta, c, Lam=LAM):
    """每个 [2jpi,(2j+2)pi] 内 net = ∫ f w, |net| 供 E|w|<Ew 核对"""
    out = []
    for j in range(5):
        a, b = 2*pi*j, 2*pi*(j+1)
        if b > Lam:
            break
        fw = quad(lambda r: np.exp(-beta*(1-c*np.cos(0.5*r)))*np.sin(0.5*r)**2/r*np.sin(0.5*r)/r,
                  a, b, epsabs=1e-15, epsrel=1e-13, limit=2000, points=[2*pi*j+pi])[0]
        fa = quad(lambda r: np.exp(-beta*(1-c*np.cos(0.5*r)))*np.sin(0.5*r)**2/r*abs(np.sin(0.5*r))/r,
                  a, b, epsabs=1e-15, epsrel=1e-13, limit=2000, points=[2*pi*j+pi])[0]
        out.append((j, fw, fa, fw-fa))
    return out


def crit(beta, th):
    c, t2 = cos(0.5*th), sin(0.5*th)**2
    mean, m2, var, absw = w_dens(beta, c)
    lhs, rhs = 2*beta*t2*var, c*mean
    return lhs, rhs, lhs/rhs if rhs != 0 else np.inf, mean, var, absw


def beta_c(th, lo=0.0, hi=40.0):
    """满足 2 beta t^2 Var <= c<w> 的 beta 最大值（二分；在 [lo,hi] 上单调性由数值观察）"""
    f = lambda b: crit(b, th)[2]-1.0
    if f(lo) > 0:
        return 0.0
    for _ in range(80):
        mid = 0.5*(lo+hi)
        if f(mid) <= 0:
            lo = mid
        else:
            hi = mid
    return lo


L = []
def say(s=""):
    print(s)
    L.append(s)


say("# R02-CRIT 数值核对日志")
say(f"Lambda = 22pi = {LAM:.15f}   py -3.12  numpy/scipy")
say("")

# [1] beta=0 闭式 vs 求积
say("[1] beta=0 闭式(Cin/Ci) vs 直接求积, Lam=22pi")
Z0, A0, B0 = exact_b0()
Zq, Aq, Bq, Cq = mom(0.0, 0.0)
say(f"  Z   exact {Z0:.15f}   quad {Zq:.15f}   diff {Z0-Zq:+.2e}")
say(f"  A   exact {A0:.15f}   quad {Aq:.15f}   diff {A0-Aq:+.2e}")
say(f"  B   exact {B0:.15f}   quad {Bq:.15f}   diff {B0-Bq:+.2e}")
m0, v0 = Aq/Zq, Bq/Zq-(Aq/Zq)**2
say(f"  beta=0: <w>={m0:.12f}  <w^2>={Bq/Zq:.12f}  Var={v0:.12f}  E|w|={Cq/Zq:.12f}")
say(f"          E|w|/<w>={Cq/Aq:.9f}   2Var/<w>={(2*v0)/m0:.9f}  R=<w>/Var={m0/v0:.9f}")
say("")

# [2] h1e 核准（Lam=2pi, beta=0）
say("[2] h1e 核准: Lam=2pi (rho<=2pi), beta=0/任意 c 等价")
Z2, A2, B2 = exact_b0(2*pi)
Zq2, Aq2, Bq2, Cq2 = mom(0.0, 0.0, 2*pi)
say(f"  <w> = {A2/Z2:.12f}  (h1e 0.3574113534)")
say(f"  Var = {B2/Z2-(A2/Z2)**2:.12f}  (h1e 0.0110592938)")
say(f"  quad 复核 <w>={Aq2/Zq2:.12f} Var={Bq2/Zq2-(Aq2/Zq2)**2:.12f}")
N1, N2, Z1 = 0.871244998311689, 0.676703157945306, 1.218826696528229
say(f"  h1e 组合式 N2/(4Z)-(N1/(2Z))^2 = {N2/(4*Z1)-(N1/(2*Z1))**2:.12f}")
say("")

# [3] Bessel 级数核对
say("[3] beta>0: Bessel 级数 Z=e^{-b}sum (bc)^k/k! M_k  vs 求积")
for beta, c in [(1.0, 0.7), (2.5, 0.5), (0.5, 0.99)]:
    Zs = bessel_series(beta, c)
    Zq3 = mom(beta, c)[0]
    say(f"  beta={beta:5.2f} c={c:4.2f}: series {Zs:.12f}  quad {Zq3:.12f}  diff {Zs-Zq3:+.2e}")
say("")

# [4] 无条件性网格
say("[4] 无条件性网格: 2Var<=<w>?  E|w|<=<w>?")
say("    beta      c      <w>         Var        E|w|      2Var/<w>   E|w|/<w>   R=<w>/Var")
worst2 = worstabs = -1e9
worst2at = worstabsat = None
for beta in [0.0, 0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.5, 4.0, 6.0, 10.0, 20.0]:
    for c in [0.0, 0.1, 0.25, 0.5, 0.7, 0.85, 0.95, 1.0]:
        mean, m2, var, absw = w_dens(beta, c)
        r2, ra = 2*var/mean, absw/mean
        if r2 > worst2:
            worst2, worst2at = r2, (beta, c)
        if ra > worstabs:
            worstabs, worstabsat = ra, (beta, c)
        say(f"  {beta:6.2f}  {c:5.2f}  {mean:.9f}  {var:.9f}  {absw:.9f}  {r2:.9f}  {ra:.9f}  {mean/var:10.4f}")
say(f"  worst 2Var/<w> = {worst2:.9f} at (beta,c)={worst2at}")
say(f"  worst E|w|/<w> = {worstabs:.9f} at (beta,c)={worstabsat}")
say("")

# [5] 配对核对
say("[5] 驼峰配对 net=∫fw>0 与 ∫fw-∫f|w|>0 (beta=0.3, c=0.5)")
for j, fw, fa, d in pair_nets(0.3, 0.5):
    say(f"  pair j={j} [{2*j}pi,{2*j+2}pi]: ∫fw={fw:+.12e}  ∫f|w|={fa:.12e}  net={d:+.12e}")
say("")

# [6] Lambda 依赖（渐近 ratio -> 2B/A）
say("[6] 截断依赖 beta=0: 2Var/<w>, E|w|/<w>, R")
for Lam in [2*pi, 6*pi, 10*pi, 22*pi, 42*pi, 82*pi]:
    Zx, Ax, Bx = exact_b0(Lam)
    mean = Ax/Zx
    var = Bx/Zx-mean*mean
    Cx = quad(lambda r: np.sin(0.5*r)**2/r*abs(np.sin(0.5*r))/r, 0, Lam,
              epsabs=1e-15, epsrel=1e-13, limit=4000,
              points=[2*pi*k for k in range(1, int(Lam/(2*pi))+1)])[0]
    say(f"  Lam={Lam/pi:5.1f}pi  <w>={mean:.9f} Var={var:.9f} 2Var/<w>={2*var/mean:.9f} E|w|/<w>={(Cx/Zx)/mean:.9f} R={mean/var:.6f}")
say(f"  渐近: A_inf=(3/8)ln3={0.375*log(3):.12f}  B_inf=(1/4)ln2={0.25*log(2):.12f}  2B/A={2*(0.25*log(2))/(0.375*log(3)):.9f}")
say("")

# [7] 判据边界 beta_c(theta)
say("[7] 判据 2 beta t^2 Var <= c<w> 的边界 beta_c(theta) vs q04 充分条件 c/t^2")
say("    theta/pi   c/t^2        beta_c       gain   (2Var/<w> at beta_c)")
rows = []
for th in [0.05*pi, 0.1*pi, 0.2*pi, 0.3*pi, 0.4*pi, 0.5*pi, 0.6*pi, 0.7*pi, 0.8*pi, 0.9*pi, 0.95*pi, 0.99*pi]:
    c, t2 = cos(0.5*th), sin(0.5*th)**2
    bc = beta_c(th)
    pred = c/t2
    _, _, _, mean, var, absw = crit(min(bc, 40.0), th)
    rows.append((th/pi, pred, bc, bc/pred if pred else np.nan, 2*var/mean))
    say(f"  {th/pi:8.3f}  {pred:11.4f}  {bc:11.4f}  {bc/pred:6.3f}   {2*var/mean:.9f}")
say("")

# [8] q04 条件成立时判据成立（抽查 + 反例区确认）
say("[8] q04 充分条件 beta t^2<=c 内：LHS/RHS 最大比 (应 <1)")
mx = -1; mxat = None
for th in np.linspace(0.02*pi, 0.98*pi, 49):
    c, t2 = cos(0.5*th), sin(0.5*th)**2
    for beta in np.linspace(0.0, 40.0, 81):
        if beta*t2 <= c:
            lhs, rhs, ratio = crit(beta, th)[:3]
            if ratio > mx:
                mx, mxat = ratio, (th/pi, beta, c, t2)
say(f"  max LHS/RHS = {mx:.9f} at (theta/pi,beta)={mxat[0]:.4f},{mxat[1]:.3f}")
say(f"  若 2Var<=<w> 无条件，则 ratio<=beta t^2/c<=1 [T]；数值 max {mx:.6f}<1 一致")
say("")
say("[9] theta=pi (c=0) 反例确认：判据必假，但 2Var<=<w> 仍真")
for beta in [0.5, 2.5, 10.0]:
    lhs, rhs, ratio, mean, var, absw = crit(beta, pi)
    say(f"  beta={beta:5.2f}: LHS={lhs:.9e} RHS={rhs:.1e} 2Var/<w>={2*var/mean:.6f}")
say("")
# [10] h1e 式配对核验证: A-2B = (1/2)∫ e sin^3 s (s-sin s)/s^3 ds >= 0
say("[10] 核心恒等式 A-2B = (1/2)∫ e sin^3 s (s-sin s)/s^3 ds 与配对非负性")
def hkern(beta, c, Lam=LAM):
    """按 s 变量: h = e sin^3 s (s-sin s)/s^3, 块 k=[2pi k,2pi k+2pi] 内配对非负"""
    out = []
    kmax = int(np.floor(Lam/(2*pi)))
    for k in range(kmax):
        a, b = 2*pi*k, 2*pi*(k+1)
        if b <= Lam:
            val = quad(lambda s: np.exp(-beta*(1-c*np.cos(s)))*sin(s)**3*(s-sin(s))/s**3,
                       a, b, epsabs=1e-16, epsrel=1e-13, limit=3000,
                       points=[a+pi])[0]
            out.append((k, a/pi, b/pi, val, "full"))
    if kmax*2*pi < Lam:
        a, b = 2*pi*kmax, Lam
        val = quad(lambda s: np.exp(-beta*(1-c*np.cos(s)))*sin(s)**3*(s-sin(s))/s**3,
                   a, b, epsabs=1e-16, epsrel=1e-13, limit=3000)[0]
        out.append((kmax, a/pi, b/pi, val, "tail"))
    return out

for beta, c in [(0.0, 0.0), (0.3, 0.5), (2.5, 0.9), (10.0, 1.0)]:
    Zx, Ax, Bx, Cx = mom(beta, c)
    hid = 0.5*quad(lambda s: np.exp(-beta*(1-c*np.cos(s)))*sin(s)**3*(s-sin(s))/s**3,
                   0, 0.5*LAM, epsabs=1e-16, epsrel=1e-13, limit=8000,
                   points=[pi*k for k in range(1, 6)])[0]
    say(f"  beta={beta:5.2f} c={c:4.2f}: A-2B={Ax-2*Bx:+.12e}  half-integral={hid:+.12e}  diff={Ax-2*Bx-hid:+.2e}  A={Ax:.9f} 2B={2*Bx:.9f} 2B/A={2*Bx/Ax:.9f}")
say("  逐块(h 在 s 变量, 块 [2kpi,2(k+1)pi], 末块为半个正驼峰; 每块应 >=0):")
for beta, c in [(0.3, 0.5), (10.0, 1.0)]:
    say(f"   beta={beta}, c={c}:")
    for k, a, b, val, kind in hkern(beta, c):
        say(f"     k={k} s/pi in [{a:.1f},{b:.1f}] {kind:4s}  ∫h={val:+.12e}")
say("")
say("[11] 2B/A 网格最大（应 <1，即 2<w^2><=<w> [T]）")
mx2 = -1; at2 = None
for beta in [0.0, 0.5, 2.5, 6.0, 20.0, 40.0, 100.0]:
    for c in [0.0, 0.3, 0.7, 0.95, 1.0]:
        Zx, Ax, Bx, Cx = mom(beta, c)
        r = 2*Bx/Ax
        if r > mx2:
            mx2, at2 = r, (beta, c)
say(f"  max 2B/A = {mx2:.9f} at (beta,c)={at2};  渐近 beta->inf (w->1/2) 与 beta=0,Lam->inf 分别给 1 与 {2*0.25*log(2)/(0.375*log(3)):.9f}")
say("")
say("DONE")

with open(r"D:\Dev\math\w7\r02_crit.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
