# -*- coding: utf-8 -*-
"""
p05_cinf.py -- P05: E[H_vv + Phi''|0] >= c * beta * ||v||^2, 显式 c
   单 plaquette SU(2), 边界 holonomy theta, 切片 xi_3=0, v=e_3, a=0, ||v||=1
   约定与 O21-O21-CURLIN.md / O20-VARBOUND.md / O15-THETAPI.md 逐项一致

精确约化 (O20 命题 1, [T]):
   u = rho/2 in (0,pi),   nu_B ∝ sin^2 u * e^{-B(1-cos u)} du/u,   B = beta*cos(theta/2)
   g(u) = sin u/(2u)  (= O15 的 A_0, 逐点 > 0, g(0)=1/2)
   Q := H_vv/beta = cos(theta/2) * g(u)                      [beta-无关, 逐点 >= 0 当 theta<=pi]
   G(B) := E_{nu_B}[g] = <A_0>      (0 < G < 1/2,  随 B 单调 up)
   Phi''_t(u) = 1/(2u^2) - cot(u)/(2u) = 1/6 + u^2/90 + O(u^4)  (点态 >= 1/6)
   F(B) := E_{nu_B}[Phi''_t] = <Phi''|0>,  >= 1/6,  随 B 单调 down  ->  1/6
因此
   E[H_vv|0]     = beta * (c/2) * G(B)                        (c := cos(theta/2))
   E[Phi''|0]    = F(B)
   c_tot(beta,theta) = (c/2) G(B) + F(B)/beta ,   B = beta*c
   c_inf(theta)  = lim_{beta->inf} c_tot = c/4  = cos(theta/2)/4 = Re tr U_p/(4N)
"""
import math
import numpy as np

LOG = []
def P(s=""):
    print(s)
    LOG.append(s)

NG = 24
XP, WP = np.polynomial.legendre.leggauss(NG)
SQ2 = math.sqrt(2.0)


def _nodes(npan):
    edges = np.linspace(0.0, SQ2, npan + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    half = 0.5 * (edges[1:] - edges[:-1])
    z = (mid[:, None] + half[:, None] * XP[None, :]).ravel()
    w = (half[:, None] * WP[None, :]).ravel()
    return z, w


_Z, _W = _nodes(600)


def _g(u):
    return np.where(u < 1e-8, 0.5, np.sin(u) / (2.0 * np.where(u < 1e-8, 1.0, u)))


def _phitt(u):
    small = u < 1e-3
    us = np.where(small, 1.0, u)
    exact = 0.5 / (us * us) - 1.0 / (2.0 * us * np.tan(us))
    ser = 1.0 / 6.0 + u * u / 90.0 + u ** 4 / 945.0
    return np.where(small, ser, exact)


def _moments(B, z=None, w=None):
    """返回 (Z0, G=<A0>, F=<Phi''>) 在 nu_B 下。z,w 可换更高分辨率做收敛检验。"""
    z = _Z if z is None else z
    w = _W if w is None else w
    wz = z * z
    u = np.arccos(np.clip(1.0 - wz, -1.0, 1.0))
    jac = 2.0 * wz * np.sqrt(np.maximum(0.0, 2.0 - wz)) / u   # phi(w) dw 在 z 下
    ex = -B * wz
    dens = jac * np.exp(ex - ex.max())   # 平移防溢出 (比值不变)
    Z0 = float(np.sum(w * dens))
    G = float(np.sum(w * dens * _g(u))) / Z0
    F = float(np.sum(w * dens * _phitt(u))) / Z0
    return Z0, G, F


def G_of(B):
    return _moments(B)[1]


def F_of(B):
    return _moments(B)[2]


def c_tot(beta, theta):
    c = math.cos(0.5 * theta)
    B = beta * c
    _, G, F = _moments(B)
    return 0.5 * c * G + F / beta, G, F


# ---------------------------------------------------------------- main
P("=" * 84)
P("p05_cinf.log -- P05: E[H_vv+Phi''] >= c*beta*||v||^2  (单 plaquette SU(2), 切片 xi_3=0, v=e_3)")
P("=" * 84)
P("约化: u=rho/2; nu_B ∝ sin^2 u e^{-B(1-cos u)} du/u;  B=beta*cos(theta/2)")
P("      Q=H_vv/beta = c*g(u), g=sin u/(2u)>0;  G=<A0>=E[g] in (0,1/2);  F=<Phi''> >= 1/6")
P("      c_tot = (c/2)G(B) + F(B)/beta ;   c_inf = c/4 = cos(theta/2)/4")
P("")

# ---- [P0] 收敛检验 + 与 o20/o15/o21 的存档对账
P("[P0] 求积收敛与存档对账 (npan 600 vs 1200, 目标: 与 O21/O15 的表面值逐位一致)")
_z2, _w2 = _nodes(1200)
rows = [(2.5, 0.0), (2.5, 0.5 * math.pi), (200.0, 0.0), (200.0, 0.75 * math.pi)]
P("   beta   theta/pi |   c_tot(npan=600)   c_tot(1200)     diff      |  E[H]=beta*(c/2)G    F=<Phi''>")
for b, th in rows:
    v1, G1, F1 = c_tot(b, th)
    c = math.cos(0.5 * th); B = b * c
    _, G2, F2 = _moments(B, _z2, _w2)
    v2 = 0.5 * c * G2 + F2 / b
    P("  %6.3f  %7.4f |  %.12f  %.12f  %8.1e |  %.9f      %.9f"
      % (b, th / math.pi, v1, v2, abs(v1 - v2), b * 0.5 * c * G1, F1))
P("")
P("   O21 存档: c(2.5,0)=0.291970  c(2.5,pi/2)=0.223273  c(2.5,3pi/4)=0.152250")
P("             c(200,0)=0.250418  c(200,3pi/4)=0.096090  c_inf(0)=0.25")
P("   O15 存档: <A0>(N=2,beta=2.5,theta=0)=0.442503 ; E[H](beta=2.5,theta=0)=0.553129 ;")
P("             E[H](beta=200,theta=0)=0.249584")
P("   本脚本同点: G(2.5)=%.9f  E[H](2.5,0)=%.9f  E[H](200,0)=%.9f"
  % (_moments(2.5)[1], c_tot(2.5, 0.0)[1] * 2.5 * 0.5, c_tot(200.0, 0.0)[1] * 200.0 * 0.5))
P("")

# ---- [P1] O21 的 c 表复算
P("[P1] c_tot(beta,theta) 复算 O21 表 (行 theta/pi, 列 beta)")
bs = [2.5, 6.0, 10.0, 50.0, 200.0]
ths = [0.0, 0.25, 0.5, 0.75]
P("   theta/pi |" + "".join("%12s" % ("b=%g" % b) for b in bs) + "   c_inf=c/4")
for t in ths:
    th = t * math.pi
    line = "   %7.4f |" % t
    for b in bs:
        line += "%12.6f" % c_tot(b, th)[0]
    line += "   %.6f" % (0.25 * math.cos(0.5 * th))
    P(line)
P("")

# ---- [P2] G(B), F(B): 单调性与量级  (B = beta*c 单参数族)
P("[P2] 单参数族 B=beta*c: G(B)=<A0> 与 F(B)=<Phi''>  (B=1.135 = beta_min*cos(0.35pi))")
P("        B        G(B)=<A0>      F(B)=<Phi''>     (1/2-G)*B      (F-1/6)*B")
Bs = [0.5, 0.8, 1.0, 1.135, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0, 25.0, 50.0,
      100.0, 200.0, 500.0, 1000.0, 1e4]
Gs, Fs = [], []
for B in Bs:
    _, G, F = _moments(B)
    Gs.append(G); Fs.append(F)
    P("   %7.4g   %.10f    %.10f    %+.6f      %+.6f"
      % (B, G, F, (0.5 - G) * B, (F - 1.0 / 6.0) * B))
dG = np.diff(Gs); dF = np.diff(Fs)
P("   单调性(G 增 / F 减):  max(dG)=%+.3e (应 >0) ; max(dF)=%+.3e (应 <0) ; G<1/2: max G=%.10f"
  % (dG.max(), dF.max(), max(Gs)))
P("   B->inf 渐近: (1/2-G)B -> ?  (F-1/6)B -> ?   (Rayleigh 预测: 1/12 与 1/45)")
P("")

# ---- [P3] c_tot 对 beta 的单调性 (窗口 [2.5, 1e4])
P("[P3] c_tot(beta,theta) 在 beta>=2.5 上单调递减? (max 增量, 应 <=0)")
bgrid = np.exp(np.linspace(math.log(2.5), math.log(1e4), 40))
worst = -1e9
for t in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9]:
    th = t * math.pi
    vals = [c_tot(b, th)[0] for b in bgrid]
    inc = max(vals[i + 1] - vals[i] for i in range(len(vals) - 1))
    worst = max(worst, inc)
    P("   theta/pi=%.2f  c(2.5)=%.8f  c(1e4)=%.8f  c_inf=%.8f  max_inc=%+.2e"
      % (t, vals[0], vals[-1], 0.25 * math.cos(0.5 * th), inc))
P("   => 窗口内 max 增量 = %+.2e  (负/零 === 单调递减, inf 在 beta=inf 取到)" % worst)
P("")

# ---- [P4] c_min over theta in [0,0.7pi], beta in [2.5, inf)
P("[P4] c_min = inf_{theta<=0.7pi} inf_{beta>=2.5} c_tot   (稠密 theta 网格)")
tg = np.linspace(0.0, 0.7 * math.pi, 141)
for bmax, tag in [(200.0, "beta in [2.5,200]"), (1e4, "beta in [2.5,1e4]"), (float("inf"), "beta->inf (解析)")]:
    best = (1e9, None, None)
    for th in tg:
        if bmax == float("inf"):
            v = 0.25 * math.cos(0.5 * th)
        else:
            v = c_tot(bmax, th)[0]
        if v < best[0]:
            best = (v, th, bmax)
    P("   %-22s min = %.9f  at theta/pi = %.4f, beta = %s" % (tag, best[0], best[1] / math.pi, best[2]))
P("   解析极限: cos(0.7pi/2)/4 = cos(0.35pi)/4 = %.9f" % (0.25 * math.cos(0.35 * math.pi)))
for t in [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]:
    P("     theta/pi=%.2f :  c_inf=c/4 = %.9f ;  c(2.5)=%.9f ;  c(200)=%.9f"
      % (t, 0.25 * math.cos(0.5 * t * math.pi), c_tot(2.5, t * math.pi)[0], c_tot(200.0, t * math.pi)[0]))
P("")

# ---- [P5] 可证 ([T]) 的 beta-一致 c: 只用 (i) 精确约化 (ii) G 对 B 单调 (iii) F>=1/6
P("[P5] [T] beta-一致下界: c_tot >= (c/2) G(beta_min*c)   [丢掉 F/beta>0, 用 G 单调]")
th_max = 0.7 * math.pi
c_max = math.cos(0.5 * th_max)
Bmin = 2.5 * c_max
G_Bmin = _moments(Bmin)[1]
P("   theta_max = 0.7pi  => c_min_theta = cos(0.35pi) = %.9f" % c_max)
P("   B_min = beta_min*c = 2.5*%.9f = %.9f  =>  G(B_min) = <A0> = %.9f" % (c_max, Bmin, G_Bmin))
P("   [T] c >= (c/2)*G(B_min) = %.9f * %.9f = %.9f"
  % (0.5 * c_max, G_Bmin, 0.5 * c_max * G_Bmin))
P("   加上 F>=1/6 的窗口项 (只在 beta<=B_win 有效): c >= (c/2)G(beta_min c) + N/(12 B_win)")
for Bw in [2.5, 6.0, 10.0, 200.0]:
    P("      B_win=%7.1f :  c >= %.9f" % (Bw, 0.5 * c_max * G_Bmin + 2.0 / (12.0 * Bw)))
P("   O21 的纯窗口界 (不用 G): c >= N/(12 B_win)  ->  %.6f (B=2.5) / %.6f (B=6)"
  % (2.0 / (12 * 2.5), 2.0 / (12 * 6)))
P("")

# ---- [P6] 符号: theta 阈值, 两支
P("[P6] c_tot 的零点 theta_c^c(beta)  (c>0 的 theta-范围; 用 Phi'' 的 1/6 底救回来的部分)")
for b in [2.5, 6.0, 200.0]:
    lo, hi = 0.0, 2 * math.pi
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if c_tot(b, mid)[0] > 0:
            lo = mid
        else:
            hi = mid
    P("   beta=%6.1f :  theta_c^c/pi = %.6f  (c=cos(theta/2)=%+.6f)"
      % (b, 0.5 * (lo + hi) / math.pi, math.cos(0.5 * 0.5 * (lo + hi))))
P("   对照 O15 的 E[H|0] 零点 (N=2 精确 theta=pi; N=3,beta=2.5 时 0.9333pi)")
P("   theta>pi 支: Q = c*g < 0 逐点  =>  E[H_vv|0] = beta*(c/2)G < 0, 随 beta 线性发散(不是更快)")
for t in [1.0, 1.25, 1.5]:
    th = t * math.pi
    P("     theta/pi=%.2f : E[H_vv]/beta = %.6f (beta=2.5) / %.6f (beta=200)   [c=%.4f]"
      % (t, c_tot(2.5, th)[1] * 0.5 * math.cos(0.5 * th),
         c_tot(200.0, th)[1] * 0.5 * math.cos(0.5 * th), math.cos(0.5 * th)))
P("")

# ---- [P7] 判据合成
P("[P7] 判据合成")
cmin_inf = 0.25 * math.cos(0.35 * math.pi)
cmin_200 = min(c_tot(200.0, th)[0] for th in tg)
P("   (a) beta-一致负下界: E[H_vv|0] = beta*(c/2)G(B) ;  theta<=pi 时 >=0 逐点(模型内 [T])")
P("       => C=0, c 的下界由 Phi'' 与 G<1/2 的亏空竞速决定: c_tot = c/4 + 1/(12 beta) + O(beta^-2)")
P("   (b) 无下界? 否. theta<=pi: 有下界 0;  theta>pi: 线性 -|c|*G/2*beta, 有显式速率, 非超线性")
P("   (c) <A0> 变负? 否: A0=g(u)>0 逐点, 与 theta,beta 无关 (O15 §4.1 同).")
P("       E[H|0] 变负由 c=cos(theta/2) 单独决定 (N=2: theta=pi 精确)")
P("   显式数: c_min(theta<=0.7pi, beta>=2.5) = %.9f ;  beta<=200: %.9f ;  [T] 粗界 %.9f"
  % (cmin_inf, cmin_200, 0.5 * c_max * G_Bmin))
P("")
# ---- [P8] theta>pi 支: 集中点从 rho=0 移到 rho=2pi (c<0, B=beta*c<0) -> 速率不是线性
P("[P8] theta>pi 支 (c=cos(theta/2)<0): E[H_vv]/beta = (c/2)G(B), B=beta*c<0;  G(B)->?")
P("     => beta->inf 时测度集中到 u->pi (rho->2pi), 那里 g->0, 故 E[H_vv] 不随 beta 线性发散")
P("        beta    theta/pi=1.25 (c=-0.38268)  theta/pi=1.50 (c=-0.70711)  theta/pi=1.75 (c=-0.92388)")
bl = [2.5, 10.0, 100.0, 1000.0, 1e4]
for b in bl:
    line = "     %8.1f " % b
    for t in [1.25, 1.5, 1.75]:
        th = t * math.pi
        cc = math.cos(0.5 * th)
        _, G, F = _moments(b * cc)
        line += "   E[H]=%+11.6f (c_tot=%+.6f)" % (b * 0.5 * cc * G, 0.5 * cc * G + F / b)
    P(line)
P("     斜率 dlog|E[H_vv]|/dlog(beta):")
for t in [1.25, 1.5, 1.75]:
    th = t * math.pi
    cc = math.cos(0.5 * th)
    e1 = 10.0 * 0.5 * cc * _moments(10.0 * cc)[1]
    e2 = 1e4 * 0.5 * cc * _moments(1e4 * cc)[1]
    P("       theta/pi=%.2f : %.6f  (线性发散会是 1.0; sqrt(beta) 是 0.5)" % (t, math.log(abs(e2 / e1)) / math.log(1e3)))
P("     c_tot 在该支上的零点 beta_c (c_tot 由 + 穿到 -, 因 (c/2)G~-C*beta^-1/2 终将压过 F/beta~1/(6beta)):")
for t in [1.05, 1.25, 1.5, 1.75]:
    th = t * math.pi
    lo, hi = 1.0, 1e4
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if c_tot(mid, th)[0] > 0:
            lo = mid
        else:
            hi = mid
    P("       theta/pi=%.2f : beta_c = %8.2f   (c_tot(2.5)=%+.6f, c_tot(1e4)=%+.3e)"
      % (t, math.sqrt(lo * hi), c_tot(2.5, th)[0], c_tot(1e4, th)[0]))
P("     => O21 §1 的 '一个 theta>pi 支上任何 c>0 都不成立' 正确; 但 'E[H+Phi''] 发散' 的速率是")
P("        ~beta^(1/2) (不是 beta), 且 c_inf = cos(theta/2)/4 这只对 theta<pi 成立 (集中点在 rho=0)")
P("")
P("done.")

with open("p05_cinf.log", "w", encoding="utf-8") as f:
    f.write("\n".join(LOG) + "\n")
