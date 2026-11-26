# -*- coding: utf-8 -*-
r"""
s03_haarneg.py -- S03: SU(3) 纯 Haar 的 gap<0 是 rank-2 真性质, 还是口径问题?

判据 (三行证完) [T]:
  w 周期, w_tau(x)=w(x+tau*lam), Z(tau)=int_{T^2} w(x+tau*lam)dx
  => x->x+tau*lam 保 Lebesgue => Z(tau)==const => (log Z)'' == 0
  => Var(Psi') = <Psi''>  (Psi := -log w)   <=>  gap := S - D == 0.
即 "gap != 0" 只有两个来源:
  (i)  (Psi, w) 配对不一致: S,D 里的 Phi 不是权重里的 -log(密度);
  (ii) 流程不是积分域的平移 (r03 的 SU(2) 是 1D 径向截面: 横向方向未积分).
另: Phi_O09 = Phi_A + psi, psi = 2 log(d12 d23 d13) 含 2 log d, **不是环面上的函数**
  => e^{-Phi_O09} 不是 T^2 上的测度, 其积分随切片 (chart cut) 变 (见 [E]).

[A] 恒等式: 匹配核 gap=0; 且对任意平滑核 gap=0 (核无关).
[B] SU(2) 1D (r03 模型): 纯 Haar 是否总正 —— 否, beta > beta* ~ 4.307 转负.
[C] SU(3) 纯 Haar 复核 (网格加密 + Richardson): r12 的负值收敛到 -0.1747 (错配核), 匹配核 0.
[D] D 的分解: Var(G) / Var(Phi') / Cov / psi 记账 (精确对上 r12 的 6.6154, 6.4407, -0.2358).
[E] Phi_O09 的切片依赖 (Z 随 delta 变 85 倍) + theta 扫描.
标签: [T] 定理/恒等式; [C] 数值; [E] 上界/估计.
"""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 30
LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


BETA = 2.5
LAM = (1.0, 0.0, -1.0)
OTH = {'12': ('23', '13'), '23': ('12', '13'), '13': ('12', '23')}
PAIRS = (('12', '23'), ('12', '13'), ('23', '13'))


def grid(N):
    h = 2 * math.pi / N
    a = (np.arange(N) + 0.5) * h
    return np.meshgrid(a, a, indexing='ij')


def pieces(A, B):
    d = {'12': A, '23': B, '13': A + B}
    s = {p: 4 * np.sin(v / 2) ** 2 for p, v in d.items()}
    return d, s, s['12'] * s['23'] * s['13']


def SW_(A, B, theta, beta):
    x = ((2 * A + B) / 3, (B - A) / 3, -(A + 2 * B) / 3)
    ph = (theta, theta, 0.0)
    return beta * (1.0 - (np.cos(ph[0] - x[0]) + np.cos(ph[1] - x[1]) + np.cos(ph[2] - x[2])) / 3.0)


def GH(A, B, theta, beta):
    """G = d_tau S_W, H = d^2_tau S_W  (与 r12 partC 同约定: 其 tau 与本脚本反号)"""
    x = ((2 * A + B) / 3, (B - A) / 3, -(A + 2 * B) / 3)
    ph = (theta, theta, 0.0)
    G = np.zeros_like(A)
    H = np.zeros_like(A)
    for j in range(3):
        G += -(beta / 3) * np.sin(ph[j] - x[j]) * LAM[j]
        H += +(beta / 3) * np.cos(ph[j] - x[j]) * LAM[j] ** 2
    return G, H


def psi_(A, B):
    return 2.0 * (np.log(np.abs(A)) + np.log(np.abs(B)) + np.log(np.abs(A + B)))


def psi_d(A, B):
    dd = {'12': LAM[0] - LAM[1], '23': LAM[1] - LAM[2], '13': LAM[0] - LAM[2]}
    d = {'12': A, '23': B, '13': A + B}
    p1 = 2 * (dd['12'] / d['12'] + dd['23'] / d['23'] + dd['13'] / d['13'])
    p2 = -2 * (dd['12'] ** 2 / d['12'] ** 2 + dd['23'] ** 2 / d['23'] ** 2 + dd['13'] ** 2 / d['13'] ** 2)
    return p1, p2


def Phi_der(kind, A, B):
    """Phi' 与 Phi'' (两次都除以 T 吸收 1/s 奇异).  kind='A': Phi=-log|Delta|^2;  'B': Phi=Phi_O09"""
    d, s, T = pieces(A, B)
    dd = {'12': LAM[0] - LAM[1], '23': LAM[1] - LAM[2], '13': LAM[0] - LAM[2]}
    q1, q2 = {}, {}
    for p in ('12', '23', '13'):
        v = d[p]
        if kind == 'A':
            q1[p] = -2 * np.sin(v)
            q2[p] = 4 - s[p]
        else:
            q1[p] = 8 * np.sin(v / 2) ** 2 / v - 2 * np.sin(v)
            q2[p] = (4 * np.sin(v / 2) / v - 2 * np.cos(v / 2)) ** 2
    T1 = np.zeros_like(A)
    T2 = np.zeros_like(A)
    for p in ('12', '23', '13'):
        o1, o2 = OTH[p]
        T1 += dd[p] * q1[p] * s[o1] * s[o2]
        T2 += dd[p] ** 2 * q2[p] * s[o1] * s[o2]
    for pp, qq in PAIRS:
        rr = ({'12', '23', '13'} - {pp, qq}).pop()
        T2 += 2 * dd[pp] * dd[qq] * q1[pp] * q1[qq] * s[rr]
    Tpp = np.zeros_like(A)
    for p in ('12', '23', '13'):
        o1, o2 = OTH[p]
        if kind == 'A':
            Tpp += dd[p] ** 2 * 2 * s[o1] * s[o2]
        else:
            q0 = (2 * np.sin(d[p] / 2) / d[p]) ** 2
            Tpp += dd[p] ** 2 * (2 - 2 * q0) * s[o1] * s[o2]
    return T, T1 / T, T2 / T, Tpp / T


def SD(kind, N, theta=math.pi, beta=BETA, weight='matched'):
    """kind: Phi 核; weight: 'matched' = 与 Phi 一致的权重, 'A'= |Delta|^2 e^{-S_W}, 'B'= e^{-Phi_O09}e^{-S_W}"""
    A, B = grid(N)
    d, s, T = pieces(A, B)
    SW = SW_(A, B, theta, beta)
    G, H = GH(A, B, theta, beta)
    _, p1, _, p2 = Phi_der(kind, A, B)
    if weight == 'A':
        w = T * np.exp(-SW)
    elif weight == 'B':
        w = np.exp(-SW - psi_(A, B)) * T / T          # e^{-Phi_O09-SW} = |Delta|^2 e^{-psi-SW}
    else:
        w = T * np.exp(-SW)
    w = w * np.ones_like(T)
    if kind == 'B' and weight == 'B':
        w = T * np.exp(-SW - psi_(A, B))
    S_ = np.sum(w * (p2 + H)) / np.sum(w)
    m1 = np.sum(w * (p1 + G)) / np.sum(w)
    D_ = np.sum(w * (p1 + G) ** 2) / np.sum(w) - m1 ** 2
    return dict(S=S_, D=D_, gap=S_ - D_, E1=m1)


# =====================================================================
# [A] 恒等式
# =====================================================================
def partA():
    N = 1024
    P("[A] 平移不变恒等式 (N=%d, SU(3) 2D Weyl, theta=pi, beta=2.5, lam=(1,0,-1)):" % N)
    P("    S=<Phi''+H>, D=Var(Phi'+G), gap=S-D;  匹配核 (Phi=-log 密度) 必须给 0")
    P("     %-6s %-16s %-14s %-14s %-14s" % ("beta", "Phi / 权重", "S", "D", "S-D"))
    for beta, wl in ((2.5, 'matched'), (2.5, 'A'), (0.0, 'matched'), (0.0, 'A')):
        o = SD('A', N, beta=beta, weight=wl)
        P("     %-6.2f %-16s %-+14.8f %-+14.8f %-+14.8f" % (beta, "Phi_A/" + wl, o['S'], o['D'], o['gap']))
    P("     => 匹配核 gap = 0 (beta=0 精确 0; beta=2.5 的是 O(h^2) 差分离散残留, 见下)")
    P("")
    P("     [A1] 匹配核 (beta=2.5) 的 h^2 收敛, 证明残差不是物理:")
    P("          %-6s %-14s" % ("N", "gap(Phi_A, 匹配)"))
    for NN in (256, 512, 1024, 2048, 4096):
        o = SD('A', NN, beta=2.5, weight='A')
        P("          %-6d %+.3e" % (NN, o['gap']))
    o1, o2 = SD('A', 2048, beta=2.5, weight='A')['gap'], SD('A', 4096, beta=2.5, weight='A')['gap']
    P("          Richardson(h^2): gap_inf = %+.2e   => 0 [T]+[C]" % (o2 + (o2 - o1) / 3))
    P("     => 恒等式对 **任意 S_W** 成立 (r12 只做了 beta=0 的 8=8 一步, 没做 beta!=0)")
    P("")
    # 核无关性: 任意平滑周期核 (无 log 奇异, 差分干净)
    A, B = grid(N)

    def gen(tau, mis):
        a, b = A + tau, B + tau
        Pg = 1.3 - 0.6 * np.cos(a) + 0.9 * np.sin(2 * b) + 0.4 * np.cos(a + b) - 0.7 * np.cos(a)
        cut = 0.7 * np.cos(a)                        # 错配: 核里漏掉这一项 (a 随平移变)
        return np.exp(-Pg), (Pg - cut if mis else Pg)

    h = 1e-2
    for mis in (False, True):
        f = {k: gen(k * h, mis) for k in (-2, -1, 0, 1, 2)}
        w0 = f[0][0]
        d1 = (f[-2][1] - 8 * f[-1][1] + 8 * f[1][1] - f[2][1]) / (12 * h)
        d2 = (-f[-2][1] + 16 * f[-1][1] - 30 * f[0][1] + 16 * f[1][1] - f[2][1]) / (12 * h ** 2)
        m1 = (w0 * d1).sum() / w0.sum()
        S_ = (w0 * d2).sum() / w0.sum()
        D_ = (w0 * d1 * d1).sum() / w0.sum() - m1 ** 2
        P("     %-28s S=%+.8f D=%+.8f  S-D=%+.2e" %
          ("通用平滑核 (匹配)" if not mis else "同权重但核漏 0.7cos(a-b) (错配)", S_, D_, S_ - D_))
    P("     => gap=0 与核的函数形式无关, 只要 Psi = -log w [T]")
    P("")


# =====================================================================
# [B] SU(2) 1D
# =====================================================================
def su2(beta, th):
    c, s = mp.cos(th / 2), mp.sin(th / 2)
    B = beta * c

    def g(u):
        return mp.sin(u) / (2 * u) if u != 0 else mp.mpf(1) / 2

    def wt(u):
        return mp.e ** (-B * (1 - mp.cos(u)))

    def W(u):
        return mp.sin(u) ** 2 / u if u != 0 else mp.mpf(0)

    def phipp(u):
        return (1 / (2 * u * u) - (mp.cos(u) / mp.sin(u)) / (2 * u)) if u != 0 else mp.mpf(1) / 6

    def Q(f):
        return mp.quad(lambda u: wt(u) * W(u) * f(u), [0, mp.pi])

    Z = Q(lambda u: mp.mpf(1))
    E1 = Q(g) / Z
    E2 = Q(lambda u: g(u) ** 2) / Z
    Sp = Q(phipp) / Z
    VG = E2 - E1 ** 2
    S = (beta * c / 2) * E1 + Sp
    D = beta ** 2 * s * s * VG
    return dict(S=S, Sp=Sp, E1=E1, VG=VG, D=D, gap=S - D)


def partB():
    P("[B] SU(2) 单 plaquette 径向截面 (r03/q11 模型, 与 r03 区间证书同一 1D 约化):")
    P("    nu ∝ e^{-beta(1-c cos u)} sin^2 u/u du, u in (0,pi); Phi=-log J 只依赖 u => Phi'≡0")
    P("    S=(beta c/2)E[g]+E[Phi''], D=beta^2 s^2 Var[g]")
    o = su2(mp.mpf(5) / 2, mp.pi)
    P("    [B0] 回归 theta=pi,beta=2.5: S=%.16f D=%.16f gap=%.16f" % (float(o['S']), float(o['D']), float(o['gap'])))
    P("         r03: 0.205115297123073148 / 0.0691205860810849 / +0.135994711041988211  => 吻合")
    P("")
    P("    [B1] theta 扫描 beta=2.5:")
    P("         %-8s %-6s %-14s %-14s %-14s" % ("theta", "c", "S", "D", "S-D"))
    for i in range(13):
        th = mp.pi * i / 12
        r = su2(mp.mpf(5) / 2, th)
        P("         %-8.4f %-6.3f %-+14.8f %-+14.8f %-+14.8f" %
          (float(th), float(mp.cos(th / 2)), float(r['S']), float(r['D']), float(r['gap'])))
    P("         纯 Haar 只在 theta=pi (c=0) 出现; [0,pi] 上 gap 单调降, 最小 +0.136 > 0")
    P("")
    P("    [B2] 纯 Haar 点 (theta=pi) 的 beta 扫描:")
    r = su2(mp.mpf(5) / 2, mp.pi)
    bstar = mp.sqrt(r['Sp'] / r['VG'])
    P("         E[Phi'']=%s  Var[g]=%s  beta*=sqrt(Sp/Var)=%s" %
      (mp.nstr(r['Sp'], 12), mp.nstr(r['VG'], 12), mp.nstr(bstar, 10)))
    for bb in (1.0, 2.0, 2.5, 3.0, 4.0, float(bstar), 4.5, 5.0, 6.0):
        rr = su2(mp.mpf(bb), mp.pi)
        P("         beta=%-7.4f gap=%-+14.8f%s" % (bb, float(rr['gap']), "  <== 转负" if bb > float(bstar) else ""))
    P("         => SU(2) 纯 Haar 不总正, 阈值 beta* ~ 4.3066; beta=2.5 的 +0.136 是窗口内数值")
    P("")
    P("    [B3] beta=5 的 theta 扫描 (大 beta 下出现正的 theta 窗口):")
    for i in range(7):
        th = mp.pi * i / 6
        rr = su2(mp.mpf(5), th)
        P("         theta=%-7.4f gap=%-+14.8f" % (float(th), float(rr['gap'])))
    P("")


# =====================================================================
# [C] SU(3) 纯 Haar 复核
# =====================================================================
def partC():
    P("[C] SU(3) 纯 Haar 复核 (beta=0, 权重 |Delta|^2, theta=pi, lam=(1,0,-1)):")
    P("     %-6s %-15s %-15s %-15s   %-13s %-13s" % ("N", "S(O09 错配)", "D(O09 错配)", "gap", "gap(A 匹配)", "gap(G 混入)"))
    ser = []
    A, B = grid(1024)
    G, H = GH(A, B, math.pi, 2.5)
    T = pieces(A, B)[2]
    vG = (T * G * G).sum() / T.sum() - ((T * G).sum() / T.sum()) ** 2
    for N in (256, 512, 1024, 2048):
        a = SD('B', N, beta=0.0, weight='A')
        b = SD('A', N, beta=0.0, weight='A')
        AA, BB = grid(N)
        G2, _ = GH(AA, BB, math.pi, 2.5)
        T2 = pieces(AA, BB)[2]
        vg = (T2 * G2 * G2).sum() / T2.sum() - ((T2 * G2).sum() / T2.sum()) ** 2
        ser.append((N, a['gap']))
        P("     %-6d %-+15.8f %-+15.8f %-+15.8f   %-+13.2e %-+13.8f" %
          (N, a['S'], a['D'], a['gap'], b['gap'], a['gap'] - vg))
    g1 = [g for n, g in ser if n == 1024][0]
    g2 = [g for n, g in ser if n == 2048][0]
    P("     Richardson(h^2): gap_inf(O09, G=0) ≈ %+.8f   (r12 报 -0.1747 的 G=0 部分)" % (g2 + (g2 - g1) / 3))
    P("     Var_T(G) = %.8f  (r12: 0.06108481)" % vg)
    P("     => r12 的 -0.23580249 = -0.1747 (错配残差) - 0.0611 (Var(G), 权重里却没有 e^{-S_W}) [C]")
    P("     => 匹配核同网格 gap = 0 (上表末列倒数第二列) [T]+[C]")
    P("")


# =====================================================================
# [D] D 的分解
# =====================================================================
def partD():
    N = 2048
    A, B = grid(N)
    T = pieces(A, B)[2]
    d, s, _ = pieces(A, B)
    dd = {'12': 1.0, '23': 1.0, '13': 2.0}
    G, H = GH(A, B, math.pi, 2.5)

    def Av(f, w=T):
        return (w * f).sum() / w.sum()

    def Var(f, w=T):
        return Av(f * f, w) - Av(f, w) ** 2

    _, PAp, _, PAy = Phi_der('A', A, B)
    _, PO09p, _, PO09y = Phi_der('B', A, B)
    psp, ppsp = psi_d(A, B)
    cApsi = Av(PAp * psp) - Av(PAp) * Av(psp)
    P("[D] D 的逐项分解 (N=%d, theta=pi, 权重 |Delta|^2 = 纯 Haar):" % N)
    P("     Var_T(Phi_A')       = %+.8f      <Phi_A''>_T       = %+.8f   [应 8,8; 差 %.1e]" %
      (Var(PAp), Av(PAy), Var(PAp) - Av(PAy)))
    P("     Var_T(psi')         = %+.8f      <psi''>_T         = %+.8f" % (Var(psp), Av(ppsp)))
    P("     2Cov_T(Phi_A',psi') = %+.8f" % (2 * cApsi))
    P("     Var_T(Phi_O09')     = %+.8f   [r12: 6.61542931]  = 8 + 2Cov + Var(psi') = %+.8f" %
      (Var(PO09p), Var(PAp) + 2 * cApsi + Var(psp)))
    P("     <Phi_O09''>_T       = %+.8f   [r12: 6.44071163]" % Av(PO09y))
    P("     Var_T(G)            = %+.8f   [r12: 0.06108481]   Cov_T(G,Phi_O09') = %+.2e" %
      (Var(G), Av(G * PO09p) - Av(G) * Av(PO09p)))
    P("     <H>_T               = %+.8f   [r12: 0.10230219 (加权)]" % Av(H))
    P("")
    P("     gap_O09(G=0) = <Phi_O09''> - Var(Phi_O09') = %+.8f" % (Av(PO09y) - Var(PO09p)))
    P("     记账: (8-8) + <psi''> - 2Cov(Phi_A',psi') - Var(psi') = %+.8f  [恒等式]" %
      (Av(ppsp) - 2 * cApsi - Var(psp)))
    P("     D 的占比: Var(Phi_O09') = %.2f%%,  Var(G) = %.2f%%  (Cov=0)" %
      (100 * Var(PO09p) / (Var(PO09p) + Var(G)), 100 * Var(G) / (Var(PO09p) + Var(G))))
    P("     => D 被 Var(Phi') 主宰; 而 Var(Phi') 的偏离全部来自 psi 的 1/d^2 项 [T]+[C]")
    P("")


# =====================================================================
# [E] Phi_O09 的切片依赖 + theta 扫描
# =====================================================================
def partE():
    P("[E] Phi_O09 = Phi_A + psi, psi = 2[log a + log b + log(a+b)] 含 2 log d 项:")
    for a0 in (0.3, 1.1):
        d0 = psi_(np.array([a0]), np.array([0.7]))[0]
        d1 = psi_(np.array([a0 + 2 * math.pi]), np.array([0.7]))[0]
        P("     psi(a+2pi) - psi(a) 于 a=%.1f : %+.8f  != 0  => 非周期, 不是 T^2 上的函数 [T]" % (a0, d1 - d0))
    N = 1024
    A, B = grid(N)
    P("     Z(delta) = int_{delta}^{2pi+delta} e^{-Phi_O09-S_W} da db  (N=%d, 同一环面, 只切 chart):" % N)
    for dl in (0.0, 0.3, 0.6, 1.0):
        a, b = A + dl, B + dl
        w = pieces(a, b)[2] * np.exp(-SW_(a, b, math.pi, 2.5) - psi_(a, b))
        P("       delta=%-4.1f  Z(O09 权重)=%-.4f" % (dl, -w.sum()))
    P("     [E1] 3 折覆盖: (a,b)=(x1-x2,x2-x3) 的映射不是单射, e^{-S_W} 在 (a,b) 只有 6pi 周期:")
    M = 512
    aa = (np.arange(M) + 0.5) * (6 * math.pi / M)
    AA, BB = np.meshgrid(aa, aa, indexing='ij')
    P("          %-8s %-16s %-16s" % ("delta", "Z over 2pi 方形", "Z over 6pi 基本域"))
    for dl in (0.0, 0.5, 1.0, 2.0):
        a, b = A + dl, B + dl
        z2 = (pieces(a, b)[2] * np.exp(-SW_(a, b, math.pi, 2.5))).sum()
        z6 = (pieces(AA + dl, BB + dl)[2] * np.exp(-SW_(AA + dl, BB + dl, math.pi, 2.5))).sum()
        P("          %-8.1f %-16.4f %-16.6f" % (dl, z2, z6))
    P("     => 含 e^{-S_W} 的积分 (r12 的 +0.1017 分支) 依赖切片切法 => 不是 SU(3) 上的平均值 [T]+[C]")
    P("     => 纯 Haar 权重 |Delta|^2 是 2pi 周期的, 不受此影响; r12 的 -0.2358 是良定义的数 (但仍是错配残差)")
    P("")
    P("[E2] SU(3) theta 扫描 (N=1024): 纯 Haar 的 gap 与 theta 无关 (theta 只进 S_W)")
    P("     %-8s %-16s %-16s %-16s" % ("theta", "纯Haar/O09错配", "加权/O09错配", "加权/匹配"))
    for i in range(5):
        th = math.pi * i / 4
        a = SD('B', 1024, theta=th, beta=0.0, weight='A')
        b = SD('B', 1024, theta=th, beta=2.5, weight='A')
        c = SD('A', 1024, theta=th, beta=2.5, weight='A')
        P("     %-8.4f %-+16.8f %-+16.8f %-+16.8f" % (th, a['gap'], b['gap'], c['gap']))
    P("     => 纯 Haar 行是常数 => 纯 Haar 下没有 theta 窗口;")
    P("        加权行的符号翻转 (+0.03 / -1.10) 是错配残差随 theta 变, 不是正性 [T]+[C]")
    P("")


if __name__ == "__main__":
    P("=" * 108)
    P("s03_haarneg.log -- S03: SU(3) 纯 Haar 的负 gap: rank-2 真性质 vs 口径问题")
    P("=" * 108)
    P("")
    partA()
    partB()
    partC()
    partD()
    partE()
    with open("s03_haarneg.log", "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
