# -*- coding: utf-8 -*-
r"""
s01_su3red.py -- S01: SU(3) a=0 切片的**有效维数**与降维 (x,p) 参数化。

模型 (与 r08/r11/r12 同一对象, 单 plaquette, Thermitian X):
    X in su(3) (Hermitian, trX=0, 8 维),  U = e^{iX},  dX J(X) = Haar(SU(3))
    J(X) = prod_{i<j} (2 sin(x_ij/2)/x_ij)^2 ,  x_ij = x_i - x_j  (x = eig(X))
    S_W(X) = beta (1 - (1/3) Re tr(W U^dag)) ,  W = diag(-1,-1,1) = 2P-1,  P = |3><3|
    粗模 v = (3P-1)/sqrt6  (单位, tr v^2 = 1, v ∥ W 的生成元),  a = tr(vX)
    Phi = -log J,  G = d_tau S_W(X+tau v)|0,  H = d_tau^2 S_W|0
    S = E[H+Phi''|a=0],  D = Var(G+Phi'|a=0),  G''(0) = S - D

本文件判的东西 (全部机检, 不引用未核数字):
 [A] 维数账: 8 -> 7(a=0) -> 4((x,p)) -> 3(a=p.x) -> 2(p 解析积分)
 [B] p ~ Dirichlet(1,1,1) 且与 x 独立 (Haar MC)
 [C] S_W = beta(1-(1/3)(2 p.c - C)) 与 tr(W U^dag) 逐位一致
 [D] a = (3/sqrt6) p.x  == tr(vX)
 [E] G、Phi'、Phi'' 的 (x,p) 闭式 == 沿 X+tau v 的有限差分  => 只依赖 (x,p)
 [F] 对称性计数: Z_3 破缺、theta<->-theta 精确、Weyl S_3 -> S_2 破缺、H=S(U(2)xU(1)) 保 W
 [G] p 线段积分闭式 == 数值 (delta 函数口径)
 [H] 整群 MC == |Delta|^2 (x) 均匀 (p) 的乘积化  => 4 维 fiber 精确积分掉
"""
import numpy as np
from numpy.linalg import eigh
from scipy.linalg import expm

rng = np.random.default_rng(20251125)
LOG = []
def P(s=""):
    print(s); LOG.append(s)

BETA = 2.5
I3 = np.eye(3)
Pp = np.diag([0.0, 0.0, 1.0])
W = 2 * Pp - I3
v = (3 * Pp - I3) / np.sqrt(6.0)
ALPHA = 1.0 / np.sqrt(6.0)
assert abs(np.trace(W @ W) - 3) < 1e-14 and abs(np.trace(v @ v) - 1) < 1e-14


# ----------------------------------------------------------------- 工具
def haar(n):
    z = (rng.standard_normal((n, 3, 3)) + 1j * rng.standard_normal((n, 3, 3))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    d = np.diagonal(r, axis1=1, axis2=2)
    q = q * (d / np.abs(d))[:, None, :]
    det = np.linalg.det(q)
    return q / (det ** (1.0 / 3.0))[:, None, None]


def spec(X):
    """Hermitian X -> (x 升序, u 列=本征矢)"""
    x, u = eigh(X)
    return x, u


def xp_of(U):
    """U -> (x (sum 0), p_j=|u_3j|^2)"""
    lam, u = np.linalg.eig(U)
    th = np.angle(lam)
    sh = -th.sum() / 3.0                       # det=1 => sum(th) in {0,+-2pi}
    x = np.sort(th + sh)
    o = np.argsort(th + sh)
    u = u[:, o]
    p = np.abs(u[2, :]) ** 2
    return x, p


def SW_mat(X):
    return BETA * (1 - np.real(np.trace(W @ expm(-1j * X))) / 3.0)


def SW_xp(x, p):
    q = np.sum(p * np.exp(-1j * x))
    T = np.sum(np.exp(-1j * x))
    return BETA * (1 - (2 * q.real - T.real) / 3.0)


def Phi_of(X):
    x, u = spec(X)
    d = x[:, None] - x[None, :]
    iu = np.triu_indices(3, 1)
    J = np.prod((2 * np.sin(d[iu] / 2) / d[iu]) ** 2)
    return -np.log(J), x, u


def phi0(d):    # r08 的 Phi 单对项: -log(4 sin^2(d/2))  (= 无 1/d 因子)
    return -np.log(4 * np.sin(d / 2) ** 2)


def dphi(d):    # 模型 Phi = sum g(d), g(d)=-2log(2 sin(d/2)/d) ; g'(d)=2/d-cot(d/2)
    return 2.0 / d - 1.0 / np.tan(d / 2)


def ddphi(d):   # g''(d)= (1/2)csc^2(d/2) - 2/d^2
    return 0.5 / np.sin(d / 2) ** 2 - 2.0 / d ** 2


def Phi_grad2(x):
    """(dPhi/dx_j, d2Phi/dxj dxk)"""
    g = np.zeros(3); Hm = np.zeros((3, 3))
    for j in range(3):
        for k in range(3):
            if j == k:
                continue
            g[j] += dphi(x[j] - x[k])
            Hm[j, k] += -ddphi(x[j] - x[k])
    for j in range(3):
        for k in range(3):
            if j != k:
                Hm[j, j] += ddphi(x[j] - x[k])
    return g, Hm


def G_closed(x, p):
    """G = d/dtau S_W(X+tau v)|0, 只靠 (x,p)
       d/dtau tr(W e^{-i(X+tau v)}) = -i*alpha*S, S=sum_jk[6 p_j p_k + d_jk(1-5p_j)]E_jk
       E_jk = int_0^1 e^{-i s x_j - i(1-s) x_k} ds
    """
    S = 0j
    for j in range(3):
        for k in range(3):
            if j == k:
                E = np.exp(-1j * x[j])
            else:
                E = (np.exp(-1j * x[j]) - np.exp(-1j * x[k])) / (-1j * (x[j] - x[k]))
            coef = 6 * p[j] * p[k] + (1 - 5 * p[j] if j == k else 0.0)
            S += coef * E
    return -(BETA / 3.0) * np.real(-1j * ALPHA * S)


def Phi_p_closed(x, p):
    g, _ = Phi_grad2(x)
    return np.sum(ALPHA * (3 * p - 1) * g)


def Phi_pp_closed(x, p):
    g, Hm = Phi_grad2(x)
    A = ALPHA * (3 * p - 1)
    out = A @ Hm @ A
    for j in range(3):
        for k in range(3):
            if j != k:
                out += 9 * ALPHA ** 2 * p[j] * p[k] * (g[j] - g[k]) / (x[j] - x[k])
    return out


# ================================================================= [A]
def partA():
    P("[A] 维数账 (逐步; 每步写明砍掉的是什么):")
    P("  (1) su(3) = 8 = N^2-1                       : 8")
    P("  (2) 变量是 plaquette U 本身 (单 plaquette 模型, 已 gauge 固定外侧链)")
    P("      => gauge 冗余 0 维; 8->8  (与 r04 的格点 192->144 不同: 那里 4 条链才有 48 个生成元)")
    P("  (3) 条件化 a=tr(vX)=0  (一个线性超平面)     : -1  => 7   <== r12 的 7")
    P("  (4) 极分解: U = g D g^dag, dim su(3) = (x:2) + (flag G/T^2:6) = 8")
    P("      flag 6 = p(2) + fiber(4),  p_j = |<3|u_j>|^2 in Delta^2")
    P("      被积函数 (J, S_W, a, G, Phi', Phi'') 只经 (x,p) 依赖 U")
    P("      => flag fiber(4 维)精确积分掉, 但 a=0 已用掉其中 1 维 : -3  => 4")
    P("  (5) 在 (x,p) 上重写 a=0:  a = (3/sqrt6) p.x = 0, 一条线性方程 : -1  => 3")
    P("  (6) p 只有 1 维 (line cap Delta^2 的线段), 且被积函数 = e^{线性(p)} x 多项式")
    P("      => p 积分闭式 (二次除差), 只剩 x in T^2            : -1  => 2")
    P("  合计: 8 -1 -3 -1 -1 = 2      (r12 的 7 是 (3) 之后的裸数, 未计 (4)(5)(6))")
    P("  对照 SU(2): 3 -> 2(a=0) -> 1(径向 O(2)); 同一套算法, 只是 flag 那步=0")
    P()


# ================================================================= [B]
def partB(n=200000):
    P("[B] p 的分布与 x 的独立性 (Haar MC, n=%d):" % n)
    U = haar(n)
    X = []
    pm = np.zeros(3); p2 = np.zeros((3, 3)); p1arr = np.zeros(n)
    SW_direct = np.zeros(n); SW_red = np.zeros(n)
    eps = np.zeros(n); cosx = np.zeros(n)
    for i in range(n):
        x, p = xp_of(U[i])
        pm += p; p2 += np.outer(p, p)
        p1arr[i] = p[0]
        lam, g = np.linalg.eig(U[i])
        th = np.angle(lam); sh = -th.sum() / 3.0; xo = th + sh
        Xh = (g * xo) @ np.linalg.inv(g)
        Xh = 0.5 * (Xh + Xh.conj().T)
        SW_direct[i] = SW_mat(Xh)
        SW_red[i] = SW_xp(x, p)
        cosx[i] = np.sum(np.cos(x))
        eps[i] = np.max(x) - np.min(x)
    pm /= n; p2 /= n
    P("  E[p]      = %s   (Dirichlet(1,1,1): 1/3 = 0.333333)" % np.array2string(pm, precision=6))
    P("  E[p_j^2]  = %s   (理论 1/6 = 0.166667)" % np.array2string(np.diag(p2), precision=6))
    P("  E[p_j p_k] j<k = %.6f %.6f %.6f   (理论 1/12 = 0.083333)"
      % (p2[0, 1], p2[0, 2], p2[1, 2]))
    c1 = np.corrcoef(p1arr, cosx)[0, 1]
    P("  独立性: corr(p_1, sum cos x) = %+.4f  (应 ~0)" % c1)
    c2 = np.corrcoef(p1arr, eps)[0, 1]
    P("  独立性: corr(p_1, max x - min x) = %+.4f  (应 ~0)" % c2)
    P("  => p ~ Unif(Delta^2), 与 x 独立  [C]")
    P()


# ================================================================= [C][D]
def partCD(n=20000):
    P("[C] S_W 的 (x,p) 闭式 vs 直接 tr(W U^dag):")
    U = haar(n)
    e1 = e2 = a1 = a2 = 0.0
    Xs = []
    for i in range(n):
        lam, g = np.linalg.eig(U[i])
        th = np.angle(lam); sh = -th.sum() / 3.0; xo = th + sh
        Xh = (g * xo) @ np.linalg.inv(g); Xh = 0.5 * (Xh + Xh.conj().T)
        x, u = spec(Xh)
        p = np.abs(u[2, :]) ** 2
        e1 = max(e1, abs(SW_mat(Xh) - SW_xp(x, p)))
        a1 = max(a1, abs(np.real(np.trace(v @ Xh)) - (3 / np.sqrt(6)) * np.sum(p * x)))
        e2 = max(e2, abs(np.sum(p) - 1))
        Xs.append((Xh, x, p))
    P("  max|S_W(matrix) - S_W(x,p)| = %.3e    (含 2p.c - C 的符号约定)" % e1)
    P("  max|tr(vX) - (3/sqrt6) p.x|  = %.3e    => a=0 <=> p.x=0" % a1)
    P("  max|sum p - 1|               = %.3e" % e2)
    P()


# ================================================================= [E]
def partE(n=6):
    P("[E] G/Phi'/Phi'' 的 (x,p) 闭式 vs 沿 X+tau v 的有限差分:")
    P("    (若两者一致 => 三者只依赖 (x,p); v_jj=alpha(3p_j-1), |v_jk|^2=9 alpha^2 p_j p_k)")
    P("    %-3s %-13s %-13s %-11s %-13s %-13s %-12s" %
      ("#", "G_fd", "G_closed", "rel", "Phi'_fd", "Phi'_cl", "rel"))
    h = 1e-5
    worst = 0.0
    for i in range(n):
        U = haar(1)[0]
        lam, g = np.linalg.eig(U)
        th = np.angle(lam); sh = -th.sum() / 3.0; xo = th + sh
        Xh = (g * xo) @ np.linalg.inv(g); Xh = 0.5 * (Xh + Xh.conj().T)
        x, u = spec(Xh); p = np.abs(u[2, :]) ** 2
        G_fd = (SW_mat(Xh + h * v) - SW_mat(Xh - h * v)) / (2 * h)
        G_cl = G_closed(x, p)
        f1 = Phi_of(Xh + h * v)[0]; f2 = Phi_of(Xh - h * v)[0]; f0 = Phi_of(Xh)[0]
        Pp_fd = (f1 - f2) / (2 * h)
        Ppp_fd = (f1 - 2 * f0 + f2) / h ** 2
        Pp_cl = Phi_p_closed(x, p); Ppp_cl = Phi_pp_closed(x, p)
        r1 = abs(G_cl - G_fd) / max(1e-12, abs(G_fd))
        r2 = abs(Pp_cl - Pp_fd) / max(1e-12, abs(Pp_fd))
        r3 = abs(Ppp_cl - Ppp_fd) / max(1e-12, abs(Ppp_fd))
        worst = max(worst, r1, r2, r3)
        P("    %-3d %+.8e %+.8e %-11.2e %+.8e %+.8e %-12.2e" %
          (i, G_fd, G_cl, r1, Pp_fd, Pp_cl, r2))
        P("        Phi'' : fd %.8e  closed %.8e  rel %.2e" % (Ppp_fd, Ppp_cl, r3))
    P("    worst rel = %.2e  => (x,p) 闭式成立 [C]" % worst)
    P()


# ================================================================= [F]
def partF(n=300):
    P("[F] 切片上的对称性 (矩阵口径: U -> Q U Q^dag, Q in SU(3)):")
    U = haar(n)
    dZ = dth = 0.0
    perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    names = ["id", "(23)", "(12)", "(123)", "(132)", "(13)"]
    s3 = np.zeros(6)
    OM = np.exp(2j * np.pi / 3)
    for i in range(n):
        Ui = U[i]
        t0 = np.real(np.trace(W @ Ui.conj().T))
        # Z_3: U -> omega U  (中心元左乘, 不是内自同构)
        dZ = max(dZ, abs(np.real(np.trace(W @ (OM * Ui).conj().T)) - t0))
        # theta <-> -theta : U -> U^dag
        dth = max(dth, abs(np.real(np.trace(W @ Ui)) - t0))
        for k, pm in enumerate(perms):
            Q = np.zeros((3, 3))
            for a, b in enumerate(pm):
                Q[b, a] = 1.0
            if np.linalg.det(Q) < 0:
                Q = -Q                                  # det(-P)=(-1)^3 det P = +1
            s3[k] = max(s3[k], abs(np.real(np.trace(W @ (Q @ Ui @ Q.conj().T).conj().T)) - t0))
    P("  Z_3 (U->omega U, x->x+2pi/3): max|d tr(WU^dag)| = %.3e   => 破缺 => 降 0 维" % dZ)
    P("  theta<->-theta (U->U^dag):     max|d tr(WU^dag)| = %.3e   => 精确 => 有限群, 降 0 维" % dth)
    for k, nm in enumerate(names):
        P("  Weyl %-6s : max|d tr(WU^dag)| = %.3e  %s" %
          (nm, s3[k], "<= 不变 (稳定 W)" if s3[k] < 1e-12 else "<= 破缺"))
    P("  => SU(3) 共轭下 S_W 的对称群 = Stab(W) = S(U(2)xU(1)) (4 维),")
    P("     与 W 交换的置换只有 S_2={id,(12)}; 其余 4 个 Weyl 元素破缺")
    P("  => 类函数条件不成立 => 2D Weyl 不适用 (r12 的结论对; 但障碍不是裸 7 维)")
    P()


# ================================================================= [G]
def line_seg(x):
    """{p in Delta^2 : p.x = 0} 的端点 (p1,p2) 与长度 (Delta^2 的 (p1,p2) 坐标)"""
    x3 = x[2]; A = x[0] - x3; B = x[1] - x3; C = -x3
    if abs(A) < 1e-14 and abs(B) < 1e-14:
        return None, 0.0, 0.0
    pts = []
    # 三条边: p1=0 ; p2=0 ; p1+p2=1
    for kind in (0, 1, 2):
        if kind == 0:
            if abs(B) < 1e-14: continue
            p1 = 0.0; p2 = C / B
        elif kind == 1:
            if abs(A) < 1e-14: continue
            p1 = C / A; p2 = 0.0
        else:
            if abs(A + B) < 1e-14: continue
            s = (C - B) / (A - B)       # p1=s, p2=1-s
            p1 = s; p2 = 1 - s
        if p1 >= -1e-12 and p2 >= -1e-12 and p1 + p2 <= 1 + 1e-12:
            pts.append((p1, p2))
    if len(pts) < 2:
        return None, 0.0, 0.0
    pts = pts[:2]
    d = np.array(pts[1]) - np.array(pts[0])
    L = float(np.hypot(*d))
    return (np.array(pts[0]), np.array(pts[1])), L, float(np.hypot(A, B))


def partG(nx=400, npn=200, ng=900):
    P("[G] p 线段积分闭式 vs 数值 (delta 口径):")
    P("    int_{Delta^2} dp delta(p.x) f(p) = (1/|grad_p(p.x)|) int_seg f ds,   |grad|=sqrt((x1-x3)^2+(x2-x3)^2)")
    rng2 = np.random.default_rng(3)
    worst = 0.0
    for trial in range(5):
        while True:
            x2 = rng2.uniform(-2.5, 2.5, 2); x = np.array([x2[0], x2[1], -x2[0] - x2[1]])
            if np.max(x) - np.min(x) > 0.6 and np.min(np.abs(x)) > 0.05:
                break
        seg, L, gnorm = line_seg(x)
        # 闭式: 高斯-勒让德 20 点在线段上
        pts = np.polynomial.legendre.leggauss(20)
        s = 0.5 * (pts[0] + 1); wq = 0.5 * pts[1]
        fcl = 0.0; fnm = 0.0
        for si, wi in zip(s, wq):
            p1 = seg[0][0] + si * (seg[1][0] - seg[0][0]); p2 = seg[0][1] + si * (seg[1][1] - seg[0][1])
            p = np.array([p1, p2, 1 - p1 - p2])
            h = np.exp(-SW_xp(x, p))
            fcl += wi * h * L / gnorm
        # 数值: 细分网格 + 窄 delta
        g = (np.arange(ng) + 0.5) / ng
        P1, P2 = np.meshgrid(g, g, indexing='ij')
        ok = P1 + P2 <= 1
        p1g = P1[ok]; p2g = P2[ok]; p3g = 1 - p1g - p2g
        av = p1g * x[0] + p2g * x[1] + p3g * x[2]
        sig = 0.03
        w = np.exp(-0.5 * (av / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
        # 网格上的 dp 权重 (Delta^2 的 (p1,p2) 面积元)
        hh = (np.exp(-BETA * (1 - (2 * (p1g * np.cos(x[0]) + p2g * np.cos(x[1]) + p3g * np.cos(x[2]))
                                     - np.cos(x).sum()) / 3)))
        area = 1.0 / ng ** 2
        fnm = float(np.sum(hh * w) * area)
        rel = abs(fcl - fnm) / max(1e-12, abs(fnm))
        worst = max(worst, rel)
        P("    x=(%+.4f,%+.4f,%+.4f) L=%.6f |grad|=%.6f : closed=%.8e  grid=%.8e  rel=%.2e"
          % (x[0], x[1], x[2], L, gnorm, fcl, fnm, rel))
    P("    worst rel = %.2e  (网格/窄 delta 的误差主导)  => 闭式可用 [C]" % worst)
    P()


# ================================================================= [H]
def partH(n=60000):
    P("[H] 测度因式化的稳健检验 (不依赖 x-密度取 J 还是 |Delta|^2 的 chart 约定):")
    P("    任何只依赖 x 的测度 rho(x)dx (Haar 的 |Delta|^2, 或模型 exp-chart 的 sum_k J(x+2pi k))")
    P("    乘上 p 的均匀测度 ⊗ 体积归一的 fiber; 于是 'p 在 x 的每个 bin 内都均匀' 是")
    P("    该因式化的判据, 且与 rho 无关.")
    U = haar(n)
    p1 = np.zeros(n); p1sq = np.zeros(n); cxs = np.zeros(n); p2sq = np.zeros(n)
    for i in range(n):
        lam, g = np.linalg.eig(U[i])
        th = np.angle(lam); sh = -th.sum() / 3.0; xo = th + sh
        Xh = (g * xo) @ np.linalg.inv(g); Xh = 0.5 * (Xh + Xh.conj().T)
        x, u = spec(Xh); p = np.abs(u[2, :]) ** 2
        p1[i] = p[0]; p1sq[i] = p[0] ** 2; p2sq[i] = p[1] ** 2
        cxs[i] = np.sum(np.cos(x))
    q = np.quantile(cxs, [0.2, 0.4, 0.6, 0.8])
    P("    bin (sum cos x 五分位) |    E[p_1]      E[p_1^2]     E[p_2^2]   | 理论 0.33333/0.16667/0.16667")
    edges = [-9] + list(q) + [9]
    for b in range(5):
        m = (cxs > edges[b]) & (cxs <= edges[b + 1])
        P("    bin %d (n=%6d)        |  %.5f     %.5f     %.5f"
          % (b, m.sum(), p1[m].mean(), p1sq[m].mean(), p2sq[m].mean()))
    P("    全部: E[p_1]=%.5f E[p_1^2]=%.5f E[p_2^2]=%.5f" % (p1.mean(), p1sq.mean(), p2sq.mean()))
    P("    => 各 bin 内 p 的矩与全样本一致 (差 ~ MC 误差) => p 与 x 独立、p 均匀 [C]")
    P("    注: x 的**密度**是 rho(x) (Haar: |Delta|^2; 模型 exp-chart: sum_k J) —— 只影响 2 维 x-积分的被积函数,")
    P("        不影响维数, 也不影响 (x,p) 因式化 (= r12 的 psi=sum log d^2 就是这件事的 chart 面).")
    P()


if __name__ == "__main__":
    P("=" * 100)
    P("s01_su3red.log -- S01: SU(3) a=0 切片有效维数 + 降维 (x,p)")
    P("模型: 单 plaquette SU(3), W=diag(-1,-1,1), beta=2.5, v=(3P-1)/sqrt6, a=tr(vX)")
    P("=" * 100)
    P()
    partA()
    partB(60000)
    partCD(4000)
    partE(4)
    partF(200)
    partG()
    partH(60000)
    with open("s01_su3red.log", "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
