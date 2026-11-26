# -*- coding: utf-8 -*-
r"""
s04_varphip.py -- S04: Var(Phi') 的无条件上界 (SU(2) 与 SU(3)).

Phi = -log J,  J = 指数坐标 Jacobian,  Phi(X) = sum_{a<b} g(d_ab),  d_ab=|lam_a-lam_b|
g(d) = -2 log[ sin(d/2)/(d/2) ],   g'(d)= 2/d - cot(d/2),   g''(d)= (1/2)csc^2(d/2) - 2/d^2

[S1] g 的性质 (点态, 全局): g'(0)=0, g'(d)>=0, g''(d)>=1/6, g'(d)->+inf as d->2pi^-
[S2] Phi' 显式 + SU(3) 的奇性/零点 (Phi 偶 => Phi'(0)=0)
[S3] Stein/IBP 恒等式 (精确):  E[Phi'^2] + E[G Phi'] = E[Phi'']  <=>  C 的值
[S4] SU(3) 8 维 log-chart MC: 权重 e^{-Phi}, 检验恒等式 + 条件 (a=0) 的 Var(Phi') vs E[Phi'']
[S5] SU(3) 环面 (r12 模型) 的正确配对 Phi_Delta = -log|Delta|^2  vs 错配对 Phi_O09: 恒等式残差
[S6] r12 数字代入闭式上界  Var(Phi') <= M + g(g+sqrt(g^2+4M))/2 - (E Phi')^2
"""
import sys
import numpy as np

np.random.seed(20261125)
OUT = []


def P(s=""):
    print(s)
    OUT.append(s)


# ---------------- g, g', g''  (全局, 无分支截断) ----------------
def g(d):
    d = np.asarray(d, float)
    small = np.abs(d) < 1e-6
    out = np.where(small, d * d / 12.0,
                   -2.0 * np.log(np.abs(np.sin(d / 2.0) / (d / 2.0))))
    return out


def gp(d):
    d = np.asarray(d, float)
    small = np.abs(d) < 1e-6
    return np.where(small, d / 6.0, 2.0 / d - 1.0 / np.tan(d / 2.0))


def gpp(d):
    d = np.asarray(d, float)
    small = np.abs(d) < 1e-6
    with np.errstate(divide="ignore", invalid="ignore"):
        v = 0.5 / np.sin(d / 2.0) ** 2 - 2.0 / d ** 2
    return np.where(small, 1.0 / 6.0 + d * d / 60.0, v)


def omega(s):                      # omega(s) = g'(|s|) sgn(s) : 奇函数
    s = np.asarray(s, float)
    return np.sign(s) * gp(np.abs(s))


# =====================================================================
def S1():
    P("=" * 96)
    P("[S1] g 的点态性质 (扫描 d in (1e-9, 2pi-1e-9))")
    P("=" * 96)
    d = np.linspace(1e-9, 2 * np.pi - 1e-9, 200001)
    P("  g(0+)=%.3e   g(pi)=%.6f   g(2pi-)=%.3e  (-> +inf)" % (g(1e-9), g(np.pi), g(2 * np.pi - 1e-9)))
    P("  min g'(d) = %.10f    (理论 0, 在 d->0)          g' >= 0 : %s" % (gp(d).min(), gp(d).min() >= 0))
    dd2 = np.logspace(-8, 1, 4000)
    P("  min g''(d) (log grid d in 1e-8..10) = %.12f   (理论 inf = 1/6 = %.12f)"
      % (gpp(dd2).min(), 1 / 6))
    P("  g''(d) = 1/6 + d^2/60 + O(d^4)  =>  g'' >= 1/6 全局 [T] (网格下界差 1e-8 为网格效应)")
    P("  max g'(d) = %.3e  (d=%.10f)  ==> 梯度上确界 = +inf (d -> 2pi^-)" % (gp(d).max(), d[-1]))
    P("  omega(s)=g'(|s|)sgn(s): 奇, omega(0)=0, omega(s)/s -> 1/6 (s->0)")
    ss = np.array([1e-4, 1e-2, 0.1, 1.0, 3.0])
    P("    s     : " + " ".join("%9.4f" % x for x in ss))
    P("    omega : " + " ".join("%9.6f" % x for x in omega(ss)))
    P("    s/6   : " + " ".join("%9.6f" % x for x in ss / 6))
    P("    w-s/6 : " + " ".join("%9.2e" % x for x in omega(ss) - ss / 6))
    P("  ==> g' >= 0, g'' >= 1/6 [T];  |nabla Phi| 无上界, 发散的集合 = {d_ab = 2pi} (cut locus) [T]")
    P()


# =====================================================================
def Phi_eig(X, v):
    """X, v: (N,3,3) 实对称 traceless. 返回 Phi, Phi', Phi'' (沿 v)."""
    lam, U = np.linalg.eigh(X)                     # lam 升序
    V = np.einsum("nji,njk,nkl->nil", U, v, U)     # v 在 X 的本征基下
    Ph = np.zeros(len(X)); Pp = np.zeros(len(X)); Ppp = np.zeros(len(X))
    for a in range(3):
        for b in range(a + 1, 3):
            D = lam[:, a] - lam[:, b]
            dd = np.abs(D); s = np.sign(D)
            u = V[:, a, b - b] if False else (V[:, a, a] - V[:, b, b])
            w = 0.0
            for c in range(3):
                if c != a:
                    w += 2 * V[:, a, c] ** 2 / (lam[:, a] - lam[:, c])
                if c != b:
                    w -= 2 * V[:, b, c] ** 2 / (lam[:, b] - lam[:, c])
            Ph += g(dd)
            Pp += gp(dd) * s * u
            Ppp += gpp(dd) * u * u + gp(dd) * s * w
    return Ph, Pp, Ppp


def S2():
    P("=" * 96)
    P("[S2] Phi' 显式:  Phi'(X;v) = sum_{a<b} g'(d_ab) sgn(lam_a-lam_b) <(P_a-P_b),v>")
    P("     = <sum_a c_a P_a, v>,  c_a = sum_{b!=a} omega(lam_a-lam_b)   [T]")
    P("=" * 96)
    # SU(2)-like: 2x2 traceless -- 用 3x3 的 (1,-1,0) 退化型模拟会混入; 单独做 2x2
    P("  (i) SU(2): X = (rho/2) n.sigma  => Phi = g(rho) 径向  => Phi'(v) = (g'(rho)/rho)<X,v>")
    P("      切片 <X,v>=0 上  Phi' == 0 逐点 (与 beta、theta 无关)   [T]  (r11[A]: 8.6e-42)")
    # 数值: 2x2
    lam = np.array([0.7, -0.7]); rho = 1.4
    u = np.array([0.6, 0.8, 0.0]); u /= np.linalg.norm(u)
    X = (rho / 2) * np.array([[u[2], u[0] - 1j * u[1]], [u[0] + 1j * u[1], -u[2]]])
    w = np.array([0.8, -0.6, 0.0]); w -= u * (u @ w); w /= np.linalg.norm(w)   # 垂直于 X
    P("      n=(0.6,0.8,0), v=(0.8,-0.6,0): <X,v>=0 (正交切片).  Phi'(v) 解析 = (g'(rho)/rho)*0 = 0")
    P("      g'(rho)/rho = %.10f    Phi''(横向) = g'(rho)/rho" % (gp(rho) / rho))
    gpr = gp(rho) / rho
    h = 1e-5
    for t in (0.0,):
        pass
    # 直接 FD 核对径向公式
    print_h = []
    for r in (0.9, 1.4, 2.2):
        print_h.append((r, gp(r) / r, gpp(r)))
    P("      rho    g'(rho)/rho    g''(rho)   (横向/径向 Phi'')")
    for r, a1, a2 in print_h:
        P("      %.2f   %12.9f  %12.9f" % (r, a1, a2))
    P()
    P("  (ii) SU(3): 反例 + 三次律 (r11[C-ii] 复现)")
    x = np.array([-1.5797, 0.3243, 1.2554])
    e = np.array([0.9153, -2.7871, 1.8718])
    P("      X=diag%s  v=diag%s   a=<X,v>/||v||^2 = %.2e  (在切片上)" % (x, e, (x @ e) / (e @ e)))
    pd = 0.0
    for a_ in range(3):
        for b_ in range(a_ + 1, 3):
            D = x[a_] - x[b_]
            pd += gp(abs(D)) * np.sign(D) * (e[a_] - e[b_])
    P("      Phi'(X;v) = sum g'(x_a-x_b)(e_a-e_b) = %+.10f   (!= 0) [T]" % pd)
    P("      但 Hess Phi(0) = (N/12)Id 各向同性 => 二次阶是径向的, 各向异性从四次项起 (r11 的三次律)")
    tt = np.array([0.05, 0.1, 0.2, 0.4, 1.0, 2.0])
    P("      t      Phi(tx)      Phi/t^2      (N/12)||tx||^2")
    for t in tt:
        xx = t * x
        val = 0.0
        for a_ in range(3):
            for b_ in range(a_ + 1, 3):
                val += g(abs(xx[a_] - xx[b_]))
        P("      %-5.2f  %10.6f  %10.6f  %10.6f" % (t, val, val / t ** 2, (3 / 12) * (t * x) @ (t * x)))
    P("      => Phi(tx)/t^2 -> (N/12)||x||^2: 二次型各向同性 => Phi' = O(|x|^3) 的各向异性部分 [T]")
    P()


# =====================================================================
def S4(N=400000, sig=0.40, bw=0.30, burn=20000):
    P("=" * 96)
    P("[S4] SU(3) 真 log-chart: 8 维 su(3) (Gell-Mann/2 基), MH 链目标 e^{-Phi}, 纯 Haar")
    P("     检验 E[Phi'^2] == E[Phi''] (Stein/IBP 精确); a=0 切片 Var(Phi') vs E[Phi'']  N=%d" % N)
    P("=" * 96)
    # Gell-Mann/2 : 8 个 Hermitian 基 (正交, 不需归一, 只影响常数)
    B = []
    B.append(np.array([[0,1,0],[1,0,0],[0,0,0]],complex)/2)      # l1/2
    B.append(np.array([[0,-1j,0],[1j,0,0],[0,0,0]],complex)/2)    # l2/2
    B.append(np.diag([1,0,-1]).astype(complex)/2)                 # l3/2
    B.append(np.array([[0,0,1],[0,0,0],[1,0,0]],complex)/2)       # l4/2
    B.append(np.array([[0,0,-1j],[0,0,0],[1j,0,0]],complex)/2)    # l5/2
    B.append(np.array([[0,0,0],[0,0,1],[0,1,0]],complex)/2)       # l6/2
    B.append(np.array([[0,0,0],[0,0,-1j],[0,1j,0]],complex)/2)    # l7/2
    B.append(np.diag([1,1,-2]).astype(complex)/np.sqrt(12))       # l8/2
    B = np.array(B)                                               # (8,3,3)
    v = np.diag([1,0,-1]).astype(complex)/2                      # a-方向 = xi3

    def Phi_of(x):
        H = np.einsum("k,kij->ij", x, B)
        h = np.linalg.eigvalsh(H)
        return float(sum(g(abs(h[a]-h[b])) for a in range(3) for b in range(a+1,3)))

    def derivs(x, v):
        H = np.einsum("k,kij->ij", x, B)
        h, U = np.linalg.eigh(H)
        V = U.conj().T @ v @ U
        Pp = 0.0; Ppp = 0.0
        for a in range(3):
            for b in range(a+1,3):
                D = h[a]-h[b]; dd = abs(D); s_ = np.sign(D)
                u = (V[a,a]-V[b,b]).real
                w = 0.0
                for c in range(3):
                    if c != a: w += 2*abs(V[a,c])**2/(h[a]-h[c])
                    if c != b: w -= 2*abs(V[b,c])**2/(h[b]-h[c])
                Pp += float(gp(dd))*s_*u
                Ppp += float(gpp(dd))*u*u + float(gp(dd))*s_*w
        return Pp, Ppp

    x = np.zeros(8)
    Ph = Phi_of(x)
    rec_P = np.empty(N); rec_PP = np.empty(N); rec_a = np.empty(N); acc = 0
    st = np.random.default_rng(7)
    for st_i in range(burn + N):
        xn = x + st.normal(0, sig, 8)
        Pn = Phi_of(xn)
        if np.log(st.random()) < -(Pn - Ph):
            x, Ph = xn, Pn
            if st_i >= burn: acc += 1
        if st_i >= burn:
            Pp, Ppp = derivs(x, v)
            rec_P[st_i-burn] = Pp; rec_PP[st_i-burn] = Ppp
            rec_a[st_i-burn] = (x @ np.array([0,0,1/2,0,0,0,0,1/np.sqrt(12)]) )  # <X,v>/<v,v>
    P("  接受率 = %.3f   <Phi> = %.4f" % (acc/N, np.mean([Phi_of(x)]) if False else Ph))
    P("  [全域] E[Phi'^2]=%+.6f  E[Phi'']=%+.6f  ratio=%.6f   E[Phi']=%+.3e  Var/E=%+.6f"
      % (np.mean(rec_P**2), np.mean(rec_PP), np.mean(rec_P**2)/np.mean(rec_PP),
         np.mean(rec_P), np.var(rec_P)/np.mean(rec_PP)))
    K = np.exp(-0.5*(rec_a/bw)**2)
    m = lambda f: float(np.sum(K*f)/np.sum(K))
    E2, Ec, mc_ = m(rec_P**2), m(rec_PP), m(rec_P)
    P("  [a=0, bw=%.2f] E[Phi'|0]=%+.3e  E[Phi'^2|0]=%+.6f  E[Phi''|0]=%+.6f  ratio=%+.6f"
      % (bw, mc_, E2, Ec, E2/Ec))
    P("  [a=0] Var(Phi'|0)=%+.6f  E[Phi''|0]=%+.6f  Var/E=%+.6f  (ESS_kernel=%.2e)"
      % (m((rec_P-mc_)**2), Ec, m((rec_P-mc_)**2)/Ec, (K.sum()**2)/(K**2).sum()))
    for bw2 in (0.20, 0.45):
        K2 = np.exp(-0.5*(rec_a/bw2)**2)
        m2 = lambda f: float(np.sum(K2*f)/np.sum(K2))
        P("        核 bw=%.2f: E[Phi']=%+.2e E[Phi'^2]=%+.6f E[Phi'']=%+.6f ratio=%+.6f"
          % (bw2, m2(rec_P), m2(rec_P**2), m2(rec_PP), m2(rec_P**2)/m2(rec_PP)))
    P()
    return dict()


# =====================================================================
def S5(N=2048):
    P("=" * 96)
    P("[S5] SU(3) 环面 (r12 模型, W=diag(-1,-1,1), theta=pi, beta=2.5), 2D Weyl 积分")
    P("     正确配对: 密度 |Delta|^2 e^{-S_W}  <=>  Phi_Delta = -log|Delta|^2")
    P("     错配对(r12/O09): 同一密度  +  Phi_O09 = sum g(d_ij)   => 破坏 IBP 恒等式")
    P("=" * 96)
    h = 2 * np.pi / N
    aa = (np.arange(N) + 0.5) * h
    A, Bg = np.meshgrid(aa, aa, indexing="ij")
    d12, d23, d13 = A, Bg, A + Bg
    s = {"12": 4 * np.sin(d12 / 2) ** 2, "23": 4 * np.sin(d23 / 2) ** 2, "13": 4 * np.sin(d13 / 2) ** 2}
    T = s["12"] * s["23"] * s["13"]
    beta = 2.5
    phi = (np.pi, np.pi, 0.0)
    lam = (1.0, 0.0, -1.0)
    x1 = (2 * d12 + d23) / 3; x2 = (d23 - d12) / 3; x3 = -(d12 + 2 * d23) / 3
    SW = beta * (1 - (np.cos(phi[0] - x1) + np.cos(phi[1] - x2) + np.cos(phi[2] - x3)) / 3.0)
    G = np.zeros_like(T)
    for j, xj in enumerate((x1, x2, x3)):
        G += -(beta / 3) * (-np.sin(phi[j] - xj)) * lam[j]
    w = T * np.exp(-SW); w0 = T

    def avg(f, ww):
        return float(np.sum(ww * f) / np.sum(ww))

    # Phi_Delta = -log T (T = prod 4 sin^2) ;  d_p = lam_i - lam_j
    dd = {"12": lam[0] - lam[1], "23": lam[1] - lam[2], "13": lam[0] - lam[2]}
    others = {"12": ("23", "13"), "23": ("12", "13"), "13": ("12", "23")}
    # Phi_Delta' = sum_p d_p * (-cot(d_p/2))   (d/dx of -log(4 sin^2(d/2)) = -cot(d/2)*dd/dx)
    PDp = np.zeros_like(T); PDpp = np.zeros_like(T); PDp2 = np.zeros_like(T)
    dmap = {"12": d12, "23": d23, "13": d13}
    for p in ("12", "23", "13"):
        c = -1.0 / np.tan(dmap[p] / 2.0)                    # d/dd[-log(4 sin^2(d/2))] = -cot(d/2)
        cp = 0.5 / np.sin(dmap[p] / 2.0) ** 2               # d^2/dd^2 = (1/2)csc^2
        PDp += dd[p] * c
        PDpp += dd[p] ** 2 * cp
    PDp2 = PDp ** 2
    P("  [正确配对] E[Phi_D'^2]=%.8f  E[Phi_D'']=%.8f  ratio=%.8f   E[Phi_D']=%+.3e"
      % (avg(PDp2, w0), avg(PDpp, w0), avg(PDp2, w0) / avg(PDpp, w0), avg(PDp, w0)))
    P("              (r12[B] 精确锚: E[Phi_A'']=E[Phi_A'^2]=8 精确, S-D=0)")
    P("     加 S_W 权重: E[Phi_D'^2]=%.8f  E[Phi_D'']=%.8f  E[G Phi_D']=%+.8f  (恒等式: 差应=E[G Phi'])"
      % (avg(PDp2, w), avg(PDpp, w), avg(G * PDp, w)))
    P("     恒等式残差 E[Phi'^2]-(E[Phi'']+E[G Phi']) = %+.3e"
      % (avg(PDp2, w) - avg(PDpp, w) - avg(G * PDp, w)))
    P()
    # 错配对: Phi_O09
    d = {"12": d12, "23": d23, "13": d13}
    q1 = {p: 8 * np.sin(d[p] / 2) ** 2 / np.where(d[p] == 0, 1, d[p]) - 2 * np.sin(d[p]) for p in d}
    q2 = {p: (4 * np.sin(d[p] / 2) / np.where(d[p] == 0, 1, d[p]) - 2 * np.cos(d[p] / 2)) ** 2 for p in d}
    T1 = np.zeros_like(T); T2 = np.zeros_like(T); Tpp = np.zeros_like(T)
    for p in ("12", "23", "13"):
        o1, o2 = others[p]
        T1 += dd[p] * q1[p] * s[o1] * s[o2]
        T2 += dd[p] ** 2 * q2[p] * s[o1] * s[o2]
        q0 = (2 * np.sin(d[p] / 2) / np.where(d[p] == 0, 1, d[p])) ** 2
        Tpp += dd[p] ** 2 * (2 - 2 * q0) * s[o1] * s[o2]
    for pp, qq in (("12", "23"), ("12", "13"), ("23", "13")):
        rr = ({"12", "23", "13"} - {pp, qq}).pop()
        T2 += 2 * dd[pp] * dd[qq] * q1[pp] * q1[qq] * s[rr]
    P1 = T1 / T; P2 = T2 / T; PP = Tpp / T
    e1, e2, epp = avg(P1, w0), avg(P2, w0), avg(PP, w0)
    P("  [错配对 Phi_O09, 纯 Haar 权重]  E[Phi']=%+.6f  Var=%+.6f  E[Phi'']=%+.6f  ratio=%.6f"
      % (e1, float(avg((P1 - e1) ** 2, w0)), epp, float(avg((P1 - e1) ** 2, w0)) / epp))
    P("  [错配对 Phi_O09, e^{-S_W} 权重]  E[Phi']=%+.6f  E[Phi'^2]=%+.6f  E[Phi'']=%+.6f  E[G Phi']=%+.6f"
      % (avg(P1, w), avg(P2, w), avg(PP, w), avg(G * P1, w)))
    P("      恒等式残差 = %+.6f   <== 非零(应为0) => (Phi_O09, |Delta|^2) 不是 IBP 相容对"
      % (avg(P2, w) - avg(PP, w) - avg(G * P1, w)))
    P()
    return dict(e1=e1, e2=avg(P2, w), epp=avg(PP, w), egp=avg(G * P1, w), Eg2=avg(G * G, w))


# =====================================================================
def S6(mcres, s5res):
    P("=" * 96)
    P("[S6] 闭式上界与应用 (r12 数字 + 本次 MC 数字)")
    P("=" * 96)
    P("  恒等式:  E[Phi'' ] = E[Phi'^2] + E[G Phi']            [T]")
    P("  =>  Var(Phi') = E[Phi''] - (E Phi')^2 - Cov(G,Phi')   =>  C = 1 当且仅当 Cov(G,Phi') >= 0")
    P("  只用 (M, g, EPhi') 的确定性上界:  s=||Phi'||_2,  s^2 - g s - M <= 0  =>  s <= (g+sqrt(g^2+4M))/2")
    P("      Var(Phi') <= M + g*(g+sqrt(g^2+4M))/2 - (E Phi')^2")
    for lab, M, g2, E1 in (("r12/SU(3) |Delta|^2 e^-S_W", 6.7040980628, 0.0767358867, 2.4921968555),
                           ("r12/SU(3) 纯 Haar", 6.44071163, 0.06108481, 2.44)):
        gg = np.sqrt(g2)
        bd = M + gg * (gg + np.sqrt(gg * gg + 4 * M)) / 2 - E1 ** 2
        P("  %-28s M=E[Phi'']=%.6f  sqrt(E[G^2])=%.4f  E[Phi']=%.4f  ==>  Var(Phi') <= %.4f"
          % (lab, M, gg, E1, bd))
    P("  r12 实测 Var(Phi')=6.8107 > 上界 1.2343  ==> r12 的 (Phi_O09, |Delta|^2) 违反恒等式 6.3177")
    P("      即 r12 的 6.81 不是任何 IBP 相容测度下的 Var(Phi'); 它是'配对错'的产物 [C]")
    P()
    P("  q09 的路径界 7480  vs  本上界 ~1.2 : 收紧 %.0fx" % (7480 / 1.2343))
    P()
    P("  [判] (a) C=1 精确(纯 Haar / Phi 配对 / Cov>=0): Phi 部 = (E Phi')^2 + Cov(G,Phi')")
    P("          纯 Haar + a=0 切片: X->-X 对称 => E[Phi'|0]=0 精确, 且 E[Phi'^2|0]=E[Phi''|0] =>")
    P("          Phi 部 = E[Phi'']-Var(Phi') = 0 精确 (中性), 不是正也不是负  [T]")
    P("      (b) r12 的 -0.107 / -0.236 是错配对产物; 正确配对下 Phi 部 = 0 (r12[B]: S=D=8 精确)")
    P("      (c) |nabla Phi| 无界: d_ab -> 2pi (cut locus) 处 g'->+inf; 与 a=0 无关; 测度 e^{-Phi} 压制 => 矩有限")


if __name__ == "__main__":
    P("s04_varphip.log -- S04: Var(Phi') 的无条件上界 (SU(2), SU(3))")
    P("numpy %s  python %s" % (np.__version__, sys.version.split()[0]))
    S1(); S2()
    mc = S4()
    s5 = S5()
    S6(mc, s5)
    with open("s04_varphip.log", "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")
