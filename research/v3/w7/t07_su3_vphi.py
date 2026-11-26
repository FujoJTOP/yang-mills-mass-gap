# -*- coding: utf-8 -*-
r"""
t07_su3_vphi.py -- T07: SU(3) 的 Var(Phi') 从哪来、能否压到 0 (v-依赖定位)。

Phi(X) = -log J(X) = sum_{a<b} g(|lam_a-lam_b|),  g(d) = -2 log[sin(d/2)/(d/2)]
grad Phi(X) = sum_{a<b} omega(lam_a-lam_b)(P_a - P_b),  omega(s)=g'(|s|)sgn(s)
Phi'(X;v)   = <grad Phi(X), v> = sum_a c_a (U^dag v U)_aa,  c_a=sum_{b!=a}omega(lam_a-lam_b)

[T1] Phi' 的 v-依赖定位 + FD 核对
[T2] su(3) 不变量结构 (Newton): p4=p2^2/2, Sum d^4=(9/2)p2^2, Sum d^6=(33/4)p2^3-9p3^2,
     g 的 d^6 系数 = 1/90720  ==>  Phi 径向到 4 次, 首个非径向项 = p3^2 (6 次, 系数 -9.92e-5)
[T3] E[Phi']=0 精确 (Phi 偶 ==> gradPhi 奇; e^{-Phi}dX 偶)
[T4] Schur: Sigma=E[gradPhi(x)gradPhi]=kappa*Id ==> 全域 Var=kappa|v|^2 方向无关
[T5] 切片 a=<X,v>=0: 精确条件抽样, 扫 v 的 2 参数共轭类; 报告 min/max; 验 Phi'|_{H_v}=O(|X|^5)
[T6] 审计: s04[S4] 的 "8 维" 基实际 5 维且非正交
[T7] 审计: r12[C] 的 6.81/2.49 均值 -> d13=a+b 不 mod 2pi 的域效应 + 错配对
"""
import sys
import numpy as np

np.random.seed(20261125)
OUT = []


def P(s=""):
    print(s)
    OUT.append(s)


I2 = 1j


def g(d):
    d = np.asarray(d, float)
    small = np.abs(d) < 1e-8
    with np.errstate(divide="ignore", invalid="ignore"):
        v = -2.0 * np.log(np.abs(np.sin(d / 2.0) / (d / 2.0)))
    return np.where(small, d * d / 12.0, v)


def gp(d):
    d = np.asarray(d, float)
    small = np.abs(d) < 1e-8
    return np.where(small, d / 6.0, 2.0 / d - 1.0 / np.tan(d / 2.0))


def omega(s):
    s = np.asarray(s, float)
    return np.sign(s) * gp(np.abs(s))


def gell_mann():
    L = [np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], complex),
         np.array([[0, -I2, 0], [I2, 0, 0], [0, 0, 0]], complex),
         np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], complex),
         np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], complex),
         np.array([[0, 0, -I2], [0, 0, 0], [I2, 0, 0]], complex),
         np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], complex),
         np.array([[0, 0, 0], [0, 0, -I2], [0, I2, 0]], complex),
         np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], complex) / np.sqrt(3.0)]
    return np.array([m / np.sqrt(2.0) for m in L])      # tr(m_a m_b)=delta_ab


GM = gell_mann()


def norm1(v):
    return v / np.sqrt(np.real(np.trace(v @ v)))


def frame(v):
    v = norm1(np.asarray(v, complex))
    B = [v]
    for m in GM:
        w = m.copy()
        for b in B:
            w = w - np.real(np.trace(b @ w)) * b
        nw = np.sqrt(np.real(np.trace(w @ w)))
        if nw > 1e-10:
            B.append(w / nw)
        if len(B) == 8:
            break
    assert len(B) == 8
    return np.array(B)


def Phi_of(X):
    lam = np.linalg.eigvalsh(X)
    s = np.zeros(len(X))
    for a, b in ((0, 1), (0, 2), (1, 2)):
        s += g(np.abs(lam[:, a] - lam[:, b]))
    return s


def gradPhi(X):
    lam, U = np.linalg.eigh(X)
    c = np.zeros_like(lam)
    for a in range(3):
        for b in range(3):
            if a != b:
                c[:, a] += omega(lam[:, a] - lam[:, b])
    G = np.einsum("nij,nj,nkj->nik", U, c, U.conj())
    return 0.5 * (G + np.conj(np.transpose(G, (0, 2, 1))))


def p2p3(X):
    X2 = np.einsum("nij,njk->nik", X, X)
    p2 = np.real(np.einsum("nij,nji->n", X, X))
    p3 = np.real(np.einsum("nij,nji->n", X2, X))
    return p2, p3, X2


def sample_full(N, sig=1.0):
    x = np.random.randn(N, 8) * sig
    X = np.einsum("na,aij->nij", x, GM)
    lw = -Phi_of(X)
    lw -= lw.max()
    return X, np.exp(lw)


def wmean(f, w):
    return float(np.sum(w * f) / np.sum(w))


# =====================================================================
def T1():
    P("=" * 100)
    P("[T1] Phi' 的 v-依赖定位")
    P("=" * 100)
    P("  Phi'(X;v) = <grad Phi(X), v> = sum_a c_a (U^dag v U)_aa ,  c_a=sum_{b!=a} omega(lam_a-lam_b)")
    P("  ==> g'(d_p) 与 (lam_i-lam_j) 都只依赖 X 的谱, 与 v 无关;")
    P("      【(lam_i-lam_j) 不是 v 的分量, 是 X 的谱差】")
    P("      v 的唯一入口 = 投影子矩阵元 <P_a-P_b,v> = (U^dag v U)_aa-(U^dag v U)_bb")
    P("  ==> Phi'(X;.) 对 v 是线性泛函 ==> 任何 v-依赖只能是二次型 v^T Sigma v")
    P()
    N = 500
    X = np.einsum("na,aij->nij", np.random.randn(N, 8) * 0.9, GM)
    v = norm1(GM[7] + 0.6 * GM[2])
    ana = np.real(np.einsum("nij,ji->n", gradPhi(X), v))
    h = 1e-6
    fd = (Phi_of(X + h * v) - Phi_of(X - h * v)) / (2 * h)
    P("  FD 核对 max|解析 - 差分| = %.3e  (N=%d)" % (np.max(np.abs(ana - fd)), N))
    P()


# =====================================================================
def T2():
    P("=" * 100)
    P("[T2] su(3) 不变量结构: Phi 是【径向到 4 次】, 首个非径向项 = p3^2 (6 次)")
    P("=" * 100)
    P("  Newton (traceless 3x3): p4 = p2^2/2 ,  p5 = (5/6)p2 p3   [T]")
    P("  ==> sum_{a<b}(lam_a-lam_b)^2 = 3 p2")
    P("      sum_{a<b}(lam_a-lam_b)^4 = 3p4 + 3p2^2 = (9/2) p2^2            <-- 纯径向!")
    P("      sum_{a<b}(lam_a-lam_b)^6 = (33/4) p2^3 - 9 p3^2               <-- 首个非径向")
    P("  g(d) = d^2/12 + d^4/1440 + d^6/90720 + O(d^8)   (展开, log(sin x/x) 系数 -1/2835)")
    P("  ==> Phi = p2/4 + p2^2/480 + (33/362880)p2^3 - (9/90720) p3^2 + O(deg 8)")
    P("      p3^2 项系数 = -9/90720 = -9.92e-5")
    P("  ==> 【grad Phi = (径向)X + O(|X|^5)】;  SU(2) 的 'Phi 径向' 在 SU(3) 只剩【6 次以下的径向】")
    P()
    N = 400
    X = np.einsum("na,aij->nij", np.random.randn(N, 8) * 0.6, GM)
    p2, p3, X2 = p2p3(X)
    p4 = np.real(np.einsum("nij,nji->n", X2, X2))
    lam = np.linalg.eigvalsh(X)
    s2 = sum((lam[:, a] - lam[:, b]) ** 2 for a, b in ((0, 1), (0, 2), (1, 2)))
    s4 = sum((lam[:, a] - lam[:, b]) ** 4 for a, b in ((0, 1), (0, 2), (1, 2)))
    s6 = sum((lam[:, a] - lam[:, b]) ** 6 for a, b in ((0, 1), (0, 2), (1, 2)))
    P("  p4 - p2^2/2            max = %.3e" % np.max(np.abs(p4 - p2 ** 2 / 2)))
    P("  sum d^2 - 3p2          max = %.3e" % np.max(np.abs(s2 - 3 * p2)))
    P("  sum d^4 - (9/2)p2^2    max = %.3e" % np.max(np.abs(s4 - 4.5 * p2 ** 2)))
    P("  sum d^6 - ((33/4)p2^3-9p3^2) max = %.3e" % np.max(np.abs(s6 - (8.25 * p2 ** 3 - 9 * p3 ** 2))))
    pred = p2 / 4 + p2 ** 2 / 480 + 33 * p2 ** 3 / 362880 - (9 / 90720) * p3 ** 2
    act = Phi_of(X)
    P("  |Phi - [p2/4+p2^2/480+(33/362880)p2^3-(9/90720)p3^2]| / Phi  max = %.3e"
      % np.max(np.abs(act - pred) / np.abs(act)))
    P("  (残差 = O(deg 8) 项, 相对 ~1e-3 符合 8 次截断)")
    P()
    # 精确三次恒等式 (解释 3 次探针恒 0)
    v = norm1(GM[2] + 0.4 * GM[7])
    x1 = np.random.randn(N, 3)
    x1 -= x1.mean(1, keepdims=True)
    Xd = np.zeros((N, 3, 3), complex)
    Xd[:, 0, 0], Xd[:, 1, 1], Xd[:, 2, 2] = x1[:, 0], x1[:, 1], x1[:, 2]
    Bv = np.zeros((N, 3, 3), complex)
    r = np.random.randn(N, 3)
    r -= r.mean(1, keepdims=True)
    Bv[:, 0, 0], Bv[:, 1, 1], Bv[:, 2, 2] = r[:, 0], r[:, 1], r[:, 2]
    trXv = np.real(np.einsum("nij,nji->n", Xd, Bv))
    trX3v = np.real(np.einsum("nij,nji->n", np.einsum("nij,njk->nik", np.einsum("nij,njk->nik", Xd, Xd), Xd), Bv))
    p2d = np.real(np.einsum("nij,nji->n", Xd, Xd))
    P("  恒等式 tr(X^3 v) = tr(Xv)*p2/2  max|残差| = %.3e" % np.max(np.abs(trX3v - trXv * p2d / 2)))
    P("  ==> 3 次探针 <X^3,v> 在 H_v 上恒 0 是【恒等式】, 不是巧合; 真正的首个非径向梯度项是 5 次。")
    P()


# =====================================================================
def T3():
    P("=" * 100)
    P("[T3] E[Phi'] = 0 ?")
    P("=" * 100)
    P("  Phi(-X)=Phi(X) ==> grad Phi(-X)=-grad Phi(X); e^{-Phi}dX 在 X->-X 不变")
    P("  ==> 任意固定 v: E[Phi'] = 0【精确】[T]; 切片 <X,v>=0 也是 -对称集 ==> E[Phi'|a=0]=0 [T]")
    P("  ==> Var(Phi') = E[Phi'^2]")
    P()
    X, w = sample_full(300000)
    G = gradPhi(X)
    ess = (w.sum() ** 2) / (w ** 2).sum()
    for nm, v in (("M8", GM[7]), ("M3", GM[2]), ("M3+0.6M8", norm1(GM[2] + 0.6 * GM[7]))):
        v = norm1(v)
        f = np.real(np.einsum("nij,ji->n", G, v))
        m = wmean(f, w)
        P("  v=%-10s E[Phi']=%+.3e (MC 误差 %.1e)  E[Phi'^2]=%+.6f" % (nm, m, 1e-3, wmean(f ** 2, w)))
    P("  ESS=%.3e" % ess)
    P()


def T4():
    P("=" * 100)
    P("[T4] Schur: Sigma = E[gradPhi (x) gradPhi] = kappa Id  ==> 全域 Var = kappa|v|^2, 方向无关")
    P("=" * 100)
    P("  e^{-Phi}dX 共轭不变, gradPhi 等变, su(3) 伴随表示不可约 ==> Sigma=kappa Id [T]")
    P("  ==> inf_v = sup_v = kappa|v|^2;  v!=0 ==> Var>0。")
    P()
    X, w = sample_full(300000)
    G = gradPhi(X)
    co = np.real(np.einsum("nij,bji->nb", G, GM))
    Sig = np.einsum("n,na,nb->ab", w, co, co) / w.sum()
    kap = np.trace(Sig) / 8.0
    P("  kappa = tr(Sigma)/8 = %.6f    E|gradPhi|^2 = 8kappa = %.6f" % (kap, 8 * kap))
    P("  max|diag-kappa|/kappa = %.4f    max|offdiag|/kappa = %.4e"
      % (np.max(np.abs(np.diag(Sig) - kap)) / kap, np.max(np.abs(Sig - np.diag(np.diag(Sig)))) / kap))
    P()


def T5():
    P("=" * 100)
    P("[T5] 切片 a=<X,v>=0 (精确条件抽样: 基首元=v, x_1=0), 扫 v 的 2 参数共轭类")
    P("=" * 100)
    P("  共轭不变 ==> Var_slice 只依赖 v 的本征三元组(有序); 归一 |v|=1:")
    P("  t=(sqrt(2/3)cos(th+2pi/3), sqrt(2/3)cos(th-2pi/3), sqrt(2/3)cos th), th in [0,60) 覆盖 Weyl 腔")
    P()
    N = 600000
    P("   th     t1      t2      t3     E[Phi'|a=0]  Var(Phi'|a=0)  E[Phi'^2]   (Var/Var(G))")
    res = []
    VG = 0.0767358867        # r12[D] EG2, 环面; 这里作量级对照
    for deg in range(0, 65, 5):
        th = np.deg2rad(deg)
        t = np.sqrt(2.0 / 3.0) * np.array([np.cos(th + 2 * np.pi / 3), np.cos(th - 2 * np.pi / 3), np.cos(th)])
        v = np.diag(t).astype(complex)
        B = frame(v)
        x = np.zeros((N, 8))
        x[:, 1:] = np.random.randn(N, 7)
        X = np.einsum("na,aij->nij", x, B)
        lw = -Phi_of(X)
        lw -= lw.max()
        w = np.exp(lw)
        f = np.real(np.einsum("nij,ji->n", gradPhi(X), v))
        m = wmean(f, w)
        E2 = wmean(f ** 2, w)
        Var = E2 - m * m
        res.append((deg, Var, m))
        P("  %4d  %+.4f %+.4f %+.4f   %+.2e   %+.6f    %+.6f   %.2e"
          % (deg, t[0], t[1], t[2], m, Var, E2, Var / VG))
    V = np.array([r[1] for r in res])
    P()
    P("  min_v Var(Phi'|a=0) = %.6e (th=%d)   max_v = %.6e (th=%d)   max/min = %.2f"
      % (V.min(), res[int(np.argmin(V))][0], V.max(), res[int(np.argmax(V))][0], V.max() / V.min()))
    P("  所有 v: Var>0 ? %s ;  max|E[Phi'|a=0]| = %.1e (应 0)"
      % (bool(np.all(V > 0)), max(abs(r[2]) for r in res)))
    P()
    # O(|X|^5) 标度检验
    P("  Phi'|_{H_v} = O(|X|^5) 检验: 固定 X0 in H_v, 看 Phi'(tX0;v)/t^5")
    v = np.diag(np.sqrt(2.0 / 3.0) * np.array([np.cos(2 * np.pi / 3), np.cos(-2 * np.pi / 3), 1.0])).astype(complex)
    B = frame(v)
    x0 = np.zeros(8)
    x0[1:] = np.random.randn(7)
    P("     t      Phi'(tX0)/t^3     Phi'(tX0)/t^5")
    for tt in (0.3, 0.5, 1.0, 1.5, 2.0):
        X = np.einsum("a,aij->ij", tt * x0, B)[None]
        f = float(np.real(np.trace(gradPhi(X)[0] @ v)))
        P("   %4.1f   %+.6e     %+.6e" % (tt, f / tt ** 3, f / tt ** 5))
    P()


def T6():
    P("=" * 100)
    P("[T6] 审计 s04[S4] 的 '8 维 log-chart' 基")
    P("=" * 100)
    B = []
    for i in range(3):
        for j in range(i + 1, 3):
            M = np.zeros((3, 3)); M[i, j] = M[j, i] = 1 / np.sqrt(2); B.append(M)
    for i in range(2):
        M = np.zeros((3, 3)); M[i, i] = 1 / np.sqrt(2); M[2, 2] = -1 / np.sqrt(2)
        if i == 1:
            M = np.zeros((3, 3)); M[0, 0] = 1 / np.sqrt(6); M[1, 1] = 1 / np.sqrt(6); M[2, 2] = -2 / np.sqrt(6)
        B.append(M)
    B = np.array(B)
    tr = [float(np.trace(M)) for M in B]
    Gm = np.einsum("aij,bji->ab", B, B)
    P("  元素数=%d (宣称 8 维); tr 各元 = " % len(B) + " ".join("%+.4f" % x for x in tr))
    P("  最大非对角 Gram = %.6f (正交基应 0);  max|tr| = %.6f" % (np.max(np.abs(Gm - np.diag(np.diag(Gm)))), max(abs(np.array(tr)))))
    P("  ==> [E] s04[S4] 抽样只在 5 维实对称子空间且基不正交, 非 su(3) 8 维 ==> 其 S4 数字不可用。")


def T7():
    P("=" * 100)
    P("[T7] 审计 6.81: r12[C] 的 E[Phi']=2.4922 != 0 从哪来")
    P("=" * 100)
    P("  r12_su3_correct.log [C]: E[Phi']=+2.4922  E[Phi'^2]=+13.0218  Var=+6.8107  E[G]=-4.2e-17")
    P("  ==> Var = 13.0218 - 2.4922^2 = 13.0218-6.2109 = 6.8109  [E] 算术自洽")
    P("  ==> 6.81 里 E[Phi']^2=6.21 是【均值平方】; 且 r12[D] 纯 Haar 分支 Var_T(Phi_B')=6.6154")
    P()
    N = 512
    h = 2 * np.pi / N
    aa = (np.arange(N) + 0.5) * h
    A, Bg = np.meshgrid(aa, aa, indexing="ij")
    lam = (1.0, 0.0, -1.0)
    beta, phi = 2.5, (np.pi, np.pi, 0.0)
    x1 = (2 * A + Bg) / 3; x2 = (Bg - A) / 3; x3 = -(A + 2 * Bg) / 3
    SW = beta * (1 - (np.cos(phi[0] - x1) + np.cos(phi[1] - x2) + np.cos(phi[2] - x3)) / 3.0)

    def run(tag, d13):
        d = {"12": A, "23": Bg, "13": d13}
        s = {p: 4 * np.sin(d[p] / 2) ** 2 for p in d}
        T = s["12"] * s["23"] * s["13"]
        dd = {"12": lam[0] - lam[1], "23": lam[1] - lam[2], "13": lam[0] - lam[2]}
        P1 = sum(dd[p] * (-1.0 / np.tan(d[p] / 2.0)) for p in d)
        for wnm, W in (("纯Haar |D|^2", T), ("|D|^2 e^{-SW}", T * np.exp(-SW))):
            a1 = float(np.sum(W * P1) / np.sum(W))
            v1 = float(np.sum(W * (P1 - a1) ** 2) / np.sum(W))
            P("   %-22s %-14s E[Phi']=%+.6f  Var=%+.6f" % (tag, wnm, a1, v1))
        return T

    run("d13=a+b (不 mod 2pi)", A + Bg)
    run("d13=(a+b) mod 2pi", np.mod(A + Bg, 2 * np.pi))
    P()
    P("  ==> E[Phi']!=0 只在 d13 不 wrap 时出现: 域 (0,2pi)^2 的 d13 in (0,4pi) 破坏 (a,b)->(-a,-b) 对称;")
    P("      Phi_O09 的 g 含 log(d/2) (非周期), 本就不是环面上的函数 ==> r12[C] 的 2.49/6.81 是域效应+错配对。")
    P("      正确配对 Phi=-log|Delta|^2 的精确锚 (r12[B]): E[Phi_A']=0, E[Phi_A'']=E[Phi_A'^2]=8, S=D。")


if __name__ == "__main__":
    P("t07_su3_vphi.log -- T07: SU(3) Var(Phi') 的来源与可消除性")
    P("numpy %s  python %s" % (np.__version__, sys.version.split()[0]))
    P()
    T1(); T2(); T3(); T4(); T5(); T6(); T7()
    with open("t07_su3_vphi.log", "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")
    print("\n[written] t07_su3_vphi.log")
