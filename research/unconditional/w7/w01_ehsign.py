# -*- coding: utf-8 -*-
r"""
w01_ehsign.py -- W01: 裁决 E[H_vv] 在 theta=pi 上的符号 (p05 的 [T] vs v04 的 [C])

固定定义 (AGENTS §2.1, 先定死再用):

  [S1] S_W = beta * (1 - (1/N) Re tr U_p),  U_p = W e^{iX},  X in Cartan (对角).
       phi = log 的本征相: W = diag(e^{i phi_j}), sum phi = 0 (mod 2pi).
  [S2] v = 方向 lambda in Cartan,  sum lam_j = 0,  ||lam||^2 = 2 (归一).
       tau 是一维流:  x_j -> x_j + tau * lam_j .
       这是【离开切片】的导数, 与 r12 §1 的 G = d_tau S_W 同一约定.
  [S3] H_vv := d_tau^2 S_W (同一 tau, 同一 X).  <-- 符号由 S2/S1 唯一确定, 不是自由选择.
       闭式:  H_vv = +(beta/N) sum_j lam_j^2 cos(phi_j - x_j)
       (检验: theta=0 (phi=0), x=0  =>  H_vv = (beta/N)||lam||^2 > 0, 真空是极小 -> 正号正确)
  [S4] 测度:
       (M1) 无条件: 2D Weyl cell (a,b) in (0,2pi)^2,  w = |Delta|^2 e^{-S_W},  |Delta|^2 = s12 s23 s13
       (M2) 条件 (p05 型): 切片 <lam,X> = 0 上的诱导测度
  [S5] SU(2): 全 su(2) (3 维), v = e_3 = sigma_3/2, 切片 a = 2 tr(E X) = x_3 = 0,
       测度 = Haar 的切片 (sin^2 u du dpsi) x e^{-S_W}

输出:
  [A] SU(2), theta=pi: E[H_vv] 精确 = 0 ? (解析恒等式 + FD 数值)
  [B] SU(2) 锚点: theta=0, beta=2.5 -> 0.553061 (p05 §2 核对)
  [C] SU(3), theta=pi, beta=2.5, 无条件测度 M1: 三个 c_j = <cos x_j>_w
      => 对【任意】lam 给 E[H_vv](lam) = (beta/3)(lam3^2 c3 - lam1^2 c1 - lam2^2 c2)  (λ 线性!)
  [D] 区间 (libmpi, PREC=200, 向外舍入) + Richardson 外推
  [E] 与 p05 / v04 / r12 / s02 的对账 (两个符号分支)
"""
import math
import time
from fractions import Fraction as Fr

import numpy as np
import mpmath as mp
import mpmath.libmp as M
import mpmath.libmp.libmpi as LI

mp.mp.dps = 40
PREC = 200
# 2^-100 外扩 (再保险), 与 s02 同口径
PAD = M.mpf_pow(M.from_int(2, PREC), M.from_int(-(PREC // 2), PREC), PREC, 'c')

LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


# ===================== 区间引擎 (libmpi, 真向外舍入) =====================
def frac(p, q):
    return (M.from_rational(p, q, PREC, 'f'), M.from_rational(p, q, PREC, 'c'))


def mkint(n):
    return (M.from_int(n, PREC), M.from_int(n, PREC))


def ipi():
    return (M.mpf_pi(PREC, 'f'), M.mpf_pi(PREC, 'c'))


def iadd(a, b):
    return LI.mpi_add(a, b, PREC)


def isub(a, b):
    return LI.mpi_sub(a, b, PREC)


def imul(a, b):
    return LI.mpi_mul(a, b, PREC)


def idiv(a, b):
    return LI.mpi_div(a, b, PREC)


def ineg(a):
    return (M.mpf_neg(a[1], PREC), M.mpf_neg(a[0], PREC))


def icos(a):
    return LI.mpi_cos(a, PREC)


def isin(a):
    return LI.mpi_sin(a, PREC)


def iexp(a):
    return LI.mpi_exp(a, PREC)


def ipad(a):
    lo, hi = a
    w = M.mpf_sub(hi, lo, PREC, 'c')
    m = M.mpf_add(M.mpf_add(M.mpf_abs(lo, PREC), M.mpf_abs(hi, PREC), PREC, 'c'), w, PREC, 'c')
    p = M.mpf_mul(m, PAD, PREC, 'c')
    return (M.mpf_sub(lo, p, PREC, 'f'), M.mpf_add(hi, p, PREC, 'c'))


def wd(a):
    return M.mpf_sub(a[1], a[0], PREC, 'c')


def num(a, nd=20):
    return M.to_str(a[0], nd), M.to_str(a[1], nd)


# ============================================================
# [A] SU(2), theta=pi:  E[H_vv] 是否精确 0
# ============================================================
def su2_Hvv_fd(x, lam, beta, phi):
    """对 v = (x . sigma/2), 沿 x -> x + tau*lam 的 d_tau^2 S_W 的中心差分."""
    def SW(xx):
        # U_p = W e^{iX}, W = diag(e^{i phi0}, e^{-i phi0})  (SU(2), phi = (phi0,-phi0))
        # e^{iX} = cos(r/2) + i sin(r/2) (xhat . sigma),  r = |x|
        r = math.sqrt(sum(t * t for t in xx))
        if r < 1e-14:
            c, s = 1.0, 0.0
            nh = (0.0, 0.0, 0.0)
        else:
            c, s = math.cos(r / 2), math.sin(r / 2)
            nh = (xx[0] / r, xx[1] / r, xx[2] / r)
        # Re tr(W U) = 2[cos(phi0) c - sin(phi0) s nh_3]
        return beta * (1.0 - (math.cos(phi[0]) * c - math.sin(phi[0]) * s * nh[2]))

    h = 1e-4
    xp = [x[i] + h * lam[i] for i in range(3)]
    xm = [x[i] - h * lam[i] for i in range(3)]
    return (SW(xp) - 2 * SW(x) + SW(xm)) / (h * h)


def partA():
    P("=" * 78)
    P("[A] SU(2), theta=pi  (phi0 = pi/2),  H_vv = d_tau^2 S_W,  v = e_3,  切片 x_3 = 0")
    P("    闭式 (解析): H_vv = beta * cos(theta/2) * sin(u)/(4u),  u = |x|/2   [T]")
    P("    theta=pi => cos(theta/2)=0 => H_vv = 0 【逐点恒等】, 与测度/与 beta 无关")
    beta = 2.5
    lam = (0.0, 0.0, 1.0)
    for th, tag in ((math.pi, "theta=pi"), (0.0, "theta=0"), (0.7 * math.pi, "theta=0.7pi")):
        phi = (th / 2, -th / 2)
        # 切片 x_3 = 0 上的网格 (u, psi); Haar: sin^2 u du dpsi, u=|x|/2 in (0,pi)
        Nu, Np = 400, 200
        du, dp = math.pi / Nu, 2 * math.pi / Np
        accW = 0.0
        accH = 0.0
        accHan = 0.0
        for iu in range(Nu):
            u = (iu + 0.5) * du
            r = 2 * u
            for ip in range(Np):
                psi = (ip + 0.5) * dp
                x = (r * math.cos(psi), r * math.sin(psi), 0.0)
                wgt = math.sin(u) ** 2
                SWv = beta * (1.0 - (math.cos(phi[0]) * math.cos(u)))
                wgt *= math.exp(-SWv)
                hnum = su2_Hvv_fd(x, lam, beta, phi)
                han = beta * math.cos(th / 2) * math.sin(u) / (4 * u)
                accW += wgt
                accH += wgt * hnum
                accHan += wgt * han
        EH = accH / accW
        EHan = accHan / accW
        P("    %-12s  E[H_vv]_FD = %+.15e   E[H_vv]_closed = %+.15e   (|diff| = %.1e)"
          % (tag, EH, EHan, abs(EH - EHan)))
    P("    => theta=pi: 0 是【精确恒等式】, 不是数值小量; theta=0 复现 p05 的 +0.553061")
    P()


# ============================================================
# [B] SU(2) 锚点 (p05 §2)
# ============================================================
def partB():
    P("=" * 78)
    P("[B] SU(2) 锚点: theta=0, beta=2.5 -> E[H_vv] 应为 +0.553061 (p05 §2 / O21)")
    beta = 2.5
    Nu, Np = 600, 60
    du, dp = math.pi / Nu, 2 * math.pi / Np
    accW = accH = 0.0
    for iu in range(Nu):
        u = (iu + 0.5) * du
        w0 = math.sin(u) ** 2 * math.exp(-beta * (1 - math.cos(u)))
        accW += w0 * Np
        accH += w0 * Np * (beta * math.sin(u) / (4 * u))
    P("    E[H_vv](theta=0, beta=2.5) = %.9f   [p05: 0.553061]  diff = %.1e"
      % (accH / accW, abs(accH / accW - 0.553061)))
    P()


# ============================================================
# [C]/[D] SU(3), theta=pi, beta=2.5 —— 无条件 Weyl cell 测度 M1
# ============================================================
LAM1 = (1.0, 0.0, -1.0)
PHI3 = (math.pi, math.pi, 0.0)


def cell_beta(h, chunk=256):
    """np.chunked 中点求和: 返回 (Z, Zc1, Zc2, Zc3, Zh) —— 全用【规范】H = +(beta/3)sum lam^2 cos."""
    N = int(round(2 * math.pi / h))
    h = 2 * math.pi / N
    g = (np.arange(N) + 0.5) * h
    beta = 2.5
    acc = np.zeros(5)
    for i0 in range(0, N, chunk):
        a = g[i0:i0 + chunk][:, None]
        b = g[None, :]
        d12 = a + 0 * b
        d23 = 0 * a + b
        d13 = a + b
        s12 = 2 - 2 * np.cos(d12)
        s23 = 2 - 2 * np.cos(d23)
        s13 = 2 - 2 * np.cos(d13)
        T = s12 * s23 * s13
        x1 = (2 * d12 + d23) / 3
        x2 = (d23 - d12) / 3
        x3 = -(d12 + 2 * d23) / 3
        cs = np.cos(PHI3[0] - x1) + np.cos(PHI3[1] - x2) + np.cos(PHI3[2] - x3)
        SW = beta * (1 - cs / 3)
        w = T * np.exp(-SW)
        acc[0] += w.sum()
        acc[1] += (w * np.cos(x1)).sum()
        acc[2] += (w * np.cos(x2)).sum()
        acc[3] += (w * np.cos(x3)).sum()
        # 规范 H_vv for lam=(1,0,-1):  +(beta/3)*(lam1^2 cos(phi1-x1)+0+lam3^2 cos(phi3-x3))
        H = (beta / 3) * ((1.0) * np.cos(PHI3[0] - x1) + (1.0) * np.cos(PHI3[2] - x3))
        acc[4] += (w * H).sum()
    return acc, N


def partCD():
    P("=" * 78)
    P("[C] SU(3), theta=pi (W=diag(-1,-1,1), phi=(pi,pi,0)), beta=2.5, 无条件 Weyl cell (M1)")
    P("    w = |Delta|^2 e^{-S_W},  S_W = beta(1-(1/3) sum_j cos(phi_j-x_j))")
    P("    规范 H_vv = +(beta/3) sum_j lam_j^2 cos(phi_j - x_j)   [= d_tau^2 S_W, 由真真空凸性定号]")
    P("    => H_vv(lam) 对 lam 是【线性】: E[H_vv] = (beta/3)( lam3^2 c3 - lam1^2 c1 - lam2^2 c2 )")
    P("       c_j := <cos x_j>_w  (phi=(pi,pi,0): cos(phi_j-x_j) = -cos x_j (j=1,2), +cos x_j (j=3))")
    P()
    rows = []
    for N in (128, 256, 512, 1024, 2048):
        t0 = time.time()
        acc, NN = cell_beta(2 * math.pi / N)
        Z = acc[0]
        c1, c2, c3 = acc[1] / Z, acc[2] / Z, acc[3] / Z
        EH_direct = acc[4] / Z
        b = 2.5
        # 三个 lam 的预测
        def EH(l):
            return (b / 3) * (l[2] ** 2 * c3 - l[1] ** 2 * c2 - l[0] ** 2 * c1)
        lA = (1.0, 0.0, -1.0)
        lB = (1.0, -1.0, 0.0)
        lC = (0.0, 1.0, -1.0)
        lD = (1 / math.sqrt(2), 1 / math.sqrt(2), -math.sqrt(2))
        rows.append((NN, c1, c2, c3, EH_direct, EH(lA), EH(lB), EH(lC), EH(lD)))
        P("    N=%-5d (%.1fs)  c1=%+.10f c2=%+.10f c3=%+.10f" % (NN, time.time() - t0, c1, c2, c3))
        P("            E[H_vv] 直接=%+.10f | lam(1,0,-1)=%+.10f  lam(1,-1,0)=%+.10f"
          % (EH_direct, EH(lA), EH(lB)))
        P("            lam(0,1,-1)=%+.10f  lam(1,1,-2)/sqrt6=%+.10f" % (EH(lC), EH(lD)))
    # Richardson (O(h^2))
    P()
    P("    —— O(h^2) Richardson 外推 (N 加倍) ——")
    for i in range(len(rows) - 1):
        N0, N1 = rows[i][0], rows[i + 1][0]
        for k, nm in ((4, "E[H]|lam=(1,0,-1)"),):
            f0, f1 = rows[i][k], rows[i + 1][k]
            Ric = f1 + (f1 - f0) / 3.0
            P("      %-20s N=%-5d %+.12f  N=%-5d %+.12f  Richardson %+.12f  (相邻差/3=%.1e)"
              % (nm, N0, f0, N1, f1, Ric, abs(f1 - f0) / 3))
    return rows


def partD_iv(N):
    """libmpi 区间: Z, Zc1..Zc3, ZH, 全部向外舍入 + 2^-100 pad."""
    pi = ipi()
    h = idiv(imul(mkint(2), pi), mkint(N))
    beta = frac(5, 2)
    three = mkint(3)
    Z = (M.fzero, M.fzero)
    Zc = [(M.fzero, M.fzero) for _ in range(3)]
    ZH = (M.fzero, M.fzero)
    ph = [idiv(pi, mkint(1)), idiv(pi, mkint(1)), (M.fzero, M.fzero)]
    lam2 = [frac(1, 1), frac(0, 1), frac(1, 1)]
    for i in range(N):
        a = imul(h, frac(2 * i + 1, 2))
        for j in range(N):
            b = imul(h, frac(2 * j + 1, 2))
            d12, d23 = a, b
            d13 = iadd(a, b)
            cs = []
            for d in (d12, d23, d13):
                cs.append(icos(d))
            s = [isub(mkint(2), imul(mkint(2), c)) for c in cs]
            T = imul(imul(s[0], s[1]), s[2])
            x1 = idiv(iadd(imul(mkint(2), d12), d23), three)
            x2 = idiv(isub(d23, d12), three)
            x3 = ineg(idiv(iadd(d12, imul(mkint(2), d23)), three))
            xs = (x1, x2, x3)
            sum_c = (M.fzero, M.fzero)
            for k in range(3):
                sum_c = iadd(sum_c, icos(isub(ph[k], xs[k])))
            SW = imul(beta, isub(mkint(1), idiv(sum_c, three)))
            w = imul(T, iexp(ineg(SW)))
            Z = iadd(Z, w)
            for k in range(3):
                Zc[k] = iadd(Zc[k], imul(w, icos(xs[k])))
            H = (M.fzero, M.fzero)
            for k in range(3):
                H = iadd(H, imul(lam2[k], icos(isub(ph[k], xs[k]))))
            ZH = iadd(ZH, imul(w, imul(idiv(beta, three), H)))
    Z, ZH = ipad(Z), ipad(ZH)
    Zc = [ipad(z) for z in Zc]
    return Z, Zc, ZH


def partD():
    P("=" * 78)
    P("[D] 区间 (libmpi PREC=%d, 中点求和, 节点取区间, 向外舍入 + 2^-100 pad)" % PREC)
    P("    只包住【中点和】; 求积截断另计 (O(h^2), 见 [C] 的 Richardson) => 总标签 [C], 不是 [T]")
    for N in (32, 64):
        t0 = time.time()
        Z, Zc, ZH = partD_iv(N)
        P("    N=%-4d (%.1fs)" % (N, time.time() - t0))
        P("      Z    [%s, %s] w=%s" % (num(Z, 16) + (M.to_str(wd(Z), 3),)))
        for k in range(3):
            c = idiv(Zc[k], Z)
            lo, hi = num(c, 16)
            P("      c_%d  [%s, %s] w=%s" % (k + 1, lo, hi, M.to_str(wd(c), 3)))
        e = idiv(ZH, Z)
        lo, hi = num(e, 18)
        P("      E[H_vv]|lam=(1,0,-1)  规范支  [%s, %s]  w=%s" % (lo, hi, M.to_str(wd(e), 3)))
        P("      E[H_vv]|lam=(1,0,-1)  s02 支 (= -规范)  [%s, %s]"
          % num(ineg(e), 18))
    P()


def partE(rows):
    P("=" * 78)
    P("[E] 对账")
    N, c1, c2, c3, EHd, EhA, EhB, EhC, EhD = rows[-1]
    P("    规范支 (H_vv = d_tau^2 S_W):  E[H_vv] | lam=(1,0,-1) = %+.11f" % EhA)
    P("    s02/v04 支 (H = -(beta/3)sum lam^2 cos = -H_vv):  %+.11f  <-> v04 记 -0.10230219" % (-EhA))
    P("    r12 §3.2 记的 E[H] = +0.1023  = 规范支 (与 v04 只差一个整体符号)   => |差| = %.1e"
      % abs(EhA - 0.1023021930))
    P("    SU(2), theta=pi: E[H_vv] = 0 精确  (v04 §4 表也这么写: '0 (精确, 中心元)')")
    P("    SU(3) 的 S_W 在 theta=pi 上【不是常数】: S_W(x=0) = beta(1-(-1)/3) = 4beta/3 = %.6f"
      % (2.5 * 4 / 3))
    P("      (p05 的 'S_W ≡ beta' 只在 N=2: phi=(pi/2,-pi/2) 时 sum_j cos phi_j = 0)")
    P("    lam 依赖 (同一测度、同一 beta、同一 theta):")
    P("      lam=(1,0,-1)      E[H_vv] = %+.11f" % EhA)
    P("      lam=(1,-1,0)      E[H_vv] = %+.11f" % EhB)
    P("      lam=(0,1,-1)      E[H_vv] = %+.11f" % EhC)
    P("      lam=(1,1,-2)/v6   E[H_vv] = %+.11f" % EhD)
    P("    => E[H_vv] 在 v 上【变号】; 极小 = (beta/3)*min(2*c3, -2*c1, -2*c2) = %.11f"
      % ((2.5 / 3) * min(2 * c3, -2 * c1, -2 * c2)))
    P()
    P("    ★ 符号的唯一裁定依据 (不靠约定): theta=0 (phi=0) 在真空 x=0 上")
    P("      H_vv = (beta/N)||lam||^2 = (+)%s ; 真空是 S_W 的极小 => d_tau^2 S_W > 0 必成立"
      % ((2.5 / 3) * 2,))
    P("      故 H_vv = +(beta/N) sum lam^2 cos(phi-x) 是【规范】支; s02 的 H = -(...) 是 -H_vv")
    P("      (s02 自己的 docstring 就写了: 'r08/r12 用 H=+...; 与它们自己的 G 不自洽')")


def partF(rows):
    """[F] 在【无迹】Cartan 上扫 v: lam(psi) = (2/sqrt3)(cos psi, cos(psi-2pi/3), cos(psi+2pi/3)),
    sum lam = 0, ||lam||^2 = 2.  E[H_vv] 对 lam 线性 => 用 c_j 精确算, 两个端点用直接求积复核."""
    N, c1, c2, c3 = rows[-1][0], rows[-1][1], rows[-1][2], rows[-1][3]
    b = 2.5
    a1, a2, a3 = -c1, -c2, c3          # E[H_vv] = (b/3) sum_j a_j lam_j^2
    P("=" * 78)
    P("[F] 无迹 Cartan 上的 v-扫描 (sum lam=0, ||lam||^2=2),  theta=pi, beta=2.5")
    P("    a_j = (-c1,-c2,c3) = (%+.10f, %+.10f, %+.10f)   E[H_vv] = (b/3) sum a_j lam_j^2"
      % (a1, a2, a3))
    best = None
    lo = (1e9, None)
    hi = (-1e9, None)
    zeros = []
    M = 3600
    prev = None
    for i in range(M + 1):
        psi = math.pi * i / M           # 0..pi 覆盖 mod Weyl 的整族
        lam = tuple(2 / math.sqrt(3) * math.cos(psi - 2 * math.pi * j / 3) for j in range(3))
        # 注意: 与 (1,0,-1) 的对应关系另测, 这里只取族
        E = (b / 3) * sum(a * l * l for a, l in zip((a1, a2, a3), lam))
        if E < lo[0]:
            lo = (E, psi, lam)
        if E > hi[0]:
            hi = (E, psi, lam)
        if prev is not None and prev[0] * E < 0:
            zeros.append((prev[1], psi))
        prev = (E, psi)
    P("    min = %+.10f  at psi=%.6f  lam=(%+.6f,%+.6f,%+.6f)" % (lo[0], lo[1], *lo[2]))
    P("    max = %+.10f  at psi=%.6f  lam=(%+.6f,%+.6f,%+.6f)" % (hi[0], hi[1], *hi[2]))
    P("    变号区间 (psi) 个数 = %d ; 若非空 => E[H_vv] 在 v 上【变号】" % len(zeros))
    # 直接求积复核两个端点 + 说明 lambda=(1,0,-1) 在族里的位置
    def EH_direct(lam):
        acc, NN = cell_beta(2 * math.pi / 4096)
        Z = acc[0]
        c = (acc[1] / Z, acc[2] / Z, acc[3] / Z)
        return (b / 3) * ((-c[0]) * lam[0] ** 2 + (-c[1]) * lam[1] ** 2 + (c[2]) * lam[2] ** 2)
    P("    复核 (直接 2D 求积 N=4096, 用 c_j 线性式):")
    for lam, nm in ((lo[2], "min-lam"), (hi[2], "max-lam"),
                    ((1.0, 0.0, -1.0), "lam=(1,0,-1) [v04]"),
                    ((2 / math.sqrt(3), -1 / math.sqrt(3), -1 / math.sqrt(3)), "lam=(2,-1,-1)/v3")):
        P("      %-22s E[H_vv] = %+.10f" % (nm, EH_direct(lam)))
    P()


def partG(rows):
    """[G] p05 型【条件测度】(切片 <lam,X>=0) 在 SU(3) 2D Weyl cell 上是否良定义."""
    b = 2.5
    P("=" * 78)
    P("[G] p05 型条件测度 (<lam,X>=0 的切片) 在 SU(3) 2D Weyl cell 上的命运")
    P("    <lam,X> 的切片:  lam=(1,0,-1) -> x1-x3 = a+b ; (1,-1,0) -> a ; (0,1,-1) -> b")
    for lam, tag, line in (((1, 0, -1), "lam=(1,0,-1)", "a+b=0"),
                           ((1, -1, 0), "lam=(1,-1,0)", "a=0"),
                           ((0, 1, -1), "lam=(0,1,-1)", "b=0")):
        mx = 0.0
        for i in range(2001):
            t = 2 * math.pi * i / 2000
            d12 = t if tag == "lam=(1,-1,0)" else (t if tag == "lam=(1,0,-1)" else 0.0)
            d23 = -t if tag == "lam=(1,0,-1)" else (0.0 if tag == "lam=(1,-1,0)" else t)
            d13 = d12 + d23
            s = lambda d: 4 * math.sin(d / 2) ** 2
            mx = max(mx, s(d12) * s(d23) * s(d13))
        P("      %-14s 切片 %-6s :  max|Delta|^2 = %.3e  => %s"
          % (tag, line, mx, "权重恒 0, 条件期望 0/0 【无定义】" if mx < 1e-12 else "良定义"))
    P("    => p05 的条件化在 2D Weyl frame 里对【简单根方向】退化到奇异轨迹上")
    P("       (SU(2) 里 a=0 是整个 2-球面, 非退化 —— 这正是 r12 判 (c) 的 '维度' 缺口)")
    P()
    # 非退化的切片: lam=(2,-1,-1)/sqrt3 => x1=0  => b=-2a ; 该线上 S_W 恒为常数
    P("    非退化例子 lam=(2,-1,-1)/v3: 切片 x1=0 (b=-2a), 线上 x2=-a, x3=a")
    P("       S_W = beta(1-(1/3)(-cos0 - cos(-a) + cos a)) = beta(1+1/3) = 4beta/3 = %.6f  【常数】"
      % (4 * b / 3))
    P("       => p05 的 'S_W ≡ const 于切片, 但 H_vv ≠ 0' 结构在 SU(3) 里确实存在 (只是常数是 4beta/3 不是 beta)")
    num_, den_ = 0.0, 0.0
    for i in range(200001):
        t = 2 * math.pi * (i + 0.5) / 200000
        a = t
        s = lambda d: 4 * math.sin(d / 2) ** 2
        w = s(a) * s(-2 * a) * s(-a)          # |Delta|^2, exp 因子常数 -> 约掉
        H = b * (4 + 2 * math.cos(a)) / 9     # 该线上 H_vv (规范支), lam=(2,-1,-1)/v3
        num_ += w * H
        den_ += w
    P("       E[H_vv | x1=0] = %+.10f  (>0, 切片上逐点 H_vv = beta(4+2cos a)/9 >= 2beta/9)"
      % (num_ / den_))
    P()


def main():
    P("w01_ehsign.py  --  E[H_vv] 在 theta=pi 上的符号: p05 [T] vs v04 [C]")
    P("定义见文件头 [S1]-[S5].   beta=2.5, theta=pi, SU(2)/SU(3).")
    P()
    partA()
    partB()
    rows = partCD()
    partD()
    partE(rows)
    partF(rows)
    partG(rows)


if __name__ == "__main__":
    main()
