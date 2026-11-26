# -*- coding: utf-8 -*-
r"""
s02_su3_2dcert.py -- S02: SU(3) 2D Weyl gap 的区间算术证书 (r12 框架)

[A] g 的严格界:  g(d) = -2 log[sin(d/2)/(d/2)] 的 Taylor 级数
      log(sin x/x) = -sum_{n>=1} zeta(2n)/(n pi^{2n}) x^{2n}      (x=d/2)
    zeta(2n) in [1, zeta(2)=pi^2/6]  =>  系数全正, 上下界 + 显式尾部:
      tail_{n>N} <= (pi^2/3) * 1/(N+1) * (d/2pi)^{2N+2} / (1-(d/2pi)^2)
[B] 关键恒等式 (消掉 log/cot 奇性): 记 s(d)=4sin^2(d/2),
      s*g'   = 4(1-cos d)/d - 2 sin d            =: A(d)   整函数
      s*g''  = 2 - 4(1-cos d)/d^2                =: B(d)   整函数
      s*g'^2 = A(d)^2 / s(d)                     =: C(d)
    => 允许 1/d 的奇性只出现在 A,B,C 的定义式里, 被 s 精确吸收; 不需要
       log/cot 的逐点界 (r12 说的障碍在 |Delta|^2 权重下不存在)。
[C] 区间: mpmath.libmp.libmpi (mpi_add/mul/div/sin/cos/exp) + 200 bit
    + 节点 2pi/N 本身取区间 => 中点和的真外包围。
[D] 量 (与 r12 同框架):
      w=|Delta|^2 e^{-S_W},  S=<H+Phi''>,  D=Var(G+Phi')
      G=(beta/3)sum lam_j sin(phi_j-x_j)   =>  H=d_tau G = -(beta/3)sum lam_j^2 cos(phi_j-x_j)
      (r08/r12 用 H=+...; 与它们自己的 G 不自洽 —— 本脚本两个分支都算)
"""
import math
import numpy as np
import mpmath as mp
import mpmath.libmp as M
import mpmath.libmp.libmpi as LI

mp.mp.dps = 40
PREC = 200
PAD = M.mpf_pow(M.from_int(2, PREC), M.from_int(-(PREC // 2), PREC), PREC, 'c')  # 2^-100 外扩
LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


# ===================== 区间引擎 (libmpi, 真向外舍入) =====================
def iv(x):                                   # raw mpf -> 点区间
    return (x, x)


def frac(p, q):                              # 有理数 -> 区间
    return (M.from_rational(p, q, PREC, 'f'), M.from_rational(p, q, PREC, 'c'))


def mkint(n):
    return (M.from_int(n, PREC), M.from_int(n, PREC))


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


def isin(a):
    return LI.mpi_sin(a, PREC)


def icos(a):
    return LI.mpi_cos(a, PREC)


def iexp(a):
    return LI.mpi_exp(a, PREC)


def iabs(a):
    return LI.mpi_abs(a, PREC)


def ipad(a):                                 # 再外扩 (双保险)
    lo, hi = a
    w = M.mpf_sub(hi, lo, PREC, 'c')
    m = M.mpf_add(M.mpf_add(M.mpf_abs(lo, PREC), M.mpf_abs(hi, PREC), PREC, 'c'),
                  w, PREC, 'c')
    p = M.mpf_mul(m, PAD, PREC, 'c')
    return (M.mpf_sub(lo, p, PREC, 'f'), M.mpf_add(hi, p, PREC, 'c'))


def sng(a):
    return M.mpf_cmp(a[0], M.fzero)


def num(a, nd=20):
    return M.to_str(a[0], nd), M.to_str(a[1], nd)


def wd(a):
    return M.mpf_sub(a[1], a[0], PREC, 'c')


# ===================== [A] g 的严格界 (zeta(2n) 级数) =====================
ZETA2_LO = frac(1, 1)                        # zeta(2n) >= 1
ZETA2_UP = frac(1, 1)                        # 占位, 下面用 pi^2/6


def partA():
    P("[A] g(d) = -2 log[sin(d/2)/(d/2)] 的严格界")
    P("    级数: log(sin x/x) = -sum_{n>=1} c_n x^{2n},  c_n = zeta(2n)/(n pi^{2n}),  x=d/2")
    P("    系数正, zeta(2n) in [1, zeta(2)] => 部分和给下界; 尾部用显式界:")
    P("      R_N(x) <= (pi^2/3) * (1/(N+1)) * (x/pi)^{2N+2} / (1-(x/pi)^2)")
    pi = (M.mpf_pi(PREC, 'f'), M.mpf_pi(PREC, 'c'))
    pi2 = imul(pi, pi)
    zeta2 = idiv(pi2, mkint(6))              # zeta(2)=pi^2/6 (区间)
    # 精确有理系数 c_n = zeta(2n)/(n pi^{2n}) = 2^{2n-1}|B_{2n}|/(n (2n)!)  (Bernoulli, Fraction)
    from fractions import Fraction as Fr
    NB = 64
    BB = [Fr(1)]
    for m in range(1, NB + 1):
        t = Fr(0)
        for k in range(m):
            t += Fr(math.comb(m + 1, k)) * BB[k]
        BB.append(-t / Fr(m + 1))

    def ccoef(n):
        return Fr(2 ** (2 * n - 1) * abs(BB[2 * n]), n * math.factorial(2 * n))
    rows = []
    for ds in ("0.05", "0.5", "1.0", "2.0", "3.0", "4.0", "5.0", "6.0"):
        d = (M.from_str(ds, PREC, 'f'), M.from_str(ds, PREC, 'c'))
        x = idiv(d, mkint(2))
        for NN in (6, 12, 24):
            pw = iv(M.fone)                  # x^{2n}
            x2 = imul(x, x)
            lo = (M.fzero, M.fzero)
            for n in range(1, NN + 1):
                pw = imul(pw, x2)
                nn, dd_ = ccoef(n).numerator, ccoef(n).denominator
                lo = iadd(lo, imul(pw, frac(nn, dd_)))
            # 尾部 (用 c_n <= zeta2/(n pi^2n), 1/n <= 1/(N+1))
            r = idiv(x, pi)
            r2 = imul(r, r)
            # (x/pi)^{2N+2}
            q = iv(M.fone)
            for _ in range(2 * NN + 2):
                q = imul(q, r)
            tail = imul(imul(idiv(pi2, mkint(3)), frac(1, NN + 1)), q)   # 2*zeta(2)=pi^2/3
            one_m = isub(iv(M.fone), r2)
            if sng(one_m) > 0:
                tail = idiv(tail, one_m)
            glo = imul(mkint(2), lo)
            gup = iadd(glo, tail)
            # 参照值 (mpmath 40 dps -> 外扩 1e-30 相对)
            gm = mp.mpf(ds)
            gref = -2 * mp.log(mp.sin(gm / 2) / (gm / 2))
            raw = gref._mpf_
            sl = M.mpf_mul(raw, M.from_str('0.999999999999999999999999999999', PREC, 'f'),
                           PREC, 'f')
            sh = M.mpf_mul(raw, M.from_str('1.000000000000000000000000000001', PREC, 'c'),
                           PREC, 'c')
            g_iv = (sl, sh)
            rows.append((ds, NN, glo, gup, g_iv, tail))
    P("    %-6s %-4s %-24s %-24s %-26s %-10s" % ("d", "N", "g_lo", "g_up", "width", "tail"))
    for ds, NN, glo, gup, gr, tail in rows:
        P("    %-6s %-4d %-24s %-24s %-26s %-10s"
          % (ds, NN, num(glo, 12)[0], num(gup, 12)[1], M.to_str(wd(isub(gup, glo)), 4),
             M.to_str(wd((M.fzero, tail[1])), 3)))
        if NN == 24:
            ok = (M.mpf_cmp(glo[0], gr[1]) <= 0) and (M.mpf_cmp(gup[1], gr[0]) >= 0)
            P("           contain mpmath ref? %s" % ok)
    P()

    P("    [B] 奇性吸收恒等式 (60 dps 机检 s*g', s*g'', s*g'^2):")
    mp.mp.dps = 60
    okb = True
    for ds in ("0.3", "1.7", "4.1", "5.9"):
        d = mp.mpf(ds)
        g1 = 2 / d - mp.cot(d / 2)
        g2 = 1 / (2 * mp.sin(d / 2) ** 2) - 2 / d ** 2
        s = 4 * mp.sin(d / 2) ** 2
        A1 = 4 * (1 - mp.cos(d)) / d - 2 * mp.sin(d)
        B1 = 2 - 4 * (1 - mp.cos(d)) / d ** 2
        e1 = abs(s * g1 - A1)
        e2 = abs(s * g2 - B1)
        e3 = abs(s * g1 ** 2 - A1 ** 2 / s)
        okb = okb and e1 < mp.mpf('1e-50') and e2 < mp.mpf('1e-50') and e3 < mp.mpf('1e-50')
        P("      d=%-4s |s g' - A|=%.1e  |s g''-B|=%.1e  |s g'^2-A^2/s|=%.1e" %
          (ds, float(e1), float(e2), float(e3)))
    P("      => A,B,C 是整函数 (无 1/d 奇性): log/cot 的逐点严格界 **不需要**. [T] ok=%s" % okb)
    mp.mp.dps = 40
    P()


# ===================== [C] 区间中点求和 =====================
BETA = frac(5, 2)                            # 2.5
LAM = (frac(1, 1), frac(0, 1), frac(-1, 1))
PHI = ((M.mpf_pi(PREC, 'f'), M.mpf_pi(PREC, 'c')),
       (M.mpf_pi(PREC, 'f'), M.mpf_pi(PREC, 'c')),
       frac(0, 1))
DD_PAIRS = (('12', '23', '13'), ('23', '12', '13'), ('13', '12', '23'))  # (p, other1, other2)
PAIRS = ('12', '23', '13')
DDF = {'12': isub(LAM[0], LAM[1]), '23': isub(LAM[1], LAM[2]), '13': isub(LAM[0], LAM[2])}


def cell(a, b):
    d = {'12': a, '23': b, '13': iadd(a, b)}
    c = {k: icos(v) for k, v in d.items()}
    sn = {k: isin(v) for k, v in d.items()}
    s = {k: imul(mkint(2), isub(iv(M.fone), c[k])) for k in PAIRS}      # 4sin^2(d/2)=2-2cos d
    Ad = {k: isub(idiv(imul(mkint(4), isub(iv(M.fone), c[k])), d[k]),
                  imul(mkint(2), sn[k])) for k in PAIRS}
    Bd = {k: isub(mkint(2), idiv(imul(mkint(4), isub(iv(M.fone), c[k])), imul(d[k], d[k])))
          for k in PAIRS}
    Cd = {}
    for k in PAIRS:
        half = idiv(d[k], mkint(2))
        # C = s g'^2 = (4 sin(d/2)/d - 2 cos(d/2))^2  -- 整函数 (无 1/s, 无 cot 极点)
        Ee = isub(idiv(imul(mkint(4), isin(half)), d[k]), imul(mkint(2), icos(half)))
        Cd[k] = imul(Ee, Ee)
    T = imul(imul(s['12'], s['23']), s['13'])
    x = (idiv(iadd(imul(mkint(2), a), b), mkint(3)),
         idiv(isub(b, a), mkint(3)),
         ineg(idiv(iadd(a, imul(mkint(2), b)), mkint(3))))
    G = (M.fzero, M.fzero)
    H = (M.fzero, M.fzero)
    csum = (M.fzero, M.fzero)
    for j in range(3):
        u = isub(PHI[j], x[j])
        cu, su = icos(u), isin(u)
        csum = iadd(csum, cu)
        G = iadd(G, imul(imul(BETA, LAM[j]), su))
        H = isub(H, imul(imul(BETA, imul(LAM[j], LAM[j])), cu))
    G = idiv(G, mkint(3))
    H = idiv(H, mkint(3))
    SW = imul(BETA, isub(iv(M.fone), idiv(csum, mkint(3))))
    e = iexp(ineg(SW))
    w = imul(T, e)
    Z1 = (M.fzero, M.fzero)      # sum dd_p A_p s_o1 s_o2
    Z2 = (M.fzero, M.fzero)
    Z11 = (M.fzero, M.fzero)
    for p, o1, o2 in DD_PAIRS:
        Z1 = iadd(Z1, imul(imul(DDF[p], Ad[p]), imul(s[o1], s[o2])))
        Z2 = iadd(Z2, imul(imul(imul(DDF[p], DDF[p]), Bd[p]), imul(s[o1], s[o2])))
        Z11 = iadd(Z11, imul(imul(imul(DDF[p], DDF[p]), Cd[p]), imul(s[o1], s[o2])))
    for pp, qq in (('12', '23'), ('12', '13'), ('23', '13')):
        rr = ({'12', '23', '13'} - {pp, qq}).pop()
        Z11 = iadd(Z11, imul(imul(imul(mkint(2), imul(DDF[pp], DDF[qq])),
                                 imul(Ad[pp], Ad[qq])), s[rr]))
    return (w, imul(w, H), imul(e, Z2), imul(e, Z1), imul(e, Z11),
            imul(w, G), imul(w, imul(G, G)), imul(imul(e, G), Z1))


def quad_int(N):
    pi = (M.mpf_pi(PREC, 'f'), M.mpf_pi(PREC, 'c'))
    h = idiv(imul(mkint(2), pi), mkint(N))
    acc = [(M.fzero, M.fzero)] * 8
    half = frac(1, 2)
    for i in range(N):
        ai = imul(iadd(mkint(i), half), h)
        for j in range(N):
            bj = imul(iadd(mkint(j), half), h)
            v = cell(ai, bj)
            acc = [iadd(acc[k], v[k]) for k in range(8)]
    return [ipad(a) for a in acc]


def quad_np(N):
    h = 2 * math.pi / N
    g = (np.arange(N) + 0.5) * h
    A, B = np.meshgrid(g, g, indexing='ij')
    d = {'12': A, '23': B, '13': A + B}
    c = {k: np.cos(v) for k, v in d.items()}
    s = {k: 2 - 2 * c[k] for k in PAIRS}
    Ad = {k: 4 * (1 - c[k]) / d[k] - 2 * np.sin(d[k]) for k in PAIRS}
    Bd = {k: 2 - 4 * (1 - c[k]) / d[k] ** 2 for k in PAIRS}
    Cd = {k: (4 * np.sin(v / 2) / v - 2 * np.cos(v / 2)) ** 2 for k, v in d.items()}
    T = s['12'] * s['23'] * s['13']
    x1 = (2 * A + B) / 3
    x2 = (B - A) / 3
    x3 = -(A + 2 * B) / 3
    lam = (1.0, 0.0, -1.0)
    phi = (math.pi, math.pi, 0.0)
    G = np.zeros_like(T)
    H = np.zeros_like(T)
    cs = np.zeros_like(T)
    for j, xj in enumerate((x1, x2, x3)):
        G += (2.5 / 3) * lam[j] * np.sin(phi[j] - xj)
        H -= (2.5 / 3) * lam[j] ** 2 * np.cos(phi[j] - xj)
        cs += np.cos(phi[j] - xj)
    SW = 2.5 * (1 - cs / 3)
    w = T * np.exp(-SW)
    e = np.exp(-SW)
    Z1 = np.zeros_like(T)
    Z2 = np.zeros_like(T)
    Z11 = np.zeros_like(T)
    dd = {'12': 1.0, '23': 1.0, '13': 2.0}
    oth = {'12': ('23', '13'), '23': ('12', '13'), '13': ('12', '23')}
    for p in PAIRS:
        o1, o2 = oth[p]
        Z1 += dd[p] * Ad[p] * s[o1] * s[o2]
        Z2 += dd[p] ** 2 * Bd[p] * s[o1] * s[o2]
        Z11 += dd[p] ** 2 * Cd[p] * s[o1] * s[o2]
    for pp, qq in (('12', '23'), ('12', '13'), ('23', '13')):
        rr = ({'12', '23', '13'} - {pp, qq}).pop()
        Z11 += 2 * dd[pp] * dd[qq] * Ad[pp] * Ad[qq] * s[rr]
    return dict(Z=w.sum(), ZH=(w * H).sum(), Z2=(e * Z2).sum(), Z1=(e * Z1).sum(),
                Z11=(e * Z11).sum(), ZG=(w * G).sum(), ZG2=(w * G * G).sum(),
                ZG1=(e * G * Z1).sum(), G2=(w * G * np.ones_like(T)).sum())


def sums_to_SD(Z, ZH, Z2, Z1, Z11, ZG, ZG2, ZG1, hsign):
    E = lambda v: idiv(v, Z)
    E2 = E(Z2)
    EH = E(ZH)
    S = iadd(E2, EH) if hsign > 0 else isub(E2, imul(mkint(-1), EH))
    # hsign=-1: H 已按 -(beta/3) 计算, 直接用; hsign=+1: 翻回 r08/r12 的符号
    if hsign > 0:
        S = iadd(E2, ineg(EH))
    else:
        S = iadd(E2, EH)
    D = iadd(iadd(E(Z11), E(ZG2)), imul(mkint(2), E(ZG1)))
    D = isub(D, imul(iadd(E(ZG), E(Z1)), iadd(E(ZG), E(Z1))))
    return S, D, E2, EH, E(Z1), E(ZG), E(ZG2), E(ZG1)


def partC():
    P("[C] 区间中点求和 (libmpi, PREC=200, 节点 2pi/N 也取区间; 再 pad 2^-100):")
    P("    %-6s %-6s %-22s %-22s" % ("N", "t(s)", "S_lo", "S_up"))
    import time
    res = {}
    for N in (64, 128):
        t0 = time.time()
        acc = quad_int(N)
        res[N] = acc
        S, D, E2, EH, E1, EG, EG2, EG1 = sums_to_SD(*acc, hsign=-1)
        S2, D2, *_ = sums_to_SD(*acc, hsign=+1)
        P("    %-6d %-6.1f %-22s %-22s" % (N, time.time() - t0, num(S, 12)[0], num(S, 12)[1]))
        P("           [H=-(b/3)..] S in [%s, %s]  D in [%s, %s]"
          % (num(S, 14)[0], num(S, 14)[1], num(D, 14)[0], num(D, 14)[1]))
        P("           E[Phi''] in [%s, %s]  E[H] in [%s, %s]"
          % (num(E2, 14)[0], num(E2, 14)[1], num(EH, 14)[0], num(EH, 14)[1]))
        P("           gap(dS-D) width = %s ;  [H=+..] S in [%s, %s]"
          % (M.to_str(wd(isub(S, D)), 4), num(S2, 12)[0], num(S2, 12)[1]))
        P("           Z interval width = %s" % M.to_str(wd(acc[0]), 3))
        P("           widths: " + " ".join(M.to_str(wd(v), 3) for v in acc))
        g_lo = isub(S, D)
        P("           S_lo - D_up = %s   SIGN=%s"
          % (M.to_str(M.mpf_sub(S[0], D[1], PREC, 'c'), 12),
             M.mpf_cmp(M.mpf_sub(S[0], D[1], PREC, 'c'), M.fzero)))
        g2 = isub(S2, D2)
        P("           [H=+..] S_lo - D_up = %s  SIGN=%s"
          % (M.to_str(M.mpf_sub(S2[0], D2[1], PREC, 'c'), 12),
             M.mpf_cmp(M.mpf_sub(S2[0], D2[1], PREC, 'c'), M.fzero)))
        P("    %.1f s total" % (time.time() - t0))
        P()
    return res


def partD():
    P("[D] 浮点收敛表 (同一被积函数, 显示 端点非周期性 => O(h^2)):")
    P("    %-6s %-13s %-13s %-13s %-13s" % ("N", "E[Phi'']", "E[H]", "S", "D"))
    prev = {}
    for N in (128, 256, 512, 1024, 2048):
        o = quad_np(N)
        Z = o['Z']
        E = lambda k: o[k] / Z
        E2, EH, E1, EG, EG2, EG1, E11 = (E('Z2'), E('ZH'), E('Z1'), E('ZG'),
                                         E('ZG2'), E('ZG1'), E('Z11'))
        S = E2 + EH
        D = E11 + EG2 + 2 * EG1 - (EG + E1) ** 2
        P("    %-6d %-13.9f %-13.9f %-13.9f %-13.9f" % (N, E2, EH, S, D))
        prev[N] = (S, D, E2, EH)
    P()
    P("    相邻差 (外推):")
    for N in (256, 512, 1024, 2048):
        m = N // 2
        P("      N=%-5d dS=%.3e  dD=%.3e  dE[Phi'']=%.3e  ratio S/D = %.2f/%.2f"
          % (N, prev[N][0] - prev[m][0], prev[N][1] - prev[m][1],
             prev[N][2] - prev[m][2],
             (prev[m][0] - prev[m // 2][0]) / (prev[N][0] - prev[m][0]) if N > 256 else float('nan'),
             (prev[m][1] - prev[m // 2][1]) / (prev[N][1] - prev[m][1]) if N > 256 else float('nan')))
    S, D = prev[2048][0], prev[2048][1]
    dS = S - prev[1024][0]
    dD = D - prev[1024][1]
    P("    Richardson (O(h^2), 用 N=1024/2048):")
    P("      S_inf = %.9f  (± %.1e)   D_inf = %.9f  (± %.1e)"
      % (S + dS / 3, abs(dS) / 3, D + dD / 3, abs(dD) / 3))
    EHInf = prev[2048][3] + (prev[2048][3] - prev[1024][3]) / 3
    E2Inf = prev[2048][2] + (prev[2048][2] - prev[1024][2]) / 3
    P("      gap_der (H=d_tau G) = %+.9f ;  gap_r12 (H=+(b/3)sum..) = %+.9f"
      % ((E2Inf + EHInf) - (D + dD / 3), (E2Inf - EHInf) - (D + dD / 3)))
    P()
    P("    所需 N (用 C h^2, C=|dS|*N^2/(3*4pi^2*...) 标定):")
    C = abs(dS) / 3 * (2048 / (2 * math.pi)) ** 2
    for tgt in (1e-3, 1e-4, 1e-6):
        P("      width %-8.0e -> N >= %.0f" % (tgt, (2 * math.pi) * math.sqrt(C / tgt)))
    P("      (上表 C=%.3e 由 N=1024/2048 的 dS/3 标定; 严格化需 Cauchy/Hessian 界)" % C)
    P()


if __name__ == "__main__":
    P("=" * 108)
    P("s02_su3_2dcert.log -- S02: SU(3) 2D Weyl gap 区间 (r12 框架)")
    P("=" * 108)
    P()
    partA()
    partC()
    partD()
    P("[E] 工具/精度: libmpi PREC=%d, pad 2^-100; 浮点对照 numpy float64" % PREC)
    with open("s02_su3_2dcert.log", "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
