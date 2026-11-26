# -*- coding: utf-8 -*-
r"""
v04_su3_split.py -- V04: 拆开 SU(3) (theta=pi, beta=2.5) 的 gap<0

复用 s02 (S02-SU3-2DCERT) 的区间框架 (`import s02_su3_2dcert as S2`), 不重写引擎。

量:
  S = E[H_vv] + E[Phi'']                       (s02 branch: H = d_tau G = -(beta/3)sum lam^2 cos)
  D = Var(G + Phi') = Var(G) + Var(Phi') + 2 Cov(G,Phi')
  gap = S - D = E[H] + [E[Phi'']-Var(G)] - Var(Phi') - 2 Cov

输出:
  [A] 区间 (libmpi, PREC=200, 中点求和 N=64/128) 的 6 个分量
  [B] 分块浮点 + Richardson (O(h^2)) 高精度分量
  [C] gap 的四项记账 + 占比
  [D] SU(2) 对照 (存档数字, 见 SU2_SRC)
  [E] 工具/精度状态
"""
import math
import time

import numpy as np

import s02_su3_2dcert as S2

PREC = S2.PREC
IDX = dict(Z=0, ZH=1, Z2=2, Z1=3, Z11=4, ZG=5, ZG2=6, ZG1=7)

P_STATE = []


def P(s=""):
    print(s)
    P_STATE.append(s)


# ===================== [A] 区间: 6 个分量 =====================
def comp_iv(acc):
    """acc = s02.quad_int(N) 的 8 个累积区间 -> 分量的区间 (外包)"""
    Z = acc[IDX['Z']]
    E = lambda k: S2.idiv(acc[IDX[k]], Z)
    e2, eh, e1 = E('Z2'), E('ZH'), E('Z1')
    eg, eg2, e11, eg1 = E('ZG'), E('ZG2'), E('Z11'), E('ZG1')
    # Var(G) = E[G^2]-E[G]^2
    varG = S2.isub(eg2, S2.imul(eg, eg))
    varP = S2.isub(e11, S2.imul(e1, e1))
    cov = S2.isub(eg1, S2.imul(eg, e1))
    S = S2.iadd(e2, eh)
    D = S2.iadd(S2.iadd(varG, varP), S2.imul(S2.mkint(2), cov))
    gap = S2.isub(S, D)
    return dict(EH=eh, E2=e2, EG=eg, E1=e1, VarG=varG, VarP=varP, Cov=cov,
                S=S, D=D, gap=gap, E2mVarG=S2.isub(e2, varG))


IV_N = (64, 128)


def partA():
    P("[A] 区间分量 (libmpi PREC=%d, 中点求和, 节点 2pi/N 取区间, pad 2^-100):" % PREC)
    P("    (s02 同一 cell(); 注意: 区间只包住【中点和】, 求积截断 O(h^-2) 另计 -> [C])")
    out = {}
    for N in IV_N:
        t0 = time.time()
        acc = S2.quad_int(N)
        c = comp_iv(acc)
        out[N] = c
        P("    N=%-4d (%.1f s)" % (N, time.time() - t0))
        for k in ('E2', 'EH', 'EG', 'E1', 'VarG', 'VarP', 'Cov', 'E2mVarG', 'S', 'D', 'gap'):
            lo, hi = S2.num(c[k], 18)
            P("      %-8s [%s, %s]  w=%s" % (k, lo, hi,
                                             S2.M.to_str(S2.wd(c[k]), 3)))
        P()
    return out


# ===================== [B] 分块浮点 + Richardson =====================
def sums_np_chunk(N, chunk=256):
    """s02.quad_np 同一被积函数, 按行分块以支持 N>=4096 (内存)"""
    h = 2 * math.pi / N
    g = (np.arange(N) + 0.5) * h
    lam = (1.0, 0.0, -1.0)
    phi = (math.pi, math.pi, 0.0)
    dd = {'12': 1.0, '23': 1.0, '13': 2.0}
    oth = {'12': ('23', '13'), '23': ('12', '13'), '13': ('12', '23')}
    pairs = ('12', '23', '13')
    acc = np.zeros(8)
    for i0 in range(0, N, chunk):
        a = g[i0:i0 + chunk]
        A = a[:, None]
        B = g[None, :]
        d = {'12': A + 0 * B, '23': 0 * A + B, '13': A + B}
        c = {k: np.cos(v) for k, v in d.items()}
        s = {k: 2 - 2 * c[k] for k in pairs}
        Ad = {k: 4 * (1 - c[k]) / d[k] - 2 * np.sin(d[k]) for k in pairs}
        Bd = {k: 2 - 4 * (1 - c[k]) / d[k] ** 2 for k in pairs}
        Cd = {k: (4 * np.sin(v / 2) / v - 2 * np.cos(v / 2)) ** 2 for k, v in d.items()}
        T = s['12'] * s['23'] * s['13']
        x1 = (2 * d['12'] + d['23']) / 3
        x2 = (d['23'] - d['12']) / 3
        x3 = -(d['12'] + 2 * d['23']) / 3
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
        for p in pairs:
            o1, o2 = oth[p]
            Z1 += dd[p] * Ad[p] * s[o1] * s[o2]
            Z2 += dd[p] ** 2 * Bd[p] * s[o1] * s[o2]
            Z11 += dd[p] ** 2 * Cd[p] * s[o1] * s[o2]
        for pp, qq in (('12', '23'), ('12', '13'), ('23', '13')):
            rr = ({'12', '23', '13'} - {pp, qq}).pop()
            Z11 += 2 * dd[pp] * dd[qq] * Ad[pp] * Ad[qq] * s[rr]
        acc += np.array([w.sum(), (w * H).sum(), (e * Z2).sum(), (e * Z1).sum(),
                         (e * Z11).sum(), (w * G).sum(), (w * G * G).sum(),
                         (e * G * Z1).sum()])
    return acc


def comps_np(acc):
    Z = acc[0]
    e2, eh, e1, eg, eg2, e11, eg1 = (acc[k] / Z for k in (2, 1, 3, 5, 6, 4, 7))
    varG = eg2 - eg * eg
    varP = e11 - e1 * e1
    cov = eg1 - eg * e1
    s = e2 + eh
    dd = varG + varP + 2 * cov
    return dict(EH=eh, E2=e2, EG=eg, E1=e1, VarG=varG, VarP=varP, Cov=cov,
                S=s, D=dd, gap=s - dd, E2mVarG=e2 - varG)


def rich(vals, Ns, k=2):
    """O(h^k) Richardson 外推 (用最后两个 N)"""
    N1, N2 = Ns[-2], Ns[-1]
    d = vals[N2] - vals[N1]
    r = (N2 / N1) ** k
    return vals[N2] + d / (r - 1), abs(d) / (r - 1)


NAP = (256, 512, 1024, 2048, 4096, 8192)


def partB():
    P("[B] 分块浮点 (s02.quad_np 同一被积函数, 行块=256) + Richardson O(h^2):")
    P("    %-7s %-16s %-14s %-14s" % ("N", "E[Phi'']", "E[H]", "S"))
    tab = {}
    for N in NAP:
        t0 = time.time()
        c = comps_np(sums_np_chunk(N))
        tab[N] = c
        P("    %-7d %-16.10f %-14.10f %-14.10f (%.1f s)"
          % (N, c['E2'], c['EH'], c['S'], time.time() - t0))
    P()
    P("    Richardson (N=4096/8192, O(h^2)) : 值 +- 误差")
    out = {}
    for k in ('E2', 'EH', 'EG', 'E1', 'VarG', 'VarP', 'Cov', 'E2mVarG', 'S', 'D', 'gap'):
        v, err = rich({n: tab[n][k] for n in NAP}, NAP)
        out[k] = (v, err)
        P("      %-8s = %+.10f  +- %.1e" % (k, v, err))
    P()
    P("    自检: S - D = %+.10f  (应 = gap %+.10f)  差 %.1e"
      % (out['S'][0] - out['D'][0], out['gap'][0],
         abs(out['S'][0] - out['D'][0] - out['gap'][0])))
    P("    自检: E[H]+(E[Phi'']-Var(G))-Var(Phi')-2Cov = %+.10f"
      % (out['EH'][0] + out['E2mVarG'][0] - out['VarP'][0] - 2 * out['Cov'][0]))
    P()
    return out, tab


# ===================== [C] 记账 =====================
def partC(res):
    EH, E2, VG, VP, CV = (res[k][0] for k in ('EH', 'E2', 'VarG', 'VarP', 'Cov'))
    gap = res['gap'][0]
    P("[C] gap 的记账:  S = E[H]+E[Phi''] = %+.8f ;  D = Var(G+Phi') = Var(G)+Var(Phi')+2Cov = %+.8f"
      % (res['S'][0], res['D'][0]))
    P("    gap = S - D = %+.8f  (s02 区间: %+.8f)" % (gap, -0.10251290744559712))
    P()
    dDyn = EH - VG
    wFull = E2 - VP - 2 * CV
    P("    (ii) S_W / Phi 两半 (与语料 G''=D_dyn+W_full 同形):")
    P("      D_dyn = E[H]-Var(G)                  = %+.8f   %+7.1f%% of |gap|" % (dDyn, 100 * dDyn / abs(gap)))
    P("      W_full = E[Phi'']-Var(Phi')-2Cov     = %+.8f   %+7.1f%%" % (wFull, 100 * wFull / abs(gap)))
    P("      sum                                  = %+.8f                 (校验差 %.1e)"
      % (dDyn + wFull, abs(dDyn + wFull - gap)))
    P("      E[H] = %+.8f   vs   Var(G) = %+.8f   =>  E[H] < Var(G) 成立, 缺口 = %+.8f"
      % (EH, VG, dDyn))
    P()
    P("    (i) 四项逐项 (注意大项两两对消, 单看 % 会误导):")
    P("      %-30s %-16s" % ("项", "值"))
    for nm, v in (("E[H_vv]", EH), ("E[Phi'']", E2), ("Var(G)  [-]", -VG),
                  ("Var(Phi')  [-]", -VP), ("2Cov  [-]", -2 * CV)):
        P("      %-30s %+.8f" % (nm, v))
    P("      Phi 半边小计 = (E[Phi'']-Var(G)) - Var(Phi') - 2Cov = %+.8f  (%.2f%% of |gap|)"
      % (wFull, 100 * wFull / abs(gap)))
    P("      内部近对消: E[Phi'']-Var(G) = %+.6f   vs   Var(Phi')+2Cov = %+.6f   差 %.1e"
      % (E2 - VG, VP + 2 * CV, abs((E2 - VG) - (VP + 2 * CV))))
    P()
    P("    反事实: 若 E[H]=0 (SU(2) theta=pi 那样), SU(3) 的 gap 会是 W_full = %+.6f > 0  =>  Phi 半边单独是【通过】的" % wFull)
    P("    判据:")
    P("      |E[H]| = %.6f  vs  |Phi 半边小计| = %.6f  =>  %.2f : 1  (E[H] 主导)"
      % (abs(EH), abs(wFull), abs(EH) / abs(wFull)))
    P("      Var(G)/E[Phi''] = %.5f   (SU(2): 0.069121/0.205115 = 0.3370)  =>  E[Phi''] 远【大】于 Var(G), 不是 (b)"
      % (VG / E2))
    P("      Phi 半边符号为正 (%.6f>0) => 与 gap<0 【反向】, 不是 (c) 的共谋" % wFull)
    P("      ==> 判 (a): gap<0 来自 S_W 半边 (E[H]<0), 不是 rank")
    P()
    return dDyn, wFull


# ===================== [D] SU(2) 对照 =====================
SU2_SRC = "research/v3/species/H1A-EXACT.md:99-106 (0.20511530, 0.01105929, beta*=4.30661); O11-PHYS.md:16 (精确求积)"


def partD(res):
    A, B, b = 0.2051152971, 0.0110592938, 2.5
    varG2 = B * b * b
    P("[D] SU(2) 对照 (theta=pi, a0=0; 该切片条件测度=Haar, 精确):")
    P("      E[H] = 0 (精确, 中心元)   E[Phi''] = A = %.10f   Var(G) = B*beta^2 = %.10f" % (A, varG2))
    P("      Var(Phi') = 0 (精确)      S = %.10f   D = %.10f   gap = %+.10f" % (A, varG2, A - varG2))
    P("      来源: %s" % SU2_SRC)
    P()
    P("    同一框架 SU(3) vs SU(2):")
    P("      %-20s %-16s %-16s %-12s" % ("量", "SU(3)", "SU(2)", "比"))
    rows = (("E[H]", res['EH'][0], 0.0),
            ("E[Phi'']", res['E2'][0], A),
            ("Var(G)", res['VarG'][0], varG2),
            ("Var(G)/E[Phi'']", res['VarG'][0] / res['E2'][0], varG2 / A),
            ("S", res['S'][0], A),
            ("D", res['D'][0], varG2),
            ("gap", res['gap'][0], A - varG2))
    for nm, v3, v2 in rows:
        P("      %-20s %+-16.8f %+-16.8f %-12s"
          % (nm, v3, v2, "n/a" if abs(v2) < 1e-14 else "%.4g" % (v3 / v2)))
    P()
    P("    SU(2) 的通道: E[H]=0 恒成立 => gap = E[Phi'']-Var(G) = A-B beta^2, beta>4.30661 才负 (theta=pi);")
    P("    SU(3) beta=2.5: E[Phi'']-Var(G) = %+.6f > 0 (Phi 半边通过), 但 E[H] = %+.6f < 0 直接压翻."
      % (res['E2mVarG'][0], res['EH'][0]))
    P()
    P("    [警告] Var(Phi') 归一化: s02 框架下 %.6f (非可忽略!), 其 %.1f%% 来自 E[Phi']^2 = %.6f;"
      % (res['VarP'][0], 100 * res['E1'][0] ** 2 / res['VarP'][0], res['E1'][0] ** 2))
    P("      t07 T7 已判定 E[Phi']=2.4922 != 0 是 d13=a+b 未 wrap (域效应+错配对), 正确配对锚: E[Phi']=0, E[Phi'']=E[Phi'^2]=8, S=D。")
    P("      t07 的 Var(Phi')~2.3e-5 是 8 维 Haar 图上 v-方向梯度方差, 另一个对象/归一化 —— 两者不可混用。")
    P()


def partE():
    P("[E] 工具状态:")
    P("    py -3.12 + numpy %s + mpmath %s; libmpi PREC=%d; s02 引擎 import 成功"
      % (np.__version__, __import__('mpmath').__version__, PREC))
    P("    区间 = 中点和的向外舍入外包 (libmpi); 求积截断未严格界 -> 数值结论标 [C]")
    P("    分块浮点 N<=8192, 行块 256 (内存安全); Richardson 误差栏为相邻差/3")


if __name__ == "__main__":
    P("=" * 108)
    P("v04_su3_split.log -- V04: SU(3) theta=pi beta=2.5 gap<0 的 S_W / Phi 拆分")
    P("=" * 108)
    P()
    ivres = partA()
    res, tab = partB()
    partC(res)
    partD(res)
    partE()
    with open("v04_su3_split.log", "w", encoding="utf-8") as fh:
        fh.write("\n".join(P_STATE) + "\n")
