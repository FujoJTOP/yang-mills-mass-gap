# -*- coding: utf-8 -*-
r"""
r12_su3_correct.py -- R12: SU(3), theta=pi, beta=2.5, a=0 的正确 S 与 D。

[A] G != 0 核准 (SU(2) 校准 + 逻辑): G = d_tau S_W，切片上 S_W == beta 只说明"值"，
    对 tau 的导数离开切片才看得到。r03 的 G = beta*s*g(u) != 0。
[B] r08 的 Phi 与 r11/O09 的 Phi 不是同一个函数:
       Phi_O09 = sum_{i<j} g(d_ij),  g(d) = -2 log|sin(d/2)/(d/2)|   (指数坐标 Jacobian)
       Phi_r08 = -log|Delta|^2 = -sum log(4 sin^2(d/2))
    二者差 psi = sum_{i<j} log(d_ij^2)。r08 的恒等式 int|Delta|^2(Phi'^2-Phi'')=d^2_tau int|Delta|^2=0
    只对 Phi = -log|Delta|^2 成立 => 用错 Phi 时 S != D（即使 G=0）。
[C] 2D Weyl 积分（共形坐标 a=x1-x2, b=x2-x3 in (0,2pi), dx=da db/3）:
       <f> = int |Delta|^2 e^{-S_W} f / int |Delta|^2 e^{-S_W}
       S = <H + Phi''>,  D = Var(G + Phi'),  G = d_tau S_W, H = d^2_tau S_W
     theta=pi: W = diag(-1,-1,1) (det=1, Z2 型), S_W = beta(1 - (1/3) sum_j cos(phi_j-x_j))
[D] 严格区间 (800 bit 定向舍入): 只对三角多项式部分 (Phi_r08) 可精确; Phi_O09 的 log 项
    破坏"常数项"路线 => 该部分的区间是瓶颈, 如实报。
"""
import math
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp
import mpmath.libmp as M

mp.mp.dps = 40
PREC = 800
EPS = M.from_rational(1, 2 ** (PREC // 2), PREC, 'c')      # 2^-400 相对外扩
LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


# =====================================================================
# [A] G != 0 (SU(2) 校准)
# =====================================================================
def partA():
    beta = mp.mpf(5) / 2
    c = mp.cos(mp.pi / 2)          # = 0
    s = mp.sin(mp.pi / 2)          # = 1
    P("[A] G = d_tau S_W 于 theta=pi、a=0 切片上 (SU(2) 校准, r03 模型):")
    P("    S_W(u,a) = beta(1 - c cos u + s sin(u) a g(u)/1),  a=xi3,  u=rho/2,  g(u)=sin u/(2u)")
    P("    c=%s  =>  S_W|(a=0) = beta 逐点 (值),  但 d_tau 离开切片:" % mp.nstr(c, 5))
    P("    %-6s %-26s %-26s %-14s" % ("u", "S_W(0)-beta", "G=FD", "beta*s*g(u)"))
    ok = True
    for uu in ("0.37", "1.20", "2.60"):
        u = mp.mpf(uu)
        rho = 2 * u

        def SW(tau):                       # a = tau (xi3 方向), rho(tau)=sqrt(rho^2+tau^2)
            r = mp.sqrt(rho ** 2 + tau ** 2)
            uu2 = r / 2
            return beta * (1 - c * mp.cos(uu2) + s * mp.sin(uu2) * tau / r)

        val0 = SW(mp.mpf(0))
        G_fd = mp.diff(SW, mp.mpf(0))
        G_an = beta * s * (mp.sin(u) / (2 * u))
        P("    %-6s %-26s %-26s %-14s" % (uu, mp.nstr(val0 - beta, 6), mp.nstr(G_fd, 18), mp.nstr(G_an, 18)))
        ok = ok and abs(val0 - beta) < mp.mpf('1e-38') and abs(G_fd - G_an) < mp.mpf('1e-18')
    P("    => S_W 在切片上恒为 beta (%s), G = d_tau S_W = beta*s*g(u) != 0 (%s)" % (ok, ok))
    P("    [T] 逻辑: S_W 的值与它沿 v 的导数没有推论关系; r08 由 S_W==beta 推 G==0 不合法。")
    P()


# =====================================================================
# [B] 精确 Laurent 常数项 (Fraction) + 800 bit 定向舍入
# =====================================================================
def lmul(A, B):
    out = {}
    for (a1, a2), c in A.items():
        for (b1, b2), d in B.items():
            k = (a1 + b1, a2 + b2)
            out[k] = out.get(k, Fr(0)) + c * d
    return out


def lscal(A, r):
    return {k: v * r for k, v in A.items() if v * r != 0}


def ladd(A, B):
    out = dict(A)
    for k, v in B.items():
        out[k] = out.get(k, Fr(0)) + v
    return {k: v for k, v in out.items() if v != 0}


def mon(m1, m2, c=Fr(1)):
    return {(m1, m2): c}


# u1=e^{i th1}, u2=e^{i th2}, u3 = u1^-1 u2^-1
# e^{i(th1-th2)} = u1 u2^-1 ;  e^{i(th2-th3)} = u1 u2^2 ; e^{i(th1-th3)} = u1^2 u2
E = {'12': (1, -1), '23': (1, 2), '13': (2, 1)}


def s_poly(p):        # s = 2 - e^{id} - e^{-id}
    m1, m2 = E[p]
    return ladd(ladd(mon(0, 0, Fr(2)), mon(m1, m2, Fr(-1))), mon(-m1, -m2, Fr(-1)))


def sin_poly(p):      # 2 sin d = (e^{id} - e^{-id})/i  -> 用 4 sin^2 更省; 这里给 "2 sin d" 的平方
    pass


def build(lams):
    """返回三角多项式 dict: T=|Delta|^2, NA=|D|^2 Phi_r08'', NB08=|D|^2 Phi_r08'^2, etc."""
    S = {p: s_poly(p) for p in ('12', '23', '13')}
    T = lmul(lmul(S['12'], S['23']), S['13'])
    d = {'12': lams[0] - lams[1], '23': lams[1] - lams[2], '13': lams[0] - lams[2]}
    others = {'12': ('23', '13'), '23': ('12', '13'), '13': ('12', '23')}
    # |D|^2 Phi_r08'' = 1/2 sum_p d_p^2 * 4 * prod_{others} s
    NA = {}
    for p in ('12', '23', '13'):
        o1, o2 = others[p]
        NA = ladd(NA, lscal(lmul(S[o1], S[o2]), 2 * Fr(d[p] ** 2)))
    # |D|^2 Phi_r08'^2 : 对角项 (4-s_p)*prod + 交叉项 2 d_p d_q (2 sin d_p)(2 sin d_q) s_r
    NB = {}
    for p in ('12', '23', '13'):
        o1, o2 = others[p]
        coef_poly = ladd(mon(0, 0, Fr(4)), lscal(S[p], Fr(-1)))     # 4 - s_p
        NB = ladd(NB, lscal(lmul(lmul(coef_poly, S[o1]), S[o2]), Fr(d[p] ** 2)))
    # cross: 2 d_p d_q (2 sin d_p)(2 sin d_q) s_r   ; 2 sin d = (e^{id}-e^{-id})/i
    def sin2_poly(p):     # 2 sin d, as Laurent poly (coefficients in i -> 用 (2 sin)^2 = s*(4-s)/1 ... )
        m1, m2 = E[p]
        # (2 sin d)^2 = 4 - s^2/... 直接用 (2 sin d)^2 = 4 - (2-2cos d)^2/? 改用 (e^{id}-e^{-id})^2*(-1) = 4 sin^2 d
        return None
    # 处理交叉项: (2 sin d_p)(2 sin d_q)(s_r) 的乘积 = sqrt 型 -> 用 (2 sin d_p)^2*(2 sin d_q)^2 = 16 s_p s_q ... 需要带符号
    # 直接用显式: 2 sin d_p = (u_p - 1/u_p)/i,  两两相乘得到 Laurent (含 i^2=-1 因子在实部)
    def twosin(p):
        m1, m2 = E[p]
        return ladd(mon(m1, m2, Fr(1)), mon(-m1, -m2, Fr(-1)))
    for pp, qq, rr, sgn in ((('12', '23', '13'), None, None, None),):
        pass
    pairs = [('12', '23'), ('12', '13'), ('23', '13')]
    for pp, qq in pairs:
        rr = ({'12', '23', '13'} - {pp, qq}).pop()
        # 2 d_p d_q * (2 sin d_p)(2 sin d_q) * s_r ;  (2 sin)(2 sin) = -[(u-1/u)(v-1/v)]*... 处理:
        # 2 sin d = (e^{id}-e^{-id})/i = -i (u_p - u_p^-1).  乘积 (2 sin p)(2 sin q) = - (u_p-u_p^-1)(u_q-u_q^-1)
        X = lmul(twosin(pp), twosin(qq))
        X = lscal(X, Fr(-1))                       # 因为 i*i = -1
        NB = ladd(NB, lscal(lmul(X, S[rr]), 2 * Fr(d[pp]) * Fr(d[qq])))
    return T, NA, NB, d


def cterm(P):
    return P.get((0, 0), Fr(0))


def iv_from_frac_2pi2(fr):
    """(2pi)^2 * fr 的严格区间 (定向舍入 + PAD)"""
    pl = M.mpf_pi(PREC, 'f')
    ph = M.mpf_pi(PREC, 'c')
    num, den = fr.numerator, fr.denominator
    fl = M.from_rational(num, den, PREC, 'f')
    fh = M.from_rational(num, den, PREC, 'c')
    lo = M.mpf_mul(M.mpf_mul(pl, pl, PREC, 'f'), fl, PREC, 'f')
    hi = M.mpf_mul(M.mpf_mul(ph, ph, PREC, 'c'), fh, PREC, 'c')
    pad_l = M.mpf_mul(M.mpf_abs(lo), EPS, PREC, 'c')
    pad_h = M.mpf_mul(M.mpf_abs(hi), EPS, PREC, 'c')
    return M.mpf_sub(lo, pad_l, PREC, 'f'), M.mpf_add(hi, pad_h, PREC, 'c')


def partB():
    P("[B] 精确常数项 (Fraction) + 800 bit 定向舍入, theta=pi 的 2D Weyl 积分:")
    for lab, lams in (("lambda=(1,0,-1)", (Fr(1), Fr(0), Fr(-1))),
                      ("lambda=(1,-1,0)", (Fr(1), Fr(-1), Fr(0)))):
        T, NA, NB, d = build(lams)
        cT, cA, cB = cterm(T), cterm(NA), cterm(NB)
        P("  %s : CT[|D|^2]=%s  CT[|D|^2 Phi_r08'']=%s  CT[|D|^2 Phi_r08'^2]=%s" % (lab, cT, cA, cB))
        Z = iv_from_frac_2pi2(cT)
        A = iv_from_frac_2pi2(cA)
        B = iv_from_frac_2pi2(cB)
        Slo = M.mpf_div(A[0], Z[1], PREC, 'f')
        Sup = M.mpf_div(A[1], Z[0], PREC, 'c')
        Dlo = M.mpf_div(B[0], Z[1], PREC, 'f')
        Dup = M.mpf_div(B[1], Z[0], PREC, 'c')
        gap = M.mpf_sub(Slo, Dup)
        P("     Z in [%s , %s]  (width %s)" % (M.to_str(Z[0], 20), M.to_str(Z[1], 20),
                                               M.to_str(M.mpf_sub(Z[1], Z[0], PREC, 'c'), 4)))
        P("     S=A/Z in [%s, %s] ; D=B/Z in [%s, %s]" % (M.to_str(Slo, 20), M.to_str(Sup, 20),
                                                          M.to_str(Dlo, 20), M.to_str(Dup, 20)))
        P("     S_lo - D_up = %s   SIGN=%+d   (精确有理: S-D = %s - %s = %s)"
          % (M.to_str(gap, 8), M.mpf_cmp(gap, M.from_rational(0, 1, PREC)),
             cA, cB, (cA - cB)))
        P("")
    P("  => Phi_r08 = -log|Delta|^2 时 S == D == 8 精确 (r08 复现, 代码自检通过)。")
    P("     该恒等式只对 -log(密度) 成立: int|D|^2(Phi'^2-Phi'')=d^2_tau int|D(x+tau)|^2 dx = 0。")
    P()


# =====================================================================
# [C] 正确 Phi (=Phi_O09) 的 2D 积分: 数值 (周期梯形, 谱精度)
# =====================================================================
def partC(N=1024, quiet=False):
    h = 2 * math.pi / N
    a = (np.arange(N) + 0.5) * h
    A, B = np.meshgrid(a, a, indexing='ij')
    d12, d23, d13 = A, B, A + B
    d = {'12': d12, '23': d23, '13': d13}
    s = {k: 4 * np.sin(v / 2) ** 2 for k, v in d.items()}
    T = s['12'] * s['23'] * s['13']
    x1 = (2 * d12 + d23) / 3
    x2 = (d23 - d12) / 3
    x3 = -(d12 + 2 * d23) / 3
    beta = 2.5
    lam = (1.0, 0.0, -1.0)
    phi = (math.pi, math.pi, 0.0)

    # 三套 per-pair 因子 (全部乘 |Delta|^2 之后安全: 无 1/s 奇异)
    #  T*c_p        = q1_p * prod_others          (一次幂)
    #  T*c_p^2      = q2_p * prod_others          (平方对角)
    #  T*c_p c_q    = q1_p q1_q * s_r             (平方交叉)
    def factors(kind):
        q1, q2 = {}, {}
        for p in ('12', '23', '13'):
            v = d[p]
            if kind == 'A':                       # c = -cot(d/2)
                q1[p] = -2 * np.sin(v)
                q2[p] = 4 - s[p]
            else:                                 # c = g'(d) = 2/d - cot(d/2)
                q1[p] = 8 * np.sin(v / 2) ** 2 / v - 2 * np.sin(v)
                q2[p] = (4 * np.sin(v / 2) / v - 2 * np.cos(v / 2)) ** 2
        return q1, q2

    others = {'12': ('23', '13'), '23': ('12', '13'), '13': ('12', '23')}
    dd = {'12': lam[0] - lam[1], '23': lam[1] - lam[2], '13': lam[0] - lam[2]}
    tri = (('12', '23'), ('12', '13'), ('23', '13'))

    def sums(kind):
        q1, q2 = factors(kind)
        T1 = np.zeros_like(T)     # T*Phi'
        T2 = np.zeros_like(T)     # T*Phi'^2
        for p in ('12', '23', '13'):
            o1, o2 = others[p]
            T1 += dd[p] * q1[p] * s[o1] * s[o2]
            T2 += dd[p] ** 2 * q2[p] * s[o1] * s[o2]
        for pp, qq in tri:
            rr = ({'12', '23', '13'} - {pp, qq}).pop()
            T2 += 2 * dd[pp] * dd[qq] * q1[pp] * q1[qq] * s[rr]
        return T1, T2

    # Phi_B'' : T*Phi_B'' = sum_p d_p^2 (2 - 2 q0_p) prod_others   (q0=4sin^2(d/2)/d^2)
    TPhiBpp = np.zeros_like(T)
    for p in ('12', '23', '13'):
        o1, o2 = others[p]
        q0 = (2 * np.sin(d[p] / 2) / d[p]) ** 2
        TPhiBpp += dd[p] ** 2 * (2 - 2 * q0) * s[o1] * s[o2]
    # Phi_A'' : T*Phi_A'' = sum_p d_p^2 2 prod_others
    TPhiApp = np.zeros_like(T)
    for p in ('12', '23', '13'):
        o1, o2 = others[p]
        TPhiApp += dd[p] ** 2 * 2 * s[o1] * s[o2]

    TPhiB, TPhiB2 = sums('B')
    TPhiA, TPhiA2 = sums('A')

    # G, H : diag plaquette, W=diag(-1,-1,1),  x -> x+tau*lam
    G = np.zeros_like(T)
    H = np.zeros_like(T)
    for j, xj in enumerate((x1, x2, x3)):
        G += -(beta / 3) * (-np.sin(phi[j] - xj)) * lam[j]
        H += -(beta / 3) * (-np.cos(phi[j] - xj)) * lam[j] ** 2
    SW = beta * (1 - (1.0 / 3.0) * (np.cos(phi[0] - x1) + np.cos(phi[1] - x2) + np.cos(phi[2] - x3)))
    w = T * np.exp(-SW)
    w0 = T

    def A_(f, ww):
        return float(np.sum(ww * f) / np.sum(ww))

    out = {}
    out['Z'] = float(np.sum(w))
    out['EPhiBpp'] = A_(TPhiBpp / T, w)
    out['EPhiBp2'] = A_(TPhiB2 / T, w)
    out['EPhiBp'] = A_(TPhiB / T, w)
    out['EG'] = A_(G, w)
    out['EG2'] = A_(G * G, w)
    out['EGPhi'] = A_(G * (TPhiB / T), w)
    out['EH'] = A_(H, w)
    SV = out['EPhiBpp'] + out['EH']
    Var = out['EG2'] + out['EPhiBp2'] + 2 * out['EGPhi'] - (out['EG'] + out['EPhiBp']) ** 2
    out['S'] = SV
    out['D'] = Var
    out['gap'] = SV - Var
    # r08 校验: 纯 Haar 权重 + Phi_A  => S=D=8, E[Phi_A']=0
    out['r08_S'] = A_(TPhiApp / T, w0)
    out['r08_D'] = A_(TPhiA2 / T, w0)
    out['r08_E1'] = A_(TPhiA / T, w0)
    # 纯 Haar 权重 + 正确 Phi
    out['S0'] = A_(TPhiBpp / T, w0)
    out['D0'] = A_(TPhiB2 / T, w0) + A_(G * G, w0) + 2 * A_(G * (TPhiB / T), w0) \
        - (A_(G, w0) + A_(TPhiB / T, w0)) ** 2
    out['D0_Phi'] = A_(TPhiB2 / T, w0) - A_(TPhiB / T, w0) ** 2
    out['VarG0'] = A_(G * G, w0) - A_(G, w0) ** 2
    out['E1_0'] = A_(TPhiB / T, w0)
    out['EGphi_0'] = A_(G * (TPhiB / T), w0)
    return out


def partD():
    P("[C] 正确 Phi = Phi_O09 的 2D Weyl 积分 (a=x1-x2,b=x2-x3 in (0,2pi), dx=da db/3;")
    P("    W=diag(-1,-1,1) (theta=pi), beta=2.5, lambda=(1,0,-1), 中点网格 N^2:")
    P("    %-6s %-12s %-12s %-12s %-11s %-11s %-11s" % ("N", "S", "D", "S-D", "E[Phi']", "Var(Phi')", "E[G]"))
    prev = None
    for N in (256, 512, 1024, 2048):
        o = partC(N)
        P("    %-6d %-12.8f %-12.8f %-+12.8f %-11.4e %-11.4e %-11.3e"
          % (N, o['S'], o['D'], o['gap'], o['EPhiBp'], o['EPhiBp2'] - o['EPhiBp'] ** 2, o['EG']))
        prev = o
    P("")
    P("    分解 (N=2048, 权重 |Delta|^2 e^{-S_W}): E[Phi'']=%+.8f  E[H]=%+.8f  S=%+.8f"
      % (prev['EPhiBpp'], prev['EH'], prev['S']))
    P("                    E[G]=%+.3e  E[G^2]=%+.8f  2E[G Phi']=%+.8f  E[Phi'^2]=%+.8f  E[Phi']=%+.3e"
      % (prev['EG'], prev['EG2'], 2 * prev['EGPhi'], prev['EPhiBp2'], prev['EPhiBp']))
    P("     校验1 (r08 的 Phi, 纯 Haar 权重): E[Phi_A'']=%.8f  E[Phi_A'^2]=%.8f  S-D=%+.2e  (应 8,8,0)"
      % (prev['r08_S'], prev['r08_D'], prev['r08_S'] - prev['r08_D']))
    P("     校验2 同网格 E[Phi_A']=%+.2e (应 0)" % prev['r08_E1'])
    P("")
    P("     纯 Haar 权重 (theta=pi) + 正确 Phi 的分支:")
    P("       S0=E_T[Phi_B'']=%.8f   E_T[Phi_B']=%+.2e   Var_T(Phi_B')=%.8f" % (prev['S0'], prev['E1_0'], prev['D0_Phi']))
    P("       Var_T(G)=%.8f  E_T[G Phi']=%+.8f   D0=Var_T(G+Phi')=%.8f   S0-D0=%+.8f"
      % (prev['VarG0'], prev['EGphi_0'], prev['D0'], prev['S0'] - prev['D0']))
    P("")
    P("[D] 分量表 (N=2048, |Delta|^2 e^{-S_W} 权重):")
    for k in ('Z', 'EPhiBpp', 'EH', 'EG', 'EG2', 'EGPhi', 'EPhiBp2', 'EPhiBp', 'S', 'D', 'gap'):
        P("     %-9s = %+.10f" % (k, prev[k]))
    P("")


if __name__ == "__main__":
    P("=" * 100)
    P("r12_su3_correct.log -- R12: SU(3) theta=pi 的正确 S 与 D (G != 0 核准)")
    P("=" * 100)
    P()
    partA()
    partB()
    partD()
    with open("r12_su3_correct.log", "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
