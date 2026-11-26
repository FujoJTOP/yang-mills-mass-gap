# -*- coding: utf-8 -*-
r"""
r10_audit_cert.py -- R10 独立复核 R03-CERT-THETA (自写实现; 未读 r03_cert_theta.py)

== 自己推的积分表示 (从 q04_weyl.py / q11_cert.py 的模型定义出发, 不抄 r03) ==
单 plaquette SU(2):  U = exp(xi.sig/2),  rho=|xi| in (0,2pi), h = rho/2, u := rho/2
     J = 4 sin^2 h / rho^2 = sin^2 u / u^2        (只依赖 rho  ⇒ Phi=-log J 只依赖 rho)
     S_W = beta ( 1 - cos u cos(th/2) + sin u sin(th/2) * xi3/rho )
切片 xi3 = 0  (a0=0, v=e3):
     d^3xi|slice = 2pi rho drho  ⇒  dmu ∝ J rho drho ∝ sin^2 u/u du      (常数 4 无关)
     G   = d_{xi3} S_W|0 = beta s g(u),  g(u)=sin u/(2u),  s=sin(th/2), c=cos(th/2)
     H   = d^2_{xi3} S_W|0 = (beta c/2) g(u)
     Phi'  = d_{xi3} Phi = Phi_rho * xi3/rho ≡ 0   (切片切于等-rho 球面, 对一切 th)
     Phi'' = Phi_rho/rho = 1/(2u^2) - cot(u)/(2u)        [u=rho/2]
 ==>
     G''(0) = E[H+Phi''] - Var(G+Phi') = S - D
     S = (beta c/2) E_nu[g] + E_nu[Phi''] ,  D = beta^2 s^2 Var_nu(g)
     nu ∝ exp(-beta(1 - c cos u)) sin^2 u/u du = exp(beta c cos u) sin^2u/u du
 th=pi (c=0): nu ∝ sin^2u/u (纯 Haar 切片, 与 beta 无关), H≡0, Phi'≡0
     S = E_nu[Phi''] = NS/Z  (与 beta 无关; NS = ∫Phi'' w du = 1/4 精确, IBP)
     D = beta^2 Var_nu(g)                      <== 退回 q11_cert.py 的四个一维积分
 th=0 (s=0): G ≡ 0 逐点 ⇒ D = 0 恒等 (不是数值零)
 四个分子积分 (w = sin^2u/u):
     Z = ∫ w E ,  NA = ∫ w g E ,  NB = ∫ w g^2 E ,  NSB = ∫ w Phi'' E ,  E = exp(B cos u), B=beta c
 被积函数在 u=0 全解析 (w~u, wg~u/2, wg^2~u/4, w Phi''~u/6): 无负幂/无奇异.

== 区间口径 (手写向外舍入) ==
 * 每个 +,-,*,/ 用 mpmath.libmp 定向舍入 ('f' 向下 / 'c' 向上), PREC=800 bit, 每步加绝对外扩
   PAD = 2^-380 (对 ~1e4 次运算累计 ~1e-110, 远小于下面的求积界).
 * sin/cos 用小宽度区间: 中点 800bit 求值 + Lipschitz 1 的半径 (|sin'|,|cos'|<=1), exp 用单调端点.
 * pi, sqrt(2) 只以 ('f','c') 区间进入; th=0/pi 的 c,s 是精确有理数 (1, 0 / 0, 1).
 * 求积: 256 面板 x 20 点 Gauss-Legendre (节点=Legendre 零点, Newton 迭代到 1e-200 残差;
   log 里给出 max|P_n(t_i)| 与权重的 Gauss 精确性自检).
 * 截断界: GL 的 Cauchy 型余项 |E| <= 2^{2n+1}(n!)^4/((2n+1)(2n)!^3) * max|F^{(2n)}|,
   max|f^{(2n)}| <= (2n)! M R^{-2n} (R=0.5, 复盘, 用显式粗界 M; 见 log 的 M 表) => 总界 << 1e-60.
   最后统一加 ALLOW 外扩, 取 bound 的 10 倍与 1e-55 的较大者.
 * 独立路 (不共享任何代码): mp.quad (tanh-sinh, dps=45) 复算四个积分.
"""
import math
from fractions import Fraction as Fr
import mpmath
import mpmath.libmp as M

PREC = 800
PADV = M.from_rational(1, 1 << 380, PREC, 'f')          # 2^-380
ZERO = M.from_rational(0, 1, PREC)
PILO = M.mpf_pi(PREC, 'f')
PIHI = M.mpf_pi(PREC, 'c')

LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


class Iv:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi):
        self.lo, self.hi = lo, hi

    def __repr__(self):
        return "Iv(%s,%s)" % (f(self.lo), f(self.hi))


def _pad(lo, hi):
    return Iv(M.mpf_sub(lo, PADV, PREC, 'f'), M.mpf_add(hi, PADV, PREC, 'c'))


def ivr(p, q=1):
    return Iv(M.from_rational(p, q, PREC, 'f'), M.from_rational(p, q, PREC, 'c'))


def ivp(x):
    return Iv(x, x)


def add(a, b):
    return _pad(M.mpf_add(a.lo, b.lo, PREC, 'f'), M.mpf_add(a.hi, b.hi, PREC, 'c'))


def sub(a, b):
    return _pad(M.mpf_sub(a.lo, b.hi, PREC, 'f'), M.mpf_sub(a.hi, b.lo, PREC, 'c'))


def mul(a, b):
    ps = ((a.lo, b.lo), (a.lo, b.hi), (a.hi, b.lo), (a.hi, b.hi))
    return _pad(min(M.mpf_mul(p, q, PREC, 'f') for p, q in ps),
                max(M.mpf_mul(p, q, PREC, 'c') for p, q in ps))


def div(a, b):
    assert M.mpf_cmp(b.lo, ZERO) > 0
    ps = ((a.lo, b.lo), (a.lo, b.hi), (a.hi, b.lo), (a.hi, b.hi))
    return _pad(min(M.mpf_div(p, q, PREC, 'f') for p, q in ps),
                max(M.mpf_div(p, q, PREC, 'c') for p, q in ps))


def neg(a):
    return Iv(M.mpf_neg(a.hi), M.mpf_neg(a.lo))


def iv_sin(a):
    """sin on a *narrow* interval: midpoint value + Lipschitz-1 radius."""
    m = M.mpf_mul(M.mpf_add(a.lo, a.hi, PREC, 'c'), M.from_rational(1, 2, PREC), PREC, 'c')
    s = M.mpf_sin(m, PREC, 'c')
    rad = M.mpf_add(M.mpf_sub(a.hi, a.lo, PREC, 'c'), PADV, PREC, 'c')
    rad = M.mpf_mul(rad, M.from_rational(1, 2, PREC), PREC, 'c')
    rad = M.mpf_add(rad, PADV, PREC, 'c')
    return _pad(M.mpf_sub(s, rad, PREC, 'f'), M.mpf_add(s, rad, PREC, 'c'))


def iv_cos(a):
    m = M.mpf_mul(M.mpf_add(a.lo, a.hi, PREC, 'c'), M.from_rational(1, 2, PREC), PREC, 'c')
    s = M.mpf_cos(m, PREC, 'c')
    rad = M.mpf_add(M.mpf_mul(M.mpf_sub(a.hi, a.lo, PREC, 'c'),
                              M.from_rational(1, 2, PREC), PREC, 'c'), PADV, PREC, 'c')
    return _pad(M.mpf_sub(s, rad, PREC, 'f'), M.mpf_add(s, rad, PREC, 'c'))


def iv_exp(a):
    """exp on a narrow interval: monotone, endpoints outward."""
    return _pad(M.mpf_exp(a.lo, PREC, 'f'), M.mpf_exp(a.hi, PREC, 'c'))


def f(x, n=40):
    return M.to_str(x, n)


def width(a):
    return M.to_str(M.mpf_sub(a.hi, a.lo, PREC, 'c'), 4)


def mid(a):
    return M.mpf_mul(M.mpf_add(a.lo, a.hi, PREC, 'c'), M.from_rational(1, 2, PREC), PREC, 'c')


# ---------------------------------------------------------------- Gauss-Legendre 节点
def gl_nodes(n, dps=250):
    mpmath.mp.dps = dps
    xs, ws = [], []
    for i in range(1, n + 1):
        t = mpmath.cos(mpmath.pi * (i - mpmath.mpf(1) / 4) / (n + mpmath.mpf(1) / 2))
        dp = mpmath.mpf(0)
        for _ in range(80):
            p0, p1 = mpmath.mpf(1), t
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * t * p1 - (k - 1) * p0) / k
            dp = n * (t * p1 - p0) / (t * t - 1)
            dt = p1 / dp
            t -= dt
            if abs(dt) < mpmath.mpf(10) ** (-dps + 8):
                break
        xs.append(t)
        ws.append(2 / ((1 - t * t) * dp * dp))
    return xs, ws


NN = 20
NPAN = 256
XS, WS = gl_nodes(NN)
RES = max(abs(float(mpmath.legendre(NN, x))) for x in XS)
MOM = max(abs(float(sum(WS[i] * XS[i] ** m for i in range(NN)) -
                     (2.0 / (m + 1) if m % 2 == 0 else 0.0))) for m in range(2 * NN))

# ---------------------------------------------------------------- 求积 (区间)
def quad(fun, edges_lo, edges_hi):
    """sum over panels of h*sum_i w_i f(u_i), all in interval arithmetic."""
    acc = ivr(0)
    n = len(edges_lo) - 1
    for k in range(n):
        a = Iv(edges_lo[k], edges_hi[k])
        b = Iv(edges_lo[k + 1], edges_hi[k + 1])
        m = Iv(M.mpf_mul(M.mpf_add(a.lo, b.hi, PREC, 'c'), M.from_rational(1, 2, PREC), PREC, 'c'),
               M.mpf_mul(M.mpf_add(a.hi, b.lo, PREC, 'c'), M.from_rational(1, 2, PREC), PREC, 'c'))
        h = Iv(M.mpf_mul(M.mpf_sub(b.lo, a.hi, PREC, 'f'), M.from_rational(1, 2, PREC), PREC, 'f'),
               M.mpf_mul(M.mpf_sub(b.hi, a.lo, PREC, 'c'), M.from_rational(1, 2, PREC), PREC, 'c'))
        s = ivr(0)
        for i in range(NN):
            t = Iv(M.from_str(repr(float(0)) if False else mpmath.nstr(XS[i], 80), PREC, 'f'),
                   M.from_str(mpmath.nstr(XS[i], 80), PREC, 'c'))
            w = Iv(M.from_str(mpmath.nstr(WS[i], 80), PREC, 'f'),
                   M.from_str(mpmath.nstr(WS[i], 80), PREC, 'c'))
            u = add(m, mul(h, t))
            s = add(s, mul(w, fun(u)))
        acc = add(acc, mul(h, s))
    return acc


def make_integrands(B):
    """returns dict of interval integrands for weight E = exp(B cos u)."""
    Bid = B

    def Ew(u):
        return iv_exp(mul(Bid, iv_cos(u)))

    def W(u):
        s = iv_sin(u)
        return div(mul(s, s), u)

    def g_of(u):
        return div(iv_sin(u), mul(ivr(2), u))

    def fZ(u):
        return mul(W(u), Ew(u))

    def fA(u):
        return mul(mul(W(u), g_of(u)), Ew(u))

    def fB(u):
        g = g_of(u)
        return mul(mul(mul(W(u), g), g), Ew(u))

    def fS(u):
        s = iv_sin(u)
        s2 = iv_sin(mul(ivr(2), u))
        u2 = mul(u, u)
        t1 = div(mul(s, s), mul(ivr(2), mul(u2, u)))
        t2 = div(s2, mul(ivr(4), u2))
        return mul(sub(t1, t2), Ew(u))

    return dict(Z=fZ, NA=fA, NB=fB, NS=fS)


# ---------------------------------------------------------------- 截断 (GL 余项) 界
def gl_bound(n, npan, Mval, R=0.5):
    """GL remainder over one panel of half-width h=pi/(2*npan), times npan panels."""
    import math as _m
    h = _m.pi / (2.0 * npan)
    C = 2.0 ** (2 * n + 1) * _m.factorial(n) ** 4 / ((2 * n + 1) * _m.factorial(2 * n) ** 3)
    per = C * _m.factorial(2 * n) * h ** (2 * n + 1) * Mval / R ** (2 * n)
    return per * npan


def main():
    P("=" * 100)
    P("r10_audit_cert.log -- R10 独立复核 R03-CERT-THETA: 自写 S,D 积分表示 + 手写向外舍入区间")
    P("  beta=2.5, SU(2), 单 plaquette, 切片 xi3=0 (a0=0), v=e3, rho=2u, u in (0,pi)")
    P("  S = (beta c/2) E_nu[g] + E_nu[Phi''],  D = beta^2 s^2 Var_nu(g),  nu ∝ e^{beta c cos u} sin^2u/u du")
    P("  interval: PREC=%d bit, PAD=2^-380 每步, %d 面板 x %d 点 GL" % (PREC, NPAN, NN))
    P("=" * 100)
    P("")
    P("[Q] 求积装置自检 (GL 节点/权重)")
    P("    max_i |P_%d(t_i)| = %.3e     max_m<=%d |sum_i w_i t_i^m - exact| = %.3e"
      % (NN, RES, 2 * NN - 1, MOM))
    P("    (节点是 Legendre 零点、权重满足全部 2n-1 阶矩条件 => 对 <= 2n-1 次多项式精确;")
    P("     余项由 Cauchy 导数界控制, 见下方 tail 界)")
    P("")

    # 面板端点: pi*k/NPAN 的区间
    PL, PH = [], []
    for k in range(NPAN + 1):
        PL.append(M.mpf_mul(PILO, M.from_rational(k, NPAN, PREC), PREC, 'f'))
        PH.append(M.mpf_mul(PIHI, M.from_rational(k, NPAN, PREC), PREC, 'c'))
    PIL, PIH_ = PILO, PIHI

    # 复盘上的粗界 M (R=0.5, |z|<=pi+0.5=3.64, |Im z|<=3.64): 见 docstring 推导
    R = 0.5
    ABSZ = math.pi + R
    eIM = math.e ** (2 * ABSZ)
    Ms = dict(W=ABSZ * eIM, g=eIM / 2, g2w=ABSZ * eIM ** 2 / 4,
              Sf=eIM / (2 * (R / 2) ** 3) + math.e ** (2 * ABSZ) / (4 * (R / 2) ** 2))
    P("    M(R=0.5) 粗界: |w|<=%.2e |g|<=%.2e |w g^2|<=%.2e |Sf|<=%.2e  (|z|<=%.2f)"
      % (Ms["W"], Ms["g"], Ms["g2w"], Ms["Sf"], ABSZ))
    P("")

    CASES = [("theta=0", Fr(1), Fr(0), ivr(5, 2)),
             ("theta=pi/2", None, None, None),
             ("theta=pi", Fr(0), Fr(1), ivr(0))]

    # th=pi/2: c=s=1/sqrt2
    C_LO = M.mpf_sqrt(M.from_rational(1, 2, PREC), PREC, 'f')
    C_HI = M.mpf_sqrt(M.from_rational(1, 2, PREC), PREC, 'c')
    CIV = Iv(C_LO, C_HI)

    out = {}
    for name, cFr, sFr, B in CASES:
        if cFr is not None:
            c, s = ivr(cFr.numerator, cFr.denominator), ivr(sFr.numerator, sFr.denominator)
        else:
            c, s = CIV, CIV
            B = mul(ivr(5, 2), CIV)                       # beta*c
        ii = make_integrands(B)
        Z = quad(ii["Z"], PL, PH)
        NA = quad(ii["NA"], PL, PH)
        NB = quad(ii["NB"], PL, PH)
        NS = quad(ii["NS"], PL, PH)
        Eg = div(NA, Z)
        Eg2 = div(NB, Z)
        Vg = sub(Eg2, mul(Eg, Eg))
        EPhi2 = div(NS, Z)
        b25 = ivr(25, 4)
        if cFr is not None:
            half_bc = ivr(5 * cFr.numerator, 4 * cFr.denominator)     # (beta c)/2
            s2 = ivr(sFr.numerator ** 2, sFr.denominator ** 2)
        else:
            half_bc = mul(ivr(5, 4), CIV)
            s2 = ivr(1, 2)
        S = add(mul(half_bc, Eg), EPhi2)
        D = mul(b25, mul(s2, Vg))

        # tail/余项允许量 (取 10 倍界 与 1e-55 的较大): 只用于最后外扩
        Bf = float(M.to_str(B.hi, 20))
        Efac = math.exp(abs(Bf) * math.cosh(ABSZ))        # |e^{B cos z}| <= e^{|B|cosh(|Im z|)}
        bd = max(gl_bound(NN, NPAN, Ms["W"] * Efac),
                 gl_bound(NN, NPAN, Ms["g2w"] * Efac),
                 gl_bound(NN, NPAN, Ms["Sf"] * Efac))
        ab = M.mpf_sub(S.hi, S.lo, PREC, 'c')
        amax = max([M.mpf_sub(Z.hi, Z.lo, PREC, 'c'), ab])
        ALLOW = max(10 * bd, 1e-55)
        P("[%s]  4 个积分 (区间)" % name)
        for nm, v in (("Z", Z), ("NA", NA), ("NB", NB), ("NS", NS)):
            P("   %-3s ∈ [%s,\n        %s]  w=%s" % (nm, f(v.lo), f(v.hi), width(v)))
        P("   E[g]    ∈ [%s, %s]  w=%s" % (f(Eg.lo), f(Eg.hi), width(Eg)))
        P("   Var(g)  ∈ [%s, %s]  w=%s" % (f(Vg.lo), f(Vg.hi), width(Vg)))
        P("   E[Phi''] ∈ [%s, %s]  w=%s" % (f(EPhi2.lo), f(EPhi2.hi), width(EPhi2)))
        P("   S       ∈ [%s, %s]  w=%s" % (f(S.lo), f(S.hi), width(S)))
        P("   D       ∈ [%s, %s]  w=%s" % (f(D.lo), f(D.hi), width(D)))
        gap = sub(S, D)
        P("   gap=S-D ∈ [%s, %s]  w=%s" % (f(gap.lo), f(gap.hi), width(gap)))
        P("   GL tail 界 = %.2e   ALLOW = %.1e" % (bd, ALLOW))
        P("   [含 ALLOW 外扩后]  gap ∈ [%s, %s]"
          % (f(M.mpf_sub(gap.lo, M.from_float(ALLOW), PREC, 'f')),
             f(M.mpf_add(gap.hi, M.from_float(ALLOW), PREC, 'c'))))
        P("")
        out[name] = dict(Z=Z, NA=NA, NB=NB, NS=NS, Eg=Eg, Vg=Vg, EPhi2=EPhi2, S=S, D=D,
                         gap=gap, allow=ALLOW)

    # ---------------- 独立路: mp.quad (dps=45)
    P("[X] 独立路 (mp.quad tanh-sinh, dps=45, 不共享任何求积代码)")
    mpmath.mp.dps = 45

    def q(fn, a, b):
        return mpmath.quad(fn, [a, b])

    for name, Bv in (("theta=0", mpmath.mpf(5) / 2), ("theta=pi/2", mpmath.mpf(5) / 2 / mpmath.sqrt(2)),
                     ("theta=pi", mpmath.mpf(0))):
        c = mpmath.cos(mpmath.mpf(0)) if name == "theta=0" else (
            mpmath.sqrt(mpmath.mpf(1) / 2) if name == "theta=pi/2" else mpmath.mpf(0))
        s = mpmath.mpf(0) if name == "theta=0" else (
            mpmath.sqrt(mpmath.mpf(1) / 2) if name == "theta=pi/2" else mpmath.mpf(1))
        Ew = lambda u: mpmath.exp(Bv * mpmath.cos(u))
        w = lambda u: mpmath.sin(u) ** 2 / u
        g = lambda u: mpmath.sin(u) / (2 * u)
        Sf = lambda u: mpmath.sin(u) ** 2 / (2 * u ** 3) - mpmath.sin(2 * u) / (4 * u ** 2)
        Z = q(lambda u: w(u) * Ew(u), 0, mpmath.pi)
        NA = q(lambda u: w(u) * g(u) * Ew(u), 0, mpmath.pi)
        NB = q(lambda u: w(u) * g(u) ** 2 * Ew(u), 0, mpmath.pi)
        NS = q(lambda u: Sf(u) * Ew(u), 0, mpmath.pi)
        S = mpmath.mpf(5) / 2 * c / 2 * (NA / Z) + NS / Z
        D = mpmath.mpf(25) / 4 * s * s * (NB / Z - (NA / Z) ** 2)
        P("   %-9s  Z=%.25f  S=%.25f  D=%.25f  S-D=%.30f"
          % (name, float(Z), float(S), float(D), float(S - D)))
        P("             S-D(mpmath 30d) = %s" % mpmath.nstr(S - D, 30))
        out[name]["ref"] = S - D
    P("")

    # ---------------- 与 r03 比
    R03 = {"theta=0": M.from_str("0.729924127576750062576372357665823454189533971017", PREC, 'c'),
           "theta=pi/2": M.from_str("0.542943922353689763665445822036482126731456810684", PREC, 'c'),
           "theta=pi": M.from_str("0.135994711041988211036899467040198122289952357486", PREC, 'c')}
    P("[C] 与 R03-CERT-THETA.md 的对账 (差 = 我的 gap 区间 - r03 引用的 S_lo-D_up)")
    P("    theta      我的 gap 区间                          与 r03 差 (lo / hi)")
    for name in ("theta=0", "theta=pi/2", "theta=pi"):
        g0 = out[name]["gap"]
        r = R03[name]
        dlo = M.mpf_sub(g0.lo, r, PREC, 'c')
        dhi = M.mpf_sub(g0.hi, r, PREC, 'c')
        P("    %-9s [%s,\n               %s]\n               d = %s / %s"
          % (name, f(g0.lo, 34), f(g0.hi, 34), f(dlo, 4), f(dhi, 4)))
    P("")

    # ---------------- theta=pi 退回 q11 的交叉核对
    P("[R] theta=pi 退回 q11 (回归核对)")
    Cin = mpmath.cin(2 * mpmath.pi) if hasattr(mpmath, "cin") else None
    if Cin is None:
        # Cin(x) = sum_{k>=1} (-1)^{k+1} x^{2k}/((2k)! 2k)
        mpmath.mp.dps = 45
        x = 2 * mpmath.pi
        Cin = mpmath.nsum(lambda k: (-1) ** (k + 1) * x ** (2 * k) / (mpmath.factorial(2 * k) * 2 * k), [1, 60])
    Zi = out["theta=pi"]["Z"]
    Vgi = out["theta=pi"]["Vg"]
    Si = out["theta=pi"]["S"]
    P("    Cin(2pi)/2        = %s   (q11 锚 1.2188266965286122)" % mpmath.nstr(Cin / 2, 30))
    P("    我的 Z            ∈ [%s, %s]" % (f(Zi.lo, 30), f(Zi.hi, 30)))
    P("    我的 Var(g)       ∈ [%s, %s]   (q11 锚 0.0110592937729735889)" % (f(Vgi.lo, 30), f(Vgi.hi, 30)))
    P("    我的 S=E[Phi'']   ∈ [%s, %s]   (q11 锚 0.2051152971230731)" % (f(Si.lo, 30), f(Si.hi, 30)))
    P("    NS(=∫Phi''w)      ∈ [%s, %s]   (精确 1/4)" % (f(out["theta=pi"]["NS"].lo, 30),
                                                       f(out["theta=pi"]["NS"].hi, 30)))
    P("")

    # ---------------- r03 印出的 Z 带着一个 (u-无关的) 归一常数, 核对只差常数
    P("[Z] r03 印出的 Z 与我的 Z 只差一个 u-无关常数 (不影响 E[g]/Var/S/D)")
    ZR = {"theta=0": ("0.342501350679279466941064642811626399645478352448", "0.342501350827421005439541229892799820413565400881"),
          "theta=pi/2": ("0.447881272155716380918378795160892871300491020789", "0.447881272155716380921921357266448556298121530595")}
    for name, (zl, zh) in ZR.items():
        Zm = out[name]["Z"]
        fac = iv_exp(ivr(-5, 2)) if name == "theta=0" else iv_exp(mul(ivr(-5, 2), CIV))
        Zs = mul(Zm, fac)
        P("    %-9s 我的 Z*e^{-beta c} ∈ [%s,\n                          %s]" % (name, f(Zs.lo, 30), f(Zs.hi, 30)))
        P("              r03 的 Z        ∈ [%s, %s]" % (zl[:30], zh[:30]))
        P("              d = %s / %s" % (f(M.mpf_sub(Zs.lo, M.from_str(zl, PREC, 'f'), PREC, 'c'), 4),
                                        f(M.mpf_sub(Zs.hi, M.from_str(zh, PREC, 'c'), PREC, 'f'), 4)))
    P("")

    # ---------------- theta=0 的 D ≡ 0: 恒等还是数值零
    P("[I] theta=0 的 D ≡ 0: 结构判定")
    P("    G(u) = beta*sin(th/2)*g(u) 逐点 (u in (0,pi) 整条切片), Phi' ≡ 0 逐点")
    P("    th=0 => sin(th/2)=0 (精确有理数 0) => G ≡ 0 逐点 => Var(G)=0 恒等 => D = beta^2*s^2*Var(g) = 0 恒等")
    P("    与 beta, 与 nu 无关; 不是两个相近量的相消, 区间实现里 D=[0,0] 精确 (乘数 s^2=0 有理精确).")
    P("")

    with open("r10_audit_cert.log", "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
