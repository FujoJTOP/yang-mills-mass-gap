# -*- coding: utf-8 -*-
r"""
r03_cert_theta.py -- R03: 把 q11 的 theta=pi 严格区间证书扩到 theta != pi

== 单 plaquette 精确约化 (与 q11_cert.py / q04 同一模型) ==
    S_W = beta(1 - cos(u) c + sin(u) s xi3/rho),  c=cos(th/2), s=sin(th/2), u=rho/2
    J   = 4 sin^2(u)/rho^2 = sin^2(u)/u^2      (只依赖 rho)
切片 xi3 = 0:
    测度  dmu ∝ e^{-beta(1-c cos u)} sin^2(u)/u du            [q04: e^{-beta(1-c cos(rho/2))} sin^2(rho/2) rho^{-1} drho]
    G   = d_{xi3} S_W = beta s g(u),        g(u)=sin u/(2u)
    H   = d^2_{xi3} S_W = (beta c/2) g(u)                     [ theta!=pi 时 H != 0 ]
    Phi = -log J 只依赖 rho  =>  Phi' = d_{xi3} Phi = Phi_rho' (xi3/rho) ≡ 0  (切于等-rho 球面)
    Phi''= Phi_rho'/rho = 1/(2u^2) - cot(u)/(2u)
故 (对一切 theta, 不只是 pi):
    S = E[H + Phi''] = (beta c/2) E_nu[g] + E_nu[Phi'']
    D = Var(G + Phi') = beta^2 s^2 Var_nu[g]         <== 无交叉项, Phi'≡0
    nu ∝ e^{-beta(1-c cos u)} sin^2(u)/u du ,  u in (0,pi)

== 严格性 ==
把 e^{-beta(1-c cos u)} =: Ebar(u) 作为整函数级数 (Ebar(0)=1, 偶函数), 与
    W=sin^2u/u, A=sin^3u/(2u^2), B=sin^4u/(4u^3), Sf=Phi''W
的整函数级数相乘, 逐项积分 int_0^pi u^k du = pi^{k+1}/(k+1).
 1) 系数在 Q[sqrt2] 里精确 (theta=0: Q; theta=pi/2: c=sqrt2/2 => 系数 a+b sqrt2, a,b in Q).
 2) Ebar 系数由 ODE  Ebar' = -beta c sin(u) Ebar, Ebar(0)=1  精确递推 (无浮点).
 3) pi, sqrt2 只以 libmp 定向舍入区间进入 (prec=800, 'f'/'c').
 4) 每个 +,-,*,/ 定向舍入 + 绝对外扩 PAD=2^-400.
 5) tail: Cauchy  |a_k| <= M(r)/r^k  (r=4.5, |u|=r),
    M(r) <= e^{beta|c|(cosh r -1)} * sup|f0|,  sup|f0| 用 sinh 的严格上界 (e^4.5<=91, e^9<=8104);
    tail <= M(r) (pi/r)^{K+1} / ((K+2)(1-pi/r)) ,  pi<3.2 => pi/r < 32/45.
"""
from fractions import Fraction as Fr
import math
import mpmath.libmp as M

PREC = 800
EPS = M.from_rational(1, 2 ** (PREC // 2), PREC, 'c')     # 2^-400 (相对外扩, 见 _w)
KMAX = 400
NK = KMAX // 2
LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


# ---------------------------------------------------------------- Q[sqrt2] 精确元 a + b*sqrt2
class Alg:
    __slots__ = ("a", "b")

    def __init__(self, a, b=Fr(0)):
        self.a = Fr(a)
        self.b = Fr(b)


def Aa(x, y):
    return Alg(x.a + y.a, x.b + y.b)


def As(x, y):
    return Alg(x.a - y.a, x.b - y.b)


def Am(x, y):
    return Alg(x.a * y.a + 2 * x.b * y.b, x.a * y.b + x.b * y.a)


def Ac(x, r):
    return Alg(x.a * r, x.b * r)


# ---------------------------------------------------------------- 区间算术 (libmp 定向舍入)
class Iv:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi):
        self.lo, self.hi = lo, hi


def iv_rat(p, q=1):
    return Iv(M.from_rational(p, q, PREC, 'f'), M.from_rational(p, q, PREC, 'c'))


def _w(lo, hi):
    """相对外扩: pad = max(|lo|,|hi|)*2^-400  (纯相对, 无绝对项).
    theta!=pi 的级数有灾难性相消(峰值项 ~1e7, 系数小到 1e-227), 绝对外扩(如 2^-400=3.9e-121)
    会把 1e-227 的系数整条吃掉, 再乘 pi^{k+1}(~1e199) 放大成 1e78 的假误差 => 必须纯相对.
    输入的区间宽度(来自 sqrt2/pi 的定向舍入)由区间运算法则自行传播, 相消时自动正确。"""
    mag = M.mpf_abs(lo)
    m2 = M.mpf_abs(hi)
    if M.mpf_cmp(m2, mag) > 0:
        mag = m2
    pad = M.mpf_mul(mag, EPS, PREC, 'c')
    return Iv(M.mpf_sub(lo, pad, PREC, 'f'), M.mpf_add(hi, pad, PREC, 'c'))


def iv_add(a, b):
    return _w(M.mpf_add(a.lo, b.lo, PREC, 'f'), M.mpf_add(a.hi, b.hi, PREC, 'c'))


def iv_sub(a, b):
    return _w(M.mpf_sub(a.lo, b.hi, PREC, 'f'), M.mpf_sub(a.hi, b.lo, PREC, 'c'))


def iv_mul(a, b):
    ps = ((a.lo, b.lo), (a.lo, b.hi), (a.hi, b.lo), (a.hi, b.hi))
    return _w(min(M.mpf_mul(p, q, PREC, 'f') for p, q in ps),
              max(M.mpf_mul(p, q, PREC, 'c') for p, q in ps))


def iv_div(a, b):
    assert M.mpf_cmp(b.lo, M.from_rational(0, 1, PREC)) > 0
    ps = ((a.lo, b.lo), (a.lo, b.hi), (a.hi, b.lo), (a.hi, b.hi))
    return _w(min(M.mpf_div(p, q, PREC, 'f') for p, q in ps),
              max(M.mpf_div(p, q, PREC, 'c') for p, q in ps))


PI = Iv(M.mpf_pi(PREC, 'f'), M.mpf_pi(PREC, 'c'))
PI2 = iv_mul(PI, PI)
SQ2 = Iv(M.mpf_sqrt(M.from_rational(2, 1, PREC), PREC, 'f'),
         M.mpf_sqrt(M.from_rational(2, 1, PREC), PREC, 'c'))


def aiv(x):
    """Alg -> 严格区间 (sqrt2 用 [SQ2.lo,SQ2.hi], 按 b 的符号定向)."""
    base = iv_rat(x.a.numerator, x.a.denominator)
    if x.b == 0:
        return base
    r = iv_rat(abs(x.b.numerator), x.b.denominator)
    pr = iv_mul(r, SQ2)
    if x.b > 0:
        return iv_add(base, pr)
    return iv_sub(base, pr)


def f(x, n=48):
    return M.to_str(x, n)


def wd(iv):
    return M.to_str(M.mpf_sub(iv.hi, iv.lo, PREC, 'c'), 4)


# ---------------------------------------------------------------- 一个 theta 的全套计算
SIN_E = [Fr((-1) ** j, math.factorial(2 * j + 1)) for j in range(NK + 2)]   # sin u 的偶/奇索引


def ebar(B):
    """Ebar(u)=e^{-B(1-cos u)} 的偶次幂系数 e[m] (u^{2m}); ODE 精确递推."""
    e = [Alg(0)] * (NK + 1)
    e[0] = Alg(1)
    for n in range(NK):
        acc = Alg(0)
        for j in range(n + 1):
            acc = Aa(acc, Ac(e[n - j], SIN_E[j]))
        e[n + 1] = Am(acc, Ac(B, Fr(-1, 2 * (n + 1))))
    return e


def f0_series():
    """四个整函数被积因子 (去掉 Ebar), 返回 dict: odd k -> rational 系数."""
    W, A, B4, Sf = {}, {}, {}, {}
    for m in range(1, NK + 1):        # k = 2m-1
        k = 2 * m - 1
        # sin^2 u / u  : 系数 (-1)^{m+1} 2^{2m-1}/(2m)!
        W[k] = Fr((-1) ** (m + 1) * 2 ** (2 * m - 1), math.factorial(2 * m))
        # sin^3 u/(2u^2) : (-1)^m (3-3^{2m+1}) / (8 (2m+1)!)
        A[k] = Fr((-1) ** m * (3 - 3 ** (2 * m + 1)), 8 * math.factorial(2 * m + 1))
        # sin^4 u/(4u^3): 幂 k=2m-1 对应 sin^4 的 u^{2n} 项 n=m+1: (-1)^n (16^n-4*4^n)/(32 (2n)!)
        n = m + 1
        B4[k] = Fr((-1) ** n * (16 ** n - 4 * 4 ** n), 32 * math.factorial(2 * n))
        # Sf = sin^2u/(2u^3) - sin 2u/(4u^2)
        #  sin^2u/(2u^3):  k=2n-3 系数 (-1)^{n+1} 2^{2n-2}/(2n)!   (n=m+1)
        #  sin2u/(4u^2):   k=2n-1 系数 (-1)^n 2^{2n-1}/(2n+1)!     (n=m)
        Sf[k] = (Fr(((-1) ** (m + 2)) * 2 ** (2 * m), math.factorial(2 * m + 2))
                 - Fr(((-1) ** m) * 2 ** (2 * m - 1), math.factorial(2 * m + 1)))
    return W, A, B4, Sf


W0, A0, B0, S0 = f0_series()


def prod(e, f0):
    """Ebar(偶) x f0(奇) -> 奇次幂系数 (Alg)."""
    out = {}
    for m in range(NK + 1):
        for k, c in f0.items():
            kk = k + 2 * m
            if kk > KMAX:
                break
            t = Ac(e[m], c)
            out[kk] = Aa(out[kk], t) if kk in out else t
    return out


# ---- tail: Cauchy 界; r = 4.2, 严格有理上界 sinh(4.2)<=34 (e^4.2<=67), sinh(8.4)<=2230,
#      cosh(4.2)-1 <= 33.4, pi < 3.2 (=> pi/r < 16/21)
R = Fr(21, 5)
SINH_R = Fr(34)
SINH_2R = Fr(2230)
COSH_M1 = Fr(167, 5)
PI_UB = Fr(16, 5)
RATIO_UB = PI_UB / R       # 16/21


def mfac(name):
    return {"W": SINH_R ** 2 / R,
            "A": SINH_R ** 3 / (2 * R ** 2),
            "B": SINH_R ** 4 / (4 * R ** 3),
            "S": SINH_R ** 2 / (2 * R ** 3) + SINH_2R / (4 * R ** 2)}[name]


def tail_ub(name, B_abs_up):
    ex = M.mpf_exp(M.mpf_mul(B_abs_up, M.from_rational(COSH_M1.numerator, COSH_M1.denominator, PREC), PREC, 'c'), PREC, 'c')
    mf = mfac(name)
    Mup = M.mpf_mul(ex, M.from_rational(mf.numerator, mf.denominator, PREC), PREC, 'c')
    geo = M.from_rational((RATIO_UB ** (KMAX + 1)).numerator, (RATIO_UB ** (KMAX + 1)).denominator, PREC, 'c')
    comb = PI_UB / (Fr(KMAX + 2) * (1 - RATIO_UB))
    t = M.mpf_mul(M.mpf_mul(Mup, geo, PREC, 'c'),
                  M.from_rational(comb.numerator, comb.denominator, PREC), PREC, 'c')
    return M.mpf_add(t, EPS, PREC, 'c')


def integ(ser, name, B_abs_up, pe=None):
    acc = iv_rat(0)
    pw = PI2
    zero = M.from_rational(0, 1, PREC)
    mx = zero
    for k in range(1, KMAX + 1, 2):
        c = ser.get(k)
        if c is not None and (c.a != 0 or c.b != 0):
            term = iv_div(iv_mul(aiv(c), pw), iv_rat(k + 1))
            acc = iv_add(acc, term)
            if M.mpf_cmp(M.mpf_abs(term.lo), mx) > 0:
                mx = M.mpf_abs(term.lo)
            if M.mpf_cmp(M.mpf_abs(term.hi), mx) > 0:
                mx = M.mpf_abs(term.hi)
        pw = iv_mul(pw, PI2)
    if pe is not None:
        pe[name] = mx
    t = tail_ub(name, B_abs_up)
    return Iv(M.mpf_sub(acc.lo, t, PREC, 'f'), M.mpf_add(acc.hi, t, PREC, 'c')), t


def certify(th_lab, ca, cb, beta=Fr(5, 2)):
    """c = ca + cb*sqrt2; s^2 = 1-c^2; B = beta c."""
    B = Alg(beta * ca, beta * cb)
    Bab = M.mpf_add(M.mpf_mul(M.from_rational(abs(B.a.numerator), abs(B.a.denominator), PREC),
                              M.from_rational(1, 1, PREC), PREC, 'c'),
                    M.mpf_mul(M.from_rational(abs(B.b.numerator), abs(B.b.denominator), PREC),
                              SQ2.hi, PREC, 'c'), PREC, 'c')
    e = ebar(B)
    PE = {}
    Z, tZ = integ(prod(e, W0), "W", Bab, PE)
    NA, tA = integ(prod(e, A0), "A", Bab, PE)
    NB, tB = integ(prod(e, B0), "B", Bab, PE)
    NS, tS = integ(prod(e, S0), "S", Bab, PE)
    E1 = iv_div(NA, Z)
    E2 = iv_div(NB, Z)
    VG = iv_sub(E2, iv_mul(E1, E1))
    S_h = iv_mul(aiv(Alg(beta * ca / 2, beta * cb / 2)), E1)     # (beta c/2) E[g]
    S_p = iv_div(NS, Z)
    S_iv = iv_add(S_h, S_p)
    s2 = Fr(1) - (ca * ca + 2 * cb * cb)
    assert s2 >= 0
    D_iv = iv_mul(iv_rat((beta * beta * s2).numerator, (beta * beta * s2).denominator), VG)
    gap = M.mpf_sub(S_iv.lo, D_iv.hi)
    return dict(lab=th_lab, c=(ca, cb), s2=s2, Z=Z, NA=NA, NB=NB, NS=NS, tZ=tZ, tA=tA, tB=tB, tS=tS,
                E1=E1, E2=E2, VG=VG, Sh=S_h, Sp=S_p, S=S_iv, D=D_iv, gap=gap, PE=PE,
                sign=M.mpf_cmp(gap, M.from_rational(0, 1, PREC)))


def pp(nm, iv):
    return "  %-9s ∈ [%s ,\n               %s]   width=%s" % (nm, f(iv.lo), f(iv.hi), wd(iv))


CASES = [("theta=0    (c=1, s=0)", Fr(1), Fr(0)),
         ("theta=pi/2 (c=s=sqrt2/2)", Fr(0), Fr(1, 2)),
         ("theta=pi   (c=0, s=1)", Fr(0), Fr(0))]

P("=" * 104)
P("r03_cert_theta.log -- R03: theta != pi 的严格区间证书   beta=2.5, SU(2), 单 plaquette, 切片 xi3=0")
P("  prec=%d bit, 相对外扩 2^-%d, KMAX=%d, Q[sqrt2] 精确系数, tail=Cauchy(r=4.2)" % (PREC, PREC // 2, KMAX))
P("=" * 104)
P("")
P("[0] 积分表示 (对一切 theta):")
P("  nu ∝ Ebar(u) sin^2(u)/u du,  Ebar=e^{-beta(1-c cos u)},  u in (0,pi)")
P("  S = (beta c/2) E_nu[g] + E_nu[Phi''],   D = beta^2 s^2 Var_nu[g],   Phi'≡0 (切片切于等-rho 球面)")
P("  归约核对: 我核 Phi' = Phi_rho'*(xi3/rho) = 0 逐点; H=(beta c/2)g(u) != 0 当 c!=0  [T]")
P("")

res = []
for lab, ca, cb in CASES:
    r = certify(lab, ca, cb)
    res.append(r)
    P("-" * 104)
    P("[%s]" % lab)
    P("  Z   tail<=%s   NA tail<=%s   NB tail<=%s   NS tail<=%s" % (f(r['tZ'], 4), f(r['tA'], 4), f(r['tB'], 4), f(r['tS'], 4)))
    P(pp("Z", r['Z']))
    P(pp("E[g]", r['E1']))
    P(pp("Var(g)", r['VG']))
    P(pp("E[Phi'']", r['Sp']))
    P(pp("(bc/2)E[g]", r['Sh']))
    P(pp("S", r['S']))
    P(pp("D", r['D']))
    P("  S_lower - D_upper = %s   SIGN=%+d" % (f(r['gap']), r['sign']))
    P("")

# 回归: theta=pi 必须重现 q11_cert 的 0.2051152971230731 / 0.0691205860810849
P("[R] 回归 (theta=pi 应与 q11_cert.log 逐位一致):")
P("  S = %s   (q11: 0.205115297123073147983329003178018512976649324657)" % f(res[2]['S'].lo))
P("  D = %s   (q11: 0.0691205860810849369464295361378203906866969638109)" % f(res[2]['D'].lo))
P("  S-D = %s  (q11: +0.135994711041988211036899467040198122289952360846)" % f(res[2]['gap']))
P("")
P("[C] 独立浮点核对 (复合 Simpson, N=20000, 与上面区间值互相对照; 仅 [C]):")


def simpson(th):
    import math
    c, s = math.cos(th / 2), math.sin(th / 2)
    B = 2.5 * c
    N = 20000
    h = math.pi / N

    def g(u):
        return math.sin(u) / (2 * u) if u else 0.5

    def w(u):
        return math.sin(u) ** 2 / u if u else 0.0

    def wt(u):
        return math.exp(-B * (1 - math.cos(u)))

    def phipp(u):
        return 1 / (2 * u * u) - (math.cos(u) / math.sin(u)) / (2 * u) if u else 1.0 / 6

    def quad(f):
        t = 0.0
        for i in range(N + 1):
            u = i * h
            wgt = 1 if i in (0, N) else (4 if i % 2 else 2)
            t += wgt * f(u)
        return t * h / 3

    Z = quad(lambda u: wt(u) * w(u))
    NA = quad(lambda u: wt(u) * g(u) * w(u))
    NB = quad(lambda u: wt(u) * g(u) ** 2 * w(u))
    NS = quad(lambda u: wt(u) * phipp(u) * w(u))
    S = 2.5 * c / 2 * NA / Z + NS / Z
    D = 6.25 * s * s * (NB / Z - (NA / Z) ** 2)
    return S, D, S - D


import math
for th, nm in ((0.0, "theta=0"), (math.pi / 2, "theta=pi/2"), (math.pi, "theta=pi")):
    Sf, Df, Gf = simpson(th)
    P("  %-11s S=%.12f  D=%.12f  S-D=%.12f" % (nm, Sf, Df, Gf))
P("  (对照: 区间证书 S-D = %s / %s / %s )"
  % (f(res[0]['gap'], 12), f(res[1]['gap'], 12), f(res[2]['gap'], 12)))
P("")
P("[P] 峰值项 (灾难性相消证据): " + ", ".join("%s=%.2e" % (k, float(M.to_str(v, 6))) for k, v in res[1]["PE"].items()))

with open("r03_cert_theta.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG) + "\n")

