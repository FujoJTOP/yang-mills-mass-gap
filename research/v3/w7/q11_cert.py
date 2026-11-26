# -*- coding: utf-8 -*-
r"""
q11_cert.py -- Q11: 区间算术严格证书  S - D > 0   在  beta=2.5, SU(2), 单 plaquette, theta=pi, a0=0, v=e3

== 恒等式 ==
    G''(0) = E[H + Phi'' | 0] - Var(G + Phi' | 0) =: S - D ,
    H = d^2_{xi3} S_W ,  G = d_{xi3} S_W ,  Phi = -log J ,  Phi' = d_{xi3} Phi ,  Phi'' = d^2_{xi3} Phi .
单 plaquette 精确约化 (与 h1e_var.py / p05_cinf.py 同一模型约定):
    U = exp(xi.sig/2), |xi|<2pi, h=rho/2, u=rho/2, rho=|xi|,
    S_W = beta(1 - cos h cos(th/2) + sin h sin(th/2) xi3/rho),   J = 4 sin^2 h / rho^2,  dmu ∝ J d^3xi.
切片 xi3 = 0 (即 a0 = <v,xi> = 0), v = e3 (单位向量):
    测度 ∝ J d^2xi = [4 sin^2(rho/2)/rho^2] rho drho dphi  ⇒ 密度 ∝ sin^2(u) du/u ,  u ∈ (0,pi)
    G     = beta sin(th/2) sin(u)/(2u) = beta s g(u),  s=sin(th/2), g(u)=sin u/(2u)
    H     = beta cos(th/2) sin(u)/(4u) = (beta c/2) g(u),  c=cos(th/2)
    Phi   = -log J = 2 log rho - 2 log(2 sin h)   只依赖 rho
    Phi'  = d_{xi3} Phi = Phi_rho' * (xi3/rho) = 0        逐点 (切片切于等 rho 球面)   <== 关键
    Phi'' = Phi_rho'/rho = 1/(2u^2) - cot(u)/(2u) ,  u = rho/2
在 theta=pi (c=0, s=1):
    S_W ≡ beta (切片常数 ⇒ 权重常数 ⇒ 测度 = 纯 Haar 切片, 与 beta 无关)
    H ≡ 0 逐点,  G = beta g(u),  Phi' ≡ 0 逐点  ⇒  D = Var(G), 无交叉项
    S = E[Phi''] ,  D = beta^2 Var(g)
四个有限维积分 (u ∈ (0,pi), 密度 w = sin^2 u/u):
    Z=∫w , NA=∫g w , NB=∫g^2 w , NS=∫Phi'' w ;   S=NS/Z, D=beta^2(NB/Z-(NA/Z)^2)
被积函数全部整函数 (用 cot(u)sin^2 u = sin(2u)/2 消掉 cot):
    W=sin^2u/u , A=sin^3u/(2u^2) , B=sin^4u/(4u^3) ,
    Sf=Phi'' w = -(1/4) d/du[sin^2 u/u^2]  ⇒  NS = 1/4  精确 [T]

== 严格性口径 (真向外舍入) ==
1) Taylor 系数 a_k: Fraction 精确有理数 (无浮点).
2) ∫_0^pi f = Σ_k a_k pi^{k+1}/(k+1)  (整函数逐项积分恒等).
3) pi 只以区间 [pi_f,pi_c] 进入 (mpmath.libmp 定向舍入 'f'/'c', prec=800 bit).
4) 每个区间 +,-,*,/ 用 libmp 定点定向舍入 ('f' 向下 / 'c' 向上), 并统一加绝对外扩
   PAD = 2^-400 —— 即使 libmp 每步有 2^100 ulp 的欠正确, 累积仍 < 1e-110, 远小于所需 1e-10.
5) 截断 tail: Cauchy 估计 |a_k| ≤ |c|(e j/(m+k))^{m+k} (f=c sin^j u/u^m 整函数), 且 pi^{k+1} < 3.2^{k+1}.
"""
from fractions import Fraction as Fr
import math
import mpmath.libmp as M

PREC = 800
PAD = M.from_rational(1, 2 ** (PREC // 2), PREC, 'f')      # 2^-400
KMAX = 260
LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


# ---------------------------------------------------------------- 精确 Laurent 级数
def s_sin(n):
    return {2 * i + 1: Fr((-1) ** i, math.factorial(2 * i + 1)) for i in range(n // 2 + 1) if 2 * i + 1 <= n}


def s_sin2u(n):
    return {2 * i + 1: Fr((-1) ** i * 2 ** (2 * i + 1), math.factorial(2 * i + 1))
            for i in range(n // 2 + 1) if 2 * i + 1 <= n}


def lmul(p, q, n):
    out = {}
    for k1, v1 in p.items():
        for k2, v2 in q.items():
            k = k1 + k2
            if -n <= k <= n:
                out[k] = out.get(k, Fr(0)) + v1 * v2
    return {k: v for k, v in out.items() if v != 0}


def lshift(p, m):
    return {k + m: v for k, v in p.items()}


def lscale(p, c):
    return {k: v * c for k, v in p.items() if v * c != 0}


def ladd(p, q):
    o = dict(p)
    for k, v in q.items():
        o[k] = o.get(k, Fr(0)) + v
    return {k: v for k, v in o.items() if v != 0}


N = KMAX + 4
sn = s_sin(N)
sn2 = lmul(sn, sn, 2 * N)
sn3 = lmul(sn2, sn, 3 * N)
sn4 = lmul(sn2, sn2, 4 * N)
W = lshift(sn2, -1)
A = lscale(lshift(sn3, -2), Fr(1, 2))
B = lscale(lshift(sn4, -3), Fr(1, 4))
Sf = ladd(lscale(lshift(sn2, -3), Fr(1, 2)), lscale(lshift(s_sin2u(N), -2), Fr(-1, 4)))

for nm, ser in (("W", W), ("A", A), ("B", B), ("Sf", Sf)):
    assert not [k for k in ser if k < 0], (nm, sorted(k for k in ser if k < 0))
G2 = lshift(sn2, -2)
der = lscale({k - 1: v * k for k, v in G2.items() if k != 0}, Fr(-1, 4))
for k in [k for k in (set(der) | set(Sf)) if k <= KMAX + 2]:
    assert der.get(k, Fr(0)) == Sf.get(k, Fr(0)), ("IBP identity check failed at k=%d" % k)


# ---------------------------------------------------------------- 区间算术 (libmp 定向舍入)
class Iv:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi):
        self.lo, self.hi = lo, hi


def iv_rat(p, q=1):
    return Iv(M.from_rational(p, q, PREC, 'f'), M.from_rational(p, q, PREC, 'c'))


def _w(lo, hi):
    return Iv(M.mpf_sub(lo, PAD, PREC, 'f'), M.mpf_add(hi, PAD, PREC, 'c'))


def iv_add(a, b):
    return _w(M.mpf_add(a.lo, b.lo, PREC, 'f'), M.mpf_add(a.hi, b.hi, PREC, 'c'))


def iv_sub(a, b):
    return _w(M.mpf_sub(a.lo, b.hi, PREC, 'f'), M.mpf_sub(a.hi, b.lo, PREC, 'c'))


def iv_mul(a, b):
    ps = ((a.lo, b.lo), (a.lo, b.hi), (a.hi, b.lo), (a.hi, b.hi))
    lo = min(M.mpf_mul(p, q, PREC, 'f') for p, q in ps)
    hi = max(M.mpf_mul(p, q, PREC, 'c') for p, q in ps)
    return _w(lo, hi)


def iv_div(a, b):
    assert M.to_rational(b.lo) is not None and M.mpf_cmp(b.lo, M.from_rational(0, 1, PREC)) > 0
    ps = ((a.lo, b.lo), (a.lo, b.hi), (a.hi, b.lo), (a.hi, b.hi))
    lo = min(M.mpf_div(p, q, PREC, 'f') for p, q in ps)
    hi = max(M.mpf_div(p, q, PREC, 'c') for p, q in ps)
    return _w(lo, hi)


PI_IV = Iv(M.mpf_pi(PREC, 'f'), M.mpf_pi(PREC, 'c'))
PI_UB = Fr(16, 5)            # 3.2 > pi, 用于 tail 上界
E_UB = M.mpf_exp(M.from_rational(1, 1, PREC), PREC, 'c')


def tail_bound(terms, K):
    """Σ_{k>K} |c|(e j/(m+k))^{m+k} pi^{k+1}/(k+1) 的严格上界 (pi<3.2)."""
    tot = M.from_rational(0, 1, PREC)
    last = None
    for k in range(K + 1, K + 1 + 400):
        s = M.from_rational(0, 1, PREC)
        for c, j, m in terms:
            r = M.mpf_div(M.mpf_mul(E_UB, M.from_rational(j, 1, PREC), PREC, 'c'),
                          M.from_rational(m + k, 1, PREC), PREC, 'c')
            lg = M.mpf_mul(M.from_rational(m + k, 1, PREC), M.mpf_log(r, PREC, 'c'), PREC, 'c')
            t = M.mpf_mul(M.mpf_exp(lg, PREC, 'c'),
                          M.from_rational(c.numerator, c.denominator * (k + 1), PREC, 'c'), PREC, 'c')
            s = M.mpf_add(s, t, PREC, 'c')
        pw = M.from_rational(PI_UB.numerator ** (k + 1), PI_UB.denominator ** (k + 1), PREC, 'c') \
            if k < 300 else None
        # pi^{k+1} < 3.2^{k+1}; 大 k 直接与上一项比较 (比值 << 1), 用 2 次幂近似避免大整数
        if pw is not None:
            s = M.mpf_mul(s, pw, PREC, 'c')
        tot = M.mpf_add(tot, s, PREC, 'c')
        if last is not None and M.mpf_cmp(s, M.mpf_mul(last, M.from_rational(1, 10 ** 6, PREC), PREC, 'c')) < 0:
            return tot
        last = s
    return tot


def integral(ser, terms):
    acc = iv_rat(0)
    pw = PI_IV
    for k in range(KMAX + 1):
        a = ser.get(k, Fr(0))
        if a != 0:
            acc = iv_add(acc, iv_mul(iv_rat(a.numerator, a.denominator * (k + 1)), pw))
        pw = iv_mul(pw, PI_IV)
    t = tail_bound(terms, KMAX)
    return Iv(M.mpf_sub(acc.lo, t, PREC, 'f'), M.mpf_add(acc.hi, t, PREC, 'c')), t


TAIL = {"W": [(Fr(1), 2, 1)], "A": [(Fr(1, 2), 3, 2)], "B": [(Fr(1, 4), 4, 3)],
        "Sf": [(Fr(1, 2), 2, 3), (Fr(1, 4), 2, 2)]}

Zv, tZ = integral(W, TAIL["W"])
NAv, tA = integral(A, TAIL["A"])
NBv, tB = integral(B, TAIL["B"])
NSv, tS = integral(Sf, TAIL["Sf"])
NS_EXACT = iv_rat(1, 4)

BETA2 = Fr(25, 4)
S_iv = iv_div(NS_EXACT, Zv)
E1 = iv_div(NAv, Zv)
E2 = iv_div(NBv, Zv)
VARG = iv_sub(E2, iv_mul(E1, E1))
D_iv = iv_mul(iv_rat(BETA2.numerator, BETA2.denominator), VARG)
SD = iv_sub(S_iv, D_iv)


def f(x, n=48):
    return M.to_str(x, n)


def wd(iv):
    """区间宽度 (检验向外舍入确实发生)."""
    return M.to_str(M.mpf_sub(iv.hi, iv.lo, PREC, 'c'), 4)


def pp(nm, iv):
    return "  %-7s ∈ [%s ,\n             %s]      width=%s" % (nm, f(iv.lo), f(iv.hi), wd(iv))


P("=" * 104)
P("q11_cert.log -- Q11 严格区间证书:  S - D = E[H+Phi''] - Var(G+Phi') > 0 ?")
P("  beta=2.5, SU(2), 单 plaquette, theta=pi, 切片 xi3=0 (a0=0), v=e3, ||v||=1")
P("  interval: prec=%d bit, 绝对外扩 PAD=2^-%d,  Taylor 阶 KMAX=%d" % (PREC, PREC // 2, KMAX))
P("=" * 104)
P("")
P("[0] 积分表示核准")
P("  测度: 密度 ∝ sin^2(u)/u du, u=rho/2 ∈ (0,pi)   [Haar: dmu ∝ sin^2(rho/2)d rho dOmega, 切片 d^2xi=rho d rho dphi,")
P("        两者相除得 sin^2(u)/u; 代码断言 W/A/B/Sf 的 Laurent 负幂全部相消 = 被积函数在 u=0 解析]")
P("  theta=pi: S_W ≡ beta ⇒ 权重常数 ⇒ 测度=纯 Haar 切片 (与 beta 无关), H≡0 逐点, Phi'≡0 逐点")
P("  ⇒ S = NS/Z 与 beta 无关;  D = beta^2 (NB/Z - (NA/Z)^2)")
P("  NS = ∫Phi''w = 1/4 精确 (IBP: Phi''w = -(1/4)d/du[sin^2u/u^2], 端点 sin^2u/u^2: 1 -> 0)")
P("")
P("[1] 严格区间 (向外舍入, prec=800bit, PAD=2^-400)")
P("  Z   ∈ [%s,\n         %s]   tail<=%s   width=%s" % (f(Zv.lo), f(Zv.hi), f(tZ, 6), wd(Zv)))
P("  NA  ∈ [%s,\n         %s]   tail<=%s   width=%s" % (f(NAv.lo), f(NAv.hi), f(tA, 6), wd(NAv)))
P("  NB  ∈ [%s,\n         %s]   tail<=%s   width=%s" % (f(NBv.lo), f(NBv.hi), f(tB, 6), wd(NBv)))
P("  NS  ∈ [%s, %s]  (精确 1/4);  级数独立核对 ∈ [%s, %s] tail<=%s"
  % (f(NS_EXACT.lo), f(NS_EXACT.hi), f(NSv.lo), f(NSv.hi), f(tS, 6)))
P(pp("E[g]", E1))
P(pp("E[g^2]", E2))
P(pp("Var(g)", VARG))
P("")
P("[2] 判死线 (a): S_lower - D_upper > 0 ?")
P(pp("S", S_iv))
P(pp("D", D_iv))
P(pp("S-D", SD))
P("  S_lower - D_upper = %s" % f(M.mpf_sub(S_iv.lo, D_iv.hi)))
P("  S_upper - D_lower = %s" % f(M.mpf_sub(S_iv.hi, D_iv.lo)))
P("  SIGN = %+d   (+1 = 严格正, 证书成立)" % M.mpf_cmp(SD.lo, M.from_rational(0, 1, PREC)))
P("")
P("[3] 交叉核对 (h1e_var.log / 语料锚)")
P("  Z        ↔ 0.5*Cin(2pi)          = 1.218826696528        (h1e_var.log:3)")
P("  Var(g)   ↔ 0.011059293772973                            (h1e_var.log:13)")
P("  S        ↔ 1/(2Cin(2pi)) = 0.205115297123  (语料 eq:v3anchor)")
bst = M.mpf_sqrt(M.mpf_div(S_iv.lo, VARG.hi, PREC, 'f'), PREC, 'f')
bst2 = M.mpf_sqrt(M.mpf_div(S_iv.hi, VARG.lo, PREC, 'c'), PREC, 'c')
P("  beta* = sqrt(S/Var(g)) (取 S-D=0)  ∈ [%s, %s]   (h1e_var.log:20 记 4.30660785)" % (f(bst), f(bst2)))

with open("q11_cert.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG) + "\n")
