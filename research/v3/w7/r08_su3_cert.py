# -*- coding: utf-8 -*-
r"""
r08_su3_cert.py -- R08: SU(3) 上 q11/r03 型严格证书的**二维 Weyl 积分**（精确有理 + 800bit 定向舍入）

== 模型 (SU(3), Weyl/环面形式, theta=pi 切片 => 权重常数 => 测度 = 纯 Haar) ==
    链接变量 U ∈ T ⊂ SU(3):  θ=(θ1,θ2,θ3), θ1+θ2+θ3 ≡ 0 (mod 2π), θ3 = -θ1-θ2
    Weyl 公式 (类函数):  ∫_{SU(3)} F dU = (1/3!)(2π)^{-2} ∫_T  ∏_{i<j}|e^{iθi}-e^{iθj}|^2 F(θ) dθ1dθ2
    Vandermonde:  J(θ) := ∏_{i<j} 4 sin^2((θi-θj)/2) = ∏_{i<j}(2-2cos(θi-θj))
    Φ := -log J ,  Φ' = Σ_k λ_k ∂_k Φ ,  Φ'' = Σ_{kl} λ_k λ_l ∂^2_{kl} Φ
    G := ∂_τ S_W ,  H := ∂^2_τ S_W   (τ: θ -> θ + τλ, λ ∈ Cartan, Σλ=0)
    S := E[H + Φ''] ,  D := Var(G + Φ') = E[(G+Φ')^2] - E[G+Φ']^2
    G''(0) = S - D

== 本文件的判定 ==
    theta=pi 切片: S_W ≡ β ⇒ 权重常数 ⇒ G ≡ 0, H ≡ 0, ν = |Δ|^2 dθ / Z (纯 Haar)
    ⇒ S = E[Φ''] = NS/Z ,  D = Var(Φ') = NP/Z   (E[Φ']=0 由 τ-平移不变性)
    恒等式 (沿平移方向的 Weyl–Schwinger / IBP):  ∫J Φ'' = ∫J Φ'^2  ⇒  S ≡ D
    ⇒ S_lo - D_up < 0  (区间仅因外扩 PAD 而重叠, 真差 = 0, 精确有理数相等)

== 严格性口径 ==
    1) |Δ|^2 Φ'' 与 |Δ|^2 Φ'^2 都是 (θ1,θ2) 的**三角多项式**:
       s_ij cot_ij = 2 sin x_ij ,  s_ij cot^2_ij = 4 - s_ij  (s_ij := 4sin^2(x_ij/2))
       所以 s_ij 的分母 sin^2 被 |Δ|^2 的分子整除, cot 全部消掉.
    2) ∫ over [0,2π)^2 of trig poly = (2π)^2 × 常数项, 常数项用 **Fraction 精确**算出 (无浮点).
    3) π 只以区间 [π_f, π_c] 进入 (libmp 定向舍入 'f'/'c', prec=800); 每个 +,-,*,/ 走定向舍入,
       统一加绝对外扩 PAD = 2^-400 (真向外舍入的机检指纹: 宽度 ≈ 操作数×PAD).
"""
from fractions import Fraction as Fr
import mpmath.libmp as M

PREC = 800
PAD = M.from_rational(1, 2 ** (PREC // 2), PREC, 'f')
LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


# ------------------------------------------------------------------ 高斯有理数
class Gq:
    __slots__ = ("a", "b")

    def __init__(self, a, b=Fr(0)):
        self.a, self.b = Fr(a), Fr(b)

    def __add__(self, o):
        return Gq(self.a + o.a, self.b + o.b)

    def __sub__(self, o):
        return Gq(self.a - o.a, self.b - o.b)

    def __mul__(self, o):
        return Gq(self.a * o.a - self.b * o.b, self.a * o.b + self.b * o.a)

    def __neg__(self):
        return Gq(-self.a, -self.b)

    def __eq__(self, o):
        return self.a == o.a and self.b == o.b

    def __hash__(self):
        return hash((self.a, self.b))

    def __repr__(self):
        return "%s%+si" % (self.a, self.b)


ZERO = Gq(0)
ONE = Gq(1)
HALF_I_INV = Gq(0, Fr(-1, 2))          # 1/(2i)


def gnum(p, q=1):
    return Gq(Fr(p, q))


# ------------------------------------------------------------------ Laurent 多项式 (u=e^{i(θ1-θ2)}, v=e^{i(θ1-θ3)})
def lmul(p, q):
    o = {}
    for (a1, b1), c1 in p.items():
        for (a2, b2), c2 in q.items():
            k = (a1 + a2, b1 + b2)
            o[k] = o.get(k, ZERO) + c1 * c2
    return {k: v for k, v in o.items() if v != ZERO}


def ladd(*ps):
    o = {}
    for p in ps:
        for k, v in p.items():
            o[k] = o.get(k, ZERO) + v
    return {k: v for k, v in o.items() if v != ZERO}


def lscale(p, c):
    return {k: v * c for k, v in p.items() if v * c != ZERO}


def lpow(p, n):
    r = {(0, 0): ONE}
    for _ in range(n):
        r = lmul(r, p)
    return r


U = {(1, 0): ONE, (-1, 0): ONE}
V = {(1, 0): ONE, (-1, 0): ONE}

# u12 = e^{i(θ1-θ2)}，u13 = e^{i(θ1-θ3)}，u23 = u13/u12
u12 = {(1, 0): ONE}
u12i = {(-1, 0): ONE}
u13 = {(0, 1): ONE}
u13i = {(0, -1): ONE}
u23 = {(-1, 1): ONE}
u23i = {(1, -1): ONE}
TWO = gnum(2)

s12 = ladd({(0, 0): TWO}, lscale(u12, gnum(-1)), lscale(u12i, gnum(-1)))
s13 = ladd({(0, 0): TWO}, lscale(u13, gnum(-1)), lscale(u13i, gnum(-1)))
s23 = ladd({(0, 0): TWO}, lscale(u23, gnum(-1)), lscale(u23i, gnum(-1)))
D2 = lmul(lmul(s12, s13), s23)                    # |Δ|^2 = J

# sin x_ij = (u - u^-1)/(2i) ;  s_ij cot_ij = 2 sin x_ij ;  s_ij cot^2_ij = 4 - s_ij
sn12 = lscale(ladd(u12, lscale(u12i, gnum(-1))), HALF_I_INV)
sn13 = lscale(ladd(u13, lscale(u13i, gnum(-1))), HALF_I_INV)
sn23 = lscale(ladd(u23, lscale(u23i, gnum(-1))), HALF_I_INV)
four_minus = lambda s: ladd({(0, 0): gnum(4)}, lscale(s, gnum(-1)))

# λ = (1, 0, -1)  (λ1+λ2+λ3 = 0, 超荷方向)
LAM = (Fr(1), Fr(0), Fr(-1))
d12, d13, d23 = LAM[0] - LAM[1], LAM[0] - LAM[2], LAM[1] - LAM[2]

# |Δ|^2 Φ'' = (1/2) Σ_{i<j} (λi-λj)^2 · 4 · s_kl s_mn      [csc^2(x/2)=4/s_ij]
def w4(w, s_a, s_b):
    return lscale(lmul(s_a, s_b), gnum(4 * w))

D2Phi2 = lscale(ladd(w4(d12 * d12, s13, s23),
                     w4(d13 * d13, s12, s23),
                     w4(d23 * d23, s12, s13)), gnum(1, 2))

# Phi' = -Σ_{i<j} (λi-λj) cot_ij   ⇒   |Δ|^2 Φ'^2
#   项 (ij)^2      : (λi-λj)^2 (4-s_ij) s_kl s_mn
#   项 (ij)(kl)    : 2(λi-λj)(λk-λl) (2 sin x_ij)(2 sin x_kl) s_mn
sq = ladd(
    lscale(lmul(lmul(four_minus(s12), s13), s23), gnum(1)),      # (λ1-λ2)^2 = 1
    lscale(lmul(lmul(four_minus(s13), s12), s23), gnum(4)),      # (λ1-λ3)^2 = 4
    lscale(lmul(lmul(four_minus(s23), s12), s13), gnum(1)),      # (λ2-λ3)^2 = 1
)
cross = ladd(
    lscale(lmul(lmul(lscale(sn12, gnum(2)), lscale(sn13, gnum(2))), s23), gnum(2 * d12 * d13)),
    lscale(lmul(lmul(lscale(sn12, gnum(2)), lscale(sn23, gnum(2))), s13), gnum(2 * d12 * d23)),
    lscale(lmul(lmul(lscale(sn13, gnum(2)), lscale(sn23, gnum(2))), s12), gnum(2 * d13 * d23)),
)
D2PhiP2 = ladd(sq, cross)
# 注: 系数 4,2,4 = 2·(λ1-λ2)(λ1-λ3)=2·1·2, 2·(λ1-λ2)(λ2-λ3)=2·1·1, 2·(λ1-λ3)(λ2-λ3)=2·2·1

# ------------------------------------------------------------------ 常数项 (精确)
def CT(p):
    c = p.get((0, 0), ZERO)
    assert c.b == 0, ("常数项非实!", c)
    return c.a


ctZ, ctA, ctB = CT(D2), CT(D2Phi2), CT(D2PhiP2)
ct_sq, ct_cross = CT(sq), CT(cross)

# ------------------------------------------------------------------ 定向舍入区间
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
    return _w(min(M.mpf_mul(x, y, PREC, 'f') for x, y in ps),
              max(M.mpf_mul(x, y, PREC, 'c') for x, y in ps))


def iv_div(a, b):
    ps = ((a.lo, b.lo), (a.lo, b.hi), (a.hi, b.lo), (a.hi, b.hi))
    return _w(min(M.mpf_div(x, y, PREC, 'f') for x, y in ps),
              max(M.mpf_div(x, y, PREC, 'c') for x, y in ps))


PI = Iv(M.mpf_pi(PREC, 'f'), M.mpf_pi(PREC, 'c'))
FOURPI2 = iv_mul(iv_rat(4), iv_mul(PI, PI))          # (2π)^2 的区间

Ziv = _w(M.mpf_mul(M.from_rational(ctZ.numerator, ctZ.denominator, PREC, 'f'), FOURPI2.lo, PREC, 'f'),
         M.mpf_mul(M.from_rational(ctZ.numerator, ctZ.denominator, PREC, 'c'), FOURPI2.hi, PREC, 'c'))
Aiv = _w(M.mpf_mul(M.from_rational(ctA.numerator, ctA.denominator, PREC, 'f'), FOURPI2.lo, PREC, 'f'),
         M.mpf_mul(M.from_rational(ctA.numerator, ctA.denominator, PREC, 'c'), FOURPI2.hi, PREC, 'c'))
Biv = _w(M.mpf_mul(M.from_rational(ctB.numerator, ctB.denominator, PREC, 'f'), FOURPI2.lo, PREC, 'f'),
         M.mpf_mul(M.from_rational(ctB.numerator, ctB.denominator, PREC, 'c'), FOURPI2.hi, PREC, 'c'))

S_iv = iv_div(Aiv, Ziv)
D_iv = iv_div(Biv, Ziv)
SD = iv_sub(S_iv, D_iv)


def f(x, n=40):
    return M.to_str(x, n)


def wd(iv):
    return M.to_str(M.mpf_sub(iv.hi, iv.lo, PREC, 'c'), 4)


# ------------------------------------------------------------------ 独立核对 (浮点, [C])
import numpy as np

def ct_float(p):
    return float(sum(v.a for v in p.values() if True) * 0)  # placeholder unused


def const_term_np(p):
    # 直接对三角多项式做周期梯形求和核对常数项 (独立于字典算法)
    N = 400
    t = 2 * np.pi * np.arange(N) / N
    T1, T2 = np.meshgrid(t, t, indexing='ij')
    val = np.zeros_like(T1, dtype=complex)
    for (a, b), c in p.items():
        val += float(c.a) * np.exp(1j * (a * T1 + b * T2))
    return val.mean().real


# ------------------------------------------------------------------ 输出
P("=" * 108)
P("r08_su3_cert.log -- SU(3) 二维 Weyl 积分严格证书: S - D > 0 ?")
P("  模型: 单 plaquette, 链接 U ∈ T ⊂ SU(3), theta=pi 切片 (S_W ≡ β ⇒ 权重常数 ⇒ 纯 Haar)")
P("  Weyl:  ∫_{SU(3)}F dU = (1/3!)(2π)^-2 ∫_T ∏_{i<j}4sin^2((θi-θj)/2) F dθ1dθ2,  θ3=-θ1-θ2")
P("  涨落方向 λ=(1,0,-1) (超荷),  S=E[Φ''],  D=Var(Φ')   [G≡H≡0 on the slice]")
P("  interval: prec=%d bit, 绝对外扩 PAD=2^-%d" % (PREC, PREC // 2))
P("=" * 108)
P("")
P("[1] 三角多项式常数项 (Fraction 精确, 无浮点)")
P("  CT[|Δ|^2]                     = %s" % ctZ)
P("  CT[|Δ|^2 Φ'']                 = %s" % ctA)
P("  CT[|Δ|^2 Φ'^2]                = %s" % ctB)
P("  检验 CT[|Δ|^2 Φ''] - CT[|Δ|^2 Φ'^2] = %s   (精确有理差; IBP 恒等式要求 = 0)" % (ctA - ctB))
P("  分解: CT[|Δ|^2Φ'^2] = CT[sq] + CT[cross] = %s + %s = %s" % (ct_sq, ct_cross, ct_sq + ct_cross))
P("  Weyl 归一化: ∫_T|Δ|^2 = (2π)^2·CT = 24π^2 ; (1/3!)(2π)^-2·24π^2 = %s  (须 = 1)"
  % (Fr(ctZ, 6)))
P("  独立核对 (周期梯形 %d×%d, [C]): CT 浮点 = %.12f / %.12f / %.12f"
  % (400, 400, const_term_np(D2), const_term_np(D2Phi2), const_term_np(D2PhiP2)))
P("")
P("[2] 严格区间 (定向舍入, prec=800bit, PAD=2^-400)")
P("  Z   = ∫_T|Δ|^2   ∈ [%s , %s]   width=%s" % (f(Ziv.lo), f(Ziv.hi), wd(Ziv)))
P("  A   = ∫_T|Δ|^2Φ''∈ [%s , %s]   width=%s" % (f(Aiv.lo), f(Aiv.hi), wd(Aiv)))
P("  B   = ∫_T|Δ|^2Φ'^2 ∈ [%s , %s] width=%s" % (f(Biv.lo), f(Biv.hi), wd(Biv)))
P("  S = A/Z ∈ [%s , %s]   width=%s" % (f(S_iv.lo), f(S_iv.hi), wd(S_iv)))
P("  D = B/Z ∈ [%s , %s]   width=%s" % (f(D_iv.lo), f(D_iv.hi), wd(D_iv)))
P("")
P("[3] 判死线: S_lower - D_upper > 0 ?")
P("  S_lower - D_upper = %s" % f(M.mpf_sub(S_iv.lo, D_iv.hi)))
P("  S_upper - D_lower = %s" % f(M.mpf_sub(S_iv.hi, D_iv.lo)))
P("  S-D (区间) ∈ [%s , %s]   width=%s" % (f(SD.lo), f(SD.hi), wd(SD)))
P("  SIGN(S_lo - D_up) = %+d" % M.mpf_cmp(S_iv.lo, D_iv.hi))
P("  ⇒ 真差 = 0 (精确有理恒等式 CT[A]-CT[B]=%s); 区间仅因 PAD 外扩而重叠." % (ctA - ctB))
P("")
P("[4] SU(2) 对照 (q11, θ=π, β=2.5)")
P("  S_SU(2)=0.205115297123073147983329003178 , D_SU(2)=0.069120586081084936946429536138")
P("  S-D    =+0.135994711041988211036899467040  (>0, 非恒等)")
P("  SU(3) Weyl 环面上 S ≡ D ⇒ 差 = 0 (无隙, 不是负)")

with open("r08_su3_cert.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG) + "\n")
