# -*- coding: utf-8 -*-
"""
gm1alg_check.py — GM1-ALGEBRA 的确定性核（无随机数、无 MC、无格点采样）
[0] 自由场边缘：ω_a ≡ m 精确（asinh/sinh 互逆）
[1] (1.2a)(1.2b) 精确（sympy 级数 + Fraction）
[2] 2D 分支：Δ_a ≡ Δ_cont ⇒ δ₁=Δ³/24, δ₂=Δ⁵/1920, c₁=c₂=0 精确
[3] 自由场多粒子扇区：K₂ = −Σω³/24 vs H³/24 = (Σω)³/24 的显式差异（精确有理数）
[4] 阿贝尔格点核：K₂ 的 (DF)² 系数与符号（闭式 1−J₀(2a sin(qa/2)) 精确到 a⁶；直接格点和核）
"""
from fractions import Fraction as F
import sympy as sp

log = []
def P(s=""):
    print(s); log.append(s)

a = sp.symbols('a')
Dv = sp.symbols('Delta', positive=True)
m = sp.symbols('m', positive=True)

# ---------------------------------------------------------------- [0]
P("=" * 80)
P("[0] 自由场边缘：Δ_a = (2/a)asinh(am/2) ⇒ ω_a = (2/a)sinh(aΔ_a/2) ≡ m （精确）")
w_free = sp.simplify(2 / a * sp.sinh(a * (2 / a * sp.asinh(a * m / 2)) / 2))
P(f"    ω_a = {sp.simplify(w_free)}   （差 = {sp.simplify(w_free - m)}）")

# ---------------------------------------------------------------- [1]
P("=" * 80)
P("[1] (1.2a)(1.2b) 精确核")
P("    Δ_a = ω_a − (a²/24)ω_a³ + (3a⁴/640)ω_a⁵ ,   ω_a = D + s1 a² + s2 a⁴")
P("    c1 = s1/D − D²/24 ;   c2 = s2/D − s1 D/8 + 3D⁴/640")
ok = True
for (D, s1, s2) in [(F(1), F(0), F(0)), (F(1), F(7, 10), F(-2, 5)),
                    (F(3, 2), F(-1, 3), F(5)), (F(5, 7), F(2), F(-11, 3))]:
    R = lambda x: sp.Rational(x.numerator, x.denominator)
    w = R(D) + R(s1) * a**2 + R(s2) * a**4
    Del = sp.series(w - sp.Rational(1, 24) * a**2 * w**3 + sp.Rational(3, 640) * a**4 * w**5,
                    a, 0, 5).removeO()
    p = sp.Poly(sp.expand(Del), a)
    c1_num = sp.nsimplify(p.coeff_monomial(a**2)) / R(D)
    c2_num = sp.nsimplify(p.coeff_monomial(a**4)) / R(D)
    c1_f = s1 / D - D * D / F(24)
    c2_f = s2 / D - s1 * D / F(8) + 3 * D**4 / F(640)
    good = (sp.simplify(c1_num - R(c1_f)) == 0) and (sp.simplify(c2_num - R(c2_f)) == 0)
    ok &= good
    P(f"    (D,s1,s2)=({D},{s1},{s2})   c1={c1_f}  c2={c2_f}   match={good}")
P(f"    ⇒ 4/4 OK = {ok}")

# ---------------------------------------------------------------- [2]
P("=" * 80)
P("[2] 2D：H_a = (g²/2)ΣĈ₂ 与 a 无关（无空间 plaquette）⇒ Δ_a ≡ Δ_cont ≡ Δ 精确")
w2 = sp.series(2 / a * sp.sinh(a * Dv / 2), a, 0, 6).removeO()
pol = sp.Poly(sp.expand(w2), a)
s1 = sp.simplify(pol.coeff_monomial(a**2)); s2 = sp.simplify(pol.coeff_monomial(a**4))
P(f"    s1 = {s1}  (=Δ³/24? {sp.simplify(s1 - Dv**3/24)==0})   s2 = {s2}  (=Δ⁵/1920? {sp.simplify(s2 - Dv**5/1920)==0})")
c1 = sp.simplify(s1 / Dv - Dv**2 / 24); c2 = sp.simplify(s2 / Dv - s1 * Dv / 8 + 3 * Dv**4 / 640)
P(f"    c1 = {c1}   c2 = {c2}")
P("    ⇒ 2D: δ₁ = Δ³/24 ≠ 0，而 c1 = 0 ≠ −Δ²/24。k=2（其实 C=0）但 δ₁≠0 ⇒ (SM-2) 的'⟺'假")
P("    数值：Δ=1 时 δ₁=+0.041667，−Δ²/24=−0.041667，c1=0")

# ---------------------------------------------------------------- [3]
P("=" * 80)
P("[3] 自由场（全离散 Wilson 转移矩阵）：K₂ 一体、H³/24 三体")
P("    E_a(k) = (2/a)asinh(aω_k/2) ⇒ K₂ = Σ_k(−ω_k³/24)n_k ;  H = Σω_k n_k ⇒ H³/24 = (Σω_k n_k)³/24")
P("    双粒子态： ⟨K₂ + H³/24⟩ = [−(ω1³+ω2³)+(ω1+ω2)³]/24 = ω1ω2(ω1+ω2)/8 ≠ 0")
allok = True
for (w1, w2) in [(F(1), F(2)), (F(1), F(1)), (F(1, 2), F(3, 2)), (F(1), F(4)), (F(3), F(5))]:
    lhs = (-(w1**3 + w2**3) + (w1 + w2)**3) / 24
    pred = w1 * w2 * (w1 + w2) / 8
    allok &= (lhs == pred)
    P(f"    (ω1,ω2)=({w1},{w2}):  K₂+H³/24 = {lhs}   ω1ω2(ω1+ω2)/8 = {pred}   match={lhs==pred}")
P(f"    ⇒ 5/5 OK = {allok}；单粒子扇区 (ω2=0) 才为 0 ⇒ '逐算子'的说法假")

# ---------------------------------------------------------------- [4]
P("=" * 80)
P("[4] 阿贝尔格点核（U(1)，A_2(x₁)=f(x₁)，其余为 0）")
P("    U_□(x)=exp(i a[f(x+a)−f(x)]) ⇒ 1−ReU_□ = 1−cos(aΔf)")
P("    连续统预言：(1/a)Σ_x(1−ReU) = ∫f'²/2 + a²∫f'f'''/24 + O(a⁴)")
P("    取 f=cos(qx), L=2πn/q：Σ_x = n[1−J₀(2a sin(qa/2))]")
q = sp.symbols('q', positive=True)
zz = 2 * a * sp.sin(q * a / 2)
ser = sp.simplify(sp.expand(sp.series(1 - sp.besselj(0, zz), a, 0, 9).removeO()))
polz = sp.Poly(sp.expand(ser), a)
z4 = sp.simplify(polz.coeff_monomial(a**4)); z6 = sp.simplify(polz.coeff_monomial(a**6))
P(f"    1−J₀(2a sin(qa/2)) = a⁴·({z4}) + a⁶·({z6}) + O(a⁸)")
P(f"      a⁴ 系数 ÷ (q²/4)     − 1 = {sp.simplify(z4/(q**2/4) - 1)}")
P(f"      a⁶ 系数 ÷ (−q⁴/48)   − 1 = {sp.simplify(z6/(-q**4/48) - 1)}")
P("    ⇒ 精确吻合：K₂ = +(1/g²)∫f'f'''/24 = −(1/g²)∫f''²/24 ≤ 0（kinetic 片负定，系数 1/24）")

import numpy as np
def P1_lattice(n, qv=1.0):
    """n₂=n₃=1 ⇒ 横向体积 = a²；L=2π/qv 一个周期；返回 (1/a)Σ_x(1−ReU_□) 与 a"""
    L = 2 * np.pi / qv
    aL = L / n
    x = np.arange(n) * aL
    U = np.exp(1j * aL * np.cos(qv * x))
    Up = U * np.conj(np.roll(U, -1))          # exp(i a[f(x+a)-f(x)])
    return aL, np.sum(1 - Up.real) / aL, L

P("    直接格点和核：n₂=n₃=1 ⇒ ∫d³x = a²L；预言 (1/a)Σ_x = a²L(q²/4 − a²q⁴/48)")
P("    故比值 R := (1/a)Σ_x / (a²Lq²/4) 应 = 1 − a²q²/12 + O(a⁴)：")
for n in [2**6, 2**8, 2**10, 2**12]:
    aL, P1, L = P1_lattice(n)
    base = aL**2 * L * 1.0 / 4.0
    R = P1 / base
    P(f"      n={n:5d}  a={aL:.4e}   R={R:.12f}   1−a²/12={1-aL**2/12:.12f}   R−(1−a²/12)={R-(1-aL**2/12):.3e}")

P("=" * 80)
P("全部确定性：无随机数、无 MC、无格点采样（格点和是确定性求和）。")
with open(r"D:\Dev\math\research\v3\species\gm1alg_check.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(log) + "\n")
