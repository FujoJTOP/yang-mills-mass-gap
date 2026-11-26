# -*- coding: utf-8 -*-
r"""
r08_ibp_check.py -- R08 的两条独立核对 (不属于区间证书本身, 只做 falsification)

[V1] 1D 酸测 (SU(2) 环面): J=4sin^2(θ/2), Φ=-logJ ⇒ JΦ''=2, JΦ'^2=2+2cosθ,
     ∫(JΦ''-JΦ'^2)=0 精确 ⇒ 关系式 s·cot^2=(4-s), s·csc^2=4 用对了.
[V2] 点态 IBP 恒等式 (SU(3), λ=(1,0,-1)): J(Φ'^2-Φ'') == d^2_τ J(θ+τλ)|_{τ=0} 逐点.
[V3] 网格均值 (= 三角多项式常数项): CT[|Δ|^2Φ''] == CT[|Δ|^2Φ'^2] == 48.
"""
import mpmath as mp
import numpy as np

mp.mp.dps = 40
LAM = (1, 0, -1)
PAIRS = ((0, 1), (0, 2), (1, 2))

print("=" * 90)
print("[V1] SU(2) 环面 1D 酸测")
th = mp.mpf('0.83')
J1 = 4 * mp.sin(th / 2) ** 2
P1 = -1 / mp.tan(th / 2)
P11 = mp.mpf(1) / 2 / mp.sin(th / 2) ** 2
print("  J Phi''        =", mp.nstr(J1 * P11, 20), " (应 = 2)")
print("  J Phi'^2       =", mp.nstr(J1 * P1 ** 2, 20), " (应 = 2+2cos θ =", mp.nstr(2 + 2 * mp.cos(th), 20), ")")
print("  ∫_0^{2π}(...)  =", mp.nstr(mp.quad(lambda t: 4 * mp.sin(t / 2) ** 2 *
      (1 / mp.tan(t / 2) ** 2 - mp.mpf(1) / 2 / mp.sin(t / 2) ** 2), [0, mp.pi, 2 * mp.pi]), 10))

print("[V2] SU(3) 点态 IBP 恒等式 (3 个随机点)")


def J_of(t):
    s = mp.mpf(1)
    for i, j in PAIRS:
        s *= 4 * mp.sin((t[i] - t[j]) / 2) ** 2
    return s


def phip(t):
    return sum(-mp.mpf(LAM[i] - LAM[j]) / mp.tan((t[i] - t[j]) / 2) for i, j in PAIRS)


def phipp(t):
    return sum(mp.mpf(1) / 2 * (LAM[i] - LAM[j]) ** 2 / mp.sin((t[i] - t[j]) / 2) ** 2 for i, j in PAIRS)


for a, b in ((mp.mpf('0.7'), mp.mpf('-1.1')), (mp.mpf('2.3'), mp.mpf('0.41')), (mp.mpf('-1.9'), mp.mpf('2.05'))):
    t = [a, b, -(a + b)]
    lhs = J_of(t) * (phip(t) ** 2 - phipp(t))
    rhs = mp.diff(lambda tau: J_of([t[i] + tau * LAM[i] for i in range(3)]), 0, 2)
    print("  θ=(%s,%s): J(Φ'^2-Φ'')=%s  d²_τJ=%s  residual=%s"
          % (mp.nstr(a, 6), mp.nstr(b, 6), mp.nstr(lhs, 18), mp.nstr(rhs, 18), mp.nstr(lhs - rhs, 3)))

print("[V3] 网格均值 (N=900, 周期梯形 ⇒ 常数项)")
N = 900
t = 2 * np.pi * np.arange(N) / N
T1, T2 = np.meshgrid(t, t, indexing="ij")
TH = [T1, T2, -(T1 + T2)]
S = {p: 4 * np.sin((TH[p[0]] - TH[p[1]]) / 2) ** 2 for p in PAIRS}
Jg = S[(0, 1)] * S[(0, 2)] * S[(1, 2)]


def prod2(ex):
    o = np.ones_like(T1)
    for p in PAIRS:
        if p != ex:
            o = o * S[p]
    return o


A = sum(mp.mpf(1) / 2 * 0 + (LAM[i] - LAM[j]) ** 2 * 2 * prod2((i, j)) for i, j in PAIRS)
sq = sum((LAM[i] - LAM[j]) ** 2 * (4 - S[(i, j)]) * prod2((i, j)) for i, j in PAIRS)
cr = 0.0
for x in range(3):
    for y in range(x + 1, 3):
        i, j = PAIRS[x]
        k, l = PAIRS[y]
        rem = [p for p in PAIRS if p not in (PAIRS[x], PAIRS[y])][0]
        cr = cr + 2 * (LAM[i] - LAM[j]) * (LAM[k] - LAM[l]) * \
            (2 * np.sin(TH[i] - TH[j])) * (2 * np.sin(TH[k] - TH[l])) * S[rem]
print("  mean |Δ|^2      =", Jg.mean())
print("  mean |Δ|^2 Φ''  =", A.mean())
print("  mean |Δ|^2 Φ'^2 =", (sq + cr).mean(), " = sq %s + cross %s" % (sq.mean(), cr.mean()))
print("  S = A/meanJ =", A.mean() / Jg.mean(), "  D =", (sq + cr).mean() / Jg.mean())
print("=" * 90)
