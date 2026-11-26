# -*- coding: utf-8 -*-
"""atom2_m6.py — M6 (谐振子 ⊗ 环上自由模, d=0+1) 与 单 plaquette：
   M6  : m=0（b-零模，C->inf）, G''=2w>0, 而 O=a 的 KL 测度 = 单个原子 (1/(2w)) at E1=w
         => 有原子无隙；且 I1 ∧ I2 在该模型上真而 m=0 => '原子性' 不是隙的充分条件
   单 plaquette : 无时间方向 => 无 T => 无谱测度 => 原子性 undefined
只写 w7/atom2_m6.log。无 git commit。mpmath（Hermite 求积 + Bessel）。
"""
import mpmath as mp

mp.mp.dps = 30
w = mp.mpf(1)          # 振子频率（w=1）

print("=" * 78)
print("[1] M6 的 a-边缘：rho_a(a)=sqrt(w/pi) e^{-w a^2}  =>  G''(0)=2w = %.6f [T]" % float(2 * w))
print("    谱: spec(H) = {(n+1/2)w} + {2 pi^2 k^2/C^2}  => m_C = min(w, (1/2)(2pi/C)^2) -> 0")

def psi0(x):
    return (w / mp.pi) ** mp.mpf(0.25) * mp.e ** (-w * x ** 2 / 2)
def psi(n, x):
    # 归一化谐振子本征函数
    c = (w / mp.pi) ** mp.mpf(0.25) / mp.sqrt(2 ** n * mp.factorial(n))
    return c * mp.hermite(n, mp.sqrt(w) * x) * mp.e ** (-w * x ** 2 / 2)

print("-" * 78)
print("[2] O = a 的 KL 测度：w_n = |<0|a|n>|^2  (数值求积，区间 [-20,20])")
tot = mp.mpf(0)
for n in range(0, 7):
    me = mp.quad(lambda x, n=n: psi0(x) * x * psi(n, x), [-20, 0, 20])
    tot += abs(me) ** 2
    print("    n=%d : <0|a|n> = %+.12f   |.|^2 = %.12f" % (n, float(me), float(abs(me) ** 2)))
print("    sum_n w_n = %.15f    精确值 1/(2w) = %.15f" % (float(tot), float(1 / (2 * w))))
print("    => nu_a = (1/(2w)) delta_{E1},  E1 = w = 1 ；单原子（宇称：a 只连 n=1）")
print("    => Lambda* = inf supp(nu - f1^2 delta) = +inf  => I2 平凡成立；I1 成立 (f1^2=1/2>0)")

print("-" * 78)
print("[3] 但 m = 0（b-模无力，周长 C -> inf）：")
for C in (2, 10, 40, 1000, mp.mpf('1e6')):
    mC = min(w, (mp.mpf(1) / 2) * (2 * mp.pi / C) ** 2)
    print("    C=%-10s : m_C = %.12e" % (str(C), float(mC)))
print("    => 原子上限与 C 无关（1/(2w)），隙 -> 0  => M6 = 有原子无隙")

print("-" * 78)
print("[4] 单 plaquette：无时间方向 => 无转移矩阵 => 原子性 undefined")
print("    唯一可算的是特征标系数 c_n(beta)=2(n+1) I_{n+1}(beta)/beta ; Z=c_0")
for beta in (mp.mpf('0.1'), mp.mpf(1), mp.mpf('2.5'), mp.mpf(20)):
    c0 = 2 * mp.besseli(1, beta) / beta
    c1 = 4 * mp.besseli(2, beta) / beta
    c2 = 6 * mp.besseli(3, beta) / beta
    print("    beta=%-5s : -log(c1/c0)=%+.6f   -log(c2/c0)=%+.6f   (两者之差非常数 => 无'谱')"
          % (str(beta), float(-mp.log(c1 / c0)), float(-mp.log(c2 / c0))))
print("    => 这些比值随 beta 单调趋于 0（beta->inf），无下界、无常数差 => 不是任何")
print("       算符的谱；单 plaquette 没有 T，nu 无定义 => 原子性 undefined [T]")

print("-" * 78)
print("[5] 四模型表（本脚本负责 M6 与 单 plaquette 两行）")
print("    模型         m          f1^2                    原子性")
print("    自由场(m=0)  0          0                       无隙无原子   [T]")
print("    自由场(m>0)  m>0        0                       有隙无原子   [T]")
print("    2D YM        3g^2/8>0   0 (nu==0, O=1/4trF^2)   有隙无原子   [T]")
print("    M6           0          1/(2w)=0.5 (E1=w)       有原子无隙   [E]/[T]")
print("    单 plaquette 无定义     无定义                  无定义       [T]")
print("=" * 78)
