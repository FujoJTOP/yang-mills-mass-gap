# -*- coding: utf-8 -*-
"""Q16b — 高精度复核 Psi(B) 与 p04 centre 支的分项（mpmath 30 dps）。"""
import math
import mpmath as mp
import numpy as np
from scipy.integrate import quad

mp.mp.dps = 30
PI = mp.pi


def phi(w):
    if w <= 0:
        return mp.mpf(1)
    if w >= 2:
        return mp.mpf(0)
    u = mp.acos(1 - w)
    return mp.sin(u) / u


def Z_mp(B):
    B = mp.mpf(B)
    return mp.quad(lambda w: phi(w) * mp.e ** (-B * w), [0, 1, 2])


def Psi_mp(B):
    B = mp.mpf(B)
    Z = Z_mp(B)
    Eh = mp.quad(lambda w: (1 - phi(w)) * mp.e ** (-B * w), [0, 1, 2]) / Z
    Eh2 = mp.quad(lambda w: (1 - phi(w)) ** 2 * mp.e ** (-B * w), [0, 0.5, 1, 2]) / Z
    return 9 * B * B * (Eh2 - Eh * Eh), Eh, Eh2


def Psi_mp_u(B):
    """u 变量独立路：dnu ∝ sin^2 u / u e^{-B(1-cos u)} du"""
    B = mp.mpf(B)
    f = lambda u: mp.sin(u) ** 2 / u * mp.e ** (-B * (1 - mp.cos(u)))
    Z = mp.quad(f, [mp.mpf('1e-12'), mp.mpf(0.1), 1, 2, PI])
    Eg = mp.quad(lambda u: (mp.sin(u) / (2 * u)) * f(u), [mp.mpf('1e-12'), mp.mpf(0.1), 1, 2, PI]) / Z
    Eg2 = mp.quad(lambda u: (mp.sin(u) / (2 * u)) ** 2 * f(u), [mp.mpf('1e-12'), mp.mpf(0.1), 1, 2, PI]) / Z
    return 36 * B * B * (Eg2 - Eg * Eg)


print("[R1] Psi(B) 高精度（w 路 / u 路 + p17 的级数 1-2/(15B)-149/(315B^2)）")
print("   B        Psi_w(mp30)          Psi_u(mp30)          series        Psi_float_quad")
for B in (0.5, 1.0, 2.0, 2.43, 3.0, 5.0, 10.0, 50.0, 200.0):
    pw, Eh, Eh2 = Psi_mp(B)
    pu = Psi_mp_u(B)
    eps = mp.mpf(1) / B
    ser = 1 - mp.mpf(2) / 15 * eps - mp.mpf(149) / 315 * eps ** 2
    # float 路（与 q16_conflict.py 相同实现）
    def phif(w):
        if w <= 0: return 1.0
        if w >= 2: return 0.0
        u = math.acos(1.0 - w); return math.sin(u) / u
    Z = quad(lambda w: phif(w) * math.exp(-B * w), 0, 2, epsabs=1e-16, epsrel=1e-14, limit=800)[0]
    e1 = quad(lambda w: (1 - phif(w)) * math.exp(-B * w), 0, 2, epsabs=1e-16, epsrel=1e-14, limit=800)[0] / Z
    e2 = quad(lambda w: (1 - phif(w)) ** 2 * math.exp(-B * w), 0, 2, epsabs=1e-16, epsrel=1e-14, limit=800)[0] / Z
    pf = 9 * B * B * (e2 - e1 * e1)
    print("  %6.2f  %s  %s  %s  %12.8f"
          % (B, mp.nstr(pw, 12).ljust(20), mp.nstr(pu, 12).ljust(20), mp.nstr(ser, 12).ljust(12), pf))

print("\n[R2] p04 §6 的 centre 支分项 @ B=2.39（M=0.406）")
B = 2.39
M = 0.406
f1 = lambda w: (1 - w / 3 - w ** 2 / 45 - w ** 3 / 189 - 0.0035 * w ** 4) * math.exp(-B * w)
Z_lb = quad(f1, 0, 1, epsabs=1e-16, epsrel=1e-14)[0] + quad(
    lambda w: (2 / math.pi) * (2 - w) * math.exp(-B * w), 1, 2, epsabs=1e-16, epsrel=1e-14)[0]


def dphi(w):
    if w <= 0: return -1 / 3
    if w >= 2: return -1e18
    a = math.sqrt(w * (2 - w)); b = math.acos(1 - w)
    return ((1 - w) * b / a - 1) / (b * b)


Z1 = quad(lambda w: (lambda z: math.sqrt(z * (2 - z)) / math.acos(1 - z) if 0 < z < 2 else 0)(w)
          * math.exp(-B * w), 0, 2, epsabs=1e-16, epsrel=1e-14)[0]
Kbar = M * quad(lambda w: w * math.exp(-B * w), 0, 1, epsrel=1e-14)[0] + math.exp(-B) * quad(
    lambda w: w * abs(dphi(w)), 1, 2, epsrel=1e-13, limit=600)[0]
K = quad(lambda w: w * abs(dphi(w)) * math.exp(-B * w), 0, 2, epsrel=1e-13, limit=600)[0]
r = Kbar / Z_lb
Wlo = (1 - r) / B
Ek = {}
for k in (2, 3, 4):
    Ek[k] = quad(lambda w, k=k: w ** k * (1 - w / 3) * math.exp(-B * w), 0, 2, epsrel=1e-14)[0] / Z_lb
print("  Z1=%.8f  Z_lb=%.8f  Z1/Z_lb=%.6f" % (Z1, Z_lb, Z1 / Z_lb))
print("  K=%.8f  Kbar=%.8f  r=Kbar/Z_lb=%.6f  (K/Z1=%.6f)" % (K, Kbar, r, K / Z1))
print("  E[w]=%.6f  W_lo=(1-r)/B=%.6f  B*W_lo=%.6f" % ((1 / B) * (1 - K / Z1), Wlo, B * Wlo))
termA = B * B * (Ek[2] - Wlo ** 2)
termB = 6 * B * B * (Ek[3] / 45 + 0.06 * Ek[4])
termC = 9 * B * B * (0.01185 / B ** 4 + 0.1654 / B ** 5 + 0.6918 / B ** 6)
print("  E2=%.8f E3=%.8f E4=%.8f" % (Ek[2], Ek[3], Ek[4]))
print("  分项: A=B^2(E2-Wlo^2)=%.4f   B=(2/3)*(...)*9B^2=%.4f   C=E[e^2]*9B^2=%.4f  合计=%.4f"
      % (termA, termB, termC, termA + termB + termC))
print("  p04 §6 报: A≈1.05  B≈0.44  C≈0.21  合计 1.699")

print("\n[R3] 若 W_lo 用 '零'（B<=5/3 的合法替代）:")
for B in (1.2, 1.5, 1.667, 2.0, 2.39):
    Zl = Z_lb0 = None
    f = lambda w: (1 - w / 3 - w ** 2 / 45 - w ** 3 / 189 - 0.0035 * w ** 4) * math.exp(-B * w)
    Zl = quad(f, 0, 1, epsabs=1e-16, epsrel=1e-14)[0] + quad(
        lambda w: (2 / math.pi) * (2 - w) * math.exp(-B * w), 1, 2, epsabs=1e-16, epsrel=1e-14)[0]
    E2 = quad(lambda w: w * w * (1 - w / 3) * math.exp(-B * w), 0, 2, epsrel=1e-14)[0] / Zl
    E3 = quad(lambda w: w ** 3 * (1 - w / 3) * math.exp(-B * w), 0, 2, epsrel=1e-14)[0] / Zl
    E4 = quad(lambda w: w ** 4 * (1 - w / 3) * math.exp(-B * w), 0, 2, epsrel=1e-14)[0] / Zl
    Kb = M * quad(lambda w: w * math.exp(-B * w), 0, 1, epsrel=1e-14)[0] + math.exp(-B) * quad(
        lambda w: w * abs(dphi(w)), 1, 2, epsrel=1e-13, limit=600)[0]
    rr = Kb / Zl
    wlo = max(0.0, (1 - rr) / B)
    tot = 9 * B * B * ((E2 - wlo ** 2) / 9 + (2 / 3) * (E3 / 45 + 0.06 * E4)
                       + 0.01185 / B ** 4 + 0.1654 / B ** 5 + 0.6918 / B ** 6)
    print("  B=%-6.3f r=%.4f W_lo=%.4f  C_centre=%.4f   (A=%.4f B=%.4f)"
          % (B, rr, wlo, tot, B * B * (E2 - wlo ** 2),
             9 * B * B * (2 / 3) * (E3 / 45 + 0.06 * E4)))
print("\n[DONE]")
