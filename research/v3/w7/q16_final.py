# -*- coding: utf-8 -*-
"""Q16c — 修复 Psi 的权重 bug 后的复核 + 正确 C 与误差分解。"""
import math
import numpy as np
from scipy.integrate import quad

PI = math.pi


def phi(w):
    if w <= 0:
        return 1.0
    if w >= 2:
        return 0.0
    return math.sin(math.acos(1.0 - w)) / math.acos(1.0 - w)


def dphi(w):
    if w <= 0:
        return -1.0 / 3.0
    if w >= 2:
        return -1e18
    a = math.sqrt(w * (2 - w)); b = math.acos(1 - w)
    return ((1 - w) * b / a - 1.0) / (b * b)


def Z1(B):
    return quad(lambda w: phi(w) * math.exp(-B * w), 0, 2, epsabs=1e-17, epsrel=1e-14, limit=600)[0]


def Psi(B):
    """Psi = 9B^2 Var_nu(phi), nu ∝ phi(w)e^{-Bw}dw（这次带权重）"""
    Z = Z1(B)
    Eh = quad(lambda w: (1 - phi(w)) * phi(w) * math.exp(-B * w), 0, 2, epsabs=1e-18, epsrel=1e-14,
              limit=600)[0] / Z
    Eh2 = quad(lambda w: (1 - phi(w)) ** 2 * phi(w) * math.exp(-B * w), 0, 2, epsabs=1e-18,
               epsrel=1e-14, limit=600)[0] / Z
    return 9 * B * B * (Eh2 - Eh * Eh)


def Z_lb(B):
    return (quad(lambda w: (1 - w / 3 - w ** 2 / 45 - w ** 3 / 189 - 0.0035 * w ** 4) * math.exp(-B * w),
                 0, 1, epsabs=1e-17, epsrel=1e-14)[0]
            + quad(lambda w: (2 / PI) * (2 - w) * math.exp(-B * w), 1, 2, epsabs=1e-17, epsrel=1e-14)[0])


def Kbar(B, M):
    return (M * quad(lambda w: w * math.exp(-B * w), 0, 1, epsrel=1e-14)[0]
            + math.exp(-B) * quad(lambda w: w * abs(dphi(w)), 1, 2, epsrel=1e-13, limit=600)[0])


def pieces(B, M):
    zl = Z_lb(B)
    r = Kbar(B, M) / zl
    wlo = max(0.0, (1 - r) / B)
    E = {}
    for k in (2, 3, 4):
        E[k] = quad(lambda w, k=k: w ** k * (1 - w / 3) * math.exp(-B * w), 0, 2,
                    epsabs=1e-17, epsrel=1e-14)[0] / zl
    A = B * B * (E[2] - wlo ** 2)
    Bt = 6 * B * B * (E[3] / 45.0 + 0.06 * E[4])
    Ct = 9 * B * B * (0.01185 / B ** 4 + 0.1654 / B ** 5 + 0.6918 / B ** 6)
    return A, Bt, Ct, A + Bt + Ct, r, wlo, E[2], E[3], E[4], zl


W1 = np.concatenate([[1e-4], np.linspace(1e-3, 1.995, 4000)])
L1 = np.array([abs(dphi(float(w))) for w in W1])


def lip(B):
    return float(np.min(18.0 * L1 ** 2 + 9.0 * B * B * np.exp(-B * W1)))


print("[S1] Psi(B)（权重已修：phi*e^{-Bw}dw）")
for B in (0.5, 1.0, 1.33, 1.5, 2.0, 2.43, 3.0, 5.0, 10.0, 50.0):
    print("   B=%-6.2f  Psi=%.8f   (p17 级数 1-2/(15B)=%.6f)" % (B, Psi(B), 1 - 2 / (15 * B)))

print("\n[S2] centre 支分项")
for B in (1.33, 2.39, 2.43):
    M = PI / 6
    A, Bt, Ct, tot, r, wlo, E2, E3, E4, zl = pieces(B, M)
    print("  B=%.2f  A=B^2(E2-Wlo^2)=%.4f  B_e=%.4f  E[e^2]=%.4f  合计=%.4f  (r=%.3f Wlo=%.4f E2=%.4f E3=%.4f E4=%.4f)"
          % (B, A, Bt, Ct, tot, r, wlo, E2, E3, E4))

print("\n[S3] 正确 C：三条支路 min，range 用 p17 的 2.25B^2")
print("   B      R=2.25B^2    L        C_centre   min")
best = (0, 0)
for B in np.exp(np.linspace(math.log(0.55), math.log(6.0), 46)):
    R = 2.25 * B * B; Lb = lip(B); C = pieces(B, PI / 6)[3]
    v = min(R, Lb, C)
    if v > best[0]:
        best = (v, B)
    print("  %5.3f  %9.4f  %8.4f  %9.4f  %8.4f" % (B, R, Lb, C, v))
print("  ⇒ max_B min = %.4f @ B=%.3f" % best)

print("\n[S4] 精化 crossing")
fine = np.exp(np.linspace(math.log(1.15), math.log(1.55), 41))
bb = [(B, min(2.25 * B * B, lip(B), pieces(B, PI / 6)[3])) for B in fine]
bb.sort(key=lambda z: -z[1])
print("  top5:", ["B=%.4f min=%.4f" % (b, v) for b, v in bb[:5]])
B = bb[0][0]
M2 = 0.406
print("  同点 M=0.406: C_centre=%.4f" % pieces(B, M2)[3])

print("\n[S5] r<=1 的下界（centre 支合法性）")
lo, hi = 0.05, 3.0
for _ in range(40):
    mid = 0.5 * (lo + hi)
    if pieces(mid, PI / 6)[4] > 1.0:
        lo = mid
    else:
        hi = mid
print("  r=1 于 B≈%.4f  ⇒ centre 支在 B>%.3f 合法（M=pi/6）" % (hi, hi))

print("\n[S6] p04 的 1.7454 是 range 支与 centre 支的交叉点？")
for B in (2.39, 2.43, 2.47):
    R_o20 = 0.2971 * B * B
    C = pieces(B, PI / 6)[3]
    print("  B=%.2f  R_o20=0.2971B^2=%.4f   C_centre=%.4f   min=%.4f" % (B, R_o20, C, min(R_o20, C)))
print("\n[DONE]")
