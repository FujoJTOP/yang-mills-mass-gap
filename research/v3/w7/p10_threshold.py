"""P10b — 判死线所需的三个量：单调性翻转角 theta_c(beta)、H1 判据 f=D+W 的箱内极值、
以及求积误差标定。同一模型/约定（见 p10_monotone.py docstring）。
"""
import numpy as np
import importlib.util, sys

spec = importlib.util.spec_from_file_location("p10", "p10_monotone.py")
# 直接复用其函数而不执行正文（避免重复打印）：重新定义最小实现
def grid(n):
    r = np.linspace(1e-9, 2 * np.pi - 1e-9, n)
    Sr, Cr = np.sin(r / 2), np.cos(r / 2)
    return r, Sr, Cr, Sr / r, Sr / (2 * r), 2 * np.pi * Sr**2 / r, 2.0 / r**2 - (Cr / Sr) / r


def _stats(g, beta, thetas):
    r, Sr, Cr, GG, HH, BASE, P2 = g
    c = np.cos(thetas / 2.0)[:, None]
    s = np.sin(thetas / 2.0)[:, None]
    S0 = 1.0 - c * Cr
    w = BASE * np.exp(-beta * S0)
    Z = np.trapezoid(w, r, axis=1)[:, None]
    P = w / Z
    Eg = np.trapezoid(P * GG, r, axis=1)
    Eh = np.trapezoid(P * HH, r, axis=1)
    Eg2 = np.trapezoid(P * GG**2, r, axis=1)
    Vg = Eg2 - Eg**2
    W = np.trapezoid(P * P2, r, axis=1)
    ES0 = np.trapezoid(P * S0, r, axis=1)
    Cg = np.trapezoid(P * GG * S0, r, axis=1) - Eg * ES0
    Cg2 = np.trapezoid(P * GG**2 * S0, r, axis=1) - Eg2 * ES0
    Ch = np.trapezoid(P * HH * S0, r, axis=1) - Eh * ES0
    Cp2 = np.trapezoid(P * P2 * S0, r, axis=1) - W * ES0
    c1, s1 = c[:, 0], s[:, 0]
    D = beta * c1 * Eh - beta**2 * s1**2 * Vg
    dD = c1 * (Eh + beta * (-Ch)) - s1**2 * (2 * beta * Vg + beta**2 * (-(Cg2 - 2 * Eg * Cg)))
    df = dD - Cp2                                        # d/db (D+W), dW/db = -Cov(P'',S0)
    return dict(D=D, dD=dD, W=W, f=D + W, df=df, Vg=Vg, Eh=Eh, c=c1, s=s1)


def one(g, beta, th):
    o = _stats(g, beta, np.array([th]))
    return {k: v[0] for k, v in o.items()}


gS, gF = grid(20001), grid(400001)

print("=" * 82)
print("[A] 求积误差标定: Delta_dyn 在 n=2e4 / 2e5 / 4e5 上的差")
g2, g3 = grid(200001), grid(800001)
for b in [1.0, 6.0]:
    for th in [0.0, np.pi / 2, np.pi]:
        a1 = one(gS, b, th)['D']
        a2 = one(g2, b, th)['D']
        a3 = one(gF, b, th)['D']
        a4 = one(g3, b, th)['D']
        print(f"  b={b:4g} theta={th / np.pi:.2f}pi  D: {a1:+.12f} {a2:+.12f} {a3:+.12f} {a4:+.12f}"
              f"   |n8e5-n4e5|={abs(a4 - a3):.2e}")

print("=" * 82)
print("[B] 单调性翻转角 theta_c(beta): dD/dbeta = 0 的 theta")
for b in [0.5, 1.0, 2.5, 4.0, 6.0, 10.0]:
    tg = np.linspace(0.5 * np.pi, np.pi, 4001)
    f = _stats(gS, b, tg)['dD']
    k = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]
    if len(k) == 0:
        print(f"  beta={b:4g}: 无翻转 (dD/db 端点 {f[0]:+.3e} / {f[-1]:+.3e})")
        continue
    a, z = tg[k[0]], tg[k[0] + 1]
    for _ in range(80):
        m = 0.5 * (a + z)
        if (one(gS, b, a)['dD']) * (one(gS, b, m)['dD']) <= 0:
            z = m
        else:
            a = m
    tc = 0.5 * (a + z)
    print(f"  beta={b:4g}: theta_c = {tc / np.pi:.6f} pi   (cos(theta_c/2) = {np.cos(tc / 2):.6f})"
          f"   D(theta_c) = {one(gS, b, tc)['D']:+.6f}")

print("=" * 82)
print("[C] H1 判据 f = D+W 在 beta in [2.5,6] x theta in [0,pi] 的箱内极值")
bg = np.linspace(2.5, 6.0, 141)
tg = np.linspace(0.0, np.pi, 181)
F = np.zeros((len(bg), len(tg)))
DF = np.zeros((len(bg), len(tg)))
for j, b in enumerate(bg):
    o = _stats(gS, b, tg)
    F[j] = o['f']
    DF[j] = o['df']
print(f"  min f = {F.min():+.6f} at beta={bg[np.unravel_index(F.argmin(), F.shape)[0]]:.3f}, "
      f"theta={tg[np.unravel_index(F.argmin(), F.shape)[1]] / np.pi:.3f}pi")
print(f"  max f = {F.max():+.6f} at beta={bg[np.unravel_index(F.argmax(), F.shape)[0]]:.3f}, "
      f"theta={tg[np.unravel_index(F.argmax(), F.shape)[1]] / np.pi:.3f}pi")
print(f"  df/dbeta > 0 比例 {np.mean(DF > 0):.4f};  f 在箱内每 theta 是否单调 beta: "
      f"{np.all((np.all(np.diff(F, axis=0) > 0, axis=0)) | (np.all(np.diff(F, axis=0) < 0, axis=0)))}")
bad = [tg[i] / np.pi for i in range(len(tg))
       if not (np.all(np.diff(F[:, i]) > 0) or np.all(np.diff(F[:, i]) < 0))]
if bad:
    print(f"  非 beta-单调的 theta/pi: min {min(bad):.3f}, max {max(bad):.3f}, count {len(bad)}")
print("  端点: f(2.5,pi)=%+.6f  f(6,pi)=%+.6f  f(2.5,pi/2)=%+.6f  f(6,pi/2)=%+.6f"
      % (one(gS, 2.5, np.pi)['f'], one(gS, 6.0, np.pi)['f'],
         one(gS, 2.5, np.pi / 2)['f'], one(gS, 6.0, np.pi / 2)['f']))

print("=" * 82)
print("[D] Delta_dyn 在 beta 上的内点极值: theta=0.9pi / pi")
for th in [0.85 * np.pi, 0.9 * np.pi, 0.95 * np.pi, np.pi]:
    bg2 = np.linspace(0.2, 14.0, 3000)
    D = np.array([one(gS, b, th)['D'] for b in bg2])
    k = int(np.argmax(D)); k2 = int(np.argmin(D))
    print(f"  theta={th / np.pi:.2f}pi : max D = {D[k]:+.6f} at beta={bg2[k]:.3f} ; "
          f"min D = {D[k2]:+.6f} at beta={bg2[k2]:.3f} ; D(2.5)={one(gS, 2.5, th)['D']:+.6f} "
          f"D(6)={one(gS, 6.0, th)['D']:+.6f}")
print("=" * 82)
