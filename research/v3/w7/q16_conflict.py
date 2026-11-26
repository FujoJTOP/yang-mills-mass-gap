# -*- coding: utf-8 -*-
"""Q16 — 判定常数冲突：P04-SHARPEN (C=1.7454 [T]) vs P17-AUDIT-VAR (C=10.05)。

只写 research/v3/w7/q16_*.py|.log。禁 git commit。
全部数值自立实现：
  u = arccos(1-w), phi(w) = sin u/u = sqrt(w(2-w))/u, dnu ∝ phi(w) e^{-Bw} dw on [0,2]
  Psi(B) = 36 B^2 Var(g) = 9 B^2 Var(phi)
分支：
  R  : range 支 (Popoviciu on range of g)
  L  : Lipschitz+tail 支 (18 L^2 + 9B^2 e^{-B w1})
  C  : p04 的 centre 支 (L4-L8 除式矩界)
"""
import math
import numpy as np
from scipy.integrate import quad

PI = math.pi
A_LIP_O20 = 9.0 * 3.327907          # 29.9512  (o20/p04 的矩常数支)
C1_O20 = 1.0 / (1.0 - math.exp(-1.0))  # 1.581977


def phi_w(w):
    if w <= 0.0:
        return 1.0
    if w >= 2.0:
        return 0.0
    u = math.acos(1.0 - w)
    return math.sin(u) / u


def g_u(u):
    if u <= 0.0:
        return 0.5
    return math.sin(u) / (2.0 * u)


def dphi_w(w):
    if w <= 0.0:
        return -1.0 / 3.0
    if w >= 2.0:
        return -math.inf
    a = math.sqrt(w * (2.0 - w))
    b = math.acos(1.0 - w)
    return ((1.0 - w) * b / a - 1.0) / (b * b)


def _q(f, a, b, lim=400):
    v, _ = quad(f, a, b, epsabs=1e-15, epsrel=1e-13, limit=lim)
    return v


def Z1(B):
    return _q(lambda w: phi_w(w) * math.exp(-B * w), 0.0, 2.0)


def Z_lb(B):
    v1 = _q(lambda w: (1.0 - w / 3.0 - w ** 2 / 45.0 - w ** 3 / 189.0 - 0.0035 * w ** 4)
            * math.exp(-B * w), 0.0, 1.0)
    v2 = _q(lambda w: (2.0 / PI) * (2.0 - w) * math.exp(-B * w), 1.0, 2.0)
    return v1 + v2


def Kbar(B, M):
    t1 = M * _q(lambda w: w * math.exp(-B * w), 0.0, 1.0)
    t2 = math.exp(-B) * _q(lambda w: w * abs(dphi_w(w)), 1.0, 2.0)
    return t1 + t2


def moments(B):
    zl = Z_lb(B)
    Ek = {}
    for k in (1, 2, 3, 4):
        num = _q(lambda w, k=k: (w ** k) * (1.0 - w / 3.0) * math.exp(-B * w), 0.0, 2.0)
        Ek[k] = num / zl
    return zl, Ek


def centre(B, M, extended):
    """p04 的 centre 支。extended=True 时把 W_lo 覆盖到 W_lo=max(0,(1-r)/B)（全 B 合法）。"""
    zl, Ek = moments(B)
    r = Kbar(B, M) / zl
    wlo = (1.0 - r) / B
    if extended:
        wlo = max(0.0, wlo)
    else:
        if r > 1.0:
            return float('inf'), r
    var = ((Ek[2] - wlo * wlo) / 9.0
           + (2.0 / 3.0) * (Ek[3] / 45.0 + 0.06 * Ek[4])
           + 0.01185 / B ** 4 + 0.1654 / B ** 5 + 0.6918 / B ** 6)
    return 9.0 * B * B * var, r


def Var_phi(B):
    """真值 Var_nu(phi)，稳定形式 h = 1-phi。"""
    Z = Z1(B)
    Eh = _q(lambda w: (1.0 - phi_w(w)) * math.exp(-B * w), 0.0, 2.0) / Z
    Eh2 = _q(lambda w: (1.0 - phi_w(w)) ** 2 * math.exp(-B * w), 0.0, 2.0) / Z
    return Eh2 - Eh * Eh


# w1 网格（Lipschitz 枝）
def _build_w1():
    head = np.array([1e-4])
    body = np.linspace(1e-3, 1.995, 4000)
    W1 = np.concatenate([head, body])
    L1 = np.array([abs(dphi_w(float(w))) for w in W1])
    return W1, L1


W1, L1 = _build_w1()


def lip_C1(B):
    """p04 §4.2 的 Lip+tail 支：18 L(w1)^2 + 9B^2 e^{-B w1}（C1=1）。"""
    return float(np.min(18.0 * L1 ** 2 + 9.0 * B * B * np.exp(-B * W1)))


def lip_o20(B):
    """o20/p04 §3 的 Lip+tail 支：A_lip L^2 + 9B^2 min(1, C1 e^{-B w1})。"""
    return float(np.min(A_LIP_O20 * L1 ** 2
                        + 9.0 * B * B * np.minimum(1.0, C1_O20 * np.exp(-B * W1))))


def scan(Bgrid, M, which):
    out = []
    for B in Bgrid:
        R = 2.25 * B * B          # 修正 range 支
        Ro = 0.2971 * B * B       # o20/p04 range 支
        Lb = lip_C1(B)
        Cc, r = centre(B, M, extended=(which == 'ext'))
        if which == 'o20rig':
            val = min(Ro, lip_o20(B))
        elif which == 'p04':
            val = min(Ro, Lb, Cc)
        else:
            val = min(R, Lb, Cc)
        out.append((B, val, Ro, R, Lb, Cc, r))
    return out


def head(title):
    print("\n" + "=" * 74)
    print(title)
    print("=" * 74)


def main():
    print("[Q0] repro / env")
    print("  phi(1e-6)=%.15f  dphi(0+)=%.12f  dphi(1)=%.10f  dphi(1.9)=%.6f"
          % (phi_w(1e-6), dphi_w(0.0), dphi_w(1.0), dphi_w(1.9)))

    # ---------------- Q1: g 的值域 ----------------
    head("[Q1] g(u)=sin u/(2u) 在 (0,pi) 上的值域（核 p17 的 range 主张）")
    us = np.linspace(1e-9, PI, 2000001)
    gs = np.sin(us) / (2.0 * us)
    print("  sup g = %.12f  at u=%.3e   (理论 0.5)" % (gs.max(), us[np.argmax(gs)]))
    print("  inf g = %.12f  at u=pi      (理论 0)" % gs.min())
    print("  g(pi/2)=%.12f   1/pi=%.12f   g(pi/2)-1/pi=%.6e"
          % (g_u(PI / 2), 1.0 / PI, g_u(PI / 2) - 1.0 / PI))
    print("  ⇒ g((0,pi)) 的值域 = (%.1f, %.1f) ；g((0,pi/2)) 的值域 = (%.6f, 0.5)"
          % (gs.min(), gs.max(), g_u(PI / 2)))
    # 值域长度的两种常数
    r_o20 = 36.0 * (0.5 - 1.0 / PI) ** 2 / 4.0
    r_corr = 36.0 * (0.5 - 0.0) ** 2 / 4.0
    print("  Popoviciu 常数: 用 (1/pi,1/2) -> 36*(1/2-1/pi)^2/4 = %.6f   (o20/p04 的 0.2971)" % r_o20)
    print("  Popoviciu 常数: 用 (0,1/2)    -> 36*(1/2)^2/4     = %.6f   (p17 的 2.25)" % r_corr)
    # 质量落在 (pi/2,pi) 的比例
    for B in (1e-6, 1e-3, 0.1, 1.0, 10.0):
        num = _q(lambda w: phi_w(w) * math.exp(-B * w), 1.0, 2.0)
        den = Z1(B)
        print("  nu_B 在 w>1 (即 u>pi/2) 上的质量: B=%-8.0e %.6f" % (B, num / den))

    # ---------------- Q2: L2 尾界 (C1=1) 数值核 ----------------
    head("[Q2] p04 的 L2：P(w>t) <= e^{-Bt}（C1=1）数值核")
    worst = 0.0
    for B in (0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 200.0):
        for t in np.linspace(0.02, 1.98, 25):
            S = _q(lambda w: phi_w(w) * math.exp(-B * w), t, 2.0)
            A = _q(lambda w: phi_w(w) * math.exp(-B * w), 0.0, t)
            if S + A <= 0:
                continue
            ratio = (S / (S + A)) / math.exp(-B * t)
            worst = max(worst, ratio)
    print("  max_{B,t} P(w>t)/e^{-Bt} = %.6f   (<=1 则 L2 成立，C1=1)" % worst)
    # 与 C1=1.582 的 o20 版本对照
    worst_o = 0.0
    for B in (0.3, 0.5, 1.0, 2.0, 5.0):
        for f in (1.0, 1.5, 2.0, 3.0, 5.0):
            t = f / B
            if t >= 2.0:
                continue
            S = _q(lambda w: phi_w(w) * math.exp(-B * w), t, 2.0)
            A = _q(lambda w: phi_w(w) * math.exp(-B * w), 0.0, t)
            worst_o = max(worst_o, (S / (S + A)) / math.exp(-B * t))
    print("  o20 版本（t=f/B 扫）max P/e^{-Bt} = %.6f" % worst_o)

    # ---------------- Q3: 两个矩常数 ----------------
    head("[Q3] 矩常数：3.3279/B^2 是否用到 g 的值域？")
    print("  p17 §4.2 复算: E[w^2] <= (1 + 4 C1/e)/B^2，C1=1/(1-e^-1)")
    print("     C1=%.6f -> (1+4C1/e) = %.6f  (o20 的 3.3279)" % (C1_O20, 1.0 + 4 * C1_O20 / math.e))
    print("     C1=1        -> (1+4/e)    = %.6f  （若在矩界里也用 L2 的 C1=1）" % (1.0 + 4.0 / math.e))
    print("  实测 E[w^2]B^2:")
    for B in (1.0, 2.0, 5.0, 20.0, 100.0):
        zl, Ek = moments(B)
        print("     B=%-6.1f E2*B^2 = %.6f" % (B, Ek[2] * B * B))
    print("  ⇒ 该常数只由 (截断 1/B)+(尾界 C1) 给出，与 sup/inf g 无关。")

    # ---------------- Q4: 三条支路，重算 p04 的 C ----------------
    head("[Q4] p04 §3 的 6.505 复算（grid 粗细之争）")
    Bg = np.exp(np.linspace(math.log(0.2), math.log(200.0), 4000))
    vals = [(B, min(0.2971 * B * B, lip_o20(B))) for B in Bg]
    v = np.array([x[1] for x in vals])
    i = int(np.argmax(v))
    print("  4000 点网格: C_rig(o20 range, C1=1.582) = %.4f  @ B=%.4f" % (v[i], vals[i][0]))
    for w1 in (1.0, 1.202, 1.224, 1.3):
        Lw = abs(dphi_w(w1))
        print("     w1=%.3f L=%.6f  A_lip L^2=%.4f" % (w1, Lw, A_LIP_O20 * Lw * Lw))
    print("  ⇒ p04 的 6.505 与 p17 的 6.696/6.70 的差 = 网格分辨率/装配细节，非 range。")

    head("[Q5] p04 的 C_centre（M=pi/6 全解析档）复算 + 支路表")
    Ms = {'M=pi/6 (p04 [T] 档)': PI / 6.0, 'M=0.406 (p04 [C] 档)': 0.406}
    res = {}
    for name, M in Ms.items():
        # p04 档：range 用 o20 的 0.2971，centre 只用 r<=1 的 B
        rows = scan(np.exp(np.linspace(math.log(0.4), math.log(12.0), 260)), M, 'p04')
        arr = np.array([r[1] for r in rows if np.isfinite(r[1])])
        vmax = float(np.max(arr))
        idx = int(np.argmax([r[1] if np.isfinite(r[1]) else -1 for r in rows]))
        Bs = rows[idx][0]
        res[name] = (vmax, Bs)
        print("  %s : max_B min(R_o20, L, C) = %.4f  @ B=%.3f   (p04 报 1.7454 / 1.6971)"
              % (name, vmax, Bs))
        print("     该 B 处: R_o20=%.4f  L=%.4f  C=%.4f  r=%.4f"
              % (rows[idx][2], rows[idx][4], rows[idx][5], rows[idx][6]))
        # r 何时越过 1（=centre 支的合法下界）
        bb = np.exp(np.linspace(math.log(0.2), math.log(6.0), 25))
        for B in bb:
            zl, Ek = moments(B)
            r = Kbar(B, M) / zl
            if r <= 1.0:
                print("      r<=1 起始: B≈%.3f (r=%.4f)  ⇒ p04 的 'B>5/3' 门槛" % (B, r))
                break

    head("[Q6] 修正 range 后：p04 支路的 uniform 结论还站得住吗")
    print("  真 range 支: Psi <= 2.25 B^2  ⇒  只到 B<=2/3 才 <=1 (p17 §4.4 对)")
    print("  p04 的 uniform 覆盖: B<=5/3 靠 range 支(R_o20<=0.825)，B>5/3 靠 centre 支")
    print("  ⇒ 修正后 B in (2/3, 5/3] 出现空档，除非把 centre 支扩到全 B (W_lo=max(0,...))")
    rows = scan(np.exp(np.linspace(math.log(0.6), math.log(6.0), 240)), PI / 6.0, 'ext')
    print("\n   B      min(R2.25,L,C_ext)   R(2.25B^2)   L         C_ext     r")
    for (B, val, Ro, R, Lb, Cc, r) in rows[::12]:
        print("  %6.3f   %12.4f      %9.4f  %8.4f  %8.4f  %6.3f" % (B, val, R, Lb, Cc, r))
    arr = np.array([r[1] for r in rows])
    i = int(np.argmax(arr))
    print("\n  max_B min(R_corr, L, C_ext) = %.4f  @ B=%.3f  (R=%.4f L=%.4f C=%.4f)"
          % (rows[i][1], rows[i][0], rows[i][3], rows[i][4], rows[i][5]))

    head("[Q7] 正确 C 的精化（两支交叉点附近加密）")
    fine = np.exp(np.linspace(math.log(1.0), math.log(3.0), 400))
    rows = scan(fine, PI / 6.0, 'ext')
    arr = np.array([r[1] for r in rows])
    i = int(np.argmax(arr))
    print("  C_correct(M=pi/6) = %.4f @ B=%.4f" % (rows[i][1], rows[i][0]))
    print("     该 B: R=%.4f  L=%.4f  C_ext=%.4f  r=%.4f" % (rows[i][3], rows[i][4], rows[i][5], rows[i][6]))
    fine2 = np.exp(np.linspace(math.log(1.0), math.log(3.0), 400))
    rows2 = scan(fine2, 0.406, 'ext')
    arr2 = np.array([r[1] for r in rows2])
    j = int(np.argmax(arr2))
    print("  C_correct(M=0.406) = %.4f @ B=%.4f   [C]" % (rows2[j][1], rows2[j][0]))

    # 完全不用 centre 支时的上限（= p17 的合成，修正 range）
    head("[Q8] p17 的合成（丢掉 p04 的 centre 支）复核")
    rows = scan(np.exp(np.linspace(math.log(0.5), math.log(50.0), 1200)), PI / 6.0, 'o20rig')
    # 这里 which='o20rig' 用的是 R_o20；再单独算修正版
    vals2 = [(B, min(2.25 * B * B, lip_C1(B))) for B in
             np.exp(np.linspace(math.log(0.5), math.log(50.0), 1200))]
    v2 = np.array([x[1] for x in vals2]); k = int(np.argmax(v2))
    print("  max_B min(R_corr=2.25B^2, L_C1) = %.4f @ B=%.4f   (p17 报 10.046，用的是 C1=1.582 与 A_lip 支)"
          % (v2[k], vals2[k][0]))
    vals3 = [(B, lip_C1(B)) for B in np.exp(np.linspace(math.log(0.5), math.log(6.0), 800))]
    v3 = np.array([x[1] for x in vals3]); m = int(np.argmax(v3))
    print("  max_B L_C1 (Lip 支单独) = %.4f @ B=%.4f   (p17 报 10.641 @1.481)" % (v3[m], vals3[m][0]))

    head("[Q9] 真值参考：Psi(B) = 9B^2 Var(phi)")
    for B in (0.5, 1.0, 1.5, 2.0, 2.43, 3.0, 5.0, 10.0, 50.0):
        print("  B=%-6.2f  Psi=%.8f   9B^2*1/16=%.3f  2.25B^2=%.3f"
              % (B, 9.0 * B * B * Var_phi(B), 9.0 * B * B / 16.0, 2.25 * B * B))

    print("\n[DONE]")


if __name__ == '__main__':
    main()
