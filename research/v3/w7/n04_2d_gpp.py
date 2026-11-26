"""n04_2d_gpp.py — N04: 2D YM 上的 G'', W_full, Δ_dyn, Δ₀²

链条（全部来自语料，逐字）：
    G''(0) = E[H|0] − Var(G|0) + W_full = Δ_dyn + W_full        (eq:split)
    W_geo  = E[Φ''] − Var(Φ'),  Φ = −log J,  dU = J dX          (eq:wgeo)
    2D 精确解 (Migdal/热核)：ρ_φ(φ) = K_t(φ) sin²φ,  J = sin²φ (class angle)
      ⇒ W_geo(φ) = −∂_φ² log sin²φ = 2/sin²φ       （精确，t-无关）
      ⇒ G''(φ)   = −∂_φ² log(K_t sin²φ)
      ⇒ Δ_dyn    = G'' − W_geo = −∂_φ² log K_t     （精确恒等，非拟合）
    ⇒ Δ₀² = (G'' − W_full)/a² = Δ_dyn/a²           （a² 只改标度，不改符号）

注：Δ_dyn 与 J 的约定无关（J 在 G'' 与 W_geo 中同时出现、相消），但 G'' 本身依赖
    粗变量坐标（这正是 T03 的 0 与 H1A:132 的 2.0101 的全部差别）。

全解析导数（无 diff、无采样）。用法:  py -3.12 n04_2d_gpp.py > n04_2d_gpp.log
"""
from mpmath import mp, mpf, sin, cos, exp, acos, sqrt, pi, nstr, fsum

mp.dps = 40
LINE = "-" * 100


def K_ders(t, phi, jmax=40):
    """(K, K', K'')，SU(2) 特征标级数的逐项解析导数。
       K  = Σ_j (2j+1) χ_j e^{−tj(j+1)},  χ_j = sin((2j+1)φ)/sinφ
       χ'  = (w cos(wφ) sinφ − sin(wφ) cosφ)/sin²φ
       χ'' = [(1−w²) sin(wφ) sin²φ − 2w cosφ cos(wφ) sinφ + 2cos²φ sin(wφ)]/sin³φ
    """
    s, c = sin(phi), cos(phi)
    s2, s3 = s * s, s * s * s
    K = K1 = K2 = mpf(0)
    for j in range(jmax):
        w = 2 * j + 1
        e = exp(-t * j * (j + 1))
        sw, cw = sin(w * phi), cos(w * phi)
        K += w * sw / s * e
        K1 += w * (w * cw * s - sw * c) / s2 * e
        K2 += w * (sw * (1 - w * w) * s2 - 2 * w * c * cw * s + 2 * c * c * sw) / s3 * e
    return K, K1, K2


def neg_d2log_K(t, phi):
    K, K1, K2 = K_ders(t, phi)
    return (K1 / K) ** 2 - K2 / K


def Gpp_phi(t, phi):
    """G''(φ) = −∂²log(K_t sin²φ) = −∂²logK_t + 2/sin²φ"""
    return neg_d2log_K(t, phi) + 2 / sin(phi) ** 2


def main():
    print(LINE)
    print("[N04] 2D YM: G'', W_full, Δ_dyn, Δ₀²   (mpmath dps=%d)" % mp.dps)
    print("[恒等] G'' = Δ_dyn + W_full (eq:split) ; W_full = W_geo + X ; 2D 上 X = 0")
    print("[2D]   ρ_φ = K_t(φ)sin²φ ⇒ W_geo = 2/sin²φ ; Δ_dyn = −∂_φ²log K_t")
    print(LINE)

    phis = [mpf("0.3"), mpf("0.6"), mpf("1.0"), mpf("1.5"), mpf("2.0"), mpf("2.6")]
    ts = [mpf("0.5"), mpf("1"), mpf("2"), mpf("5"), mpf("40")]
    H1A = {  # H1A-EXACT.md:124–130 的存档值（t × φ），用于独立交叉验证
        (0, 0): "26.5616", (0, 1): "9.9143", (0, 2): "6.4123", (0, 3): "5.4495",
        (0, 4): "5.4594", (0, 5): "7.8812",
        (1, 0): "24.5616", (1, 1): "7.9143", (1, 2): "4.4123", (1, 3): "3.4495",
        (1, 4): "3.4596", (1, 5): "6.1344",
        (2, 0): "23.5617", (2, 1): "6.9144", (2, 2): "3.4125", (2, 3): "2.4526",
        (2, 4): "2.5105", (2, 5): "6.5586",
        (3, 0): "22.9850", (3, 1): "6.3480", (3, 2): "2.8783", (3, 3): "2.0243",
        (3, 4): "2.3853", (3, 5): "7.4417",
        (4, 0): "22.9011", (4, 1): "6.2731", (4, 2): "2.8246", (4, 3): "2.0101",
        (4, 4): "2.4189", (4, 5): "7.5261",
    }

    # ================= [A] 逐点 + 与 H1A 存档对账 =================
    print("\n[A] 逐点精确值（t = g²A/2，A = 面积）。W_geo(φ)=2/sin²φ（t-无关）")
    print("    末列 = 与 H1A-EXACT.md:124–130 存档值之差（独立交叉验证）\n")
    hdr = "%-6s %-6s %-14s %-12s %-14s %-11s %-12s" % (
        "t", "phi", "G''", "W_geo", "Δ_dyn=G''-W", "残差", "H1A差")
    print(hdr); print("-" * len(hdr))
    worst = mpf(0); worst_h1a = mpf(0)
    for kt, t in enumerate(ts):
        for kp, phi in enumerate(phis):
            K, K1, K2 = K_ders(t, phi)
            Gpp = Gpp_phi(t, phi)
            W = 2 / sin(phi) ** 2
            D = Gpp - W
            res = abs(D - neg_d2log_K(t, phi)); worst = max(worst, res)
            ref = mpf(H1A[(kt, kp)])
            dh = abs(Gpp - ref); worst_h1a = max(worst_h1a, dh)
            print("%-6s %-6s %-14s %-12s %-14s %-11s %-12s" % (
                nstr(t, 4), nstr(phi, 4), nstr(Gpp, 9), nstr(W, 8),
                nstr(D, 9), nstr(res, 3), nstr(Gpp - ref, 4)))
        print("-" * len(hdr))
    print("max |Δ_dyn(G''−W) − (−∂²logK_t)| = %s   ⇒ 2D 上 G'' = Δ_dyn + W_full 精确" % nstr(worst, 3))
    print("max |G'' − H1A 存档|              = %s   ⇒ 与语料表逐位一致（< 1e−4，存档舍入）" % nstr(worst_h1a, 3))

    # ================= [B] Δ_dyn 的符号与幅度（= Δ₀²·a²） =================
    print("\n[B] Δ₀²·a² = Δ_dyn：φ 网格 (0,π) 128 格。Δ_dyn = −∂_φ²log K_t")
    print("%-7s %-16s %-16s %-9s %-16s %-12s" % (
        "t", "sup_φ Δ_dyn", "inf_φ Δ_dyn", "argmaxφ", "max|Δ_dyn|", "max/8e^{-2t}"))
    print("-" * 86)
    for t in [mpf("0.5"), mpf("1"), mpf("2"), mpf("5"), mpf("10"), mpf("20"), mpf("40")]:
        grid = [mpf(k) * pi / 128 for k in range(1, 128)]
        vals = [neg_d2log_K(t, p) for p in grid]
        imax = max(range(len(vals)), key=lambda i: vals[i])
        mx = max(abs(v) for v in vals)
        print("%-7s %-16s %-16s %-9s %-16s %-12s" % (
            nstr(t, 3), nstr(vals[imax], 6), nstr(min(vals), 6), nstr(grid[imax], 4),
            nstr(mx, 6), nstr(mx / (8 * exp(-2 * t)), 5)))
    print("⇒ Δ_dyn 随面积指数衰减（比值列 → const），且符号随 φ 变（sup>0、inf<0）")
    print("  ⇒ 2D 上 sup_φ Δ₀²·a² → 0⁺ 且 inf_φ Δ₀²·a² → 0⁻：既无正下界，也无定号。")

    print("\n[B2] 三个代表点 sign(Δ_dyn)（t=5 / t=40）")
    for t in [mpf("5"), mpf("40")]:
        for phi in [mpf("0.6"), mpf("1.2"), pi / 2, mpf("2.0"), mpf("2.8")]:
            v = neg_d2log_K(t, phi)
            print("   t=%-4s φ=%-9s Δ_dyn=%-16s sign=%s" % (
                nstr(t, 3), nstr(phi, 6), nstr(v, 8),
                "正" if v > 0 else ("零" if v == 0 else "负")))

    # ================= [C] G'' 的坐标依赖 =================
    print("\n[C] 同一测度、换粗变量坐标 ⇒ G'' 变（W_geo 与之同步变），这是 T03 与 H1A 的全部差别")
    print("    平坦/线性坐标 (class angle 的 cos 坐标, J_u=sqrt(1−u²)):")
    print("%-6s %-8s %-14s %-13s %-14s %-13s" % ("t", "φ", "G''_φ(J=sin²φ)", "W_geo(φ)", "G''_u(J=sinφ)", "W_geo(u)"))
    print("-" * 82)
    for t in [mpf("0.5"), mpf("2"), mpf("40")]:
        for phi in [mpf("0.6"), mpf("1.0"), mpf("2.0")]:
            u0 = cos(phi)
            s, c = sin(phi), cos(phi)
            K, K1, K2 = K_ders(t, phi)
            G_phi = (K1 / K) ** 2 - K2 / K + 2 / s ** 2
            W_phi = 2 / s ** 2
            # u 坐标：φ=acos u, dφ/du=−1/s, d²φ/du²=−u/s³
            dphi, d2phi = -1 / s, -u0 / s ** 3
            g1 = K1 / K + c / s
            g2 = K2 / K - (K1 / K) ** 2 - 1 / s ** 2
            G_u = -(g2 * dphi ** 2 + g1 * d2phi)
            h1, h2 = c / s, -1 / s ** 2
            W_u = -(h2 * dphi ** 2 + h1 * d2phi)
            print("%-6s %-8s %-14s %-13s %-14s %-13s" % (
                nstr(t, 3), nstr(phi, 4), nstr(G_phi, 8), nstr(W_phi, 8),
                nstr(G_u, 8), nstr(W_u, 8)))
    print("⇒ 同一 (t,φ)：G''_φ ≠ G''_u（差 O(1)~O(10)），W_geo 同步变。")
    print("  ⇒ '2D 的 G''=0'(T03, 平坦坐标/纯 Haar 极限) 与 '2D 的 δ=+2.0101'")
    print("    (H1A:132, class-angle 坐标 t→∞, 恰是 2/sin²(1.5)=2.01006) 差的是整整一块 W_full；")
    print("    H1A 的 2.0101 = W_geo，与动力学无关（t→∞ 时 Δ_dyn→0）。")

    # ================= [D] 地板 + m =================
    print("\n[D] Haar/Jacobian 地板：W_geo(φ)=2/sin²φ ≥ 2 ≫ N/12=1/6（SU(2)），t-无关")
    print("    最小点 φ=π/2：W_geo = %s ⇒ 地板成立（[T] eq:wgeo）" % nstr(2 / sin(pi / 2) ** 2, 10))
    print("\n    2D 精确谱隙 m=(g²/2)C₂(fund)=3g²/8，与面积 A 无关；")
    print("    Δ_dyn=O(e^{−2t})=O(e^{−g²A}) ⇒ A→∞ 时 Δ_dyn→0 而 m 不动（彻底解耦）")
    for g2 in [mpf("1"), mpf("0.5"), mpf("2")]:
        m = 3 * g2 / 8
        print("      g²=%-5s m=3g²/8=%-9s a·m(a=0.15fm)=%-9s max_φΔ_dyn(t=40)=%s" % (
            nstr(g2, 3), nstr(m, 8), nstr(m * mpf("0.15"), 8),
            nstr(max(abs(neg_d2log_K(mpf(40), mpf(k) * pi / 128)) for k in range(1, 128)), 5)))
    print("\n[E] Δ₀² = Δ_dyn/a² ：2D 上 sup→0⁺、inf→0⁻，符号与面积都使其不能当正下界：")
    print("    'Δ₀>0' 在 2D 上被 m=3g²/8>0 与 Δ₀²→0（且变号）同时否定。")
    print(LINE)


if __name__ == "__main__":
    main()
