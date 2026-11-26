"""P2-CHAT 补充：按类分裂后 '通用熵就够用' 的折密度阈值 epsilon*，以及极值族的多项式计数。"""
import mpmath as mp
mp.mp.dps = 40
OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)

eD = mp.e * 20
N = 2
P("eΔ =", mp.nstr(eD, 16), " Cbar = N =", N, " C^_S = 4")
# 类 (E4/F = eps) 若仍用通用熵 eΔ：条件 4*0.625^eps*Cbar*eΔ*u < 1
for beta in (mp.mpf(1) / 48, mp.mpf(1) / 24, mp.mpf(1) / 12):
    u = beta / 8
    thr = 1 / (4 * N * eD * u)          # 需要 (5/8)^eps < thr
    eps = mp.log(thr) / mp.log(mp.mpf('0.625'))
    P(f" beta={mp.nstr(beta,8)} (u={mp.nstr(u,8)}): (5/8)^eps < {mp.nstr(thr,10)}"
      f"  => eps > {mp.nstr(eps,10)}  (即 E4/F 超过这个折密度才有通用熵可用)")
# 极值族（轴对齐 L^3 立方盒表面）过固定 plaquette 的个数：2 个第三方向 * 2 个面 * L^2 位置
P("")
P("轴对齐 L 立方盒表面过固定 plaquette 的个数 = 2*2*L^2 = 4L^2，F = 6L^2 => 个数 = (2/3)F（多项式）")
for L in (1, 2, 3, 4, 5):
    P(f"  L={L}: F={6*L*L:>4} 个数={4*L*L:>4}  (4u)^F 比值={mp.nstr((4*mp.mpf(1)/384)**(6*L*L),6)}")
P("")
P("== 需要的两条数（判死线的数值）==")
P("  Ĉ<4 的目标值 = 7.06328527049/N =", mp.nstr(mp.mpf('7.063285270491693') / 2, 12), "（从 4 砍 11.7%）")
P("  按类分裂：A_S < 1/beta =", 48, "(1/48) /", 24, "(1/24) /", 12, "(1/12)；通用 eΔ =", mp.nstr(eD, 10))
P("  ⇒ 1/48 只需砍", mp.nstr(eD / 48, 8), "x；2x 窗需砍", mp.nstr(eD / 24, 8), "x")
with open(__file__.replace(".py", ".log"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
print("[log written]")
