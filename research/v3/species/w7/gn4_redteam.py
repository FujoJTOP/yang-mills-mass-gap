# -*- coding: utf-8 -*-
"""GN4-REDTEAM witness: all numbers quoted in research/v3/species/GN4-REDTEAM.md.

Pure stdlib (math only). No MC, no numpy, no I/O besides stdout.  <= 1 s.
py -3.12 research/v3/species/w7/gn4_redteam.py > research/v3/species/w7/gn4_redteam.log
"""
import math

sep = "-" * 78


def pr(tag, *a):
    print(tag, *a)


print("GN4-REDTEAM witness  (py -3.12, stdlib only)")
print(sep)

# ---------------------------------------------------------------- [0] ratios
print("[0] 口径：beta=2.5 相对 SU(2) 去禁闭 beta_c(N_tau)")
bc4, bc6 = 2.2986, 2.43
pr("  beta_c(N_tau=4) =", bc4, " beta=2.5 ->", round(2.5 / bc4, 4), "(+%0.2f%%)" % (100 * (2.5 / bc4 - 1)))
pr("  beta_c(N_tau=6) =", bc6, " beta=2.5 ->", round(2.5 / bc6, 4))
print("  (A6-TARGET.md:93 记 1.088 / +8.7%%；A4-SCALE-SEP.md:50 同)")

# ------------------------------------------------- [1] finite-T mass ratios
print(sep)
print("[1] 同一 beta 下换成 L_tau：m(N_tau=4)/m(N_tau=6)（= 同 a 下的有限温质量比）")
# research/unconditional/tm_gap_results.csv  (gevp48_tL, err)
rows = {  # beta: (ma_tau4, err4, ma_tau6, err6)
    2.2: (1.1506, 0.0406, 1.0800, 0.6001),
    2.5: (1.2562, 0.1008, 0.8467, 0.1481),
    2.8: (1.3975, 0.0666, 0.7153, 0.2257),
}
wsx = wsum = 0.0
for b in sorted(rows):
    m4, e4, m6, e6 = rows[b]
    r = m4 / m6
    s = r * math.hypot(e4 / m4, e6 / m6)
    lr, ls = math.log(r), math.hypot(e4 / m4, e6 / m6)
    pr("  beta=%.1f : r=%.4f +- %.4f   ln r=%.4f +- %.4f" % (b, r, s, lr, ls))
    wsx += lr / ls ** 2
    wsum += 1.0 / ls ** 2
lr = wsx / wsum
ls = 1.0 / math.sqrt(wsum)
pr("  weighted mean ln r = %.4f +- %.4f  ->  r = %.4f  (x/ %.4f)  z(below 1) = %.2f sigma"
   % (lr, ls, math.exp(lr), math.exp(ls), -lr / ls))
pr("  判据：凹陷要求 r<1（质量随 T 下降）。实测 3/3 点 r>1。")
print("  出处 tm_gap_results.csv:3,4,5,6,7,8 ; A6-TARGET.md:19-20 只引了 beta=2.5 两行。")

# 同一 N_tau 下 beta 斜率（C2-MEASURED-SLOPE.md:53-56 复核）
print(sep)
print("[2] 同一 N_tau 下 d ln(ma)/d beta（复核 C2-MEASURED-SLOPE.md:53-56）")
for tag, b1, m1, b2, m2 in (
        ("N_tau=4  2.20->2.50", 2.2, 1.1506, 2.5, 1.2562),
        ("N_tau=4  2.50->2.80", 2.5, 1.2562, 2.8, 1.3975),
        ("N_tau=6  2.20->2.50", 2.2, 1.0800, 2.5, 0.8467),
        ("N_tau=6  2.50->2.80", 2.5, 0.8467, 2.8, 0.7153)):
    pr("  %s : dln(ma)/dbeta = %+.3f" % (tag, math.log(m2 / m1) / (b2 - b1)))
dln_a = -2.523  # 两圈 dln(aL)/dbeta at beta=2.5（GW7-FIXEDA-SCAN.md:82）
pr("  物理质量斜率 dln(m)/dbeta = dln(ma)/dbeta - dln a/dbeta，取 dln a/dbeta = %.3f" % dln_a)
for tag, s in (("N_tau=4 2.20->2.50", 0.293), ("N_tau=6 2.50->2.80", -0.561)):
    pr("  %s : dln m/dbeta = %+.3f" % (tag, s - dln_a))
print("  => 有限温线上物理质量随 beta(=随 T 升) 上升，4/4 点为正。")

# ------------------------------------------------ [3] lattice vs continuum
print(sep)
print("[3] beta=2.5 (SU(2)) 格点隙 vs 两圈连续统外推")
ma_lat = 1.707           # RUN3.md:51 m_fit（S/N(2)=5.0<10，B04.md:169 限定）
ma_2loop = 0.023         # B04.md:64 两圈
pr("  ma(lat, 实测 m_fit) = %.3f   ma(2-loop) = %.3f   ratio = %.1fx  (格点【高于】外推)"
   % (ma_lat, ma_2loop, ma_lat / ma_2loop))
pr("  xi/a = 1/%.3f = %.3f  ->  关联长度 < 1 格距（不在标度区）" % (ma_lat, 1 / ma_lat))
pr("  凹陷要求 ratio < 1 且 a->0 时保持 < 1；此处方向相反。")

# ------------------------------------------------------- [4] 2D YM: exact
print(sep)
print("[4] 2D 纯 YM：H_a = -(1/a)log T_a 的 a-消去（精确，机器核）")
g2 = 1.0
C2 = 0.75          # C2(fund), SU(2)
sigma = 0.5 * g2 * C2
pr("  sigma = (g^2/2)C2(fund) = %.6f = 3g^2/8" % sigma)
worst = 0.0
for L in (1, 2, 3):
    for a in (1e-3, 1e-2, 1e-1, 1.0, 10.0):
        lam2 = math.exp(-(g2 * a / 2.0) * C2 * L)     # T_a 的第二本征值
        Delta_a = -(1.0 / a) * math.log(lam2)
        Delta_c = L * sigma
        worst = max(worst, abs(Delta_a - Delta_c))
pr("  max_{L in 1..3, a in 1e-3..10} |Delta_a - Delta_cont| = %.3e" % worst)
pr("  gap(L) = L*3g^2/8 : L=1,2,3 -> %.4f, %.4f, %.4f  (ATOM2.md:69-71 同数)"
   % tuple(L * sigma for L in (1, 2, 3)))
print("  => 2D 上 (NEW) 是【恒等】，不是凹陷实例。")

# --------------------------------------------- [5] free field: no O(a^2) dent
print(sep)
print("[5] 自由场：零模隙（naive 与 Wilson 均精确 = m；Wilson 项在 k=0 为零）")
m = 1.0
for a in (1.0, 0.5, 0.1, 0.01):
    for r in (0.0, 1.0):
        # k=0: 4/a^2 sum sin^2(0)=0; Wilson term 2r/a sum (1-cos 0)=0
        m_lat = math.hypot(m + 2 * r / a * 0.0, 0.0)
        assert abs(m_lat - m) < 1e-15
pr("  Delta_a - Delta_cont = 0 对 a in {1,0.5,0.1,0.01}, r in {0,1} 全部精确成立")
print(sep)
print("[6] 有限体积（min-max 单调性见证）：自由标量 E_1(L) = sqrt(m^2+(2pi/(La))^2) >= m")
for L in (2, 4, 8, 16, 64):
    E = math.hypot(m, 2 * math.pi / (L * 1.0))
    pr("  L=%3d : E_1 = %.6f   E_1 - m = %.6f" % (L, E, E - m))
print("  => 体积只能把底抬高；L->inf 才降到 m。同型给出 m_scr(T) >= m_gap。")

# ------------------------------- [7] REPORT:281 correction at the T level
print(sep)
print("[7] 复核 glue-b-blocking/REPORT.md:281 的严格性例子")
pr("  H=diag(1,3), H_a=diag(1-eps_a,3), eps_a->1/2:")
for eps in (0.0, 0.1, 0.4, 0.49, 0.4999):
    pr("    eps=%.4f : ||H_a-H|| = %.4f   Delta_a = %.4f" % (eps, eps, 1 - eps))
pr("  => ||H_a-H|| -> 1/2 != 0，故 ':281' 的 '||H_a-H||->0（强预解必成立）' 不成立。")
pr("  正确的例子在 T 层（= ':283' 的家族）：T_a=exp(-a H_a)，H_a=diag(1-eps_a,3)")
for a in (1.0, 1e-1, 1e-2, 1e-3):
    eps = 0.5 - a                      # eps_a -> 1/2
    lam = [math.exp(-a * (1 - eps)), math.exp(-a * 3.0)]
    pr("    a=%.0e : ||T_a-I|| = %.3e   Delta_a = %.6f" % (a, max(abs(x - 1) for x in lam), 1 - eps))
pr("  => T 层一致收敛、隙 -> 1/2 < 1；极限 T=I 无 H（Delta_cont=+inf）。")
print(sep)
print("[8] 固定物理体积线（GW7）斜率复核")
r25, r266 = 2.6075, 2.7904
pr("  dln(Rbar)/dbeta (2.50->2.66) = %+.4f   (GW7-FIXEDA-SCAN.md:33 记 +0.424)"
   % (math.log(r266 / r25) / (2.66 - 2.50)))
pr("  dln(Rbar)/dbeta (2.59->2.66) = %+.4f" % (math.log(2.7904 / 2.7537) / (2.66 - 2.59)))
pr("  物理预测 dln a/dbeta = %.3f -> 差距 %.1f sigma (GW7 :33)"
   % (dln_a, (0.424 + 2.523) / 0.119))
print(sep)
print("done. 无 QED；Clay = 0%。")
