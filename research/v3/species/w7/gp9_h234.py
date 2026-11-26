# -*- coding: utf-8 -*-
"""GP9-H234：H2/H3/H4 的精确算术核验（无 MC、无随机）。
只消费仓内精确整数 B_n 与 count_n（GV2:120-125），外加 ln/root 的浮点。
输出：research/v3/species/w7/gp9_h234.log
"""
from fractions import Fraction as F
from math import log, exp

out = []
def p(s=""):
    out.append(s)
    print(s)

# ---------- 0. 输入（GV2:120-127 逐位） ----------
B = [None, 1, 20, 514, 14748, 450255, 14305816]          # B_n = sum |mu|
CNT = [None, 1, 20, 458, 11132, 280175, 7218768]          # 簇数（每 plaquette）
p("== §0 输入 ==")
p("B_n            = " + str(B[1:]))
p("count_n        = " + str(CNT[1:]))
p("sum count_n    = %d   (GV2:126 = 7,510,554)" % sum(CNT[1:]))

# ---------- 1. H4 的三个常数 ----------
ln3 = log(3.0)
p()
p("== §1 H4 常数（精确到 1e-15） ==")
p("ln3                          = %.15f" % ln3)
p("gamma_ent + ln3      (g=1)   = %.15f   (GV2:91 '≈2.10')" % (1.0 + ln3))
p("0.5 + gamma_ent + ln3 (g=1)  = %.15f   (GV2/GO3 '2.5986')" % (0.5 + 1.0 + ln3))
p("round(.,4)                   = %.4f" % (0.5 + 1.0 + ln3))
p("rho at c0q = 2.598612       = %.9f   (GO3:138 '0.202179')" % exp(1.0 - (0.5 + 1.0 + ln3)))
p("1 - rho                      = %.9f   (GO3:138 '0.797821')" % (1.0 - exp(1.0 - (0.5 + 1.0 + ln3))))

# ---------- 2. O3 的 epsilon（GO3 建议值） ----------
p()
p("== §2 O3 'epsilon >= ?' （AS 单独关 CPB-W 所需） ==")
for n in (1, 2, 3, 4, 5, 6):
    eps = (1.0 + 0.81 ** n) / 2.0
    p("  n=%d : R_n >= 1-2rho = 2eps-1 >= 0.81^n = %.9f  =>  eps >= %.9f" % (n, 0.81 ** n, eps))
p("  => n=6 最紧: eps >= (1+0.81^6)/2 = %.15f  (GO3 '0.641215')" % ((1.0 + 0.81 ** 6) / 2.0))
p("  0.81^6 = %.15f" % (0.81 ** 6))

# ---------- 3. gamma_ent 的数值：可证下界 vs '≈1' ----------
p()
p("== §3 gamma_ent 的数值状态 ==")
p("  自由边界（每 plaquette）连通片集计数率，用仓内精确 count_n：")
for n in range(1, 7):
    p("    n=%d  count_n=%9d  (1/n)ln count = %.6f   (1/n)ln(count/6) = %.6f"
      % (n, CNT[n], log(CNT[n]) / n, log(CNT[n] / 6.0) / n))
p("  => 任何满足 N(n)<=e^{gamma n} (all n) 的 gamma 必须 >= (1/6)ln 7218768 = %.6f" % (log(CNT[6]) / 6.0))
p("     若按 'per site' 归一 (/6): gamma >= %.6f" % (log(CNT[6] / 6.0) / 6.0))
p("  壳质量 m_n = u^n B_n 的 n 次根（= kappa_c 表的口径）：")
for u, name in ((F(5, 18), "SU3 u=5/18"), (F(5, 8), "SU2 u=5/8")):
    row = []
    for n in range(1, 7):
        m = float(u) ** n * B[n]
        row.append("%.4f" % (m ** (1.0 / n)))
    p("    %s : m_n^{1/n} = %s" % (name, " ".join(row)))
    rowB = ["%.4f" % (B[n] ** (1.0 / n)) for n in range(1, 7)]
    p("    %s : B_n^{1/n} = %s" % (name, " ".join(rowB)))

# ---------- 4. H4 的 q-门槛（窗内自动性） ----------
p()
p("== §4 H4 自动性：c0 >= 0.5+g+ln3, c0*q >= g+ln3  =>  q >= (g+ln3)/(0.5+g+ln3) ==")
for g in (1.0, 1.5, 2.0, log(CNT[6] / 6.0) / 6.0, log(CNT[6]) / 6.0, log(64.0)):
    win = 0.5 + g + ln3
    thr = g + ln3
    p("  gamma=%8.6f : 窗 c0>=%.6f ; H4 阈值 c0q>=%.6f ; 自动 iff q >= %.9f"
      % (g, win, thr, thr / win))
p("  => gamma=1 时 q>=%.9f 才有自动性；q 是块单位质量，CPB-S 的窗在 q≈1（BL-5）" % ((1.0 + ln3) / (0.5 + 1.0 + ln3)))

# ---------- 5. AS 的 kappa_c 下界（常数 1/2 的代价） ----------
p()
p("== §5 AS 的 1/2 => kappa_c >= 2^{-1/|C0|}（GO3 O4） ==")
for n in (1, 2, 3, 4, 5, 6, 12):
    v = 2.0 ** (-1.0 / n)
    p("  |C0|=%2d : 2^{-1/|C0|} = %.9f   >=0.80? %s   >=0.81? %s"
      % (n, v, v >= 0.80, v >= 0.81))
p("  3^{-1/6} = %.9f  (rho<=1/3 的更紧版 1-rho/(1-rho)=1/2 之外，GO3 用过 (2/3)^{1/n})"
  % (3.0 ** (-1.0 / 6.0)))

# ---------- 6. 两条引用勘误的自核 ----------
p()
p("== §6 引用自核 ==")
p("  GV2:29 引 BianchiLemma.md:224 -> 该行是自编号行（见 §5 文字），c0>kappa 检索 0 命中")
p("  GV2:91 引 BianchiLemma.md:117 -> 实为 :111(公式)/:116(表行 '1.0 | 2.60')")
p("  gamma=1 表行算值 = %.4f (与 :116 的 2.60 一致)" % (0.5 + 1.0 + ln3))

with open(r"D:\Dev\math\research\v3\species\w7\gp9_h234.log", "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")
p()
p("[log written] research/v3/species/w7/gp9_h234.log")
