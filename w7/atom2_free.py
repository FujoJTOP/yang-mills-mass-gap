# -*- coding: utf-8 -*-
"""atom2_free.py — 自由场上的"原子性"：
   O = :phi^2: (项目里 O=1/4 trF^2 的自由场类似物)
   (1) m=0  : 阈值 0      => inf supp nu = 0  => 无隙；nu 绝对连续 => 无原子
   (2) m>0  : 阈值 2m     => inf supp nu = 2m > 0 => 有隙，但仍无原子（第二个'有隙无原子'见证）
   (3) 有限 L: 每个二粒子能级都是原子（紧转移矩阵 => 谱纯离散）=> 任何'格点上 f1^2>0'
       的读数对自由场也报正 => 空洞（本脚本给 w_max/S 的 1/L 标度）
只写 w7/atom2_free.log。无 git commit。
"""
import mpmath as mp

mp.mp.dps = 25
PI = mp.pi

def levels(L, m):
    """环上 L 格点、格点色散 omega_k=sqrt(m^2+4 sin^2(k/2))；O=:phi^2: 的谱测度
    只在二粒子态上有权重（<0|:phi^2:|1particle>=0，Wick），能级 E=2 omega_k，
    权重 w_k = (1/(2 omega_k))^2  （归一无关，只用比值）。"""
    out = []
    for n in range(L):
        k = 2 * PI * n / L
        w = mp.sqrt(m ** 2 + 4 * mp.sin(k / 2) ** 2)
        if w == 0:
            continue          # m=0 的 k=0 零模：:phi^2: 上是 c-数平移，不进连通测度
        out.append((2 * w, (1 / (2 * w)) ** 2))
    return out

print("=" * 78)
print("[1] 有限 L：'格点上的原子'是空洞的（紧 T => 谱纯离散 => 处处是原子）")
print("    m=1.0 （有质量，有隙）")
print("      L      E_min=2m(数值)     w_max/S          #levels")
for L in (8, 16, 32, 64, 128, 256, 512):
    lv = levels(L, mp.mpf(1))
    Emin = min(e for e, w in lv)
    S = sum(w for e, w in lv)
    wmax = max(w for e, w in lv)
    print("    %5d   %.12f   %.6e   %d" % (L, float(Emin), float(wmax / S), len(lv)))
print("    => w_max/S ~ 1/L -> 0：有限 L 上 f1^2>0 恒真，但极限里没有原子")
print("    => 任何'格点测得 f1^2>0'都不构成原子证据（对自由场同样报正）")

print("-" * 78)
print("[2] 无质量 m=0（酸测靶：必须什么都不给）")
for L in (8, 32, 128):
    lv = levels(L, mp.mpf(0))
    Emin = min(e for e, w in lv)
    S = sum(w for e, w in lv)
    wmax = max(w for e, w in lv)
    print("    L=%4d : E_min(L) = %.10f  (L->inf 时 ->0 = 阈值)  w_max/S = %.3e"
          % (L, float(Emin), float(wmax / S)))
print("    => inf supp nu = 0 => 无隙；nu 绝对连续 => nu({E})=0 对每点 E  => 无隙无原子")

print("-" * 78)
print("[3] 判词（自由场行）")
print("    m=0 : 隙 = 0           原子 f1^2 = 0   => 无隙无原子   [T]")
print("    m>0 : 隙 = m > 0       原子 f1^2 = 0   => 有隙无原子   [T]  <- 第二个见证")
print("         理由：<0|:phi^2:|1particle> = 0（Wick），谱从二粒子阈 2m 起且为割")
print("    => '有隙无原子' 不是 2D YM 的特产；2D 的独特点是 nu 恒等于 0（算符=守恒荷）")
print("=" * 78)
