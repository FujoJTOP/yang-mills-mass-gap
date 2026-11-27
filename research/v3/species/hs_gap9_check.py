"""HS-GAP-9 price-surface check. Exact rationals only. No sampling, no MC."""
from fractions import Fraction as Fr
import math

def show(tag, *xs):
    print(tag, *[str(x) for x in xs])

# ---- 1. 价格面 w0 = F / Mmax  (F = 分子地板, Mmax = M0 天花板) --------------
for N in (2, 3):
    r0   = Fr(1, 4) * Fr(1, (2 * N) ** 2)     # ¼(2N)^{-2}   (Z01:204)
    perL = Fr(1, 2) * Fr(1, (2 * N) ** 2)     # ½(2N)^{-2}   (Z01:203, 逐 L)
    Mp   = Fr(N * N, 4)                       # Popoviciu: (b-a)^2/4, supp [0,N]
    Mc   = 2 * N * N                          # HZ1:47 粗界
    show(f"N={N}", f"r0={r0}", f"perL={perL}", f"perL/r0={perL/r0}")
    show(f"N={N}", f"w0_Pop(r0)={r0/Mp}", f"w0_Pop(perL)={perL/Mp}",
         f"w0_crude(r0)={r0/Mc}")
    # 要 w0 = 1/(4N^4) 所需的分子地板
    tgt = Fr(1, 4 * N ** 4)
    show(f"N={N}", f"target={tgt}", f"F_need(Mp)={tgt*Mp}", f"F_need(Mc)={tgt*Mc}",
         f"F_need(Mp)/r0={(tgt*Mp)/r0}", f"F_need(Mp)/perL={(tgt*Mp)/perL}")
    # 1/(4N^4) 是否等于 r0 （仅 N=2 处作为数）
    show(f"N={N}", f"target==r0 ? {tgt == r0}", f"target/r0={tgt/r0}")

# ---- 2. Popoviciu 与 Bhatia-Davis 的差 -----------------------------------
for N in (2, 3):
    Mp = Fr(N * N, 4)
    for m in (Fr(N, 10), Fr(N, 4), Fr(N, 2), Fr(4 * N, 5)):
        bd = (N - m) * m
        show(f"N={N} m={m}", f"BD={(N-m)*m}", f"Pop={Mp}",
             f"N^2/4-(m-N/2)^2={Mp - (m - Fr(N,2))**2}", f"gain={Mp/bd}")

# ---- 3. 买回 F 倍分母所需的均值偏离门槛 ----------------------------------
for gain in (Fr(3, 2), Fr(2), Fr(5, 2), Fr(4)):
    t = math.sqrt(1 - float(1 / gain))
    show(f"gain={gain}", f"need|2m/N-1|>={t:.6f}", f"need m/N<={(1-t)/2:.6f}")

# ---- 4. 窗内外距离 -------------------------------------------------------
show("dist", f"2.5/(4*3.777028e-3)={2.5/(4*0.003777028):.6f}")
show("neg", "r0=1/2*a0*(1-ubar) at a0=1/16:",
     Fr(1, 2) * Fr(1, 16) * (1 - Fr(1)), ">0 requires ubar<1")
