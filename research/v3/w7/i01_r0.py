# i01_r0.py — r₀(ū) 的精确式 / 变号点 / β=2.5 的位置。只算，不跑格点。
# 复现：cd D:\Dev\math; py -3.12 research\v3\w7\i01_r0.py
import mpmath as mp
mp.mp.dps = 50

LOG = []
def p(s=""):
    LOG.append(str(s)); print(s)

N   = mp.mpf(2)
a0  = 1/(2*N)**2          # (2N)^{-2} = 1/16   [T]  P5-R0:12,17
eD  = 20*mp.e             # eΔ = 20e = 54.3656365691809...   [T]  ADJ2:43
eD6 = mp.mpf('54.3656')   # 语料 6 位写法（P5-R0:35）
b1  = 2*N**2/eD           # 8/ē
b16 = 2*N**2/eD6

def u(b):   return b/(2*N**2)
def ub(b, e=eD): return e*u(b)
def r0(b, e=eD, half=True):
    c = mp.mpf(1)/2 if half else mp.mpf(1)
    return c*a0*(1-ub(b, e))

p("§0 常数（精确）")
p("  N                 = %s" % N)
p("  a₀ = (2N)^{-2}    = %s   [1/16]" % mp.nstr(a0, 20))
p("  ½a₀ = 1/32        = %s" % mp.nstr(a0/2, 20))
p("  ¼a₀ = 1/64        = %s" % mp.nstr(a0/4, 20))
p("  ē = eΔ = 20e      = %s" % mp.nstr(eD, 28))
p("  β₁ = 2N²/ē = 8/ē  = %s" % mp.nstr(b1, 28))
p("  同式用 6 位 ē      = %s   （P5-R0:13 记 0.147151875451）" % mp.nstr(b16, 28))
p("  β₁/(1/48)         = %s   （= 7.0633 语境）" % mp.nstr(b1*48, 20))
p("  β*(ū=1/2) = 4/ē   = %s   （P5-R0:38 的 0.073576）" % mp.nstr(2*N**2/eD/2, 28))

p()
p("§1 恒等式：ū ≡ β/β₁（不是近似，是同一式两写）")
for s in ('2.5',):
    b = mp.mpf(s)
    p("  β=2.5: ū=ēβ/8 = %s ; β/β₁ = %s ; 差 = %s"
      % (mp.nstr(ub(b), 22), mp.nstr(b/b1, 22), mp.nstr(ub(b)-b/b1, 5)))

p()
p("§2 ★ 变号点：解 1−ū = 0")
p("  ū* = 1  （系数 K=1，逐字 P5-R0:37 的 K=1 行）")
root = mp.findroot(lambda b: 1-ub(b), mp.mpf('0.147'))
p("  findroot ⇒ β* = %s ; 残差 1−ū = %s" % (mp.nstr(root, 25), mp.nstr(1-ub(root), 5)))
p("  β* = β₁ = 2N²/ē 逐位同值 ⇒ 变号点 = 窗边本身，不是独立数")
p("  第二窗边（D_∞≥½ ⇒ ū≤½）：β ≤ 4/ē = %s" % mp.nstr(4/eD, 25))

p()
p("§3 表（r₀ 两种前系数：½a₀ 与 a₀）")
p("  %-14s %-14s %-14s %-16s %-16s" % ("β", "ū=β/β₁", "1−ū", "r₀=½a₀(1−ū)", "r₀=a₀(1−ū)"))
for b in (mp.mpf(0), mp.mpf(1)/48, 4/eD, b1, mp.mpf(1)/12, mp.mpf(1), mp.mpf('2.2'), mp.mpf('2.5')):
    p("  %-14s %-14s %-14s %-16s %-16s"
      % (mp.nstr(b, 12), mp.nstr(ub(b), 10), mp.nstr(1-ub(b), 10),
         mp.nstr(r0(b), 10), mp.nstr(r0(b, half=False), 10)))

p()
p("§4 β=2.5 逐数")
b = mp.mpf('2.5')
p("  ū(2.5)        = 6.25e = %s" % mp.nstr(ub(b), 22))
p("  1−ū           = %s" % mp.nstr(1-ub(b), 22))
p("  r₀(2.5)       = ½a₀(1−ū) = %s" % mp.nstr(r0(b), 22))
p("  r₀(2.5)/(½a₀) = %s  （纯 1−ū）" % mp.nstr(1-ub(b), 18))
p("  a₀→4096·a₀ : r₀ = %s  （P5-R0:43 逐字 −2046.6）" % mp.nstr(r0(b)*4096, 12))
p("  检查 8/54.3656 口径: r₀ = %s" % mp.nstr(r0(b, e=eD6), 12))

p()
p("§5 三个 r₀ 口径在 β=2.5 上的值（互不等价）")
p("  C1  r₀ = ½a₀(1−O(ū))            = %s  （F05:53,65 · U15:34 · UNIF2:82）" % mp.nstr(r0(b), 12))
p("  C2  r₀ = ¼(2N)^{-2}（常数）      = %s  （Z01:204 · X02:111；其推导要 ū≤½）" % mp.nstr(a0/4, 12))
p("  C3  r₀ = ½c₁u^p（c₁=2,p=1 ⇒ u）  = %s  （PROOF-P3:31；u=β/8）" % mp.nstr(mp.mpf(1)*u(b), 12))

p()
p("§6 u 的三种归一 ⇒ ū(2.5) 三个值（都 ≫1）")
for name, uu in (("u=β/(2N²)=β/8（species）", mp.mpf('0.3125')),
                 ("u=β/N²=β/4（SU(2) 精确 c_f/N）", mp.mpf('0.625')),
                 ("u=1.014（T1:164 的 c_f 口径）", mp.mpf('1.014'))):
    p("  %-34s ū=ēu=%s ; 1−ū=%s ; r₀=½a₀(1−ū)=%s"
      % (name, mp.nstr(eD*uu, 12), mp.nstr(1-eD*uu, 12), mp.nstr(a0/2*(1-eD*uu), 12)))

open(r"D:\Dev\math\research\v3\w7\i01_r0.log", "w", encoding="utf-8").write("\n".join(LOG) + "\n")
print("\n[log written]")
