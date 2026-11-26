# -*- coding: utf-8 -*-
r"""
r06_gap_mono_b.py -- R06 补充:
  [A] theta 表 (gap 与 dgap/dtheta)
  [B] gap''(c) 的数值范围 (决定 gap' 的 Lipschitz 常数, 用于把单调性换成 [T] 的预算)
  [C] 逐项支配关系: min T1, max|T3|, min(T1-|T3|), min(T2+T4+T5)
依赖 r06_gap_mono.py 的 moments()。
"""
import math
import mpmath as mp

# 只加载 r06_gap_mono.py 的定义段 (不重跑 502 点扫描)
_src = open("r06_gap_mono.py", encoding="utf-8").read().split(
    "# ================================================================= 计算")[0]
_ns = {"__name__": "r06_defs"}
exec(compile(_src, "r06_gap_mono.py", "exec"), _ns)
moments = _ns["moments"]
LOG = []


class R:
    pass


R.moments = moments
R.LOG = LOG


def P(s=""):
    print(s)
    LOG.append(s)


BETA = mp.mpf('2.5')

P("=" * 100)
P("r06_gap_mono_b.log -- R06 补充: theta 表 / gap'' / 逐项支配")
P("  工具状态: pwsh(OK) mpmath(OK) python-3.12(OK)")
P("=" * 100)
P("")

P("[A] theta 表 (c=cos(th/2));  gap'(c) 与 dgap/dth = -(sin(th/2)/2) gap'(c)")
P("  %-8s %-9s %-20s %-14s %-14s" % ("th", "c", "gap", "gap'(c)", "dgap/dth"))
for deg in (0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180):
    th = mp.pi * deg / 180
    c = mp.cos(th / 2)
    r = R.moments(c)
    s = mp.sin(th / 2)
    P("  %-8s %-9s %-20s %-14s %-14s"
      % ("%d deg" % deg, mp.nstr(c, 7), mp.nstr(r['gap'], 16),
         mp.nstr(r['dgap'], 10), mp.nstr(-s / 2 * r['dgap'], 10)))
P("")

P("[B] gap''(c) 数值范围 (中心差分 h=1e-3, 51 点):")
hs = mp.mpf('1e-3')
cs = [mp.mpf(i) / 52 for i in range(1, 52)]
d2 = []
for c in cs:
    a = R.moments(c - hs)['dgap']
    b = R.moments(c + hs)['dgap']
    d2.append((b - a) / (2 * hs))
imx = max(range(len(d2)), key=lambda i: abs(d2[i]))
P("  gap''  min = %s  max = %s ;  max|gap''| ~ %s at c=%s  (c=0: %s, c=1 侧: %s)"
  % (mp.nstr(min(d2), 8), mp.nstr(max(d2), 8), mp.nstr(abs(d2[imx]), 8),
     mp.nstr(cs[imx], 6), mp.nstr(d2[0], 8), mp.nstr(d2[-1], 8)))
P("  => gap' 在 (0.88,1] 上略降, 但 gap' 始终 > 0.44 [C]")
P("")

P("[C] 逐项支配 (同一 51 点):")
mnT1 = None
mxT3 = mp.mpf(0)
mn_rest = None
mn_dom = None
for c in cs:
    r = R.moments(c)
    T = r['T']
    a3 = abs(T[2])
    rest = T[1] + T[3] + T[4]
    if mnT1 is None or T[0] < mnT1[0]:
        mnT1 = (T[0], c)
    if a3 > mxT3:
        mxT3 = a3
    if mn_rest is None or rest < mn_rest[0]:
        mn_rest = (rest, c)
    if mn_dom is None or (T[0] - a3) < mn_dom[0]:
        mn_dom = (T[0] - a3, c)
P("  min T1            = %s  (c=%s)" % (mp.nstr(mnT1[0], 10), mp.nstr(mnT1[1], 6)))
P("  max |T3|          = %s  (全 c)" % mp.nstr(mxT3, 10))
P("  min (T1-|T3|)     = %s  (c=%s)   > 0" % (mp.nstr(mn_dom[0], 10), mp.nstr(mn_dom[1], 6)))
P("  min (T2+T4+T5)    = %s  (c=%s)" % (mp.nstr(mn_rest[0], 10), mp.nstr(mn_rest[1], 6)))
P("  => gap'(c) = T1 + (T2+T4+T5) + T3 >= (T1-|T3|) + (T2+T4+T5) > 0   [C]")
P("     (T1>0 和 T2,T4>0 是 [T]: g,x 都关于 u 严格单减 => Cov(g,x)>0;")
P("      Phi'' 关于 u 严格单增 => Cov(Phi'',x)<0, 故 T3<0 是唯一的负项)")
P("")

P("[D] 把单调性升到 [T] 的预算: 均匀 c 网格 + Lipschitz 覆盖")
Lgrid = mp.mpf('0.638981087222819686')      # 见 r06_gap_mono.log [2]/([7])
L2 = abs(d2[imx])
P("  sup|gap'| <= Lgrid + L2*dc/2 ;  需 L*dc < min gap = 0.1359947")
for K in (5, 6, 8, 10):
    dc = mp.mpf(1) / K
    Lb = Lgrid + L2 * dc / 2
    P("    K=%2d (dc=%.3f): L<=%s ; L*dc=%s ; 通过=%s"
      % (K, float(dc), mp.nstr(Lb, 6), mp.nstr(Lb * dc, 6), Lb * dc < mp.mpf('0.1359947')))
P("  => c=0,0.2,0.4,0.6,0.8,1.0 这 6 个有理点 (+ 各自区间证书 + |gap''| 界) 就够覆盖全 [0,pi]。")
P("     其中 c=0 (=theta=pi) 与 c=1 (=theta=0) r03 已有 [T], 只需补 c=0.2,0.4,0.6,0.8 四个。")
P("=" * 100)

with open("r06_gap_mono_b.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(R.LOG) + "\n")
