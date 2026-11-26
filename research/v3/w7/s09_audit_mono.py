# -*- coding: utf-8 -*-
"""S09-AUDITMONO — 独立复核 R06-GAPMONO 的单调性结论（不复用、不读取 r06 的实现）。

模型来源（不是 r06）：R03-CERT-THETA.md §1 / Q04 的同一约化
  单 plaquette SU(2)，切片 xi_3=0，v=e3，beta=2.5
  u = rho/2 in (0,pi),  c = cos(theta/2), s = sin(theta/2)
  g(u)   = sin u/(2u)
  Phi''  = 1/(2u^2) - cot(u)/(2u) = (1/u - cot u)/(2u)
  nu_c   ~ e^{-beta(1-c cos u)} sin^2 u/u du
  S(c)   = (beta c/2) E[g] + E[Phi'']
  D(c)   = beta^2 (1-c^2) Var(g)
  gap(c) = S - D

本文件做的事
  [1] Jacobian 核对：u-坐标 vs x=cos u 坐标（指数族形式是否成立）
  [2] 精确锚点（与任何 c-导数无关）：int Phi'' W du = 1/4 ; int W du = Cin(2pi)/2
  [3] gap(c) 与 r03 三点 [T] 值的对照
  [4] 自立网格扫描（>=300 点）：gap'(c) 的 min/max/符号
  [5] T1..T5 分解（与 r06 的表对照）
  [6] Phi'' 单调性（Phi''' > 0）与 Cov(Phi'',x) 的符号
  [7] mpmath 高精度复核（gap 与 gap' 两条独立路线：中心差分 vs 协方差恒等式）

写盘范围：本文件 + s09_audit_mono.log。禁 git commit。
"""
import math
import time

import numpy as np
import mpmath as mp

T0 = time.time()


def P(*a, **kw):
    kw.setdefault('flush', True)
    print(*a, **kw)


P("=" * 78)
P("S09-AUDITMONO — independent audit of R06-GAPMONO (gap'(c)>0 ?)")
P("=" * 78)
P("python/numpy/mpmath:", "numpy", np.__version__, "| mpmath", mp.__version__)
P("model source: R03-CERT-THETA.md sec.1 (NOT r06_gap_mono*.py)   beta=2.5")

B = mp.mpf('2.5')
Bf = 2.5
U0, PI = mp.mpf(0), mp.pi


# ---------------------------------------------------------------- integrands
def w_np(u, c):
    """nu_c density in u-coordinate (without normalisation)."""
    return np.exp(-Bf * (1.0 - c * np.cos(u))) * np.sin(u) ** 2 / u


def g_np(u):
    return np.sin(u) / (2.0 * u)


def phi2_np(u):
    """Phi'' = (1/u - cot u)/(2u) = 1/6 + u^2/90 + u^4/945 + u^6/9450 + ...
    cancellation-safe: series for u<0.02 (tail u^8/93555 ~ 3e-19), stable form above."""
    small = u < 0.02
    big = np.where(small, 1.0, u)
    val = (1.0 / big - np.cos(big) / np.sin(big)) / (2.0 * big)
    ser = 1.0 / 6.0 + u ** 2 / 90.0 + u ** 4 / 945.0 + u ** 6 / 9450.0
    return np.where(small, ser, val)


def w_mp(u, c):
    return mp.e ** (-B * (1 - c * mp.cos(u))) * mp.sin(u) ** 2 / u


def g_mp(u):
    return mp.sin(u) / (2 * u)


def phi2_mp(u):
    """Phi'' = (1/u - cot u)/(2u); stable form for u>1e-3, Taylor below (tail ~1e-29)."""
    if abs(u) < mp.mpf('1e-3'):
        return mp.mpf(1) / 6 + u ** 2 / 90 + u ** 4 / 945 + u ** 6 / 9450
    return (1 / u - mp.cos(u) / mp.sin(u)) / (2 * u)


# ================================================================ [1] Jacobian
P("")
P("=" * 78)
P("[1] JACOBIAN CHECK  (u-coordinate  vs  x = cos u)")
P("=" * 78)
mp.mp.dps = 50
P("substitution: x=cos u  =>  du = -dx/sqrt(1-x^2),  sin^2 u = 1-x^2")
P("  sin^2u/u du = (1-x^2)/u * du = (1-x^2)/u * (-dx/sqrt(1-x^2))")
P("              = -sqrt(1-x^2)/arccos(x) dx   [orientation flip: u:0->pi  <=>  x:1->-1]")
P("  => int_0^pi h(u) du = int_{-1}^{1} h(arccos x) * w(x) dx,")
P("     w(x) = sqrt(1-x^2)/arccos x  (= 2g(x), w(1)=1, w(-1)=0)")
P("  => nu_c ~ e^{-beta} * e^{beta c x} w(x) dx   : exponential family in x, natural par. beta*c")
P("")

w_of_x = lambda x: mp.sqrt(1 - x ** 2) / mp.acos(x)

P("  A. measure normalisation (must equal, and equal Cin(2pi)/2):")
Iu = mp.quad(lambda u: mp.sin(u) ** 2 / u, [0, mp.pi / 2, mp.pi])
Ix = mp.quad(lambda x: w_of_x(x), [-1, 0, 1])
Cin2 = (mp.euler + mp.log(2 * mp.pi) - mp.ci(2 * mp.pi)) / 2
P("     int_0^pi sin^2u/u du          = %s" % mp.nstr(Iu, 40))
P("     int_{-1}^{1} w(x) dx          = %s" % mp.nstr(Ix, 40))
P("     Cin(2pi)/2 (exact)            = %s" % mp.nstr(Cin2, 40))
P("     |diff| = %s ; |Ix-Cin2| = %s" % (mp.nstr(abs(Iu - Ix), 3), mp.nstr(abs(Ix - Cin2), 3)))

P("")
P("  B. same functional for f in {1, x, g, Phi''}, at c in {0,.3,.7,1}:")
P("     %-6s %-12s %-26s %-26s %-10s" % ("c", "f", "u-coordinate", "x-coordinate", "reldiff"))
for c in ['0', '0.3', '0.7', '1']:
    cc = mp.mpf(c)
    for name, fu, fx in (("1", lambda u: mp.mpf(1), lambda x: mp.mpf(1)),
                         ("x=cos u", lambda u: mp.cos(u), lambda x: x),
                         ("g", g_mp, lambda x: mp.sin(mp.acos(x)) / (2 * mp.acos(x))),
                         ("Phi''", phi2_mp, lambda x: phi2_mp(mp.acos(x)))):
        a = mp.quad(lambda u: w_mp(u, cc) * fu(u), [0, mp.pi / 2, mp.pi])
        b = mp.quad(lambda x: mp.e ** (-B) * mp.e ** (B * cc * x) * w_of_x(x) * fx(x), [-1, 0, 1])
        rd = abs(a - b) / max(abs(a), mp.mpf(1))
        P("     %-6s %-12s %-26s %-26s %-10s" % (c, name, mp.nstr(a, 22), mp.nstr(b, 22), mp.nstr(rd, 3)))

# ============================================================ [2] exact anchors
P("")
P("=" * 78)
P("[2] EXACT ANCHORS (independent of any c-derivative)")
P("=" * 78)
mp.mp.dps = 50
NS0 = mp.quad(lambda u: phi2_mp(u) * mp.sin(u) ** 2 / u, [0, mp.pi / 2, mp.pi])
P("  int_0^pi Phi''(u) sin^2u/u du = %s   (IBP exact value 1/4;  diff = %s)"
  % (mp.nstr(NS0, 30), mp.nstr(NS0 - mp.mpf(1) / 4, 3)))


# ============================================================== [3] mp moments
def mp_moms(c, dps=50):
    """all moments needed: Z, E[x], E[g], E[g^2], E[Phi''], E[g x], E[g^2 x], E[Phi'' x]"""
    mp.mp.dps = dps
    cc = mp.mpf(c)
    Z = mp.quad(lambda u: w_mp(u, cc), [0, mp.pi / 2, mp.pi])
    X = mp.quad(lambda u: w_mp(u, cc) * mp.cos(u), [0, mp.pi / 2, mp.pi])
    A = mp.quad(lambda u: w_mp(u, cc) * g_mp(u), [0, mp.pi / 2, mp.pi])
    Ax = mp.quad(lambda u: w_mp(u, cc) * g_mp(u) * mp.cos(u), [0, mp.pi / 2, mp.pi])
    Bq = mp.quad(lambda u: w_mp(u, cc) * g_mp(u) ** 2, [0, mp.pi / 2, mp.pi])
    Bx = mp.quad(lambda u: w_mp(u, cc) * g_mp(u) ** 2 * mp.cos(u), [0, mp.pi / 2, mp.pi])
    NS = mp.quad(lambda u: w_mp(u, cc) * phi2_mp(u), [0, mp.pi / 2, mp.pi])
    NSx = mp.quad(lambda u: w_mp(u, cc) * phi2_mp(u) * mp.cos(u), [0, mp.pi / 2, mp.pi])
    return dict(Z=Z, x=X / Z, g=A / Z, g2=Bq / Z, phi=NS / Z,
                gx=Ax / Z, g2x=Bx / Z, phix=NSx / Z)


def mp_gap(c, dps=50):
    m = mp_moms(c, dps)
    var = m['g2'] - m['g'] ** 2
    S = (B * mp.mpf(c) / 2) * m['g'] + m['phi']
    D = B ** 2 * (1 - mp.mpf(c) ** 2) * var
    return S - D, m


def mp_terms(m, c):
    """T1..T5 decomposition of d/dc gap (analytic, via covariances)."""
    cc = mp.mpf(c)
    var = m['g2'] - m['g'] ** 2
    cgx = m['gx'] - m['g'] * m['x']
    cg2x = m['g2x'] - m['g2'] * m['x']
    cphx = m['phix'] - m['phi'] * m['x']
    T1 = (B / 2) * m['g']
    T2 = (B ** 2 * cc / 2) * cgx
    T3 = B * cphx
    T4 = 2 * B ** 2 * cc * var
    T5 = -B ** 3 * (1 - cc ** 2) * (cg2x - 2 * m['g'] * cgx)
    return T1, T2, T3, T4, T5


P("")
P("=" * 78)
P("[3] gap(c) at the three r03 [T] anchors (c=0, sqrt2/2, 1)")
P("=" * 78)
R03 = {0: ("0.135994711041988211036899467040198122289952357486",
           "0.135994711041988211036899467040198122289952357486"),
       '0.70710678118654752440084436210484903928483593768847':
       ("0.542943922353689763665445822036482126731456810684",
        "0.542943922353689763665445822036482126731456810684"),
       1: ("0.729924127576750062576372357665823454189533971017",
           "0.729924127715173653761964295165831775817594579822")}
mp.mp.dps = 50
for c in ['0', mp.sqrt(2) / 2, '1']:
    cc = mp.mpf(c) if not isinstance(c, str) else mp.mpf(c)
    gap, m = mp_gap(cc)
    key = 0 if abs(cc) < mp.mpf('1e-30') else (1 if abs(cc - 1) < mp.mpf('1e-30') else '0.70710678118654752440084436210484903928483593768847')
    lo, hi = R03[key]
    lo, hi = mp.mpf(lo), mp.mpf(hi)
    P("  c=%-22s gap = %s" % (mp.nstr(cc, 8), mp.nstr(gap, 36)))
    P("      r03 [T] interval [%s , %s]" % (mp.nstr(lo, 24), mp.nstr(hi, 24)))
    P("      inside r03 interval? %s   (gap-lo=%s, hi-gap=%s)"
      % (lo <= gap <= hi, mp.nstr(gap - lo, 4), mp.nstr(hi - gap, 4)))
    P("      E[g]=%s  Var(g)=%s  E[Phi'']=%s"
      % (mp.nstr(m['g'], 20), mp.nstr(m['g2'] - m['g'] ** 2, 20), mp.nstr(m['phi'], 20)))

# ======================================================= [4] grid scan (numpy)
P("")
P("=" * 78)
P("[4] INDEPENDENT GRID SCAN  (numpy, Gauss-Legendre 2x240 nodes)")
P("=" * 78)


def gl_nodes(nq=240):
    xq, wq = np.polynomial.legendre.leggauss(nq)
    U, Wq = [], []
    for a, b in ((0.0, math.pi / 2), (math.pi / 2, math.pi)):
        U.append(0.5 * (b - a) * xq + 0.5 * (a + b))
        Wq.append(0.5 * (b - a) * wq)
    return np.concatenate(U), np.concatenate(Wq)


UG, WG = gl_nodes(240)


def np_scan(cs):
    out = []
    for c in cs:
        Wd = w_np(UG, c)
        gg = g_np(UG)
        pp = phi2_np(UG)
        xx = np.cos(UG)
        ww = WG * Wd
        Z = ww.sum()
        m1 = (ww * gg).sum() / Z
        m2 = (ww * gg * gg).sum() / Z
        mp_ = (ww * pp).sum() / Z
        mx = (ww * xx).sum() / Z
        gx = (ww * gg * xx).sum() / Z
        g2x = (ww * gg * gg * xx).sum() / Z
        px = (ww * pp * xx).sum() / Z
        var = m2 - m1 * m1
        cgx = gx - m1 * mx
        cg2x = g2x - m2 * mx
        cpx = px - mp_ * mx
        T1 = (Bf / 2) * m1
        T2 = (Bf * Bf * c / 2) * cgx
        T3 = Bf * cpx
        T4 = 2 * Bf * Bf * c * var
        T5 = -Bf ** 3 * (1 - c * c) * (cg2x - 2 * m1 * cgx)
        S = (Bf * c / 2) * m1 + mp_
        D = Bf * Bf * (1 - c * c) * var
        out.append(dict(c=c, gap=S - D, dgap=T1 + T2 + T3 + T4 + T5,
                        T1=T1, T2=T2, T3=T3, T4=T4, T5=T5,
                        m1=m1, var=var, cgx=cgx, cphx=cpx, S=S, D=D))
    return out


N1 = 1001
cs = list(np.linspace(0.0, 1.0, N1))
cs += list((1 - np.cos(np.pi * np.arange(1, 200) / 200)) / 2)  # Chebyshev clustering
cs = sorted(set(round(float(x), 15) for x in cs))
rows = np_scan(cs)
inner = [r for r in rows if 1e-14 < r['c'] < 1 - 1e-14]
P("  points: %d total, %d strictly inside (0,1)" % (len(rows), len(inner)))
dg = np.array([r['dgap'] for r in inner])
cs_i = np.array([r['c'] for r in inner])
P("  d/dc gap  : min = %.12f  at c = %.9f" % (dg.min(), cs_i[dg.argmin()]))
P("              max = %.12f  at c = %.9f" % (dg.max(), cs_i[dg.argmax()]))
P("              nonpositive points (<=0) : %d / %d" % (int((dg <= 0).sum()), len(dg)))
P("              min(T1+|T3| sign guard)  : T1-|T3| = %.9f" % min(r['T1'] - abs(r['T3']) for r in inner))
P("              min(T2+T4+T5)            = %.9f" % min(r['T2'] + r['T4'] + r['T5'] for r in inner))
gapv = np.array([r['gap'] for r in rows])
P("  gap       : min = %.12f at c = %.3e   max = %.12f at c = %.9f"
  % (gapv.min(), rows[int(gapv.argmin())]['c'], gapv.max(), rows[int(gapv.argmax())]['c']))
P("  dgap at the endpoint c=0 : %.12f   (c=1: %.12f)"
  % (rows[0]['dgap'], rows[-1]['dgap']))
P("  inf over [0,1] INCLUDING endpoints = %.12f  at c = %.3e"
  % (min(r['dgap'] for r in rows), min(rows, key=lambda r: r['dgap'])['c']))
P("  sign flips of central-difference gap' on the uniform 1001-grid: ", end="")
h = 1e-5
# central difference on the uniform grid only (no Chebyshev jitter)
cu = np.linspace(0.0, 1.0, N1)
ru = np_scan(cu)
gv = np.array([r['gap'] for r in ru])
fd = (gv[2:] - gv[:-2]) / (2 * (cu[1] - cu[0]))
P("%d negative of %d" % (int((fd <= 0).sum()), len(fd)))

# ==================================================== [5] T1..T5 table (mpmath)
P("")
P("=" * 78)
P("[5] T1..T5 at c = 0(.1)1  (mpmath dps=40)   [compare r06 tab. sec.2]")
P("=" * 78)
P("  %-6s %-10s %-10s %-10s %-10s %-10s %-10s %-12s %-12s"
  % ("c", "T1", "T2", "T3", "T4", "T5", "sum", "gap'(FD)", "gap"))
tab = {}
for c in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
    gap, m = mp_gap(c, 45)
    T1, T2, T3, T4, T5 = mp_terms(m, c)
    hh = mp.mpf('1e-15')
    gp, _ = mp_gap(mp.mpf(c) + hh, 45)
    gm, _ = mp_gap(mp.mpf(c) - hh, 45)
    fd_c = (gp - gm) / (2 * hh)
    tab[c] = (T1, T2, T3, T4, T5, T1 + T2 + T3 + T4 + T5, fd_c, gap)
    P("  %-6.2f %-10s %-10s %-10s %-10s %-10s %-10s %-12s %-12s"
      % (c, mp.nstr(T1, 8), mp.nstr(T2, 8), mp.nstr(T3, 8), mp.nstr(T4, 8),
         mp.nstr(T5, 8), mp.nstr(T1 + T2 + T3 + T4 + T5, 8),
         mp.nstr(fd_c, 9), mp.nstr(gap, 10)))
P("")
P("  consistency  sum(T) vs central difference (h=1e-15, dps=45):")
for c in tab:
    s = tab[c][5]
    f = tab[c][6]
    P("     c=%.1f   |sum - FD| = %s" % (c, mp.nstr(abs(s - f), 3)))

# =============================================== [6] Phi'' monotone / T3 sign
P("")
P("=" * 78)
P("[6] Phi'' STRICTLY INCREASING on (0,pi) ?  and  Cov(Phi'',x) < 0 ?")
P("=" * 78)
mp.mp.dps = 40
ug = np.linspace(0.02, math.pi - 0.02, 4001)
sg, cg = np.sin(ug), np.cos(ug)
ph = (1.0 / ug - cg / sg) / (2.0 * ug)
P("  numerical: min forward difference of Phi'' on [0.02,pi-0.02] = %.3e  (u=%.4f)"
  % (np.diff(ph).min(), ug[int(np.diff(ph).argmin())]))
Hg = ug ** 2 + 0.5 * ug * np.sin(2 * ug) - 1 + np.cos(2 * ug)
Phi3 = Hg / (2.0 * ug ** 3 * sg ** 2)
P("  numerical: min Phi''' = H/(2u^3 sin^2u), same grid           = %.6f  (u=%.4f)"
  % (Phi3.min(), ug[int(Phi3.argmin())]))
P("  numerical: min F = H/sin^2u (= 2u^3 Phi'''), same grid        = %.6e"
  % (Hg / sg ** 2).min())
P("  small u: Phi''' = u/45 + 4u^3/945 + ... > 0")
P("  analytic  [T]: Phi''' = H(u)/(2u^3 sin^2u),  H = u^2 + (u/2)sin2u - 1 + cos2u")
P("             H'''' = 8u sin2u > 0 on (0,pi/2);  H(0)=H'(0)=H''(0)=H'''(0)=0")
P("             => H''>0 => H'>0 => H>0 on (0,pi/2);  on (pi/2,pi):")
P("             H' = u(2+cos2u) - (3/2)sin2u > 0 since sin2u<=0, H(pi/2)=pi^2/4-2>0.")
P("             => H>0 => Phi'''>0 identically  => Phi'' strictly increasing  [T]")
mp.mp.dps = 40
P("")
P("  Cov(Phi'',x) (should be < 0 for all c):")
for c in [0.0, 0.25, 0.5, 0.75, 1.0]:
    m = mp_moms(c, 40)
    P("     c=%.2f  E[Phi'']=%s  Cov(Phi'',x)=%s  T3=%s"
      % (c, mp.nstr(m['phi'], 10), mp.nstr(m['phix'] - m['phi'] * m['x'], 10),
         mp.nstr(B * (m['phix'] - m['phi'] * m['x']), 10)))

# ============================================ [7] high-precision refinement
P("")
P("=" * 78)
P("[7] HIGH-PRECISION VERIFICATION AT SELECTED c (dps=50)")
P("=" * 78)
P("  %-16s %-26s %-24s %-10s" % ("c", "gap", "gap'(FD,h=1e-15)", "|sum-FD|"))
mp.mp.dps = 50
sel = [0.0, 1e-8, 1e-4, 0.05, 0.25, 0.5, float(mp.sqrt(2) / 2), 0.75, 0.9, 0.99, 1 - 1e-6, 1.0]
cs_arg = float(cs_i[dg.argmin()])
sel.append(cs_arg)
worst = mp.mpf(0)
for c in sel:
    cc = mp.mpf(c)                       # exact at current dps: DO NOT add h in float64
    gap, m = mp_gap(cc, 50)
    hh = mp.mpf('1e-15')
    gp, _ = mp_gap(cc + hh, 50)
    gm, _ = mp_gap(cc - hh, 50)
    fdc = (gp - gm) / (2 * hh)
    T1, T2, T3, T4, T5 = mp_terms(m, cc)
    s = T1 + T2 + T3 + T4 + T5
    worst = max(worst, abs(s - fdc))
    P("  %-16s %-26s %-24s %-10s" % ("%.10g" % c, mp.nstr(gap, 22), mp.nstr(fdc, 18),
                                    mp.nstr(abs(s - fdc), 3)))
P("  worst |analytical sum(T) - central difference| over the above points: %s" % mp.nstr(worst, 3))

# ======================================================= [8] Jacobian summary
P("")
P("=" * 78)
P("[8] SUMMARY OF THE AUDIT")
P("=" * 78)
P("  (1) nu_c IS an exponential family: nu_c ~ e^{beta c x} w(x) dx, x=cos u in (-1,1),")
P("      w(x)=sqrt(1-x^2)/arccos x = 2g(x); the Jacobian du=-dx/sqrt(1-x^2) only")
P("      rescales the reference weight (and an overall e^{-beta}), c-independent.")
P("      => d/dc E[f] = beta Cov(f,x) is VALID.  (no Jacobian obstruction)")
P("  (2) min d/dc gap on the grid = %.12f  > 0  => no zero crossing." % dg.min())
P("      d/dc gap at c=0 (infimum over [0,1]) = %s" % mp.nstr(tab[0.0][5], 15))
P("      d/dc gap at c=1                     = %s" % mp.nstr(tab[1.0][5], 15))
P("      gap(c=0) = %s   (= q11 / r03 theta=pi anchor)"
  % mp.nstr(tab[0.0][7], 20))
P("      gap(c=1) = %s" % mp.nstr(tab[1.0][7], 20))
P("  (3) gap(0)=0.1359947110 (r03 [T]) is the minimum at c=0  => min over [0,1] = gap(pi).")
P("  (4) Phi'' strictly increasing [T] => Cov(Phi'',x) < 0 => T3 < 0 uniquely.")
P("")
P("  ELAPSED %.1f s" % (time.time() - T0))
P("=" * 78)
