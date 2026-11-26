# -*- coding: utf-8 -*-
"""b05_invariant.py -- B05 物种：三个 𝔐 候选在四模型试金石上的判定
py -3.12, numpy only, deterministic, wall < 5 s.
纪律：数值只做能精确算的；不能精确算的标 [M]（模型给定）或 [T]（标准恒等式）。
"""
import numpy as np
from math import pi, log, exp, gamma, sqrt

np.set_printoptions(precision=6, suppress=True)
TOL = 1e-9


def hdr(s):
    print("\n" + "=" * 88)
    print(s)
    print("=" * 88)


# ============================================================
# 0. 四模型给定数据（试金石输入，[M]）
# ============================================================
# m / G'' / f1^2 per the touchstone (task-given; single-plaquette G'' from
# w7/o02_plaquette.log archived value 0.135995 at beta=2.5, th=pi, matched to 2.9e-7)
TOUCH = {
    "free":      dict(m=0.0,    Gpp=0.0,    f1sq=0.0,  transfer=True),
    "2D-YM":     dict(m=0.375,  Gpp=0.0,    f1sq=0.0,  transfer=True),   # 3g^2/8, g=1
    "M6":        dict(m=0.0,    Gpp=2.0,    f1sq=None, transfer=True),   # 2w, w=1, L->inf
    "plaquette": dict(m=None,   Gpp=0.135995, f1sq=None, transfer=False),
}
TARGET = {"free": 0.0, "2D-YM": +1.0, "M6": 0.0, "plaquette": None}  # None = 无定义
# +1.0 = "必须 >0"，0.0 = "必须 =0"


def verdict(value, target):
    if target is None:
        return "PASS(undef)" if value is None else "FAIL"
    if value is None:
        return "FAIL(undef)"
    if target == 0.0:
        return "PASS" if abs(value) < TOL else "FAIL"
    return "PASS" if value > TOL else "FAIL"


# ============================================================
# 1. 2D YM (SU(2), g^2=1)：转移谱 E_j = (g^2/2) j(j+1)
# ============================================================
def e2dym(j, g2=1.0):
    return 0.5 * g2 * j * (j + 1.0)


def levels_2dym(g2=1.0, N=800):
    j = 0.5 * np.arange(1, N + 1)
    return j, e2dym(j, g2), 2.0 * j + 1.0


hdr("[1] 2D YM SU(2), g^2=1 : transfer spectrum E_j = (g^2/2) j(j+1),  vacuum j=0")
j2, E2, d2 = levels_2dym()
print("   first levels:", np.round(E2[:5], 6), " (j=1/2,1,3/2,2,5/2)")
m_2dym = E2[0]
spacings = np.diff(E2)
r_iso_2dym = float(spacings.min())
print(f"   m = E_(1/2)                        = {m_2dym:.6f}   (=3g^2/8)")
print(f"   min spacing of connected spectrum  = {r_iso_2dym:.6f}   (=5g^2/8, between j=1/2 and j=1)")
print(f"   ratio r_iso/m                      = {r_iso_2dym/m_2dym:.6f}")

# M3 (HS rate) = 2 * lambda_bot  -- 数值确认 -log||T~^t||_HS / t -> 2 E_1
print("   [M3] HS rate  -log( sum_i d_i^2 e^{-2E_i t} ) / t :")
for t in (5.0, 10.0, 20.0, 40.0):
    S = float(np.sum(d2 ** 2 * np.exp(-2.0 * E2 * t)))
    print(f"        t={t:5.1f}   rate = {-(1.0/t)*log(S):.6f}   (limit 2E_1 = {2*m_2dym:.6f})")

# SU(3) cross-check: r_iso/m is NOT a universal rescaling
c2f_su3, c2adj_su3 = 4.0 / 3.0, 3.0
m_su3 = 0.5 * c2f_su3
r_iso_su3 = 0.5 * (c2adj_su3 - c2f_su3)
print(f"   [E] SU(3): m = {m_su3:.6f}, r_iso = {r_iso_su3:.6f}, ratio = {r_iso_su3/m_su3:.6f}"
      f"   vs SU(2) ratio {r_iso_2dym/m_2dym:.6f}  -> r_iso != const* m")

# ============================================================
# 2. M6 (w=1)：E = n*w + 2 pi^2 k^2 / L^2  (E0 = w/2 removed)
# ============================================================
hdr("[2] M6 (w=1): H = (p_a^2+p_b^2)/2 + a^2/2, b on circle of length L")
print("   spectrum = { n*w } + { 2 pi^2 k^2 / L^2 } ; E1-E0 = min(w, 2 pi^2/L^2) -> 0")


def m6_min_gap(L, w=1.0):
    return min(w, 2.0 * pi ** 2 / L ** 2)


def m6_min_spacing(L, w=1.0):
    """min nonzero |E_i - E_j| over the discrete spectrum at finite L"""
    n = np.arange(0, 40)
    k = np.arange(0, 400)
    E = np.unique(np.round(np.outer(n, np.ones_like(k)) * w
                           + np.outer(np.ones_like(n), 2.0 * pi ** 2 * k ** 2 / L ** 2), 12))
    E = E[E > 1e-12]
    return float(np.diff(E).min())


for L in (10.0, 40.0, 160.0):
    print(f"   L={L:6.1f}   m_L = min(w,2pi^2/L^2) = {m6_min_gap(L):.8f}"
          f"   min spacing = {m6_min_spacing(L):.8f}")
print("   both -> 0 as L -> inf :  the a-tower atoms are buried in the b-continuum closure")

# M3 for M6 at finite L : rate = 2 * min{E>0}
for L in (10.0, 40.0, 160.0):
    print(f"   [M3] L={L:6.1f}   HS rate (per vacuum-subtracted bottom) = "
          f"{2.0*m6_min_gap(L):.8f} -> 0")

# ============================================================
# 3. Free field : rho(E) = c E^{1/2}, C(t) = Gamma(3/2) t^{-3/2}
# ============================================================
hdr("[3] Free field (massless): mu = c E^{1/2} dE, C(t) = Gamma(3/2) t^{-3/2}")
print("   no atoms, continuum starts at 0  ->  min spacing = 0 exactly [T]")
print("   [M3] HS rate:  -log( C_HS(t) ) / t  with C_HS ~ t^{-3/2} :")
for t in (10.0, 100.0, 1000.0):
    print(f"        t={t:7.1f}   rate = {1.5*log(t)/t:.8f}  -> 0")
print("   contact carrier exists in the SAME theory: O = box(phi)  =>  C_O(t) = 0 (t>0) [T]")
print("     <box phi(x) box phi(0)> = int d^2k (k^2)^2/(k^2) e^{ikx} ~ box^2 delta(x) = contact")

# ============================================================
# 4. Single plaquette : no transfer direction
# ============================================================
hdr("[4] Single plaquette (d=0): one SU(2) matrix, weight e^{(beta/2) tr U}")
print("   G'' = 0.135995 at beta=2.5, theta=pi  [M/V: w7/o02_plaquette.log, archived 2.9e-7]")
print("   no time direction  =>  no t, no spectrum, no extensive area")
print("   => M1 (area law) : UNDEFINED   M2 (spacing) : UNDEFINED   M3 (HS rate) : UNDEFINED")

# ============================================================
# 5. 三个候选的四模型值
# ============================================================
hdr("[5] M1/M2/M3 candidates x 4 models")
VAL = {
    "M1 = sigma (charged-sector gap / area law)": {
        "free": 0.0, "2D-YM": m_2dym, "M6": 0.0, "plaquette": None},
    "M2 = isolation radius (bottom atom <-> rest)": {
        "free": 0.0, "2D-YM": r_iso_2dym, "M6": 0.0, "plaquette": None},
    "M3 = HS rate (= 2*lambda_bot, control)": {
        "free": 0.0, "2D-YM": 2.0 * m_2dym, "M6": 0.0, "plaquette": None},
    "m (reference, banned)": {
        "free": 0.0, "2D-YM": m_2dym, "M6": 0.0, "plaquette": None},
}
rows = ["free", "2D-YM", "M6", "plaquette"]
print(f"{'candidate':46s}" + "".join(f"{r:>13s}" for r in rows) + "   verdict")
allpass = {}
for name, v in VAL.items():
    cells, ok = "", True
    for r in rows:
        x = v[r]
        cells += f"{'undef' if x is None else f'{x:.6f}':>13s}"
        ok &= (verdict(x, TARGET[r]).startswith("PASS"))
    allpass[name] = ok
    print(f"{name:46s}{cells}   {'4/4 PASS' if ok else 'FAIL'}")

# ============================================================
# 6. signature collision : every passer has m's own signature
# ============================================================
hdr("[6] signature collision (why the 4-model test cannot enforce `not a rename`)")
def pat(x):
    """the only thing the 4 rows can pin: undef / exact 0 / strictly positive"""
    if x is None:
        return "undef"
    return "0" if abs(x) < TOL else "+"


print("   pattern = tuple over (free, 2D-YM, M6, plaquette) of (undef | 0 | +):")
sig_m = tuple(pat(VAL["m (reference, banned)"][r]) for r in rows)
print("   pat(m) =", sig_m)
for name in VAL:
    s = tuple(pat(VAL[name][r]) for r in rows)
    print(f"   {name:46s} pat = {s}   {'== pat(m)  COLLISION' if s == sig_m else '!= pat(m)'}")
print("   => every 4/4 passer has pattern (0,+,0,undef) == pat(m); the numeric value is NOT pinned.")
print("   => M1 is numerically identical to m on 2D-YM (sigma == 3g^2/8, degenerate identity [T]).")
print("   => the test is a (transfer-exists? , m>0?) bit test; pass-set = monotone relabelings of m.")

# ============================================================
# 7. 5th/6th model : what WOULD certify a non-rename
# ============================================================
hdr("[7] extra rows needed to certify a non-rename")
print("   (5) massive free scalar, mass w : m=w>0, G''=2w>0, f1^2=1/(2w)>0, sigma=0, r_iso=w")
print("   (6) free field + contact carrier : C(t)=0 (t>0), m=0  <-- kills every functional of one C")
print("   consequence: sigma is NOT a function of m (m=w>0, sigma=0) -> sigma passes (5);")
print("                but sigma is not a gap witness -> it witnesses a DIFFERENT sector.")
print("   consequence: r_iso(w)=w and M3(w)=2w -> both collapse to m on any relativistic target.")

# ============================================================
# 8. U13 fixture audit : its 2D-YM carrier has f1^2 > 0
# ============================================================
hdr("[8] audit of the existing U13 `pass` (c1/c3/c4)")
print("   U13 used C(t) = sum_j w_j e^{-E_j t} with c1 -> 0.375 = E_(1/2).")
print("   that limit requires w_(1/2) > 0  =>  f1^2 > 0  in the fixture.")
print("   the touchstone's 2D-YM row has f1^2 == 0 (trF^2 = H conserved => C is a contact term).")
print("   => U13's c1/c3/c4 passes are carrier artefacts, not properties of 2D YM.")

hdr("SUMMARY")
for name, ok in allpass.items():
    print(f"   {'PASS 4/4' if ok else 'FAIL      '}  {name}")
print("   -> (a) does NOT hold for any non-rename; (b) holds (pass-set = relabelings of m);")
print("      (c) each candidate's exact failure cause is printed above.")

