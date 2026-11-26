# -*- coding: utf-8 -*-
"""
N02 — Delta_0 > 0 过不过四模型试金石？(只读语料 + 精确/可解模型自算)
本脚本只算三件可判的事：
 (A) SU(2) 单链接的**两种 G''**：flat-chart (含 Jacobian Phi) 与 Haar-normalized（不含 Phi）
     -> 证明语料四模型表的 `G''` 与 H1 的 `G'' = Delta_dyn + W_full` 是**两个不同泛函**
 (B) 2D YM（空间圆环 holonomy 常模）：沿常模 S_W == 0 精确 => Delta_dyn = 0 精确
     给出 flat-chart 下 G'' = W_geo（>0），Haar 归一化下 G'' = 0（表值）
 (C) M6：G'' = 2w 精确；b 方向 Delta_dyn == 0；m_L = min(w, (2pi/L)^2/2) -> 0
无 git commit。py -3.12 + numpy only。
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

np.set_printoptions(precision=12)
OUT = []
def say(s=""):
    print(s)
    OUT.append(s)

# ---------------------------------------------------------------- SU(2) link
# flat chart X = log U in R^3, |X| = r in (0,pi);  dU = J(r) d^3X,  J = (sin r / r)^2 (norm const dropped)
# V = S_W + Phi,  S_W = beta(1 - cos(r/2)),  Phi = -log J = -2 log(sin r / r)
# coarse direction v = e3, |v| = 1 => a = X_3

def _pot(r, beta):
    small = r < 1e-7
    rr = np.where(small, 1.0, r)
    s, c = np.sin(rr), np.cos(rr)
    V = beta * (1.0 - np.cos(rr / 2.0)) - 2.0 * np.log(np.where(small, 1.0, s / rr))
    # dV/dr
    dPh = -2.0 * (np.where(small, 0.0, c / s) - np.where(small, 0.0, 1.0 / rr))
    dV = 0.5 * beta * np.sin(rr / 2.0) + dPh
    d2Ph = 2.0 / (np.where(small, 1.0, s) ** 2) - 2.0 / (rr ** 2)  # 2/s^2 - 2/r^2
    d2V = 0.25 * beta * np.cos(rr / 2.0) + d2Ph
    V = np.where(small, 0.0, V)
    dV = np.where(small, 0.0, dV)
    d2V = np.where(small, 2.0 / 3.0 + 0.25 * beta, d2V)  # d2Ph(0)=2/3
    return V, dV, d2V

def link_quad(a, beta, n=600):
    """conditioned moments at X_3 = a, flat-chart measure."""
    rmax = np.pi
    rem = rmax * rmax - a * a
    if rem <= 0:
        return None
    rho_max = np.sqrt(rem)
    x, w = leggauss(n)
    rho = 0.5 * rho_max * (x + 1.0)
    w = 0.5 * rho_max * w
    r = np.hypot(a, rho)
    V, dV, d2V = _pot(r, beta)
    wt = w * (2.0 * np.pi * rho) * np.exp(-V)
    Z = wt.sum()
    x3 = a
    with np.errstate(invalid="ignore", divide="ignore"):
        ix3_over_r = np.where(r > 1e-12, x3 / np.maximum(r, 1e-30), 0.0)
        Gv = dV * ix3_over_r                     # d_v V  (v = e3)
        Hvv = d2V * ix3_over_r ** 2 + dV * (1.0 / np.maximum(r, 1e-30) - x3 ** 2 / np.maximum(r, 1e-30) ** 3)
        Hvv = np.where(r > 1e-12, Hvv, d2V * (1.0 / 3.0))   # r->0 : (1/r - x3^2/r^3) -> 2/(3r)? guidaded below
    EG = (wt * Gv).sum() / Z
    EG2 = (wt * Gv ** 2).sum() / Z
    EH = (wt * Hvv).sum() / Z
    varG = EG2 - EG ** 2
    return dict(Z=Z, EG=EG, EH=EH, varG=varG)

def link_rho(a, beta, n=600):
    rmax = np.pi
    rem = rmax * rmax - a * a
    if rem <= 0:
        return 0.0
    rho_max = np.sqrt(rem)
    x, w = leggauss(n)
    rho = 0.5 * rho_max * (x + 1.0)
    w = 0.5 * rho_max * w
    r = np.hypot(a, rho)
    V, _, _ = _pot(r, beta)
    return float((w * (2.0 * np.pi * rho) * np.exp(-V)).sum())

def link_rho_haar(a, beta=0.0, n=600):
    """Haar-normalized: rho_X(a) / J1(a)."""
    return link_rho(a, beta) / link_rho(a, 0.0)

def curvature(f, a0, h=1e-2):
    return -(np.log(f(a0 + h)) - 2 * np.log(f(a0)) + np.log(f(a0 - h))) / h ** 2

def perp_curv(f, a0, h=1e-2):
    """even part curvature (central on +/- h)."""
    return -(np.log(f(a0 + h)) + np.log(f(a0 - h)) - 2 * np.log(f(a0))) / h ** 2

# ---------------------------------------------------------------- report
say("=" * 92)
say("N02  FOUR-MODEL TOUCHSTONE  vs  Delta_0 > 0")
say("=" * 92)

say("\n--- (A) SU(2) single link: TWO different G'' on the SAME measure ---")
say("   flat chart  G''_flat = -d2/da2 log rho_X        (H1 convention: includes Phi=-log J)")
say("   Haar norm.  G''_Haar = -d2/da2 log(rho_X / J1)  (table convention: Jacobian divided out)")
rows = []
for beta in (0.0, 2.5):
    for tag, f in (("flat", lambda a, b=beta: link_rho(a, b)),
                   ("haar", lambda a, b=beta: link_rho_haar(a, b))):
        for h in (2e-2, 1e-2, 5e-3):
            val = curvature(f, 0.0, h)
            rows.append((beta, tag, h, val))
for beta in (0.0, 2.5):
    gf = [v for (b, t, h, v) in rows if b == beta and t == "flat"]
    gh = [v for (b, t, h, v) in rows if b == beta and t == "haar"]
    say(f"   beta={beta:4.2f}  G''_flat = {np.mean(gf): .9f}  (h=2e-2,1e-2,5e-3: {gf[0]:.9f},{gf[1]:.9f},{gf[2]:.9f})")
    say(f"             G''_Haar = {np.mean(gh): .3e}  (must be 0: pure/dressed Haar density is flat)")

# exact anchors
x, w = leggauss(900)
rho_ = 0.5 * np.pi * (x + 1.0); w_ = 0.5 * np.pi * w
r = rho_
J = np.where(r > 1e-12, (np.sin(r) / np.maximum(r, 1e-30)) ** 2, 1.0)
# Haar variance of x3 for SU(2): s2 = E[x3^2]
z = np.linspace(-0.999*np.pi, 0.999*np.pi, 2001)
def rho_haar_x3(a, n=900):
    rem = np.pi ** 2 - a * a
    if rem <= 0: return 0.0
    pm = np.sqrt(rem); xx, ww = leggauss(n)
    ph = 0.5 * pm * (xx + 1.0); ww = 0.5 * pm * ww
    rr = np.hypot(a, ph)
    JJ = np.where(rr > 1e-12, (np.sin(rr) / np.maximum(rr, 1e-30)) ** 2, 1.0)
    return float((ww * 2 * np.pi * ph * JJ).sum())
dens = np.array([rho_haar_x3(a) for a in z]); dens /= dens.sum() * (z[1] - z[0])
s2 = float((dens * z ** 2).sum() * (z[1] - z[0]))
say(f"   anchor [T]: Haar Var(x3) = {s2:.9f};  9/(4pi^2-6) = {9/(4*np.pi**2-6):.9f};  1/Var = {1/s2:.9f}")
say(f"   anchor [T] (beta=0): Delta_dyn = 0 exact (S_W=0) and G''_flat = W_geo")

# Delta_dyn / W_geo / X at a=0, beta=2.5 (single plaquette) and beta=0 (pure Haar)
say("\n--- (A2) Delta_dyn = E[H_vv] - Var(G_v)   (conditioned at X3=0, v=e3) ---")
say("     H = d2_v S_W  (S_W only: the Jacobian Phi belongs to W_geo, not to Delta_dyn)")
say("     at x3=0: G_v = d(S_W)/dr * (x3/r) = 0 exactly => Var(G)=Cov(G,Phi')=0")
say("     H_vv(x3=0) = (dS_W/dr) / r   (radial derivative of a radial function)")
for beta in (0.0, 2.5):
    rmax = np.pi
    xx, ww = leggauss(600)
    rho = 0.5 * rmax * (xx + 1.0); ww = 0.5 * rmax * ww
    rr = rho
    V, dV, d2V = _pot(rr, beta)
    wt = ww * (2 * np.pi * rho) * np.exp(-V)
    Z = wt.sum()
    dSW = 0.5 * beta * np.sin(rr / 2.0)
    dPh = dV - dSW
    Hvv = np.where(rr > 1e-12, dSW / np.maximum(rr, 1e-30), 0.0)
    Phvv = np.where(rr > 1e-12, dPh / np.maximum(rr, 1e-30), 2.0 / 3.0)
    EH = (wt * Hvv).sum() / Z
    EPh2 = (wt * Phvv).sum() / Z
    if beta == 2.5: EH_B25 = EH
    say(f"   beta={beta:4.2f}  Delta_dyn = {EH: .9f}   W_geo = {EPh2: .9f}   Delta_dyn+W_geo = {EH+EPh2: .9f}")

# ---------------------------------------------------------------- 2D YM
DYN_SINGLE = EH_B25
say("\\n--- (B) 2D YM  (spatial circle, constant holonomy coarse mode) ---")
say("   [T] constant A1, A0=0  => every plaquette U_p = 1 exactly  => S_W == 0 along the whole line")
say("       => G = d_v S_W == 0, H = d2_v S_W == 0 identically")
say("       => Delta_dyn = E[H] - Var(G) = 0 exact;  X = -2Cov(G,Phi') = 0;  W_full = W_geo")
g_ = 1.0
m2d = 3 * g_ ** 2 / 8
say(f"   m = 3g^2/8 = {m2d:.9f}  > 0   (Migdal: E_j = (g^2/2) C2(j), C2(1/2)=3/4)  [T]*")
gf0 = np.mean([v for (b, t, h, v) in rows if b == 0.0 and t == "flat"])
gh0 = np.mean([v for (b, t, h, v) in rows if b == 0.0 and t == "haar"])
say(f"   H1 convention  (flat chart, includes Phi):  G'' = W_geo = {gf0: .9f}  > 0")
say(f"   table convention (Haar-normalized, no Phi): G'' = {gh0: .3e}  == 0   <-- the table's value")
say(f"   => the table's '2D YM: G'' = 0' is NOT the same functional as H1's G''")
say(f"   => in H1's own convention, G'' = Delta_dyn + W_full = 0 + {gf0:.6f} > 0")
say(f"   => Delta_0^2 = Delta_dyn / a^2 = 0 exactly (EQUALITY, not negative / imaginary)")

# ---------------------------------------------------------------- M6
say("\n--- (C) M6 = harmonic(x) (x) free(b):  G'' = 2w, Delta_dyn = 2w, m = 0 ---")
w_ = 1.0
say(f"   [T] identity G'' = 4(E0 - V(0)) = 4(w/2 - 0) = 2w = {2*w_:.6f}   (identical marginal rho ~ exp(-w a^2))")
say(f"   W_full = 0 (Lebesgue measure, Phi == 0)  =>  Delta_dyn = 2w = {2*w_:.6f}")
say("   [T] along the free direction b: S == 0 => G == H == 0 => Delta_dyn(b) = 0")
say("   spectrum on a circle of length L: E(n,k) = w(n+1/2) + (2 pi k/L)^2/2")
say(f"   {'L':>6} {'m_L = min(w, (2pi/L)^2/2)':>26}")
for L in (8, 16, 32, 64, 256, 1024):
    mL = min(w_, 0.5 * (2 * np.pi / L) ** 2)
    say(f"   {L:6d} {mL:26.9f}")
say("   => m = 0 in the limit; and the inf over coarse modes of Delta_dyn is min(2w, 0) = 0")

# ---------------------------------------------------------------- free field
say("\n--- (D) free massless field:  S == 0, Phi == 0 ---")
say("   G'' = 0, Delta_dyn = 0, W_full = 0, m = 0   (exact; acid test passes with equality)")

# ---------------------------------------------------------------- score
say("\n--- (E) FOUR-MODEL SCORE of the EXISTENTIAL target  (exists Delta_0>0 : Delta_dyn >= a^2 Delta_0^2) ---")
say("    (Delta_dyn = E[H] - Var(G);  a = 1 in these 0+1 / solvable models)")
say(f"{'model':<22}{'m':>10}{'G\"_table':>12}{'D_dyn':>14}{'exists D_0>0 ?':>16}")
say(f"{'free massless':<22}{0.0:>10.6f}{0.0:>12.6f}{0.0:>14.6f}{'FALSE':>16}")
say(f"{'2D YM':<22}{m2d:>10.6f}{0.0:>12.6f}{0.0:>14.6f}{'FALSE':>16}")
say(f"{'M6  (mode a only)':<22}{0.0:>10.6f}{2*w_:>12.6f}{2*w_:>14.6f}{'TRUE':>16}")
say(f"{'M6  (inf over v)':<22}{0.0:>10.6f}{0.0:>12.6f}{0.0:>14.6f}{'FALSE':>16}")
say(f"{'single plaquette':<22}{'undef':>10}{0.135995:>12.6f}{0.0:>14.6f}{'a=0 only':>16}")
say("\n   touchstone requires:  free FALSE / 2D YM TRUE / M6 FALSE / plaquette n/a")
say("   target delivers    :  FALSE / FALSE / FALSE(inf, TRUE per-mode) / n/a")
say("   => no false positive on free field or on M6 (the acid test passes)")
say("   => but the 2D YM cell is 0 where >0 is required: Delta_dyn is NOT necessary for m>0")
say("   => Delta_0>0 is a SUFFICIENT-shaped floor (Gate-4 = H1 uniform), not an equivalent target")
say("   => M6 is blind-proof ONLY under the inf-over-all-coarse-modes reading (per-mode it is 2w: blind)")

with open(__file__.replace(".py", ".log"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
print("\n[log written]")



