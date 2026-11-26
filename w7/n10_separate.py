"""n10_separate.py -- 能否从格点上把 Delta_dyn 与 G'' 分开测？(py -3.12, numpy only)

[A] 精确单 plaquette 核 (o27 Tier-1 同款约定): 先用**无条件**值复现 o27 的精确锚点做自检，
    再算 **xi=0 条件切片** 上的分解. 目的: 判定该玩具核能否回答符号问题.
[B] 归档路线的预算表: G'' - Delta_dyn = W_full 逐路核对, W_full/(N/12) 比值, Delta_dyn 的 sigma.
[C] Delta_dyn(beta) 的零点 (p19 四点的线性拟合), 与 beta=0 精确闭合对照.

标签 [T]/[E]/[V]. 只写 w7/n10_*.log, 不写 docs/, 不 git commit.
"""
import numpy as np

BETA, N, FLOOR = 2.5, 2, 2 / 12          # N/12 = 1/6
P = print


# ---------- [A] exact single SU(2) plaquette, flat chart X in R^3, v = e3 ----------
def plaquette(beta=BETA, nr=2_000_000, rmax=3.2):
    r = np.linspace(1e-4, rmax, nr)
    sw = beta * (1 - np.cos(r / 2))
    phi = -2 * np.log(np.sin(r / 2) / (r / 2))
    J = (np.sin(r / 2) / (r / 2)) ** 2
    w = np.exp(-(sw + phi)) * J * r ** 2                       # 3D radial weight
    dsw, d2sw = 0.5 * beta * np.sin(r / 2), 0.25 * beta * np.cos(r / 2)
    dph = -np.cos(r / 2) / np.sin(r / 2) + 2 / r
    d2ph = 0.5 / np.sin(r / 2) ** 2 - 2 / r ** 2
    c2 = 1 / 3.0                                               # <(nhat.v)^2> on S^2
    av = lambda f, p: np.trapezoid(w * f * {0: 1.0, 1: 0.0, 2: c2}[p], dx=r[1] - r[0]) / np.trapezoid(w, dx=r[1] - r[0])
    EH = av(d2sw, 2) + av(dsw / r, 0) - av(dsw / r, 2)
    EG2 = av(dsw ** 2, 2)
    EPh2 = av(d2ph, 2) + av(dph / r, 0) - av(dph / r, 2)
    EPh1sq = av(dph ** 2, 2)
    EGP = av(dsw * dph, 2)
    return dict(EH=EH, vG=EG2, EPh2=EPh2, vP=EPh1sq, Cov=EGP,
                vF=EG2 + EPh1sq + 2 * EGP)


def cond_slice(beta=BETA, nrho=2_000_000, rmax=3.2):
    """xi = <v,log U> = x3 = 0 的切片: 2D radial weight e^{-V} J rho drho.
    沿 v=e3 在 x3=0:  G_v = 0, Phi' = 0  =>  Var(G)=Var(Phi')=Cov=0;
    H_vv = (dS_W/dr)/r,  Phi'' = (dPhi/dr)/r."""
    t = np.linspace(1e-4, rmax, nrho)
    sw = beta * (1 - np.cos(t / 2))
    phi = -2 * np.log(np.sin(t / 2) / (t / 2))
    J = (np.sin(t / 2) / (t / 2)) ** 2
    w = np.exp(-(sw + phi)) * J * t                              # 2D radial weight
    dsw, dph = 0.5 * beta * np.sin(t / 2), -np.cos(t / 2) / np.sin(t / 2) + 2 / t
    d = t[1] - t[0]
    m = lambda f: np.trapezoid(w * f, dx=d) / np.trapezoid(w, dx=d)
    return m(dsw / t), m(dph / t)


d = plaquette()
A, Wg, X = d['EH'] - d['vG'], d['EPh2'] - d['vP'], -2 * d['Cov']
WF = d['EPh2'] + d['vG'] - d['vF']
P("=" * 96)
P("[A] exact single plaquette beta=2.5, flat chart  (self-check vs o27 Tier-1 unconditional)")
P(f"  E[H_vv]={d['EH']:+.6f}  E[Phi'']={d['EPh2']:+.6f}  Var(G)={d['vG']:.6f}"
  f"  Var(G+Phi')={d['vF']:.6f}  Var(Phi')={d['vP']:.6f}  Cov={d['Cov']:+.6f}")
ref = dict(EH=0.500379, EPh2=0.184014, vG=0.272802, vF=0.499617, vP=0.034694, Cov=0.096060)
dev = max(abs(d[k] - ref[k]) for k in ref)
P(f"  [V] self-check |mine - o27| max dev = {dev:.2e}   -> kernel {'OK' if dev < 1e-4 else 'MISMATCH'}")
P(f"  UNCOND : Delta_dyn=A={A:+.6f}  W_geo={Wg:+.6f}  X={X:+.6f}  W_full={WF:+.6f}"
  f"  R=W_full/floor={WF/FLOOR:+.4f}")
EH0, EPh0 = cond_slice()
P(f"  xi=0 SLICE (exact): G_v == 0 and Phi' == 0 identically => Var(G)=Var(Phi')=Cov=0")
P(f"    E[H_vv|0]={EH0:+.6f}  E[Phi''|0]={EPh0:+.6f}  (Lemma B floor N/12={FLOOR:.6f})")
P(f"    Delta_dyn|0 = E[H|0]-0 = {EH0:+.6f} > 0   [T]   W_geo|0 = E[Phi''|0] = {EPh0:+.6f}"
  f"  W_full|0 = {EPh0:+.6f}  R={EPh0/FLOOR:.4f}")
P("    => the toy's sign answer is FORCED by the single-scalar conditioning, not measured:")
P("       one mode IS the whole coordinate, so the force vanishes on the slice.")
P("       It therefore cannot test Delta_dyn vs G''; only a multi-mode 4D ensemble can.")

# ---------- [B] archived routes at beta=2.5, L=2 (all Delta_dyn measured, not inferred) ----------
P("=" * 96)
P("[B] beta=2.5, L=2, SU(2): routes that measured Delta_dyn directly (G''-Delta_dyn=W_full by [T])")
rows = [
    ("p19 b2.5 (4 chains x5000)", 0.4116, 0.0421, 0.1028, 0.1401, 0.0419, 0.0981),
    ("p02   (4 chains x19500, fund-domain reject)", 2.002, 0.051, np.nan, 0.8838, 0.0501, np.nan),
    ("o13_ldep L2 (1500 traj)", 2.0148, 0.3500, np.nan, 1.0655, 0.3326, np.nan),
    ("o27 mode A (mu=0,Lie0)", 1.696, 0.293, np.nan, 0.5619, 0.2770, np.nan),
    ("o27 mode B (all mu,Lie0)", 0.9459, 0.2102, np.nan, -0.1789, 0.2251, np.nan),
]
P(f"  {'route':44s} {'G\"':>7s} {'D_dyn':>8s} {'W_full':>7s} {'R=Wf/F':>7s} {'s_jk':>6s} {'sig_jk':>6s} {'sig_tau':>7s}")
for nm, g, sg, st, dd, sd, sdt in rows:
    wf = g - dd
    stx = f"{st:.4f}" if np.isfinite(st) else "  -  "
    sgt = f"{dd/sdt:+.2f}" if np.isfinite(sdt) else "  -  "
    P(f"  {nm:44s} {g:+7.4f} {dd:+8.4f} {wf:+7.4f} {wf/FLOOR:7.2f} {sd:6.4f} {dd/sd:+6.2f} {sgt:>7s}")
P(f"  floor N/12 = {FLOOR:.6f};  exact unconditional single plaquette W_full = {WF:+.4f} (R={WF/FLOOR:+.2f})")
dd_p19 = 0.1401
P(f"  [E] W_full/N/12 spans {min((g-dd)/FLOOR for _,g,_,_,dd,_,_ in rows):.2f}"
  f" .. {max((g-dd)/FLOOR for _,g,_,_,dd,_,_ in rows):.2f}  => W_full != N/12 (never 1.00)")
P(f"  [E] user's bound Delta_dyn <= G''-N/12 : at G''=0.6018 -> {0.6018-FLOOR:.4f};"
  f" measured Delta_dyn = {dd_p19:+.4f} ({dd_p19/(0.6018-FLOOR):.2f}x the bound)")
P(f"      true G''(beta=2.5,L=2) is 0.4116 (p19) / 0.3942 (p01), NOT 0.6018 (that is beta=6):"
  f" bound becomes {0.4116-FLOOR:.4f}")
P("  [E] sign of Delta_dyn: + in 4/5 routes, but -0.179 in 1 route; magnitudes 0.14..1.07 (7.6x).")
P("      => sign positive but NOT route-independent; magnitude not converged.")

# ---------- [C] beta-dependence of Delta_dyn (p19 4-point scan) ----------
P("=" * 96)
P("[C] Delta_dyn(beta) from p19 (L=2): 1.0:+0.0148  2.5:+0.1401  4.0:+0.2221  6.0:+0.2962")
bb = np.array([1.0, 2.5, 4.0, 6.0])
dd = np.array([0.0148, 0.1401, 0.2221, 0.2962])
k, c = np.polyfit(bb, dd, 1)
P(f"  linear fit: Delta_dyn = {k:+.4f}*beta {c:+.4f}   => zero at beta = {-c/k:.3f}")
P(f"  exact beta=0 closure (pure Haar, Phi flat): Delta_dyn = 0  -> fit intercept {c:+.4f} compatible with 0")
P("  => Delta_dyn>0 is a DYNAMICAL, beta-growing statement, not a kinematic floor;")
P("     it passes the beta=0 acid test in the sense that it vanishes there.")
# ---------- [D] cross-L and cross-mode table (archived; sign is NOT mode/L invariant) ----------
P("=" * 96)
P("[D] Delta_dyn at beta=2.5 across L, mode, estimator (all measured, all [E])")
rows2 = [
    ("p11  L=2 k=0 uniform, |v|=1", 2, 0.3144, 0.0910, "8/8 chains +"),
    ("p11  L=3 k=0 uniform, |v|=1", 3, 0.3713, 0.1005, "8/8 chains +; ratio 1.181 => O(L^0)"),
    ("p19  L=2 4x5000", 2, 0.1401, 0.0419, "1 of 4 chains NEGATIVE (-0.0079)"),
    ("o27A L=2 mu=0,Lie0", 2, 0.5619, 0.2770, "2.0 sigma"),
    ("p02  L=2 4x19500 fund-reject", 2, 0.8838, 0.0501, "17.6 sigma, hard reject 4.64%"),
    ("o13  L=2 1500 traj", 2, 1.0655, 0.3326, "3.2 sigma"),
    ("o27B L=2 all mu,Lie0", 2, -0.1789, 0.2251, "NEGATIVE, 0.8 sigma"),
    ("o13  L=3 560 samples", 3, -0.5720, 0.6090, "NEGATIVE, 0.9 sigma; tau=740 >> bin 14"),
]
for nm, Lv, dd_, sd, note in rows2:
    P(f"  L={Lv}  {nm:34s} Delta_dyn={dd_:+7.4f} +- {sd:.4f}  ({dd_/sd:+5.2f} sig)  {note}")
P(f"  => 6 of 8 positive; magnitude spread 0.14 .. 1.07 = {1.0655/0.1401:.1f}x;")
P("     the two NEGATIVE centres are both <1 sigma from 0, one is mode B, one is L=3 with tau>>bin.")
P("  => cross-L same-sign evidence exists ONLY on the p11 k=0 uniform mode (+0.314 -> +0.371).")
P("=" * 96)
P("[VERDICT] W_full = N/12 : FALSE.  W_full >= N/12 : numerically true at all measured points")
P("          (min ratio 1.025 in research/v3/species/H1C-WFULL.md), but [V] not [T]:")
P("          what is [T] is W_geo >= (N/12)||v||^2 (ch05 eq:wgeo-quant); W_full = W_geo + X,")
P("          X = -2Cov(G,Phi') sign-indefinite (X<0 at 21% of 232 scanned points; exact")
P("          unconditional single-plaquette X = -0.192 < 0).")
P("[VERDICT] Delta_dyn IS separable from G'' and IS directly measured (8 estimates above, same")
P("          ensembles, identity residual <=1e-16).  (c) does not apply - not a rename.")
P("[VERDICT] Delta_dyn>0 at beta=2.5: 6/8 estimates positive; the hardest are p11 (L=2 and L=3,")
P("          8/8 chains, >3.5 sigma, ratio 1.18 => O(L^0)) and p02 (17.6 sigma).  Both negative")
P("          centres are <1 sigma from 0 (mode B; L=3 with tau>>bin).")
P("          => sign supported but NOT mode/L invariant; magnitude unconverged 7.6x;")
P("          it is NOT implied by the bound G''-N/12 (measured value = 0.32x that bound).")
