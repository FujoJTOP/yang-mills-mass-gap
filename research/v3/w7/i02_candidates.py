# -*- coding: utf-8 -*-
"""i02_candidates.py — beta=2.5 上六个候选常数 c（f1^2 >= c）的数值。
只读语料数字，不做新格点作业。标签: [T] 严格 / [C] 条件 / [E] 实测。
锚: RUN3.md:11 (C(0)=19.888, w1=0.3427+-0.0119, L=8 raw, beta=2.5, m1=1.546) [E]
    T1-POLYMER-RG.md:164,352 (u(2.5)=1.01439) [T]+[V]; :396 (KP 阈值 u<0.01839)
    T2-TUBE-OVERLAP.md:329 (u(2.5)=0.719，另一口径) [T]
    A01-CBAR.md:39 (eDelta=lambda_a=20e=54.3656, ubar=0.5654 @ beta_W=1/12)
    u14_w1L.log (beta=2.5 L=4 fresh: M0=2.474, f1^2=1.309, w1=0.529) [E]
"""
import math
from scipy.special import ive

N = 2
BETA = 2.5
EDELTA = 20.0 * math.e
A0 = (2 * N) ** -2.0          # 单片下界 a0 >= 1/16 (N=2)
R0_MAG = 0.5 * A0             # = 1/32

# ---------- 0. 口径核验: SU(2) I2/I1 ----------
b_rat = lambda x: ive(2, x) / ive(1, x)
u_conv = {"I2(2b)/I1(2b)": b_rat(2 * BETA), "I2(b)/I1(b)": b_rat(BETA),
          "I2(b/2)/I1(b/2)": b_rat(BETA / 2)}
U_CORPUS = 1.01439            # T1
U_T2 = 0.719                  # T2
UALL = [("corpus-T1", U_CORPUS), ("corpus-T2", U_T2)] + list(u_conv.items())

print("=" * 72)
print("0. u(beta=2.5) 口径核验  [T]")
for k, v in UALL:
    print(f"   u[{k:16s}] = {v:.6f}   uK=u*eD*Cbar(Cbar=1) = {v*EDELTA:8.3f}   "
          f"{'KPOUT' if v > 0.018394 else 'in'}")
print(f"   KP 阈值 1/(eDelta) = {1/EDELTA:.6f}  (=0.01839, T1:396)")

# ---------- 1. M0 / f1^2 ----------
M0 = 19.888; W1 = 0.3427; DW1 = 0.0119
F1 = W1 * M0
DF1 = DW1 * M0
print("=" * 72)
print("1. M0(beta=2.5) 与 f1^2 = w1*M0   [E]")
print(f"   (L=8, raw)      M0=C(0)={M0}   w1={W1}+-{DW1}  =>  f1^2 = {F1:.4f} +- {DF1:.4f}")
for tag, m, w in [("L=4 fresh", 2.474, 0.529), ("L=8 APE1", 45.527, 0.0168),
                  ("L=8 TAPE2", 39.428, 0.0107), ("L=8 w1(A/C0)", 19.888, 0.4058)]:
    print(f"   ({tag:12s}) M0={m:8.3f}  w1={w:.4f}  f1^2={m*w:8.4f}")

# ---------- 2. 六候选 ----------
print("=" * 72)
print("2. 六候选在 beta=2.5 的值")
# (1)
print(f" (1) c = f1^2 实测 = w1*M0          = {F1:.4f}  (+-{DF1:.3f})  [E]  L=4 口径 {2.474*0.529:.4f}")
# (2)
print(f" (2) c = |r0| = 1/32                = {R0_MAG:.6f}  (正数, 但 signed r0 在下)")
for k, u in UALL:
    r0 = R0_MAG * (1.0 - u)          # C_corr=1
    print(f"      signed r0(2.5)[{k:16s}, C=1] = {r0:+.6f}    (C=2: {R0_MAG*(1-2*u):+.6f})")
# (3)
c3 = {}
c3["c'=1"] = math.exp(-1.0 * BETA)
c3["锚 r0(1/48)=0.005"] = math.exp(-(math.log(1 / 0.005) / (1 / 48.0)) * BETA)
c3["c'=ln(1/u_corpus)"] = math.exp(-math.log(1 / U_CORPUS) * BETA)
for k, v in c3.items():
    print(f" (3) c=e^(-c'b): {k:20s} = {v:.6e}")
# (4)
for n in (1.0, 1.26, 2.0):
    print(f" (4) c=c0/b^n: c0=1,n={n:<4}      = {BETA**-n:.6f}   c0=r0(1/48): {0.005*BETA**-n:.6f}")
# (5)
print(f" (5) c = 1/(eDelta)                 = {1/EDELTA:.6f}   (前提 u<0.018394 在 2.5: "
      f"超 {U_CORPUS/(1/EDELTA):.0f}x / {U_T2/(1/EDELTA):.0f}x)")
# (6)
for k, u in UALL:
    print(f" (6) ubar(2.5)[{k:16s}] Cbar=1: {EDELTA*u:8.2f}   Cbar=2: {EDELTA*u*2:8.2f}")
print(f"      u(1/48)=0.0104 -> u(2.5)=1.01439 : x{U_CORPUS/0.0104:.0f} (单调增, ubar 无界)")

# ---------- 3. 相容性判 ----------
print("=" * 72)
print("3. 相容性: c <= f1^2_meas ?  (+ 是否 2.5 上有推导)")
print(f"   f1^2_meas = {F1:.4f};  1/32={R0_MAG:.4f}(x{F1/R0_MAG:.0f} 倍偏小);  "
      f"1/(eD)={1/EDELTA:.4f}(x{F1*EDELTA:.0f})")
print(f"   M0 离散度 (L=4 2.474 / APE1 45.527) = x{45.527/2.474:.1f}  => 绝对 c 无 scheme/L 锚 [T]")
print("=" * 72)
