# -*- coding: utf-8 -*-
"""
m3c_tail.py — 道盲余项：把 (★-gap) 的第三条（尾部界）从 KP 机器里取出来并核对常数。

任务（L3 第三关）：`a_c` 的定义 + `a_c>=2 => m_c >= (1/2)ln(1/u)` 的推导 + "道盲"的含义
             + 这条尾部界是否已由 I1 的 KP 机器给出。

结论（数值部分；判词见 research/v3/species/M3C-TAIL.md）：
  T1  u 的两套归一化精确相差 2：u_P3 = c1/c0 = 2*u_T5（u_T5 = I2/I1）。
      质量式 m_c = a_c ln(1/u) 吃的是 u_T5（m1a=4ln(1/u_T5)=18.257 @1/48），
      KP 判据式 rho = lam_a*Phi 吃的是 u_P3。混用就会把窗口边读错 2 倍。
  T2  道盲余项 = ln(Cbar*e*Delta) = 4.689（只吃动物常数；零 J^{PC} 输入）。
      尾部界 m_c >= a_c*L - ln(Cbar*e*Delta)，L = ln(1/u_T5)。
  T3  该余项在**差式**里精确对消：m_c - m_{0++} >= (a_c-1)*L，零余项。
      => a_c>=2 给 m_c >= m_{0++} + L >= (1/2)L，无需任何 "C <= 2.691" 算术。
      => L3C §4/P3 的 "C <= 2.691" 是**幻影任务**（见 §4）。
  T4  但在窗边 delta 的数值依 Cbar 约定而变：Cbar=2 -> delta=ln2=0.693；
      双倍约定（UC 的 n=1 项 2*lam_a*u_P3）-> delta=0（窗边恰是 KP 关闭点）。
      故第三项必须写成 delta = ln(1/ubar)/a，**不能**写成 (1/2)ln(1/u)=2.691。
  T5  自由场酸测：u->1 时 ubar->1, delta->0，判据自行失效 => 非空洞。

无 git commit。确定性，无 MC，无下载。
"""
import json
import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import ive

OUT = {}


def u_t5(beta):
    """u_T5 = I2/I1 —— 质量式 m_c = a_c ln(1/u) 用的那套（T23/T5 归一化）。"""
    return ive(2, beta) / ive(1, beta)


def u_p3(beta):
    """u_P3 = c1/c0 = (I1-I3)/(I0-I2) —— KP 判据 rho=lam_a*Phi 用的那套（P3 归一化）。"""
    return (ive(1, beta) - ive(3, beta)) / (ive(0, beta) - ive(2, beta))


def Phi0(beta):
    """Phi(0,beta) = beta*e^beta/(2 I1) - 1 = beta/(2*ive1) - 1  （UC §1.2 (★)）"""
    return beta / (2.0 * ive(1, beta)) - 1.0


# 动物常数（语料值，逐字）：Delta=20 (F07:39/Z12:124)，Cbar=2 (Z12:124)，lam_a=e*Delta=54.366 (UC-COMPLETE:50)
DELTA = 20.0
CBAR = 2.0
LAM_A = math.e * DELTA                 # 54.3656...
CEFF = math.log(CBAR * math.e * DELTA)  # ln(Cbar*e*Delta) = ln(108.73) = 4.6890...

BETA_EDGE = 0.018394231   # L3C §4 的窗边（= 1/lam_a 的数值巧合，见 §1 T1）
BETA_SZZ = 1.0 / 48.0

print("=" * 78)
print("[0] 常数台账（逐字取自语料，非本轮新定）")
print(f"    Delta  = {DELTA:g}                     (F07:39 / Z12:124)")
print(f"    Cbar   = {CBAR:g}                     (Z12:124 / A02:60)")
print(f"    lam_a  = e*Delta = {LAM_A:.6f}       (UC-COMPLETE:50)")
print(f"    C_eff  = ln(Cbar*e*Delta) = {CEFF:.6f}   <-- 这就是“道盲余项”的数值")
print(f"    Cbar*e*Delta = {CBAR*math.e*DELTA:.6f}")
OUT["constants"] = dict(Delta=DELTA, Cbar=CBAR, lam_a=LAM_A, C_eff=CEFF)

print("=" * 78)
print("[1] u 的两套归一化：u_P3 = c1/c0 = 2*u_T5（精确因子 2）")
print("    beta        u_T5=2I2/I1?  u_T5=I2/I1      u_P3=c1/c0    ratio u_P3/u_T5")
for b in (1.0 / 48.0, 0.5, 1.0, 2.5):
    a, c = u_t5(b), u_p3(b)
    print(f"    {b:<10.6f}  -            {a:.10f}   {c:.10f}   {c/a:.10f}")
    OUT.setdefault("dict", []).append(dict(beta=b, u_t5=a, u_p3=c, ratio=c / a))
ratio_dev = max(abs(d["ratio"] - 2.0) for d in OUT["dict"])
print(f"    max |ratio-2| = {ratio_dev:.3e}   -> 精确字典 u_P3 = 2 u_T5  [T]+[V]")
OUT["ratio_dev"] = ratio_dev

print("=" * 78)
print("[2] 三个根（交叉验证 UC-COMPLETE §2.3 的台账，验证归一化字典没读错）")
b_p3_need = brentq(lambda b: u_p3(b) - 1.0 / (2.0 * LAM_A), 1e-4, 1.0)   # 必要条件 u_P3 < 1/(2 lam_a)
b_kp_true = brentq(lambda b: LAM_A * Phi0(b) - 1.0, 1e-4, 1.0)          # 真判据 lam_a*Phi(0)<1
b_tube = brentq(lambda b: u_p3(b) - 1.0, 1e-4, 3.0)                    # u_P3<1
print(f"    u_P3 = 1/(2*lam_a)  必要性根 beta = {b_p3_need:.9f}   [UC: 0.018394231]")
print(f"    lam_a*Phi(0,beta)=1 真判据根 beta = {b_kp_true:.9f}   [UC: 0.018268566]")
print(f"    u_P3 = 1            四壁根 beta = {b_tube:.9f}   [UC: 2.446918311]")
OUT["roots"] = dict(b_p3_need=b_p3_need, b_kp_true=b_kp_true, b_tube=b_tube)
print(f"    => 判据根比必要性根更严 {100*(1-b_kp_true/b_p3_need):.3f}%（L3C 用的窗边 0.018394 是**必要条件**根）")
b_ubar = brentq(lambda b: CBAR * math.e * DELTA * u_t5(b) - 1.0, 1e-4, 1.0)  # ubar=1 in u_T5 units
print(f"    ubar = Cbar*e*Delta*u_T5 = 1  =>  beta = {b_ubar:.9f}  （u_T5 归一化下的尾部关闭点）")
OUT["b_ubar"] = b_ubar

print("=" * 78)
print("[3] 道盲余项与尾部界：L = ln(1/u_T5), delta = L - C_eff = ln(1/ubar)")
print("    beta        u_T5         L=ln(1/u)   ubar        delta       delta/L")
grid = [1e-3, 2e-3, 5e-3, 0.01, 0.018394231, 1.0 / 48.0, 0.0368, 2.5, 10.0, 100.0]
rows = []
for b in grid:
    u = u_t5(b)
    L = math.log(1.0 / u)
    ubar = CBAR * math.e * DELTA * u
    d = L - CEFF
    d_chk = math.log(1.0 / ubar)
    rows.append(dict(beta=b, u_t5=u, L=L, ubar=ubar, delta=d, delta_over_L=d / L, chk=d_chk))
    print(f"    {b:<10.6f}  {u:.6e}  {L:10.6f}   {ubar:10.6f}  {d:+10.6f}  {100*d/L:+8.2f}%")
print(f"    max |delta - ln(1/ubar)| = {max(abs(r['delta']-r['chk']) for r in rows):.3e}  -> 恒等式 [T]")
OUT["tail_grid"] = rows

print("=" * 78)
print("[4] 三个算术靶：L3C 的 (1/2)L 靶 vs 带 a_c>=2 的正确靶 vs 纯正性靶")
L_edge = math.log(1.0 / u_t5(BETA_EDGE))
print(f"    窗边 beta={BETA_EDGE}: L = {L_edge:.6f}, C_eff = {CEFF:.6f}")
print(f"    靶 A（L3C §4/P3 写的）: C <= (1/2)L = {0.5*L_edge:.6f}   -> C_eff 通过？ {CEFF <= 0.5*L_edge}"
      f"   (差 {CEFF-0.5*L_edge:+.6f})")
print(f"    靶 B（a_c>=2 的正确靶）: C <= (3/2)L = {1.5*L_edge:.6f}   -> 通过？ {CEFF <= 1.5*L_edge}")
print(f"    靶 C（只要 delta>0）  : C <  L      = {L_edge:.6f}   -> 通过？ {CEFF < L_edge}"
      f"   余量 {100*(1-CEFF/L_edge):.2f}%")
print("    ==> 靶 A 通不过（C_eff=4.689 > 2.691），但靶 A 用不着：见 §5 差式对消。")
OUT["targets"] = dict(L_edge=L_edge, target_halfL=0.5 * L_edge, target_1p5L=1.5 * L_edge,
                      pass_halfL=CEFF <= 0.5 * L_edge, pass_1p5L=CEFF <= 1.5 * L_edge,
                      pass_L=CEFF < L_edge)

print("=" * 78)
print("[5] 差式对消：余项 c-盲 => m_c - m_{0++} >= (a_c-1)*L（零余项）")
print("    a_c   m_c >= a_c*L - C_eff   m_c - m_0 >= (a_c-1)*L   (a_c-1)*L >= (1/2)L ?")
for ac in (1, 2, 3):
    absb = ac * L_edge - CEFF
    diffb = (ac - 1) * L_edge
    print(f"    {ac}     {absb:+18.6f}   {diffb:+18.6f}            {diffb >= 0.5*L_edge}")
OUT["diff_form"] = dict(L_edge=L_edge, halfL=0.5 * L_edge,
                        rows=[dict(a_c=ac,
                                   abs_bound=ac * L_edge - CEFF,
                                   diff_bound=(ac - 1) * L_edge) for ac in (1, 2, 3)])

print("=" * 78)
print("[6] (★-gap) 的第三项：delta，不是 (1/2)ln(1/u)")
print(f"    L3C 写的第三项 : (1/2)L          = {0.5*L_edge:.6f}")
print(f"    机器实际交付   : delta = ln(1/ubar) = {L_edge-CEFF:.6f}")
print(f"    比             : (1/2)L / delta  = {0.5*L_edge/(L_edge-CEFF):.3f}x")
print(f"    Cbar 双倍约定下 : delta = L - ln(2*Cbar*e*Delta) = {L_edge - math.log(2*CBAR*math.e*DELTA):+.6f}")
print("    => 窗边 delta in {0, ln2}：约定决定窗是“恰好关闭”还是“留 2 倍余量”。")
OUT["third_term"] = dict(halfL=0.5 * L_edge, delta=L_edge - CEFF,
                         ratio=0.5 * L_edge / (L_edge - CEFF),
                         delta_doubled_convention=L_edge - math.log(2 * CBAR * math.e * DELTA))

print("=" * 78)
print("[7] 道盲审计：每个常数吃什么（逐条；这是“道盲”的可核形式）")
audit = [
    ("kappa/4 = ln(1/u)", "u = I2/I1(β)", "NO", "单 plaquette 的 4 条边 + Wilson 单片成本"),
    ("w(B_T) = 4", "源 plaquette 的 4 条边", "NO", "几何"),
    ("u = c_f/N", "β", "NO", "Wilson 测度数据"),
    ("C_eff = ln(Cbar*e*Delta)", "Delta=20, Cbar=2", "NO", "KP 动物常数"),
    ("B, ubar", "w=4, u, KP 盘", "NO", "Y04:137 开管"),
    ("delta = L - C_eff", "以上四件", "NO", "=> 同一数对全道"),
    ("a_c", "管的每单位时间 plaquette 数", "YES", "唯一吃道的量，且只以乘子出现"),
]
print("    量                          输入                             吃 J^{PC}？  备注")
for name, inp, eats, note in audit:
    print(f"    {name:<27} {inp:<32} {eats:<12} {note}")
OUT["audit"] = [dict(q=n, inputs=i, eats_jpc=e, note=t) for n, i, e, t in audit]

print("=" * 78)
print("[8] 自由场酸测：u -> 1 时 ubar -> 1, delta -> 0（判据自行失效 => 非空洞）")
for b in (2.5, 10.0, 100.0):
    u = u_t5(b)
    d = math.log(1.0 / u) - CEFF
    print(f"    beta={b:<8g} u_T5={u:.6f}  ubar={CBAR*math.e*DELTA*u:.3f}  delta={d:+.6f}  {'空' if d<0 else '非空'}")
print("    自由无质量场 m=0 => 尾部界给 m_c >= 负/零 => 不存在假定理；YM 上给 delta>0 => 有内容。")

print("=" * 78)
print("[9] 判词（数值侧；文字判见 M3C-TAIL.md）")
print("    (a) 尾部界**已由** I1 的五环 KP 机器（环2/环3）给出：delta = ln(1/ubar)/a。")
print("        其 c-盲性 = C_eff 吃 (Delta,Cbar) 不吃 J^{PC}；差式里对消 => 零余项。")
print("    (b) 不需要新论证；但 L3C 的第三项 (1/2)ln(1/u) 必须改写成 delta（差 3.9x）。")
print("    (c) a_c 的定义不清：它由**管/算符**定义，不是谱道的“级”；m_c=a_c L+O(1) 是引入的。")
print("    窗：beta < min(1/48, 0.018269[判据根], 0.0368[尾部根])；L3 的 0.018394 是必要条件根。")

with open("m3c_tail.json", "w", encoding="utf-8") as f:
    json.dump(OUT, f, indent=1, ensure_ascii=False)
print("\nexit 0 ; witness -> w7/m3c_tail.json")
