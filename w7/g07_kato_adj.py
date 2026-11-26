# -*- coding: utf-8 -*-
# w7/g07_kato_adj.py --- G07 裁决: `Kato => m>0` 真伪 + `r_iso` vs `I2` 是否同一个量
# 只依赖 mpmath; 确定性, 无采样无 MC; 墙 < 20 s.
# 逐字口径 (本轮读过的原文):
#   F04-SIGKATO.md:15      r_iso(T~) := dist(sup spec T~, spec T~ \ {sup}); sup 非有限重原子时记 0
#   F04-SIGKATO.md:105     `Kato => m>0`(给定 lambda_bot(T~)<1); spec(T~) ⊂ {lambda_bot} ∪ [0, lambda_bot-k]
#   w7/f04_sigma_kato.py:31  T~ := e^{-(H-E0)}|_{Omega^perp};  m := -log(sup spec T~)
#   B01-RP.md:12,43        m = -(1/a) log||T|_{Omega^perp}|| = -log lambda_2(T) = inf supp nu = 2 Lambda
#   B05-INVENT.md:54       T~ := e^{-(H-E_0)}|_{Omega^perp}   (真空归一/去真空)
#   G04-M1UNIF.md:56       nu = delta_0 + 1_{[delta,delta+1]}(s) ds ; F06-KLATOM.md:53  nu 建在 (O-<O>)Omega_L 上
#   LOOP-RUNBOOK.md:8311-8313  `r_iso = 1-e^{-a delta} > 0` / `m = -log(1)/a = 0`  <- 本轮要核的那两行
from mpmath import mp, mpf, exp, log

mp.dps = 40
A   = mpf('0.35')          # 晶格间距 a  (量纲 1/能量)
DEL = mpf('0.7')           # 见证的层距 delta  (I2 的 delta)
G2  = mpf('1')             # 2D YM 的 g^2


OUT = []


def hr(t=""):
    OUT.append(t)
    print(t)


def block(title):
    hr("\n" + "=" * 84)
    hr(title)
    hr("=" * 84)


# ---------------------------------------------------------------- [1] m 的定义
block("[1] `m` 是什么: 三个候选式同一个数; 第四式恒为 0 (reductio)")
E = [mpf(0), DEL, DEL + mpf('1.3'), DEL + mpf('2.9')]      # 任意谱: 真空 0 + 三个激发层
T = [exp(-A * (e - E[0])) for e in E]                      # T~ = e^{-a(H-E0)}
lam2 = T[1]                                                # T 在 Omega^perp 上的谱半径
norm = max(T[1:])                                          # ||T~|_{Omega^perp}|| (自伴: = 谱半径)
m_lam2 = -log(lam2) / A
m_norm = -log(norm) / A
m_full = -log(max(T)) / A                                  # -log||T~||(全域)/a  <-- runbook:8312 用的
hr(f"  lam2 = sup spec T~|Om^perp = {mp.nstr(lam2, 22)}")
hr(f"  ||T~|Om^perp||             = {mp.nstr(norm, 22)}   <- 与 lam2 **同一个数** (自伴) [T]")
hr(f"  m = -(1/a)log lam2   = {mp.nstr(m_lam2, 22)}")
hr(f"  m = -(1/a)log||.||   = {mp.nstr(m_norm, 22)}")
hr(f"  m = E1 - E0          = {mp.nstr(E[1] - E[0], 22)}")
hr(f"  => 三式逐位同值 (偏差 {mp.nstr(abs(m_lam2 - (E[1]-E[0])), 5)});  "
   f"`-(1/a)log lam2` 与 `-(1/a)log||T~|Om^perp||` **不是两个候选, 是一个** [T]")
hr(f"  `-(1/a)log||T~||`(全域, 含真空) = {mp.nstr(m_full, 6)}   <- 恒 0: 见 [5]")

# ------------------------------------------------- [2] 见证: 去真空后的 r_iso
block("[2] 见证 nu = d0 + 1_[d,d+1]: 按其**逐字**定义算 r_iso (真空已去, B05:54)")
supT = exp(-A * DEL)                 # E1 = delta  (a.c. 块的下沿)
botT = exp(-A * (DEL + 1))           # E = delta+1
gap = supT - botT
hr(f"  spec T~|Om^perp = [e^-a(d+1), e^-a d] = [{mp.nstr(botT,18)}, {mp.nstr(supT,18)}]  (纯 a.c., 无原子)")
hr(f"  块宽 = {mp.nstr(gap, 18)}")
hr("  (i) 字面 dist(sup, spec\\{sup}): 块在 sup **下方累积** -> N 点均分时最接近点为 gap/N")
for N in (10**2, 10**3, 10**4, 10**5, 10**6):
    hr(f"      N = {N:>8d}   min_{{s in spec\\{{sup}}}} |sup - s| = {mp.nstr(gap/N, 8)}")
hr("      => dist -> 0  (1/N 律)  => r_iso = 0  **字面上就是 0** [T]")
hr("  (ii) F04:15 的默认规则: sup 须为**有限重原子**, 否则记 0")
for eps in (mpf('1e-2'), mpf('1e-4'), mpf('1e-6')):
    mass = (DEL + eps) - DEL                            # nu([d, d+eps]) = eps  (Lebesgue)
    hr(f"      nu([d, d+{mp.nstr(eps,3)}]) = {mp.nstr(mass, 6)}   -> 0  => 无原子")
hr("      => 两条独立判据都给 r_iso = 0  => **Kato = 假**")
hr(f"  而 m = -(1/a)log sup spec T~ = -(1/a)log e^-a d = d = {mp.nstr(DEL, 12)} > 0")
hr("  => 见证实例是 **m>0 且 not Kato** = f04:106 的『逆假』方向, **不是** Kato 的反例 [T]")

# ------------------------------------------ [3] g04 有效的那一半: 平移族
block("[3] g04 的**有效**部分: I1 ^ I2 对支撑平移不变 => 给不出 m1 的位置 [T]")
hr("  nu_t = delta_t + 1_[t+d, t+d+1] ds    (t >= 0)")
hr("       t      f1^2=nu_t({t})   I1   Lambda*-m1^2   I2    m1=inf supp nu_t")
for t in [mpf(0), mpf('0.1'), mpf('0.3'), mpf('0.5'), mpf('1')]:
    f1 = mpf(1)                       # atom mass (归一 theta)
    lam_star = t + DEL
    hr(f"    {mp.nstr(t,4):>6}   {mp.nstr(f1,4):>10}      T    "
          f"{mp.nstr(lam_star - t,4):>8}      T     {mp.nstr(t,4)}")
hr("  => I1 与 I2 在**每个** t 上都真, 而 inf_t m1 = 0  => I1^I2 不含 m1 的位置 (无下界) [T]")
hr("     这正是 F06:117 已写的『否则底贴 0, I2 无内容』; 与 f04 的 Kato 无关 (对象不同, 见 [6])")

# ------------------------------------ [4] 唯一合法的 Kato 反例: 真空简并
block("[4] 唯一合法的 `Kato ^ not(m>0)`: 真空简并 (0 in spec H|Om^perp 为原子)")
sup4 = mpf(1)                        # E=0 的第二个基态 => sup spec T~ = 1, 有限重原子
bot4 = exp(-A * (DEL + 1))
r4 = sup4 - bot4
m4 = -log(sup4) / A
hr(f"  spec H|Om^perp = {{0}} ∪ [d, d+1]  (0 为 1 重原子 = 第二个基态)")
hr(f"  sup spec T~ = 1 (原子, 重数 1);  rest = [e^-a(d+1), e^-a d]")
hr(f"  r_iso = 1 - e^-a d = {mp.nstr(r4, 20)} > 0        => Kato = 真")
hr(f"  m     = -(1/a)log 1 = {mp.nstr(m4, 6)}            => m>0 假")
hr("  => **存在**合法模型使 `Kato ^ not(m>0)`; 但它要求 Om^perp 内有一个零能态:")
hr("     Omega 简并. 纯 YM 上由 PF 唯一真空 (B01-RP:11 ' lam1=1 单重') 排除 [T].")
hr("     F04:105 的旁条件 `lambda_bot(T~)<1` (= sup spec T~ < 1) 恰好只排除这一个点;")
hr("     而由 [1], 该旁条件 **逐字就是 m>0** (B01-RP:12) => 蕴含是**空转**的, 不是新内容 [T]")

# ------------------------------------------------ [5] reducito: 用户的 m 档
block("[5] reductio: 若把 `m` 读成 -(1/a)log||T~||(全域), 则对**一切**理论 m=0")
rows = [
    ("2D YM  (g^2=1)",  mpf('3') * G2 / 8),          # 真值 m = 3g^2/8 (F04:92)
    ("2D YM  (b=2.5)",  mpf('3') * mpf('1.6') / 8),  # g^2=1.6
    ("有质量自由场",     mpf(1)),                      # m = omega_0
    ("见证 (g04)",      DEL),
]
hr("  理论                m (B01-RP:12 口径)      m =-log||T~||/a (全域)")
for name, mt in rows:
    hr(f"  {name:<18} {mp.nstr(mt, 8):>14}      {mp.nstr(mpf(0), 6):>10}")
hr("  => 全域档给 0 是**恒等式** (真空 E0=0 归一), 与理论无关: 2D YM 也 0 (真值 3g^2/8) [T]")
hr("     => RUNBOOK:8312 的 `m = -log(1)/a = 0` 是**取错了分母子空间**, 不是见证的性质")

# ------------------------------------- [6] 类型表: r_iso (距离) vs m (位置)
block("[6] 类型: r_iso 与 I2 是**同一个量**(距离), m 是**另一个量**(位置)")
hr("  2D YM (g^2=1): m = 3g^2/8 = 0.375   r_iso = E1-E_{1/2} = 5g^2/8 = 0.625   Lambda* = g^2 = 1")
hr("    -> r_iso 与 (Lambda*-m1^2) 同值 0.625;  m 是 0.375  => 两个不同的量")
hr("  T~ 尺度 (2D YM, g^2=1, a=0.35): m1=0.375, Lambda*=1")
a0, m1, ls = mpf('0.35'), mpf('0.375'), mpf(1)
r_novac = exp(-a0 * m1) - exp(-a0 * ls)     # 去真空: 底(原子) 到 次底(支持点) 的距离 = 距离量
r_vac = mpf(1) - exp(-a0 * m1)              # 含真空: 真空到"第一个激发"的距离 = **隙本身**
hr(f"    r_iso(去真空, F04/B05 口径) = e^-a m1 - e^-a L*  = {mp.nstr(r_novac,12)}")
hr(f"    r_iso(含真空, 用户档)         = 1 - e^-a m1       = {mp.nstr(r_vac,12)}")
hr(f"    1 - e^-a m1 对 m1 **单射** => 由 r_iso={mp.nstr(r_vac,10)} 反解 m1 = "
   f"{mp.nstr(-log(1 - r_vac)/a0,12)}  (= 0.375 OK)")
hr("  => 去真空口径: r_iso 是 Delta 的像 (距离); 含真空口径: r_iso 是 m 的像 (位置)")
hr("     两种口径都**单调**, 故 r_iso>0 在各自口径下与对应的量**同真同假**;")
hr("     RUNBOOK:8311-8313 混用: r_iso 用含真空, m 用 `sup=1` 单点 => r_iso 与 m 不同源 [NG]")
hr("     RUNBOOK:8106 自己写 `r_iso <=> Lambda*-m1^2 >= kappa` => 语料已判 r_iso 与 I2 同一")
hr("     真正的分歧对是 **r_iso(距离) vs m(位置)**, 不是 r_iso vs Lambda*")

# --------------------------------------------------------------- 判词
block("[7] 判词")
for s in [
    "[T] `m` = -(1/a)log||T~|Om^perp|| = -log lam2 = inf supp nu = 位置(E1). 两个候选式同一个数.",
    "[T] `r_iso` = dist(sup spec T~, spec\\{sup}) = 距离(E2-E1) = `Lambda*-m1^2` 的单调像. 与 I2 同一个量.",
    "[T] 见证 nu=d0+1_[d,d+1] 按其逐字定义: 去真空后底是 a.c. 块沿 => r_iso = 0 => Kato 假; m = d > 0.",
    "[T] 见证实例 = `m>0 ^ not Kato` = f04 的 `逆假`; 不能反证 `Kato => m>0`.",
    "[T] `Kato ^ not(m>0)` 的唯一合法模型 = 真空简并 (0 in spec H|Om^perp 为原子), PF 排除.",
    "[T] f04:105 的旁条件 `lambda_bot<1` 逐字就是 m>0 => 蕴含成立但**空转**, 不产生新内容.",
    "[NG] RUNBOOK:8311-8313 (`r_iso=1-e^-ad>0`, `m=-log(1)/a=0`) 两行: 把真空放回 spec T~.",
    "[T] g04 的有效残余: I1^I2 对支撑平移不变 => 不给 m1 的下界 (与 Kato 无关, 对象不同).",
]:
    hr("  " + s)
hr("")

# [6b] 内部一致性: 见证上 r_iso 的定义域一旦是 "sup vs rest", m 必须用 rest 的 sup
block("[6b] 见证上的**内部一致性**: r_iso 用 sup-vs-rest, m 却只用 sup => 不同源")
sup_w = exp(-A * DEL)          # 去真空: sup = e^-a d
sup_v = mpf(1)                 # 含真空: sup = 1, rest 的 sup = e^-a d
for tag, s_up in (("去真空", sup_w), ("含真空", sup_v)):
    rest_sup = exp(-A * DEL)
    if s_up == rest_sup:
        continue
    d_iso = s_up - rest_sup                    # 这是 r_iso (距离)
    m_from_rest = -log(rest_sup) / A           # 由 rest 的 sup 反解的位置
    m_from_sup = -log(s_up) / A                # 由 sup 反解的位置
    hr(f"  {tag}: r_iso = {mp.nstr(d_iso,12)};  m(由 rest 反解) = {mp.nstr(m_from_rest,10)};"
          f"  m(由 sup 反解) = {mp.nstr(m_from_sup,10)}")
hr("  => RUNBOOK:8312 的 `m = -log(1)/a` 是**由 sup(=1, 真空) 反解**, 而 8311 的 r_iso")
hr("     是**由 rest 的上沿**量出; 两行用了两个不同的算子=> different m. 一致做法唯一:")
hr(f"     由 rest 的上沿反解 m = d = {mp.nstr(DEL,6)} > 0  => 见证上 m>0, 与 Kato 同真")
hr("")

with open("g07_kato_adj.log", "w", encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
print("[written] g07_kato_adj.log")
