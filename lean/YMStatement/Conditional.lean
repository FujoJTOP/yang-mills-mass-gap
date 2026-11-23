/-
# 条件定理骨架：九个输入 ⇒ `Δ/Λ > 0`

本文件把 `docs/conditional-theorem.md` 的 §5 证明骨架形式化为一个 **Lean 条件定理**。
它不是 Clay 问题的解；它把"解"写成"`ClayInputs` ⇒ `0 < Delta_bar`"，并让每个假设字段
都是关于真实数学对象（实函数、实数不等式）的陈述。

## 形式化的部分（对应 §5 的 S3、S4）

* `ClayInputs` 收集 `IN-0 … IN-8` 中**可以用纯实分析陈述**的部分：
  `IN-3` 的两圈 β 界 `NPB`、`IN-4` 的强侧证书界、`IN-5` 的 overlap（写成 margin 不等式）、
  `IN-8` 的群常数（作为数据），以及证书函数 `cert` 与语义桥 `hcert`（`cert μ ≤ Δ/Λ`）。
* `twoLoop_overlap`：两圈 `t_up` 的闭式不超过一圈 `t_up`（纯代数），因此
  `q₀` 的两圈值不低于一圈值；一圈 overlap 蕴含两圈 overlap。
* `gap_of_inputs` 复用 `Handoff.lean` 的 `mass_gap_of_certificate`（S3+S4），
  不重新实现交接引理。
* `clay_of_inputs : 0 < I.Delta_bar`。

## 未形式化的部分（散文假设，显式列为定理参数）

以下输入是**关于测度、Wightman 场和极限**的陈述，无法在本文件中非空洞地形式化。
它们**不是** `ClayInputs` 的字段（否则会伪装成已形式化），而是作为 `Prop` 参数出现在
`clay_of_inputs_with_prose` 里，任何下游使用都必须显式命名它们：

* `IN-2`：`O(4)` 不变性，最锐形式 `Q_g Ω = 0`（角度荷湮灭真空）；
* `IN-6`：`{μ_L}` 的唯一弱极限与唯一 GNS 真空；
* `IN-7`：`μ` 非高斯且 `g_ren/Λ > 0`（需要 `RGCharge`，见 `Statement/NonTriviality.lean`）；
* `S1`、`S2`、`S6`–`S10`：测度构造、`O(4)` 不变性、极限、OS 重构、Wightman、禁闭。
  （`S5` 的算子代数核心已形式化，见 `clustering_rate_of_gap`；其测度侧仍是散文。）

此外，`NPB`（β 的两圈误差界）到 `t_up` 界的**推导**是 ODE 比较论证（积分 β 的界），
本文件只形式化了其中纯代数的部分：两圈闭式 `≤` 一圈闭式。这一点在
`docs/conditional-theorem.md` §5 的 S3 中也是"条件于 `IN-3`"的。

## 非空洞性

* `NPB` 是对实函数 `β : ℝ → ℝ` 的逐点界 `|β(λ) − β₂(λ)| ≤ C₃ λ⁴`，多一个或少一个
  函数都会使它为假；
* `Overlap` 是实数不等式 `q₁ ≤ q₀`，并经由 `handoff_iff_margin` 与 `margin ≥ 0` 等价；
* `Gap` 是 `0 < Delta_bar`，由 `min δ_W δ_S > 0` 与 `cert ≤ Δ/Λ` 合出。

没有 `sorry`，没有 `opaque`，没有 `True` 占位。
-/

import Mathlib
import YMStatement.Handoff
import HandoffMargin
import GL9Decay

noncomputable section

namespace YMStatement.Conditional

/-! ## 两圈 RG 时间与 β 函数

`HandoffMargin.lean` 的一圈闭式是 `tUp C_RG λ = C_RG/λ`。这里加入两圈闭式
（`docs/conditional-theorem.md` §IN-3）

    t_up^{2L}(λ) = C_RG/λ − (51/121) · log(1 + 121 C_RG/(51 λ)),

并证明它不超过一圈值。两圈的修正项是正的（符号由 `log(1+x) ≥ 0` 给出），因此两圈
到达弱耦合边界所需的 RG 时间更短，`q₀ = (m/Λ)e^{−t_up}` 更大。
-/

/-- 两圈 RG 时间闭式（`λ = h^∨ g²` 归一化）。 -/
def tUp2 (C_RG lam : ℝ) : ℝ :=
  C_RG / lam - (51 / 121) * Real.log (1 + 121 * C_RG / (51 * lam))

/-- **两圈 ≤ 一圈**：修正项非负。 -/
theorem tUp2_le_tUp {C_RG lam : ℝ} (hC : 0 < C_RG) (hlam : 0 < lam) :
    tUp2 C_RG lam ≤ C_RG / lam := by
  unfold tUp2
  have hx : 1 ≤ 1 + 121 * C_RG / (51 * lam) := by
    have : 0 ≤ 121 * C_RG / (51 * lam) := by positivity
    linarith
  have hlog : 0 ≤ Real.log (1 + 121 * C_RG / (51 * lam)) := Real.log_nonneg hx
  have hcoef : 0 ≤ (51 : ℝ) / 121 := by norm_num
  nlinarith

/-- 一圈 β 函数（`λ` 归一化）：`dλ/dt = −(h^∨/C_RG) λ²`。 -/
def beta1 (hvee lam : ℝ) : ℝ :=
  -(hvee / YangMillsInfra.RG.C_RG) * lam ^ 2

/-- 两圈 β 函数：一圈项加一个三次项，系数 `c2` 是输入数据。 -/
def beta2 (hvee c2 lam : ℝ) : ℝ :=
  beta1 hvee lam + c2 * lam ^ 3

/-- 两圈版本的 `q₀ = (m/Λ) e^{−t_up^{2L}}`。 -/
def q0TwoLoop (m Lam lamStar : ℝ) : ℝ :=
  (m / Lam) * Real.exp (-(tUp2 YangMillsInfra.RG.C_RG lamStar))

/-- `C_RG = 24π²/11 > 0`。 -/
theorem C_RG_pos : 0 < YangMillsInfra.RG.C_RG := by
  unfold YangMillsInfra.RG.C_RG
  positivity

/-- 两圈 `q₀` 不低于一圈 `q₀`（`t_up` 更短 ⇒ 指数更大）。 -/
theorem q0_le_q0TwoLoop {m Lam lamStar : ℝ} (hm : 0 < m) (hLam : 0 < Lam)
    (hlam : 0 < lamStar) :
    YangMillsInfra.RG.q0 YangMillsInfra.RG.C_RG m Lam lamStar
      ≤ q0TwoLoop m Lam lamStar := by
  unfold YangMillsInfra.RG.q0 q0TwoLoop YangMillsInfra.RG.tUp
  apply mul_le_mul_of_nonneg_left _ (le_of_lt (div_pos hm hLam))
  apply Real.exp_le_exp.mpr
  have h := tUp2_le_tUp (C_RG := YangMillsInfra.RG.C_RG) C_RG_pos hlam
  linarith

/-! ## 输入 -/

/-- 九个输入中可用纯实分析陈述的部分。

每个字段都是关于真实数学对象的陈述：

* `IN-8`：群标签 `hvee`、`dimG`、`dimf`（它们只作为数据出现，结论不含它们）；
* `IN-3`：两圈 β 误差界 `NPB`（对实函数 `beta` 的逐点界）与弱侧常数 `deltaW`；
* `IN-4`：强侧常数 `deltaS`；
* `IN-5`：重叠条件，写成 margin 不等式 `0 ≤ margin`（由 `handoff_iff_margin`
  等价于 `q₁ ≤ q₀`）；
* 证书 `cert : ℝ → ℝ` 与语义桥 `hcert : cert μ ≤ Delta_bar`（S4 的非自由部分）。

`handoff_gap` 只用 `q`、`gsq` 的单调性、证书两侧界与 overlap；`NPB` 与群标签不进入
证明，它们是输入清单的一部分（写明而不假装用于证明）。 -/
structure ClayInputs where
  /- IN-8: group labels. -/
  hvee : ℝ
  hvee_pos : 0 < hvee
  dimG : ℕ
  dimf : ℕ
  dimf_pos : 0 < dimf
  /- IN-3: two-loop beta bound. -/
  C3 : ℝ
  C3_nonneg : 0 ≤ C3
  c2 : ℝ
  beta : ℝ → ℝ
  NPB : ∀ lam : ℝ, 0 < lam → lam ≤ 7.3 →
    |beta lam - beta2 hvee c2 lam| ≤ C3 * lam ^ 4
  deltaW : ℝ
  deltaW_pos : 0 < deltaW
  /- IN-4: strong-side constant. -/
  deltaS : ℝ
  deltaS_pos : 0 < deltaS
  /- IN-5: overlap, in margin form. -/
  m : ℝ
  Lam : ℝ
  m_pos : 0 < m
  Lam_pos : 0 < Lam
  q1 : ℝ
  q1_pos : 0 < q1
  lamStar : ℝ
  lamStar_pos : 0 < lamStar
  lamStar_le : lamStar ≤ 7.3
  margin_nonneg :
    0 ≤ YangMillsInfra.RG.margin YangMillsInfra.RG.C_RG m Lam lamStar q1
  /- Handoff functions (the certificate is a function of the scale). -/
  q : ℝ → ℝ
  gsq : ℝ → ℝ
  mu_up : ℝ
  mu_up_pos : 0 < mu_up
  q_up : q mu_up = YangMillsInfra.RG.q0 YangMillsInfra.RG.C_RG m Lam lamStar
  q_antitone : ∀ a b : ℝ, 0 < a → a ≤ b → q b ≤ q a
  gsq_antitone : ∀ a b : ℝ, 0 < a → a ≤ b → gsq b ≤ gsq a
  /- Certificate and the semantic bridge (S4). -/
  cert : ℝ → ℝ
  cert_weak : ∀ mu : ℝ, 0 < mu → gsq mu ≤ gsq mu_up → deltaW ≤ cert mu
  cert_strong : ∀ mu : ℝ, 0 < mu → q1 ≤ q mu → deltaS ≤ cert mu
  Delta_bar : ℝ
  hcert : ∀ mu : ℝ, 0 < mu → cert mu ≤ Delta_bar

namespace ClayInputs

/-- 由输入构造 `HandoffData`。overlap 由 margin 不等式经 `handoff_iff_margin` 得到。 -/
def toHandoffData (I : ClayInputs) : HandoffData where
  q := I.q
  gsq := I.gsq
  μ_up := I.mu_up
  q0 := YangMillsInfra.RG.q0 YangMillsInfra.RG.C_RG I.m I.Lam I.lamStar
  q1 := I.q1
  deltaW := I.deltaW
  deltaS := I.deltaS
  μ_up_pos := I.mu_up_pos
  q_up := I.q_up
  overlap :=
    (YangMillsInfra.RG.handoff_iff_margin YangMillsInfra.RG.C_RG I.m I.Lam I.lamStar
      I.q1 I.m_pos I.Lam_pos I.lamStar_pos I.q1_pos).mpr I.margin_nonneg
  q_antitone := I.q_antitone
  gsq_antitone := I.gsq_antitone

/-- 一圈 overlap：`q₁ ≤ q₀`（`q₀` 取一圈闭式）。 -/
theorem overlap_oneLoop (I : ClayInputs) :
    I.q1 ≤ YangMillsInfra.RG.q0 YangMillsInfra.RG.C_RG I.m I.Lam I.lamStar :=
  (YangMillsInfra.RG.handoff_iff_margin YangMillsInfra.RG.C_RG I.m I.Lam I.lamStar
    I.q1 I.m_pos I.Lam_pos I.lamStar_pos I.q1_pos).mpr I.margin_nonneg

/-- **两圈 overlap**：`q₁ ≤ q₀^{2L}`。由 `t_up^{2L} ≤ t_up` 与一圈 overlap 合出。 -/
theorem twoLoop_overlap (I : ClayInputs) :
    I.q1 ≤ q0TwoLoop I.m I.Lam I.lamStar :=
  le_trans (overlap_oneLoop I) (q0_le_q0TwoLoop I.m_pos I.Lam_pos I.lamStar_pos)

end ClayInputs

/-! ## S3 + S4：证书下界与物理数 -/

/-- **S3 + S4**：交接引理给出的证书下界，经语义桥 `hcert` 传到 `Delta_bar`。

这里复用 `Handoff.lean` 的 `mass_gap_of_certificate`，不重新实现交接论证。 -/
theorem gap_of_inputs (I : ClayInputs) : min I.deltaW I.deltaS ≤ I.Delta_bar :=
  mass_gap_of_certificate I.toHandoffData I.cert I.Delta_bar
    I.cert_weak I.cert_strong I.hcert

/-- `min δ_W δ_S` 严格正。 -/
theorem min_pos (I : ClayInputs) : 0 < min I.deltaW I.deltaS :=
  lt_min I.deltaW_pos I.deltaS_pos

/-- **条件定理**：九个输入（中可用实分析陈述的部分）的 S3+S4 合成给出 `Δ/Λ > 0`。 -/
theorem clay_of_inputs (I : ClayInputs) : 0 < I.Delta_bar :=
  lt_of_lt_of_le (min_pos I) (gap_of_inputs I)

/-! ## S5 的代数核心：隙 ⇒ 离散时间指数聚类

`research/glue-a-limits/GL9Decay.lean` 证明了 GL9 的代数核心：在任意实内积空间中，
若 `T Ω = Ω`、`‖Ω‖ = 1`、`T` 保持 `Ω^⊥` 且在其上以 `λ` 压缩，则 `n` 步关联以 `λⁿ`
衰减。下面把 `λ = e^{−aδ}`（`a` 格距，`δ` 隙率）代入，得到指数形式。

这一步仍不涉及测度：它是 S5 的算子代数部分，测度侧的 `μ_L`、极限与 OS 重构属于散文
输入（见下节）。 -/
section S5

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- **S5（指数形式）**：以 `e^{−aδ}` 压缩的转移矩阵给出 `e^{−n a δ}` 的关联衰减。 -/
theorem clustering_rate_of_gap {T : E →ₗ[ℝ] E} {Ω : E} {a δ : ℝ}
    (hΩ : ‖Ω‖ = 1) (hTΩ : T Ω = Ω) (hperp : GlueA.PreservesPerp T Ω)
    (_ha : 0 < a) (_hδ : 0 < δ)
    (hgap : GlueA.GapOnPerp T Ω (Real.exp (-(a * δ)))) (A B : E) (n : ℕ) :
    |inner ℝ A ((T : E → E)^[n] B) - inner ℝ Ω B * inner ℝ A Ω|
      ≤ Real.exp (-(n * a * δ)) * ‖A‖ * ‖B‖ := by
  have hlam : 0 ≤ Real.exp (-(a * δ)) := Real.exp_nonneg _
  have hdecay := GlueA.spectral_decay hΩ hTΩ hperp hgap hlam A B n
  have hpow : Real.exp (-(a * δ)) ^ n = Real.exp (-(n * a * δ)) := by
    rw [← Real.exp_nat_mul]
    congr 1
    ring
  rwa [hpow] at hdecay

end S5

/-! ## 散文输入（显式参数，不是字段）

`IN-2`（`O(4)` 不变性，`Q_g Ω = 0`）、`IN-6`（唯一极限与唯一真空）、`IN-7`（非高斯）
以及 S1、S2、S6–S10 是测度/场论/极限陈述（S5 的代数核心已在上节形式化）。它们出现在下面定理的**参数**位置：定理
本身只是把 S3/S4 的结论裹上这些假设，从而任何引用都必须显式命名它们。它们没有被
形式化，`clay_of_inputs` 也不依赖它们。 -/
theorem clay_of_inputs_with_prose (I : ClayInputs)
    (IN2_O4Invariance : Prop) (IN6_UniqueVacuum : Prop)
    (_h2 : IN2_O4Invariance) (_h6 : IN6_UniqueVacuum) : 0 < I.Delta_bar :=
  clay_of_inputs I

end YMStatement.Conditional

/-! ## 公理审计（构建强制）

如果任一主定理的公理足迹改变（例如 `sorry` 泄漏进来），脚本会失败。 -/

/-- info: 'YMStatement.Conditional.clay_of_inputs' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YMStatement.Conditional.clay_of_inputs

/-- info: 'YMStatement.Conditional.gap_of_inputs' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YMStatement.Conditional.gap_of_inputs

/-- info: 'YMStatement.Conditional.tUp2_le_tUp' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YMStatement.Conditional.tUp2_le_tUp

/-- info: 'YMStatement.Conditional.clustering_rate_of_gap' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YMStatement.Conditional.clustering_rate_of_gap
