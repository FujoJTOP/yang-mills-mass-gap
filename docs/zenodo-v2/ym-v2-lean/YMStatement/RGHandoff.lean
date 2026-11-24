/-
# RG 交接骨架：由 RG 不变标度比 q(μ)=m(μ)/μ 控制弱→O(1) 耦合

本文件把 `research/physics_handoff.md` 第 2、3 节的数学内容形式化为一个**纯实分析骨架**。
约定 `t = ln μ`（RG time），所有对象都是 `ℝ → ℝ` 上的函数；`μ = exp t ∈ (0,∞)` 与 `t ∈ ℝ`
通过 `exp` 一一对应，因此"覆盖 `(0,∞)`"与"覆盖所有 `t ∈ ℝ`"等价。

包含六部分：
1. 一圈流 `g' = -b g³` 的显式解及其导数恒等式（`g²` 形式）；
2. RG 不变量 `Λ(t) = e^t · e^{-1/(2b g²)}` 的常数性与导数 `Λ'=0`；
3. 方案协变：`g'² = g² + c g⁴` 下 `Λ` 的相对因子精确表达，以及 `g²=1` 穿越时刻随 `c` 平移；
4. 质量与 `q = m/μ`：严格单调、`q(t₀)=1` 归一化与唯一穿越；
5. 交接胶合引理：IVT 存在唯一 + 集合覆盖 + `gap ≥ min(δ_W,δ_S)`；
6. 与 `Basic.lean` 的衔接说明（见文件末注释）。
-/

import Mathlib

noncomputable section

open Filter
open scoped Topology

namespace YMStatement.RG

/-! ## 1. 一圈流 `g' = -b g³` 的显式解 -/

/-- 一圈流的**耦合平方**显式解。设 `A = g₀²`，则
`g²(t) = A / (1 + 2 b A (t - t₀))`。该函数满足 ODE `(g²)' = -2b (g²)²`，
它正是 `g' = -b g³`（因 `(g²)' = 2 g g' = -2b g⁴ = -2b (g²)²`）的平方形式。 -/
def gSq (b A t0 : ℝ) : ℝ → ℝ := fun t => A / (1 + 2*b*A*(t - t0))

/-- 初值：`g²(t₀) = A = g₀²`。 -/
theorem gSq_at_t0 (b A t0 : ℝ) : gSq b A t0 t0 = A := by
  unfold gSq; simp

/-- 显式解满足一圈流的 ODE（`g²` 形式）：在分母非零处，
`d/dt g²(t) = -2b · (g²(t))²`。这是一条代数导数恒等式。 -/
theorem gSq_deriv (b A t0 t : ℝ)
    (hd : 1 + 2*b*A*(t - t0) ≠ 0) :
    deriv (gSq b A t0) t = -(2*b) * (gSq b A t0 t)^2 := by
  have h : HasDerivAt (fun s : ℝ => A / (1 + 2*b*A*(s - t0)))
      (-(2*b) * (A/(1 + 2*b*A*(t - t0)))^2) t := by
    have hnum : HasDerivAt (fun _ : ℝ => A) 0 t := hasDerivAt_const t A
    have hden : HasDerivAt (fun s : ℝ => 1 + 2*b*A*(s - t0)) (2*b*A) t := by
      have h1 : HasDerivAt (fun s : ℝ => s - t0) 1 t := (hasDerivAt_id t).sub_const t0
      have h2 : HasDerivAt (fun s : ℝ => (2*b*A) * (s - t0)) (2*b*A) t := by
        simpa using h1.const_mul (2*b*A)
      have h3 : HasDerivAt (fun s : ℝ => (2*b*A) * (s - t0) + 1) (2*b*A) t := by
        simpa using h2.const_add 1
      simpa [add_comm] using h3
    have hquot := hnum.div hden hd
    have hval : ((0 * (1 + 2*b*A*(t - t0)) - A * (2*b*A)) / (1 + 2*b*A*(t - t0))^2)
        = -(2*b) * (A/(1 + 2*b*A*(t - t0)))^2 := by
      field_simp
      ring
    rwa [hval] at hquot
  have hd' := h.deriv
  rw [show gSq b A t0 = (fun s : ℝ => A / (1 + 2*b*A*(s - t0))) from rfl]
  simpa using hd'

/-! ## 2. RG 不变量 `Λ(t) = e^t · e^{-1/(2 b g²)}` 与 `Λ' = 0` -/

/-- RG 不变量。`Λ(t) = exp(t - 1/(2 b g²(t)))`（与 `e^t · e^{-1/(2b g²)}` 相等）。 -/
def lam (b A t0 : ℝ) : ℝ → ℝ := fun t => Real.exp (t - 1/(2*b*(gSq b A t0 t)))

/-- **不变量核心事实**：由 1 的显式解，`Λ(t)` 与 `t` 无关；"两个时刻算出的 `Λ` 相等"。
显式地 `t - 1/(2b g²(t)) = t₀ - 1/(2b A)`。 -/
theorem lam_const (b A t0 t : ℝ) (hb : b ≠ 0) (hA : A ≠ 0)
    (_hd : 1 + 2*b*A*(t - t0) ≠ 0) :
    lam b A t0 t = lam b A t0 t0 := by
  unfold lam
  congr 1
  unfold gSq
  field_simp
  ring

/-- 由常数性得 `Λ' = 0`：`Λ` 在局部是常函数，故导数恒为零。 -/
theorem lam_deriv (b A t0 t : ℝ) (hb : b ≠ 0) (hA : A ≠ 0)
    (hd : 1 + 2*b*A*(t - t0) ≠ 0) :
    deriv (lam b A t0) t = 0 := by
  have hcont : ContinuousAt (fun s : ℝ => 1 + 2*b*A*(s - t0)) t :=
    (by fun_prop : Continuous (fun s : ℝ => 1 + 2*b*A*(s - t0))).continuousAt
  have hev : ∀ᶠ s in 𝓝 t, 1 + 2*b*A*(s - t0) ≠ 0 := hcont.eventually_ne hd
  have hfun : lam b A t0 =ᶠ[𝓝 t] (fun _ => lam b A t0 t0) := by
    filter_upwards [hev] with s hs
    exact lam_const b A t0 s hb hA hs
  rw [hfun.deriv_eq]
  exact deriv_const t (lam b A t0 t0)

/-! ## 3. 方案协变：`g'² = g² + c g⁴` -/

/-- 方案平移后的耦合平方：`g'² = g² + c g⁴ = g² (1 + c g²)`。 -/
def gSqShift (b A t0 c : ℝ) : ℝ → ℝ :=
  fun t => gSq b A t0 t * (1 + c * gSq b A t0 t)

/-- 方案平移后的不变量（用 `g'²` 构造）。 -/
def lamShift (b A t0 c : ℝ) : ℝ → ℝ :=
  fun t => Real.exp (t - 1/(2*b*(gSqShift b A t0 c t)))

/-- 纯代数恒等式（方案协变的核）：对 `u ≠ 0`、`1 + c u ≠ 0`，
`-1/(2b·u(1+cu)) = -1/(2b u) + c/(2b(1+cu))`。 -/
theorem exp_factor_identity (b c u : ℝ) (hb : b ≠ 0) (hu : u ≠ 0)
    (h1cu : 1 + c*u ≠ 0) :
    -(1 : ℝ)/(2*b*(u*(1 + c*u))) = -(1/(2*b*u)) + c/(2*b*(1 + c*u)) := by
  have h1cu' : 1 + u*c ≠ 0 := by simpa [mul_comm] using h1cu
  field_simp
  ring

/-- **方案协变的精确相对因子**：`Λ'(t) / Λ(t) = exp(c / (2b(1 + c g²(t))))`。
在一圈意义下，这正是"`Λ'` 与 `Λ` 差一个乘性因子"的精确形式；当 `g² → 0`（UV 弱耦合）
时该因子趋于常数 `e^{c/(2b)}`，即报告中的线性化 `Λ' = Λ · e^{c/(2b)}`。
（精确因子随 `g²` 缓变，故严格的"常数因子"只在 UV 极限/线性化下成立——诚实起见在此点明。） -/
theorem lamShift_ratio (b A t0 c t : ℝ) (hb : b ≠ 0)
    (hu : gSq b A t0 t ≠ 0) (h1c : 1 + c * gSq b A t0 t ≠ 0) :
    lamShift b A t0 c t
      = lam b A t0 t * Real.exp (c / (2*b*(1 + c*gSq b A t0 t))) := by
  unfold lamShift lam
  rw [← Real.exp_add]
  have hinner : t - 1/(2*b*(gSqShift b A t0 c t))
      = (t - 1/(2*b*(gSq b A t0 t))) + c/(2*b*(1 + c*gSq b A t0 t)) := by
    unfold gSqShift
    have h := exp_factor_identity b c (gSq b A t0 t) hb hu h1c
    linear_combination h
  rw [hinner]

/-- UV 极限（`g² → 0`）下相对因子退化为常数 `exp(c/(2b))`。 -/
theorem lamShift_uv (b c : ℝ) :
    c / (2*b*(1 + c*0)) = c/(2*b) := by
  simp

/-- **`g²=1` 穿越时刻随 `c` 平移**（代数事实）：在 `g² = 1` 的尺度上
`g'² = 1 + c·1⁴ = 1 + c ≠ 1`（当 `c ≠ 0`）。因此 `g'² = 1` 的穿越发生在不同的耦合值
（即 `x + c x² = 1` 的非平凡根），故"`g²=1` 的尺度"不是 scheme-不变的。 -/
theorem crossing_not_invariant (c : ℝ) (hc : c ≠ 0) : (1 : ℝ) + c * (1 : ℝ)^4 ≠ 1 := by
  simp only [one_pow, mul_one]
  intro h
  exact hc (by linarith)

/-! ## 4. 质量 `m(t)` 与不变量 `q(t) = m(t)/e^t = m/μ` -/

/-- 有效质量 `m(t) = m₀ e^{γ(t-t₀)}`，取归一化 `m₀ = e^{t₀}`（使 `q(t₀)=1`）。 -/
def mass (γ t0 : ℝ) : ℝ → ℝ := fun t => Real.exp t0 * Real.exp (γ*(t - t0))

/-- RG 不变标度比 `q(t) = m(t)/μ = mass(t)/e^t`。 -/
def qfun (γ t0 : ℝ) : ℝ → ℝ := fun t => mass γ t0 t / Real.exp t

/-- `q` 的指数闭式：`q(t) = exp(t₀ + γ(t-t₀) - t)`。 -/
theorem qfun_eq_exp (γ t0 t : ℝ) :
    qfun γ t0 t = Real.exp (t0 + γ*(t - t0) - t) := by
  unfold qfun mass
  rw [← Real.exp_add, ← Real.exp_sub]

/-- 归一化：`q(t₀) = 1`。 -/
theorem qfun_t0 (γ t0 : ℝ) : qfun γ t0 t0 = 1 := by
  rw [qfun_eq_exp]; simp

/-- **严格单调**：当 `0 < γ < 1` 时 `q` 关于 `t` 严格递减。
等价地：`t` 减小（IR 方向）时 `q` 增大。 -/
theorem qfun_strictAnti (γ t0 : ℝ) (hγ : γ < 1) : StrictAnti (qfun γ t0) := by
  intro a b hab
  rw [qfun_eq_exp, qfun_eq_exp]
  apply Real.exp_lt_exp.mpr
  nlinarith [mul_pos (sub_pos.mpr hab) (sub_pos.mpr hγ)]

/-- `q` 连续。 -/
theorem qfun_continuous (γ t0 : ℝ) : Continuous (qfun γ t0) := by
  have h : qfun γ t0 = fun t => Real.exp (t0 + γ*(t - t0) - t) := by
    funext t; exact qfun_eq_exp γ t0 t
  rw [h]
  fun_prop

/-- UV 极限：`q(t) → 0`（`t → +∞`）。 -/
theorem qfun_tendsto_zero (γ t0 : ℝ) (hγ : γ < 1) :
    Tendsto (qfun γ t0) atTop (𝓝 0) := by
  have hc : 0 < 1 - γ := by linarith
  have hinner : Tendsto (fun t : ℝ => t0 + γ*(t - t0) - t) atTop atBot := by
    rw [tendsto_atBot]
    intro a
    rw [eventually_atTop]
    refine ⟨((t0 - γ*t0) - a) / (1 - γ), ?_⟩
    intro t ht
    have hmul : (1 - γ) * (((t0 - γ*t0) - a) / (1 - γ)) ≤ (1 - γ) * t :=
      mul_le_mul_of_nonneg_left ht (le_of_lt hc)
    have heq : (1 - γ) * (((t0 - γ*t0) - a) / (1 - γ)) = (t0 - γ*t0) - a := by
      field_simp
    rw [heq] at hmul
    nlinarith [hmul]
  have key : Tendsto (fun t : ℝ => Real.exp (t0 + γ*(t - t0) - t)) atTop (𝓝 0) :=
    Real.tendsto_exp_atBot.comp hinner
  rw [show qfun γ t0 = (fun t => Real.exp (t0 + γ*(t - t0) - t)) from by
    funext t; exact qfun_eq_exp γ t0 t]
  exact key

/-- **`q = 1` 的唯一穿越**：在归一化下，`q(t) = 1` 当且仅当 `t = t₀`。
显式 `t* = t₀`（即 `μ* = e^{t₀} = m(t₀)`，与"交接点 = 质量隙标度"一致）。 -/
theorem qfun_unique_crossing (γ t0 : ℝ) (hγ : γ < 1) :
    ∃! t : ℝ, qfun γ t0 t = 1 := by
  refine ⟨t0, qfun_t0 γ t0, ?_⟩
  intro t ht
  rw [qfun_eq_exp] at ht
  have hzero : t0 + γ*(t - t0) - t = 0 := by
    apply Real.exp_injective
    rwa [Real.exp_zero]
  have hfac : (γ - 1) * (t - t0) = 0 := by linear_combination hzero
  have hγ1 : γ - 1 ≠ 0 := by linarith
  have ht0 : t - t0 = 0 := (mul_eq_zero.mp hfac).resolve_left hγ1
  linarith

/-! ## 5. 交接胶合引理（核心，纯实分析） -/

/-- **IVT：穿越点存在唯一**。设 `q` 连续、严格递减、`q(+∞) = 0`，且 `q₁ > 0`
落在 `q` 的值域内（某个 `t` 上 `q ≥ q₁`）。则 `q(t) = q₁` 有唯一解 `t*`。 -/
theorem crossing_exists_unique (q : ℝ → ℝ)
    (hcont : Continuous q) (hmono : StrictAnti q)
    (h0 : Tendsto q atTop (𝓝 0)) (q1 : ℝ) (hq1 : 0 < q1)
    (hEx : ∃ t, q1 ≤ q t) :
    ∃! t, q t = q1 := by
  obtain ⟨T, hT⟩ := eventually_atTop.mp
    (h0.eventually (Iio_mem_nhds hq1))
  have hTlt : q T < q1 := hT T le_rfl
  obtain ⟨t0, ht0⟩ := hEx
  have ht0T : t0 < T := by
    rcases lt_or_ge t0 T with h | h
    · exact h
    · exfalso
      have := hT t0 h
      linarith
  have hy : q1 ∈ Set.Icc (q T) (q t0) := ⟨le_of_lt hTlt, ht0⟩
  have himg : q1 ∈ q '' Set.Icc t0 T :=
    intermediate_value_Icc' ht0T.le hcont.continuousOn hy
  obtain ⟨c, _hc, hcEq⟩ := (Set.mem_image q (Set.Icc t0 T) q1).mp himg
  refine ⟨c, hcEq, ?_⟩
  intro t ht
  exact hmono.injective (ht.trans hcEq.symm)

/-- **集合覆盖**：弱域 `W = {t | g²(t) ≤ g*²}`（记 `{t | t ≥ t_up}`）与
强域 `S = {t | q(t) ≥ q₁}`（记 `{t | t ≤ t*}`）。若重叠条件 `t_up ≤ t*` 成立，
则每个尺度 `t ∈ ℝ` 都落在 `W ∪ S` 中（对应 μ = e^t ∈ (0,∞) 全被覆盖）。 -/
theorem coverage (q gc : ℝ → ℝ)
    (hmono_q : StrictAnti q) (hmono_g : StrictAnti gc)
    (tUp tStar gs q1 : ℝ)
    (hg_at : gc tUp = gs) (hq_at : q tStar = q1)
    (hoverlap : tUp ≤ tStar) :
    ∀ t : ℝ, gc t ≤ gs ∨ q t ≥ q1 := by
  intro t
  rcases lt_or_ge t tUp with h | h
  · right
    have h2 : t < tStar := lt_of_lt_of_le h hoverlap
    have : q tStar < q t := hmono_q h2
    linarith [hq_at]
  · left
    have hle : gc t ≤ gc tUp := by
      rcases lt_or_eq_of_le h with hlt | heq
      · exact le_of_lt (hmono_g hlt)
      · rw [heq]
    rwa [hg_at] at hle

/-- **每个尺度上的证书下界**（抽象陈述，`cert` 为给定函数，不构造测度）。

`cert t` 是**同一个固定数** `Δbar`（`Δ/Λ` 或 `√σ/Λ`）由尺度-`t` 分析给出的
**下界证书**，不是尺度 `t` 处的物理隙（后者会是 `m/μ`，在 `μ→∞` 时趋于 0，
写成 `δ_W ≤ cert t` 会不满足）。若弱域上 `cert ≥ δ_W`、强域上 `cert ≥ δ_S`，
且重叠条件成立，则所有尺度上 `cert t ≥ min(δ_W, δ_S)`。
证书到物理数的桥见 `certificate_transfer`。 -/
theorem gap_lower_bound (q gc cert : ℝ → ℝ)
    (hmono_q : StrictAnti q) (hmono_g : StrictAnti gc)
    (tUp tStar gs q1 : ℝ) (hg_at : gc tUp = gs) (hq_at : q tStar = q1)
    (hoverlap : tUp ≤ tStar)
    (δW δS : ℝ)
    (hW : ∀ t, gc t ≤ gs → cert t ≥ δW)
    (hS : ∀ t, q t ≥ q1 → cert t ≥ δS) :
    ∀ t, cert t ≥ min δW δS := by
  intro t
  rcases coverage q gc hmono_q hmono_g tUp tStar gs q1 hg_at hq_at hoverlap t with h | h
  · exact le_trans (min_le_left δW δS) (hW t h)
  · exact le_trans (min_le_right δW δS) (hS t h)

/-- **证书到物理数的桥**。若每个尺度上的证书都是固定数 `Δbar` 的下界
（这是 RG 的 gap-transfer 引理，非免费），则交接引理给出
`min(δ_W, δ_S) ≤ Δbar`。 -/
theorem certificate_transfer (cert : ℝ → ℝ) (Δbar δW δS : ℝ)
    (hbound : ∀ t, cert t ≥ min δW δS)
    (hcert : ∀ t, cert t ≤ Δbar) :
    min δW δS ≤ Δbar :=
  le_trans (hbound 0) (hcert 0)

/-! ## 6. 与 `Basic.lean` 的 `MassGapStatement` 如何衔接

`gap_lower_bound` 给出的是**证书的均匀正下界** `m_* := min δ_W δ_S > 0`（在所有尺度 `t`
上）；`certificate_transfer` 再把 `cert t ≤ Δbar` 转成 `Δbar ≥ m_*`。
`Basic.lean` 的 `MassGapStatement` 中，`ExponentialClustering μ m` 要求一个速率 `m > 0`。
本文件的交接引理提供的正是"该速率可由两侧拼接出的一个正数"这一**候选下界**；要真正得到
`MassGapStatement` 尚需：(i) 把 `cert`/`Δbar` 对应到 `ExponentialClustering` 的 `m`；(ii) 在两侧分别
证明 `δ_W, δ_S` 的定量界（H1,H2，属各侧分析）；(iii) 证明重叠条件 (H3) `q₀ ≥ q₁`——这正是公开难点。
故本文件**不是** `MassGapStatement` 的证明，只是其纯实分析骨架与坐标 `q = m/μ` 的严格化，
且不夸大：它把"交接"约化为一个 scheme-无关的**纯数不等式** `q₀ ≥ q₁`。
-/

end YMStatement.RG
