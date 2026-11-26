/-
  GW1-bridge.lean
  ===============
  第 216 轮（GW1）：单侧 Grönwall 桥引理 + D04-ODE 命题 2.1 的等价性半边 + 加性（Hellmann–Feynman）版。

  三个定理 + 一个助手，全部 `sorry = 0`（文末 `#print axioms` 自检）：

  * `GW1.gronwall_floor` —— 单侧桥（D04-ODE 命题 3.1 的单侧形式；
    证明 = 命题 2.1 (i)⇒(ii) 的首零点论证，未用 mathlib 的 ODE/Grönwall API）：
    `f ∈ C¹([b₀,b₁])`、`f b₀ ≥ ε > 0`、微分不等式在 `f ≥ 0` 的点上成立
    （`f' x ≥ −C x · f x`）、`C` 连续且有整体原函数 `G`（`G' = C`；对连续 `C`
    这正是 "∫C 有限" 的可行形式：`G b₁ − G b₀ = ∫_{b₀}^{b₁} C`）
    ⇒ 全区间 `f x ≥ ε · exp(−(G x − G b₀)) > 0`。
  * `GW1.pos_C_exists` —— 等价性的另一半（D04 命题 2.1 (ii)⇒(i)）：
    `f ∈ C¹` 且在紧区间上恒正 ⇒ ∃常数 `C`：`f' ≥ −C·f`。
    两条合并 = "可证 C 的存在 ⟺ f 无零点"，即 D04 的"假设逐字等价于结论"。
  * `GW1.additive_floor` —— 加性版（Hellmann–Feynman 货币，(F-m) 载体实际能给的形状）：
    `f' ≥ −A`（`A` 连续、与 `f` 无关）⇒ `f x ≥ f b₀ − (G x − G b₀)`。
    取 `A ≡ 12V`（D04 (1.4)）即正文 §1 的 `seed − 29.8·L⁴` 之死。

  环境（第 216 轮实测）：lean 4.33.0-rc1 + mathlib badcache
  （`vendor/OSforGFF/.lake/packages/mathlib/.lake/build.badcache/lib/lean`，
  `LEAN_PATH` 依 `research/lean-infra/check.sh` 的方案拼
  `batteries;aesop;Qq;proofwidgets;LeanSearchClient;importGraph;plausible`）。

  对应正文：`research/v3/species/GW1-FLOW-TARGET.md`（§2 桥引理；§1 差分诊断）。
  无 QED；Clay 0%。
-/

import Mathlib

namespace GW1

open Set Filter MeasureTheory intervalIntegral

/-! ## 0 · 助手：连续性把单侧缺口闭合 -/

/-- 若 `g` 在 `[a,b]`（`a < b`）上连续，且 `c ≤ g t` 对一切 `t ∈ [a,b)` 成立，
则 `c ≤ g b`。首零点论证里"端点由左侧单调性接管"的一步。 -/
theorem left_le_of_leOn_Ico {g : ℝ → ℝ} {a b c : ℝ} (hab : a < b)
    (hg : ContinuousOn g (Icc a b)) (h : ∀ t ∈ Ico a b, c ≤ g t) : c ≤ g b := by
  by_contra hcon
  have hgb : g b < c := by linarith
  have hCW : ContinuousWithinAt g (Icc a b) b :=
    hg.continuousWithinAt (mem_Icc.mpr ⟨hab.le, le_refl b⟩)
  have hpre : g ⁻¹' (Iio c) ∈ nhdsWithin b (Icc a b) :=
    hCW.preimage_mem_nhdsWithin (isOpen_Iio.mem_nhds hgb)
  rw [mem_nhdsWithin] at hpre
  obtain ⟨U, hUopen, hbU, hsub⟩ := hpre
  obtain ⟨δ, hδ, hball⟩ := Metric.isOpen_iff.mp hUopen b hbU
  set δ' := min (δ / 2) ((b - a) / 2) with hδ'def
  have hδ'pos : 0 < δ' := lt_min (by linarith) (by positivity)
  have hδ'lt : δ' < δ := by
    rw [hδ'def]
    exact lt_of_le_of_lt (min_le_left _ _) (half_lt_self hδ)
  set z := b - δ' with hzdef
  have hz_lt : z < b := by linarith
  have hz_ge : a ≤ z := by
    have h1 : δ' ≤ (b - a) / 2 := min_le_right _ _
    linarith
  have hzinU : z ∈ U := by
    have hdist : dist z b = δ' := by
      rw [Real.dist_eq, hzdef]
      have hzz : b - δ' - b = -δ' := by ring
      rw [hzz, abs_neg, abs_of_pos hδ'pos]
    exact hball ((Metric.mem_ball).mpr (by rw [hdist]; exact hδ'lt))
  have hzgt : g z < c := hsub (mem_inter hzinU (mem_Icc.mpr ⟨hz_ge, le_of_lt hz_lt⟩))
  have hleft := h z (mem_Ico.mpr ⟨hz_ge, hz_lt⟩)
  linarith

/-! ## 1 · 主定理：单侧 Grönwall 桥（首零点论证，可积 C 版） -/

/-- **单侧 Grönwall 桥**。假设：`f` 在 `[b₀,b₁]` 上可导且导数连续（C¹）、
`f b₀ ≥ ε > 0`（种子）、微分不等式在 `f ≥ 0` 的点上成立（`f' x ≥ −C x · f x`）、
`C` 连续、`G` 是 `C` 在区间上的原函数。结论：对一切 `x ∈ [b₀,b₁]`，
`f x ≥ ε · exp(−(G x − G b₀)) > 0`。

对连续 `C`，取 `G x := ∫_{b₀}^{x} C` 即得任务书要求的 "C 可积" 形式
（此时 `G b₁ − G b₀ = ∫ C`；该 FTC 步未形式化，见正文 §2 注）。
指数上的量无任何限制——`exp` 恒正；**正性由种子单独承担**，
可积性只决定地板的数量级。 -/
theorem gronwall_floor
    (f C G : ℝ → ℝ) (b₀ b₁ ε : ℝ) (hb : b₀ < b₁) (hε : 0 < ε)
    (hfd : ∀ x ∈ Icc b₀ b₁, DifferentiableAt ℝ f x)
    (hC1 : ContinuousOn (deriv f) (Icc b₀ b₁))
    (hCc : ContinuousOn C (Icc b₀ b₁))
    (hG : ∀ x ∈ Icc b₀ b₁, HasDerivAt G (C x) x)
    (hseed : ε ≤ f b₀)
    (hineq : ∀ x ∈ Icc b₀ b₁, f x ≥ 0 → deriv f x ≥ -(C x) * f x) :
    ∀ x ∈ Icc b₀ b₁, ε * Real.exp (G b₀ - G x) ≤ f x := by
  classical
  have heG : ContinuousOn (fun t => Real.exp (G t - G b₀)) (Icc b₀ b₁) := by
    intro x hx
    have h1 : ContinuousAt (fun t => G t - G b₀) x :=
      (hG x hx).continuousAt.sub continuousAt_const
    exact (Real.continuous_exp.continuousAt.comp h1).continuousWithinAt
  have hfe : ContinuousOn
      (fun t => f t * Real.exp (G t - G b₀)) (Icc b₀ b₁) := by
    refine ContinuousOn.mul (fun x hx => (hfd x hx).continuousAt.continuousWithinAt) heG
  have hDcont : ContinuousOn
      (fun t => deriv f t * Real.exp (G t - G b₀) +
        f t * (Real.exp (G t - G b₀) * C t)) (Icc b₀ b₁) := by
    refine ContinuousOn.add (ContinuousOn.mul hC1 heG)
      (ContinuousOn.mul (fun x hx => (hfd x hx).continuousAt.continuousWithinAt) (heG.mul hCc))
  have hPhiDeriv : ∀ t ∈ Icc b₀ b₁,
      HasDerivAt (fun u => f u * Real.exp (G u - G b₀))
        (deriv f t * Real.exp (G t - G b₀) +
          f t * (Real.exp (G t - G b₀) * C t)) t := by
    intro t ht
    have h1 : HasDerivAt f (deriv f t) t := (hfd t ht).hasDerivAt
    have h2 : HasDerivAt (fun u => G u - G b₀) (C t) t :=
      HasDerivAt.sub_const (G b₀) (hG t ht)
    have h3 : HasDerivAt (fun u => Real.exp (G u - G b₀))
        (Real.exp (G t - G b₀) * C t) t := h2.exp
    exact h1.mul h3
  -- 步一：在 "f > 0 于 [b₀,y)" 的伸长上，Φ = f·e^{G−G(b₀)} 单调 ⇒ 结论在 y 处成立
  have step : ∀ y ∈ Icc b₀ b₁, (∀ t ∈ Ico b₀ y, f t > 0) →
      ε * Real.exp (G b₀ - G y) ≤ f y := by
    intro y hy hpos
    rcases lt_or_eq_of_le hy.1 with hby | hby
    · -- 主情形 b₀ < y：先证 ∀ y' ∈ [b₀, y): Φ y' ≥ ε，再用助手闭合端点
      have hPhi : ∀ y' ∈ Ico b₀ y, ε ≤ f y' * Real.exp (G y' - G b₀) := by
        intro y' hy'
        have hy'Icc : y' ∈ Icc b₀ b₁ :=
          ⟨hy'.1, le_trans (le_of_lt hy'.2) hy.2⟩
        have hy'le : y' ≤ b₁ := hy'Icc.2
        have hy'le₀ : b₀ ≤ y' := hy'.1
        have hposy' : f y' > 0 := hpos y' hy'
        have hDpos : ∀ t ∈ Icc b₀ y', 0 ≤ deriv f t * Real.exp (G t - G b₀) +
            f t * (Real.exp (G t - G b₀) * C t) := by
          intro t ht
          have htIcc : t ∈ Icc b₀ b₁ :=
            ⟨ht.1, le_trans ht.2 hy'le⟩
          have htfpos : f t > 0 := by
            rcases lt_or_eq_of_le ht.2 with hLt | hEq
            · exact hpos t ⟨ht.1, lt_trans hLt hy'.2⟩
            · rw [hEq]; exact hposy'
          have hsign : 0 ≤ deriv f t + C t * f t := by
            have h' := hineq t htIcc (le_of_lt htfpos)
            linarith
          have hEq : deriv f t * Real.exp (G t - G b₀) +
              f t * (Real.exp (G t - G b₀) * C t) =
              Real.exp (G t - G b₀) * (deriv f t + C t * f t) := by ring
          rw [hEq]
          exact mul_nonneg (Real.exp_nonneg _) hsign
        have hDmono : ContinuousOn
            (fun t => deriv f t * Real.exp (G t - G b₀) +
              f t * (Real.exp (G t - G b₀) * C t)) (uIcc b₀ y') := by
          rw [uIcc_of_le hy'le₀]
          exact hDcont.mono (Icc_subset_Icc_right hy'le)
        have hint : IntervalIntegrable
            (fun t => deriv f t * Real.exp (G t - G b₀) +
              f t * (Real.exp (G t - G b₀) * C t)) volume b₀ y' :=
          hDmono.intervalIntegrable
        have hPhiD : ∀ t ∈ uIcc b₀ y',
            HasDerivAt (fun u => f u * Real.exp (G u - G b₀))
              (deriv f t * Real.exp (G t - G b₀) +
                f t * (Real.exp (G t - G b₀) * C t)) t := by
          intro t ht
          rw [uIcc_of_le hy'le₀] at ht
          exact hPhiDeriv t (Icc_subset_Icc_right hy'le ht)
        have hFTC : ∫ t in b₀..y',
            deriv f t * Real.exp (G t - G b₀) +
              f t * (Real.exp (G t - G b₀) * C t) =
            f y' * Real.exp (G y' - G b₀) -
              f b₀ * Real.exp (G b₀ - G b₀) :=
          integral_eq_sub_of_hasDerivAt hPhiD hint
        have h2 : 0 ≤ f y' * Real.exp (G y' - G b₀) -
            f b₀ * Real.exp (G b₀ - G b₀) := by
          rw [← hFTC]
          exact intervalIntegral.integral_nonneg hy'le₀ hDpos
        have h3 : f b₀ * Real.exp (G b₀ - G b₀) ≤
            f y' * Real.exp (G y' - G b₀) := by linarith
        simp only [sub_self, Real.exp_zero, mul_one] at h3
        exact le_trans hseed h3
      -- 用助手把 y' ↑ y 的下界搬到 y
      have hclose : ε ≤ f y * Real.exp (G y - G b₀) :=
        left_le_of_leOn_Ico (g := fun t => f t * Real.exp (G t - G b₀))
          (a := b₀) (b := y) (c := ε) hby
          (hfe.mono (Icc_subset_Icc_right hy.2)) hPhi
      have hepos : 0 < Real.exp (G y - G b₀) := Real.exp_pos _
      calc ε * Real.exp (G b₀ - G y)
          ≤ f y * Real.exp (G y - G b₀) * Real.exp (G b₀ - G y) :=
            mul_le_mul_of_nonneg_right hclose (Real.exp_nonneg _)
        _ = f y := by
            rw [mul_assoc, ← Real.exp_add,
              show (G y - G b₀) + (G b₀ - G y) = 0 from by ring,
              Real.exp_zero, mul_one]
    · -- y = b₀：结论即种子
      rw [hby, sub_self, Real.exp_zero, mul_one, ← hby]
      exact hseed
  -- 步二：首零点论证 ⇒ f 在全区间为正
  have hmain : ∀ x ∈ Icc b₀ b₁, f x > 0 := by
    intro x hx
    by_contra hcon
    have hfx0 : f x ≤ 0 := by linarith
    -- 零集 Z 非空、有下界；s := sInf Z
    have hZne : ∃ z, z ∈ {t | t ∈ Icc b₀ b₁ ∧ f t ≤ 0} := ⟨x, hx, hfx0⟩
    obtain ⟨x₁, hx₁Icc, hx₁le⟩ := hZne
    have hZbdd : BddBelow {t | t ∈ Icc b₀ b₁ ∧ f t ≤ 0} :=
      ⟨b₀, fun z hz => hz.1.1⟩
    set s := sInf {t | t ∈ Icc b₀ b₁ ∧ f t ≤ 0} with hsdef
    have hx₁mem : x₁ ∈ {t | t ∈ Icc b₀ b₁ ∧ f t ≤ 0} :=
      Set.mem_setOf.mpr ⟨hx₁Icc, hx₁le⟩
    have hs_ge : b₀ ≤ s :=
      le_csInf ⟨x₁, hx₁mem⟩ (fun z hz => hz.1.1)
    have hs_le : s ≤ b₁ := le_trans (csInf_le hZbdd hx₁mem) hx₁Icc.2
    have hs_mem : s ∈ Icc b₀ b₁ := ⟨hs_ge, hs_le⟩
    -- (i) s 左侧全正
    have hposleft : ∀ t ∈ Ico b₀ s, f t > 0 := by
      intro t ht
      have htIcc : t ∈ Icc b₀ b₁ :=
        ⟨ht.1, le_of_lt (lt_of_lt_of_le ht.2 hs_le)⟩
      by_contra hcon'
      have ht2 : t < s := ht.2
      have hmemZ : t ∈ {u | u ∈ Icc b₀ b₁ ∧ f u ≤ 0} :=
        Set.mem_setOf.mpr ⟨htIcc, by linarith⟩
      have h' : s ≤ t := csInf_le hZbdd hmemZ
      linarith
    -- (ii) step 在 s 处 ⇒ f s > 0
    have hsstep := step s hs_mem hposleft
    have hexp : 0 < Real.exp (G b₀ - G s) := Real.exp_pos _
    have hfspos : 0 < f s := by
      have : 0 < ε * Real.exp (G b₀ - G s) := mul_pos hε hexp
      linarith
    -- (iii) 连续性：s 附近（在 Icc 内）f > 0
    obtain ⟨U, hUopen, hsU, hsubU⟩ : ∃ U : Set ℝ, IsOpen U ∧ s ∈ U ∧
        ∀ t ∈ U ∩ Icc b₀ b₁, f t > 0 := by
      have hCW : ContinuousWithinAt f (Icc b₀ b₁) s :=
        (hfd s hs_mem).continuousAt.continuousWithinAt
      have hmem : f ⁻¹' (Ioi 0) ∈ nhdsWithin s (Icc b₀ b₁) :=
        hCW.preimage_mem_nhdsWithin (isOpen_Ioi.mem_nhds hfspos)
      rw [mem_nhdsWithin] at hmem
      obtain ⟨U, hUopen, hsU, hsubU⟩ := hmem
      exact ⟨U, hUopen, hsU, fun t ht => hsubU ht⟩
    -- (iv) 下界稠密性 + s 的正性邻域 ⇒ 矛盾
    obtain ⟨δ, hδ, hball⟩ := Metric.isOpen_iff.mp hUopen s hsU
    have hA : ∀ d > 0, ∃ z ∈ {t | t ∈ Icc b₀ b₁ ∧ f t ≤ 0}, z < s + d := by
      intro d hd
      by_contra hcon2
      push_neg at hcon2
      have h1 : s + d ≤ s := le_csInf ⟨x₁, hx₁mem⟩ hcon2
      linarith
    obtain ⟨z, hzZ, hzlt⟩ := hA δ hδ
    have hzge : s ≤ z := csInf_le hZbdd hzZ
    have hdist : dist z s = z - s := by
      rw [Real.dist_eq, abs_of_nonneg (by linarith)]
    have hzinU : z ∈ U :=
      hball ((Metric.mem_ball).mpr (by rw [hdist]; linarith))
    have hzinIcc : z ∈ Icc b₀ b₁ := hzZ.1
    have hzle : f z ≤ 0 := hzZ.2
    have hzpos : 0 < f z := @hsubU z (mem_inter hzinU hzinIcc)
    linarith
  -- 步三：收尾
  intro x hx
  refine step x hx (fun t ht => hmain t ⟨ht.1, le_of_lt (lt_of_lt_of_le ht.2 hx.2)⟩)

/-! ## 2 · 等价性的另一半：恒正 ⇒ 常数 C 存在（D04 命题 2.1 (ii)⇒(i)） -/

/-- `f ∈ C¹` 且在紧区间上恒正 ⇒ 存在**常数** `C` 使 `f' ≥ −C·f` 在全区间成立。
与 `gronwall_floor` 合并：可行 `C` 的存在 ⟺ `f` 在区间上无零点
（D04-ODE 命题 2.1 的"假设逐字等价于结论"）。 -/
theorem pos_C_exists
    (f : ℝ → ℝ) (b₀ b₁ : ℝ) (hb : b₀ < b₁)
    (hfd : ∀ x ∈ Icc b₀ b₁, DifferentiableAt ℝ f x)
    (hC1 : ContinuousOn (deriv f) (Icc b₀ b₁))
    (hpos : ∀ x ∈ Icc b₀ b₁, f x > 0) :
    ∃ C : ℝ → ℝ, ContinuousOn C (Icc b₀ b₁) ∧
      ∀ x ∈ Icc b₀ b₁, deriv f x ≥ -(C x) * f x := by
  classical
  have hgcont : ContinuousOn (fun x => deriv f x / f x) (Icc b₀ b₁) :=
    ContinuousOn.div hC1
      (fun x hx => (hfd x hx).continuousAt.continuousWithinAt)
      (fun x hx => ne_of_gt (hpos x hx))
  obtain ⟨M, hM⟩ := isCompact_Icc.exists_bound_of_continuousOn hgcont
  refine ⟨fun _ => M, continuousOn_const, fun x hx => ?_⟩
  have hfx0 : f x ≠ 0 := ne_of_gt (hpos x hx)
  have hid : deriv f x = deriv f x / f x * f x := by field_simp
  have h2 : |deriv f x| ≤ M * f x := by
    rw [hid, abs_mul, abs_of_pos (hpos x hx)]
    exact mul_le_mul_of_nonneg_right (hM x hx) (le_of_lt (hpos x hx))
  have h3 : -(M * f x) ≤ deriv f x := (abs_le.mp h2).1
  show deriv f x ≥ -(M) * f x
  rw [neg_mul]
  exact h3

/-! ## 3 · 加性版：Hellmann–Feynman 货币（与 f 无关的 A） -/

/-- **加性地板**（(F-m) 载体的实际形状）：`f ∈ C¹`、`f' ≥ −A` 逐点成立
（`A` 连续、**与 f 无关**——这正是 Hellmann–Feynman 的分子形状：
`dm/dβ = ⟨1|S_W|1⟩ − ⟨0|S_W|0⟩`，无 `1/Δ` 也无 `1/f`）、
`G` 是 `A` 的原函数 ⇒ `f x ≥ f b₀ − (G x − G b₀)`。
取 `A ≡ 12V` 即得 D04 §4.3 的 `seed − 29.8·L⁴` 之死。 -/
theorem additive_floor
    (f A G : ℝ → ℝ) (b₀ b₁ : ℝ) (hb : b₀ < b₁)
    (hfd : ∀ x ∈ Icc b₀ b₁, DifferentiableAt ℝ f x)
    (hC1 : ContinuousOn (deriv f) (Icc b₀ b₁))
    (hAc : ContinuousOn A (Icc b₀ b₁))
    (hG : ∀ x ∈ Icc b₀ b₁, HasDerivAt G (A x) x)
    (hineq : ∀ x ∈ Icc b₀ b₁, deriv f x ≥ -(A x)) :
    ∀ x ∈ Icc b₀ b₁, f x ≥ f b₀ - (G x - G b₀) := by
  intro x hx
  have hxle : b₀ ≤ x := hx.1
  have hΨ : ∀ t ∈ uIcc b₀ x, HasDerivAt (fun u => f u + G u) (deriv f t + A t) t := by
    intro t ht
    rw [uIcc_of_le hxle] at ht
    exact (hfd t (Icc_subset_Icc_right hx.2 ht)).hasDerivAt.add
      (hG t (Icc_subset_Icc_right hx.2 ht))
  have hcont : ContinuousOn (fun t => deriv f t + A t) (uIcc b₀ x) := by
    rw [uIcc_of_le hxle]
    exact ContinuousOn.add (hC1.mono (Icc_subset_Icc_right hx.2))
      (hAc.mono (Icc_subset_Icc_right hx.2))
  have hint : IntervalIntegrable (fun t => deriv f t + A t) volume b₀ x :=
    hcont.intervalIntegrable
  have hFTC : ∫ t in b₀..x, deriv f t + A t = f x + G x - (f b₀ + G b₀) :=
    integral_eq_sub_of_hasDerivAt hΨ hint
  have hintpos : 0 ≤ ∫ t in b₀..x, deriv f t + A t := by
    refine intervalIntegral.integral_nonneg hxle (fun t ht => ?_)
    have h' := hineq t (Icc_subset_Icc_right hx.2 ht)
    linarith
  rw [hFTC] at hintpos
  linarith

-- 自检：应只列 [propext, Classical.choice, Quot.sound]（无 sorryAx）
#print axioms GW1.left_le_of_leOn_Ico
#print axioms GW1.gronwall_floor
#print axioms GW1.pos_C_exists
#print axioms GW1.additive_floor

end GW1
