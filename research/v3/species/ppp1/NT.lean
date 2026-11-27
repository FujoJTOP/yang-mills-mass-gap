/-
# PPP1 / NT.lean — `NT` 的陈述层 + 引理 C + 「`κ₃ ≠ 0 ⇒ ¬Gaussian`」的已证一半

本轮（PPP1）判词：`NT` 保持 `[NG]`（`research/v3/species/PPP1-NT.md`）。本文件只固定
**陈述层**，并把直路唯一真能用的那条推理——引理 C 的 `n = 3` 逆否——机检化。
本文件 **`sorry` 为零**，因此它证明的东西不依赖任何缺口。

## 依赖（逐行标注；本轮**未编译**，按任务书只写代码 + 标依赖）

* `import YangMillsInfra.Statement.NonTriviality`
  - 经 `YangMillsInfra.Basic` 提供 `EuclideanTheory`（`Basic.lean:52`，字段 `Z`）、
    `TestFunction d = SchwartzMap (SpaceTime d) ℝ`（`Basic.lean:31`）。
  - 提供 `IsGaussianTheory`（`Statement/NonTriviality.lean:51`）、
    `cumulantGF`（`:63`）、`nthCumulant T J n = iteratedDeriv n (cumulantGF T J) 0`（`:67`）。
  - 提供**引理 C 充分方向** `isGaussianTheory_cumulants_vanish`（`:78`，已证，无公理）、
    必要方向接口 `CumulantCharacterization`（`:112`）、
    条件等价 `isGaussianTheory_iff_cumulants_vanish`（`:120`）、
    `GaugeNonTrivial`（`:198`）与其 `gaugeNonTrivial_implies_not_gaussian`（`:202`）。
  - 包在 `research/lean-infra`（`lakefile.toml`: `yangMillsInfra`, mathlib `v4.33.0-rc1`）。
    核验命令（本轮**未执行**）：`cd research/lean-infra && lake build YangMillsInfra.Statement.NonTriviality`。
* measure 侧的逐字目标 `¬ YMStatement.IsGaussianMeasure μ` 在 `lean/YMStatement/Basic.lean:31,52`，
  依赖 `OSforGFF`（本机未构建）⇒ 本文件用测度无关的同形 `IsGaussianTheory` 陈述；
  两者在支撑测度的理论上相同（`lean-infra/REPORT.md:417`）。
-/
import YangMillsInfra.Statement.NonTriviality

noncomputable section

namespace PPP1

open YangMillsInfra

variable {d : ℕ}

/-- `NT`：与 `YMStatement.NonTrivial μ = ¬ YMStatement.IsGaussianMeasure μ`
（`lean/YMStatement/Basic.lean:52–53`）同形的 `EuclideanTheory` 侧陈述。

**它不是新的判据**：`Basic.lean:84–87` 的 `MassGapStatement.nonTrivial` 字段就是这个谓词，
即 `¬Gaussian` **恰好**是 Clay 的 `NT`（不缺更强的东西；Clay 的另一合取项是 `GAP`）。 -/
def NT (T : EuclideanTheory d) : Prop := ¬ IsGaussianTheory T

theorem nt_iff_not_gaussian (T : EuclideanTheory d) : NT T ↔ ¬ IsGaussianTheory T :=
  Iff.rfl

/-- **引理 C**（`research/nontriviality/PROOF.md:52–53` 逐字）：
`IsGaussianMeasure μ ↔ ∀ n ≥ 3, ∀ f₁…f_n, κₙ(f₁,…,fₙ) = 0`。

此处以 `CumulantCharacterization`（充分方向已证 + 必要方向为接口）为前提引用；
必要方向需要 OS0/OS1 型解析性与矩唯一性（`PROOF.md:59`），本文件不重证。 -/
theorem lemmaC {T : EuclideanTheory d} (h : CumulantCharacterization T) :
    IsGaussianTheory T ↔ ∀ (n : ℕ), 3 ≤ n → ∀ J : TestFunction d, nthCumulant T J n = 0 :=
  isGaussianTheory_iff_cumulants_vanish h

/-- **引理 C 的 `n = 3` 逆否（主定理，`sorry` 为零）**：
`κ₃(⟨φ,J⟩) ≠ 0 ⇒ ¬Gaussian`。

* 只用**充分方向**取逆否，故**不需要**解析性/唯一性（`NonTriviality.lean:78`）。
  ⇒ 直路所需的**推理**是免费的，缺口全在**前提**（见 `NT-C3.lean`）。
* 辖域：`J : TestFunction d` 是**线性**测试泛函，`κ₃` 是随机变量 `⟨φ,J⟩` 的累积量。
  **`q = ¼tr F²` 之类的复合算子不在辖域内**（`nontriviality/REPORT.md:55–57` 命题 1.2）。 -/
theorem nthCumulant_three_ne_zero_implies_NT {T : EuclideanTheory d}
    {J : TestFunction d} (h : nthCumulant T J 3 ≠ 0) : NT T :=
  fun hG => h (isGaussianTheory_cumulants_vanish hG 3 (by decide) J)

/-- 合法见证形（`NonTriviality.lean:202` 已证）：`g_ren ≠ 0 ⇒ NT`。

与本文件上面的定理并列，标出**两条路的分工**：线性累积量（推理免费，但无规范不变的载体）
vs `g_ren`（唯一同时满足「规范不变」与「自由值 = 0」的载体，`REPORT.md:55–57`）。 -/
theorem nt_of_gaugeNonTrivial {T : EuclideanTheory d} (h : GaugeNonTrivial T) : NT T :=
  gaugeNonTrivial_implies_not_gaussian h

end PPP1
