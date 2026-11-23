# RG 交接形式化 · lean_handoff

**文件**：`lean/YMStatement/RGHandoff.lean`（`import Mathlib`）。最终编译 **EXIT_CODE=0，无 `sorry`，无警告**。

## 编译命令（关键：题给的命令不可用）
在 `vendor/OSforGFF` 内 elan 选 `v4.33.0-rc1`（与 olean 匹配；全局是 4.36.0）。但 `mathlib/.lake/build/lib/lean` **不完整**（仅 1014 个 olean，缺整个 `Analysis`、无 `Mathlib.olean`），故 `lake env lean` 对 `import Mathlib` 报 `Mathlib.olean does not exist`。
完整构建在 `mathlib/.lake/build.badcache`（6.3 GB，8279 olean，含所需全部模块）。**未改 vendor 任何文件**，用 `LEAN_PATH` 前置该目录：

```
cd /d/Dev/math/vendor/OSforGFF && BAD='D:\Dev\math\vendor\OSforGFF\.lake\packages\mathlib\.lake\build.badcache\lib\lean' && export LEAN_PATH="$BAD;$(lake env printenv LEAN_PATH)" && lean /d/Dev/math/lean/YMStatement/RGHandoff.lean
```

冒烟测试（`import Mathlib`+`#check Real.exp_pos`）用此命令 37 s 通过；整个文件约 20 s。

## 通过的定理（0 sorry）
1. **一圈流**：`gSq_at_t0`、`gSq_deriv`——`d/dt g² = -2b(g²)²`，即 `g'=-bg³` 的平方形式。
2. **不变量**：`lam_const`（两时刻 Λ 相等）、`lam_deriv`（`Λ'=0`）。
3. **方案协变**：`exp_factor_identity`、`lamShift_ratio`（精确相对因子 `e^{c/(2b(1+cg²))}`）、`lamShift_uv`（`g²→0` 退化为常数 `e^{c/2b}`）、`crossing_not_invariant`（`g²=1` 的穿越随 `c` 平移）。
4. **质量 q**：`qfun_eq_exp`、`qfun_t0`（`q(t₀)=1`）、`qfun_strictAnti`（严格递减 ⇒ IR 方向 `q` 增）、`qfun_continuous`、`qfun_tendsto_zero`（`q→0`）、`qfun_unique_crossing`（唯一穿越，显式 `t*=t₀`）。
5. **交接**：`crossing_exists_unique`（IVT：`μ*` 存在唯一）、`coverage`（重叠 `t_up≤t*` ⇒ 全尺度覆盖）、`gap_lower_bound`（`gap≥min(δ_W,δ_S)`）。

## sorry 清单
无。

## 与 `q=m/μ` 框架、`Basic.lean` 的对应
`Λ=e^t·e^{-1/(2b g²)}` 与 `q=m/μ` 同一框架：覆盖引理把每个尺度归入弱域 `{μ≥μ_up}`∪强域 `{μ≤μ*}`，重叠条件化为 **scheme-无关的纯数不等式 `q₀≥q₁`**。`gap_lower_bound` 给出均匀正下界 `m_*=min(δ_W,δ_S)`，是 `Basic.lean` 中 `ExponentialClustering` 速率 `m` 的候选。要真正得到 `MassGapStatement` 尚需：(i) `gap↔m` 对应；(ii) 两侧 `δ_W,δ_S` 的定量界（H1,H2）；(iii) `q₀≥q₁`（公开难点）。本文件仅为纯实分析骨架，**不声称已证 MG**。

## 备选方案
无需纯 core 版本——badcache 路线可用。若将来 badcache 失效，可退回纯 Lean core（抽象有序域类或 ℚ 上显式公式），但 `Real.exp/log/deriv` 需重建，成本高。
