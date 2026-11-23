# IN-2 — O(4) 不变性：状态、障碍、推进

> 源：`conditional-theorem.md` §3 IN-2/诚实清单 9；`os2-rotation/{REPORT,ROUTE-D,NPSYM}.md`；`glue-b-blocking/REPORT.md` §4；`lean/SYMDim4.lean`。

## 1 精确陈述
**IN-2**：`μ=lim_{L→∞}lim_{a→0}μ_{L,a}`（先 `a` 后 `L`）在 `O(4)` 下不变。等价于 (i) ℝ⁴ 上 Wilson 环/Schwinger 期望在 `E(4)=ℝ⁴⋊O(4)` 下协变；(ii) `∀g∈so(4)`，角荷湮灭真空 `Q_gΩ=0`，`Q_g=∫d³x(x_iT_{4j}−x_jT_{4i})`；(iii) 连通关联无 `ℓ=4` 各向异性。限定：有限环面 `𝕋_L` 只有 `W=(ℤ/2)⁴⋊S₄`（`|W|=384`），`O(4)/W` 仅在 `L→∞` 出现（REPORT:13）。

## 2 真实状态：开放，非 RP 推论
- 2026 构造将其作假设 `Selection.Certificate.sym`（`Mixture.lean:136–148`）；`SourceScalar.lean:34` 只证"字段⇒E1"。IN-2 即抽出该字段。
- **RP 与之正交**：有质量 `ℓ=4` 各向异性高斯 `Ĉ=1/(p²+m²+εΛ⁻²Σp_i⁴)` 满足 OS0/1/3/4、平移、`W`、Ward，却破 `U2`（`ROUTE-D §14`，`check-d14.py`）。故 RP+构造 ⇏ O(4)。
- 倾斜不能给：`B_ij aℤ⁴=R_ij(5aℤ⁴)` ⇒ 极限自动 `(R_ij)_*μ`，循环（R28）。

## 3 推进：KIN-2（已证；张量部分 Lean）
**命题 KIN-2**（`glue-b §4.2`，`SYMDim4.lean`）：各向同性点无维 `≤4` 的 `O(4)`-破坏规范不变算符。`dim<4` 只常数；`dim=4` 为 `{tr F²,tr F F̃}`，均 `O(4)`-不变。唯一候选 `Σ c_{μν}tr F_{μν}²`：`c` 是 `W`-不变对称张量 ⇒ `c∝δ`（`invariant_sym2_eq_smul_one`，0 sorry），本身 `O(4)`-不变。
**刀刃**：点群若仅 `(ℤ/2)⁴⋊S₃`（时间被区分），两族 `{(01,02,03)}`、`{(12,13,23)}` 不互换，`dim-4` 各向异性**存在**（各向异性晶格 YM，`β₀≠βᵢ`，极限有优先标架）；`S₄` 时空换位杀死它——**承重，故 OS2 非空洞**。
**推论**：**无旋转破坏的"相关"(relevant)算符**；`dim-4` 为 **marginal**（仅各向异性正则化、须调参）；残留从 `dim-6` 起，**无关、`O(a²)`**。题干"相关算符"前提为假——related-op 论证在此**满足**（GL5-KIN ✅），这正是预期恢复之由。

## 4 精确障碍（KIN-2 后仅两处）
1. **`a→0`（固定 `L`）**：`SYM-4D+SYM-IR`——把"`dim-6` 无关⇒`a²`"升为非微扰、体积一致界（`a²` 系数 = 全空间连通积分的体积一致控制）。Symanzik 系弱耦合；Balaban RG 停在 `λ~1`；4D 无非微扰 Symanzik 定理（`NPSYM §7.3`）。
2. **`L→∞`**：`Q-CONS`（=`NoSpacetimeGoldstone`）。局部 Ward 只给 `[H,Q_g]=iB_g`，`B_g=lim_R∮x_{[i}T_{j]k}n_k`；要 `Q_gΩ=0` 需边界角流湮灭。

## 5 翻转条件
- **(F-flip，最省)** 证 `Q-CONS`：`[H,Q_g]=0`+`kerH=ℂΩ`（三行谱定理）⇒`Q_gΩ=0`；配聚类(MG) 得 `IN-2⟸IN-3/4`，独立输入 6→5。`ℓ=4` 高斯 `REVERSAL-OPEN`（无对称守恒 `T`，非反例），可信未证。
- **(F-Sym)** 或证非微扰 `SYM-4D` + `L→∞` 无自发破坏。

二者均**非 RP 所给**。
