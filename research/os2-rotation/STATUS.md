# STATUS — os2-rotation

> 2026-11-22 · OS2 / Euclidean invariance for the 4D Yang–Mills continuum measure.

## 当前状态

**调研完成，REPORT.md 已交付。** 无代码改动，未运行 `lake build`（主 session 在重建 Mathlib）。
本目录只包含本报告与状态文件。

## 已完成（可验证）

- 读 `docs/yang-mills-prior-work.md`、`lean/YMStatement/Basic.lean`、
  `vendor/OSforGFF/OSforGFF/OS/Axioms.lean`、`OS2_Invariance.lean`、`Spacetime/Euclidean.lean`。
  确认 OS2 的 Lean 陈述是全 `E(d) = ℝᵈ ⋊ O(d)`，且 GFF 的 OS2 是**识别式**证明
  （协方差显式计算 + 变量替换）。
- 读 OpenAI 2D O(3) 预印本源文件 `build/sections/uniqueness.tex`（经 GitHub raw），
  核对其 (3.7)/(3.8)、Lemma 3.8、Niven 与 `d/γ ∉ ℚ` 论证。
- 新证：**命题 D** —— 三个 (3,4,5) 倾斜分块 `B₁₂, B₁₃, B₂₄ = 5·R_ij(ϑ) ∈ M₄(ℤ)`
  生成的旋转群稠密于 `SO(4)`（李代数括号张成 `so(4)`；附录 A 有完整计算）。
  **这是纯密度引理，不是 O(4) 证据**（gluing R28）。
- 明确否证了一个说法：不存在关于"格点极限恢复 Lorentz 不变性"的严格
  "Weinberg 定理"；真实内容是 Symanzik 有效作用量 + 无关性（检索确认）。
- 与兄弟 session 交叉核对（REPORT §8.3）：与 gluing 的"OS2 = 输出 / GL5"
  一致；补齐顺序极限次序与 `L`-一致量词；发现并记录倾斜分块与标准时间反射
  不对易（`θ_fine` 不保 `B₁₂ℤ⁴`）。
- **gluing R28 已采纳（重要更正）**：`B_ij aℤ⁴ = R_ij(5aℤ⁴)`，所以倾斜正规化的
  极限**自动**是 `(R_ij)_*μ`；"倾斜极限 = μ"与结论等价，是循环的。因此
  **命题 D 降为纯密度引理**，路线 A 从"路线"降为"归约"；O(4) 只能由 B
  （无关算符）或 C（唯一性）得到。原 `(REG-IND)` 假设已删除，改为 `(Γ-INV)`
  （三个旋转的不变性），并在 §1.3 加备注说明它不能由倾斜格点推出。
  见 REPORT §0.2、§1.3 备注、§3.3、§3.4、§6、§7、§8.2–8.3。

## 路线 B（NP-SYM）判定 — 见 NPSYM.md

- **NP-SYM ⟹ O(4) 不成立。** 两重理由：(i) ρ-形变 + `H(4)` 只生成单项矩阵，与
  `O(4)` 的交是 `W`；(ii) 每个 `μ_{a,ρ}` 的点群是 `B₃×ℤ₂`（与 ρ 无关），等极限
  只继承有限群。另加第三重：NP-SYM 只覆盖对角/度规型的 dim-6 方向。
- **与 L7**：`L7 ⟹ NP-SYM`；`NP-SYM ⇏ O(4)`。所以 NP-SYM 是 L7 的弱推论，
  不是路线。（`L7 ⟹ O(4)` 需类在 `O(4)` 下封闭。）
- **B 的 NP-SYM 形式：refuted。** 已写入 REPORT §0.2 与 §8.3(7)；保留的活形式是
  `SYM-4D`。
- **正确目标**：`SYM-4D`（有效作用量的 `O(4)`-破坏部分一致 `O(a²)`）。它 ⟹ OS2，
  但**强于** OS2，且与 L7 不可比。另有 `SYM-IR`：`a²` 系数是全空间连通关联积分，
  体积一致需要聚类型输入——**但该纠缠只对 SYM-4D 成立**；D 的 U2 与 MG 逻辑独立。
- **路线表按 uniqueness/R32 的映射更正**：C = O1(a)（两轨迹差估计，`B ⊂ C`）、
  C′ = 唯一不动点（独立更强，2D 不用，⇏ B）、B = SYM-4D、D 独立。取舍：证 C 得到
  U1+U2 但最难；证 B 只解旋转方向。
- **§14 决定性计算：`REVERSAL-OPEN`。** ℓ=4 massive Gauss 的 Noether 应力张量不对称
  （`A_{μν}=λp_μp_ν(p_μ²−p_ν²)`），且 Belinfante 改进像为 `{0}`（两个允许结构的散度
  全为零；`check-d14.py` 精确整数检验）。故它不满足 `Conserved.symmetric`/`StressLaws`，
  **不是** `MG+QA ⇏ U2` 的反例。纪律：这不证明 `(Q-CONS)`，只说明反转未被此模型否定；
  自由理论反例已被排除，反例须来自"有对称守恒 T 却破坏 U2"的交互理论。
- **精度修正（NG 的强弱）**：`NG` 的强形式是 `Q_gΩ=0`（∀g ∈ so(4)），
  等价于 U2；`ℓ=4` 只是**领头检验**（必要不充分：`W`-不变通道还有 ℓ=6,8,…）。
  对 T–T 关联，spin-2×spin-2 ⟹ ℓ≤4，故 `(NG-ℓ4^T)` 可能更锐，但
  "`Q_gΩ≠0 ⇒ T–T 有 ℓ=4`" 未证。层级：`(NG-ℓ4) ⟸ (NG-ℓ4^T) ⟸? (NG) ⟺ U2`。
- **ROUTE-D.md §续 2（谱论证）**：谱定理（Q 自伴 + [H,Q]=0 + 唯一零能 ⟹ QΩ=0）正确；
  但 `[H,Q_g]=0` 只到"边界项"：局部 Ward 给 `[H,Q_g]=i·边界`。故 `MG ⇒ U2` 是
  **条件性反转**：`MG + QA + (Q-CONS) ⟹ U2`，缺口是边界算子 `(Q-CONS)`。ℓ=4 massive
  Gauss 是决定性检验（若满足 `StressLaws`，则它是 `MG ⇏ U2` 的反例）。§9.3 已修正：
  二次各向异性不 `W`-不变，须用 ℓ=4。
- **ROUTE-D.md §续**：`NoSpacetimeGoldstone` 写成 ℓ=4 命题（三等价形式：谐波投影/
  `Q_gΩ=0`/积分 Ward 边界项）；**可证伪**——各向异性 Gauss 满足平移+W+RP+StressLaws
  却违反它，故它不是这些公理的推论。`kills` 路线的承重项 (K2)（轨道导数恒等式）
  等价于无自发破缺，不是捷径。Lean 冗余性测试**被阻塞**（缺上游 olean，整仓构建
  超出范围，未做）。
- **路线 D 已加入并判定**（`ROUTE-D.md`）：RP + Ward + 应力衰减 ⇒ U2，缺口
  `NoSpacetimeGoldstone`。判定：(i) "无矢量序参量"不足，真正的序参量是 `ℓ=4`
  各向异性；(ii) 证书字段 `sym` **非冗余**（`Oriented` 有限协变 vs `flux` 无穷小
  Ward，各向异性 Gauss 可分离）；(iii) D **不绕开聚类**，消费 `StressTensorDecay`
  （边界项），但通道比 SYM-4D 窄。
- **§7.4 已按源码更正**：2D 的致命项是 (3.8) 的**精确性**（一维流/超可重整），
  不是秩/Niven；Lemma 3.8（`lem:dense-periods`）结论是"旋转不变 **且** 与 s 无关"，
  前提是精确形状关系 + 连续性，**不是**唯一不动点。
- **源级事实**（据 uniqueness session 核对）：`Selection.Certificate` 的
  `sym i : R.Oriented (C.ensemble i)` 是旋转不变性的 certificate 字段；U2 = 证它。
- **已交付 Lean 工件**：`lean/SYMDim4.lean` + `check-symdim4.sh`。三个定理
  （`invariant_eq_trace_div_four`、`invariant_sym2_eq_smul_one`、
  `average_eq_trace_div_four`）**已实际编译通过**，0 sorry，公理
  `[propext, Classical.choice, Quot.sound]`。用与 `lean/check-handoff.sh` 相同的
  pin 与 vendored olean cache，只跑单文件 `lean`，不跑 `lake`。
- **已形式化**：`SYM-DIM4`（dim-4 无各向异性），284 行，见上。注意其精确措辞是
  "系数张量 `= cδ`"，不是"dim-4 算子空间一维"（后者是 `{tr F², tr F̃F}`，两者均
  `O(4)`-不变）；"系数张量 ⟹ 无各向异性"那一步是物理，非 Lean。
- **审计**：4D 无非微扰 Symanzik 定理；"Bauerschmidt–Hairer 3D" 查无此文（应为
  CCHS, Invent. Math. 2024）。2D 的证明在秩/形状维数/精确分块/超可重整四处
  本质用 `d=2`。

## 结论（一句话）

OS2 不能由格点对称性 + 连续性得到（`W` 是有限闭群），也不能由倾斜正规化得到
（R28）；它归约为三个特定旋转的不变性，而这三个不变性只能由无关算符（L6）或
唯一性（L7）证明。形式骨架（L1–L4）可证；物理内核与公开难题同级。

## 未完成 / 未验证

- L6（Symanzik 一致 O(a²)）、L7（不动点/极限唯一）：**开放**。
  L5（"倾斜关系"）在字面表述下是恒真推前，**不构成输入**。
- 4D 构造（16 输入）是否含"形状参数 / 规范化唯一性"结构：**未审计**（无本地副本）。
- Zenodo 18825033（Kirk）与 Zenodo 18644280：**未审计，仅作线索，不入依赖链**。
- 所有"可证"项都未实际形式化（本 session 不开工 Lean）。

## 下一步（按价值排序）

1. **审计 4D 构造的 `RG1–RG5` / `Y1–Y4`**：有效作用量类是否在 `O(4)` 下封闭？
   是否跟踪 `O(4)`-破坏方向（L6 是否已隐含）？→ 决定路线 B 的具体形式。
   需要主 session 或用户提供仓库/论文访问。
2. **Lean 引理 L4** `os2_of_dense_subgroup`（~100–200 行，0 sorry，与 YM 无关）。
3. **Lean 引理 L4′**（OS0+OS1 ⟹ `Z_μ` 连续，Zorn）；GFF 情形先做模板。
4. **Lean 模板**：2D `⟨Rot(2ϑ)⟩` 稠密于 `SO(2)`（`AddCircle.denseRange`）。
5. **数值**：验证命题 D 的秩条件；并用 2D 高斯自由场自检 R28（倾斜极限 =
   旋转后的极限，应当自动成立）。
6. **陈述层**：OS2 的 YM 版本改写成 loop expectation 形式
   （当前 `MassGapStatement` 是标量场配置空间，见 REPORT §1.2 备注）。

## 风险 / 诚实说明

- 倾斜分块不能用来证明 O(4)（R28）；不要再把它当作路线。若 L6 与 L7 都不可得，
  OS2 就是开放的，没有第三条路。
- 不要误读命题 D：它是密度引理，不是不变性证据。
- 不要把"OS2 是技术细节"当作前提：它和唯一性同级，可能是比 MG 更难的形式化
  前置条件。但 OS2 不需要质量间隙。

## 文件

- `REPORT.md` — 主交付物（定理陈述、引理链、两条路线、可证/开放、下一步）。
- `NPSYM.md` — 路线 B 的判定：NP-SYM 陈述、Duhamel 化归、`NP-SYM ⇏ O(4)` 的证明、
  与 L7 的关系、正确的 `SYM-4D`、`SYM-DIM4` 引理、4D/2D 审计。
- `ROUTE-D.md` — 路线 D 的源级判定：`Oriented`/`StressLaws` 逐字引用、冗余性判定、
  `NoSpacetimeGoldstone` 的精确形式、与 B/C 对比。
- `check-d14.py` — §14 的精确整数检验（改进像为零、`A≠0`）。
- `lean/SYMDim4.lean` — 284 行，0 sorry；三个定理，标准公理。
- `check-symdim4.sh` — 一条命令编译 + 验公理。
- `STATUS.md` — 本文件。
