# STATUS — gauge-invariance 子项目

> 2026-11-22 · 状态快照。详细内容见 [REPORT.md](REPORT.md)。

## 一句话

规范不变性不是 4D 质量间隙的解析障碍；它是可分离的结构/记账问题。框架 = 轨道空间测度 + 等变块自旋 + 不变控制泛函 `q`。真正的开放点（Balaban 定量类保持、`q₀≥q₁` 交接、4D counterterm）与规范无关。

## 已确定（含依据）

| 结论 | 依据 | 状态 |
|---|---|---|
| Wilson 测度 `𝒢_Λ`-不变，下降到 `𝒜/𝒢` | Haar + 迹不变 | 已证/初等 |
| 面平均块自旋 `Φ` 精确等变（`δ≡1`） | REPORT §2.6 的显式计算 | 本文，初等 |
| `Φ` 等变 ⟹ pushforward 不变 ⟹ 轨道空间 RG 良定义 | REPORT L2.2/L2.4 | 本文，初等 |
| OS 反射正性 + 转移矩阵在不变子空间 | Osterwalder–Seiler 1978 | 文献已证 |
| 强耦合体积一致聚类（MG 的原型） | OS 1978；SZZ 2023 | 文献已证（仅强耦合） |
| `q = m/μ` 规范不变 | 定义在不变代数上 | 本文 |
| `q` 对耦合 scheme 变换不变；对 blocking 有 `O((a/ξ)²)` 修正 | REPORT §5.2 | 本文精确化 |
| 有限格点规范固定独立性 | REPORT §3.5 | 本文，初等 |
| 桥：`q` 可在任意规范固定下计算 | 由上一行 | 本文推论 |
| `Basic.lean` 是标量陈述，非 YM | 文件诊断 | 本文 |

## 已提出、未验证

- Balaban 的"精确定量等变界"在文献中的确切形式（需读 CMP 98/99/109 原文）；`【待核实】`
- Bauerschmidt–Dagallier U(1) 4D LSI 的 arXiv 号与确切陈述；`【待核实】`
- Giles 1981 的精确范围（U(N)/SU(N)/一般紧群）；`【待核实】`

## 开放（与规范无关）

1. **O1** Balaban 型块自旋的定量类保持（4D SU(N)）——本框架下唯一"规范相关"的硬点，但难在正则性/类保持，不难在等变性。
2. **O2** `q₀ ≥ q₁` 弱→强交接（交接报告的范围；本报告只保证 `q` 的定义与不变性）。
3. **O3** 4D 无穷 counterterm / 非超可重整（CCHS 3D 的唯一 gauge-covariant 重整化无 4D 对应）。
4. **O9** 2026 构造的标度极限规范不变性是否被证明（审计项）。

## 下一步（按优先级）

| # | 工件 | 类型 | 验收 |
|---|---|---|---|
| N1 | 规范不变陈述层 `GaugeStatement.lean` | Lean | 编译；β=0 测度过 `os`/`gap` 但不过 `nonTrivial`；自由 Maxwell 过 `os` 不过 `gap` |
| N2 | 有限格点规范固定独立性 | Lean | 0 sorry；小例子数值验证 |
| N3 | 等变 RG 结构引理（L2.1–L2.4） | 文档 | 可作为主 session 的引用 |
| N4 | Balaban CMP 98/99/109 精确陈述审计 | 文献 | 提取正则类、误差界、缺陷界 |
| N5 | 与主 session 的 `q` 接口 | 交接 | 不重复解析估计 |

**建议先做 N2**：它是本子项目唯一完全自足、可编译、且直接回应"把规范固定构造接到不变可观测量"的 Lean 定理。

## 文件

- `REPORT.md` — 完整框架、引理链、障碍、下一步、文献核实表。
- `STATUS.md` — 本文件。
- 共享 hand-over：`C:\Users\hooya\AppData\Local\Claude-3p\side-session-notes\local_893d10be-22f3-4514-9901-4f61e62e9204\gauge-invariance-handoff.md`

## 与其他 session 的边界

- **不重复**主 session 的 RG 交接控制量分析（`research/physics_handoff.md`）。
- **不重复** OS2/OS3、唯一性、非平凡性、group-uniformity、mass-spectrum、confinement、lean-infra 各 session 的工作。
- 本子项目的输出是：规范不变的**陈述层**、**等变 RG 结构引理**、**桥定理**，以及 `q` 的不变性论证。
