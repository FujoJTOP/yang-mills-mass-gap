# STATUS — glue-a (GL3 / GL4 / GL9)

> 2026-11-22 · 4D pure Yang–Mills · full report: `REPORT.md`
> 记号逐字引自 `research/gluing/REPORT.md` §1；本目录只增加引理层。
> **分类约定**：`✅` = 不消费开放输入；`✅-蕴含` = 消费**点名**输入并在陈述中写出。

## 0. 一句话

**GL4 = ✅（作为引理，无条件，全锥机器检验）**。组合层面经 GL3 仍依赖 `UVS`（弱极限存在）。**GL3 = ✅-蕴含**（消费 `UVS(L) ∧ UVS-L`）；**GL9 = ✅-蕴含**（消费 `GAP = CPB-W`）。

```
GL4（作为引理）：无条件 ✅（全锥已验证）
GL4（作为组合步骤）：消费 GL3 的 UVS（弱极限存在）
```
GL3 的有界类常数 = 1，是**空洞**的，不算成果；GL9 的核心谱不等式已 Lean 机器检验。

## 1. 状态表（可记账）

| 引理 | 分类 | 陈述（简） | 机器证据 | 点名消费 |
|---|---|---|---|---|
| **GL3** | **✅-蕴含** | `UVS(L) ∧ UVS-L ⟹ {S^{(n)}_{L,a}}` 相对紧、子列弱-* 收敛、界 a-一致；`UVS-L` 下阶段 2 子列紧 | `exists_subseq_tendsto_pi`（Tychonoff `isCompact_univ_pi` + `IsCompact.tendsto_subseq` + `tendsto_pi_nhds`）；`gl3_imp_with_minlos`（**`hMinlos` 显式**，`P` 正性由 GL4） | **`UVS(L)`, `UVS-L`**；Minlos 为显式外部假设 |
| **GL4** | **✅（全锥）** | 弱极限保持 RP 锥（有界连续 + L² 截断闭包）与指数聚类界（`(U_a)` 显式） | `rp_transfer`（真实弱收敛）、`rp_l2_closure_trunc`（L² 闭包）、`clustering_transfer`（`hU`）；反例脚本 | 无开放输入（消费 GL3 极限、OS'78）；唯一显式假设 `trunc A n ∈ L²` |
| **GL9** | **✅-蕴含** | 一致隙 ⇒ 体积一致指数聚类（速率 `δ/2`） | `spectral_decay`（核心谱不等式） | **`GAP` = `CPB-W`** |

按 `gluing/REPORT.md` §3.4 的"Lean 或可复现计算"标准：**四个 Lean 定理 + 一个可复现反例**，全部可一键复现。

## 2. 三条可供其它 session 直接使用的结论

1. **GL3 的紧性对有界可观测是空洞的**：归一化 Wilson 环 `|Φ|≤1` ⇒ 矩 `∈[-1,1]`，紧性是 Tychonoff 的免费推论，常数 = 1。GL3 的**全部**内容在 `UVS` 的 `(U1)-(U3)`（a-一致矩界、正确重正化减除、Minlos 等度连续）。**不得把有界类记作成果。**
2. **弱极限不保持聚类速率**（R28）：massive Gaussian `mₙ→0` ⇒ 弱极限 = massless GFF，RP 与弱收敛都成立但无指数聚类。**`MG` 必须"体积一致"**；`clustering_transfer` 把 `(U_a)` 写成显式假设。
3. **GL9 的常数与 L 无关是概率测度的免费推论**（有界局部算符 `‖O‖_{L²(μ)} ≤ ‖O‖_∞`）。Lieb–Robinson 只在 (i) 无界 smeared 算子、(ii) "Hamiltonian 隙 ⇒ 转移矩阵隙"处进入；一般欧氏方向速率是 **`δ/2`**，不是 `δ`（R29）。

## 3. 依赖（与 §3.0 供给矩阵）

- **GL3 消费**：A0/GL1、GL4-RP（正性）、**`UVS(L) ∧ UVS-L`**；**被消费**：I-OS2(A2)、I-OS3(极限)、I-UNIQ(TIGHT/子列)。
- **GL4 消费**：GL3（极限存在）、GL2a/OS'78（有限 a RP）、（GL4-C）a-一致聚类界；**被消费**：I-OS3(B3)、I-UNIQ（MG 传递）、I-OS2。
- **GL9 消费**：GL2a/OS3（`T≥0`）、OS'78（逐方向 RP）、**`GAP`**、可选 `LR-G`；**被消费**：I-UNIQ(W1_0，经 MG)、I-OS2、I-SPEC、GL3-L/GL4-C。
- **新增内部边（应补进 §3.0 下一版）**：`GL4-C → GL9`；`GL3-L → MG`。

## 4. 证伪测试（详见 REPORT §2.4/§3.5/§4.5）

| 测试 | 对象 | 做法 | 失败含义 |
|---|---|---|---|
| **F3-1** | `UVS` | SU(3) 沿 RT 测 plaquette 与重正化 2 点；看 `ξ/a` 是否增长 | 无 scaling window ⇒ a-一致界假 |
| **F3-2** | RT 归一化 | SU(3) `g²=1` 的 `t=7.178`（`24π²`） | 匹配 `48π²` ⇒ scheme 常数错，差 `e^{3.59}` |
| **F3-3** | GL3 证明结构 | 若证明不消费 `UVS`，应用到固定 β 序列 | 会"证明"不存在的连续极限 |
| **F3-5** | `UVS-L` | massless GFF on `T⁴_L`，`f=1_B/|B|`：`S^{(2)}_L(f,f)~L²` | 证伪"`UVS-L` 无条件成立" |
| **F4-1** | GL4-C | `gl4_gaussian_counterexample.py`（去掉 `hU`） | 极限无指数聚类 ⇒ `(U_a)` 必需（已复现） |
| **F4-2** | GL4-RP | 2D Ising `β_c`：RP、唯一、但幂律 | 只用 RP 的聚类证明错 |
| **F4-3** | GL4-RP 证明 | 检查"L² 截断 + 支撑闭包"步骤 | 只控有限测试族的证明错 |
| **F4-4** | GL4-RP | 反射面/分块不对齐（违反 §1.7 约束 1） | `(θ)_*μ=μ` 失效 |
| **F9-1** | GL9 | 2D Ising `T<T_c` 周期 bc：`Δ_L→0` | 只用"`Δ_L>0 ∀L`"的证明错 |
| **F9-2** | GL9 | `T⁴_L` massless GFF：`λ₂=e^{-2π/L}→1` | 同上 |
| **F9-3** | GL9 | SU(3) `β≈2.6`，`L=16,20,24` 的 `0^{++}`：`m_eff` 平台与 `C(L)` | 平台随 L 下降或 `C(L)` 增长 ⇒ 体积一致假 |
| **F9-4** | GL9 扇区 | 用 Wilson 环检验中性隙不控制弦 | 量词反例（用错隙） |

## 5. 下一步（按价值排序）

1. **给 GL9 补整条链**：`λ=e^{-aδ}`、`t=na`、四方向合并进 `GL9Decay.lean`；~200 行，机械。
2. **Minlos 仍是外部假设**（Mathlib 无）。若将来有 Mathlib 的 Minlos/Prokhorov 核形式，直接替换 `hMinlos` 即可，其余模块不动。

## 6. 团队层面（如实计数）

- 按本轮约定：**GL4 = ✅（全锥机器检验）**；**GL3 = ✅-蕴含**（`UVS ∧ UVS-L`）；**GL9 = ✅-蕴含**（`GAP`）。
- 按 §3.4 的"可复现证据"标准：四条 Lean 定理 + 一个反例脚本均可一键复现；但 GL3/GL9 的**结论本身**依赖点名输入，不得写成无条件 ✅。
- 加上 glue-b 的 GL2a（os3 定理，待机器化）、glue-c 的 GL1（合规表）与 GL6（有限格点初等），团队 ≥3/9 可达；仍需各 session 交可复现证据。
- **最大的诚实提醒**：九条里没有一条会自己产生质量隙。GL3/GL4/GL9 合起来只把"若有一致隙，则 MG"变成定理；`CPB-W` 仍是唯一的物理墙。

## 7. 复现

```bash
bash research/glue-a-limits/check-glue-a.sh
```

编译 `GL9Decay.lean`、`GL4Transfer.lean`、`GL4L2.lean`、`GL3Limits.lean`（公理审计：`[propext, Classical.choice, Quot.sound]`，无 `sorryAx`；后两者用 `#guard_msgs` 冻结），并运行 `gl4_gaussian_counterexample.py`。末行应为 `ALL GLUE-A CHECKS PASS`。
