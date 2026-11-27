# GGG5-RECON — OS 重构给的是哪种"存在"，`Δ>0` 是不是它的推论

> 2026-11-27 · 靶：`OS0–OS3/OS4 ⇒ ℋ_phys, H_op ≥ 0, Ω, φ` 这条桥，对 Clay 供给多少。
> **净判（三句）**：**(1)** 是 **(a) Wightman 场**，不是 (b) 只给矩、也不是 (c) 只给代数 —— 代数是重构的中间物，输出比它多；但带三条限定：**所选相 / 子列 / 规范不变场**。**(2)** Clay 的"存在" = `AX ∧ NT ∧ GAP`（`research/unconditional/B2-bypass-gap.md:6`），重构只兑掉 `AX`；差 `NT`、差"理论唯一"、差一句场内容陈述。**(3)** `Δ>0` 是 `H_op` 的谱性质，**不是**重构的推论：无质量自由标量是 `OS0–4 ∧ 唯一真空 ∧ 聚类 ∧ Δ=0` 的精确反例。**零新 MC，无 commit。**

## §1 重构的产出（逐条）

| 产出 | 内容 | 出处 |
|---|---|---|
| Hilbert 空间 | 反射正内积下商零模后完备化，`ℋ_phys` | `research/mass-spectrum/REPORT.md:48` |
| `H_op` | 自伴、`H ≥ 0`、`HΩ = 0` | 同上 `:48` |
| 真空唯一 | `ker H = ℂΩ` ⟸ 聚类/遍历（OS4 ≡ AN1-E4） | `research/uniqueness/REPORT.md:66`；`research/handoff-audit/16-input-audit.md:59` |
| Poincaré 表示 + 谱条件 + 局域性 + 场为算子值分布 | 需 OS0/OS1（解析性、增长）+ 解析延拓（管域、Laplace、`O(1,3)`） | `research/lean-infra/REPORT.md:486`；`docs/yang-mills-prior-work.md:64` |
| **场内容** | 只能是**规范不变**局域算子代数；规范固定极限不供给 | `research/unconditional/T3-construction.md:11` |
| 输入范围 | 只对**所选相**的关联给 E0–E4；子列极限逐点满足（含 E4） | `16-input-audit.md:31`；`research/uniqueness/REPORT.md:78` |
| **不在产出内** | 非平凡、理论唯一、`Δ>0` | `16-input-audit.md:40,59`；`B2-bypass-gap.md:10` |

`[T]` 结论：重构是**跨维一致**的翻译定理，对 4D 无特设附加公理；4D 缺的是"存在满足 OS0–OS4 的测度"这件**输入**，不是重构多要一条条件。
`[T]` 结论（子问题二）：OS4 簇聚 ⟹ 该理论内真空唯一（up to 相位），但这只到"相内"，**不**到"理论唯一"——`uniqueness` 判：OS4 是每个极限点自身的性质，不能排除另一子列落到另一 Gibbs 态（`research/uniqueness/REPORT.md:80–87`）。
`[NG]` 结论：**`I-RECON: OS0–OS4 ⇒ Wightman + Δ>0` 这一格不成立**（`research/gluing/REPORT.md:1348–1350,1620`）；正确形态是该报告自己的 `docs/conditional-theorem.md:301`："OS0–OS4 **+ 隙** ⇒ Wightman + 谱隙"。
》下一手：把 `gluing:1348–1350/1620` 的盒标签拆成"重构（无隙）"与"隙（额外输入）"两项，并由 `Lead` 同步 `docs/position-current.html`。

## §2 与 Clay 要求的差

1. **`NT`（非平凡）。** Gaussian 测度同样过 OS0–OS4，故重构不排除自由场；合格见证只能取 `g_ren ≠ 0` / `Λ > 0`（`research/nontriviality/REPORT.md:11–13,23`）。`[T]`：差一项独立输入。
2. **理论唯一。** Clay 要"a theory"，构造给"每相一个 + 子列"；补法是 Gibbs 态唯一（W1/Dobrushin，`uniqueness` 定理 A）。`[T]`：差一项独立输入。
3. **场内容。** `[T]` 重构输出为规范不变场代数；写成 Wightman 理论时须显式声明"场 = 规范不变局域算子"，否则第 (a) 项退化为 (c)。
4. **`Δ>0` 不在其中。** 见 §3.3。`[T]`
》下一手：在重构接口处补三条并列字段 `NT`、`UNI`、`fieldContent`，每条指向`16-input-audit`的具名输入编号，不塞进 `AX`。

## §3 判定

1. **精确产出 = (a)，带三限定**（所选相 / 子列 / 规范不变场）。`[T]`。下一手：把这三限定写成接口签名的三个参数，使"存在"一词在 Lean 层无法被读宽。
2. **Clay 的"存在" ≠ 重构的产出**：Clay = `AX ∧ NT ∧ GAP ∧ 构造`（`B2:6`），重构兑 `AX`；非平凡与唯一性是另外两项，`Δ>0` 是第三项。`[T]`。下一手：条件定理里把三者并列，禁止用"已由 OS 重构保证"合并陈述。
3. **`Δ>0` 是 `H_op` 的谱性质，但只能由额外的场论陈述供给**：Euclidean 侧"对 `μ` 一致速率的指数聚类"（`16-input-audit.md:59`；`docs/conditional-theorem.md:301` 已把它列在**前提**列），重构后经 Källén–Lehmann 转成 `spec(H_op) ∩ (0,Δ) = ∅`。`[T]` 反例：`B2-bypass-gap.md:10`；`research/unconditional/BASE-PHASE-RIGOROUS.md:25`（RP 不是隙证书）。`[NG]`："由 OS 公理自动得隙"这一读法。下一手：交接不等式只对"一致速率"这个量写，绝不写成"OS⟹隙"。

## §附 工具状态（AGENTS §1）

| 工具 | 状态 | 检索式 | 命中 |
|---|---|---|---|
| `pwsh` | OK（exit 0，`Get-Location` 回 `D:\Dev\math`） | `Test-Path research\v3\species`；`Get-ChildItem` | 目录存在；292 项 |
| `grep` | OK（250/940 截断，全表落盘） | `OS\s*重构\|Osterwalder\|Schrader\|Wightman\|reconstruction` @ `D:\Dev\math`，`*.md` | 940 |
| `grep` | OK | `重构.*规范不变\|规范不变.*重构\|gauge-invariant.*reconstruct` | 15 |
| `glob` | OK | `research/v3/species/GGG*` → 0；`GG*` → 3（`GG1-BIT/GG2-KILL/GG3-CLASS`） | 0 / 3 |
| `read` | OK | `H01-OSSCHRO.md`、`yang-mills-prior-work.md:28–137`、`uniqueness/REPORT.md:56–115`、`16-input-audit.md:1–80`、`nontriviality/REPORT.md:10–69`、`gluing/REPORT.md:1335–1384,1605–1649`、`conditional-theorem.md:292–321`、`B2-bypass-gap.md:1–20`、`lean-infra/REPORT.md:475–504`、`README.md:8–27` | — |

标签纪律：只用 `[T]`（重构定理/文件内明文）与 `[NG]`（`I-RECON` 含 `Δ>0` 那一格）。无 `TOOL-FAIL`。未改既有文件，未 commit。
