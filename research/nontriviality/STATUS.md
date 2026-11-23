# STATUS — 非平凡性（非自由场）

> 2026-11-22 · 一页速览。详细论证见同目录 `REPORT.md`。

## 一句话

4D 纯 YM 连续极限的"非平凡"= **非 Gaussian（非准自由）**；其实质内容是**连续极限来自 `Λ>0` 的非平凡渐近自由轨迹，而非 `g₀≡0` 的 Gaussian 平凡轨迹**。这是**审计/定义层**结论，不是独立 UV 定理；一旦 `Λ>0` 成立，"短程 `g_ren≠0` 传到极限"是 UV 侧可控的，**不需 MG**。

> **测度来源（GL1/E2 已补）**：全部结论针对共享格点测度 `μ_{L,a,β} ∝ exp((β/N_c)Σ_p Re tr U_p)∏_e dU_e`，`β=2N_c/g²`、`N_c:=dim f`、`γ:=β/N_c=2/g²` 是测度 `exp(·)` 里的系数（`SU(3),g²=1 → β=6, γ=2`；勿把 `γ` 记作 `β`）。逐字引用 `research/gluing/REPORT.md` §1.2、`research/glue-c-gauge-group/REPORT.md` (1.1)–(1.2)；显式见 `PROOF.md §1.1`。

## 已建立的（确定）

- **正确的抽象层**：`NonTrivial := ¬IsGaussianMeasure`（`lean/YMStatement/Basic.lean`）方向正确；自由场 = Gaussian 测度。
- **核心陷阱（诊断完成）**：规范不变局域算子都是 `F` 的二次及以上（`tr F²`…），而自由规范场的 `tr F²` 是 Wick 平方，**连通 4 点非零**。故"规范不变连通 4 点 ≠ 0"作非平凡见证是**空洞**的（自由理论也满足）。命题 1.1/1.2，REPORT §1.2。
- **唯一合格见证**：**重整化耦合 `g_ren`**（规范不变的顶角/传播子比；自由 ⟹ `≡0`；相互作用 ⟹ `≠0`）。`g_ren≠0` 是 **scheme-无关谓词**（`g'²=g²+cg⁴` 保 `g=0`）；对照 `g²=1` 的穿越点**随 scheme 平移**（见 `numerics/t2_scheme_shift.txt`）。REPORT §3.3、§6.3。
- **平凡性对比**：`φ⁴₄` 平凡（Aizenman–Duminil-Copin, Ann. Math. 2021, arXiv:1912.07973）因四次耦合 irrelevant、Gaussian FP **IR-吸引**；YM 的 Gaussian FP **UV-吸引**，相互作用在 IR 端生成。平凡性定理**不适用**。REPORT §2。
- **OPE 不判别**：规范不变 OPE 系数自由也不为零；判别需非零**凝聚 `⟨tr F²⟩`** 或 `g_ren≠0`（皆 IR 量）。REPORT §4。
- **机器检验（新）**：`NonTrivialCore.lean`（core Lean，`#print axioms` **无公理**）证：Rademacher `κ₄=−2≠0`；**偶测度可有 `κ₄=0` 而 `κ₆≠0`**（`four_point_insufficient`）。⟹ **"连通 4 点单独见证"是错的**，须控制 `n≥3` 全体。已据此修正 REPORT §6.1/§6.4。完整证明稿：`PROOF.md`。

## 引理链状态

| 引理 | 状态 |
|---|---|
| C：Gaussian ⟺ 高阶累积全零（偶性消奇阶；**每偶阶 4,6,… 都需**；4 点单独不足） | **Lean 部分已证**（`NonTrivialCore.lean`，无公理）；解析方向待 Mathlib |
| C′：`NonTrivial ↔ ∃ n≥3` 非零累积见证（标量） | 代数骨架已 Lean；解析方向待 Mathlib |
| G：`g_ren≡0 ⟺` Gaussian（纯规范） | 可证（定义层） |
| S：`g=0` 在 scheme 变换下不变 | 平凡（已验证） |
| AF：`β=−b₀g³+…`，`b₀>0` | 2026 输入 |
| DT：`Λ=μe^{−1/(2b₀g²)}`，`Λ>0 ⟺ g_ren≠0 ∀μ>Λ` | 定义层 |
| **N1**：短程 `g_ren(p)` 一致收敛、匹配一圈 | **待审计（2026 UV 输入）** |
| **N2**：极限沿 `Λ>0` 非平凡轨迹、非自由子列 | **待审计（轨迹定义）** |
| **T**：定理 N（非平凡性） | N1+N2 ⇒ 成立 |

## 开放 / 假设 / 风险

- **真正的缺口全在 N1、N2**（审计问题），不在新数学。三问（REPORT §5.2）：
  1. 迹反常是**结构方程**还是**定标**（`β≠0` 且 `⟨tr F²⟩≠0`/`g_ren≠0`）？自由场也满足结构方程。
  2. 是否显式要求 `Λ>0` / 非零重整化耦合，从而排除 `g₀≡0` 平凡（Gaussian）轨迹？
  3. **子列**极限 + 未排除自由塌缩 ⟹ 陈述可能被自由理论满足（**vacuity 风险**，同 `formal-conjectures` 教训）。
- **诚实判断**：非平凡性**弱于 MG**（只需一个非零尺度，不需一致指数速率），但若只拿到结构性方程，则与 MG 同源于 IR 障碍。
- **未核实**：2026 原始论文在知识截止之后；§5 全部为据 `docs/yang-mills-prior-work.md` 的转述，arXiv 编号未独立确认。

## 下一步（最小动作）

1. **P0（Lean）✅ 核心完成**：`NonTrivialCore.lean` 已编译、无公理（`conn4_eq_zero_iff`、`radem_K4_ne_zero`、`four_point_insufficient`）。**待补**：接 Mathlib 证引理 C 解析方向（`κₙ=0 ∀n≥3 ⟹` Gaussian）；`NonTrivial` 与"`n≥3` 全体非零累积"等价化（**勿收紧到 n=4**）；自由 GFF `¬NonTrivial` 的 falsification（需 OSforGFF 构建）。
2. **P1（陈述层）**：见证从标量 `φ` 换成 `g_ren`；证 `IsGaussian ⟹ g_ren=0`（引理 G）。→ 纯规范的**非空洞**陈述。
3. **P2（审计，读不建）**：读 `yang-mills-group/YM` 陈述签名 + 16 输入，回答 §5.2 三问（对应 `docs` M0）。
4. **P3/P4（协同）**：把 N1（短程收敛匹配一圈）推成引理（UV 侧）；用 `numerics/` 的梯度流耦合做非严格数值 certificate。

## 与其他 session 的关系

- **正交**于主 session 的弱→强 RG 交接（`research/physics_handoff.md`，控制量 `q=m/μ`）：那攻 **IR/MG**；本文的非平凡性若 P2 判为"轨迹非平凡 + UV 连续性"，则**在 UV 侧、不需 MG**。二者共享 `g_ren`/scheme-无关性的语言（本文 §3.3 已对齐 handoff §2）。
