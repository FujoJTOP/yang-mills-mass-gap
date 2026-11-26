# S01-ANATOMY — `(★)` 的九个符号逐个钉死

> 对象：`research/v3/C2-BALABAN-UNIFORM.md:214`（逐字）
> `(★)  Var_{μ(·|U′)}(⟨Pᵀ∇S_j, v⟩)  ≤  E_{μ(·|U′)}⟨PᵀH_jP v, v⟩ + W_j(v) − δ_W      (δ_W > 0),`
> 以及 `:218` `G″_j(v) = E[H_j] − Var_j + W_j ≥ δ_W`、`:211` "uniform in `j,U′,a`"。
> 本级判词：**(a) 为主 + (b) 一处（`δ_W` 在一般理论中只有模型级定义）+ (c) 量纲齐**。

## 0. 引文核验（工具状态）

| 检索式 | 工具 | 命中 | 判定 |
|---|---|---|---|
| `Test-Path research/v3/species/C2-BALABAN-UNIFORM.md` | pwsh | 0（文件实存 `research/v3/C2-BALABAN-UNIFORM.md`） | 题给路径少一级 |
| `read C2-BALABAN-UNIFORM.md` 180–259 | read | 有效 | 214 行逐字命中所引 |
| `grep 'W_j\|δ_W\|S_j\|H_j\|U′\|Var_j'` on `G1-L3-UV-FORCE-BOUND.md` | grep | 3 | 有效 |
| `grep 'Pᵀ\|Hessian\|collective mode\|transverse projector'` 全 `research/**.md` | grep | 1271（截断 250） | 有效 |
| `grep 'δ_W'` 全 `research/**.md` | grep | 483（截断 250） | 有效 |

无 TOOL-FAIL。以下"未定义"均为**检索执行成功后的 0 命中或仅模型级命中**。

## 1. 符号表（逐字定义 + 行号 + 量纲）

| 符号 | 逐字定义 | 行号 | 量纲 | 标签 |
|---|---|---|---|---|
| `j` | "At shell `j` (block spacing `a_j`)…"；`β_j = 2N/g_j²`，"one-loop `β_j = 4Nβ₀ j`"；`(★)` "uniform in `j,U′,a`" | `S6-DOBRUSHIN.md:64,69-70`；`C2:211` | 无量纲（标号） | **[T]** |
| `U′` | "the **conditional block measure** given the boundary configuration `U′ = U_{Λ\Λ_B}`"；亦作 coarse field："`U′` enters the physical collective mode `v`" | `S6:65-66`；`S3-UV-FESHBACH.md:71` | 群值（无量纲） | **[T]** |
| `S_j` | "`S_j = Wilson action = β_j Σ_{p ∩ B≠∅} (1 − (1/N)Re tr U_p)`"；`∇S_j` 对**块内链路** `U_B` 求导，`H_j=Hess S_j` 同变量 | `S6:67,70` | `S_j` 无量纲；`∇S_j` 无量纲（`dU` 为 Haar）| **[T]** |
| `P` | "`P` is the block-averaging (coarse-graining) map, `Pᵀ` its adjoint"；线性（"the campaign's `⟨Pᵀ∇S,v⟩` already is"）；"harmless" **仅对常值 `P`** | `S6:70`；`G2-TWO-PLAQUETTE.md:328`；`G1-BLOCK-IDENTITY.md:165-167` | 无量纲线性算子 | **[T]** |
| `H_j` | "`H_j = Hess S_j` is the **flat (Wilson) Hessian**"；`⟨PᵀH_jP v,v⟩ = vᵀPᵀH_jPv` 为 coarse 方向上的 Hessian 双线性型 | `S6:70`；`G1:93-98` | 无量纲 | **[T]** |
| `W_j(v)` | "`W_j` is the **Haar/Jacobian curvature** budget"；`W_j = E[Φ″] − Var(Φ′)`，`Φ=−log J`；且 "`W_j` means `W_j = W_geo + X`，`X = −2Cov(G,Φ′)`"（`X` 号不定） | `S6:71,83`；`RECONCILIATION-STAR.md:50-53` | 无量纲（曲率/单位方向） | **[T]**（定义）/[C]（`W_geo≥N/12`，`S3:152`） |
| `δ_W` | "`δ_W(β) := 2·inf_{U′,v} W_geo > 0`"；"`δ_W(β)→0` as `β→0`"；"the worst shell is the IR head `β≈2.5`, where `δ_W ≈ 0.13`"；"increases along the tower" | `S2-BOUNDARY-COUNTEREXAMPLE.md:201-206,239,36` | 无量纲；量级 `O(0.1)`@head，深 UV `≳N/12` | **[C]** |
| `μ(·|U′)` | "`μ_j(·|U′) ∝ exp(−S_j(U_B;U′)) ∏_{ℓ∈Λ_B} dU_ℓ`"——条件化在**块外/边界场 `U′`** 上 | `S6:66-67` | 概率测度 | **[T]** |
| `v` | "along a **unit** collective mode `v`"；"all directional, along a field `v = v(ξ)` in the tangent space `TΩ`" | `S6:76`；`G1:90,98` | 无量纲，`‖v‖=1` | **[T]** |

## 2. 逐条回答（题面 1–9）

1. **`j` = 尺度/壳标号，不是块标号。** 块用 `B`、`Λ_B`（`S6:64-65`；`S4-DECIMATION.md:118` 对块用 `B`）。`j` 同时给出该壳的耦合 `β_j = 2N/g_j²`。**[T]** 语料内 `j` 从未作块标号使用。
2. **`U′` = 已被积掉/留在块外的其余链路构型 `U_{Λ\Λ_B}`**（边界条件），同时就是 RG 的 coarse 场。它只作为 `S_j` 的参数进入（"The boundary enters **only** as a parameter of `S`"，`G1:88`）。**[T]**
3. **`S_j` = 壳 `j` 上、限于块 `B` 的 Wilson 作用量**（`S6:67`）。`∇S_j` 对**块内链路 `U_B`**（fibre 变量）求导，`U′` 固定。**[T]**
4. **`P` = 块平均/粗粒化（线性）算子，`Pᵀ` 其伴随**（`S6:70`）。`w := Pv` 是 coarse 方向（`G1:165`）。**[T]**
5. **`H_j` = `S_j` 的 flat(Wilson) Hessian**（`S6:70`）。关于 `P` 两侧看似不对称：**不是笔误，是张量类型使然**——`⟨Pᵀ∇S_j,v⟩` 里 `∇S_j` 是 **covector（1 个 `P`）**，`⟨PᵀH_jP v,v⟩` 里 `H_j` 是 **双线性型（2 个 `P`）**；两者都等于 `⟨H_j w,w⟩` 型的同一个 coarse 二次型，`G1:165` 明写 `⟨Pᵀ∇S,v⟩ = ⟨∇S,Pv⟩ = ∂_{Pv}S`。**[T]**
6. **`W_j(v)` = Haar/Jacobian 曲率预算**，且是**完整版** `W_j = W_geo + X`，`X` 符号不定、不可被吸收（`RECON:50-53`）。**符号：正为主**——`W_geo ≥ N/12`（`S3:152`）、Lemma B 给 `W≥δ_W`（`G1-L3-UV-FORCE-BOUND.md:3`）；但含 `X` 后 `W_j` 本身**不保证正**。**[T]/[C]**
7. **`δ_W`：无量纲、`O(0.1)`@`β=2.5`，是 `β` 的函数（`δ_W(β_j)`），不是与 `β` 无关的常数**；`→0` 当 `β→0`（自由场酸测）；沿 tower 增大（`S2:201-206,239`）。深 UV 由 `Φ″⪰N/12` 给 `δ_W ≳ N/12`（`S3:268`）。**[C]**
8. **`μ(·|U′)` = 块内条件测度**，密度 `∝ e^{−S_j(U_B;U′)}∏dU_ℓ`，条件化在块外/边界场 `U′` 上（`S6:66-67`）。**[T]**
9. **`v` = coarse/collective 方向，`‖v‖=1`**（`S6:76` "unit collective mode"）；`v` 可与场有关 `v=v(ξ)`，那样 `div v ≢ 0` 并产生 `Δ`（`G1:173-175`）；`(★)` 的"常值 `P`"一栏要求 `v` 常值（`G1:162,166`）。**[T]**

## 3. 判死线

### (a) 可钉死者（7/9 为完全 `[T]`）
`j, U′, S_j, P, H_j, μ(·|U′), v` — 七项均有语料内**逐字定义 + 行号**，且与 `(★)` 的用法**一一吻合**（无一处需外推）。

### (b) "未定义"者（2 项，附最可能读法）
1. **`δ_W` 在一般理论中无定义。** 检索式 `grep 'δ_W'`（工具 OK，483 命中截断）显示：**唯一的显式定义**出现在**可解的单 plaquette 模型**——`δ_W(β)=2 inf_{U′,v} W_geo`（`S2:201`），且同文件自标 `[C]`（"+[V]"）。在完整 4D 块上，`δ_W` 只作为"`(★)` 的严格余量"**被使用**（`S4:270` "`δ_W=O(1)`（Lemma B, scale-independent）"、`S3:37` 作为待证目标）。**最可能读法**：`δ_W := 2 inf_{U′,v} W_geo(β_j)` 的块版推广（即"Haar 曲率的一致最坏值"），量级 `O(0.1)`。
   ⚠ **同名不同物警告**：`research/gluing/**` 与 `physics_handoff.md:45` 的 `δ_W` 是**弱侧证书常数**（`cert μ ≥ δ_W`），与 `(★)` 的余量**不是同一对象**，不得互相代入。
2. **`P` 的"协变 vs 常值"两读。** `G1:165-167` 只证 `P` **常值**时 `(★)` 无害；`RECON:237` 另列"covariant `P` 的 `O(g_j)` transport-defect 界"为**未证**项。**最可能读法**：`P` 线性常值（`G2:328` 断言 `(★)` 本就是线性类）。

### (c) 量纲：**齐**（不指向"式子写错"）
四项全为**无量纲**（格点单位 `a=1`；`S_j` 无量纲、`v` 单位化、`Φ=−log J` 无量纲）：`Var(力)`, `E[H]`, `W_j`, `δ_W` 同量纲。故**无** `Var` vs `E[H]` vs `W_j` 的量纲冲突。
但两点**非量纲**的结构问题需记账（`[C]`）：
- **`a` 不出现在式子里**：`uniform in `a`` 是靠 `β_j(a)`、per-block 归一隐式承载，不是 `(★)` 自身显式携带——`a`-均匀性因此**不可从该式读出**。
- **`W_j` 是"另一种对象"**（Haar 曲率预算，含不定 `X`）被加进 `Var` 与 `E[H]` 的比较，这正是 `A1-twosided-UV.md:3` 与 `A2` 判 `(c)` 为"空/需 G2 修形"的那一点；不是量纲病，是**类型病**。

## 4. 一句话总结

九个符号里七个 `[T]` 钉死（`j`=壳、`U′`=块外场、`S_j`=块 Wilson 作用量、`P`=块平均、`H_j`=flat Hessian、`W_j`=Haar 预算 `W_geo+X`、`μ(·|U′)`=块内条件测度、`v`=单位 coarse 方向）；`δ_W` 只有**模型级**定义（`[C]`，β-函数、`O(0.1)`、`→0` 当 `β→0`），并与 gluing 的证书常数**同名不同物**；`P` 有协变/常值两读。**量纲齐**，`(★)` 无书写量纲错；代价在 `a`-均匀性未显式化与 `W_j` 的类型混入。

---

*纪律：所有断言标 `[T]/[C]`；未标 `[T]` 者未证。检索执行状态见 §0。*
