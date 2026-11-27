# NNN3-RECON — 把「重构」闭成 `[T]` 并钉死它的范围：不夹带 `Δ>0`

> 轮 NNN3 · 子代理 `nnn3-recon` · `D:\Dev\math`。只交 `[T]`/`[NG]`。编号**按题面**：OS0 温性/增长、OS1 平移、OS2 反射正性、OS3 `O(4)`、OS4 聚类（仓内 `docs/conditional-theorem.md:45` 与 G1 反号，跨文件必带号）。**窗**：`L` 固定、`a→0`；`L`-一致是 `UVS-L`，另账；OS4 的定量路线在 `β≥2.5`。零新 MC、未编译 Lean、未改既有文件、未 commit。
> 题面「OS0–OS3 + 给定 UVS 的 OS0」读作**同一条 OS0 的两种记法**：前件集 = `OS0–OS3`，其中 OS0 的供给形 = 固定 `L` 的 `UVS(L)`。

## §0 判词（四条，各带下一手）

1. **[T] 重构可闭，但须闭成两条**：`Re-A = OS0–OS3 ⇒ Wightman`；`Re-B = Re-A + OS4 ⇒ ker H=ℂΩ`。题面的六项输出（`ℋ_phys / H≥0 / HΩ=0 / ker H=ℂΩ / Poincaré 表示 / 谱条件`）**恰是 A∪B**——只给 OS0–OS3 时 `ker H=ℂΩ` **不在**结论内（反例：两纯相混合，OS0–OS3 成立、`ker H` 二维）。→ **下一手**：AX 接口签名分 `Re-A`/`Re-B` 两行，`ker H=ℂΩ` 只许挂 B。
2. **[T] `Δ>0` 不是重构的推论，任何一条 OS 公理都不产它**：无质量自由标量满足 OS0–OS4（含唯一真空、幂律聚类）而 `Δ=0`；故 `AX ∧ OS4 ∧ 唯一真空 ⇏ GAP`。→ **下一手**：`Δ>0` 只写进 `GAP` 合取项，禁出现在 `Re` 输出列；`Clay = AX ∧ NT ∧ GAP` 三合取不动。
3. **[NG] 仓内两处错位成立，逐字复核**：(A) `research/gluing/REPORT.md:1348–1350` 盒逐字「`I-RECON: OS0–OS4 ⇒ Wightman + 谱条件 + 唯一真空 + Δ>0`」——输出列夹带 `Δ>0`，该盒为 `[NG]`；(B) `research/lean-infra/REPORT.md:486` 逐字「**Level B（不做）**」⇒ 第 4 项产出非机器。→ **下一手**：盒拆「`Re`（无隙，六项）」+「`GAP`（独立输入）」两格；接口注「文献级 / 未机器化」。
4. **[NG] Mathlib 侧不可得，须自建，量级 = 上千行**（非几十/几百）；陈述骨架可几百行但只是定义。→ **下一手**：若机器化，先出「文件-行」报价与 `axiom` 清单，不先铺 `sorry`。

## §1 精确陈述（前件逐条，含温性）

**定理（Osterwalder–Schrader，CMP 31 (1973) 83；CMP 42 (1975) 281；`[T]` 文献级、未机器化）。** 设 `S_0=1`、`S_n∈𝒮′(ℝ^{4n})` 且置换对称，并满足下表 `OS0–OS3`；则存在：可分 `ℋ`；单位向量 `Ω`；`E(4)=ℝ⁴⋊O(4)` 的强连续酉表示 `U`；自伴 `H≥0`（`U(t,1)=e^{−tH}`，`HΩ=0`）；平移联合谱 `spec(P)⊂V̄₊`；算子值温分布场 `φ`（共同稠密 `U`-不变域 `D∋Ω`，实 `f` 自伴），使 `S_n(f_1,…,f_n)=⟨Ω,φ(f_1)⋯φ(f_n)Ω⟩`；经 `t↦it` 解析延拓，`φ` 为 Wightman 场（Poincaré 协变、谱条件、局域性、温分布、真空循环）。**加 OS4 ⇒ 唯一真空 `ker H=ℂΩ`。**

| 前件 | 逐字内容 | 本项目状态 | 出处 |
|---|---|---|---|
| **OS0**（唯一 `∀n` 条） | 温性 + `n!` 增长：`S_n∈𝒮′(ℝ^{4n})` 且 `∃C,p: \|S_n(f_1⊗…⊗f_n)\| ≤ C^n n! ∏‖f_i‖_{𝒮^p}`（⇔ `Z(f)=Σ_n i^n/n! S_n(f^{⊗n})` 在 `0` 邻域解析） | **`[NG]`**；供给形 = 固定 `L` 的 `UVS(L)`：`(U1)` `a`-一致半范界、`(U2)` 同方案减除、`(U3)` `n!` 增长 + Minlos；**温性不由 `(C_c^k)′` 弱收敛自动继承**（`Σe^kδ_k`） | `AA1-OS0.md:16,24–27`；`glue-a-limits/REPORT.md:58–62` |
| **OS1** | 平移：`S_n(·+a)=S_n(·)`，`a∈ℝ⁴` | **`[T]`**（群像 `aℤ⁴` 稠密 + 弱闭） | `GGG3-OSCLOSED.md:18` |
| **OS2** | 反射正性：`Σ_{i,j≤k} S_{i+j}(θf_i,f_j) ≥ 0`，`f_i∈𝒮` 支于 `t>0`，`θ(t,x⃗)=(−t,x⃗)` | **`[T]`**（RP 锥弱闭；Lean `rp_transfer`、`rp_l2_closure_trunc`） | `GGG3:19`；`glue-a-limits:145–177` |
| **OS3** | `O(4)` 不变：`S_n(R·)=S_n(·)` | **`[NG]`**（= `SYM-4D`/`IN-2`；`ℓ=4` 各向异性 Gaussian 过其余四条仍破它） | `AA2-OS3.md:11,20` |
| **OS4** | 聚类（`⟨A T_t B⟩→⟨A⟩⟨B⟩`） | **`[T]`⟸`GAP`**（`GlueA.spectral_decay`）；仅 `ker H=ℂΩ` 需要 | `conditional-theorem.md:294` |

**不在结论内**：`Δ>0`／`spec(H)∩(0,Δ)=∅`／非平凡／理论唯一（子列）——四条全是并列合取项，不是重构的输出。

## §2 能否 Lean 化（Mathlib 逐条核）

| 对象 | Mathlib 现状（`mathlib4` rev `79d0395`，v4.33.0-rc1） | 判 |
|---|---|---|
| Wightman 场 / OS 重构 | `Wightman\|Osterwalder\|Schwinger` → **0** | 无 |
| Poincaré 群 / Lorentz / `O(1,3)` | `Lorentz\|LorentzGroup` → **0**；`Poincar` **26** 条全是递推定理/引理/猜想/度量 | 无 |
| 正定正交群 | 仅 `Matrix.orthogonalGroup`（矩阵子幺半群；无 Lie/酉表示） | 不够 |
| 温分布 | `TemperedDistribution` **109** 条（Fourier/导数/支撑齐备） | 有 |
| 谱测度 / PVM / Stone | `SpectralMeasure\|ProjectionValued` → **0** | 无 |
| Laplace 变换 | `[Ll]aplace` **3** 条（Laplacian 记号、行列式展开、文献名） | 无 |
| 管域定理 / 边值分布 / 多复变 | `tubeDomain\|TubeDomain` → **0** | 无 |

**判**：`[NG]`（Mathlib 侧）⇒ 须自建，**量级上千行**。分层：①**陈述骨架**（`structure WightmanTheory` + 场作 `𝓢(ℝ⁴)→` 线性映射 + 五条公理作字段）约几百行，但**只是定义**、不构成机器化；②真证明卡三处：`H`/`P` 的无界自伴谱定理 + Stone（仓内 `G-4/G-5`）、`E(4)`/`O(1,3)` 的强连续酉表示、管域定理 + 分布边值 + Laplace 解析延拓（**局域性也依赖它**：谱延拓上置换对称 = Jost 型论证）。**量级锚**：仓内 Euclidean 半边 `vendor/OSforGFF/OSforGFF/` 已 58 文件 / 28,017 行（本轮机读；`lean-infra:558` 记 31,429 行）——Level B 不低于此，且三处各自上千行。

## §3 正确的表述句（可直接入链）

> **`Re`（`[T]`，文献级、未机器化）**：`OS0–OS3`（OS0 = 温性 + `n!` 增长，本项目由固定 `L` 的 `UVS(L)` 供）⇒ 存在 `ℋ`、真空 `Ω`、`E(4)` 的强连续酉表示、`H≥0`（`HΩ=0`）、谱条件 `spec(P)⊂V̄₊`、算子值温分布场 `φ`，`S_n` 为其真空期望，解析延拓给 Wightman 场；**`ker H=ℂΩ` 只在加 OS4 后成立**；**`Δ>0` 不在其中**——它与 `Re` 并列，属 `GAP`。

**禁写**：「`OS0–OS4 ⇒ Wightman + Δ>0`」（`gluing:1348–1350` 原盒，`[NG]`）；「重构已 Lean 化」（`lean-infra:486` 逐字「Level B（不做）」）。
**接口三条**：①`Re` 输出列止于六项，`Δ` 单列 `GAP` 行；②`ker H=ℂΩ` 与 `OS4` 同现，禁单独引用；③OS0 一律写「`UVS(L)` 的 Minlos 面 + 温性」，禁写「弱收敛已给温性」。→ **下一手**：位置图 `AX` 行按 `Re-A / Re-B / GAP` 三分改写。

## §附 AGENTS §1（工具状态 / 检索式 / 命中 / 失败标记）

| 工具 | 状态 | 检索式（逐字） | 命中 |
|---|---|---|---|
| `grep` | **OK**（有效响应，无 `provider failure`/`ENOENT`/空响应） | ①`UVS` @ `research/v3` ②`I-RECON\|AXFINAL\|Wightman` @ `research/v3` ③`Wightman\|Osterwalder\|Schwinger` @ Mathlib ④`Poincar\|LorentzGroup\|Lorentz` @ Mathlib ⑤`OrthogonalGroup` @ Mathlib ⑥`TemperedDistribution\|tubeDomain\|SpectralMeasure\|ProjectionValued` @ Mathlib ⑦`[Ll]aplace` @ Mathlib ⑧`Wightman\|Poincar\|Minkowski` @ `lean/` ⑨同式 @ `OSforGFF/OSforGFF` | ①176 ②44 ③**0** ④26 ⑤16 ⑥109 ⑦3 ⑧8 ⑨6 |
| `glob` | **OK** | `research/v3/species/*.md`；`research/v3/species/NNN*.md` | 994（显示 100）/ **0** |
| `read` | **OK**（带行号、非空、未截断） | `GGG5-RECON`(47)·`AA7-AXFINAL`(79)·`gluing:1335–1364`·`lean-infra:470–559`·`conditional-theorem:285–324`·`G1-CONTINUUM-OS`(15)·`B2-bypass-gap:1–14` | — |
| `pwsh` | **OK**（exit 0） | Mathlib 定位、`lake-manifest.json` 的 `rev`、`OSforGFF/OSforGFF` 行数 | rev `79d0395`；58 文件 / 28,017 行 |
| `web_search` | **OK**（返回有效源） | OS I/II 的 DOI | `10.1007/BF01645738`、`10.1007/BF01608978`（OS II 摘要逐字 "analytic continuations to a relativistic field theory"） |

**失败标记 0**（无 `TOOL-FAIL`；上表 `0` 均为工具正常下的真实 0 命中）。**直读/转引**：§1 定理按 OS73/75 标准形（OS I 正文未逐字核、只核 DOI；OS II 读到摘要）；仓内 `file:line` 全部直读上列文件，无转引未标者。
