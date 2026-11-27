# GO1-RPPROOF — 撤回任务的存根（stub）

> **声明**：本文件是任务「对具体 4D 规范改进作用量（LW/Iwasaki/DBW2）证 RP 或证 no-go」的**存根**。该任务已被 lead 从关键路径撤下（依据 `GN3-NEARCRIT.md` 定理 B + `gluing/REPORT.md:79`：RP 是本项目所用 Wilson 测度的直接性质，不是项目欠的账）。**本文件不含判词、不含 `[T]`/`[NG]`、不含任何"4D gauge RP 可/不可证"的结论；不得引用为结论。**
> 本轮：未改任何既有文件；未 `git commit`；无 MC。Clay = 0%。

## §1 已核事实（只摘录，不判级）

| # | 摘录 | 出处 | 读入级别 |
|---|---|---|---|
| 1 | 「作用量取 Wilson（基本表示），不取 Manton/adjoint/improved。理由：Wilson 作用量是**已知反射正性**的最简选择（Osterwalder–Seiler 1978），而 RP 是承重的（§3 GL2）。改进作用量（Symanzik）虽有更好的 `a²` 行为，但其 RP 需要单独证明，**不在第一阶段**」 | `research/gluing/REPORT.md:79` | 逐字 |
| 2 | R16 风险行：「Wilson vs **改进作用量：某 session 用 Symanzik** … 低 ｜ **RP 未证**」 | `research/gluing/REPORT.md:1504` | 逐字 |
| 3 | 「`𝒞₊ ⊆ {RP}` 是包含，不是相等。负系数作用量**可能**仍然 RP（Lüscher–Weisz 只说 "in general"）」 | `research/v3/species/GL2-RPORDER.md:284` | 逐字 |
| 4 | `T = e^{−aH}`、`0 ≤ T ≤ 1`、`H ≥ 0`；`spec T ⊂ [0,1]`（R1/R2） | `research/v3/species/B01-RP.md:26–28` | 逐字 |
| 5 | 「When operators of **dimension 6** are added to the standard Wilson action in lattice gauge theories, **physical positivity is lost in general**. We show that a transfer matrix can nevertheless be defined. Its properties are, however, unusual: **complex eigenvalues may occur** …, and **there are always contributions in the spectral decomposition of two-point functions that come with a negative weight.**」 | Lüscher–Weisz, NPB **240** (1984) 349, DOI `10.1016/0550-3213(84)90270-0` | **摘要逐字**（SciDirect+MathSciNet MR0762195 互证）；正文未核 |
| 6 | App. B.2.1 标题「**Unphysical poles in the propagator**」；正文片段「for the Iwasaki and DBW2 action, due to **the violation of physical positivity** described in appendix B, one expects also negative contributions」 | Necco, DESY-THESIS-2003-016 = arXiv:`hep-lat/0306005` | 目录+片段；正文未逐字核 |
| 7 | 「**for the DBW2 action** … we have calculated the transfer matrix in the Euclidean time direction for **free domain wall fermions** and shown that **negative eigenmodes** arise for `M>1` and produce oscillatory behavior」 | arXiv:`0710.0425`「Oscillatory terms in the domain wall transfer matrix」 | 摘要+正文片段（**范围：费米子转移矩阵，非规范转移矩阵**——与仓内既有转引的口径差在此） |
| 8 | 「In general, adding operators of dimension 6 to the standard Wilson action violates the physical positivity at finite lattice spacing. A transfer matrix …」；`c_j` 表：Iwasaki `0.331`、DBW2 `1.4088` | arXiv:`hep-lat/0208052` | 片段；**符号约定未逐字核** |
| 9 | 「physical positivity holds in **Wilson's** lattice gauge theories」；`Osterwalder–Seiler 1978` 仅仓内转引 | Lüscher CMP **54** (1977) 283, DOI `10.1007/BF01614090`；Ann. Phys. **110** (1978) 440, DOI `10.1016/0003-4916(78)90039-8` | 摘要逐字 / 转引 |

## §2 交接

- 若有后续：「**负本征模 ⇒ 谱非负性失败**」这条仍**可能**是一条独立有效的论证路线；但它现在**不在关键路径上**，未定级、不得作为判词引用。
- 未定级工作（**不得引用**）：确定性 free-field 脚本 `%TEMP%\go1_rp\{fy.py, fy2.py, run1.py … run6b.py}`（numpy 2.4.6 / scipy 1.18.1 / sympy 1.14.0，零随机、零 MC）。仅留档，未成文、未定级。
