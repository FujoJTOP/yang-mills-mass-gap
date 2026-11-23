# LOWER-AMPLITUDE-MECHANISMS — can any correlation/spectral-support LOWER-bound machine reach the edge and give `m₁>0`?

> 2026-11-23 · 八个候选（GKS / FKG / Simon–Lieb / reverse-BL / Carathéodory / RP-cone / Ginibre / GCI）逐一核验。**判词：无一 5/5。** 它们是"下界机器瞄准上界目标"——质量隙需要 connected 关联的**上界**（衰减率）或支撑底 `inf supp ν>0`，而它们的原生输出是关联/能量的**下界**，物理上推向**轻模 = anti-gap**。唯一"类型对"的下界是**残数权重** `f₁²≥r₀`（经夹逼 `f₁²/χ≤m₁²`）配 `χ≤1/λ` 得 `m₁²≥r₀λ`；但由 Schur 凹性，正性/矩/积分数据可证地**只给 `f₁²` 的上界**（`f₁²≤M₀` 原子≤总质量；`f₁²≤χm₁²`，由 `χ≥f₁²/m₁²`），八个机制**无一**供应下界。

## §1 三分：关联积分下界 ≠ 原子权重下界 ≠ 边位置下界
任务的关键。质量隙 = `ν([0,Δ²))=0`（**位置**）。三条逻辑路线：
- **(A) 经权重**：`m₁² ≥ f₁²/χ`（G2A-INFSUPP）。需权重地板 `f₁²≥r₀` **且** `χ≤1/λ`。这是**唯一**"下界 → `m₁>0`"的正道。
- **(B) 直接支撑**：证 `supp ν` 离 0。**任何**正性/连续性数据都做不到：自由无质量场 `ν` 支在 `[0,∞)`、`f₁²=0` 是合法 OS 理论——凡正性可证的支撑下界必在 `g=0` 失效（master obstruction，SCHUR-CONVEX / GOLDS-POS-MASS）。故 (B) 的正性路线**被自由场封死**。
- **(C) 关联积分**：`M_k=∫s^k dν`、`∫φdν`（φ≥0）全是 `ν` 的 Schur-凹泛函（矩引理 GAP-G2-edge），给原子权重**上界**、且 `δ₀+δ₂` 与 `2δ₁` 同矩异边——积分**永不见**原子/边。下界关联 = 慢衰减 = `m₁` 上界 = **anti-gap**。

所以缺口是「**残数权重地板 `f₁²≥r₀`**」（= `H*`，RESIDUE-CONDITIONAL），而它是**非矩、非积分**的独立数据。下表逐机制问：给不给下界？在 4D gauge 成立否？触边否？过 5 否？

## §2 逐机制
| 机制 | 原生下界？ | 4D SU(N)+Wilson 成立？ | 触原子/边？ | 5/5 |
|---|---|---|---|---|
| **GKS**（Griffiths'67 / Kelly–Sherman'68） | 是：`⟨σ_A⟩≥0`、`⟨σ_Aσ_B⟩≥⟨σ_A⟩⟨σ_B⟩` | **否**：铁磁/XY 结构；SU(N) 受挫 `cos` 和、torsion 无格序（PROOF-COMPARISON-FKG） | 否（关联积分，且方向反） | 否 |
| **FKG / 正关联**（FK G'71） | 是：`Cov≥0` | **否**：无格序；且**错号**（`Cov≥0⇒X≤0`，降 `G″`，C6） | 否 | 否 |
| **Simon–Lieb** | 否（**方差上界**，`Var≤ΣE[∂²]`） | 是但 sign-blind | 否：`=Var(a)≤1/λ=H1`（C6 Bin 2） | 否（#2/#3/#5 败） |
| **reverse/inverse BL** | 是：log-凹测度上 `E[fg]≥…` | **否**：需 `Hess V⪰0`；`λ_min=−g\vert B\vert<0`（Savvidy 墙） | 否（关联，且方向反） | 否 |
| **Carathéodory / 正实函数** | 是：`Re f(z)≥0`、Pick PSD | 是（`f=∫dν/(s−z)` 恒真） | 否：仅重编码 `ν≥0`；自由场盲；矩/Schur 全给上界 | 否 |
| **RP-cone** | 仅 PSD：`⟨θF,F⟩≥0` | 是（OS 定义） | 否：PSD ⇒ C-S **上界**单矩阵元，造不出个体 `≥c>0`；RP 尺度协变、任意小隙都有 RP 模型 | 否 |
| **Ginibre / 高自旋关联** | 是：铁磁 `⟨S_A·S_B⟩≥0` | **否**：需正（铁磁）耦合 + 格序；受挫非吸引 | 否（关联积分，方向反） | 否 |
| **Gaussian correlation**（Royen'14） | 是：`P(∩Cᵢ)≥∏P(Cᵢ)` | **否**：仅高斯测度；Wilson 非高斯 | 否：凸集联合概率，非谱陈述 | 否 |

逐条精确死法：**GKS/Ginibre/FKG** 给的是"关联不衰减"的**正**下界——`0≤Cov` 即 `X=−2Cov≤0`，把 `G″` 往**下**压（C6 已测 C-S `0.40>W_geo 0.32`，`G″≳−0.08`），**方向反**，且需的铁磁/格序 SU(N) 没有。**Simon–Lieb / BL** 是方差**上界**，逐字 `Var(a)≤1/λ`：这不是新的下界，是 H1 的另一件外衣。**reverse-BL** 唯一的"下界"需 log-凹，而 log-凹性恰是 Savvidy 墙（`Hess` 不定）；即便成立，`E[fg]≥…` 仍是 anti-decay、触不到原子。**Carathéodory** 的 `Re f≥0` 是平凡正性（自由场饱和），Pick/Schur 算法与矩同阶——**全上界**。**RP-cone** 是 PSD 锥，C-S 只把**非对角元从上方**压住，**绝不可能**生产单元素 `≥c>0`（GAP-G2-edge §"RP is a PSD condition"）；且 RP 尺度协变 ⇒ 不能强制 `m>0`。**GCI** 是高斯专属的联合概率下界，测度错、方向错、非谱。

## §3 唯一的"类型对"下界，以及为何八者都到不了
缝隙的正道 (A) 要 `f₁²≥r₀>0`。但所有 `∫φdν`（φ≥0，含任意关联积分、任意 Carathéodory 值）对 `f₁²` 都是 Schur-凹 ⇒ 给**上界** `f₁²≤χm₁²`、`f₁²≤M₀`；下界 `f₁²≥r₀` 是**非积分**数据（GAP-G2-edge："edge 下界是独立非矩数据"），且必在 `g=0` 失效（自由场 `f₁²=0`）。八个机制产出的是积分/方差/PSD——**全是 Schur-凹或 sign-blind**，故按 SCHUR-CONVEX / INEQ-T3 的类定理，**没有任何一个**能给 `f₁²` 下界；更遑论位置 `m₁`。

**关于"位置而非权重"**：除经 (A) 的 `f₁²/χ` 夹逼外，无任何直接位置路线。(B) 被自由场封死；(C) 位置是支撑底、非矩可定（`δ₀+δ₂` vs `2δ₁`）。位置下界**只能**从权重下界 + `χ` 天花板合成（G2A-INFSUPP），而权重下界缺席。

## §4 裁决
**无 `5/5 ACHIEVED`。** 八个是八种方言的一堵墙：下界机器、上界目标。最诚实的一句——任务设想的"下界振幅"在**物理方向上是 anti-gap**：更强的关联 = 更轻的模。唯一能被下界推动的坐标是残数**权重** `f₁²`，而它恰是所有正性机制被 Schur-凹性禁止触碰的 Schur-凹对象。精确缺失步仍是 `H*`：`f₁²≥r₀` 在 β=2.5、`(L,a)` 一致地成立。无机制，无捷径，无捏造。

*Refs.* SCHUR-CONVEX; GAP-G2-edge（矩引理、残数–能级–隙、edge 下界=非矩数据）; G2A-INFSUPP-DIRECT（`f₁²/χ≤m₁²≤M₀/χ`）; GOLDS-POS-MASS（唯一正性→原子桥钉 0 质量）; REDTEAM-SPEC; C6-DEPENDENCE-CLASS（FKG 错号、Simon–Lieb=H1、BL 需 `Hess⪰0`）; PROOF-COMPARISON-FKG（torsion 无格序）; C8-SPECIES-NOGO; RESIDUE-CONDITIONAL (`H*`). External: Griffiths, J. Math. Phys. **8** (1967) 478; Kelly–Sherman, J. Math. Phys. **9** (1968) 466; Fortuin–Kasteleyn–Ginibre, Physica **57** (1972) 536; Ginibre, Commun. Math. Phys. **16** (1970) 310; Simon–Lieb 条件方差不等式（cf. Efron–Stein, Ann. Statist. **9** (1981) 586）; Brascamp–Lieb, JFA **22** (1976) 366; Prékopa 1971; Royen, Far East J. Theor. Stat. **48** (2014) 139; Osterwalder–Schrader, CMP **31** (1973) 83; Savvidy, PLB **71** (1977) 133.
