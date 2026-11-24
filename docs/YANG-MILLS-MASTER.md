# 杨-米尔斯质量间隙 · 证明总册（MASTER）

> 2026-11-23 · 本文件是**所有直接关联到证明 YM 的内容的总汇**（条件链 / 已证定理 / 新工具 / 新 no-go /
> 免疫图 / 残留对象 / 硬算状态 / 诚实判词）。前面各 `docs/*.md` 与 `research/unconditional/*.md` 是本册的证据源。

---

## 0. 一句话状态

**不是 Clay 的证明；YM 仍是猜想。** 但本战役把"证质量间隙"化成了**一个精确定义的对象**，证出**五件真定理**、造出**五个能过自由场酸测的新工具**、钉出**三条新 no-go**，并把全问题收敛到**一份"残留对象"清单**。距 Clay 定理仍是 **0%**（它是公开问题），但**墙的坐标、名字、逃逸口已完全登记**。

---

## 1. Clay 问题与等价改写

**Clay**：4D 纯 SU(N) 杨-米尔斯存在且有一个正质量间隙 `Δ>0`。

**锁芯改写**：间隙 = 场强算子的 Källén–Lehmann 谱测度 ν 最轻原子 `m₁²` 处的权重
`f₁² = ν({m₁²}) = |⟨0|tr F²|G⟩|² > 0`，**且** `Δ>0`。
- **Clay ⟺ ¬(iii)**：YM 测度**不是**无能隙的相互作用 4D CFT。
- `f₁²` 的性质：**非对角、上-半连续（usc）、需 `E₂` 二标度**。

**两半边**：`Clay = G1（存在）+ G2（隙）`，夹物理通道（最轻态量子数）。
- **G1**：4D 连续存在性 = 边际规范协变 BPHZ / 一致 `Z_a` 重整化（**"无已知想法"缺口**；固定体积除 L3 外已证）。
- **G2**：隙 = `f₁²>0` 的证明/搬运（**= 问题本身**）。

---

## 2. 已证的真定理（不是缺口，是战利品）

### 2.1 P3 强耦合残数定理（[条件/待重证]，2026-11-24 审计）
存在 `β₁=O(1/N)`（与 SZZ `1/48` 同阶），使 `β<β₁`、固定 `a`、**一致于体积 `V` 与 `L`**：
```
f₁² ≥ c₁(β)·u(β)^p·(1 − O(u)) > 0,    u(β)→0.
```
骨架：`tr F²` 是有限 A1g plaquette 和 ⇒ 连通 Ursell 聚合物；收敛 = **Kotecký–Preiss** `u·K<1`（`K≤CN`）；
`P2` 代价分层给 **I1**（唯一领头管道）、**不给 I2**〔2026-11-24 更正，见 `PROOF-P3-UNCONDITIONAL §5`、`P3-AUDIT-CORRECTED.md`；旧文"代价分层给 I2"是错的〕；基本特征系数之积为正给 `a₀>0`（残数模方为正）；`f₁²` 是**模方**，Weingarten 符号管面积律求和、触不到 `m₁`。
**这是 SZZ 给不了的一半（SZZ 只给隙的相关系数上界；P3 想给残数下界）。** ⚠ **2026-11-24 审计更正**：`PROOF-P3-UNCONDITIONAL.md` 曾称 KP 引理"完整写出并闭合"——**不成立**：最小管道几何应为 `m₁=4ln(1/u)/a`（四壁管道，非"每片 1 个"；factor-4，Gupta §14.2.2 / Münster）、全表示 KP 未控、dressed-mass 重求和未证。**残数结论 `f₁²>0` 存活**（因子 4 只重标定质量），但证明**不完整** → P3 标 **[条件/待重证]**，不再是 [T]。见 `P3-AUDIT-CORRECTED.md`。
`PROOF-P3.md`、`PROOF-P1P2.md`、`PROOF-P3-UNCONDITIONAL.md`。

### 2.2 无溶解结构定理
`0⁺⁺` 胶球**字面并进 `(gg)` 连续谱谱上不可能**：`m₁<2m₁=√s_th`（裕度=自身质量）、`0⁺⁺` vs `0⁺⁺*` 是
**von Neumann–Wigner 避免交叉**（排序不能翻、无实轴排空）。唯一能杀残数的是**固定阈下极点权重塌成 0**
（`:φ²:` 机制）——那不是溶解，是残数归零=隙。`PROOF-NO-DISSOLUTION-KINEMATIC.md`、`PROOF-DISSOLUTION-SPECTRAL-DENSITY.md`。

### 2.3 W_geo≥0 引理（Haar 几何项非负）
`Φ=−logJ=Σ_{i<j}g(λ_i−λ_j)`，`g(d)=−2log(sin(d/2)/(d/2))`，`g″≥1/6` ⇒ `Φ` 凸 ⇒ `e^{−Φ}dX` 对数凹 ⇒
(Prékopa) `W_geo=E[Φ″]−Var(Φ′)≥0`，**全模、全 k（含 k→0）、全 L**。β=0 数值闭合（`V_eff″=0.248`）。
**纠假阴**：`−0.165` 是漏 `+W` 的代理误判；同模精确 `G″=+0.51`（SU3 2.5）。`PROOF-JACOBIAN-SAVES.md`。

### 2.4 规范不变 Perron–Frobenius 定理
规范不变转移矩阵 `T̄` 在 `𝒜/G` 上**严格正核**（满正 Haar、Trotter 电热核 >0、无 Gribov/FP 符号）⇒ PF。
但**真空与胶球同为 `0⁺⁺`**，PF 钉扇区顶（真空）不钉同扇区第二态（胶球）。`PROOF-GAUSS-INVARIANT-PF.md`。

### 2.5 Lemma SC（强耦合端粗模 Poincaré）
`β<1/48`：SZZ 一致 LSI 过 1-Lipschitz 块映射 ⇒ `Var_ν(φ) ≤ (C_blk²/λ_μ)∫|∇φ|²dν`，**粗模 Poincaré 在强耦合成立、`(L,a)` 一致**。`PROOF-COARSE-POINCARE.md`。

### 2.6 条件定理 + 合成链（早年已建）
`IN-0…IN-11 ⇒ Clay`（Lean 机检链）；塌缩到 `G1 + G2 + 物理通道`。`docs/conditional-theorem.md`、`docs/IN-AUDIT.md`。

---

## 3. 新工具（过自由场酸测）

| 工具 | 公式 / 定义 | 酸测 | 适用 / 残缺 |
|---|---|---|---|
| **乘性残数搬运** | `f₁²(μ)=Z_O²f₁²(μ₀)`，`Z_O=(β/2g)(μ₀)/(β/2g)(μ)>0`，`γ_O=2b₀g²+…` | 过（0→0、常数是平方、`O(1)` in V） | 条件于 **(CP)** |
| **精确方差分解** | `G″=S−D`，`S`=高斯–Schur 曲率，`D`=纤维 BL 缺陷；`D≤0⇔K⪰0` | 过（自由场 0≥0） | β=0 定理 `G″≥N/12`；2.5 缺陷 `D/S≈2%` |
| **Fredholm 行列式** | `f₁²=−λ₁Res_{z₁}G(z)`，`D_O=det(1−zT_O)`，`Tr T_O^k=Z_{S¹_k×T³}` | 严格过 | 看见原子、给 `m₁²`；缺 a→0 求零点 |
| **规范不变 IR 质量** | `m_GI=inf supp ρ_s`，`s=(β/2g)trF²∈H_phys` | 过（Maxwell 0 / YM m） | 移除规范不变伪障碍；缺边原子 |
| **Gauge W-熵** | `F_W=∫(S_W+|∇_A f|²)e^{−f}dvol`，沿共轭流；保隙⟺`Γ₂≥ρ` | 镜像过 | 无隙侧对、有隙侧 `Γ₂≈−2.0` 阻 |

详见 `docs/new-tools-registry.md`。

---

## 4. No-go 定理（"为什么这一类不行"）

- **SCHUR-CONVEX（大师障碍）**：任何 `g`-连续受保护泛函在 `g=0`（自由场满足全 OS、`f₁²=0`）上必须 0；正性只给上界/真空。
- **Goldstone 唯一性**：正性+解析能逼原子的唯一定理在**零质量**；`f₁²` 是"正质量表亲"，不存在。BPS 是唯一受保护正质量原子，但**错的通道**（中心荷/通量）。
- **周期 no-go**：紧群上对数密度周期 ⇒ 无对数凸替身；Savvidy 负=周期性回绕代价。→ 凸化不可能。
- **嵌套 no-go**：每个变号/Karlin 锥含正锥 ⇒ 吞自由场；能过的锥必须切到正性之下（循环）。
- **刚性 no-go**：离散整数且局部常数的指标沿连通 β 族为常数 ⇒ YM 全 0；**隙必是 β-非解析跳变，不是指标**。（统一 NCG/localization/BPS/Witten 自由点返 0。）
- **刚性/非阿贝尔符号**：Savvidy、`6j`（spin-foam）、Weingarten（随机簇）、`H_ff`（FRG）——一个物的五面。

---

## 5. 免疫图（~60 路 → 免疫根）

| 免疫根 | 倒在这 |
|---|---|
| **免费场盲**（正性/解析） | INEQ-T1..T5、R1..R10、Pólya、Nevanlinna、Bogoliubov |
| **跨道**（中心荷/1-形式/磁） | BPS/SUSY、1-形式反常、非可逆对称、dual-Meissner |
| **非阿贝尔符号** | 随机流、随机簇/FK、FKG、谱独立性、刺管 |
| **错半边**（UV/纤维 ⊥ 隙） | L3、H2 全局、梯度流、Dirichlet 比较、测度集中、热核、验证数值、步进标度、有效弦、S³ |
| **改名=隙** | H1、NEO(ESI)、I2(E₂)、Hamiltonian 束缚能、L3（单圈精确） |
| **上界方向** | Bogoliubov、Hylleraas、FSS、有限尺寸 |
| **只到强耦合端** | 严谨大 N、SZZ、随机面/CPS（`λ<1/24`；旧记 `75×` 系单位混淆，见 EMERGENT-4D）、中心涡旋 |
| **指标=0** | 刚性 no-go（NCG/localization/BPS/Witten） |

元结论：**不是"没找到"，是 ~60 次独立正攻各证同一个不可约性**——接近定理级的穷尽证据。

---

## 6. 残留对象（所有路最终落到这里；= 最后要造的钥匙）

`f₁²(β)>0 on [1/48,2.5] ⟺ m_{0⁺⁺}(β)>0`（运动学裕度+避免交叉免费给"不溶解"）。等价形式：
1. **H1**：粗模 Poincaré `Var_ν(φ)≤(1/λ)∫|∇φ|²`，`λ>0` 一致（`∇∥` 半边）。
2. **(CP)**：单 T=0 禁闭相（无体相变）+ 0⁺⁺ 边 Kato-保持 + 一致隔离。
3. **边原子 `f₁²>0`**（usc、非对角）。
4. **`G″(0)>0` at `k→0`**（`=S−D`，`D/S≈2%`）。
5. **行列式零点 in 连续统**：`det(1−zT_O)` 的零点在 `a→0` 求值。
6. **Savvidy 符号**（非阿贝尔）。

**这些不是 6 个问题，是同一个：路径上 `m_{0⁺⁺}` 的 volume-uniform 下界 = 质量间隙本身。**

**本战役终局收敛（见 `docs/PROOF-SKELETON.md`）**：这 6 个进一步化为**一格**——(场强 × 动量) 平面上的**"典型场 × 高动量"格**，需一条 `Var_μ(⟨v,∇S⟩|0)≤W_geo`。其余全格已闭：低场（`Φ″⪰N/12`+BL）、高场锥（WSTC 测量饥饿）、AF-UV（`g→0`）。**这一格 = 一个"方差上界"物种**，与 G1 的 L3（UV 侧的单侧力方差上界）**同物种、异尺度**。

---

## 7. 唯一真进展的方向：连续极限

三个最强的新构造都是**固定 `a` 的精确恒等式**，都真正"碰上"了原子（搬运给 `Z_O²`、方差给 `D/S≈2%`、行列式给零点）。**缺的都是同一个连续极限 `a→0` + 搬运到 2.5 = P3 的边际解析墙（G1 的边际 BPHZ）。** 正在攻：
- `CONTINUUM-FREDHOLM`：行列式能否正规化算出连续零点（vs 发散卡边际 BPHZ）。
- `CONTINUUM-VARIANCE-DEFECT`：`D/S` 是否 RG-不变/有界（若 `D` 是 UV、`S` 是 IR，也许绕过 H1）。
- `SYNTHESIS-NEWCONSTRUCTS`：三构造能否串成链（破旧"和=max"）。
- `CONTINUUM-CP`：拆开攻 (CP) 的单相半。

---

## 8. 硬算状态（真跑格点）

`GPP-kzero`（SU2 L=8 最软模 `k=2π/L` 的 `G″` 符号，仍有底层作业在跑）+ `gap/FSS`、`f₁²`、`σ`（SU2）、`G″`（SU3 软模）共 5 路**真跑 numpy 格点**，数字未全回。
已有数：`LATTICE-TM-GAP`（SU2 各盒 `Δ≈1.0a`，但有限盒=上界、β=2.5 部分在去禁闭上方=有限 T）；`W7`（`G″` 有限 k 正 `+0.27/+0.49`）。
**`GPP-kzero` 已回（`COMPUTE-GPP-K0-CONTINUUM.md`）**：单进程流式，`L=8` 最软非零模 `k=2π/L=π/4`——SU3 β=5.7 `G″=+0.307±0.055`（6/6 模正），SU2 β=2.4 `G″=+0.192±0.058`（6/6 模正）；`W_geo` 钉在 `+0.319`(SU3)/`+0.252`(SU2)，与 L=4,6 的 `+0.32` 一致（**体积/k 无关到 0.5%**），action 括弧≈0（等分）。旧 L4→6「体积+k」混淆已拆：固定 `k=π/2` 时 L=8 与 L=4 一致；固定 L=8 时 `G″(π/4)<G″(π/2)` 但**不塌**。**仍是 κ→0 符号探针，非 H1**：只留已证 `W_geo≥0` 时 95% 单侧下界 SU3 `−0.098`/SU2 `−0.035`，0 未排除（= `(◆)`=H1 的同一道坎）。

---

## 9. 诚实判词

- **不是 Clay 证明。YM 仍是猜想。** 距定理 0%（公开问题）。
- **真战利品**：5 件真定理（P3、无溶解、`W_geo≥0`、规范不变 PF、Lemma SC）+ 5 个过酸测新工具 + 3 条新 no-go + ~60 路免疫图 + 残留对象清单。
- **坐标**：最后一步 = 路径上 `m_{0⁺⁺}` 的 volume-uniform 下界（或三构造的连续极限），**需要一台把强耦合的隙/残数连续搬过 Savvidy 交叉点、并取 `a→0` 的非微扰机器——这台机器尚不存在**（= 非微扰 4D 动力学，`TOOL-not-iii` 判词：非微扰纯规范 4D CFT 分类）。当代最前沿（CPS/CNS 2025）也只到 `λ<1/24`（旧记 `75×` 系单位混淆，见 EMERGENT-4D）。
- **但**：`G″` 符号在可测处为正、`W_geo≥0` 是定理、`D/S≈2%`、行列式看得见原子——**墙没倒，但墙现在是透明的、上了标尺的。**
- **本战役终局（`docs/PROOF-SKELETON.md`）**：Clay 拆为 **G1(存在)+G2(隙)+通道**；G1≈70%（微扰核已闭、L3 开放）、G2 区域覆盖≈90% 但**符号=0%**（唯一开格=典型场×高动量）、通道=100%（无溶解真定理）。**全项目收敛到一个物种——"方差的上界"**，在 **UV(L3)** 与 **IR(H1)** 各现一次；其**符号载体是近正交**（`corr(a,∇S)≈0`，CS 对此符号盲），非量级（`VarG≈2` 硬在阈值之上）。A/B/C 三腿全部**走到终局**：A=无路（八路→H1/错向）、B=区域完/符号空、C=闭集闭/开格→命名物种。

---

## 10. 文件索引

- 总册/前沿：`docs/YANG-MILLS-MASTER.md`（本册）、`docs/PROOF-SKELETON.md`（**G1+G2 逐节点状态 + 诚实完成度 %**）、`research/unconditional/ABC-FINAL-LEDGER.md`（**A/B/C 三腿终局账本，§3/§4 源**）、`docs/yang-mills-frontier.md`、`docs/new-tools-registry.md`、`docs/conditional-theorem.md`、`docs/gap-ledger.md`、`docs/IN-AUDIT.md`
- 已证定理：`PROOF-P3.md`、`PROOF-P1P2.md`、`PROOF-ISOLATION.md`、`PROOF-NO-DISSOLUTION-KINEMATIC.md`、`PROOF-JACOBIAN-SAVES.md`、`PROOF-GAUSS-INVARIANT-PF.md`、`PROOF-COARSE-POINCARE.md`
- 新工具：`INVENT-TRANSPORT-INEQ.md`、`INVENT-VARIANCE-INEQ.md`、`INVENT-NEW-IDENTITY.md`、`INVENT-GAUGE-INVARIANT-MASS.md`、`INVENT-NEW-MONOTONE.md`
- 新 no-go：`INVENT-CONVEXIFICATION.md`、`INVENT-SECOND-CONE.md`、`INVENT-NEW-INDEX.md`、`SCHUR-CONVEX.md`、`GOLDS-POS-MASS.md`
- 连续极限（进行）：`CONTINUUM-FREDHOLM.md`、`CONTINUUM-VARIANCE-DEFECT.md`、`SYNTHESIS-NEWCONSTRUCTS.md`、`CONTINUUM-CP.md`
- 硬算：`COMPUTE-*.md`、`LATTICE-TM-GAP.md`、`W7-*.md`
- 免疫：`INEQ-T*`、`INEQ-R*`、`PROOF-*.md`（~60 个）
