# TRANSPORT-RIGOROUS-NUMERICS — 严格随机量子化 / 严格区律 / 严格数值 三条 Douglas 路线的诚实裁决

> 2026-11-24 · 任务：评估 Douglas 2026 (Nature Rev. Phys.) 点名的三条当前活跃严格路线，是否对**交叉点墙** (β≈2.5, α_s~O(1)) 提供 `(L,a)`-一致的 **G2 隙**新攻击。逐条查一手源（arXiv/DOI）。**一律标注** [T] 真定理 / [C] 条件 / [O] 精确缺失步 / [NG] 具名类 no-go / [V-here] 本轮新算 / [R] 外部据报（已核）。**不在未证上敲 QED。**

---

## 0. VERDICT — 一句话

**这条线不是新的交叉点攻击，它就是同一堵墙换语言（rigorous-SPDE / rigorous-RG / validated-numerics 版）。三路全部只到 G1 或强耦/大 N/低维/小体积；无一条在质量生成尺度 α_s~O(1) 给出小控制参数。** 具体：
- **随机量子化 = G1 技术**（构造 + OS）。4D 连 G1 都开（无构造）；它在 G2 层面的**唯一**输出是强耦 LSI（SZZ `β<1/48`）——**恰好就是本项目的 P3**。Langevin 生成元的谱隙 = 虚时间松弛率 ≠ 物理质量隙（H1-T1 已钉）。桥 `LSI⇒clustering` 需 Bakry–Émery ⇒ 在 Savvidy `β=2.5` 死。
- **严格区律 = 禁闭/强耦/大 N**，非 0⁺⁺ 隙。CNS 2505.16585 是**格点 U(N)、`β≤β₀(d)=10^{-10d}/d`（均匀于 N）**的强耦结果；物理 SU(3) 对应 `β≈0.95`，与证明区差 **~10⁴¹ 倍**（[V-here]）。't Hooft regime (2509.04688 / Nissim 2510.22788) 的"质量隙"**依赖 DF80 的 σ-模型聚类作为独立输入**。且逻辑方向是 **隙 ⇒ 区律**（Chatterjee / DF80），区律**下游**，证区律不给隙。
- **严格数值 = 无 `(L,a)`-一致隙**。有限体积 `β=5.7` 的证书（A5）用现有 validated-numerics 不可行，且**有限 k 酸测就败**（自由无质量场同样 `G″(k>0)>0`）。

**G1/G2 分割如今异常锐利**：G1 已有一份（有争议但已过 Gross 委员会 Lean 审查的）构造 [R]；它**明确不含质量隙、唯一性、区律**。故 G2 = 整个剩余问题，而三路**无一触及**。**Clay(G2) 仍 0%。** `∀class ≠ ∀possible`：这是具名三族的裁决，不是全局不可能定理。

---

## 1. 逐条裁决表

| # | 方法 | 已证 (标签) | regime | G1/G2 | 跨 crossover? | 酸测（自由无质量场） |
|---|---|---|---|---|---|---|
| 1 | **Hairer–Mattingly 强 Feller** (AIHP 2018) | 奇异 SPDE 半群强 Feller ⇒ 不变测度唯一/遍历 [R] | fixed-cutoff, KPZ/Φ⁴₃ | 既非 G1 构造亦非 G2；定性 | **否** | **败**：强 Feller 对无质量场**也成立**（不区分）⇒ 作隙工具 **vacuous** |
| 2 | **Tsatsoulis–Weber Φ⁴₂ 谱隙** (AIHP 2018) | Φ⁴₂ Langevin 指数遍历/Poincaré [R] | 2D 超可重正, 固定体积 | G1+（2D 隙另证） | **否**（2D 特殊） | **过**（∞ 体积无质量 `λ=0`），但仅 2D |
| 3 | **CCHS 3D YM–Higgs** (Invent. Math. 2024) | T³ 上规范协变 Markov 过程（局部时间，可爆） [R] | 3D, T³ | **G1**（构造） | **否**（非 4D） | N/A（无隙声明） |
| 4 | **4D 随机量子化** | **开**：无构造（[g]=0 临界、无规范协变 BPHZ） [R] | 4D | G1 开；G2 未触 | **否** | — |
| 5 | **SZZ 强耦** (CMP 400, 2023) = **本项目 P3** | LSI/Poincaré, `f₁²≥r₀>0`, ∞体积唯一+隙 [R] | `β<1/48`(=1/(16(d-1))), 任意 d | **G2**（强耦） | **否** | **过**（自由无质量 `f₁²=0`），但强耦窗口 |
| 6 | **CNS 区律** (arXiv:2505.16585) | 格点 U(N) 区律, `d≥2`, `β≤β₀(d)=10^{-10d}/d` 均匀于 N [R] | **强耦**（β 极小） | 非 G1 非 G2（禁闭） | **否**（差 ~10⁴¹×） | **过**（自由场 perimeter 律、无区律），但**错通道**+下游 |
| 7 | **CNS/Nissim 't Hooft** (2509.04688; 2510.22788) | 't Hooft regime 区律；U(N) 质量隙+∞体积唯一，**经检查 DF80 质量隙条件** [R] | `β~cN`（大 N） | G2 in 't Hooft | **否**（大 N；DF80 独立输入） | **部分**：条件于 DF80 的 σ-模型聚类（独立硬输入） |
| 8 | **Chatterjee 区律 / 强隙⇒禁闭** (PTRF 2019; CMP 2019,2021) | 强耦区律；"强质量隙 ⇒ 禁闭"（条件）[R] | 强耦 / 条件 | 隙⇒区律（**下游**） | **否** | **条件**（强隙假说未证于弱耦） |
| 9 | **Chatterjee SU(2) YM–Higgs 标度极限** (PMP 2026) | >2D 首个非阿贝尔标度极限；**经 Higgs 机制**得质量 ⇒ 大质量**高斯**场 [R] | d>2, 有 Higgs | 质量但**外置** | **否** | **败**：硬质量（Higgs VEV），非纯 YM 动力学生成 |
| 10 | **严格有限体积数值** (本项目 A5) | 有限 k `G″>0` 的条件证书；`β=5.7` 不可证 [V-here/O] | 固定 (β,L,k>0) | G2 尝试 | **否** | **败**：有限 k 自由场 `G″=4sin²(k/2)>0` 同号 |
| 11 | **OpenAI/Lean G1 构造** (2026-11-10) | 4D Wilson 格点 → Wightman QFT 标度极限，任意单连通紧群 [R] | **连续 4D** | **仅 G1** | G1 only | 明确**排除**隙/唯一性/区律 |

> 表 6 的 `β≤β₀(d)` 数值与"`λ<1/24`"互斥——见 §5 订正（本项目 AREALAW-GAP 的 1/24 记法在源文中不存在，SOTA-VERIFIED 已正确标记为 unverified）。

---

## 2. 决定性分析：为什么是同一堵墙换语言

### 2.1 随机量子化 = G1 技术，G2 面 = 虚时间 ≠ 物理时间
- **构造层**：CCHS 在 2D (Publ. IHÉS 2022) 与 3D YM–Higgs (Invent. Math. **237**, 2024) 用 regularity structures 造出规范协变 Markov 过程 [R]。**4D 是临界**（`[g]=0`）：Hairer 的次临界度条件失败、无规范协变 BPHZ、无良定性 [R]（H1-T1 §1）。故**连 G1 在 4D 都未由该法给出**；且 4D Φ⁴ 是唯一解出的临界情形而它**平凡**(Aizenman–Duminil-Copin)，故即便良定性也不给非平凡测度。
- **隙层**：`L=Δ−∇S·∇` 在 `L²(e^{−S})` 自伴，其谱隙 `λ₁` 是**虚时间 L² 松弛率**，**不是**转移矩阵隙 `e^{−aH}`（H1-T1 §2，NA-1 v2 的 object-mismatch）。致命证据：NA-3 实测模拟时间混合 `τ_P≈O(1)`（有限体积 `λ₁=O(1)`）而质量隙仍开——**解耦**。唯一桥 `LSI(μ)⇒Euclid clustering`（Zegarlinski）需 Bakry–Émery，而 `Γ₂=Ric+Hess S` 在 `β=2.5` 的 Savvidy 方向 `≈−1.1β`（SU(2) L=4：min −2.81，~22% 负）⇒ **`Γ₂≥ρ>0` 假**（W6；H1-T1 §3）。
- **唯一真实隙输出**：SZZ (CMP 400, 2023) `β<1/48` 的 LSI + 指数遍历 + Wilson-loop clustering ⇒ ∞体积正质量隙 [R]——**这正是本项目 P3**（PROOF-P3-UNCONDITIONAL）。**所谓"正交的随机量子化线"，在隙层面给出的东西，就是我们已经拥有的 P3。**
- **酸测**：强 Feller 对**无质量**场也成立（如有限体积自由场），作"隙工具"**vacuous**；Poincaré/LSI 版本在 ∞体积自由无质量场 `λ=0` ⇒ 技术上**过**酸测，但它的 `λ` 一均匀化就等于 `H1`（隙本身），**循环**。

### 2.2 区律 = 强耦/大 N + 错通道 + 逻辑下游
- **regime（[R] 已核一手）**：CNS 2505.16585 是格点 U(N) on `[-L,L]^d⊂ℤ^d`, `d≥2`, 阈值 `β≤β₀(d)`, 显式 `β₀(d)=10^{-10d}/d`，**均匀于 N**，改进 Osterwalder–Seiler 1978 的 `β<c/N`；**只证区律，不证隙** [R]。**不**是 't Hooft regime，也**不是**连续统。
- **regime 分离（[V-here]）**：论文密度 `∏_p exp(2β ReTr U_p)dU`（未归一化迹）↔ 标准 Wilson `exp((β_W/N)ReTr)`，故 `β=β_W/(2N)`。SU(3) 物理 `β_W≈5.7` ⇒ `β≈0.95`；而 `β₀(4)=2.5×10^{-41}`。**物理值 / 证明阈值 ≈ 3.8×10⁴⁰**——差 **~41 个数量级**。（即便 SZZ 的 `1/48≈0.021` 也在物理 `0.95` 之下 ~46×。）**证明区与质量生成尺度之间有 40+ 个数量级的鸿沟：这正是"无小参数"墙的定量化身。**
- **'t Hooft 侧**：2509.04688 用 SZZ 动力学方法在 't Hooft regime 得区律，**经检查 DF80 (Dürhuus–Fröhlich, CMP 75 (1980) 103) 的质量隙条件**后"区律直接推出" [R]。DF80 的内容是：**(ν−1) 维 σ-模型的指数聚类 ⇒ ν 维 YM 区律**；聚类只在 1D σ-模型、2D U(n)×U(n)、或大耦合 `g²≳O(ν)` 证得 [R]。即：**"隙"是 σ-模型的、且作为独立输入**，不是纯 YM 自生成，且只在强耦/低维证。
- **错通道 + 逻辑下游**：区律的转移矩阵内容是**弦/流管 1-形式**道（能量 `V(R)≈σR`, `m_string=√σ`）；`G″(k→0)` 是**中性 0⁺⁺ 0-形式** plaquette 关联。`m_{0⁺⁺}=c√σ` 是**实测透射（c≈3.5–4），不是定理**（AREALAW-GAP §channel）。且 Chatterjee："强质量隙 + 中心对称 ⇒ 区律"（CMP 385, 2021）[R]——**方向是 隙⇒区律**，区律**下游**，证它**不给**隙。
- **酸测**：区律对自由无质量场**失败**（自由场 perimeter 律、无区律）⇒ 作为**判别量**过酸测（非 vacuous），但它**不 lower-bound `G″`**，故**不是隙工具**。

### 2.3 严格数值：无 `(L,a)`-一致，且有限 k 酸测就败
- **不可能性（本项目 A5，[V-here] 已算）**：`β=5.7, L∈{4,6}` 的精确有限区证书架构存在、且用已证的 `E[Φ″|0]≥N/12=¼` 给出**正的逐点余量** `+0.092/+0.044`；但证书**不是定理**，三面独立墙：(i) MC 的 block-bootstrap σ 是**统计**误差、不是**严格区间**（无证的 sampler/mixing ⇒ 正是隙）；(ii) 唯一**先验**上界 Cauchy–Schwarz `|X|≤1.06` 是整个预算 `0.044` 的 **24×**；(iii) 严格计算是 `G^{4L⁴}`（L=4 = **8192 实维**）上的非光滑（branch cut + eigenvalue collision）、符号不定（Savvidy）比值积分——超出现有 interval/TT/转移矩阵技术（A5 §3–4）。
- **即使闭合也非隙结果**：在固定 `k=2π/L`，**自由无质量场已有 `G″=4sin²(k/2)=2(L4)/1(L6)>0`**，同号且更大 ⇒ 有限-k 正性**不是** YM-vs-free 判别量，**不触及** `k→0`（自由场 `G″(0)=0`）处的酸测对象（A5 §5）。
- **诚实 SOTA**：经同行评审的**严格 (L,a)-一致 4D 格点隙不存在**；严格的格点隙只在强耦（cluster/KP 展开）。严格的连续统需 `(L,a)` 一致的双端界——UV = 本项目的 L3，IR = H1，皆开。2026 年若干 Zenodo/GitHub/hexagonmath 自出版 "interval-arithmetic + Lean 解 Clay" PDF **非同行评审、未核实**，**不传播为结果**（宪章 §2）。

### 2.4 结论：为什么打不破墙
三路的控制参数**只在** 强耦 (`β≪O(1)`)、大 N、低维 (2D/3D)、或小体积 处小——**没有一条**在质量生成尺度（`α_s~O(1)`，`ρ~1/Λ`）有小参数。这与 BRIEF §3 的"构造性 RG 到交叉区停"、"大 N/对偶只到强耦/给上界"完全一致。且本轮外部事实把结构照得更亮：**G1 已有一份（有争议的）Lean 构造，明确不含隙**；故剩余问题是**纯 G2**，而三路**无一**在 G2 的交叉点处说话。**同一堵墙，SPDE/RG/validated-numerics 三种新方言。**

---

## 3. 每条的精确机制（为何死不在于名字，在于缺失对象）

| 线 | 死亡机制 | 落点 |
|---|---|---|
| 随机量子化 | 4D 临界 ⇒ G1 开；G2 面 = 虚时间隙 ≠ 物理隙；桥需 BE ⇒ Savvidy 负曲率 | `H1`（Bakry–Émery 面）+ 4D 构造 |
| 区律 (CNS/Chatterjee) | 强耦/大 N；1-形式弦道；DF80 σ-模型聚类独立输入；隙⇒区律下游 | `H1`/`I2` 的**产物的下界**，非隙本身 |
| 严格数值 (A5) | 8192 维非光滑符号不定积分；先验界 24× 预算；有限-k 酸测败 | `H1`（方差上界 (◆)） |

三路在正确坐标里都定义**同一缺失物种**——一个**单侧方差上界**，在 UV (L3) 与 IR (H1) 各现一次（PROOF-SKELETON §1）。**本轮没有改变这一点，只是用 Douglas 的语言重新确认了一遍。**

---

## 4. 唯一值得追的子方向 + 精确缺失步（诚实：它仍落 H1）

**先给否定的诚实**：三路**均非**新攻击；随机量子化与耦合/Doeblin（D2-COUPLING-DOEBLIN 已把耦合物种磨到 `λ=−b′(0)=G″` 恒等，von Renesse–Sturm 锐 ⇒ 零松弛）**都已耗尽**到 H1。故**没有一条能破墙的新物种**。

**若必须点名一个"探索最少、且非本项目已试技术"的子方向**：**严格的张量网络 RG（TRG/HOTRG 型）带认证截断误差，做有限体积隙证书 + 用精确乘性残数恒等式 `f₁²=Z_O²f₁²` 搬运**。理由：(a) 这正是 **Douglas 明确建议**物理学家发展的"更强大的计算 RG（如 tensor RG）"[R]，是活跃领域方向而非我们臆造；(b) 本项目**未严格试过**它——A5 试的是 **interval 积分**（8192 维），PROOF-COMPUTER-ASSISTED 试的是**认证转移矩阵/Peter–Weyl Feshbach**，二者与 TRG 是**不同技术物种**；(c) 它是一个真正的离散数值物种。

**精确缺失步 [O]**：**对截断误差的、键维-uniform 的严格界，使 TRG 映射在交叉点 `β≈2.5` 处有认证的压缩常数**（等价的"discarded-weight 可控的 RG 映射收缩"）。这正是非微扰控制 = `H1` 在 numerical-RG 语言里的重述。**酸测**：自由无质量场的 TRG 是不动点（临界、不收缩），故"TRG 隙证书"对自由无质量场**失败** ⇒ 原则上过酸测；但**认证其收缩常数跨交叉点均匀**，与现有一切路线同难。**它买到的是"用第三种技术再测一次 H1"，不是削弱它。**

**一个具体的、项目未审计的可证伪子目标 [O]**：**在固定 N=3、`β≈2.5` 检查 DF80 的质量隙条件**（即 3D 非线性 σ-模型的指数聚类）。这是 CNS/Nissim 't Hooft 结果的**唯一外部硬输入**，项目尚未单独审计其内容。**精确缺失步**：DF80 条件本身是一个定量的 3D σ-模型质量隙下界；在 `β≈2.5` 固定 N 证它 = 一个已知-open 的独立硬问题（σ-模型聚类现仅大耦合/低维证）——**仍非新钥匙，但把"区律线的隙从哪来"钉死到一个可命名、可单独攻的对象上**。

---

## 5. 诚实订正（宪章 §2：引用逐个核对）

1. **CNS 阈值**：AREALAW-GAP 记的 `λ<λ₀=1/(8(d−1))=1/24` **在源文中不存在**。2505.16585 的真实阈值是 `β≤β₀(d)=10^{-10d}/d`（均匀于 N，未归一化迹）[R]。SOTA-VERIFIED §3 已正确标 `1/24` 为 unverified；**本文件据原文订正**。
2. **"Chandra–Chevyrev–Hairer–Shen, *…in 3D via random geometries*"**：查无此确切标题。已核实的 CCHS 3D 文为 **"Stochastic quantisation of Yang–Mills–Higgs in 3D"**, Invent. Math. **237**(2), 541–696 (2024), arXiv:2201.03487 [R]。若他处确存在"random geometries"篇，未检索到——**不宜据错题名引用**。
3. **"Chatterjee, *Rigorous QFT functional integrals over the p-adics*"**：**归属错误**。该 p-adic 系列实为 **Abdesselam–Chandra–Guadagni** (arXiv:1210.7717; 1302.5971)，且**不涉 YM/质量隙**（p-adic BMS 标量）[R]。Chatterjee 的严格 YM 工作实为 强耦区律 (PTRF **175**, 2019; CMP **370**, 2019)、"probabilistic mechanism for quark confinement" (CMP **385**, 2021)、"Yang–Mills for probabilists" (arXiv:1803.01950)、及 SU(2) YM–Higgs 标度极限 (arXiv:2401.10507; PMP **7**, 2026) [R]。**不传播错误归属。**
4. **G1 状态更新**：2026-11-10，Gross 委员会（Douglas, Jaffe, Klebanov, Nekrasov）审查通过 OpenAI 的 Lean 4 构造：4D Wilson 格点 → Wightman QFT 标度极限，任意单连通紧群；**明确不含** 隙/唯一性/区律，**非 Clay 裁决、未期刊评审** [R]。（旁支：Gowers/Hairer/Witten 顾问组曾要求立即开源。）**不影响本裁决：G1 不是 G2 的钥匙。**

---

## 6. 外部引用（本轮逐一核实）

- **[R]** Douglas, M. R., "The Yang–Mills Millennium problem," *Nature Rev. Phys.* **8**, 86–97 (2026), DOI 10.1038/s42254-025-00909-2（12 Jan 2026；审稿含 Ünsal, Chatterjee）。点名：rigorous stochastic quantization、rigorous strong-coupling expansion（概率论）、computational/tensor RG；"no clear path known"。
- **[R]** Hairer–Mattingly, "The strong Feller property for singular stochastic PDEs," *Ann. IHP B* **54**(3), 1314–1340 (2018), DOI 10.1214/17-AIHP840, arXiv:1610.03415（KPZ, Φ⁴₃）。
- **[R]** Tsatsoulis–Weber, "Spectral gap for the stochastic quantization equation on the 2-dimensional torus," *Ann. IHP B* **54**(3), 1204–1249 (2018), arXiv:1609.08447（Φ⁴₂ 指数遍历）。
- **[R]** Chandra–Chevyrev–Hairer–Shen, "Stochastic quantisation of Yang–Mills–Higgs in 3D," *Invent. Math.* **237**(2), 541–696 (2024), DOI 10.1007/s00222-024-01264-2, arXiv:2201.03487。
- **[R]** Chandra–Chevyrev–Hairer–Shen, "Langevin dynamic for the 2D Yang–Mills measure," *Publ. Math. IHÉS* (2022), arXiv:2006.04987。
- **[R]** Shen–Zhu–Zhu（= P3 源）, "A stochastic analysis approach to lattice Yang–Mills at strong coupling," *CMP* **400**(2), 805–851 (2023), arXiv:2204.12737（`β<1/48`, d=4）。
- **[R]** Cao–Nissim–Sheffield, "Expanded regimes of area law for lattice Yang–Mills theories," arXiv:2505.16585（`d≥2`, 格点 U(N), `β≤β₀(d)=10^{-10d}/d` 均匀于 N；只证区律）。Cao–Nissim, "Dynamical approach to area law," arXiv:2509.04688（'t Hooft regime；经 DF80 条件）。Nissim, arXiv:2510.22788（U(N) 't Hooft regime 质量隙）。
- **[R]** Osterwalder–Seiler, *Ann. Phys.* **110**(2), 440–471 (1978), DOI 10.1016/0003-4916(78)90039-8（强耦 `β<c/N` 区律）。
- **[R]** Dürhuus–Fröhlich, "A connection between ν-dim Yang–Mills and (ν−1)-dim non-linear σ-models," *CMP* **75**(2), 103–151 (1980), DOI 10.1007/BF01222514（σ-模型聚类 ⇒ 区律）。
- **[R]** Chatterjee, PTRF **175** (2019) 强耦区律; CMP **385**, 1007–1039 (2021) probabilistic confinement; arXiv:1803.01950 (Yang–Mills for probabilists); arXiv:2401.10507 / PMP **7**, 339–381 (2026)（SU(2) YM–Higgs 标度极限 → 大质量高斯）。
- **[R]** Abdesselam–Chandra–Guadagni, arXiv:1210.7717 / 1302.5971（p-adic；**非** Chatterjee，**非** YM）。
- **[R]** OpenAI/Lean 4 G1（2026-11-10；Gross/Jaffe/Klebanov/Nekrasov/Douglas 审查；明确排除隙）；repo `github.com/yang-mills-group/YM`。
- **内部**：PROOF-SKELETON、BRIEF-FOR-MODEL、H1-T1-stochastic-quant、D2-COUPLING-DOEBLIN、SOTA-VERIFIED、A5-RIGOROUS-NUMERICS、AREALAW-GAP、G2E-RIGOROUS-NU、PROOF-P3-UNCONDITIONAL、PROOF-COMPUTER-ASSISTED、W6-BE-at-2.5、TRANSPORT-MACHINE、TRANSPORT-FIVE-WAY-VERDICT。
- **[V-here]** `β=β_W/(2N)`；SU(3) `β_W=5.7→β=0.95`，`β₀(4)=2.5×10^{-41}`，比值 `3.8×10⁴⁰`；SZZ `1/48≈0.021` 相对物理 `0.95` 差 ~46×。

**Sharpest line.** 三路名不同、方言不同，坐标相同：随机量子化把 G2 约到**虚时间隙 ≠ 物理隙**、并在 4D 连 G1 都开；区律把控制参数压在 `β≤10^{-40}`（差物理 41 个量级）或大 N、且"隙⇒区律"方向使区律**下游**；严格数值在 `β=5.7` 不可证、且有限-k 酸测就与自由场同号。**G1 已被一份（有争议的）Lean 构造拿走，剩余的纯 G2 三路皆不触。这叫同一堵墙，不叫新攻击。** 唯一可追的 TRG 或 DF80-条件审计，精确缺失步仍=`H1` 的又一次重述。**Clay(G2) 0%，不在未证上敲 QED。**
