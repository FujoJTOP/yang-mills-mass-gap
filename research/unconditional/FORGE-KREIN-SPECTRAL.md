# FORGE-KREIN-SPECTRAL — definitizable/Pontryagin spectral theory on the indefinite Savvidy Hessian: the negative index is a *count* and its floor is the *UV* scale `g|B|`; the gap lives at the `k→0` point where the definitizing split degenerates

> 2026-11-23 · FORGE（Krein/Pontryagin 谱论）. 判词 (verdict)：**不 5/5。** Krein 空间谱论确是 `H=−D²+2adF_B` 这个不定算子的**精确、正确**语言——它真正非 Hilbert，且 definitizable 算子的谱函数确实存在——但它只交付两样：(a) **负指标** `r_−`，是一个**计数**（`r_−/V=c_N(gB)²/8π²` 广延，**不是**一致的 Pontryagin `κ<∞`）；(b) 一个 **J-正下界 `=Θ(g|B)`**，即 **Landau 能级/UV 尺度**，**不是**规范不变 `m₁`。`计数≠底`、`UV≠IR`。到 `inf supp ν>0` 的桥仍是同一条 `(H)`/H1。**#3 PASS（其主场）；#1 部分；#2,#4,#5 FAIL。**

## 1 对象——为什么不定 Hessian 天生是 Krein 对象
在常阿贝尔（Abelian）`F₁₂=B` 上，横模 Hessian `H=−D²+2adF_B` 的谱是 Landau 阶梯 `λ_{n,s}=k_∥²+(2n+1)g|B|+2s\,gB`，最低阶 `λ_{0,−}=k_∥²−g|B|`，底 `−g|B|`（`k_∥=0`），负模数 `N_neg/V=c_N(gB)²/8π²`（`c_{SU2}=1,c_{SU3}=3`），深度 `Θ(g|B)`，正阶在 `k=0` 处也落回 `Θ(g|B)`。`H` 对普通的 **Hilbert** 内积仍自伴、有下界（`H⪰−2‖F‖_∞`）——所以它的“不定”是**值的符号**，不是内积。自然地，`J:=sign H` 让 `H` 在 `⟨·,J·⟩` 上正定，得到一个 Krein 空间，负指标 `κ=r_−`。Krein 谱论（Pontryagin Π_κ、Langer definitizable、谱函数）正是处理“自伴但不定”的完备框架——这是它对其它正性工具（全部要 Hilbert 正性、全部只给上界）的**唯一结构性优势**。

## 2 锻造的工具：definitizable + 谱函数 + 负指标
若 `H` 在 Π_{r_−} 中自伴且有非空预解集，`r_−<∞` ⇒ **definitizable**（Langer）：存在多项式 `p` 使 `p(H)⪰0`，并有谱函数 `E(Δ)` 落在 `ℝ` 上、仅在有限“临界点”退化；`J`-正交分解把空间拆成极大正子空间 `H_+` 与有限维负子空间 `H_−`（`dim=r_−`）。**锻造目标**：`inf spec(H|_{H_+})≥κ_{IR}>0` 且 `(L,a)`-一致 ⇒ 直接是间隙。负指标 `r_−` 控制不定部分的“秩”，`definitizing` 多项式 `p` 把不定性吸收成正性。这是 SPEC 想要的“对凹/不定区仍给下界”的那台机器。

## 3 它到底交付什么——两处盲点
**(a) 地板是 `Θ(g|B)`，不是 `m₁`。** `H_+` 的底 `=λ_{0,+}(k=0)=Θ(g|B)`（与 tachyon 同尺度——Landau 能级劈裂）。这是**微扰胶子磁质量/场强尺度**（质量²，`g|B|∼1/a²`，被钉在 cutoff），**不**经量纲透射成 IR `m₁²`：`m₁²/(g|B|)=O(a²)→0`。它是**背景依赖、规范非不变量**的量，不出现在胶球谱里。
**(b) 负指标是计数，间隙是极小。** Pontryagin 全套（负指标、符号差、谱移动）是**数值域/计数**理论；`m₁=inf(supp ν)` 是**底部/尺度**理论。正如 D-BBS-LOWRANK-KATO §a：“秩小是计数命题，而间隙是 min——一个负本征值即决定”。Krein 计数对 min 结构性地盲。**临界/软化点 `k_∥=0` 正是 Savvidy 带底部与粗模重合处**（D-BBS §c：`PᵀHP` min `≈−0.3`）；definitizing 劈裂恰在此**失去一致裕度**——谱函数给出的是 `g|B|` 级 J-正底，而物理（Hilbert）正底被自由零模 `h_+(k→0)→0` 拽到 0。Krein 把“UV 塔锥空”精确证成形（`λ_j>k_⋆` 壳上 `κ=0`、`Tr H_−≲∫|F|³`，D-BBS 的 UV 半定理），但间隙坐落在**唯一非空的那层 `k→0`**，正是临界点。

## 4 分叉 no-go：不存在一致有限 `κ`；且过不了酸测
要让 Krein 地板变成物理间隙，需在**物理侧**有 `inf spec(H)≥κ_{IR}>0`（`(L,a)`-一致）。两条路互斥封死：
- **取 `J=sign H`（使地板物理）⇒ `κ=dim H_+=∞`**，且 `sign H` 需已知正谱——**循环**。Π_κ 的 `κ<∞`（Pontryagin）在热力学极限破裂：`r_−∝∫|F|²∝V→∞`。
- **取先验 Krein 度量（BRST/Gupta–Bleuler 不定度规）⇒ `κ=O(1)/模<∞`，但与 β、与相互作用无关**——自由无质量场有**逐字相同**的不定度规结构（timelike 胶子 + ghost）。故任何由此读出的谱结论**自由/相互作用同结论**，酸测**空洞通过或恒无隙**。
- **自由无质量场 + 外部背景**：`F≠0` 时 `g|B|>0` 的“地板”照常出现，理论**仍无隙**（自由！）⇒ 这是 §3(a) 的 `g|B|` 对错判据的正解——它**不**区分有隙/无隙。至于“signed amplitude = 负指标符号差”：它是**整数**，沿连通 β 族（`β=∞` 达自由高斯）局部常值 ⇒ 必等于自由值 `0`（**B3 刚性 no-go**，INVENT-NEW-INDEX）。整数逃逸也被封。

## 5 五点记分（对 SPEC i–v）
1. **非微扰** — **部分**：算子是精确算子、谱论无反 `1/N,g` 展开；但对象 `−D²+2adF_B` 是二阶变分/单圈曲率（半经典）。
2. **support 下界 或 signed amplitude** — **FAIL**：给的是 `g|B|`（UV）与计数 `r_−`，非 `inf supp ν=m₁`；signed 整数死在 B3。
3. **活过 Savvidy 不定 Hessian** — **PASS**（这是它唯一的真主场：definitizable 谱函数对不定区成立）。
4. **自由无质量场失效** — **FAIL/空洞**：BRST 度规自由场逐字相同；`g|B|` 对有背景的无隙自由场照给“地板”。
5. **`(L,a)` 一致** — **FAIL**：`r_−∝V→∞`，地板 `~1/a²`（UV 钉死）。

**计分：仅 #3 干净 PASS、#1 部分，~2/5。非 5/5。**

## 6 精确卡点
```
Krein: H=−D²+2adF_B 在 Π_{r_−} 中 definitizable, 谱函数 E 存在     [真]
J-正底 = λ_{0,+}(k=0) = Θ(g|B|)  = Landau/UV 尺度  ≠ m₁           [真, 非隙]
r_−/V = c_N(gB)²/8π² → r_− ∝ V（广延）⇒ 无一致有限 Pontryagin κ   [真]
先验 BRST κ=O(1) 但自由场逐字相同 ⇒ 酸测空洞                      [真]
J=sign H 才有物理底 ⇒ κ=∞ 且循环（需正谱）                        [循环]
signed 负指标 = 整数 ⇒ B3 刚性 ⇒ 恒等自由值 0                      [no-go]
```
**精确缺失的一步**：`inf spec(H|_{物理侧}) ≥ κ_{IR}>0`，`κ_{IR}` 为 IR 尺度、`(L,a)`-一致——但物理侧要求 `J=sign H`（循环），先验 `J` 自由场盲；且 `r_−` 广延使任何一致 Pontryagin `κ<∞` 不存在。这一缺失不等式**逐字就是** BLIND-LOWER §6 的 `(H)`，等价 H1。

**Sharpest line.** Krein/Pontryagin 谱论是这台不定算子的**正确且不完备**的钥匙：它精确证明 UV 塔锥空、把负指标与 `Tr H_−≲∫|F|³` 严格界住（D-BBS 的 UV 半定理，记成谱函数语言）——**真收获**。但它的一切不变量是**计数/符号**（负指标、符号差）或**UV 尺度**（`g|B|` Landau 底），而间隙是 `k→0` 底部的 **min**、是**经量纲透射的 IR 尺度**；definitizing 劈裂恰在 `k=0`（Savvidy 带底）失去一致裕度，而那里正是粗模与最负模重合处、正谱落回自由零模——**Krein 的临界点就是间隙要印数的那一点**。取 `J=sign H` 循环（`κ=∞`+需正谱），取先验 BRST 度规自由场盲（酸测空洞），signed 整数死在 B3。**故 Krein 空间给的是 H1 的第五张脸（计数/谱函数脸），不是它的削弱；#3 真、其余是 Clay 内容。5/5 未达，Clay 仍 0%。**

*Refs.* Pontryagin, *Izv. Akad. Nauk SSSR Ser. Mat.* **8** (1944) 243; Krein 1950; Iokhvidov–Krein, *Spectral Theory of Operators in Spaces with an Indefinite Metric* (AMS 1982); Langer, Habilitationsschrift Dresden 1965, and “Spectral functions of definitizable operators in Krein spaces,” *Lecture Notes in Math.* **948** (1982) 1–46; Jonas, *Beiträge Anal.* **16** (1981) 121; Čurgus, *Integral Equations Operator Theory* **8** (1985) 462; Azizov–Iokhvidov, *Linear Operators in Spaces with an Indefinite Metric* (Wiley 1989); Jonas–Langer, *J. Operator Theory* **2** (1979) 63. Internal: PROOF-SAVVIDY-EXACT-SIGN (阶梯、`λ_min=−g|B|`、`N_neg/V=c_N(gB)²/8π²`、`Tr H_−≤8C_N∫|F|³`); D-BBS-LOWRANK-KATO (UV 塔锥空、`k→0` 软底、`PᵀHP` min `−0.3`、计数≠min); TRANSPORT-MACHINE (Stage 3=H1/I2、三面 no-go); BLIND-LOWER §6 (`(H)`、唯一下界引擎); INVENT-NEW-INDEX (B3 刚性 no-go); F2-KREIN-RUTMAN-QUANT / G2G-TRANSFER-SPECTRAL (Krein–Rutman=锥，非 Krein 空间；rate=G″=H1); PROOF-EDGE-E2 (bottom-functional 盲). 无虚构界。
