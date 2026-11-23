# 连续统 Wilson 环面积律（禁闭）：定理、引理链、与质量隙的关系、障碍

> 2026-11-22 · 范围：4D 纯 Yang–Mills，`G` 为紧致连通单 Lie 群（Clay 要求任意紧单群；`SU(N)` 为特例）。存在性已条件化于 16 输入、仅子列极限；质量隙（MG）等价于体积一致的指数聚类速率。
> 本报告不重复主 session 的弱→强耦合交接控制量分析（见 `research/physics_handoff.md`），而是处理**禁闭/面积律**这一独立子问题。
> 结论先行：**连续统面积律是公开问题，与 MG 在逻辑上独立、在技术上共享同一个 RG 墙**；可证的核心是格点强耦合面积律及其精确曲面表示。

---

## 0. 结论摘要（≤300 字）

连续统面积律 = 弦张量 `σ > 0`：对光滑闭曲线 C，`|⟨W(C)⟩| ≤ C₁ e^{−σ·Area(C)+c·Perim(C)}`。它**不在** Clay 陈述内，也**不蕴含、不被蕴含于** MG：
- 面积律管的是**带基本电荷的弦扇区**（Wilson 环 = 外源夸克探针的响应），其转移矩阵基态能量 `E₀(L) ≥ σL`（大 L 线性增长）；
- MG 管的是**中性胶球扇区**（局域规范不变算子的第一激发），Clay 只要求后者 `Δ>0`。

二者是同一 `E₀(L)` / 转移矩阵谱的不同端，**逻辑独立**；完整解需两者（"胶球有隙 + 弦有张力"）。

已证：格点强耦合面积律（Osterwalder–Seiler 1978、Seiler 1982），随机曲面表示（Cao–Park–Sheffield 2023），扩展 β 区间（Cao–Nissim–Sheffield 2025），3D U(1) 单极凝聚面积律（Polyakov/Göpfert–Mack），4D Z₂ 自对偶面积律（Wegner）。
未证：**4D 非阿贝尔连续统面积律**——卡在 `β→∞` 的 UV→IR RG 交接，与 MG 同一堵墙。弦张量的连续统值必为 `σ = c_σ Λ²`（量纲 transmutation），`c_σ>0` 是 scheme-无关纯数；证明 `c_σ>0` 正是交接。

---

## 1. 精确定理

### 1.1 记号

- **群范围（G）**：`G` 为紧致连通单 Lie 群（Clay 要求任意紧单群；以下定理对任意此类 `G` 成立，`SU(N)` 为特例）。取基本（定义）表示 `f`、`d_f := dim f`；**`SU(N)` 时 `d_f = N`**（注：`N²−1` 是伴随维数/群维数，**不是** `d_f`）。
- **盒**：`Λ_{L,a} ⊂ aℤ⁴`（边长 `L`、格距 `a`）。QCD 标度另记 `Λ := Λ_QCD`（仅出现于 §4.4/§6 的比率 `√σ/Λ`、`m/Λ`，与盒记号不冲突）。
- **边界/通量扇区**：周期边界；若为环面则**固定平凡 't Hooft 通量扇区 `θ = 0`**（否则环面上的 Wilson 环期望依扇区而变）。
- **周长因子**：定理 1 的界为 `C_β e^{−σA}`（**无**周长因子 `K^{|C|}`）。这是**已证的加强**（`PROOF.md` §3.3：连通分支分解 + 闭曲面/体积精确相消 + 连通动物界），**非新假设**；GL4/GL9 的周长版 `K^{|C|}e^{−σA}` 可由之推出（取 `K≥2`）。建议台账记为"已证加强，GL4/GL9 可引用"。
- 连接变量 `U_e ∈ G`（有向边），Haar 测度 `∏_e dU_e`（`∫_G dU = 1`）。
- Wilson 作用量权重 `∏_{p} exp((β/d_f) Re tr_f U_p)`，**`β := 2 d_f/g²`**（`SU(N)` 退化为 `(β/N) Re tr U_p`、`β=2N/g²`），`U_p` 为 plaquette 的沿边有序积。
- 矩形 Wilson 环 `C = L×T`：`W(C) = (1/d_f) tr_f ∏_{e∈C} U_e`（`SU(N)` 为 `(1/N)tr`）。一般光滑闭曲线 C 同定义。
- 配分函数 `Z_{Λ_{L,a}}(β)`，期望 `⟨·⟩_{Λ_{L,a},β} = Z_{Λ_{L,a}}^{-1} ∫ (·) e^{−S} ∏ dU`。

### 1.2 定理 1（格点强耦合面积律；**已证**）

存在 `β₀ > 0`、`C_β < ∞`、`σ: (0,β₀) → (0,∞)`，使得对所有 `β ∈ (0,β₀)` 与所有可收缩闭曲线 `C`（含矩形 `L×T`）：

```
|⟨W(C)⟩_{Λ_{L,a},β}|  ≤  C_β · exp(−σ(β) · Area(C)),      σ(β) = −log(c_f(β)/c₀(β)) − γ > 0,
```

其中 `Area(C)` 为以 `C` 为边界的 2-链最小 plaquette 数（矩形 `L×T` 即 `L·T`），`c_f(β)` 为基本表示的特征标系数（见 §3 L1）。
**常数已钉**：`C_β := (1 − e^{γ} c_f(β)/c₀(β))^{-1}`，`γ = log(e Δ_p) + log d_f + 1`，`Δ_p ≤ 8`（4D plaquette 邻接度；见 `PROOF.md` 引理 3 与定理 A/B）。因 `c_f/c₀ → 0`（`β→0`），存在 `β₀` 使 `β<β₀` 时 `C_β ≤ 2`；`σ(β) := −log(c_f/c₀) − γ > 0`（上界速率，真实弦张力 ≥ 此）。
**周长因子（已证可去，非新假设）**：主界 `C_β e^{−σA}` **无**周长因子——`PROOF.md` §3.3 用连通分支分解 + 闭曲面/体积精确相消 + 连通动物界证得，是对文献周长版 `K^{|C|}e^{−σA}` 的**严格加强**；取 `K := max(e^{γ},2)` 即得周长版（`K^{|C|} ≥ 2 ≥ C_β`）。GL4/GL9 可直接引用本加强。

> 出处：Wilson 1974（物理）；Osterwalder–Seiler, *Ann. Phys.* 110 (1978) 440（公理化 + 面积律 + 质量隙）；Seiler 专著 LNP 159 (1982)。随机曲面严格版本：Cao–Park–Sheffield, arXiv:2307.06790（任意亏格嵌入曲面，全强耦合区，全 N）。扩展 β 区间：Cao–Nissim–Sheffield, arXiv:2505.16585 / 2509.04688。

### 1.3 定理 2（连续统面积律；**公开问题**）

设 `μ` 为 4D 纯 Yang–Mills（紧单群 `G`，`SU(N)` 为特例）的连续统测度（存在性见 `docs/yang-mills-prior-work.md` §3.3 的 4D 构造；满足 OS 公理，子列无限体积极限）。**弦张量**定义为双极限

```
σ  :=  − lim_{L→∞}  (1/L) · lim_{T→∞} (1/T) · log |⟨W(L×T)⟩_μ| .
```

则断言：

1. `σ > 0`（弦张量为正 ⇒ 线性势 ⇒ 禁闭）；
2. 存在常数 `C₁ < ∞`、`c < ∞`，使得对每个光滑简单闭曲线 `C`，
   `|⟨W(C)⟩_μ| ≤ C₁ · exp(−σ·Area(C) + c·Perim(C))`。

**等价表述（面积率）**：对求长矩形族，`σ(L) := −lim_T (1/T) log|⟨W(L×T)⟩|`，则 `σ = lim_{L→∞} σ(L)/L > 0`，且对每个固定 L，`σ(L) ≥ σ·L`（大 L 线性）。

`c·Perim` 项来自外源静夸克的自能（连续统中线性 UV 发散，需质量重整化后为有限 `c`）；**普适（universal）的禁闭内容只在面积项**，故 σ 的定义取 `L,T→∞` 双极限以孤立面积项。

### 1.4 定理 2′（条件化归约；可形式化的"若…则…"）

若下列交接判据 `(H3_σ)` 成立，则由定理 1 + 弱耦合 RG 传输推出定理 2：

> **(H3_σ)** 存在弱侧有效域边界 `μ_up` 与强侧边界，使 `√σ(μ)/Λ` 的不变不等式在整个尺度上闭合（详见 §4.4，与 `physics_handoff.md` 的 `q = m/μ` 判据同构，替换 `m → √σ`）。

即：**面积律的连续统问题 = 一个 RG 交接问题**，其结构与 MG 的交接**逐字同构**，只是不变量 certificate 由 `m/Λ` 换成 `√σ/Λ`。这是本报告最可操作的结构性发现。

---

## 2. 与 Clay 陈述 / 质量隙（MG）的关系

### 2.1 逻辑独立性（核心）

对反射正性（OS3）的格点规范理论，转移矩阵 `T = e^{−aH}` 存在。对 `T×L` 环，谱表示

```
⟨W(L×T)⟩  =  Σ_n |a_n(L)|² · e^{−E_n(L)·T},     a_n(L) = ⟨n|W_L|Ω⟩,
```

`W_L` 为沿空间 L-环的（闭合 ⇒ 规范不变）Wilson 算子，`E_n(L)` 为它生成的转移矩阵谱。定义**弦基态能量**

```
E₀(L)  :=  − lim_{T→∞} (1/T) log ⟨W(L×T)⟩  =  inf_n E_n(L)（主导项）.
```

则：

- **面积律 ⇔ `E₀(L) ≥ σ·L`（大 L 端）**：弦能量随夸克间距线性增长。
- **质量隙 `Δ` = 胶球（局域、规范不变算子）第一激发**（小 L / 局域端；等价地 `Δ = inf{ E_n(0) : n, 非真空 }`，`E_n(0)` 为局域算子的能级）。

二者是**同一转移矩阵谱的不同扇区 / 不同端**：

| | 面积律 | 质量隙（MG） |
|---|---|---|
| 扇区 | 带基本电荷的**弦/通量管** | 中性**胶球** |
| 谱对象 | `E₀(L) ≥ σL`（大 L） | `Δ = inf E_n(0) > 0`（局域） |
| 物理 | 外源夸克线性势（禁闭） | 胶球质量下界 |
| Clay | **不要求** | **要求** |

**结论：面积律与 MG 互不蕴含。** 排斥例（逻辑上允许，物理上被通量管图像排除）：
- 有隙无面积律：希格斯型屏蔽相（规范玻色子有质量 ⇒ `Δ>0`，但基本电荷被屏蔽 ⇒ 周长律）。→ `Δ>0` 不推出 `σ>0`。
- 有面积律无隙：禁闭但存在无质量胶球的相（弦仍线性，中性扇区无隙）。→ `σ>0` 不推出 `Δ>0`。

因此**完整解必须同时给 MG（中性隙）与面积律（弦张力）**："整个激发谱有隙"。仅给 MG 不覆盖 Clay 之外的禁闭物理；仅给面积律不覆盖 Clay。

### 2.2 已证的联合出现

在格点强耦合区，**同一簇展开**同时给出面积律（定理 1）与质量隙（Osterwalder–Seiler；SZZ 的 Langevin/Bakry–Émery 路线给出无限体积 + 正质量隙）。即：强耦合是已证的"胶球有隙 + 弦有张力"微扰区。连续统/弱耦合区二者都缺失，且缺的原因**相同**（§4.4、§6 O1）。

### 2.3 对 4D 构造的接入

4D 存在性构造给的是 `μ`（OS 测度）+ 中性扇区结构。面积律需要**额外**在 `μ` 上证明 `σ>0`，这**不在 16 个输入之内**，也不能由 MG 推出。故本条与 `nontriviality`/`mass-spectrum`/`uniqueness` 等 sibling 任务平行、不重叠。

---

## 3. 引理链 L0–L8

标记：✅已证 / 🟡半严格或条件 / ❌开放。依赖用箭头。

```
L0 ── L1 ── L2 ── L3 ── L4 ──┐
                             ├── 定理1（格点面积律）✅
L0 ── L5 ─────────────────────┘
                             
L4 + L5 ── 面积律谱表述（§2.1）✅/🟡
L4 + L7 ── 定理2（连续统面积律）❌
```

**L0（格点测度与环）✅** — 有限盒上的积测度 `∏ dU_e`（紧群 Haar，概率测度）、`Z_{Λ_{L,a}}`、`W(C)`。形式化友好度最高：纯有限维紧群积分，无测度论障碍。`Mathlib` 有紧群 Haar（`MeasureTheory.Measure.haar`）与 `PeterWeyl` 基础。

**L1（特征标展开）✅** — 由 `Peter–Weyl` 定理，
`exp((β/N) Re tr U) = Σ_R c_R(β) χ_R(U)`，`c_R(β) = d_R^{-1} ∫ exp((β/N)Re tr U) χ_R(U) dU`。
小 `β`：`c_R(β) ~ β^{|R|}/|R|!`（`|R|` 为表示阶），`c_f(β) = β/(2N) + O(β²)`，`c_f ≥ 0`，`c_f(β) ≤ Kβ`。基本表示 `f` 主导。

**L2（曲面表示 / 群积分正交）✅** — 由 `∫_G χ_R(U)χ_{R'}(U†) dU = δ_{RR'}`（`Peter–Weyl` 正交），带 Wilson 环的配分函数化为**二维曲面**求和：
`⟨W(C)⟩ = Σ_{S : ∂S = C} w(S),   w(S) ≥ 0,   w(S) ≤ ∏_{p∈S} c_{f}(β)`（`= (c_f/c₀)^{|S|}c₀^{|S|}`，`c₀^{|S|}` 并入 `γ`），
`S` 遍历以 `C` 为边界的 plaquette 集合（允许亏格、自交重数）。严格化（任意亏格、全 N、全强耦合区）：**Cao–Park–Sheffield, arXiv:2307.06790**。

**L3（曲面熵界）✅(强耦合)** — 以 `C` 为边界、面积 `s` 的曲面计数满足
`#{S : |S| = s, ∂S = C} ≤ exp(γ·(s − A(C)) + O(Perim(C)))`（平坦面主导，每多一个 plaquette 的"超额面积"被指数压制）。工具：Penrose 树图界 / Combes–Thomas / CPS 的曲面簇展开。**这是全部组合难点的所在**：熵界必须**对 loop 尺寸 L 与格距 a 一致**。

**L4（格点面积律）✅** — L2+L3 相加（连通分支分解 + 闭曲面/体积相消，见 `PROOF.md` §3.3）：`|⟨W(C)⟩| ≤ Σ_s e^{γ(s−A)}(c_f/c₀)^s ≤ e^{−σA}`，`σ = −log(c_f/c₀) − γ`（无周长因子）。即定理 1。

**L5（转移矩阵谱表示 / `E₀(L)`）✅** — 反射正性（OS3）⇒ `T = e^{−aH} ≥ 0` 存在；`T×L` 环的谱分解给出 `E₀(L)`（§2.1）。形式化友好度中（需 OS3 ⇒ 转移矩阵；`OSforGFF` 已有 OS 陈述层可复用）。

**L6（强耦合质量隙）✅** — 同一簇展开作用于规范不变关联：`|⟨O(x)O(y)⟩_c| ≤ C e^{−Δ|x−y|}`，`Δ(β)>0`（Osterwalder–Seiler）。或 SZZ 的 Langevin 指数遍历性（`arXiv:2204.12737`, CMP 400）给出无限体积唯一测度 + 正质量隙。

**L7（连续统标度 `c_σ > 0`）❌** — 量纲 transmutation：纯 YM 唯一尺度 `Λ`，故 `σ = c_σ Λ²`，`c_σ` 为 scheme-无关纯数。需证 `lim_{a→0} a² σ_lat(β(a)) = c_σ Λ²` 存在且 `c_σ>0`，其中 `β(a) = 2d_f/g(a)²`（`SU(N)` 时退化为 `2N/g²`，值不变），`g(a)² ~ 1/(b₀ log(1/(aΛ)))`（渐近自由，`b₀ = 11 C₂(adj)/(48π²)`，`SU(N)` 时 `=11N/(48π²)`）。**方案约定**：scheme 变换 `g'²=g²+C g⁴` 的常数 `C` **与 `A0 §1.4` 的 `RT` 同一 RT、同一 `C`**（不另立；`Λ'=Λ e^{C/(2b₀)}` 的 scheme 协变随之与 A0 一致）。**要求把 L4 的强耦合结果沿 RG 流传输到 `β→∞`**——正是交接墙。

**L8（定理 2）❌** — L4 + L7 + L5。连续统面积律。

**证明骨架**：`L0 → L1 → L2 → L3 → L4` 给出强耦合面积律（干净、可形式化）；`L0 → L5` 给出扇区分离（§2.1，可形式化）；**唯一缺口是 L7**，与 MG 的交接缺口同址。

---

## 4. 五条路线综述

### 4.1 强耦合格点展开（经典、已证、到不了连续统）

`e^{β tr U_p}` Taylor/特征标展开，逐 plaquette 展开后群积分正交把项组织成**二维弦世界面**。收敛半径 `β₀>0`。产物：面积律（定理 1）、质量隙、OS 公理（Osterwalder–Seiler 1978、Seiler 1982、Münster–Weisz 综述）。
**卡点**：连续统要求 `β → ∞`（渐近自由 `g→0`）。强耦合展开的收敛半径有限，**结构上到不了连续统**。但它给出**最干净的可形式化核心**（L0–L4）。

### 4.2 Cao–Nissim–Sheffield：随机曲面与扩展面积律区

- **Cao–Park–Sheffield**（arXiv:2307.06790）：任意亏格嵌入曲面的**严格求和表示**，全强耦合区、全 N——把 L2 做成了定理。
- **Cao–Nissim–Sheffield**（arXiv:2505.16585, 2509.04688）：把面积律推到**更大的 β 区**（含 't Hooft 大 N 区），用"随机曲面的表面张力"与**动力学/熵**论证（2509 为 dynamical approach）。
**价值**：扩大了已证的 `β` 区间，是当前"把强耦合结果往弱耦合推"的最前沿严格工作。**卡点**：仍 `β` 有上界，未到连续统；具体 β 阈值与大 N 细节本报告未逐条核实。

### 4.3 Shen–Zhu–Zhu：随机分析 / Langevin

- SZZ（arXiv:2204.12737, CMP 400 (2023) 805–851）：任意维格点 YM 强耦合，Langevin 动力学的**指数遍历性**（Bakry–Émery）⇒ 唯一不变测度、无限体积极限、**正质量隙**。
**价值**：质量隙侧最强的概率路线；与面积律的"禁闭侧"互补。**卡点**：强耦合；4D 弱耦合的 Langevin 需无穷多 counterterm，regularity structures 不能直接用。

### 4.4 RG（Balaban；Bauerschmidt–Bodineau–Dagallier）

- **Balaban**（CMP 1980s）：4D 格点 YM 的 **UV 稳定性** + 小盒子有效作用量控制。**未**到 IR/面积律。
- **Bauerschmidt–Bodineau–Dagallier**（arXiv:2307.07619 综述，Probab. Surveys 2024；Polchinski flow / multiscale Bakry–Émery）：U(1)₄ 格点规范**弱耦合**的一致 log-Sobolev 不等式——目前最接近的**弱耦合**严格动力学控制。**非阿贝尔 4D 开放**。
- **交接判据 (H3_σ) 的结构**：与 `physics_handoff.md` 的 `q = m/μ` 完全同构。取 `Q(μ) := √σ(μ)/μ = √σ/Λ · Λ/μ`，则：
  - 弱侧 `μ≫Λ`：`√σ ≪ μ`（`Q→0`，RG 展开）；
  - 强侧 `μ→0`：`Q → ∞`（`√σ` 固定、`μ→0`），或在格点上 `Q ~ √σ_lat · a`；
  - 交接 ⇔ 两域重叠 ⇔ 若 `√σ/Λ = c_σ^{1/2}` 足够大以容下两侧半径。
  **与 MG 的差别仅在不变量：`m/Λ` → `√σ/Λ`。** 故"禁闭 + 质量隙"是**同一个交接、两个 certificate**。

### 4.5 中心涡旋与 't Hooft 环（严格的是二择一，不是证明）

- **'t Hooft 1978**：Wilson 环 `W(C)` 与 't Hooft 磁环 `W'(C')` 满足中心对易代数 `⟨W W'⟩ = e^{2πi/N}⟨W'W⟩`。由此**严格**得：`W` 面积律与 `W'` 面积律**不能同时**成立（至多一个面积、另一个周长）。这是**二择一约束**，不是"哪个成立"的证明。
- **中心对称判据**：`SU(N)` 的格点测度有中心 `Z_N` 对称；**若中心对称不破**，则 `W` 取面积律（`W'` 取周长律）。这是**条件性判据**——证明"连续统中中心不破"是开放问题。
- 严格成立的具体模型：**4D `Z_2` 规范理论**（Wegner 1971）——**自对偶**（4D 中 plaquette 的对偶仍为 plaquette），相变（设唯一）位于自对偶点 `sinh(2β_c)=1 ⇒ β_c = ½ln(1+√2) ≈ 0.4407`；`β<β_c` 面积律、`β>β_c` 周长律。**这是最有价值的严格"面积律相变"范本**，但 `Z_N`/`SU(N)` 无精确自对偶。
  - **勿混淆**：与 **3D Ising** 对偶的是 **3D** `Z_2` 规范理论（Kogut 1979），其 `β_c = arctanh(e^{−2K_c^{Ising}}) ≈ 0.761`（`K_c^{Ising}≈0.2217`），**不是** `0.4407`；4D 与 3D 的 `β_c` 互不相同。
  - **数值旁证（非证明，L=4/6，Loment）**：4D `Z_2` 模拟（`D:/Dev/math/loment-exp/z2b/`）在 `β≈0.51–0.52` 出现有限体积伪相变——L=4：plaquette 由 `0.53`(β=0.49) 跳到 `0.84`(β=0.52)，`χ(2,2)` 由 `0.486` 塌到 `0.011`；L=6：`β=0.50→0.54` 间 plaquette `0.56→0.87`、`χ(2,2)` `0.45→0.006`。`β=0.4407` 处 L=6 仍 `plaquette≈β=0.44`、`χ≈0.76`（强耦合/禁闭侧）。伪相变点高于 `β_c≈0.44`，符合一阶相变的大有限体积移位与滞后。**数值证据，非定理。**
- **涡旋凝聚图像**（'t Hooft、Engelhardt–Reinhardt、Greensite）：物理上中心涡旋的随机凝聚给出面积律；**严格化仅限 `Z_N` 强耦合**。

### 4.6 单极凝聚（dual superconductor）

- **3D 紧致 U(1)**：**严格**——Polyakov 1977（单极凝聚）、**Göpfert–Mack 1982**（严格面积律 + 质量隙）。这是唯一**完全严格**的 4D 以下"凝聚 ⇒ 面积律"范例。
- **4D 纯 YM**：单极凝聚（'t Hooft abelian projection、dual superconductor）是 **conjectural**；Seiberg–Witten 在软破缺 `N=2` 已解，纯 YM 未严格化（Douglas 称"迄今最好的解析处理"）。
**价值**：给出**机制**与 3D 的严格范本；4D 是图像而非定理。

### 4.7 其他：2D 精确、对偶、条件判据

- **2D YM**：Migdal 1975 / Gross–Witten 1980 **精确**面积律，Driver 1989（连续统期望、格点收敛）、Lévy 2003、Sengupta 1997。**形式化的理想第一个目标**（测度显式、面积律精确、OS 可证）；但 2D 无传播自由度、质量隙退化，不直接推 4D。
- **对偶**：4D `Z_2` 规范理论**自对偶**（Wegner 1971）；3D `Z_2` 规范理论 ↔ **3D Ising**（Kogut 1979，`β_c≈0.761`）。非阿贝尔无精确对偶。
- **Chatterjee**（*Gauge Theory* 书 / *Yang–Mills for probabilists*）：把面积律表述为**条件判据**（单调性/对偶型充分条件），是框架而非 4D 证明。

---

## 5. 已证 / 未证状态表

| 陈述 | 状态 | 出处 |
|---|---|---|
| 格点强耦合面积律（全 N，`β<β₀`）| ✅ 已证 | Osterwalder–Seiler 1978；Seiler 1982；CPS 2307.06790 |
| 随机曲面表示（任意亏格，强耦合）| ✅ 已证 | Cao–Park–Sheffield 2307.06790 |
| 扩展 β 区间面积律（含大 N/ 't Hooft）| ✅ 已证（区间接有上界）| Cao–Nissim–Sheffield 2505.16585, 2509.04688 |
| 格点强耦合质量隙 + 无限体积 | ✅ 已证 | Osterwalder–Seiler；SZZ 2204.12737 |
| 3D U(1) 面积律 + 隙（单极凝聚）| ✅ 已证 | Polyakov 1977；Göpfert–Mack 1982 |
| 4D `Z_2` 强耦合面积律 + 自对偶定相变 | ✅ 强耦合严格已证；相变由自对偶定位（数值 L=4/6 见 §4.5）| Wegner 1971 |
| 2D YM 精确面积律 | ✅ 已证 | Migdal 1975；Gross–Witten 1980；Driver 1989 |
| 弱耦合 U(1)₄ 一致 LSI（abelian）| ✅ 已证 | Bauerschmidt–Bodineau–Dagallier 2307.07619 及其 U(1) 文 |
| `W`/`W'` 面积律二择一 | ✅ 已证（约束，非选择）| 't Hooft 1978 |
| **4D 非阿贝尔连续统面积律** | ❌ **公开** | — |
| 连续统中心对称不破 | ❌ 公开 | — |
| 4D 纯 YM 单极凝聚 | ❌ conjectural | 't Hooft；Seiberg–Witten（软破缺） |
| 弱耦合非阿贝尔 4D 一致 LSI | ❌ 公开 | — |

**非空洞性测试**（陈述层正确性，用户要求的"反例测试"）：
- ✅ `Z_2`/`SU(N)` 强耦合**应通过**（定理 1）；
- ✅ **`U(1)₄` 应失败**面积律——它有无质量光子（Guth 1980、Fröhlich–Spencer 1983 的 Coulomb 相），面积律为假。陈述必须写明"非阿贝尔"或 `Z_N`，否则被 `U(1)` 反例击穿；
- ✅ 自由场（无相互作用）**应失败**（`σ=0`，周长/平凡律）。

---

## 6. 障碍（按致命性排序）

**O1（致命，与 MG 共享）—— UV→IR RG 交接。** 连续统要求 `β→∞`，强耦合展开 `β<β₀` 结构上到不了；Balaban 只做 UV 稳定性，不做 IR；非阿贝尔 4D 弱耦合无一致动力学控制（`U(1)` 有、非阿贝尔无）。**面积律的墙 = MG 的墙（`physics_handoff.md` 的 `q₀ ≥ q₁`），只是把 `m` 换成 `√σ`。**

**O2 —— 弦张量的连续统标度与 UV 减除。** 需 `σ = c_σ Λ²`，`c_σ>0` 且极限存在、有限；周长项 `c·Perim` 来自外源自能（连续统线性 UV 发散），须质量重整化后孤立面积项。`c_σ` 是 scheme-无关纯数（transmutation），故问题是"纯数是否为正"。

**O3 —— 中心对称在连续统中是否不破。** 中心判据（`W` 面积 ⇔ `Z_N` 不破）是**条件性**的；证明"连续统里中心不破"本身开放。't Hooft 代数只给二择一、不给选择。

**O4 —— 拓扑扇区不受簇展开控制。** 瞬子、单极、涡旋是"大场"构型，polymer/簇展开在本质上是**小场**工具；4D 的拓扑扇区如何进入 `μ` 并影响 `σ` 未受控。

**O5 —— 曲面熵界对 L 与 a 的一致性。** L3 的熵界须同时均匀于 loop 尺寸 `L` 与格距 `a`；CPS 在强耦合做到，弱耦合无。这是"把面积律沿 RG 传输"的技术前提。

**O6 —— 与 MG 的逻辑独立性。** 即使 MG 解决，也不自动给面积律（§2.1）；完整解需**两个**独立证明（中性隙 + 弦张力），虽共享 O1。

---

## 7. 下一具体步骤

按"可验收、可形式化"排序（与 M0–M6 计划对齐，`docs/yang-mills-prior-work.md` §5.3）：

**N1（2–4 周，可证，Lean 目标；对应 M2 的面积律版）。** 形式化**格点强耦合面积律**（定理 1）。路线：先 **`Z_2` 规范理论、有限格点**（群 `Z_2` 离散、Haar = 均匀，`Peter–Weyl` 退化为简单傅里叶），再 **`SU(2)`**（紧群 Haar + `Peter–Weyl`）。交付：

```lean
-- 目标签名（有限格点 Λbox，固定）
theorem strong_coupling_area_law
    (N : ℕ) (hN : 2 ≤ N) (Λbox : Finset (Fin 4 → ℤ)) :
    ∃ β₀ > 0, ∃ σ : ℝ → ℝ, ∃ Cβ : ℝ → ℝ, (∀ β ∈ Ioo 0 β₀, 0 < σ β) ∧
      ∀ (C : RectLoop Λbox) (β : ℝ), β ∈ Ioo 0 β₀ →
        |wilsonLoopExpectation Λbox β C| ≤
          Cβ β * Real.exp (-(σ β) * area C)
```

复用 M2 的 polymer/character 机器（同一簇展开给面积律与质量隙，故 `confinement` 与 `mass-spectrum` 两个 session 应共享 `character expansion` 基础库）。

**N2（2–3 周，半严格，Lean 目标）。** 形式化**扇区分离命题**（L5，§2.1）：`OS3 ⇒ 转移矩阵存在`，`⟨W(L×T)⟩ = Σ |a_n|² e^{−E_n T}`，并从面积律推出 `E₀(L) ≥ σL`。这是面积律与 MG 的**逻辑桥**，直接服务 `mass-spectrum` session 的 `E₀(L)` 对象。

**N3（长期，概念 + 陈述）。** 把定理 2 写成**条件定理**：`(H3_σ)（§4.4 的 √σ/Λ 不变不等式）+ 定理 1 + RG 传输 ⇒ 连续统面积律`。结构上复制 `physics_handoff.md` 的 `q` 判据，把 `m/Λ` 换成 `√σ/Λ`。

**N4（协调）。** 交 `physics_handoff`：`√σ/Λ`（或 `σ/Λ²`）作为与 `m/Λ` **并列**的不变量 certificate，共用 O1 的交接判据。交 `mass-spectrum`：`E₀(L)`（弦基态能量）与其 `Δ`、`m_glueball` 的关系（§2.1 双端）。交 `uniqueness`：面积律是"带电荷扇区"的聚类，与中性扇区聚类（MG）独立。

**可走的更粗但严格的旁路**：`Z_2`/`SU(N)` 强耦合 + CNS 扩展区 + `Z_2` 自对偶，可形式化出一条**完整的严格面积律链**（不需 RG），作为"形式化规范理论"的示范工件；连续统部分如实标为开放。

---

## 8. 文献

**禁闭 / 面积律（经典）**
- K. Wilson, *Confinement of quarks*, Phys. Rev. D 10 (1974) 2445.
- K. Osterwalder, E. Seiler, *Gauge field theories on a lattice*, Ann. Phys. 110 (1978) 440–471.
- E. Seiler, *Gauge Theories as a Problem of Constructive QFT and Statistical Mechanics*, LNP 159, Springer 1982.
- G. Münster, P. Weisz, *On the convergence of the strong coupling expansion*, 及后续（综述性，卷号/年份本报告未逐条核实）。

**随机曲面 / 扩展面积律（CNS）**
- S. Cao, M. Park, S. Sheffield, *Random surfaces and lattice Yang–Mills*, arXiv:2307.06790, Comm. AMS (2025).
- S. Cao, R. Nissim, S. Sheffield, *Expanded regimes of area law for lattice Yang–Mills theories*, arXiv:2505.16585 (2025).
- S. Cao, R. Nissim, S. Sheffield, *Dynamical approach to area law for lattice Yang–Mills*, arXiv:2509.04688 (2025).

**随机分析 / Langevin（SZZ）**
- H. Shen, R. Zhu, X. Zhu, *A stochastic analysis approach to lattice Yang–Mills at strong coupling*, arXiv:2204.12737, CMP 400 (2023) 805–851.

**RG / 弱耦合**
- T. Balaban, *Renormalization group approach to lattice gauge field theories* I–IV, CMP 1980s（卷号部分未核实）。
- R. Bauerschmidt, T. Bodineau, B. Dagallier, *Stochastic dynamics and the Polchinski equation: an introduction*, arXiv:2307.07619, Probab. Surveys 21 (2024).
- （U(1)₄ 弱耦合一致 LSI 的 BBD 系列文，arXiv 号本报告未核实。）

**中心涡旋 / 't Hooft 环 / 对偶**
- G. 't Hooft, *On the phase transition towards permanent quark confinement*, Nucl. Phys. B 138 (1978) 1.
- F. Wegner, *Duality in generalized Ising models and phase transitions without local order parameters*, J. Math. Phys. 12 (1971) 2259.
- J. Kogut, *An introduction to lattice gauge theory and spin systems*, Rev. Mod. Phys. 51 (1979) 659.
- M. Engelhardt, H. Reinhardt；J. Greensite, *An Introduction to the Confinement Problem*（涡旋凝聚，综述）。

**单极凝聚**
- A. Polyakov, *Quark confinement and topology of gauge theories*, Nucl. Phys. B 120 (1977) 429.
- M. Göpfert, G. Mack, *Proof of confinement of static quarks in 3-dimensional U(1) lattice gauge theory for all values of the coupling constant*, CMP 82 (1982) 545.
- N. Seiberg, E. Witten, *Monopole condensation and confinement in N=2 SUSY QCD*, 1994.

**2D 精确**
- A. A. Migdal, *Recursion equations in gauge theories*, Sov. Phys. JETP 42 (1975) 413.
- D. Gross, E. Witten, Phys. Rev. D 21 (1980) 446.
- B. Driver, *YM₂: continuum expectations, lattice convergence, and lassos*, CMP 123 (1989) 575.
- T. Lévy, Mem. AMS (2003)；A. Sengupta, Mem. AMS (1997).

**U(1) 相与反例**
- A. Guth, *Existence proof of a nonconfining phase in four-dimensional U(1) lattice gauge theory*, Phys. Rev. D 21 (1980) 2291.
- J. Fröhlich, T. Spencer, *Massless phases and symmetry restoration in abelian gauge theories and spin systems*, CMP 83 (1982) 411.

**条件判据 / 综述**
- S. Chatterjee, *Gauge Theory*（书）；*Yang–Mills for probabilists*, arXiv:1803.01950.
- M. R. Douglas, *The Yang–Mills Millennium problem*, Nature Rev. Phys. 8 (2026) 86–97.
- 项目内：`docs/yang-mills-prior-work.md`；`research/physics_handoff.md`（交接控制量 `q`）。

---

## 附录 A：一页速查

```
面积律 = σ>0 弦张力: |<W(C)>| ≤ C1 e^{-σA+c·Perim}
与MG：独立（弦扇区 E0(L)≥σL vs 中性胶球隙 Δ>0）；完整解需两者；共享 RG 墙
已证：格点强耦合面积律（OS78/Seiler82/CPS23）、扩展区（CNS25）、3D U(1)（GM82）、4D Z2（Wegner79）、2D精确
未证：4D 非阿贝尔连续统面积律（公开）；中心对称不破（公开）；4D 单极凝聚（conjectural）
引理链：L0 测度 → L1 特征标 → L2 曲面表示 → L3 熵界 → L4 格点面积律（全 ✅）
         L0 → L5 转移矩阵谱（✅）；L4+L5 → 扇区分离（✅/🟡）
         缺口 = L7 连续统标度 c_σ>0（❌ = 交接墙）
不变量：σ=c_σ Λ²，√σ/Λ 与 m/Λ 并列（结构性发现）
反例测试：U(1)4 必须失败（有光子）；自由场必须失败（σ=0）
下一：N1 形式化格点强耦合面积律（Z2→SU2）；N2 形式化扇区分离；N3 条件定理
```
