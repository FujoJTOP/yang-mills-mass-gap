# W95C-ORDER — 语料里的质量序：`0⁺⁺/2⁺⁺/0⁻⁺/1⁻⁻`，`a_c=1`，以及 `m>0` 需要哪几个道

> **判词（先给）**：**(b)。** 语料**有**定论的质量序只有一条，且它是**简并的**：强耦合 `β<β₁` 下 `a_{0⁺⁺}=a_{2⁺⁺}=1`、同一条四壁管 ⇒ `m_{2⁺⁺}/m_{0⁺⁺}→1` `[T]`（`PROOF-IN11.md:23`），其它一切 `J^{PC}` 严格更重且差距发散 `[T]`（`:13`）。因此 `inf spec` 由 `{0⁺⁺,2⁺⁺}` **共同**决定，`m>0` **至少要 `0⁺⁺` + `2⁺⁺` 两道**。`(a)` 的前提「`0⁺⁺` 严格最轻」在语料里**不是定理、且领头阶为假**；`(c)` 不成立（强耦合区有序）——无定论的只是 `β≈2.5` 侧，那里只有 `[J]`。
>
> 标 `[T]` 真定理 · `[C]` 条件 · `[J]` 格点/文献经验值 · `[M]` 模型 · `[NG]` no-go · `[O]` 精确缺失步 · `[E]` 精确/数值复核。Clay 0%。

---

## 1 · 质量序表（逐条给出处）

### 1.1 强耦合侧（`β<β₁`，`β₁=min(1/48,c/N)`，`λ>96`）：**唯一有定理的序**

| 道 | `a_c`（每单位时间最小 plaquette 数） | 相对 `0⁺⁺` 的质量 | 标签 | 出处（逐字） |
|---|---|---|---|---|
| **`0⁺⁺`** | **1** | 基准 | `[T]` | `PROOF-IN11.md:11` `a_{0⁺⁺}=a_{2⁺⁺}=1`；`:9` `m_c=a_c\ln(1/u)+O(1)` |
| **`2⁺⁺`** | **1** | **`m_{2⁺⁺}/m_{0⁺⁺}→1`（领头阶精确简并）** | `[T]` | `PROOF-IN11.md:23` 逐字 `The two anomaly-protected channels share the same minimal flux tube (one plaquette per time slice), so m_{2⁺⁺}/m_{0⁺⁺} → 1` |
| **`0⁻⁺`** | `≥2` | `m_{0⁻⁺}−m_{0⁺⁺} ≥ \ln(1/u)+O(1) → +∞` | `[T]` | `PROOF-IN11.md:11` `0⁻⁺ requires the ε-contraction tr FF̃ of two field strengths`；`:13` 不等式逐字 |
| **`1⁻⁻`** | **`≥3`（C-奇需 ≥3 胶子）`[T]`/`[C]`**（语料未给显式数值，见下） | `m_{1⁻⁻}−m_{0⁺⁺} ≥ 2\ln(1/u)+O(1) → +∞`（**严格更重**） | `[T]`+`[C]` | `PROOF-IN11.md:11` 逐字 `every other J^{PC} needs a_c ≥ 2 (… C-odd needs ≥3 gluons …)`；`1⁻⁻` 是 C-奇、`J=1`（`K0-TRANSVERSE-DECOUPLING.md:15` 逐字 `the physical gluon is J^{PC}=1⁻⁻`、`:19` `The tachyon is 1⁻⁻`）⇒ 落该子句（**该归属是本文件的推断，标 `[C]`**）。**但语料从未给 `1⁻⁻` 的通道专有陈述** |

**⇒ 强耦合序**：`0⁺⁺ = 2⁺⁺`（简并，最轻）`≺ 0⁻⁺ ≺ 1⁻⁻`（后者由 `a_c≥2` 的统一界 + C-奇 `≥3` 胶子给出 `[T]/[C]`）。`PROOF-IN11.md:17` 给出这条序的用途：`inf_c Δ_c = min(Δ_{0⁺⁺}, Δ_{2⁺⁺})` `[T]`。
**通道穷举的原文**：`PROOF-IN11.md:11` —— 单 plaquette 算符在 `O_h` 下 `A1g⊕Eg⊕T2g = 0⁺⁺⊕2⁺⁺`，**其它一切 `J^{PC}` 需 `a_c≥2`**（`0⁻⁺` 要 `ε`-收缩，C-奇要 ≥3 胶子，`J≥3` 要横向宽度）；同文 `Y02-CHANNELS.md:43–44` 逐字转引。⇒ 简并**不是**偶然丢了一项，而是**同一算符**的 `O_h` 分解，`2⁺⁺` 无法靠「多一张 plaquette」被推开。

### 1.2 格点/物理点 `β≈2.5`：**只有 `[J]`，不是定理**

| 道 | 值 | 标签 | 出处 |
|---|---|---|---|
| `0⁺⁺` | 1710 MeV，**最轻** | `[J]` | `docs/conditional-theorem.md:208`；`PROOF-IN11.md:27`；`KT-RHO2.md:65`；`R6-ATOM-ADVERSARIAL.md:300` |
| `2⁺⁺` | 2390 MeV，次之（`m₂/m₀=1.39`） | `[J]` | 同上（`conditional-theorem.md:208` 逐字 `0⁺⁺ 1710 MeV 最轻、2⁺⁺ 2390 MeV 次之（与 0⁺⁺ 之比 1.39）`） |
| `0⁻⁺` | 2560 MeV，第三 | `[J]` | `KT-RHO2.md:65` 逐字 `0⁺⁺ 1710 < 2⁺⁺ 2390 < 0⁻⁺ 2560 MeV 只能作 [J] 引用（GAP-IN11）` |
| `1⁻⁻` | **语料无格点值**（位置由强耦合侧 §1.1 定：严格更重） | — | `Y02-CHANNELS.md:55`；全库 `1⁻⁻` 命中 27，无一条给质量数值或位次 |

- **为什么不升级成定理**：`PROOF-IN11.md:27` 逐字 —— 强耦合序只解析到 `β₁`（`λ>96`），`β≈2.5` 的序是 **controlled extrapolation, not a theorem**，`No rigorous continuation exists`。
- **无对称性保护**：`docs/campaign-final.md:31` 逐字 `无可对称/反常论证能锁死最轻态量子数（T1：迹反常只定插值算子、不定质量序）`；`conditional-theorem.md:208` 逐字 `不保证最轻态是自旋-0；无 Ward 恒等式/anomaly 固定最轻 channel`；`PROOF-IN11.md:31` 逐字 `None orders the degenerate pair.`（1-形式中心盲、KL 道内、大 `N` 同迹同号）。
- **量级差与判据的方向无关**：1.39（格点）vs 1.00（强耦合领头）——**两者不同 regime**，不可拼接；`1.39` 不能用来替代 `m>0` 需要的下界（下界不是比值）。

### 1.3 各道**能给出什么**（这决定了第 4 节的判）

| 道 | `f₁²` 有定义？ | 有下界定理？ | `I2` 型隔离？ | 出处 |
|---|---|---|---|---|
| `0⁺⁺` | 有（`ν`=`tr F²` 的 KL 测度） | `f₁²_L ≥ r₀(β)>0`，`β<1/48`，`L`-一致 `[T] mod OS/KL`（审计已降级 `[C]`） | 有形式、**未闭合**（`I2` 仍是开的精确步） | `Y02-CHANNELS.md:52`；`V03-I2GAP.md:127` |
| `2⁺⁺` | **无**（只出现 `m_{2⁺⁺}`，从不出现 `f₁²_{2⁺⁺}`） | **无残数**；只有质量序（简并）`[T]` | **无**；语料明写 `relocates the demand to H1 (2⁺⁺ isolation)` | `Y02-CHANNELS.md:53`；`BLAST-EXCLUDE-BUMP.md:40` |
| `0⁻⁺` | 无 | **只有上界** `m²_{0⁻⁺} ≤ χ_t/f₁²` `[T]/[NG]` | 无 | `Y02-CHANNELS.md:54`；`BREAK-IN-11.md:20,27` |
| `1⁻⁻` | 无 | **无** | 无 | `Y02-CHANNELS.md:55` |

`m_{0⁻⁺} ≥ m_{0⁺⁺}`（`0⁻⁺` 的下界）**不可证**：`BREAK-IN-11.md:3–4` 逐字 `Every route gives an upper bound on m_{0⁻⁺} … m_{0⁻⁺} ≥ m_{0⁺⁺} stays physical input.`

---

## 2 · `a_c = 1` 是什么；`m₂/m₀→1` 在什么极限

**`a_c` 的定义（逐字）**：`PROOF-IN11.md:11` —— `a_c` = **“the minimal plaquette content per unit *time* of the lightest tube interpolating channel `c`”**，质量展开为

```
m_c(β) = a_c · ln(1/u(β)) + O(1),   u(β) → 0⁺  as β → 0
```

（`:9`。`u` 是 fundamental character coefficient，`u ≍ β`，`U01-F1AUDIT.md:111` 逐字；数值 `u(1/48)=5.208e-3`，`T23-TEMPLATE.md:63–64` `[V]`。）

**⇒ `a_c = 1` 不是临界 `β`，是**「**一根时间步内只含一张 plaquette**」这个**组合计数**，即 `m_c` 的**对数发散系数**。`β` 侧的临界量是 `β₁`（解析半径，`λ>96 ⇒ β<1/48`；`PROOF-IN11.md:15,27`），与 `a_c` **无关**。把 `a_c` 读成 `β_c` 是坐标错配。
（注：`P3` 审计把 `m₁` 由 `ln(1/u)/a` 改成 `4ln(1/u)/a` ⇒ `a_c=1→4`；`P3-AUDIT-CORRECTED.md:128` 逐字 `the factor 4 is common to both channels`，故**简并结论不变**，只改系数 `[O]`。）

**`m₂/m₀ → 1` 的极限**：**强耦合 `β→0`（`u→0`）的领头阶**，且**只是因为两条道共用同一条四壁管**（`a₀=a₂=1`）。逐字 `PROOF-IN11.md:23`：`The two anomaly-protected channels share the same minimal flux tube (one plaquette per time slice), so m_{2⁺⁺}/m_{0⁺⁺} → 1 and continuity/analyticity deliver nothing.`
**极限的精度**：`O(1)` 余项**未被压住**——`ln(1/u)→∞` 只让**比值**趋 1，`m₂−m₀` 的符号仍由 `O(1)` 的次领头分裂决定，`PROOF-IN11.md:23` 逐字 `The ordering and isolation of 0⁺⁺ live entirely in the O(1) subleading splitting … not carried out … and confined to β<β₁`。⇒ **简并是领头阶的、符号未定的**：`0⁺⁺` 略轻 / `2⁺⁺` 略轻 **两者都与语料相容**。
**数值复核（模型量，`[E]`，`w7/w95c_order.py|.log`）**：先对齐语料锚点 —— `u(β)=I₂(β)/I₁(β)` 给 `u(1/48)=5.208239e-3`，与 `T23-TEMPLATE.md:64` 的 `5.208e-3` 相对偏差 `4.6e-5` ✔。再取**示意性**的对称 `O(1)` 分裂 `c₀=−C, c₂=+C`（余项符号**未知**，见 `PROOF-IN11.md:23`）。**定理窗口内（`β≤1/48`）**：`C=1` 时 `m₂/m₀` = `1.47`（`β=1/48`）→ `1.27`（`β=10⁻³`）；`C=0.5` 时 `1.21` → `1.13`。`a_c=2` 的道（`0⁻⁺`）则 `m_c−m₀ = \ln(1/u)` 从 `1.43`（`β=1`）单调升到 `5.26`（`β=1/48`）、`8.29`（`β=10⁻³`）；`a_c=3`（`1⁻⁻`）在同一列乘 2。**⇒ 三点**：(i) `m₂/m₀→1` 是**对数慢**的渐近陈述——整个定理窗内比值仍在 `1.1–1.5`，**简并不是窗口内的数值事实，是 `u→0` 的极限事实**，故不能拿 `1.39`（格点）或任何有限 `β` 比值把它"验掉"；(ii) `a_c≥2` 的道**严格发散**，`0⁻⁺` 永远不参与竞争；(iii) `β≳1/48` 的列（`β=1` 给 `5.69`）只作趋势参考：`u=I₂/I₁` 的强耦合领头式在窗外已不可信，本脚本不对窗外下任何判。

---

## 3 · 判：「`m₂/m₀→1` ⇒ 简并 ⇒ `inf spec` 由两者共同决定 ⇒ `m>0` 要两者都有」——**成立**

三步，逐条核：

1. **`inf spec` 是「道底取 min」**：`PROOF-IN11.md:17` `inf_c Δ_c = min(Δ_{0⁺⁺}, Δ_{2⁺⁺})` `[T]`（强耦合）；`conditional-theorem.md:202` 逐字 `cert` 取两个 channel 界的 `min`。**min 的数学**：`min(x,y)>0 ⟺ x>0 ∧ y>0`（`w95c_order.log §3` 的 `[E]`：`inf` 只在两者同正时为正）。
2. **简并 ⇒ 不能靠「谁更轻」把另一道扔掉**：`0⁺⁺` 的隔离/排序是 `O(1)` 次领头动力学，**无对称性可定**（`PROOF-IN11.md:23,31`；`campaign-final.md:31`）。故「`0⁺⁺` 严格最轻」这个前提在语料里**既非定理、又在领头阶被否定**（领头阶是简并）。⇒ 用 `0⁺⁺` 单道覆盖 `inf spec` 的那一步**无依据**。
3. **但「共同决定」到「必须两条证书」之间还差一步，语料明写它没做**：若简并且**两个插值算符都看见同一态**（两个 `f₁²>0`），则 `Δ₀=Δ₂=m`，证一个就够。语料**没有**这条：`Y02-CHANNELS.md:53` 逐字 `2⁺⁺ … 无（只出现 m_{2⁺⁺} 质量，从不出现 f₁²_{2⁺⁺}）`；`:69` 逐字 `(2⁺⁺…(f_1|residue|atom) → 12 命中，无一条给 2⁺⁺ 的 f₁²>0`。⇒ 「`2⁺⁺` 也看见该态」本身是 `[O]`（`Y02-CHANNELS.md:95–97`：`2⁺⁺` 的 `I1` 只差把同一四壁管算一遍，**结构可行、语料未写**）。
   **⇒ 保守且语料可支持的正解：两条都要**；「证明一条即可」需要额外一条 `f₁²_{2⁺⁺}>0`（`[O]`），而它**不改变** `cert=min(cert₀,cert₂)` 的记账（`IN-10-attack.md:7` 逐字 `cert μ … defined as min(cert_{0⁺⁺}(μ), cert_{2⁺⁺}(μ))`）。

**反证式见证（`[E]`，`w95c_order.log §3`）**：取 `Δ₀=δ>0`、`Δ₂=0` ⇒ `inf spec=0`：**`0⁺⁺` 被证书覆盖而隙仍为 0**。这正是「只证 `0⁺⁺` 不够」的最小反例，且与 `T6-5(i)` 的 `m₁=0, I=T, G=F` 见证同构（`V03-I2GAP.md:120`）。反向：`Δ₀=0, Δ₂=δ` ⇒ 同样 `0`。⇒ **`min` 的两个因子任一为 0，隙即失**。

**「道盲证书」这条逃生口也没有**：`C2-bypass-INb.md:13` 逐字 —— 需要 `cert_blind ≤ inf_c Δ_c = Δ`（零松弛的隙强度），`The known channel-blind mechanisms … are exactly the IN-3 wall: rigorous only below the crossover. Not constructible today`。⇒ 不存在「用一条不分道的证书同时覆盖两道」的替代。

---

## 4 · 判：`m>0` 到底需要哪几个道

**答：`0⁺⁺` + `2⁺⁺` 两道（且两道都要 `≥δΛ`），加 `IN-11` 的成员性 `(b′)` 排除其余道。**

| 需要什么 | 内容 | 状态 |
|---|---|---|
| **(i) 成员性 `(b′)`** | `δ` 之上无 `0⁻⁺`/奇 C/自旋≥3 的 bulk 态（`IN-B.md:12,23`；`conditional-theorem.md:207`） | 强耦合 `[T]`（`PROOF-IN11.md:19` `This proves (b′) at strong coupling`）；`β≈2.5` 侧 `[C]` 物理输入 |
| **(ii) `0⁺⁺` 道底 `>δ`** | `Δ_{0⁺⁺} ≥ δ`（残数 `I1` 已有 `[C]`，隔离 `I2` 开） | `[C]/[O]`（`V03-I2GAP.md:127–129,152`） |
| **(iii) `2⁺⁺` 道底 `>δ`** | `Δ_{2⁺⁺} ≥ δ`：**无残数、无 `I2`、无定理**，与 `0⁺⁺` 简并 ⇒ 与 (ii) **同等承重** | `[O]`（`Y02-CHANNELS.md:53,95–97`；`BLAST-EXCLUDE-BUMP.md:40` 逐字 `relocates the demand to H1 (2⁺⁺ isolation)`） |
| (iv) 其它道（`0⁻⁺`、`1⁻⁻`、…） | **不需要各自的证书**——只需要 (i)；`0⁻⁺` 的下界不可证（`BREAK-IN-11.md:3–4`），`1⁻⁻` 只有统一界（§1.1）、通道专有陈述 0 条 | `[T]`（强耦合）`/ [C]` |

**精确形式**：`m>0` 需要的不是「一个道」，是 `min(Δ_{0⁺⁺},Δ_{2⁺⁺}) ≥ δΛ` + `(b′)`。逐字可引的桥是 `V03-I2GAP.md:149`：`I2 ∧ (m₁²>0) ∧ IN-11 ⇒ m>0（全道）`——**但它的 `IN-11` 必须按简并读**：`IN-11` 不是「最轻在 `0⁺⁺`」，而是「最轻在 `{0⁺⁺,2⁺⁺}`」（`conditional-theorem.md:205,207` 逐字 `the lightest bulk state lies in 0⁺⁺ or 2⁺⁺`）。

---

## 5 · 判死线逐条

- **(a)「`0⁺⁺` 严格最轻 ⇒ `m>0` 只需 `0⁺⁺`（`I2` 就够）」——前提假，`ch05:468` 对。**
  - 逻辑形式**有效**：若 `0⁺⁺` 是全谱严格最轻且 `Δ_{0⁺⁺}≥δ`，则 `inf spec=Δ_{0⁺⁺}≥δ`（其余道更重）。
  - **前提假**：`PROOF-IN11.md:3` 逐字 `the two lightest channels are exactly degenerate at leading order … the lightest object is the degenerate pair`；`:23` 同。语料**没有一处**把 `0⁺⁺` 的严格最轻写成定理/条件（`campaign-final.md:31` 明写无对称论证能锁）。
  - **`ch05:468` 的读法**：原文 `one physical input (IN-11, the lightest bulk state in $0^{++}/2^{++}$)` ⇒ 语料自己把 `2⁺⁺` 列为可能的最轻道。`V03-I2GAP.md:110` 已下同一判：`语料自己承认最轻态可能在 2⁺⁺——那正是 I2 管不到的道`。**故 `(a)` 撤回，(b) 立。**
- **(b)「简并或 `2⁺⁺` 更轻 ⇒ 至少要 `0⁺⁺`+`2⁺⁺`」——成立，逐字记：**
  > `m>0 ⟺ inf spec(H_a|_{Ω^⊥}) ≥ δΛ`（`ch05:455–456`）需要 **`inf_c Δ_c ≥ δΛ`**，而 `inf_c Δ_c = min(Δ_{0⁺⁺}, Δ_{2⁺⁺})`（`PROOF-IN11.md:17` `[T]`）；两位因子**领头阶简并**（`:23`）且排序/隔离是**无对称可定的 `O(1)` 次领头**（`:31`）。⇒ **`m>0` 至少需要 `0⁺⁺` 与 `2⁺⁺` 两个道的底都 `>δ`**；证书侧即 `cert = min(cert_{0⁺⁺}, cert_{2⁺⁺}) > 0`，**要求两条都正**（`IN-10-attack.md:7,16`）。`2⁺⁺` 一侧的缺件恰被语料标为 `H1`（`BLAST-EXCLUDE-BUMP.md:40`；`PROOF-HAMILTONIAN-VARIATIONAL.md:13` 逐字 `If instead a lighter 2⁺⁺ exists (leading-order degenerate, PROOF-IN11), B is set against 2m_{2⁺⁺}, itself uncontrolled`）。
- **(c)「质量序在语料里无定论」——部分成立，但不是全称，且无定论处已定位。**
  - **有定论**：`β<β₁` 强耦合（`[T]`）：`{0⁺⁺,2⁺⁺}` 简并最轻，`0⁻⁺` 及 C-奇/高自旋 ≫；`1⁻⁻` 连 `a_c` 都未给。
  - **无定论**：① `0⁺⁺` vs `2⁺⁺` 的**次领头序**（无对称、`β<β₁` 内未算，`PROOF-IN11.md:23`）；② `β≈2.5` 侧整条序（`[J]` 外推，`PROOF-IN11.md:27`）；③ `1⁻⁻` 的**通道专有**位置（0 条）——但它的**位次**已被 `PROOF-IN11.md:11` 的 "every other `J^{PC}` / C-odd ≥3 gluons" **统一界**覆盖 ⟨本文件推断 `[C]`⟩，故不是"无定论"，只是"未单列"。
  - **判它需要什么**：`(α)` `0⁺⁺/2⁺⁺` 次领头 `O(1)` 分裂的有限阶强耦合计算 + OS 余项（`PROOF-IN11.md:23`，`[O]`）；`(β)` `1⁻⁻` 的**显式** `a_c` 与残数（语料未做；`Y02-CHANNELS.md:55` 逐字"无"）；`(γ)` 跨 `β₁→2.5` 的搬运（= `IN-3` 墙，`ch05:469–471`，`[O]`）。**注意**：即使拿到 `(α)`，`m>0` 仍要 (ii)+(iii)，因为**简并是领头阶、不是全部**——只要两道的底可以独立地掉到 0，`min` 就必须两条都控。

---

## 6 · 一句话总表

| 问题 | 判 |
|---|---|
| 四个道的序（强耦合 `[T]`） | `0⁺⁺ = 2⁺⁺ ≺ 0⁻⁺ ≺ 1⁻⁻`（后两级由 `a_c≥2` / C-奇 `≥3` 胶子给出） |
| 四个道的序（`β≈2.5`，`[J]`） | `0⁺⁺ 1710 < 2⁺⁺ 2390 < 0⁻⁺ 2560 MeV`；`1⁻⁻` 语料无值 |
| `a_c=1` | 单 plaquette / 单位时间的**最小 plaquette 计数**（`m_c=a_c\ln(1/u)+O(1)` 的系数），**不是**临界 `β` |
| `m₂/m₀→1` 的极限 | 强耦合 `β→0`（`u→0`）**领头阶**；同一四壁管；`O(1)` 余项未压 ⇒ 序**符号未定** |
| `m>0` 需要哪几个道 | **`0⁺⁺` + `2⁺⁺`**（`min(Δ₀,Δ₂)≥δΛ`）**+ (b′) 成员性**；`0⁻⁺/1⁻⁻` 不需要各自证书 |
| (a)/(b)/(c) | **(a) 撤回**（前提假，`ch05:468` 对）· **(b) 立** · **(c) 部分**（仅次领头序 / `β≈2.5` / `1⁻⁻` 无定论） |

---

## 7 · 检索与工具记账（`AGENTS.md §1`）

| # | 检索式 / 动作 | 工具 | 状态 | 命中 |
|---|---|---|---|---|
| 1 | `质量序\|mass order\|mass ordering\|m_?2/m_?0\|a_c`（全仓） | `grep` | OK | 6070（显示 250，全量落 spill 文件） |
| 2 | `2⁺⁺\|2\+\+\|0⁻⁺\|1⁻⁻\|inf spec\|最轻态`（`*.md`） | `grep` | OK | 327（显示 250） |
| 3 | `*ch05*`、`*y02*`（全仓文件名） | `pwsh Get-ChildItem` | OK | `docs/paper-v2/ch05-theorems.tex`、`research/v3/species/Y02-CHANNELS.md` |
| 4 | `ch05-theorems.tex` 440–489 逐行读 | `read` | OK | 455–456、468 逐字核到 |
| 5 | `PROOF-IN11.md` 全文（35 行） | `read` | OK | `a_c`、简并、`β₁` 全部逐字核到 |
| 6 | `V03-I2GAP.md` 80–169、`Y02-CHANNELS.md` 25–104 | `read` | OK | `ch05:468` 读法、逐道有无 `f₁²/I2` |
| 7 | `u\(β\)\|c_f\(β\)\|1/48`（`research/v3/species/*.md`） | `grep` | OK | 356（显示 250）；`u` 定义 = `U01-F1AUDIT.md:111` |
| 8 | `read D:\Dev\math\w7`（当文件读目录） | `read` | **操作错**（`not a regular file`，非 `TOOL-FAIL`） | 改 `pwsh Get-ChildItem` → OK，112 文件 |
| 9 | `**/[Ww]95[Cc]*` | `glob` | OK | **0**（无命名冲突） |
| 10 | `py -3.12` + numpy/scipy | `pwsh` | OK | 1.18.1 / 2.4.6 |

**未查到的条目（按 §1 记法）**：
- `1⁻⁻` 的**通道专有**质量数值/序：**未查到**（检索式 `1⁻⁻`，工具 OK，命中 27 条中 0 条给序；`0⁻⁻`、`3⁻⁻` 全库 0 命中，`Y02-CHANNELS.md:36–37`）。**但**它在统一界里被覆盖：`PROOF-IN11.md:11` `every other J^{PC} needs a_c ≥ 2` ⇒ 位次由该句给出 ⟨本文件推断 `[C]`⟩。
- `0⁺⁺/2⁺⁺` 的**次领头**序：**未查到**（检索式 `m_2/m_0`、`subleading splitting`，工具 OK，`PROOF-IN11.md:23` 明写"未做"）。
- 未找到任何一处把 `0⁺⁺` 的**严格最轻/隔离**写成 `[T]` 或 `[C]`（检索式 `0\^\{?\+\+\}?.{0,80}(strictly lighter|isolated|唯一最轻)`，工具 OK，命中 0）⇒ `(a)` 的前提在语料里**无支撑**，不是"未检索到"。

---

## 8 · 复现

```
py -3.12 D:\Dev\math\w7\w95c_order.py  >  D:\Dev\math\w7\w95c_order.log
```
`w95c_order.py` 三节：`§1` 用 `u(β)=I₂/I₁` 复算 `u(1/48)=5.208e-3`（对齐 `T23-TEMPLATE.md:64` `[V]`）；`§2` 打印 `m₂/m₀` 与 `(m_c−m₀)` 随 `β→0` 的两条走向（简并慢收敛 vs `a_c=2` 发散）；`§3` 打印 `min(Δ₀,Δ₂)` 网格与两个 `Δ=0` 反例见证。**`[E]` 数值，不增定理。**

*Refs（内部）：`docs/paper-v2/ch05-theorems.tex:455–456,468`；`docs/conditional-theorem.md:202,205–209,280,303`；`docs/campaign-final.md:31`；`research/unconditional/PROOF-IN11.md:3,9,11,13,15–19,23,27,31,33`；`research/unconditional/IN-10-attack.md:7,16`；`research/unconditional/C2-bypass-INb.md:13`；`research/unconditional/BREAK-IN-11.md:3–4,20,27`；`research/unconditional/BLAST-EXCLUDE-BUMP.md:37–40`；`research/unconditional/PROOF-HAMILTONIAN-VARIATIONAL.md:13`；`research/unconditional/P3-AUDIT-CORRECTED.md:128`；`research/v3/species/Y02-CHANNELS.md:33,43–44,52–55,69,95–97`；`research/v3/species/V03-I2GAP.md:110,120,127–129,149,152–154`；`research/v3/species/KT-RHO2.md:65`；`research/v3/species/R6-ATOM-ADVERSARIAL.md:300`；`research/v3/species/T23-TEMPLATE.md:63–64`；`research/v3/species/U01-F1AUDIT.md:111`。外部：Morningstar–Peardon, *PRD* **60** (1999) 034509（`hep-lat/9901004`）——`[J]` 格点值，读到的是转引（未逐字核原表）。*
