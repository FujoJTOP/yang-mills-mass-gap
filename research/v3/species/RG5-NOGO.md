# RG5-NOGO — 核 `l4a` 的关键一句：`C2:373/376` 是否真的把"体积增长内在"钉成 `[NG] refuted`

> 标签（charter §2）：**[T]** 已证定理 · **[C]** 条件定理 · **[NG]** 对**指名类别**的 no-go · **[O]** 开放/未排除 · **[V]** 书目核实 · **[M]** 未逐字核。
> 本轮任务：判定构造性 RG 是否**在原理上**被排除。只读，不改其他文件。

---

## §0 工具状态与检索（AGENTS.md §1 四项）

| 项 | 值 |
|---|---|
| 工具状态 | `pwsh` = **OK**（有效响应，无 `ENOENT`）；`glob` = **OK**；`read` = **OK**；`Select-String` = **OK**。无 `TOOL-FAIL`。 |
| 检索式 1 | `Get-ChildItem -Recurse research -Include *.md \| Select-String "intrinsic" -List` → 命中文件列表（30+），用于定位 |
| 检索式 2 | `glob "**/C2*.md"`（path=`research`）→ 命中 **3**：`research/unconditional/C2-bypass-INb.md`、`research/v3/C2-BALABAN-UNIFORM.md`、`research/v3/species/C2-MEASURED-SLOPE.md` |
| 检索式 3 | 逐个计数 `intrinsic` 与行数 → `C2-bypass-INb.md` 32 行 / `intrinsic`=0；`C2-MEASURED-SLOPE.md` 104 行 / 0；**`C2-BALABAN-UNIFORM.md` 411 行 / `intrinsic`=5** ⇒ **`C2` 唯一候选**（`:373`、`:376` 存在且可读） |
| 检索式 4 | `Select-String research\v3\C2-BALABAN-UNIFORM.md -Pattern "intrinsic\|no-go\|volume-free\|Lemma B\|c\*"` → **32** 行命中（复算确认） |
| 检索式 5 | `Select-String research\v3\species\L4A-RG.md -Pattern "intrinsic\|C2\|no-go"` → **24** 行命中（复算确认） |
| 检索式 6 | `glob "**/*[ClL]4[aA]*.md"` → 命中 **1**：`research/v3/species/L4A-RG.md` |

**`C2` = `research/v3/C2-BALABAN-UNIFORM.md`（411 行）。** 结论：`l4a` 指的 `C2` 唯一、可判定。

---

## §1 `C2:373` 与 `C2:376` 逐字原文

**`C2-BALABAN-UNIFORM.md:373`**（Status ledger 表行，整行逐字）：

```
| "volume growth intrinsic to the method" | **not supported** — entropy channel volume-free | **[NG] refuted** |
```

**`C2-BALABAN-UNIFORM.md:376`**（整行逐字）：

```
**Sharpest line.** The (a,Λ)-uniform Balaban step is **alive**: its combinatorial half — the large-field
```

**另需一行才完整** —— `:380` 逐字：

```
`(★)` with `(★dec)`, identical to the campaign's L3 — the UV face of the one species. **No intrinsic no-go.**
```

### §1.1 对 `l4a` 引用的精度核（非实质错误，不触发判死线 (c)）

`l4a`（`research/v3/species/L4A-RG.md:19–20`）逐字写：

> 证据：`C2-BALABAN-UNIFORM.md:376` 逐字 —— *"The (a,Λ)-uniform Balaban step is **alive** … **No intrinsic no-go**."*；
> `:373` 逐字 —— "volume growth intrinsic to the method" = **`[NG] refuted`**（最自然的组合学 no-go **已被反驳**）。

裁决：**实质正确，行号跨了 4 行**。`:376` 承载 `alive`；*"No intrinsic no-go"* 逐字在 **`:380`**（`:377–379` 是中间的解释句）。`l4a` 用 `…` 省略号合并了 `:376→:380`，属引文拼接；**其断言的语义与 `C2` 一致**，故 **不构成 (c)**，只记一条引用精度注记：引用时应写 `:376–380`。

---

## §2 "体积增长内在"这个 no-go 的**精确陈述**

它不是一个裸句子，而是 `C2` 从文献 (N1) 反推出来的**两条机制**之一直觉（`:89–96` 逐字）：

```
The literature states (N1): "all constants depend on the volume `Λ` … no control of `Λ→∞`." Two *mechanisms*
are logically possible, and they must be separated because only one is real:

**(Channel U — UV large-field entropy.)** More volume ⇒ more blocks ⇒ more ways to place large-field
regions ⇒ a union bound over blocks could carry a factor `(Λ/a)^d`. **This channel is `FALSE` as a source of
`Λ`-growth** — its per-unit-volume form is volume-independent; convergence is a constant inequality
(**Lemma B, [T]**, §2). Witness confirms: the per-unit-volume large-field density depends on the number of
blocks nowhere (§2.3).
```

**精确陈述（我按 `:373` + `:89–96` 复原，不是引用）**：

> **(INT)** 若 `(a,Λ)`-一致 Balaban 步的成功必须付出一份"随总体积 `Λ` 增长"的代价 —— 因为体积越大、大场区块越多、区块并集界可带因子 `(Λ/a)^d` —— 那么该步的常数不可能对 `Λ→∞` 一致；于是"用 Balaban 型逐壳展开构造 4D YM"这条路 **在原理上** 被封死（组合学熵爆炸 = 方法内禀的 no-go）。

`C2` 把 **(INT)** 拆成两半判：**Channel U**（UV 大场熵）与 **Channel I**（IR 塔底）。**(INT) 的全部"组合学"力量落在 Channel U 上**；Channel I 是另一种 `Λ`-依赖（`:98–106`），**不是熵**。

---

## §3 谁反驳、怎么反驳、等级

### §3.1 反驳者与反驳内容 —— `C2` §2.2 **Lemma B**

`C2:154` 逐字：`### 2.2 Lemma B (large-field entropy is volume-free per unit volume) — **[T]**`

`C2:156–170` 逐字（Lemma 陈述的核心）：

```
> **Lemma B.** Fix `d = 4`, blocking factor `L`, and shells `j = 1,…,J`. Suppose the shell fibre at scale
> `j` satisfies the Gaussian-type large-deviation bound `μ(SF_j(B)^c) ≤ exp(−c/g_j²)` for a universal
> …
>   ρ(c) := Σ_{j=1}^{J} a_j^{-d} · exp(−c/g_j²) = a_0^{-d} Σ_{j≥1} exp( j ( d ln L − 2 c β₀ ) ) ,
> ```
> which is **finite and independent of the total number of blocks** `(Λ/a_0)^d`, provided the **constant
> inequality**
> ```
>   2 c β₀ > d ln L      ⟺      c > c* := d ln L / (2β₀)
> ```
> holds. The **total** expected large-field block count is `ρ(c)·(Λ/a_0)^d` and grows like `Λ^d` — but that is
> the *extensivity of a density*, not an obstruction: an action-**density** bound needs `ρ(c) < ∞`, not the total.
```

**反驳的机制（怎么反驳的）**，三步（`C2:172–178` proof sketch）：
1. 逐壳大偏差 `μ(‖F‖ ≥ t_j) ≤ exp(−c t_j²/g_j²)`（Chernoff/Borell–TIS + Holley–Stroock 处理紧致色因子）；
2. 单位体积内标度 `j` 的区块数 `a_j^{-d} = a_0^{-d} L^{dj}`；
3. 几何求和 → 比值 `exp(d ln L − 2cβ₀)`，**收敛 ⟺ 常数不等式** `2cβ₀ > d ln L`；**整个论证不引用粗胞数 `(Λ/a_0)^d`**，只用 `a_0` 与常数。

**关键一击（Kotecký–Preiss 加强，`C2:180–184` 逐字）**：把"按区块并集界"换成"按连通聚合物求和"，`c*` 只会**变小**（"`κ` an O(1) lattice-animal growth constant, not `d ln L`"），而 **"volume-free, constant-level" 的定性结论不变**。

**结论行**（`C2:186–189` 逐字）：

```
**Consequence.** **Channel U is closed [T].** "The volume growth is intrinsic to the method" cannot be
sourced from large-field entropy. The residual `Λ`-dependence, if any, is Channel I (IR tower bottom) and the
analytic per-shell bound `(★)`. This is the *positive* structural result of this attack: it removes the most
natural candidate for an intrinsic combinatorial no-go and *redirects* the obstruction to `(★)`.
```

**证人数值见证（`C2:191–206`）**：`research/v3/w7/c2_largefield_count.py`，`d=4,L=2`：`c*=29.85 (N=2)`、`c*=19.90 (N=3)`；`ρ(c)=0.638`（`c=40,N=2`）**固定**，而总量 `Λ^d` 增长（163 → 2615 → 4.2×10⁴）。**诚实保留（同节逐字）**：粗并集界需 `c* ≈ 20–30`，物理大场权重常数是否超过此阈 **"a normalisation question I did **not** settle from primary sources"**。

### §3.2 等级裁决 —— **分裂，不可一概标 `[T]`**

| 子命题 | 等级 | 出处（逐字/行号） |
|---|---|---|
| **Channel U 关闭**：大场熵**每单位体积**与总区块数无关，(INT) 的熵来源**不成立** | **[T]**（"proof sketch + witness"） | `C2:367` 逐字 `\| **Lemma B** — large-field entropy volume-free per unit volume; const inequality `2cβ₀>d ln L` \| proof sketch + witness \| **[T]** \|`；`:186` "Channel U is closed [T]"；`:93` "This channel is `FALSE` as a source of `Λ`-growth" |
| **(INT) 作为"方法内禀 no-go"被反驳** | **[NG] refuted** | `C2:373` 逐字 `\| "volume growth intrinsic to the method" \| **not supported** — entropy channel volume-free \| **[NG] refuted** \|` |
| **剩余 `Λ`-依赖 = Channel I（IR 塔底）** | **[C]**（作者自陈"inferred"） | `C2:366` 逐字 `\| mechanism of `Λ`-growth = IR tower bottom (Channel I), not UV entropy (Channel U) \| adjudicated \| **[C]** (U settled by Lemma B; I inferred) \|`；`:108–111` 逐字 `[C]` — this adjudication is my analysis … the *mechanism* attribution U-vs-I is **inferred**` |
| **整个 (a,Λ)-一致步** | **[C]** 条件于 `(★)∧(★dec)` | `C2:371` 逐字 `\| **the uniform step** \| `[C]` conditional on `(★)∧(★dec)` \| **[C]** \|` |
| **`(★)∧(★dec)` = L3** | **[O]** | `C2:372` 逐字 `\| **`(★)∧(★dec)` = L3** \| open, one-sided variance bound \| **[O]** \|` |

**所以反驳的等级是：在"熵/组合学通道"上 `[T]`（结构性，非定量完整——`c*` 阈值未从原始文献核定）；在"这条路的剩余障碍"上 `[C]/[O]`。**

一处必须点明的**范围限制**（AGENTS.md §2.1，不得升级为普适）：`C2:353–356` 逐字划定 —— 该结果是 **UV** 的、**固定物理体积 `Λ`**、弱耦合 `g_j→0`、`a→0`，并且 *"makes no claim about `s→0` of the spectral measure `ν`; a tool that reproduces Lemma B but says nothing about `ν` near `s=0` has **not** addressed the mass gap"*。**Lemma B 反驳的是"熵型 no-go"，不是"无隙"。**

---

## §4 判：构造性 RG 在原理上是否被排除？

**判死线 (a) 成立 —— 构造性 RG 是 `[O]`，不是 `[NG]`。**

逐字写清：

1. **(INT) 这条"体积增长内在于方法"的候选 no-go 已被反驳，反驳的支柱是 `[T]`（Lemma B，`C2:367`；`Channel U is closed [T]`，`C2:186`），不是 `[C]`。** 因此"构造性 RG 因组合学熵随体积爆炸而在原理上不可行"这一排除理由 **不成立** ⇒ **构造性 RG 在原理上未被排除 = `[O]`，标签不是 `[NG]`**（`C2:380` 逐字：`**No intrinsic no-go.**`；`C2:17` 逐字：`**The (a,Λ)-uniform step is ALIVE, and it is REDUCED to a single named inequality — it is not a no-go.**`）。
2. **但 `[O]` 不等于 `[T]`/`[O]-1 已完成`。** 该步本身仍是 `[C]` 条件于 `(★)∧(★dec)`（`C2:371`），而 `(★)∧(★dec)=L3` 是 `[O]`（`C2:372`）。所以精确状态是 **`[O]`：缺一条具名不等式，不是缺一条路线**。
3. **残余的 `[C]` 项：** Channel I 归因（`C2:366`）是推断而非证明；即便它成立，Channel I 就是隙本身（`:103` "the *same wall seen from two ends*"），因此**不能被用作新的 no-go**——它恰恰是靶。
4. **(b) 不触发、(c) 不触发。** (b)：反驳并非"仅条件性"——熵通道的关闭是 `[T]`；(c)：`C2` 存在、`:373/:376` 确实是那个意思，`l4a` 仅行号跨 4 行（`:376–380`），实质无误。

### §4.1 一句话判词

> **`[NG]` 不是这条路的标签；`[O]` 才是。** 构造性 RG 被从"原理性不可能"降级/升级为"一条具名不等式 `(★)` 未证"，而**不是**被判死。**唯一仍活的路仍是活的**——本轮的判定**不改变**全项目"唯一未判死的路"这一地位，相反，它**确认**了该地位。

### §4.2 本轮不声称的东西（诚实线）

- 不声称 Lemma B 的**定量阈值** `c*≈20–30` 已从原始文献核定（`C2:201–206` 自陈未定）。
- 不声称 `(★)` 可解。
- 不声称 `[O]` ⇒ Clay 有进展（`C2:382` 逐字：`Clay 0%.`）。

---

## §5 引用清单（供上层核对）

- `research/v3/C2-BALABAN-UNIFORM.md` — `:17`、`:25–29`、`:89–96`、`:98–111`、`:154–189`、`:191–206`、`:214–219`、`:353–356`、`:366–373`、`:376–382`、`:390–399`。（**正文逐字读**，本轮已读 `:17–46`、`:87–114`、`:150–219`、`:350–411`）
- `research/v3/species/L4A-RG.md` — `:19–20`、`:32–36`、`:78–79`、`:143–144`、`:170–171`、`:261`（**逐字读**其引用句）
- 脚本证据：`research/v3/w7/c2_largefield_count.py`（提出而未重跑；`C2:193` 称 "deterministic; run quoted"）
- 文献：Balaban CMP 109/116/119/122 系列、Magnen–Rivasseau–Sénéor CMP 155(2) 325–383 — **`[V]` 二级**（`C2:390–399` 转录，原文本轮**未逐字核**，付费墙）；Holley–Stroock Ann. Probab. 15 (1987) 1541 — `[V]`（carried）。
