# Z03-DM — `P3-DRESSED-MASS` 的 (a)(b)(d) 逐条核 + `(3.3)` 反例复核：**(a)(b)(d) 拼得出"重求和本体"（= 审计 §1.3 items 1/3/4 + item 2 的 mass 值），但拼不出 `δ`；`(3.3)` 的 `δ ≥ κ/3` 是假的（`w05` 对，且可再简化 ⇒ 真值 `δ = f_dec = O(u/a)`）⇒ `y06` 的"整环 `[T]`"不成立**

> 2026-11-25 · 子代理 `z03-dressed-mass` · 只写本文件；未 git commit。
> 逐字对象：`research/unconditional/P3-DRESSED-MASS.md`（174 行，全读）、`research/v3/species/W05-TDECAY.md`（136 行，全读）、`research/v3/species/Y06-P3DOWN.md`、`docs/YANG-MILLS-MASTER.md:42`、`research/unconditional/P3-AUDIT-CORRECTED.md:53–61`。
> 标签 `[T]/[C]/[NG]/[M]/[O]/[R]/[J]`。No QED。

---

## §0 检索执行记录（`AGENTS.md §1` 四要素）

| # | 检索式 | 工具 | 状态 | 命中 |
|---|---|---|---|---|
| 1 | `glob **/*P3-DRESSED-MASS*`、`**/*MASTER*` | `glob` | OK | 1 / 6 |
| 2 | `read` `P3-DRESSED-MASS.md` 全 174 行 | `read` | OK | 逐字取 `:13,14,16,50,51,55,64,67–68,70,78,86–88,93–95,115–116,121–124,170` |
| 3 | `grep "P3-DRESSED-MASS"` @ `.` | `grep` | OK | 97（关键 `Y06:88,99`、`LOOP-RUNBOOK:5567,5581,6090`、`W05:18,20,53,65,66`） |
| 4 | `grep "w05\|y06"` @ `*.md` | `grep` | OK | 16 |
| 5 | `read` `W05-TDECAY.md` 全 136 行 | `read` | OK | `:30`（判别标准）、`:54`、`:65–72`、`:123–125` |
| 6 | `read` `Y06-P3DOWN.md:60–129`、`MASTER:30–54`、`P3-AUDIT-CORRECTED:50–69` | `read` | OK | `:42` 逐字取出 |

**无 TOOL-FAIL。** 未逐字核：`P3-FULLREP-KP.md` 的 `K≤CN` 证明主体（只经 `P3-AUDIT-CORRECTED:57` 转引 `[T]`）；Kotecký–Preiss CMP 103 (1986) 491 原文（语料内转引）。

---

## §1 (a)(b)(d) 各是什么、各证什么（逐字 + 行号）

先钉一句**编号**：`P3-DRESSED-MASS(3.3)` = **方程 (3.3)**（`:94`），与 **§3.3**（`:97`，"is `m₁^{dressed}` a spectral-atom energy?"）**不是同一个东西**。下文 `(3.3)` 一律指方程。

### 1.1 (a) = 装饰气的 Mayer/Ursell 重求和**存在**

逐字（`:13`，裁决表行）：

> | **(a) decoration-gas resummation exists** | **`[T]`** | The decoration gas attached to the tube is an abstract polymer gas with per-attachment activity `O(u)`; Kotecký–Preiss on the **full** representation sum (Seiler: factorial domination) gives absolute convergence of the Mayer/Ursell series, uniformly in tube length `T=t/a` and volume `L`. |

条文（`:63–68`）：`(2.1)` KP 判据 `Σ_{γ∋x}|w(γ)|e^{α|γ|} ≤ α`，`α=1`；`(2.2)` `log Z_dec(B_T) = Σ_{connected Γ attached to B_T} φ^U(Γ)`，**绝对收敛**，在 `uK<1` 上解析。

**证的是什么命题**：*装饰气*（`B_T` 外、与之邻接、两两不交的连通聚合物集合）的**对数**可以写成连通簇（Ursell）绝对收敛级数。这正是 `P3-AUDIT-CORRECTED §1.3` 的 **item 1** 逐字（`:56`："existence of a resummed/exponentiated weight (need the linked-cluster/Mayer expansion of the decoration gas along the tube)"）。**依赖**：全表示 KP `(2.1)` 本身（`P3-FULLREP-KP`，审计 `:57` 判 `[T]`）。

### 1.2 (b) = 有效质量存在、`L`-一致

逐字（`:14`）：

> | **(b) effective mass exists, `L`-uniform** | **`[T]`** | `\log Z_dec = φ_dec·T + O(1) = f_dec·t + O(1)` with `φ_dec = O(u)`, `f_dec=O(u/a)` and `|defect| ≤ C u`; existence of the free-energy density is the standard KP thermodynamic-limit theorem for the (signed) pressure; `L`-dependence enters only through wrapped decoration polymers of weight `O(u^{L/a})=e^{−cL}`, hence exponentially uniform. The rate is isolated from the next level (`≥1.5κ`) by `δ ≥ κ/3`. |

条文：`(2.3)`（`:70`）`φ_dec := lim_T (1/T)log Z_dec = O(u)`、`f_dec = φ_dec/a`；`(2.4)`（`:74`）`φ_dec = C_dec u + O(u²)`；`(2.5)`（`:78`）`f_dec(L) = f_dec(∞) + O(e^{−cL/a})`；`(3.2)`（`:88`）`m₁^{dressed} := κ − f_dec`。

**证的是什么命题**：自由能密度极限存在（KP 压力定理）、`f_dec = O(u/a)`、`m₁^{dr} = κ − f_dec`、以及空间环绕装饰 `e^{−cL}` ⇒ **`L`-一致**。= 审计 **item 2 的 mass 值** + **item 3**（volume-uniformity）。

⚠ **但该行后半句把 §3.2 的 `(3.3)` 塞进了 (b)**（"The rate is isolated from the next level (`≥1.5κ`) by `δ ≥ κ/3`"）。这句**不是** (b) 的证明内容，是 `(3.3)` 的复述；它使 (b) 的 `[T]` 标签携带一条假命题（见 §3）。另注：该行的 `|defect| ≤ Cu` 其实出自 (d) 的 `(4.1)`（`:116`），交叉引用。

### 1.3 (d) = 打扮幅度的**均匀正下界**

逐字（`:16`）：

> | **(d) uniform positive amplitude floor** | **`[T]`** | Decorations are `O(u)`-small in **absolute value** (`|D_fin−1| ≤ Cu` by absolute convergence), so `A = a₀(1±Cu) ≥ (1−Cu)a₀ > 0`; `a₀>0` uses only `c_f>0` and the singlet recoupling `f⊗f̄→1` (multiplicity `1`). **No blanket `c_r≥0` is used.** |

条文：`(4.1)`（`:116`）`|log D_fin(T)| ≤ C u`（只用**绝对收敛**，不用逐项非负）；`(4.2)`（`:122`）`A^{dressed} = a₀D_fin ≥ a₀(1−Cu) ≥ ½a₀ ≥ ½(2N)^{-2} > 0`；`(4.3)`（`:124`）`inf_L A^{dressed}(L) ≥ ¼(2N)^{-2} > 0`。

**证的是什么命题**：`log D_fin` 的**双边界** ⇒ `D_fin = e^{O(u)}` 的**真下界** `1−Cu`（审计算术要害：`|1+r| ≤ 1+O(u)` 只给上界，(4.1) 给的是 log 的双边 ⇒ 下界合法），加上 `a₀>0`（Fact 0.1 重耦合多重数 + `c_f>0`）⇒ 正地板。= 审计 **item 4**，逐字兑现。**且这是全文件唯一写出"跨切割簇绝对和 ≤ Cu"的一处**——即 §2 所说"换序"的实际落点。

**域**（三条共用）：`β<β₁`、固定 `a`、柱面 `T³_L ×_a ℤ`、固定有限 `L`、`t→∞` 先于 `L→∞`。**(c)**（`:15`）是另一件事，`[T] mod [KL]`，不在本题三条内。

---

## §2 重求和需要什么？(a)(b)(d) 里有吗？

把单管权重 `u^{4T}=e^{−κt}` 变成 `A^{dr}e^{−m₁^{dr}t}(1+O(u))`，需要四件：

| 需要 | 在哪 | 在 (a)(b)(d) 里？ |
|---|---|---|
| **(i) 精确分解** `C_L(t)=a₀u^{4T}Z_dec+R^{hi}` (1.2) | `:50`（§1.3 骨架引理） | **不在**——第 4 件，`y06` 未单列 |
| **(ii) 绝对收敛** (Mayer/Ursell) (2.1)(2.2) | `:64,:67` | **在 (a)** |
| **(iii) log↔线性 的换序**（压力极限 `log Z_dec = φ_dec T+O(1)`） | (2.2)→(2.3) `:68,:70` | **在 (a)+(b)** |
| **(iv) 缺陷的双边界** `|log D_fin| ≤ Cu` ⇒ `D_fin→D_∞=1+O(u)` | (4.1) `:116` | **在 (d)** |

**⇒ (a)(b)(d) 覆盖 (ii)(iii)(iv)：收敛性与交换求和齐备，可拼。拼法（5 步，逐条落行号）：**

1. `:50` (1.2)：`C_L(t) = a₀u^{4T}Z_dec(B_T) + R^{hi}_L(t)`（精确恒等，`w` 在不相交片上乘性）。
2. `:67–68` (2.2)：`log Z_dec = Σ_{Γ} φ^U(Γ)` **绝对收敛** ⇒ 可取 `exp`，`log` 与 `Σ` 可换（KP，`uK<1`，与 `L,T` 无关）。
3. `:70–71` (2.3)：KP 压力定理 ⇒ `log Z_dec = φ_dec T + O(1)`，`φ_dec=O(u)`。**这一步就是"单管权重 → `e^{−φ_dec T}`"的全部内容**：`u^{4T}Z_dec = e^{−κt}e^{φ_dec t}D_fin = e^{−(κ−f_dec)t}D_fin`。
4. `:88` (3.2)：`m₁^{dr} := κ − f_dec = [4ln(1/u) − C_dec u + O(u²)]/a`。
5. `:116–124` (4.1)(4.2)(4.3)：`|log D_fin| ≤ Cu` ⇒ `A^{dr}=a₀D_∞ ≥ ½a₀ > 0`，`L`-一致。

**⇒ 结论（对"重求和"这个词）：能拼，且这就是 `MASTER:42` 的"dressed-mass 重求和未证"所指的 item 1/3/4 + item 2 的 mass 值。** 标签 `[T]`（在 (i) 已给、`(2.1)` 全表示 KP `[T]` 的前提下）。

**⇒ 但 (a)(b)(d) 里没有 `δ`。** 余项的**率**只在 §3.2/`(3.3)`（`:93–95`）——`y06` 的"整环"把它算进去了，`(b)` 的行也把它写进去了。下一节核它。

---

## §3 ★ `P3-DRESSED-MASS(3.3)` 逐字 + `w05` 反例复核

### 3.1 逐字（`:93–95`）

> **3.2 Rate isolation from the next level. `[T]`** The remainder `R^{hi}_L(t)` in (1.2) is `𝔓_{\text{rest}}`, whose lowest rate is the perimeter-6 wider tube (`|P|≥6T`, rate `6\ln(1/u)/a=1.5κ`), décor-renormalized to `1.5κ−O(u/a)`; or two tubes / two-particle at `2κ` (rate `2κ−O(u/a)`). Hence, resumming decorations in `𝔓_{\text{rest}}` exactly as above,
> `C_L(t) = A^{dressed}e^{−m₁^{dressed}t}\bigl(1+O(u)\bigr) + O\bigl(e^{−m₁^{dressed}t}e^{−δ t}\bigr)`, `δ ≥ \tfrac12κ − O(u/a) ≥ \tfrac13κ\ (u \text{ small})`, `(3.3)`
> **uniformly in `L`.** I.e. `e^{m₁^{dressed}t}C_L(t) → A^{dressed}` with an exponentially small error at rate `δ≥κ/3`. `[T]`

它**断言**：`e^{m₁^{dr}t}C_L(t)` 收敛到 `A^{dr}`，且**误差的率** `δ` 有一个 `κ/3` 的**与分割无关的正下界**。算术来源：周长 6 宽管 `6T ⇒ 1.5κ`，减 `m₁^{dr}≈κ` ⇒ `δ≈½κ`。

**关键前提（两处，均标 `[T]`）**：

- `:51`：「and every `P∈𝔓_{\text{rest}}` satisfies `|P| ≥ 6T + O(1)` (perimeter-6 wider tube, the lowest non-tube) or `|P| ≥ 8T+O(1)`」
- `:55`（§1.4）：「"long" decorations (`ℓ` comparable to `T`) carry cost `≳ T`, i.e. belong to `𝔓_{\text{rest}}` (a parallel-tube / two-tube object, rate `≥ 2κ`), **not** to the leading class」

### 3.2 `w05` 的反例（`W05-TDECAY.md:66` 逐字）

> 2. **但"危险族"不是多管，是单管 + 亚线性长装饰**：`P3-DRESSED-MASS:55` 逐字只说"`ℓ` comparable to `T` 的装饰属 `𝔓_rest`（rate `≥2κ`）"，并断言 `:51`"every `P∈𝔓_rest` satisfies `|P| ≥ 6T+O(1)`"。**这两句只对"整管型"装饰成立**：一条**单壁**平行装饰（每片 1 个 plaquette、`dr=ℓ`）给 `|P| = 4T + ℓ`，`ℓ` 可一直取到 `T` ⇒ 率 `κ(1+ℓ/T)`，`ℓ/T→0` ⇒ **`𝔓_rest` 的率下确界 = `κ`，不是 `1.5κ`**。⇒ `P3-DRESSED-MASS:94` 的 `δ ≥ κ/3` **不成立**（勘误级，见 §7 E1）。

### 3.3 我的裁决：**成立**，而且比 `w05` 说的更简单（不需要 `ℓ/T→0`）

按文件**自己的归类**（`:55` 把 `ℓ` comparable to `T` 的长装饰推入 `𝔓_rest`），取 `ℓ=T`（最不利、且完全在 `:55` 的"comparable to T"内）：

- `|P| = 4T + T = 5T` ⇒ 率 `5\ln(1/u)/a = \tfrac54κ = 1.25κ`；
- `m₁^{dr} = κ − f_dec` ⇒ 相对率 `δ = 1.25κ − (κ − f_dec) = \tfrac14κ + f_dec`；
- 而 `f_dec = O(u/a) ≪ κ` ⇒ **`δ = ¼κ + O(u/a) < κ/3`**。⇒ `(3.3)` 的 `δ ≥ κ/3` **假**（不需要任何极限过程，一条 `5T` 的账就够）。

`w05` 的极限版更狠，且同样成立：对 `ℓ ≤ T` 全族求和 `Σ_ℓ (Cu)^ℓ ≤ Cu`（几何），故整族给 `≤ Cu·e^{−κt}` ⇒ 相对 `e^{−m₁^{dr}t}` 的率 `δ = κ − m₁^{dr} = f_dec = C_dec u/a = O(u/a)`。⇒ **`δ` 的真值（本分解能给的最好值）就是 `f_dec`**，与 (d) 的 `Cu` **同阶同源**——这正是 `w05` 的 E2（`:124`）：`(3.3)` 的真形状是 `|η| ≤ Cu + Be^{−δt}`，那个 `Cu` 不可去。

**⇒ 连带作废**：`:51` 的 `|P| ≥ 6T+O(1)`（假：`4T+ℓ`，`ℓ≤T`）；`:55` 的 "long decorations `rate ≥ 2κ`"（假：`ℓ=T` 给 `1.25κ`）；`:93` 的 "lowest rate is 6T ⇒ 1.5κ"（假）。`(3.3)` 的 `δ ≥ ½κ` 与 `δ ≥ κ/3` 同时死。

### 3.4 一处必须写清的保留（不许编）

**字面上**，`P = B_T ∪ {单壁}` 与 `B_T` 共链 ⇒ `P ⊇ B_T` ⇒ 按 `(1.1)` 的定义它在 `𝔓_{\text{tube}}` 里（进 `Z_dec`），**不是** `𝔓_rest`。`w05` 的例子严格说必须按 `:55` 的"推入"读（文件自己把长装饰划给 `𝔓_rest`）。两种读法都杀 `(3.3)`：

- **读 A（按 `:55` 推入）**：`P∈𝔓_rest`，率 `1.25κ` ⇒ `δ=¼κ+f_dec` ⇒ `(3.3)` 式子**假**。
- **读 B（按 `(1.1)` 留在 `Z_dec`）**：`:55` 为假，且 `(3.3)` 的**证明**用的正是 `:51` 的 `𝔓_rest ⊂ {|P|≥6T}` ⇒ 写下的证明**无效**（式子或可活，但无证明）。

⇒ 无论哪种读法："`(3.3)` 标 `[T]`"这件事**不成立**。另：若还有人主张把残差类换成"缺一段 `B_T` 的位移四壁管"（`|P|=4T+O(1)`，`P⊉B_T`，字面即属 `𝔓_rest`），率 `→κ`，`δ→f_dec`，结论同样只到 `f_dec`。

---

## §4 判死线

### (a) `(a)(b)(d)` 是否足以拼出重求和？—— **只对"重求和本体"够，对 `y06` 的"整环"不够**

- **拼得出的**（§2 五步）：`log Z_dec = φ_dec T+O(1)` ⇒ `m₁^{dr}=κ−f_dec`、`A^{dr}=a₀D_∞>0`、`L`-一致。= `MASTER:42` 的"dressed-mass 重求和"所指的**全部四项**（审计 `:56–59` 的 items 1–4），其中 item 2 的"positive gap to the next level"这句**只拿到"正"（`δ=f_dec>0`），拿不到 `κ/3`**。
- **拼不出的**：`(3.3)` 的 `δ ≥ κ/3`。它**不在** (a)(b)(d) 的陈述内，只挂在 §3.2，而 §3.2 是假的（§3.3）。
- ⇒ `y06:99` 的「**dressed-mass 重求和未证 → (a)(b)(d) `[T]`；(c) ⇒ 单一命名输入 `[KL]`=`[T]`**」这句：**对象引对了，等级给高了**。整环的正确标签是 **重求和本体 `[T]`（mod `(2.1)` 全表示 KP）+ 余项率 `[C]`（须重写，真值 `f_dec=O(u/a)`）**。

### (b) `(3.3)` 是否确实不成立？—— **是（`w05` 对）**

`δ ≥ κ/3` = **`[NG]`**（假，`¼κ+f_dec < κ/3`；`ℓ≤T` 全族给 `δ=f_dec`）。`δ ≥ ½κ` 同死。
**但降级落点要说准**：死的是**"率隔离/下一能级 ≥1.5κ−O(u/a)"**（那是 `I2`/`(L1)` 的载重物，`w05:69–74`），**不是**重求和本体（items 1/3/4 + mass 值）。所以：

- **`I1`（`f₁² ≥ r₀>0`）不被 `(3.3)` 击穿**：Lemma 3.4（`:99`）只需 `C(t)e^{mt}→A>0`，即余项 `= o(e^{−m₁^{dr}t})`；`δ=f_dec>0` 就够（`e^{−κt}=e^{−m₁^{dr}t}e^{−f_dec t}`）。更省：Fact 0.1 的正性（逐项 `≥0`）直接给原子权重 `≥ A^{dr} > 0`。**但两者都要一行**：`f_dec>0`。文件 `:74` 只写 `φ_dec = C_dec u+O(u²)`，**未声明 `C_dec>0`**；正性（Fact 0.1）只给 `φ_dec ≥ 0`，若 `φ_dec≡0` 则余项率 `=0`，Lemma 3.4 的 `o(e^{−m₁^{dr}t})` 不成立 ⇒ 需补一句"`Z_dec ≢ 1`（存在非零装饰，例如管壁的位移段，其重耦合多重数 `≥1`）"。记 **`[R]`（一行，本轮可补）**。
- **但"写下的"不行**：文件里唯一的余项陈述是错的 `(3.3)`，`I1` 的现行证明走 `(3.3)+Lemma 3.4`。⇒ **现状 `I1` 只能 `[C]`**（与 `y06` 口径一致），**且这一 `[C]` 的成因是"引用了错的 `(3.3)`"，不是"重求和没证"**。落一条勘误（§5）后，`I1` 可回 `[T] mod [KL]`。
- **`I2` 的 `δ`**：本路彻底给不出 `κ` 阶（`δ=f_dec=O(u/a)` 只是"有正率"，够 `I1` 不够 `I2`）；`w05` 的 `(L1″)` 仍是空位。

### (c) `(a)(b)(d)` 是否证的是别的东西（引错）？—— **不触发**

(a)(b)(d) 证的正是审计 §1.3 的 items 1/3/4 + item 2 的 mass 值；不是引错。唯一"引错"是 `(b)` 行**多带了**一句 `(3.3)`（`:14` 末句），那是**超额引用**，不是对象错位。另：`(c)`（§3.3 谱原子）本就不在 (a)(b)(d) 内，`y06` 也没把它算进来（它单列为 `[KL]`）——这一处 `y06` 是对的。

---

## §5 路径（最小、本轮可完）

1. **勘误（改上游文件的可引用性，`[T]` 级簿记）**：`P3-DRESSED-MASS:94` 撤 `δ ≥ ½κ − O(u/a) ≥ κ/3`，改为 `δ ≥ f_dec − O(u) = O(u/a)`；或按 `w05:113` 写成 `|η| ≤ Cu + Be^{−δt}`。删 `:95` 的 "error at rate `δ≥κ/3`"。
2. **限定**：`:51` 的 `|P| ≥ 6T+O(1)` 加"**整管型**"；`:55` 的 `rate ≥ 2κ` 同加；补 `|P| = 4T+ℓ`（`ℓ≤T`）这一族的成本账与"`δ=¼κ+f_dec` 上界"。
3. **`(b)` 行删末句**（`P3-DRESSED-MASS:14` 的 "The rate is isolated … `δ ≥ κ/3`"），该句属 §3.2 且假。
4. **补一行**：`φ_dec>0`（`Z_dec ≢ 1`，正性 + 存在非零装饰）——`I1` 的原子抽取靠它，而不是靠 `δ≥κ/3`。
5. **之后**：重求和本体 `[T]`（mod 全表示 KP）+ `[KL]` ⇒ `I1` 可回 `[T] mod [KL]`；`MASTER:42` 的"dressed-mass 重求和未证"可回填（**但须同时记 `δ` 只到 `O(u/a)`**，否则又越级）。`I2` 的 `δ` 仍缺 = `w05` 的 `(L1″)`。

---

## §6 标签表

| 对象 | 陈述 | 等级（本次核后） |
|---|---|---|
| (a) (2.2) Ursell/Mayer log 绝对收敛 | 审计 item 1 | **`[T]`**（mod 全表示 KP `[T]`） |
| (b) (2.3)(2.5)(3.2) `φ_dec=O(u)`、`m₁^{dr}=κ−f_dec`、`L`-一致 | item 2 值 + item 3 | **`[T]`** |
| (d) (4.1)(4.2)(4.3) `|log D_fin|≤Cu` ⇒ `A^{dr}≥¼(2N)^{-2}` | item 4 | **`[T]`** |
| (3.3) `δ ≥ ½κ−O(u/a) ≥ κ/3` | 率隔离 | **`[NG]`**（假；真值 `δ=f_dec=O(u/a)`） |
| `:51` `𝔓_rest ⊂ {|P|≥6T}`；`:55` 长装饰 `rate≥2κ` | 成本分级 | **`[NG]`**（如陈述；仅对整管型成立） |
| 重求和本体（拼法 §2 五步） | `MASTER:42` 所指 | **`[T]`** |
| `y06` 的"dressed-mass 重求和 `[T]`（整环）" | `:99` | **`[C]`**（含假的 `(3.3)`；勘误后可回 `[T]`） |
| `I1` `f₁²≥r₀>0`（现状，依 `(3.3)`） | — | **`[C]`**；勘误后 **`[T] mod [KL]`** |
| `I2` 的 `δ` | — | **`[R]`**（`(L1″)`，本路给不出 `κ` 阶） |
| `w05` 的反例 | `:66` | **成立 `[T]`**（结论；例子需按 `:55` 的"推入"读，见 §3.4） |

*Refs.* `research/unconditional/P3-DRESSED-MASS.md:13,14,16,50,51,55,64,67–68,70,74,78,86–88,93–95,99–101,115–116,119,121–124,170`；`research/unconditional/P3-AUDIT-CORRECTED.md:53–61`；`docs/YANG-MILLS-MASTER.md:42`；`research/v3/species/W05-TDECAY.md:30,54,65–72,113,123–125`；`research/v3/species/Y06-P3DOWN.md:88,99,131`；`research/v3/species/W02-L1.md:65,81,104`；`LOOP-RUNBOOK.md:5567,5581,6090`；`SUB2-I2-STRONG.md:15`（`Cu+Be^{−δt}` ⇒ 无 `δ`）。外部（仅转引，未逐字核）：Kotecký–Preiss CMP **103** (1986) 491；Osterwalder–Seiler Ann. Phys. **110** (1978) 440；Seiler LNP **159** (1982)；Münster PLB **95** (1980) 59 / NPB **190** (1981) 439；Gupta arXiv:hep-lat/9807028 §14.2.2。
