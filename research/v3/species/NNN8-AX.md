# NNN8-AX — `T16` 的 axiom 账本逐条判：3 → 2 → 0 能不能走通

> 轮 NNN8 · 子代理 `nnn8-ax` · `D:\Dev\math`。写域仅本文件；既有文件一字未改；未 `commit`；零新 MC；**零 Lean 编译**。进展单位 = 轮。
> 口径：只 `[T]`/`[NG]`；无百分数；每结论带「下一手」。`AGENTS §1` 四要素见 §4。

## §0 判词

| # | 断言 | 判 | 依据 | 下一手 |
|---|---|---|---|---|
| **0.1** | ★ **仓内不存在 `#print axioms hy6_clay_of_three` 的输出**：`hy6_*`/`IotaL2Seam`/`L3Bound`/`BFloor`/`ClayStatement` 在**任何 `.lean`** 里 **0 命中** ⇒ 该常量不存在，Lean 侧 `#print axioms` 只会报 unknown constant。「3 条」是**记账数**，不是机检结果 | **`[T]`**（记账） | §1；§4 #3 | 若要真输出，先落 D0 与 `ClayStatement` 定义 |
| **0.2** | 记账数**自相矛盾**：新账 = **3**（`KK1:16`，`QQ12:18` 承）；旧图 = **2**（`position-current:186–187`、`HZX:14` 的「2 条」）；题面第三条名 `hy6_N4_b_floor` 已由 `KK1:0.4` 退役（`R5` 代替）⇒ **三处指代不同物** | **`[T]`**（记账） | §1 | 以 `KK1:0.5` 的 3 条为准，旧图作废；交 `lead` 回填（`AGENTS §6`） |
| **0.3** | **xx3-seam 的修形未落地**：`XX3-SEAM` 全文**无任何文件引用**（`grep` 2 命中，含自身 2 处）；其 §4 自记 `IotaL2Seam\|HaUI\|L2Dom` @ `*.lean` = **0** ⇒「`axiom` 退场」是**条件式纸面**，账本一条未动 | **`[T]`** | §2 A1；§4 #2 | 落 D0 三定义（`L2Dom`/`HaUI`/`HlimUI`）+ 4 行 body |
| **0.4** | 三条里只有 **A1 可降级**（且是条件式、要 `T16` 换签名）；A2 = `L3` 本体，A3 = `GAP` 侧的 β-模数 ⇒ **两条不可降**，它们就是缺件的名字 | **`[T]`**（判） | §2 | 见 §2 各行 |
| **0.5** | 「最终恰好 = `L3`∧`GAP`∧`NT` 三件」：**否**。三条 `axiom` 与三件**不同位**（`NT` 在 Lean 台账里零槽，被未定义的 `ClayStatement` 吞掉）；且 `Clay ⟸ AX∧NT∧GAP∧[Wilson]`（`AA7:5` 逐字），`AX` 的开义务最小集是 `{L3,GAP}`（`AA7:0.2`）⇒ 数量 3=3 只是巧合 | **`[NG]`** | §3 | 把 `ClayStatement` 的纸面框（`KK2:6.1` 五项）拆成 Lean 槽，`NT` 才有位置 |

## §1 当前「输出」的逐字出处（无一条是机检）

| # | file:line | 逐字 | 性质 |
|---|---|---|---|
| 1 | `research/v3/species/KK1-LEAN.md:16` | 「仓内 `axiom` 台账 = **3 条**：`N1`(C，`HV5:40`)·`N2`(A，本条)·`N4`(B，本条 R5 形)。`HZX:0.4` 的「2 条」按此登记」 | **记账（最新）** |
| 2 | `research/v3/species/QQ12-LEAN.md:18` | 「`#print axioms hy6_clay_of_three` 现状 **3** 条（`KK1:0.5`）；**只改名 ⇒ 仍 3 条**」 | 转引 `KK1` |
| 3 | `research/v3/species/HV5-LEAN.md:70` | 「`T16` **一旦引用**，`#print axioms hy6_clay_of_three` **会**列出 `hy6_N1_*`/`hy6_N4_*`」 | **未来式** ⇒ 从未跑过 |
| 4 | `docs/position-current.html:186–187` | 只列 2 条：`hy6_N1_iota_l2_seam`（柱 1 出口）、`hy6_N4_b_floor`（柱 2 出口） | 旧图（2 条） |
| 5 | `research/v3/species/HY6-ASSEMBLE.md:101` | `theorem clay_of_three (hC : IotaL2Seam) (hA : L3Bound) (hB : BFloor) : ClayStatement` | 纸面签名 |

## §2 三条 axiom 逐条判

| # | axiom（逐字类型） | file:line | 能否 → `theorem` | 需要什么输入 | 判 | 下一手 |
|---|---|---|---|---|---|---|
| **A1** | `hy6_N1_iota_l2_seam : IotaL2Seam μ ι Ha Hlim`（`memLp` + `strong` 两字段） | `HV5:34–40`；修形 `XX3:33–38` | **能，仅条件式**：`theorem hy6_N1 (hL3 : L3Bound varpi a0) (hd : L2Dom …) (hW1 : …) : IotaL2Seam …`（`XX3:59`） | ① `L3`（= A2，账内）；② `L2Dom`；③ `W1`；④ `ι_a→ι₀` 强收敛；⑤ Laplace 1 行；⑥ **`T16` 换签名**（`Ha`/`Hlim` 不再是输入） | **`[T]`**（条件式）/ **`[NG]`**（无条件） | 落 D0 三定义；把 ⑥ 与 `lead` 结账（`QQ12:114`） |
| **A2** | `kk1_N2_varpi_bound (ha0 : 0<a0) : L3Bound varpi a0`，`L3Bound := ∃C, 4<C ∧ ∀a∈Ioc 0 a0, varpi a ≤ C` | `KK1:47`；def `:37` | **不能**，除非证 `L3` 本身；无更弱形态（`C=4` 逐字给 `m≥∞`，`IAE-M:0.2`） | 一个具体 `C>4` + `a`-一致 `ϖ` 上界 = `L3`（= `AA7` 的 `{L3}` 件） | **`[NG]`** | 反解 `m₀(C)=2arccoth((C/4)^{1/4})`（`IAE-M:0.1`），须先有非 log 型率（`HZU-M:60`） |
| **A3** | `kk1_N4_r5_floor (eps) (N) (hN : 2≤N) : R5 eps N`（`∃β₀∀β≥β₀∀L≥2, w₀≤eps`） | `KK1:73`；def `:54` | **降的是消费面，不是账本**：`R5 ⇒ BFloor` 桥已写（`KK1:64`，`c:=w₀`），但 `R5` 本体仍是 `axiom` | β-一致模数（`inf_β` 那层量词，`HZV-R5:41`）+ 夹逼左腿（`Reduction.lean:100–124`，机检 `[T]`） | **`[NG]`** | 或把 `T16` 第 3 参直接写成 `R5`（省一层桥），或攻 β-模数 |
| — | ~~`hy6_N4_b_floor`~~ | `HV5:51` | **已退役**：`KK1:0.4` 由 `R5` 导出（`c` 不再自由） | — | — | 题面沿用旧名 ⇒ 记账须改 |

## §3 最终几条

1. **修形落地（A1 → theorem，`T16` 换签名）⇒ 3 → 2**（剩 A2、A3）。**代价**：`T16` 至少加 3 参（`L2Dom`/`W1`/`ι→ι₀`）+ Laplace 1 行 ⇒ **`axiom` 减 1、前提加 ≥4**，是**换算不是清零**。
2. **三格全闭 ⇒ 0 条**（`QQ12:0.7`）：但那要求 `L3`（A2）与 `R5`（A3）**本身被证**；此时 `ClayStatement` 仍是纸面框（`KK2:6.1` 五项）⇒ 账本归零 ≠ Clay 归零。
3. **「恰好 = `L3`∧`GAP`∧`NT`」`[NG]`**，三条硬理由：(i) 映射错位——A2 = `L3`（同位）、A3 = `GAP` 的**一根柱**（`R5⇒BFloor⇒(B)⇒Δ>0`，非 `GAP` 本身）、A1 = seam（`AX→Re` 的极限唯一，`AA7` 的 E6），**没有一条对应 `NT`**；(ii) `NT` 在 Lean 台账里**零槽**（`NT` 只活在未定义的 `ClayStatement` 里）；(iii) `Clay ⟸ AX∧NT∧GAP∧[Wilson]` 是**四件**（`AA7:5`），`AX` 的最小开义务是 `{L3,GAP}`（`AA7:0.2`）⇒ 三件的正确写法是 `{L3,GAP}` + `NT` + `[Wilson]`。
4. **下一手（按成本）**：① 只改名 ⇒ 0 成本、0 收益，账本逐字不变；② 落 D0 三定义 + seam 4 行 ⇒ 3→2（可执行）；③ 拆 `ClayStatement` 五项成 Lean 槽 ⇒ `NT` 才有位置（这是「三件」唯一合法入口）。

## §4 台账（`AGENTS §1` 四要素）

| # | 工具 | 状态 | 检索式 / 调用（逐字） | 命中 | 失败 |
|---|---|---|---|---|---|
| 1 | `grep` | **OK**（有效响应；无 provider failure／ENOENT／空响应） | ①`hy6_iota`→`hy6_N1_iota_l2_seam\|hy6_N4_b_floor` @ 全仓 ②`hy6_clay_of_three` @ 全仓 ③`#print axioms` @ 全仓 ④`IotaL2Seam\|HaUI\|L2Dom` 类 →`hy6_\|IotaL2Seam\|L3Bound\|BFloor\|ClayStatement` @ `*.lean` ⑤`XX3-SEAM\|xx3-seam` @ `species` ⑥`ClayStatement\|def Clay` @ `species/*.md` | ①**23** ②**11** ③**221**（`hy6_*` **0**）④**0**（真 0）⑤**2**（皆自身）⑥**10** | **TOOL-FAIL 0** |
| 2 | `read` | **OK**（带行号、非空） | `XX3-SEAM`(74)·`KK1-LEAN`(109)·`AA7-AXFINAL`(79)·`HV5-LEAN:25–74`·`position-current:175–199` | 5/5 | 0 |
| 3 | `glob` | **OK** | `**/*.lean`（19962 路径，落 spill）· `species/*.md`（995） | 见左 | 0 |
| 4 | `pwsh` | **OK**（exit 0，表格非空） | `Get-ChildItem research\v3\species -Filter *.md \| Sort-Object LastWriteTime -Descending \| Select -First 12` | 12 行；`XX3-SEAM.md` **不在**最新 12 ⇒ 无下游消费 | 0 |

**直读 vs 转引**：直读 = 上列 5 份 + `KK1:16,47,73`、`HV5:40,51,70`、`HY6:101`、`position-current:186–187`。转引级（只经命中行）：`KK2:6.1`、`IAE-M:0.1`、`HZU-M:60`、`HZV-R5:41`、`Reduction.lean:100–124`。
**限定**：未编译 Lean、未核 Mathlib 名；「3→0」路径按 `QQ12:0.7` 的记账口径，未独立验证 `ClayStatement` 五项是否可拆。零新 MC；既有文件零改动；未 `commit`。

*NNN8-AX · `#print axioms hy6_clay_of_three` **无机检输出**（`hy6_*` 在 `.lean` 真 0）⇒「3 条」= `KK1:0.5` 记账数；A1 可条件式降为 theorem（3→2，代价 ≥4 前提）、A2=`L3`、A3=`R5` 不可降；「恰好 `L3`∧`GAP`∧`NT`」`[NG]`：`NT` 零槽、`Clay` 实为四件。*
