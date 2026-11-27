# IAB-N1 — `N1 ← L3` 的四行归约

> **轮** IAB · 子代理 `iab-n1` · `D:\Dev\math`。只新建本文件；既有文件原样；未编译 Lean；未 `git commit`；零 MC。
> **靶（题给，不复推）**：`hy6_N1_iota_l2_seam : IotaL2Seam μ ι Ha Hlim`（`HV5-LEAN:34–40`）← `L3`（`HQ1-L3:36–45`）；`HY3-C` 的 `(C-n)` 是 `[T]`（条件 L3），不用 `ι_a`、不用强收敛。
> **规矩**：只交 `[T]`/`[NG]`；不用百分数；每条结论附「下一手」；进展单位是轮。

## §0 判词

| # | 断言 | 判 | 一行依据 | 下一手 |
|---|---|---|---|---|
| **0.1** | 归约**只到 `ι_a`-free 形**：L3 供 `μ`、全 `n` 矩收敛、一致模；`IotaL2Seam` 里**显式 `ι_a` 的 `Tendsto`** 不是 L3 的推论 | `[T]`（归约账）/ `[NG]`（显式形） | §1 R1–R4；`HY3-C:32,44` | 把 `strong` 改写成矩-强形（ι_a-free） |
| **0.2** | `memLp` **的 `∀a>0` 只有在 `ι` 随 `a` 变时才有内容**：取 `ι a O := 同名嵌入`（a-无关、纤维 rank-one）⇒ 量词塌缩、`a>a₀` 段自动满足、`[T]`（条件 L3）；一般 `Ω→H` 纤维 ⇒ `[NG]` | `[T]` / `[NG]` | §1 R2；`HV5-LEAN:36` | 钉 `ι` 的一行定义（乘子-真空嵌入） |
| **0.3** | L3 供不到 `(L²-dom)`：它**同时**是 C 的缺件与 N1 显式形的缺件 ⇒ **N1 不是第二堵墙**，是同一缺件的第二个出口 | `[T]` | `HY3-C:33,47`；§2-1 | 与 `G1-UNIFORM-ZA` 两侧界核名 |
| **0.4** | **`Tendsto` 是全极限形，而 GL3 只给子列** ⇒ N1 还借 `W1`（唯一性；柱 2 的 `(B)`+GL9 买） | `[T]` | `glue-a-limits:70–75`；`HQ1-L3:190` | 核 T1 的 `liminf` 对子列形是否仍立 |
| **0.5** | `(A-abs)` 不吃 N1：走 `(C-n)`（`[T]`，条件 L3）⇒ N1 只服务 `T1`/`T14`，**降为半承重** | `[T]` | `HY3-C:15,66–81` | 位置图按此记，不按"三格并列" |

**一句话。** 四行里 R1/R2 可归约（R2 需一个定义选择），R3 是 `[T]`（模一处逐行），R4 是 `[NG]`；归约后 N1 **不是一条 axiom，是 `(L²-dom)` + `W1` + 三条 ≤1 行的具名步**。

## §1 四行归约

| 行 | 供给的字段 | L3 的哪一半 | 怎么供 | 剩什么 |
|---|---|---|---|---|
| **R1** | `μ`（`MemLp` 的第三变元；`Hlim` 的出生地） | `(U1)`+`(U3)`＝`UVS(L)` | `GL3-IMP (i)(ii)(iii)`：Banach–Alaoglu ⇒ 相对紧；可数对角 ⇒ 全 `n` 同时 weak-* 收敛；**极限满足同一界**（范数下半连续）；等度连续 ⇒ Minlos 给概率测度 `μ_L` | 唯一性 `W1`（R1 只给子列）；`L→∞` = `UVS-L`，**不属 L3**（`glue-a-limits:64`） |
| **R2** | `memLp` | `(U1)` 在 `n=2` | `abs(S^{(2)}_a(ḡ,g)) ≤ C₂‖g‖²_{C^k}`；`n=2` 的收敛 ⇒ `‖g‖²_μ = lim_a S^{(2)}_a = S^{(2)}_μ(ḡ,g) ≤ C₂‖g‖²`；取 `ι a O := g` ⇒ `MemLp` 对一切 `a>0` 由**一条**事实满足 | ①"同名/rank-one"是**定义选择**，不是定理；②`ι` 随 `a` 变 ⇒ 同一界不搬到 `μ`（要 `(L²-dom)`）；③一般纤维 = `μ⊗μ` 二重对象，L3 只供单重 |
| **R3** | `strong` 的**模**（Trotter–Kato (ii)） | `(U2)` + F2 读法 `sup_a Z_a^{-1}‖H_a‖_D<∞` | 命名同一（`HQ1-L3:58`：仓库断言 F1⟺F2） | **F1⇒F2 的逐行推导未做**（仓库只写 "equivalently"）⇒ `[NG]` |
| **R4** | `strong` 的**收敛**（Trotter–Kato (i)+(iii)） | `(U1)` 全 `n` + transfer 恒等 | `⟨g,T_a^{j a₀/a}g⟩_a → ⟨g,e^{−j a₀H}g⟩_μ`（=`HY3-C:70` 的 (n1)，L3 特例）⇒ `t∈a₀ℕ` 上弱收敛且 `‖T_a^j g‖²_a`＝`2j` 矩 convergent；+R3 的模 ⇒ 弱+模⇒强（2 行，μ_a 侧） | ①把 μ_a 侧改写成 `ι_a` 向量形 = `(L²-dom)`；②`t` 连续化（1 行未写）；③`z.im>0` 半平面（Laplace+解析延拓，转引未机检）；④全极限 ⇒ `W1` |

**四行为什么不是五行**：R1 供"对象"、R2 供"对象属于 L²"、R3 供"模"、R4 供"收敛"——`IotaL2Seam` 的两个字段恰好各吃两行（`memLp ← R1+R2`，`strong ← R3+R4`）。

## §2 归约后 N1 还需要什么（逐条）

1. **`(L²-dom)`** —— `[NG]`（承重）。`C^{-1}‖g‖_{μ_a} ≤ ‖g‖_μ ≤ C‖g‖_{μ_a}`，`C` 与 `a,g` 无关（`HY3-C:33`）。它是 R2② 与 R4① 的**共同**前件，也是 C 的唯一缺件。**下一手**：与 `G1-UNIFORM-ZA` 的两侧界核名，写成一条 `[NG]`→有名件。
2. **`W1` 唯一性** —— `[NG]`，不属 L3。R1/R4 只给子列；N1 的 `Tendsto` 是全极限。**下一手**：核 T1 的 `liminf` 对"子列-子列同极限"形是否仍立——若立，N1 可降为子列形，`W1` 不必由 N1 承。
3. **F1⇒F2 逐行** —— `[NG]`。`sup_a Z_a^{-1}‖H_a‖_D<∞` 与 `(U1)–(U3)` 的等价目前只有"命名同一"。**下一手**：用 `log T_a = −aH_a` + 二阶矩写 2 行。
4. **`t` 连续化** —— `[NG]`（未写）。`t∈a₀ℕ` → `t≥0`：`ε/3` + 半群 + 一致压缩，1 行。**下一手**：写这一行；若 T1 只用 `liminf`，先查是否需要。
5. **`z` 半平面** —— `[NG]`（形式面）。`strong` 要对所有 `0<z.im` 成立；L3 只到实时间半群。**下一手**：Laplace 表示（`Re z<0`）+ 解析延拓到 `Im z>0`（Trotter–Kato，转引）。
6. **`ι` 的物理识别**（=`HQ1-L3:191` 第 8 格是否定义格）—— `[NG]`。R2 的"同名"是选择，物理 `ι_a` 是另一件；不做也能走 `(C-n)`。**下一手**：不做；`(A-phys)` 立项时再做。
7. **`L` 固定** —— 窗口纪律：R1–R4 全部只对固定 `L`；`L→∞` 是 `UVS-L`（`[NG]`，`glue-a-limits:137` 的 F3-5 证其**无条件假**）。**下一手**：不并入本账。

**净账**：归约后 `N1 = (L²-dom) ∧ W1 ∧ [F1⇒F2] ∧ [t 连续化]`；`memLp` 面在"同名"选择下清空，`strong` 面留在 `ι_a`-free 形。`(C-n)` 已把 `(A-abs)` 从这张账上摘走。

## §3 台账

### 3.1 `AGENTS §1` 四要素

| # | 工具 | 状态 | 实际调用（逐字） | 命中 |
|---|---|---|---|---|
| 1 | `pwsh` | **OK**（exit 0） | `Get-ChildItem -Force D:\Dev\math \| Select Mode,Name`；`Get-ChildItem D:\Dev\math\research\v3\species \| Select Name,Length,LastWriteTime` | 根 **22** 项；species 全列（输出截断落盘） |
| 2 | `glob` | **OK**（无 provider failure） | ①`**/*IotaL2Seam*` @repo ②`research/v3/species/IAB*` ③`**/BREAK-IN-9.md` ④`**/glue-a-limits/REPORT.md` | ①**0**（`IotaL2Seam` 只作 `HV5-LEAN` 内 Lean 文本，非文件名）②**0**（本文件为首次创建）③**1** ④**1** |
| 3 | `grep` | **OK**（无 ENOENT / 空响应） | ①`N1\s*←\s*L3\|4 行归约\|四行归约\|L3 ⇒ memLp\|memLp.*L3` @ `species/*.md` ②`MemLp\|rp_l2_closure\|rp_transfer` @ 全仓 `*.lean` | ①**4**（`HV5-LEAN:15,67,95`；`HZX-COMPLETE:14`）②**44**（`GL4L2:54–58,86–89,210–216`；`GL4Transfer:48,97`；zenodo 副本同） |
| 4 | `read` | **OK**（带行号正文） | `HV5-LEAN`(97 全)·`HY3-C`(112 全)·`HZX-COMPLETE`(107 全)·`HY6-ASSEMBLE`(130 全)·`HQ1-L3:1–200`·`BREAK-IN-9`(22 全)·`REPORT.md:50–89,125–144`·`GL4Transfer.lean:40–64` | 8 份全返回 |
| 5 | 失败标记 | **TOOL-FAIL 0；空响应 0；PATH-MISS 1** | PATH-MISS = `research/v3/species/BREAK-IN-9.md`（`read` 报 not found）⇒ 经 ③`glob` 定位到 `research/unconditional/BREAK-IN-9.md` 后重读成功 | — |
| 6 | `web_*` | 未调用 | — | **0**（本任务不引新文献） |

### 3.2 直读 vs 转引（`AGENTS §3`）

- **直读（全文/行号级）**：`HV5-LEAN`·`HY3-C`·`HZX-COMPLETE`·`HY6-ASSEMBLE`·`HQ1-L3:1–200`·`BREAK-IN-9`(22)·`glue-a-limits/REPORT.md:50–89`·`GL4Transfer.lean:40–64`。
- **只经 `grep` 命中行采信（未开全文）**：`GL4L2.lean:54–58,86–89,210–216`（只证 `MemLp` 是**输入**、L² 闭包是 RP 传递，不是 `MemLp` 的产出）。
- **经本文件转引、本轮未开**：`G1-UNIFORM-ZA`·`G1-CONTINUUM-OS`·`GN3-SPECTRAL:97–107,187–199`·`HZ8-APHYS:52–57,63–67`·`HZG-BLOCK` ⇒ 记转引级，只作行内容证据。
- **外部文献**：零新增；Trotter–Kato 经 `BREAK-IN-9:9` 转引 Reed–Simon I Thm VIII.25（**未逐字核正文**）。

### 3.3 边界

1. 本文件**不证 `N1`、不证 `L3`**：四行是归约（蕴含骨架），`L3` 本身仍 `[NG]`（`HQ1-L3:0.2`）。
2. R2 的"同名/rank-one"是**定义选择**，写成 `[T]` 只因"在该选择下 `MemLp` 由一条 L3 事实给出"；选择未落 Lean。
3. R4 的"弱+模⇒强"是 `HY3-C:44` 的 2 行，本轮**逐字沿用**，未重证。
4. 未编译 Lean；签名未过 elaboration；未 `git commit`；零 MC；无百分数；无 QED。
5. **总下一手**：先钉 `ι` 的定义（一行，决定 `memLp` 是否清空），再把 `(L²-dom)` 与 `G1-UNIFORM-ZA` 核名——它同时销 C 与 N1 的一半。

*IAB-N1 · `N1 ← L3` = R1 对象 + R2 范数 + R3 模 + R4 收敛；剩 `(L²-dom)` + `W1` + 三条 ≤1 行。*
