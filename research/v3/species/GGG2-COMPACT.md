# GGG2-COMPACT — 「2 阶矩有界 ⇒ 子列收敛」：紧性只翻译界，不产生界

> 2026-11-27 · 本轮 · 靶：`(i) 连续极限存在 ⟸ L3（紧性）+ W1（唯一性）` 的紧性半步，与 `L3` 的标量形 `O_B ≤ C − Var(q)` 能否喂它。
> 纪律：只新建本文件；不改既有文件、不 commit、零 MC；只交 `[T]`/`[NG]`，每结论带「下一手」。

## §0 判词

| # | 断言 | 判 | 依据（`file:line`） | 下一手 |
|---|---|---|---|---|
| 0.1 | 紧性 + 可数对角 + 极限泛函这一步是 `[T]`：Banach–Alaoglu / Tychonoff + `IsCompact.tendsto_subseq`，Lean 机检（公理 `propext/Classical.choice/Quot.sound`） | `[T]` | `research/glue-a-limits/GL3Limits.lean:58–69`；`research/glue-a-limits/REPORT.md:83–84`；`research/unconditional/CONTINUUM-TIGHTNESS.md:5` | 这半步不再立项 |
| 0.2 | 缺的只是**输入界** `hS : ∀k j, abs(S k j) ≤ C j`（= `UVS(L)` 的 `(3.U1)`，n-线性泛函的范数界）。标量 `O_B≤C−Var(q)` 只住一个坐标，给不出半范数族 | `[NG]` | `REPORT.md:58–62`；`research/v3/species/HQ1-L3.md:37–45` | 把 `L3` 验收句换成 `(U1)(U2)(U3)` 三子句（改陈述，不是补估计） |
| 0.3 | **2 阶 ⇏ n 阶**（无界类）：`μ_ε=(1−ε)δ₀+εδ_{ε^{-1/2}}` 有 `E[X²]≡1`、`E[X³]=ε^{-1/2}→∞`；靶算符 `¼trF²` 正是无界类 | `[NG]` | `HQ1-L3.md:78`；反例=本文件两行算术 | 不用 2 阶推 n 阶；用 `(U1)` 的 `∀p` 或 `(U3)` |
| 0.4 | Minlos 是外部假设，且要**核空间**：`C_c^k` 的 sup 范数不核 ⇒ 检验类须换成 `S`/Sobolev，或把界换成 `S`-半范数族 | `[T]`（标准，正文未核）/ `[NG]`（仓内无证书） | `GL3Limits.lean:19–27, 82–94`；`REPORT.md:85,119` | 写清 `hMinlos` 签名（空间 + 正性 + `0` 处连续） |
| 0.5 | 阶段 2（`L→∞`）在 `(C_c^k)'` 上的字面相对紧**是假命题**：`‖S^{(2)}_L‖ ≳ χ₂L⁴` | `[T]`（对字面式） | `research/v3/species/HS-L3-5.md:15,48–54` | stage 2 改走定额化/固定支撑（`HS-L3-5:19`）或 `MG` |
| 0.6 | OS（编号按仓内：OS2=`O(4)`、OS3=RP；`GGG3` 用靶编号 OS2=RP/OS3=对称，跨文件必带号）：OS1、OS3=RP 过弱极限 `[T]`；OS0 只给定 `UVS` 后 `[T]`；OS2、OS4 与 W1 `[NG]` | `[T]`/`[NG]` | `REPORT.md:145–176,218–225`；`G1-CONTINUUM-OS.md:5,13`；`GGG3-OSCLOSED.md:17–21`；`GGG4-W1.md:43` | OS4 只由柱 2 `(B)_gap`+GL9 买；`W1` 的 `(B)+L3` 买法已判循环+扇区缺口（`GGG4:43`） |

**一句话。** 五步里紧性、对角、极限泛函全 `[T]`（Lean 已证）；OS 的**闭性半边**（OS1/OS3=RP，及给定 `UVS` 后的 OS0）也 `[T]`；两处 `[NG]` 是**输入界**（步 2）与**唯一性**（步 5），OS2/OS4 另计。故 `L3 ⇒ 子列收敛` 可给 `[T]`，但它只是「消费 L3 的紧性包装」，且只在**矩坐标 + 固定 L** 上；`S'(ℝ⁴)`/`L→∞` 上不给（`χ₂L⁴`）。→ 下一手：位置图把「连续极限」拆成「盒中 `μ_L`」与「`ℝ⁴` 极限」两行。

## §1 五步论证（固定 `L`，`a→0`；要证 `∃ a_k↓0, φ, μ_L`）

**步 1 · 空间/拓扑 `[T]`。** 环境空间 = 矩坐标族 `S_k ∈ ∏_{j∈ℕ} ℝ`，`j` 枚举 `(n, f₁,…,f_n)`（`f` 取 `C_c^k(B(0,L/4))` 的可数稠密族），拓扑 = 逐坐标收敛，即 `(C_c^k)'` 弱-* 在稠密族上的坐标化（`GL3Limits.lean:5–8`；`HQ1-L3.md:49`）。`S'(ℝ⁴)` **不是**紧性空间：`L3` 只控盒 `B(0,L/4)`，盒外零信息；`S'` 只能是 Minlos 的输出。→ 下一手：所有引用处写死「盒 + `C^k` + `k`」，否则「空间」无定义。

**步 2 · 输入界：半范数族 + 各阶矩 `[NG]`。** 需要的是 `|S^{(n)}_{L,a}(f)| ≤ C_{n,k}(L)∏‖f_i‖_{C^k}`（`(3.U1)`），即 n-线性泛函的**范数**界（`REPORT.md:58–62`）；标量 `O_B≤C−Var(q)` 给的是 `|S^{(2)}(f₀)|` 一个值，对 `sup_{‖f‖≤1}` 零信息，Banach–Alaoglu 要范数球（`HS-L3-5.md:15`）。各阶矩也不能从 2 阶得到：无界类上 2 阶不控 3 阶（`μ_ε` 反例）；正确来源是 `(U1)` 的 `∀p`（多项式算符的 n 阶矩 = `p=n` 的 `L^p` 界）**加上** `(U3)` 的 n-一致因子增长 `C_n(L) ≤ A^n n! P(L)`——后者是 Minlos 连续性的输入，严格强于 `(U1)`（`REPORT.md:91–92`）。注意 `(3.U1)` 写的是 `n≤n₀`：**有限 `n₀` 只给有限维矩问题，全 `n` 的那一半在 `(U3)`（生成函数在 `0` 连续）**，不在 `(U1)`。→ 下一手：`(U3)` 是这条链上唯一仍无任何供给的格子。

**步 3 · tightness 判据 `[T]`。** 逐 `a`：`G^{#links}` 紧 + `μ` 概率 ⇒ Prokhorov 自动、无内容。跨 `a`：**没有**公共紧空间；推前到有界类 `∏[-1,1]` 得 Tychonoff 紧、常数 = 1（空洞，`REPORT.md:77,94`）；物理类的 tightness **就是** `UVS(L)`（= `L3`），「无可分离的 tightness 捷径」（`CONTINUUM-TIGHTNESS.md:3`）。真正要的判据不是 Prokhorov 而是 Banach–Alaoglu（对偶范数球弱-* 紧）+ 可数对角；Prokhorov 的核形式（Minlos）只在末端把极限矩泛函变成测度。**逐坐标有界已使族在 `K=∏_j[-C_j,C_j]` 上 tight**（Tychonoff），但 `K` 太大——它的点不是矩泛函，缩不到测度，缩的手续正是 `(U3)` 的连续性。→ 下一手：不再把 tightness 当独立靶（它 = `L3` = 边际 BPHZ）。

**步 4 · 子列 + 极限泛函 `[T]`（条件于步 2 与 `hMinlos`）。** `S_k ∈ K` ⇒ 子列 `φ`、极限 `L∈K`（`GL3Limits.lean:58–69`）；正性/对称在逐点极限下闭；范数界过极限（弱-* 下半连续），再用稠密族把 `L` 延拓到全部 `f`（`REPORT.md:84`）——**这一步也吃步 2 的界**。`hMinlos` + 正性/归一 ⇒ `μ_L`（`GL3Limits.lean:82–94`；Mathlib 无 Minlos，`REPORT.md:119`）。→ 下一手：无（除 `hMinlos` 的签名写全）。

**步 5 · 极限的 OS 与唯一性 `[T]`/`[NG]`。** 编号按仓内（OS2=`O(4)`、OS3=RP；`GGG3-OSCLOSED.md:4` 按靶编号 OS2=RP/OS3=对称，跨文件必带号）。OS1（平移，群像稠密）与 OS3=RP 弱闭 `[T]`（GL4；Lean `rp_l2_closure_trunc`，唯一显式假设 `MemLp A 2 μ`，`REPORT.md:145–176,215`；`GGG3:18–19`）；OS0 单靠弱收敛 `[NG]`、给定 `UVS`/`L3` 后 `[T]`（`GGG3:17`）；靶算符 `¼trF²` 的 `A∈L²(μ)` 由 tightness/`UVS` 供、不由 RP 供（`GK5-CONTINUUM.md:227`）。OS2=`O(4)` `[NG]`（`cl(W)=W⊊O(4)`、各向异性 Gaussian 锐性，`GGG3:20`）；OS4/聚类 `[NG]`（massive Gaussian `m_n→0` → massless GFF，`REPORT.md:218–225`）；`W1` 全收敛无条件 `[NG]`（`REPORT.md:75`），且 `(B)+L3 ⇒ W1` 的买法**不是 `[T]`**（循环 + 带电扇区缺口，`GGG4-W1.md:43`）。→ 下一手：`W1` 与 `L3` **并联**（`HZ8-APHYS.md:215` 只是形状），不串联。

## §2 瓶颈

**判：拓扑的选择。** `tightness` 不是瓶颈——它是 `[T]`（给定界后免费；有界类常数 = 1 空洞）。矩的阶数也不是独立瓶颈——它是拓扑的**坐标影子**：同一个缺件，在矩坐标拓扑里读成「缺 `(U3)`/`n≥3`」，在 `(C_c^k)'` 范数拓扑里读成「范数球无界」，在局部代数里读成「无」（免费空洞）。

| 空间 | 紧性 | 需要的界 | 判 |
|---|---|---|---|
| `∏_j ℝ`（矩坐标） | Tychonoff | 全族 `C_j` + `(U3)` | `[T]`（输入 = `L3`） |
| `(C_c^k)'`，固定 `L` | Banach–Alaoglu | 范数球 `C_{n,k}(L)` | `[T]`（输入 = `UVS(L)`） |
| `(C_c^k)'`，`L`-一致 | 同上 | `sup_L C_{n,k}(L)` | 假（`χ₂L⁴`，`HS-L3-5:15`） |
| 局部代数（有界类） | Tychonoff | 常数 1 | `[T]` 空洞 |
| `S'(ℝ⁴)` | — | 盒外控制 + 核空间 | `[NG]`（不可达，只作输出） |

**为什么是拓扑而不是阶数**：三处读数不同、缺的是同一个「界」，说明缺口不在紧性判据的强弱，而在**要求哪个半范数族有界**——这正是拓扑选择。标量 `O_B≤C−Var(q)` 住在一个坐标上，子列收敛住在全族坐标的积拓扑上，故「2 阶矩 ⇒ 子列」这个问法本身就是拓扑错位。→ 下一手：`L3` 的验收句固定为 `(U1)(U2)(U3)`；`L→∞` 侧另立（定额化/固定支撑或 `MG`）。

## §3 判定：`L3 ⇒ 子列收敛` 有没有 `[T]`

**(A) `L3` 读作 `UVS(L)`（全族界）+ 矩坐标 + 固定 `L`：能，`[T]`。** 骨架（= Lean 签名，`GL3Limits.lean:58–69`）：

1. `j` 枚举 `(n,f₁,…,f_n)`，`f` 取可数稠密族；记 `S k j = S^{(n_j)}_{L,a_k}(f^{(j)}_1,…,f^{(j)}_{n_j})`。
2. 输入 `hS : ∀k j, |S k j| ≤ C j`（= 步 2 的 `(3.U1)`）。
3. `S k ∈ K := ∏_j[-C_j,C_j]`，`K` 紧（Tychonoff）。
4. `K` 紧 + `ℕ→ℝ` 第一可数 ⇒ `IsCompact.tendsto_subseq`：子列 `φ`、极限 `L∈K`。
5. `tendsto_pi_nhds` ⇒ 全 `j` 同时收敛；正性/对称闭；`hMinlos` ⇒ `μ_L`。

**(B) `L3` 读作标量 `O_B≤C−Var(q)`：`[NG]`。** 精确障碍三条：①一个坐标 ↛ 半范数族（Banach–Alaoglu 要范数球，一个 `f₀` 的值给不出 `sup_{‖f‖≤1}`）；②`n=2` ↛ `n≥3`（`μ_ε` 反例，无界类）；③无 `(U3)` ⇒ 生成函数在 `0` 处不连续 ⇒ 无测度，即使有子列极限也只是 `ℝ^{ℕ}` 里的一个点。

**(C) `S'(ℝ⁴)` 或 `L→∞`：`[NG]`。** 障碍 = 盒外无控制 + 字面 `UVS-L` 为假（`χ₂L⁴`，`HS-L3-5.md:15`）；修法只有定额化/固定支撑，修好后内容落到 `MG`/`W1`（`HS-L3-5.md:19`）。

**净判。** `L3 ⇒ 子列收敛` 的 `[T]` **存在**，但它是「消费 `L3` 的紧性包装」，不减 `L3` 一分——`hS` 就是 `L3` 本身。故 `(i)` 应改写为两行：**盒中 `μ_L` 存在 ⟸ `L3`**；**`ℝ⁴` 上存在 ⟸ `L3` + stage-2 紧性（字面假，须改陈述）+ `W1`**。→ 下一手：位置图按两行记，`Clay` 不加格（`L3` 与 `W1` 并联，`HQ1-L3.md:190`）。

## §附 工具状态（`AGENTS §1` 四要素）

| 工具 | 状态 | 实际调用 | 命中 |
|---|---|---|---|
| `pwsh` | OK（exit 0） | `Get-Location`；`Get-ChildItem` 顶层 / `species` / `GGG*` / `*COMPACT*` / `GG*.md` / `glue-a-limits` 递归 | 顶层 23 项；`species` 1064 项；`GGG*`→4（GGG2/3/4/5）；`*COMPACT*`→1（本文件）；`GG*.md`→7 |
| `glob` | OK | `research/v3/species/GG*.md`；`**/glue-a-limits/**/*.{md,tex,lean}`；`research/unconditional/CONTINUUM-TIGHTNESS.md` | 3（创建前）/ 7（现在）；10；1 |
| `grep` | OK（无 provider failure / ENOENT / 空响应） | ①`O_B\|Var\(q\)\|\\mathrm\{Var\}` @ `research` `*.md` ②`Prokhorov\|紧性\|tightness\|子列\|S'\(\|弱紧\|预紧` @ `research/v3` `*.md` | ①250 显示 / 796 全量（spill 落盘）②250 显示 / 422 全量 |
| `read` | OK（带行号正文） | `GL3Limits.lean`(104)·`CONTINUUM-TIGHTNESS.md`(13)·`glue-a-limits/REPORT.md`(405)·`HS-L3-5.md`(128)·`HQ1-L3.md`(286)·`G1-CONTINUUM-OS.md`(15)·`GGG5-RECON.md`(47)·`GGG3-OSCLOSED.md`(53)·`GGG4-W1.md`(57) | 9 文件、1108 行 |

`TOOL-FAIL`：0；空响应：0；`PATH-MISS`：0。直读 vs 转引：上列 9 份为直读；`HZ8-APHYS.md:215`、`GK5-CONTINUUM.md:227`、`UNIF3.md:151` 及 `HS-L3-5/HQ1-L3` 外的 species 行只按 grep 命中行采信（上下文未核，不升任何判）。外部标准（Minlos/Bochner、Banach–Alaoglu、Prokhorov、OS）全部未逐字核正文，`[T]` 只对「标准定理 + 仓内陈述一致」这一层。

*GGG2-COMPACT · 只新建本文件；既有一字未改；未 `git commit`；零 MC、零 Lean 编译；位置图更新交 `lead`（`AGENTS §6`）。*
