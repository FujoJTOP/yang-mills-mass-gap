# Z06 · 洞 A：`u = c_f` vs `u = c_f/N` 的归一化冲突

> 2026-11-25 · 子代理 `z06-hole-a` · 判官位（≤15 min 补洞）
> 判据（任务书指定）：**Gupta (14.20) `(β/2N²)^{4T}`**。不改任何源文件，只写本文件。
> 标签 `[T]/[C]/[NG]`；`[O]` = 仅转引（本会话未重读 PDF）。

---

## §0 判决（先给结论）

**(a)** —— **`u := c_f(β)/N` 对**。`P3-DRESSED-MASS.md:28,36` 必须改；`A^{dressed}` 一侧**不动**。

| 对象 | 值 | 标签 |
|---|---|---|
| `c_f(β)`（基本表示**字符系数**，`=` 未归一化单片迹期望 `⟨tr U_p⟩`） | `β/(2N)+O(β³)` | `[T]` |
| `u`（**每 tile 的归一化权重**，`=` `⟨(1/N)tr U_p⟩`） | `c_f/N = β/(2N²)+O(β³)` | `[T]` |
| 四壁管权重 | `w(B_T) = a₀·(c_f/N)^{4T} = a₀u^{4T}` | `[T]` |
| `A^{dressed}` | `= a₀D_∞ = a₀(1+O(u))`，`a₀=(2N)^{-p}s`，`s∈ℤ_{≥1}` | `[T]` |
| `f₁²` 与 `A^{dressed}` | **`f₁²_{,L} = A^{dressed}(L)`（mod `[KL]`）**，即 `A^{dressed}` **就是** `f₁²` | `[T]` mod `[KL]` |

`DR:87` 的 `A^{dressed}=a₀` **不需要改**；Z05 的推论"`A^{dressed}=a₀`（DR:87）**不是** `f₁²`" **撤回**——那是混用归一化造成的假象，不是残数本身的问题。

---

## §1 `c_f` 的显式定义：它是**每片的活动度**，不是"每片每色"

`P3-DRESSED-MASS.md:28` 逐字：

> `e^{(β/N)Re tr U_p} = Σ_r c_r(β)χ_r(U_p)`. The **fundamental** coefficient is `u:=c_f(β)∈(0,1)`, `u≍β`

**`c_f` 的正确定义（`[T]`，由 Weyl 正交一步得出，不依赖任何外部文献）：**

`c_f(β) = ∫_SU(N) dU\, χ_f(U)\, e^{(β/N)\mathrm{Re\,tr}\,U} = ⟨tr U_p⟩ = β/(2N)+O(β³)`.

它**没有**任何"每色"的东西：`χ_f = tr` 已经是 `N` 维指标求和，`∫χ_fχ_f̄ = 1` 把这个 `N` 吃掉了，剩下的就是 `β/(2N)`。所以

> **`c_f` = 每片（未归一化迹）的活动度。它不是"每片每色"，也不是"每片"的最终权重。**

**"每片"最终权重是它再乘一个 `1/N = 1/dim f`：**

`⟨(1/N)tr U_p⟩ = c_f/N = β/(2N²)+O(β³)` `[T]`（同一行积分：`(1/N)∫dU trU·(β/N)RetrU = (1/N)(β/N)(1/2) = β/(2N²)`，用 `∫|trU|²=1`、`∫(trU)²=0`）。

**这就是"每片带一个 `1/N`"的确切来源**：不是 `SU(N)` 的"每色权重"这种含糊说法，而是**闭合迹归一化 `(1/N)tr` 在每一个 tile 上的 `1/dim f`**。Gupta 把它叫作 "another factor of `1/N` from the integration"（§3 逐字）。

**精确形式（`[T]`，可选的 `c₀` 修正）**：`u = c_f/(Nc₀)`，`c₀ = ∫dU e^{(β/N)RetrU} = 1+O(β²)` 是真空片系数；SU(2) 精确地 `u = I₂/I₁`（`R4-PERCENT:19` 的 `u_tile`，`R4-PERCENT:116–118`；`AUDIT-R5B:222` 给 `u_P3=c_f/c₀=2I₂/I₁=1.0144`，`I₂/I₁=u/2` ✔）。领头阶 `c₀=1` ⇒ `u=c_f/N=β/(2N²)`。**洞 A 是 `1/N` 的归属问题，不是 `c₀` 问题**；`c₀` 已由 `Q4-LOGCONVEX:23` 折进 `γ`。

---

## §2 `N^{4T}` 从哪来：`4T` = `T` 片 × 4 壁，**每个 tile 一个 `1/N`**

1. **`4T` = 管的 tile 数** `[T]`（`P3-FOURWALL-DERIVATION §1.2–1.5`、`DR:34–35`）：源片 4 条边 ⇒ 4 条通量链；只有 temporal plaquette 推进时间 ⇒ 每条链跨 `T=t/a` 片，每片每壁 1 个 temporal plaquette ⇒ **每时间片 4 个、共 `4T` 个**。`4T` 是 tile 数，不是长度本身；`4` = 壁数。
2. **每 tile 一个 `1/N`** `[T]`：见 §1，`1/N = 1/dim f`，来自该 tile 的闭合迹归一化。`4T` 个 tile ⇒ **`N^{-4T}`**。
3. **`a₀` 不能接这个因子** `[T]`：`DR:36`、`DR:119`、`T2-TUBE-OVERLAP:114`、`Z04-AMPFLOOR:25–28,101` 四处**全写** `a₀=(2N)^{-p}s`，`s∈ℤ_{≥1}` 是 link 重耦合重数的**整数**乘积（`s=1`，`Z04 §2.3` 已证）。整数**装不下** `N^{-4T}`；而 `a₀` 又被定义为 `T`-无关的振幅（`DR:36` "the `T`-independent amplitude"）。⇒ **`N^{-4T}` 必须由 `u^{4T}` 携带 ⇒ `u = c_f/N`。**
4. **`(2N)^{-p}` 不是重复计入** `[T]`（回答 `Z04-AMPFLOOR:134` 的"双重计入"疑点）：`p=2` 来自 `O` 的**算符**归一 `O=(2N)^{-1}Re tr U_p`（`DR:25` 的定义），它属于源/汇两端；`c_f/N` 属于中间 `4T` 个 tile。两者是不同的一次性因子，无重叠。

**三种互相独立的确证（同一结论：tile 权重 `= c_f/N`）：**

| # | 证 | 内容 | 标签 |
|---|---|---|---|
| 1 | 单片**精确** | `⟨W₁₁⟩ ≡ ∫dU (N^{-1}tr U)e^{(β/N)RetrU} = c_f/N + O(c_f³)`；同时 `⟨tr U_p⟩ = c_f` ⇒ 未归一化给 `c_f`，归一化给 `c_f/N`。**`1/N` 在迹归一上，一次一片。** | `[T]` |
| 2 | 2D **精确**（面积律按片因子化） | `⟨W(A)⟩ = (w₁)^A` 是**纯 `A` 次幂** ⇒ `A` 个 tile 就带 `A` 个 `1/N`（不是每个 link、也不是整体一个）。若 `1/N` 只在整体迹上，`(1/N)c_f^A N^{k}`（`k` 与 `A` 无关）**不可能**化成 `(β/2N²)^A`。 | `[T]`（经典） |
| 3 | Gupta 逐字 | (14.18) `⟨W_RT⟩=(β/2N²)^{RT}`、(14.20) `(β/2N²)^{4T}`、且**明说**每 tile 一个 `1/N`（§3） | `[O]`（转引 `FOURWALL:150–153`，该文件称 "read from PDF"；**本会话未重读 PDF**） |

---

## §3 Gupta 的逐字（源：`P3-FOURWALL-DERIVATION.md:150–153`，标 `[O]`，称读自 PDF）

- **(14.17)** `W_{RT} = (1/Z)∫dU\,W_{RT}[1+(\beta/2N)(W₁₊₁+W₁₋₁)+…]` —— **每片展开因子 `β/2N`**。
- **(14.18)** 每个 tiling plaquette "brings a factor of `β/2N` from the expansion and another factor of **`1/N` from the integration** … combined factor **`β/2N²`**"（附注："SU(2) is special and the combined factor is `β/N²`"）；`⟨W_{RT}⟩=(β/2N²)^{RT}(1+…)`。
- **(14.20)** "The minimum tiling … is the rectangular tube consisting of four sides of size `1×T` … `W₁₁(T)W₁₁(0) ∼ e^{−M₀₊₊T} = (β/2N²)^{4T}(1+…)`."
- **(14.21)** `M₀₊₊ = −4log(β/2N²)+…`.

**`β/2N²` 是哪来的** `[T]`：`(β/2N)`：`e^{(β/N)RetrU_p}` 展开的首项系数，**这恰好就是 `c_f`**（含 `Re tr` 的 `1/2` 与 `tr` 的 `1/dim`）。(`1/N`)：该 tile 的 `1/dim f` 归一。两者之积 `= c_f/N = β/(2N²)`。

**与标准 `SU(N)` 单片展开一致吗** `[T]`：**一致，且它就是标准值**。`⟨W₁₁⟩ = c_f/N = β/(2N²)+O(β³)` 与 Gupta 完全同式（SU(2) 特例：两定向相等 ⇒ `c_f` 加倍为 `β/2` ⇒ `c_f/N=β/4=β/N²` ✔ `FOURWALL:167–171`、`R4-PERCENT:118`）。**不是约定差别，是可算的同一个数。**

---

## §4 判：**(a)**。`u = c_f/N` 对，`DR:28,36` 要改

**(b) 不成立。** 若坚持 `u=c_f`，则须把 Gupta 的 `(β/2N²)` 重新解释成"`β/2N` × 某个非物理解释的 `1/N`"。但 §2 表第 1、2 行把 `1/N` 定死为 `1/dim f`/每 tile（精确单片 + 2D 精确幂律），没有任何解释空间。`[T]` ⇒ (b) 死。

**(c) 不成立（但有一处限定要说清）。** 不是"polymer 用 `c_f`、管用 `c_f/N`"的分工：`c_f` **从来不是任何权重**，它只是**字符系数/未归一化迹期望** `⟨tr U_p⟩`。唯一"`c_f` 单独出现"的合法场合是**上界**——`LOOP-RUNBOOK:1827` 逐字 `|w_β(S)| ≤ ∏_{p∈S}c_f(β)`（`PROOF.md` 引理 2′）。这是**界**，比真权重松 `N^{|S|}`；KP 判据用现成的界是合法的，但**判据半径必须与权重同归一**：取 `u=c_f/N`（`Z02-SEILER §4 表 C`，`uK<1` ⇒ `β_W<0.0368N²`）；若取 `u=c_f`（表 A）则把 `β₁` 收紧 `N` 倍，与管子本身的权重不自洽。⇒ 形式上像 (c)，**实质是 (a) + 一条"界 ≤ ∏c_f"的旁注**，不是两条并存的正确定义。

**`P3-DRESSED-MASS.md` 的两处最小修改（逐字给出）：**

- **`:28`** 现文 "The **fundamental** coefficient is `u:=c_f(β)∈(0,1)`, `u≍β`"
  → "The fundamental **character** coefficient is `c_f(β)∈(0,1)`, `c_f=β/(2N)+O(β³)`; the **per-tile normalized weight** is `u:=c_f(β)/N∈(0,1)`, `u≍β/(2N²)` (Gupta 14.18: `β/2N` from the expansion × `1/N` from the integration). `[T]`"
- **`:36`** 现文 "where `∏_{p∈B_T}c_f = u^{4T}` (each of the `4T` temporal plaquettes carries `c_f=u`)"
  → "where `∏_{p∈B_T}(c_f/N) = u^{4T}` (each of the `4T` temporal plaquettes carries the **normalized** tile weight `c_f/N=u`; the `1/N=1/dim f` per tile is Gupta's `1/N` from the integration). `[T]`"

`:35` 的 `κ:=4\ln(1/u)/a` **形式不变**，只是 `u` 换成 `c_f/N`，于是 `κ = 4\ln(N/c_f)/a = −4\ln(β/2N²)/a`，与 (14.21) 逐字相符（`[T]` 代数 + `[O]` 源）。

---

## §5 统一后 `A^{dressed}` 与 `f₁²` 的关系（逐字）

统一取 `u=c_f/N` 后，`DR:86–89` 的链条**一步不改、且现在自洽**：

`Σ_{P∈𝔓_tube}w(P) = a₀u^{4T}Z_dec(B_T) = a₀u^{4T}e^{f_dec t}D_{fin}(T) = A^{dressed}e^{−m₁^{dr}t}D_{fin}(T)`

其中 `u^{4T}=e^{−κt}`、`A^{dressed}=a₀`、`m₁^{dr}=κ−f_dec`。文字化：

> **`A^{dressed} = a₀ = (2N)^{-p}s`（`T`-、`L`-无关，`>0`），`D_{fin}(T)=1+O(u)→D_∞=1+O(u)`，
> 故 `f₁²_{,L} := dμ_O^{(L)}({m₁^{dr}(L)}) = A^{dressed}(L) = a₀D_∞ = a₀(1+O(u))`。
> 即 **`f₁² = A^{dressed}·D_∞`，`A^{dressed}` 与 `f₁²` 是同一个量**（差 `D_∞=1+O(u)`），
> 下界 `f₁² ≥ ½a₀ ≥ ½(2N)^{-2}`，`inf_L f₁² ≥ ¼(2N)^{-2}`（`DR:122–125` 原封不动）。**`[T]` mod `[KL]`**

**这正是 Z05 的"`A^{dressed}=a₀` 不是 `f₁²`"被撤回的机制**：混用时 `w(B_T)^{DR:36}=a₀c_f^{4T}=a₀N^{4T}(c_f/N)^{4T}`，于是残数被高估 `N^{4T}`、且若仍写 `A^{dressed}=a₀` 则 `C(t)e^{m₁t}` 的极限不是常数 ⇒ Lemma 3.4（`DR:99–101`）前提失效。**改 `u` 之后 `A^{dressed}=a₀` 恢复 `T`-无关，Lemma 3.4 前提自动满足，`A^{dressed}=f₁²` 成立。**

**会动 / 不会动的清单：**

| 条目 | 变化 | 标签 |
|---|---|---|
| `a₀=(2N)^{-p}s>0`、`s=1`、`f₁²` 的下界 `r₀=¼(2N)^{-2}` | **不变** | `[T]` |
| `A^{dressed}=a₀`（`DR:87`）、`(4.2)(4.3)`（`DR:122–125`） | **不变** | `[T]` |
| `m₁ = 4\ln(1/u)/a` | 形式不变；数值 **上移 `4\ln N/a`**（`u:c_f→c_f/N`）⇒ 与 (14.21) 对齐 | `[T]` |
| 管权重 `∏c_f` | **改**为 `∏(c_f/N)`（`:36`） | `[T]` |
| `u≍β` | 改为 `u≍β/(2N²)` | `[T]` |
| KP 半径 `uK<1` | 同一 `u=c_f/N`（`Z02-SEILER` 表 C） | `[T]` |
| `DR:51` 的 `≥6T` 二分法 | **不在本洞**（洞 B，`Z05 §B`） | `[NG]` 本文件不处理 |

**酸测**：改后结论仍对 YM 正（`a₀>0` 来自 `f⊗f̄→1` 的 `1>0` 重耦合）、对自由无质量场失效（`a₀=0` ⇒ `f₁²=0`），`DR:127` 的 PASS 不因归一化修正而变。`[T]`

---

## §6 检索核验（AGENTS §1）

| # | 检索式 | 工具 | 状态 | 命中 |
|---|---|---|---|---|
| 1 | `c_f\|P3-DRESSED-MASS\|P3-FOURWALL\|14\.20` @ `research/v3/**/*.md` | `grep` | **OK** | **245** |
| 2 | `Get-ChildItem species\Z0*` | `pwsh` | **OK** | 3（Z02/Z04/Z05，**无 Z06**，本文件新建） |
| 3 | `read P3-DRESSED-MASS.md:20–94` | `read` | **OK** | 全段 |
| 4 | `read P3-DRESSED-MASS.md:95–139`（Lemma 3.4、`D_fin`、`(4.1)–(4.3)`） | `read` | **OK** | 全段 |
| 5 | `read P3-FOURWALL-DERIVATION.md:70–179`（§2.1、§3.1–3.5、§4.1–4.5 Gupta） | `read` | **OK** | 全段 |
| 6 | `read PROOF-P3-UNCONDITIONAL.md:14–36` | `read` | **OK** | 全段 |

**无 TOOL-FAIL。** 未逐字核：**Gupta arXiv:hep-lat/9807028 的 PDF**（`[O]`，按 `FOURWALL:150–153` 的"读自 PDF"记录转引，本会话未重读）；Münster 两篇（同）。

*Refs.* 内部：`P3-DRESSED-MASS.md:25,28,34–36,86–89,99–101,106,119,122–127`；`P3-FOURWALL-DERIVATION.md:79–86,150–155,159–171`；`PROOF-P3-UNCONDITIONAL.md:22,25,28,31`；`Z05-AUDITI1.md:62,88`；`Z04-AMPFLOOR.md:25–28,101,134`；`X01-L1C.md:15–16,92`；`W1-STRONG-ALIGN.md:83,196–199`；`R4-PERCENT.md:19,116–118,130`；`AUDIT-R5B.md:222`；`Z02-SEILER.md:46,82–84,114,121`；`Q4-LOGCONVEX.md:23`；`LOOP-RUNBOOK.md:1827,5829,6208–6212`。外部（转引，未重读）：Gupta hep-lat/9807028 §14.2.2 (14.17)–(14.21)；Münster PLB 95 (1980) 59 / NPB 190 (1981) 439；Seiler LNP 159 (1982) §II；Kotecký–Preiss CMP 103 (1986) 491。

**No QED on anything beyond `(a)`.** 洞 A 本身（`u` 的归属）判死为 `(a)`，`[T]`：判据是 Gupta 逐字 + 两个精确一步计算（单片、2D 幂律）+ `a₀` 的整性与 `T`-无关性。洞 B–F（`Z05:89–93`）本轮不动。Clay 0%。
