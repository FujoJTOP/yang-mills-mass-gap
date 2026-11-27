# HW2-A6 — `PROOF-L3:6` 的两行 `N_□` 论证：用的是**未减除**矩 ⇒ A6 旁路**未被关闭**（活）；但 A6 的链条有一格终止在 `(★)` 之外、且 `A2:6` 逐字把它认作 `L3` 本身 ⇒ A6 是 `L3` 的**等价重写（finiteness 面）**，不是捷径

> **轮**：HW2 · 子代理 `hw2-a6` · 工作目录 `D:\Dev\math`。**写域**：只新建本文件；既有文件一字未改；未 `git commit`；零新 MC；零 Lean；未动 `docs/position-current.html`。
> **任务（只做两件）**：(1) 把 `PROOF-L3:6`（=`research/unconditional/PROOF-L3.md:6`）的 `N_□` 论证**逐字**抄出并判：用未减除还是减除后矩？(2) 据此判 A6：活则给逐格链条，死则给具名理由。
> **纪律**：只交 `[T]`/`[NG]`；每个数带 `file:line`；每条结论带「下一手」；`AGENTS §1` 四要素见 §3。不用百分数。

---

## §0 判词（先给；每条带「下一手」）

| # | 断言 | 判 | 逐字依据 / `file:line` | 下一手 |
|---|---|---|---|---|
| **0.1** | **`N_□` 论证用的是【未减除】矩**：`N_□=6(L/a)⁴`（`L` 固定）里唯一的 `a`-依赖就是 `a^{-4}`，而 `a^{-4}` 正是 `(U2)` 要减掉的那一幂 ⇒ 该句等价于「`p_bare=−4`」，不是「减除后 `p=0`」 | **`[T]`** | `PROOF-L3:6`（逐字见 §1）；`PROOF-CONCENTRATION:6`（`N_□=6(L/a)⁴`）；`HS-L3-6:0.2`（R1 `p_bare=−4n`、R2 `n≥2` 最坏 `−4`、**R3 减除后 `p=0`**，全 `[T]`）；`glue-a-limits:90` | 凡引 `PROOF-L3:6` 此句，一律注「讲的是裸矩」 |
| **0.2** | **它是"对象偏离"，不是定理**：主语名词写「renormalized moments (`UVS(L)`)」，支点却是「compactness gives `O(N_□)`」——紧性给的是 `osc≤2βN_□`（`PROOF-CONCENTRATION:7`），**加性减除改不动振荡** ⇒ 支点天生在减除之前 | **`[T]`** | `PROOF-L3:6`；`PROOF-CONCENTRATION:7`；`glue-a-limits:90`（`(U2)` 的用途逐字） | 把该句降级为「裸界备忘」，不进 `L3` 的必要条件表 |
| **0.3** | 该句结论 **"a-uniform tails require log-concavity"**（及同段的 `L3 ⟺ (★)`）判 **`[NG]`**：`(★)` **严格强于** `L3`、只是充分条件；自由场上 `(★)` 精确假而 `L3` 真 | **`[NG]`** | `HQ1-L3:59`（逐字「**严格强于 `L3`**…充分条件」）·`HQ1-L3:63`（自由场 `G″=0` ⇒ `(★)` 假，`L3` 必真）·`HS-L3-10:0.4`（`¬(★) ⇏ ¬L3`） | 删掉「require」；`(★)` 只留在 `𝒮_t` 上的充分条件位 |
| **0.4** | **⇒ A6 活**：拦阻（`HS-L3-7:49` 全部拦阻 = 这一句）不成立；A6 的链条**没有一格消费 `(★)`、软模或 `U′` 量词** | **`[T]`** | 0.1+0.2+0.3；§2 逐格 | 见 §2 三条下一手 |
| **0.5** | **但 A6 不给捷径**：链条唯一实质格 `A6-1`（行 1 的非微扰面）= `A2` 的 `L3` 层，`A2:6` 逐字「The additive `a⁻⁴` is benign **provided** the fluctuation is controlled — exactly L3」 | **`[NG]`** | `A2-scheme-match:3,6,12`；`G1-EXISTENCE-THEOREM:21`（"the nonperturbative a-uniform bound **is exactly L3**"） | 把 A6 记账为 `L3` 的等价重写，不记「旁路/缩短」 |

**一句话。** `PROOF-L3:6` 那两句的计数 `N_□∼L⁴/a⁴` 与 `(U2)` 的 `a^{-4}` 是**同一幂**：`N_□=6(L/a)⁴`，`L` 固定时唯一的 `a`-依赖是 `a^{-4}`；而 `(U2)` 的全部用途就是减掉这一幂（`glue-a-limits:90`），减除后 `p=0` 是仓内 `[T]`（`HS-L3-6` R3）。故该句讲的是**裸矩**，把裸幂安到「renormalized moments」上——对象偏离，不构成关闭 A6 的论证。⇒ **A6 活**。但逐格写下去，`(★)` 腿被合法判掉（`HQ1-L3:59` 只充分），债全落到一格：`sup_a[ρ_a^{ren}𝒢_a^{ren}]<∞`（`HS-L3-9:2.3`），而这一格就是 `A2:6`/`G1-EXISTENCE:21` 认领的 `L3`。

---

## §1 `PROOF-L3:6` 逐字（只切句，不改字）

| 段 | 逐字 | 它在说什么（本文件判） | 判 |
|---|---|---|---|
| **6a** | 「With the exact identity `G″_j=E_{μ(·|U′)}[PᵀH_jP]−Var_{μ(·|U′)}(Pᵀ∇S_j)+W_j` (PROVE-Gpp), L3 is: the normalized generator is bounded, `sup_a Z_a^{-1}‖H_a‖_D<∞` (SPEC-RES′), **and** the renormalized moments are a-uniform (UVS(L)).」 | 目标定式：两面并列（生成元面 + 矩面）。与 `HQ1-L3:37–44` 同 | `[T]`（定式） |
| **6b** | 「The second half is the substance: `UVS(L)` moments scale with `N_□∼L⁴/a⁴→∞`,」 | **本核**。主语名词是「`UVS(L)` moments」（前句刚说是 renormalized），但给的数是 `N_□=6(L/a)⁴` ⇒ `a^{-4}`，即 R1/R2 的**裸幂** | `[T]`（该数为真）/ **`[NG]`**（安错对象） |
| **6c** | 「so operator-norm *finiteness* (compactness gives `O(N_□)`) is not a-uniform」 | 支点 = 紧性界 `osc(βH_fib)≤2βN_□`；`Z_a^{-1}` 未入、减除不入（加性常数不改振荡）⇒ 天生是**未减除/未归一**的计数 | `[T]`（该界为真）/ **`[NG]`**（作"减除后仍不 a-一致"的理由） |
| **6d** | 「— a-uniform tails require log-concavity.」 | 结论句。"require" 是必要性断言 | **`[NG]`** |
| **6e** | 「Hence L3 ⟺ the **missing inequality** (★) plus, for the fixed-`L`, `a→0` count `N_{E′}∼1/a⁴`, the decay `|Cov(s_B,s_B′)|≤Ce^{−κ|B−B′|/a}` (`★dec`).」 | 升成 `⟺`；仓内后判：只 `⟸`（`(★)` 充分、且 `∀U′` 版为假） | **`[NG]`**（`⟺` 的那一向） |

**判（问题 1）：用的是【未减除】矩。** 三条独立理由，逐条可核：

1. **幂对齐**：`N_□=6(L/a)⁴`（`PROOF-CONCENTRATION:6`）。`L` 固定 ⇒ `N_□∝a^{-4}`；`a^{-4}` 就是 `tr F²` 的加性发散幂（`glue-a-limits:90`）。减除后 `p=0`（`HS-L3-6` R3，`[T]`）⇒ 「scales with `N_□`」只可能是减除前的陈述。
2. **支点对齐**：6c 的理由是紧性，而紧性给的是 `osc`（`2βN_□`）——减常数不改振荡 ⇒ 该 `O(N_□)` 无论怎么读都在减除之前。
3. **文本对齐**：同段 6e 自带 `N_{E′}∼1/a⁴` 的裸格点计数；该段的"计数语言"自始至终是裸的，只在 6a 出现过一次「renormalized」。

---

## §2 A6 判定（活）+ 逐格链条

**A6 的陈述**（沿用 `HS-L3-7:49,87`）：行 1（`(U2)` 减除，非微扰一致）+ 行 7（头部 finiteness）⇒ 尾部 finiteness，不经软模、不经 `(★)`。

| 格 | 内容 | 判 | 依据 | 残余 / 下一手 |
|---|---|---|---|---|
| **A6-1** | 行 1 = `(U2)`：`a^{-4}` 用**同一** scheme C 减除。`A2` 的三分解：`L1`[定义] + `L2`[全阶定理] + `L3`[开] | **`[T]`**（L1+L2）/ **`[NG]`**（非微扰） | `A2:3,9,12`；`G1-EXISTENCE:21` | `A2:6` 逐字「benign **provided** the fluctuation is controlled — exactly L3」⇒ 这一格就是 `L3`，不是它的替代 |
| **A6-2** | 计数被吃：`N_□` 的 `a^{-4}` = `(U2)` 减掉的那一幂；减除后 `p=0`（+`log^k`）；`L⁴` 部分只是 `C_n(L)` 允许的常数 | **`[T]`** | `HS-L3-6:0.2`（R1–R3）；`glue-a-limits:90` | 此格闭合，无残余 |
| **A6-3** | `(U3)` 不独立：树界给 `C_n≤A(a)^n n!P(L)`；故 `(U1)∧(U2)∧(U3)` 的**全部** `a`-内容 = 一个标量对 `sup_a[ρ_a^{ren}𝒢_a^{ren}]<∞`（同方案减除后） | **`[T]`**（归约） | `HS-L3-9:0.5,2.3` | 该标量对本身 `[NG]`；这是 A6 唯一实质格 |
| **A6-4** | 行 7 = 头部 finiteness：`E′` 固定 + 有界（负）Hessian ⇒ 有界自伴 ⇒ `C_head(L)<∞`，**不求正性** | **`[T]`** | `L3-scale-uniform:10` | 块局部版已具名：`Λ_UV(b)=e^{−2βN_p(B)}/C_G`，`U′,L,a` 一致（`PROVE-H2-uvgap:12–18`）——只喂有限性，不喂地板 |
| **A6-5** | 塔的归一：`Z_j=Θ(β_j)` ⇒ `Z_j^{-1}H_j=O(Q)` 有界（免费） | **`[T]`**（自由块）/ **`[NG]`**（相互作用 `(log)^γ`） | `V1:78,237,320`；`HS-L3-6:89`（`Z_a` 块和/块平均未钉，**一行可闭**） | 先钉 `Z_a` 的块和 vs 块平均约定 |
| **A6-6** | 矩面 ⟺ 生成元面（F1⟺F2）的逐行等价 | **`[NG]`** | `HQ1-L3:58`（仓内只断言 "equivalently"，`log T_a` 展开未落文） | 不做也不伤 A6：A6 全程只在矩面用 |
| **A6-7** | `(★)` 腿**不必要**：`(★)` 严格强于 `L3`（只充分）；自由场 `(★)` 假、`L3` 真 | **`[T]`** | `HQ1-L3:59,63`；`HS-L3-10:0.4`；`PROOF-L3:22`（"L3 is **true** but not proven"） | 把 `(★)` 从本链的必要条件表删除（这是 A6 的真产出） |

**净账（`[T]`）**：A6 的七格里 **五格 `[T]`**、无一处消费 `(★)`/软模/`U′`；`[NG]` 只有三处：**A6-1 的非微扰面**（= `L3` 层）、**A6-3 的标量对**（= 同一物的矩面写法）、**A6-5 的 log 余项**（一行约定 + 一个 `(log)^γ`）。⇒ **A6 活**（未被 `PROOF-L3:6` 关闭，也未被任何具名 no-go 拦），但它的产出是**重写**：`L3` ⟺ `sup_a[ρ_a^{ren}𝒢_a^{ren}]<∞` ∧ `sup_a Z_a^{-1}‖H_a‖_D<∞`，两头同格，且都不是 `(★)`。

**下一手（按性价比）**

| 序 | 动作 | 类型 | 代价 |
|---|---|---|---|
| **1** | 把 A6 记账为「`L3` 的 finiteness 面重写」，`HS-L3-7:49` 的 A6 行从 `[NG]` 改 `[T]`（活）但去掉「旁路/缩短」字样 | 记账 | 零 |
| **2** | 钉 `Z_a` 的块和 vs 块平均（`HS-L3-6:89` 已指出它决定 `Z_a` 的精确幂）；钉后 A6-5 的 `(log)^γ` 变成可算 | 约定 + 幂计数 | 半轮 |
| **3** | 攻 A6-3 的标量对：`sup_a[ρ_a^{ren}𝒢_a^{ren}]<∞`——它是 `L3` 的全部 `a`-内容，且在自由场上可精确标定（`HS-L3-9:3.3` 已给 `χ_2^{ren}`） | 解析（唯一真格） | 多轮 |

---

## §3 台账（`AGENTS §1` 四要素）

| # | 工具 | 状态 | 实际调用（逐条） | 命中 |
|---|---|---|---|---|
| 1 | `read` | **OK**（全部返回带行号正文） | `PROOF-L3.md`(全 25)·`HS-L3-7.md`(全 129)·`HS-L3-6.md:1–60`·`HS-L3-9.md:36–65`·`HS-L3-10.md`(全 118)·`HW3-K6.md`(全 90)·`HQ1-L3.md:28–107`·`L3-scale-uniform.md`(19)·`G1-L3.md`(29)·`G1-UNIFORM-ZA.md`(33)·`PROVE-H2-uvgap.md`(29)·`G1-L3-VIA-FLOW.md`(25)·`G1-EXISTENCE-THEOREM.md:1–30`·`A2-scheme-match.md`(全 19)·`V1-UV-NORMALIZATION.md`(全 358)·`glue-a-limits/REPORT.md:30–149` | 16 份 |
| 2 | `glob` | **OK** | `**/PROOF-L3*`·`**/G1-L3-VIA-FLOW.md` | 1·1 |
| 3 | `grep` | **OK**（无 provider failure、无 ENOENT、无空响应） | ①`HS-L3-7` @ 仓 ②`N_□` @ `*.md` ③`A6` @ `research` ④`减除` @ `research` ⑤`renormaliz\|additive divergence` @ 仓 | ①**5** ②**125** ③**87** ④**135** ⑤**524**（显示 250，spill 落盘） |
| 4 | `pwsh` | **OK** | `Get-Location`·`Get-ChildItem`×3·`glob` 定位 `A2-scheme-match.md`/`G1-EXISTENCE-THEOREM.md` | 全部有效响应 |
| 5 | `write` | **OK** | 新建 `research/v3/species/HW2-A6.md` | 本文件 |
| 6 | `py`/MC/Lean/`web_*`/`git`/`skill` | **0 次** | — | — |

**`TOOL-FAIL`：0 次。空响应：0 次。`PATH-MISS`：0 次。** 承重行全部直读原文（`PROOF-L3:6`、`A2:6,12`、`G1-EXISTENCE:21`、`HS-L3-6:0.2`、`HS-L3-9:2.3`、`HQ1-L3:58–59,63`、`glue-a-limits:90`、`V1:237`、`L3-scale-uniform:10`）；转引仅 `ch15:34–47`（经 `HQ1-L3:59`）与 `PROOF-CONCENTRATION:7` 的 `osc` 式（经 `PROVE-H2-uvgap:12` 同式复核）。

### 诚实边界（4 条）

1. **本文件不证 `L3`、不证其否定、不产估计。** 两个产物是：`N_□` 句的对象判定（0.1–0.3）与 A6 的逐格账（§2）。
2. **0.4 的「活」只裁「未被 `PROOF-L3:6` 关闭」**，不裁「A6 可交付」；0.5 已把它的价格写明（行 1 = `A2` 的 `L3` 层）。
3. **§1 的切句是我的切法**：`PROOF-L3:6` 是一整段，原文无 6a–6e 标号；逐字未改一字，只加切点。
4. **`N_□` 与 `(U2)` 的幂对齐用了 `L` 固定**：这是 `L3` 的窗口（`HQ1-L3:72`）；若读成 `L=c·a` 则 `N_□` 不随 `a` 变，那时变的是量词结构（`HS-L3-10` P2 `[T]`），本判不适用。

*HW2-A6 · 只新建本文件；既有文件一字未改；未 `git commit`；零 MC；零 Lean；未动 `docs/position-current.html`。*
