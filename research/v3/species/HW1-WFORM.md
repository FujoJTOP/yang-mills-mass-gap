# HW1-WFORM — `W` 的精确陈述：绕开 `Δ₀` 的那条路，价格是靶本身

> 第 HW1 轮 · 子代理 `hw1-wform` · 工作目录 `D:\Dev\math`。**只新建本文件**；既有文件原样；未 `commit`；零 MC；未动 `docs/position-current.html`（`AGENTS §6` 交 `lead`）。
> 靶（逐字，不重推）：`(B1^rate)` 的产出 —— "连续极限谱隙 >0（= 靶）"（`HY1-B1.md:42`）。
> 正文 `research/v3/species/_ht7/faizal.txt`，`fax:` = 其行号（arXiv v1 抽文）。标签只 `[T]`／`[NG]`。
> 本轮的账：`HZ3-VRATE.md:44` 的 `W` 行（"极限聚簇 + 不经 `Δ₀`"），以及 `:46` 的"新格 W"。

## §0 判词

| # | 断言 | 判 | 依据 | 下一手 |
|---|---|---|---|---|
| **0.1** | **`W` 的陈述逐字在 `fax:5114–5120`**；纸面自述口径"只用极限聚簇、不用有限截断谱信息"在 `fax:5110–5112` | `[T]` | §1 | 引 `W` 一律引 `fax:5114–5120`，不引二手转述 |
| **0.2** | **一行式 = `∃m>0 ∀(A,B) ∃K(A,B) ∀t≥0` 的聚簇 ⇒ `spec(H)⊆{0}∪[m,∞)`**：`m` 全局一个、`K` 逐对，产出是**下界**不是谱隙值 | `[T]` | §2 | 台账只留 `m`、`K(A,B)` 两槽，不加第三输入 |
| **0.3** | **`W` 与靶等价（ge3 换名格），不是加强** | `[T]` | §3.2 两向：`⟸` 两行（`(8.37)` + 谱范数），`⟹` Thm 8.6 | `W` 记"同难"，不记"更强"、不记"更弱" |
| **0.4** | **"充分且严格弱"这一格在靶层是空集**，故 `W` 只能落"换名/更强"；它落换名 | `[T]` | `HZ6-SYNTH.md:248`（ge3 第三格 = ∅，逻辑空） | 停止在 `W` 上找"比靶便宜"的版本 |
| **0.5** | **把 `K` 提成全局一致也不买强度**：`K(A,B)≤4‖A‖‖B‖` 由谱隙直接给 | `[T]` | §3.3 三行 | 强度不在 `K` 的一致性轴上，停止加工该轴 |
| **0.6** | **相对 `B1^rate_{a-uni}` 它严格弱，但那是推论身份**：`a`-一致 ⇒ `W` 的假设；`W` ⇏ `a`-一致 | `[T]` | §3.4；`HZ3-VRATE.md:11` | `B1^rate_{a-uni}` 单列一格；`W` 不得当它的路 |
| **0.7** | **纸面自己的供给方就是第 (5)–(7) 节**（阶跃标度 + 局部性），与 `Δ` 路同一批供给 ⇒ 换账不改价 | `[T]` | `fax:5195–5199` | `W` 省下的只是第 (3) 节的 `Δ₀^{SC}`，省不下第 (5)–(7) 节 |

**一句话。** `W` 是靶的换名：在 OS 重构环境里它的假设与结论互相推出（`⟸` 两行、`⟹` 是 Thm 8.6 本体），它唯一"新"的是**记账位置**（唯一不经 `Δ₀` 的入口），不是强度；实践上它的供给方仍是 (5)–(7) 节，故"绕开 `Δ₀`"省掉塔底那一格，没省掉靶。

## §1 逐字（`fax:5114–5120`）

```text
Theorem 8.6 (Gap from exponential clustering). Assume that, for some m > 0, there exists
for each pair of bounded, gauge-invariant observables A, B with compact, disjoint supports and
positive time separation a constant K(A, B) such that the connected Euclidean correlator

                        CA,B (t) := ⟨A τt (B)⟩ − ⟨A⟩ ⟨B⟩                       (8.35)

satisfies |CA,B (t)| ≤ K(A, B)e−mt for all t ≥ 0. Then spec(H) ⊂ {0} ∪ [m, ∞).
```

（纯文本抽文，`τt`=`τ_t`，`e−mt`=`e^{−mt}`，`CA,B`=`C_{A,B}`；抽文的减号是 U+2212。）

**口径逐字**（`fax:5110–5112`）："We next present a Tauberian mechanism that uses only exponential clustering of connected Euclidean correlations in the limit and does not rely on finite-cutoff spectral information."

**兼容性逐字**（`fax:5194–5199`）："the spectral gap of H is bounded below by min(∆, m) in our setting. … Sections (5)-(7) established step-scaling and locality estimates that imply exponential clustering in the scaling window and hence m > 0 in the limit."

## §2 一行式

```text
W :   ∃m>0  ∀(A,B)∈𝒪  ∃K(A,B)<∞  ∀t≥0 :
          | ⟨A τ_t B⟩_lim − ⟨A⟩_lim ⟨B⟩_lim |  ≤  K(A,B) · e^{−mt}
      ⟹   spec(H) ⊆ {0} ∪ [m,∞)            （质量隙 ≥ m）

𝒪        = 有界、规范不变、紧支撑且时间分离为正的观测量对
⟨·⟩_lim  = σ↓0 的极限 Schwinger 函数
H        = OS 重构的 Hamiltonian（HΩ=0，H≥0）
量词     = m 在 ∀(A,B) 之前（一个速率管所有对）；K(A,B) 在每一对之后（前置因子逐对）
产出     = 谱隙的下界 m（不是"谱隙的值"：数值层只有 m ≤ Δ_真）
```

读量词：`m` 全局 ⇒ 堵掉"逐对速率"的弱化读法；`K` 逐对 ⇒ 看起来是强度所在，`§3.3` 证明连这一条也可以去掉而不改强度 —— **`W` 的强度不落在任何一个量词上，而落在"假设即结论"上。**

## §3 强度（ge3 二分）

### 3.1 二分的口径

```text
ge3 二分（GE3-GAPLB.md:24；本仓已用形 HZ6-SYNTH.md:242–249）：
    任一 R 满足 R ⟹ T  ⇒  要么 T ⟹ R（换名）   要么 R 严格强于 T
    第三格（充分且严格弱）= ∅（逻辑空：充分 ∧ 更弱 ⇒ 等价）
```

### 3.2 两向（判据本体）

| 方向 | 论证 | 判 |
|---|---|---|
| `H_W ⟹ 靶` | Thm 8.6 本体（`fax:5122–5189`）：取 `B=A*` 得 `µ_A` 的 Laplace 界 ⇒ `µ_A(0,m)=0`；Reeh–Schlieder 密度把"逐 `A`"抬到 `H^⊥` ⇒ `E(0,m)=0` | `[T]` |
| `靶 ⟹ H_W` | `ψ_A := AΩ−⟨Ω,AΩ⟩Ω`；`C_{A,B}(t)=⟨ψ_A,e^{−tH}ψ_B⟩` 是**恒等**（`fax:5134`，不在假设里）；谱隙 ⇒ `‖e^{−tH}|_{Ω^⊥}‖≤e^{−mt}`；取 `K(A,B)=‖ψ_A‖‖ψ_B‖` | `[T]` |

⇒ **`W ⟺ 靶`：第一格（换名）。不是第二格（严格强），也不是"仅充分"。** 它的"绕开 `Δ₀`"是把债从塔底挪到聚簇，而聚簇与谱隙是同一句话的两个坐标。

### 3.3 `K` 的一致化不买强度

`‖ψ_A‖ ≤ 2‖A‖` ⇒ `K(A,B) ≤ 4‖A‖‖B‖`。把"逐对 `K`"换成"`‖A‖,‖B‖≤1` 上的全局常数 4"，假设仍由靶逐字给出 ⇒ `W_uni ⟺ W ⟺ 靶`。**强度不在 `K` 的这条轴上**：想靠收紧前置因子造第二格，只会造出换名。

### 3.4 与 `B1^rate_{a-uni}` 的关系（不同层，别混账）

```text
B1^rate_{a-uni}  ⟹  W 的假设（丢掉 a-一致性，只剩极限）
W  ⟹  B1^rate_{a-uni}  ✗（W 不产 a-一致、不产格上速率）
```

⇒ 相对 `B1^rate_{a-uni}`，`W` **严格弱**；但此时 `W` 是它的**推论**，不是它的充分条件 ⇒ 不进 ge3 二分（二分只排"充分者"），也不能记作"残余"或"路"。**判：`W` 只在极限层有格位。**

### 3.5 数值层（唯一的非等价处）

`W` 给的下界是假设里的 `m`；靶的谱隙是 `Δ`。`W` 只保证 `m ≤ Δ`（`fax:5194` 的 `min(∆,m)` 逐字同向）。⇒ 等价只在"`∃m>0`"层成立，不在"`m=Δ`"层。记 `[T]`，不据此升降级。

**下一手（§3）**：`W` 行只写"与靶同难（换名格）"；`B1^rate_{a-uni}` 行写"严格强于 `W`，且是 `W` 的唯一已知供给"；不要在 `K` 的轴、`m` 的数值层上找残余。

## §4 台账（`AGENTS §1` 四要素）

| # | 工具 | 状态 | 实际检索式 / 调用 | 命中 | 失败标记 |
|---|---|---|---|---|---|
| 1 | `read` | **OK**（带行号返回） | `HZ3-VRATE.md`(62)；`fax` 区间 `5040–5239`；`GE3-GAPLB.md`(1–100)；`HY1-B1.md`(28–57)；`HZ6-SYNTH.md`(224–267) | 全部返回 | 0 |
| 2 | `grep` | **OK**（ripgrep 正常） | 路径 `research/v3`、`include=*.md`、式 `ge3` | 25 | 0 |
| 3 | `glob` | **OK** | `**/*fax*` → 0（文件名是 `faizal.txt`，不含 `fax` 子串，属真 0 而非工具故障） | 0 | 0 |
| 4 | `glob` | **OK** | `research/v3/species/HW*.md` → `HW1-CLOSURE`、`HW2-A6`、`HW3-K6`、`HW4-BIAS`、`HW5-RNWEIGHT`（本文件为 `HW1-WFORM`，不撞名） | 5 | 0 |
| 5 | MC / Lean / `web_search` / `git` / `skill` | **0**（未调用） | 本靶是仓内逐字裁决 + 两行谱论；catalog 无匹配技能 | — | — |

**诚实边界**：①`靶 ⟹ H_W` 是本轮的两行推导，纸面**没有**写这条逆向（纸面只把 (5)–(7) ⇒ 聚簇当供给，`fax:5197–5198`）；②Thm 8.6 的 `⟹` 用 Reeh–Schlieder 密度，属纸面证明内部，未独立核 RS 在规范不变子代数上的适用性；③行号以 arXiv v1 抽文为准，未核出版版页码；④`W` 的"新格"（`HZ3-VRATE.md:46`）是**记账格**，本轮的"换名"是**强度格**，两根轴不冲突，本文件不改 HZ3 的账；⑤无 QED，未证之物未标 `[T]`。

*HW1-WFORM · 只新建本文件；既有文件原样；未 `commit`；零 MC。*
