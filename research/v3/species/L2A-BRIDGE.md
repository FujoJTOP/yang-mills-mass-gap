# L2A-BRIDGE — `残数 ⇒ 曲率缺陷` 的精确陈述（L2 第一步）

> 2026-11-25 · 第 L2A 轮 · 物种笔记 · **只写本文件；未 git commit**
> 任务源（逐字）：`N11-FALLBACK.md:10128` 一行 ——「L2 桥定理：`β<1/48` 的 `f₁² ≥ c`、`G″ = S−D` 已 `[T]`，缺"**残数 ⇒ 曲率缺陷**"引理；难度中」。
> 逐字读过：`docs/paper-v2/ch05-theorems.tex:230–249`（`eq:split`，label 在 `:239`）· `docs/paper-v2/ch18-scalar-h1.tex:1–90`（`eq:scalar-identity`/`eq:scalar-rule`）·
> `docs/paper-v2/ch02-squeeze.tex:42–130,179–228`（squeeze、`prop:thinatom`）· `docs/paper-v2/app.tex` 的 `Δ₀` 块（`:182–198`）与 `D/S` 对象勘误（`:236–238`，`:681–683`）·
> `research/unconditional/INVENT-VARIANCE-INEQ.md`（全文 34 行）· `research/unconditional/H1-DIRECT-SPLIT.md:60–130` · `research/unconditional/KEY-DIRECT-DERIVATION.md` ·
> `research/v3/species/DS-RATIO.md`（判官文，D14/D17/D10）· `research/v3/H2-DEFECT-BOUND.md:105–110` · `research/v3/Q3-SCHUR-DEFECT.md:55–80,210–250` ·
> `research/v3/species/T01-IMPLY.md:17–19,64–123,150` · `research/unconditional/PROOF-P3-UNCONDITIONAL.md:1–45` · `docs/YANG-MILLS-MASTER.md:40` · `lean/YMStatement/Reduction.lean:58–80`
> 标签：`[T]` 定理/定义级 · `[C]` 条件 · `[NG]` no-go 或缺步 · `[O]` 精确缺失步 · `[E]` 数值 · `[J]` 判断 · `[V]` 已存档读数 · `[X]` 语料内冲突
> **No QED；Clay 0%。**

---

## §0 判词（先结论）

| 问 | 判 | 标签 |
|---|---|---|
| ① `S`、`D` 各是什么 | **两套读法，`G″` 同、`D` 不同**（§1）。`eq:split` 第一行：`S=E[h_n|0]`、`D=Var(f|0)`；SCBL 第二行：`S=E[h_n−⟨K⁻¹c,c⟩|0]`、`D=Var(f|0)−E[⟨K⁻¹c,c⟩|0]` | `[T]` |
| ② `D` 是 `Var(G)` 吗 | **只有"未做 Schur 移项"的读法是**。语料默认的 `D`（`D≤0⟺K⪰0` 那条符号律所属的读法）是 `Var(f|0) − E⟨K⁻¹c,c⟩`，**不是方差** | `[T]` |
| ③ "曲率缺陷"是什么 | 语料的 `D/S≈2%` **不是** `D/S`，是 `1−D/S = G″/S`（`DS-RATIO:25` 逐字）。β=2.5 处真值是 `D/S≈0.97`、`G″/S≈2–3%`；且 `D/S` 是**分裂商**，随分裂取 `0.02/0.34/0.97/∞` | `[T]`+`[V]` |
| ④ (i)/(ii)/(iii) 哪个对 | **三个是同一个**（由恒等式等价，(i)≡(ii)≡(iii)），且 **`app.tex` 逐字判它"insufficient"** | `[T]` |
| ⑤ 判死线 | **(c) 主判 + (b) 的精确陈述**；不是 (a)。见 §5 | `[T]` |

**一句话**：`残数` 是**分子**，`曲率缺陷` 是**分母的倒数**；squeeze 把它们**乘起来**，没有把前者**蕴含**成后者。所以"残数 ⇒ 曲率缺陷"作为**蕴含**不成立（三个反例，§6），作为**乘法桥**成立且已 `[T]`（§4）。

---

## §1 逐字核 `G″ = S − D`

### 1.1 论文原式（`[T]`，`ch05-theorems.tex:230–240`，`eq:split` 逐字）

```
G″(0) = E[H_vv + Φ″ | ξ=0] − Var(G_v + Φ′ | ξ=0) = Δ_dyn + W_geo + X ,
Δ_dyn = E[H] − Var(G) ,      X = −2Cov(G, Φ′) .
```

记号陷阱（先钉）：`G_v := ∂_v S_W` 是**力**，`G(a)` 是**自由能**，同一字母两物（`GPP-STRUCT.md:54`）。

### 1.2 SCBL 分裂（`[T]`，`INVENT-VARIANCE-INEQ.md:6–8` 逐字）

`μ ∝ e^{−V}dX`，`V=S_W+Φ`，`a=⟨v,X⟩`，`h_n=vᵀHess V v`，`f=∂_nV=G_v+Φ′`，`Π₀=I−vvᵀ`，`c=Π₀ Hess V v`，`K=Π₀ Hess V Π₀`：

```
G″_0 = S(v) − D(v)
S := E_μ[h_n − ⟨K⁻¹c,c⟩ | 0]        （Gaussian–Schur 曲率；Gaussian 下 S = 1/Var_μ(a) 精确）
D := Var_μ(f|0) − E_μ[⟨K⁻¹c,c⟩|0]   （条件化纤维 μ₀ 的 Brascamp–Lieb 缺陷）
符号律：D ≤ 0 ⟺ K ⪰ 0；D = 0 ⟺ μ₀ Gaussian；D > 0 ⟺ K 不定（Savvidy）
```

### 1.3 三种读法（同一 `G″`，不同 `D`）

| 读法 | `S` | `D` | `D` 是方差吗 | 出处 |
|---|---|---|---|---|
| A（`eq:split` 第一行，未 Schur） | `E[H_vv+Φ″\|0]` | `Var(G_v+Φ′\|0)` = `Var(力\|0)` | **是**（且 `≥0`） | `ch05:233–234` |
| B（SCBL，`∂_n` 归一） | `E[h_n−⟨K⁻¹c,c⟩\|0]` | `Var(f\|0) − E⟨K⁻¹c,c⟩\|0` | **不是** | `INVENT:8`；`H1-DIRECT-SPLIT:62` |
| C（`Δ_dyn` 块 / Lean） | `E[H]` | `Var − W`（`W=W_geo+X`） | **不是** | `Reduction.lean:71`（`gDouble_eq_S_sub_D`，`ring`）+ `ch05:236–237` |

**`D` 是 `Var(G)` 吗 —— 答案**：`G″ = E[H] − Var + W`（Lean `gDouble` def，`:66`）⇒ 只有当把 `W_full` 整个塞进 `S` 时 `D = Var`；语料默认读法 B/C 里 `D = Var − (Schur 损失或 W_full)`，**不是方差**。三种读法的 `S−D` 逐字相同（Schur 项抵消，`KEY-DIRECT-DERIVATION:121` 逐字 "the discrepancy decomposition is an identity that cancels its `K⁻¹` term"）⇒ **`D` 是分裂商（split quotient），不是坐标**（`DS-RATIO` 判词；`:234` 逐字 "`D/S` 的数值由分裂决定（`0.02/0.34/0.97/∞` 皆可），而 `S−D=G″` 与分裂无关"）。

---

## §2 "曲率缺陷"是什么（必须先把两个数分开）

| 量 | β=2.5 的语料值 | 性质 | 出处 |
|---|---|---|---|
| `D/S`（SCBL/H2 分裂、胶子模、β=2.4、L=8） | **0.9651**（`0.8311–1.1028`，`5/12` 模 `D>S`） | 分裂商 | `H2-DEFECT-BOUND:108,292`；`SM2-RANK-ONE:324` |
| `1 − D/S = G″/S` | **≈0.035**（另一路给 **0.02**） | **这才是语料说的"2%"** | `DS-RATIO:25,233` 逐字；`app.tex:236–238` |
| 绝对值 `G″` | `+0.049`（mean）；`+0.27…+0.51`（另三处） | **与分裂无关** | `SM2:324`；`W7`/`H1-NUMERIC-RIGHT:31` |

- **`D/S≈2%` 是记账事故**：它是 `1−D/S`（`DS-RATIO:25` 逐字："`2%` 不是 `D/S`，是 `1−D/S`；`D/S≈0.97` 意味着 `D` 与 `S` 近抵消"）。⇒ 所谓"因子 50 余量"应为 **1.03×**。
- 精确可解模型给闭式 `D/S = β²/β*²`，`β* = 4.30661`（`DS-RATIO:18–21`，`[T]`+`[E]`）：**单调、无上界、在 `β*` 精确穿越 1**，而 `β*` 同时是 `G″` 的零点 ⇒ `D/S<1 ⟺ G″>0`（`S>0` 时）**是同义反复**（D17）。
- ⇒ **判**：「曲率缺陷」= `D`（`G″` 相对 `S` 的亏量），可用；但 **`D/S` 与 `1−D/S` 会互相冒充**，所以任何引理必须写在**绝对值 `G″ ≥ λ`** 上，不许写比值。

---

## §3 三个候选：它们不是一个问题，是一个不等式

由恒等式 `G″ = S − D`（任一读法）：

```
(i)  D ≤ (1−c′)S   ⟺  G″ ≥ c′S
(ii) D/S ≤ φ(c)    ⟺  G″ ≥ (1−φ(c))S ,  即 φ ≡ 1−c′ 时与 (i) 逐字相同
(iii) G″ ≥ c″S     ⟺  D ≤ (1−c″)S ,     即 c″ ≡ c′ 时与 (i) 逐字相同
```

**⇒ (i) ≡ (ii) ≡ (iii)。** 三个候选不是一个三选一，是同一个**相对**不等式换三个字母。

而语料对这一个不等式的判词是**逐字的**（`app.tex` 的 `Δ₀` 块）：

> "Equivalently `D ≤ S − λ` **together with** `S ≥ λ` (the bare defect form `D ≤ (1−ε)S` is **insufficient**)"
> —— 独立同判：`H1-DIRECT-SPLIT:109`：「`D ≤ (1−ε)S` … *necessary-looking, insufficient.* … closing needs the **pair** `{S≥λ_S, D≤S−λ}` … There is **no factorization of H1** into "denominator large" + "defect small"」

缺的正是**绝对地板 `S ≥ λ`**：`S` 无一致下界（自由场 `S=0`；`DS-RATIO:21,253` 给 `θ=π,a₀=1,β=4.30661` 处 `S=−0.0339`；`k→0` 高斯块 `S=0.198→0.058→0.016→0.004`，`S1-CONDITIONAL-POINCARE:344`）。⇒ 相对形式在 `S→0` 处**恒真而无用**。

---

## §4 ★ 精确陈述（正确读法：乘法桥，不是蕴含桥）

**设定**：有限格点 Λ，`β>0`，粗方向 `v`（`‖v‖=1`），`a=⟨v,X⟩`，`μ ∝ e^{−V}dX`，`V=S_W+Φ`；`S`、`D` 如 §1.2；`χ := Var_μ(a_{k→0})`；`f₁² := ν({m₁²})`，`ν` 为 `trF²` 的 Källén–Lehmann 测度。

**(BL-0) `[T]` 恒等式。** `G″(0) = S − D = E_μ[h_n|0] − Var_μ(f|0)`。（`ch05` `eq:split` + `INVENT §1`）

**(BL-1) `[T]` 双端 squeeze**（`ch02-squeeze.tex:51` 逐字）：
```
f₁²/χ  ≤  inf supp ν  ≤  M₀/χ .
```
（下界证明两行：`χ m₁² = f₁² + m₁²∫_{s>m₁²}s⁻¹dν ≥ f₁²`，`ch02:55–56`。)

**(BL-2) `[T]` 装配**（`ch02:199–201` 逐字）：
```
(f₁² ≥ r₀) ∧ (χ ≤ 1/λ)  ⟹  m₁² ≥ f₁²/χ ≥ r₀λ > 0 .
```

**(BL-3) `[T]` 分母 = 曲率**（`D3-SUSCEPTIBILITY-VARIATIONAL` 三 face；`ch18` `eq:scalar-rule`）：
```
χ ≤ 1/λ  ⟺  G″(0) ≥ λ  （固定 L,a；与 H1(λ) 逐字同一句）
```

**⇒ 唯一成立的"桥"是三件套（不是两件）：**
```
L2A:   f₁² ≥ c   ∧   G″ ≥ λ      ⟹   m₁² ≥ cλ > 0 .        [T]（BL-1+BL-2+BL-3）
```
**而「`f₁² ≥ c ⟹ G″ ≥ λ`」——即 (i)/(ii)/(iii)/N11-L2 的翻译引理本身——不在上述任一条的推论里**，理由不是难，是**结构**：`f₁²` 进的是**分子**，`χ`（`=1/G″`）进的是**分母**，squeeze 是乘法的，等价于 `m₁² ≥ f₁²/χ`；它对 `χ` 只给**下界方向的约束**（要 `χ≤1/λ`，即 H1），对 `G″` 一个字都不给。

**N11-L2 写下的目标式 `S − D ≥ a²Δ₀(c)²` 本身不成立为陈述**：`Δ₀` 是 `app.tex:182–194` 里**那条不等式自身的见证变量**（`∃Δ₀>0 s.t. Var ≤ E[H] − a²Δ₀²`），`N11 §0` ① 已自判"`Δ₀>0 ⟺ H1 ⟺ (M)`，是同一句的换名"；且该命名**已被撤回**（`docs/paper-v6/sec07-assembly.tex:70`、`sec06-kl.tex:120`）。⇒ 不存在映射 `c ↦ Δ₀(c)`；桥的正确书写形式只有 `G″_{k=0}/a² ≥ δ² > 0` 这一种。

---

## §5 判死线

### (a) 不成立 —— 不由 `I1` 直接给出

三条独立证据（任一足够）：

1. **P3 自己逐字否认**（`PROOF-P3-UNCONDITIONAL.md §5`）："The polymer expansion controls a **weight** (I1), not **eigenvalue spacing** (I2)"。桥正是 weight→spacing。
2. **thin atom `[T]`（机检）**：`ch02:112–129` `prop:thinatom`——对任意 `λ,M₀` 存在 `ν = wδ_ε + (1−w)δ_{s₂}`，`w = ε/(2λ)`、`inf supp = ε`。取 `w ≥ c` 固定、`λ→∞` ⇒ `f₁² ≥ c` 成立而 `G″ = 1/χ → 0`。Lean：`H1_blind_to_thin_atom`（编译通过，`ch02:164–170`）。
3. **`χ` 无界即死**（`ch02:226–228` 逐字）："If `f₁² ≥ r₀` but `χ` is unbounded, then `m₁²` can be arbitrarily small: the squeeze only gives `m₁² ≥ f₁²/χ`, and `χ→∞` kills the bound."

### (b) 成立 —— 需要的新输入**不是**一条新不等式

新输入是 **(BL-3)**：`G″ ≥ λ`，或等价的 `{S ≥ λ, D ≤ S − λ}` 对。而这两者与 `H1(λ)` **逐字同一句**（`ch18` `eq:scalar-rule`；`H1-DIRECT-SPLIT:109` 的 "no factorization"）。⇒ 若坚持"桥"这个词，它是 **`I1 ∧ H1 ⟹ m₁² ≥ cλ`**（`[T]`，`ch02:199`），**不含任何关于 `D` 的新估计**。

### (c) 主判 —— "残数 ⇒ 曲率缺陷"真实说的是什么

**它说的是两件 `[T]` 的事，都不是蕴含：**
1. **同源**：`f₁²` 与 `G″=1/χ` 是同一 KL 测度 `ν` 的两端（原子权重 vs 逆磁化率），`D3` 三 face 逐字给出这条绑定；
2. **单向乘法**：`m₁² ≥ f₁²/χ`，即残数给隙一个**只在分母有上界时才有效**的下界。
**它不说的**：残数给曲率任何正下界。⇒ 「残数 ⇒ 曲率缺陷」应改写为 **「残数 × 磁化率上界 ⇒ 隙」**；把 `残数` 单独放在左边是一个**类型错误**（分子当分母）。

---

## §6 三个反例：残数与曲率**逻辑独立**（两向都断）

| # | 模型 | `f₁²` | `G″` | 杀了哪个方向 | 标签/出处 |
|---|---|---|---|---|---|
| 1 | thin atom `wδ_ε+(1−w)δ_{s₂}` | `=w>0` 可固定 | `→0`（`χ→∞`） | **`R ⇒ C` 假** | `[T]` `ch02:112–129`；`T01-IMPLY:77` |
| 2 | 有质量自由场 `:φ²:` | `=0`（连续谱，无原子） | `=m²>0`（`D=0,S=m²`） | **`C ⇒ R` 假** | `[T]` `INVENT §4`；`SYNTHESIS-CROSSROUTE:12` 逐字 |
| 3 | `β=0` 纯 Haar | `=0` | `≥S≥N/12>0` | **`C ⇒ R` 假**（且**过不了纯 Haar 酸测**） | `[T]` `INVENT §4`；`TCI-COMPENSATION:314`；`T01-IMPLY:150` |

自由场（`f₁²=0` 且 `G″=0`）过酸测 `[T]`（`INVENT §4`），但**它过得太容易**——真正有分辨力的是 #2/#3：**曲率为正而无残数**。⇒ 任何"桥"必须解释为什么 #2/#3 不算反例，而它们正是 `[T]`。

---

## §7 为什么 L2 的"难度：中"估价不成立（诚实线）

N11-L2 的承诺是"**只承诺窗口内**"（`N11:67`）。但在该窗口内：

1. **结论已有 `[T]`，不是本桥给的**：强耦合 `m>0` 是 OS'78 / SZZ'23（P3 §5 逐字："**Does NOT** give … the strong-coupling *gap* itself (that is Osterwalder–Seiler '78 / SZZ '23, not P3)"）；
2. **`β=0` 的曲率地板已 `[T]`**：`V=Φ` 凸 ⇒ `D≤0` ⇒ `G″ ≥ S ≥ N/12`（`INVENT §4`、`ch05:406`）；
3. **强耦合组合的判词已写在语料里**（`SYNTHESIS-CROSSROUTE:12` 逐字）："… yields one genuine **small theorem** — a-uniform marginal log-concavity at `β<1/48` — but it produces the **gap** (already SZZ), **never the residue**"；
4. 而窗口外（`β≈2.5`），`D>0 ⟺ K` 不定，`D<S` 逐字 = H1（`INVENT §5`）。

⇒ L2 只有两种结局：**(α)** 在窗内真但空转（结论已被 SZZ/β=0 地板覆盖）；**(β)** 一旦要求搬运到窗外，要证的正是 H1，不是"初等估计"。**"中难度"是把 (β) 记成了 (α)。**

---

## §8 检索记录（`AGENTS.md §1` 四项）

| # | 工具 | 检索式 | 范围 | 命中 | 判 |
|---|---|---|---|---|---|
| 1 | `grep` | `eq:split\|G″ = S\|G''=S` | 全仓 | 250+（截断，存盘 `dsh-spill … b2b4f72eac6b`） | OK |
| 2 | `grep` | `D/S\|η ?= ?D\|eta ?= ?D/S` | `*.md` | 250+（含 `DS-RATIO`、`H2-DEFECT-BOUND`、`app.tex` 引文） | OK |
| 3 | `grep` | `S ?≥ ?N/12\|S ?= ?1/Var` | 全仓 | 162 | OK |
| 4 | `grep` | `1/48.{0,80}G\|G.{0,80}1/48` | 全仓 | 227（**未见**"`β<1/48` 上 `G″≥c` 的定理"；窗口内曲率只有 `β=0` 地板与 KL 改名） | OK |
| 5 | `grep` | `residue 与 Hess 同现`（承 `N11 §5` #2–#4 的 0 命中） | 本文件未重跑 | — | 承 `N11:144–146` |
| 6 | `read` | 逐字核 8 个文件（见文件头） | — | — | OK |
| 7 | `write` | 本文件 | — | 1 | OK；**未 commit** |

**读法**：#4 的"未见"是**检索成功后的 0 命中语义**（工具 OK），不是 TOOL-FAIL。**未做**：`web_*`（本轮不需要外部文献）。

---

## §9 标签汇总 · 冲突 · 自否证

**标签**：§1 三读法 `[T]`（逐字，行号已给）；§2 `D/S≈0.97`、`1−D/S≈2–3%` `[T]`(身份)+`[V]`(读数)；§2 闭式 `D/S=β²/β*²` `[T]`+`[E]`；§3 (i)≡(ii)≡(iii) `[T]`（代数）、"insufficient" `[T]`（`app.tex` 逐字）；§4 BL-0–BL-3 `[T]`；§4 L2A `[T]`；§5(a) 三条 `[T]`；§5(b) "no factorization" `[T]`；§6 三反例 `[T]`；§7 `[J]`（判，非定理）。

**冲突登记（不自行改判）**：
- **(X1) 窗口数字**：任务书 ① 写 `β < 0.018394`；`docs/YANG-MILLS-MASTER.md:40` 写窗口 `β < 3.53×(1/48) = 0.073576`（`1/48` 只是基半径）。两说并存，本文件不改判，只记 `[X]`。
- **(X2) `D/S` 的 β-稳健性**：`H1-NUMERIC-RIGHT:45` 逐字 "`D/S` is NOT a clean `2%` … `−0.45…+0.60` at strong coupling … the `2%` is a *weak-coupling artifact*"，与 `SM2`/`H2` 的 `0.97` 并存。⇒ §2 的读数**只在 `β=2.4–2.5`、胶子模、`L≤8` 上成立**（窗口限定，`AGENTS §2.1`）。
- **(X3)** `P3` 的 `[T]` 曾在 `MASTER:42` 被降级、2026-11-24 已解除降级（`C̄≡N`）。当前状态 `[T]`；本文件不依赖其窗口以外的任何东西。

**自否证条件（本文件被打到即作废）**：
1. 若语料里存在形如 `f₁² ≥ c ⟹ D ≤ (1−c′)S` 或 `⟹ G″≥c″S` 的**已证**命题 ⇒ §5 的 (a) 撤销，L2 变"整理"。
2. 若 `D/S` 能被证明是**分裂无关**的坐标（与 `DS-RATIO:234` 相反）⇒ §2/§3 重写。
3. 若 thin-atom 命题或 `H1_blind_to_thin_atom` 的机检被推翻 ⇒ §6 #1 撤销。
4. 若 `β<1/48` 上出现 `G″ ≥ c > 0` 且 `c` **不**依赖 β→0 地板（即真搬运了曲率而非换名）⇒ §7 的"空转"判撤销。

**未做/未写**：未改 `position-current.html`（本轮指令为"只写一个文件"，故 `AGENTS §6` 的更新义务本轮**不履行**，只在此记账）；未 git commit；未跑数值。
