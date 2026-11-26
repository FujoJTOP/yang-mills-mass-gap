# N05-SPLIT — `G″ = Δ_dyn + W_full` 是恒等式还是不等式？（独立核验 T01-IMPLY）

**日期**：2026-11-25 · **只写** 本文件 + `w7/n05_split_audit.py|.log` · **禁 git commit**（遵守）
**被核对象**：`research/v3/species/T01-IMPLY.md:18, 111–131, 208`（判死线 (b)/§4.1/§4.2）里的分解
`G″ = Δ_dyn + W_full`，`W_full ≥ N/12 > 0` ⇒ `G″>0 ⇏ Δ_dyn>0`。
**标签**：`[T]` 真定理 / `[C]` 条件定理 / `[NG]` 对具名类 no-go / `[M]` 精确缺失步 / `[E]` 数值。

---

## §0 判词（先给）

| 问 | 判 | 等级 |
|---|---|---|
| 分解是恒等式还是不等式？ | **恒等式**（两行双线性代数），逐字见 `ch05:230` "The identity is"、`ch18:31` "an identity with no content beyond that display" | `[T]` |
| `G″(0)` 是哪个量的二阶导？ | 粗场条件化**边缘自由能** `V_eff(a)=−log ρ_v(a)` 对粗幅 `a` 的二阶导（逐粗模 `v`、在 `a=a₀`/`ξ=0`）；**不是** `log⟨Ω|O(t)O(0)|Ω⟩`，语料从无以时间关联函数定义 `G″` 之处 | `[T]` |
| `W_full ≥ N/12 > 0` 是定理吗？ | **不是**。已证的只有 `W_geo ≥ (N/12)‖v‖²`（`ch05:214–220`）；`W_full = W_geo + X` 需 `X ≥ 0` ⟺ `Cov(G,Φ′) ≤ 0` ＝ `[NA]` ＝ `ch04` G2.6 **0% 单元**。机器反例：`Φ″≡N/12` 下 **144/144** 点 `W_full < N/12`（min `−0.157`） | `[M]`（T01 标 `[T]` **错标**） |
| 分解有内容吗（AGENTS §2.3）？ | **零内容**。它是**定义**；`G″>0 ⇏ Δ_dyn>0` 是"一个方程两个未知量"的维数计数，**同一义反复** | `[T]` |
| 判死线 | **(a) 触发（恒等式支）**；**(b) 对恒等式不适用，对 `W_full ≥ N/12` 触发（无证明 + 反例）**；**(c) 恒等式：无类限制；地板式：类受限** | — |

---

## §1 `G″(0)` 的逐字定义（哪个函数的二阶导）

`docs/paper-v2/ch05-theorems.tex:230–240`（`eq:split`，label 在 `:239`）逐字：

```
230: \textbf{The exact split and what fails.} The identity is
233: G''(0)&=E[H_{vv}+\Phi''\,|\,\xi=0]
234:         -\operatorname{Var}(G_v+\Phi'\,|\,\xi=0)
235:         =\Delta_{\mathrm{dyn}}+W_{\mathrm{geo}}+X,\\
236: \Delta_{\mathrm{dyn}}&=E[H]-\operatorname{Var}(G),\qquad
237: X=-2\operatorname{Cov}(G,\Phi').
```

`docs/paper-v2/ch15-uvface.tex:20–28`（`eq:v3object`，label 在 `:26`）逐字：`G''(v;a_0)=E_a[∂_v²V]−Var_a(∂_vV)`，
`V=S_W+Φ`，`Φ=−log J`，**conditioned at the coarse field `a=a_0`**。
`docs/paper-v2/ch04-gap.tex:23`（G2.0 恒等层）：`G''(k→0)=E[H|0]−Var(∇S|0)=V_eff''=m_eff²(k→0)`。

**⇒ 定义**：`G″(0)` ＝ **逐粗模 `v` 的粗场条件化边缘自由能曲率**
`G″ = ∂_a²[−log ρ_v(a)]`，`ρ_v` = `⟨v,·⟩` 方向的边缘密度，取 `a=a₀`（`ξ=0`）。
**不是**某时间关联函数：`G = ∂_v S_W` 是**作用力**，`G(a)` 是自由能（`GPP-STRUCT.md:52` 逐字 "可核的定义是 `G″`，不是 `G`"）。
`μ` 上的二阶导与方差取自**同一个条件化测度**，这是后面全部代数的来源。

## §2 `Δ_dyn` 与 `W_full` 的逐字定义

| 量 | 逐字定义 | 出处 |
|---|---|---|
| `Δ_dyn` | `Δ_dyn = E[H] − Var(G)`，`H = Hess S_W` 的 `vv` 元，`G = ∂_v S_W`（条件化于 `ξ=0`/`a=a₀`） | `ch05:236`；`ch15:23–25` |
| `W_geo` | `W_geo := E[Φ″] − Var(Φ′)`，`Φ = −log J`（Haar Jacobian） | `ch05:203`（`eq:wgeo`） |
| `X` | `X = −2Cov(G,Φ′)` | `ch05:237`；`ch15:25` |
| `W_full` | `W_full := W_geo + X`（"the full two-point remainder"） | **`ch08-immune.tex:18`**；`app.tex:796` |
| `(M)` | `Var_μ(⟨v,∇S⟩\|0) ≤ E_μ[H_vv\|0] + W_full − δ`，`(L,U′)`-一致 ⇒ `ch08:19` "the display is then `G″ ≥ δ`" | `ch08:11–19`（`eq:species`） |

**注意 `(M)` 本身是不等式**（`≤`），而 `G″ = Δ_dyn + W_full` 是**另一行**——两者相差一次移项：

```
ch18-scalar-h1.tex:30  Subtracting E[H]+W_full from both sides, W_full := W_geo + X,
       :31             gives an identity with no content beyond that display:
       :33-35  Var_μ(⟨v,∇S⟩|0) − E_μ[H_vv|0] − W_full  =  − G''(0).
```

## §3 ★ 恒等式还是不等式？（核心）

**恒等式。两行证明（双线性）：**

```
Var(G_v+Φ′|0) = Var(G|0) + Var(Φ′|0) + 2Cov(G,Φ′|0)          (方差的二次型恒等式)
E[∂_v²V|0]    = E[H|0] + E[Φ″|0]
⇒  G″ = E[H]+E[Φ″] − Var(G) − Var(Φ′) − 2Cov
      = (E[H]−Var(G)) + (E[Φ″]−Var(Φ′)) + (−2Cov) = Δ_dyn + W_geo + X.  ∎
```

不用模型、不用一致性、不用任何"存在隙"的假设（`ch18:38` 逐字："It uses no model, no uniformity, and no hypothesis that a gap exists"）。
**机器核准**（`w7/n05_split_audit.py` §[1]；12 个随机光滑 `V`，1D 边缘核，`Φ=a²/12`）：
`max |G″ − (Δ_dyn + W_geo + X)| = 6.7e−16` ⇒ 恒等式，**且不需要 bounded 势**（`[E]`）。

**⇒ 判死线 (a) 触发**：分解是**定义**。于是

- `G″ > 0 ⇏ Δ_dyn > 0` 是**直接推论**（一个方程解不出两个符号）：`G″>0` 等价于 `Δ_dyn > −W_full`。
- 机器的**独立见证**（不借用任何 `W_full` 正性）：双井核 `V = −1.100/2 a² + 0.132/4 a⁴`、`Φ=a²/12`（N=2）：
  `G″ = +0.200109`、`Δ_dyn = −0.551151`、`W_geo = +0.027530`、`X = +0.723730`、`W_full = +0.751260`、`resid = −1.1e−16`。
  **`G″>0` 与 `Δ_dyn<0` 精确并存** ⇒ 蕴含假，无需地板式（`[E]`）。

### §3.1 但 T01 的**前提**`W_full ≥ N/12 > 0` 是错的标（这是本轮唯一实质新发现）

- 已证 `[T]`：`W_geo ≥ λ(v) = inf_ξ⟨v, Hess Φ(ξ) v⟩ ≥ (N/12)‖v‖² > 0`，对每个 `a,β,L`（`ch05:214–220`，`eq:wgeo-quant`，Prékopa + `NS-NOSUBTRACTION` Prop. N1）。
- `W_full = W_geo + X`，`X` **符号不定**。故 `W_full ≥ N/12` 需 `X ≥ 0` ⟺ `Cov(G,Φ′|0) ≤ 0` ＝ `ch08:24–27` 的 `(eq:neqassoc2)` ＝ `ch04:37–39` 的 **G2.6，0%**，`ch08:29–30` 逐字 "absent from every known dependence class"。
- **机器反例**（`w7/n05_split_audit.py` §[3]；`Φ″≡1/6=N/12` 逐点成立）：**144/144** 点 `W_full < N/12`，`min W_full = −0.156991`，`X<0` 144/144 ⇒ **`W_full ≥ N/12` 在一般测度下假**，`[M]`。
- 语料已判过同一件事：`H1C-WFULL.md:3–4` 逐字 "`W_full ≥ W_geo` 不成立（扫描点上 21% 失效）⇒ **(a) 的免费 `W_full ≥ (N/12)‖v‖²` 拿不到**"；`:89` 判 `[V]`（无证明）；`:76` 仅 `[E]`：232 点 `min R = W_full/((N/12)‖v‖²) = 1.025`。
- **⇒ T01:129,208 把 `ch05:203–219` 当作 `W_full ≥ N/12` 的 `[T]` 出处，是一次 `W_geo → W_full` 的静默替换**（正是 T01 自己在 §4.2 警示的那类替换，`LOOP-RUNBOOK:4712` 亦带同一错标）。**结论方向不变**（恒等式 + §3 见证已够），但 `[T]` 标签必须降为 `[M]`/`[E]`。
- 附带：T01:129 的措辞"与测度、`β`、条件点都无关的正数"只对 `W_geo` 的 `λ(v)` 成立，对 `W_full` 为假（`X` 依赖测度）。

## §4 判：这个分解有内容吗？（AGENTS §2.3）

**没有。** 三条独立理由，全部语料内逐字：

1. `ch18:31` —— 它是 "an identity with no content beyond that display"。
2. `ch18:59–63`（claimbox 第 2 条）—— "`Var`, `E[H]` and `W_full` are **not separately constrained** by `H1` beyond their appearance in `eq:scalar-identity`; **the sign lives in the curvature**"。
3. 分解**约束 0 个方向**：`Δ_dyn` 与 `W_full` 由 `G″` 与 `X` 反解，恒等式对任意 `C²` 势成立，故对任何符号组合都自洽（§[1] 12/12 残余机器零，符号组合遍布）。

**⇒ `G″>0 ⇏ Δ_dyn>0` 是"`δ` 族成员严格性"的换名，不是新信息。**

## §5 判死线 (a)：逐字写清"`Δ₀>0` 不可能是比 `I2`/`H*` 更强的条件"

`app.tex:182–191`（`app:sharp`）逐字：

```
182: \LM{} ... There is $\Delta_0>0$ such that, uniformly in $(L,a)$ ...,
185: Var_μ(⟨v,∇V⟩|0) ≤ E_μ[H_vv|0] − a²Δ_0²,
190: This is the (M) form ... with W_full = W_geo + X and
191: δ := a²Δ_0² + W_full; the two displays are the same inequality.
```

`ch18:42–47`（`eq:scalar-rule`）：`H1(δ) ⟺ G″(0) ≥ δ`（逐点、逐模、每个 `(β,U′)`），`δ_max = G″(0)`；
`ch18:55–58`：族 `{H1(δ)}` **全序**，`H1(δ) ⇒ H1(δ′)` 对 `δ′ ≤ δ`；"`H1`" without a `δ` **is** `G″(0)>0`。

**⇒ 逐字结论**：

1. `Δ₀>0` 那个靶的算式恰好是 `G″ − W_full ≥ a²Δ₀²`（把 `δ := a²Δ₀²+W_full` 代进 `H1(δ)`），即**把 `G″` 减去一个已知项**——`P03-EQUIV.md:253–258` 已判 `[T]`：`H1 ⟺ G″ ≥ a²Δ₀²+W_full ⟺ Δ_dyn ≥ a²Δ₀²`，`W_full` 两端**同一个对象、直接相消**，"零新估计，只是换名"。
2. 故 `Δ₀>0` **就是** `H1(δ_phys)`，而 `H1`/`I2`/`(E-I2)`/`H*`/`Δ₀>0` 是**同一命题的五个名字**（`paper-v6/sec07-assembly.tex:131`；`sec07:70`）。
3. `G″>0` 是 `δ=0` 成员，由全序（`ch18:55–58`）**严格弱于** `δ_phys>0`；反向不成立。
   ⇒ **`G″>0` 不充分"这句话不含任何超出 `H1(δ_phys) ⇒ H1(0)` 的内容**，因此 `Δ₀>0` **不可能**比 `I2`/`H*` 更强——它是**同一句**，只换了记号与量纲（`G″_a = a²m_phys²`，`ch18:81–83`）。

## §6 判死线 (b)/(c)

- **(b) 对分解不适用**（恒等式无"证明"，只有两行代数；不可反转，因为它是等式）。
  **(b) 对 `W_full ≥ N/12` 触发**：无证明，且**一般测度下假**（§3.1：144/144）。语料只在**格点扫描点**上 `[E]` 成立（`H1C-WFULL:76`，`min R = 1.025`）。**可反转性**：`X` 不定 ⇒ 该式在 `X<0` 大到吃掉 `W_geo` 时反转，已由机器给出。
- **(c) 范围**：
  - **恒等式：无类限制**（任意测度、任意 `C²` 势、只要条件二阶矩有限；1D 边缘核机器核到 `6.7e−16`；不需要 bounded 势、不需要紧群、不需要格点）。
  - **`W_geo ≥ N/12`：类受限** —— `SU(N)` Haar 平图 `X_ℓ=log U_ℓ`、`|d|<2π`（`g″≥1/6`）、Prékopa + 强对数凹在边缘化下不衰减（`ch05:196–228`）。这是**唯一**带"类"的部件。
  - **`W_full ≥ N/12`：无已证类**（`[M]`）。

## §7 附带记账（不属判死线，但须留痕）

1. **过期标签**：T01:18,118 引 `eq:h1st-dyn`；T01:113 引 `eq:h1st-rule`。全仓 `label{eq:` 检索：`ch15-uvface.tex` 只有 `eq:v3object/v3good/...`（`:26,:59,...`），`ch18` 只有 `eq:scalar-identity/eq:scalar-rule`（`:36,:47`）。**`eq:h1st-dyn`/`eq:h1st-rule` 在现行 v2 文件里不存在**（`paper-v5/sec02:100` 仍引 `ch15:60–68 (eq:h1st-dyn)`，同类过期）。指向的数学内容无争议，但标签须改 `eq:scalar-rule` + `app:191`。
2. **T01:139/208 的 `[T]`** 应为 `W_full ≥ N/12` 的 `[M]`；`G″=Δ_dyn+W_full` 的 `[T]` **成立**（恒等式）。

## §8 检索执行核验（AGENTS.md §1）

| 项 | 记录 |
|---|---|
| 工具状态 | `read` **OK**（`ch05/ch15/ch08/ch18/app/ch04` 逐行读到）；`glob` **OK**（`**/T01-IMPLY.md` → 1）；`pwsh` **OK**；`py -3.12` **OK**（脚本正常退出）。无 `ENOENT`、无 `ripgrep provider failure`、无超时。 |
| 检索式 1 | `Select-String -SimpleMatch 'W_{\mathrm{full}}','W_full',...`（全仓 `*.tex,*.md`）→ **命中大量（输出截断，完整落盘 `…dsh-spill-…/027b3336d075-pwsh.txt`）**，**记 TRUNCATED 不记 0** |
| 检索式 2 | `Select-String -SimpleMatch 'h1st-dyn','h1st-rule','Delta_0'` → `Delta_0` 命中 40+ 行；`h1st-dyn` **0 条定义处**（仅 `paper-v5/sec02:100` 引用） |
| 检索式 3 | `Select-String -SimpleMatch 'eq:h1st','eq:scalar-rule','label{eq:'`（ch15/ch18/ch13）→ **ch15 15 条 / ch18 10 条 / ch13 8 条 label**，其中 `eq:h1st*` **0 条** |
| 未查到（采信） | **"`W_full ≥ N/12` 的证明"**：检索式 1 的并集 + `H1C-WFULL.md` 通读，**0 处证明**，只有 `[E]`（232 点）与 `[V]` 判词（工具 OK） |

## §9 出处（逐条可核）

`docs/paper-v2/ch05-theorems.tex:196–228, 230–249` · `ch15-uvface.tex:20–32, 55–70` · `ch08-immune.tex:11–30` ·
`ch18-scalar-h1.tex:28–48, 50–70, 77–83` · `ch04-gap.tex:23, 33–39` · `app.tex:179–202, 796` ·
`research/v3/species/P03-EQUIV.md:11–14, 165–193, 253–258` · `H1C-WFULL.md:3–5, 30–46, 76, 89–101` ·
`GPP-STRUCT.md:52, 144` · `LOOP-RUNBOOK.md:4712–4723` · `paper-v6/sec07-assembly.tex:70, 131–132`。
**本轮机器证据**：`research/v3/w7/n05_split_audit.py` → `n05_split_audit.log` §[1]/§[2]/§[3]。

---

*本文件不证新定理、不改任何数值、不 QED。唯一新增的**判决性陈述**是：分解是**恒等式**（故 `G″>0 ⇏ Δ_dyn>0` 为零内容同义反复），
而其伴随不等式 `W_full ≥ N/12` 是 **`[M]` 且一般测度下假**（T01 的 `[T]` 为错标）。**Clay 仍 0%。***
