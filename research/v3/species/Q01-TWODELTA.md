# Q01-TWODELTA — `δ` 到底有几个？`δ_curv` 与 `δ_I2` 是不是同一个符号？

**日期**：2026-11-25 · **只写** 本文件 · **禁 git commit**（遵守）
**被核对象**：`N05-SPLIT.md:119–122`（"`H1/I2/(E-I2)/H*/Δ₀>0` 是同一命题的五个名字"）vs `N12-VERDICT.md:67–74`（"两个 `δ` 同名不同物 `[T]`"）vs `N06-RELATIONS.md:13,50,145`（"`Δ₀` 与 `I2` 不可比"）
**标签**：`[T]` 真定理 / `[C]` 条件 / `[NG]` 对具名类 no-go / `[M]` 精确缺失步 / `[E]` 数值

---

## §0 判词（先给）

| 问 | 判 | 等级 |
|---|---|---|
| `ch18:43` 的 `H1(δ)⟺G″≥δ` 里的 `δ` 是哪个？ | **`δ_curv = a²Δ₀²+W_full`**。`ch18:40` 同一子句逐字给出 `Var ≤ E[H]+W_full−δ`，与 `app:191` 的 `δ` 定义移项后逐字同式 | `[T]` |
| `app:789` 的 I2 `Λ★−m₁² ≥ δ` 里的 `δ` 是哪个？ | **质量² 槽的层距**，与 `m₁²`/`Λ★` 同量纲；**与 `:191` 的 `δ` 不是同一物** | `[T]` |
| 两个 `δ` 是同一个吗？ | **不是**。量纲判据：`app:191` 的 `δ` 是**格点单位无量纲曲率地板**（`ch18:81–83` 逐字 "every δ in this chapter is a dimensionless curvature floor and **none is a physical mass bound**"）；`app:789` 的 `δ` 坐落在 `Λ★−m₁²` 里 ⇒ **质量²** | `[T]` |
| 语料混用同一符号吗？ | **混用**（同一 `app.tex` 内 `δ` 三义：`:191` 曲率 / `:787` Dirac / `:789` 层距；`ch05` 内四义）。符号表 **无 `δ` 条目** ⇒ 无全局定义可依 | `[T]` |
| `n05` 的"⟺ I2"是否把两个 `δ` 当一个？ | **是**。`N05:120` 的桥只有一句"五个名字"，且**该引文在现行语料里已被撤回**（`sec07:70` "That identification is withdrawn"、`sec07:132` "同名不同物"）⇒ **误并** | `[T]` |
| 判死线 | **(c) 成立，其内容正是 (b) 对 `(δ_curv, δ_I2)` 这一对的判定**：两 `δ` 确为不同物 ⇒ `n05` 的 `⟺ I2` **无效**；`n12` 的"同名不同物" **对但不全**（全仓 `δ` 至少 6 义） | — |

**一句话**：`P03:253` 的 `H1 ⟺ G″ ≥ a²Δ₀²+W_full ⟺ Δ_dyn ≥ a²Δ₀²` 在**曲率槽内**是 `[T]`（`W_full` 两端同一对象相消），到 `I2` 的那一步**没有推导，只有一个共用字母**。

---

## §1 `δ` 的全部出现处（逐行，三文件，共 52 行）

检索执行：`grep -n "delta"`（工具 `grep`，**OK，返回带行号正文，无 provider failure / 无超时**）
- `docs/paper-v2/ch18-scalar-h1.tex` → **41 行**
- `docs/paper-v2/app.tex` → **3 行**
- `docs/paper-v2/ch05-theorems.tex` → **8 行**
- 合计 **52 行**（部分行含 2 个 `δ` token，如 `ch18:43`）

### 1.1 `ch18`：41 行，**全部** `δ_curv`（无量纲曲率地板）

| 行 | 原文片段 | 归类 |
|---|---|---|
| 10 | `$H1(\delta)$ \emph{is} the scalar positivity $G''(0)\ge\delta$` | 曲率 |
| 11 | `$\delta_{\max}=G''(0)$` | 曲率 |
| 15 | `$\delta$ \emph{increases} with $\beta$` | 曲率 |
| **40** | `$\operatorname{Var}\le E[H]+W_{\mathrm{full}}-\delta$` | **曲率（定义式）** |
| **43** | `H1(\delta)\iff G''(0)\ge\delta` | **曲率（被问的那一行）** |
| 46 | `\boxed{\ \delta_{\max}=G''(0).\ }` | 曲率 |
| 55–57 | `the family $\{H1(\delta)\}_\delta$ is totally ordered … $H1(\delta)\Rightarrow H1(\delta')$ for $\delta'\le\delta$ … ``$H1$'' without a $\delta$ is the statement $G''(0)>0$` | 曲率（全序族） |
| 65 | `must produce $G''(0)\ge\delta$ on the whole $(\beta,\text{vacuum})$ plane` | 曲率 |
| 73–74 | `the display is then $G''\ge\delta$ in curvature form`；`H1A-EXACT.md \S1.3 (identity, $\delta_{\max}$…)` | 曲率 + 引文 |
| **81–83** | `the curvature floor is a lattice-unit object, $G''_a=a^2m_{\mathrm{phys}}^2$, so **every $\delta$ in this chapter is a dimensionless curvature floor and none is a physical mass bound**` | **★ 语料自证** |
| 114 | `\delta(\beta)=G''(0)=W_0-\beta^2V_0` | 曲率（单 plaquette 实例） |
| 118 | `\label{eq:scalar-plaq-delta}` | 标签 |
| 122 | `\delta(\beta)<0\iff\beta>\beta^*=4.30661` | 曲率 |
| 140 | `$\delta_{\max}=\operatorname{Var}-E[H]-W_{\mathrm{full}}$` | 曲率 |
| 183 | `$\delta<0\iff$ the vacuum does not sit at the conditioned point` | 曲率 |
| 192 | `$\beta<\beta^*$ and $\delta=+0.136>0$` | 曲率（`[E]`） |
| **195–196** | `the $\delta$ of \eqref{eq:scalar-plaq-delta} is a **lattice-unit curvature with no $a\to0$ scaling performed**` | **★ 语料自证** |
| 216, 220, 224 | `$\delta$ is read off the identity`；表头 `Model & $\delta$ & …` | 曲率（表） |
| 239 | `exact zero, $\delta_{\max}=0$` | 曲率 |
| 253 | `exactly, so $\delta=W_{\mathrm{full}}$` | 曲率 |
| 271 | `\emph{a positive $\delta$ is a ceiling, not an atom}` | 曲率 |
| 281, 323 | `\eqref{eq:scalar-plaq-delta}` | 交叉引用 |
| 308 | `$\delta=+0.136$` | 曲率（`[E]`） |
| 329 | 表头 `$\beta$ & $\delta$ & $E[H]$ & $\operatorname{Var}$` | 曲率（表） |
| 337–339 | `$\delta$ is monotone \emph{increasing} in $\beta$ … $V_0(3)<0$ in the notation of \eqref{eq:scalar-plaq-delta}` | 曲率 |
| 342 | `$\delta$ from` [SU(4) lift] | 曲率 |
| 367–368 | `$\delta$ grows with $\beta$ … $\delta>0$` | 曲率 |
| 605, 609, 614 | `$\delta_{\max}=G''(0)$`；`$\delta(\beta)=0.2051153-0.0110593\beta^2<0\iff\beta>\beta^*$`；`$\delta$ \emph{increasing} in $\beta$` | 曲率 |
| 630 | `No $a\to0$ scaling is performed on any $\delta$ of` … | 曲率（限定） |

**⇒ `ch18` 41/41 行是同一个 `δ_curv`，且 `ch18` 自己在 `:81–83` 声明它是无量纲曲率、**不是**质量界。**

### 1.2 `app`：3 行，**三义并列**（这是记账错误现场）

```
191: $\delta:=a^{2}\Delta_0^{2}+W_{\mathrm{full}}$; the two displays are the same inequality.
787: \item[$\Lambda_\star$] the isolation threshold $\inf\operatorname{supp}(\nu-f_1^2\delta_{m_1^2})$.
789: \item[I2] the isolation statement $\Lambda_\star-m_1^2\ge\delta>0$.
```

- `:191` → `δ_curv`（无量纲曲率地板）
- `:787` → **`δ_{m_1^2}` = 质量² 处的 Dirac（测度原子）**，第三义
- `:789` → `δ_I2`（质量² 层距），第二义

**★ 符号表核对**：`app.tex:780–797` 逐条列出 `ν, m_1², f_1², χ, M_k, Λ_★, I1, I2, H1, H^*, G'', S/D, K, W_geo, W_full, L3` —— **没有 `δ` 条目**。故 `:191` 与 `:789` 两个 `δ` 在符号表里**都没有全局定义**，相隔 598 行同字母不同物。[T]

### 1.3 `ch05`：8 行，**四义**

| 行 | 原文片段 | 义 |
|---|---|---|
| 74 | `an $O(e^{-\delta t})$` | **衰减率**（量纲 = 质量¹） |
| **151** | `the correct I2 is a gap to the excited glueball of size $\delta\approx0.5\,m_1$, fixed by power counting` | **层距的"质量"版**（`0.5m_1`，量纲 = 质量¹；与 `app:789` 的 `Λ★−m₁²` 差一次方） |
| 222 | `$g(x)=\int f\,\delta(x-\langle v,\xi\rangle)\,d\xi$` | **Dirac** |
| **423** | `$\Lambda_{\star}:=\inf\operatorname{supp}(\nu-f_1^2\delta_{m_1^2})>m_1^2$` | **Dirac** + 层距（`>m_1^2`）同屏 |
| 436 | `$|\eta(t)|\le Be^{-\delta t}$` | **衰减率** |
| 440 | `Write $\nu=f_1^2\delta_{m_1^2}+\mu$` | **Dirac** |
| 454 | `$\Delta/\Lambda\ge\min(\delta_W,\delta_S)>0$` | **O(1) 常数对** |
| 456 | `$\inf\operatorname{spec}(H_a|_{\Omega^\perp})\ge\delta\Lambda$` | **无量纲 × Λ** |

**⇒ `δ` 在现行三文件里至少 6 义**：`①` 曲率地板（无量纲，`app:191`+`ch18` 全 41 行）、`②` 层距-质量²（`app:789`）、`③` 层距-质量¹（`ch05:151`）、`④` Dirac 原子（`app:787`，`ch05:222,423,440`）、`⑤` 衰减率（`ch05:74,436`）、`⑥` O(1) 常数/无量纲（`ch05:454,456`）。

---

## §2 `ch18:43` 的 `δ` 是 `δ_curv` —— 逐字三步闭合

1. `ch18:38–41` 逐字：`Substituting the curvature form \eqref{eq:opencell-curv} of Chapter~\ref{ch:gap} --- which is $\operatorname{Var}\le E[H]+W_{\mathrm{full}}-\delta$ --- and moving terms again yields the decision rule`
2. `app:185–191` 逐字：`Var_μ(⟨v,∇V⟩|0) ≤ E_μ[H_vv|0] − a²Δ_0²` … `with $W_{\mathrm{full}}=W_{\mathrm{geo}}+X$ and $\delta:=a^{2}\Delta_0^{2}+W_{\mathrm{full}}$; the two displays are the same inequality.`
   ⇒ 把 `δ` 代回：`Var ≤ E[H] − a²Δ₀² = E[H] + W_full − δ`。**与 `ch18:40` 同一个不等式**。
3. `ch18:81–83` 逐字自我限定：`every $\delta$ in this chapter is a dimensionless curvature floor and none is a physical mass bound`。

**⇒ `ch18:43` 的 `δ` = `δ_curv = a²Δ₀²+W_full`，`[T]`，无歧义。**

---

## §3 `I2` 的 `δ`：`app:789` 逐字，且 `δ_I2` 这个记号**不存在于语料**

- 逐字：`\item[I2] the isolation statement $\Lambda_\star-m_1^2\ge\delta>0$.`（`app:789`）
- 全仓检索 `delta_\{\\mathrm\{I2\}\}|delta_I2|delta^\{I2\}` → **0 命中**（工具 `grep` OK）。**`δ_I2` 是 `N12` 造的区分符，不是语料记号**；语料里它就是裸 `δ`。
- `I2` 的 `δ` 属**质量²** 槽：由 `ch05:426–428` 的核 `f_1^2 ≥ (Λ_★M_0−M_1)/(Λ_★−m_1^2) = M_0(Λ_★−R̄)/(Λ_★−m_1^2)` 的齐次性独立核出 —— 使 `Λ_★M_0−M_1` 成立，`Λ_★` 必须与 `m_1^2`、`M_1/M_0` 同量纲 ⇒ `Λ_★−m_1^2` 为**质量²**，其下界 `δ` 亦然。[T]
- `ch05:423` 那一屏里同时出现 `δ_{m_1^2}`（Dirac）与 `m_1^2`（槽位）——**Dirac 与层距在同一行共用字母**，是误读 `Λ_★−m_1^2≥δ` 的最短路径。

**⇒ 即便只谈"层距"这一义，`ch05:151`（`δ≈0.5m_1`，质量¹）与 `app:789`（`Λ_★−m_1^2≥δ`，质量²）也不自洽**：`ch05:149` 给的 `Λ_★≈1.56\,m_{0^{++}}` ⇒ `Λ_★−m_1^2≈1.43\,m_1^2`，而 `:151` 的 `0.5m_1` 是**质量差** `(1.56−1)m_1`。同一 `I2` 的两种写法相差一次方。[T]（`1.56^2−1=1.4336`，算术）

---

## §4 量纲是判据 —— 两个 `δ` 不可能是同一个

| | `δ_curv` | `δ_I2` |
|---|---|---|
| 定义处 | `app:191`：`δ := a²Δ₀²+W_full` | `app:789`：`Λ_★−m_1^2 ≥ δ > 0` |
| 槽 | `G″` 的**地板**（曲率） | `ν` 两支撑点之间的**距离**（层距） |
| 量纲 | **无量纲**（`ch18:81–83` 逐字：`G''_a=a^2m_{\mathrm{phys}}^2` ⇒ `δ_curv` 是 `a²·质量²` 组合后的无量纲曲率地板；`ch18:630` 逐字声明未做 `a→0` 标度） | **质量²**（与 `m_1^2`、`Λ_★`、`M_1/M_0` 同槽，§3 齐次性） |
| 载体 | 条件化路径测度 `μ_a(·|ξ=0)` 上的 `G″` | 谱测度 `ν` 的 `supp` |
| 是否随 `a→0` 变 | 是（`a²` 前置，`ch18:196,630`） | 否（物理质量²，`a`-无关；`N12:57` 引 `A02:34`/`F05:85`） |
| 平移协变 | 不适用（曲率无平移槽） | 是（`N12:72` 引 `G07:52–53`） |

**量纲不同 ⇒ 不能是同一个 `δ`。** 这一步不依赖任何解释，只依赖 `ch18:81–83` 与 `ch05:426–428` 两处语料内部的自证。[T]

---

## §5 判：`n05:120` 的"⟺ I2"是误并 —— 逐字

**`N05-SPLIT.md` 逐字（`:119–122`）**：

```
119: 1. `Δ₀>0` 那个靶的算式恰好是 `G″ − W_full ≥ a²Δ₀²`（把 `δ := a²Δ₀²+W_full` 代进 `H1(δ)`）…
     ——`P03-EQUIV.md:253–258` 已判 `[T]`：`H1 ⟺ G″ ≥ a²Δ₀²+W_full ⟺ Δ_dyn ≥ a²Δ₀²`，`W_full` 两端**同一个对象、直接相消**，"零新估计，只是换名"。
120: 2. 故 `Δ₀>0` **就是** `H1(δ_phys)`，而 `H1`/`I2`/`(E-I2)`/`H*`/`Δ₀>0` 是**同一命题的五个名字**（`paper-v6/sec07-assembly.tex:131`；`sec07:70`）。
122:    ⇒ … 因此 `Δ₀>0` **不可能**比 `I2`/`H*` 更强——它是**同一句**，只换了记号与量纲（`G″_a = a²m_phys²`，`ch18:81–83`）。
```

**逐条核**：

1. **`:119` 成立**，但只覆盖 `δ_curv` 内部：`P03:253` 的 `⟺` 链全程在 `δ_curv` 一个符号上（`W_full` 两端同一对象相消），`[T]`。**它没有跨到 `ν`/`I2`。**
2. **`:120` 第一句成立**：`Δ₀>0` 就是 `H1(δ_phys)`，此处 `δ_phys = δ_curv`。
3. **`:120` 第二句（`I2` 挂进"五个名字"）无推导**：句中唯一的支撑是 `paper-v6/sec07:131` + `:70` 的引文 —— **而这两处在现行语料里正是撤回该识别的地方**：
   - `sec07:70` 逐字：`An earlier draft named this gap $\Delta_0>0$ and identified it with $H1$, $I2$, $(E\text{-}I2)$ and $H^*$. **That identification is withdrawn**, on three independent grounds: (i) \code{n12} … (ii) \code{n06}: $\Delta_0$ and $\Lambda_\star$ live on different carriers ($G''$ vs $\operatorname{supp}\nu$) … **The correct equivalence class is three, not five**`
   - `sec07:132` 逐字：`$\delta_{\mathrm{curv}}:=a^2\Delta_0^2+W_{\mathrm{full}}$（\code{m05}）是**曲率槽**的量，与 $I2$ 的 $\Lambda_\star-m_1^2$ **同名不同物**`
   ⇒ `N05:120` 的引文**与它所引的现行正文相反**：引文被用来证"五个名字"，正文逐字写"撤回"+"三不是五"+"同名不同物"。[T]
4. **`:122` 自证其反**：逐字写"只换了记号**与量纲**"。按量纲判据，换了量纲就**不是**同一个量；这句话本身就是误并的自供。且 `ch18:81–83`（被 `:122` 引来支撑）恰恰是**否认**曲率地板是质量界的同一段。
5. **独立旁证（两条，均 `[T]`）**：`N06:50` `载体 ｜ 条件化路径测度 μ_a(·|ξ=0) ｜ 算符 ¼trF² 的谱测度 ν（质量²轴） ｜ 否`；`N06:145` `假：+ Δ₀>0 … 换了载体（G″ vs supp ν）、换了量词域；M05 [T] 判不可比`。`N06:148` 逐字把 `g03` 的并入判为 `误并`。

**⇒ 判**：`n05` 的 `⟺ I2` **无效** —— 它不是被算出来的，是被**同一个字母**接上去的。`n05` 的曲率侧结论（恒等式、`G″>0 ⇏ Δ_dyn>0`、`W_full ≥ N/12` 为 `[M]`）**不受影响**，仍成立；死掉的只有跨槽那一步。

---

## §6 判死线 (a)/(b)/(c)

- **(a) 不触发**：两个 `δ` 不是同一个（§2–§4）。语料无任何一处把 `a²Δ₀²+W_full` 与 `Λ_★−m_1^2` 写成同一个量；相反，`ch18:81–83` 与 `sec07:132` 逐字否认。
- **(b) 触发**：量纲不同（无量纲曲率地板 vs 质量²层距），载体不同（`μ_a` vs `ν`），`a`-依赖不同（`a²` 前置 vs `a`-无关）⇒ `n05` 的"⟺ I2"无效。
- **(c) 成立，且是本次记账的**主判**：语料本身混用同一符号。全部出现处见 §1；混用清单（须在论文里改名）：
  - `app:191` `δ_curv` ↔ `app:789` `δ_I2` ↔ `app:787` `δ_{m_1^2}`（Dirac）—— **同一文件、598 行内三义，符号表 `δ` 条目 0**。
  - `ch05:423` 同屏 `δ_{m_1^2}`（Dirac）与层距槽 `m_1^2`；`ch05:436` `e^{−δt}`（率）；`ch05:454` `δ_W,δ_S`（常数）；`ch05:456` `δΛ`（无量纲）；`ch05:151` `δ≈0.5m_1`（质量¹）vs `app:789`（质量²）。
  - 建议改名（`N12:74` 已提案一半）：`δ_curv` / `δ_I2`（或 `d_I2 := Λ_★−m_1^2`，直接避开 `δ`）/ `δ_{m_1^2} → \mathbf 1_{m_1^2}` 或 `\mu_{\mathrm{at}}` / `δt` 的率改 `γt` / `δ_W,δ_S → c_W,c_S`。

**保留的残余不一致（不属本问，但须记）**：`sec07:131` 逐字仍称"主力陈述 `G'_{k=0}/a²\ge\delta²>0` 的 `(E\text{-}I2)` 读法**属同一等价类**"。这是**换了一条路**（取比值使量纲对上）把曲率槽并进 `(E-I2)`；`N06:145` 已判该并入 `[NG]`（换载体）。改 `δ` 的名字**不能**修好这一句 —— 它是命题级误并，不是记号级。`[C]`

---

## §7 限定与不做的

- 本文件**不**证任何新定理，**不**改任何数值，**不**判 `H1`/`I2` 的真假，**不** QED。
- 结论只在**符号层 + 量纲层**：`δ_curv ≠ δ_I2`，`n05:120` 跨槽那一步无效。**曲率侧的 `[T]`（`P03:253`）原样保留**，`I2` 侧的 `[T]`（`app:789` 的声明形状）也原样保留；死的是二者的**桥**。
- `δ_curv` 的 `[T]` **不是**"`δ_curv>0` 已证"：它由 `Δ_dyn` 与 `W_full` 拼出，`W_full ≥ N/12` 仍是 `[M]`（`N05:88`、`N05` §3.1）。
- 未做：`a→0` 标度（`ch18:630` 同状）。

---

## §8 检索执行核验（AGENTS.md §1）

| 项 | 记录 |
|---|---|
| 工具状态 | `grep` **OK ×6**（全部返回带行号正文，无 `ripgrep provider failure`）；`glob` **OK ×2**（`docs/paper-v2/*.tex`→29；`**/{N05,N12,N06,P03}*.md`→4）；`read` **OK ×9**（`ch18:30–89`、`app:180–199/780–797`、`ch05:145–158/415–464`、`N05` 全文 159 行、`N12` 全文 144 行、`P03:235–284`）；`write` **OK**。**无 `ENOENT`、无超时、无空响应。** |
| 检索式 1 | `delta` @ `ch18-scalar-h1.tex` → **41 命中**（行号见 §1.1） |
| 检索式 2 | `delta` @ `app.tex` → **3 命中**（`:191, :787, :789`） |
| 检索式 3 | `delta` @ `ch05-theorems.tex` → **8 命中**（`:74, :151, :222, :423, :436, :440, :454, :456`） |
| 检索式 4 | `delta_\{\\mathrm\{I2\}\}\|delta_I2\|delta^\{I2\}` @ 全仓 → **0 命中** ⇒ `δ_I2` 非语料记号（采信，工具 OK） |
| 检索式 5 | `delta` @ `ch19-h1-status.tex` → **3 命中**（`:7, :8` 均 `δ_curv`；`:638` 为 `β²V₀+δ ≤ inf W_full`，同槽） |
| 检索式 6 | `delta` @ `ch08-immune.tex` → **2 命中**（`:14, :19`，均 `δ_curv`，即 `G″≥δ` 的 `(M)` 原形） |
| 检索式 7 | `五个名字\|delta\|Delta\|H\^\*` @ `docs/paper-v6/sec07-assembly.tex` → **13 命中**（关键：`:70` 撤回句、`:131` 主力陈述、`:132` `δ_curv`↔`I2` "同名不同物"） |
| 未查到（采信） | **"`δ` 的全局定义/符号表条目"**：`app.tex:780–797` 逐条通读，**0 条 `δ` 条目**（工具 `read` OK） |

---

## §9 出处（逐条可核）

**论文（行号已 `read` 核过）**：`docs/paper-v2/app.tex:185–191, 194, 780–797`；`ch18-scalar-h1.tex:10–15, 30–48, 55–58, 73–83, 114–122, 140, 183–196, 216–224, 239–253, 271–281, 308–342, 367–368, 605–630`；`ch05-theorems.tex:74, 145–158, 222, 415–463`；`ch08-immune.tex:14, 19`；`ch19-h1-status.tex:7–8, 638`；`docs/paper-v6/sec07-assembly.tex:70, 131–132`。

**species（行号已核）**：`N05-SPLIT.md:103, 119–122`；`N12-VERDICT.md:15, 56–61, 67–74, 105`；`N06-RELATIONS.md:13, 17, 37–41, 50, 80–85, 143–148`；`P03-EQUIV.md:253–258`。

---

*状态：`ch18:43` 的 `δ` = `δ_curv`（`[T]`，`ch18:40` + `app:191` + `ch18:81–83` 三重闭合）；`app:789` 的 `δ` = 质量²层距（`[T]`）。两者**非同一物**（量纲/载体/`a`-依赖三项全不同）⇒ 判死线 **(b) 成立、(c) 成立（同一 `δ` 六义，`app` 内三义并列且符号表无 `δ` 条目）**；`n05:120` 的"⟺ I2"是**误并**，其引文（`sec07:131, :70`）在现行语料里恰是撤回该识别处。**曲率侧 `[T]` 保留，桥作废。No QED；Clay 0%。***
