# V03-I2GAP — `I2`（`Λ★>m₁²`）与 `m>0` 的关系：**不可比**，不是「严格弱于」

> 2026-11-25 · 任务：判 `I2` 是否**严格弱于** `m>0`（`t02`「距离 `m>0` 恰好一步」）。
> 只读 + 只写本文件；未 git commit。**No QED；Clay 0%。**
> 标签：`[T]` 定理 · `[C]` 条件 · `[NG]` no-go（对某命名类）· `[M]` 模型级 · `[O]` 语料原标签「精确缺失步」。

---

## 0 · 判词（先给结论）

1. **`I2` 与 `m>0` 不可比（正交）**：`I2 ⇏ m>0`（`[T]`，T6-5(i) 见证 `m₁=0, I=T, G=F`），
   `m>0 ⇏ I2`（`[T]`，T6-4 见证 `m₁=2, 原子=T, 隙=T, I=F`）。
2. **所以「严格弱于」为假**——它要求 `m>0 ⇒ I2`，而这一半就被 T6-4 打死。
   **等价也为假**，**必要/充分皆假** ⇒ 判死线 **(a)(b)(c) 三支的前提全不成立**。
3. **你的 §3 论证方向要翻过来**：`I2` 只管道（`ν`= `tr F²` 的 `0⁺⁺` 测度），`m>0` 要全道 ⇒
   只管道**不足以**推出全道 ⇒ 那是「**范围更窄、命题不可比**」，不是「严格弱」。
   窄范围 + 不可比 = 连「A 比 B 弱」都不能说，只能说 **A 换了轴**。
4. **若硬按「还差什么」记账**，`I2 → m>0` 的缺件是 `m₁²>0` **与**「其它道底部 `>0`（或 `IN-11`）」，
   而 `m₁²>0` 正是目标本身的 `0⁺⁺` 投影 ⇒ **循环**。故 `t02` 的「一步」不是步数少算，
   是**方向错**：`I2` 不在 `m>0` 的下游，它在一个正交轴上更深一级（`I2 ⇒ I1`，`R2:118`）。

---

## 1 · `I2` 的精确陈述（逐字）

### 1.1 `docs/paper-v2/ch05-theorems.tex:421–437`（`I2` 的定义现场）

```
421: \begin{claimbox}{cT}
422: \claimhead{Lemma}{(residue--level-gap)} \emph{Suppose the lightest state is an isolated atom at
423: $m_1^2$, and let $\Lambda_{\star}:=\inf\operatorname{supp}(\nu-f_1^2\delta_{m_1^2})>m_1^2$.
424: Then}
425: \begin{equation}
426: f_1^2\ \ge\ \frac{\Lambda_{\star} M_0-M_1}{\Lambda_{\star}-m_1^2}
427: = M_0\,\frac{\Lambda_{\star}-\bar R}{\Lambda_{\star}-m_1^2},
428: \qquad \bar R:=\frac{M_1}{M_0}.
429: \label{eq:reslevel}
430: \end{equation}
431: \emph{The right side is positive \textbf{iff} $\Lambda_{\star}>\bar R$, and $\Lambda_{\star}>\bar
432: R$ suffices (if $f_1^2=0$ then $\Lambda_{\star}=m_1^2$ and $\bar R\ge\Lambda_{\star}$). The
433: converse fails (a positive atom with a heavy tail can have $\bar R>\Lambda_{\star}$). So this is
434: a self-certifying sufficient certificate, not an equivalence. By this lemma and Theorem A of
435: \code{mass-spectrum}, $\text{atom(I1)}\wedge\text{I2}\iff$ sharp two-exponential decay
436: $C(t)=f_1^2e^{-m_1t}(1+\eta(t))$, $|\eta(t)|\le Be^{-\delta t}$.}
437: \end{claimbox}
```

### 1.2 `ch05-theorems.tex:69–77`（`I2` 不是 P3 的推论——自我勘误）

```
69: \textbf{What it is not --- the correction that must be preserved.} Earlier drafts claimed the
70: integer cost grading ``additionally yields I2 (isolation).'' \textbf{That claim is wrong} and is
71: corrected in \code{PROOF-P3-UNCONDITIONAL.md} \S5 and \code{PROOF-P1P2.md} \S3. The cost grading
72: grades \textbf{polymer amplitudes}, not \textbf{transfer-matrix eigenvalue spacing}. It makes the
73: single tube the unique leading term (I1), but the KP/OS remainder is a $t$-uniform $O(u)$, never
74: an $O(e^{-\delta t})$; so it cannot separate a second atom from one at $\varepsilon\to0$, nor
75: place the excited $0^{++*}$ (at $\approx1.5\,m_1$) below the $2m_1$ continuum. \textbf{I2 remains
76: exactly the open step.} Hence P3 is a \textbf{residue theorem, not a gap theorem.} The
77: strong-coupling gap itself is classical Osterwalder--Seiler (1978) and is not what P3 proves.
```

### 1.3 `app.tex:787–789`（术语表，"绝对 `δ` 形式"）

```
787: \item[$\Lambda_\star$] the isolation threshold $\inf\operatorname{supp}(\nu-f_1^2\delta_{m_1^2})$.
788: \item[I1] the atom statement $f_1^2>0$.
789: \item[I2] the isolation statement $\Lambda_\star-m_1^2\ge\delta>0$.
```

⇒ **`I2` 逐字两种形式**：谓词式 `Λ★ := inf supp(ν − f₁²δ_{m₁²}) > m₁²`（`ch05:423`）；
定量式 `Λ★ − m₁² ≥ δ > 0`（`app:789`）。**`ν` = `tr F²` 的 Källén–Lehmann 测度，
`m₁` = `0⁺⁺` 道最轻质量**（`ch01:76` 逐字："where $m_1$ is the lightest mass in the $0^{++}$
channel"；`ch01:79` 的原子 `f₁²=ν({m₁²})=|⟨0|tr F²|𝔾⟩|²`）。**`I2` 的全部内容活在一个道里。**

---

## 2 · `m>0` 的精确陈述（逐字，取 `ch05` 的条件定理 + `ch01` 的读法）

```
451: \textbf{The conditional theorem.} Under the twelve named inputs \code{IN-0} \ldots{} \code{IN-11}
452: (\code{docs/conditional-theorem.md}), there exists a gauge-invariant measure $\mu$ on
453: $\mathbb{R}^4$ satisfying OS0--OS4, non-Gaussian with $g_{\mathrm{ren}}/\Lambda>0$, with
454: $\Delta/\Lambda\ge\min(\delta_W,\delta_S)>0$, unique vacuum, area law, and group-uniform
455: constants; and the gap is equivalent to
456: $\inf\operatorname{spec}(H_a|_{\Omega^\perp})\ge\delta\Lambda$ uniformly in $a$.
```

```
ch01:27: exponentially at a rate $m>0$ in every channel carrying the lightest state.
ch01:16: exists on $\mathbb{R}^4$ and has a mass gap $\Delta>0$.
```

⇒ **`m>0` 逐字 = `inf spec(H_a|_{Ω^⊥}) ≥ δΛ`（`a` 一致）+ `Δ/Λ ≥ min(δ_W,δ_S) > 0`**
（`ch05:454–456`）= **全 Hamiltonian 在真空正交补上的底 `>0`，跨所有道**；
`ch01:27` 的读法是「每个载最轻态的道都以速率 `m>0` 指数聚簇」。

### 2.1 与 `m₁` 的关系（关键，且语料自己写明了两者不同）

`ch05:271–272` 逐字（PF 定理的诚实话）：

```
271: of ``a positivity gauge-fixed PF lacked.'' But (i) $\lambda_1$ is the lightest state over
272: \emph{all} channels, not the $0^{++}$ specifically; (ii) $\inf_V\Delta_V>0$ is the mass gap and
```

⇒ **语料明说：全道的最轻态 ≠ `0⁺⁺` 的 `m₁`**。二者只在 `IN-11` 下被认同：
`ch05:468–469` 逐字 `one physical input (\code{IN-11}, the lightest bulk state in $0^{++}/2^{++}$)`；
`ch13-machine.tex:49–50` 逐字 `These are not six problems; they are one --- the volume-uniform lower
bound on $m_{0^{++}}$ along the path, which is the mass gap itself.`
**注意 `ch05:468` 的 `0⁺⁺/2⁺⁺`：语料自己承认最轻态可能在 `2⁺⁺`——那正是 `I2` 管不到的道。**

---

## 3 · 判：`I2` vs `m>0`

记 `I := (Λ★>m₁²)`，`G := (m>0)`。四问：

| 问题 | 判 | 依据 |
|---|---|---|
| `I ⇒ G`（`I` 充分） | **✘** `[T]` | `ν=f₁²δ₀+1_{[3,∞)}ρ`：`Λ★=3>0=m₁²` 而底在 `0`（`T6-BINDING-ADVERSARIAL.md:163–165`，T6-5(i)，witness `m₁=0, I=T, G=F`）；且 `ν` 只约束 `0⁺⁺`（`ch01:76`），其它道全空 |
| `G ⇒ I`（`I` 必要） | **✘** `[T]` | `ν=f₁²δ_{m₁²}+1_{(m₁²,∞)}ρ`，`m₁²>0`：原子+真隙但 `Λ★=m₁²`（`T6:158–161`，T6-4，witness `m₁=2, A=T, G=T, I=F`） |
| `I ⟺ G`（等价） | **✘** | 上面两行 |
| `I` **严格弱于** `G` | **✘** | 「严格弱」= `G ⇒ I` 且 `I ⇏ G`；**前半就假**（T6-4） |

### 3.1 道范围（你的 §3）：方向要翻

- **范围**：`I2 ⊂ {0⁺⁺}`（`ν` 是 `tr F²` 的谱测度，`ch01:79`、`R2-ISOLATION.md:4–5`、`ch15:318`
  "the atom on the $0^{++}$ Källén--Lehmann measure"）；`G` 要**全道**（`ch05:271–272`、`:455–456`）。
- **推论**：`I` 只覆盖一个道 ⇒ **`I ⇏ G`**（缺其它道的底）。
  **这不是「弱」，这是「不可比」**：`G` 也不覆盖 `I` 的信息（隔离/第二支撑点），T6-4 已证。
- ⇒ **「A 只管一个道、B 要所有道」⇒ A 不蕴含 B；而"严格弱"要求 A 被 B 蕴含。故前提矛盾。**
  **严格弱的正确替代词是：`I` 在范围上更窄、在强度上不可比。**

### 3.2 为什么 `Λ★` 的下界反不过来（唯一 `[T]` 界是上界）

- `R2:95–101` Lemma C：`Λ★ ≤ Λ_max(f₁²,M₀,M₁)`——**只有上界**；`ch05:550–551` 汇总表末行明写
  `$\Lambda_{\star}>\bar R$ & equivalence; $\Lambda$ lower bound` 排在 **"does not give"** 列。
- 正性/矩/Hankel 全部**不能**下界 `Λ★−m₁²`（`R2:130–151` W4，thin-gap `9→0.01`）`[NG]`。
⇒ `m>0` 无法机械地「免费」给出 `I2`；这也是 T6-4 的结构原因，不只是巧合。

---

## 4 · 若坚持「`I2` 到 `m>0` 还差什么」：两项，且第一项即目标

唯一精确的桥（逐字可核）：

```
I2 ∧ (m₁²>0)              ⇒ 0⁺⁺ 道有隙（底在 m₁²>0）
I2 ∧ (m₁²>0) ∧ IN-11      ⇒ m>0（全道；IN-11 把 0⁺⁺/2⁺⁺ 认成全道最轻）
```

- **缺件 1：`m₁²>0`。** `I2` 允许底在 `0`（T6-5(i)）。而 `m₁²>0` **就是目标在 `0⁺⁺` 上的投影**。
- **缺件 2：其它道底部 `>0`**（`2⁺⁺`、`0⁻⁺`、`1⁻⁺`…）。语料只在 `IN-11`（`ch05:468–469`，**物理输入
  `[C]`**）下把这一项消掉；**没有**任何定理说「`0⁺⁺` 隙 ⇒ 全道隙」。
- **缺件 3（若要求 `a→0`/`β≈2.5` 一致）：一致化**，即 `R2:276–282` 的 `(E-I2)` `ρ≥1+c`，`[O]`。
  `β<1/48` 的 KP/OS 余项是 `t`-一致 `O(u)`、**never** `O(e^{−δt})`（`ch05:73–74`）⇒ 强耦合不给。

⇒ **两项都未被语料覆盖为定理**；其中缺件 1 是**循环**（`U03-I2.md:234–236` 已记同一循环）。

---

## 5 · 判死线逐条落判

- **(a) `I2` 严格弱于 `m>0` 且缺件已被覆盖 ⇒「一步」成立**
  **✘ 不成立。** 前提假：`I2` 不蕴含 `m>0`（T6-5(i) `[T]`）。连带「缺件已被覆盖」也无从谈。
- **(b) `I2` 严格弱于 `m>0` 但缺件未覆盖 ⇒「一步」是低估，其实两步**
  **✘ 不成立**（前提同一处假）。**但它是三支中最接近的一支**，改写后才真：
  「`I2` 与 `m>0` **不可比**；缺件 = `m₁²>0` ∧ 其它道（或 `IN-11`），未覆盖；第一件即目标 ⇒ 循环」。
- **(c) `I2` 等价于 `m>0` ⇒「一步」是换名，`[NG]`**
  **✘ 不成立。** `m>0 ⇏ I2`（T6-4 `[T]`）。
  **唯一像换名的地方**是绝对 `δ` 形式 + SSC（`Λ★=ρ²m₁²`，`ρ` 有界）：`Λ★−m₁²≥δ ⟺ m₁²≥δ/(ρ²−1)`
  （`R2:248–250`、`:256–257`）——**但那是 `0⁺⁺` 单道内部的换名**，仍不含其它道；且 SSC 在
  Savvidy 双井上失效（`T6:166–169`）。故即使在这一形式下，`⟺` 也只是「`0⁺⁺` 隙 ⟺ `I2`」，**不是** `⟺ m>0`。

### 5.1 对 `t02` 那句话

`t02`（`LOOP-RUNBOOK.md:4565`、`:4681`）：「距离 `m>0` **恰好一步**：`Λ★>m₁²`（条件 I2）」。

- **作为「目标只剩 `I2` 一个前件」**：`[T]` 成立的方向只有 `I2 ⇒ I1`（`R2:118`，Lemma A）——
  即 `I2` 是**原子的父**，不是隙的儿子。
- **作为「离 `m>0` 一步」**：**✘**。不是少算步数，是**方向错**：`I2` 与 `m>0` 正交，
  二者共享的只有 `m₁²>0`，而那正是要证的东西。
- **正确说法**：`I2` 是「同一道内更高一级」的对象（`I2 ⇒ I1`），它与 `m>0` 的交集为空：
  单道 vs 全道，隔离 vs 正性。**「一步」应改写成「正交轴上更深一级」**。

---

## 6 · 引用的核验记录（AGENTS.md §1）

| 项 | 工具状态 | 实际检索式 | 命中 |
|---|---|---|---|
| `I2` 定义现场 | `read` OK（`ch05:410–469`、`:50–94`、`:120–149`、`:258–312`） | 行段直读 | 命中 |
| `I2`/`Λ★`/`m₁` in `ch05` | `grep` OK | `I2\|I1\|Lambda_\{\\star\}\|Lambda_\\star\|m_1` @ `ch05-theorems.tex` | 28 |
| `I2` 全 `paper-v2` | `grep` OK | 同上 @ `docs/paper-v2/*.tex` | 51 |
| `m>0` 逐字 | `grep` OK | `m>0\|m\\gt 0\|mass gap.*:=\|m:=\\inf\|\\Delta:=\\inf` @ `paper-v2` | 3（`ch01:27`、`ch01:120`、`ch18:462`） |
| `m>0` 的算符形式 | `read` OK | `ch05:451–462` | 命中（`:455–456`） |
| 全道 vs 单道 | `grep` OK | `2\^\{\+\+\}\|0\^\{\+\+\}\|other channel\|all channels` @ `paper-v2` | 104；关键 `ch05:271–272`、`:468` |
| 术语表 | `read` OK | `app.tex:783–796` | 命中（`:787–789`） |
| `ch01` 的 `ν`/`m₁` 定义 | `read` OK | `ch01-intro.tex:74–85` | 命中（`:76`、`:79`） |
| 反例 T6-4/T6-5 | `read` OK | `research/v3/T6-BINDING-ADVERSARIAL.md:140–174`（注：**不在** `species/` 下，首读 `species/T6-…` 报 not found，换路径后 OK） | 命中 |
| SSC / `(E-I2)` | `read` OK | `research/v3/R2-ISOLATION.md:238–282` | 命中 |
| `t02` 原句 | `grep` OK | `t02\|T02` @ `research/v3` | 37（含 `LOOP-RUNBOOK.md:4545/4589–4590/4678–4691/4835`） |
| 「其它道/全谱隙」是否已有语料论断 | `grep` OK | `其它道\|全谱隙\|all-channel\|every channel` @ `research/v3/*.md` | 1（`ICA-TEST.md:230`，信息论义，**非**道义）⇒ **0 实质命中** |
| 同义英文式 | `grep` OK | `single channel\|one channel\|only the $0^{++}$\|other sectors\|其它道\|全道` @ `paper-v2/*.tex` | 0 |

**未做**：未核 `PROOF-P3-UNCONDITIONAL.md §5`、`PROOF-P1P2.md §3`（仅经 `ch05:71` 转引，`[T]`-转引）；
未核 `G2D-ISOLATION-I2.md`、`STATE-DENSITY-I2.md`、`SUB2-I2-STRONG.md` 正文（经 `R2` 转引）。
外部文献未直接读（`ch05` 内的 Osterwalder–Seiler 1978 仅作被引记录）。

---

*本文件只写 `research/v3/species/V03-I2GAP.md`；未 git commit；无 QED；Clay 0%。*
