# U01-F1AUDIT — `f₁²` 定义与 t02 四条"已知蕴含"逐条核准

**判词（先给结论）**
- **① `Var_Ω(O) ≥ f₁²` 是 `f₁²` 的【上界】**（`f₁² ≤ Var_Ω(O)`）⇒ 不能用来证 `f₁²>0`。**【方向反，判死线 (a) 触发】**
- **② `C(t) ≥ f₁²e^{−E₁t}` 对 `C` 是下界，对 `f₁²` 仍是【上界】**（等价 `f₁² ≤ e^{E₁t}C(t) ∀t>0`，且该族随 `t→∞` 单调降到 `f₁²`）⇒ 同样不能证 `f₁²>0`。**【(a) 触发】**
- **④ `r₀(β)=½c₁(β)u(β)^p>0` 是【真下界】，`[T]` mod OS/KL，窗口 `β<β₁=min(1/48,c/N)`、固定 `a`、柱面、有限 `L`。** 从 `r₀` 到 `f₁²>0` 的链**完整闭合**（§5）。
- **③ `+I2 ⟺ 尖锐双指数` 的 `⟺` 是【定义等价】`[T]`，不是新定理；`I2` 本身 `[NG]`；`δ≈0.5m₁` 依赖 `[C]`+`[E]`。**【(c) 触发】**
- **t02 的"距离一步"必须改口径**：缺的不是"到 `f₁²>0`"（`β<1/48` 已证），而是 `d₂`（隔离）与 `d₅`（搬运）——同一堵墙的两个坐标，**两步**。逐字改写见 §6。

---

## 1. `f₁²` 的逐字定义 + 三个对象（`C`、`𝒢`、`ν_O`）

### 1.1 `ch05-theorems.tex:32`（eq.`p3`，claimbox `cT` = `\LT` = `[T]`，`main.tex:38,128`）

> `f_{1,L}^2:=\nu_L(\{E_{1,L}\})=\lim_{t\to\infty}e^{E_{1,L}t}\,C_L(t)\ \ge\ r_0(\beta)\ >\ 0,`

同行上文的窗口（`:29–35` 逐字）：`on the cylinder $\mathcal C_L=T^3_L\times\mathbb R$ and for every spatial volume $L$`，`(limit $t\to\infty$ at fixed $L$ first, then $L$-uniformity; $L\to\infty$ is not taken)`。

### 1.2 `ch05-theorems.tex:60`（逐字）

> `$f_1^2=|\langle\Omega|O|\mathbb{G}\rangle|^2$ is a modulus-square, so no sign can undo it.`

`:52` 给出 `O`：`$\operatorname{tr}F^2$ is a finite A1g sum of local plaquette observables`；`O = ¼trF²`（`T02-F1SQ.md:17` 表述）；`𝔾` = 最轻 `0^{++}` 胶球态。

### 1.3 `ch06-tools.tex:65–72`（eq.`fredholm`，逐字）

> ```
> f_1^2=-\lambda_1\operatorname*{Res}_{z=z_1}\mathcal{G}(z)
> =-\frac{1}{z_1}\frac{N_O(z_1)}{\mathcal{D}_O'(z_1)},
> ```
> `with $\mathcal{D}_O(z)=\det(1-zT_O)$, $\log \mathcal{D}_O=-\sum_k z^k\operatorname{Tr}(T_O^k)/k$, $z_1=e^{am_1}$ the smallest-modulus zero (so $\lambda_1=e^{-am_1}=1/z_1$ is the leading transfer-matrix eigenvalue), and "periodic orbits" $\operatorname{Tr}T_O^k=Z_{S^1_k\times T^3}$ (closed-time flux tubes).`

**`𝒢(z)` 是什么：[M]——同段只给出它的残数，未定义它。** 语料内唯一描述在 `ch06-tools.tex:219–220`（逐字）：`$G$ is the boundary density (the resolvent/ Stieltjes transform of $\nu$ near the real axis)`；`ch09-weakstar.tex:23–26` 给出 `\Pi(z)=\int_{[0,\infty)}d\nu(s)/(s-z)`。
**唯一与 `λ₁` 因子自洽的读法**：`𝒢(z)=⟨Ω|O(1−zT_O)^{−1}O|Ω⟩`（对**转移矩阵**的预解式）。此时 `1−ze^{−E₁}=−e^{−E₁}(z−z₁)`，`Res_{z=z₁}𝒢=−f₁²/λ₁`，故 `−λ₁Res𝒢=f₁²` ✔。若把 `𝒢` 读成能量变量 Stieltjes 变换 `∫dν/(s−z)`，则 `−λ₁Res𝒢=f₁²/z₁≠f₁²`（仅 `z₁=1` 相符）✘ ⇒ **两种读法有一个不相容，eq.`fredholm` 定义缺失，必须补**（本项记账 `[M]`，非定理问题）。

### 1.4 `C(t)` 是什么关联函数（逐字两处，同一对象）

- `research/unconditional/PROOF-P3-UNCONDITIONAL.md:19`：`Then $C_L(t)=\sum_{n\ge1}f_{n,L}^2e^{-E_{n,L}t}$, $f^2=|\langle n|O|\Omega\rangle|^2\ge0$; $f_1^2_L=\nu_L(\{E_{1,L}\})=\lim_t e^{E_{1,L}t}C_L(t)$`；同段 `0=E_0<E_1<E_2\le\ldots` 离散谱。
- `research/mass-spectrum/REPORT.md:51,57`：`A:=O−ω(O)·1`，`W_O(t):=ω(A(t)A(0))=∫e^{−iEt}dν_O(E)`；`:65` 零空间动量 `C_O(t):=W_O^{p=0}(t)=∫₀^∞e^{−m|t|}·dρ_O(m²)/(2m)`。

⇒ **`C(t)` = `O=¼trF²` 在零空间动量、连通（真空减除后）的时间关联函数，`C(t)=∫e^{−Et}dν_O(E)`；不是胶子传播子、不是 Wilson 环。**

### 1.5 `ν_O` 是什么测度（逐字）

- `research/mass-spectrum/REPORT.md:54`：`ν_O(B):=⟨AΩ|P_H(B)AΩ⟩，B∈𝔅([0,∞))`，`[0,∞)` 上**有限正测度**，`ν_O([0,∞))=‖AΩ‖²`。
- `docs/paper-v2/ch02-squeeze.tex:10–12`：`Let $\nu\ge0$ be the Källén–Lehmann measure of $\operatorname{tr}F^2$, with $\operatorname{supp}\nu\subset[m_1^2,\infty)$, atom weight $f_1^2=\nu(\{m_1^2\})$, and moments $M_k=\int s^k\,d\nu(s)$.`
- `docs/paper-v2/ch05-theorems.tex:423`：`$\Lambda_{\star}:=\inf\operatorname{supp}(\nu-f_1^2\delta_{m_1^2})>m_1^2$`；`:440`：`$\nu=f_1^2\delta_{m_1^2}+\mu$`，`$\mu\ge0$`。

⇒ **`ν_O` = `O` 的 Källén–Lehmann 谱测度（= 谱定理作用于向量 `AΩ` 的谱测度），正、有限、`supp⊂[m₁²,∞)`，`f₁²=ν_O({m₁²})`。** 记号漂移（记账 `[M]`）：`PROOF-P3`/`ch05:32` 用能量 `E_{1,L}`，`ch02`/`ch05:423` 用质量平方 `m₁²`；`E₁=m₁`（`a·m₁` 出 `z₁`），非新对象。

---

## 2. ① `Var_Ω(O) ≥ f₁²`——**方向：上界**

**证明一句（一行恒等，`[T]`）**：`ν_O≥0` 且 `Var_Ω(O)=‖AΩ‖²=ν_O([0,∞))=M₀`（`REPORT.md:54`；`ch02:12` "M₀=∫dν is the total spectral weight"），故
`Var_Ω(O)=M₀ ≥ ν_O({E₁})=f₁²`。

**它给的是上界**：`Var_Ω(O) ≥ f₁²` ⟺ `f₁² ≤ Var_Ω(O)`。**`Var` 再大也不产生 `f₁²` 的正下界**——要证 `f₁²>0` 需要 `f₁² ≥ 正数`，与 ① 符号相反。**⇒ 不能用来证 `f₁²>0`，用户判断正确。**

**同文件紧邻处即是反向反例**（`ch05-theorems.tex:284–287` 逐字）：`The decisive check is on the physical space itself: free Maxwell is PF and strictly positive with a finite-volume gap, yet its $0^{++}$ spectrum is a two-photon cut from zero, $f_1^2=\nu(\{0\})=0$. Gauge-invariant PF plus strict positivity plus a finite-volume gap \textbf{does not imply} a glueball atom.` ⇒ 自由场 `Var>0`、`f₁²=0`，①**不可逆**。

**引文错配（记账 `[M]`）**：`T02-F1SQ.md:40` 把 I-1 的出处写作 `ch05-theorems.tex:284–287`；但那三行讲的是 **`Var>0 ⇏ f₁²>0` 的反例**，不是 `Var ≥ f₁²`。全仓检索（§7）确认：**①在 `ch05`/`ch06`/`PROOF-P3-UNCONDITIONAL` 中无逐字出处，`Var_Ω` 在论文侧亦未定义**——它是 `T02-F1SQ.md:40,142` 自造的记号。作为 `f₁²>0 ⇒ Var ≥ f₁²>0` 的前向蕴含它自洽（前件已含结论）；一旦被摘录成独立"已知蕴含"就变成上界工具，**必须按上界记账**。

---

## 3. ② `C(t) ≥ f₁²e^{−E₁t}`——**方向：对 C 下界，对 f₁² 仍是上界**

**逐字来源**：`PROOF-P3-UNCONDITIONAL.md:19`（§1 已引）`C_L(t)=Σ_{n≥1}f_{n,L}²e^{−E_{n,L}t}`，`f²≥0`。一行推论（`[T]`）：
`C_L(t) ≥ f₁²e^{−E₁t}`（`∀t>0`；`E_n≥E₁`，`f²≥0`）。

**方向核准**：上式 ⟺ `f₁² ≤ e^{E₁t}C_L(t)`。且
`e^{E₁t}C_L(t)=f₁²+\sum_{n\ge2}f_{n,L}²e^{−(E_{n,L}−E_1)t}` 对 `t` **单调降**、`≥f₁²`、`t→∞` 极限 `=f₁²`。
⇒ **② 恰是把 `ch05:32` 的等号定义 `f₁²=lim_t e^{E₁t}C_L(t)` 改写成不等式族 `f₁²=inf_{t>0}e^{E₁t}C_L(t)`，零新增信息；它对 `f₁²` 是上界族，不能给 `f₁²>0`。**【(a) 触发】

**附带记账**：`ch05:32` 原文是**等号定义**（`=lim`），t02 的 `≥` 版本是其弱化，不是独立定理；`T02-F1SQ.md:41` 的"不给上界"仅指**对 `C`** 不给上界——对 `f₁²` 恰恰给的是上界，措辞需加限定。

---

## 4. ③ `+I2` 的精确陈述 + "尖锐双指数"

### 4.1 `ch05-theorems.tex:422–436`（逐字，claimbox `cT`）

> `Suppose the lightest state is an isolated atom at $m_1^2$, and let $\Lambda_{\star}:=\inf\operatorname{supp}(\nu-f_1^2\delta_{m_1^2})>m_1^2$. Then`
> `f_1^2\ \ge\ \frac{\Lambda_{\star} M_0-M_1}{\Lambda_{\star}-m_1^2}= M_0\,\frac{\Lambda_{\star}-\bar R}{\Lambda_{\star}-m_1^2},\qquad \bar R:=\frac{M_1}{M_0}.`
> `The right side is positive \textbf{iff} $\Lambda_{\star}>\bar R$ ... So this is a self-certifying sufficient certificate, not an equivalence. By this lemma and Theorem A of mass-spectrum, $\text{atom(I1)}\wedge\text{I2}\iff$ sharp two-exponential decay $C(t)=f_1^2e^{-m_1t}(1+\eta(t))$, $|\eta(t)|\le Be^{-\delta t}$.`

`:443`（逐字）：`i.e.\ $f_1^2\le M_0-M_1/\Lambda_{\star}$ (an upper bound, not the claim)` —— 同一引理里同形状的上界被明确标为"不是结论"，可作 ①②的旁证。

### 4.2 `I2` 是什么

`I2` = **隔离（上限隙）**：`Λ★>m₁²`，即 `ν` 在 `(m₁², Λ★)` 无谱（`ch05:423`）。等价谱表述（`research/mass-spectrum/REPORT.md:77–81` 逐字）：`IsolatedOneParticle(O;m₁,m₂) :⇔ 0<m₁<m₂ ∧ ν_O({m₁})>0 ∧ ν_O((0,m₂)\{m₁})=0`，`:83`：`"上限隙"就是 $m_1$ 在 $\operatorname{supp}\nu_O$ 中孤立`。
**它不是 `f₁²>0`**：`ch05:69–77` 逐字 `The cost grading grades \textbf{polymer amplitudes}, not \textbf{transfer-matrix eigenvalue spacing} ... \textbf{I2 remains exactly the open step}. Hence P3 is a \textbf{residue theorem, not a gap theorem}.`；`PROOF-P3:16`：`\textbf{NOT obtained:} isolation I2`.

### 4.3 "尖锐双指数"是什么

`research/mass-spectrum/REPORT.md:88–90`（逐字定义）：`SharpDecay(O;m₁,Z,m₂,C) :⇔ 0<Z ∧ m₁<m₂ ∧ ∀t≥0, |W_O(−it)−Ze^{−m₁t}| ≤ Ce^{−m₂t}`。与 `ch05:436` 的 `C(t)=f₁²e^{−m₁t}(1+η(t))`、`|η(t)|≤Be^{−δt}` 一致，`δ=m₂−m₁`。

### 4.4 **(c) 触发：`⟺` 依赖什么**

1. **`⟺` 是定义等价，`[T]`（测度论），但前件自带 `m₂`。** `REPORT.md:93–98` 定理 A：`设 ν 是 [0,∞) 上有限正测度，0<m₁<m₂ 则 IsolatedAtom(ν;m₁,m₂) ⇔ ∃Z C, SharpDecay(ν;m₁,Z,m₂,C)`。其中 `IsolatedAtom`（`:77–81`）**已含** `ν_O((0,m₂)\{m₁})=0`，即已含 `I2`。故 `ch05:435` 的 `⟺` **不产生 I2**，只把"I2 成立"翻译成"双指数"——**把缺的东西留在前件**。`REPORT.md:110` 亦写明方向：`⇐ 需要正性`（`ν≥0`，本体满足）。
2. **`ch05:435` 的 "By this lemma and Theorem A" 是多余引用**（`[M]`）：定理 A 单独即给 `⟺`；残数–level-gap 引理（`Λ★>R̄ ⇒ f₁²>0`，`ch05:426–434`，`IMP-ARROWS.md:66` 行 14 标 `[T]`）不是 `⟺` 的一环，它的用途是**另一条**路（前件 `Λ★` 比 Clay 更强）。
3. **`δ≈0.5m₁` 依赖未证/条件物**：`ch05:149–151` 逐字 `the \emph{operative} threshold $\Lambda_{\star}$ for isolation is the excited $0^{++*}$ at $\approx1.5\,m_g$ (physically $1.56\,m_{0^{++}}$) ... the correct I2 is a gap to the excited glueball of size $\delta\approx0.5\,m_1$, fixed by power counting, not symmetry`。其中"`0^{++}` 严格低于两体阈"由 **no-dissolution** 定理给，而该定理 `ch05:127–135` 是 `[C]`（假设 `(K-A)` 一致隙、`(K-B) m_{0++}<2m_lightest`）；`1.56` 是数值 `[E]`。⇒ **`δ` 的数值不是 `[T]`。**

---

## 5. ④ `r₀` 的等级 + 从 `r₀` 到 `f₁²>0` 的完整链

### 5.1 逐字（`ch05-theorems.tex:25–36`，claimbox `cT`=`[T]`）

> `There exist $\beta_1=\min(1/48,c/N)>0$ and, for each $\beta<\beta_1$, a constant $r_0(\beta)=\tfrac12 c_1(\beta)u(\beta)^p>0$ ($c_1$ a positive $\mathrm{SU}(N)$ group constant, $u(\beta)$ the fundamental character coefficient, $u\asymp\beta$ ...), such that, on the cylinder $\mathcal C_L=T^3_L\times\mathbb R$ and for every spatial volume $L$,` [`eq.`p3`]

`PROOF-P3-UNCONDITIONAL.md:31`（逐字）：`f₁²_L ≥ r₀(β)=½c₁u^p>0, uniform in L, on 𝒞_L, for β<β₁ — [T] modulo the standard OS/KL representation (§1)`。表行 `ch05:531–532`：`P3 residue (I1) & \LT{} ($+$OS/KP) & fixed $a$, $\beta<1/48$, $T^4_L$ & $f_1^2\ge r_0(\beta)>0$, uniform in $V,L$ & I2 isolation; $a\to0$; $\beta\approx2.5$`。

**等级：`[T]`（mod 标准 OS/KL 表示；`[T]` 的定义见 `main.tex:128`：`\LT{} theorem — proved (any imported standard result is named)`）。不是 `[C]`。**

### 5.2 `f₁²>0` 的完整链（每步都有出处，全 `[T]` mod OS/KL）

1. **KP 收敛**：`β<1/48 ⌐ uK<1` ⇒ `Σ_{γ∩γ₀≠∅}|w(γ)|e^{a|γ|}≤a` 绝对收敛（`ch05:94–108` eq.`kp`；`PROOF-P3:25` `K≤6.1×10²N`，邻接 `Δ=20`）。
2. **四壁管唯一极小**：连接时 `0` 与 `t` 的极小聚合物 `|P|≥4t/a`（`PROOF-P3:22` 逐字 `The minimal connected, gauge-invariant, nonvanishing polymer linking the two single-plaquette insertions at times 0,t has |P| ≥ 4t/a`；`ch05:52–63,111–116`），等号 = 直四壁管，权 `a₀u^{4t/a}=a₀e^{−m₁t}`，`m₁=4ln(1/u)/a`，`a₀>0`（`ch05:57–59`：正基本特征标系数 × 正单态重耦合）。
3. **谱桥**：`f₁²_L=ν_L({E_{1,L}})=lim_{t→∞}e^{E_{1,L}t}C_L(t)`（`PROOF-P3:19`；`ch05:32`）。
4. **振幅地板**：`inf_L A^dressed ≥ ¼(2N)^{-2}>0`（`PROOF-P3:28(d)`）。
5. ⇒ `f₁²_L ≥ r₀(β)=½c₁u^p>0`，`L`-一致（`PROOF-P3:31`）；`f₁²=|⟨Ω|O|𝔾⟩|²` 是模方，符号不可被抹（`ch05:60`）。

**关键限定（`PROOF-P3:16,19` 逐字）**：`\textbf{Isolation I2 is NOT needed for the finite-$L$ floor} (compactness isolates every atom); it is needed only for the infinite-volume "atom vs continuum threshold" and $f_1^2_\infty$ — the open rung.`
⇒ **在 `β<1/48`、固定 `a`、柱面、有限 `L` 上，`f₁²>0` 已闭合（`[T]`），不需要 `I2`、不需要 ①②。** `I2` 只在 `L→∞`/一致化时才成为缺口。

### 5.3 窗口核验（AGENTS §2.1）

P3 的窗口（逐字，`ch05:29–35` + `PROOF-P3:11,15`）：固定 `a>0`、`N`、`β₁=min(1/48,c/N)`、`β<β₁`、柱面 `𝒞_L=T³_L×ℝ`、`t→∞` 先在固定 `L` 取、`L`-一致、**不取 `L→∞`**。
**勘误 `[M]`（记法漂移，两处）**：`ch05:531` 表行与 `IMP-ARROWS.md:69`（行 17）把窗口写成 `T⁴_L`，与 `ch05:30`（`𝒞_L`）、`ch05:35` 及 `PROOF-P3:19`（逐字 `\textbf{$T^4_L$ is ill-posed} (periodic propagator bounded ⇒ limit does not exist) — hence the cylinder`）**矛盾**。`T02-F1SQ.md:67` 已记同一漂移；本次复核确认是**论文侧也有的**错误，不只是 IMP-ARROWS 的。

---

## 6. 判死线回执：(a) 触发后的逐字重述

### (a) 结论

**①、② 方向都是 `f₁²` 的上界**（§2、§3）⇒ **`f₁²>0` 不能从它们推出**。原文 `T02-F1SQ.md:142` 的句子必须重写：

> **原文（作废）**：`不是"不衰减"，是**单侧地板** `C(t)≥f₁²e^{−E₁t}`（I-2）+ `Var(O)≥f₁²`（I-1）。**与 `m>0` 的距离 = 恰好一个对象：`Λ★>m₁²`（I2）**，即把"地板"升级成"率"所需的隔离`

> **改写（逐字采用）**：
> "**①②不是通往 `f₁²>0` 的路。** ① 给 `f₁² ≤ Var_Ω(O)`（`ν≥0` 一行恒等，`Var_Ω(O)=M₀=ν_O([0,∞))`），② 给 `f₁² ≤ e^{E₁t}C(t)`（`∀t>0`，该族随 `t→∞` 单调降到 `f₁²`）——**两者都是 `f₁²` 的上界族，是 `f₁²>0` 的推论，不是它的来源**；把它们当工具只会给出"原子不能太大"（反例方向：自由 Maxwell `Var>0`、`f₁²=0`，`ch05:284–287`），永远不会给出 `f₁²>0`。
> **`f₁²>0` 本身在 `β<1/48`、固定 `a`、柱面 `𝒞_L`、有限 `L` 上已经由 P3 的 `r₀(β)=½c₁u^p>0` 无条件给出（`[T]` mod OS/KL；`ch05:24–36`、`PROOF-P3:31`），且有限 `L` 不需要 `I2`（紧算符使每个原子孤立，`PROOF-P3:19`）。**
> **从 `f₁²>0` 到"隙"缺两步，不是一步**：`d₂` = 隔离 `Λ★>m₁²`（= `I2`，只在 `L→∞`/一致化与 `f₁²_∞` 处缺），`d₅` = 搬运（`a→0` 且 `[1/48,2.5]`，条件于 `(CP)`）。二者是同一堵墙的两个坐标（`IMP-ARROWS.md:126–135`，链 II：长度 5、已证 3、缺 2）。"

### (b) 部分成立

`④` 确实给**下界**；从 `r₀` 到 `f₁²>0` 的链见 §5.2（5 步全 `[T]` mod OS/KL）。但该链**只覆盖** `β<1/48`、固定 `a`、有限 `L`；它**不给** `L→∞` 的 `f₁²_∞`（需 `I2`）、不给 `a→0`、不给 `β≈2.5`。

### (c) 触发

`③` 的 `⟺`：`[T]` 但为**定义等价**（定理 A 的前件 `IsolatedAtom` 已含 `m₂`=`I2`，`REPORT.md:77–81,93–98`）；`I2` 本身 `[NG]`；`ch05:435` 对残数–level-gap 引理的引用是**多余**（`[M]`）；`δ≈0.5m₁` 依赖 `[C]`（no-dissolution，`(K-A)+(K-B)`，`ch05:127–135`）+ `[E]`（`1.56`，`ch05:149–151`）。

---

## 7. AGENTS §1 检索/工具状态

| # | 检索式 | 工具 | 状态 | 命中 |
|---|---|---|---|---|
| 1 | `Get-ChildItem D:\Dev\math -Filter "ch0*.tex" -Recurse` | pwsh | **OK**（有效输出，含 `ch05-theorems.tex` 36587B、`ch06-tools.tex` 15957B） | 10 |
| 2 | `grep "\s*pattern=C_L\(t\)\s*:=\|C\(t\)\s*:=\|\\nu_L\|two-exponential\|sharp two"` on `docs/paper-v2/*.tex` | grep | **OK** | 6 |
| 3 | `grep "C\(t\)\|C_L\|mathcal\{G\}\|spectral measure\|nu_L"` on `ch05-theorems.tex` | grep | **OK** | 6 |
| 4 | `grep "K[aä]ll[eé]n\|spectral measure\|resolvent\|T\)\^\{-1\}"` on `docs/paper-v2` | grep | **OK** | 15 |
| 5 | `Set-Content`/`Select-String "Var_\\Omega\|Var_Ω"`，`Get-ChildItem -Recurse -Include *.md,*.tex` 全仓 | pwsh | **OK** | 10（全在 `research/unconditional/{AQFT-CAUSAL,BLAST-TRANSFER,JONES-SUBFACTOR}.md`，**均非 `f₁²` 语境**；`Var_Ω(O)` 在论文侧 0 命中） |
| 6 | `Select-String "mass-spectrum"` 全仓 | pwsh | **OK** | 20（含 `ch05:435`；**`research/mass-spectrum/REPORT.md` 实存 41272B**） |
| 7 | `Select-String "Theorem A\|Thm A"` on `research/mass-spectrum/REPORT.md` | pwsh | **OK** | **0**（`Theorem` 为英文串；中文 `定理 A` 命中，见 #8）⇒ 非"不存在"，是字符串语言不同 |
| 8 | `Select-String "Theorem\|定理"` on `REPORT.md` | pwsh | **OK** | 25（定位 `§1.3 定理 A`、`§1.4 定理 B/B′`、`§1.6 定理 D`） |
| 9 | `read`：`ch05:1–110,111–180,275–304,420–459,515–544`；`ch06:1–40,40–109,205–232`；`ch02:1–80`；`ch09:1–30`；`PROOF-P3-UNCONDITIONAL.md`（全 36 行）；`T02-F1SQ.md`（全 176 行）；`IMP-ARROWS.md:50–139`；`REPORT.md:48–132`；`ch15:470–499` | read | **OK** | — |
| 10 | `grep "I2\|(I2)"` on `ch05-theorems.tex` | grep | **OK** | 6 |

**TOOL-FAIL：无。** 未重跑任何数值；Lean/`Handoff*.lean` 未打开；`refs.tex` 未核（本节不引入新文献，仅用内部行号）。

---

## 8. 一句话总账

`f₁²` = `ν_O({E₁})` = `lim_t e^{E₁t}C(t)` = `−λ₁Res_{z₁}𝒢(z)` = `|⟨Ω|O|𝔾⟩|²`（`ch05:32,60`；`ch06:65–72`；`𝒢` 定义缺 `[M]`）。
**①=上界、②=上界 ⇒ (a) 触发**；`f₁²>0` 的真正来源是**④ `r₀`（`[T]`，`β<1/48`，有限 `L`）**，缺口是 `d₂`(I2)+`d₅`(搬运) **两步**而非"恰好一步"；`③` 的 `⟺` 是定义等价、`I2` 仍 `[NG]` ⇒ **(c) 触发**。

*本文件只写 `research/v3/species/U01-F1AUDIT.md`；未 git commit。*
