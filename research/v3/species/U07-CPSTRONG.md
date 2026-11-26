# U07-CPSTRONG — 强耦合侧 `(CP)`/`f₁²>0`：`β<1/48` 上 **`f₁²>0` 已是 `[T]`**（判 (a)）；`(CP)` 只拿到 `I1` 一半，`I2`（隔离）仍缺；`Z_O` **不是**强耦合侧的一环

> 2026-11-25 · 任务：核准 `t02` 的 `[T]`（`β<1/48` 有 `r₀>0`）到底证了什么，以及 `r₀>0 → f₁²>0` 是否已闭合。
> **判词**：**(a) 成立** —— `f₁²>0` 在 `β<1/48`（固定 `a`、柱面、`L`-一致、有限 `L`）**已经 `[T]`**，且这一环**既不需要 `(CP)` 也不需要 `Z_O`**：`r₀>0` 与 `f₁²>0` 不是两件事，而是**同一条不等式的两端**。
> 但若靶子是 `(CP)` 本身，则**只拿到 `I1` 半**；`I2`（`Λ★>m₁²`）在 `β<1/48` 上同样未证（`IMP-ARROWS.md:128` `d₂ ✗`）。⇒ **强耦合侧新靶 = `I2`，不是 `f₁²>0`。**

---

## §1 检索/核验记录（`AGENTS.md §1` 四要素）

| # | 工具状态 | 检索式 | 命中数 |
|---|---|---|---|
| 1 | `glob` OK | `**/PROOF-P3-UNCONDITIONAL.md` | 1（`research/unconditional/PROOF-P3-UNCONDITIONAL.md`） |
| 2 | `grep` OK（正常返回，非 ENOENT / 非 provider failure） | `f₁²\(2\.5\)\|Z_O\|1/48` in `research/v3/species/T02*.md` | 14 |
| 3 | `grep` OK | `1/48\|48` in `research/unconditional/P3*.md` | 7（`P3-AUDIT-CORRECTED.md` 4、`P3-FULLREP-KP.md` 3） |
| 4 | `grep` OK | `1/48\|beta_1\|β₁` in `docs/paper-v2/ch05-theorems.tex` | 11 |
| 5 | `grep` OK | `1/48\|I2\|isolat\|d₂` in `research/v3/species/IMP-ARROWS.md` | 15 |
| 6 | `grep` OK | `\(CP\)\|1/48\|I2` in `research/v3/species/DF-CONTINUITY.md` | 10 |
| 7 | `read` OK | `PROOF-P3-UNCONDITIONAL.md` 全 36 行；`ch05-theorems.tex:20–89`；`ch06-tools.tex:28–53`；`T02-F1SQ.md` 全 176 行；`IMP-ARROWS.md:100–144`；`U02-CP.md:1–40` | — |

**无 TOOL-FAIL。** 未逐字核：`P3-FULLREP-KP.md` / `P3-FOURWALL-DERIVATION.md` / `P3-DRESSED-MASS.md` / `P3-CYLINDER-RESIDUE.md` 的**证明主体**（只读了 `P3-FULLREP-KP.md` 的 `48` 命中行 120/124/136），按 `PROOF-P3-UNCONDITIONAL.md` §1–§5 的汇编采信。

---

## §2 `P3` 的精确陈述（逐字）

`research/unconditional/PROOF-P3-UNCONDITIONAL.md:10–16`：

> ## Theorem (corrected, honest scope)
> Fix `a>0`, `N`, `β₁=min(1/48, c/N)`. For `β<β₁`, on the **cylinder** `𝒞_L=T³_L×ℝ` (spatial torus, infinite time), for every spatial volume `L`, the lightest `0⁺⁺` atom of `O=¼tr F²` obeys
> ```
>   f₁²_L := ν_L({E_{1,L}}) = lim_{t→∞} e^{E_{1,L}t} C_L(t)  ≥  r₀(β) > 0,
> ```
> **uniformly in `L`** (limit `t→∞` at fixed `L` first; `L→∞` is not taken — `L` is only the uniformity axis). This is the **residue (I1)**.
> **Label: `[T]` modulo the standard OS/Källén–Lehmann transfer-matrix representation on the cylinder.** **NOT obtained:** isolation `I2`, continuum `a→0`, transport to `β≈2.5`.

同一文件 `:31`：

> `f₁²_L ≥ r₀(β)=½c₁u^p>0`, uniform in `L`, on `𝒞_L`, for `β<β₁` — **`[T]` modulo the standard OS/KL representation** (§1).

`docs/paper-v2/ch05-theorems.tex:25–35`（Theorem P3）：

> `Fix $a>0$. There exist $\beta_1=\min(1/48,c/N)>0$ and, for each $\beta<\beta_1$, a constant $r_0(\beta)=\tfrac12 c_1(\beta)u(\beta)^p>0$ ($c_1$ a positive $\mathrm{SU}(N)$ group constant, $u(\beta)$ the fundamental character coefficient, $u\asymp\beta$; the exponent $p$ is fixed by the character normalization --- the corpus uses $p\ge1$ --- and its value only shifts the prefactor), such that, on the cylinder $\mathcal C_L=T^3_L\times\mathbb R$ and for every spatial volume $L$,`
> `f_{1,L}^2:=\nu_L(\{E_{1,L}\})=\lim_{t\to\infty}e^{E_{1,L}t}\,C_L(t) \ \ge\ r_0(\beta)\ >\ 0,`
> `(limit $t\to\infty$ at fixed $L$ first, then $L$-uniformity; $L\to\infty$ is not taken).`

### §2.1 `r₀` 是什么

`r₀(β)=½c₁(β)u(β)^p>0`：一个 **`β`（因而 `N`,`a`）依赖的正下界**，`c₁>0` 是 `SU(N)` 群常数，
`u(β)` 是基本特征标系数（`u≍β`），`p` 是算符内容指数（`ch05:28` 写 "the corpus uses $p\ge1$"；`PROOF-P3` 汇编在 `:31` 只写 `u^p`，`P3-AUDIT-CORRECTED.md:87` 写 `p=0` for the normalized A1g plaquette）。
**`p` 的取值只移动前因子，不改 `r₀>0`**（`ch05:28–29`）——这是三处文件的唯一不齐（**记账**，非缺陷）。

`r₀` **不是**质量、**不是**隙、**不是**隔离率：它是 `O=¼tr F²` 的 `0⁺⁺` 道**最低原子重量的下界**。
机制（`ch05:54–62`）：连接时刻 `0` 与 `t` 的极小聚合物是**四壁管**，`|P|≥4t/a`，`m₁=4ln(1/u)/a`，
振幅 `a₀(β)>0`（正基本特征标系数 × 正单态重耦合），其余聚合物各多付一个 `u`，KP 判据与骨架长度无关 ⇒ 余项 `O(u)` 对 `t,V,L` 一致。

### §2.2 `β<1/48` 的条件来自哪

三层，**必须分开**：

1. **真前提**是 **KP/Ursell 收敛条件** `u(β)K_KP<1`（`ch05:45–50`，逐字：*"The hypothesis is the per-site condition $u(\beta)K_{\mathrm{KP}}<1$ with explicit $K_{\mathrm{KP}}\le CN$ (the plaquette adjacency degree is $\le20$ …)"*），来自 `P3-FULLREP-KP.md:120`：`K ≤ C_0e^βN`，`C_0=4λ_a e`，`λ_a=O(20)`，*"in particular $K\le10^3N$ $(\beta\le1/48)$"*。
2. **数值 `1/48` 不是从 KP 解出来的**：`P3-FULLREP-KP.md:124` 逐字：*"the KP radius is $\beta_1=O(1)$ — **N-independent**, consistent with Shen–Zhu–Zhu's $\beta<1/48$. If instead $\beta$ is Hooft-rescaled so that $u\asymp\beta$ … the radius reads $\beta_1=O(1/N)$. Either way $K\le CN$ and the theorem hold; only the stated numerical $\beta_1$ shifts."* ⇒ `β₁=min(1/48, c/N)` 里的 `1/48` 是**保守取值 / 借用 SZZ 的半径**，且 `c₁`–`C_g` 常数**未化成十进制**（`ch05:48–50`）。
3. **`1/48` 的语义是方法半径，不是谱墙**（`ch05:79–84`，逐字）：*"The analyticity radius $1/48$ is the Ursell series radius, not a wall in the spectrum … so $1/48$ is a limitation of the *method*, not of the object. But the method stops there, and no absolute strong-coupling method can cross its analyticity radius."*
   同源的 `1/48` 出处：Shen–Zhu–Zhu, *CMP* **400** (2023) 805, arXiv:2204.12737，`|β|<1/(16(d−1))`，`d=4` ⇒ `1/48`（`P3-AUDIT-CORRECTED.md:51`；**该文证的是 LSI/Poincaré/指数聚类，不是残数、也不含 KP 全表示估计**，`:51` 明写 *"it cannot be cited for the KP full-rep estimate"*）。

---

## §3 `r₀>0` 是 `[T]` 吗？——是，`[T]` **mod OS/KL**，一句话证明

**`[T]`（mod OS/KL 表示）**。一句：四壁管是该道的**唯一**主导聚合物、其振幅 `a₀(β)>0` 是两个正因子之积，且 KP 余项是 `t,V,L`-一致的 `O(u)`（与骨架长度无关），故
`lim_{t→∞}e^{m₁t}C_L(t) ≥ a₀(1−O(u)) ≥ ½c₁u^p = r₀(β)>0`，`L`-一致——**不**需要任何谱假设超出已进口的 `ν_O≥0`（`ch05:52–62`、`PROOF-P3:§2–§5`）。

严格等级（**不得升格**）：
- `r₀>0` 的**算术**部分（`c₁>0`, `u>0`）是 `[T]` 恒等，平凡；
- **有内容的**是 `f₁² ≥ r₀` 这条不等式，等级 `[T] mod OS/KL`（`PROOF-P3:§1` 逐字 *"The one imported line: RP ⇒ `H≥0` self-adjoint, `C(t)=∫e^{−Et}dμ_O`, `dμ_O≥0`"*）；
- 域：**固定 `a`**、**柱面 `𝒞_L`**（`T⁴_L` 上 ill-posed：`PROOF-P3:19` *"`T⁴_L` is ill-posed (periodic propagator bounded ⇒ limit does not exist) — hence the cylinder"*）、**不取 `L→∞`**。
- ⇒ **`IMP-ARROWS.md:69` 行 17 与 `ch05:531` 表里把平台写成 `T⁴_L` 是记法漂移（勘误级）**，`T02-F1SQ.md:67` 已记账。

---

## §4 从 `r₀>0` 到 `f₁²>0`：需要什么？——**什么也不需要**（逐条）

这条"链"其实**没有箭头**，两端是同一式：

| # | 需要的环节 | 在 `β<1/48` 上是否需要 | 依据 |
|---|---|---|---|
| 0 | `f₁² ≥ r₀(β)` | **就是 P3 的结论本身**（不是一步推理） | `PROOF-P3:13,31`；`ch05:32`（eq.`p3`） |
| 1 | `ν_O≥0`（KL/OS） | 需要，但**已在 P3 的进口**，非新增缺口 | `PROOF-P3:19` |
| 2 | `r₀>0`（即 `c₁>0`、`u>0`、`p` 有限） | 是 P3 的**参数**，非独立命题 | `ch05:26–29` |
| 3 | `(CP)` | **不需要**。`(CP)` 只在 `β=2.5` 的**搬运**处出现 | `ch06-tools.tex:45–48`（`[LC]` on `(CP)`）；`T02-F1SQ.md:44` |
| 4 | `Z_O` | **不需要**。`Z_O²` 是 `2.5` 侧的乘性因子 | `ch06-tools.tex:34–48` |
| 5 | `I2`（`Λ★>m₁²`） | **不需要** `f₁²>0`；**需要**才谈隙 | `ch05:69–77` |
| 6 | 连续极限 `a→0`、`L→∞` | 不需要 `f₁²>0`；需要才谈物理残数 | `PROOF-P3:16`；`ch05:79` |

**一句话**：`β<1/48` 上 `f₁²>0` 的**充分条件集合 = {OS/KL 表示 + KP 半径}**，两者都已在 `[T]` 里；
`(CP)` 与 `Z_O` 属于 **`2.5` 侧的搬运链**，不属强耦合侧。

---

## §5 判死线

### (a) `β<1/48` 上 `f₁²>0` **已 `[T]`** ⇒ 强耦合侧**有基础**（完整链，逐字）

```
[进口] RP(自伴,0≤T_L≤1) + KL/OS(ν_O≥0, C(t)=∫e^{-Et}dμ_O)          [T] mod OS/KL     PROOF-P3:19
   ↓
[几何] 四壁管 |P|≥4t/a, m₁=4ln(1/u)/a, 振幅 a₀(β)>0                 [T]              PROOF-P3:§2
   ↓
[KP]   全表示 K≤10³N (β≤1/48), 邻接度 Δ=20, 余项 O(u) 对 t,V,L 一致  [T]              P3-FULLREP-KP:120; PROOF-P3:§3
   ↓
[重求和] 骨架-装饰, m₁^dressed, inf_L A^dressed ≥ ¼(2N)^{-2}>0       [T] mod OS/KL     PROOF-P3:§4
   ↓
[结论]  f₁²_L = lim_{t→∞} e^{E_{1,L}t}C_L(t) ≥ r₀(β)=½c₁u^p > 0, L-一致   [T] mod OS/KL  PROOF-P3:13,31
   ↓   （无箭头：同一式）
[f₁²>0 在 β<1/48]        ⟸ 这已是 (CP) 的 I1 半                      [T]              本例 §4
```

**但**：`(CP) = [CP-a] ∧ [CP-b] = I1 ∧ I2`（`DF-CONTINUITY.md:39` 逐字：`(CP) = [CP-a] ∧ [CP-b] = I1 ∧ I2 , 且 [CP-b] ⇒ [CP-a]（单向）`；`U02-CP.md:40` 同）。
⇒ 强耦合侧拿到的是 **I1 半**；`I2` 在 `β<1/48` 上**同样未证**（`IMP-ARROWS.md:128`：`d₂ ✗ ★隔离★ Λ★ > m₁²（= I2 = Clay 上一级）`）。
⇒ **"强耦合侧的 `(CP)` 有基础"只能读到"`I1` 有基础"**；`(CP)` 整句在 `β<1/48` 上**仍是 `[C]`/`[NG]`，不是 `[T]`**。

### (b) 差一环？——`β<1/48` 的 `f₁²>0` **不缺环**；缺环只在 `(CP)` 与 `2.5` 侧

| 靶子 | 缺哪环 | 依据 |
|---|---|---|
| `f₁²>0` @ `β<1/48` | **无**（已 `[T]` mod OS/KL） | 本例 §4 |
| `(CP)` @ `β<1/48` | **`I2`（`Λ★>m₁²`，隔离）** | `IMP-ARROWS.md:128`；`ch05:69–77` |
| `f₁²(2.5) ≥ Z_O²r₀` | **`(CP)` + `Z_O` 的两圈以上 + `β_eff` 匹配 + `a→0`** | `ch06-tools.tex:45–48`；`T02-F1SQ.md:80–84` |

⇒ **新靶 = `I2`（隔离），不是 `f₁²>0`、不是 `Z_O`。** `Z_O` 是 `[T]` 的**恒等**（`ch06:34–48`：*"It transports the sign and positivity of the residue; it transports neither the size nor the isolation."*），把它当"缺环"是**降级错误**：`Z_O²>0` 是平方，永远不会"造出"正性，也永远不会破坏正性。

### (c) `r₀` 与 `f₁²` 无关？——**否，`t02` 的链不错**

`r₀` **逐字就是 `f₁²` 的下界**（`PROOF-P3:13` 同一行），不是另一个量。`T02-F1SQ.md:43` 的 I-4 标注 `[T] mod OS/KL`、域 `β<β₁=min(1/48,c/N)`、固定 `a` —— **与原文一致**（本例 §2 逐字比对）。
**唯一的记账偏差**（不影响判词）：`T02-F1SQ.md:43` 的出处写 `PROOF-P3-UNCONDITIONAL.md:31`，`:31` 确实有 `r₀(β)=½c₁u^p>0`；但 `[T]` 的**主陈述**在 `:10–16`（`Theorem`）——引用指向 `:31` 无误，只是**不是首现场**。

---

## §6 路径（下一步该打哪）

1. **不要**再证 `β<1/48` 上的 `f₁²>0` —— 已 `[T]`，零新地盘（`IMP-ARROWS.md:135` 同判：*"链 II 在 `(a 固定, β<1/48)` 内闭合——代价是那儿的隙早被 OS/SZZ 独立证过 ⇒ 零新地盘"*）。
2. **新靶 = `I2`，且必须在 `β<1/48` 上先打**：这是 `(CP)` 的**另一半**，也是唯一能在**已进口最少**（无 `(CP)`、无 `Z_O`、无 `2.5`）的域里被攻的隔离命题。若在 `β<1/48` 上 `I2` 也立不住，则 `(CP)` 的强耦合锚点本身就是空的，`2.5` 侧搬运**无处可搬**。
3. **`I2` 的精确障碍已定位**（`ch05:74–76` 逐字）：*"the KP/OS remainder is a $t$-uniform $O(u)$, never an $O(e^{-\delta t})$; so it cannot separate a second atom from one at $\varepsilon\to0$, nor place the excited $0^{++*}$ (at $\approx1.5\,m_1$) below the $2m_1$ continuum."* ⇒ 需要**第二指数尺度**（`T02-F1SQ.md:167–168`），这是换物种问题，不是加常数问题。
4. **`T⁴_L` 记法漂移**（`ch05:531`、`IMP-ARROWS.md:69` 行 17）应回填为柱面 `𝒞_L`——**勘误级**，已在本例与 `T02-F1SQ.md:67` 双记账。

---

*本文件只写 `research/v3/species/U07-CPSTRONG.md`；未 git commit。*
