# T01-IMPLY — `G″(0) > 0` 到底能推出什么：四条候选路径的逐条判死，与"第 4 级台阶"的精确形态

**日期**：2026-11-25 · **只写** 本文件 · **禁 git commit**（遵守）
**输入（逐字核过，行号在下）**：`docs/paper-v2/ch05-theorems.tex:230–240`（`eq:split`）、`ch15-uvface.tex:20–28`（`eq:v3object`）、
`ch18-scalar-h1.tex:32–48/50–70/77–83`、`ch19-h1-status.tex:32–47/56–58/82–93/125–149`、`ch08-immune.tex:11–21`（`(M)`）、
`ch04-gap.tex:23`（G2.0 恒等层）、`ch06-tools.tex:11–16`（酸测）、`research/v3/species/Q11-CERT.md`、`R06-GAPMONO.md:11–20`、
`R10-AUDITCERT.md`、`IMP-ARROWS.md:53–95`（行 5(iv)/行 8/§2.1(a′)）、`MU5-SIGMA-M.md:14–18`（Chatterjee 方向）。

**标签**：`[T]` 真定理 / `[C]` 条件定理 / `[NG]` 对具名类的 no-go / `[M]` 精确缺失步。

---

## §0 判词（先给，后逐句核）

| 判死线 | 判定 |
|---|---|
| **(a)** 存在已被证明的链？ | **存在，但只有一条，且它是改名链，内容长度 = 0**：`G″>0 ⟺ H1(0)`（`[T]`，`eq:scalar-rule`）`⟺ m_eff²(k→0)>0`（`[T]`，`PROOF-IR-H1` 自标 "rename"）。它证的是 `H1 ⟺ H1`。**这不是"今天最重要的发现"，是"今天最重要的记账"：`G″(0)>0` 不是通往 H1 的一级台阶，它就是 H1 的 `δ=0` 切片。** |
| **(b)** 链条断在哪？ | **断在"符号 → 带余量的 `(L,a)`-一致下界"这一步**。逐字：`G″ = Δ_dyn + W_full` 且 `W_full ≥ N/12 > 0` ⇒ **`G″>0` 只给 `Δ_dyn > −W_full`，不给 `Δ_dyn > 0`**；而物理 H1 要的是 `Δ_dyn ≥ a²Δ_0²`（`eq:h1st-dyn`）。差的那一项正是 `W_full + a²Δ_0²`。**这就是第 4 级台阶。** |
| **(c)** `G″>0` 与 `m>0` 无蕴含（有`G″>0`而`m=0`的模型）？ | **在恒等式适用域内不可能**（`G″` 与 `m_eff²(k→0)` 是同一个数，模型无关代数，举不出反例）。**但裸蕴含在自然类内为假**（`G_YM × U(1)_weak`，§4.3）⇒ H1 路线**不必重新定位，但必须把量词写进对象**：对象不是"`G″>0`"，而是"`G″` 的一致正下界 + 纯性"。 |

---

## §1 `G″(0)` 是什么：哪个量的二阶导

**逐字核准（不是记忆）。**

`ch15-uvface.tex:20–28`（`eq:v3object`，标签在 `:26`）：

> "The object in $(\star)$ is the *coarse-field-conditioned marginal curvature*
> $G''(v;a_0)=E_a[\partial_v^2V]-\operatorname{Var}_a(\partial_vV)=\Delta_{\mathrm{dyn}}+W_{\mathrm{geo}}+X$, $X=-2\operatorname{Cov}(G,\Phi')$,
> with $V=S_W+\Phi$, $\Phi=-\log J$ the Haar Jacobian, **conditioned at the coarse field $a=a_0$**."

`ch05-theorems.tex:230–240`（`eq:split`，标签在 `:239`）：

> "$G''(0)=E[H_{vv}+\Phi''\,|\,\xi=0]-\operatorname{Var}(G_v+\Phi'\,|\,\xi=0)=\Delta_{\mathrm{dyn}}+W_{\mathrm{geo}}+X$,
> $\Delta_{\mathrm{dyn}}=E[H]-\operatorname{Var}(G)$, $X=-2\operatorname{Cov}(G,\Phi')$."

**⇒ 它是哪个量的二阶导：粗场条件化边缘自由能（有效势）`V_eff(a)` 对粗幅 `a` 的二阶导，在 `a=0` 处取值、逐粗模。**
`V = S_W + Φ`，`∂_v := ⟨v,∇⟩`（未归一；`‖v‖=1` 时即上面两式，`‖v‖²=1/L` 时差 `L²`，见 `ch05:214–228`）。

**它有四个同义名（这是全部麻烦的来源）**，`ch08-immune.tex:11–21` 逐字：

```
(M): Var_μ(⟨v,∇S⟩|0) ≤ E_μ[H_vv|0] + W_full − δ ,  uniform in (L,U')
"the display is then G''≥δ in curvature form.
 Equivalently Var(ξ_{k→0}) ≤ 1/λ (H1), equivalently m_eff²(k→0) ≥ λ, equivalently G''(0) ≥ λ."
```

再加 `ch04-gap.tex:23`（G2.0 恒等层）：`G''(k→0) = E[H|0] − Var(∇S|0) = V_eff'' = m_eff²(k→0)`，
以及 Ward 恒等式 `G″·Var(a) = 1`（⇔ 逆感受率）。**⇒ `G″(0)` = 逆磁化率 = `m_eff²(k→0)`（格点单位）。**

**两条必须钉死的量纲/记账事实**：

1. **它没有"物理质量"的量纲**：`ch18-scalar-h1.tex:81–83` 逐字 —— "the curvature floor is a lattice-unit object,
   `G″_a = a²m_phys²`, so **every δ in this chapter is a dimensionless curvature floor and none is a physical mass bound**" `[T]`。
2. **`G(a)` 本身从未写出**（只有 `G″`）：`GPP-STRUCT.md:52` —— "可核的定义是 `G″`，不是 `G`"。故本文件所有陈述都以 `G″` 为对象。

---

## §2 四条候选路径：每条"需要什么额外假设"

| # | 形式链 | 除 `G″(0)>0` 外还**需要**什么 | 出处 |
|---|---|---|---|
| **P1** | `G″>0 ⇒ f₁²>0 ⇒` 原子 ⇒ 隔离 ⇒ 双指数 ⇒ 隙 | 需要"次阈值**原子**存在"= `Λ★>m₁²` = I2 = **Clay 上一级**（`IMP-ARROWS:66` 行 14、`:128` `d₂`） | `ch06-tools.tex:31–48`（`f₁²(μ)=Z_O²f₁²(μ₀)` 乘性残数） |
| **P2** | `G″>0 ⇒ 某关联函数不衰减` | 需要"**该道**就是 Clay 要的道"：`G″` 活在 `a`/`S_W`-力道，Clay 要 `trF²`（0⁺⁺ 中性）道 ⇒ 需要 **e₃ 扇区转移** | `IMP-ARROWS:60`（行 8 KL 双向）、`:113`（`e₃`） |
| **P3** | `G″>0 ⇒ 转移矩阵谱隙` | 需要 **`inf_{(L,a)}` 一致**（有限 `V` 的 `Δ_V>0` 是 PF/Jentzsch 平凡；`inf_V Δ_V>0` 逐字 = 隙） | `ch08:19–21`、`IMP-ARROWS:140–145`（链 III） |
| **P4** | `G″>0 ⇒ ⟨¼trF²⟩>0 ⇒ ?` | 需要 `⟨¼trF²⟩` 与 `ν` 原子之间的**已证通道**（OPE 系数？减除？）—— 语料**一处也没写过** | `ch06-tools.tex:74–75`；`IMP-ARROWS:57`(iv) |

---

## §3 逐条判

### P1 `G″>0 ⇒ f₁²>0` —— **`[NG]`：这条不是"缺证"，是"反了"**

- **反例是 `[T]`**：2D 纯 YM 有精确面积律 `σ_R ∝ C₂(R)>0`、谱离散有隙（⇒ H1 成立 ⇒ 由 §1 的等价 `G″>0`），
  而 **`f₁²(trF²) = 0` 精确**（`trF² = trE² = H` 守恒 + 特征标正交）。`IMP-ARROWS:92–95` 逐字；底料 `R5-EXACT-ANCHORS.md` §1.1。
- **机检第二反例**：`H1 ⇏ atom`（thin atom `wδ_ε+(1−w)δ_{s₂}`，`w=ε/2λ`，`χ=1/λ`，`inf supp=ε→0`），`IMP-ARROWS:33–35`。
- **方向也错级**：`f₁²>0` + 隔离 = `B_atom` **严格强于** `B_gap`（`IMP-ARROWS:33`）。`G″>0`（= `B_gap` 级）
  去推 `B_atom` 级，是把结论升级当前提。
- **⇒ 判：`[NG]`（对"RP 格点规范理论 + 曲率正性"类），且需要的那条假设本身就是 Clay 上一级。**

### P2 `G″>0 ⇒ 某关联函数不衰减` —— **`[T]` 同义反复 + `[O]` 跨道**

- **同道内是双向 `[T]`**：Källén–Lehmann `m_eff²(k→0)≥λ ⟺ ⟨OO⟩_c ≤ Ce^{−λ|x|}`（逐道，`IMP-ARROWS:60` 行 8）。
  ⇒ 这半边**零增益**：`G″>0` 与"该道关联不衰减"是同一句话。
- **跨道是 `[O]`**：`G″` 的条件化对象是粗幅 `a`/`∂_v S_W`；Clay 要 `trF²` 道。
  运动学桥被 `Z_N^{(1)}` 荷分级**排除**（`[T]`，`IMP-ARROWS:76`）；动力桥唯一实例是 `u`-同一性（同域锁定，`IMP-ARROWS:117–122`）。
- **"根本没写过"的那句**：**没有任何文件写过"粗模 `a`-道的正曲率蕴含 `trF²`-道的谱测度在 0 以上有支撑"**。这是 `e₃`，不是新缺口，但它在这里第一次以"蕴含"形式被点名。

### P3 `G″>0 ⇒ 转移矩阵谱隙` —— **`[T]` 改名 + `[NG]` 循环**

- 改名是 `[T]`：`ch08:19–21` 的四个名字逐字等价；`ch18:50–70` claimbox 第 1 条 "``H1'' without a δ **is** the statement `G″(0)>0`"。
- 循环是 `[NG]`：链 III（`IMP-ARROWS:137–145`）缺口 `k₂ = inf_V Δ_V > 0`，**逐字就是隙**。
- 剩下的真内容只有一处：`(L,a)` 一致。**这正是 (b) 的断点。**

### P4 `G″>0 ⇒ ⟨¼trF²⟩>0 ⇒ ?` —— **前半 `[NG]`/未写，后半 `[O]`**

- **2D 已否证"凝聚 ⇒ 原子"**：`trF² = H`，`⟨H⟩>0` 对有隙的 2D YM 成立，而 `f₁²(trF²)=0`（`IMP-ARROWS:92–95`）。
- **4D 的 `⟨¼trF²⟩` 是 OPE 量**（需减除、`[C]` 级），与 `ν` 原子之间**无已证通道**；语料零处写过后半箭头。
- **"同一个对象的多个名字"是纲领措辞，不是定理**：`CLAUDE.md §5` 把 `f₁²>0` / H1 / `(CP)` / `⟨¼trF²⟩` 并列成"同一个对象"，
  但**已证的只有前两者之间的恒等（§1）**；`⟨¼trF²⟩` 与它们的同一性是 `[C]`（OPE 图像），不可当 `[T]` 用。
- **⇒ 判：P4 无可判的链。它既不是 `[T]` 也不是 `[C]`，是 `[O]` + 一个 `[NG]` 的 2D 否定。**

---

## §4 最短可行链（每步给等级）+ 断点

### §4.1 最短链：**3 步，两步是改名，一步是墙**

```
(0) [T] 定义层   G″(0) := V_eff″(a=0) = E[H+Φ″|0] − Var(G+Φ′|0) = Δ_dyn + W_geo + X
                  = 1/Var(a) = m_eff²(k→0)（格点单位）        ch05:230–240 / ch15:23–28 / ch04:23 / Ward
(1) [T] 改名     G″(0) > 0  ⟺  H1(0)                          eq:scalar-rule / eq:h1st-rule（纯代数）
(2) [T] 改名     H1(0) ⟺ m_eff²(k→0) > 0（固定 L,a）          ch08:19–21 + PROOF-IR-H1（自标 rename）
─── 到此为止内容长度 = 0：证的是 H1 ⟺ H1 ───
(3) [M] ★墙★     符号 → 带余量的 (L,a)-一致下界：
                  inf_{L,a} [ G″_{L,a} − W_full,{L,a} − a²Δ_0² ] > 0
                  ⟺ Δ_dyn ≥ a²Δ_0² 一致（eq:h1st-dyn）
                  ⟹ m_eff² ≥ λ > 0 一致 ⟹ Clay 的隙
```

**判**：**最短链长度 = 1 步非平凡**（(3)）。`G″(0)>0` 不是走向 H1 的第 0.5 步——**它与 H1(δ=0) 是同一句**。
把 `G″>0` 记账成"H1 路线的进展"是把 `δ_max>0` 冒充 `δ_max ≥ a²Δ_0²`。

### §4.2 断点逐字（这是第 4 级台阶的精确形态）

断在 **(2) → (3)**，一句话：

> **`W_full ≥ N/12 > 0` 把 `G″` 与 `Δ_dyn` 隔开了一个与测度、`β`、条件点都无关的正数。**
> `[T]`：`ch05:203–206`（`eq:wgeo`：`W_geo = E[Φ″] − Var(Φ′) ≥ 0`）、`ch05:214–219`（`eq:wgeo-quant`：`W_geo ≥ λ(v) ≥ (N/12)‖v‖² > 0`，Prékopa + `NS-NOSUBTRACTION` Prop. N1）。
> 于是 `G″>0 ⇔ Δ_dyn > −W_full`，**`Δ_dyn` 可以是负的**，而 `Δ_dyn ≤ 0` 的物理来源是已知的：
> Savvidy 大场凹性 `−D² + 2adF`，`λ_min = −g|B| < 0`（`PROOF-SAVVIDY-EXACT-SIGN`；`W4` 规范商后 7.1% 负）。

**数值把这条记账钉死**（`ch19:125–149`，L=2, β=2.5, SU(2)）：

| 量 | 值 | 读数 |
|---|---|---|
| `G″(0)` | `+1.696 ± 0.293`（5.8σ） | 符号已 `[T]` 级（Q11/R06/R10） |
| `W_full` | `+1.134 ± 0.063` | **= 67% of `G″`**；其中 `E[Φ″] = +1.2325 ± 0.0189` 占 **109%** |
| `Δ_dyn` | `+0.562 ± 0.277` | **只有 2.0σ** |

⇒ **`G″>0` 的 5.8σ 里有 67% 是运动学（Jacobian）；真正的动力学余量 `Δ_dyn` 只有 2.0σ，且未被任何一致化证明覆盖。**
`ch19:142–149` 已把这条列为 trap 3（"lattice gives H1 at 5.8σ" = 悄悄用 `eq:h1st-rule` 替换 `eq:h1st-dyn`、丢掉 `W_full`）。

**酸测的位置（本文件把它接到蕴含链上）**：

- **过**：自由无质量场 `G″=0`（`ch06-tools.tex:11–16`：`f₁²=0, G″=0`；`CAVITY-MEAN-FIELD:21`：`Φ` 平、Gaussian 等分 ⇒ `Δ_dyn=0`）。⇒ 符号**不是**平凡真。
- **不过**：`β=0` 的纯 Haar（无任何 YM 动力学）上 `G″ = 0.248 > 0`，`E[Φ″]−Var(Φ′) = 0.246`（`ch05:209–211` 逐字 "an exact closure"）。
  ⇒ **符号本身不携带动力学信息**；`CAVITY-MEAN-FIELD:11` 已写 "it is the **Jacobian**, present before any YM dynamics … its positivity is non-dynamical"。
- **⇒ 判**：`G″>0` **通过"自由场"酸测、失败"纯 Haar"酸测**。因此它可以是**判据**（`H1` 的等价形式），不能是**路标**（不能当"已过动力学关"记）。

### §4.3 `(c)` 的正面回答：类内反例存在，且它不改目标、只改量词

**取 `G_YM × U(1)_weak`（U(1) 在 Coulomb 相、耦合足够弱）**（`IMP-ARROWS:83–90` §2.1(a′)）：
`G_YM` 因子有面积律、RP 仍成立、其 `a`-道 `G″>0`（YM 因子内部有隙 ⇒ 由 §1 等价 H1 成立），
但**全理论的谱无隙**（光子扇区）⇒ **`G″>0` 与 `m=0` 并存**。

- **这意味着**：裸蕴含 `G″>0 ⇒ m>0` **在"RP 格点规范理论"类内为假**；任何证明都必须用到 **"无退耦中性无质量扇区"（纯性）**，
  这一条与 `e₃`（扇区转移）是同一堵墙的两个面。
- **这并不意味着** H1 路线重新定位：对**纯** 4D YM，恒等式层（§1）使 `G″>0 ⇔ m_eff²(k→0)>0` 为 `[T]`，
  所以 `G″>0` 而 `m=0` 在纯 YM 内需要 `inf_{(L,a)}` 失效（有限体积有隙、无穷体积无隙），**那是 (b) 而不是 (c)**。
- **⇒ 结论**：**(c) 的形态不成立（举不出纯 YM 内反例），(b) 的形态成立（断在一致化 + 余量 + 扇区）。**

---

## §5 对"第 4 级台阶"的定义（可直接当任务书用）

> **Gate-4（精确缺失步 `[M]`）**：证明存在 `λ>0`，使对物理点附近的一切 `(L,a)` 与一切粗模 `k→0`，
> `G″_{L,a}(0) − W_full,{L,a} ≥ a²Δ_0² + λ`，即 `Δ_dyn ≥ a²Δ_0²` 一致；且该下界经 `a→0` 后以 `m_phys² = G″_a/a²` 的形式存活。
> **等价改写（可判定的失败面）**：`inf_{L,a} Δ_dyn ≤ 0` 或 `sup_{L,a}[W_full/Δ_dyn] = ∞`。
> **不得使用**：`G″>0`（符号）、`Var(G)≤C`（等分律下 `Δ_dyn≈0`，`C7-SPECIES-CONSTRUCT:6/43` 已判"恒等式形状的死路"）。

**两条配套警示（防止把 Gate-4 又降级成 (a) 的改名）**：

1. **不要把 `W_geo` 当余量**：`W_geo ≥ N/12` 是 `[T]`，但它与 `Δ_dyn` 同时出现在 `G″` 里，**加号方向相反于需要**（它掩盖 `Δ_dyn<0`）。
2. **不要把"L=1 单 plaquette 的 `[T]`"当理论结果**：Q11 的证书对象逐字是"单 plaquette、`SU(2)`、`θ=π`、切片 `ξ₃=0`、`v=e₃`、`‖v‖=1`"（`Q11-CERT.md` §0）；
   `ch18:17` 逐字："What stays open is the step from that one-plaquette (patch) model to the theory `[M]`"。

---

## §6 检索执行核验（AGENTS.md §1）

| 项 | 记录 |
|---|---|
| 工具状态 | `pwsh`（PowerShell `Select-String`/`Get-Content`）**OK**；`read` **OK**；`glob` **OK**（`**/*ch05*` → 1 命中）。**无 `ENOENT`、无 `ripgrep provider failure`、无超时。** |
| 检索式 1 | `Select-String -SimpleMatch 'eq:split'`（全仓 `*.tex,*.md`）→ **16 行 / 7 文件** |
| 检索式 2 | `Select-String -SimpleMatch 'eq:v3object'` → **9 行 / 5 文件** |
| 检索式 3 | `Select-String -Pattern 'Δ_dyn|Delta_\{\\mathrm\{dyn\}\}'`（`research/**/*.md`）→ **大量命中，输出被截断**（完整日志 `…\Temp\dsh-spill-…\1ee76d76059a-pwsh.txt`）；**记 TRUNCATED，不记 0** |
| 检索式 4 | `Select-String -Pattern "G''\(0\) *> *0|G″\(0\) *> *0"`（`research/v3/species/*.md`）→ **25 行 / 17 文件** |
| 检索式 5 | `Get-ChildItem -Filter 'p21*/q11*/r06*/r10*'` → 定位 `P21-REBUILD / Q11-CERT / R06-GAPMONO / R10-AUDITCERT` + `w7/*.py/.log`，**命中** |
| 未查到（采信） | **"`G″(0)>0 ⇒ f₁²>0`"、"`G″(0)>0 ⇒ ⟨¼trF²⟩>0 ⇒ ·`"、"粗场道 `G″>0` ⇒ `trF²` 道在 0 以上有支撑"这三条蕴含**（检索式 3+4 的并集上人工通读，命中 0 条以蕴含形式写出的句子；工具 OK） |

---

## §7 出处（逐条可核）

**内部（行号已核）**：`docs/paper-v2/ch05-theorems.tex:203–206, 209–211, 214–228, 230–240` · `ch15-uvface.tex:20–28` ·
`ch18-scalar-h1.tex:17, 32–48, 50–70, 77–83` · `ch19-h1-status.tex:32–47, 82–93, 125–149` · `ch08-immune.tex:11–21` ·
`ch04-gap.tex:23` · `ch06-tools.tex:11–16, 31–48, 74–75` · `research/v3/species/Q11-CERT.md §0` · `R06-GAPMONO.md:11–20` ·
`R10-AUDITCERT.md §0` · `IMP-ARROWS.md:33–35, 57, 60, 66, 76, 83–95, 113, 117–122, 137–145` · `MU5-SIGMA-M.md:14–18` ·
`C7-SPECIES-CONSTRUCT.md:6, 43` · `CAVITY-MEAN-FIELD.md:11, 21` · `GPP-STRUCT.md:52`.

**外部（本次未开卷，按语料著录，`[B]` 级）**：Källén 1952 / Lehmann 1954（KL 表示）· Prékopa 1973 · Brascamp–Lieb 1976 ·
Kotecký–Preiss 1986 · Chatterjee, arXiv:2006.16229（`MU5` 已逐字核摘要：方向为 `m>0 ⇒ σ>0`）· Savvidy, PLB **71** (1977) 133。

---

*本文件不证任何新定理、不改任何数值、不 QED。唯一新增的"蕴含级"陈述是 §4.2 的两句：`G″>0` 不蕴含 `Δ_dyn>0`（`[T]`，由 `W_full ≥ N/12` 直接得出），
以及与之配套的 Gate-4 定义。**Clay 仍 0%。***
