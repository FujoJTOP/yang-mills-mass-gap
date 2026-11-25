# AUDIT-CH16 — `docs/paper-v2/ch16-sixth-species.tex` 逐数字溯源审计

> 审计员 **ch16-audit**，2026-11-24。**只读**（未改任何 `.tex`，未 git commit），本文件是唯一新建文件。
> 对象：`docs/paper-v2/ch16-sixth-species.tex`（960 行，64 979 B，mtime 23:40:12）。
> 一手源（唯一合法出处，共 13 份 + `w7/*.log` 24 份）：
> `ROUND1-VERDICTS.md`、`AUDIT-R1.md`、`AUDIT-C0.md`、`RB-RATIO-BOUND.md`、`WA-ALIGNMENT.md`、
> `QI-BRACKET.md`、`TCI-COMPENSATION.md`、`NA2-NONDEGENERATE.md`、`BN-RESIDUAL.md`、`LIT-CHECK.md`、
> `LEAD-NOTE-NA-SU3.md`、`SIXTH-SPECIES-CHARTER.md`（实际路径 `research/v3/`）、`SIXTH-SPECIES-SYNTHESIS.md`（仅作旁证，**不计入**判定）、
> 以及 `research/v3/species/w7/*.log`。
> 方法：`\times10^{k}` 归一成 `ek`、源侧 `−`(U+2212)→`-`、指数零填充归一后，对**每个**数字串做逐字 grep；
> 三态判定 OK / SHIFTED / MISSING，另设“算术复核”与“OMISSION”两类。

---

## §0 结论（先看这里）

- **提取**：291 个数字 token（去重；`\to` 后的数已修复抽取），其中 **148 条**为实质数值断言（其余是表宽 `2.05/6.35/7.6/4.1cm`、`N` 取值、网格/样本数与编号）。
- **计数**：**OK 140 ／ SHIFTED 4 ／ MISSING 2 ／ 算术复核失败 2 ／ OMISSION 3 ／ 标签问题 10 ／ 引用问题 5**；另有 **逻辑/范围失败 2**（§3 第 6 条“every anchor”量词、第 9 条“SU(N≥3) 未测”），不计入上面 148 条数值口径。
- **追加高危点（lead 指令，详见 §8）**：①TR 混模型 —— **ch16 未犯**（`0.6444↔3.013`、`0.267884↔10.5` 配对正确；错的是 `SYNTHESIS §1`）；②`e^{−0.4595N}` 的 `\LE` —— **未逐处降标**（两处无标、一处合标，且邻近 `\LE` 错挂在“min≡0”上）；③`floor⇒ratio` —— **犯了**：引了 RB-6(iv) 的等价并写“given (ii) … yields ρ²≤C²”，**全章无 DT 限制** ⇒ OMISSION（重）；④SW 两句 —— **未犯**（“At fixed x”“At the mode”分开写）。
- **最严重三条**
  1. **`1.56` / `2.43` 在一手源里不存在**（§1 第 1、2 行）：chapter 把它当作“C=2 的物理决定性论据”，并把出处挂到 `Z1-HOLO.md`；该文件只写 lattice `0^{++*}/0^{++}≈1.6–1.8`（模型自身给 `1.639`），**从无 1.56**。1.56 的真实来源在别处（`RATIO-RHO.md` / `E2-FLIP-KINEMATIC.md` 的 `ρ_lat≈1.54–1.56`），故属**错引 + 走形**，非凭空造数。
  2. **`175.5757 / 175.5762` 挂错了公式**（§2-S1）：一手源是 **`q=−c`** 的 `(4π²−c)/(π²−c)`（`dt_dirichlet.log:71–73` 逐字 `lam_k = k^2 pi^2 - c`）；ch16 印的是 `(4π²+c)/(π²+c)`，`c=9.7` 给 `2.513`，**不是 175.58** ⇒ “FD 确认该式”这句自证不成立。
  3. **QI 精确恒等式缺失**（§3 第 2 条）：`QI-BRACKET.md` 的核心是恒等式 `ϑ=w₁−Δ/ε`，即 `w₁=(s₂−Ē)/(s₂−s₁)+Δ/(s₂−s₁)`，`Δ=∫_{[s₂,∞)}(s−s₂)dν≥0`；ch16 只留不等式 `w₁≥(Λ−Ē)/(Λ−s₁)`，slack 的 UV 内容整段丢失 ⇒ **OMISSION（中）**。
- 另：**末段“SU(N≥3) 未测”被一手源正相反地否证**（§3 第 9 条）：`NA2-NONDEGENERATE.md`（23:41，晚于本章 23:40）已在一 plaquette SU(3) 上测出 `X(a₀=0)≠0` 且随边界 holonomy 变号。
- 正面：五个物种的死因数字、数据模缺陷、`C=2`、`1/48`、(TCI) 全套数字**逐字可核**，未发现伪造量级。

---

## §1 MISSING 清单（一手源里根本没有这个数）

| 数字 | ch16 位置 | 一手源里最接近的是什么 |
|---|---|---|
| `1.56` | 行 582–583 §“Why the correction is forced by the physical spectrum”：*“$m_{0^{++*}}/m_{0^{++}}\approx1.56$ (the corpus records the lattice range $1.6$–$1.8$ in `Z1-HOLO.md`; $1.56$ is the conservative end)”* | `Z1-HOLO.md:78` 只有 lattice `0^{++*}/0^{++}≈1.6–1.8` 与硬壁模型预测比值 `{1, 1.639, 2.263, 2.881}`；**无 1.56**。1.56 只出现在 `SIXTH-SPECIES-SYNTHESIS.md:146`、`research/unconditional/RATIO-RHO.md:3,21`、`E2-FLIP-KINEMATIC.md:45`（`ρ_lat≈1.54–1.56`）——**均不在本章一手源内**。且 `1.56` 低于被引区间下端 `1.6`，把它称“conservative end”与所引文本不符（区间下端是 1.6，模型给 1.639） |
| `2.43` | 行 585 eq. `\eqref{eq:physwindow}`：*“$\rho^2_{\mathrm{phys}}\approx1.56^2=2.43>2.25$”* | 纯由 1.56 平方得来（算术对：`1.56²=2.4336`），**但被乘数无一手出处** ⇒ 该行 eq 与行 588 “the theory violated the window”的**定量依据缺失**。（若改用一手源给的 1.6：`1.6²=2.56>2.25`，结论仍立、数字须换。） |

> 口径说明：这不是“凭空编数”，是**把语料别处的数挂到了本章一手源名下**，并顺手把区间读反。审计按任务口径记 MISSING。

---

## §2 SHIFTED 清单（数字在源里，但被用到另一个量/另一个模型上）

| 数字 | 源里的含义 | ch16 里的含义 | 严重度 |
|---|---|---|---|
| **S1** `175.5757` / `175.5762` | `rb_ratio_bound.log:71–73`：`--- constant q = -c : exact lam_k = k^2 pi^2 - c ---`，`c= 9.70 exact rho^2= 175.5757 FD= 175.5762`；即 **`q=−c`** 下的 `ρ²=(4π²−c)/(π²−c)`（`c<π²`，`c↗π²` 发散） | 行 552–563：*“For $-u''+cu=\lambda u$ … $\mu_k=k^2\pi^2+c$, hence $\rho^2(c)=(4\pi^2+c)/(\pi^2+c)$ … The finite-difference check at $N=4000$ **confirms** (eq. dirichletratio) ($c=9.7$: exact $175.5757$, FD $175.5762$)”*。**`+c` 版在 `c=9.7` 给 `2.513`**，与 175.58 差 70 倍；`ρ²≤4⟺c≥0` 在两版都成立，所以错的是**符号约定 + “confirms”** 这一句 | **高**（被印成定理的恒等式，且用它“确认”自己） |
| **S2** `4.48 → 807.36`，`4.34 → 20.05`，`log ρ²≈0.558g` | `RB-RATIO-BOUND.md:275–285` 的 **SUSY 双井** `H=−∂_x²+W`，`W=V'²/4−V''/2`，`V=(g/2)(x²−1)²`，FD `N=6000`，`g=0.5` 处 `μ₁=0.79209`；`dt_dirichlet.log:139` 独立复现（`g=2: 0.74868`） | 行 811–813：*“**The same model** kills any upper-bound-data machine … ρ²=4.48→807.36 … T₂ grows only 4.34→20.05”*，承接的上表（行 794–800）是 `H=−½d²/dx²+(g/2)(x²−1)²`、`g=0.5` 处 `s₁=0.7250` 的**另一个**双井。`0.7250≠0.79209` ⇒ 两个模型 | **中** |
| **S3** `2 of 320 self-loops tie` | `tr_maxplus.log:17`：`critical graph: 2 of 320 self-loops tie at rho`；但**同一日志 27 行**：`=> every self-loop ties at rho_inf; the tropical spectrum is n-fold degenerate`；`AUDIT-R1.md:139`：`64/64 self-loops tie`，且 `best 2-cycle mean = 2.487962 < 2.5`（`deficit=β(1−cos2π/n)`） | 行 180：*“the critical graph is the whole diagonal: two of 320 self-loops tie at the same cycle mean, the max-plus spectrum is n-fold degenerate”* — “whole diagonal”与“two of 320”**自相矛盾**；“2”其实是**最佳 2-圈**（`2.487962`），不是“两个自环并列” | **中** |
| **S4** `10.5` | `AUDIT-R1.md:136–137`：`width/Δ = 4.2552`（`n=64`）、`10.4651`（`n=4096`）；`ROUND1-VERDICTS.md:76` 只记 `10.5` | 行 190：*“an exact Δ=0.267884 that is n-independent (n=64 and n=4096 give the same value), **with width/Δ=10.5**”* — 在强调“Δ 与 n 无关”的同一句里，引了唯一与 n 有关的量却**不标 n**（实为 `n=4096`；`n=64` 是 4.255） | **低** |

---

## §3 不等式方向与恒等式的核对（逐条 OK / FAIL）

| # | ch16 的断言 | 一手源 | 判 |
|---|---|---|---|
| 1 | 行 713 eq.`l2`：`w_1 ≥ (Λ−Ē)/(Λ−s_1)`，前提“`ν` 在 `(s_1,Λ)` 无质量”，两点测度取等 | `QI-BRACKET.md:127`（`L2`：`w₁≥(Λ−Ē)/(Λ−s₁)`，任意 `Λ≤s₂`，**本文不改不等式**）+ `:27`（两点族 `ϑ≡w₁` 精确）。方向、前提、锐性**全对** | **OK** |
| 2 | **精确恒等式** `w_1 = (s_2−Ē)/(s_2−s_1) + Δ/(s_2−s_1)`，`Δ=∫_{[s₂,∞)}(s−s₂)dν≥0` | `QI-BRACKET.md:18–23`、`:111–119`：`ϑ := (s₂−Ē)/(s₂−s₁) = w₁ − Δ/ε`（“**恒等式，非不等式**”，数值 `max|ϑ−(w₁−Δ/ε)|=3.6e−15`）。**ch16 全文无 `ϑ`、无 `Δ`（传输隙的 Δ=0.6444 与地板的 Δ 是另一个符号）、无 `Δ/ε`** ⇒ 不等式在、恒等式**不在** | **FAIL（OMISSION，中）** |
| 3 | 判据形态 `ϑ≥w₀ ⟺ Δ ≤ (w₁−w₀)ε`（“slack 是 UV 敏感的超额一阶矩”） | `QI-BRACKET.md:159–163`、`:291`；ch16 只写 `s₂ ≥ (Ē−w₀s₁)/(1−w₀)`（行 722–725），**没有任何一处把 slack 与 Δ 挂钩** | **FAIL（OMISSION，轻）** |
| 4 | 行 767–772 eq.`eta`：`ρ²≤C² ⟺ μ_1 ≥ t/(1+C²) ⟺ η≥1` | `RB-RATIO-BOUND.md:186–192`（RB-2″，`[T]`）、`:90–96`（D3 定义 `η`）。**方向、等价、阈值的可算性全对** | **OK** |
| 5 | 行 516：`C=2`，不是 `C≈1.5`；表 5 行违反 `2.25`、2 行饱和 `4` | `RB-RATIO-BOUND.md:308–323`（RB-8：Mehler 2、plaquette `8/3`、2D YM `SU(2)` `8/3`、path `N=2` 3、盒/转子/弦 4、2D `U(1)` 4）+ `:419`（PART F：`SU(3) 9/4`）。逐行对齐 | **OK** |
| 6 | 行 526 的**全称**：*“That value is violated by **every** exact anchor large enough to test it”* | 同表自身就有 Mehler `2`（过）与 2D `SU(3)` `9/4=2.25`（**boundary, pass**）；`RB-RATIO-BOUND.md:321` 只说“被 **5 个**锚点违反”。全称量词与自己的表冲突 | **FAIL（量词过强，中）** |
| 7 | 行 573–576：代价 `μ_1 ≥ Δ/(C²−1)`，`Δ/1.25 → Δ/3` | `RB-RATIO-BOUND.md:323`（`Δ/1.25` 降到 `Δ/3`）、`:219–221`（`t/5=Δ/3`，`C=1.5→0.800` vs `C=2→0.333`）。方向（**变弱**）与数都对 | **OK** |
| 8 | 行 460–484：`X≡0`（`a₀=0`）只是 **SU(2)** 径向 Jacobian 的偶然；`(NA)` 的内容只在 `SU(N≥3)`；`∂_vI_3=3tr(T_v(ξ_⊥·T)²)≠0` | `LEAD-NOTE-NA-SU3.md:10–17`、`:36–40`；`NA2-NONDEGENERATE.md:13–15`。**范围限定与公式都对**（记法见 §5-R5） | **OK** |
| 9 | 行 934–936：*“the $\mathrm{SU}(N\ge3)$ version of $(NA)$ **was not tested**, because no solvable model with a Weyl-symmetric (non-radial) Jacobian was built”* | `NA2-NONDEGENERATE.md:18–22`：*“实测（800k Haar 样本 + 1.2M 拟合交叉验证）：**SU(3) 一 plaquette 上 `X(a₀=0) ≠ 0`**，量级 `1e−4…1.7e−2`，且**符号随边界 `B` 变**”*；`:3` 标题即“`a₀=0` 的恒零是 rank-1（SU(2)）专属”。**该句被一手源直接否证** | **FAIL（与一手源相反，高）** |
| 10 | 行 659–662：自由场上 `ℋ₋={0}` ⇒ `∀φ∈ℋ₋` **空真**，不是“返回 0” | `TCI-COMPENSATION.md:313`（`ℋ₋={0}` ⇒ **EMPTY** ⇒ holds **vacuously**；`does not "return 0"`）、`:440` | **OK** |
| 11 | 行 307–311：`d_N~(e^{−0.4595N})`，末三点回升，不可外推；`[T]` 只覆盖零方向 | `dm_floor_defect.log`（`fitted log-slope -0.4595`、`N=12/13/14: 1.0728e-2/1.1192e-2/1.3634e-2`、`strictly positive…True`）、`AUDIT-R1.md:272`（C-5） | **OK** |
| 12 | 行 178–183：常模切片 `S_W≡0` ⇒ `Δ_trop≡0`（每个 β），对 `Δ=0.6444` | `ROUND1-VERDICTS.md:71`、`tr_maxplus.log:17,27`、`AUDIT-R1.md:143–145` | **OK**（“2 of 320”见 §2-S3） |
| 13 | 行 631–632：`(TCI) ⇏ Σ⪰μ`，Schur 更强，差距恰为 `‖B‖²/γ` | `TCI-COMPENSATION.md:164–168`、`:202`（`COMP ⇒ ELIM` 恰付 `‖B‖²/γ`） | **OK** |
| 14 | 行 555–557：`ρ²(c)=(4π²+c)/(π²+c)`，`ρ²≤4 ⟺ c≥0`（sharp） | 等价式本身对两版都成立；但**被它“确认”的 175.5757 属 `q=−c` 版** ⇒ 见 §2-S1 | **FAIL（符号约定）** |
| 15 | 行 722–725 eq.`bridge`：`w_1≥w₀ ⟸ s_2 ≥ (Ē−w₀s_1)/(1−w₀)` | 代数复核：`(s₂−Ē)/(s₂−s₁)≥w₀ ⟺ s₂(1−w₀)≥Ē−w₀s₁` ✓；与 `WA-ALIGNMENT.md:223–226` 一致 | **OK** |
| 16 | 行 700–701：`w_1 ≤ Ē/s_1`、`Σ_{n≥2}w_n ≤ Ē/s_2`（都是**上**界） | `WA-ALIGNMENT.md`（方向判决 §2.2；`Ē` 为一阶矩 ⇒ 上界机）；代数复核 ✓ | **OK** |
| 17 | 行 750–753 eq.`ratio`：`ρ² ≤ T₂(V)/μ_1 − 1`（对任意测试对） | `RB-RATIO-BOUND.md:160–167`（RB-1′ 2.2）、PART H（`2000 random V`、`min slack 2.53`、`1e-9` 取等） | **OK** |
| 18 | 行 732–735：自由场上桥“不弱而是**未定义**” | `QI-BRACKET.md:199`（`ν=δ_∞` ⇒ 桥未定义 `∞/∞`）、`:66` | **OK** |
| 19 | `×5.6`：`T₂: 4.34→20.05` 是 `×5.6` | 源亦写 `×5.6`（`RB-RATIO-BOUND.md:293,420`），但**算术复核**：`RB-RATIO-BOUND.md:279,284` 的 `T₂` 列 `4.34099 → 20.05462` ⇒ 比值 **`4.62`**（`dt_dirichlet.log:139` 同模型 `g=12: T₂=20.005`）⇒ 应为 **`×4.6`** | **FAIL（算术，低—中）** |
| 20 | 行 440：*“magnitudes are 10^{-2} to 10^{-1}, **nine** orders of magnitude above the quadrature agreement（1.5e-13）”* | 派生量：`1e−2/1.5e−13 ≈ 6.7e10`，`1e−1/1.5e−13 ≈ 6.7e11` ⇒ **10–12 个量级**，不是 9 | **FAIL（算术，低）** |
| 21 | 行 191：`βδ_2=8.7e−4` vs `Δ=0.6444`，*“short by 750×”* | `ROUND1-VERDICTS.md:77`（差 750×）。复核 `0.6444/8.7e−4 = 740.7` ⇒ 750× 是一位有效数字的取整，可接受 | **OK（≈7.4e2）** |
| 22 | 行 643–645：`2.72×` short、`P(μ_3≥0)=0.070`、`β_c≈0.76` | `TCI-COMPENSATION.md:332`（`2.723`）、`:339–343`；`tci_comp.log` | **OK** |
| 23 | 行 602–604：`β≈2.5` 是 `120×` 修正半径；`λ_tH=1` 时 `48/N`，即 `24×`(N=2)/`16×`(N=3) | 算术：`2.5/(1/48)=120` ✓；`(1/N)/(1/48)=48/N` ✓ | **OK** |
| 24 | 行 564：215 形状电池 `sup ρ²=3.999999`、`inf ρ²=1.00000` | `RB-RATIO-BOUND.md:252–257`（`15+200=215`）、`dt_dirichlet.log`（`[PASS] A2 … 3.999999`、`1.000000`） | **OK** |
| 25 | 行 582–583：*“$1.56$ is the conservative end”*（of the cited `1.6–1.8`） | 被引区间下端是 `1.6`；`1.56` 在区间外 ⇒ 方向的读法错 | **FAIL（读法，中）** |
| 26 | 行 562：*“$c\ge0$ means … i.e.\ $e_1\ge\pi^2/L_s^2$ --- a floor, not a degree condition”* ＋ 行 859–861：*“given (ii), the rank-two trace (ratio) yields $\rho^2\le C^2$ with $C=2$”* | `DT-DIRICHLET.md:39–40`、`:161`、`:270–272`：**floor ⇏ ratio**——“`e₁=18.0963≥π²` 时优化器给 `ρ²=27.1846`”（局部变号见证 B）与“非局部正项 `ρ²=4+M/π²=4.5,9.0`、rank-r `4,9,16,24,…`”（见证 A）；`DT-3` 的合法形式是**双侧**比较 `c₋L₀≤L≤c₊L₀ ⇒ C=2√(c₊/c₋)`。ch16 既引了 RB-6(iv) 的等价、又下了 `(ii)⇒ρ²≤C²` 的结论，却**一处也未提 DT 的限制** | **FAIL（OMISSION，重）** |

---

## §4 标签纪律（位置 | 现标签 | 应标签 | 理由）

> 标签语义（`main.tex:37–42`）：`\LT`=T 定理/恒等式，`\LC`=C 条件，`\LM`=M 精确缺失步，`\LN`=N 具名类 no-go，`\LE`=E 数值证据。
> 全书**无 `[J]`（判断）宏**，判断类请落到 `\LC` 并写明是读法。

| # | 位置 | 现标签 | 应标签 | 理由 |
|---|---|---|---|---|
| **A** | 行 199–200：*“Holds … \LT{} for the max-plus spectrum statement, \LM{} for the compactness reading \LC”* | **`\LT`+`\LM`+`\LC`（一句三标）** | `\LT`（“max-plus 圈均值谱对隙盲”）+ `\LC`（“紧致性是测度性质、不是作用量性质”= 读法/判断） | `AUDIT-R1.md:293`（C-8）明确要求：`[T]` 只覆盖 max-plus 谱那一句，“紧致性是测度性质”标 `[J]`；`\LM`=精确缺失步，此处**不是**缺失步 |
| **B** | 行 903 台账 “Radius”：*“1/48 (SZZ, N-independent), `λ_tH<1/24` absent from the source”* | `\LT` | `\LT` 只覆盖 SZZ 的 `1/48`（被引定理）；“`1/24` 不在源文”是**书目核查**→ `\LC`（`SIXTH-SPECIES-SYNTHESIS.md:18` 用的是「`[NG]` 于源文检索面 / `[X]` 于语料」） | 把“某串数在已核源文里查不到”标成定理，会把否证面升格 |
| **C** | 行 625–629：*“The **only self-consistent completion** is `Q=3λ_4φ_0²=3cF` … under which `(TCI)⟺λ_min(A_K|H_-)≥cF+μ` \LT”* | `\LT` | 等价式 `(TCI)⟺λ_min≥cF+μ` → `\LT`；**“唯一自洽完成”+七元枚举 → `\LC`** | `TCI-COMPENSATION.md:122` 把 `Q=3λ₄φ₀²=3cF` 标为 `[MINE]`（作者选择），`:128` 称 `(D3)` 是“唯一非显然的选择”——**是选择，不是定理**；`SM5-TCI-HEAD.md:98` 的 `[T]` 只覆盖“该量为正” |
| **D** | 行 307：*“the minimum is $0.0$ at every $N$ and every window \LE”* | `\LE` | `\LT` | 零方向是初等的（同向量的 `δ₀` 竞争者），`AUDIT-R1.md:272` 明确 `[T]` 覆盖这一半；ch16 自己在行 314 写“the `\LT`-grade part … is the elementary zero direction” ⇒ **同章自相矛盾** |
| **E** | 行 896 台账 “Data-mode theorem” | `\LT/\LE`（合标） | 拆三条：`d_N>0`（`N≤14`）= `\LE`；`min≡0` = `\LT`；`~e^{−0.4595N}` = `\LE`（拟合，且源已撤“`N=O(1/ε)`”） | 合标把“有限计算”和“初等定理”混成一格；`AUDIT-R1.md:272`（C-5）要求逐条分开 |
| **F** | 行 813：*“a factor $180$ against a factor $5.6$ \LT+\LE”* | `\LT+\LE` | `\LE`（并把 5.6 改 4.6） | 纯数值比值；且该比值算错（§3-19） |
| **G** | 行 545：*“Five exact anchors violate $C=1.5$; two of them saturate $C=2$ exactly.”* | **无标签** | `\LT`（精确锚点） | 下一句才挂 `\LT/\LE`；本句是全章“常数更正”的承重句，缺标 |
| **H** | 行 646：*“The compensated well … sits at $a_*\approx9$–$14$, outside any small-block amplitude”* | **无标签** | `\LE` | `TCI-COMPENSATION.md:496` 台账标 `[V]`（`tci_comp.log [D]`）；整段其余句都带 `\LE` |
| **I** | 行 287–297 `claimbox{cT}`“Theorem（two-sided blindness…）” | 盒型 `cT`，正文却写“the positive direction is an explicit extremal construction, **computed for all $N\le14$**” | 保留 `cT` 给“双向不可判”的定理壳，**盒内补 `\LE`** 限定“computed for $N\le14$”（或降为 `cC`） | 定理盒里装有限计算，正是 `CLAUDE.md §2` 禁止的那类升级；源 `ROUND1-VERDICTS.md:119` 把拟合标 `[E]` |
| **J** | 行 307–308：*“A log-slope fit over $N=5,\ldots,14$ gives $d_N\sim e^{-0.4595\,N}$”*；行 358：*“(moments $e^{-0.46N}$ **by fit**; …)”*；行 896 台账 | **两处无标签**；台账处 `\LT/\LE` 合标；而紧邻的 `\LE`（行 307 句末）挂在 **“the minimum is 0.0”** 上 | 拟合律两处各补 `\LE`；行 307 句末的 `\LE` 改 `\LT`（零方向）；台账拆成 `\LE`（拟合）／`\LT`（`d_N>0` 与 `min≡0` 的双向性） | lead 指令第 2 条：`d_N~e^{−0.4595N}` 是 `N≤14` 的**拟合**，只能 `\LE`；作者自称“拟合降 `\LE`”，但**逐处核对未落实**，且 `\LE` 被挂到了方向相反的句子上 |

> 另（低）：行 662 *“which is the opposite of a pass \LT”* —— 空真性本身是定理（`TCI` 台账 `[T]` defect），但“opposite of a pass”是判断，建议 `\LT` 只留前半句。

---

## §5 引用与 label 核实（缺失的文件 / 标签列表）

**`\code{...}` 目标逐一 `Test-Path`（结论：全部存在，但多数不在 `research/v3/species/` 下，读者按裸文件名不可解析）**

| # | ch16 引用 | 实际路径 | 判 |
|---|---|---|---|
| R1 | `research/v3/SIXTH-SPECIES-CHARTER.md`（行 7） | ✅ 同名同位（**注意：不在 `species/` 下，ch16 写对了**） | OK |
| R2 | `research/unconditional/SPECIES-EXISTENCE.md`（行 264） | ✅ | OK |
| R3 | `Z1-HOLO.md`（行 583） | ✅ `research/v3/Z1-HOLO.md`，但**内容不含 `1.56`**（只 `1.6–1.8`、`1.639`） | **错引**（§1 第 1 行） |
| R4 | `SM5-TCI-HEAD.md:72,74,91`（行 621–623） | ✅ `research/v3/SM5-TCI-HEAD.md`；**逐行核过**：72=`Q_Sav=−D²+2adF_B`、74=`cF:=g|B|=:c`、91=`cF=g|B|`、98=`Q=3λ₄φ₀²=3c` | OK（路径未给，建议补 `research/v3/`） |
| R5 | `U6-BINDING-LOGIC.md`、`T6-BINDING-ADVERSARIAL.md`、`V3-COUNTEREXAMPLE-SCOPE.md`、`FOSP-A` | ✅ 均在 `research/v3/`（**不在** `species/`） | OK（建议统一前缀） |
| R6 | `IN-11`（行 525） | ⚠️ **无同名文件**；最接近 `research/unconditional/BREAK-IN-11.md`、`IN-11-attack.md`、`GAP-IN11.md` | 记为“登记项名，不是文件名”，建议写 `GAP-IN11.md` |
| R7 | `w7/dm_floor_defect.py`、`ROUND1-VERDICTS.md`、`AUDIT-R1.md`、`LIT-CHECK.md`、`RB-RATIO-BOUND.md`、`TCI-COMPENSATION.md`、`WA-ALIGNMENT.md`、`QI-BRACKET.md`、`SIXTH-SPECIES-SYNTHESIS.md` | ✅ 全部存在（`species/` 下） | OK |

**`\ref`/`\eqref`：全部可解析，无悬挂 label**

- `ch:sixth`（自定，行 4）被 `app.tex:652` 引用 ✅；`tab:species16`、`eq:{hccrosswalk,ps404,namodel,naobj,nacounter,swlogcurv,ray,window,dirichletratio,physwindow,dcom,barE,l2,bridge,kyfan,ratio,eta,nasun}` 均自定自用 ✅；
- `ch:immune` → `ch08-immune.tex:3` ✅；`app:erratum-radius` → `app.tex:611` ✅；
- 盒/宏：`claimbox`、`claimhead`、`qbox`、`\LT/\LC/\LM/\LN/\LE` 均在 `main.tex:37–56` 定义 ✅（`ch16` 已 `\input`）。
- **未被任何章节 `\ref` 的 label**：`eq:dcom`、`eq:kyfan`、`eq:naobj`、`eq:nasun`、`eq:nacounter`、`eq:physwindow`、`eq:ray`、`eq:swlogcurv`（共 8 个只定义不使用——不是错误，但 `eq:nacounter`/`eq:nasun` 是承重式，正文却是纯文字引用，建议至少 `eq:nasun` 加 `\eqref`）。

**文字/记法问题（非数字）**

- R8：行 472 `\partial_v I_3=3\operatorname{tr}\bigl(T_v(\xi\perp T)^2\bigr)` —— 源为 `(ξ_\perp\cdot T)²`（`LEAD-NOTE-NA-SU3.md:16`），`\perp` 在此被排成**关系符**而非下标，读者会读成 “ξ ⊥ T”。应为 `(\xi_\perp\!\cdot T)^2`。
- R9：行 409–410 引号句 *“vacuous at the only point where it was evaluated”* —— 一手源**无逐字对应**（`ROUND1-VERDICTS.md:152` 作“`(NA)` 的字面形式在自然条件化点上是空的”、`:155` 作“作为 `a₀=0` 的逐点陈述 vacuous”）。属**转述加引号**，建议去引号或改成逐字。
- R10：行 924 *“from thirty-three root families … to forty”* —— `SIXTH-SPECIES-CHARTER.md:§4` 的旧口径确为“33 族”，但 `ch08-immune.tex:92` **现表已是 40 行**（33+7）；两章并读会让读者以为 ch08 只有 33 族。建议写“ch08 的 40 族 = 本章新增的 7 族 + 原有 33 族”。（`SIXTH-SPECIES-SYNTHESIS.md §3` 同用 33→40，故这是**全局口径**，不是 ch16 独创的错。）

---

## §6 建议的最小修补清单（15 条；只动 ch16 文本，不动事实源）

| # | 哪一句（行号） | 改成什么 |
|---|---|---|
| 1 | 行 582–583 `m_{0^{++*}}/m_{0^{++}}\approx1.56`（“$1.56$ is the conservative end”） | `the corpus records the lattice range $1.6$--$1.8$ in \code{research/v3/Z1-HOLO.md} (§2), and the operative lattice ratio $\rho_{\mathrm{lat}}\approx1.54$--$1.56$ in \code{research/unconditional/RATIO-RHO.md} §5; taking the conservative end $1.6$` ⇒ `\rho^2_{\mathrm{phys}}\approx1.6^2=2.56>2.25`。同步改行 585（`eq:physwindow`）与行 901–902 台账里的 `2.43`。 |
| 2 | 行 552–563 `-u''+cu=\lambda u` / `(4π²+c)/(π²+c)` / “confirms … $c=9.7$: exact 175.5757”；**及行 859–861 “given (ii), … yields $\rho^2\le C^2$”** | (a) 写成**源的真实约定**：`q=-c` ⇒ `μ_k=k^2π^2-c`，`ρ²(c)=(4π²-c)/(π²-c)`（`0≤c<π²`），`ρ²≤4\iff c≥0`（sharp），`c=9.7` ⇒ exact `175.5757`、FD `175.5762`（`rb_ratio_bound.log` PART D）；或保留 `+c` 版但删掉 175.5757 半句。(b) **补 DT 限制**：`e_1\ge\pi^2/L_s^2` **不蕴含** `ρ²≤4`（`DT-DIRICHLET.md` DT-2(b)/DT-3：见证 B `ρ²=27.1846$ @ $e_1=18.0963\ge\pi^2`、见证 A `ρ²=4+M/π^2=4.5/9.0`、rank-$r$ `4,9,16,24,\dots` 无界）；合法形式是**双侧** `c_-L_0\le L\le c_+L_0 ⇒ C=2\sqrt{c_+/c_-}`。(c) 行 859–861 改成 *“given (ii) **and** `T_2(V)\le 5e_1` (i.e. `t\le(1+C^2)\pi^2/L_s^2`) …”*——单靠 (ii) 不给窗口。 |
| 3 | 行 713 `eq:l2` 之后 | 补 **QI 精确恒等式**（`QI-BRACKET.md` QI-1）：`w_1=\dfrac{s_2-\bar E}{s_2-s_1}+\dfrac{\Delta}{s_2-s_1}`，`\Delta:=\int_{[s_2,\infty)}(s-s_2)\,d\nu\ge0`；并加一句“`slack=\Delta/\varepsilon` 是 UV 敏感的超额一阶矩 ⇒ 桥是**少能级证书**，不是机制”；同处可加 `\vartheta\ge w_0\iff\Delta\le(w_1-w_0)\varepsilon`。 |
| 4 | 行 934–936 *“the SU(N≥3) version of (NA) was not tested, because no solvable model … was built”* | **撤回该句**：`\code{NA2-NONDEGENERATE.md} §0` 已在一 plaquette $SU(3)$ 上测出 $X(a_0=0)\neq0$（量级 $10^{-4}$–$1.7\times10^{-2}$）且**随边界 holonomy 变号** ⇒ $(NA)$ 在自然条件化点上被更强地否证；仍开放的是 **4D 的 $\partial_vI_3$ 源**，不是“模型没造”。 |
| 5 | 行 526 *“violated by every exact anchor large enough to test it”* | 改为 *“violated by five of the exact anchors below (all but Mehler $2$ and 2D $SU(3)$ $9/4$, which sits exactly on the bound)”*（`RB-RATIO-BOUND.md:321` 原文就是“5 个”）。 |
| 6 | 行 811–813 *“The same model kills …”* + `×5.6` | 加模型限定：*“the SUSY double well of \code{RB-RATIO-BOUND.md} §RB-7 ($H=-\partial_x^2+W$, $W=V'^2/4-V''/2$, FD $N=6000$) — not the $H=-\tfrac12d^2/dx^2+(g/2)(x^2-1)^2$ of the table above”*；`×5.6` → `×4.6`（`20.05/4.34=4.62`）。 |
| 7 | 行 180 *“two of 320 self-loops tie”* | 改为 *“all $n=320$ self-loops tie at $\rho_\infty$ (the critical graph is the whole diagonal); the best 2-cycle mean is strictly lower, $2.487962<2.5$, with deficit $\beta(1-\cos 2\pi/n)$”*（`tr_maxplus.log:17` 与 `:27` 自相矛盾，取 `:27`+`AUDIT-R1.md:139`）。 |
| 8 | 行 199–200 `\LT{} … \LM{} … \LC` | 删 `\LM`，改成两句两标：max-plus 谱句 `\LT`；“compactness is a property of the measure” → `\LC`（读法）。 |
| 9 | 行 307 `… \LE` ＋ 行 307–308、行 358 的拟合律 | 行 307 句末 `\LE` 改 `\LT`（零方向是初等的，与行 314 自述一致）；**行 307–308 与行 358 两处 `d_N\sim e^{-0.4595N}` 各补 `\LE`**（拟合，`N≤14`），台账行 896 拆标（lead 指令 H2）。 |
| 10 | 行 625–629 `\LT` | 拆：`(TCI)\iff\lambda_{\min}(A_{\mathbb K}\mid\mathcal H_-)\ge cF+\mu` 留 `\LT`；“the **only self-consistent** completion + seven-element enumeration” 改 `\LC`（`Q=3\lambda_4\varphi_0^2=3cF` 在源里是 `[MINE]`/`(D3)`）。 |
| 11 | 行 440 *“nine orders of magnitude”* | `ten to twelve orders`（`10^{-2}/1.5\times10^{-13}\approx7\times10^{10}$，$10^{-1}\sim7\times10^{11}$）。 |
| 12 | 行 190 *“with width/Δ=10.5”* | `width/\Delta=10.465 at n_z=4096 (4.255 at n_z=64)`（`AUDIT-R1.md:136–137`）。 |
| 13 | 行 903 台账 “Radius” `\LT` | 拆：`1/48`（SZZ 被引定理）`\LT`；*“`\lambda_{\mathrm{tH}}<1/24` absent from the source”* → `\LC`（书目核查，非定理）。 |
| 14 | 行 472 `3\operatorname{tr}(T_v(\xi\perp T)^2)` | `3\operatorname{tr}\bigl(T_v\,(\xi_\perp\!\cdot T)^2\bigr)`（源 `LEAD-NOTE-NA-SU3.md:16`）。 |
| 15 | 行 545（无标）与行 646（无标） | 行 545 补 `\LT`（五个精确锚点违反 1.5、两个饱和 4）；行 646 补 `\LE`（`a_*\approx9`--`14` 出自 `tci_comp.log [D]`）。 |

---

## §7 附：本轮**逐字核销**的承重数字（供复核，非缺陷）

`8.5e-16`、`20.72→36.84`、`(0,2,6,12)`、`2.5e-11`、`0.0796`、`1/6,1/4,1/3`、`4.3e-9`（§1 HC）；
`0.7871`、`0.001944`、`404.8`、`2.47e-3`、`12.18`、`54.6`、`2.9e215`、`5.2e4`、`0.8871`、`V=8/16/256/65536`（§2 PS）；
`0.6444`、`3.013`、`0.94`、`1.718→3.443`、`2.415`、`4.85→3.89`、`0.267884`、`8.7e-4`、`0.9327`、`0.0697`、`9.2×`、`750×`（§3 TR）；
`1891.27`、`−1`、`2.7e-7`、`−1.0000/−1.7071/−3.6180/−15.5495/−36.5311`、`84.89`、`0.5034`（§4 SW/§5 LX）；
`0.25/0.4/301`、`5.8219e-1/3.4479e-1/5.3758e-2/1.3634e-2`、`5.601e-1/5.717e-1`、`1.0728e-2/1.1192e-2`、`−0.4595`、`4.879e-1/2.183e-4`、射线表 21 个值 + `0.183h`（§6 DM）；
`1.5e-13`、`1e-8`、`−0.0158229063`、`+0.0316458126`、`0.20511530`、`0.01105929`、`4.30661`、`0.13599`、`−0.19302`、`0.0691`、`0.3981`、符号表全表（§7 NA）；
`1.5/2.25/2/2.6667/9/4/3/4`、`1.25/3`、`1/48/1/24/120/48/24/16`、`215/3.999999/1.00000`（§常量）；
`0.1667→0.7280`、`0.0109`、`0.52`、`4c²`、`1.00000/0.97980/0.91652/0.60000`、`6e4`、`dim12`、`0.166680/0.166667`、`1.70/0.851`、`−0.1796/−0.6745`、`0.070/2.72/0.76`、`a_*9–14`、七元 `\mathbb K`（§TCI）；
`4.3e-5/1.3e-15`、`0.9935…0.9778` 七对过报、`0.023439/0.052609/0.076048305`、`2.529/5.420`、双井表 42 个值、`180`、`4.7e-3`、`0.5308/0.5467/0.5629`、`0.5288/0.5475/0.5633`、`206/1170/3894`（§bridge）。
→ 以上在 13 份一手文件或 `w7/*.log` 中**逐字命中且用法一致**（`OK=140` 的构成）；§1/§2/§3 已判缺陷的数字**不在**此列。

## §8 追加高危点专项对账（lead 指令；一手源之间的已知矛盾，逐条查 ch16 是否继承）

> 对账时点：`SIXTH-SPECIES-SYNTHESIS.md` 仍为 23:43:12 版，**尚无 `## §-1 对账记录`**（`Select-String '^## '` 只出 §0–§7）。故下列对账基于 `SYNTHESIS §1/§3` 现文；§-1 落盘后第 1 项应按其新口径复核一次。

| # | 一手源之间的矛盾 | ch16 的位置与原句 | 判 |
|---|---|---|---|
| **H1** | **TR 混模型**：`SYNTHESIS §1` TR 行写“`Δ_trop≡0` … vs 精确 `Δ=0.267884`（n=64 与 n=4096 同值）**；夹逼 `width/Δ=3.013`（L=1）→ `2.415`（L=6）**”——把 **flat** 的 `0.267884` 与 **P4 锚**的 `3.013` 并列。正确配对：`Δ=0.267884 ↔ width/Δ=4.255（flat, n=64）` 或 `10.465（flat, n=4096）`；`Δ=0.6444 ↔ 3.013（P4 / tr_maxplus）` | ch16 **分两处写，未混配**：行 181–187 先给 `Δ=0.6444`（periodic anchor），再给 `width/Δ=3.013`（β=2.5）、`2.415`（L=6）、`s_L=4.85→3.89`；行 189–190 另起一句给 **flat** 模型 `Δ=0.267884`，配 `width/Δ=10.5`（= `tr_flat`/`AUDIT-R1.md:137` 的 `10.465`，n=4096） | **ch16 未犯** ⇒ **不记 SHIFTED**。唯一残留是 `10.5` **未标 n**（`n=64` 给 `4.255`）——已记 §2-S4（低）。**建议顺手加 `(n_z=4096)`，以免下一次回填再把这两个 Δ 配错** |
| **H2** | **指数律标签**：`d_N~e^{−0.4595N}` 是 `N≤14` 的**窗口拟合**，只能 `\LE`；`\LT` 只覆盖 `d_N>0` 与 `min≡0` 的双向性 | 逐处：①行 307–308 *“A log-slope fit over $N=5,\ldots,14$ gives $d_N\sim e^{-0.4595\,N}$”* → **无标**；②行 358 *“(moments $e^{-0.46N}$ **by fit**; …)”* → **无标**；③行 896 台账 → `\LT/\LE` **合标**。且行 307 句末唯一那个 `\LE` 挂在 **“the minimum is 0.0 at every $N$ and every window”**（应 `\LT`）上 | **FAIL**：作者自称“拟合降 `\LE`”，但**两处根本没标**、一处合标，且 `\LE`/`\LT` 在相邻两句上**互换了方向** ⇒ 记 §4-J（标签问题 +1，不计数值口径） |
| **H3** | **`floor ⇒ ratio`**：`DT-DIRICHLET.md` DT-2(b)/DT-3 用两个见证否证“单侧 floor ⇒ `ρ²≤4`”：见证 B `ρ²=27.1846 @ e₁=18.0963≥π²`（局部变号，`dt_dirichlet.log:139` PART C3）、见证 A `ρ²=4+M/π²=4.5/9.0`、rank-r `4,9,16,24,…`（无界）；合法形式是**双侧** `c₋L₀≤L≤c₊L₀ ⇒ C=2√(c₊/c₋)`。这与 `RB-RATIO-BOUND.md` RB-6(iv) 的“`q≥0 ⟺ e₁≥π²/L²`”读法冲突（DT 称 RB 原句“容易被读成 floor 就够”） | ch16 **引了 RB 的等价**（行 562：*“$c\ge0$ means … i.e.\ $e_1\ge\pi^2/L_s^2$”*），并在行 859–861 下结论 *“given (ii), the rank-two trace (ratio) yields $\rho^2\le C^2$ with $C=2$”*；**全章无 `DT`、无 `27.18`、无 `18.0963`、无“floor ⇏ ratio”**（grep `27\.1\|18\.09\|DT-` 零命中） | **FAIL（OMISSION，重）**：按 lead 口径必须同时给出 DT 的限制，否则 `(ii)⇒ρ²≤C²` 这条推导被印成合法而实际不成立。见 §6 补丁 2（并修 patch 2 的措辞） |
| **H4** | **SW 两句必须分开**：固定 `x` 时 `(logρ)″` 与 β 无关；模态处随 β 变（`−1.0 → −36.5`） | ch16 行 219–230 **主动照做**：*“The audit requires two separate sentences here, and both are true. **At fixed $x$**, the right-hand side … does not involve $\beta$ at all … a finite-difference check at $x=0.2,1.0,1.8$ for $\beta=0,2.5,4$ agrees to $\le2.7\times10^{-7}$. **At the mode**, the value does move with $\beta$: … $-1.0000$ ($\beta=0$), $-1.7071$, $-3.6180$, $-15.5495$, $-36.5311$ ($4.0$) \LE.”*，并保留 `ρ″` 不必处处负的 caveat（`β=10, x=0.5: 84.89>0`） | **ch16 未犯** ⇒ **不记 OMISSION**。两句齐、分开、数值全对（`AUDIT-R1.md:217–220`） |

**H 组净判**：4 条里 **2 条 ch16 未继承**（H1、H4），**1 条继承为标签换向**（H2），**1 条继承为缺失限制**（H3，重）。§0 计数据此更新：`OMISSION 2→3`、`标签 9→10`、`SHIFTED 仍 4`。

---

*ch16-audit · 2026-11-24 · 未改任何 `.tex` · 未 git commit · 全部判定可由 `grep` 复现*
