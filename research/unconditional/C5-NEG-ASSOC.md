# C5-NEG-ASSOC — negative association is not a species for SU(N)+Wilson: BBBL/strong-Rayleigh wants a real-stable/determinantal lattice, log-submodularity wants an attractive partial order, Chebyshev–TP₂ wants a common monotone scalar, and Stein–Helffer–Sjöstrand's constant *is* the spectral gap = H1; the β=2.5 `Cov<0` is the chart-local log-cut shadow of proven `Φ`-convexity, and a tight `|X|` is `Var(G|0)`=H1

> **判词 (verdict).** 最后一块是 `Cov(G,Φ′|0)≤0`（anti-FKG）：若立，`X=−2Cov≥0`，**免方差直穿 H1**。四种物种**无一存活**。(a) Borcea–Brändén–Liggett 的 strong-Rayleigh⇒NA 是 **Boolean 格/实稳定多项式/行列式点过程** 的性质；SU(N)+Wilson 是 **连续、受挫（`e^{βRe tr U_□}` 是 cos 之和）、非吸引** 的 Gibbs 测度，既非实稳定也非行列式，无承载格/核。(b) log-次模⇒FKG 要 **偏序+超模** 密度；`|F|` 与 `d(·,z)` 非乘积格坐标，`S_W` 非单调，且 FKG 给**正**关联（错号）。(c) Chebyshev/TP₂ 要 **共同单调标量**，而 `G`、`Φ′` 在 `2(N²−1)V` 维 fiber 上**无公共全序**；唯一可用 TP 在谱层 `≡ν≥0`、对自由场失效。(d) Stein/Helffer–Sjöstrand `Cov(f,g)=⟨∇f,(−L)^{−1}∇g⟩` 的常数 `C_p=1/λ`,粗模 `λ=Hess V_eff=G″`——即 `C_p=Var(a)=H1`，且只给**量级**、永不给**符号**。更致命：符号前提非普遍律——`β=5.7` 逐模 `Cov=−X/2` **变号**；`β=2.5` 的 `Cov≈−0.11` 只是 `Φ` 凸性（Prékopa）在 log-cut 的图卡投影。**裁决：物种不存在；任何"紧 `|X|`"经 Ward `G″Var(a)=1` 逐字即 `Var(G|0)≤C`=H1。五点：过 1、4（第 4 对 Stein 为 vacuous），败 2、3、5。**

## §1 四种方言 = 两种死法
(a)(b) **结构不适用**：NA / log-次模都是*格上*命题（Boolean/乘积格 + 偏序 + 吸引）。Wilson 受挫（`Re tr U_□=cos` 之和，非逐坐标超/次模）、连续、无偏序；`|F|` 与"近中心"共享全部 link，非独立坐标；强耦合团簇展开只能把 `|Cov|` 控成量级，收敛半径 `β≲1/48`（P3 已闭区），到不了 2.5。(c)(d) **循环**：Chebyshev/TP₂ 要一条把两变量拉到同序的标量（不存在，C4 §2.2）；Stein 常数 `=‖(−L)^{−1}‖=1/λ=`谱隙，粗模 `λ=G″`，于是 `C_p=Var(a)`——求 `Cov` 的界就是求谱隙 = H1。故四方言非四希望：两种 by-structure 不适用，两种 by-construction 循环。

## §2 "β=2.5 的 anti-FKG" 是影子，不是律
`X=−2Cov`，`β=5.7`（crossover）逐模 `Cov={+0.044,−0.079,+0.009,+0.030,−0.036,−0.024}`——**符号逐模翻转**，universal `Cov≤0` 为假。`β=2.5`（强耦合）测得 `Cov≈−0.11`（`X≈+0.22`，W7 Part B），恰是已证 **`Φ` 凸性→Prékopa→`W_geo≥N/12`** 压过（等分平衡的）Wilson 括号之**同一件事**：log-cut 处 `d→2π`、link→中心元 `z`、`U_□≈1`、**`F≈0`**，`Φ″` 在此暴增，与 `∝|F|` 的 Savvidy 驱动**反定域**（PROOF-SAVVIDY-WGEO-DOMINATION §2）——**单坐标奇性**，非测度级负相关。`Z_N` 心对称不强制符号（`U→zU` 只平移奇点，图卡 artifact）。把它当律 ≈ 把"Haar 凸性主导"当定理 = H1。

## §3 精确短差
```
需:      Cov_μ(G,Φ′|0) ≤ 0.11   (⟺ X ≥ −0.22，取 λ≈0.1)
C-S（测度无关界）: Cov ≤ √(VarG·VarΦ′) = √(2·0.02) = 0.20      (差 1.8×)
NA 若存在会给 Cov ≤ 0 < 0.11 —— 但:
  (a) BBBL 要 Boolean/实稳定 (无)   (b) FKG 错号 + 无偏序
  (c) Chebyshev 要共同标量 (无;TP≡ν≥0) (d) Stein C_p=1/λ=Var(a)=H1（循环,仅量级）
universal Cov≤0 假: β=5.7 逐模变号; 即便给符号, 余项 Δ_dyn=E[H|0]−Var(G|0)≥0 仍是 conditioning=H1.
紧 |X|<0.22 ⟺ Var(G|0)<W_geo²/(4VarΦ′)≈0.8–1.3 (vs 实测 ≈2) = H1.
```

## §4 五点裁决
1. **非微扰 — PASS.** 四种都是精确测度恒等式/不等式，无圈展开。
2. **控号下界 `m_eff²(k→0)≥λ` — FAIL.** NA 能定号却不存在；(c)(d) 只给量级/同义反复；β=5.7 逐模变号。
3. **非循环 — FAIL.** (c)(d) 常数即 `1/λ=H1`；(a)(b) 结构不适用；唯一非循环对象（真 NA 定理）缺席。
4. **非-SUSY + 败自由场 — PASS（对 (d) vacuous）.** `F≡0⇒G=Φ′=0⇒X=0`，不造伪正性；Stein 对自由无质量 `Var(a)=∞` 而**盲目**（SUB1-COV-BOUND），其过是空的。
5. **均匀 `≥λ` — FAIL.** 无 `λ`（=gap）；β=5.7 变号毁掉普遍性。

**Sharpest line.** 缺的"关联符号物种"不是缺一条不等式，而是 **H1 穿上了受挫测度的几何**：负关联是**吸引格与实稳定多项式**的性质，SU(N)+Wilson 都没有；两条**能在测度上写出来**的路（Chebyshev–TP、Stein–Helffer–Sjöstrand）所携带的算子范数*就是*谱隙——故它们以 H1 为代价只证**量级**、永证不出所需**符号**；而 β=2.5 的反相关是已证 Haar 凸性的图卡局部影子，非新律。**紧 `|X|`=Var(G|0)=H1；符号路既无不平等支撑、又无真普遍前提——本轮再次落在 H1 上，精确缺失的一步被点名：一条对连续、受挫、非吸引 SU(N) 测度的确定性负依赖定理，四种框架内均未知存在。**

*Refs.* C4-COV-LASTMILE（C-S `0.40>0.32`，逐模 `X` 变号，紧 `|X|`=H1）；B3-TYPICAL-COERCIVITY（`ρ` 不变，`Cov≤0.11` vs C-S `0.20`，1.8×）；B3-FUNCTIONAL-BE（`Γ₂=G″=H1`，conditioning）；OTHER-H1-POINCARE（conditioning=100% gap）；PROOF-SAVVIDY-WGEO-DOMINATION（log-cut 反定域、图卡局部、不可升级）；PROOF-POLYA-FREQUENCY（TP≡`ν≥0`，自由场失效）；SUB1-COV-BOUND（Stein 盲目）；ABC-NEW-INEQ（Ward `G″Var(a)=1`）；COMPUTE-GPP-SU3 / w7/su3_L4.log（逐模 `X`，`β=5.7` 变号；W7 Part B `+0.51`、`X≈+0.22`）；Fortuin–Kasteleyn–Ginibre 1971；Holley 1974；Borcea–Brändén–Liggett, *JAMS* **32** (2019) 205；Karlin, *Total Positivity* (1968)；Helffer–Sjöstrand 2000；Stein 1972；Prékopa 1971；Savvidy PLB 71 (1977) 133.
