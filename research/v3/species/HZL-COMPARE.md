# HZL-COMPARE — 物理尺度的正比较算子：定义、Cheeger 关系、三路判定

> 轮 HZL · 子代理 `hzl-compare` · `D:\Dev\math`。只新建本文件；既有文件原样；未 `commit`；零 MC；未动 `docs/position-current.html`（`AGENTS §6` 交 `lead`）。
> 靶（承 `HZA-LOWER §0.5/§2 末行`）：18 条工具里唯一未覆盖的形状 = **物理尺度的正比较算子**；最便宜检验 = Cheeger `h(L)` 的 `L`-标度（零 MC）。
> 窗口：`β=2.5`、`SU(2)`、Wilson、固定 `a`、`L→∞`；`ξ/a=1.0±0.3` 是 `A4-SCALE-SEP:29` 的测量输入，不是定理。标签只 `[T]`／`[NG]`。

## §0 判词

| # | 断言 | 判 | 依据 | 下一手 |
|---|---|---|---|---|
| 0.1 | **定义必须排除电项比较**：`K_a=(g²/2a)Σ_l Ê_l²` 的谱已知，`γ_K=3g²/(8a)`（SU(2) 基本表示）；但 `K_a` 与 `H_a` 的比较有两种写法，都不满足物理尺度——带磁项常数 `C∝N_p/a` 只给 `Δ_a≥γ_K−C`，`L→∞` 变负；把常数吸收进 c-数则只剩 `O(g²/a)` 的 UV Casimir 标度 | `[T]` | `HZA-LOWER:27,44`；`GM1-ALGEBRA:61`；`HX6-DBETAH:27`；`GO2-SOFTMODE:162` | 把 (P3) 写成“非 `K_a−C_{N_p}`、非 UV 标度”，再检验任何候选 |
| 0.2 | **Cheeger 关系是唯一精确的一行**：对任何满足 (P1)-(P2) 的 `B_a`，`Δ(B_a)≥h(B_a)²/4`；若 `B_a` 是粗模边际 `ν_L` 的生成元，则 `h(B_a)=h(L)`，`h(L)≥h_0>0 ⇒ Δ≥h_0²/4`；`h(L)∼c/L ⇒` Cheeger 证书只有 `O(L^{-2})`，无 `L`-一致地板。方向是**充分**证书，不是等价（细刺域 `λ₁>0` 而 `h→0`） | `[T]` | Cheeger 1970（DOI `10.1515/9781400869312-013`，书目级；`O18-CHEEGER:48`）；`O18-CHEEGER:44,54`；`CHEEGER-BLIND:121`；`HZA-LOWER:46` | 在既有 `L=4,6,8` 配置上量 `h(L)`（零 MC），看 `h(L)·L` 是否平台 |
| 0.3 | **三路全 `[NG]`**：① 由 `H_a` 显式形 → 只有电项 `K_a`（UV/广延）与磁项 `V_a`（正但平坦、广延）；② 由 Savvidy 背景 → `Hess=−D²+2adF_B` 有 `λ_min=−g|B|<0`，任何 `B_a≤Hess` 必负；③ 由 `F²` 逐点正性 → 只给 `V_a≥0` 的符号，零集含真空/平坦方向，无正隙 | `[NG]` | `GM1-ALGEBRA:61`；`H1E-SAVVIDY:104–113`；`O18-CHEEGER:130`；`H1A-EXACT:102` | 三条都不再试；把“造 `B_a`”改成“量 `h(L)`” |
| 0.4 | **净判**：三路都造不出物理尺度 `B_a`；唯一非空的形状是 Cheeger 生成元，而 `h(ν)≍1/√χ`、`h>0` 一致 ⟺ H1 ⇒ 这个“未覆盖形状”在精确化后**就是 H1 的瓶颈形状**，不是第六种机制 | `[NG]`（射程：三路） | `CHEEGER-BUILD:3,101`；`B3-ISOPERIMETRIC:32`；`HZA-LOWER:15,66` | 若 `h(L)·L` 有平台 → 交 `λ₃`/H*；若 `h(L)∝1/L` → 关掉 Cheeger 这一格，回 `λ₃` 与 `f₁²≥r₀` |

**一句话**：正比较算子不是“再找一个正项”，而是“序 + 正隙 + 物理尺度”三件同时成立；电项有隙但无尺度，磁项有尺度但无隙，Savvidy 有尺度但有负模。`[NG]` 下一手：见 0.2/0.4 行末。

## §1 定义

### 1.1 对象与归一

`H_a := −(1/a)\log T_a` 于 `1^\perp`，`H_a\Omega=0`（真空能已扣），`Δ_a := \inf\sigma(H_a)\cap(0,\infty)`（`GN3-SPECTRAL:87–88`；`B01-RP:11–12`）。所有比较都在 `\Omega^\perp` 的二次型意义下做；`H_a` 的 c-数平移不改变 `Δ_a`，但改变“谁是正算子”，故 (P1) 必须对**已扣真空能**的 `H_a` 写。

### 1.2 定义（物理尺度的正比较算子）

称 `(B_a,c)` 是**物理尺度的正比较算子**，当且仅当：

- **(P1) 序**：`B_a` 自伴，`Q(B_a)\supseteq Q(H_a)`，且 `B_a\le H_a`（`\Omega^\perp` 上的二次型）。于是 Weyl 单调性给 `\lambda_k(H_a)\ge\lambda_k(B_a)`。
- **(P2) 正隙**：`B_a\ge0`，`\Delta(B_a):=\inf\sigma(B_a)\cap(0,\infty)\ge c>0`，`c` 与 `L` 无关。
- **(P3) 物理尺度**：`B_a` 的二次型只依赖粗模 `k\lesssim1/\xi`（不是单链 Casimir 的和），且 `\Delta(B_a)` 在 `a\to0` 时保持 `O(1)`（物理单位下的质量），不随 `g²(a)\to0` 消失，也不是 `O(g²/a)` 的 UV 标度。

**一句读法**：`K_a` 满足 (P1)(P2) 但违反 (P3)；`V_a` 满足 (P1) 但违反 (P2)（平坦方向/广延）；Savvidy Hessian 违反 (P2)（负模），且是 UV 二次型、不满足 (P3)。三者都不是 `B_a`。

### 1.3 与 Cheeger `h(L)` 的一行式

对任意 `B_a` 满足 (P1)-(P2)，取 `\mu_{B_a}` 为其 Gibbs 测度、`h(B_a)` 为其 Cheeger 常数：

\[
\Delta(B_a)\ \ge\ \frac{h(B_a)^2}{4},\qquad
h(B_a)=\inf_A\frac{|\partial A|_{\mu_{B_a}}}{\min(\mu_{B_a}(A),\mu_{B_a}(A^c))}.
\]

若 `B_a` 是粗模边际 `\nu_L` 的生成元，则 `h(B_a)=h(L)`，于是

\[
h(L)\ge h_0>0\ \Longrightarrow\ \Delta(B_a)\ge h_0^2/4;\qquad
h(L)\sim c/L\ \Longrightarrow\ \text{Cheeger 证书}=O(L^{-2})\text{，无 }L\text{-一致地板}.
\]

**方向**：这是**充分**证书，不是等价；`h(L)` 小不排除别的 `B_a` 有大隙（细刺域：`λ₁>0` 而 `h→0`，`O18-CHEEGER:33`）。但 Cheeger 是唯一已知的、把“序 + 正性”同时写出来的通道，故它是最便宜的**可证伪**检验。库内 `h(\nu)\asymp1/\sqrt{\chi}`，`h>0` 一致 ⟺ H1（`CHEEGER-BUILD:3,101`；`CHEEGER-BLIND:121`）。

## §2 三路判定

### 2.1 由 `H_a` 的显式形：`[NG]`

`H_a=(g²/2a)\sum_l \hat C_2(l)+(1/(g²a))\sum_p(N-\mathrm{Re\,tr}\,U_p)`（`GM1-ALGEBRA:61`；`HX6-DBETAH:27`）。显式形里只有两块正项：

- **电项 `K_a`**：谱已知，`\gamma_K=3g²/(8a)`（SU(2) 基本表示）。`K_a` 与 `H_a` 的比较绕不开磁项常数：扣常数后 `K_a-C≤H_a` 只给 `\Delta_a\ge\gamma_K-C`（`C\propto N_p/a`，`L\to\infty` 变负，`HZA-LOWER:44`），第一种写法在 `L\to\infty` 时连 (P2) 都保不住；不扣常数则 `K_a` 与 `H_a` 差一个 c-数，`\gamma_K` 只是 `O(g²/a)` 的 UV 标度，违反 (P3)。
- **磁项 `V_a`**：逐点 `\ge0`，但在真空/环面平坦方向上为零，且广延（`\propto N_p`）⇒ 无正隙，违反 (P2)。

显式形没有第三块正项；任何“把 `V_a` 下界成 `-C` 再与 `K_a` 相加”的做法都回到广延减除。**下一手**：不再从显式形找 `B_a`；把 `H_a` 的显式形只用于**验证**候选 `B_a` 的 (P1)。

### 2.2 由 Savvidy 背景：`[NG]`

Savvidy 二次型 `H_{\rm Sav}=-D²+2adF_B`（`H1E-SAVVIDY:104`；`BianchiLemma.md:212`；`REPORT.md:465`）。其最低模为

\[
\lambda_{\min}(H_{\rm Sav})=-g|B|<0
\]

（横向磁负模；规范商后仍在，W4 实测负模比 0.071，`P^THP` min `≈-0.3`，`O18-CHEEGER:130`）。若 `B_a\le H_{\rm Sav}`，则 `\lambda_{\min}(B_a)\le -g|B|<0` ⇒ 违反 (P2)。平移 `+g|B|` 使 `\lambda_{\min}=0`，但 `B_a=H_{\rm Sav}+g|B|\ge H_{\rm Sav}`，违反 (P1)；平移 `-g|B|` 则 `\lambda_{\min}=-2g|B|<0`，仍违反 (P2)。**下一手**：Savvidy 背景只用于**否定**候选（找负模），不用于造正比较。

### 2.3 由 `F²` 的逐点正性：`[NG]`

`\frac14\mathrm{tr}F²=\frac18\sum_a(F^a_{\mu\nu})²\ge0` 逐点（`\mathrm{tr}T^aT^b=\frac12\delta^{ab}`），故 `N-\mathrm{Re\,tr}\,U_p\ge0`，`V_a\ge0`。但这是**符号信息**，不是**谱信息**：

- `V_a` 的零集包含真空/平坦方向（`U_p=1`），故 `\inf\sigma(V_a)=0`，无正隙；
- `V_a` 是广延的（`\propto N_p`），其下界 0 不给 `L`-一致的正数；
- 逐点正性不能产生粗模 `k\lesssim1/\xi` 的算子；它只能复现 2.1 的电项比较或平凡 `H_a\ge0`。

**下一手**：把 `F²` 正性只当**符号输入**（例如 `V_a\ge0` 用于 (P1) 的验证），不当谱输入。

## §3 台账（`AGENTS §1` 四要素）

| # | 工具 | 状态 | 实际检索式 / 读入 | 命中 | 失败标记 |
|---|---|---|---|---|---|
| 1 | `pwsh` | **OK**（exit 0） | `Get-ChildItem research\v3`；`research\v3\species` 清点 | 2 目录 | 0 |
| 2 | `read` | **OK** | `HZA-LOWER.md`(90)·`O18-CHEEGER.md`(263)·`H1A-EXACT.md`(234)·`H1E-SAVVIDY.md`(186)·`B01-RP.md:1–60`·`RUN3.md:1–80` | 6 份 | 0 |
| 3 | `grep` | **OK** | `HZA-LOWER\|HZL\|正比较` @ 全仓 | **9** | 0 |
| 4 | `grep` | **OK** | `H_a\s*=\|H_a:=\|K_a\s*=` @ `research/v3/species/*.md` | **250 显示 / 549 全表**（截断，落 spill） | 0 |
| 5 | `grep` | **OK** | `物理尺度\|physical scale\|正比较\|positive comparison` @ 全仓 `*.md` | **116** | 0 |
| 6 | `grep` | **OK** | `h\(L\)\|h_L\|Cheeger.*L\|L-标度` @ `research/**/*.md` | **136** | 0 |
| 7 | `grep` | **OK** | `Kogut\|Kogut–Sitter\|electric term\|电项\|Ĉ₂\|Ê_l` @ `research/v3/species/*.md` | **78** | 0 |
| 8 | `glob` | **OK** | `**/*SAVVIDY*.md` | **17** | 0 |
| 9 | `web_*` | **0**（未调用） | 本文件无新外部断言；Cheeger 1970 按书目级转引 | — | — |
| 10 | MC / Lean / `git` / `skill` | **0**（未调用） | 零 MC；catalog 无匹配技能 | — | — |

**诚实边界（`AGENTS §1`、`§2.1`）**

1. 本文件**未新算**任何数；`γ_K=3g²/(8a)` 与 `λ_min=-g|B|` 是 SU(2) 基本表示/单模的标准值，承重判词落在仓内 `file:line`。
2. 外部标准件只引 Cheeger 1970（DOI `10.1515/9781400869312-013`，书目级，`O18-CHEEGER:48`）；Weyl 单调性与 Savvidy 二次型按标准形式使用，**未逐字核原文**，故不列具体文献条目。
3. 窗口：固定 `β=2.5`、`SU(2)`、Wilson、固定 `a`、`L→∞`；`SU(3)`/`β≈7.6` 不混入。
4. 本文件判的是“三路能否造出物理尺度 `B_a`”，不是“隙不存在”；`[NG]` 的射程只到这三路 + Cheeger 生成元。
5. `h(L)` 的 `L`-标度是**下一步的测量**，本文件未做；本文件只给它的判据与方向。

*HZL-COMPARE · 只新建本文件；既有文件原样；未 `commit`；零 MC。*
