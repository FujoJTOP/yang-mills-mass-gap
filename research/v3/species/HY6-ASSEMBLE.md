# HY6-ASSEMBLE — 最终证明骨架（Lean 依赖序 · 逐定理 · 每条引机检）

> **轮** HY6 · 子代理 `hy6-assemble` · `D:\Dev\math`。只新建本文件；既有文件保持原样；未 `git commit`；未编译 Lean；零 MC。
> **范围（题给）**：三格 C（`ι_a`/L² 缝）→ A（`sup_a[ρ^ren 𝒢^ren]<∞`）→ B（`s₁ ≥ c·M₀/χ`）+ 五条已机检 + 六条本批（HV9/HV7/HV8/HU1/HU2/HU3）。HX7 的其余记账不重排。
> **标签**：只交 `[T]`/`[NG]`；`[T]`=已证/恒等/闭式算术/量词演算/机检蕴含，`[NG]`=具名缺步或具名 no-go。不用百分数；每条结论附「下一手」。
> **五条机检（只此五条）**：`gap_ge_liminf`（`GL2bPrime.lean:128`）· `susceptibility_squeeze`（`Reduction.lean:101`）· `spectral_decay`（`GL9Decay.lean:85`）· `exists_subseq_tendsto_pi` / `gl3_imp_with_minlos`（`GL3Limits.lean:58,82`）。
## §0 判词

| # | 断言 | 判 | 一行依据 | 下一手 |
|---|---|---|---|---|
| **0.1** | **骨架完备**：五机检 + 13 条 `[T]` 派生 + 7 条 `[NG]`（6 义务 + 1 反例）+ 3 条装配，全部有位，每条边有类型，三个开输入具名 | **`[T]`**（记账） | §1–§3 | 先落 D0 定义文件 |
| **0.2** | **★ Lean 数据序 ≠ 题给箭头**：`A→T4→μ→C→桥→gap`；题给 `C→A→B` 是 HX7 的**成本序**（C 一轮、A 多轮、B 新机器） | **`[T]`** | `HX7-PATH:91`；§2 | 两序分写，不混用 |
| **0.3** | **★ HV7 的 `s₁ ≥ 4arccoth(ρ_β)²` 不成立**：`φ_β=coth(√s/2)` 减 ⇒ `⟨φ_β⟩ ≤ φ_β(s₁)` ⇒ 正确方向 `≤`（cap）；其日志 [10] 两原子 `1.7385<2.164` 即反例 | **`[NG]`**（`≥` 作地板）/ **`[T]`**（`≤` cap） | `HV7:0.5,3.2(B4)` 与其 `:72` 的 [10] 数冲突 | 改 B4 方向；只登记 cap/位置 |
| **0.4** | **★ E1 只覆盖有界+范数收敛**；物理 `H_a` 无界、范数收敛不可能 ⇒ 物理桥是强预解 4 行版（`[T]` 非机检） | **`[T]`**（引用审计） | `GN2:52–54`；`GR1:47`；`GP9:131` | 桥分 T1/T1′；C 内含 SPEC-RES |
| **0.5** | **HV9 撤 β-积分路**（预算为零）· **HV8 地板改 `(1−ε)c_L`** · **HU3 删 M3**：三条都是删前提/删路，不是新估计 | **`[T]`** | `HV9:0.3`；`HV8:0.2`；`HU3:0.4` | 按此重写 B 的签名 |
| **0.6** | **HU1 欠账 = ④ 的 4D 面**（`M^4D`），③ 在切片类上免费；**HU2 承重 = 条件方差**（平均=Fisher `[T]`，点值无恒等 `[NG]`） | **`[T]`**/`[NG]` | `HU1:0.2–0.3`；`HU2:0.3,0.5` | L3 拆两面、M2 拆两截 |
| **0.7** | **Clay `[NG]`**：C/A/B 三输入开；装配 `[T]` 只是蕴含 | **`[NG]`** | §1 T16 | 攻 C 的定义层 |

**一句话。** 骨架是真骨架：五机检各自有位，六条本批各删一条路或改一条签名，三格各具名，装配是 `[T]` 蕴含。两处 ★ 勘误（HV7 方向、E1 类型）对任何格的判无影响，但挡掉两条误引。
## §1 定理清单（Lean 依赖序：导入 → 派生 → 义务 → 装配）

| # | 定理 | 一行陈述 | 判 | 引机检 | 下一手 |
|---|---|---|---|---|---|
| E1 | `YMGlue.gap_ge_liminf` | 有界 `H_k`、`‖H_k−H‖→0` ⇒ `gap H Ω ≥ liminf_k gap (H_k) Ω` | `[T]` | — | 只作有界特例（T1′） |
| E2 | `susceptibility_squeeze` | `0<χ ∧ f₂≤χm ∧ χm≤M₀ ⇒ f₂/χ≤m≤M₀/χ` | `[T]` | — | 供 T2 |
| E3 | `GlueA.spectral_decay` | `TΩ=Ω`、`T` 保 `Ω^⊥`、`‖Tx‖≤λ‖x‖` ⇒ 真空扣除后 `λⁿ` 衰减 | `[T]` | — | 供 T5 |
| E4 | `exists_subseq_tendsto_pi` | 坐标一致有界 ⇒ 全坐标收敛子列 | `[T]` | — | 供 T4 |
| E5 | `gl3_imp_with_minlos` | E4+`hMinlos` ⇒ 子列+`HasRealization` | `[T]` | E4 | 供 T4 |
| **T1** | `hy6_bridge_strong` | `H_a→H` 强预解、自伴下有界 ⇒ `Δ_cont ≥ liminf_a Δ_a` | `[T]`（4 行，非机检） | —（E1 仅特例） | 与 C 的类型对齐 |
| **T1′** | `hy6_bridge_norm` | 有界+范数收敛特例，即 E1 逐字 | `[T]` | E1 | 物理族不可用，只作形式核 |
| **T2** | `hy6_factor_le_one` | `χs₁≤M₀ ⇒ factor:=s₁χ/M₀≤1` | `[T]` | E2 | 记 factor 的窗 |
| **T3** | `hy6_b_iff_factor` | `χ,M₀>0 ⇒ (s₁≥cM₀/χ ⟺ factor≥c)` | `[T]` | — | B 的两写法同义 |
| **T4** | `hy6_moments_subseq` | A 的界 ⇒ 对角子列 + μ（Minlos 显式假设） | `[T]` | E4,E5 | μ 交 N1 |
| **T5** | `hy6_cluster_from_gap` | `GapOnPerp` ⇒ 体积一致指数聚类（速率 λ） | `[T]` | E3 | 供 T15 |
| **T6** | `hy6_hv9_budget_zero` | `(∀wA≥w₀, F·wA≥w₀) ⇒ 1≤F`；而 Grönwall `F=e^{−∫C}<1`（`∫C>0`）⇒ 预算为零 | `[T]` | — | 停投 β-积分路 |
| **T7** | `hy6_hv7_rho_cap` | `0≤s₁ ∧ 1<ρ ∧ ρ≤coth(√s₁/2) ⇒ s₁≤4arccoth²ρ` | `[T]` | — | 只作 cap/位置读数 |
| **T8** | `hy6_hv7_rho_floor` | `s₁≥4arccoth²ρ` | `[NG]`（反例=HV7[10]） | — | 退役，不作 B 供给 |
| **T9** | `hy6_hv8_floor_gate` | 自由块 `c_L+π²/(11Nj)` ⇒ `∀δ>0` 逐 j 不可；地板只能 `(1−ε)c_L` | `[T]` | — | 改 N5 签名 |
| **T10** | `hy6_hu1_slice_merge` | (F.51) 切片界 ⇒ `(M^sl)`，`C_n` 与 `(Λ,a)` 无关 | `[T]`（文本） | — | ③ 记账为免费 |
| **T11** | `hy6_hu2_fisher_avg` | `∫Var(∂V \mid a₀)dρ_{a₀}=I(ρ)` | `[T]`（文本） | — | 只作平均，不作点值 |
| **T12** | `hy6_hu3_drop_dec` | `(★)` mode-wise；`(★dec)⇒Var<δ_W⇒(★)` 单向（194.7×） | `[T]`（文本） | — | M3 不入 B 的前置 |
| **T13** | `hy6_b_to_gap` | `factor≥c ⇒ (B)`（HZ7-B 格 1–5 逐格 `[T]`） | `[T]`（文本） | — | B 的两层同义 |
| **N1** | `hy6_iota_l2_seam` | C：`ι_a`、`Ha`、`MemLp`、`(Ha−z)⁻¹→(H−z)⁻¹` 强 | `[NG]` | — | 一轮定义层：判是否定义格 |
| **N2** | `hy6_l3_scalar_pair` | A：`∃C ∀a≤a₀, ρ^ren_a 𝒢^ren_a≤C` | `[NG]` | — | 先钉 `Z_a` 块和/块平均 |
| **N3** | `hy6_m4d_face` | A 的 ④ 面：`(M^4D)` 的 a-一致 4D 矩界 | `[NG]` | — | 拆出切片行后单攻 |
| **N4** | `hy6_b_floor` | B：`∃c>0 ∀L ∀β≥β_A, s₁≥cM₀/χ` | `[NG]` | —（词法用 E2） | 按 T9/T12 重写 |
| **N5** | `hy6_starsf_uv` | B 塔面：`Z_j^{-1}G″≥(1−ε)c_L`（j-一致，UV 壳） | `[NG]` | — | 五缺件 D-1..D-5 |
| **N6** | `hy6_condvar_pointwise` | B 核：`Z_j^{-1}[E∂²V−Var(∂V \mid a₀)]≥δ̂` 逐 U′ | `[NG]` | — | 只攻"交互修正保号" |
| **T14** | `hy6_gap_positive` | C∧B ⇒ `∃δ>0, Δ_cont≥δ` | `[T]`（装配） | —（T1；E1 特例） | 等 N1 |
| **T15** | `hy6_unique_limit` | C∧A∧B ⇒ 唯一极限测度 | `[T]`（装配） | E3,E4,E5 | 等 N1/N4 |
| **T16** | `hy6_clay_of_three` | C∧A∧B ⇒ Clay（OS 重构 + 两柱） | `[T]`（装配） | E2–E5 | 等三格 |
| — | `Clay` | 三格未证 ⇒ 不成立 | `[NG]` | — | C → A → B |
## §2 依赖序

**2.1 Lean 声明 DAG（不是时间序）**
```
D0 定义层（无判，一个文件）: iota, μ, Ha, Delta, DeltaCont, factor, rhoBeta
D1 导入: E1 E2 E3 E4 E5（唯一边 E5→E4）
D2 [T] 派生（互不依赖，全并行）: T1←∅ · T1′←E1 · T2←E2 · T3←∅ · T4←{E4,E5} · T5←E3 · T6–T13←∅
D3 [NG] 义务（声明互不依赖，签名可同时写）: N1 N2 N3 N4 N5 N6
D4 装配: T14←{T1,T13} · T15←{T4,T5,N4} · T16←{T14,T15}+OS
数据边: {N2,N3}→T4→μ→N1→T1→T14 ; N4→T13→T14 ; N4→T5→T15
```

**2.2 三格的真序**
- **声明序**：A、B 可在未构造 μ 时先写（纯测度语言）；C 的实例化要 A 的 μ；桥要 C。
- **攻击序**（HX7 成本）：C（定义层，一轮）< A（解析）< B（新机器+机检）；题给 `C→A→B` 读作成本序，Lean 数据序是 `A→C`。两序不冲突，但 brief 必须写清是哪一序。

**2.3 可并行 / 必须串行**

- 可并行：A 两面（N2,N3）· B 三面（N4,N5,N6）· 六条本批（T6–T13）· C 与 B 的签名 · E1–E4。
- 唯一串行链（关键路径）：`{N2,N3}→T4→N1→T1→T14`；汇合点 T16 要 `N2∧N3`（存在）+ `N4`（地板/唯一性）+ `N1`（对象）。

## §3 ★ 必须新写的 Lean 声明（名字 + 类型签名草案；仅签名）

**3.1 定义层（6 条，一个文件，无判）**
```lean
iota : ℝ → Obs → E                        -- ι_a[O] ∈ L²(μ)
Ha   : ℝ → E →L[ℝ] E                      -- H_a := −(1/a)·log T_a
Delta : ℝ → ℝ := fun a => gap (Ha a) Ω    DeltaCont : ℝ := gap H Ω
factor : ℝ := s1 * chi / M0               rhoBeta : ℝ := |dBetaA| / M0
```

**3.2 ★ 三格主签名（C / A / B）**
```lean
theorem iota_l2_seam (a : ℝ) (O : Obs) : MemLp (iota a O) 2 μ
theorem ha_strong_resolvent (z : ℂ) (hz : 0 < z.im) (ψ : E) : Tendsto (fun a => (Ha a - z • 1)⁻¹ ψ) (𝓝[>] 0) (𝓝 ((H - z • 1)⁻¹ ψ))
theorem l3_scalar_pair : ∃ C, ∀ a ∈ Set.Ioc 0 a0, rhoRen a * calG a ≤ C
theorem m4d_face : ∀ n, ∃ C, ∀ a ∈ Set.Ioc 0 a0, ∀ f, |S a n f| ≤ C * ∏ i, ‖f i‖
theorem b_floor : ∃ c > 0, ∀ L β, βA ≤ β → s1 L β ≥ c * M0 L β / chi L β
theorem starsf_uv : ∃ ε ∈ Set.Ioo 0 1, ∃ t0, ∀ a ≤ a0, ∀ j ≥ jUV, ∀ U' ∈ S t0,
    Zinv j * Gpp j U' ≥ (1 - ε) * cL
theorem condvar_pointwise : ∃ δ > 0, ∀ j ≥ j0, ∀ U', Zinv j * (Econd2 j U' - VarCond j U') ≥ δ
```

**3.3 ★ 桥与装配（`[T]` 蕴含；前件即三格）**

```lean
theorem bridge_strong (hres : StrongResolvent Ha H) : DeltaCont ≥ Filter.liminf (fun a => Delta a) (𝓝[>] 0)
theorem gap_positive (hres : StrongResolvent Ha H) (hB : ∀ a, Delta a ≥ δB) : DeltaCont ≥ δB
theorem clay_of_three (hC : IotaL2Seam) (hA : L3Bound) (hB : BFloor) : ClayStatement
```

**3.4 `[T]` 派生的写法**：T1′–T5 = E1–E5 的逐字实例化（签名 = E 的签名 + `hy6_` 前缀）；T6–T13 = §1 陈述逐字形式化，不另立类型。名字以 Mathlib 实际为准（`arccoth` 可能需 `Real.artanh` 组合），签名未过 elaboration，不据此宣称新机检。

## §4 台账

### 4.1 `AGENTS §1` 四要素

| # | 工具 | 状态 | 实际调用（逐字） | 命中 |
|---|---|---|---|---|
| 1 | `glob` | **OK** | `research/v3/species/*.md` | **764**（显示 100，全量落盘 `…677ebfc1b6a6-glob-results.txt`） |
| 2 | `grep` | **OK**（无 provider failure / ENOENT / 空响应） | ①`gap_ge_liminf\|susceptibility_squeeze\|spectral_decay\|gl3_imp_with_minlos\|exists_subseq_tendsto_pi` @ 全仓 ②`HV9\|HV7\|HV8\|HU1\|HU2\|HU3` @ `research/v3/**/*.md` ③`形式等级\|签名外\|有界算子.*范数收敛\|强预解` @ `species/*.md` | ①**157** ②**41** ③**103** |
| 3 | `read` | **OK**（全部带行号正文） | `HX7-PATH`(121)·`HV9:1–100`·`HV7:1–120`·`HV8`(123)·`HU1:1–70`·`HU2:1–75`·`HU3:1–80`·`GL2bPrime:40–139`·`Reduction:90–144`·`GL9Decay:1–107`·`GL3Limits:40–99` | 11 份全返回 |
| 4 | `pwsh` | **OK**（exit 0） | `Get-ChildItem research\v3 -Recurse -File` | 输出截断，全量落盘 `…8673b7f6846a-pwsh.txt` |
| 5 | 失败标记 | **TOOL-FAIL 0；空响应 0；PATH-MISS 0；provider failure 0** | — | — |

### 4.2 直读 vs 转引

- **直读（行号级）**：#3 的 11 份（六份本批文件各读到所引段；五份 Lean 只读签名区）。
- **只经命中行采信（未开全文）**：`GN2:52–54,309`·`GR1:47`·`GM1:52`·`HZ8:0.1–0.3,20`·`GP9:131`·`GT1B:43,130` ⇒ 记为转引级；**纯转引（经六份本批文件）**：`V1`（`Z_j` 归一）·`HZ7-B:15,57`（`(B)⟺`靶）·`HS-L3-9:15`·`HQ1-L3:191`·`HT1:133`·`S1:287`·`HZ5:168–172` ⇒ 记为转引级。

### 4.3 边界

1. 本文件不证任何格；`[T]` 全来自被引文件或两行算术。**唯一新算术 = HV7 方向**（`⟨coth(√s/2)⟩ ≤ coth(√s₁/2)`，用 `HV7` 自己的 [10] 数）。
2. 装配 `[T]` 是蕴含，不是 Clay；C/A/B 三输入开 ⇒ Clay `[NG]`。
3. E1 的类型限制是引用审计，E1 的机检状态不受影响；物理桥另立 T1（非机检）。
4. 签名草案未过 elaboration；不据此宣称任何新机检；既有文件保持原样、未 `git commit`、零 MC、未编译 Lean；无百分数、无 QED。

*HY6-ASSEMBLE · 五机检 + 13 `[T]` + 7 `[NG]` + 3 装配；★ 两处勘误：HV7 方向、E1 类型。*
