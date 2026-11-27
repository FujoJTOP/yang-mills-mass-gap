# HV5-LEAN — `N1` / `N4` 的 Lean 声明草案（两件缺件 · `axiom`）

> **轮** HV5（`HY6-ASSEMBLE` 收尾批）· 子代理 `hv5-lean` · `D:\Dev\math`。只新建本文件；既有文件原样；**未编译 Lean**；未 `git commit`；零 MC。
> **靶（题给，不复推）**：`T16 ← {T14,T15} + OS`；`T14 ← {T1 ← N1, T13}`；`T15 ← {T4,T5,N4}`；`N1 = C`（可归约到 `L3`），`N4 = B`。
> **标签**：只交 `[T]`/`[NG]`。`[T]` = 类型/量词/恒等/机检蕴含可核；`[NG]` = 具名缺件。不用百分数；每条结论带「下一手」。

## §0 判词

| # | 断言 | 判 | 一行依据 | 下一手 |
|---|---|---|---|---|
| **0.1** | `N1` 修类型后可写成合法声明，**内容仍是缺件** | `[T]`（形式）/ `[NG]`（内容） | §1.1；`MemLp` 要 `Ω → H` | 先落 D0 的 `ι`/`Ha`/`Obs` |
| **0.2** | `N4` 的 `∃c>0` 形式完整，缺的只是 `c` 本身 | `[T]` / `[NG]` | §1.2；原稿缺 `χ>0`、`L≥2` | 按 `T9`/`T12` 重写后接 `T13` |
| **0.3** | **HY6 §3.2 原稿三处类型错**：`ι a O : E`（应 `Ω → H`）；`z : ℂ` 作用在 `ℝ`-线性 `Ha` 上；无界 `H` 的 `(H−z)⁻¹` 无对象 | `[T]` | `GL2bPrime:54`（`gap` 只吃 `→L[ℝ]`）；`Spectral/Unbounded:1–20` | 三处已在本文件修 |
| **0.4** | Mathlib v4.33.0-rc1 **有** `resolvent`（有界代元）、`MemLp`、`ProbabilityMeasure.tendsto_iff_forall_integral_tendsto`；**无** 无界自伴预解/谱定理 | `[T]` | `Algebra/Algebra/Spectrum/Basic:79`；`Unbounded:1–20` 自述 | 物理无界面改 `HasGapAbove`（`Unbounded:82`） |
| **0.5** | Clay `[NG]`：两条只是缺件签名，三格仍开 | `[NG]` | §2 | 攻 `N1` 的 `memLp` 侧（`L3` 归约） |

**一句话**：两件都写得出来；写不出来的不是类型，是内容。

## §1 两个声明（`axiom`，不是 `theorem`）

**1.0 前置（D0 引用面；本草案不定义，也不新建 D0 文件）**

```lean
variable {Ω H Obs : Type*} [MeasurableSpace Ω]
variable [NormedAddCommGroup H] [InnerProductSpace ℂ H] [CompleteSpace H]
variable (μ : Measure Ω) [IsProbabilityMeasure μ]
variable (ι : ℝ → Obs → Ω → H)                    -- ι_a[O] : Ω → H（原稿缺 `Ω →`）
variable (Ha : ℝ → H →L[ℂ] H) (Hlim : H →L[ℂ] H)  -- 有界代表；无界面见 §2(c)
```

**1.1 `N1` = `hy6_N1_iota_l2_seam`**

```lean
structure IotaL2Seam (μ : Measure Ω) (ι : ℝ → Obs → Ω → H)
    (Ha : ℝ → H →L[ℂ] H) (Hlim : H →L[ℂ] H) : Prop where
  memLp  : ∀ (a : ℝ), 0 < a → ∀ O : Obs, MemLp (ι a O) 2 μ
  strong : ∀ (z : ℂ), 0 < z.im → ∀ ψ : H,
    Tendsto (fun a : ℝ => resolvent (Ha a) z ψ) (𝓝[>] (0 : ℝ)) (𝓝 (resolvent Hlim z ψ))

axiom hy6_N1_iota_l2_seam : IotaL2Seam μ ι Ha Hlim
```

- **名字**：`hy6_N1_iota_l2_seam`（即 HY6 §1 的 `N1`；结构名 `IotaL2Seam`）。
- **假设**：`memLp`（C 的对象面，供 `T4` 的 OS 重构）；`strong`（C 的预解面，供 `T1` 的 `hres : StrongResolvent Ha H`）。
- **结论**：`IotaL2Seam μ ι Ha Hlim` 成立 —— C 的两面同时成立。
- **符号**：Mathlib `resolvent a z = (algebraMap z − a)⁻¹ʳ`（`Basic.lean:79`），故本式 = `−((Ha a − z•1)⁻¹)`；整体负号不影响 `Tendsto`。

**1.2 `N4` = `hy6_N4_b_floor`**

```lean
axiom hy6_N4_b_floor (s1 M0 chi : ℕ → ℝ → ℝ) (βA : ℝ) (hβA : 0 < βA)
    (hχ : ∀ (L : ℕ), 2 ≤ L → ∀ β : ℝ, βA ≤ β → 0 < chi L β)
    (hM : ∀ (L : ℕ), 2 ≤ L → ∀ β : ℝ, βA ≤ β → 0 < M0 L β) :
    ∃ c : ℝ, 0 < c ∧ ∀ (L : ℕ), 2 ≤ L → ∀ β : ℝ, βA ≤ β →
      s1 L β ≥ c * (M0 L β / chi L β)
```

- **名字**：`hy6_N4_b_floor`（即 HY6 §1 的 `N4`）。
- **假设** `hβA`：地板声明在 `β ≥ βA` 上，`βA > 0` ⇒ `β = 0` 不在域内，`hχ` 可满足（非空）。
- **假设** `hχ`/`hM`：`χ > 0`、`M₀ > 0`，正是 `T3` 把 `s₁ ≥ cM₀/χ` 翻成 `factor ≥ c` 所需的两前件。
- **结论**：`∃ c > 0`，`c` 与 `L, β` 无关；`2 ≤ L` 排除 `L = 0,1`（无空间 plaquette ⇒ `χ = 0` ⇒ Lean 的 `x/0 = 0` 使不等式失去意义）。

## §2 三判

| 声明 | (a) 陈述完整（量词齐、类型对） | (b) 与 `T14`/`T15` 消费形式匹配 | (c) Mathlib 现成定义 | 下一手 |
|---|---|---|---|---|
| **N1** | 修后 ✓：`a>0`、`z.im>0`、`ψ` 三量词齐；`ι : ℝ→Obs→Ω→H`、`μ` 概率、`Ha/Hlim : H→L[ℂ]H` 自洽。原稿 ✗：`ι a O : E` 与 `MemLp` 冲突；`z : ℂ` 无 ℂ-模作用于 `→L[ℝ]`。 | `T14` ✓：`T1` 吃 `hres : StrongResolvent Ha H`，`strong` 逐字是其 `Tendsto` 面。`T15` 半匹配：OS 重构要 `S^{(n)}` 的弱极限，`memLp` 只给对象、不给 `ι_a → ι`（缺 `GL4Transfer:52` 那一型的 `Tendsto`）。 | `MemLp`（`GL4L2:54`）；`ProbabilityMeasure.tendsto_iff_forall_integral_tendsto`（`GL4Transfer:57`）；`Filter.Tendsto`、`nhdsWithin`、`Filter.liminf`（`Order/LiminfLimsup:69`）；`resolvent`/`resolventSet`（`Basic:79`，仅**有界**代元，需 `[Ring A][Algebra ℂ A]`，`A = H→L[ℂ]H` 的实例未核）。**无** `LinearPMap` 预解（`Unbounded:1–20`） | 若 `N1` 真可归约到 `L3`，诚实形是 `theorem … (hL3 : L3Bound) : IotaL2Seam …`，`axiom` 只在 `L3` 未闭时占位 —— 先写那 4 行归约 |
| **N4** | 原稿 ✗：`L`/`chi` 无类型、`∃c>0` 缺 `χ>0`、`L=0,1` 分母为 0。修后 ✓：`L : ℕ`、`β : ℝ`、`chi/M0/s1 : ℕ→ℝ→ℝ` 全显式。 | `T14` ✓ 经 `T3`+`T13`（`hχ,hM` 恰是 `T3` 的两前件）；`T15` ✓ 直接吃实数不等式 `s₁ ≥ cM₀/χ`。 | 无 Mathlib 对象；只有代数重写 `div_le_iff₀`/`le_div_iff₀`（`Reduction:93,98`；`susceptibility_squeeze` 在 `:101`）与 `mul_le_mul_of_nonneg_left`。`s1/M0/chi` 属 D0，无一物可代。 | 把 `c` 与 `T9` 的 `(1−ε)c_L` 对齐；改 `L→∞` 序列写法时重判 |

- 两 `axiom` 都进 `#print axioms` 审计：机检 `E1–E5` 不依赖它们；`T16` 一旦引用，`#print axioms hy6_clay_of_three` 会列出 `hy6_N1_*`/`hy6_N4_*` —— 缺件可见，这是要的。
- 签名未过 elaboration ⇒ 不据此宣称任何新机检。

## §3 台账

**3.1 `AGENTS §1` 四要素**

| # | 工具 | 状态 | 实际调用（逐字） | 命中 |
|---|---|---|---|---|
| 1 | `pwsh` | **OK**（exit 0） | `Get-ChildItem research\v3\species -Name`；`Get-Content` 定段（6 份 Lean、1 份 md）；`Select-String`（Mathlib 定义与行号、两副本 SHA256） | 目录全列；`resolvent` 定义 1、`liminf` 2、`HasGapAbove` 2；两 `Unbounded.lean` 哈希相同 |
| 2 | `grep` | **OK**（非 provider failure） | ①`HY6-ASSEMBLE\|T16\|N4 = B\|N1 = C` @ `*.md` ②`MemLp\|resolvent\|Tendsto\|liminf\|def gap` @ `*.lean` ③`N1\|N4\|b_floor\|iota_l2` @ `species/*.md` | ①**8** ②**100**（显示上限） ③**250 / 1398**（截断，spill 落盘；对判词贡献 0 —— 槽名不是文件名） |
| 3 | `read` | **OK**（带行号正文） | `HY6-ASSEMBLE`(130 全)·`Reduction:60–190`·`GL2bPrime:1–70`·`Basic:1–80`·`Unbounded:1–140`·`GL4L2:200–245`·`GL4Transfer:40–70` | 全返回 |
| 4 | `web_*` | 未调用 | — | 0（本任务不需外部文献） |
| 5 | 失败标记 | **TOOL-FAIL 0；空响应 0；PATH-MISS 0** | — | — |

**3.2 直读 vs 转引**

- 直读（行号级）：`HY6` 全 130 行；`#3` 的六份 Lean 定段；Mathlib `Spectrum/Basic:55–85`（`resolvent` 的 section 变量）与 `Rayleigh.lean` 命中行（只为确认“只有有界 `resolvent`”）。
- 只经命中行/转引：`HX7-PATH:33,89`（C = 格 5 及其下一手，未开全文）·`HQ1-L3:191`（经 `HX7` 引）·`T9`/`T12` 只在 `HY6` 表内（未开 `HV8`/`HU3`）⇒ 记转引级，不进 (a)(c)。

**3.3 边界**

1. 本文件不证 `N1`/`N4`：两条都是 `axiom`，`[NG]`。
2. 签名未过 elaboration：`resolvent` 的 `Ring`/`Algebra ℂ` 实例、`MemLp` 的 `2 : ℝ≥0∞` 强制、`IsProbabilityMeasure` 均按仓库同名用法写，未机检。
3. 未编译、未 `git commit`、零 MC；无百分数；无 QED。
4. **总下一手**：先写 `L3 ⇒ memLp` 的 4 行，再按 `T9` 把 `c` 与 `(1−ε)c_L` 对齐；两件都闭后才轮到 `T13`/`T14` 的 `hB`。

*HV5-LEAN · `N1`/`N4` 各一条 `axiom` 声明；三处原稿类型错已修；Mathlib 无无界预解 ⇒ 物理面另立 `HasGapAbove`。*
