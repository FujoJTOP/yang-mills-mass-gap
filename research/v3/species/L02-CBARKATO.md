# L02 · C̄ ≤ 1.054（c07）与 Kato ⇒ m>0（f04）

> 物种 L02 交付。文件：`research/v3/lean/L02-CBARKATO.lean`。
> 标签：`[T]` 真定理（机器验证）/ `[C]` 条件定理或待形式化 / `[NG]` no-go。

---

## 0 工具状态（先报，按 AGENTS §1）

| 项 | 值 | 状态 |
|---|---|---|
| `lean` | `Lean (version 4.36.0, x86_64-w64-windows-gnu, commit 611c117…)` | OK |
| `lake` | `Lake version 5.0.0-src+611c117 (Lean version 4.36.0)` | OK |
| Mathlib | `import Mathlib` → `error: unknown module prefix 'Mathlib'`；`No directory 'Mathlib' or file 'Mathlib.olean' in the search path entries` | **TOOL-FAIL** |
| 全仓 Mathlib 探测 | `Get-ChildItem D:\AI,D:\Dev -Filter Mathlib.olean -Recurse`（工具 OK）→ 命中 **0** | 0 命中 |
| 已有 Mathlib 工程 | `docs/zenodo-v2/ym-v2-lean/lakefile.lean` require mathlib `79d0395` (toolchain `v4.33.0-rc1`)，但 `.lake` 目录不存在 ⇒ 从未构建 | 不可用 |
| Mathlib master 工具链 | `lean-toolchain` = `leanprover/lean4:v4.37.0-rc1`，与本机 4.36.0 不匹配 | 不匹配 |
| `Rat`（core ℚ） | `OfNat ℚ 288` → `synthInstanceFailed`；`HDiv ℚ Nat` → 同 | core 无 ℚ 算术实例 |

**结论**：`Real.exp` / `Real.log` / `spectrum` / `Metric.infDist` **在本机一律不可调用**。
因此交付走 **`Nat` 上的精确整数证书**（内核可验）+ 解析转移步的手写证明（标「待形式化」）。

编译命令与实输出（`lean L02-CBARKATO.lean`，无参数直接编译）：

```
1.053946
1.054191
1.053464
'L02.delta_over_96_eq_five_over_24' does not depend on any axioms
'L02.worst_case_at_N_two' depends on axioms: [propext, Quot.sound]
'L02.exp_bound_certificate' does not depend on any axioms
'L02.exp_bound_strict' does not depend on any axioms
'L02.n_two_insufficient' does not depend on any axioms
'L02.cbar_product_bound' does not depend on any axioms
'L02.cbar_product_exact' does not depend on any axioms
'L02.gap_ratio_in_unit_interval' depends on axioms: [propext, Classical.choice, Quot.sound]
'L02.gap_ratio_collapses_at_boundary' depends on axioms: [propext, Quot.sound]
'L02.kato_instance_half_quarter' does not depend on any axioms
=== EXIT=0 ===
```

首行 `1.053946` 是 `#eval (288/283)^3`（数值对账，非证明）。
无 error、无 warning、**无 `sorry`**；`#print axioms` 全部不含 `sorryAx`。

---

## 1 ① `C̄ ≤ e^{Δu} = 1.054`（c07）

### 1.1 链条（逐步，每步给出类型）

1. `Δ = 20`（纯数）；`u = c_f/N² = β/(2N²)`。
2. `β < 1/48 ⇒ u < 1/(96N²)`。
3. `Δu < 20/(96N²) = 5/(24N²)`。 **← 机器已验（`delta_over_96_eq_five_over_24`）**
4. `N ≥ 2 ⇒ 5/(24N²) ≤ 5/96`（最坏在 `N = 2`，且 `N = 2` 取等）。 **← 机器已验（`worst_case_at_N_two`, `N_two_sharp`）**
5. `exp` 单调 ⇒ `e^{Δu} ≤ e^{5/96}`。 **← 解析，待形式化**
6. `e^{5/96} ≤ (288/283)³`（初等上界，见 1.3）。 **← 解析，待形式化**
7. `(288/283)³ ≤ 1054/1000 = 1.054`。 **← 机器已验（`exp_bound_certificate`，严格版 `exp_bound_strict`）**
8. `Ĉ·C̄ ≤ 4 × 1.054 = 4.216 < 7.0633`。 **← 机器已验（`cbar_product_bound`, `cbar_product_exact`）**

### 1.2 机器验证的六条（`Nat`，交叉相乘，无除法）

| 定理 | 整数命题 | 依据 |
|---|---|---|
| `delta_over_96_eq_five_over_24` | `20·24 = 5·96` | `decide` |
| `worst_case_at_N_two` | `N ≥ 2 ⇒ 96 ≤ 24·N²` | `Nat.mul_le_mul` + `omega` |
| `N_two_sharp` | `96 = 24·(2·2)` | `decide` |
| `exp_bound_certificate` | `1000·288³ ≤ 1054·283³` | `decide` |
| `exp_bound_strict` | `1000·288³ < 1054·283³` | `decide` |
| `cbar_product_bound` | `4·1054·10⁴ < 70633·10³`，即 `4.216 < 7.0633` | `decide` |

**负对照（诚实性）**：`n_two_insufficient : 1054·187² < 1000·192²`，即 `(192/187)² = 1.054191 > 1.054`
⇒ `e^x ≤ (1 − x/2)^{-2}` 这条更弱上界**在此失败**，必须用 `n = 3`。已机器验证。

### 1.3 解析转移步（手写，待形式化）

**引理 A（初等 exp 上界）**：对 `x ∈ [0, n)`，
```
e^x ≤ (1 − x/n)^{-n}.
```
证：对 `t > −1` 有 `e^t ≥ 1 + t`。取 `t = −x/n ∈ (−1, 0]`：`e^{−x/n} ≥ 1 − x/n > 0`，
取倒数 `e^{x/n} ≤ 1/(1 − x/n)`，两边 `n` 次幂。 ∎
（Mathlib 对应件：`Real.add_one_le_exp : x + 1 ≤ Real.exp x`。**本机不可调用**。）

取 `x = 5/96, n = 3`：`e^{5/96} ≤ (283/288)^{-3} = (288/283)³ = 1.053946…`。
再由 `e^{5/96} = 1.053464…`（`#eval Float.exp` 对账）可见余量 `1.054 − 1.0535 ≈ 5.4×10⁻⁴ > 0`。

**引理 B（单调）**：`x ≤ y ⇒ e^x ≤ e^y`（`Real.exp_le_exp`）。本机不可调用。

### 1.4 ① 的等级裁定

- **算术核心（步 3、4、7、8）：`[T]`，机器验证，`sorry = 0`。**
- **整句 `Real.exp (5/96) ≤ 1.054`：`[C]`，待形式化**（缺 `Real.exp` 与引理 A/B；
  一旦 Mathlib 可用，`Real.add_one_le_exp` + `Real.exp_le_exp` + `nlinarith` 即可闭合，无实质障碍）。
- 严格说，机器证书证明的是 **`e^{5/96} ≤ 1.054` 在初等上界下的精确有理化形态**，不是 `exp` 本身；
  这一点不缩水表述：**`C̄ ≤ 1.054` 的数值内容已机器验证，`exp` 的接口未验证。**

---

## 2 ② `Kato ⇒ m > 0`（f04）

### 2.1 `spec` / `dist` 在 Mathlib 里有吗？

- **`spectrum`**：Mathlib 有（`Algebra.Spectrum.spectrum`，记 `spectrum 𝕜 a`），
  配 `spectrum.mem_iff`、`spectrum.isCompact`（紧自伴时）等。
- **`Metric.infDist`**：Mathlib 有（`Metric.infDist : Set α → Set α → ℝ`，另有 `Metric.infEdist`）。
- **本机核验状态**：**TOOL-FAIL**。`lean` 可执行但 `import Mathlib` 失败（`unknown module prefix`），
  无法 `#check spectrum` / `#check Metric.infDist`。
  按 AGENTS §1：**这不等于「不存在」，记为「未核验」**。上文名称来自训练记忆，**标 `[?]`，未逐字核**。
- 更关键：`spec T ∖ {sup spec T}` 不是闭集，`dist` 到它需 `Metric.infDist`（下确界，非 `dist`）；
  且 `sup` 不达时的退化正是 f04 的边界情形（见 2.4）。

### 2.2 纯实分析等价形式（**建议采用的陈述**）

> **定理 G（Gap ⇒ Mass）。**
> `H` 实 Hilbert，`T ∈ B(H)` 紧自伴，`λ₁ := ‖T‖ > 0`。
> 设
> (H1) `E₁ := ker(λ₁ − T)`，`dim E₁ = 1`（`λ₁` 单重）；
> (H2) `δ := λ₁ − ‖T|_{E₁^⊥}‖ > 0`。
> 令 `λ₂ := ‖T|_{E₁^⊥}‖`（`>0` 时即第二特征值，按 `λ₁ ≥ λ₂ ≥ …` 排序）。
> 则
> (i) `λ₂ ≤ λ₁ − δ`；
> (ii) 若 `λ₂ > 0`，则 `m := −log(λ₂/λ₁) ≥ −log(1 − δ/λ₁) > 0`。

### 2.3 证明（手写，可直接翻译成 Lean）

1. **谱定理**：`T` 紧自伴 ⇒ `H = ⊕̂_k ker(λ_k − T)`，`spec T = {λ_k} ∪ {0}`，`λ_k → 0`；
   且 `sup spec T = ‖T‖ = λ₁`。因 `T` 紧自伴，`‖T‖` 被达到，故 `λ₁ ∈ spec T` 是特征值，
   `E₁ ≠ 0`，`E₁` 有限维且闭（`dim = 1`）。
2. `E₁^⊥ = ⊕̂_{k≥2} ker(λ_k − T)`（Hilbert 直和的正交补）。
3. `T|_{E₁^⊥}` 仍紧自伴，`spec(T|_{E₁^⊥}) = {λ_k}_{k≥2} ∪ {0}`，
   故 `‖T|_{E₁^⊥}‖ = sup_{k≥2} |λ_k|`。按 `λ₁ ≥ λ₂ ≥ λ₃ ≥ …` 排序且 `λ₂ > 0` 时，
   `|λ_k| ≤ λ₂ (k ≥ 2)`，于是 `‖T|_{E₁^⊥}‖ = λ₂`。由 (H2)：`λ₂ ≤ λ₁ − δ`。 — (i)
4. `δ > 0` 且 `λ₂ > 0 ⇒ δ ≤ λ₁ − λ₂ < λ₁`，故 `0 < δ/λ₁ < 1`，
   即 `0 < 1 − δ/λ₁ < 1`。 — **算术核，已机器验证（`gap_ratio_in_unit_interval`）**
5. `λ₂/λ₁ ≤ 1 − δ/λ₁`，`log` 在 `(0,∞)` 严格增 ⇒
   `log(λ₂/λ₁) ≤ log(1 − δ/λ₁) < 0` ⇒ `m = −log(λ₂/λ₁) ≥ −log(1 − δ/λ₁) > 0`。 ∎

**这条链里没有 Kato。** 需要的只是谱定理 + 紧性。
`Kato` 真正的位置在**反方向**：由耦合参数的解析性/预解式连续性去**得到** `r_iso > 0`（或 `δ > 0`），
不是由 `r_iso > 0` 推出 `m > 0`。f04 的内容应拆成
`[Kato ⇒ r_iso>0]`（分析，难）∘ `[r_iso>0 ⇒ m>0]`（谱论，易，本文给出）。
**标 `[C]`：作为「给定 δ>0 的谱论推论」成立；作为「Kato 推出 m>0」则缺少 `Kato ⇒ δ>0` 一环。**

### 2.4 不可逆性与反例的两个读法（**必须钉住**）

反例 `T = mult((1+x)/4)` on `L²(0,1)`：`spec T = [1/4, 1/2]`，`sup = 1/2` 不达，
`dist(sup, spec ∖ {sup}) = dist(1/2, [1/4,1/2)) = 0`，故 `r_iso = 0`。**这一半成立。**

但 `m` 的数值依赖 `λ₂` 的定义，两个读法给出**不同**答案：

| 读法 | `λ₂` | `m` | 与 `r_iso = 0` 相容？ |
|---|---|---|---|
| `λ₂ := sup spec(T|_{E₁^⊥})` | `E₁ = 0` ⇒ `E₁^⊥ = H` ⇒ `λ₂ = 1/2` | `m = log(1/2 ÷ 1/2) = 0` | 相容（`m = 0`） |
| `λ₂ := inf spec(T|_{E₁^⊥})` | `1/4` | `m = log(1/2 ÷ 1/4) = log 2` | 相容（`m = log 2`，**与题述一致**） |

⇒ 题述「`m = log2`、`r_iso = 0`」只在**第二读法**（`λ₂ := inf`，即谱宽/谱底）下成立。
若 `m` 定义为 `−log(λ₂/λ₁)` 且 `λ₂ = sup(…)`，则该反例给 `m = 0`，**不能**充当
「`m>0` 而 `r_iso=0`」的反例。
**结论：`Kato ⇒ m>0` 的严格性成立，但反例的数值须按 `λ₂ := inf spec(T|_{E₁^⊥})` 记账。**
建议在正式陈述里把 `λ₂` 的定义写死（这是本条唯一的口径风险）。

### 2.5 结构旁注（`[T]`，易证，但决定了 f04 的价值）

`T` 紧自伴、`λ₁ = ‖T‖` 且 `λ₁` 单重 ⇒ `spec T` 在 `λ₁` 附近离散
（`λ_k → 0`，`λ₁` 是孤立点）⇒ **`r_iso > 0` 自动成立**。
所以在**紧**框架下 `Kato` 是多余的；f04 只有在 `T` **非紧**（`λ₁` 落在本质谱/连续谱里，
或 `r_iso = 0`）时才有内容——而那时 `δ > 0` 恰恰不可能，`m > 0` 必须换别的来源（如谱隙的
一致估计、或 IR 正则化后的粗模界）。这是 f04 的真正边界，**标 `[C]`，需在下轮明确框架**。

---

## 3 ② 的可行性判定

| 层次 | 可用性 | 判定 |
|---|---|---|
| `Real.exp`/`Real.log`/`spectrum`/`Metric.infDist` | 本机**无** Mathlib | **不可形式化**（TOOL-FAIL） |
| `δ>0 ⇒ m>0` 的算术核 | `Nat` | **已机器验证**（`gap_ratio_in_unit_interval` 等 3 条，`sorry = 0`） |
| 定理 G 的谱论部分 | 需 `spectrum` + 紧自伴谱定理 | Mathlib 有谱定理，但本机不可调用 ⇒ **待形式化** |
| `Kato ⇒ δ>0` | — | **缺环**，f04 的实质难点在此，标 `[C]` |

**可行性结论**：定理 G **可以**写成纯实分析/谱论形式且**证明是初等的**（2.3 的 5 步），
但其中「谱定理 + 紧自伴特征值展开」在 Mathlib 里是重活（`spectrum` + `CompactOperator` 谱定理），
**本机 25 分钟内不可能形式化**。建议形式化的最小切口：
把定理 G 降为 **`H = ℓ²(ℕ)` 上的对角算子**（`T e_k = λ_k e_k`，`λ_k ↓ 0`），
此时谱定理变成逐坐标，`λ₂ ≤ λ₁ − δ` 就是 `λ₂ = ‖T|_{E₁^⊥}‖` 的直接计算，
可在 core/mathlib 上用 `Finset.sup` 完成——这是下一轮的形式化靶。

---

## 4 死线三项

**(a) ① 是否完全证出（`sorry` = 0）？**
→ **算术核心是**：步 3/4/7/8 全部机器验证，`sorry = 0`，`#print axioms` 无 `sorryAx`。
→ **但「`Real.exp (5/96) ≤ 1.054`」整句未证**（无 `Real.exp`）。故**不宣布「`C̄` 的界升为机器验证」**；
   诚实等级：**算术机器验证 `[T]` + 解析接口 `[C]`**。

**(b) ② 是否可写成纯实分析形式且有证明？**
→ **是**：定理 G（2.2）+ 证明（2.3），初等谱论，**手写已完成**。
→ 但 (i) 本机不能形式化；(ii) `Kato ⇒ δ>0` 一环缺失，故整体标 `[C]`。

**(c) Lean/Mathlib 不可用 ⇒ 精确陈述 + 手写证明，标「待形式化」**：**部分适用**。
`lean` 可用、Mathlib 不可用 ⇒ 交付 = 机器验证的算术核 + 待形式化的解析/谱论陈述与手写证明。

---

## 5 文件与产物

- `research/v3/lean/L02-CBARKATO.lean` — 10 条定理，`decide`/`omega`/`Nat.mul_le_mul` 关闭，`sorry = 0`。
- 本文件。
- 未提交 git（遵守禁令）。

## 6 弃用/待办

1. 装 Mathlib（toolchain 需 `v4.37.0-rc1` 或 mathlib 对应 `v4.36.0` 的 commit + `lake exe cache get`）⇒ 方可证 `Real.exp` 整句。
2. 定理 G 的最小形式化切口：`ℓ²` 对角算子版本（§3 末）。
3. 钉住 `λ₂ := sup` vs `λ₂ := inf`（§2.4）——影响 f04 的反例记账。
