# M02-LEANC BAR — `C̄ ≡ N` 的算术核：mathlib 跑通 + 25 条定理、`sorry` = 0

> 第 N 轮 · 物种 M02。口径逐字取 `f07` 的裁决（`research/v3/species/F07-ADJCBAR.md` §0 行 3 / 5 / 7）。
> 只落两件：本文件 + `research/v3/lean/M02-CBAR.lean`。**未 git commit。**
> 标签 `[T]/[C]/[NG]/[M]/[E]`。**No QED；Clay 0%。**

---

## 0 · 判词（一表）

| # | 判决 | 等级 |
|---|---|---|
| 1 | **mathlib 跑通**：`import Mathlib` + `#check Real.exp_pos` → `Real.exp_pos (x : ℝ) : 0 < Real.exp x`，`exit=0` | `[T]` |
| 2 | **路径不是 `vendor\OSforGFF\.lake`**：完整的 8,279 个 `.olean` 在 `vendor\OSforGFF\.lake\packages\mathlib\.lake\**build.badcache**\lib\lean`；同包的 `…\.lake\build\lib\lean` 只有 2,136 个且**无根 `Mathlib.olean`** ⇒ 在那里 `import Mathlib` 必失败 | `[T]` |
| 3 | **toolchain 必须 `v4.33.0-rc1`**：`v4.36.0` 读同一份 `.olean` 报 `failed to read file …Mathlib.olean, incompatible header`，`exit=1`。（`l01` 报 `4.36.0` 读数成功，本机复现不了。） | `[T]` |
| 4 | **`LEAN_PATH` 必须含 toolchain 的 `lib\lean`**：一旦设置 `LEAN_PATH` 就顶替默认搜索路径，漏了它 `import Mathlib` 在 `Std.Http.Data.Body.Basic` 上报 `failed to read file`。 | `[T]` |
| 5 | **任务书的 `thresh` 式子对**：`2N²·48/eΔ = 19.2/e = 7.0632853`（`eΔ = 20e = 54.3656366`）。但 **`7.0633` 是上取整**：`T = 7.06328… < 7.0633`（机器证：`7.063 < T < 7.0633`）。 | `[T]` |
| 6 | **(c) 任务书的 `window` 式子错**：`48/(8/(2·2²·48/(e·20))) = 42.37`，不是 `3.53`。正确窗口比 = `48·2N²/(eΔ·C̄) = 48/(5e) = 3.5316426`；窗口 `β_w = 1/(5e) = 0.073575888`（`f07` 的 `0.073576 = 3.532×(1/48)` 一致：`3.53164` 三位小数就是 `3.532`）。 | `[T]` |
| 7 | **核心链 `sorry`-free 机器验证**：25 条定理，`lean` `exit=0`，`#print axioms` 12 条全为 `[propext, Classical.choice, Quot.sound]`（**无 `sorryAx`**）。 | `[T]` |
| 8 | **判死线 = (a) + (c)**：mathlib 可用 ⇒ 走 (a)；窗口式子的**精确形式**与任务书不同 ⇒ 同时触发 (c)，以 `f07` 为准逐字纠正。 | `[T]` |

---

## 1 · 环境（实测，原样贴）

### 1.1 定位几何

```
D:\Dev\math\vendor\OSforGFF\.lake\packages\mathlib\.lake\build\lib\lean
     → 2136 个 .olean，**没有 Mathlib.olean**       （import Mathlib 在此必挂）
D:\Dev\math\vendor\OSforGFF\.lake\packages\mathlib\.lake\build.badcache\lib\lean
     → 8279 个 .olean，含 Mathlib.olean (764,408 B)  ← **真正可用的那份**
lean-toolchain: leanprover/lean4:v4.33.0-rc1        lake-manifest rev: 79d0395a1825a6264ad5d269e35e60537518955e
```

`elan toolchain list`：

```
leanprover/lean4:v4.33.0-rc1
leanprover/lean4:v4.36.0        (default, resolved from 'stable')
```

### 1.2 有效 `LEAN_PATH`（顺序无关，但必须**同时**含 toolchain 的 `lib\lean`）

```
C:\Users\hooya\.elan\toolchains\leanprover--lean4---v4.33.0-rc1\lib\lean
D:\Dev\math\vendor\OSforGFF\.lake\packages\mathlib\.lake\build.badcache\lib\lean
…\.lake\packages\{aesop,batteries,importGraph,Qq,proofwidgets,plausible,
                  LeanSearchClient,BochnerMinlos,GaussianField,kolmogorov_extension4}\.lake\build\lib\lean
```

### 1.3 冒烟测试（逐字输出）

```
$ lean --version
Lean (version 4.33.0-rc1, x86_64-w64-windows-gnu, commit 62eed1db4d67327ec8120be05f1a1b0847d74561, Release)

$ lean env_check.lean     # 内容：import Mathlib ; #check Real.exp_pos ; #check Real.exp_one_gt_d9
Real.exp_pos (x : ℝ) : 0 < Real.exp x
Real.exp_one_gt_d9 : 2.7182818283 < Real.exp 1
exit=0
```

反例（**为什么不是 4.36.0**）：

```
$ <v4.36.0>\lean.exe env_check.lean
error: failed to read file '…\build.badcache\lib\lean\Mathlib.olean', incompatible header
exit=1
```

反例（**为什么 `LEAN_PATH` 要带 toolchain 的 lib**）：

```
$ <v4.33.0-rc1>\lean.exe env_check.lean      # LEAN_PATH 未含 toolchain lib
error: failed to read file 'c:\Users\hooya\.elan\toolchains\…v4.33.0-rc1\lib\lean\Std\Http\Data\Body\Basic.olean'
exit=1
```

---

## 2 · 式子核准（先判精确形式，再形式化）

设 `N = 2`（`SU(2)`）、`Δ = 20`、`Ĉ = 4`、`C̄ = c_f/u = N`、`eΔ = 20e`。

| 量 | 任务书写的 | 机器核出的精确值 | 判 |
|---|---|---|---|
| `c_f/u` | `= N`（`N=2`） | `β/(2N) ÷ (β/(2N²)) = N`，**恒等式**，`N=2 ⇒ 2` | **对** |
| `Ĉ·C̄` | `4N = 8` | `4·2 = 8` | **对** |
| 阈值 | `2N²·48/eΔ` | `384/(20e) = 19.2/e = **7.0632853**` | **对**；`7.0633` 是**上取整**（`T < 7.0633`） |
| `eΔ` | `= 20e = 54.3656` | `54.3656366`，机器夹在 `(54.3656, 54.3657)` | **对** |
| 比值 | `8 / 阈值` | `= 5e/12 = **1.1326174**`（`= f07 §1.4 的 eΔ·N·u` 在 `β=1/12`） | **对** |
| 窗口 | `48/(8/(2N²·48/(eΔ))) = 3.53` | 上式算得 **42.37**：把「比值」又除了 48。正确：`48·2N²/(eΔ·C̄) = 48/(5e) = **3.5316426**` | **错 ⇒ 更正** |
| 窗口 `β` | `β < 3.53×(1/48)` | `1/(5e) = **0.073575888**`；`×48 = 3.53164` ✓ | **对（数值）**，式形需如上更正 |

**关键结构**（机器证 `thresh_eq_eight_div_x`）：`阈值 = (Ĉ·C̄) / (eΔ·N·u |_{β=1/12})`，
于是 `8 > 阈值 ⟺ eΔ·N·u > 1 ⟺ e > 2.4`（`thresh_iff`）。`f07` 的 `7.0633` 不是孤立数 ——
它是「余量 = `e/(12/5) = 5e/12`」的另一种写法。

---

## 3 · `M02-CBAR.lean`：25 条定理（全 `sorry`-free）

| 组 | 定理 | 内容 |
|---|---|---|
| 归一 | `cf_SUN2` / `u_SUN2` | `c_f(1/12,2) = 1/48`、`u = 1/96`（`f07` 实测 `⟨(1/N)trU⟩ = 0.0208273`，差 0.03%） |
| `C̄ ≡ N` | `cbar_eq_N` / `cbar_eq_N_two` | `c_f/u = N`（一般 `N`、`β≠0`）；`N=2 ⇒ 2` |
| `2N²` | `phi0_over_u` / `phi0_over_u_two` | `Φ(0)/u = 2N²`；`N=2 ⇒ 8`（与 `Ĉ·C̄` **同名不同物**） |
| `Ĉ·C̄` | `chat_times_cbar` / `chat_times_cbar_N` | `4·2 = 8`；通式 `4N` |
| 阈值 | `thresh_eq` / `thresh` | `= 192/(10e)`；**`8 > 2·2²·48/(e·20)`** |
| 等价/锐性 | `thresh_iff` / `sharp_at_e_2_4` / `thresh_eq_eight_div_x` | `8 > T ⟺ e > 2.4`；`e=12/5` 时 `T = 8`；`T = 8/(5e/12)` |
| `x` | `x_at_beta_twelfth` / `ratio` / `x_gt_one` | `eΔ·N·u = 5e/12`；`8/T = 5e/12`；**`x = 1.1326 > 1`（`f07` 的 ★ 缺口）** |
| 窗口 | `window_general` / `window` / `window_beta` | `48·2N²/(eΔN) = 48·2N/(20e)`；`N=2` ⇒ `48/(5e)`；`β_w = 1/(5e)` |
| 数值夹逼 | `window_bounds` / `window_beta_bounds` / `thresh_bounds` / `eDelta_bounds` | `3.53 < 48/(5e) < 3.54`；`0.073575 < β_w < 0.073576`；`7.063 < T < 7.0633`；`54.3656 < 20e < 54.3657` |
| 有理/整数层备份 | `rational_certificate` / `integer_certificate` | `7063/1000 < 384/(20·(27182818283/10¹⁰))`；`7063·20·27182818283 < 384·10¹⁴` |

**`sorry` 数 = 0**。位置：无（文件内 `sorry` 三个字符只出现在注释/文档串的行 28、100、234，
由 `Select-String '\bsorry\b'` 逐条核过，均为 `` `sorry` 数 = 0`` 这类说明文字）。

**编译（末次，逐字）**：

```
$ lean M02-CBAR.lean
'M02.cbar_eq_N' depends on axioms: [propext, Classical.choice, Quot.sound]
'M02.thresh' depends on axioms: [propext, Classical.choice, Quot.sound]
'M02.thresh_iff' depends on axioms: [propext, Classical.choice, Quot.sound]
'M02.ratio' depends on axioms: [propext, Classical.choice, Quot.sound]
'M02.x_gt_one' depends on axioms: [propext, Classical.choice, Quot.sound]
'M02.window' depends on axioms: [propext, Classical.choice, Quot.sound]
'M02.window_bounds' depends on axioms: [propext, Classical.choice, Quot.sound]
'M02.window_beta_bounds' depends on axioms: [propext, Classical.choice, Quot.sound]
'M02.thresh_bounds' depends on axioms: [propext, Classical.choice, Quot.sound]
'M02.eDelta_bounds' depends on axioms: [propext, Classical.choice, Quot.sound]
'M02.rational_certificate' depends on axioms: [propext, Classical.choice, Quot.sound]
'M02.integer_certificate' depends on axioms: [propext, Classical.choice, Quot.sound]
EXIT=0   （12 行 = 全部输出；无 error、无 warning；墙钟 19.5 s）
```

`#print axioms` 干净（无 `sorryAx`）是「`sorry`-free」的**机器证据**，不是自述。

---

## 4 · 复现（一条命令）

```powershell
$lean="C:\Users\hooya\.elan\toolchains\leanprover--lean4---v4.33.0-rc1\bin\lean.exe"
$tc  ="C:\Users\hooya\.elan\toolchains\leanprover--lean4---v4.33.0-rc1\lib\lean"
$base="D:\Dev\math\vendor\OSforGFF\.lake\packages"
$dirs=@($tc,"$base\mathlib\.lake\build.badcache\lib\lean")
foreach($n in @("aesop","batteries","importGraph","Qq","proofwidgets","plausible",
                "LeanSearchClient","BochnerMinlos","GaussianField","kolmogorov_extension4")){
  $p="$base\$n\.lake\build\lib\lean"; if(Test-Path $p){ $dirs+=$p } }
$env:LEAN_PATH=($dirs -join ";")
& $lean "D:\Dev\math\research\v3\lean\M02-CBAR.lean"; "exit=$LASTEXITCODE"
```

`lake env` 亦可（在 `vendor\OSforGFF` 下），但只要 `LEAN_PATH` 漏 toolchain lib 就挂 —— 见 §1.3。

---

## 5 · 与 `f07` 的逐字对照 / 限定

| `f07` 数 | 本文件机器值 | 判 |
|---|---|---|
| `c_f/u = N = 2.0000` | `cbar_eq_N_two`（恒等式） | 一致 |
| `Φ(0)/u = 2N²`，`N=2 → 8` | `phi0_over_u` / `phi0_over_u_two` | 一致 |
| `eΔ = 54.3656` | `54.3656 < 20e < 54.3657` | 一致 |
| `eΔ·N·u = 1.132617` | `ratio = 5e/12 = 1.1326174` | 一致（`x_at_beta_twelfth` 同） |
| 窗口 `0.073576 = 3.532×(1/48)` | `β_w = 0.073575888`，`48β_w = 3.5316` | 一致（三位小数取整） |
| 阈值 `7.0633`（任务书） | `T = 7.06328…`，机器夹逼 `7.063 < T < 7.0633` | **上取整**，方向无害（`8 > 7.0633 > T`） |
| 任务书 `window` 式 | `48/(8/(2·2²·48/(e·20))) = 42.37` | **式错**，已更正为 `48·2N²/(eΔ·C̄)` |

**限定**：以上全是 `SU(2)`、`Δ = 20`、`β` 归一取 `u = c_f/N = β/(2N²)` 口径下的算术；
**不含**任何跨 `β` 的一致性或连续极限断言。窗口 `β < 1/(5e)` 是**该口径**的窗口，
精确 `SU(2)` 口径（`f07` 的 `0.0735925`）本文件不涉及。

**No QED；Clay 0%。**
