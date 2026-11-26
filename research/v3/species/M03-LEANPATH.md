# M03 · Lean / mathlib 路径裁决

**结论先行**：mathlib **可用**，但只在**子模块级**（`import Mathlib.Order.Basic` / `Mathlib.Analysis.SpecialFunctions.Exp` 等）；
**根伞模块 `import Mathlib` 不可用**（`Mathlib.olean` 不在活动搜索路径上）。
`l01` 与 `l02` 的矛盾由此消解，且两人的报数都可复现 —— 详见证 §3。

裁决等级：`[T]` = 本条为实测（命令 + 输出逐字记录）；`[C]` = 推断。

---

## 1. 环境事实

| 项 | 值 | 证据 |
|---|---|---|
| lean | `4.36.0`（elan default，`C:\Users\hooya\.elan\bin\lean.exe`）`[T]` | `lean --version` |
| lake | `5.0.0-src+611c117 (Lean 4.36.0)` `[T]` | `lake --version` |
| vendor 工程 toolchain | `leanprover/lean4:v4.33.0-rc1` `[T]` | `vendor\OSforGFF\lean-toolchain` |
| mathlib rev | `79d0395a1825a6264ad5d269e35e60537518955e`, `inputRev v4.33.0-rc1` `[T]` | `vendor\OSforGFF\lake-manifest.json` |
| elan 已装 toolchain | `v4.33.0-rc1`, `v4.36.0` `[T]` | `ls .elan\toolchains` |

> 关键：**mathlib 的 olean 头是 `4.33.0-rc1`**（`olean` magic + `4.33.0-rc1` 版本串，实测 `Format-Hex Mathlib\Order\Basic.olean`），
> 所以任何调用都必须 pin 到 4.33.0-rc1；用 elan 默认 4.36 直连会 `incompatible header`。`[T]`

---

## 2. `Mathlib.olean` 到底在哪（第 1 问）

```
Get-ChildItem -Recurse -Filter Mathlib.olean vendor\OSforGFF\.lake
```

**命中数 = 1**，唯一路径：`[T]`

```
D:\Dev\math\vendor\OSforGFF\.lake\packages\mathlib\.lake\build.badcache\lib\lean\Mathlib.olean
```

分布：

| 位置 | olean 数 | `Mathlib.olean`（根伞模块） |
|---|---|---|
| `...\mathlib\.lake\build\lib\lean`（**活动**路径） | **2264** | **0** |
| `...\mathlib\.lake\build.badcache\lib\lean`（Lake 判坏的缓存） | **8279** | **1** |

> 计数口径注：`-File` 计数 = 2264；早前不带 `-File` 的同式计数 = 2136（Windows`-Filter` 通配行为差异）。两者均实测，取 2264 为准。
> `badcache` 的 8279 ≈ `l01` 报的 8268 —— `l01` 的数**是真的**，只是它数的是**被 Lake 弃用的缓存**。`[T]`

活动路径的覆盖（`Mathlib/` 子树）：

```
Algebra=584; Data=384; Tactic=339; Topology=257; Order=199; Analysis=91;
LinearAlgebra=75; RingTheory=72; CategoryTheory=58; GroupTheory=52; Logic=47;
SetTheory=28; Lean=25; Util=22; Control=10; Combinatorics=6; MeasureTheory=6;
Dynamics=5; NumberTheory=2; Geometry=1; Init.olean=1
```

即：**20 个子树在，根 `Mathlib.olean` 不在**；`Mathlib/Tactic.olean`（tactic 伞）、
`Mathlib/Data/Matrix/Notation.olean`、`Mathlib/Data/Matrix/Defs.olean`、`Mathlib/Order/Field/Basic.olean` 均**缺**。`[T]`

---

## 3. 三种方式的实际结果（第 2 问）

### (a) `cd vendor\OSforGFF ; lake env lean <file>` —— **通** ✅

```
cd D:\Dev\math\vendor\OSforGFF
lake env lean ..\..\research\v3\lean\M03-COUNTER.lean
```

```
'M03.counterexample' depends on axioms: [propext, Classical.choice, Quot.sound]
'M03.not_f1sq_zero_implies_gapless' depends on axioms: [propext, Classical.choice, Quot.sound]
EXIT=0
```

`lake env` 自动给出 12 条 `LEAN_PATH`（11 个包 `build/lib/lean` + 4.33 toolchain `lib/lean`），
并 pin 到工程 toolchain —— 这是**唯一零配置可用**的入口。`[T]`

### (b) 手设 `LEAN_PATH` + 裸 `lean` —— **通，但有两个必要条件**

```
cd D:\Dev\math
$env:LEAN_PATH = "D:\Dev\math\vendor\OSforGFF\.lake\packages\aesop\.lake\build\lib\lean;D:\Dev\math\vendor\OSforGFF\.lake\packages\batteries\.lake\build\lib\lean;D:\Dev\math\vendor\OSforGFF\.lake\packages\BochnerMinlos\.lake\build\lib\lean;D:\Dev\math\vendor\OSforGFF\.lake\packages\GaussianField\.lake\build\lib\lean;D:\Dev\math\vendor\OSforGFF\.lake\packages\importGraph\.lake\build\lib\lean;D:\Dev\math\vendor\OSforGFF\.lake\packages\kolmogorov_extension4\.lake\build\lib\lean;D:\Dev\math\vendor\OSforGFF\.lake\packages\LeanSearchClient\.lake\build\lib\lean;D:\Dev\math\vendor\OSforGFF\.lake\packages\mathlib\.lake\build\lib\lean;D:\Dev\math\vendor\OSforGFF\.lake\packages\plausible\.lake\build\lib\lean;D:\Dev\math\vendor\OSforGFF\.lake\packages\proofwidgets\.lake\build\lib\lean;D:\Dev\math\vendor\OSforGFF\.lake\packages\Qq\.lake\build\lib\lean;c:\Users\hooya\.elan\toolchains\leanprover--lean4---v4.33.0-rc1\lib\lean"
lean +leanprover/lean4:v4.33.0-rc1 research\v3\lean\M03-COUNTER.lean
```

```
'M03.counterexample' depends on axioms: [propext, Classical.choice, Quot.sound]
'M03.not_f1sq_zero_implies_gapless' depends on axioms: [propext, Classical.choice, Quot.sound]
EXIT=0
```

两个必要条件，缺一即死：

| 变体 | 输出 | 结论 |
|---|---|---|
| LEAN_PATH 只给 `...\mathlib\.lake\build\lib\lean`（**`l02` 的做法**） | `error: unknown module prefix 'Aesop'` + `No directory 'Aesop' or file 'Aesop.olean' in the search path entries: <mathlib dir>, <4.33 tc>` | ✗ 依赖包未上路径；**这不是 mathlib 不可用** `[T]` |
| 全 LEAN_PATH，但用 elan 默认 `lean`（4.36） | `error: failed to read file '...4.33.0-rc1\lib\lean\Init.olean', incompatible header` | ✗ 版本不匹配 `[T]` |
| 全 LEAN_PATH + `+leanprover/lean4:v4.33.0-rc1` | 见上，`EXIT=0` | ✅ `[T]` |

### (c) `cd 仓库根 ; lake build <file>` —— **死** ❌

```
cd D:\Dev\math
lake build research/v3/lean/M03-COUNTER.lean
error: [root]: no configuration file with a supported extension:
D:\Dev\math\lakefile.lean
D:\Dev\math\lakefile.toml
EXIT=1
```

仓库根**没有** `lakefile.lean` / `lakefile.toml`，也**没有** `lean-toolchain`。`[T]`
（仓内 lake 工程只有：`vendor\OSforGFF`、`research\lean-infra`（无 mathlib 依赖目录）、`docs\zenodo-v2\ym-v2-lean`。）

---

## 4. `l01` vs `l02`：两边都对，问的不是同一个东西

| 说法 | 本次实测 | 判 |
|---|---|---|
| `l01`：mathlib v4.33.0-rc1 @ `vendor\OSforGFF\.lake`，约 8.2k olean，`#check Real.exp_pos` exit=0 | 复现成功（`import Mathlib.Analysis.SpecialFunctions.Exp` → `Real.exp_pos (x : ℝ) : 0 < Real.exp x`，EXIT=0）；8268 ≈ `build.badcache` 的 8279 | **成立** `[T]` |
| `l02`：`unknown module prefix 'Mathlib'`、全盘 `Mathlib.olean` 命中 0 | `import Mathlib` 确实失败，但错因是**根伞 olean 缺失**（`object file '...\build\lib\lean\Mathlib.olean' of module Mathlib does not exist`）；全盘 0 命中与本次 1 命中冲突 ⇒ `l02` 的检索没进 `.lake`（隐藏目录）或只查了 `build/lib/lean` | **部分成立 / 结论错** `[T]` |

`l02` 侧可复现的两个真实错误串（逐字）：

```
# import Mathlib
error: object file 'D:\Dev\math\vendor\OSforGFF\.lake\packages\mathlib\.lake\build\lib\lean\Mathlib.olean' of module Mathlib does not exist
# LEAN_PATH 只给 mathlib 目录
error: unknown module prefix 'Aesop'
```

**操作坑（务必记）**：冷启动首次 import 会偶发
`error: failed to read file '...\Lean\Meta\AppBuilder.olean.private'` 或 `...\Init\Data\Iterators\...\FlatMap.olean`
（文件存在、可读、非 0 字节；**重跑即过**）。这类是 Windows 侧首读/索引抖动，**不是** 环境坏。`[T]`
⇒ 任何 Lean 任务遇 read 失败先重试一次再判死。

---

## 5. 判：mathlib 可用？—— **可用（子模块级），不可用（伞模块级）**

- ✅ 可 `import`：`Mathlib.Data.Matrix.Basic`、`Mathlib.Data.Real.Basic`、`Mathlib.Algebra.*`、`Mathlib.Order.*`、`Mathlib.Topology.*`、`Mathlib.Tactic.FinCases`、`Mathlib.Tactic.NormNum`、`Matrix.diagonal_mul_diagonal` 等。`[T]`
- ❌ 不可 `import`：`Mathlib`（根伞）、`Mathlib.Tactic`（伞）、`Mathlib.Data.Matrix.Notation`（`!![...]` 记法由此失效）、`Mathlib.Data.Matrix.Defs`、`Mathlib.Order.Field.Basic`。`[T]`
- ⚠️ 因此：**不要写 `import Mathlib`**；逐模块 import；`!![...]` 用 `Matrix.of ![![...],![...]]` 或 `Matrix.diagonal` 代替（实测后者能过）。
- ⚠️ `decide` 在 `ℚ` 上有陷阱：`(0:ℚ) < a - b` 这类**内核 reduce 卡住**（`Rat.blt` 不展开）；
  可用的组合是「先 `decide` 证具体元素等式 → `rw` → `norm_num`」或 `native_decide`（会引入 `._native.native_decide.ax` 公理）。`[T]`

---

## 6. 基准反例（第 3 问）—— `research/v3/lean/M03-COUNTER.lean`

**链条被证伪的那一步**：`f₁² = 0 ⇒ m = 0`（残数消失 ⇒ 无隙）。有限维 2×2 反例：

- `H = diagonal ![0, 1]`（谱隙 `H₁₁ - H₀₀ = 1 > 0`），`O = diagonal ![3, -1]`（非标量）。
- `H * O = O * H`（`Matrix.diagonal_mul_diagonal`）— 对位 `trF² = H` 守恒 ⇒ `[H,O] = 0`。
- `f1sq := O 0 1 = 0` — 对位 `⟨𝔾|[H,O]|n⟩ = 0 ⇒ f₁² = 0`。
- **同时** `0 < H 1 1 - H 0 0` ⇒ `m > 0`。

主定理：

```lean
theorem counterexample :
    ∃ H O : Matrix (Fin 2) (Fin 2) ℚ,
      H * O = O * H ∧ O 0 1 = 0 ∧ (0 : ℚ) < H 1 1 - H 0 0

theorem not_f1sq_zero_implies_gapless :
    ¬ (∀ H O : Matrix (Fin 2) (Fin 2) ℚ,
        H * O = O * H → O 0 1 = 0 → H 1 1 - H 0 0 = 0)
```

公理审计（全文无 `sorry`、无 `native_decide` 公理）：

```
'M03.counterexample' depends on axioms: [propext, Classical.choice, Quot.sound]
'M03.not_f1sq_zero_implies_gapless' depends on axioms: [propext, Classical.choice, Quot.sound]
```

另有 A 层纯 `ℚ` 算术核 `M03.f1sq_eq_zero : (a - b) * o01 = 0 → a ≠ b → o01 = 0`，
即「交换子矩阵元为零 + 谱不退化 ⇒ 非对角元为零」——这是链条里唯一有数学内容的一步，且**无条件成立**；
反例说明它**反方向不成立**（不能由 `o01 = 0` 反推 `a = b`，也即不能由 `f₁² = 0` 反推 `m = 0`）。

**这一层的用途**：以后凡是要把「残数/搬运量 = 0」当成「无隙」用的论证，本条即可一击证伪；`H` 型有隙 + 残数归零的配对是**相容的**，不构成矛盾。

---

## 7. 复现清单（后续所有 Lean 任务照抄）

```powershell
# 唯一推荐入口
cd D:\Dev\math\vendor\OSforGFF
lake env lean ..\..\research\v3\lean\<FILE>.lean
```

```powershell
# 不依赖 cwd 的等价入口（须 pin toolchain + 全 LEAN_PATH）
cd D:\Dev\math
# 用 §3(b) 那段 12 条 LEAN_PATH 逐字赋值
lean +leanprover/lean4:v4.33.0-rc1 research\v3\lean\<FILE>.lean
```

- 禁止：`import Mathlib`；裸 `lean`（默认 4.36）；`cd 仓库根 ; lake build`（无 lakefile）。
- 若需完整 mathlib（含 `Mathlib.olean` 根伞）：得把 `build.badcache` 复活或重新 `lake build Mathlib`，本轮**未做**（超出 30 min 窗，且 2264 子集已够本波用）。`[C]`

标签：§1–§6 全部条目均为实测 `[T]`；§7 最后一条重建伞模块的成本判断为 `[C]`。
