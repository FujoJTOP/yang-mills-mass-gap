# RRR1-LAKE · OS1 最小子集的 lakefile / toolchain 建立与逐件编译

轮次：本轮（第 1 轮 lake 化）。仓库：`D:\Dev\math`。**未 commit。零新 MC。**

---

## §0 判词

| 项 | 裁决 | 证据 |
|---|---|---|
| `lean/os1` 是否有可复现编译目标 | **`[T]`** — 有 | `lake build OS1Tight` → `Build completed successfully (2604 jobs)`，exit 0，产物 `lean/os1/.lake/build/lib/lean/OS1Tight.olean`（46920 B） |
| mathlib 是否复用 vendor（不重拉） | **`[T]`** — 复用 | `lean/os1/.lake/packages` 为 NTFS junction → `vendor/OSforGFF/.lake/packages`；mathlib `rev-parse` = `79d0395a18`，与 `vendor/OSforGFF/lake-manifest.json` 的 `79d0395a1825a6264ad5d269e35e60537518955e` / `inputRev v4.33.0-rc1` 逐字符一致 |
| mathlib olean 是否齐 | **`[T]`** — 齐 | `lake exe cache get`：`Decompressing 6005 already-cached file(s) (2638 already decompressed)` / `No files to download` / `Completed successfully in 31724 ms`；`Mathlib` 下 `.olean` 计数 **8268**（= 源码 8268） |
| 四件全编译 | **`[NG]`** — 3 件不过 | 见 §2；首个错误逐件给 |

**裁决一句话**：**tower 建成了，最小可编译子集 = `{OS1Tight}` 单件（0 sorry）**，其余三件的首错都不是「环境缺件」而是**源文件级首错**（scoped notation 未开 / 未知类名 / 未标 noncomputable）。

**下一手**：按 §2 的首错行，给三件各补 1–2 行 header（`open scoped BoundedContinuousFunction` + `import Mathlib.Topology.Metrizable.Basic` + `noncomputable`），重编一次；`OS1JointShift:100` 与 `OS1Schwartz:95` 的 `sorry` 是内容缺件，不属 header 修复范围。

---

## §1 lakefile / toolchain

### 命令（按序，全部已执行）

```powershell
cd D:\Dev\math\lean\os1
New-Item -ItemType Directory -Force .lake
cmd /c mklink /J ".lake\packages" "D:\Dev\math\vendor\OSforGFF\.lake\packages"   # 复用，不重拉
Copy-Item ..\..\vendor\OSforGFF\lake-manifest.json .\lake-manifest.json -Force     # 锁 rev
Set-Content lean-toolchain "leanprover/lean4:v4.33.0-rc1"
lake exe cache get                                                                # 解压已有 ltar，不联网下载
lake build OS1Tight
```

### 结果

- `lake-toolchain`：`leanprover/lean4:v4.33.0-rc1`（`elan toolchain list` 命中：`leanprover/lean4:v4.33.0-rc1` + `v4.36.0`，已装，无需下载）。
- `lake --version` = `Lake 5.0.0-src+611c117 (Lean version 4.36.0)`（宿主默认；项目内由 `lean-toolchain` 切到 4.33.0-rc1）。
- `lake exe cache get`：`No files to download` —— vendor 已存 `.ltar`，本步只解压，**0 网络下载**。
- 失败/异常：无（首个 turn 曾把 `autoImplicit := false` 写进 lakefile，造成 3 件假错；已删，见 §2 备注）。

### 文件

`lean/os1/lakefile.lean`

```lean
import Lake
open Lake DSL

package «os1»

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "v4.33.0-rc1"

@[default_target]
lean_lib «OS1» where
  roots := #[`OS1Tight, `OS1JointShift, `OS1Statement, `OS1Schwartz]
```

`lean/os1/lean-toolchain` → `leanprover/lean4:v4.33.0-rc1`

---

## §2 逐个编译（`lake build <module>`）

| 文件 | 退出码 | 裁决 | 第一个错误（原文） |
|---|---|---|---|
| `OS1Tight.lean` | 0 | **`[T]`** | 无。`✔ [2604/2604] Built OS1Tight (24s)` / `Build completed successfully (2604 jobs)` |
| `OS1JointShift.lean` | 1 | **`[NG]`** | `OS1JointShift.lean:53:38: invalid binder annotation, type is not a class instance / ?m.7` —— `PseudoMetrizableSpace` 未在作用域（类定义在 `Mathlib/Topology/Metrizable/Basic.lean:40`，本文件未 import）。其后 `61:41 expected token` ⇒ `→ᵇ` 未开 scoped notation |
| `OS1Statement.lean` | 1 | **`[NG]`** | `OS1Statement.lean:38:18: expected token` —— `(f : E →ᵇ ℝ)` 的 `→ᵇ` 需 `open scoped BoundedContinuousFunction`；本文件只 `open scoped Topology ENNReal` |
| `OS1Schwartz.lean` | 1 | **`[NG]`** | `OS1Schwartz.lean:42:0: failed to compile definition, consider marking it as 'noncomputable' because it depends on 'Complex.instNormedField', which is 'noncomputable'`（`example : AddCommGroup S' := inferInstance`）。另 `102:14` `ProbabilityMeasure` 未知（缺 import） |

**备注（必读）**：这三个首错是**真源码首错**，不是 toolchain 假错。依据：删掉 `autoImplicit := false` 前后重编，首错行/文案不变（`53:38 / 38:18 / 42:0`），只是错误类型从 `Unknown identifier` 收敛为 `invalid binder annotation` / `expected token` —— 说明 `autoImplicit=false` 是我的 lakefile 注入的假象，已移除。

---

## §3 可编译子集（0 sorry）

### 编译 × sorry 二维表

| 文件 | 编译 | 真 `sorry` 处 | 0-sorry | 备注 |
|---|---|---|---|---|
| `OS1Tight.lean` | `[T]` | 0（第 7 行 `sorry` 仅出现在文件头注释「本文件不含 `sorry`」） | ✅ | **唯一入子集** |
| `OS1JointShift.lean` | `[NG]` | 1（`100: sorry -- [GLUE-1] uniform_prod_eq_comap_prod + mem_prod_iff`，第 9 行为注释） | ❌ | ②b 带 sorry 不入 |
| `OS1Statement.lean` | `[NG]` | 0（第 4 行为注释） | （编译不过） | 补 header 后可望入；已探针确认**不止 header** |
| `OS1Schwartz.lean` | `[NG]` | 1（`95: sorry`） | ❌ | 另有 `42/44/83` noncomputable 与 `102` 缺 import |

### 裁决

**最小可编译子集 = `{OS1Tight}`，0 sorry，`[T]`。**

### 探针（不是修复，只用来判「补 header 够不够」）

`lean/os1/probe/OS1Statement.lean` = 原文件 + 两行 header（`import Mathlib.Topology.ContinuousMap.Bounded.Basic`、`open scoped BoundedContinuousFunction`），`lake env lean probe/OS1Statement.lean` → exit 1：

```
probe\OS1Statement.lean:44:48: error: Application type mismatch: The argument h has type
  ∀ (f : E →ᵇ ℝ), ∫ x, f x ∂↑μ = ∫ x, f x ∂↑ν
but is expected to have type ∀ (f : ?m.32 →ᵇ ℝ), ...
in the application FiniteMeasure.ext_of_forall_integral_eq h
probe\OS1Statement.lean:60:40: error: invalid binder annotation, type is not a class instance / ?m.7
```

⇒ `OS1Statement` **不是「差一行 header」**：还有 `ext_of_forall_integral_eq` 的实例合成失配（`E` 的 `PseudoMetrizableSpace`/`BorelSpace` 未被 `inferInstance` 落下）与第二个未知类名。补 header 后仍 `[NG]`。

### 下一手

1. `OS1Statement`：补 `import Mathlib.Topology.Metrizable.Basic` + `open scoped BoundedContinuousFunction`，再修 `44:48`（给 `FiniteMeasure.ext_of_forall_integral_eq` 显式喂实例，或改用 `Measure.ext_of_forall_integral_eq`）。
2. `OS1JointShift`：同一套 header；`100` 的 `[GLUE-1]` 是 4 行 `uniformity_prod_eq_comap_prod` 搬运，属内容缺件而非环境。
3. `OS1Schwartz`：`example` 三行加 `noncomputable`，`102` 补 `import Mathlib.MeasureTheory.Measure.ProbabilityMeasure`；`95` 的 `continuous_shift_pair` 仍是内容缺件。
4. 子集达标后，把 `lean/os1` 的 `default_target` 收缩为 `roots := #[`OS1Tight]`，让 `lake build` 默认绿。

---

## §4 工具执行核验（AGENTS §1）

| 检索式 | 工具 | 状态 | 命中 |
|---|---|---|---|
| `Get-ChildItem lean\os1` | pwsh | OK | 4 个 `.lean`（4277/6250/5128/4619 B） |
| `lake --version` | pwsh | OK | `Lake 5.0.0-src+611c117 (Lean 4.36.0)` |
| `elan toolchain list` | pwsh | OK | `v4.33.0-rc1`、`v4.36.0` |
| `Test-Path .lake\packages\mathlib\.lake\build\lib\lean\Mathlib\MeasureTheory\Measure\Prokhorov.olean` 等 8 条 | pwsh | OK（首次全 False 因路径少了 `lib\lean`，已更正重测） | cache 后 8/8 True |
| `Get-ChildItem Mathlib -Recurse -Filter *.olean \| Measure` | pwsh | OK | 8268 |
| `Select-String -Pattern "class PseudoMetrizableSpace"`（mathlib 全树） | pwsh | OK | 1（`Mathlib/Topology/Metrizable/Basic.lean:40`） |
| `Select-String -Pattern '\bsorry\b'` × 4 文件 | pwsh | OK | Tight 1（注释）/ JointShift 2（1 注释 1 真）/ Statement 1（注释）/ Schwartz 2（1 注释 1 真） |
| `lake build <module>` × 4 | pwsh | OK | 见 §2 |

TOOL-FAIL：无。
