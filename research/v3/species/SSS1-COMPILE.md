# SSS1-COMPILE —— OS2 四文件真编译裁决

轮次：第 1 轮（本会话）。执行者：`sss1-compile`。
对象：`research/v3/species/ooo2/{OS2Clip,OS2Topology,OS2Theta,OS2Wiring}.lean`。
**不 commit。零新 MC（无一处新未机检断言）。**

---

## §0 判词

**`[T]` —— 四个文件全部 0 error，`lake build` 整包通过（2909 jobs）。**

`ooo2` 自报「四文件机检 `[NG]`，唯一缺口 = 跑一次 `lean`」**只对了一半**：
「从未跑过」为真，但「只差跑一次」为假 —— 每份文件都有**真错误**，共 4 类 11 处：

| 类 | 实例 | 性质 |
|---|---|---|
| 缺 import | `abs_le`、`Mathlib.Topology.Metrizable.Basic` | 缺引理（准确说是引入缺失） |
| 常量改名 | `le_or_lt` → `le_or_gt`（mathlib4 v4.33 已更名） | 未定义常量 |
| 定义顺序 | `OS2Clip` 第 85 行先用 `clip_eq_self`，该引理第 97 行才定义 | 未定义常量（前向引用） |
| 命名空间 / 宇宙 / 依赖 | `PseudoMetrizableSpace` 实际在 `TopologicalSpace` 名下；`LoopSpace (κ : Type*) : Type` 宇宙不匹配；`GL4Transfer` 不在搜索路径 | 类型不符 / 未定义模块 |

外加**一处数学错误**（唯一实质错误，非机械错误）：
`OS2Topology` 原写「任意 `κ`，`κ → [−1,1]` 是 `PseudoMetrizableSpace`」——**假**。
不可数积不拟可度量化。已加 `[Countable κ]`（同文件 `SecondCountableTopology` / `BorelSpace` /
`OpensMeasurableSpace` 三行本来就带该假设，故这是原稿漏写而非结论削弱）。
`κ` 任意时仍成立的：`TopologicalSpace`、`CompactSpace`（Tychonoff）、`T2Space`——三行保持无条件。

**公理审计**（`#print axioms`，落盘 `lean/os2/axcheck.lean`）：

```
rp_full_cone              : [propext, Classical.choice, Quot.sound]
measurePreserving_of_weak_limit : [propext, Classical.choice, Quot.sound]
map_eq_of_weak_limit      : [propext, Classical.choice, Quot.sound]
clip_memLp                : [propext, Classical.choice, Quot.sound]
```

`sorry` / `admit` 扫描：四文件 **0 命中**。无 `sorryAx` ⇒ 主定理 `rp_full_cone` 是**真机检定理**，
不是带洞骨架。

**下一手**：
1. `docs/position-current.html` 需按 AGENTS §6 记一次等级变更（OS2 四文件 `[NG]` → `[T]`）——本轮未动，
   留给 lead 合并（避免与其它代理的写冲突）。
2. `OS2Wiring` 的 `rp_full_cone` 仍有两条**显式假设** `hRP`（μ_a 侧 clip 层 RP）与 `hA2 : MemLp A 2 μ`；
   下一手是把 `hRP` 从签名降为推论（或证明它是 RP 在 `C_b` 读法下的免费项），而不是继续加假设。
3. `GL4Transfer` / `GL4L2` 本轮顺手机检也过了（见 §2.5），`research/glue-a-limits/check-glue-a.sh`
   的 `LEAN_PATH` 可以撤掉 `build.badcache` 那一段（原脚本指向不存在目录）。

---

## §1 lakefile / toolchain

**复用 `rrr1-lake` 已建好的 `lean/os1`**（不重建、不重拉 mathlib）：

```
lean/os1/lakefile.lean    package «os1», require mathlib @ v4.33.0-rc1
lean/os1/lean-toolchain   leanprover/lean4:v4.33.0-rc1
lean/os1/.lake/packages   JUNCTION → D:\Dev\math\vendor\OSforGFF\.lake\packages
```

工具链与缓存实测：

```
lean --version  → Lean (version 4.33.0-rc1, x86_64-w64-windows-gnu, commit 62eed1db…)
lake --version  → Lake version 5.0.0-src+62eed1d (Lean version 4.33.0-rc1)
mathlib olean   → 8279 个 / 源码 8268 个 .lean（缓存完整，无缺口）
junction 实测   → packages [D:\Dev\math\vendor\OSforGFF\.lake\packages]
```

新建 `lean/os2`（同一 junction，不下载任何东西）：

```powershell
New-Item -ItemType Directory -Force -Path lean\os2
Copy-Item lean\os1\lean-toolchain, lean\os1\lake-manifest.json lean\os2\ -Force
New-Item -ItemType Directory -Force -Path lean\os2\.lake
cmd /c mklink /J lean\os2\.lake\packages D:\Dev\math\vendor\OSforGFF\.lake\packages
Copy-Item research\v3\species\ooo2\OS2*.lean lean\os2\ -Force
```

`lean/os2/lakefile.lean`（本文件即最终版）：

```lean
import Lake
open Lake DSL

package «os2» where
  leanOptions := #[⟨`autoImplicit, false⟩]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "v4.33.0-rc1"

@[default_target]
lean_lib «OS2» where
  roots := #[`OS2Clip, `OS2Topology, `OS2Theta, `OS2Wiring]

/-- GL4 依赖文件（OS2Wiring 的 `rp_transfer` / `rp_l2_closure` 出处），随包一起编译。 -/
lean_lib «GL4» where
  roots := #[`GL4Transfer, `GL4L2]
```

编译命令（逐个 = `lake build <module>`，等价于对单文件跑 `lean`，但带完整 `LEAN_PATH`）：

```powershell
cd D:\Dev\math\lean\os2
lake build OS2Clip; lake build OS2Topology; lake build OS2Theta; lake build OS2Wiring
lake build            # 整包，末行 Build completed successfully (2909 jobs)
```

逐文件复核（直接 `lake env lean`，看的是源码级 0 error，不看 lake 缓存）：

```powershell
foreach ($f in "OS2Clip","OS2Topology","OS2Theta","OS2Wiring") { lake env lean "$f.lean" }
# → OS2Clip : OK (0 error) / OS2Topology : OK / OS2Theta : OK / OS2Wiring : OK
```

---

## §2 逐个编译结果

### 2.1 `OS2Clip.lean` —— `[T]`（修复前 `[NG]`）

修复前第一个错误（原始状态，`lake build OS2Clip`）：

```
error: OS2Clip.lean:61:11: Unknown identifier `le_or_lt`
error: OS2Clip.lean:61:37: Tactic `rcases` failed: `x✝ : ?m.26` is not an inductive datatype
error: OS2Clip.lean:70:13: Unknown identifier `le_or_lt`
error: OS2Clip.lean:85:10: Unknown identifier `clip_eq_self`
```

三类错：
1. **常量改名**：mathlib4 v4.33 中 `le_or_lt` → `le_or_gt`（`@le_or_gt : a ≤ b ∨ b < a` 实测在，
   `@le_or_lt` 实测 `Unknown identifier`；`lt_or_le` → `lt_or_ge` 同理）。
2. **前向引用**：第 85 行用 `clip_eq_self`，但该引理在第 97 行才定义 —— 文件顺序错，非引理缺。
3. **缺 import**：`abs_le` 在 `Mathlib.Algebra.Order.Ring.Abs`，原 import 段没有它
   （实测：原 import 下 `#check @abs_le` → `Unknown identifier`；加该 import 后
   `@abs_le : |a| ≤ b ↔ -b ≤ a ∧ a ≤ b`）。

修复（3 处，均为机械修复，不改陈述）：
- import 段加 `import Mathlib.Algebra.Order.Ring.Abs`；
- `le_or_lt` → `le_or_gt`（61、70 两处）；
- `clip_eq_self` 整块前移到 `clip_tendsto` 之前。

结果：`Build completed successfully (2887 jobs)`；`lake env lean OS2Clip.lean` → 0 error。
证明内容（`continuous_clip` / `clip_abs_le` / `clip_tendsto` / `clip_abs_le_n` / `clip_eq_self` /
`clipBCF` / `clip_memLp`）**逐字未改**。

### 2.2 `OS2Topology.lean` —— `[T]`（修复前 `[NG]`）

修复前第一个错误：

```
error: OS2Topology.lean:52:42: Unknown identifier `PseudoMetrizableSpace`
  Note: … `autoImplicit` is set to `false`
error: OS2Topology.lean:69:39: Type mismatch
  κ → LoopCoord has type Type u_1 … but is expected to have type Type
```

四类错：
1. **命名空间**：`class PseudoMetrizableSpace` 定义在 `namespace TopologicalSpace` 内
   （`Mathlib/Topology/Metrizable/Basic.lean:40`），根名不解析 ⇒ 全文件 6 处改
   `TopologicalSpace.PseudoMetrizableSpace`。
2. **宇宙不匹配（类型不符）**：`abbrev LoopSpace (κ : Type*) : Type` —— `κ → LoopCoord` 是
   `Type u_1` 却声明成 `Type`。改为 `universe u` + `(κ : Type u) : Type u`。
3. **缺 import**：`PseudoMetrizableSpace ℝ` / 子类型实例需要 `Mathlib.Topology.Metrizable.Real`；
   注意本 mathlib **没有** `Mathlib.Topology.Instances.Real`（`import` 实测报
   `object file … Instances/Real.olean does not exist`，源码也不存在）——原依赖清单里那条已失效。
4. **数学错（唯一实质错）**：`example (κ : Type*) : PseudoMetrizableSpace (LoopSpace κ)`
   无条件为**假**：不可数积非拟度量化。实测 `import Mathlib` 全量下同样
   `failed to synthesize PseudoMetrizableSpace (κ → ↑(Set.Icc (-1) 1))`。
   加 `[Countable κ]` 后实测四个实例全过（`PseudoMetrizableSpace` / `SecondCountableTopology` /
   `BorelSpace` / `CompactSpace (ProbabilityMeasure …)`）。
   `CompactSpace`（Tychonoff）、`T2Space`、`TopologicalSpace` 三行**保持任意 `κ` 无条件**。

结果：`Build completed successfully (2604 jobs)`；`lake env lean` → 0 error。

### 2.3 `OS2Theta.lean` —— `[T]`（修复前 `[NG]`）

修复前第一个错误（连锁型：变量段第一行挂 ⇒ 后面 40+ 条 `Ω` 全变未定义）：

```
error: OS2Theta.lean:53:3: Unknown identifier `PseudoMetrizableSpace`
error: OS2Theta.lean:59:41: Unknown identifier `Ω`
…
error: OS2Theta.lean:93:37: cannot coerce to function θe
error: OS2Theta.lean:99:18: Unknown identifier `map_eq_of_weak_limit`
```

根因只有一条：`variable … [PseudoMetrizableSpace Ω]` 不解析 ⇒ `Ω` 根本没进上下文 ⇒ 后面每条都是
级联噪声（`⇑θe`、`map_eq_of_weak_limit` 均为假警报）。

修复（2 处）：import 段加 `Mathlib.Topology.Metrizable.Basic`；类名限定为
`TopologicalSpace.PseudoMetrizableSpace Ω`。

结果：`Build completed successfully (2546 jobs)`；`lake env lean` → 0 error。
`map_eq_of_weak_limit` / `measurePreserving_of_weak_limit` **证明未被改动**，仅公理
`propext, Classical.choice, Quot.sound`。

### 2.4 `OS2Wiring.lean` —— `[T]`（修复前 `[NG]`：依赖模块缺失）

修复前第一个错误：

```
error: OS2Wiring.lean:26:0: unknown module prefix 'GL4Transfer'
  No directory 'GL4Transfer' or file 'GL4Transfer.olean' in the search path entries: …
```

**未定义模块**，不是证明错。处理：把 `research/glue-a-limits/{GL4Transfer,GL4L2}.lean` 拷入
`lean/os2/`（保持模块名），并在 lakefile 里加 `lean_lib «GL4»`。
这两个依赖本身**先编过了**（见 2.5），于是 `import GL4Transfer` / `import GL4L2` 成立。

另修 2 处：`import Mathlib.Topology.Metrizable.Basic`（供类名解析）；
`[PseudoMetrizableSpace Ω]` → `[TopologicalSpace.PseudoMetrizableSpace Ω]`。

结果：`lake env lean OS2Wiring.lean` → 0 error；`rp_full_cone` 无 `sorryAx`，公理三条标准公理。
**主定理陈述与证明脚本逐字未改**（`hRP`、`hA2` 两条显式假设仍在签名里，见 §0 下一手 2）。

### 2.5 顺手核（依赖侧，非本任务四文件）

| 文件 | 结果 | 依据 |
|---|---|---|
| `research/glue-a-limits/GL4Transfer.lean` | `[T]` | `Build completed successfully (2588 jobs)` |
| `research/glue-a-limits/GL4L2.lean` | `[T]` | `Build completed successfully (2869 jobs)` |

即：`ooo2` 交接链上游不是断的 —— 断的只是「没进 `LEAN_PATH`」。

---

## §3 可编译子集

**`{OS2Clip, OS2Topology, OS2Theta, OS2Wiring} = 4/4 全在可编译子集内`**（加依赖
`{GL4Transfer, GL4L2}` 共 6/6）。

产物落盘：

```
lean/os2/lakefile.lean        新建（package «os2» + lean_lib «OS2»/«GL4»）
lean/os2/lean-toolchain       leanprover/lean4:v4.33.0-rc1
lean/os2/lake-manifest.json   复用 os1
lean/os2/.lake/packages       JUNCTION → vendor/OSforGFF/.lake/packages
lean/os2/{OS2Clip,OS2Topology,OS2Theta,OS2Wiring}.lean   已修，0 error
lean/os2/{GL4Transfer,GL4L2}.lean                        依赖，0 error
lean/os2/.lake/build/lib/lean/OS2*.olean                 OS2Clip 90264 / OS2Topology 40896
                                                         OS2Theta 41096 / OS2Wiring 45296
```

修好的四文件已回写 `research/v3/species/ooo2/`（`SHA256` 与 `lean/os2/` 逐一相符，
`OS2Wiring.lean` 实测同为 `16BA9E47…52EB0A`）。

**判**：0 error = 4/4；有错 = 0/4。原 `[NG]` 的错误分布：
缺引理 2 处（`abs_le`、`Metrizable.Real`）、未定义常量 4 处（`le_or_lt`×2、`clip_eq_self`、
`PseudoMetrizableSpace`）、类型不符 2 处（宇宙、`⇑θe` 级联）、未定义模块 1 处（`GL4Transfer`）、
数学错 1 处（不可数积拟度量化，已加 `[Countable κ]`）。
