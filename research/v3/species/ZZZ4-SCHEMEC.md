# ZZZ4-SCHEMEC — 钉死「方案 C」：它是**重正化方案**（RT 上以有限常数 `C` 标号的 scheme），不是仓内 A/B/C 编号；它的 `[T]` 只到定义/逐阶，不到值

> **轮** ZZZ4 · 子代理 `zzz4-schemec` · `D:\Dev\math`。只新建本文件；既有文件一字未改；未 `commit`；**零新 MC**；未动 `docs/position-current.html`（`AGENTS §6` 交 `lead`）。
> **靶**：`YYY1-U123` 判「缺 `Z_O(k)` 沿方案 C 的 `k`-极限」⇒ 先钉方案 C 的所指，再判仓内有无它的 `[T]`。
> **纪律**：只 `[T]`／`[NG]`；无百分数；每条带「下一手」；时间单位是轮。

## §0 判词（先行）

| # | 断言 | 判 | 依据 `file:line` | 下一手 |
|---|---|---|---|---|
| **0.1** | **方案 C = 重正化方案**，不是仓内三方案编号之一 | `[T]` | `research/unconditional/A2-scheme-match.md:1,5–6`·`research/gluing/REPORT.md:106,112,180`；反证 `方案 A`／`方案 B` 各 0 命中（§1 第 7 行） | 引用统一写 `scheme C`／方案 C，禁写"方案三" |
| **0.2** | **定义**：`C` 是 RT 上的一个 scheme 常数，`1/g²(a)=2b₀log(1/aΛ)+(b₁/b₀)loglog(1/aΛ)+C+o(1)`；不同 `C` = 不同 RT，九个接口必须同 `C` | `[T]` | `A2-scheme-match.md:6`·`gluing/REPORT.md:106,112`·`research/glue-a-limits/REPORT.md:90`（`(U2)` 逐字「同一 scheme 常数 `C`」） | 引 `(U2)` 必带常数 `C`，不带"某个 scheme" |
| **0.3** | **它挂在这一族的 scheme 依赖上**：换 scheme `g'²=g²+cg⁴` 保 `b₀`、平移 `Λ'=Λe^{c/2b₀}`，`m/Λ,√σ/Λ,q=a/ξ` 不变 ⇒ 方案 C 是**成员**，不变量不由它产生 | `[T]` | `A2-scheme-match.md:6`·`gluing/REPORT.md:114`；`b₀` 归一/方案依赖见 `research/v3/species/YY5-ANOM.md:0.4`·`research/v3/species/VVV3-CRIT.md:14,43` | 凡写"方案 C 钉死 X"，先核 X 是否 scheme-不变 |
| **0.4** | **仓内有方案 C 的 `[T]`，只在定义层与逐阶匹配层**：L1（改 scheme = 有限重参数化，definitional）+ L2（`Z_a`↔scheme C 全阶微扰匹配，定理） | `[T]`（L1+L2） | `research/unconditional/G1-ASSEMBLE-100.md:69`·`research/unconditional/G1-UNIFORM-ZA.md:6,22`·`research/v3/species/HQ1-L3.md:100`·`research/v3/species/HS-L3-6.md:14`·`research/v3/species/HZ6-ZDEF.md:0.1,45` | `L3` 台账把 `[T]` 限定写死成「L1+L2」，不加词 |
| **0.5** | **方案 C 的"值"无 `[T]`**：`Z_O(k)` 沿方案 C 的 `k`-极限 + 与 `C` 的有限常数一致 = `(Z-last)`，判 `[NG]` | `[NG]` | `research/v3/species/HZ9-FINAL.md:35`（行 10）·`HZ6-ZDEF.md:16,60`（D4）·`research/v3/species/HZ1-UNBOUND.md:65`·`research/v3/species/YYY1-U123.md:37` | 只攻 `ϖ(k)=Θ(1)` 双侧 + `lim_k ϖ(k_*)/ϖ(k)`；不再问"方案 C 是什么" |

## §1 grep（`AGENTS §1` 四要素）

| # | 实际检索式 | 工具 | 状态 | 命中 | 承重 `file:line` |
|---|---|---|---|---|---|
| 1 | `方案\s*C` | `grep` | **OK**（无 provider failure／ENOENT／空响应） | **44** | `HQ1-L3.md:100,184`·`HS-L3-6.md:14`·`HZ9-FINAL.md:35,55,58`·`YYY1-U123.md:4,37` |
| 2 | `方案C`（无空格） | `grep` | **OK** | **0** | — |
| 3 | `scheme C`（小写） | `grep` | **OK** | **28** | `A2-scheme-match.md:1,3,5,6,14`·`G1-EXISTENCE-THEOREM.md:3,7,12`·`L3-scale-uniform.md:3,6` |
| 4 | `Scheme C`（首字母大写） | `grep` | **OK** | **0** | — |
| 5 | `scheme-C` | `grep` | **OK** | **5** | `A2-scheme-match.md:1`·`CONTINUUM-TIGHTNESS.md:3`·`L3-scale-uniform.md:3,6`·`PROOF-G1-BPHZ.md:5`（`scheme_C` **0**） |
| 6 | `\(C\)` | `grep` | **OK** | **729**（显示 250，全表落盘） | **无一处指方案 C**：`research/v3/G4-GENERAL-BOUND.md:74`（对象 (C)=粗粒曲率）、`docs/LEAN-STATUS.md:39`（(C) 归约定理）、`research/confinement/PROOF.md:34`（Wilson 环 `C`） |
| 7 | `方案\s*A` ／ `方案\s*B` | `grep` | **OK** | **0** ／ **0** | ⇒ 无 A/B/C 编号三件套 |

`TOOL-FAIL` **0**；空响应 **0**；`PATH-MISS` **0**。`pwsh`／`web_*`／MC／Lean／`git`／`skill` **0** 调用（本靶是仓内命名裁决）。

## §2 方案 C 是什么（一行 + 三条支撑）

**一行**：`scheme C` = 沿二维重正化轨迹（RT）以有限常数 `C` 标号的那一族重正化方案中的**一个成员**（`A2-scheme-match.md:5` 逐字 "L1 — scheme C is a choice, not an obstruction"）：同一 `C` 的 RT 才进同一套九接口（`gluing/REPORT.md:112`「不同 `C` 是不同的 RT」）。

1. **不是仓内编号**：`方案 A`／`方案 B` 各 0 命中；`C` 在方程里是有限常数（`A2:6`），并在符号表登记为「RT 的 scheme 常数，不是 `C_RP`」（`gluing/REPORT.md:180`）。
2. **它是 scheme 名，与别的 scheme 并列**：`HS-L3-8.md:60` 逐字「B4 的 counterterm 递推是 **Balaban 方案，不是方案 C**」⇒ 方案 C 与 Balaban 方案同级，属 scheme 命名，不是"第三种策略"。
3. **与 β 的关系**：换 scheme 保 `b₀`（`A2:6`·`gluing:114`），但 `b₀=11N/48π²` 的**归一**带方案/耦合约定，唯一不变量是比值 `γ₀/b₀=2`（`YY5-ANOM.md:0.4`·`VVV3-CRIT.md:14,43`）；三圈系数亦 scheme 依赖（`gluing/REPORT.md:681`）。⇒ `YYY1-U123` 的缺件不在 scheme 层，在**值层**。
4. **同名警告（必带）**：`Z_O` 是三个同名量（`HZ9-FINAL.md:17`·`HZ6-ZDEF.md:1.2`·`docs/CLAY-PROOF-CHAIN.html:262`）；本靶的 `Z_O(k)` 是锚定方差比 `[Z_a(k_*)/Z_a(k)]^{1/2}`（`HZ6-ZDEF.md:41`），**不是**迹反常比 `(β/2g)(μ₀)/(β/2g)(μ)`。

## §3 仓内有没有方案 C 的 `[T]`

| 层 | 内容 | 判 | `file:line` |
|---|---|---|---|
| L1 定义 | 改 scheme = 有限重参数化 `g'²=g²+cg⁴`，`Λ'=Λe^{c/2b₀}`，不变量 `m/Λ,√σ/Λ,a/ξ`；"Bookkeeping, not an analytic wall" | **`[T]`** | `A2-scheme-match.md:5–6`·`G1-EXISTENCE-THEOREM.md:12`·`G1-ASSEMBLE-100.md:69`·`gluing/REPORT.md:114` |
| `(U2)` 措辞 | 「同一 scheme 常数 `C`」下的 `a^{-4}` 减除；RT 进入 GL3 的唯一位置 | **`[T]`**（逐字） | `glue-a-limits/REPORT.md:90`·`HS-L3-10.md:40`·`HW1-CLOSURE.md:18` |
| L2 匹配 | `Z_a`↔scheme C：全阶微扰 lattice↔RT（Reisz BPHZ + Symanzik），dim-4 类闭于 `{TrF²,TrFF̃}`（`θ` 全导数） | **`[T]`** | `A2:9`·`G1-UNIFORM-ZA.md:6,9,22`·`HQ1-L3.md:100,184`·`HS-L3-6.md:14`·`L3-scale-uniform.md:3,6` |
| `Z_O(k)` 定义层 | `Z_O(k):=[Z_a(k_*)/Z_a(k)]^{1/2}`，裸耦合"由方案 C 定"；存在/正性/幂=记账（D1–D3） | **`[T]`** | `HZ6-ZDEF.md:0.1,12,41–45,56–58` |
| `(Z-last)` 值层 | `Z_O(k)` 沿方案 C 的 `k`-极限存在且等于 `C` 钉死的有限常数；等价 `ϖ(k)=Θ(1)` 双侧 `(a,k)`-一致 | **`[NG]`** | `HZ9-FINAL.md:35,55–63`·`HZ6-ZDEF.md:16,60`·`HZ1-UNBOUND.md:65,72`·`YYY1-U123.md:37` |

**裁决**：**有 `[T]`**（4 行：L1／`(U2)` 措辞／L2／定义层），**值层无 `[T]`**。方案 C 本身已钉死到"可选但必须声明、且九接口一致"；`YYY1-U123` 的「缺 `Z_O(k)` 沿方案 C 的 `k`-极限」因此**不是方案未定**，而是**值未定**——该判与本文件不冲突，也不被本文件升级。**下一手**：把 §3 前 4 行抄进 `L3` 台账并加「`[T]` 只到 L1+L2+定义层」一句；值格仍只攻 `(Z-last)`。

**诚实边界（3 条）**：① 四条 `[T]` 全在**上一版语料已标的级**上，本文件只做命名归并，未给任何新证。② `𝒟_C`（方案 C 的反常维因子）只到 1-loop 形（`HZ9-FINAL.md:58,88③`），高阶级数未核 ⇒ 值格是"精确形"而非"已解形"。③ 零 MC、零 Lean、未改既有文件、未 `commit`、未动位置图。**No QED.**
