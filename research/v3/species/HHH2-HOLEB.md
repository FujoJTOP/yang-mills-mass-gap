# HHH2-HOLEB — 洞 B：`≥6T` 与 `4t/a+4d` 是同一量的两个速率吗

> 2026-11-27 · 靶：`research/v3/species/EEE1-P3CHAN.md:13,54`（洞 B）。逐字读过：`P3-FOURWALL-DERIVATION.md:44,52,69,71,79,86`；`P3-DRESSED-MASS.md:34–38,44,51,53,65,93–95`；`Z07-HOLEB.md:55–129`。只新建本文件；零新 MC；未 commit。

## §0 判词

**混记，不是两个速率。`[NG]`（作为「同一量两个速率」）。** `4t/a+4d` 是**计数**：`4t/a` = 截面周长 4 × 片数（率 `κ`），`4d` = **四壁 × 横移格距**（`d≤L/2` 有界 ⇒ 率贡献 0、只进振幅）；`6T` = **另一个族的率界**（1×2 截面周长 6 × 片数，率 `1.5κ`）。冲突不在率层，在**分类层**：`𝔓_rest` 按「不含固定 `B_T`」定义、却按「`|P|≥6T`」收账，而横移管恰是「不含 `B_T`、成本 `4T+O(1)`」⇒ 全称句假 `[NG]`，假步在 `DRESSED:53`。

- ★ 下一手：`𝔓_rest` 改按**率**分层（不按含不含 `B_T`），撤 `DRESSED:51` 的全称句。

## §1 两个速率式（逐字 + 各是哪个量的指数）

| 式 | 逐字 | 出处 | 是什么量的指数 |
|---|---|---|---|
| `|P|=4(t/a)+4d` | "A transverse **displacement** `d=|x_source−x_sink|` of the whole sink adds `4d` plaquettes (all four walls become slanted staircases of `N_t+d` plaquettes), so `|P|=4(t/a)+4d` — on `T³_L` the transverse distance is bounded, `d≤L/2=O(1)` as `t→∞`, so this is still `4t/a+O(1)`, **same rate**" | `FOURWALL:69` | **计数（片数）**。`4t/a` = 周长 4 × 片数 ⇒ 衰减率 `κ`；`4d` = 四壁 × d ⇒ **有界附加片数**，不是率 |
| `≥6T+O(1)` | "every `P∈𝔓_rest` satisfies `|P| ≥ 6T + O(1)` (perimeter-6 wider tube, the lowest non-tube) or `|P| ≥ 8T+O(1)` (two tubes / two-particle)"；"lowest rate is the perimeter-6 wider tube (`|P|≥6T`, rate `6\ln(1/u)/a=1.5κ`)" | `DRESSED:51`／`DRESSED:93` | **率界**：周长 6 × 片数 ⇒ 率 `1.5κ` |

三个数字的角色：`4t/a` 的 `4` = 截面周长（单 plaquette 边界的链数，`FOURWALL:44,52`）；`4d` 的 `4` = **壁数**；`6` = 另一截面周长。**同一个 `4` 出现两次、含义不同** —— 混记的入口。`d` 是源汇横截面之间的**平移格距**（不是豪斯多夫距离），`d≤L/2` 与 `t` 无关 ⇒ `4d/T→0`；它进振幅 `a₀`（`FOURWALL:79` 的 "Summing over `y` (all transverse positions)"、`:86` 的 "finite positive transverse sum"）。

- ★ 下一手：全仓凡引 `4t/a+4d` 处一律写 `4t/a+O(1)`，禁止把它当「率族」的列名。

## §2 矛盾性质（裁决）

1. **不是「同量两率」`[NG]`。** 可比的率只有两个：`κ=4ln(1/u)/a`（`4T` 族，含横移）与 `1.5κ=6ln(1/u)/a`（`6T` 族），差 `2ln(1/u)/a`，彼此不冲突。`4d` 不产生率 ⇒ `4t/a+4d` 不能与 `6T` 放在同一坐标轴上比。
2. **真形态是「一个对象、两个文件、两次分派」`[NG]`。** 横移管在 `FOURWALL:79,86` 已折进 `a₀`（率差 0，振幅）；在 `DRESSED:44,51`（`B_T` 固定为 `p₀` 的直管）又被判入 `𝔓_rest`（率 `≥1.5κ`）。同一批对象记了两次账、落在两个率类。
3. **假步定位在 `DRESSED:53` 的逐字句**："If `|P| ≤ 4T + Λ₀` with `Λ₀` the fixed source/sink overhead, `P` is the straight four-wall tube `B_T` plus `≤Λ₀` extra plaquettes" —— 这一步把「成本有界」推成「含固定 `B_T`」。`d≤Λ₀/4` 的横移管成本合格、却不含 `B_T`，正是反例。`[NG]`
4. **两式各自都对**：`6T` = 1×2 截面周长 6 的成本、率 `1.5κ`（`FOURWALL:69(III)`）`[T]`；`4t/a+4d` 作为计数 `[T]`。错的只有全称量词（`:51` 的 "every"）。
5. `FOURWALL:18` 自承残数结论仍 conditional on 缀饰重求和，与洞 B **无关**：那条已由 `P3-DRESSED-MASS:16,116–125` 付清（`GT2-GRADES.md:135` 判回收）`[T]` ⇒ 它不额外恶化洞 B。
6. `DRESSED:38` 的「横移 ⇒ 至少 1 个空间片 ⇒ decoration」是**第三种混记**：刚体横移是领头族的平移分量（率 0），wobble 才是 `O(1)` 装饰；两者都落振幅，故不改变率结论，但措辞须改，否则 `4d` 会被读成「增率」。

- ★ 下一手：撤 `:51` 全称句；`DRESSED:38,53` 改写成「横移 = 领头族平移分量（率 0，进 `a₀`）、wobble = `O(1)` 装饰」。

## §3 正确速率（推导骨架）

换算：每片一个 `c_f=u` ⇒ `w(P)=u^{|P|}`（`FOURWALL:81`），故 `rate(F)=lim_T |P_F|/T·ln(1/u)/a`。记 `ρ:=ln(1/u)/a`、`κ=4ρ`。

| 族 | `|P|` | 未扣熵率 | 相对 `m₁^{dr}` |
|---|---|---|---|
| 直管 + 横移（`d≤L/2`） | `4T+4d` | `κ` | `0` ⇒ 振幅 |
| 单壁平行长装饰（`ℓ≤T`） | `4T+ℓ` | `(4+ℓ/T)ρ ∈ [κ,1.25κ]` | `→0.25κ`（binding） |
| 周长 6 宽管 | `6T` | `1.5κ` | `0.5κ` |
| 双管 / 双粒子 | `8T` | `2κ` | `κ` |

余项**未扣熵**的率下确界 = `1.25κ`（不是 `1.5κ`）；扣每片壁/横移计数熵 `ē`、`ū=ēu` 后，实际第二指数率

`δ = ln(1/ū)/a = ρ − ln(ē)/a`，`δ>0 ⟺ ū<1 ⟺ KP(uK<1)`（`DRESSED:65`），且 `δ<ρ=κ/4<κ/3` 严格 `[T]`（承 `X01-L1C.md:39–82`；`Z07:76–105` 复核）。

**EEE1 的下一手「把 `4t/a+4d` 并进 `𝔓_tube`」= 对但不充分**：并完二分恢复穷尽（横移族归位、振幅正确，`[T]`），但 `δ≥κ/3` 仍不立 —— 抽掉它的是上表第二行（`4T+ℓ` 长装饰，1.25κ 扣熵），与横移无关。

- ★ 下一手：①`𝔓_rest` 定义成「率 ≥ 1.25κ 的族」；②`δ` 一律写 `ln(1/ū)/a`；③回填 `X02-THEOREMA.md:84`（B2 那行挂 `:51` 的判据）与 `W03` 的 `Λ★≥(5/4)²m₁²`（承 `Z07:111`）。

## §4 核验（AGENTS §1 四要素）

| # | 工具 | 状态 | 检索式 / 读入 | 命中 |
|---|---|---|---|---|
| 1 | grep | OK | `4t/a|4t/a \+ 4d|FOURWALL` @ 全仓 | 535（截取 250，余存 spill） |
| 2 | grep | OK | `6T|≥ 6T` @ `research/v3` | 218 |
| 3 | read | OK | `FOURWALL:40–89`、`DRESSED:25–99`、`Z07:55–129`、`EEE1:13,54` | 全返回正文 |
| 失败标记 | — | — | — | 无 TOOL-FAIL |

`docs/position-current.html` 本轮未动（硬规矩「不改既有文件」）；是否回填由 lead 判。
