# LC2 — 4D 格上动物（lattice animals）增长常数的已知严格上界

> 靶：`λ_* ≤ 48 = 16(d−1)`（`eΔ = 54.3656` 为 "unoptimized" 上界，`Δ = 20`）。
> 已给数据：`a_n = 1, 20, 458, 11132, 280175, 7218768`。
> 本轮结论：**文献里没有 `16(d−1)=48` 这个数**；对 `Δ=20` 的类，可搬过来的严格上界最好只到 `Δ^Δ/(Δ−1)^{Δ−1} = 53.0007`（略优于 `eΔ = 54.3656`，但差 48 有 ~10%）；而标准 `Z^4` 多立方体类（`Δ=2d=8`）文献上界是 `6e = 16.31`，远低于 48 —— 两者**不是同一个类**，必须先定类。

## 1. 工具状态（AGENTS §1）

| 工具 | 调用 | 状态 |
|---|---|---|
| `web_search` | 5 次（四检索式各 1 次 + 1 次追 arXiv/DOI） | **OK**，均返回有效结果列表（每条 ~10 源） |
| `web_fetch` | 4 次 | 3 次 **OK**（arXiv abs ×2、Cambridge 文章页 200 但正文被截断）；1 次 **TOOL-FAIL**（`core-cms.cambridgecore.org` → `TypeError: fetch failed`） |
| `pwsh` | 2 次 | **OK**（常数精确计算、目录检查） |
| `write` | 1 次 | **OK** |

## 2. 四个检索式与命中数

| # | 检索式 | 命中（返回源数） | 有无有用命中 |
|---|---|---|---|
| 1 | `lattice animals growth constant rigorous upper bound` | 10 | 有（Barequet–Shalah、Georgakopoulos–Panagiotis、Klarner/Penrose 综述） |
| 2 | `"lattice animals" "growth rate" bound d=4` | 10 | 有（同两篇 + Luther–Mertens 数值表 `λ_4 ≈ 13.28–13.30`） |
| 3 | `Klarner animal growth constant bound` | 10 | 有（Klarner–Rivest 系、`λ_2 ≤ 4.649551`、Shalah 论文摘要） |
| 4 | `site animals hypercubic lattice upper bound` | 10 | 有（Slade 综述、Gaunt 系列展开，均为背景） |
| 补 | `"On the exponential growth rates of lattice animals and interfaces" arXiv DOI …` | 10 | 有（**DOI + arXiv ID 定死**） |

**关于 `16(d−1) = 48` 的专门检索：0 命中。**（检索式 1–4 及补充检索的全部返回源中，无任何 `λ ≤ 16(d−1)` / `λ ≤ 48` / `λ_4 ≤ 48` 形式的定理；数字 48 在我的命中里**从未**与 lambda 同现。）

## 3. 一手文献（读到的层级逐条标注）

1. **Barequet, Shalah — *Improved Upper Bounds on the Growth Constants of Polyominoes and Polycubes***
   arXiv:**1906.11447**（cs.DM，2019-06-27；DOI 10.48550/arXiv.1906.11447）。**读摘要（全文页 fetched，HTTP 200）**。
   原文：`λ_d ≤ (2d−2)e + o(1)` for any `d ≥ 2`；`λ_2 ≤ 4.5252`；`d=3`：`12.2071 → 9.8073 → 9.3835`。
   ⇒ **d = 4：`λ_4 ≤ 6e = 16.3097`**（`[T]`，仅摘要层：定理陈述可直接读，证明未逐字核）。
   适用窗口：`Z^d` 上**面连通**多立方体（site animals），最大度 `Δ = 2d = 8`。

2. **Georgakopoulos, Panagiotis — *On the exponential growth rates of lattice animals and interfaces***
   *Combinatorics, Probability and Computing* **32**(6), 912–955 (2023)；DOI **10.1017/S0963548323000214**；arXiv:**1908.03426**。
   摘要**未逐字核**（Cambridge 正文页 200 但内容被截断；core-cms 镜像 fetch 失败）。
   从检索索引片段读到其引述的界：`ȧ(Z^d) ≤ 2de − 2e + 1/(2d−2)`（归功于 Barequet–Shalah 并行工作）⇒ d=4：`8e−2e+1/6 = 16.4764`。标 `[C]`（片段级）。

3. **Luther, Mertens — *Counting Lattice Animals in High Dimensions***
   arXiv:**1106.1078**；*J. Stat. Mech.* **P09026** (2011)；DOI **10.1088/1742-5468/2011/09/P09026**。**读摘要（fetched，HTTP 200）**。
   内容：`3 ≤ d ≤ 10` 的枚举扩展（perimeter polynomials、`n ≤ 14` 的显式公式）+ 数值增长率估计。
   其表格片段（HTML 索引）：`d=3: 8.34`、`d=4: 13.28 / 13.29（MC 13.30）`。⇒ **数值** `λ_4 ≈ 13.3`，标 `[M]`，**不是界**。
   ⚠️ 检索中出现的 Zenodo 文档（`su2-lattice-mass-gap-strong-coupling.pdf`、`Yang Mills Mass Gap.PP.pdf`）把 "`λ_4 ≈ 8.35`" 归给 Mertens —— **8.34/8.35 是 `d=3` 的值**，属误引；这两份 Zenodo 文本按不可信数据对待，未采信。

4. **Klarner → Penrose（一般正则图度上界）**
   *Electron. Commun. Probab.* 系 PDF（`emis.de/journals/EJP-ECP/.../getdoc07a7.pdf`）检索片段原文：
   > "Klarner proves that for 2-dimensional nearest-neighbour site animals the growth constant is at most `27/4` and Penrose states that this can be generalised to the upper bound `α ≤ K^K/(K−1)^{K−1} ~ K e` for bond animals on an arbitrary regular graph."
   标 `[C]`（片段级，PDF 未 opened 逐字核；无 arXiv/DOI 取到，仅有 URL）。
   ⇒ 这正是 `eΔ` 家族的确切出处。**在 `Δ = 20` 处精确取值**（`pwsh` 算）：
   `20^20/19^19 = 53.000687 < 20e = 54.365637`。

5. （背景）Miranda & Slade, arXiv:**1102.3682**（*growth constants of lattice trees and lattice animals in high dimensions*，`1/d` 展开）——只对大 `d`，非 `d=4` 上界。

## 4. 判（按判死线）

**先定类，否则一切界都无效：**

- 项目的 `a_n` 比值：`20.0000, 22.9000, 24.3057, 25.1684, 25.7652`（`pwsh` 精确算）——**递增且已 > 25.7**。
- 而标准 `Z^4` 多立方体序列的 `λ_4 ≈ 13.3`（文献 3）。⇒ **项目的动物类 ≠ `Z^4` 面连通 site animals。**
- `eΔ = 54.365637 = 20·e` 精确成立 ⇒ 项目的类按 **`Δ = 20`（20-正则度参数）**处理。

于是：

- **(a) 不成立（对 `Δ=20` 的类）**：未查到任何 `λ_* ≤ 48` 或 `≤ 16(d−1)` 的定理。检索式 1–4 + 补检索，工具 OK，命中中 `48` 与 `λ` 同现 **0 次**。→ 缺口**未关闭**。
- **(b) 部分成立，但 `eΔ` 不是最好**：可搬的、对任意正则图成立的最佳**已证**形式上界是 Klarner–Penrose 的 `Δ^Δ/(Δ−1)^{Δ−1}`，`Δ=20` 给 **53.0007 < 54.3656**（`[C]`）。所以"`eΔ` 已是最好"这句**不准确**：经典精确形式就已经好 2.5%。
- **(c) 不适用**：非 0 命中。
- **附加强结论（须先核类）**：**若** `λ_*` 实际是标准 `Z^4` 多立方体增长常数（`Δ=2d=8`），则文献已有 `λ_4 ≤ 6e = 16.31`（arXiv:1906.11447）与 `≤ 16.4764`（arXiv:1908.03426）——**48 这个靶本身太松，问题会自动消失**。所以 LC2 的第一件事不是造界，而是**证明项目的 20-正则动物类与 `Z^4` 多立方体类是同一/可嵌入类**；若不可嵌入（我目前的读数指向"不同类"），才需要自造 `Δ`-参数化的 `(Δ−2)e` 型定理。

## 5. 待办（下一轮交接）

1. `[T/C]` 类识别：把 `a_n = 1,20,458,11132,280175,7218768` 与 OEIS 对（检索式 `1, 20, 458, 11132` 尚未跑）——这是判定"能否直接借用 16.31"的唯一钥匙。
2. `[C]` 若能定位 Barequet–Shalah 证明的**度参数依赖**（其构造是 Klarner–Rivest twigs + 主导有理生成函数），检查 `2d−2` 换成 `Δ−2` 是否成立；若成立，`Δ=20` ⇒ `18e = 48.929`，**仍差 48 一点点**（`48.93 > 48`），这意味着 `16(d−1)=48` 这个靶在 `Δ−2` 路线下不可达，需要更强的 twig 权重。
3. 若 `Δ`-参数化只能到 `Δ^Δ/(Δ−1)^{Δ−1}`，则 LC2 的诚实结论是：**上界 53.00，48 不可由文献搬来**。

## 6. 可信度表

| 断言 | 标 | 依据 |
|---|---|---|
| `λ_d ≤ (2d−2)e+o(1)`，故 `λ_4 ≤ 16.3097` | `[T]`（摘要层） | arXiv:1906.11447 |
| `ȧ(Z^d) ≤ 2de−2e+1/(2d−2)` | `[C]`（片段层） | arXiv:1908.03426 / DOI 10.1017/S0963548323000214 |
| `λ_4 ≈ 13.3`（数值，非界） | `[M]` | arXiv:1106.1078 |
| 任意正则图 `α ≤ Δ^Δ/(Δ−1)^{Δ−1}` | `[C]`（片段层） | EJP-ECP PDF URL（无 DOI） |
| `Δ=20 ⇒ 该式 = 53.0007` | `[T]` | `pwsh` 精确计算 |
| 无 `λ ≤ 48` 的定理 | `[NG]` | 检索式 1–4 + 补，工具 OK，命中 0 |
| 项目类 ≠ `Z^4` 多立方体类 | `[C]` | 比值 25.77 vs `λ_4≈13.3`（推论，待 OEIS 核） |
