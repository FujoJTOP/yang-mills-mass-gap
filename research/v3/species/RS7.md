# RS7 — 活度 `ū > 1`（或耦合不经小）而团簇/聚合物展开仍收敛：文献核对

轮次：本轮（第 N 轮，2026-11-25，时间戳仅作追溯）
任务：查"`Σ ū^T` 在 `ū > 1` 时仍收敛 / 仍给指数小界"的已知机制。
约束：只写本文件；无 git commit。

---

## §1 检索执行核验（AGENTS §1 四项）

| # | 检索式 | 工具 | 状态 | 命中数（返回源数） |
|---|---|---|---|---|
| 1 | `cluster expansion "activity" "greater than one" convergence` | `web_search` | OK（有效响应，无 ENOENT / provider failure） | 10 |
| 2 | `polymer expansion large activity convergence hard core` | `web_search` | OK | 10 |
| 3 | `"Kotecky-Preiss" criterion activity large` | `web_search` | OK | 10 |
| 4 | `Mayer expansion convergence radius activity` | `web_search` | OK | 10 |
| 5 | `high density expansion inverse activity 1/z convergence lattice gas particle-hole duality` + `cluster expansion activity larger than one polymer weight complement holes` | `web_search` | OK | 10 |

（`web_search` 单次上限返回约 10 源，故"命中数"记为返回源数，非全库计数。）
`web_fetch` 两次均 HTTP 200：arXiv:2112.13134、arXiv:math-ph/0605041（读到**摘要**，正文未逐字核）。

---

## §2 最相关的 3 篇（读摘要）

1. **S. Jansen, L. Kolesnikov, "Cluster expansions: Necessary and sufficient convergence conditions"**,
   arXiv:**2112.13134** [math-ph]；J. Stat. Phys. **189**, 33 (2022)；
   DOI **10.1007/s10955-022-02992-6**（arXiv DOI 10.48550/arXiv.2112.13134）。
   摘要逐字："We prove a new convergence condition for the activity expansion of correlation functions … with possibly negative pair potentials. **For non-negative pair potentials, the criterion is an if and only if condition.** The condition is formulated with a sign-flipped Kirkwood-Salsburg operator and known conditions such as Kotecký-Preiss and Fernández-Procacci are easily recovered. In addition, we deduce new sufficient convergence conditions for hard-core systems in R^d and Z^d as well as for abstract polymer systems. The latter improves on the Fernández-Procacci criterion." — 这是**精确刻画收敛域**的一篇：说明 KP/FP 都不是 sharp，但**给出的收敛域仍以"活度小"表述**。`[C]`

2. **R. Fernández, A. Procacci, "Cluster expansion for abstract polymer models. New bounds from an old approach"**,
   arXiv:**math-ph/0605041**；Comm. Math. Phys. (2007)，DOI **10.1007/s00220-007-0279-2**。
   摘要逐字："We revisit the classical approach to cluster expansions, based on tree graphs, and establish **a new convergence condition that improves those by Kotecky-Preiss and Dobrushin**, as we show in some examples. The two ingredients … (i) a careful consideration of the **Penrose identity** for truncated functions, and (ii) the use of **iterated transformations** to bound tree-graph expansions." — 改进的是**判据的锐度**（常数/阈值），不是"解除 `ū<1`"。`[C]`

3. **M. Ushcats, "Evidence for a first-order phase transition at the divergence region of activity expansions"**,
   arXiv:**1808.06058**；Phys. Rev. E **98**, 042127 (2018)。
   摘要要点："a general analytical expression is derived for the **phase-transition activity (the convergence radius of activity expansions)** … Mayer's expansion with the constant (volume-independent) cluster integrals remains correct up to the condensation beginning." 且"in the example of a lattice gas with **'hole-particle' symmetry**"给出高密度（`z` 大）状态方程 —— **活度展开的收敛半径 = 凝聚活度**，即大 `ū` 侧一般**发散**而非收敛。`[C]`

---

## §3 ★ 判：有没有 "`ū > 1` 仍收敛" 的机制？

**结论：(a) 有——但只是"换变量型"，不是"同一 `ū > 1` 型"。**

- **机制 A（正解，唯一可搬运的）**：**对偶 / 补（"hole"）聚合物系统 + 以 `ū^{-1}` 为展开变量**。
  格气有 particle–hole 对称时，`z → 1/z`，把"空位"当聚合物，其权重 `≈ ū^{-1} < 1`；KP/FP 判据施加在**对偶权重**上即得绝对收敛与指数小界。文献对应 Ushcats 的 hole–particle 高密度展开（arXiv:1808.06058），以及 CNRS 讲义中"Equivalent polymer models; generalizations"一节。**要点：`ū>1` 从来没有被"压小"，被换掉的是变量。** `[C]`
- **机制 B（几何型，非判据型）**：一维 / 有限程无相变体系，压强对**任意** `z>0` 解析 —— 此处"活度任意大仍收敛"成立，但靠的是**无相变**这一几何事实，不是任何 `ū>1` 判据。`[C]`
- **机制 C（半径数值型）**：Kalmykov, "Estimating the convergence radius of Mayer expansions: the nonnegative potential case", Theor. Math. Phys. **116**, 417 (1998)，DOI **10.1007/BF02557147** —— 给出收敛到真半径的上界序列；半径是**具体数字**，某些势下可 `>1`，但**无统一 `ū>1` 定理**。`[C]`

**反面（判死的部分）**：KP 判据 `Σ_{γ'≁γ}|a(γ')| e^{a(γ')} ≤ a(γ)` 在有界度 `D` 下强制 `|a| e^{|a|} ≲ a/D`，即 `|a| ≲ 1/D < 1`；FP 与 Jansen–Kolesnikov 的 iff 判据只改**锐度/常数**，都**未**给出"同一范数下 `ū>1` 收敛"的定理。检索中**没有任何一篇**声称 uniform `ū>1` 收敛。

⇒ 对本路（YM 聚合物展开）的意义：**不能指望"活度大仍收敛"的现成判据**；可搬运的只有**对偶/补系统换变量**（把 `ū` 换成 `ū^{-1}`）与**无相变的几何情形**。`[NG]`

---

## §4 断言等级

- `[C]` 上列三篇摘要所述内容（读摘要，正文未逐字核）。
- `[C]` "KP 在有界度下强制 `|a|<1`" —— 由判据形式直接代数推出，非文献逐字。
- `[NG]` "存在 uniform `ū>1` 收敛判据" —— 检索 5 式 / 50 源，**0 篇**如此声称（工具 OK）。

路径：`research/v3/species/RS7.md`
