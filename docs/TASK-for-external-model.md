# 任务书（第六版）——附原始数组

> 与第五版的差别：**「`(TP)` 被数据违反」已查明是**观测量不匹配**（见 `docs/ERRATA-35-observable-mismatch.md`）；`(TP)` 对**固定取向** plaquette 恢复为 `[T]`。**
>
> **本轮附上原始数组**：`docs/data-raw/Ck_L{6,8,10,12}.csv`、`Cz_L{6,8,10,12}.csv`（C-order 展平，`%.17g`），以及 `docs/data-raw/W-table.md`（`W_μ(n;κ)` 全表 + 负值清单）。

---

## §0 最重要的一条：数据测的不是论文的观测量

**Monte Carlo（`block_mc.py:88–95`）用的是**

$$\boxed{\ q^{\rm data}_x=1-\bar q_x,\qquad \bar q_x:=\frac16\sum_{0\le\mu\lt\nu\le3}\tfrac12\operatorname{tr}U_{P_{\mu\nu}(x)}\ }$$

**⇒ 六个取向的平均（clover），带一个 `1−` 平移（对协方差无影响）。**

**而论文 `§14` 的定理是对**

$$q_x=\tfrac12\operatorname{tr}U_{P_{12}(x)},\qquad \textbf{固定单取向}$$

**证明的。**

$$\boxed{\ \text{两者不是同一对象}\ \Longrightarrow\ \text{此前报的「矛盾」不成立}\ }$$

**为什么：单取向时 `A_κ` 要么是**单时间片乘法算子**（transfer 方向 ⊥ plaquette 平面），要么是**单方向插入**（transfer 方向 ∈ plaquette 平面）——两种都被 Lemma 14.1 / 14.2 覆盖。而 clover 给出**

$$\langle A_{-\kappa}(0)A_\kappa(n)\rangle=\underbrace{\langle A^s_{-\kappa}A^s_\kappa\rangle}_{\ge0}+\underbrace{\langle B_{-\kappa}B_\kappa\rangle}_{\ge0}+\underbrace{\langle A^s_{-\kappa}B_\kappa\rangle+\langle B_{-\kappa}A^s_\kappa\rangle}_{\textbf{交叉项，无符号控制}}$$

---

## §1 原始数据（本轮新增）

| 文件 | 内容 | 大小 |
|---|---|---|
`docs/data-raw/Ck_L6.csv` … `Ck_L12.csv` | `Ĉ(k)` 的完整 `L^4` 数组，C-order 展平 | `28`–`454` KB |
`docs/data-raw/Cz_L6.csv` … `Cz_L12.csv` | `Γ_L(0,z)` 的完整 `L^4` 数组 | `31`–`505` KB |
`docs/data-raw/W-table.md` | `W_μ(n;κ)` 全表 + 每体积的负值清单（最负 20 个） | `10` KB |
`research/unconditional/w7/export_raw.py` | 导出脚本（可重跑） | — |

**用法（`Cz` 展平为 `L^4`，最后一维最快）：**

```python
import numpy as np
L = 12
Cz = np.loadtxt(f"docs/data-raw/Cz_L{L}.csv").reshape(L,L,L,L)
Ck = np.loadtxt(f"docs/data-raw/Ck_L{L}.csv").reshape(L,L,L,L)
print("fftn(Cz) - Ck 最大绝对差 =", np.max(np.abs(np.fft.fftn(Cz) - Ck)))
```

**⚠️ 我已在本机跑过：`max|fftn(Cz)−Ck| = 9.5e−16 … 1.1e−15`（`L=6…12`）⇒ 归一化无误。**

---

## §2 已闭合（`[T]`）

| 编号 | 内容 | 对象 |
|---|---|---|
(P1) | `∫F²dμ=⟨dF,ℋ^{-1}dF⟩`，`ℋ≥0`，无零模 | — |
(P2) | `‖d‖²≤24` | — |
(P3) | `T_β≥0`（character expansion，`I_{2j+1}>0`） | — |
(P4) | `B_κ^*=B_{−κ}`（时间方向插入在标准 Wilson kernel 下自伴） | — |
(P5) | Lemma 14.1/14.2：`G(n)≥0` | **单取向** |
**(P6)** | **(TP)：`W_μ(n;κ)≥0`** | **单取向** |
**(P7)** | **(Q)：`Ĉ(k)≤Ĉ(0)` 全动量** | **单取向** |
**(P8)** | **`‖Γ_L‖_{2\to2}=χ_L`** | **单取向** |
(P9) | `w(n)≥0` 且 `χ=Σ_nw(n)` **无抵消** | 单取向 |
(P10) | `β≤1/30 ⇒ ℋ≥I ⇒ supĈ≤24` | 两者 |

---

## §3 任务

### 任务 1（最高优先，需重跑）：把观测量换成固定取向

**要做的：**在 `block_mc.py` 里加一个固定取向场

```python
def qfield12(U):
    # 固定 (mu,nu) = (0,1) 的 plaquette，取 q = 0.5*tr，不做 1- 平移
    P = qmul(qmul(qmul(U[0], shift(U[1], 0, -1)), qdag(shift(U[0], 1, -1))), qdag(U[1]))
    return 0.5 * qtr2(P)
```

**同时保存 clover 版**，`L=8, 12`，`β=2.5`。

**要报告的：**

1. **固定取向的 `W_μ(n;κ)`** —— 与 `[T]` 的 `(TP)` 逐点比对。**这是 `(TP)` 的第一次真检验。**
2. **clover 与固定取向的 `χ` 之比**（预期 `O(1)`，但因交叉项不一定相等）。
3. **`q_fields` 也保存**，供 bootstrap。

**判据：**若固定取向的 `W_μ(n;κ)` 全非负 ⇒ `(TP)` 获数据支持；若有负值 ⇒ 那条路真的要重审。

**⚠️ 一个可先做的诊断（无需重跑）：**`W-table.md` 显示负值集中在 **`n=L/2`**（自反点）。`L=6` 时最负是 `(μ=1,n=3)`，而 `μ=0,n=3` 是 `+8.5e−5`。**⇒ `n=L/2` 处的别名结构与 clover 的平均可能共同造成符号。**任务 3 应先把它分开。

### 任务 2（解析）：把 `(TP)` 推广到 clover

$$\boxed{\ \textbf{要证：}\ \langle A^s_{-\kappa}(0)B^{(1)}_\kappa(n)\rangle+\langle B^{(1)}_{-\kappa}(0)A^s_\kappa(n)\rangle\ \ge\ -\,\langle A^sA^s\rangle-\langle BB\rangle\ }$$

**即交叉项不能超过两个对角项之和。**

**可用的结构：**

- `A^s` 与 `B` 都是**规范不变、自伴**的；
- `A^s` 是**乘法算子**，`B` 是**插入**；
- `T≥0` 且 `B_κ^*=B_{−κ}`；
- **`A^s` 与 `B` 在时间平移下互相转移**（`A^s` 作用在同一片，`B` 跨两片）。

**⚠️ 若做不到，退而求其次：证明 clover 的 `Ĉ(k)≤Ĉ(0)`。**这条比 `(TP)` 弱：由 `(TP)` 单取向 + 六取向平均的**交叉项对称性**，`Ĉ^{(clover)}(k)≤Ĉ^{(clover)}(0)` 可能不经过逐项 `W≥0` 而直接成立。

### 任务 3（数值，零成本）：分离 `n=L/2` 的别名结构

**要做的：**用 `W-table.md` 与 `Cz_L*.csv`，把 `W_μ(n;κ)` 按 `n ∈ {0, L/2, 其他}` 分组，报告三组各自的符号分布。

**理由：**负值在 `L=6` 集中在 `n=L/2`；若三组行为不同，说明问题在别名而非物理。

### 任务 4（解析，独立）：`sup_Lχ_L<∞`

**与上面全部独立。**见 `docs/TASK-link8.md`。**注意 `χ` 现在有两个版本（单取向 / clover），二者都要控。**

---

## §4 禁止清单

| 禁止 | 理由 |
|---|---|
**把数据的 `Cz`/`Ck`/`χ` 当作论文 `§14` 定理的对象** | **数据是 clover，定理是单取向** |
用数据否证 `(TP)` | 观测量不匹配（勘误 35） |
用数据支持 `(TP)` | 同上 |
`fftn(Cz)≠Ck` 作解释 | 已核：`1e−15` |
`W` 计算有 bug 作解释 | 已核：`Σ_ne^{ikn}W=Ĉ` 差 `6e−17` |
`Cz` 有常数偏移作解释 | 已核：偏移 `≪1e−6` |
`T≥0` 不足（需插入自伴性）作解释 | 已显式核验 `B_κ^*=B_{−κ}` |
`t_a` 简并作解释 | 任一正交本征基下成立 |
`χ=O(1) ⟹ ‖Γ‖=O(1)`（无 `(Q)` 时） | 反例 `q_x=aξ(−1)^{x_1}+bη_x` |
`m≈2.48` 作稳健质量 | `η=0` 给 `1.90`；`η` 自由给 `4.4–5.3` |
行和 `0.9953`（`β≤0.419`）外推到 `2.5` | 即使 `p=0.02` 也只到 `β_c=0.53` |

---

## §5 优先级

$$\boxed{\ \textbf{任务 3（零成本，先分离别名）}\ \succ\ \textbf{任务 1（重跑固定取向）}\ \succ\ \textbf{任务 2（clover 推广）}\ \succ\ \textbf{任务 4（}χ\le C\text{）}\ }$$
