# NA-3-T 交接卡（给 `lattice-cpb`）— 一页

> 2026-11-23 · 配套详版 `NA-3-T-spec.md`。本卡是**操作摘要**，可直接照做。
> 目标：在现有 SU(3) CM 热浴上**只加观测量**，估 Dobrushin / SDPI 系数 `η`，测 crossover 是否 `η<1` 且 `−ln η` 一致。**不改动力学。**

## 0. 一句话（为何这不撞 OB-1）

**`η` 是局部系数；`η<1` 由 Dobrushin 定理导出聚类、不假定它；全局磁化率 `χ` 才是循环输入。** 故这是**决策输入**，不是证明；不建立 uniform 收缩。

## 1. 脚本改动点（**只加观测量**）

- `su3.py`：暴露**单 link 热浴条件分布**的采样器（给定其余 link，采样 `U_e`）。**不改**热浴/更新逻辑。
- 新增 `dobi.py`：
  - **path coupling（主）**：复制构型，仅在一条 link 改动；两份各跑**一步**热浴，用**最大耦合**（共享均匀随机数/拒绝）更新该 link；记期望 Hamming 距离变化 ⟹ 收缩系数 `λ`（`η ≤ λ`）。
  - **逐对影响（辅）**：邻 link 对，采样 `P_e(·|U)` 与 `P_e(·|U')`（`U'` 仅在 `e'` 改），两样本/耦合估 `‖·‖_TV`；行和得 `η`。
  - **禁 cluster**：不复用任何 cluster / 活动度界。
- 新增 `dobi_run.py`（MC + 上述量）、`dobi_analyze.py`（拟合/出图）。
- `m_eff`：复用 `corr_run.py` / `corr_analyze.py` 的**时间方向连通关联**衰减；`tau_hb`：现有热浴自关联。

## 2. CSV 列（`results_eta.csv`）

```
L, beta, op, eta, eta_err, lambda_path, cert_est, tau_hb, m_eff, ratio, n_cfg, notes
```
（`op` = 核：`link`/`block`；`cert_est ∈ {cert, est, na}`；`ratio = (1-eta)/m_eff`；见证对落盘 `witness_*.npz`：两构型 + 被改 link + 单步耦合严格值。）

## 3. 最小可执行版本（按现状）

- **默认最小版 = SU(3)、`L=4`、`β=2.2 / 3.0 / 6.0` 三点**（保守强耦合 / crossover / 弱耦合对照 `η→1`）。
- 若 SU(2) 代码可用：加 `L=4,6`、`β=2.0,2.2,2.4`（**可选，非阻塞**）。
- 跑两个核：`link` 与 `b=2 block`。

## 4. 复现命令（示例）

```
python dobi_run.py --group SU3 --L 4 --beta 2.2 --kernel link  --ncfg 2000
python dobi_run.py --group SU3 --L 4 --beta 3.0 --kernel link  --ncfg 2000
python dobi_run.py --group SU3 --L 4 --beta 6.0 --kernel link  --ncfg 2000
python dobi_run.py --group SU3 --L 4 --beta 3.0 --kernel block --ncfg 2000
python dobi_analyze.py results_eta.csv
```

## 5. 预注册死亡判据（原文，先登记后跑；读数前置规则见 spec §4.0）

> **读数前置（摘要）：** `β=6.0` 是对照、**永不判死**；`β=2.2` 是 BL-6 冗余点、**不算存活**；**只有 `β=3.0` 判 D3**；判死需 **path-coupling 严格界 `λ≥1`**（`η̂≥1` 点估计不判死）；D4 不成立读"对照无信息"，**非死**；报全表 `η(β,L)`。

- **D1**：`η` 的**唯一严格上界**只能由 cluster expansion 给出 ⟹ NA-3 塌为 BL-6，**关闭**。
- **D1′**：`β=3.0` 处任何非 cluster 工具都认证不了严格 `λ<1`、且唯一严格 `η`-界是完整解析性/cluster ⟹ **实际塌成 BL-6**（预登记塌陷路径）。
- **D2**：数值 `η<1` 但 `−ln η` 随 `L`（4→6）或随 `β`→crossover **→0** ⟹ **关闭**。
- **D3**（**只在 `β=3.0` 判**）——四结果：**D3-DEATH** 严格见证：一构型上 Dobrushin 行和 ≥1（`cert`，**非 λ**）⟹ 关闭、触发 §8；**D3-EVIDENCE** `est` 且 `η̂≫1` ⟹ **已死（待证）**、暂不触发 §8（*当前数据即此档：`η̂=8.07±0.26`*）；**D3-NO-EVIDENCE** `est` 且 `η̂<1` ⟹ 不关不开、移交 NA-2；**D3-NA** 跑不出 ⟹ 存疑。**`β=6.0` 为对照、不判死**；`β=2.2` BL-6 冗余；D4 不判死。
- **D4**：`η<1` 但 `1−η` 与 `m_eff` 无比例 ⟹ **对照无信息**，**不据此判死**。
- **通过（认证）**：严格 `λ<1` **只能来自解析（NA-2）**；格点 `λ̂<1` 仅 `est`，**不是存活认证**。

## 6. 边界（**必写进结果**）

**通过只是决策输入，不建立 uniform 收缩，不进 Lean 链，不替代 `CPB-W/S`。** 判死则触发 `NEW-ARCHITECTURE.md` §8 的结构性结论（触发须是**实测事实**：D1/D2/D3 之一）。
