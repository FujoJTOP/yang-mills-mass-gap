# 缺陷-bootstrap 路线的条件定理 — 一页摘要

> 2026-11-23 · 汇总 `research/defect-bootstrap/` 的诚实端点。
> 配套主文档：`docs/conditional-theorem.md`。
> **一句话：本路线产出的是条件定理（消费 `IN-0 … IN-11`）+ 部分排除界 `δ_bootstrap`，不是 `(iii)` 或质量隙的无条件证明。**

---

## 1. 路线想做什么

缺陷-bootstrap（DB）路线（`ROUTE.md`）要绕开中间耦合的 RG 墙，改证

> `(iii)` 纯规范、无物质、无 SUSY 的 4D CFT ⇒ 自由。

它的杠杆是 **1-形式 `Z_N` 对称性**：该对称性不作用在局域算符上（`tr F²` 中性），只作用在线缺陷（Wilson / 't Hooft）。所以必须对含 `Z_N`-荷线缺陷的谱做**水平 bootstrap**，局域 bootstrap 看不见它。

若 `(iii)` 成立，则"无隙 IR ⟹ 4D CFT + 非平凡 `IN-7`"给出矛盾，从而**有隙**。这就是这条路线与 RG 路线的分工：DB 不直接算 `β`，而是想从**唯一性**侧证 gap。

## 2. 为什么终点是条件定理

`(iii)` **证不了**（`III-UNIQUENESS.md`、主文档 §3ter）。三条理由：

1. **4D 无 CFT 唯一性定理。** 4D bootstrap 只产**岛**和**解析界**（extremal functional），**没有任何 4D CFT 被证成唯一**；离散 `Z_N` 是选择定则，不是唯一性引擎。
2. **SDPB 只给逐点严格界，不给全局排除。** 岛是对每个采样 `(Δ_i,j_i)` 的区间算术界，**不是**唯一性陈述。排除所有相互作用 CFT 需：(i) 解析 extremal functional；或 (ii) 先验紧致性归约 + 严格扫描；或 (iii) bootstrap+定理混合。有限扫描只给**逐点不可行**，**没有"多撒点"通往唯一性的路**。
3. **方向是反的。** 若已有 bulk gap，bootstrap 反而去 bound 缺陷谱。

**能诚实声称的最强陈述（部分排除）：**

> 不存在具有纯规范 `Z_N` 结构、且 channel `X` 中最轻算符低于可计算 `δ_bootstrap` 的相互作用 CFT。

`δ_bootstrap` **之上仍有相互作用窗口存活**，故严格弱于 `(iii)`。

## 3. 条件定理（精确陈述）

主文档 `docs/conditional-theorem.md` §2。设 `G` 为任意紧、连通、单连通单李群，假设 `IN-0 … IN-11`，则存在 ℝ⁴ 上规范不变的 `μ` 满足 OS0–OS4、非平凡、禁闭、群一致，且**转移矩阵有质量隙**：

```
∃ δ = min(δ_W, δ_S) > 0 :   inf spec( H_a|_{1^⊥} ) ≥ δ·Λ   (所有足够小的 a) ,
                         ⇒   inf spec H|_{1^⊥} = Δ ≥ δ·Λ ,   次大本征值 e^{−Δ} < 1 .
```

两条**新增**显式前提：

- **`IN-10`（`hcert`）**：量级桥 `∀μ, cert μ ≤ Δ/Λ`。写进 Lean 签名 `Handoff.lean:77`，旧记"第五条声明"。把"证书为正"接到"真实谱隙为正"。
- **`IN-11`（`IN-b`）**：最轻 bulk 态在 `0⁺⁺`（`tr F²`）或 `2⁺⁺`（`T_{μν}`）；弱形式 `(b′)` 要求 `δ` 之上无其他量子数的 bulk 态。**物理输入，不可由 bootstrap 导出**：格点支持（Morningstar–Peardon 1999：`0⁺⁺` 1710 MeV 最轻、`2⁺⁺` 2390 MeV，比值 1.39），**非对称性保护**（迹反常只绑自旋-0 扇区，不保证最轻态是自旋-0）。

从 bulk-channel 下界到真实 `Δ` 的转移**恰好**需要 `IN-11`；随后 `IN-10` 把 `cert` 接到 `Δ/Λ`。这是条件性的全部来源。

## 4. 缺陷-bootstrap 机器（诚实摘要）

- **分块形状**（`BLOCK-CORRECTED.md`，已定）：线为 codim-3，保留群 `SO(2,1)×SO(3)`；`𝔽_j(u,θ)=f(u)P_j(cos θ)`，含**一个** `K_z`-不变 H² 交比 `u`；共形 Casimir 是 `SO(2,1)` 二次 Casimir，即 `u` 上的**超几何型**算子（非朴素 Laplacian）；`Δ_φ` 只经 `u→0` bulk-OPE 归一化进入。**未定**：同时满足 `u→0,1` 两端 OPE 的显式 `u`-形式（`f` 在 `u→1` 错、`Q` 在 `u→0` 错）。
- **crossing/junction 的 SDP**（`SDPB-FORM.md`、`SDPB-LAYER.md`）：**未扭曲**反常**不是新求解器**——**标准 SDPB + 每算符相位旋转 + 一条线性等式**。相位 `φ_i` 由 `ω(q_i,q_j)=e^{2πi k q_i q_j/N}` 固定；约束 `Re(e^{−2iφ_i}Λ(F_i))≥0`、`Λ(G_O)≥0`、`Re/Im Λ(cross)=0`（每 `Z_N` 扇区）、`Λ(F_id)=1`。代价约 `N×` 标准 SDPB（SU(3) ~3×）。**真正新**的是**扭曲（变线）扇区**，在标准 SDPB 之外。
- **SDPB 输入组装**（`SDPB-INPUT.md`）：`SO(3)` 精确（`P_j`），投到横向分波 `j`，拆成 `(j,q)` 的 1D 问题。`N_u=60` 映射-Chebyshev（`w=log((1−u)/u)`，`u∈[10⁻⁴,1−10⁻⁴]`）、`j=0…J`（`J≈10`）、`q∈Z_N`。缺陷块 `f_Δ(u_k)`；bulk 块 `G_{O,j}` = Dolan–Osborn 分波（identity、`T_{μν}`(ℓ=2)、`tr F²`(ℓ=0)）。**未验证/未跑代码**：bulk 分波依赖缺陷运动学约定（`z,z̄↔u` 映射是弱点）；`φ_i` 依赖 level `k`；`u`-节点须分辨 `u=1` 附近 `(1−u)^{Δ_j}`。
- **更底层缺口**（`RESOLUTION.md`）：真正 bulk–defect 块有 **3 个交比**；文献里闭式超几何块是**缺陷 4-点**块（不同对象）；非 SUSY 纯规范 `Z_N` 设置**未建立**。即形状已定，3-交比块本身**缺失**。

**净结论**：机器**可组装但未建成**——形状已定、SDP 可复用标准 SDPB、输入已逐项写出；显式块、`Z_N` 纯规范设置、3-交比块、代码全部缺失。

## 5. 诚实清单

1. `(iii)` **无法**被 bootstrap 关闭：4D 无 CFT 唯一性定理；SDPB 只给部分排除 `δ_bootstrap`，其上相互作用窗口存活。
2. `IN-b` 是**物理输入**，格点事实、非对称性保护；bootstrap 只 bound `0⁺⁺/2⁺⁺`，其他通道无信息。
3. `hcert`（`IN-10`）是显式量级桥，不是机器证明。
4. DB 的块/设置/代码**未完成、未验证**；它给的是 `IN-3/IN-4` 的证书**候选**，不是它们本身。
5. **出口**：条件定理（消费 `IN-0 … IN-11`）+ 部分界 `δ_bootstrap`。**不是无条件证明。** 翻转条件：解析 extremal functional、或已证紧致性归约 + 严格扫描、或 4D CFT 唯一性定理——三者目前都不存在。

---

来源索引：`ROUTE.md` · `III-UNIQUENESS.md` · `IN-B.md` · `GAP-DICTIONARY.md` · `B-CHECK.md` · `BLOCK-CORRECTED.md` · `SDPB-FORM.md` · `SDPB-LAYER.md` · `SDPB-INPUT.md` · `RESOLUTION.md`；主文档 `docs/conditional-theorem.md` §2、§3bis、§3ter。
