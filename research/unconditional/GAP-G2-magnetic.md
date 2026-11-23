# GAP-G2-magnetic — 尝试证 `ρ_θ→0`（磁面/'t Hooft 周长律）与判词

> 2026-11-23 · 承接 `BREAK-IN-6.md` §3、`IN-6-attack.md` §3、`confinement/REPORT.md` §2.1/§4.5、`uniqueness/REPORT.md` §5.5。
> **判词：未证。** 得到三件真东西：(i) **可证部分**——强耦合处 `W1_1` 是无条件定理，平凡界为 `O(1)`；(ii) **精确等价**——`W1_1 ⟺` 't Hooft 磁环**周长律** `⟺` 磁屏蔽（对偶 Meissner/单极凝聚），即**禁闭的磁面**；(iii) 与电隙 `W1_0` 比较：**证明上不更易，卡同一堵 RG 墙**。

## 1. 陈述

`T⁴_L` 上取体基本 Wilson 环 `W(C)=(1/N)tr∏U_e`，'t Hooft 通量 `θ∈H²(T⁴_L,Z_N)`：
$$
\rho_\theta(C):=\langle W(C)\rangle_\theta-\langle W(C)\rangle_{0},\qquad
W1_1\iff \sup_{C\subset\lambda_0,\ \theta}|\rho_\theta(C)|\to0\ (\operatorname{dist}(\lambda_0,\partial)\to\infty).
$$

## 2. 能证的（无条件 / 条件）

- **平凡界**：`W=(1/N)tr` 时 `|⟨W⟩|≤1 ⇒ |ρ_θ(C)|≤2`（若取 `W=tr` 无归一则为 `2N`，即 `BREAK-IN-6` 记号）。`O(1)`，**不衰减**——弱耦合只有它。
- **强耦合 `β<β₀`**：Wilson 簇展开含**扭曲 bc**（Osterwalder–Seiler 1978；SZZ 2204.12737）给 `|ρ_θ(C)|≤C e^{−m\,dist(C,θ\text{-twist})}`，**均匀于 `θ` 与 bc** ⇒ `W1_1` 成立。但 `β→∞` 到不了连续统（本质奇点）。此即**唯一**无条件定理。
- **条件一步**：若磁环 `V` 周长律（`|⟨V(C')⟩|≤K^{|C'|}e^{−m\,\mathrm{Perim}(C')}`，`m>0`），配 OS 正性把 `W` 对扭曲的响应按磁涡面自由能控制 ⇒ `ρ_θ(C)≤C e^{−m'\,d(C,\theta\text{-面})}→0`。即"`V` 周长律 ⇒ `W1_1`"。

## 3. 精确等价（核心）

`⟨W(C)⟩_θ` = `W` 在 't Hooft 磁通 `θ` 存在下的期望；`ρ_θ` 是 `W(C)` 与远处 't Hooft/涡面插入的**连通关联**。't Hooft 代数 `⟨WV⟩=ω^{link}⟨VW⟩`（`ω=e^{2πi/N}`）蕴含：
- `V` 面积律与 `W` 面积律**不能同真**（二择一）；
- 故 `ρ_θ→0 ⟺ θ` 磁通被局域化为定域磁管 `⟺ V` **周长律**（`∝L`，非面积 `∝L²`）`⟺` **磁屏蔽 = 对偶 Meissner = 对偶单极凝聚**；
- **失效** = `V` 面积律：磁通铺满盒，`ρ_θ` 活到 `L→∞`（磁退禁闭），两扇区给不同 Gibbs 态。

所以 `W1_1` 不是新问题：它就是**电禁闭在 't Hooft 对偶侧的表述**（`V` 周长 `⇔` `W` 面积的磁面）。

## 4. 与电隙比较、精确缺失输入

- `W1_0`（中性 `0⁺⁺` 隙）`=` Clay (B)。`W1_1` 住**弦/磁扇区**，中性隙**不管它**（`BREAK-IN-6` §2 已证 `W1_0 ⇏ W1_1`）。且禁闭（`σ>0`）逻辑上独立于、不被蕴含于隙（`confinement` §2.1）。
- **更易吗？否。** 同一 UV→IR RG 墙（`confinement O1`）：簇展开半径 `β₀<∞`，连续统要 `β→∞`；非阿贝尔**无 FKG/正关联**，'t Hooft 代数只给二择一、不给选择；4D **无严格对偶**；弱耦合唯一 LSI 严格结果仅 `U(1)₄`（Coulomb，无禁闭），非阿贝尔无类比。
- **精确缺失输入**：均匀于 `L,β→∞,θ` 的连通界
  $$|\langle W(C)V_\theta\rangle-\langle W(C)\rangle\langle V_\theta\rangle|\le\varepsilon(L)\to0,$$
  等价于证明**连续统 `Z_N` 中心不破**（'t Hooft 二择一的"选择"）——即一条**4D 严格对偶 Meissner/单极凝聚定理**，或一条跨过 Savvidy 点（`β≈2.47`）且保 't Hooft 代数的 RG。**与 G2-电场恰是同一件缺失工具。**

严格范本仅 3D 紧致 `U(1)`（Polyakov 1977；Göpfert–Mack 1982）；4D 纯 YM 只有 Seiberg–Witten 软破缺 `N=2` 与 abelian-projection 图像，皆 conjectural。**判词：`W1_1 ≈` 禁闭本身（对偶面），不弱于电隙；无条件部分止于 `β<β₀`。**
