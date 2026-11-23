# ABC2-DUALITY — SHAPE A 第二轮：有限-β 动量壳大-N 因子化 / 精确环方程闭合 / 同通道可积转移 —— 皆不闭合 #5，逐字归约为 H1

> 2026-11-23 · 承接 `ABC-DUALITY`（第一轮精确可解：large-N master field / MM / 2D Migdal / 可积扇区 = 错向/错道/错区）。本轮靶已由 WSTC（高场锥 measure-starved）、浅带 `{2g|F|≤N/12}`、AF-UV 深锥三个定理**钉死在一个对象**：β=2.5（g≈1.5）的 **typical-field、high-momentum 粗模方差** `Var_μ(⟨v,∇S⟩|0)`。**判词：A2 不闭合 #5——它就是 H1。** 三结构死法不同、墙同一张：(a) 因子化把方差精确改写成**正定平面谱密度**，低-q 上限 = 平面隙 = H1；(b) 环方程闭合 = 截断累积量塔，被弃三阶累积量正是 conditioned 方差 = H1；(c) 同通道可积 sub-sector 的 min–max **方向反了**，且 β=2.5 无非-SUSY 同通道可积结构。唯一真战利品仍是 #4：自由场返回 λ=0，判别非空过。

## 0. 靶（精确）
`G″(0)[v]=E_μ[PᵀHP|0][v,v]−Var_μ(⟨Pᵀ∇S,v⟩|0)+W`（k0-CONSTRAINT-EFFPOT；实测 `G″≈1/Var(a)`）。`G″≥λ ⟺ Var_μ(a_{k→0})≤1/λ`（PROVE-Gpp §5）。WSTC 证高场锥贡献 `e^{−cβ(κ²−κ_typ²)N_□}→0`，浅带与 AF-UV 已闭（ABC-CONSTRUCTIVE §2）。残余落在 **κ_typ、high-momentum**：W6-BE 实测 SU3 42%、SU2 21.7% 方向负、90–97% 在 fluctuation/UV 块、`−1.12β`。

## 1. A2 = 有限-β 动量壳大-N 因子化
取壳 `{|k|≈q}`、`v` 局域于壳。因子化 `⟨O₁O₂⟩=⟨O₁⟩⟨O₂⟩+O(1/N²)` 使
```
Var_μ^(q)(⟨v,∇S⟩) = Σ_n |⟨0|(∇S)_v|n⟩|²/(q²+m_n²) =: ρ_∞(q) ≥ 0，  (谱表示，精确)
G″(k→0) = E_μ[H_vv] − Σ_n |⟨0|(∇S)_v|n⟩|²/m_n² + W,
Var(a)≤1/λ  ⟺  Var(a_{k→0})=Σ_n |⟨0|(∇S)_v|n⟩|²/m_n² ≤ 1/λ.
```
**这是 A2 的真新意**：动量壳因子化把 finite-β 方差精确写成**正定平面谱密度**，于是 H1 被改写成 `Σ_n |⟨0|(∇S)_v|n⟩|²/m_n² ≤ 1/λ`——**逐字就是平面 `0⁺⁺` 质量隙**（`m₁≥√λ` 且谱和收敛）。正性是自动的（ρ≥0），**符号落在 `E[H]−Σ`**；要上界低-q 权重，sum rule（`∫ρ q⁰ = ⟨[f,[H,f]]⟩`）只钉**总**权重、不钉低-q 分布——**低-q 上限必须载波有隙**（谱范式之墙）。且 `N=∞→3` 不传输：1/N² 修正 `O(1/9)` 非小、`ρ_3` 可被 finite-N 增强，恰恰**可闭合隙的修正正被丢弃**（PROOF-LARGE-N-RIGOROUS：`c/9` 在 N=3 无界；4D finite-β master field 本身未构造，严格者仅强耦 asymptotic）。**A2 ≡ H1。**

## 2. A2′：有限-N 精确 4D 环方程闭合
Makeenko–Migdal 方程 finite-N 是精确 Schwinger–Dyson（#1 真非微扰）。但「闭合」粗模二点 = 用二点表出三点；被弃的连通三阶累积量**恰是** `Var(G|0)` 的粗模分量。且方程+正性不逼衰减（自由无质量场满足环方程却无隙）；环变量是电弦 `√σ`（1-形式），跨道。不动点求解仍须一个正性判据 = H1 换名。

## 3. A2″：同通道（0⁺⁺ 标量）可积转移
非-SUSY 4D β=2.5 同通道无非平凡可积扇区（N=4/N=2 是 SUSY，排除 #4；2D 't Hooft 与各向异性 `S¹` 是弦道）。即便造出，`spec(T|_P)⊂spec(T) ⇒ gap(全)≤gap(子)` 给**上界**；要下界须把 4D **嵌进更大理论**（更难、无隙性不明）。与第一轮同一死法（`FORGE5-SOLVABLE-CORNER`）。

## 4. 酸测（唯一真过，非空过）
自由无质量场：`∇S` 的作用是零质量传播子，`ρ_∞` 在 `q=0` 为**无质量 δ** ⇒ `Var(a_{k→0})=∞`、`G″=0`、`f₁²=0`，装置**无法返回 λ>0**（`Var≤1/λ` 在任意 λ 处失败）。对 YM，ρ 是真密度、方差有限 ⟺ 隙，故判别非空。**#4 干净过**——判别量正是谱密度在 q→0 的质量位置。

## 5. 五点裁决
1. **非微扰 — 形式过、可用性败。** 因子化/环方程真非微扰，但严格大-N 在 4D finite-β 未构造；可用者仅强耦。
2. **控符号 — 败。** `ρ_∞≥0` 是正定量，符号在 `E[H]−Σ`；负号需低-q 权重上界 = 隙。
3. **非循环 — 败。** 因子化→平面隙；环方程闭合→条件方差；sum rule→要载波隙。三处皆 H1。
4. **非-SUSY + 败自由场 — 干净过。** 自由场 `Var(a)=∞`、`G″=0`、返回 λ=0。
5. **均匀下界 — 败。** 低-q 积分无 λ 上界；1/N² 与 `N=∞→3` 非均匀；min–max 反向。**A2 恰好 H1。**

**最锐一句。** 动量壳因子化把 β=2.5 的方差精确写成**正定平面谱密度**——这是把墙算得最准的一步，它证明符号完全落在「低-q 谱权重是否可积」上；而「ρ 在 q→0 有统一质量」逐字就是 `Var(a_{k→0})≤1/λ` = H1。A2 没有闭合 #5，它把 #5 写成了自己最清楚的形式；墙不在精确性，而在「低-q 上界必须自带载波隙」。

*Refs.* ABC-DUALITY；ABC-NEW-INEQ（WSTC）；ABC-CONSTRUCTIVE（浅带/AF-UV）；FORGE5-SL-UV-CONE-EMPTY（K1 sup measure-blind, K2 high-momentum cone）；FORGE5-CONE-VARIANCE（SD `E[H]=E[f²]` 无符号）；PROOF-COARSE-POINCARE（H1 精确形式）；PROOF-IR-H1（H1=隙，非循环=改名）；PROVE-Gpp-positive §5（`G″≈1/Var(a)`）；k0-CONSTRAINT-EFFPOT；SELFOPT-5OF5（(BT′)=H1）；PROOF-LARGE-N-RIGOROUS；FORGE5-SOLVABLE-CORNER；W6-BE-at-2.5；Makeenko–Migdal PLB 88 (1979) 135；Migdal Sov. Phys. JETP 42 (1976) 413；Chatterjee CMP 366 (2019) 203。
