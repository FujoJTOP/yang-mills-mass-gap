# ABC-DUALITY — SHAPE A：4D 非-SUSY **精确可解 / 严格对偶**能否给 `0⁺⁺` 隙一个非循环下界？

> 2026-11-23 · 承接 `FORGE5-DUALITY`（精确电磁对偶=同构，跨道）与 `FORGE5-DUAL-SUPERCONDUCTOR`（直算单极气体，循环+跨道）。**本探针对准"精确可解"这一支**：不问"是否存在对偶"，而问"任何非-SUSY **精确结构**（large-N master field、't Hooft/MM 环方程精确解、2D Migdal 的常模 4D 类比、可积扇区）能否**非循环**地给出同通道 `0⁺⁺` 下界"。**判词：确切失败——不是"没找到"，是"精确性本身锁死在错道+错向+错区"。四条新候选各自死法不同，合起来恰好是已登记的三面墙：跨道（`√σ`≠`0⁺⁺`）、错向（精确约化=上界）、错区（强耦/2D/β=0，非 β≈2.5 连续统）。**

## 0. 最强形态
把强耦有隙 YM 精确映到弱耦**可解**对偶 `𝒯_D`，带**同通道**字典 `trF² ↔ 对偶中性标量`（不是 't Hooft 电/磁 `W↔V` 那种 1-形式共轭），显式给出 `m_{0⁺⁺}=m_D≥λ`。要害在"精确可解"必须**独立**控制弱耦端——否则对偶=同构=循环（`FORGE5-DUALITY` S2）。

## 1. 四条新候选，逐条死法
- **(i) large-N master field（Chatterjee 强耦精确）。** 真实严格：`CMP 366 (2019) 203`（Jafarov `1610.03821`；`1604.04777` 渐近）证任意 `d` 存在 `β₀(d)`，`|β|≤β₀ ⇒` Wilson 环严格收敛为格点弦轨迹和——一条**已证 gauge/string 对偶**、确定性 master field、面积律**上界**、`N=∞` 闭 Makeenko–Migdal。但 **β 固定小=强耦**：此处隙已是 Osterwalder–Seiler/P3 的（`P3: β<β₁` `f₁²>0`），master field 本身**不声称隙**；环方程`+`正性**不逼**衰减（自由无质量场反射正且无隙）；`N=∞→3` 是渐近级数，`c/9` 在 `N=3` **无界**（`PROOF-LARGE-N-RIGOROUS`）。母题还是 Wilson 环（1-形式/`√σ`）——**跨道 #2/#5**，且强耦**错区 #3**。
- **(ii) 't Hooft / Makeenko–Migdal 环方程精确解。** 有限 N 线性、`N=∞` 二次（因子化），**4D 无已知精确解**；2D 才有 Migdal 解且**谱无隙**（自由标度，`σ∝β`，`G″=0`）。方程+正性不逼衰减（同上）。且环变量是**电弦/1-形式**，不是电标量 `0⁺⁺`——跨道。
- **(iii) 2D Migdal 的常模 4D 类比。** 常模在 2D/β=0 Haar/cactus/U(1) **确可精确解**，但**每个精确模型都是高斯切片 `D=0`**：`Var_red=1/S ≤ 1/(S−D)=Var_true ⇒ gap_red≥gap_true`，给的是**上界**——恰是要避免的方向（`k0-EXACT-ONEMODE`）。4D 更狠：**无常模不变子空间**（plaquette 是动量卷积，三阶把 `k=0` 打到一切 `k`），树势在 Cartan 谷上**恒为 0**，而消元所得曲率 `V_eff″(0)=G″(0)=m_eff²(k→0)` **就是目标换名**，一阶值 Savvidy 不定（`k0-CONSTANT-MODE-TRANSFER`）。规范不变的"常模"其实是 holonomy/中心模，与载隙的通量扇区**正交**。
- **(iv) 可积扇区。** 角落可解（2D Migdal、`d→∞` 平均场、强耦、各向异性 `S¹`），但唯一的精确搬运——**不变子空间上的 min–max**——把 4D 隙**上界**（`spec(T|_P)⊂spec(T) ⇒ gap(T|_P)≥gap(T)`），方向反了；要下界须把 4D **嵌进更大**理论（更难、无隙性不明）。2D 给的是 `0`，零搬不动正号；`1/d` 余项在 `d=4` 不受控；SU(N) **无格点 FKG** ⇒ 无 `β`-单调搬运（`FORGE5-SOLVABLE-CORNER`）。

## 2. 三面一墙
四条死法是**同一堵墙的三张脸**：**(跨道)** 精确结构皆在 Wilson/环/弦（1-形式/`√σ`），局部中性 `trF²` 证明性失明（`PROOF-1FORM-ANOMALY`；实测 `m_{0⁺⁺}/√σ≈3.9`，`GL7`）；**(错向)** 任何精确约化=高斯切片 `D=0`，给**上界**，欠的 `D` 恰是 H1；**(错区)** 严格者均在强耦/2D/β=0，而目标在 β≈2.5 连续统。唯一的同通道精确对象 `V_eff″=G″=m_eff²(k→0)` 在 β=0 是真定理（`N/12`），在 β=2.5 需 `Var(G|0)≤…`=H1（`k0-CONSTRAINT-EFFPOT`）。

## 3. 五点裁决
1. **非微扰 — 形式过、可用性败。** 精确结构真非微扰，但严格的都在强耦/β=0/错维；到 β≈2.5 无一受控。
2. **控符号 — 败。** 无一触及 Savvidy 电 Hessian 号 `G″=S−D`；master field 值是 1-形式 `√σ`，2D `G″=0`，`V_eff″` 的β=0 部分（`W_geo≥0`）是周期/Haar 项、非动力学号。
3. **非循环 — 败。** 至 `0⁺⁺` 的桥要么是未造的弱耦 4D 平面闭环（MM 不闭），要么是"最轻弦模=最轻中性态"（假：`√σ≠3.9√σ`），要么是 `V_eff″`=H1 改名。
4. **非-SUSY + 败自由场 — 干净过（仅常模/Haar 支）。** 纯 YM；自由无质量零模平 ⇒ `V_eff″=0`、`f₁²=0`，β=0 SU(2) 得 `0.248`——真判别（非空过）。但它只在 β=0 是定理。
5. **均匀下界 `m_eff²(k→0)≥λ` — 结构败。** 精确 ⇒ 上界（`gap_red≥gap_true`）；large-N 强耦隙→连续统非均匀；min–max 上界；`N→∞→3` 无界。无 `λ`。

## 4. 精确失败
```
(i) β 固定小=强耦 ⇒ 隙是 P3/OS 的，master field 不造隙；N=∞→3 无界。
(ii)/(iii) 精确 ⇒ Var_red=1/S ≤ 1/(S−D)=Var_true ⇒ 上界；D=H1。
(iv) min–max 限制 ⇒ gap(4D) ≤ gap(子扇区)；2D 给 0。
同通道 V_eff″=G″=m_eff²(0)：β=0 定理，β=2.5=H1(=Var(G|0)≤·)换名。
```
**精确可解不是新杠杆——它是错向界的构造性证明**：一旦"精确",即钉死高斯切片 `D=0`，而被丢弃的 `D` 恰是 H1/纤维 Savvidy 方差。同通道只剩 `V_eff″` 在 β=2.5 的同一估计。

## 5. 最锐一句
精确性把量算**准**，却算的是**错道（`√σ` 非 `0⁺⁺`）上的错值、以上界形式**：2D Migdal 的 `G″=0`、large-N 的强耦 `√σ`、常模的 `1/S≥1/(S−D)` 全指向同一缺失 `D=H1`。SHAPE A **确切失败于 #2（错道值非符号）、#3（精确=循环/换名）、#5（上界非下界）**；#1 形式过、#4 仅常模支干净过，皆不救。#4 的自由场判别是真的——它是唯一真战利品，却只在 β=0。

*Refs.* Chatterjee CMP 366 (2019) 203 / arXiv:1502.07719；Jafarov arXiv:1610.03821；Chatterjee–Jafarov arXiv:1604.04777；Makeenko–Migdal PLB 88 (1979) 135；Migdal Sov. Phys. JETP 42 (1976) 413；Witten CMP 141 (1991) 153；Osterwalder–Seiler Ann. Phys. 110 (1978) 440；Shen–Zhu–Zhu CMP 400 (2023) 805；Morningstar–Peardon PRD 60 (1999) 034509。内部：FORGE5-DUALITY；FORGE5-DUAL-SUPERCONDUCTOR；FORGE5-SOLVABLE-CORNER；PROOF-LARGE-N-RIGOROUS；k0-EXACT-ONEMODE；k0-CONSTANT-MODE-TRANSFER；k0-CONSTRAINT-EFFPOT；PROOF-1FORM-ANOMALY；PROOF-P3-UNCONDITIONAL；GAP-G2-magnetic；H1-T10-dual。
