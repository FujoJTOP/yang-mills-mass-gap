# STATUS — uniqueness (infinite-volume limit & vacuum)

> 2026-11-22 · 维护：uniqueness 子会话

## 已完成

- 读出 OSforGFF `OS/Axioms.lean` 的精确陈述（OS0–OS4，含 `OS4_Clustering` 与 `OS4_Ergodicity` 两形式），对齐 `YMStatement/Basic.lean`（标量场版陈述层）。
- 建立唯一性/真空唯一/普适性的精确层级与反例（REPORT.md §5，例 1、2）。
- 给出定理 A（W1 ⇒ 唯一 + 真空唯一）、B（MG-spec ⇒ 一致速率与 a-一致）、C（普适性需 RG）及引理链 L1–L7。
- 列出"超出 MG"的逐项输入表（REPORT.md §6）与障碍（§8）、Lean 计划（§9）、下一步（§10）。

## 关键结论（四行）

1. **"真空唯一"≈免费**（OS4/遍历 ⇒ 一维不变子空间，已在 16 输入）；**真缺口是"极限唯一"**。
2. **极限唯一的最小充分输入是按扇区的 W1**（中性 W1$_0$ + 基本 W1$_1$，速率任意），不是 MG；MG-单态推不出 W1（2D Ising 两相，例 1）。
3. **扇区（对接 confinement §2.1/§4.5）**：中心荷序参量由 **OS4 + 中心不变 bc 免费**给出（中心不破，§5.5.2）；残余电荷输入是"**弦扇区唯一**"（非 $\sigma>0$）。唯一性弱于质量隙**且**弱于禁闭。
4. **普适性需 RG 稳定性，独立于 MG**；MG 的可观测类（是否含 Wilson/Polyakov）须核实。

## 未决 / 待核实

- [ ] **hexagonmath 2611.00492 的 MG 量化形式**：(i) 是否对**边界条件**一致（否则唯一性失败，例 1）；(ii) 可观测类是否含 **Wilson/Polyakov 电荷算符**（若仅局域中性，则不含 W1$_1$）。**保守默认：按独立需求处理"MG for Wilson/Polyakov"**，直至审计。（最高优先）
- [ ] **(MIX) 源**：构造的 OS4 是否全强度（一切 $f,g$、一切方向 $a$）？若是，(MIX) 免费（§5.5.2 四步）；若仅局域/单向，须显式假设 (MIX)。
- [ ] BC 是否中心不变；连续层 $P$ 的重整化有界性/可积性。
- [ ] **弦扇区唯一性**（§5.5.4）：能否由"中性唯一 + Reeh–Schlieder/OPE"推出（路 (b)），从而免去 W1$_1$？格点侧未严格化。
- [ ] 构造所用 **bc 是否中心不变**（周期/自由）；若为扭曲 bc，需另证极限中心不破。
- [ ] OSforGFF 的 `OS4_Clustering`（实测试函数）是否足以推出 $\mathcal A_{\rm gi}$ 上的遍历性；需规范不变版陈述层确认。
- [ ] W1 是否有比 MG 更易的 RG 证明路径（无 FKG，暂无捷径）。
- [ ] (RG) 假设是否已隐含于 RG1–RG5/UV1；"跨正则化一致"作为独立命题的现状。

## 下一步

1. 核对 16 输入中 MG 的量化形式（读 2611.00492 §MG / §E4）。
2. 形式化 VAC 概率层：`exponentialClustering ⇒ ergodic ⇒ 唯一真空`（`YMStatement` 旁，0 sorry）。
3. 与主会话对表：把"极限唯一"从 MG 目标中拆出，登记 W1 为独立里程碑。

## 文件

- `REPORT.md` — 主报告（W1/MG/唯一性，扇区，typed contract）。
- `R32.md` — 主 session 的 R32 任务：U-FULL vs U-ROT、2D 对照、收缩路线、与 MG、Lean 接口。
- `STATUS.md` — 本文件。
- 共享备注（仓库外）：`C:\Users\hooya\AppData\Local\Claude-3p\side-session-notes\local_893d10be-22f3-4514-9901-4f61e62e9204\uniqueness-handoff.md`

## R32 判定（2026-11-22，新增）

- **U1 = U-FULL（全唯一）= hexagonmath O1**；**U2 = U-ROT（逐极限旋转不变）= OS2/E1**。U1 ⇒ U2，U2 ⇏ U1（`R32.md` §1）。
- **4D 对应物 = (b) U2⁺ = U2 + 形状无关性**（Lemma 3.8 字面：固定裸耦合下 $s$-无关 + 旋转不变），**不是 U1**。强度 `U2 ≤ U2⁺ ≤ U1`；"R32 = 母输入 ⇒ U1 ∧ U2"**作废**；os2 的"唯一不动点"（C′）是**独立更强的另一思路、非 2D 路线**。4D 的 U2⁺ 陈述可搬、方法退化为 SYM-4D+SYM-IR（路线 B）。
- **U2 是结构字段 `Spec.Final.phase.phase.sym`，不是定理**：源级链 `E1_rotation_invariant ⇐ scalarKernel_rotation ⇐ F.phase.phase.sym`（§1.5）。
- **收缩路线 = hexagonmath O1(a)**（RG2 升级为 volume-uniform difference estimates），是**缺失输入**，不在 16 内。
- **2D 不可替代处 = 一维流/超可重整**（使各向异性算符不存在、Lemma 3.8 精确）；4D 因无穷维流必含 $D_6$，只剩 $O(a^2)$ = SYM-4D + SYM-IR。
- **os2 已回源码逐字**（commit `adc7f12…`）：Lemma 3.8 = `lem:dense-periods` 是**刚性引理**（两**精确**关系 + 连续性 + 重叠 ⇒ 旋转不变 + $s$-无关），**无不动点前提**。载荷步 = (3.8) `G_{s+d}=(O±⁻¹)_*G_s` 的**精确性**（一个**两轨迹恒等式**）。⇒ **U2 是正确 4D 目标**；os2 定理 1′（唯一不动点）过强、非 2D 所用；4D 缺的正是同一"两轨迹"可控估计（与收缩路线 O1(a) 同对象）。
- **路线 D（flux decay，hexagonmath part iii）最可能弱**，缺口是 `NoSpacetimeGoldstone`（时空间对称破缺不受 Goldstone 直接约束）。
- **路线 C 须拆分**（gluing 表拆两行）：**C(a) = O1(a)**（两轨迹差估计，**带速率** $\le A\kappa^N$ ⇒ 含 B=SYM-4D，且 ⇒ U1）；**C(b) = C′**（唯一不动点，**无速率**，⇒ U1 但 **⇏ B**）。故 `B ⊂ O1(a)`，而 **`B` 与 `C′` 不可比**。

## 归一化（GL1/A0，已修）

- 测度段（§1.1）：`β=2N_c/g²`、`γ=β/N_c=2/g²`，`N_c:=dim f`。
- RT 段（§1.1）：`b₀=11h^∨/(48π²)`、`λ_G:=h^∨g²`、`C_RG=24π²/11`（配 λ_G）。**两套变量勿混**，仅 SU(N) 时 h^∨=N_c；反例 SO(3)≅SU(2)/ℤ₂。
