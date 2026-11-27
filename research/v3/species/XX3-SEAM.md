# XX3-SEAM — `ι`-显式修形的 4 行缝（`hC` 槽的唯一可行形）

> 只交 `[T]`/`[NG]` · 无百分数 · 每条带「下一手」· 未改既有文件 · 未 `commit` · 未编译 · 零新 MC（`py` 0 次）。只读 `RR10-N1SEAM.md`、`HV5-LEAN.md`。
> 口径：`ι` 的物理内容住**定义**（`HaUI`/`HlimUI`），不住假设；`hC` 由「乘积」变回「缝」的唯一手段 = 让两字段共用 `ι`。

## §0 判词

| # | 断言 | 判 | 依据 | 下一手 |
|---|---|---|---|---|
| **0.1** | 修形后两字段**共同变元 = `ι`（+ `a` 索引 + 同一 `H`）**；旧签名 **0**（`RR10:14`）。缝有了对应物 | **`[T]`** | §1；`HV5:36–38` 逐字 | 用 §1 覆盖 `HV5:34–38` 的副本 |
| **0.2** | `RR10 §2` 的反例（`Ha≡0`,`Hlim:=1`,`z=i`,`2^{−1/2}`）在修形下**不可构造**：`Ha`、`Hlim` 皆为 `ι` 的共轭像，`Hlim` 不再是自由参数 | **`[T]`** | §1 `HaUI`/`HlimUI`；`RR10:36–37` | 反例判别器改抓「自由 `ι`」——已由 `L2Dom` 封 |
| **0.3** | 4 行逐行：L1 `[T|L3+Minlos]` · L2 `[T]`（L3-free，空转告警） · L3 `[T|(L²-dom)]` · L4 `[T|(L²-dom)∧W1∧ι→ι₀∧Laplace]` ⇒ **4/4 `[T]` 条件式** | **`[T]`** | §2 | 落 D0 三定义 |
| **0.4** | **`hC` 从此可申请**（`theorem` 形，`axiom` 退场），但 `T16` 的参数必须换成修形版（`Ha` 不再自由） | **`[T]`**（条件式） | §3 | 改 `QQ12 §2` 计数，交 `lead` |
| **0.5** | 无条件版仍 `[NG]`：四个条件件本身未证 | **`[NG]`** | §3 残余账 | 先攻 `(L²-dom)`（`HY3-C:33`） |

## §1 修形后的 structure（`ι`-显式）

```lean
variable {Ω H Obs : Type*} [MeasurableSpace Ω]
variable [NormedAddCommGroup H] [InnerProductSpace ℂ H] [CompleteSpace H]
variable (μ : ℝ → Measure Ω) (O : Obs → H)
variable (ι : ℝ → H ≃L[ℂ] H) (ι₀ : H ≃L[ℂ] H)          -- ★共同变元：ι 与 ι₀
variable (A : ℝ → H →L[ℂ] H) (Alim : H →L[ℂ] H)

/-- `≤`-界两侧一致：ι_a 及其逆一致有界（`(L²-dom)` 的 Lean 形） -/
def L2Dom : Prop :=
  ∃ C, 1 ≤ C ∧ ∀ a, 0 < a → ‖(ι a : H →L[ℂ] H)‖ ≤ C ∧ ‖((ι a).symm : H →L[ℂ] H)‖ ≤ C

/-- ★定义，不是假设：格点算子搬到 ℋ 上 -/
noncomputable def HaUI (a : ℝ) : H →L[ℂ] H := (ι a) ∘L (A a) ∘L (ι a).symm
noncomputable def HlimUI : H →L[ℂ] H := ι₀ ∘L Alim ∘L ι₀.symm

structure IotaL2Seam (μ : ℝ → Measure Ω) (ι : ℝ → H ≃L[ℂ] H) (ι₀ : H ≃L[ℂ] H)
    (A : ℝ → H →L[ℂ] H) (Alim : H →L[ℂ] H) (O : Obs → H) : Prop where
  memLp  : ∀ a, 0 < a → ∀ O : Obs, MemLp (fun ω => ι a (O ω)) 2 (μ a)
  strong : ∀ z : ℂ, 0 < z.im → ∀ ψ : H,
    Tendsto (fun a : ℝ => resid (HaUI ι A a) z ψ) (𝓝[>] (0 : ℝ))
      (𝓝 (resid (HlimUI ι A ι₀ Alim) z ψ))
```

`resid := resolvent`（`Spectrum/Basic:79`）。**共同变元表**：`ι` 出现 2/2 字段；`a` 索引 2/2；`H` 2/2；旧签名对应列 = `∅ / ∅ / H`。⇒ `IotaL2Seam ≡ A(μ,O) ∧ B(ι,A,Alim)`，**有交**。

## §2 那 4 行（逐行注明消费）

| 行 | Lean 声明（骨架） | 消费什么 | 判 | 下一手 |
|---|---|---|---|---|
| **L1** | `theorem seam_memLp_fin (hL3 : L3Bound varpi a0) (hμ : ∀ a, IsProbabilityMeasure (μ a)) : ∀ a, IsFiniteMeasure (μ a) := fun a => inferInstance` | `L3` F1 的 a-一致矩界 + Minlos（`HY6:33`）产 `hμ`；本行只吃 `hμ` 翻成 `IsFiniteMeasure` | **`[T]`**（条件于 `L3`+Minlos） | `W1` 只给子列（`IAB:33`） |
| **L2** | `theorem seam_memLp (hfin : ∀ a, IsFiniteMeasure (μ a)) (hmeas : ∀ a O, AEStronglyMeasurable (fun ω => ι a (O ω)) (μ a)) (hbd : ∀ a O, ∀ᵐ ω ∂(μ a), ‖ι a (O ω)‖ ≤ M) (hM : 0 ≤ M) : ∀ a, 0 < a → ∀ O, MemLp (fun ω => ι a (O ω)) 2 (μ a) := fun a _ O => MemLp.of_bound (hmeas a O) M (hbd a O)` | `MemLp.of_bound`（`LpSeminorm/Basic:549`）一行 + L1 的有限性；**内容 L3-free ⇒ 空转**，`hC` 里 `ι` 的物理内容不在此行 | **`[T]`**（空转） | 别把 L2 记作 ι-消费点 |
| **L3** | `theorem seam_strong_mod (hd : L2Dom μ O ι) (hA : ∃ B, ∀ a, 0 < a → ‖A a‖ ≤ B) : ∃ Z, ∀ a, 0 < a → ‖HaUI ι A a‖ ≤ Z` | `L2Dom`（`‖ι‖·‖A‖·‖ι⁻¹‖ ≤ C²B`，`opNorm_comp_le` 两行）⇒ Trotter–Kato (ii) 的**模输入** | **`[T]`**（条件于 `(L²-dom)`） | 与 `G1-UNIFORM-ZA` 核名（`HY3-C:48`） |
| **L4** | `theorem seam_strong (hd : L2Dom μ O ι) (hSA : ∀ a, IsSelfAdjoint (A a)) (hTK : ∀ z (hz : 0 < z.im) ψ, Tendsto (fun a => resid (A a) z ψ) (𝓝[>] 0) (𝓝 (resid Alim z ψ))) (hι : ∀ g, Tendsto (fun a => ι a g) (𝓝[>] 0) (𝓝 (ι₀ g))) (hW1 : …) : ∀ z, 0 < z.im → ∀ ψ, Tendsto (fun a => resid (HaUI ι A a) z ψ) (𝓝[>] 0) (𝓝 (resid (HlimUI ι A ι₀ Alim) z ψ))` | ① 共轭恒等式 `resid (φ∘A∘φ⁻¹) z = φ∘resid A z∘φ⁻¹`（`hSA` 保证 `0<z.im ⊆` 预解集）；② `hTK`（= `L3` 的 F2 生成元面）；③ `L2Dom` ⇒ `hι` 的极限仍是同构；④ `W1`（子列→全极限）；⑤ Laplace 1 行（半群→预解） | **`[T]`**（条件于四件） | Laplace 那行转引未逐行（`RR10:60`） |

**读法**：L1/L2 是测度面，L3/L4 是算子面；**四行全部经 `ι`**，故 `hC` 不再是两条互不指涉的义务的乘积。L2 仍空转（`ι` 只出现在 a.e. 界里），这**不是**缺陷而是事实：`memLp` 的专属内容 = `O` 的可测性，与 `ι` 的识别无关 ⇒ 识别住 L4。

## §3 判定（贴回）

**判：4 行 `[T]`，条件式。** 条件集 = `L3 ∧ (L²-dom)`（L1–L3 全部吃这两件）＋ `W1 ∧ ι_a→ι₀ 强 ∧ [Laplace 1 行]`（仅 L4）。**不是**无条件 `[T]`：四件本身未证 ⇒ 整体断言 `[NG]`（§0.5），但**条件 schema 是 `[T]`**，且不再可被 `Ha≡0/Hlim≡1` 一击驳回——`Hlim` 已是 `ι₀∘Alim∘ι₀⁻¹`，自由度为 0。
**`hC` 从此可申请**：`hy6_N1_iota_l2_seam` 写成
```lean
theorem hy6_N1 (hL3 : L3Bound varpi a0) (hd : L2Dom …) (hW1 : …) : IotaL2Seam μ ι ι₀ A Alim O
```
即 `axiom` 退场；代价是 `T16` 的参数换签名（`Ha`/`Hlim` 不再是输入）。**下一手**：① D0 落 `L2Dom`/`HaUI`/`HlimUI`（3 定义，零证明）；② 把 `ι_a → ι₀` 强收敛独立立项（它比 `W1` 更弱、比 `(L²-dom)` 更强，是新的一格）；③ 改 `QQ12 §2` 计数交 `lead`。

## §4 台账（`AGENTS §1` 四要素）

| # | 工具 | 状态 | 逐字调用／检索式 | 命中 | 失败 |
|---|---|---|---|---|---|
| 1 | `read` | **OK**（带行号正文，非空） | `research/v3/species/RR10-N1SEAM.md`(102)·`HV5-LEAN.md`(97) | 2/2 | 0 |
| 2 | `grep` | **OK**（有效响应，非 provider failure） | `IotaL2Seam\|HaUI\|L2Dom` @ `D:\Dev\math` `*.lean` | **0**（真 0 ⇒ 修形从未落 Lean，本文件为首写） | 0 |
| 3 | `glob` | **OK** | `research/v3/species/XX*.md` | **0**（无同名文件，未覆盖） | 0 |
| 4 | `pwsh`/Lean 编译/`git`/`web_*` | **0**（未调用） | 本靶 = 结构判别 + 逐行消费账，无新数可算 | — | — |

**`TOOL-FAIL` 0；空响应 0；`ENOENT` 0。零新 MC；既有文件零改动；未 `commit`；未编译（§1/§2 的 Lean 未过 elaboration，`resid`/`≃L`/`inferInstance` 按 Mathlib 同名用法写）。**

*XX3-SEAM · `ι`-显式修形：`Ha := ι∘A∘ι⁻¹`、`Hlim := ι₀∘Alim∘ι₀⁻¹`（定义），两字段共同变元 = `ι`（旧 = ∅）；4 行 `[T]` 条件式，条件 `L3 ∧ (L²-dom)` + L4 三件；`hC` 可作 `theorem` 申请，`axiom` 退场；残余 = `(L²-dom)`∧`W1`∧`ι→ι₀`∧Laplace。*
