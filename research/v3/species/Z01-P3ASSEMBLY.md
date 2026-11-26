# Z01-P3ASSEMBLY — `PROOF-P3-UNCONDITIONAL` §3–§4 的**装配**：四条伴文 + `[KL]` 拼成 `f₁²_L ≥ r₀(β)>0`（`I1` 从 `[C]` 回 `[T] mod [KL]`）

> 2026-11-25 · 任务：按逻辑顺序把 `P3-FOURWALL-DERIVATION` / `P3-FULLREP-KP` / `P3-DRESSED-MASS` / `P3-CYLINDER-RESIDUE` / `X01-L1C` / `Y01-KL` 装成 `PROOF-P3-UNCONDITIONAL §3–§4`。
> 只写本文件；未 `git commit`；标签 `[T]/[C]/[NG]/[M]/[O]/[R]/[J]`。**No QED on Clay；Clay 0%。**
> 上游判词：`Y06-P3DOWN.md:42`（`:42` 降的正是 `f₁²≥r₀` 本身，当下 `I1=[C]`）；`U07-CPSTRONG.md:71`（`[T] mod OS/KL`）；`X02-THEOREMA.md:11`（`I2` 侧的口径）。

---

## §0 判词（up front）

**(a) 触发 —— 四步能拼成，`I1` 回 `[T] mod [KL]`（且 `[KL]=[T]`，`Y01-KL.md:55`）。**
**但不是逐字照抄：装配必须带四处 `[M]` 级修补 + 一处 `[T]` 撤销（§3）。**

**(b) 触发一处，但不落在 `I1` 上。** 真断点是 `P3-DRESSED-MASS.md:94` 的 `δ ≥ κ/3`：`W05-TDECAY.md:66` 已判其**不成立**（单壁平行装饰给 `𝔓_rest` 率下确界 `κ`，非 `1.5κ`；`:69` `δ` 与 `Cu` 是同一族的两种分割，`:72` no-go）。**但 `I1` 只消费 `∃δ>0`（分割固定后），不消费 `κ/3`** —— 残数取的是**权重**不是**间距**（`P3-CYLINDER-RESIDUE.md:151`、`:168`）。⇒ 该断点属 `I2` 的两指数形式，**不进 `I1` 的链**。

**(c) 触发三处引错/过期引用（§3 逐条给行号）。**

**层级必须分开（`Y06 §4` 口径不变）**：本文件是**证据层**装配 ⇒ `I1` = `[T] mod [KL]`。`MASTER:42` 的 `[条件/待重证]` **仍未回填**，总册口径仍 `[C]`；本文件不代替回填（`Y06:140` 逐字："不许做的一步：在 `MASTER:42` 与 `PROOF-P3-UNCONDITIONAL.md` 未回填前，把 `I1` 写回 `[T]`"）。本文件的产出正是**回填所需的 §3–§4 文本**。

---

## §0.1 目标陈述（逐字，不加码）

```
在柱面 𝒞_L = T³_L ×_a ℤ 上，β < β₁ = min(1/48, c/N)，每个有限 L，
对 O = ¼tr F² 的最轻 0⁺⁺ 原子：   f₁²_L ≥ r₀(β) > 0 ，L-一致。
```

**范围（逐字，不许扩）**：`t→∞` 在固定 `L` **先取**；`L` 只作一致性轴，**不取 `L→∞`**（`P3-CYLINDER-RESIDUE.md:12`、`:91–94`）；固定 `a`，**不取 `a→0`**，**不碰 `β→2.5`**（`PROOF-P3-UNCONDITIONAL.md:16`、`:32`）；**不给 `I2`**（`:33`）。

---

## §1 装配链（一图）

```
A0 全表示 KP   u·K<1, K≤6.1×10²N          [T]   FULLREP-KP:96–97,:120–121
      │  （给出：所有聚合物/装饰气体的绝对收敛 + 半径 L,T-无关）
      ├──────────────────────────────┬───────────────────────────────┐
      ▼                              ▼                               ▼
A1 四壁几何 |P|≥4t/a, m₁=4ln(1/u)/a  A2 装饰气体重求和 (a)(b)(d)      A3 谱侧 KL/离散谱
   [T] FOURWALL:52,:86               [T] DRESSED:63–89,:115–125      [T] CYL-RES:52–57,:67–69
      │                              │                               │
      └────────► A4 骨架–装饰精确分解 ◄┘                               │
                 C_L=a₀u^{4T}Z_dec+R_hi  [T] DRESSED:50               │
                        │                                             │
                        ▼                                             │
                 A5 端点缺陷 → 极限幅度：D_fin→D_∞=1+O(u)≥½             │
                        [T] X01:40–43,:79–80（+Y04:84 复核）           │
                        │                                             │
                        ▼                                             ▼
                 A6 原子从指数引理（用 ν≥0）  [T] DRESSED:99–101  ◄──[KL] 载体
                        │
                        ▼
                 f₁²_L = ν_L({m₁^{dr}})= a₀D_∞ ≥ ½a₀ ≥ ¼(2N)^{-2} > 0
                        [T] mod [KL]   ⇒ I1
```

---

## §2 五步逐条装配

### 步 1 · §3.1 四壁管主导：`|P| ≥ 4t/a` 是哪一步、`factor-4` 怎么出

**陈述（逐字引 `P3-FOURWALL-DERIVATION.md:52`）**

> **(1.3) Lemma (four-wall minimality).** Every connected, gauge-invariant, nonvanishing polymer linking the insertions satisfies `|P| ≥ 4N_t = 4t/a`. `[T]`

**来源 + 行号（三段拼出，缺一不成）**

| 段 | 内容（逐字） | 行号 | 等级 |
|---|---|---|---|
| **(1.1)** | "Each of its four boundary links `e₁,…,e₄` therefore carries **one open fundamental index pair**. … if `e_i` is covered by the source alone, the factor is `∫dU (U)_{ij}=0`" ⇒ "**every one of the four source boundary links must be shared with at least one other charged plaquette**" | `:28`–`:30` | `[T]` |
| **(1.1′)** | "Two distinct unit plaquettes … share **at most one link** … **no single partner plaquette can cover two source links**" ⇒ "**≥ 4 distinct partner plaquettes** — four independent 'flux strands'" | `:32`–`:34` | `[T]` |
| **(1.2)** | "**any plaquette-adjacency path from a plaquette in slice `0` to one in slice `N_t` crosses each of the `N_t` time-boundaries `j→j+1`, and each crossing requires `≥1` temporal plaquette** … a fundamental strand therefore contains `≥N_t` temporal plaquettes" | `:40`–`:42` | `[T]` |
| **(1.3)** | 合成：`4` 股 × `N_t` 片、合并无益（"higher-representation temporal plaquette … costs `|r|`, exactly what it carries"） | `:46`, `:50`–`:52` | `[T]` |
| **(1.4)/(1.5)** | 等号 ⇒ 直管唯一；非管 `|P| ≥ 4t/a+1` | `:54`, `:69`–`:71` | `[T]` |

**`factor-4` 如何从 `4t/a` 出（这是 `MASTER:42` 的第一条降级理由）**

`|P|=4N_t` ⇒ 每片 4 个基本 temporal plaquette，每个带一个 `c_f`：

- `P3-DRESSED-MASS.md:35`（逐字）：`w(B_T) = a₀(β)\,u^{4T} = a₀\,e^{−κ t}`，`κ := 4\ln(1/u)/a`；`:34` 逐字 "The bare tube has `|P|=4T+O(1)`"。
- 故 `u^{4T}=e^{−4T\ln(1/u)}=e^{−(4\ln(1/u)/a)t}` ⇒ `m₁ := 4\ln(1/u)/a`（`P3-FOURWALL-DERIVATION.md:86`，逐字 `w_tube(t)=a₀e^{−m₁t}`, `m₁ := 4\ln(1/u)/a`）。
- 外部核对（`[O]`，逐字读到 Gupta 源）：Gupta (14.20) "rectangular tube consisting of four sides of size 1×T"、(14.21) `M_{0⁺⁺}=−4\log(\beta/2N^2)`（`:152`–`:153`）；代入 `u:=β/2N²` 得 `4\ln(1/u)/a = M_{0⁺⁺}`（`:159`–`:161`，`[T]/[O]`）。
- `u` 的**唯一自洽赋值**：`u := c_f/N = β/(2N²)+O(β³)`（`P3-FOURWALL-DERIVATION.md:81`、`:159`；`X01-L1C.md:15`–`16`）。若取 `u=c_f`，则 `a₀` 须携带 `N^{−4T}`，`T`-无关性坏（`X01:16` `[T]` 勘误）——**见 §3 修补 R1**。

**接口（它接上谁）**：`|P|≥4t/a` 只给**几何/代价**；它把"最轻管道"钉成 `4×T`，于是 (i) 领头指数 = `κ=4\ln(1/u)/a` 有了定义，(ii) 其余一切连接聚合物要么 `+O(1)`（并入振幅 `a₀`），要么 `≥6T`（率分离）。**它是步 3 的 `Z_dec(B_T)` 的骨架定义域**（`P3-DRESSED-MASS.md:44` 逐字 "Fix … a **backbone** `B_T`: the straight four-wall tube"）。

**不计入本步的东西**：`a₀>0` 不在本步（属步 3(d)/步 4）；`m₁` 的**数值**不是本步结论（只有 `4` 是）。

*等级：`[T]`（`P3-FOURWALL-DERIVATION.md §1.1–1.5`；`Y06:97` 同判 "`[T]` 已补"）。*

---

### 步 2 · §3.2 全表示 KP：`ū<1` 在 `β<β₁` 下成立吗、`ū` 是什么

**先纠一处误称。** 任务里的 "Seiler 的 KP 判据 `ū<1`" **不是** P3 链的判据出处：

| 对象 | 逐字定义 | 出处 | 等级 |
|---|---|---|---|
| **`ū`（配对活动度）** | "取 `\bar u:=eΔ\bar C u`" | `X01-L1C.md:71`；`Y04-AUDITX01.md:152` 逐字 "`ū=eΔ\,\bar C\,u`" | `[T]`（记号） |
| `Δ` | "plaquette-adjacency degree `Δ = 4(6−1) = 20`" | `P3-FULLREP-KP.md:116` | `[T]` |
| `\bar C` | Ursell 树图常数（`|φ^U(Γ)|≤\bar C^{|Γ|}u^{|Γ|}`）；数值上 `eΔ\bar C≈54.4`（隐含 `\bar C≈1`） | `X01:71`；`Y04:152` | `[T]` |
| **P3 的真判据** | "KP criterion `u(β)\,K<1` ⟺ `u < 1/(2λ_aC_ge)`" | `P3-FULLREP-KP.md:121`（`:96`–`:97` Lemma 2.3；`MASTER:40` 同式） | `[T]` |
| `K` | `K := 2λ_aC_g e ≤ 80e^{2+β}N ≤ 6.1×10²N`（`β≤1/48`），`λ_a=eΔ=20e` | `P3-FULLREP-KP.md:120`、`:116` | `[T]` |
| Seiler 的位置 | "Seiler, LNP **159** (1982), §II … via exactly the mechanism above … the displayed constants are **proved in §1.2** and are **not quoted verbatim** from a numbered formula" | `P3-FULLREP-KP.md:54` | **`[O]`（先例，非出处）** |

**答：`ū<1` 在 `β<β₁` 下成立，但它是被 `uK<1` 蕴含的、更弱的条件，且它不是 Seiler 的判据。**

$$uK<1\ \Longrightarrow\ \bar u=\frac{e\Delta\bar C}{K}\,(uK)\le \frac{54.4}{6.1\times10^{2}N}\cdot 1=\frac{0.089}{N}<1 .$$

即 `uK<1 ⇒ ū<1`（`N≥2` 更强）。反向不成立：`ū<1` 只给 `u<1/54.4`，比 `u<1/(6.1×10²N)` 宽得多。`X01-L1C.md:75` 逐字亦只声称 "`\bar u<1` ⟺ KP (2.1)（`uK<1`，`:64`）——**L1-c 不引入新的 `β` 限制**"，其意在"同一盘内"，严格方向是 `uK<1 ⇒ ū<1`。
数值旁证（`X01:97`，`[M]`）：SU(2) `β_W=1/12` ⇒ `u=0.0104`，`ū≈0.57<1`，`δ=ln(1/ū)/a≈0.57/a`。

**KP 的内容（三步，逐字）**

1. 系数阶乘衰减：`0 ≤ c_r(β) ≤ \frac{β^{|r|}}{|r|!}\frac{e^{β}}{\dim r}` ⇒ `c_r\dim r ≤ β^{|r|}/|r|!·e^β`（`P3-FULLREP-KP.md:29`–`:31`，`:49`–`:50` 自证）；副产品 `c_r≥0` 全 `r`（`:44` `[T]`）。
2. 全局收缩界：`|I(P)| ≤ ∏_{p}\dim r_p`（`:70`–`:72`，Haar 是概率测度 + `|tr D^r|≤\dim r`）⇒ `|w(P)| ≤ (C_gu)^{|P|}`，`C_g ≤ 2Ne^β`（`:76`–`:78`）。
3. 动物计数 + 求和：`Φ(a) ≤ e^β(e^{βe^a}−1)`（`:85`），`S_x(a) ≤ λ_aΦ/(1−λ_aΦ)`，`λ_a=eΔ`（`:89`–`:90`）⇒ `λ_aΦ ≤ ½` 时 `Σ_{P∋x}|w(P)|e^{a|P|} ≤ 1 = a`（`:96`–`:97`）。

**接口**：`Σ_{P∋x}|w(P)|e^{α|P|}≤α` 逐 plaquette ⇒ **在每个 `Λ_{L,T}` 上同一个条件** ⇒ 半径 `L,T`-无关（`P3-DRESSED-MASS.md:65` 逐字 "a per-plaquette condition, identical on every `Λ_{L,T}`, hence independent of `L` and `T`"）。**这一句是整条链的 `L`-一致性总闸**：步 3 的装饰气体、步 1 的熵计数都吃它。

**`β₁` 的数值必须诚实（见 §3 修补 R4）**：`λ_aΦ(a)≤½` 配显式常数给 `β ≲ 3.4×10^{-3} ≈ 1/300`（`|r|=1` 项占优 ⇒ `Φ(1)≈2·(β/2)·e^{1}=eβ=2.718β`（`FULLREP-KP:94`）⇒ `λ_aΦ≈54.4·2.718β=148β≤½` ⇒ `β≤3.4×10^{-3}`），**小于 `1/48`**。`min(1/48, c/N)` 只是**形状**（`P3-FULLREP-KP.md:124` 逐字：literal 归一 ⇒ `β₁=O(1)`、与 `N` 无关、与 SZZ `1/48` 同阶；Hooft 归一 `u≍β` ⇒ `β₁=O(1/N)`）。结论不受影响：`∃β₁>0` 即可，`r₀(β)` 随 `β` 变。

*等级：`[T]`（数学），`[M]`（`β₁` 数值口径）。*

---

### 步 3 · §3.3 dressed-mass：(a)(b)(d) 各证什么、重求和如何拼出

**三个子缺口逐字（`P3-DRESSED-MASS.md`）**

| 子缺口 | 证什么（逐字） | 行号 | 等级 |
|---|---|---|---|
| **(a) 重求和存在** | "The decoration gas attached to the tube is an abstract polymer gas with per-attachment activity `O(u)`; Kotecký–Preiss on the **full** representation sum … gives absolute convergence of the Mayer/Ursell series, uniformly in tube length `T=t/a` and volume `L`." | `:13`；证明 `§2.1–2.2` `:63`–`:71`；`log Z_dec(B_T)=Σ_{connected Γ}φ^U(Γ)` 绝对收敛、在 `uK<1` 多圆盘上解析 | `[T]` |
| **(b) 有效质量存在、`L`-一致** | "`\log Z_dec=φ_dec·T+O(1)=f_dec·t+O(1)` with `φ_dec=O(u)`，`f_dec=O(u/a)`"；`m₁^{dressed}=κ−f_dec`；`L`-依赖只经绕空间周期的装饰 `O(u^{L/a})=e^{−cL}` | `:14`；`§2.2`–`§2.4` `:67`–`:79`；`§3.1` `:85`–`:89` | `[T]` |
| **(d) 一致正振幅地板** | "Decorations are `O(u)`-small in **absolute value** … so `A = a₀(1±Cu) ≥ (1−Cu)a₀ > 0`；`a₀>0` uses only `c_f>0` and the singlet recoupling `f⊗f̄→1` (multiplicity `1`). **No blanket `c_r≥0` is used.**" ⇒ `\inf_L A^{dressed}(L) ≥ ¼(2N)^{-2}` | `:16`；`§4.1`–`§4.3` `:115`–`:125` | `[T]` |
| **(c) 仍是谱原子** | 输入 = `[KL]`（步 5）；本文件自行证的是 Lemma 3.4（原子从指数） | `:15`；`§3.3` `:99`–`:107` | `[T] mod [KL]` |

**重求和如何从 (a)(b)(d) 拼出（四行代数，每行给行号）**

```
(1) 精确分解   C_L(t) = a₀u^{4T}Z_dec(B_T) + R^{hi}_L(t)            DRESSED:50  (1.2)      [T]
(2) 取出指数   a₀u^{4T}Z_dec = a₀u^{4T}e^{f_dec t}D_fin(T)          DRESSED:86  (3.1)      [T]
                          = A^{dressed} e^{−m₁^{dr}t} D_fin(T),  A^{dressed}=a₀ D_∞  ← 修补 R2
               m₁^{dr} := κ − f_dec = [4ln(1/u) − C_dec u + O(u²)]/a  DRESSED:88  (3.2)      [T]
(3) 端点缺陷有界 |log D_fin(T)| ≤ Cu（只用绝对收敛）                  DRESSED:116 (4.1)      [T]
               ⇒ D_fin = e^{O(u)} = 1+O(u) ⇒ D_∞ ∈ [1−Cu,1+Cu]，u<1/(2C) ⇒ D_∞ ≥ ½  DRESSED:123-124
(4) 地板组装   A^{dressed} = a₀D_fin ≥ a₀(1−Cu) ≥ ½a₀ ≥ ½(2N)^{-2}
               inf_L A^{dressed}(L) ≥ ¼(2N)^{-2}                    DRESSED:122,124  (4.2),(4.3)  [T]
```

**(2) 的指数必须穿衣（否则整条链假）**：用裸 `κ` 则 `η(t)~e^{(m₁^{dr}−κ)t}=e^{−f_dec t}` 不衰减 —— 等价于强推 `E₁=4ln(1/u)/a` 精确，**假**（`X02-THEOREMA.md:83` 的 B1 行）。

**(3) 有一处 `P3-DRESSED-MASS` 没证的缝（见 §3 修补 R3）**：`:71` 只给 `log Z_dec − f_dec t = O(1)`，`:116` 只给 `|log D_fin(T)|≤Cu`，**都是"有界"不是"收敛"**。步 4 要把 `D_fin(T)` 取到 `T→∞` 才有 `D_∞`。补齐有两条独立路，任取其一即闭：
- **组合路（首选，`[T]` 且不依赖 `[KL]`）**：`X01-L1C.md` 的精确恒等式 `log Z_dec(B_T)=Σ_{α≥0}(T−α+1)_+A(α)=φ_dec T+log D_∞+ε_T`（`:40`–`:43`），配 `|ε_T| ≤ Bū^T`、`B=(ū/(1−ū))³`（`:79`–`:80`；`Y04-AUDITX01.md:84` 逐字复核 `(★)`、`:92`–`:94` 数值 `ε_T/ū^T` 在 `T=0..12` 恒为 `2.329267`、`:128` 判 `[T]`）。⇒ `D_∞` 存在（`:47` "右端是绝对收敛级数，不先取极限"），`|log D_∞| ≤ ū/(1−ū)²=O(u)`（同 `:47`）。
- **谱路（用 `[KL]`）**：谱侧 `e^{E₁t}C_L(t)→f₁²` 必收敛（`P3-CYLINDER-RESIDUE.md:78`–`:83`）；与 (1)(2) 比对即强迫 `D_fin(T)` 收敛。

**`:94` 的 `δ ≥ κ/3` 不能用（`W05-TDECAY.md:66` `[NG]`）**：单壁平行装饰 `|P|=4T+ℓ`（`ℓ` 可到 `T`）⇒ `𝔓_rest` 率下确界 `=κ`，不是 `1.5κ`；`:69` 逐字 "`δ` 与那个 `Cu` 是**同一族装饰在分割点两侧的两种分配**… `δ` 不是分解不变量"；`:72` no-go：该路线只到 `Cu+Be^{−δt}`。**装配里只许写**：分割固定后 `∃δ>0`（`|R^{hi}| ≤ B e^{−(m₁^{dr}+δ)t}`），**不许写 `δ≥κ/3`**。`I1` 不受损（§3 (b)）。

**`a₀` 的正性从哪来（步 3(d) 的真内容，步 4 会再用）**

`a₀(β) = (2N)^{-p}·s`，`s ∈ ℤ_{≥1}`（`DRESSED:36`）；`s` = 链路重耦合多重数之积（`DRESSED:30`–`:32` Fact 0.1：`∫_{SU(N)}dU∏_iχ_{r_i}(U)=\dim\mathrm{Hom}(1,⊗_ir_i)=mult(1∈⊗_ir_i)∈ℤ_{≥0}`），其中端点 `f⊗f̄→1` 的多重数 **恰为 1**（`DRESSED:16`、`:119` 逐字 "singlet recoupling has multiplicity `1>0`"）。⇒ `a₀ ≥ (2N)^{-2}`，`T`-、`L`-无关。

**接口**：本步输出 `A^{dressed}(L) ≥ ¼(2N)^{-2}` 这个**数**，且它是同一个 `C_L(t)` 的领头指数的**系数**；把它读成"谱原子重量"是步 4/5 的事，**不是本步的事**。

*等级：`[T]`（(a)(b)(d) 已证）；`[T] mod [KL]`（(c)）；`[NG]`（`:94` 的 `κ/3`）。*

---

### 步 4 · §4 谱桥与振幅地板：`m₁` 怎么进 `f₁² ≥ r₀`、`a₀` 的 `s≠0` 从哪来

**(i) `m₁` 进 `f₁²` 的方式：它**不进**不等式，只**定位**被读的原子。**

- `m₁` 的全部作用：把管子重量写成 `w_tube=a₀e^{−m₁t}`（`FOURWALL:86`）⇒ 装饰后 `m₁^{dr}=κ−f_dec`（`DRESSED:88`）。
- 谱侧：`C_L(t)=∫_0^∞e^{−st}dν_L(s)=Σ_{n≥1}f_{n,L}²e^{−E_{n,L}t}`，`ν_L≥0`（`P3-CYLINDER-RESIDUE.md:67`–`:69`），`m_{1,L}:=E_{1,L}=\inf\mathrm{supp}\,ν_L`，`f₁²_L:=ν_L({m_{1,L}})`（`:73`–`:75`）。
- 桥：两个**精确**表示是同一个 `C_L(t)`（KP 绝对收敛 ⇒ 聚合物展开是精确表示，`CYL-RES §2.3` 第 3 条 `:140`）⇒ 领头指数必须相等、系数必须相等。`CYL-RES:140` 逐字："Since it computes the *same* function as the exact spectral sum (a), its leading exponential must be `E_{1,L}` and its coefficient must be `f₁²_L`."
- **故 `m₁` 只提供"读哪一个系数"的坐标；`f₁² ≥ r₀` 的正性来自 `a₀>0`（代数/重耦合）+ `D_∞≥½`（绝对收敛），与 `m₁` 的**数值**无关**（换道只换 `w(O)` 与 `a₀(O)`，`Y05-P3GEN.md:51`、`:60`–`:69`）。这也是 `factor-4` 只重标定质量、不动正性的原因（`FOURWALL:98` 逐字 §2.3）。

**(ii) 原子抽取（Lemma 3.4，`DRESSED:99`–`:101`，`[T]`）**

> *Let `ν ≥ 0` be a finite Borel measure on `[0,∞)` and `C(t):=∫e^{−Et}dν(E)`. If `C(t)e^{mt}→A` then `dν({m})=A`, `ν` has no mass on `[0,m)`, and `A≥0`.*

证明只用：`ν≥0`（否则 `C e^{mt}` 可任意快衰减 ⇒ 无下界结论）+ 控制收敛（`ν'({m})=0`，`1∈L¹(ν')`）。**这是 `ν≥0` 在 `I1` 里的第一处、也是唯一一处用法**（见步 5）。

代入 `C=C_L`、`m=m₁^{dr}`、`A=a₀D_∞`：`f₁²_L = ν_L({m₁^{dr}}) = a₀D_∞`（`DRESSED:107` (3.5)，`[T] mod [KL]`）。

**(iii) `a₀=(2N)^{-p}s` 的 `s≠0` 从哪来 —— 三层，逐层给行号**

| 层 | 内容 | 出处 | 等级 |
|---|---|---|---|
| L1 | `s` 是**非负整数**：`∫dU∏χ_{r_i}(U)=mult(1∈⊗_ir_i)∈ℤ_{≥0}` | `DRESSED:30`–`:32`（Fact 0.1）；`:36`（`s∈ℤ_{≥1}`） | `[T]` |
| L2 | **端点**单态 `f⊗f̄→1` 的多重数 **=1**（Schur；`f̄≅f*`） | `DRESSED:16`、`:119` | `[T]` |
| L3 | **四壁几何**强迫四条源边各有配对（否则 `∫dU U=0`）⇒ 参与收缩的那几个链路多重数 `≥1` | `FOURWALL:28`–`:34`（(1.1)/(1.1′)） | `[T]` |

⇒ `s ≥ 1 > 0`，`a₀ ≥ (2N)^{-2} > 0`。**注意两点，不许含糊**：
1. `s>0` **不是**来自 `c_r≥0`（`DRESSED:16` 明写不用；`DRESSED:117` 亦明写不用），只来自**整数重耦合多重数**与端点单态；
2. `s` **可以等于 0** —— 对别的 `O` 会给 `a₀=0`、地板消失（显式反例 `O_Δ:=Σ_i(∂_i^{lat})²(¼trF²)`，`Y05-P3GEN.md:48` `[T]`）。对 `O=¼trF²`，L2+L3 保证了 `s≥1`。这就是本定理**依赖 `O` 的具体形式**的唯一硬槽位（连同 `w=4`）。

**(iv) 组装出目标不等式**

```
f₁²_L = ν_L({m₁^{dr}}) = a₀D_∞ ≥ a₀(1−Cu) ≥ ½a₀ ≥ ½(2N)^{-2} > 0      （逐 L，用 D_∞≥½）
inf_L f₁²_L ≥ ¼(2N)^{-2} =: r₀(β) > 0                                  （L-一致，用 D_∞ 的 L-一致 + a₀ 的 L-无关）
```

r₀ 的两种口径（`[M]`，见修补 R6）：`PROOF-P3-UNCONDITIONAL.md:31` 写 `r₀(β)=½c₁u^p`，`DRESSED:124`/`X02:111` 写 `r₀=¼(2N)^{-2}`。二者都 `>0`，**取一个进定理**即可，不得并列。

**(v) 反酸测（`[O]`，章程 §1）**：自由无质量场 `C(t)≍t^{−2}`、无原子、`a₀=0`（无单态重耦合）⇒ 构造给 `f₁²=0`，**不 vacuous**（`DRESSED:127` 逐字 §4.4 "Acid test PASSED"）。

*等级：`[T] mod [KL]`。*

---

### 步 5 · `[KL]`：`ν ≥ 0`（`P3-CYLINDER-RESIDUE:67–69`）用在哪一步

**逐字（`P3-CYLINDER-RESIDUE.md:67`–`:69`）**

> ```
> C_L(t) = ∫_0^∞ e^{−s t} dν_L(s) = Σ_{n≥1} f_{n,L}² e^{−E_{n,L} t},   ν_L := Σ_{n≥1} f_{n,L}² δ_{E_{n,L}} ,  ν_L ≥ 0.
> ```
> "`f_{n,L}²≥0` is the spectral theorem …; `ν_L≥0` is Osterwalder–Schrader reflection positivity."

**三处用法，逐处定位**

| # | 用在哪一步 | 为什么非它不可 | 出处 |
|---|---|---|---|
| **U1** | **步 4(ii) 的 Lemma 3.4**：由 `C_L(t)e^{m t}→A` 反推 `dν({m})=A`、`ν([0,m))=0` | 证明第一句 "Boundedness of `C(t)e^{mt}` forces no mass below `m` (else `C(t)e^{mt}≥δe^{εt}→∞`)" 的反证只在 `ν≥0` 下成立；符号测度可任意快衰减（`REPORT.md:110` 逐字，经 `Y01-KL.md:77` 转引） | `DRESSED:99`–`:101` |
| **U2** | **步 4(iv) 的"模方"读法**：`f₁²_L=ν_L({m₁^{dr}})` 是**非负权重**（退化时是同级各态之和） | `f_{n,L}²=|⟨n,L|O|Ω_L⟩|²≥0` 是谱定理；`ν_L({m_{1,L}}) ≥ |⟨G_L|O|Ω_L⟩|²` 对任一该级本征态成立 | `CYL-RES:67`、`:76` |
| **U3** | **步 4(i) 的"原子 vs 阈"**：有限 `L` 谱**纯离散** ⇒ 每个谱点都是原子，故领头指数**必然**是最轻原子 | 紧核算子（`X_L` 紧 + 连续核 ⇒ HS ⇒ 紧）⇒ 无连续谱 | `CYL-RES:52`–`:57` |

**不用在哪（必须并列，防越级）**：
- **不用 Theorem A 的 `⇐`**。`X02-THEOREMA.md:110`、`Y01-KL.md:77` 把 `[KL]` 定位在 **Theorem A 的 `⇐`**（`IsolatedAtom ⇔ SharpDecay` 的"有率⇒有隔离"），那是 **`I2`** 的环；`I1` 不需要它（`CYL-RES:168` 逐字 "The residue theorem (I1) stands without `I2`"；`X02:113` 逐字 "`f₁²=ν({m₁})` **不是**假设，是 Theorem A `⇐` 的**产出**" —— 对 `I1` 我们用**较弱**的 Lemma 3.4 即可，不必借 Theorem A）。
- **不用 `L=∞` 的 `ν≥0`**。`[KL]` 的载体是**有限 `L`**（`CYL-RES:12`、`Y01:108` 限定 2）；把 `[KL]` 用在 `L=∞`（"原子不是连续阈、survives `L→∞`"）**恰是 `I2` 本身** ⇒ 会循环。
- **`[KL]` 的地位**：`[T]`（OS'78 + 有限 `L` 紧性 + 谱定理），既非框架假设、也非未解缺口（`Y01-KL.md:55`、`:94` 逐字）；它只是"当时未写进 P3 文件"（`DRESSED:133` 逐字 "Everything in §1–§4 is `[T]` **except one line**"）。

*等级：`[T]`（`P3-CYLINDER-RESIDUE.md:67–69`）；它是本装配的**唯一命名输入**，且 `[KL]=[T]`。*

---

## §3 死线判定

### (a) **触发** —— 四步能拼成：完整链（标 `[T] mod [KL]`）

```
S0 全表示 KP  uK<1, K≤6.1×10²N（逐 plaquette ⇒ L,T-无关）        [T]  FULLREP-KP:96–97,:120–121
S1 四壁几何   |P|≥4N_t=4t/a，等号=直管唯一                         [T]  FOURWALL:52,:54,:71
S2 领头类重量 w(B_T)=a₀u^{4T}=a₀e^{−m₁t}，m₁=4ln(1/u)/a，
   a₀=(2N)^{-p}s>0（s≥1：Fact 0.1 + f⊗f̄→1 多重数 1 + 四壁配对）  [T]  FOURWALL:86；DRESSED:34–36,:119
S3 骨架–装饰精确分解 C_L=a₀u^{4T}Z_dec(B_T)+R^{hi}_L              [T]  DRESSED:50
S4 装饰气体 KP 可和（全表示、逐片半径）                            [T]  DRESSED:63–65
S5 自由能密度/穿衣质量 m₁^{dr}=κ−f_dec，f_dec=O(u/a)，L-一致       [T]  DRESSED:70–79,:88
S6 端点缺陷：D_fin→D_∞，|log D_∞|≤O(u) ⇒ D_∞≥½（u<1/(2C)）        [T]  X01:40–43,:47,:79–80（Y04:84 复核）
S7 余项：∃δ>0，|R^{hi}_L(t)|≤B e^{−(m₁^{dr}+δ)t}（分割固定后）     [T]  DRESSED:93–95（去掉 κ/3，见 (b)）
S8 谱桥：C_L(t)=∫e^{−st}dν_L，ν_L≥0，有限 L 谱离散               [T]  CYL-RES:52–57,:67–69   ← [KL]
S9 原子从指数（用 ν≥0）：dν({m₁^{dr}})=a₀D_∞                       [T]  DRESSED:99–101
S10 地板：f₁²_L=a₀D_∞ ≥ ½(2N)^{-2} ⇒ inf_L f₁²_L ≥ ¼(2N)^{-2}>0   [T]  DRESSED:122–124
⇒ I1：在 𝒞_L=T³_L×_aℤ 上，β<β₁，每个有限 L，f₁²_L ≥ r₀(β)>0，L-一致
   【T】mod【KL】（且 [KL]=[T]，Y01:55）
```

**为何这是"装配"而不是"再证一遍"**：每步都是**已入库的**结论，本文件只做三件事——(i) 判定接口相容（下面三张票）；(ii) 补齐 `S6` 这条 `P3-DRESSED-MASS` **没写**的缝（从 `X01` 进口）；(iii) 撤掉 `S7` 里那条**假**的 `κ/3`，换成 `∃δ>0`。

**接口相容性三票（任务点名的怀疑处）**

| 接口 | 判 | 依据 |
|---|---|---|
| 步 1 的 `u<1` ↔ 步 2 的 `uK<1` | **相容**（后者更强） | `Σ_x u^{|x|}=O(1)` 只需 `u<1`（`DRESSED:38`）；`uK<1 ⇒ u<1`（`K≥1`） |
| 步 2 的 `uK<1` ↔ `X01` 的 `ū<1` | **相容**（`uK<1 ⇒ ū<1`，见步 2 的算术：`ū≤0.089/N`） | `X01:71,:75`；`Y04:152`；`FULLREP-KP:120` |
| 步 3 的 `u<1/(2C)` ↔ 步 2 的 `uK<1` | **相容**（同一 `C~K` 类常数，缩 `β₁` 即可） | `DRESSED:123` 逐字 "guaranteed on `β<β₁` (shrink `β₁` by an `N`-dependent constant if needed)" |
| 步 3(s7) 的 `δ≥κ/3` ↔ 步 4 的原子抽取 | **不相容，但可替换** | `W05:66` `[NG]`；替换为 `∃δ>0` 后 `I1` 不需更强（`CYL-RES:151`：系数是值、隙是另一对象） |

### (b) **触发 —— 但断点不在 `I1` 的链上**

**断点（真缺口，`[NG]`）**：`P3-DRESSED-MASS.md:94` 的 `δ ≥ ½κ−O(u/a) ≥ κ/3`，以及 `:51`/`:55` 的 `𝔓_rest ⇒ |P|≥6T+O(1)`。`W05-TDECAY.md:66` 逐字：一条**单壁**平行装饰（每片 1 个 plaquette、`dr=ℓ`）给 `|P|=4T+ℓ`，`ℓ` 可到 `T` ⇒ "**`𝔓_rest` 的率下确界 `=κ`，不是 `1.5κ`**" ⇒ "`:94` 的 `δ ≥ κ/3` **不成立**（勘误级）"；`:69` 逐字 "`δ` 不是分解不变量"；`:72` no-go "该路线对 `δ` 没有一致正下界"。

**为何它不进 `I1` 的链（三条，逐条给出处）**
1. `I1` 要的是**权重**、不是**间距**：`CYL-RES:151` 逐字 "the `O(u)` bound fixes the **coefficient** `f₁²_L` … but is **compatible with** `δ_L→0` … the coefficient is a value, the gap is a separate object"；`:168` 逐字 "The residue theorem (I1) stands without `I2`"。
2. 有限 `L` 的原子抽取只用 `δ_L>0`**（`E_{2,L}>E_{1,L}` 由离散性给），**不用**一致 `δ`（`CYL-RES:128` 逐字 "(Note: `δ_L>0` at each finite `L` is all (b) uses; *uniformity* … is **not** used and is **not** available — that is `I2`)"）。
3. 聚合物侧的余项只需 `e^{m₁^{dr}t}R^{hi}→0`；分割固定后 `δ>0`（`≥½κ−O(u/a)`，因为把 `o(T)`-超额族并入 `f_dec` 后，剩者只有周长 6 `1.5κ` 与双管 `2κ`；`f_dec=O(u/a)` 由 `S4` 的绝对收敛给），**不需要 `κ/3` 这个数值**。

⇒ **判**：`I1` 的链上**没有断点**；`δ` 的 no-go 属 `I2`（两指数/一致隙）那一环，与本装配正交。**但装配文本里若照抄 `δ≥κ/3`，就把一条假命题标成了 `[T]`** ⇒ 必须改（修补 R5）。

### (c) **触发 —— 三处引错/过期引用（逐条给行号）**

| # | 引错 | 逐字 | 正确处理 |
|---|---|---|---|
| **C1** | `P3-DRESSED-MASS.md:28` 把 `u` 定成 `u:=c_f(β)∈(0,1)`、`u≍β` | 与**同文件** `:35` 的 `w(B_T)=a₀u^{4T}`（要求 `a₀` `T`-无关）**不相容**：若 `u=c_f`，`a₀` 须携带 `N^{−4T}` | 取 `u:=c_f/N=β/(2N²)+O(β³)`（`FOURWALL:81,:159`；`X01:15–16` `[T]` 勘误；数值交叉 `X01:16` 引 `BETA-CONVENTION:187`）。**这条不改则 `m₁` 差 `4ln N/a`、`a₀` 的 `T`-无关性坏** ⇒ 装配里**逐字采用 `u=c_f/N`** |
| **C2** | 任务把 `ū<1` 称 "Seiler 的 KP 判据" | `ū` 是 `X01:71` 的**自制配对活动度** `ū=eΔ\bar C u`；Seiler 在 `FULLREP-KP:54` 只是 `[O]` 级**先例**（"not quoted verbatim from a numbered formula"） | P3 链的判据是 `uK<1`（`FULLREP-KP:121` `[T]`）；`ū<1` 是它蕴含的**弱条件**，用于 `X01` 的 `(L1-c)`，**不是** Seiler 的、也**不是** P3 的判据 |
| **C3** | `P3-DRESSED-MASS.md:117` 说 "the (audit-flagged, not-established for `N≥3`) blanket `c_r(β)≥0`" | 该缺口**已不存在**：`P3-FULLREP-KP.md:44` 逐字 "`c_r(β) ≥ 0` … *This also settles the small loose end … `c_r≥0` for all r is in fact a by-product of `(★)`*" `[T]` | 保留 `:117` 的**结论**（地板不用 `c_r≥0`，只用绝对收敛）但删掉"not-established"的时态——否则装配会凭空多带一条已解决的依赖 |

**另附一条 `[M]`（非引错，是数值口径）**：`β₁=min(1/48, c/N)` 的**显式常数**问题——见修补 R4。

---

## §4 修补清单（回填 `PROOF-P3-UNCONDITIONAL.md §3–§4` 时必须同时做）

| # | 修补 | 位置 | 等级 |
|---|---|---|---|
| **R1** | `u := c_f/N = β/(2N²)+O(β³)`（**不写** `u=c_f`） | `PROOF §3` 起首；替换 `DRESSED:28` | `[T]` 记号 |
| **R2** | `A^{dressed} := a₀D_∞`（**不写** `A^{dressed}=a₀`）；`f₁²_L=a₀D_∞` | `PROOF §4`；`DRESSED:87,:106` | `[M]`（`X02:110` 已判记号级） |
| **R3** | 在 `§4` 内**显式进口** `(L1-c.1)+(L1-c.2)`（`X01:40–43,:79–80`）以取得 `D_fin→D_∞`；`P3-DRESSED-MASS` 只有"有界"（`:71,:116`） | `PROOF §4` 新一行 | `[T]` 进口 |
| **R4** | `β₁`：写"`β₁=min(β₁^{KP},1/48)`，`β₁^{KP}` 由 `λ_aΦ(a)≤½` 显式给出（literal 归一 `β≲3.4×10^{-3}≈1/300`）；`c/N` 形状只在 Hooft 归一 `u≍β` 下成立" | `PROOF` 定理陈述 | `[M]`（`FULLREP-KP:124`） |
| **R5** | **删** `δ≥κ/3`；`§4` 余项写 `∃δ>0`（分割固定后），并注明 `δ` 非分解不变量、`δ≥κ/3` 已在 `W05:66` 判假 | `PROOF §3`(=`DRESSED:94`)、`§4` | `[NG]` 撤销 |
| **R6** | `r₀`：只留一个定义（建议 `r₀(β):=½a₀ ≥ ½(2N)^{-2}`，`inf_L` 再打对折 ⇒ `¼(2N)^{-2}`），不得与 `½c₁u^p` 并列 | `PROOF §5`；`MASTER:38` | `[M]` |
| **R7** | 平台记法 `T⁴_L → 𝒞_L=T³_L×_aℤ`；`t→∞` 先于 `inf_L`；`L→∞` 不取 | 全文 | `[M]`（`Y06:149`、`CYL-RES:12`） |

---

## §5 本装配**不**主张的（逐条，防加码）

1. **不改总册**：`MASTER:42` 的 `[条件/待重证]` 仍有效，直到 `PROOF-P3-UNCONDITIONAL.md` 与 `MASTER` 同轮按 §4 回填（`Y06:146`、`:140`）。本文件只提供**回填文本**，不是回填本身。
2. **不给 `I2`**：一致两指数 / `Λ★−m₁²≥δ` 仍未证；`δ` 连一致下界都没有（`W05:72` `[NG]`）。`I2` 走 `X01→X02→Theorem A` 那条（`X01:112–117`），与本装配正交。
3. **不给 Clay**：`a→0`、`β→2.5`（解析性墙 `1/48`）、`β_c≈2.30` 交叉点全未碰（`PROOF:32`）。Clay 0%。
4. **不给"每一道"**：`O=¼trF²` 这一道；换 `O` 须重填两个有限槽位 `(w(O), a₀(O)≠0)`（`Y05:60–69`），且 `0⁻⁺` 的 `a₀` 可能对消（`Y05:44–45` `[R]`）。
5. **不主张 `m₁` 精确等于 `4ln(1/u)/a`**：只有穿衣值 `m₁^{dr}=κ−f_dec`，且 `f_dec` 只证到 `O(u/a)`（`DRESSED:88`、`:156` 第 2 条）。

---

## §6 逐字索引（装配文本可直接引的行号）

| 断言 | 文件:行 | 等级 |
|---|---|---|
| `\|P\| ≥ 4N_t = 4t/a` | `P3-FOURWALL-DERIVATION.md:52` | `[T]` |
| 4 条源边各需配对（`∫dU U=0`） | `:30` | `[T]` |
| `≥4` 独立股（两片至多共 1 链） | `:34` | `[T]` |
| 只 temporal 推进时间 ⇒ 每股 `≥N_t` | `:42` | `[T]` |
| 等号 ⇒ 直四壁管唯一 | `:54` | `[T]` |
| 非管 `≥4t/a+1` | `:71` | `[T]` |
| `w_tube=a₀u^{4N_t}`、`m₁=4ln(1/u)/a` | `:86` | `[T]` |
| `u=β/2N²` ⇒ `m₁=M_{0⁺⁺}`（Gupta 14.21） | `:159`–`:161` | `[T]/[O]` |
| `c_rdim r ≤ β^{\|r\|}/\|r\|!e^β` | `P3-FULLREP-KP.md:29` | `[T]` |
| `c_r≥0`（副产品） | `:44` | `[T]` |
| `\|I(P)\|≤∏dim` | `:72` | `[T]` |
| `K≤6.1×10²N`、`uK<1`、逐片 ⇒ `L,T`-无关 | `:120`–`:121` | `[T]` |
| `β₁` 的归一口径 | `:124` | `[M]` |
| `ū:=eΔ\bar C u`、`ū<1` 不新增限制 | `X01-L1C.md:71`,`:75` | `[T]` |
| `(L1-c.1)` 恒等式 | `X01:40`–`:43` | `[T]` |
| `\|ε_T\|≤Bū^T`、`B=(ū/(1−ū))³`、无 `T` 前因子 | `X01:79`–`:80`；`Y04:84`,`:92`–`:94`,`:128` | `[T]` |
| `u:=c_f/N` 勘误 | `X01:16` | `[T]` |
| (a)(b)(d) 判词 | `P3-DRESSED-MASS.md:13`,`:14`,`:16` | `[T]` |
| KP 逐片、`uK<1`、`K≤C'N` | `:63`–`:65` | `[T]` |
| `φ_dec=O(u)`、`f_dec=O(u/a)` | `:70`–`:75` | `[T]` |
| `L`-一致 `O(e^{−cL/a})` | `:78`–`:79` | `[T]` |
| `m₁^{dr}=κ−f_dec` | `:88` | `[T]` |
| 精确分解 (1.2) | `:50` | `[T]` |
| `\|log D_fin\|≤Cu` | `:116` | `[T]` |
| `a₀=(2N)^{-p}s`、`s∈ℤ_{≥1}` | `:36` | `[T]` |
| `s>0` 之源：`f⊗f̄→1` 多重数 1 | `:16`,`:119` | `[T]` |
| Fact 0.1 多重数 `∈ℤ_{≥0}` | `:30`–`:32` | `[T]` |
| `inf_L A^{dressed}≥¼(2N)^{-2}` | `:124` | `[T]` |
| 自由场酸测 | `:127` | `[O]` |
| Lemma 3.4 原子从指数 | `:99`–`:101` | `[T]` |
| `[KL]` 命名输入 | `:135`–`:137` | `[T]`（载体在 CYL-RES） |
| 有限 `L` 谱离散 | `P3-CYLINDER-RESIDUE.md:52`–`:57` | `[T]` |
| KL 表示、`ν_L≥0` | `:67`–`:69` | `[T]` |
| 残数定义（谱式 = 抽取式） | `:73`–`:83` | `[T]` |
| 领头管 = 领头原子系数 | `:140` | `[T]`（几何） |
| 系数是值、隙是另一对象 | `:151` | `[T]` |
| `I1` 不需 `I2` | `:168` | `[T]/[C]` |
| `δ≥κ/3` **不成立** | `W05-TDECAY.md:66`,`:69`,`:72` | `[NG]` |
| `A^{dressed}=a₀` 应改 `a₀D_∞` | `X02-THEOREMA.md:110` | `[M]` |
| `¼(2N)^{-2}` 是 `r₀` 非 `a₀` | `X02:111` | `[T]` |
| `[KL]=[T]`（非假设非缺口） | `Y01-KL.md:55`,`:94` | `[T]` |
| 有限 `L` 限定（`L=∞` 会循环） | `Y01:108` | `[T]/[R]` |
| `O_Δ` 反例（`a₀=0` 可能） | `Y05-P3GEN.md:48` | `[T]` |
| 三缺口归属与等级 | `Y06-P3DOWN.md:97`–`:99`,`:131`–`:135` | `[J]` |

---

## §7 检索执行记录（`AGENTS.md §1` 四要素）

| # | 工具 | 工具状态 | 实际检索式 / 读取 | 命中 |
|---|---|---|---|---|
| 1 | `glob` | OK | `**/{Y06-P3DOWN,…,P3-CYLINDER-RESIDUE}.md` | 8（`P3-FOURWALL*.md` **未命中** ⇒ 见 #2） |
| 2 | `pwsh Get-ChildItem -Filter *FOURWALL*` | OK | `D:\Dev\math` 递归 | **1**（真名 `P3-FOURWALL-DERIVATION.md`，非 `P3-FOURWALL.md`） |
| 3 | `pwsh Get-ChildItem research\v3 -Filter *.md` | OK | 列 species 目录 | 232 |
| 4 | `read` | OK | `Y06-P3DOWN.md` 全 153 行；`P3-FOURWALL-DERIVATION.md` 全 202 行；`P3-DRESSED-MASS.md` 全 174 行；`P3-CYLINDER-RESIDUE.md` 全 186 行；`P3-FULLREP-KP.md` 全 136 行；`PROOF-P3-UNCONDITIONAL.md` 全 36 行；`X01-L1C.md` 全 157 行；`X02-THEOREMA.md` 全 195 行；`Y01-KL.md` 全 124 行；`Y04-AUDITX01.md:1–80`；`Y05-P3GEN.md:1–97`；`W05-TDECAY.md:55–75` | — |
| 5 | `grep` | OK | `ū\|bar\{u\}\|ū\|\\bar u` @ `research/**/*.md` | **56**（`ū` 定义处 `X01:71`、`Y04:152` 逐字取出） |
| 6 | `pwsh Get-Content docs/YANG-MILLS-MASTER.md` | OK | `[29..56]` 带行号 | `:35`,`:40`,`:42` 逐字 |

**无 `TOOL-FAIL`。** 未逐字核（沿用上游，未新查）：Gupta `arXiv:hep-lat/9807028 §14.2.2`（经 `FOURWALL:148–153` 转引，标 `[O]`）；Münster PLB 95/NPB 190（同上 `:179`）；`P3-AUDIT-CORRECTED.md`、`P3-FULLREP-KP` 以外的审计原文（经 `Y06` 转引）。

---

*Refs.* 内部：`P3-FOURWALL-DERIVATION.md`；`P3-FULLREP-KP.md`；`P3-DRESSED-MASS.md`；`P3-CYLINDER-RESIDUE.md`；`PROOF-P3-UNCONDITIONAL.md`；`X01-L1C.md`；`X02-THEOREMA.md`；`Y01-KL.md`；`Y04-AUDITX01.md`；`Y05-P3GEN.md`；`Y06-P3DOWN.md`；`W05-TDECAY.md`。外部（仅转引）：Gupta hep-lat/9807028 §14.2.2 (14.17)–(14.21)；Münster PLB 95 (1980) 59 / NPB 190 (1981) 439；Osterwalder–Seiler Ann. Phys. 110 (1978) 440；Seiler LNP 159 (1982)；Kotecký–Preiss CMP 103 (1986) 491；Källén (1952) / Lehmann (1954)；Reed–Simon IV §XIII。

*Status：**装配完成 ⇒ 证据层 `I1` = `[T] mod [KL]`（`[KL]=[T]`）；总册层仍 `[C]`，待 §4 七条修补回填。** (a) 触发；(b) 触发但断点在 `δ≥κ/3`（`[NG]`，属 `I2`，不落 `I1`）；(c) 触发三处（C1/C2/C3）。No QED；Clay 0%。*

*本文件只写 `research/v3/species/Z01-P3ASSEMBLY.md`；未 `git commit`。*
