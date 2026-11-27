# SU(2) 弱耦合区块平均的**复源倾斜稳定性**猜想（TS）

> 第 252 轮立项书（**取代第 251 轮的 SCE 主猜想版**）。
> **这不是 Clay 的证明**。目标是**存在性机器**，不是隙机器。
> 标号：`[T]` 已证 ／ `[G]` 猜想 ／ `[Δ]` 缺失输入。

---

## §0 三个 `β` 规范与一条引用纪律

| 规范 | 作用量 | `Re Tr` 系数 |
|---|---|---|
SZZ | `S = N β_SZZ Σ_p Re Tr Q_p` | `N β_SZZ` |
Wilson | `S = (β_W/N) Σ_p Re Tr Q_p` | `β_W/N` |

$$\beta_W=N^2\beta_{\rm SZZ}\qquad\Longrightarrow\qquad \text{SU(2):}\ \beta_{\rm SZZ}\lt\tfrac1{48}\iff\beta_W\lt\tfrac1{12}=0.08333$$

| 量 | 值 |
|---|---|
目标点 | `β_W = 2.5` |
SZZ 窗（Wilson 单位） | `1/12 = 0.08333` |
缺口（`β` 坐标） | **`30×`** |

**⚠️ 引用纪律：任何 `1/48`、`55×`、`β*` 必标规范。早前 `120×` 是跨规范比较。**

---

## §1 旧靶的判决：**开**，不是「假」

**本轮**自己读文件**的结果（这是第一轮不派出而得出的判决）：**

**`research/v3/species/HZF-COV.md:55` 逐字：**

> **自由场：Wick 给 `C_q(z)=2(⟨AA⟩_{0z})² ~ g⁴|z|^{-8}`（格点单位、`dim q=4`）⇒ `O_B=O(g⁴)`，靶真**

**⇒ 自由场给 `p = 8`，而 `p=8` 时 `Σ|z|^{-8}` 收敛 ⇒ `O_B` 有界 ⇒ **靶成立**。**

**而同一行还写：**

> **2D YM：plaquette 变量精确因子化 ⇒ `C_q(z)=0 (z≠0)` ⇒ `O_B≡0`：判据恒真、判别力 0（须报）**

**⇒ 所以靶 `ϖ≤C` 在自由场上**真**、在 2D YM 上**平凡真** ⇒ **对隙零判别力**（这与 `L3` 必须接受自由场一致）。**

**而早前那条 `O_B ≍ 2π²λ ln(b̂/a)`（`uuu3` 机检）的**前提**是 `C_q(z)=λ|z|^{-4}`，即 `p = 4` **恰好边际**——那是**一个 ansatz 的后果，不是定理**。**

$$\boxed{\ L3\ \text{的成立}\iff p\gt4;\quad p=4\ \text{给 } \ln;\quad \text{可证只到 }p=2\ (\text{KL})}$$

**⇒ 所以旧靶**真假未定**，取决于 `p`。**判决表：**

| `p` | `O_B` | 靶 | 出处 |
|---|---|---|---|
`8`（自由场 Wick） | 有界 | **真** | `HZF-COV:55` `[T]` |
`4.626`（`NN4-P4` 读数，误差 `±1.49`） | 有界 | **可能真** | 读数，非定理 |
`4`（临界假设） | `≍ ln` | **假** | `uuu3` ansatz |
`2`（KL 可证） | `≍ L²` | **假** | `YY3:23` `[T]` |

---

## §2 统一临界核与对数预算（`[T]`，可立即完成）

### 2.1 单位球体积与对数系数

设 `ρ` 是 `R^4` 上选定的范数，`v_ρ := Vol{ρ(x)≤1}`。对

$$g_c^\rho(z):=(1+\rho(z))^{-4}$$

**四维壳层体积 `dV_ρ(r)=4v_ρ r^3 dr` ⇒ `r^{-4}dV_ρ = 4v_ρ\,dr/r` ⇒**

$$\boxed{\ \sum_{0\lt\rho(z)\le L}g_c^\rho(z)=4v_\rho\ln L+O(1)\ }\tag{2.1}$$

| 范数 | `v_ρ` | `α_ρ = 4v_ρ` |
|---|---|---|
欧氏 `ρ_2=\|z\|_2` | `π²/2` | **`α_2 = 2π²`** |
`ℓ¹` `ρ_1=\|z\|_1` | `2^4/4! = 2/3` | **`α_1 = 8/3`** |

$$\frac{\alpha_2}{\alpha_1}=\frac{3\pi^2}{4}=7.402203$$

### 2.2 **关键澄清（早前记账错处）**

$$\boxed{\ 7.402203\ \text{不是同一核换壳层计数的换算因子，而是两个不同径向核的系数比}\ }$$

**逐字对照：**

$$\sum_{\|z\|_1\le L}|z|_2^{-4}=2\pi^2\ln L+O(1)\qquad\text{但}\qquad \sum_{\|z\|_1\le L}\|z\|_1^{-4}=\tfrac83\ln L+O(1)$$

**⇒ 同一个 `ℓ¹` 球截断下，**被求和的核不同**，主系数不同。**

**⇒ 单纯改变壳层划分**不改变**同一齐次核的对数主系数；改变核的**角向形状**才改变它。**

**⇒ 若要两路匹配同一主系数：`(8/3)a = 2π²λ` ⇒ `a = (3π²/4)λ = 7.4022λ`。这就是 `7.40` 该出现的位置。**

### 2.3 统一的对数预算（此后只用这一个符号）

**选定一个固定的临界核 `g_c`，定义**

$$\boxed{\ \Lambda_L^{(c)}\ :=\ 1+\sup_{x\in B_L}\sum_{y\in B_L}g_c(x-y)\ =\ \alpha_c\ln L+O(1)\ }\tag{2.2}$$

**此后所有估计统一用它：**

$$O_B\le C_B\Lambda_L^{(c)},\qquad \|\Gamma_L\|_{2\to2}\le C_2\Lambda_L^{(c)},\qquad |\kappa_{L,n}|\le n!A^n\Lambda_L^{(c)}$$

**⚠️ 若选欧氏核且 `O_B` 主项确为 `2π²λ\ln L`，则取 `g_c(z)=λ(1+|z|_2)^{-4}`，此时 `Λ_L^{(c)}=2π²λ\ln L+O(1)`。**

**⚠️ 而 `4.63` 是**固定求和方案下的有限部分**：`Σ_{0<|z|≤L}|z|^{-4}=2π^2\ln L+c_{\rm reg}+o(1)`，**依赖截断与核的定义，不是不同范数间的普适不变量**。**

**⇒ ②（三 log 统一）到此**已解决**；剩的只是把 `O_B` 与家族无序化两条路线的**原始求和核**放在同一 `g_c` 定义下比对。**

---

## §3 修正：三点界是**六项**，不是三项乘 2

**第 251 轮 §4 的 `[Δ-3]` 公式**错了**。正确推导：**

记 `c_{ij}=g_c(x_i−x_j)`、`s_{ij}=g_s(x_i−x_j)`。`T_3` 的**三棵**标号生成树、每棵**两条边**、**各选一条标记**：

| 树 | 标记 `12` | 标记另一条 |
|---|---|---|
`T_1={12,13}` | `c_{12}s_{13}` | `c_{13}s_{12}` |
`T_2={12,23}` | `c_{12}s_{23}` | `c_{23}s_{12}` |
`T_3={13,23}` | `c_{13}s_{23}` | `c_{23}s_{13}` |

$$\boxed{\ |G_{L,3}|\le B^3\Big[c_{12}(s_{13}+s_{23})+c_{13}(s_{12}+s_{23})+c_{23}(s_{12}+s_{13})\Big]\ }\tag{SCE-3}$$

**⚠️ 不是三个循环项乘 2。两者一般不等价——不能保证同一条临界边对应的两个短程距离有固定大小关系。**

**⇒ 今后检验三点界必须保留全部**六个**树边标记组合。**

### 3.1 三点的算子形式（`[T]`）

$$\mathcal K_{L,3}(h_1,h_2,h_3):=\sum_{x,y,z}G_{L,3}(x,y,z)\,h_1(x)h_2(y)h_3(z)$$

**临界边用 Young：`Σ_{x,y}|h_i(x)||h_j(y)|g_c(x−y) ≤ ‖g_c‖_1‖h_i‖_2‖h_j‖_2`；短程叶边用 C–S：`Σ_x|h_i(x)|g_s(x−y) ≤ ‖h_i‖_2‖g_s‖_2`。**

$$\boxed{\ |\mathcal K_{L,3}(h_1,h_2,h_3)|\le 6B^3\|g_c\|_1\|g_s\|_2\prod_{i=1}^3\|h_i\|_2\ }\tag{3.1}$$

**⇒ `‖K_{L,3}‖_{2,2,2} ≤ C_3Λ_L^{(c)}`。**

**⇒ **好消息：三点检验不需要逐点算构型，只需多线性算子范数。****

---

## §4 主猜想 · TS（复源倾斜稳定性）

### 4.1 设定

$$L=\frac{\hat b}{a},\quad N_B=L^4,\quad \beta_0=2.5,\quad M_L(J)=\mathbb E_{\mu_{\beta_0,L}}\exp\Big(\sum_x J_x q_x^\circ\Big),\quad W_L=\log M_L$$

### 4.2 [G-TS] 主猜想

$$\boxed{\ \begin{aligned} &\text{(i)}\quad M_L(J)\ne0,\qquad \|J\|_2\lt r_0\\ &\text{(ii)}\quad \sup_{\|J\|_2\lt r_0}\big\|D^2W_L(J)\big\|_{2\to2}\ \le\ C\,\Lambda_L^{(c)} \end{aligned}\ }\tag{TS}$$

**注意 `D²W_L(0)[u,v] = Cov(Σu_xq_x, Σv_yq_y)`；对一般 `J` 它是**复源倾斜后**的二阶响应，不是原 Gibbs 测度下的协方差。**

### 4.3 TS ⟹ GF（`[T]`）

**中心化给 `W_L(0)=0`、`DW_L(0)=0`。沿直线 `tJ` 积分两次：**

$$W_L(J)=\int_0^1(1-t)\,D^2W_L(tJ)[J,J]\,dt$$

$$\boxed{\ |W_L(J)|\le\frac C2\Lambda_L^{(c)}\|J\|_2^2\ }\tag{GF-b}$$

**再由 §3/§4 的全纯论证（有限体积 `q_x` 有界 ⇒ `M_L` 整；级数在球内收敛 ⇒ `W_L` 全纯且 `M_L=e^{W_L}` 在 0 的邻域；**恒等定理** ⇒ 全球成立）：**

$$\boxed{\ M_L(J)\ne0,\quad \|J\|_2\lt r_0\ }\tag{GF-a}$$

**⇒ `[G-TS] ⟹ (GF-a) + (GF-b)`，且 Cauchy 给全部阶：**

$$|\kappa_n(\langle h,Q_L\rangle)|\ \le\ n!A^n\Lambda_L^{(c)}\|h\|_2^n,\qquad A=A(C,r_0)$$

### 4.4 块平均

$$|\kappa_n(Z_L)|\le n!A^n\frac{\Lambda_L^{(c)}}{N_B},\qquad \mathrm{Var}(Z_L)\le 2A^2\frac{\Lambda_L^{(c)}}{N_B}$$

**原点矩（Fubini，`ΣF_nz^n/n!=1/(2−e^z)`，奇点 `z=ln2`）：**

$$|m_n|\le C_F\,n!\Big(\frac{A}{\ln2}\Big)^n\frac{\Lambda_L^{(c)}}{N_B}$$

**⚠️ `|E Z^n| ≠ E|Z|^n`。**

---

## §5 三条证明路线（严格分工）

### 路线 A · 二阶临界控制（`J=0`）：**壳预算，不拟合指数**

**⚠️ 先做一次目标区分：确定数值上的有效指数 `p_eff`，与证明足以构造机器的尾界，是两件事。点态 `p` 只是建立 F2 的一条路线。**

#### A.1 真正的充分条件：dyadic 壳预算 (SB)

**固定 §2 选定的 `ρ`，并令 `g_c` 与 `ρ` 完全一致。对二进制壳定义**

$$S_j(L):=\sum_{\substack{z\in B_L\\ 2^j\le\rho(z)\lt2^{j+1}}}\big|C_L(z)\big|,\qquad \sum_{\rho(z)\le L}|C_L(z)|=\sum_{j=0}^{\lfloor\log_2L\rfloor}S_j(L)+O(1)$$

$$\boxed{\ \sup_{L,j}S_j(L)\ \lt\ \infty\ }\tag{SB}$$

**⇒ (SB) 直接给 `sup_xΣ_y|C_L(x,y)| ≤ C(1+\ln L)` ⇒ Schur ⇒**

$$\boxed{\ \|\Gamma_L\|_{2\to2}\ \le\ C\Lambda_L^{(c)}\ }\tag{F2}$$

**⭐⭐ 关键优化：** **(SB) 不需要 `p>4`。只要每层壳的绝对协方差总量一致有界，`p=4` 完全可以接受。**

#### A.2 三种情形必须被区分

**若 `|C(z)| ~ A(θ)|z|^{-p}(\ln|z|)^\sigma`，则 `S_j ≍ 2^{(4-p)j}j^\sigma`：**

| 情形 | `S_j` | SB | F2 |
|---|---|---|---|
`p\gt4` | 指数衰减 | 成立 | **成立** |
`p=4,\ \sigma=0` | 一致有界 | 成立 | **成立** |
`p=4,\ \sigma\lt0` | 衰减 | 成立 | **成立** |
**`p=4,\ \sigma\gt0`** | **多项式增长 `j^σ`** | **不成立** | **不能由 SB 推出** |
`p\lt4` | 指数增长 | 不成立 | 不能由 SB 推出 |

$$\boxed{\ \text{真正危险的不是 }p=4\text{，而是 }p=4\ \wedge\ \sigma\gt0\ }$$

**⚠️ 所以必须把**乘法对数修正 `σ` 一并拟合**，不能只拟合幂指数。**

#### A.3 拟合模型与局部指数

$$\boxed{\ \log S_j=a+(4-p)j\ln2+\sigma\ln(j+j_0)+b\,2^{-\omega j}\ }\tag{FIT}$$

**比 `log|C(z)|=-p\log|z|+a` 合理，因为它直接拟合四维壳求和后的量。**

**局部指数（不混尺度）：**

$$\boxed{\ p_{\rm eff}(j)=4-\log_2\frac{S_{j+1}}{S_j}\ }$$

**若存在对数修正：`p_eff(j)=p-\sigma\log_2(1+1/j)+o(j^{-1})`。**

**⚠️ 拟合三项控制：**
1. **用 MC 独立配置块做 bootstrap，不把相邻距离壳当独立样本；**
2. **改变拟合区间，检查 `p,σ` 稳定性；**
3. **用已知自由场 `p=8` 作基准，并用合成的 `p=4`/`p=2` 数据检查程序能否区分临界情形。**

#### A.4 本轮我自己算的（`O_B` 的壳增量，条件性诊断）

**`O_B = 14.26, 24.55, 36.27, 48.89`（`R=4,8,16,32`，`uuu3` 机检）⇒ 增量 `10.29, 11.72, 12.62`。**

$$p_{\rm eff}^{(1)}=4-\log_2\frac{11.72}{10.29}=3.8123,\qquad p_{\rm eff}^{(2)}=4-\log_2\frac{12.62}{11.72}=3.8933$$

**⇒ 两个读数**朝 `4` 靠近**。**

**若按 `S_j≍Aj^σ` 反解：`σ(2→3)=0.3209`、`σ(3→4)=0.2572` —— **在下降**。**

**而 `O_B` 的二次拟合（`t=\ln R`）：**

$$O_B\approx-3.06+10.80\,t+1.21\,t^2\qquad(t=\ln32\ \text{时二次项占线性项的}\ 39\%)$$

**⚠️ 三条限制，必须写清：**
1. **只有**三个**壳增量，不能证明 `p=4`，也不能证明 `p>4`；**
2. **`O_B` 是**有号累计和**，`S_j` 是**每壳绝对和**；只有在 `C_q(z)\gt0`（`z≠0`，`II6` `[T]`）时才相等；**
3. **`O_B` 与 `C_L` 是否同一核、同一归一化**尚未确认**，不得合并成同一组拟合数据。**

**⇒ 现有数据与「`σ\gt0` **且** `σ` 在下降」相容，也与「`σ=0` 加有限尺度修正」相容。**三壳不能判。****

#### A.5 统计门槛（若坚持用 `p_eff`）

**`p_eff = 4.626 ± 1.49`（`NN4-P4` 读数）。要在双侧 95% 使下界越过 4：**

$$\sigma_p\lt\frac{4.626-4}{1.96}=0.3194\qquad(\text{单侧 95\%：}\ \sigma_p\lt0.3805)$$

**⇒ 标准差需降 `4–5×`。**

**⚠️ 但这只是**统计**误差要求；有限尺寸偏差、模型选择偏差、对数修正须另行计入。**

**⭐⭐ 更建议**直接检验 SB**，而不是把资源全花在压窄 `4.626` 的误差条上——**因为 SB 正是所需充分条件，不必区分 `p=4` 与 `p>4`。**

#### A.6 更直接的一道检验：谱

**周期边界 + 平移不变 ⇒ Fourier 对角化，`Ĉ_L(k)≥0`（协方差半正定）⇒**

$$\boxed{\ \|\Gamma_L\|_{2\to2}=\max_k\widehat C_L(k)\ }$$

**⇒ F2 的直接数值判据：**

$$\boxed{\ R_L:=\frac{\max_k\widehat C_L(k)}{\Lambda_L^{(c)}},\qquad \sup_L R_L\lt\infty\ }\tag{SPEC}$$

**优点：不依赖 `C_L` 在实空间是否处处非负。**

**⚠️ 两处早前判断要纠正：**
1. **`χ_q=Ĉ_L(0)` 只检测**零动量**，不能代替 `max_k`**——不同动量方向的涨落可有不同的有限尺寸行为；
2. **`p=2` **不必然**使 F2 失败**——若长程尾有振荡、在动量空间抵消，需查整个谱才能下结论。

| 输入 | 给出 | 判 |
|---|---|---|
**P4**：`\|C_L(z)\|≤C(1+\|z\|)^{-4}` | ⇒ (SB) ⇒ F2 | `[Δ]`，**充分不必需** |
**(SB)** | ⇒ F2 | **`[Δ]`，比 P4 弱，优先** |
**(SPEC)** | = F2 | **`[Δ]`，数值可直接查** |
`p=8`（自由场，`HZF-COV:55`） | `‖C‖_1\lt∞` ⇒ F2 | **`[T]` 自由场上成立** |
`p=2`（KL，`YY3:23`） | `‖C‖_1≍L²` ⇒ **F2 在非负尾下假** | `[T]` |

### 路线 B · 高阶连通结构

**两层，必须分开判：**

| 层次 | 命题 | 关系 |
|---|---|---|
**A（弱）** | `\|K_{L,3}\|_{2,2,2} ≤ C_3Λ_L^{(c)}`（UCB-3） | 机器需要 |
**B（强）** | 六项 SCE-3 | `B ⟹ A`，反之不然 |

$$\boxed{\ \text{SCE-3}\Longrightarrow\text{UCB-3}\qquad\text{UCB-3}\not\Longrightarrow\text{SCE-3}\ }$$

**⇒ **若 UCB-3 成立而 SCE-3 失败，机器没失败，只是点态树分解不是正确的证明机制。****

**⇒ 这避免了把过强的空间结构假设误当成数学目标本身。**

### 路线 C · Ward / Schwinger–Dyson

**目标型估计（不是现有恒等式的结论）：**

$$\big\|D_J\mathcal H_L(J)[h]\big\|\ \le\ C_1\Lambda_L^{(c)}\frac{\|h\|_2}{1-A_1\|J\|_2},\qquad \mathcal H_L:=D^2W_L$$

**⇒ 沿路径积分即得固定复源球内的响应控制。**

**SD 形式：`E[D_ℓF]=E[F\,D_ℓS_W]`（Haar 群方向分部积分）。真正任务不是写恒等式，而是证明**右侧插入算子在 `ℓ²` 源空间中一致有界**，并处理接触项、链路重叠、规范不变性约束。**

---

## §5′ 三个命题的真逻辑关系与 RTC 闭合引理

### 5′.1 三个命题

$$\begin{aligned}
\text{F2}:&\quad \|H_L(0)\|_{2\to2}\le C_0\Lambda_L^{(c)}\\
\text{TF2}:&\quad \sup_{h\in\mathbb R^{N_B},\ \|h\|_2\lt r_0}\|H_L(h)\|_{2\to2}\le C_1\Lambda_L^{(c)}\\
\text{TS}:&\quad M_L(J)\ne0\ \text{且}\ \sup_{\|J\|_2\lt r_0}\|H_L(J)\|_{2\to2}\le C_2\Lambda_L^{(c)},\quad J\in\mathbb C^{N_B}
\end{aligned}$$

**其中 `H_L(J) := D²W_L(J)`。显然 `TS ⟹ TF2 ⟹ F2`。反向皆有**严格反例**。**

### 5′.2 反例 A：`F2 ⇏ TF2` `[T]`

**取 `N=N_B`、`Λ=1+\ln L`、`p_N=Λ/N`，令 `q_1=…=q_N=X`，`P(X=1)=p_N`。**

$$\Gamma_L=p_N(1-p_N)\mathbf 1\mathbf 1^{\mathsf T},\qquad \|\Gamma_L\|=Np_N(1-p_N)\le\Lambda$$

**⇒ F2 成立。而取 `h_N=t_Ne`、`t_N=N^{-1/2}\ln\frac{1-p_N}{p_N}`：**

$$t_N\to0\ (\text{故}\ \|h_N\|\to0),\qquad P_{\mu_{L,h_N}}(X=1)=\tfrac12\ \Longrightarrow\ \|H_L(h_N)\|=\frac N4$$

**数值核（我算）：**

| `N` | `Λ` | `‖H(0)‖` | `‖H(h)‖` | 比 | `t_N` |
|---|---|---|---|---|---|
`100` | `5` | `4.75` | `25` | `5.0` | `0.294` |
`1000` | `8` | `7.94` | `250` | `31.2` | `0.152` |
`10^4` | `10` | `9.99` | `2500` | `250` | `0.069` |

$$\boxed{\ \text{F2}\not\Longrightarrow\text{TF2}\ }$$

**机制：极少数宏观配置（概率 `Λ/N`）可被**极小的实源**提升为 `50–50` 混合。**

### 5′.3 反例 B：`TF2 ⇏ TS` `[T]`

**取 `S=±1` 各 `1/2`，`Q_L=\sqrt\Lambda\,S\,e`（故 `q_x=\sqrt{\Lambda/N}S` 一致有界）。**

$$M_L(J)=\cosh\big(\sqrt\Lambda\langle J,e\rangle\big),\qquad H_L(h)=\Lambda\,\mathrm{sech}^2\big(\sqrt\Lambda\langle h,e\rangle\big)\,e\otimes e$$

$$\sup_{h\in\mathbb R^N}\|H_L(h)\|=\Lambda\qquad\Longrightarrow\qquad\text{TF2 对【一切】实源成立}$$

**但零点满足 `\sqrt\Lambda\langle J,e\rangle=i\pi(k+\tfrac12)` ⇒ 最近零点距**

$$r_L=\frac{\pi}{2\sqrt{\Lambda_L}}\ \longrightarrow\ 0$$

$$\boxed{\ \text{TF2}\not\Longrightarrow\text{TS}\ }$$

**⚠️ 两例都是**抽象有界随机变量模型**，不是 SU(2) YM 的反例。它们只严格证明：不能仅由 F2、TF2 或有界性推出 TS。**

### 5′.4 [G-RTC-$r_0$] 局部相对三阶响应闭合引理

**⚠️ 本节为第 255 轮**紧急修正**：早前写成的「在**整个** `Ω_L` 上」的全域版本**不成立**，已撤回。原因见 §5′.4a。**

#### 5′.4a 全域版本不可能的证明（`[T]`）

**设 `q_x` 有界实观测，存在方向 `e` 使 `Var(⟨e,Q_L⟩)>0`。记 `X=⟨e,Q_L⟩−E⟨e,Q_L⟩`、`F(z)=E[e^{zX}]`。**

**`X` 有界 ⇒ `F` 整且指数型：`|F(z)|≤e^{R_L|z|}`。】

**若 `F` 全平面无零点：Hadamard 因子分解（零自由 + 有限阶 `≤1`）⇒ `F(z)=e^{az+b}`。由 `F(0)=1`、`F'(0)=EX=0` 得 `a=b=0` ⇒ `F≡1`，与 `F''(0)=Var(X)>0` 矛盾。**

$$\boxed{\ \text{有限体积有界非退化观测量的矩生成函数【必有】复零点}\ }\tag{Z}$$

**在重数为 `m` 的零点 `z_*` 附近：**

$$F=(z-z_*)^m g\ \Longrightarrow\ W''=-\frac{m}{(z-z_*)^2}+O(1),\qquad W'''=+\frac{2m}{(z-z_*)^3}+O(1)$$

**⇒ 逼近零点时 `‖DH_L‖ ~ d^{−3}` 而 `‖H_L‖ ~ d^{−2}` ⇒ 比值 `~d^{−1}→∞`。**

$$\boxed{\ \text{任何在【整个零点分支】上带 }L\text{-一致常数的三阶相对界都不可能}\ }$$

**多变量情形同理：光滑零点的 Hessian 是二阶极点，其导数是三阶极点。**

#### 5′.4b 修正后的命题（保留为正式猜想）

**把源限制在**半径与体积无关**的球内 `B_{r_0}={J:‖J‖₂<r_0}`。设 `Ω_L=Comp_0{J:M_L(J)≠0}`。**

$$\text{(RTC-}r_0\text{)}\qquad \|DH_L(J)[u]\|_{2\to2}\ \le\ \|u\|_2\big(a\|H_L(J)\|_{2\to2}+b\Lambda_L^{(c)}\big)$$

**对所有 `J∈Ω_L∩B_{r_0}`、所有复方向 `u`、所有充分大的 `L` 成立。`a,b,r_0` 皆与 `L` 无关。**

**⚠️ 这个版本没有零点附近的矛盾，因为只要求在**统一有限源球**内成立。**

#### 5′.4c `F2 + RTC-r_0 ⟹ TS`（`[T]`，我逐行核过）

1. 路径 `J(t)=tJ`，`0≤t≤1`。**因 `‖tJ‖≤‖J‖<r_0`，整条路径在 `B_{r_0}` 内** ✓
2. `f(t)=‖H_L(tJ)‖`，`D^+f(t)≤‖J‖(af(t)+bΛ)`（RTC-r_0 在首个零点之前可用）
3. Grönwall：`f(1)≤Λ[C_0e^{ar_0}+\frac ba(e^{ar_0}−1)]`，取 `C_2` 为该值（`a=0` 取 `C_0+br_0`）
4. 若球内有首个零点 `J_*=t_*J`：`t<t_*` 时 `W_L(tJ)=∫_0^t(t−s)H_L(sJ)[J,J]ds` ⇒ `|W|≤½C_2Λt²‖J‖²`
5. `t↑t_*` 有界 ⇒ `W_L` 有限极限；`M_L` 整 ⇒ `M_L(J_*)=lim e^{W_L}≠0`，与零点矛盾 ∎

$$\boxed{\ \text{F2}+\text{RTC-}r_0\ \Longrightarrow\ \text{TS}\ }$$

**⭐ 而 `a=0` 是一个**更干净的目标**（我提）：取 `a=0` 即**

$$\boxed{\ \sup_{\|J\|_2\lt r_0}\big\|D^3W_L(J)\big\|_{2\to2}\ \le\ b\,\Lambda_L^{(c)}\ }$$

**⇒ 此时 `C_2=C_0+br_0` 有界，且难度集中到**一个常数 `b`**；而高斯测度给 `D³W≡0` ⇒ `b=0`。**

### 5′.4d 文献证据状态（`[T]`/`[U]`/撤回）

| 项 | 证据 | 状态 |
|---|---|---|
`D²W=Cov` | 指数族微分恒等式 | **`[T]`** |
`D³W=κ₃` | 指数族微分恒等式 | **`[T]`** |
实源三阶相对界 `\|D³f[u,v,v]\|≤R\|u\|D²f[v,v]` | **广义（2-）自协调性**（Bach） | **`[T]`** |
实源系数 `R_L≤2q_*\sqrt{N_B}` | 有界特征直径 | **`[T]`**，但 **`\sqrt{N_B}` 损失** |
Haar–Schwinger–Dyson 恒等式（紧群变量） | SU(N) 有限格点精确 | **`[T]`** |
**全域 RTC（原版本）** | 由 `(Z)` | **撤回** |
局部实源 RTC，常数一致 | 目标模型 | **`[U]`** |
局部复源 RTC，常数一致 | 目标模型 | **`[U]`** |
`β=2.5` 下 RTC-`r_0` | 目标模型 | **`[U]`** |
`RTC-r_0 ⟹ TS` | Grönwall + 零点延拓 | **`[T]`** |

**⇒ **唯一的核心缺口重新表述为：****

$$\boxed{\ \text{把实源估计的 }\sqrt{N_B}\ \text{损失降到 }O(1)\text{，并建立固定复源球内的统一响应界}\ }$$

**⚠️ SZZ（`|β_SZZ|<1/48`）只能作弱耦合基准：原测度的 Poincaré **不**推出所有倾斜测度的统一 Poincaré，更不推出复源零点自由。**

### 5′.4e 精确起点：Haar–SD（`[T]`）

**从 `0=∫D_ℓ^a[F·e^{−S_W+⟨J,Q⟩}]dU` 得精确恒等式**

$$\mathbb E_J\big[D_\ell^aF\big]=\mathbb E_J\Big[F\big(D_\ell^aS_W-\langle J,D_\ell^aQ\rangle\big)\Big]\tag{SD}$$

**对源求导产生带观测量插入的 SD 关系。**恒等式是精确的；未完成的是控制右侧。**

**三个可验证估计：**

| 估计 | 内容 |
|---|---|
**A** | 局域插入：`D_ℓ^a` 的精确 SD 关系 |
**B** | 插入算子的统一界：`‖I_L(J)[u]‖≤C_1‖u‖(‖H_L(J)‖+Λ)`，**须含重叠链路与接触项，光靠单片有界不够** |
**C** | 闭合：`D³W_L=A_L(J)D²W_L+R_L(J)`，`‖A_L[u]‖≤a‖u‖`、`‖R_L[u]‖≤bΛ‖u‖` |

**⚠️ 这三步目前都**没有**从 SU(2) Wilson 作用量推出。**

---

### 5′.5 RTC 的内容，以及它比反例 A 强在哪

**RTC 恰是宣称：反例 A 的机制在目标模型中**不发生**。**

$$\text{反例 A 的机制}=\text{概率}\ \Lambda/N\ \text{的罕见宏观配置被极小源激活}$$

**⇒ 所以 RTC 是关于**测度尾部**的陈述：**没有罕见宏观配置被小源放大**。**

**⇒ 这比 F2 强得多，也比 TF2 强（反例 B 显示 TF2 可对一切实源成立而零点半径仍趋零）。**

**⚠️ RTC 本身尚未在 `β=2.5` 的 SU(2) 中证明；也**尚无证据**。**

### 5′.6 RTC 让 Ward / SD 有了精确任务

**从 Haar 分部积分 `E[D_ℓF]=E[F\,D_ℓS_W]`，选带外部源的插入，试图得到**

$$D^3W_L(J)=\mathcal A_L(J)\,D^2W_L(J)+\mathcal R_L(J)$$

**需证：**

$$\boxed{\ \|\mathcal A_L(J)[u]\|\le a\|u\|_2,\qquad \|\mathcal R_L(J)[u]\|\le b\Lambda_L^{(c)}\|u\|_2\ }$$

**三个必须处理的实事：**
1. **接触项**：不同局部插入可共享链路/plaquette，不能按独立插入计；
2. **规范约束**：源插入与群方向导数须保持规范不变，不能用普通欧氏坐标的独立变量假设；
3. **复源**：估计须在复源区域成立，不能只证真实概率测度下的 Ward 恒等式。

**⇒ 这三处是目前**真正没有闭合**的地方。**

### 5′.7 与 Dobrushin–Shlosman 完全解析性理论的关系

**该理论把「复扰动下无零点」「高阶截断关联的统一控制」「带空间衰减的树图估计」联系起来，与本机器的目标高度相关。**

**⚠️ 但不能直接套用：其典型树图条件要求衰减函数**可和**，而四维 `g_c(z)\sim|z|^{-4}` 恰在**对数发散的边际点**。**

**⇒ 需要的是**临界完全解析性**：允许每个连通关联簇支付**一次** `Λ_L^{(c)}`，而不是要求全部传播核绝对可和。**

$$\boxed{\ \text{RTC 的价值正在于此：它允许二阶响应保留临界对数，而不要求消灭这根 log。}\ }$$

---

## §5″ 【主结果】链变量张成族：由局部中心对称性得到的 `[T]` 机器

> **第 256 轮。此节是整份文件唯一**完全闭合**的有限体积结果。**
> **它不需要 `F2` 猜想、不需要 `RTC` 猜想、不需要 `p=4`、不需要 `1/48` 窗、不需要 `β` 搬运、不需要 Ward/SD。**
> **代价：它只对**链变量**族成立，不对 plaquette 块平均 `q` 成立。**

### 5″.1 精确对称性与选择规则

**未做规范固定。在任意顶点 `v` 施加 `g_v=−I₂`、其余 `I₂`。**

$$\text{所有与 }v\text{ 相接的链变量同时变号：}\ U_\ell\longmapsto-U_\ell$$

**每个 plaquette 与 `v` 相接的链路数是 `0` 或 `2`（`v` 为其角点时恰 `2`）⇒ 两次变号抵消 ⇒ `S_W` 不变 ⇒ 测度不变。**`[T]`（精确，任意 `β`）**

**四元数坐标 `U_\ell=x_{\ell,0}I+iΣ_{a=1}^3x_{\ell,a}σ_a`，`Σx²=1`。`−U` ⇒ **四个坐标全变号**。**

$$\boxed{\ \mathbb E\Big[\prod_\ell\prod_j x_{\ell,j}\Big]=0\quad\text{除非}\quad \sum_{\ell\ni v}n_\ell\equiv0\ (2)\ \ \forall v\ }\tag{Z2}$$

**⇒ 非零矩的链路多重图必须是**偶子图**。**

### 5″.2 张成族的坐标数（改正）

**`Re U^{ij}, Im U^{ij}` 共 `8` 个但线性相关。改用四元数四坐标：**

$$X_\ell=(x_{\ell,0},\dots,x_{\ell,3})\in S^3,\qquad D_L=4|E_L|=16L^4$$

**⚠️ `SU(2)` 流形维数 `3`，但**全局嵌入需要 `4` 个坐标**（满足一个二次约束）⇒ 不能因维数 `3` 就删掉一个全局生成坐标。**
**四个坐标已足以由多项式生成全部连续函数（Stone–Weierstrass：分离点 + 含常数）⇒ Wilson loop 的矩阵乘积与迹都是它们的多项式。**`[T]`

### 5″.3 偶子图估计（`[T]`，逐步）

**令 `s=‖J‖₂`、`w_ℓ=Σ_a|J_{\ell a}|`。Cauchy–Schwarz：`Σ_ℓ w_ℓ² ≤ 4s²`。**

$$M_L(J)=\mathbb E\exp\Big(\sum_{\ell,a}J_{\ell a}x_{\ell a}\Big),\qquad |M_L(J)-1|\ \le\ \mathcal Z_{\rm even}(w)-1\tag{2}$$

**（用 `(Z2)` 与 `|x|≤1`；**不用任何关联衰减、不用强耦合展开**。）**

**按 `n_ℓ` 的奇偶分组（`t_ℓ=\tanh w_\ell`）：**

$$\mathcal Z_{\rm even}(w)=\prod_\ell\cosh w_\ell\ \cdot\!\!\sum_{\substack{F\subset E_L\\ \partial F=0}}\prod_{\ell\in F}t_\ell\tag{3}$$

**偶子图分解为边不相交的简单环 ⇒**

$$\sum_{F:\partial F=0}\prod_{\ell\in F}t_\ell\ \le\ \exp\Big(\sum_{C\in\mathcal C_L}\prod_{\ell\in C}t_\ell\Big)\tag{4}$$

**⚠️ 严格写时环应带 `1/|C|` 权（吸收同一个偶子图的多种环分解）；因我们要**上界**，用不带权的更大者即可。**

**加权邻接矩阵 `A_{xy}=t_{\{x,y\}}`、`R=‖A‖_F`。长 `m` 的简单环在 `tr(A^m)` 中贡献恰 `2m` 次：**

$$\sum_C\prod_{\ell\in C}t_\ell\ \le\ \frac12\sum_{m\ge3}\frac{\operatorname{tr}(A^m)}{m}\ \le\ \frac{R^3}{6(1-R)}\quad(R\lt1)\tag{5}$$

$$R^2=2\sum_\ell t_\ell^2\ \le\ 2\sum_\ell w_\ell^2\ \le\ 8s^2\qquad\Longrightarrow\qquad R\le\sqrt8\,s\tag{6}$$

$$\boxed{\ \mathcal Z_{\rm even}(w)\ \le\ \exp\Big[2s^2+\frac{(\sqrt8 s)^3}{6(1-\sqrt8 s)}\Big]\ }\tag{7}$$

**⭐⭐ 关键：`(7)` 里**没有 `L⁴`、没有 `L²`、没有 `ln L`**。原因不是更强的物理关联衰减，而是中心对称性把非零项限制在偶图上，而偶图加权和由**有限的 Frobenius 范数**控制。**

### 5″.4 显式一致解析半径（`[T]`，常数我逐位核过）

$$r_*:=\frac1{4\sqrt8}=0.088388347648$$

$$s\le r_*\ \Longrightarrow\ 2s^2+R^3/(6(1-R))\le\frac1{64}+\frac1{288}=\frac{11}{576}$$

$$|M_L(J)-1|\ \le\ e^{11/576}-1=0.019280741\ \lt\ 0.02\tag{8}$$

$$\boxed{\ M_L(J)\ne0,\qquad \|J\|_2\lt r_*\ }\tag{GF-a-link}$$

$$|W_L(J)|\ \le\ 3\|J\|_2^2,\qquad \|J\|_2\lt r_*\tag{GF-b-link}$$

**（用 `|log(1+z)|≤−log(1−|z|)`；`3` 是保守绝对常数。）**

**⚠️ 此结论对**所有有限体积 `L`**、且对**任意满足该局部中心对称性的规范不变概率测度**成立 ⇒ **`β` 可以是 `2.5`。**1

### 5″.5 `F2-link` 是**精确恒等式**（`[T]`）

$$\mathbb E[x_{\ell,a}]=0\quad(\text{中心对称})$$
$$\mathbb E[x_{\ell,a}x_{\ell',b}]=0\ (\ell\ne\ell')\quad(\text{取 }v\in\ell\setminus\ell'\text{，只 }\ell\text{ 变号})$$
$$\text{单链路边际}=\text{Haar}\ (\text{两端左右平移不变})\ \Longrightarrow\ \mathbb E[x_ax_b]=\tfrac14\delta_{ab}$$

$$\boxed{\ \Gamma^{\rm link}_L=\tfrac14 I_{16L^4},\qquad \|\Gamma^{\rm link}_L\|_{2\to2}=\frac14\ }\tag{F2-link}$$

**⇒ **精确，不是拟合；`β`-无关，`L`-无关。****

**⚠️ 必须保留的区别：`Γ^{link} ≠ Γ^q`。此结果证明的是链变量族的 F2，**不是** plaquette 场的磁化率界。**

### 5″.6 `RTC-link`，取 `a=0`（`[T]`）

**取 `‖J‖₂ ≤ r_*/2`，三个单位方向 `u,v,w`，`f(z_1,z_2,z_3)=W_L(J+z_1u+z_2v+z_3w)`。**
**当 `|z_j| < r_*/8` 时 `‖J+…‖ ≤ 7r_*/8 < r_*` ⇒ `|f| ≤ 3r_*²`。**

**三变量 Cauchy（**注意不要 `3!`**：`∂³f/∂z₁∂z₂∂z₃|₀ = 1!1!1!·a₁₁₁ = a₁₁₁ ≤ M/(ρ₁ρ₂ρ₃)`）：**

$$\boxed{\ \big|D^3W_L(J)[u,v,w]\big|\ \le\ \frac{1536}{r_*}\|u\|\|v\|\|w\|\ }\tag{RTC-link}$$

$$b=\frac{1536}{r_*}=6144\sqrt8=17377.856254,\qquad \|J\|_2\lt r_*/2=0.044194174$$

$$\boxed{\ \|DH_L(J)[u]\|_{2\to2}\ \le\ b\|u\|_2\qquad(a=0)\ }$$

### 5″.7 组装出的常数（`[T]`，我逐位算过）

| 量 | 值 |
|---|---|
`r_*` | `0.088388347648` |
`r_0 = r_*/2` | `0.044194174` |
`b = 1536/r_*` | `17377.856254` |
`C_0 = ‖Γ_link‖` | `0.25` |
**`C_2 = C_0 + b r_0`** | **`768.25`**（`L`-无关） |
`GF` | `\|W_L(J)\| ≤ 384.125\|J\|_2^2` |

**⚠️ `b` 很松（`C_2/‖Γ‖ = 3072×`），但 **`L`-一致**——关键的是一致性，不是数值小。**

### 5″.8 这条结论**不给**什么（结构性限制）

**论证的五个前提全部与 `β` 无关：**

```
① SU(2) 中心含 −I         （群论）
② plaquette 在 g_v=−I 下不变（−I 在中心，两次变号抵消）
③ 测度的规范不变性         （任意 β）
④ |x_{ℓ,a}| ≤ 1           （U ∈ SU(2)）
⑤ 偶子图 + Frobenius       （纯图论）
```

$$\boxed{\ \text{该结论对任意 }\beta\text{（含 }\beta=0\text{）成立}\ \Longrightarrow\ \text{判别力为零}\ }$$

**⇒ `β=0` 时测度是**乘积 Haar**，平凡无隙；`β=2.5` 时不是。该界不能区分两者。**

**⇒ 所以它给：`L`-一致的矩界 + **子列**极限存在。**
**⇒ 它不给：聚簇性、唯一性、隙。**

**⇒ 这是 `mmm5` 那条筛子的直接后果，**不是缺陷**——存在性不该依赖耦合强度。**

---

## §5‴ 【关键】A 路线不足以给聚簇性：一个决定性的反例

> **第 257 轮。此节把 A 路线的边界钉死，并把缺件从「唯一性」拆成三个独立子任务。**

### 5‴.1 反例：混合测度（`[T]`，我核过每一步）

$$\mu_*=\tfrac12\mu_F+\tfrac12\mu_H$$

- **`μ_H = ⊗_ℓ dU_ℓ`（独立链变量 Haar）**
- **`μ_F`：平坦规范场 `U_{xy}=g_xg_y^{-1}`，`g_x` iid Haar**

**两者都规范不变：`g_x→Λ_xg_x` 给 `U_{xy}→Λ_xU_{xy}Λ_y^{-1}`；`g_x→−Ig_x` 给同样的 `U_ℓ→−U_ℓ`。**

$$\Longrightarrow\quad \mu_*\ \text{有完全相同的局部中心对称性}$$

**⇒ 因此 `(GF)`、`F2-link`、`RTC-link` **全部**成立，甚至 `Γ^{link}=(1/4)I`。**

**而平坦场上闭合环 `U(C)=g_{x_0}g_{x_{\rm end}}^{-1}=I` ⇒ `P(C)=½Re tr U(C)=1`：**

$$\mathbb E_*P_A=\tfrac12,\qquad \mathbb E_*(P_AP_B)=\tfrac12\qquad(\text{两环边不相交})$$

$$\boxed{\ \mathrm{Cov}_{\mu_*}(P_A,P_B)=\tfrac14\qquad\text{与距离【无关】}\ }$$

$$\boxed{\ \mathrm{Cov}_{\mu_*}(q_A,q_B)=\tfrac14\ }\qquad(\text{同式，}q=1-P)$$

**⇒ 而该测度下 `χ_q ≍ N/4`，**随体积线性增长**。**

**⚠️ `μ_*` 不是 `β=2.5` 的 Wilson Gibbs 测度 ⇒ **不是对目标理论的反例**；它是对「由 `β`-无关的 A 路线推出唯一性」这一推理的严格反例。**

### 5‴.2 结构教训（这是本轮最重要的一句）

$$\text{中心对称是【单点/局部】对称}\ \Longrightarrow\ \text{只约束链变量奇次矩}$$

$$\text{聚簇性是关于【整个测度】的性质；凸混合可破坏它而保留全部局部对称}$$

$$\boxed{\ \text{任何纯局部对称都不能给聚簇性；聚簇性必须来自【权重】（}β\text{-有关的 Wilson 作用量）}\ }$$

**⇒ 于是「A 路线给存在性、另一台机器给相互作用」不是失败，而是**正确的模块分离**：**

$$\text{第一台机器不需要知道 }\beta\qquad\text{第二台机器必须看见 }\beta{=}2.5\text{ 的 Wilson 相互作用}$$

### 5‴.3 偶子图 ≠ 连通偶图

**选择规则只给 `∂F=0 (2)`。而独立链变量 Haar 测度下两个远隔的 Wilson 环**正是**偶子图。**

$$\boxed{\ \text{偶图}=\text{一个或多个闭合回路之并}\ \ne\ \text{连通闭合图}\ }$$

**⇒ 所以原来的 `(GF)` 里**没有距离** `R`：**

$$\big|\mathcal I_L(J_A,J_B)\big|\ \le\ 2C\big(\|J_A\|^2+\|J_B\|^2\big)\tag{9}$$

**⚠️ `D²W_L(0)[J_A,J_B]=0`（不同链路源）是对的，但它只控制**二阶**；高阶项仍允许两个不相连的闭合回路。**

### 5‴.4 新主任务：连通双块生成泛函（CTB-GF）

**取 `F_A∈𝒜(A)`、`F_B∈𝒜(B)`、`d(A,B)=R`，`𝒜(A)` = 支撑在 `A` 内的 Wilson 环及其多项式。定义**

$$\mathcal I_L(s,t)=\log M_L(s,t)-\log M_L(s,0)-\log M_L(0,t),\qquad \partial_s\partial_t\mathcal I_L\big|_{0}=\mathrm{Cov}(F_A,F_B)$$

$$\boxed{\ \big|\mathcal I_L(s,t)\big|\ \le\ C\,|s||t|\,\|F_A\|_{\mathcal A}\|F_B\|_{\mathcal A}\,\phi(R),\qquad \phi(R)\to0\ }\tag{CTB}$$

**Cauchy ⇒**

$$\boxed{\ \big|\mathrm{Cov}(F_A,F_B)\big|\ \le\ C\|F_A\|_{\mathcal A}\|F_B\|_{\mathcal A}\phi(R)\ }\tag{CL}$$

**⇒ 对源求高阶导 ⇒ 控制所有混合连通累积量。**

$$\boxed{\ \text{CTB 是新的、依赖 Gibbs 测度的非微扰输入，【不能】由偶子图估计推出}\ }$$

### 5‴.5 从偶子图走向 CTB：连接桥估计

$$\big|\mathrm{Cov}(F_A,F_B)\big|\ \le\ C\|F_A\|\|F_B\|\sum_{\gamma:A\leftrightarrow B}\mathfrak w_{2.5}(\gamma)\tag{BR}$$

$$\boxed{\ \sup_{A,B:\,d(A,B)\ge R}\ \sum_{\gamma:A\leftrightarrow B}\mathfrak w_{2.5}(\gamma)\ \le\ C e^{-mR}\ }\tag{BR-E}$$

**⇒ 这才是真正要从 Wilson 作用量里找的东西：**

$$\text{不是「每个闭环存在」，而是「连接两个远距离闭环的【桥】有可求和的总权重」}$$

**与早前的 SCE 有联系，但**比 SCE 更具体地保留了两个区域之间的距离**。**

### 5‴.6 逻辑修正：聚簇性 **不** 等于唯一性

**⚠️ 我早前写的「聚簇性 ⟹ 唯一性」是**不严的**：相共存时，每个纯相各自聚簇，而多个极端 Gibbs 态并存。**

**⇒ 唯一性还需边界稳定性：**

$$\boxed{\ \sup_{\eta,\eta'}\big|\mu_{\Lambda,\eta}(F_A)-\mu_{\Lambda,\eta'}(F_A)\big|\ \le\ C\|F_A\|_{\mathcal A}\,\phi\big(d(A,\partial\Lambda)\big)\ }\tag{BI}$$

**⇒ 故原缺件「唯一性」必须拆成三条，不得混称：**

| 子任务 | 目标 | 判 |
|---|---|---|
**CTB / CL** | 相关函数聚簇 | `[U]` |
**BI** | 边界影响消失 | `[U]` |
**极限相容性** | 得到唯一无限体积 Gibbs 态 | `[U]` |

### 5‴.7 对 §5.1 的两处状态修正

**① `max_k Ĉ_q(k) = Ĉ_q(0)` 在 `L=6` 的测量中**观测到**，但：**

- **它不是所有 `L` 上的定理；**
- **它**不**证明 `C_q(z)≥0`。**

**⇒ 我早前的推导「`C_q≥0 ⇒ max_k Ĉ = Ĉ(0)`」作为**蕴含**成立；但前提 `C_q≥0` 只在**窗内**（FKG/GKS）是 `[T]`，在 `β=2.5` 是**观测**。`[T]` 应降为 `[M]`（观测）。**

**② `Ĉ(0)=χ_q` 与 `Σ_zC_q(z)/Ĉ(0)=1` 是**定义性恒等式**，属实现正确性检查，**不是物理结果**。**

**⇒ 真正有信息量的只有：`max_k Ĉ = Ĉ(0)`（有限体积观测），以及 `χ_q` 的 `L`-依赖。**

--- 
## §5⁗ 双块聚簇机器的三路方案（第 258 轮）

> **核验纪律：本节每个数字我都独立算过；算不出的标 `[Δ]`，不写成 `[T]`。**

### 5⁗.1 精确起点：双块自由能差

$$Z_\Lambda(s,t)=\int e^{-S_W}e^{sF_A+tF_B}dU,\qquad \mathcal I_\Lambda(s,t)=\log\frac{Z_\Lambda(s,t)Z_\Lambda(0,0)}{Z_\Lambda(s,0)Z_\Lambda(0,t)}$$

$$\boxed{\ \partial_s\partial_t\mathcal I_\Lambda(0,0)=\mathrm{Cov}_{\mu_{\beta,\Lambda}}(F_A,F_B)\ }\tag{1}$$

$$\partial_s^p\partial_t^q\mathcal I_\Lambda(0,0)=\kappa\big(\underbrace{F_A\cdots F_A}_p,\underbrace{F_B\cdots F_B}_q\big)\qquad(p,q\ge1)$$

**⇒ 不必控制所有局部量的所有矩；只控制**这一个双块自由能差**。**`[T]`

### 5⁗.2 SU(2) 自旋泡沫聚合物表示（`[T]`）

$$e^{\frac\beta2\mathrm{Tr}U_p}=F_0(\beta)\Big[1+\sum_{j>0}d_jc_j(\beta)\chi_j(U_p)\Big],\qquad d_j=2j+1,\ c_j=\frac{I_{2j+1}(\beta)}{I_1(\beta)}$$

$$Z_\Lambda=F_0(\beta)^{|P_\Lambda|}\sum_{\{j_p\}}\Big(\prod_{p:j_p>0}d_{j_p}c_{j_p}\Big)\mathcal A_\Lambda(\{j_p\})\tag{2}$$

**链路 Haar 积分 = 投影到该链路所载表示的**不变量子空间** `P_inv` ⇒ 非 singlet 的标号自动消失。**

$$\boxed{\ Z_\Lambda/Z_{\rm triv}=\Xi_\Lambda=\sum_{\Gamma\ \rm compatible}\prod_{\gamma\in\Gamma}z_\beta(\gamma)\ }\tag{3}$$

**⚠️ 这一步**完全精确，不是强耦合近似**；而 `z_β(γ)` **必须保留 Haar 投影与 intertwiner 收缩**，不是 `∏_p d_jc_j`。**

### 5⁗.3 加权 KP 条件（`[T]` 结构）

$$\boxed{\ \sup_\gamma\frac1{|\gamma|}\sum_{\gamma'\not\sim\gamma}|z_\beta(\gamma')|e^{a|\gamma'|+\kappa\,\mathrm{diam}(\gamma')}\le a_0,\qquad a\gt a_0\gt0\ }\tag{KP-B}$$

**⇒ 标准 KP 保证 `log Ξ` 有绝对收敛的连通簇级数；**连接 `A,B` 的簇必跨过 `R` ⇒ 支付指数权重 ⇒**

$$|\mathcal I_\Lambda(s,t)|\le C_{\rm CTB}|st|\|F_A\|\|F_B\|e^{-mR}\tag{CTB}\ \Longrightarrow\ \text{CL}$$

### 5⁗.4 【本轮硬结果】裸 spin-foam KP 在 `β=2.5` **失败**（`[T]`，我逐位核过）

**最小非零聚合物 = 三维单位立方体的边界（6 面，全 `j=½`）。**

$$u(\beta)=2\frac{I_2(\beta)}{I_1(\beta)},\qquad u(2.5)=1.01439020009404187$$

**12 条边各用 `∫U_{ab}\overline{U_{cd}}\,dU=\tfrac12\delta_{ac}\delta_{bd}` ⇒ `2^{−12}`；8 个顶点指标收缩 ⇒ `2^8`：**

$$\boxed{\ z_\square=\frac{u(2.5)^6}{16}=0.068094225825737\ }\tag{11}$$

*（我算：`u^6=1.089507613211789`，`/16 = 0.068094225825737` — 与你的 `0.06809422583` 逐位相符。）*

**冲突数 `N_inc(C)`：**

| 量 | 值 | 判 |
|---|---|---|
你的值 | `48` | **`[Δ]`** |
**我的直接枚举（`L=4,5`，含/不含环绕）** | **`24`（去重）/ `51`（按共享边计重数）** | **`[T]`（枚举）** |
**⇒ 差一个因子 `2`，需一个定义核对** | — | **`[Δ]`** |

**⇒ 但两个读数都给出同一结论：**

$**⚠️ 第 259 轮勘误：`N_inc` 不是 `48`，也**不是**我上一轮枚举出的 `24`（那是脚本过滤掉边界立方体的 bug）。**

**重写枚举（`L=6` 周期、去重、按方向组分类）逐格复现：**

| 方向组 | `shareLink` | `multLink` | `sharePlaq` |
|---|---|---|---|
`{1,2,3}` | **`18`** | `36` | `6` |
`{1,2,4}` | **`20`** | `32` | `4` |
`{1,3,4}` | **`20`** | `32` | `4` |
`{2,3,4}` | **`20`** | `32` | `4` |
**合计** | **`78`** | `132` | `18` |

**⇒ `N_inc=78` `[T]`（我复现）。而 `sharePlaq=18` 说明：只有**共享 link** 才使 Haar 积分不因子化 ⇒ 正确的是 `78`。**

| `N` | `Nz_□` | 超 `1/e` |
|---|---|---|
`18`（误用 sharePlaq） | `1.2257` | `3.33×` |
`24`（我上轮的 bug 值） | `1.6343` | `4.44×` |
**`78`（正确）** | **`5.31135`** | **`14.44×`** |

$$78z_\square=5.31135\gt\tfrac1e=0.3678794412\quad(\text{大 }14.44\times)$$$
$$48z_\square=3.2685228396\gt\tfrac1e\ (\text{大 }8.88\times)$$

**而 KP 要求 `N z_\square e^{a}\le a`，即 `a e^{−a}\ge Nz_\square`，而 `max_{a\gt0}a e^{−a}=1/e`：**

$$\boxed{\ \text{裸立方体聚合物族【不可能】满足标准均匀 KP}\ }\tag{13}$$

### 5⁗.5 而这条硬结果与今晚的 `1/48 → 2.5` 墙是**同一堵墙**

$$u(2.5)=1.014390200094=\boxed{u_{P3}}\quad\text{（今晚 `mmmm1`/`GGG` 那条 `u\lt1/\lambda_a` 门里的同一个数）}$$

**⇒ 这里它已在 `1` 之上，而 KO 需要它小** ⇒ 所以单 plaquette 系数本身就不小。**

$$\boxed{\ \text{裸展开失败的根因}=\text{那个 `u_{P3}\gt1`}\ }$$

**⇒ 这与 §1 的 `p` 问题是**两个不同的入口**，指向同一件事：`β=2.5` 处没有一个小参数。**

### 5⁗.5b 【新的主判据】加权 Dobrushin 收缩 `MASTER`（第 259 轮）

**不再把验收条件绑在某一种聚合物展开上；改用**带空间权重的块 Dobrushin 收缩**，它一次给三件。**

**粗粒化后的**规范轨道**块变量 `X,Y`，距离 `d(X,Y)`；**

$$c^{\eta,s,t}_{XY}=\sup_{\xi,\xi'}\big\|\gamma_X(\cdot|\xi)-\gamma_X(\cdot|\xi')\big\|_{\rm TV}\quad(\xi,\xi'\ \text{只在块 }Y\ \text{上不同})$$

$$\boxed{\ \text{(MASTER)}\quad q_m:=\sup_{\Lambda,\eta,\ |s|,|t|\le s_0,\ X}\ \sum_Y c^{\eta,s,t}_{XY}\,e^{m\,d(X,Y)}\ \le\ q\lt1\ }$$

**定理 1（`[T]`，我核过）**：`D=(I−C)^{-1}=Σ_nC^n`；由三角不等式 `d(X,Z)≤d(X,Y)+d(Y,Z)` 归纳得 `‖C^n‖_m ≤ q^n`：

$$\boxed{\ |\mathrm{Cov}(F_A,F_B)|\ \le\ \frac{e^{-mR}}{4(1-q)}\|\delta F_A\|_1\|\delta F_B\|_1\ }\tag{CL}$$

**⚠️ 常数修正：严格算出来是 `e^{−m(R−2r)}`（因 `δ_XF_A=0` 当 `d(X,A)>r`）；`e^{2mr}` 吸进常数。**

**定理 2（`[T]`）**：

$$\big|\mu_\Lambda^\eta(F_A)-\mu_\Lambda^{\eta'}(F_A)\big|\ \le\ \frac{e^{-m\,d(A,\partial\Lambda)}}{1-q}\|\delta F_A\|_1\tag{BI}$$

$$\boxed{\ \text{(MASTER)}\Longrightarrow\text{(CL)}+\text{(BI)}\Longrightarrow\text{唯一无限体积 Gibbs 测度}\ }$$

**⚠️ 只要求**实源**；复源还需零点自由/聚合物解析性，不能自动推出。**

### 5⁗.5c `MASTER` 有 `β` 阈值，而这是新的缺口位置

$$\beta=0\ (\text{乘积 Haar})\ \Rightarrow\ c_{XY}=0\,(X\ne Y)\ \Rightarrow\ q_m=0\ \Rightarrow\ \text{成立}$$
$$\beta\to\infty\ \Rightarrow\ c_{XY}\to1\ \Rightarrow\ q_m\to\infty\ \Rightarrow\ \text{失败}$$

$$\boxed{\ \exists\,\beta_{\rm BD}:\ \text{MASTER}\iff\beta\lt\beta_{\rm BD};\qquad\text{靶是 }\beta_{\rm BD}\gt2.5\ }$$

**文献最好值（Rajabu 2026 预印本，**未同行评审**，我未逐字核）：线规+单链路 Dobrushin 到 `β≤0.419`；两次分解后到 `β≤0.511`。**

$$\frac{2.5}{0.511}=4.89\times$$

**⇒ 换更强的机器（唯一性）后，缺口从 `30×` 降到 **`4.89×`**，**但没有换掉那根 `β` 轴**。**

**更早的 Osterwalder–Seiler 已证充分强耦合区的无限体积极限与解析性，但不覆盖 `β=2.5`。**

### 5⁗.5d 缺陷敏感 RG 的接口（`[U]`）

$$S^{\rm eff}_k,\quad F^{\rm eff}_{A,k},\quad F^{\rm eff}_{B,k};\qquad \mathcal I_k(s,t)=\log Z_k(s,t)-\log Z_k(s,0)-\log Z_k(0,t)+\log Z_k(0,0)$$

$$\text{(RG-CTB)}\quad|\mathcal I_k(s,t)|\le C_0|st|e^{-m_kR_k}$$
$$\text{(RG-E)}\quad|E_k(s,t)|:=|\mathcal I_{k-1}(s,t)-\mathcal I_k(s_k,t_k)|\le\epsilon_k|st|e^{-m_kR_k},\qquad \sum_k\epsilon_k\lt\infty$$

**若 `m_kR_k ≥ m_0R − c_k`、`sup_kc_k<∞`，则 `Σ_k|E_k| ≤ C|st|e^{−m_0R}` ⇒ `|I_0| ≤ C'|st|e^{−m_0R}`。**

**⚠️ `RG-CTB` 与 `RG-E` 必须对**实际有效作用量**证明，不由 `Z` 的普通上下界自动获得。而**每一步必须保留两个缺陷**，不许只控无源 `Z`。**
### 5⁗.6 三路方案

| 路线 | 优点 | 最大障碍 |
|---|---|---|
**spin-foam KP** | 完全保留表示结构；直接给复解析性 | **裸聚合物已判失败（§5⁗.4）** |
**块 Dobrushin** | **直接给 CTB + 边界稳定性 + 唯一性** | 需**规范轨道块**条件核的统一收缩 |
**缺陷敏感 RG** | 有机会把目标点带到有效活动小的尺度 | 必须**同步控制两个缺陷**，不能只控 `Z` |

**⇒ 建议采用：主路线 = 缺陷敏感 RG → 粗尺度 spin-foam KP；并行验证 = 规范轨道块 Dobrushin。**

### 5⁗.7 块 Dobrushin 的精确条件（`[T]` 结构，我核过判别式）

**块变量取**规范轨道**（不是裸链变量）。`c_{XY}=\sup_{\eta,\eta'}\|\gamma_X(\cdot|\eta)-\gamma_X(\cdot|\eta')\|_{\rm TV}`，`η,η'` 只在块 `Y` 上不同；`c_{XY}=0` 当 `d(X,Y)\gt r_b`。**

$$\boxed{\ \sup_X\sum_Yc_{XY}\le q\lt1\ }\tag{BD}$$

**Dobrushin 比较不等式：**

$$|\mathrm{Cov}(F_A,F_B)|\le\tfrac14\sum_{X,Y}\delta_XF_A\,[(I-C)^{-1}]_{XY}\,\delta_YF_B$$

**⇒ 因 `[(I−C)^{-1}]=Σ_nC^n`，且连接相距 `R` 需 `n≥⌈R/r_b⌉`：**

$$\sum_{Y}[(I-C)^{-1}]_{XY}\le\frac{q^{\lceil R/r_b\rceil}}{1-q}$$

$$\boxed{\ |\mathrm{Cov}(F_A,F_B)|\le\frac{\|\delta F_A\|_1\|\delta F_B\|_1}{4(1-q)}\,q^{\lceil R/r_b\rceil}\ }\tag{BD-CL}$$

**⇒ 边界稳定性同一条件：**

$$\big|\mu_\Lambda^\eta(F_A)-\mu_\Lambda^{\eta'}(F_A)\big|\le C\|\delta F_A\|_1\,q^{\,d(A,\partial\Lambda)/r_b}$$

**⭐⭐ 这条路线的**关键优势**：`q<1` 直接给**唯一无限体积 Gibbs 态**（不只是聚簇），不需要 `BI` 作为额外假设。**

**⚠️ 但 `q<1` 对 `β=2.5` 是**新非微扰输入**，未证。**

### 5⁗.8 三个工作包（正式验收）

| 包 | 必须交付 | 状态 |
|---|---|---|
**WP1 自旋泡沫** | 完整 SU(2) 活动 + 精确 Haar 投影系数 | 表示基础 `[T]`；**裸 KP 已判失败** `[T]` |
**WP2 缺陷响应** | 两缺陷的统一连通活动界（**须保留距离**） | **`[U]`** |
**WP3 边界唯一性** | `RG-KP` 或 `BD-T` 的体积一致严格证书 | **`[U]`** |

### 5⁗.9 数值任务随之改变

**不再只测 `Cov(q_0,q_z)`，而测两个固定形状、逐渐拉开的 Wilson 环：**

$$\big\{R,\ C_{AB}(R;L),\ C_{AB}(R;L)/C_{AB}(R_0;L)\big\}\qquad(W_{1\times1},W_{1\times2},W_{2\times2})$$

**⚠️ 数值只能筛选假设，不能代替全边界条件、全体积的统一证明。**

---
## §5⁗′ 双缺陷 RG：定理层闭合 + 三阶段分界（第 260 轮）

### 5⁗′.1 定理 A（`[T]`）：空间加权双缺陷收缩

**每层 `k` 的精确有效作用量 `S_k`；四个配分函数 `Z_k(s,t)`、`Z_k(s,0)`、`Z_k(0,t)`、`Z_k(0,0)`：**

$$\mathcal I_k(s,t)=\log\frac{Z_k(s,t)Z_k(0,0)}{Z_k(s,0)Z_k(0,t)}$$

**⚠️ 必须追踪**四个**配分函数。只控 `|Z_k(0,0)|≤C` 不给任何两区间的解耦信息。**

$$\text{(K)}\quad \sup_\gamma\frac1{|\gamma|}\sum_{\gamma'\nsim\gamma}|z_k(\gamma';s,t)|\,e^{a|\gamma'|+\mu\,d_k(\gamma,\gamma')}\ \le\ a-\delta,\qquad 0\lt\delta\lt a$$
$$\text{(D)}\quad \mathcal D_{A,k}(s)\le K_A|s|,\qquad \mathcal D_{B,k}(t)\le K_B|t|$$

**证明骨架**：四个对数之比把**所有不同时接触 `A,B` 的连通簇全部抵消**：

$$\mathcal I_k(s,t)=\sum_{\substack{\mathcal C\ \text{connected}\\ \mathcal C\cap A\ne\varnothing,\ \mathcal C\cap B\ne\varnothing}}w_k(\mathcal C;s,t)$$

**⇒ 每个这样的簇必跨 `A,B`。用树图界 + `e^{μR} ≤ Π_{树边}e^{μd_k}·e^{μC_{\rm geom}Σ|\gamma|}`，沿标记间的树路径反复用 `(K)`，两标记首次插入用 `(D)`：**

$$\boxed{\ |\mathcal I_k(s,t)|\ \le\ C_{\rm KP}K_AK_B|st|\,e^{-\mu R}\ }\tag{A}$$

**⚠️ `m` 必须是**物理单位**的衰减率，不是未换算的格点单位。**

### 5⁗′.2 定理 B（`[T]`）：RG 组合

$$\text{(RG)}\quad \mathcal I_0(s,t)=\mathcal I_K(s_K,t_K)+\sum_{k=0}^{K-1}E_k(s,t)\quad(\text{exact RG；局部归一在双缺陷比值中抵消})$$
$$\text{(E)}\quad|E_k(s,t)|\le C\varepsilon_k|st|e^{-mR}\qquad\text{(S)}\quad\sum_k\varepsilon_k\le E_*\lt\infty\qquad\text{(T)}\quad|\mathcal I_K|\le C_*|st|e^{-mR}$$

$$\boxed{\ |\mathcal I_0(s,t)|\ \le\ (C_*+CE_*)|st|e^{-mR}\ }\tag{RG-CTB}$$

**⇒ `K→∞` 时误差**不累积成无穷大**——这正是 RG 能否闭合的数学条件。若 `E_k=O(1)`，迭代毫无意义。**

### 5⁗′.3 四个待证量（`[U]`）

$$\mathfrak K_k(a)=\sup_\gamma\tfrac1{|\gamma|}\sum_{\gamma'\nsim\gamma}|\zeta_{k,\beta}(\gamma')|e^{a|\gamma'|}$$
$$\mathfrak D_{A,k}=\tfrac1{|s|}\sum_{\gamma\cap A\ne\varnothing}|\partial_s\zeta_{k,\beta}(\gamma;s,0)|e^{a|\gamma|}\qquad(B\ \text{同})$$
$$\mathfrak E_k=\sup_{s,t}|E_k(s,t)|/\big(|st|e^{-mR}\big)$$

$$\boxed{\ \text{(W)}\quad\sup_k\mathfrak K_k\le a-\delta,\quad \sup_k\mathfrak D_{A,k}\le K_A,\quad \sup_k\mathfrak D_{B,k}\le K_B,\quad \sum_k\mathfrak E_k\lt\infty\ }$$

**⚠️ 四项不是四个独立猜想，是**一套证明的四个检查项**：第一项控无源有效作用量；二三控缺陷在 RG 中传播；四控粗粒化的额外连接。**缺二三，即使无源 `RG` 收敛也推不出双块衰减。**

### 5⁗′.4 【勘误】渐近自由系数：`11/(6π²)`，不是 `11/(3π²)`（`[T]`，我两条路核过）

$$\text{1-loop: }g_0^{-2}=\frac{b_0}{16\pi^2}\ln\frac1{a\Lambda},\ \ b_0=\frac{11N}{3}=\frac{22}{3}\ \Longrightarrow\ \beta_W=\frac{4}{g_0^2}=\frac{11}{6\pi^2}\ln\frac1{a\Lambda}$$

$$\frac{11}{6\pi^2}=0.1857555033\qquad\frac{11}{3\pi^2}=0.3715110067\ (\text{原写法，}2\times)$$

**独立复核：`aΛ=exp(−1/(2b_0g²))`、`b_0=11/(48π²)` ⇒ `ln(1/(aΛ))=6π²β/11` ⇒ 同式。✓**

### 5⁗′.5 **而修正后出现结构性好消息：收缩区正是 IR**

$$\beta_{\rm eff}(x)=\frac{11}{6\pi^2}\ln\frac1{x\Lambda}\qquad(\text{随 }x\ \text{下降})$$

| `x` [`Λ^{-1}`] | `β_eff` |
|---|---|
`0.01` | `0.855` |
`0.03` | `0.651` |
**`0.064`** | **`0.511`**（临界） |
`0.1` | `0.428` |
`0.3` | `0.224` |
`1.0` | `0` |

$$\boxed{\ \beta_{\rm eff}\lt\beta_{\rm BD}\iff x\gt0.0639\,\Lambda^{-1}\ }$$

**⇒ 收缩**在长距离成立、在短距离失效**；而长距离正是要聚簇的地方，短距离本来弱耦合、不需要收缩。**

$$\boxed{\ \text{⇒ RG 方向是对的：blocking 把 }\beta_{\rm eff}\text{ 往下推，推进收缩区}\ }$$

**而步数很少（`b=2`）：**

| `a_0` | 到 `x_crit` 的步数 |
|---|---|
`0.01 fm` | **`2.68`** |
`0.02 fm` | **`1.68`** |
`0.05 fm` | **`0.35`** |

$$\boxed{\ \text{(S) 条件的 }\sum_k\varepsilon_k\ \text{只有三个词——不是无穷多尺度}\ }$$

### 5⁗′.6 三阶段分界（不得混称）

| 阶段 | 内容 | 判 |
|---|---|---|
**一 · 数学定理层** | 定理 A（KP ⟹ 双块）+ 定理 B（RG 组合）+ `BD ⟹ CL+BI ⟹ 唯一` | **`[T]`** |
**二 · 目标模型** | `(W)` 的四项在 `β=2.5`（或沿 RG 后的有效 `β`） | **`[U]`** |
**三 · 连续极限** | `β_W(a)→∞` 上物理尺度一致的 `C-CTB`；`S^{\rm ren}_{n,a}→S_n`；OS 条件与非平凡性 | **`[U]`** |

**⚠️ 文献（我**未逐字核**，标 `[M]`）：**
- **Rajabu 2026 预印本（未同行评审）：SU(2) 两次分解后 `β≤0.511`；**
- **另有 2026-10 预印本报 plaquette-block Dobrushin 覆盖 `β_W≤0.345`；**
- **Osterwalder–Seiler 已证充分强耦合区的无限体积极限与解析性，不覆盖 `β=2.5`。**

$$\frac{2.5}{0.511}=4.89\times\qquad\frac{2.5}{0.345}=7.25\times$$

**⇒ 固定 `β=2.5` 的收缩**未被认证**，但**这不等于"在目标点不可能成立"**——它只是说现有 plaquette-block 常数不能直接用。**

### 5⁗′.7 三阶段之间的一道硬边界

$$\text{固定 }\beta=2.5\text{ 的无限体积极限}\ \ne\ \text{四维连续极限}$$

**因为连续极限沿 `β_W(a)→∞` 的轨迹，而 `β=2.5` 是一个**固定的格距点**（而且它靠近 `β_c≈2.7`，属**交叉区**，1-loop 公式在那里定量不可靠）。**

**⇒ 所以 `β=2.5` 的结论**不能直接接到连续极限之后**，必须沿 AF 轨迹重做物理尺度估计。**

--- 
## §5⁗″ 有效活动收缩：尺度 8 块积分 + 三层预算（第 261 轮）

### 5⁗″.1 【勘误 28】上一轮那 `2.68` 步是错的

**上一轮我用 1-loop `β_eff(x)=(11/(6π²))ln(1/(xΛ))` 换算距离并算出 `2.68` 步。核对 `a(2.5)`：**

| 位置 | 1-loop 给 `β_eff` | 真实 `β` |
|---|---|---|
`x=a(2.5)=0.1881\,Λ^{−1}` | **`0.3104`** | **`2.5`** |

$$\frac{2.5}{0.3104}=8.06\times\qquad\Longrightarrow\qquad \text{1-loop 在 }\beta=2.5\text{ 处不可用于换算格距}$$

**⇒ 我上一轮 §5⁗′.5 那条“收缩区就是 IR”的**方向结论仍对**（blocking 降 `β_eff`），但**步数是错的**。**

**⇒ 正确的校准是**弦张力两端匹配**（下面 §5⁗″.2）。**

### 5⁗″.2 可信的 landing 尺度：`b=8.08 ≈ 2³`

$$\beta_W=\frac4{g^2}\ \Longrightarrow\ g^2(2.5)=1.6;\qquad a\sqrt\sigma\,(2.5)=0.1881(28)\ \text{（格点数据）}$$

$$\beta'=0.40:\quad u(\beta')=\frac{I_2}{I_1}=0.099339929651,\qquad a'\sqrt\sigma\simeq\sqrt{-\ln u}=1.519607738$$

$$b=\frac{1.519607738}{0.1881}=8.0787226902\approx2^3\ \checkmark\qquad(8\times0.1881=1.5048)$$

$$\boxed{\ \beta=2.5\ \xrightarrow{Q_2}\ \beta_1\ \xrightarrow{Q_2}\ \beta_2\ \xrightarrow{Q_2}\ \beta_3\simeq0.40\ }$$

**⚠️ 弦张力匹配只给**粗粒化尺度的候选值**，**不证明**粗化后作用量真有 `β'=0.40`，更不证明其余活动小。它是**校准工具，不是 landing 定理**。**

### 5⁗″.3 walk-gas 阈值：`β < 2/3`（`[T]`，我核过）

$$c=\frac\beta4,\qquad 6c\lt1\iff\beta\lt\frac23=0.666667$$

| `β` | `c` | `6c` | 判 |
|---|---|---|---|
`0.40` | `0.100` | `0.600` | OK |
`0.419` | `0.105` | `0.629` | OK |
`0.511` | `0.128` | `0.767` | OK |
**`0.711`** | `0.178` | **`1.067`** | **FAIL** |
`2.5` | `0.625` | `3.750` | FAIL |

**⚠️ `β=0.711` 已过不了朴素的 `6c<1` ⇒ 那个结果用的是**精细化的界**，不是此判据。记账须分清。**

$$\boxed{\ \text{原来的高温游走级数【不能】承担 }\beta=2.5\text{ 的收缩证明}\ }$$

**⚠️ 这**不**说精确积分不存在，也**不**说真实有效活动必发散——只说这一种展开的绝对收敛证明不行。**

### 5⁗″.4 各方法的上限（到目标的倍数）

| 方法 | 上限 | `2.5/上限` |
|---|---|---|
基础强耦合 Dobrushin | `0.345` | `7.25×` |
单层 decimation | `0.419` | `5.97×` |
fiber–base 精确展开 | `0.511` | `4.89×` |
strand fiber–base（后续预印本，未独立审查） | `0.711` | `3.52×` |
**pinning obstruction（类硬顶）** | **`1.01092`** | **`2.47×`** |

$$\boxed{\ \text{pinning 顶 }1.01092\ \text{封的是【一整类】pinning-inherited 判据，而靶 }2.5\ \text{在其外 }2.47\times\ }$$

**⚠️ 它**不**排除真正不同的粗粒化规格 ⇒ 必须换受控对象。**

### 5⁗″.5 精确 RG 的入口与求导恒等式（`[T]`）

$$S_{\rm eff}(V)=S_{\rm coarse}(V)-\log\int e^{-S_{\rm fine}(V,\xi)}d\xi\tag{18}$$

$$F(V):=-\log\int e^{-S_{\rm fine}}d\xi\quad\Longrightarrow\quad \partial_iF=\mathbb E_V[\partial_iS_{\rm fine}]$$
$$\boxed{\ \partial_i\partial_jF=\mathbb E_V[\partial_i\partial_jS_{\rm fine}]-\mathrm{Cov}_V\big(\partial_iS_{\rm fine},\partial_jS_{\rm fine}\big)\ }\tag{19}$$

**⇒ 若两粗链变量不在同一细 plaquette 项中直接耦合，`∂i∂jS_fine=0`：**

$$\boxed{\ \partial_i\partial_jF=-\mathrm{Cov}_V\big(\partial_iS_{\rm fine},\partial_jS_{\rm fine}\big)\ }\tag{20}$$

$$F_{ijk}=\mathbb E[S_{ijk}]-\mathrm{Cov}(S_{ij},S_k)-\mathrm{Cov}(S_{ik},S_j)-\mathrm{Cov}(S_{jk},S_i)+\kappa_3(S_i,S_j,S_k)\tag{21}$$

$$\boxed{\ \text{远距离有效相互作用}=\text{细变量条件测度下两局部作用量插入的【连通关联】}\ }$$

**⚠️ **循环性**（你已指出）：若直接假设 `(23)` 的指数衰减，就已经用了细尺度的聚簇性。**
**⇒ 但这是**两尺度论证**，不是循环：**

$$\underbrace{\beta=2.5\ \text{细}}_{\text{弱耦合}\Rightarrow\text{聚簇较易}}\xrightarrow{3\ \text{步 }Q_8}\underbrace{\beta'\approx0.40\ \text{粗}}_{\text{强耦合}\Rightarrow\text{收缩较易}}$$

**⇒ 两半各在**自己擅长的耦合端**；难的是**搬运**，而搬运方向 `2.5→0.40` **向下**，是 RG 的自然流向。**

### 5⁗″.6 三层预算与 landing 目标（`[U]`）

$$S_{\rm eff}^{(8)}=S_W(0.40;V)+R^{(8)},\qquad R^{(8)}=R_{\rm local}+R_{\rm bridge}+R_{\rm tail}$$

$$\|R\|_{r,m}=\sup_\ell\sum_{X\ni\ell}|X|^2e^{m\,\mathrm{diam}X}H_X(r),\qquad H_X(r)=\sup_{U\in T_r^X}|R_X(U)|$$

$$\boxed{\ \|R^{(8)}\|_{r,1/30}\ \le\ 6.33\times10^{-3}\qquad(r\ge1,\ m=1/30)\ }\tag{26}$$

**三项各自的预算：**

| 项 | 内容 | 需要的界 |
|---|---|---|
`R_local` | `Σ_{|C|≤6}r_C ReTr V(C)` | 局部重正化，吸收主量级 |
| `R_bridge` | `Σ_{diam≥2}R_X` | `Σ_{X∋e}|X|^2e^{m\,diam}\|R_X\|_{T_1}\le\eta_{\rm bridge}` |
| `R_tail` | `Σ_{diam≥N}R_X` | `\|R_{\rm tail}\|_{1,1/30}\le Ce^{−\mu N}` |

**⚠️ `(26)` 的三个上界**必须从 `(19)–(22)` 的条件累积量推出，**不能从 plaquette 平均值推出来**；也**不能**由弦张力匹配或 MC 余量代替。**

**⚠️ 无限范围修正：真正的块积分生成无限范围 `R`，须用**准局部**版本（`β'≤0.405` 下加权混合差分容差 ≤ `0.01`）。不能把有限范围范数证书不加说明地用于无限范围作用量。**

**⚠️ **Gibbs 性**：块映射的推前测度**一般不自动**是 Gibbs 的（条件分布可能对远处粗变量不连续依赖）⇒ 须独立证明，否则 `S_eff` 根本不存在。**

### 5⁗″.7 数值该测什么（并说明它不能替代什么）

$$D(R)=\sup_{A,B:\,d(A,B)\ge R}\big|\Delta_A\Delta_B S_{\rm eff}^{(8)}\big|,\qquad \text{对照 }D(R)\le Ce^{-mR}$$

**⚠️ 有限样本的实轴估计**不能**证明 `sup_{V∈T_1^X}|R_X(V)|` 这类**复解析**范数。**
**⇒ MC 只用于**筛选候选 RG 映射与定位支撑范围**；最终证书须用区间算术或严格解析估计。**

--- 
## §5⁗‴ 三层分解的定量审计：`Q` 与 `R_0` 方向相反（第 262 轮）

### 5⁗‴.1 定理 1（`[T]`，我数值核过）

**ENV：`‖R_X‖_r ≤ Aϑ^{|X|−1}e^{−μτ(X)}`，`0<ϑ<1`、`μ>m`。`c:=μ−m`、`θ:=Δ²ϑe^{−c/2}`、`F(θ):=(1+θ)/(1−θ)³`。**

$$\boxed{\ \sup_e\sum_{\substack{X\ni e\\ \mathrm{diam}X\gt R}}|X|^2e^{m\,\mathrm{diam}X}\|R_X\|_r\ \le\ AF(\theta)\,e^{-cR/2}\ }$$

**证明要点（我核）：`diam ≤ τ` ⇒ `e^{m·diam−μτ} ≤ e^{−cτ}`；`τ ≥ max(R,n−1) ≥ (R+n−1)/2`；`N_n(e) ≤ Δ^{2(n−1)}`；最后 `Σ_{n≥1}n²θ^{n−1} = (1+θ)/(1−θ)³`。**

| `θ` | 数值和 | 公式 |
|---|---|---|
`0.3` | `3.790087` | `3.790087` ✓ |
`0.5` | `12.000000` | `12.000000` ✓ |
`0.8` | `225.000000` | `225.000000` ✓ |

### 5⁗‴.2 `θ<1` 比 `Q<1` 严 `~40×`（`[T]`）

$$\theta=\Delta^2\vartheta e^{-c/2},\ \ \vartheta=Q^{1/s_0}\ \Longrightarrow\ \theta\lt1\iff Q\lt e^{c/2}/\Delta^2$$

**四维超立方粗图 `Δ = 2d = 8` ⇒ `Δ² = 64`：**

| `μ` | `c` | `θ<1` 要求的 `Q <` |
|---|---|---|
`0.10` | `0.0667` | **`0.01615`** |
`0.50` | `0.4667` | `0.01973` |
`1.00` | `0.9667` | `0.02534` |

$$\boxed{\ \text{实际小参数}Q\lesssim0.025\text{，不是 }Q\lt1\ }$$

### 5⁗‴.3 尾预算对 `R_0` 的下界

$$R_0\ge\frac2c\ln\frac{AF(\theta)}{\delta_{\rm tail}},\qquad F(0.5)=12,\ A=1,\ \delta_{\rm tail}=2.11\times10^{-3}$$

| `μ` | `R_0` 下界（粗格点单位） |
|---|---|
`0.10` | **`259.4`** |
`0.50` | **`37.1`** |
`1.00` | **`17.9`** |
`1.15` | `15.5` |

### 5⁗‴.4 【阻塞】两项条件方向相反

$$\text{`Q` 小}\iff\mu\ \text{小}\ (\text{因 }e^{2\mu d_0})\qquad\text{而}\ R_0\ \text{小}\iff c=\mu-m\ \text{大}\iff\mu\ \text{大}$$

**⇒ 联合解只由 `eC₀D₀κ` 决定（`d₀=1`）：**

| `eC₀D₀κ` | `μ_max`（`Q<1`） | 需要的 `R₀` |
|---|---|---|
`0.5` | `0.347` | **`55.2`** |
`0.1` | `1.151` | **`15.5`** |
`0.01` | `2.303` | **`7.6`** |
`1.2\times10^{-3}` | `3.363` | **`5.2`** |

**而若要求 `R₀ ≤ 3`（真正的"局部"）：**

$$c\ge\frac{2\ln(AF/6.33\times10^{-3})}{3}\approx5.04\ \Longrightarrow\ \mu\gtrsim5.07$$
$$Q=eC_0D_0\kappa\,e^{2\mu d_0}\approx2.5\times10^4\,eC_0D_0\kappa\lt0.025$$

$$\boxed{\ eC_0D_0\kappa\ \lesssim\ 10^{-6}\qquad\Longrightarrow\qquad C_0\kappa\ \lesssim\ 4\times10^{-7}\ }$$

### 5⁗‴.5 判词

$$\boxed{\ \text{三层分解【逻辑】正确，但【定量】要求细尺度近乎独立}\ }$$

$$\underbrace{m=\tfrac1{30}}_{\text{允许极慢衰减}\Rightarrow\text{尾要 }R_0\text{ 大}}\quad\Longleftrightarrow\quad\underbrace{Q\lt0.025}_{\text{要求 }\mu\text{ 小}}\quad\Big|\quad c=\mu-m\ \text{被夹在中间}$$

**⇒ 而 `β=2.5` 是交叉区，`C₀κ ≲ 4\times10^{-7}` 在那里不可能被供给。**

### 5⁗‴.6 三条修法与收益（`[T]` 的估算）

| 修法 | 效果 |
|---|---|
① 换更弱的 landing 范数（`m` 更小 / 不同加权） | 直接放松 `R₀`，**最有效** |
② 不按 `η_*/3` 均分、按最优重分配 | 因 `ln` 慢，收益约 `1.3×`：`R₀` 从 `18` 降到 `15`，`Q` 只放松 `2×` ⇒ **不够** |
③ **换受控对象**（不强求 `‖R‖ ≤ 6.33e−3`） | **根本出路** |

**⇒ ① 值得先查：landing 范数的 `m=1/30` 是否可换；② 明确不够；③ 与"必须换对象"的结论一致。**

--- 
## §5⁗⁗ Uniform Fiber-KP：三条定理闭合 + landing 点必须更深（第 263 轮）

### 5⁗⁗.1 精确块映射与 Haar 分解（`[T]`）

$$V_E=U_{\ell_1}\cdots U_{\ell_b}\qquad(\text{各粗路径的细链路互不重叠，四维立方块直线边可实现})$$

**取路径末端为 pivot：`U_{\ell_b}=(U_{\ell_1}\cdots U_{\ell_{b-1}})^{−1}V_E`。由 Haar 左右不变性：**

$$\boxed{\ \prod_{\ell\in E_f}dH(U_\ell)=\prod_{E\in E_c}dH(V_E)\,dH_{\rm fib}(\xi)\ }\tag{2}$$

$$\frac{d\nu^{(b)}_\Lambda}{dH_c}(V)=\frac1{Z_\Lambda}\int e^{-S_W(U(V,\xi))}dH_{\rm fib}(\xi)\tag{3}$$

**⇒ 无近似、无截断、不丢任何细自由度。**

### 5⁗⁗.2 共同前提 UFBP

$$\boxed{\ \sum_{\gamma'\nsim\gamma}\|z_{\gamma'}\|_r\,e^{(a+\lambda)n_{\gamma'}}\ \le\ (a-\delta)n_\gamma\ }\qquad a,\lambda,\delta\gt0,\ \|z_\gamma\|_r:=\sup_{V\in T_r}|z_\gamma(V)|$$

**⚠️ 必须对**所有**粗边界配置一致成立（不只对典型配置）——否则坏配置可使被约束的细系统隐藏相共存，导致推前测度不连续。**

### 5⁗⁗.3 三条推论（`[T]`，我核过结构）

**WARN: [errata 29, r264] the three corollaries below are CONDITIONAL on V (fiber level only). They do NOT compose into UFBP => CTB => uniqueness. See section 5-quad-quad-quad-quad.**

**（一）连接桥：`N(Γ) ≥ (R−D_A−D_B)/c_*` ⇒ `Σ_{Γ:A↔B}|ω(Γ)| ≤ C_1K_AK_B|st|e^{−m_1R}`，`m_1=λ/c_*`。**
**核心不是给活动加 `e^{−λdiam}`，而是**每个连通簇的总尺寸 ≥ 连接两块所需的距离**。**

**（二）解析尾：`|X|≤N`、`diam(X)≤c_*N`（`c_*=1+r_inc`）；`δ:=λ−mc_*>0`；**

$$\sup_e\sum_{\substack{X\ni e\\ \mathrm{diam}X\gt R}}|X|^2e^{m\,\mathrm{diam}X}\|\Phi_X\|_r\ \le\ C_\delta C_{\rm root}\,e^{-\delta R/(2c_*)}\tag{TAIL}$$

$$\nu=\frac{\lambda-mc_*}{2c_*}\qquad C_{\rm tail}=C_\delta C_{\rm root}$$

**（三）Gibbs 性：三步——精确 Haar 分解给条件核 `(9)`；UFBP 给一致 quasilocality `(10)`；有界能量给 non-nullness `(11)`。**

$$\boxed{\ e^{-12\beta}\ \le\ \frac{d\gamma^b_E}{dH}\ \le\ e^{12\beta}\qquad(\beta=2.5:\ e^{-30}\le\cdot\le e^{30})\ }\tag{12}$$

**推导（我核）：改一个粗链 = 改一条细链；该细链在 `2(d−1)=6` 个 plaquette 中；每项 `(β/2)Re tr`，`|Re tr|≤2` ⇒ 每项范围 `2β` ⇒ 总变化 `≤12β` ✓**

**⇒ 三者统一来自同一个条件细变量的复解析加权聚合物展开。**

### 5⁗⁗.4 【阻塞三】`ν>0` 与小参数之间的闭式张力

**KP 要求 `Σ‖z‖e^{(a+λ)n} ≤ (a−δ)n`。以 `‖z_γ‖~θ^n`、含给定聚合物的形状数 `~μ^n`：**

$$\mu\theta e^{\lambda}\lt1\ \Longrightarrow\ \lambda\le\ln\frac1{\mu\theta}$$

$$\boxed{\ \nu\ \le\ \frac{|\ln(\mu\theta)|}{2c_*}\ \Longrightarrow\ R_1\ \gtrsim\ \frac{2c_*\,\ln(C_{\rm tail}/\delta_{\rm tail})}{|\ln(\mu\theta)|}\ }$$

**取 `c_*=2`、`ln(1/2.11e−3)=6.16`：**

| `μθ` | `R₁` 下界 |
|---|---|
`0.9` | **`278`** |
`0.5` | **`36.4`** |
`0.1` | **`10.8`** |
`0.01` | **`5.4`** |
`0.001` | `3.6` |

### 5⁗⁗.5 【关键】`β'=0.40` 处空间权重**不存在**

**以 `μθ = 10u`、`u=I₂/I₁`：**

| `β'` | `u` | `μθ` | 判 |
|---|---|---|---|
`0.05` | `0.012499` | `0.125` | `R₁=11.9` |
`0.1` | `0.024990` | `0.250` | `R₁=18.0` |
`0.2` | `0.049917` | `0.499` | `R₁=36.3` |
`0.3` | `0.074720` | `0.747` | `R₁=89.7` |
**`0.40`** | `0.099340` | **`0.993`** | **`ν<0` ⇒ 不可能** |

$$\boxed{\ \text{在 } \beta'=0.40\ \text{处 } \nu\lt0\text{——正是那里 KP 才勉强成立（}\mu\theta\approx0.99\text{）}\ }$$

### 5⁗⁗.6 因此 landing 点必须**更深**

$$R_1\le5\ \Longleftarrow\ |\ln(\mu\theta)|\ge4.92\ \Longleftarrow\ \mu\theta\le7.3\times10^{-3}\ \Longleftarrow\ u\le7.0\times10^{-4}$$

$$b=\frac{\sqrt{-\ln u}}{0.1881}\ \Big|_{u=7\times10^{-4}}=\frac{2.69}{0.1881}=14.3\approx2^{3.84}$$

$$\boxed{\ \text{landing 应为 }\beta'\approx0.02,\ b\approx2^4,\ \text{【四步】而非三步；}R_1\approx5\ \text{粗}={\sim}80\ \text{细}\ }$$

**⚠️ 而 `m` 的杠杆很小：我算了 `m: 1/30 → 1/100` 只把 `R₁` 从 `36` 降到 `35`。真正的杠杆是 `μθ`（有效活动本身）。**

**⇒ 所以"`β'` 更深"不是选择，是**必然**。**

### 5⁗⁗.7 不能把什么当成证明

**⚠️ 以下**都不能**替代 MASTER-UFBP：**

- **`L=6` 的 plaquette 数据、`χ_q/Λ=0.0114`；**
- **链变量的中心对称性结果（§5″）；**
- **弦张力匹配得到的 `b≈8`（现在应为 `b≈2⁴`）。**

**⇒ MASTER-UFBP 必须来自**固定任意粗链变量后的细变量条件测度**，且上确界必须遍历**所有**粗配置（只算典型配置会漏掉 Gibbs 性最关键的 bad configurations）。**

$$\boxed{\ \text{条件细变量的谱/混合控制}\ \Longrightarrow\ \text{Fiber-KP}\ \Longrightarrow\ \text{三重闭合}\ }$$

--- 
## §5⁗⁗⁗ 【勘误 29】纤维收缩 ≠ 粗变量聚簇：五模块审计（第 264 轮）

### 5⁗⁗⁗.1 关键恒等式（`[T]`）

$$\mathrm{Cov}(F_A,F_B)=\underbrace{\mathbb E[\mathrm{Cov}(F_A,F_B\mid\mathcal V)]}_{\text{UFKP 控这一项}}+\underbrace{\mathrm{Cov}(\mathbb E[F_A\mid\mathcal V],\mathbb E[F_B\mid\mathcal V])}_{\text{粗变量自己的协方差}}$$

**⇒ 全协方差律（恒等式）。UFBP 只控**第一项**。**

$$\boxed{\ \text{细纤维解耦}\ \ne\ \text{粗变量解耦}\ }$$

**⚠️ 我 §5⁗⁗.3 写「UFBP ⟹ CTB ⟹ 唯一」是**越权**：那里的 `I(s,t)` 是**条件在 `V` 上**的纤维双缺陷自由能，界对 `V` 一致；那只给第一项。**

**反例形状（一般概率模型）：给定全局隐变量 `S` 后各局部纤维完全独立，但 `S` 有两个取值且同时控制两远隔区域的局部均值 ⇒ 无条件协方差仍是常数。**

### 5⁗⁗⁗.2 修正后的逻辑

$$\text{UFKP}\Longrightarrow\text{有效纤维积分解析且准局部}$$
$$\text{BASE}:\quad \sup_X\sum_Yc^{\rm base}_{XY}e^{m\,d(X,Y)}\le q\lt1\ \Longrightarrow\ \text{粗变量聚簇 + 边界稳定}$$
$$\boxed{\ \text{UFKP}+\text{BASE}\ \Longrightarrow\ \text{原模型双块聚簇}\ }$$

### 5⁗⁗⁗.3 好场/坏场版 UFBP（结构 `[T]`）

$$\text{(UFKP-G)}\quad \sup_{V\in\mathcal G_k}\sup_\gamma\frac1{|S(\gamma)|}\sum_{\gamma'\nsim\gamma}\|z_{\gamma'}(V)\|_r\,e^{(a+\lambda)|S(\gamma')|}\ \le\ a-\delta$$
$$\text{(BAD)}\quad \mathbb P(\exists\ \text{坏场连通簇连接 }A,B)\ \le\ C e^{-m_{\rm bad}d(A,B)}$$

$$\boxed{\ |\mathrm{Cov}(F_A,F_B)|\ \le\ C\big(e^{-m_{\rm good}R}+e^{-m_{\rm bad}R}\big)\ }$$

**⚠️ `BAD` 仍需独立证明：单点坏场概率小**不足以**控制连接两区的坏场簇。**

### 5⁗⁗⁗.4 【本轮定量阻塞】`β=2.5` 处没有小曲率好场

**我从自己的 `L=6`、`β=2.5` 实测取数：**

$$\langle q\rangle=0.34727537,\quad \mathrm{Var}(q)=1.4200\times10^{-2},\quad sd=0.1192,\quad \tfrac{sd}{\langle q\rangle}=0.343$$

$$q\in[0,2]\ \Longrightarrow\ \text{典型曲率是最大值的 }17.4\%$$

**若好场要求 `q_p<ε`（高斯近似，仅量级）：**

| `ε` | `P(q<ε)` |
|---|---|
| `0.02` | `0.003` |
| `0.05` | `0.006` |
| `0.10` | `0.019` |
| `0.20` | `0.108` |
| **`0.35`** | **`0.509`** |
| `0.70` | `0.999` |

$$\boxed{\ \text{要有可观测度的好场须 }\varepsilon\gtrsim0.35\text{，而那时「小曲率」已不是小量}\ }$$

**⇒ `BAD` 前因子极大 ⇒ 指数界被吃掉。**

**⇒ 与今晚其他全部诊断同源：**

| 量 | `β=2.5` 的值 | 判 |
|---|---|---|
| `u(2.5)=I₂/I₁` | `0.507195` | `u_P3=2u=1.014\gt1` |
| `⟨q⟩` | `0.3473` | 曲率 `O(1)`，不小 |
| `μθ`（`β'=0.40`） | `0.993` | `ν\lt0` |
| `χ_q/Λ` | `0.011` | **这一项有余量** ✓ |

$$\boxed{\ \beta=2.5\ \text{在四个入口上都没有小参数（只有 F2 例外）}\ }$$

### 5⁗⁗⁗.5 五模块审计（正式记账）

| 模块 | 内容 | 判 |
|---|---|---|
| **A** | 链变量有限体积生成泛函（中心对称 → 偶子图 → GF + RTC） | **`[T]`**（`β`-无关，**不给聚簇**） |
| **B** | 精确块映射 `dν_b(V)∝Z_fib(V)dH_c(V)`；单粗链 non-nullness `e^{−30}` | **`[T]`** |
| **C** | 条件纤维收缩 `UFKP-G + BAD` | **`[U]`**（`β=2.5` 无证书） |
| **D** | 粗变量聚簇 + 唯一性 `BASE` | **`[U]`**（**不能从 C 推出**） |
| **E** | 连续极限 | **`[U]`** |

### 5⁗⁗⁗.6 文献窗口（我**未逐字核**，标 `[M]`）

| 结果 | `β_W` 上限 | `2.5/上限` |
|---|---|---|
| 两级 Fiber–Base 唯一性 | `0.511` | `4.89×` |
| 两级聚类/传递矩阵隙 | `0.51` | `4.90×` |
| 两级复耦合簇 | `0.465` | `5.38×` |
| strand 展开（2026-11-21 后续预印本） | `0.711` | `3.52×` |

**⚠️ 教训：不能只增加聚合物项再期待 KP 自动成立，必须改变对活动的**组织与估计方式**。**

**⚠️ 另：2026-10 的好场/坏场方案也未完成坏场在每一条件纤维上的统一控制、跨尺度组合与连续极限。**

---

## §5⁗⁗⁗′ 原 UFBP 两条收缩路线的严格排除 + Block-UFBP 的新障碍（第 265 轮）

### 5⁗⁗⁗′.1 两条障碍是**同一个数**（`[T]`，我核过）

$$\text{walk-gas: }6c=6\cdot\frac{\beta_W}{4}=\frac{3\beta_W}{2}\qquad\text{Dobrushin: }\sum_{w}C_{vw}\ \ge\ 6\cdot\frac{\beta_W}{4}=\frac{3\beta_W}{2}$$

$$\boxed{\ \text{两路是同一表达式 }\frac{3\beta_W}{2}\text{；}\beta_W=2.5\ \text{时}=3.75\gt1\ }$$

**⇒ 不是两个独立障碍，是一个障碍的两个面。**

### 5⁗⁗⁗′.2 线性响应下界（`[T]`）

**`ν_h(x)∝e^{h·x}dσ` on `S³`；`Cov_σ(x)=\frac14 I_4`（因 `Σx_a²=1` 且对称）。改一个邻居 `x_w→x_w'`：**

$$h\to h+\beta_W R_{vw}(x_w'-x_w)\ \Longrightarrow\ \big|\delta\mathbb E[x]\big|=\frac{\beta_W}{4}\delta$$

$$W_1^{\rm ch}\ \ge\ \big|u\cdot(\mathbb E_\nu x-\mathbb E_{\nu'}x)\big|=\frac{\beta_W}{4}\delta\ \Longrightarrow\ C_{vw}\ \ge\ \frac{\beta_W}{4}$$

**⇒ `Σ_{w}C_{vw} ≥ 6β_W/4 = 3.75` at `β_W=2.5`。而这是**下界** ⇒ 真实行和更大 ⇒ 障碍比陈述的更强。**

**⚠️ 且障碍在 `h=0` 的**线性响应**处，**不依赖**数值或 von Mises–Fisher 积分精度 ⇒ 提高精度无用。**

**⚠️ 且 `Σ_wC_{vw}` 的值域是 `[3β_W/2,\ 6]`（小 `β` 到铁磁有序）⇒ 单点 Dobrushin 永不收缩除非 `β_W\lt2/3`。**

### 5⁗⁗⁗′.3 pinning 障碍（`[T]`）

$$G_3(0)=\int_{[-\pi,\pi]^3}\frac{d^3p}{(2\pi)^3}\frac{1}{2\sum_i(1-\cos p_i)}=0.2527310098589096\ (\text{Watson})$$

$$\boxed{\ \beta_{\rm pin}=4G_3(0)=1.0109240394\ }$$

**⇒ `2.5/1.01092 = 2.473×` ⇒ 任何**既推唯一性、又继承任意链路固定**的原始 Wilson 收缩条件都不可能覆盖 `2.5`。**

**⚠️ 该障碍**不**直接适用于线规与特定 decimation 之后的条件规格。**

### 5⁗⁗⁗′.4 平坦扇区与 Poincaré（`[T]` 恒等式）

**沿生成树递推 `O_u^{T}R_{uv}O_v=I`，令 `y_v=O_v^{T}x_v` ⇒ 树边 `R̃=I`，非树边保留全部 holonomy。**

$$\beta\sum_{\langle u,v\rangle}(1-y_u\cdot y_v)=\frac\beta2\sum_{\langle u,v\rangle}|y_u-y_v|^2\qquad(\text{因 }|y|=1)$$

$$\sum_{v\in B}|y_v-\bar y_B|^2\ \le\ \frac1{\lambda_1(B)}\sum_{\langle u,v\rangle\subset B}|y_u-y_v|^2,\qquad \lambda_1(B)\asymp\ell^{-2}$$

### 5⁗⁗⁗′.5 【新障碍，我算的】Block-UFBP 的 Poincaré 控制在 `β=2.5` 是**空**的

**记 `q := ⟨1−y_u·y_v⟩`。块 `ℓ^4` 位点、`~2ℓ^4` 条边：**

$$\langle E\rangle\sim\beta\cdot2\ell^4 q\ \Longrightarrow\ \Big\langle\sum|\nabla y|^2\Big\rangle\sim\frac2\beta\langle E\rangle=4q\ell^4$$

$$\text{平凡界：}v\le4\ell^4\qquad\text{Poincaré 界：}v\le\ell^2\cdot4q\ell^4=4q\ell^6$$

$$\boxed{\ \text{Poincaré 有用}\iff 4q\ell^6\lt4\ell^4\iff q\,\ell^2\lt1\ }$$

| `q` | `ℓ=2` | `ℓ=3` | `ℓ=4` | `ℓ=6` |
|---|---|---|---|---|
| `0.05` | OK | OK | OK | **空** |
| `0.10` | OK | OK | **空** | **空** |
| `0.20` | OK | **空** | **空** | **空** |
| **`0.347`** | **空** | **空** | **空** | **空** |

**而 `β=2.5` 实测 `⟨q⟩=0.3473` ⇒ `qℓ²=1.39` at `ℓ=2`、`3.13` at `ℓ=3`：**

$$\boxed{\ \ell\ge2\ \text{即已空}（除非 }q\text{ 很小）\ \Longrightarrow\ \text{块化无法改善快速模控制}\ }$$

**⇒ 根因：`λ_1(B)≍ℓ^{-2}` 使 `1/λ_1` 随块**变大**，而能量只按 `ℓ^4` 涨 ⇒ 净损失 `ℓ^2`。**

$$\boxed{\ \text{要有用须 }\beta_W\gg q\ell^2\text{——即弱耦合；而 }2.5\ \text{不是}\ }$$

### 5⁗⁗⁗′.6 判词

| 路线 | 判 | 原因 |
|---|---|---|
原彩色闭合游走 UFBP | **排除** | `6c=3.75\gt1` |
单点 Dobrushin UFBP | **排除** | `Σ_wC_{vw}\ge3.75\gt1`（与上同一表达式） |
任意 pinning 继承型 | **排除** | `β_{\rm pin}=1.01092` |
**Block-UFBP** | **新障碍** | **`qℓ^2\gt1` ⇒ 快速模 Poincaré 界空** |
多尺度递推 `RG-2` | 无递推常数 | 未从 Wilson 作用量导出 |

**⇒ 真正值得保留的方向**：慢模式（`M_B`)必须留在下一层，但**快速模的控制不能靠 Poincaré**——须用块内的**能量-熵**或**平均场型**估计。**

---

## §5⁗⁗⁗″ 分隔集构造（`[T]`）+ `ℓ^{10}` 把 `ℓ` 选择排除（第 266 轮）

### 5⁗⁗⁗″.1 三个 `[T]` 结果（我逐步核过）

**① 精确因子化**：`S_ℓ` 分隔集 ⇒ `F∖S_ℓ = ⊔_B C_B` 无跨块边 ⇒

$$Q_\zeta(dx_{F\setminus S_\ell}\mid x_S)=\bigotimes_B Q_{\zeta,B}(dx_{C_B}\mid x_{\partial C_B})$$

**⇒ 不同快块之间连通累积量严格为零。不依赖 `β` 大小。**

**② 快块 Poincaré**：`λ₁(S³)=3`（`−Δ` 本征值 `l(l+2)`，`l=1` 给 `3`）；乘积 Haar 仍 `3`。`osc(V_C) ≤ 2βM_C ≤ 12βn_C`；Holley–Stroock：

$$\boxed{\ \lambda_1(Q_C)\ \ge\ 3e^{-12\beta n_C}\ \overset{\beta=2.5}{=}\ 3e^{-30n_C}\ }$$

**⇒ 对一切 `ζ`、一切慢边界条件、且不依赖全局体积。**

**③ 一致复解析半径**：`‖X_C‖²=n_C` ⇒ `‖X_C−EX_C‖ ≤ 2√n_C`；

$$r_\ell=\frac{\ln(3/2)}{2\sqrt{n_\ell}}\ \Longrightarrow\ |M_C(J)-1|\le\tfrac12\ \Longrightarrow\ M_C\ne0,\ \ |\log M_C|\le\ln2$$

$$\boxed{\ \|\log Z_C\|_{r_C}\ \le\ \beta M_C+\ln2\ \le\ 15n_C+\ln2\ }$$

### 5⁗⁗⁗″.2 【本轮关键】`‖S_eff‖ ≲ βℓ^{10}e^{mℓ}` 把 `ℓ` 选择排除

$$\|\log Z_C\|=O(\beta\ell^4),\quad |\partial C|=O(\ell^3),\quad \mathrm{diam}(\partial C)=O(\ell)$$

$$\boxed{\ \|S^{\rm eff}_{\rm slow}\|_{m,r}\ \lesssim\ \beta\,\ell^{10}e^{m\ell}\ }\tag{19}$$

**而 landing 目标是 `6.33×10⁻³`，`β=2.5`、`m=1/30`：**

| `ℓ` | `βℓ^{10}e^{mℓ}` | 超目标 |
|---|---|---|
| **`1`（不块化）** | **`2.585`** | **`408×`** |
| `2` | `2.74×10³` | `4.3×10⁵×` |
| `3` | `1.63×10⁵` | `2.6×10⁷×` |
| `4` | `3.00×10⁶` | `4.7×10⁸×` |
| `6` | `1.85×10⁸` | `2.9×10¹⁰×` |

$$\boxed{\ \ell=1\ \text{已超 }408×\text{，而 }\ell^{10}\ \text{使增大 }\ell\ \text{【更差】}\ \Longrightarrow\ \text{无任何 }\ell\ \text{可选}\ }$$

**⇒ 这与 `6c=3.75\gt1` 一致（`ℓ=1` 就是未块化的情形）。**

### 5⁗⁗⁗″.3 元凶是范数里的 `|X|²` 因子（我算的）

$$\|\Phi_X\|\sim\beta\ell^4,\qquad |X|\sim\ell^3,\qquad e^{m\,\mathrm{diam}}\sim e^{m\ell}$$

| 权重 | 增长 | `ℓ=2` 的值 |
|---|---|---|
| **`|X|²`（现行）** | **`βℓ^{10}`** | **`2560`** |
| `|X|` | `βℓ^7` | `320` |
| 不带 | `βℓ^4` | `40` |

$$\boxed{\ |X|^2\ \text{贡献了 }10^3\text{ 量级爆炸；换 }|X|\ \text{只降到 }10^2\ }$$

**⚠️ 但换范数须与 landing 定理的证明配套——不能单方面改。**

### 5⁗⁗⁗″.4 Gaussian/spin-wave 一致性检验（`[T]`）——它揭示**两层要求不同**

$$\|\Gamma_{\rm fast}\|\ \le\ \frac1{\beta\lambda_*}\ \le\ \frac{\ell^2}{c_0\beta}$$

| `ℓ` | `‖Γ_fast‖` |
|---|---|
| `2` | `1.60` |
| `3` | `3.60` |
| `4` | `6.40` |

**⇒ 且 spin-wave 的 `κ_n=0`（`n≥3`）⇒ `W=½JΓJ`，`‖Γ‖=O(1)` ⇒ **GF 层成立**。**

$$\boxed{\ \textbf{GF 层只需 }\Gamma\text{ 有界}\Rightarrow O(1)\text{ 够；而 KP/landing 层因 }|X|^2\ \text{权重爆炸}\ }$$

**⇒ **两层是不同强度的要求**：GF（生成泛函）易，KP（landing/唯一性）难。**

**⇒ 这也解释了为何 §5″ 的链变量 GF 那般干净，而 landing 一路卡住。**

### 5⁗⁗⁗″.5 状态表（本轮更新）

| 项 | 判 |
|---|---|
| 精确分隔集构造 / 快变量因子化 | **`[T]`** |
| 快块 Poincaré `3e^{−30n_C}` | **`[T]`** |
| 固定块尺度复解析半径 `r_ℓ` | **`[T]`** |
| 全阶快条件累积量界 | **`[T]`** |
| Gaussian spin-wave Block-UFBP | **`[T]`** |
| **`‖S_eff‖ ≲ βℓ^{10}e^{mℓ}` ⇒ 无 `ℓ` 可选** | **`[T]`（本轮排除）** |
| 非线性好块 Hessian 界 `(23)` | `[U]` |
| 联合测度坏块指数尾 `(BAD)` | `[U]` |
| 慢变量块影响矩阵收缩 `(BLOCK-SPEC)` | `[U]` |
| `β=2.5` 完整 Block-UFBP | **`[U]`** |

### 5⁗⁗⁗″.6 下一步该证的**唯一**一条（我同意你的选择，并加一个范数条件）

$$\boxed{\ \inf_{\|\pi\|_\infty\le\rho,\ \epsilon_B\le\epsilon}\lambda_{\min}\big(D^2H_{\rm fast}\big)\ \ge\ c\,\frac{\beta}{\ell^2}\ }$$

**⚠️ 但必须先解决 `(19)` 的 `ℓ^{10}`：若范数不改，即便 `(23)` 成立，`‖S_eff‖` 仍进不了 landing 球。**

$$\boxed{\ \text{完整的一步 = }(23)\ +\ \text{一个不随 }\ell\ \text{爆炸的范数（或 }|X|\ \text{替代 }|X|^2\text{）}\ }$$

### 5⁗⁗⁗″.7 作业记账

**`L=16` 的 `SU(2)`、`β=2.5` 生产跑已超单作业预算，按 `CLAUDE.md` 纪律 kill（一次一个、不炸机）。已有 `L=6/8/10/12` 四份数据足够支撑 B 路线的 `(SB)` 结论；`L=16` 只用于加固 `p`，不是必需。**

---

## §5⁗⁗⁗‴ 非线性好块 Hessian 升为 `[T]`；全基配置慢收缩降为 `[F]`（第 267 轮）

### 5⁗⁗⁗‴.1 `ℓ^{10}` 的**准确措辞**（采纳你的修正）

$$B(\ell)=2.5\,\ell^{10}e^{\ell/30}\ \Longrightarrow\ \frac{d}{d\ell}\log B=\frac{10}{\ell}+\frac1{30}\gt0\ \Longrightarrow\ B\ \text{递增}$$

$$\inf_{\ell\ge1}B(\ell)=B(1)=2.584738\ \Longrightarrow\ \frac{B(1)}{6.33\times10^{-3}}=408.3$$

$$\boxed{\ \text{正确表述：现有 }\ell^{10}\ \textbf{上界}\ \text{无法为任何块尺度提供 landing 证书}\ }$$

**⚠️ 它**不能**证明真实的 `‖S_eff‖` 一定超阈——得到的是**上界**，不是**下界**。**

### 5⁗⁗⁗‴.2 `|X|²` 不能删（`[T]`，理由核过）

**P4 要控固定链路 `ℓ` 下、另两条不同链路的三阶混合差分。一个支撑 `X` 中有序链路对最多 `(|X|−1)(|X|−2) ≤ |X|²`；Cauchy 条带估计给 `ε₃(R_X;ℓ,ℓ',ℓ'') ≤ κ₃(r)H_X(r)`：**

$$(|X|-1)(|X|-2)H_X(r)\ \le\ |X|^2H_X(r)$$

**⇒ `|X|²` 是**从统一复解析振幅推 P4** 时用的计数上界，**不是**实际三阶差分必达的规模。**

**⇒ 所以把 `|X|²` 换成 `|X|` 再引用原 landing 定理是**不合法的**。**

### 5⁗⁗⁗‴.3 唯一合法的降权路线：直接估计 `P1–P4`

$$\mathfrak P_4(R)=\sup_\ell\sum_{\ell'\ne\ell}\ \sum_{\substack{\ell''\ne\ell\\ \ell''\ne\ell'}}w_m(\ell,\ell',\ell'')\!\!\sum_{X\supset\{\ell,\ell',\ell''\}}\!\!\epsilon_3(R_X;\ell,\ell',\ell'')$$

$$\boxed{\ \mathfrak P_4(R)\ \le\ \kappa_3(r)\|R\|_{r,m}\qquad\text{而反向不等式不成立}\ }$$

**⇒ 可以尝试直接控 `𝔓₁,…,𝔓₄`；若满足底层 robust fibre–base 定理的四项阈值，就可走对应 landing 证明，而不必强求通用解析振幅范数低于同一门槛。**

**⚠️ **必须重新核对底层定理的四项阈值——不能直接套用 `0.00633`。****

### 5⁗⁗⁗‴.4 【升级为 `[T]`】非线性 `O(4)` 好块点态 Hessian

**在 `ℳ=(S³)^{V(F)}` 上用乘积黎曼度量与 Bakry–Émery。`s_e=x_u^TR_ex_v`、`V=−βΣs_e`。**

**好块条件取**协变自旋对齐**（不只是 base holonomy 小）：**

$$\mathcal G_\epsilon=\{x:\ s_e\ge1-\epsilon,\ \forall e\in E(F_B)\}$$

$$\boxed{\ \mathrm{Ric}+\mathrm{Hess}V\ \succeq\ \big(2-\beta\Delta\epsilon\big)I+\beta L_R,\qquad \Delta\le6\ }$$

**其中 `L_R` 是协变图 Laplacian：`⟨ξ,L_Rξ⟩=Σ_{e=(u,v)}‖ξ_u−R_eξ_v‖²`。**

**证明（我逐步核过）：**
1. **`Ric_{S³}=(n−2)g=2I`（`n=4`）✓**
2. **测地二阶变分：`Hess V(ξ,ξ)=βΣ_e[s_e(|ξ_u|²+|ξ_v|²)−2ξ_u^TR_eξ_v]` ✓**
3. **恒等式 `s(|ξ_u|²+|ξ_v|²)−2ξ_u^TRξ_v = ‖ξ_u−Rξ_v‖²−(1−s)(|ξ_u|²+|ξ_v|²)`——**我数值抽检 2000 组，IDENTITY OK****
4. **`0≤1−s_e≤ε`；`Σ_e(|ξ_u|²+|ξ_v|²) ≤ 6Σ_v|ξ_v|²`（每顶点至多 6 边）⇒ 结论。∎**

**目标点：`β=2.5`、`Δ=6`、`ε=1/15` ⇒ `βΔε = 2.5·6/15 = 1`：**

$$\boxed{\ \mathrm{Ric}+\mathrm{Hess}V\ \succeq\ I+\tfrac52 L_R\qquad\text{on }\mathcal G_{1/15}\ }$$

**⚠️ 数学边界（不可省）：**
- **这是 good set 内的**点态** Hessian 下界，不是球面乘积上的全局凸性；**
- **要得条件 LSI 还须该区域的测地凸性或另行建立带边界的 Poincaré/LSI；**
- **不能把上两式直接宣称为完整条件块 LSI；**
- **只有 base holonomy 小、而纤维自旋任意时**不**成立：即便 `R_e=I`，一条边两端取反向对齐也可使该边 Hessian 贡献为负。**

### 5⁗⁗⁗‴.5 【`[F]`】`H_C=−I` 强制坏边 ⇒ 全基配置慢收缩**不成立**

$$\theta=\arccos(14/15)=0.367208,\qquad 4\theta=1.468832\lt\pi=3.141593$$

**若四边全 good，沿回路传输给 `d(x,H_Cx)≤4θ`；而 `H_C=−I` 给 `d(x,−x)=π`：**

$$\boxed{\ H_C=-I\ \Longrightarrow\ \text{至少一条边必定为 bad}\ }$$

**⇒ 对 `H_C=−I` 的基配置，`Q_ζ(B bad)=1` ⇒ **任何 `sup_ζ Q_ζ(B bad) ≤ p\lt1` 的统一条件概率估计都不成立**。**

$$\boxed{\ \text{Uniform-all-base slow Dobrushin}\ =\ [F]\ }$$

**⚠️ `[F]` 的限定：它否定的是**会继承任意 pinning 并推出唯一性**的强版本，**不是**否定完整未固定 Wilson 测度下的唯一性。**

### 5⁗⁗⁗‴.6 【我的定量补充】坏基配置在**联合测度**下有多稀有

**`q=2 ⟺ ½Re tr U_□=−1 ⟺ U_□=−I`，正是 `H_C=−I` 那种。用我自己 `L=6`、`β=2.5` 的 MC：**

$$\langle q\rangle=0.347275,\quad sd=0.119164,\quad \frac{2-\langle q\rangle}{sd}=13.87\sigma$$

| 阈值 | `P(q>阈值)`（高斯外推） |
|---|---|
`1.933` | `e^{-88.5}=3.5\times10^{-39}` |
`1.867` | `e^{-81.3}=4.8\times10^{-36}` |
`1.500` | `e^{-46.8}=4.8\times10^{-21}` |

**⇒ 坏 plaquette 密度 `~e^{-88}` 级 ⇒ `τ ≳ 88` 每块；而 4D 格点动物数 `μ_an~7` ⇒ `lnμ_an=1.95`：**

$$\boxed{\ \frac{\tau}{\ln\mu_{\rm an}}\ \approx\ 49\ \gg\ 1\ \Longrightarrow\ \text{联合坏块簇尾求和余量极大}\ }$$

**⚠️ 这是**高斯外推**（`e^{-88}` 远在采样范围外）⇒ 标 `[M]`，**不**标 `[T]`。**

### 5⁗⁗⁗‴.7 改名为 **Annealed Good–Bad Block RG**

$$\text{(A1)}\quad \mu_\beta(E^{\rm bad}_{A,B})\le C_1e^{-m_1R}$$
$$\text{(A2)}\quad \big|\mathrm{Cov}(F_A,F_B\mid\text{good fast sector})\big|\le C_2e^{-m_2R}$$
$$\text{(A3)}\quad \big|\mathrm{Cov}(\bar F_A,\bar F_B)\big|\le C_3e^{-m_3R},\qquad \bar F_A=\mathbb E[F_A\mid\mathcal S]$$

$$\mathrm{Cov}(F_A,F_B)=\mathbb E[\mathrm{Cov}(F_A,F_B\mid\mathcal S)]+\mathrm{Cov}(\bar F_A,\bar F_B)$$
$$\Longrightarrow\quad|\mathrm{Cov}(F_A,F_B)|\le(C_1'+C_2+C_3)e^{-mR},\qquad m=\min(m_1,m_2,m_3)$$

**⚠️ `(A1)` 与 `(A3)` 是新的**依赖 Wilson 相互作用**的估计，仍未证明。**

**⚠️ 关键结构：**慢变量保留为相互作用的场**，**不要求每一个固定基配置下都收缩**；改证**联合测度**下的慢变量聚簇。**

### 5⁗⁗⁗‴.8 `ℓ^{10}` 在新策略里的位置：两条并行预算

**① 快速变量条件响应：由 Hessian 得 `κ_fast = 2 − βΔε ≥ 1`（`ε=1/15`）。若另建 good 域几何凸性，可得条件 LSI/混合差分界。**

**② 直接 landing tolerance：不用 `|X|²H_X` 粗替三阶差分，而从有效作用量导数直接算——**

$$F_i=\mathbb E[S_i],\qquad F_{ij}=\mathbb E[S_{ij}]-\mathrm{Cov}(S_i,S_j)$$
$$F_{ijk}=\mathbb E[S_{ijk}]-\mathrm{Cov}(S_{ij},S_k)-\mathrm{Cov}(S_{ik},S_j)-\mathrm{Cov}(S_{jk},S_i)+\kappa_3(S_i,S_j,S_k)$$

**⇒ 远距离三阶混合项由局部插入之间的**连通关联**控制；支撑不重叠时某些接触项严格为零；若其余协方差与三点累积量有可求和衰减，`𝔓₄` 可能远小于 `|X|²H_X`。**

$$\boxed{\ \text{这条路不改 landing 范数、不改 landing 定理，直接瞄准 }P1–P4\ }$$

### 5⁗⁗⁗‴.9 审计（哪些真的 U→T）

| 项 | 判决 | 说明 |
|---|---|---|
**非线性 `O(4)` 好块 Hessian** | **`[T]`** | `Ric+Hess V ⪰ I` 精确条件（我核过恒等式与常数） |
good-domain 条件 LSI | `[U]` | 还需凸性 + 边界条件 |
**联合测度坏块指数尾** | **`[M]`/`[U]`** | 高斯外推给 `τ/lnμ_an≈49`，但外推不可当证明 |
**全基配置慢收缩** | **`[F]`** | 与 pinning 后 3D `O(4)` 长程有序冲突 |
联合测度慢变量聚簇 `(A3)` | `[U]` | annealed 条件 |
完整 Annealed Block-UFBP | `[U]` | 等 `A1–A3` |
`ℓ^{10}` landing 证书 | **上界无法认证** | 所有 `ℓ≥1` 超门槛 `≥408×` |
直接 `P1–P4` 容差路线 | **可行替代** | 不改 landing 定理；须重核四项阈值 |

---

## §5⁗⁗⁗⁗ EFDS：快变量协方差**精确归零** + `SLOW-H` **不可能**（第 268 轮）

### 5⁗⁗⁗⁗.1 【升级】EFDS：快变量双块协方差恒为零（`[T]`）

**给定 `(x_S,ζ)`，快变量按分量 `C_i` 张量分解。`F_A` 只依赖 `𝒞(A)`、`F_B` 只依赖 `𝒞(B)`。若**

$$\mathcal C(A)\cap\mathcal C(B)=\varnothing\qquad(\text{两区相距 }R\gg\ell\ \text{时成立})$$

**⇒ 条件独立 ⇒ `Cov(F_A,F_B\mid x_S,ζ)=0`。全协方差公式第一项**恒为零**：**

$$\boxed{\ \mathrm{Cov}_\mu(F_A,F_B)=\mathrm{Cov}_{\nu_{\rm slow}}(\bar F_A,\bar F_B)\ }\tag{EFDS}$$

**⇒ 这比上一轮的 `(A2)`（good fast sector 协方差**小**）强：现在是**精确零**。**

**⇒ 快变量之间的传播尾**不需要估计**，也不需要对它们做近似独立假设。**

$$\boxed{\ \text{(SLOW-CL)}\Longrightarrow\text{(FULL-CL)}\qquad(\text{因 }\|\bar F_A\|_\infty\le\|F_A\|_\infty\text{，条件期望是 }L^\infty\text{ 收缩})\ }$$

### 5⁗⁗⁗⁗.2 精确慢作用量（`[T]`）

$$d\nu_{\rm slow}(x_S,\zeta)=\frac1Z e^{H_{SS}(x_S,\zeta)}\prod_iZ_i(x_{\partial C_i},\zeta)\,d\sigma(x_S)d\zeta$$
$$S_{\rm slow}=-H_{SS}-\sum_i\log Z_i,\qquad \mathrm{diam}(\partial C_i)\le c\ell$$

**⇒ 有限体积中的**精确边缘测度**，不是近似 RG 作用量；固定 `ℓ` 下有限范围、不依赖 `L`。但**耦合强度未必小**。**

### 5⁗⁗⁗⁗.3 【本轮决定性】`SLOW-H` 不可能：协方差项符号相反且量级压倒

$$\partial_a\partial_b(-\log Z_i)=-\underbrace{\mathbb E[\partial_{ab}H_i]}_{O(\beta)}-\underbrace{\mathrm{Cov}(\partial_aH_i,\partial_bH_i)}_{\succeq0\ \text{带负号}}$$

**⇒ 协方差项**半负定** ⇒ 破坏凸性。**

**量级（`β=2.5`，边界链 `M_∂=O(ℓ^3)`）：**

| `ℓ` | `O(β²ℓ^3)`（Fisher 项） | `O(β)`（Hessian 项） | 判 |
|---|---|---|---|
`1` | `6.25` | `2.50` | **凹** |
`2` | `50.00` | `2.50` | **凹** |
`3` | `168.75` | `2.50` | **凹** |
`4` | `400.00` | `2.50` | **凹** |

$$\boxed{\ \mathcal S_{\rm slow}\ \text{在慢变量上【强凹】}\ \Longrightarrow\ \text{(SLOW-H)}\ \lambda_{\min}(\mathcal H_{\rm slow})\ge\kappa_s\gt0\ \text{不可能成立}\ }$$

**⇒ 物理上一致：`−log Z` 是自由能，Legendre 结构下对耦合本就**凹**。**

**⇒ **所以曲率/LSI 路线到 `SECC` 是关闭的**——不是“未证”，是“不可能”。**

**⇒ 好块 Hessian 的正曲率（§5⁗⁗⁗‴.4 的 `Ric+Hess V ⪰ I`）**不能**传递到慢变量有效作用量。**

### 5⁗⁗⁗⁗.4 坏块动物计数与指数尾（`[T]` 条件定理）

$$N_n\le\Delta_s^{2(n-1)};\qquad \mu_{\rm slow}(K\subset\mathcal B)\le C_0e^{-\tau|K|}\ \text{(BAD-AN)}$$

$$\tau\gt2\ln\Delta_s\ \Longrightarrow\ P(\text{坏簇大小}\ge n)\le C_1e^{-(\tau-2\ln\Delta_s)n}$$
$$\Longrightarrow\quad P(\exists\ \text{坏簇连接 }A,B)\le C_1e^{-m_bR},\qquad m_b=\frac{\tau-2\ln\Delta_s}{c_s}$$

**⚠️ **`[T]` 只是 `BAD-AN ⇒ BAD-TAIL` 这一步**；`BAD-AN` 本身在 `β=2.5` 的慢变量测度下仍未证。**

### 5⁗⁗⁗⁗.5 【勘误 30】我上一轮那个 `e^{−88}` 高估了——Cantelli 直接给弱得多的界

**只用均值与方差（不假设高斯），Cantelli：**

$$P(q-\mathbb E q\ge t)\le\frac{\mathrm{Var}(q)}{\mathrm{Var}(q)+t^2}$$

| 阈值 | Cantelli 界（plug-in） | 我上轮的高斯数 |
|---|---|---|
`q>1.500` | **`0.010574`** | `4.8\times10^{-21}` |
`q>1.933` | **`0.005616`** | `3.5\times10^{-39}` |

$$\boxed{\ \text{高斯外推比 Cantelli 严格界弱了 }10^{17}\text{–}10^{34}\ \text{倍}\ }$$

**⇒ 所以上一轮那个 `τ/lnμ_an ≈ 49` 的“余量极大”**不成立**——Cantelli 只给 `P ≈ 10^{-2}`，即 `τ ≲ 4.6`。而 `2lnΔ_s = 2ln8 = 4.16`：**

$$\frac{\tau}{2\ln\Delta_s}\ \lesssim\ \frac{4.6}{4.16}\ =\ 1.11$$

**⇒ 只余 `1.11×`——**几乎没有余量**。**

**⚠️ 而且三重缺口仍在（你已指出）：**
1. **输入是 MC 估计，无自相关与置信区间；**
2. **`P(q_P>·)` 是**单点**尾概率，而 `BAD-AN` 要的是**整块集合**的联合概率；**
3. **`q_P` 是 plaquette holonomy 的函数，而 `s_e≥14/15` 是**条件纤维图上相邻自旋对齐**——须通过具体纤维映射证明联系，不能视为同一随机变量。**

**⇒ 所以 `e^{−88}` 那一行**撤回**，改标 `[M]` 且注明是**弱于 Cantelli 的外推**。**

### 5⁗⁗⁗⁗.6 `‖Γ‖<∞` 不蕴含聚簇（`[T]` 反例）

$$C_L(x,y)=L^{-4}\ \Longrightarrow\ \sum_yC_L(x,y)=1,\ \ \chi_q=1,\ \ \|\Gamma\|=1$$

**⇒ 二阶量全部 `O(1)`，而 `d(x,y)≍L` 时 `C_L` **不衰减**。**

$$\boxed{\ \text{F2 是必要的二阶诊断，【不是】完整聚簇的替代}\ }$$

### 5⁗⁗⁗⁗.7 SECC：唯一的目标命题

$$\boxed{\ \text{(SECC)}\quad \big|\mathrm{Cov}_{\nu_{\rm slow}}(\bar F_A,\bar F_B)\big|\le C_s\|\bar F_A\|_\infty\|\bar F_B\|_\infty e^{-m_sR}\ }$$

**⇒ 加上统一边界影响估计 `|ν_slow^η(\bar F_A)−ν_slow^{η'}(\bar F_A)| ≤ Ce^{−m d(A,∂Λ)}`，即可由紧群紧性得唯一无限体积极限。**

**⚠️ `[M]` 文献：4D `SU(2)` 强耦合公开结果仍在 `β_W≤0.711`（2026-11-21 预印本，未独立审查），中间 crossover 明确列为未解决。**

### 5⁗⁗⁗⁗.8 下一步该攻哪里（我加上本轮的符号障碍）

**你想攻 `D²S_slow`。我把本轮的结论写清：**

$$D^2(-\log Z_C)=-\mathbb E[D^2H_C]-\mathrm{Cov}(DH_C,DH_C)$$

$$\boxed{\ \text{协方差项带负号且为 }O(\beta^2\ell^3)\ \Longrightarrow\ \text{不能用曲率/LSI 得 SECC}\ }$$

**⇒ 所以 `SECC` 必须走**关联衰减**路线（`Cov(DH_C, DH_C)` 本身的距离衰减），**不是**曲率路线。**

**⇒ 而这正好回到 §1 的 `p` 与 §5⁗ 的 `BR`/`BR-E`：**

$$\big|\mathrm{Cov}(F_A,F_B)\big|\ \le\ C\|F_A\|\|F_B\|\sum_{\gamma:A\leftrightarrow B}\mathfrak w_{2.5}(\gamma)$$

**⇒ 慢变量的桥权重。这就是 `SECC` 的真内容。**

### 5⁗⁗⁗⁗.9 审计（本轮）

| 项 | 判决 |
|---|---|
精确快变量因子化 / **双块快协方差恒为零（EFDS）** | **`[T]`** |
精确慢变量有效 Gibbs 作用量 | **`[T]`** |
固定快块条件 Poincaré `3e^{−30n_C}` | **`[T]`** |
`SLOW-CL ⟹ FULL-CL` | **`[T]`** |
坏块动物计数 ⇒ 指数尾（条件定理） | **`[T]`** |
**`SLOW-H` 曲率路线** | **`[F]`（本轮：符号相反且量级压倒）** |
`β=2.5` 的 `BAD-AN` 联合概率 | `[U]`（Cantelli 只给 `P≈10^{−2}`，余量 `1.11×`） |
**`SECC`（慢变量指数聚簇）** | **`[U]` 唯一目标** |
慢变量边界稳定性 | `[U]` |
完整 Wilson 双块聚簇 | `[U]` |

---

## §6 最终构造（`[Δ-C]` / `[Δ-M]`，独立）

$$\boxed{\ \lim_{\delta\downarrow0}\sup_L\sup_{p_s(f)\le\delta}\big|\chi_L(f)-1\big|=0\ }\tag{M}$$

**还需**：有限体积特征泛函的正定性；共同核测试空间上的定义；极限特征泛函的点态收敛；零点处连续性。

**⚠️ `Λ_L^{(c)} ~ ln L` **不自动**给出上述一致连续性 ⇒ **不能把 TS 与 Minlos 混为一个命题。**

---

## §7 五台机器（严格分开）

| 机器 | 内容 | 核心输入 | 判 |
|---|---|---|---|
**A · 链变量存在性** | `(Z2)` → 偶子图 → GF → F2-link → RTC(`a=0`) → TS | **局部中心对称性（`β`-无关）** | **`[T]`**（§5″） |
**B · 复合观测量** | `q` 族的 F2、`p,σ`、全方向源响应 | Wilson 关联 | `[Δ]`（`L=6/8` 有诊断） |
**C · 双块聚簇** | **CTB-GF、桥估计 `BR`/`BR-E`、边界稳定性 `BI`** | **Wilson 作用量的桥权重（`β`-有关）** | **`[U]` 新主任务** |
**D · 连续极限** | 共同拓扑空间、尺度一致紧性、重构 | `Δ-C` | `[U]` |
**E · 质量隙** | 谱隙 `∃δ>0` | — | **独立项目** |

**⚠️ 分工是结构性要求，不是懒惰：**

$$\text{A 不需要知道 }\beta\qquad\text{C 必须看见 }\beta=2.5$$

**⇒ §5‴.1 的混合测度反例证明：A 全成立而 C 全缺。**

**⚠️ 逻辑边界：**
1. **A 不推出 C/E**；`L3` 类命题**必须接受自由无质量场** ⇒ 对隙零判别力。
2. **A 不推出完整连续极限**（仍需正性、相容性、测试空间、紧性）。
3. **B 只给二阶尺度检验，不能替代 A 的全方向全阶控制。**

---

## §8 三路线总图（第 256 轮重写）

### A 路线 · 链变量张成族（**已闭合**，§5″）

$$\underbrace{g_v=-I}_{[T]}\Rightarrow\underbrace{(Z2)+\text{偶子图}}_{[T]}\Rightarrow\underbrace{|W_L(J)|\le384\|J\|_2^2}_{[T]}\Rightarrow\underbrace{\text{RTC}(a{=}0)}_{[T]}\Rightarrow\underbrace{\text{TS}}_{[T]}\Rightarrow\underbrace{\text{全阶}}_{[T]}$$

| 项 | 内容 | 判 |
|---|---|---|
**GF 对张成族** | `M_L(J)≠0`，`\|J\|<r_*=0.088388`；`\|W_L\|\le3\|J\|^2` | **`[T]`** |
**F2** | `Γ_link=(1/4)I_{16L^4}`，`‖Γ‖=1/4` | **`[T]` 精确** |
**RTC** | `‖DH_L(J)[u]‖≤b‖u‖`，`b=1536/r_*=17377.856`，`a=0` | **`[T]`** |
**`Λ`** | **`Λ^{link}=1`（常数，无 `ln L`）** | **`[T]`** |
**组装** | `C_2=C_0+br_0=768.25`；`\|W_L\|\le384.125\|J\|^2` | **`[T]`，`L`-一致** |

**⚠️ 但见 §5″.8：五个前提全与 `β` 无关 ⇒ 判别力为零 ⇒ **不给**聚簇性、唯一性、隙。**

### B 路线 · plaquette 观测族（`q_x`）

$$\text{需独立建立}\quad\|\Gamma^q_L\|_{2\to2}\le C\Lambda_L^{(q)},\qquad \Lambda_L^{(q)}=1+\sup_x\sum_y g_c(x-y)$$

**⚠️ `Γ^{link} \ne Γ^{q}`——两条路线不可混记。**

**`L=6`、`β=2.5`、`M=20000` 实测（`bmc_L6_b2.500.npz`，本轮跑）：**

| 量 | 值 |
|---|---|
`⟨1−½Re tr U_□⟩` | `0.34727537` ⇒ `⟨½Re tr U_□⟩=0.6527246`（与仓内 `L=8/10/12` 的 `0.652402/0.652099/0.651984` 同趋势） |
`Var(q)` | `1.42001639e−02` |
`χ_q = L⁴Var(q̄)` | **`3.19255284e−02`** |
`max_k Ĉ = Ĉ(0)` | **✓**（峰在零动量，与 `C_q≥0` 一致） |
`Σ_z Cz/Ĉ(0)` | `1.000000000` ✓ |
**`χ_q/Λ`** | **`0.0114`（余 `87×`）** |

**⇒ F2 对 `q` 族在 `L=6`、`β=2.5` 有 `87×` 余量；但 `L=6` 只有 `3` 个壳，判不了 `p,σ`。**

### C 路线 · 极限测度

**固定 `a`**：`SU(2)^{E_∞}` 紧（Tychonoff）⇒ Prokhorov ⇒ **子列**极限存在 `[T]`（**不需要 Minlos**）。
**⚠️ 唯一性需聚簇性，而 A 路线给不出（`β`-无关）。**
**连续极限 `a→0`**：`Δ-C`——共同拓扑空间、插值重构、尺度一致紧性。**`[U]`**

$$\boxed{\ \text{唯一真的非微扰缺口}=\Delta\text{-C}\ \ (\text{连续极限})\ }$$

### 三个必须记住的非蕴含

$$\text{F2}\not\Longrightarrow\text{TF2}\ (\text{反例 A})\qquad \text{TF2}\not\Longrightarrow\text{TS}\ (\text{反例 B})\qquad \Gamma^{\rm link}\ne\Gamma^{q}$$
## §8′ 任务优先级（第 257 轮）

| # | 任务 | 机器 | 判 | 状态 |
|---|---|---|---|---|
**1** | 统一 `g_c` 与 `Λ_L^{(c)}`（§2） | B | **`[T]`** | 已完成 |
**2** | 链变量 GF / F2 / RTC(`a=0`)（§5″） | A | **`[T]`** | **已完成** |
**3** | **CTB 双块混合累积量 + 桥权重 `BR`（§5‴.4–.5）** | **C** | **`[U]`** | **新主任务** |
**4** | **边界稳定性 `BI`（§5‴.6）** | C | **`[U]`** | **与 3 并列，不得混称** |
**5** | 极限相容性 ⇒ 唯一无限体积 Gibbs 态 | C | `[U]` | 依赖 3+4 |
**6** | `Δ-C`：连续极限 | D | `[U]` | 独立 |
**7** | `q` 族 `(SB)`/`(SPEC)`/`κ₃` | B | `[Δ]` | 诊断 |
**8** | 机器 E | — | — | 独立项目 |

**⇒ 下一轮研究重点：从「继续加大 `L` 的单点拟合」转到「CTB 的双块混合累积量与桥权重」。**

**⚠️ 诚实：第 3–6 项是仅剩的真缺口；第 7 项是诊断，不再承担 GF 的证明任务。**
## §9 一句话规格

$$\boxed{\ \text{核心问题：在 }\beta=2.5\text{ 下，}\ \sup_L\sup_{\|J\|_2\lt r_0}\frac{\|D^2\log M_L(J)\|_{2\to2}}{\Lambda_L^{(c)}}\lt\infty\ \text{？}\ }$$

**若成立 ⇒ TS ⇒ GF ⇒ 全阶矩与累积量控制 ⇒ 再加 `[Δ-M]` 即可研究极限测度。**

**SCE 的定位：**重要的结构性证明工具，但**不再限制整台机器只能通过一种空间点态结构来建造**。**

---

## §10 不声称什么

1. **不声称质量隙**；且靶 `ϖ≤C` 在自由场与 2D YM 上都成立 ⇒ **对隙零判别力**。
2. **不声称 SCE 或 TS 是 `(GF)` 的必要条件**（都是充分）。
3. **不声称已构造极限测度**（缺 `[Δ-M]`）。
4. **不声称 `O_B` 的主系数是定理**（`2π²λ` 与 `c_reg` 均依赖核与截断）。
5. **不声称矩 = 绝对矩。**
6. **不使用 SZZ 窗作推导起点、不使用 `β` 轴搬运、不使用分段数 `N`。**
7. **不声称 A 路线有判别力**：§5″ 的五个前提全与 `β` 无关（含 `β=0`）⇒ 该界**不能区分**平凡无隙的乘积 Haar 与 `β=2.5` 的 YM。它给的是**存在性/矩界**，不是隙。
8. **不声称把 `Γ^{link}` 的结果搬到 `Γ^{q}`**。两条路线的核与常数不同，禁止混记。
9. **不声称无限体积极限（固定 `a`）是唯一的**——只给子列极限。

---

*第 252 轮。TS 版取代第 251 轮 SCE 主猜想版；SCE 降为证明路线 B 的强形式。
本轮两处实质修正：① `Λ_L` 的核与系数统一（§2，含 `7.40×` 的正确位置）；② 三点界由三项×2 改为六项（§3）。*
