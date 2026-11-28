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

## §5⁗⁗⁗⁗′ 【勘误 31】`SLOW-H` 的 `[F]` 撤回 → `[U]`；`K` 矩阵估计 `[T]`（第 269 轮）

### 5⁗⁗⁗⁗′.1 【勘误 31】我上一轮的 `[F]` 是**错的**

**上一轮我写：`Fisher` 项 `O(β²ℓ³)`、`Hessian` 项 `O(β)` ⇒ `S_slow` 强凹 ⇒ `SLOW-H` 不可能 `[F]`。**

**⚠️ 这条推理**不成立**：`O(β²ℓ³)` 是**协方差的上界量级**，不是**统一正下界**；而球面变量上两项可以**精确抵消**。**

**你的精确反例（我数值核了）：**

$$H(y,x)=\beta y^{\mathsf T}Rx,\qquad Z(y)=\int_{S^3}e^{\beta y^{\sf T}Rx}d\sigma(x)$$

**`R∈O(4)` 保测度 ⇒ `Z(y)` 与 `y` 无关 ⇒ `D²(−log Z)=0`。**

**数值：三组随机 `R`，每组 6 个随机 `y`、`4×10⁵` 个球面样本，`Z(y)/mean` 的相对散布 `2.3e−3 / 2.6e−3 / 1.2e−3` ⇒ **与 `y` 无关 CONFIRMED**。**

$$\boxed{\ Z(y)\ \text{常数}\ \Longrightarrow\ D^2(-\log Z)=0\ \Longrightarrow\ \text{Fisher 项与球面 Hessian 项【精确抵消】}\ }$$

**⇒ 所以 `SLOW-H` 回到 `[U]`，而我上一轮的 `[F]` **撤回**。**

**⚠️ 保留为 `[F]` 的只是**那条推理**（量级比较 ⇒ 定符号），不是命题本身。**

### 5⁗⁗⁗⁗′.2 `K` 矩阵估计（`[T]`，我逐行核过）

**设定（球面协变 Hessian，非 `R⁴` 自由变量的普通 Hessian）：**

$$H_C(x_C;y)=H_C^{(0)}(x_C;\zeta)+\beta\!\!\sum_{a\in\partial C}\sum_{\substack{e=(a,c)\\ c\in C}}\!\!y_a^{\mathsf T}R_ex_c,\qquad w_C(a)=\#\{e=(a,c)\}$$

**逐块估计：**

$$|D_aH_C[u_a]|\ \le\ \beta w_C(a)|u_a|\tag{1}$$
$$a\ne b:\quad \mathrm{Hess}_{ab}V_C[u_a,u_b]=-\mathrm{Cov}_{\mu_C}(D_aH_C[u_a],D_bH_C[u_b])$$
$$\boxed{\ \|\mathrm{Hess}_{ab}V_C\|_{\rm op}\ \le\ \beta^2w_C(a)w_C(b)\ }\tag{2}$$
$$\boxed{\ \|\mathrm{Hess}_{aa}V_C\|_{\rm op}\ \le\ \beta w_C(a)+\beta^2w_C(a)^2\ }\tag{3}$$

**（`a≠b` 时 `D²H_C` 的混合项为零——因 `H_C` 无慢—慢边；对角用 `y_a''=−|u_a|²y_a`。）**

$$W_{ab}=\sum_Cw_C(a)w_C(b),\quad D_a=\sum_Cw_C(a),\quad m_{ab}=\text{慢—慢边重数}$$

$$\boxed{\ K_{ab}=\begin{cases}\beta m_{ab}+\beta^2W_{ab},&a\ne b\\ \beta d_{SS}(a)+\beta D_a+\beta^2W_{aa},&a=b\end{cases}\ }\tag{4}$$
$$\big|\mathrm{Hess}S_{\rm slow}(u,u)\big|\ \le\ \sum_{a,b}K_{ab}|u_a||u_b|\tag{5}$$
$$\boxed{\ \|\mathrm{Hess}S_{\rm slow}\|_{\rm op}\ \le\ \|K\|_{2\to2}\ \le\ \sup_a\sum_bK_{ab}\ }$$

### 5⁗⁗⁗⁗′.3 带空间权重的有限范围估计（`[T]`）

$$\|K\|_{\alpha,\rm row}=\sup_a\sum_be^{\alpha d(a,b)}K_{ab}$$
$$\mathrm{diam}(\partial C)\le R_\ell,\quad B_C\le c_0\Delta\ell^3,\quad D_a\le q\Delta$$

$$\boxed{\ \|K\|_{\alpha,\rm row}\ \le\ 2\beta\Delta_{SS}+\beta q\Delta+c_0\beta^2q\Delta^2\ell^3e^{\alpha R_\ell}\ }\tag{7}$$
$$\boxed{\ \alpha=0:\quad \|\mathrm{Hess}S_{\rm slow}\|_{\rm op}\ \le\ 2\beta\Delta_{SS}+\beta q\Delta+c_0\beta^2q\Delta^2\ell^3\ }\tag{8}$$

**数值（`β=2.5`、`Δ_SS=Δ=8`、`q=4`、`c₀=1`）：**

| `ℓ` | `2βΔ_SS` | `βqΔ` | `c₀β²qΔ²ℓ³` | `‖K‖_∞` |
|---|---|---|---|---|
`1` | `40` | `80` | **`1600`** | **`1720`** |
`2` | `40` | `80` | `12800` | `12920` |
`3` | `40` | `80` | `43200` | `43320` |
`4` | `40` | `80` | `102400` | `102520` |

**⇒ 有效作用量的非对角耦合只出现在**共享快分量边界**的慢变量之间 ⇒ 作用范围被块直径控制。**

**⚠️ 但 `‖Hess S_slow‖ ≤ O(ℓ³)` 是**模长上界**，`ℓ=1` 时已 `~1.7×10³`。它**不含符号信息**，**不能**用来判凹/凸；也**不是下界**。**

### 5⁗⁗⁗⁗′.4 而 `(7)` 与关联衰减之间的距离

$$\text{(9)}\quad \big|\mathrm{Cov}_{\nu_{\rm slow}}(F_A,F_B)\big|\ \le\ C\|F_A\|_\infty\|F_B\|_\infty e^{-\alpha d(A,B)}$$

**⚠️ `(7)` 只控有效作用量的**导数耦合**，**不自动推出** `(9)`。**

**⚠️ 且 `(7)` 的粗上界随 `ℓ³` **增长**，**不能**作为对所有块尺度统一成立的小耦合条件。**

$$\boxed{\ \text{本轮完成：精确 }S_{\rm slow}\ \to\ \text{显式带权 Hessian 矩阵 }K\ ;\qquad \text{未完成：}K\to(9)\ }$$

### 5⁗⁗⁗⁗′.5 审计（本轮）

| 项 | 判决 |
|---|---|
`EFDS` 快协方差恒为零 | **`[T]`** |
精确慢作用量 | **`[T]`** |
**`K` 逐块估计 `(2)–(3)`** | **`[T]`** |
**带权块矩阵估计 `(5)–(8)`** | **`[T]`** |
**`SLOW-H` 的 `[F]`** | **撤回 → `[U]`（勘误 31）** |
由 `K` 推出 `(9)` / `SECC` | `[U]` |
由 `K` 推出慢测度唯一性 | `[U]` |
**`[F]`（仅那条推理：量级比较 ⇒ 定符号）** | `[F]` |

**⚠️ 记账纪律：`[F]` 只能挂在**推理**上，不能挂在**命题**上——除非给出反例或矛盾。**

--- 

## §5⁗⁗⁗⁗″ `K → tanh → Dobrushin` 闭合为定理，但阈值被压到 `β ≲ 0.05`（第 270 轮）

### 5⁗⁗⁗⁗″.1 核心引理（`[T]`，我数值核了 `tanh` 转换）

$$\sup_y\big\|\mathrm{Hess}_{ab}S_{\rm slow}\big\|_{\rm op}\le K_{ab},\qquad K=K^{\sf T}\succeq0\ (\text{只需非对角块})$$

$$\boxed{\ c_{ab}\ \le\ \tanh\!\Big(\frac{\pi^2}{4}K_{ab}\Big)\ }\tag{3}$$

**证明骨架**：`φ(u)=S_slow(u,y_b,·)−S_slow(u,y'_b,·)`；沿两条最短测地线积分两次得 `|φ(u)−φ(v)| ≤ K_ab d(u,v)d(y_b,y'_b) ≤ π²K_ab`（用 `diam(S³)=π`）；对数密度比振幅 `≤π²K_ab`；再用两正密度 `osc log(q/p)≤R ⇒ ‖p−q‖_TV ≤ tanh(R/4)`。∎

**`tanh(R/4)` 我核了（两点分布精确比，`p=(½,½)`、`q` 取 logistic）：**

| `t` | 真实 `TV` | `tanh(t/4)` |
|---|---|---|
`0.5` | `0.122459` | `0.124353` ✓ |
`1.0` | `0.231059` | `0.244919` ✓ |
`4.0` | `0.482014` | `0.761594` ✓ |
`10.0` | `0.499955` | `0.986614` ✓ |

**⇒ 合法（偏松）界；`R→∞` 时趋于 `1` ✓**

**⭐ 优点：不要求 `K_ab` 小、不对 `β` 做级数展开 ⇒ **非微扰**影响估计。**

### 5⁗⁗⁗⁗″.2 路径和与协方差（`[T]`，条件 `(DOB)` 下）

$$\mathcal C_{ab}=\tanh(\pi^2K_{ab}/4)\ (a\ne b),\quad \mathcal C_{aa}=0;\qquad \mathcal D=(I-\mathcal C)^{-1}=\sum_{n\ge0}\mathcal C^n$$
$$\rho_\alpha(K)=\sup_a\sum_{b\ne a}e^{\alpha d(a,b)}\tanh(\pi^2K_{ab}/4)\tag{8}$$
$$\text{(DOB)}\quad \exists\alpha\gt0:\ \rho_\alpha(K)\lt1$$
$$\sup_a\sum_be^{\alpha d(a,b)}\mathcal D_{ab}\ \le\ \frac1{1-\rho_\alpha(K)}\tag{9}$$

**（次乘性：`e^{αd(a,c)} ≤ e^{αd(a,b)}e^{αd(b,c)}` ⇒ `sup_aΣ_c e^{αd(a,c)}(𝒞^n)_{ac} ≤ ρ_α^n`。）**

$$|\mathrm{Cov}_\nu(F_A,G_B)|\ \le\ \frac14\sum_{a,b}\delta_a(F_A)\mathcal D_{ab}\delta_b(G_B)\tag{10}$$
$$\boxed{\ |\mathrm{Cov}_\nu(F_A,G_B)|\ \le\ \frac{\min(|A|,|B|)}{1-\rho_\alpha(K)}\|F_A\|_\infty\|G_B\|_\infty e^{-\alpha d(A,B)}\ }\tag{12}$$

**同一条件给边界稳定性与唯一性 `(15)`。**

$$\text{桥权重被严格实现：}\quad \boxed{\ \mathfrak w_K(\gamma)=\prod_{(a,b)\in\gamma}\tanh\Big(\frac{\pi^2}{4}K_{ab}\Big)\ }\tag{14}$$

**⇒ 这正是 §5⁗ 的 `BR`/`BR-E` 的一个**严格实现**，而不是抽象待证命题。**

### 5⁗⁗⁗⁗″.3 【本轮定量判决】`tanh` 上界在 `β=2.5` 处把阈值压到 `β ≲ 0.05`

$$\frac{\pi^2\beta}{4}\Big|_{\beta=2.5}=6.1685\ \Longrightarrow\ \tanh=0.9999912273$$
$$\frac{\pi^2\beta^2}{4}\Big|_{\beta=2.5}=15.4213\ \Longrightarrow\ \tanh=1.000000000000$$

**⇒ 每个非零耦合系数已 `≈1` ⇒ 行和只要有两项就 `>1`。**

**而要 `ρ_α<1`（最省情形：一条慢—慢边 `K_ab=β`、`2d=8` 个邻居）：**

| `α` | 需 `tanh(π²β/4) <` | `β` 上界 |
|---|---|---|
`0.00` | `0.12500` | **`< 0.0509`** |
`0.05` | `0.11890` | `< 0.0484` |
`0.10` | `0.11310` | `< 0.0460` |
`0.20` | `0.10234` | `< 0.0416` |

$$\boxed{\ \text{tanh 界给 }\beta\lesssim0.05\qquad\text{而 }\beta_{\rm pin}=1.01092,\ \text{靶 }2.5\ }$$

| 阈值 | `2.5/阈值` |
|---|---|
**`tanh` 界 `0.051`** | **`49×`** |
`β_pin` `1.01092` | `2.47×` |

**⇒ **`tanh` 上界比 pinning 障碍还紧 `20×`** ⇒ 这条上界路线**损失巨大**。**

**⚠️ 但这**只说明这条上界不足以过收缩判据**，**不能**据此断言真实慢测度不聚簇——`K` 是混合 Hessian 的**绝对值**上界，可能丢掉了 `D²(−log Z_C)` 里的抵消结构。**

### 5⁗⁗⁗⁗″.4 【`[T]` 反例】有限范围 + 有界混合 Hessian `⇏` 指数聚簇

$$S_J(y)=-J\sum_{\langle a,b\rangle}y_a^{\sf T}y_b\ \Longrightarrow\ \|\mathrm{Hess}_{ab}S_J\|_{\rm op}\le J\ (\text{邻点}),\ 0\ (\text{非邻点})$$

**而 4D `O(4)` 在 `J` 足够大时**长程有序**（相变）。**

$$\boxed{\ \text{有限范围}\ +\ \sup_{a,b}\|\mathrm{Hess}_{ab}S\|\lt\infty\ \not\Longrightarrow\ \text{指数聚簇}\ }$$

**⇒ 所以**收缩条件 `(DOB)` 是逻辑上必需**的——不是技术性保守。`(DOB)` 是**充分条件，不是必要条件**。**

### 5⁗⁗⁗⁗″.5 【记账】`ζ` 必须被积分进来，否则只到条件版本

**若 `ζ` 只作固定外部参数：只得到 `Cov(F_A,F_B\mid\zeta)`。**
**要控 `Cov(𝔼[F_A\mid\zeta],𝔼[G_B\mid\zeta])`，**必须把 `ζ` 纳入慢变量图与影响矩阵 `𝒞`**；否则 `SECC` 只是**条件**版本，不是完整慢测度聚簇。**

### 5⁗⁗⁗⁗″.6 本轮状态

| 命题 | 判 |
|---|---|
`K_ab` 控单点总变差影响 `(3)` | **`[T]`** |
影响矩阵加权路径和有界 `(9)`（条件 `ρ_α<1`） | **`[T]`** |
慢测度指数协方差衰减 `(12)`（同条件） | **`[T]`** |
慢测度边界稳定性 / 唯一性 `(15)`（同条件） | **`[T]`** |
**有限范围 + 有界混合 Hessian ⇏ 聚簇** | **`[T]`（反例）** |
**`β=2.5` 时粗 `K` 满足 `(DOB)`** | **`[U]`**（`tanh` 界给 `β≲0.05`，差 `49×`） |
**`β=2.5` 时真实慢测度指数聚簇** | **`[U]`** |
`ζ` 纳入与否的记账 | **`[T]`** |

### 5⁗⁗⁗⁗″.7 下一步（我同意你的方向，并把损失量化了）

**不用 `tanh(π²K/4)` 这个粗上界，直接估 `D²(−log Z_C)` 的条件分布影响系数。**

**⇒ 收益估计：`tanh` 界损失约 `49/2.47 ≈ 20×`（相对 pinning 障碍）。若能把这 `20×` 收回来，阈值可回到 `β ~ 1` 量级；而靶是 `2.5`，仍差 `2.5×`。**

**⇒ 所以即使 `K → tanh` 这条链完全优化，也到不了 `2.5`——这意味着**必须换判据**，而不是继续精化 `K`。**

--- 

## §5⁗⁗⁗⁗‴ 命题 A `[F]`（差 `O(ℓ³)`，非边缘）+ 命题 B 仍 `[U]`（第 271 轮）

### 5⁗⁗⁗⁗‴.1 【`[F]`】粗 `K` 的 Dobrushin 收缩在 `β=2.5` 失败

$$K_{ab}=\beta m_{ab}+\beta^2W_{ab},\quad m_{ab},W_{ab}\in\mathbb Z_{\ge0}\ \text{不同时为零}\ \Longrightarrow\ K_{ab}\ge\min(\beta,\beta^2)$$

$$\beta=2.5:\quad \mathcal C^K_{ab}\ \ge\ \tanh\Big(\frac{\pi^2}{4}\cdot2.5\Big)=0.9999912273\tag{1}$$

$$\sum_{b\ne a}e^{\alpha d(a,b)}\mathcal C^K_{ab}\ \ge\ 0.9999912273\,n_a\qquad(e^{\alpha d}\ge1)$$

$$\boxed{\ n_a\ge2\ \Longrightarrow\ \rho_\alpha(K)\gt1\quad\forall\alpha\ge0\ }\tag{2}$$

**三边界接触判据：若某快分量 `C` 同时接触 `a,b,c`，则 `W_ab,W_ac≥1` ⇒ `K≥β²=6.25`：**

$$1-\tanh\Big(\frac{\pi^2}{4}\cdot6.25\Big)=\frac{2}{e^{2x}+1}=8.06\times10^{-14}\ \checkmark\ (\text{与你给的 }8.0\times10^{-14}\text{ 相符})$$

$$\boxed{\ \sum_{b'\ne a}\mathcal C^K_{ab'}\ \ge\ 2-1.6\times10^{-13}\ \gt\ 1\ }\tag{3}$$

### 5⁗⁗⁗⁗‴.2 【我的补算】`n_a` 不是 `2`，而是 `O(ℓ³)`

$$\boxed{\ \text{机制：}-\log Z_C\ \text{依赖【整个边界 }\partial C\text{】}\ \Longrightarrow\ \text{对 }\partial C\text{ 上任意两慢变量都产生耦合}\ }$$

**⇒ `W_ab>0 ⟺ a,b` 共享某个快分量 `C`。分隔集构造下，慢变量 `a` 落在块边界上，接触 `≥2` 个块；而每个块边界 `|∂C|=O(ℓ³)` 个点，**全部与 `a` 耦合**：**

| `ℓ` | `\|∂C\|` | 行和下界 `≈0.99999·n_a` |
|---|---|---|
`1` | `~8` | **`8`** |
`2` | `~64` | **`6\times10^{1}`** |
`3` | `~216` | **`2\times10^{2}`** |
`4` | `~512` | **`5\times10^{2}`** |

$$\boxed{\ \text{行和超 }1\text{ 的倍数}\ \sim\ O(\ell^3)\text{——不是边缘性的 }2\times\ }$$

### 5⁗⁗⁗⁗‴.3 结构表述：块积分把**块边界耦合成全连接**

$$\text{块积分}\ \Longrightarrow\ -\log Z_C\ \text{在}\ \partial C\ \text{上是【全连接】}\ \Longrightarrow\ \text{有效慢相互作用【稠密】}$$

**而 KP/Dobrushin 恰好要求**稀疏 + 小活动**。**

$$\boxed{\ \text{这与前几轮的 }\ell^{10}\text{（}W\text{ 的稠密 }\times|X|^2\text{ 权重）是同一现象的两个面}\ }$$

**⇒ 所以命题 A 的失败不是常数不够细，是**块积分本身产生的稠密耦合**与 KP 要求不相容。**

### 5⁗⁗⁗⁗‴.4 退火版 `(6)` 能不能救？—— 也要典型 `K` 小

$$\text{(6)}\quad |\mathrm{Cov}|\ \le\ C\|F_A\|\|F_B\|\sum_{\Gamma:A\leftrightarrow B}\mathbb E_{\nu_{\rm slow}}\Big[\prod_{e\in\Gamma}Q_e(y)\Big]$$

**关键：`c_ab` 是邻居对 `a` 的条件分布的影响。局部场 `h=βΣ_wR_e x_w`，随机方向下 `|h|~β√(2d)=β√8`：**

| `β` | `|h|` | `tanh(π²β/4)` |
|---|---|---|
`0.05` | `0.14` | `0.122748` |
`1.00` | `2.83` | `0.985719` |
**`2.50`** | **`7.07`** | **`0.999991`** |

**⇒ `β=2.5` 时 `|h|~7.1`，条件分布强峰 ⇒ `c_ab≈1` ⇒ **退火平均也 `≈1`。****

$$\boxed{\ \text{退火只有在【典型构型的局部曲率小】时才帮忙；}\beta=2.5\text{ 的典型曲率是 }O(1)\ }$$

**⇒ `(6)` **不是自动更好**——它同样要求典型 `K` 小。⚠️ 且 `(6)` 是**新目标命题**，不是已证的协方差不等式：还须证逐步条件更新的影响传播界，并控制路径权重的**联合期望**（单点平均小不够）。**

### 5⁗⁗⁗⁗‴.5 两个命题的逻辑关系（不得合并）

$$\boxed{\ A\Longrightarrow B,\qquad \neg A\ \not\Longrightarrow\ \neg B\ }$$

| 命题 | 状态 | 原因 |
|---|---|---|
**A：粗 `K` 在 `β=2.5` 满足 Dobrushin** | **`[F]`** | 式 `(2)` 直接反证；且 `n_a=O(ℓ³)` ⇒ 差 `O(ℓ³)` |
**B：真实慢测度在 `β=2.5` 指数聚簇** | **`[U]`** | 既未由粗 `K` 证明，也未被粗 `K` 失败否定 |

**⇒ `¬A` 只否定**用当前粗 `K` 走 Dobrushin**这条路线，**不**否定慢测度本身聚簇。**

**★ 可能的其它机制（B 可经由它们成立）：随机块收缩 / 多步或块条件分布收缩 / 非均匀 Dobrushin / 连通坏块概率衰减 + 好块条件混合 / 慢测度谱隙或传输不等式。**

### 5⁗⁗⁗⁗‴.6 `β_pin` 的逻辑作用必须说准（`[T]` 记账）

**若真能构造一组固定慢变量的边界配置、使条件模型中出现**有序的三维 `O(4)`**，则它排除的是：**

$$\boxed{\ \text{对所有固定慢边界配置【一致】成立的条件收缩估计}\ }$$

**⚠️ 但完整慢测度是对慢变量构型积分后的**退火**测度，**不等同**于任意人为钉扎的条件测度：**

$$\text{存在坏的 pinning 配置}\ \not\Longrightarrow\ \text{退火慢测度长程关联}$$

**⚠️ 且要把 `β_pin` 用作**严格障碍**，还须证明该三维模型确实是原有效模型的**合法条件子系统**、且相关边界条件确实属于条件测度的**支持范围**——不能只凭物理类比把数值阈值当定理。**

### 5⁗⁗⁗⁗‴.7 下一步：换判据，不再优化最坏情形行和

**⇒ 目标不再设为 `ρ_α(𝒞^K)<1`（该路线已出现结构性障碍）。**

**⇒ 改为证明**退火的、分块的影响传播界** `(6)`，即把**

$$\prod_{e\in\Gamma}\sup_yQ_e(y)\quad\longrightarrow\quad\mathbb E_{\nu_{\rm slow}}\Big[\prod_{e\in\Gamma}Q_e(y)\Big]$$

**⚠️ 但本轮我算的 `|h|~7.1` 说明：在 `β=2.5` 的典型构型上 `Q_e≈1`，所以**退火本身不足以翻盘**；还需要一个「典型局部曲率小」的机制（例如把 `ζ` 的挫折一起退火掉，或换观测量）。**

--- 

## §5⁗⁗⁗⁗⁗ ADP：退火路径权重的严格版本 + `(18)` 的定量衰减率（第 272 轮）

### 5⁗⁗⁗⁗⁗.1 三个环节的判决（`[T]`/`[T]`/`[U]`/`[T]`）

**⚠️ 先记一条纪律：**式 `(6)` 不能仅靠把确定性 `c_ab` 换成随机 `Q_{ab}(y)` 得到**——路径乘积必须对应**真实的条件影响传播过程**；否则把不同更新步骤的随机影响写成同一个构型 `y` 上的乘积，是未经证明的替换。**

**① 局部影响权重（`[T]`）：**

$$\Delta_{a\leftarrow b}(z)=\sup_{\substack{u,v\in S^3\\ s,s'\in S^3}}\big|S(u,s,z)-S(u,s',z)-S(v,s,z)+S(v,s',z)\big|$$

$$\boxed{\ Q_{b\to a}(y)=\tanh\Big(\frac{\Delta_{a\leftarrow b}(y_{\widehat{a,b}})}4\Big)\ }\tag{4}$$

$$\Big\|\gamma_a(\cdot\mid y_b{=}s,z)-\gamma_a(\cdot\mid y_b{=}s',z)\Big\|_{\rm TV}\ \le\ Q_{b\to a}(y)\tag{5}$$

**（对数比振幅 `≤Δ`，再用 `‖p−q‖_TV ≤ tanh(osc log(p/q)/4)`。）**

**与 `K` 的关系：沿两条最短测地线积分 ⇒ `Δ ≤ π²K_ab` ⇒ `Q ≤ tanh(π²K_ab/4)` `(6)`。**

**⭐ `Q` 保留了 Hessian 的**局部环境**信息；原来的 `K` 只保留**全局最坏情形**上界。**

**② 退火影响图（`[T]`）：`P(ω_{b→a}=1|Y)=Q_{b→a}(Y)`，条件独立；**

$$\boxed{\ \mathbb P_{\rm ann}(\gamma\ \text{全开})=\mathbb E_{\nu_{\rm slow}}\Big[\prod_{e\in\gamma}Q_e(Y)\Big]\ }\tag{9}$$
$$\mathbb P_{\rm ann}(A\leftrightarrow B)\ \le\ \sum_{\gamma:A\leftrightarrow B\ \rm simple}\mathbb E_{\nu_{\rm slow}}\Big[\prod_{e\in\gamma}Q_e(Y)\Big]\tag{10}$$

**③ ADP（`[U]`，真正的缺口）：**存在条件耦合 `Π_{b,b'}`（边缘为 `ν_slow(·|Y_B=b)`、`ν_slow(·|Y_B=b')`），使 `X_A≠X'_A ⇒ A↔B`，且**

$$\int\nu_B(db)\nu_B(db')\Pi_{b,b'}(A\leftrightarrow B)\ \le\ \mathbb P_{\rm ann}(A\leftrightarrow B)\tag{ADP}$$

**④ ADP ⟹ 协方差（`[T]`，我核过）：**

$$|\mathrm{Cov}(F_A,G_B)|\ \le\ \|G_B\|_\infty\mathbb E|m(Y_B)-\mathbb E m|\ \le\ \|G_B\|_\infty\mathbb E|m(Y_B)-m(Y'_B)|$$
$$\le\ 2\|F_A\|_\infty\|G_B\|_\infty\,\mathbb P_{\rm ann}(A\leftrightarrow B)$$

$$\boxed{\ |\mathrm{Cov}_{\nu_{\rm slow}}(F_A,F_B)|\ \le\ 2\|F_A\|_\infty\|F_B\|_\infty\!\!\sum_{\gamma:A\leftrightarrow B\ \rm simple}\!\!\mathbb E_{\nu_{\rm slow}}\Big[\prod_{e\in\gamma}Q_e(Y)\Big]\ }\tag{15}$$

**⇒ `C=2`。**

### 5⁗⁗⁗⁗⁗.2 为何不能直接用单点热浴最大耦合冒充 ADP（`[T]` 记账）

**逐坐标望远镜展开只给**

$$\Big\|\gamma_a(\cdot\mid x_{-a})-\gamma_a(\cdot\mid x'_{-a})\Big\|_{\rm TV}\ \le\ \sum_{b:x_b\ne x'_b}Q_{b\to a}(z^{(b)})\tag{16}$$

**⚠️ 而 `z^{(b)}` 是**依赖望远镜展开次序的混合环境**，**不一定服从 `ν_slow`**。**

$$\boxed{\ \text{不能把 (16) 各步的 }Q\text{ 直接替换成同一个 }Y\sim\nu_{\rm slow}\text{ 上的 }Q_e(Y)\text{ 再取期望}\ }$$

**⇒ 这是此前推导不能跳过的一步。单靠 Hessian 上界、`E[Q_e]` 单点平均、甚至单点 TV 估计，**都不足以**保证它。**

### 5⁗⁗⁗⁗⁗.3 【我的定量补充】`(18)` 在 `β=2.5` 的衰减率

**先算 `H_SS` 那一项：`S(u,s)−S(u,s')=−βu^TR(s−s')`，对 `u∈S³` 的振荡 `=2β|s−s'|` ⇒ `Δ/4=β|s−s'|/2`。**

**两个独立 `S³` 样本的弦距（`2×10⁶` 样本）：`mean = 1.3584`、`median = 1.4143`。**

$$Q_e=\tanh(\beta|s-s'|/2)\Big|_{\beta=2.5}:\quad \mathbb E[Q_e]=0.8991,\quad |\ln\mathbb E Q_e|=0.1064\ \text{每边}$$

**路径乘积（模拟，长度 `2/4/8`）：**

| 路径长 | `E[Π Q_e]` | `(E[Q])^n` | 比 |
|---|---|---|---|
`2` | `8.080e−1` | `8.083e−1` | **`1.00`** |
`4` | `6.532e−1` | `6.533e−1` | **`1.00`** |
`8` | `4.266e−1` | `4.269e−1` | **`1.00`** |

$$\Longrightarrow\quad k:=-\frac{\ln\mathbb E[\prod Q_e]}{n}=0.1064\ \text{每边}\quad(\text{关联长度}\approx9.4\ \text{块})$$

**⭐⭐ 反直觉的一点：在**乘积（均匀独立顶点）**测度下，路径上相邻边权**恰好独立**——因给定共享顶点 `v_{i+1}` 后 `v_i` 与 `v_{i+2}` 独立，而 `|v_i−v_{i+1}|` 的边缘与 `v_{i+1}` 无关。所以 `E[ΠQ]=(EQ)^n` **精确成立**。**

**⚠️ 而 `ν_slow` **不是**乘积测度（`H_SS` 使相邻顶点铁磁关联）⇒ 边权**有**相关。但相关方向是**铁磁**：相邻顶点更对齐 ⇒ `|s−s'|` 更小 ⇒ `Q` 更小 ⇒ **对我们更有利**。**

$$\boxed{\ \text{所以 }k\ \text{的真值}\ \ge\ 0.106\text{（本估计是保守的）}\ }$$

**⚠️ 但本估计**只用 `H_SS`**；`−Σlog Z_C` 也贡献 `Δ`，而它只会**增大** `Q`、**减小** `k`。**

$$\boxed{\ \text{故 }k\ \le\ 0.106\quad(\text{关联长度}\ \gtrsim\ 9.4\ \text{块，可能更大})\ }$$

**⇒ `(18)` 可满足 `α < k`，取 `α ≈ 0.07` ⇒ 关联长度 `≈14` 块 `≈14ℓ` 格距。**

**⇒ 有限关联长度**存在**，但是**弱指数**——与「`β=2.5` 无小参数」的诊断一致。**

### 5⁗⁗⁗⁗⁗.4 下一步的唯一目标

$$\boxed{\ \text{(17)}\quad \text{相邻慢变量的热浴更新，其不一致传播由 }Q_{b\to a}(Y)\text{ 控制，且沿传播链的联合权重由 }\nu_{\rm slow}\text{ 下的【同构型乘积期望】控制}\ }$$

**⚠️ 必须控**联合乘积** `E[Q_{e_1}⋯Q_{e_n}]`，不能只证 `E Q_e ≪ 1`——相邻边权可能相关，单点平均小不意味着路径乘积按长度指数衰减。**

**⇒ 本轮我核到：在本模型里相关方向是**铁磁（有利）**，而最大风险来自 `log Z_C` 的贡献，它可能把 `k` 压小。**

### 5⁗⁗⁗⁗⁗.5 审计

| 项 | 判 |
|---|---|
从 `S_slow` 定义 `Q_e(y)`，式 `(3)–(5)` | **`[T]`** |
`Q ≤ tanh(π²K/4)` | **`[T]`** |
路径权重 = 退火图开路径概率，`(8)–(10)` | **`[T]`** |
**ADP：构造与退火图匹配的条件 Gibbs 耦合** | **`[U]`（唯一的实质缺口）** |
ADP ⟹ 协方差路径界，`(11)–(15)`，`C=2` | **`[T]`** |
**(16) 的混合环境不受 `ν_slow` 支配** | **`[T]`（记账）** |
`H_SS` 项给的 `k ≈ 0.106`（关联长度 `≈9.4` 块） | **`[M]`（模拟）** |
含 `log Z_C` 后 `k ≤ 0.106` | **`[T]`（上界推理）** |
**`β=2.5` 下 ADP** | **`[U]`** |

---

## §5⁗⁗⁗⁗⁗′ 【`[F]`】ADP 对一般 Gibbs 测度为假：三变量反例（第 273 轮）

### 5⁗⁗⁗⁗⁗′.1 我撤回上一轮对 ADP 的处理

**上一轮我把 ADP 标成 `[U]`（未证但可能真）。**实际它**错**：按上一轮的定义，ADP 对一般有限范围 Gibbs 测度**不成立**。**

$$\boxed{\ \text{ADP}\ =\ [F]\ \text{（一般 Gibbs 测度）；对当前 Wilson }S_{\rm slow}\ \text{则仍为 }[U]\ }$$

### 5⁗⁗⁗⁗⁗′.2 必要条件（`[T]`）

**任何耦合都满足 `Π(X_A≠X'_A) ≥ ‖ν(X_A∈·|b)−ν(X_A∈·|b')‖_TV`，故 ADP 至少要求**

$$\boxed{\ \mathbb E_{b,b'}\big\|\nu_A(\cdot\mid b)-\nu_A(\cdot\mid b')\big\|_{\rm TV}\ \le\ \mathbb P_{\rm ann}(A\leftrightarrow B)\ }\tag{19}$$

### 5⁗⁗⁗⁗⁗′.3 三变量反例（我逐项数值核过，全部相符）

$$a,b,c\in\{0,1\};\qquad S(a,b,c)=Hb+M(b-c)^2+La-Jabc;\qquad L=10,\ J=20,\ M=20$$

**交叉振幅（我修正了自己的枚举 bug 后核）：**

| 量 | 值 |
|---|---|
`Δ_{a←b}(c=0)` | **`0`** |
`Δ_{a←b}(c=1)` | **`20 = J`** |
`Δ_{b←c}(a)` | **`60 = 2M+Ja`** |

$$\Longrightarrow\quad Q_{ab}=\tanh(Jc/4),\quad Q_{ac}=\tanh(Jb/4),\quad Q_{bc}=\tanh((2M+Ja)/4)$$

$$\boxed{\ b=c=0\ \Longrightarrow\ Q_{ab}=Q_{ac}=0\ \Longrightarrow\ \text{退火图里 }a\text{ 是孤立点}\ }$$

**⇒ `P_ann(a↔b) ≤ ν(b=1 或 c=1)`。**

**配 `e^{−H}=(1+r)(1+s)/\{99[s(1+r)+1+t]\}`、`r=e^{−L}`、`s=e^{−M}`、`t=e^{J−L}`，**严格得 `ν(b=1)=0.01`**：**

$$\text{我算：}e^{-H}=4.585852\times10^{-7},\ H=14.595120,\ \nu(b=1)=0.010000\ \checkmark$$

$$\boxed{\ \mathbb P_{\rm ann}(a\leftrightarrow b)\ \le\ 0.01+0.99\tfrac{s}{1+s}=0.010000002041\ }\tag{25}$$

**而条件分布：**

$$p_0=\nu(a{=}1|b{=}0)=\tfrac{r}{1+r}=0.000045398,\qquad p_1=\nu(a{=}1|b{=}1)=0.999954602$$
$$\big|p_1-p_0\big|=0.999909204\qquad(\text{我逐位相符})$$

$$\boxed{\ \mathbb E_{b,b'}\Pi(a\ne a')\ \ge\ 2(0.01)(0.99)|p_1-p_0|=0.019798202\ }\tag{27}$$

**⇒ 与 `(25)` 矛盾：**

$$\frac{\text{需要}}{\text{可用}}=\frac{0.019798202}{0.010000002}=\mathbf{1.980}$$

$$\boxed{\ \text{需要超出 }98\%\ \Longrightarrow\ \text{ADP}\ =\ [F]\ }$$

### 5⁗⁗⁗⁗⁗′.4 反例揭示的确切错误（机制）

$$\text{单点估计 }(29)\ \text{没错};\qquad \text{错在把混合环境 }z^{(b)}\ \text{都替换成同一个 }Y\sim\nu_{\rm slow}$$

**机制：`b=1` 时 `c=1` 几乎必然（`M(b−c)²` 强耦合）；于是 `a` 与 `b` 的条件分布几乎相反；但**典型** Gibbs 构型有 `b=c=0` 的概率 `≈0.99`，而在那些构型上 `a` 的两条边**都关**。**

$$\boxed{\ \text{一份【典型】构型上的退火图，低估了两种【不同边界条件】之间的影响}\ }$$

### 5⁗⁗⁗⁗⁗′.5 该机制在 Wilson 慢作用量上的**结构对应**（`[M]` 推断，非定理）

| 反例中的角色 | Wilson `S_slow` 中的对应 |
|---|---|
`M(b−c)²` | `H_SS` 的 `β(1−y_a^TRy_b)` 型双线性耦合 |
`−Jabc` | **`−log Z_C`** 的**多边界**耦合 |
`L a`、`Hb` | 各边界变量的局部场 |

**⇒ 所以反例的**结构**（"强局部耦合使条件分布翻转，而典型构型上边权为 0"）在 `S_slow` 上有对应物。**

**⚠️ 但这不是证明**Wilson 一定违反 ADP**——`S_slow` 的精确积分结构可能提供额外性质。**目前只是说：不能假定该机制不存在。**

### 5⁗⁗⁗⁗⁗′.6 两个合法方向

**方向一（保确定性上界）：`Q̄_{b→a} = sup_z Q_{b→a}(z)`。它控制**所有**混合环境，可进标准 Dobrushin–Föllmer；**缺点：重新遇到最坏情形障碍（§5⁗⁗⁗⁗‴ 的 `[F]`）。**

**方向二（keep annealed，但重证耦合估计）：真正要控的是混合环境乘积**

$$\boxed{\ \mathbb E_{\Pi_{b,b'}}\Big[\prod_{j=1}^{n}Q_{e_j}(Z^{(j)})\Big]\ }\tag{31}$$

**要把它化回 `E_{\nu_{\rm slow}}[Π_jQ_{e_j}(Y)]`，**必须对 `S_slow` 另证一个「混合环境支配不等式」**。这才是真正要研究的命题。**

### 5⁗⁗⁗⁗⁗′.7 审计（本轮）

| 命题 | 判 |
|---|---|
单坐标条件分布影响由 `Q_{b→a}` 控制，式 `(5)` | **`[T]`**（仍成立） |
**一般 Gibbs 测度满足上一轮的 ADP** | **`[F]`（三变量反例）** |
**当前 Wilson `S_slow` 满足 ADP** | **`[U]`**（需用其特殊结构另证） |
ADP ⟹ 协方差路径界 `(11)–(15)` | **`[T]`**（条件命题，推导本身成立） |
反例机制在 `S_slow` 上的结构对应 | **`[M]`**（推断） |
`(31)` 的混合环境支配不等式 | **`[U]`（新的唯一目标）** |

--- 

## §5⁗⁗⁗⁗⁗″ Wilson 因子层展开：`(4)` `[T]`、`(8)` vs `(9)` 的缺口 `[T]`、环境依赖 `[T]`（第 274 轮）

### 5⁗⁗⁗⁗⁗″.1 因子层精确恒等式（`[T]`）

$$H_C(x_C;y)=H_C^{(0)}(x_C;\zeta)+\beta_W\!\!\sum_{a\in\partial C}\!\!y_a^{\sf T}J_{C,a}(x_C),\qquad J_{C,a}=\!\!\sum_{\substack{e=(a,c)\\ c\in C}}\!\!R_ex_c$$

$$\boxed{\ Z_C(y)=Z_C^{(0)}\,\mathbb E_{\mu_C^{(0)}}\exp\Big(\beta_W\!\!\sum_{a\in\partial C}\!\!y_a^{\sf T}J_{C,a}\Big)\ }\tag{3}$$

$$\boxed{\ \mathrm{Hess}_{ab}S_{\rm slow}(y)=-\beta_W\!\!\sum_{e\in E_{SS}(a,b)}\!\!R_e\ -\ \beta_W^2\!\!\sum_{C:\,a,b\in\partial C}\!\!\mathrm{Cov}_{\mu_{C,y}}(J_{C,a},J_{C,b})\ }\tag{4}$$

**⇒ 有效慢变量的混合 Hessian = **快变量边界力的协方差之和 + 慢—慢边的 `−βR_e`**。**

### 5⁗⁗⁗⁗⁗″.2 【我的核验】单个快顶点 + `R_e=I` 是 Bessel 型，且 `−log Z_C` **凹**

$$\boxed{\ Z_C(y)=\frac{2I_1(\beta_W|h|)}{\beta_W|h|},\qquad h=\sum_{a=1}^m y_a\ }\tag{6}$$

**（`∫_{S³}e^{z·x}dσ = (2/z)I_1(z)`。）我数值核了 `F(t)=\log(2I_1(t)/t)` 的凸性：**

| `t` | `F'(t)` | `F''(t)` |
|---|---|---|
`0.5` | `0.123718` | **`0.242384`** |
`1.0` | `0.240194` | **`0.221716`** |
`2.5` | `0.507195` | **`0.134123`** |
`10` | `0.854185` | **`0.014104`** |
`20` | `0.925988` | **`0.003659`** |

**⇒ `F''>0` 恒成立 ⇒ `F` 凸 ⇒ `log Z_C=F(β|h|)` 在 `|h|` 凸；`|h|` 在 `h` 上凸；`h=Σ_ay_a` 对 `y` 线性：**

$$\boxed{\ \log Z_C\ \text{在 }y\text{ 上凸}\ \Longrightarrow\ -\log Z_C\ \text{在 }y\text{ 上【凹】}\ }$$

**（小 `t` 渐近：`F''→0.25`，故二阶导项的系数 `β²F''≈0.25β²`。）**

### 5⁗⁗⁗⁗⁗″.3 【关键】你的精确抵消例子是 `m=1` 的**退化**情形

**`m=1` 时 `h=y_a`、`|h|=1` **常数** ⇒ `F(β|h|)` 是常数 ⇒ Hessian `=0`。这就是你上一轮那个 `Z(y)` 与 `y` 无关的例子。**

**⇒ 而 `m≥2` 时 `|h|` **变化** ⇒ `F` 严格凸 ⇒ 二阶导**严格正**。**

$$\boxed{\ \text{你的反例【不】推广到 }m\ge2\text{；两点【不】矛盾——它们是 }m=1\text{ 与 }m\ge2\ }$$

**⚠️ 所以我上一轮撤回的 `[F]`（`S_slow` 凹 ⇒ `SLOW-H` 不可能）**在结论上是对的**，只是**理由**要用 Bessel 凸性，而不是那个无效的量级比较。**

**⇒ 记账：`SLOW-H` = `[F]`（新理由：`log Z_C` 凸 ⇒ `−log Z_C` 凹；由 `(4)` 的协方差项带负号）；而我上一轮那条**量级比较的推理** = `[F]`。**

**⚠️ 但环境依赖是真的：`μ_{C,y}` 依赖整个 `y_{∂C}`。`m=2` 时 `h=y_a+y_b`，Hessian 只依赖 `y_a,y_b`（**无**环境依赖）；`m≥3` 时 `h_0=\sum_{c'\notin\{a,b\}}y_{c'}\neq0`，交叉响应才有真正的环境依赖。**

$$\boxed{\ S_{\rm slow}\ \text{一般【不是】慢变量的纯两体作用量}\ }$$

### 5⁗⁗⁗⁗⁗″.4 【`[T]`】`(8)` vs `(9)`：独立边图**不**支配多邻居联合影响

$$\text{(8)}\quad\|\gamma_a(\cdot|\eta)-\gamma_a(\cdot|\eta')\|_{\rm TV}\le\tanh\Big(\tfrac14\!\!\sum_{b\in D_a}\!\!\Delta_{ab}\Big)$$
$$\text{(9)}\quad1-\prod_{b\in D_a}\big(1-\tanh(\Delta_{ab}/4)\big)$$

**反例（两条边各 `tanh(Δ/4)=½`，即 `Δ=2.197225`）：**

$$(8)=\tanh(D/2)=0.800000\qquad(9)=1-(1-\tfrac12)^2=0.75$$

$$\boxed{\ (8)\gt(9)\ \Longrightarrow\ \text{独立边图不支配多邻居联合影响}\ }$$

### 5⁗⁗⁗⁗⁗″.5 修正：确定性边概率 `(10)–(12)`（`[T]`，我核了每条不等式）

$$p_{ab}=\tanh(\bar\Delta_{ab}/2)\tag{10}$$

**`1−tanh(Δ/2)=2/(e^Δ+1) ≤ e^{−Δ/2}` ⟺ `2≤2cosh(Δ/2)` ⟺ `cosh≥1` ✓（我核了 `Δ=1/4/10`）**

**`1−e^{−x/2} ≥ tanh(x/4)` ✓（我核了 `x=0.1/1/5/10`）**

$$\boxed{\ 1-\prod_{b\in D_a}(1-p_{ab})\ \ge\ 1-e^{-\frac12\sum\bar\Delta_{ab}}\ \ge\ \tanh\Big(\tfrac14\sum\bar\Delta_{ab}\Big)\ }\tag{11}$$

**配 `Δ̄_ab ≤ π²K_ab`：**

$$p^{\rm safe}_{ab}=\tanh(\pi^2K_{ab}/2)\tag{12}$$

$$\beta=2.5:\quad \frac{\pi^2\beta}{2}=12.3370\ \Longrightarrow\ \tanh=1-3.85\times10^{-11}$$

**⇒ 极接近 `1` ⇒ **不能**给指数衰减。**⇒ 这是**正确但对目标耦合不够强**的结论。**

**⚠️ 且 `(10)` 是确定性最坏环境概率，大于原来的 `Q`，**不能替代退火平均权重**。**

### 5⁗⁗⁗⁗⁗″.6 经典构造的独立性条件不能删（`[T]` 记账）

**disagreement percolation 的构造依赖**独立随机场**支配不一致传播（van den Berg–Maes；Häggström 明确指出独立性条件不能直接删除而不增加假设）。**

**⇒ **不能**只把这种支配关系换成随机环境平均而不另证条件。**

### 5⁗⁗⁗⁗⁗″.7 真正要证的 Wilson 专属命题

$$\text{(14)}\quad \mathbb E_{\rm coupling}\Big[\prod_{j=1}^nQ_{e_j}(Z^{(j)})\Big]\ \le\ \mathbb E_{\nu_{\rm slow}}\Big[\prod_{j=1}^nQ_{e_j}(Y)\Big]$$

$$\boxed{\ \text{需要研究：快分量边界力协方差 }\mathrm{Cov}_{\mu_{C,y}}(J_{C,a},J_{C,b})\text{ 的【混合环境可控性】}\ }$$

**⚠️ `(14)` **不能**由 `(3)` 或 `(4)` 自动推出。**

**⚠️ `ζ` 的记账：若 `ζ` 固定为外部参数，须对所有允许的 `ζ` 一致；若 `ζ` 也随机，**必须纳入耦合与退火平均的概率空间**。**

### 5⁗⁗⁗⁗⁗″.8 审计（本轮）

| 命题 | 判 |
|---|---|
Wilson 因子精确表示 `(3)` | **`[T]`** |
混合 Hessian 的协方差表示 `(4)` | **`[T]`** |
**`F(t)=log(2I₁(t)/t)` 凸 ⇒ `−log Z_C` 凹** | **`[T]`（数值核，`F''>0` 于 `t=0.5…20`）** |
**`SLOW-H`** | **`[F]`（正确理由：Bessel 凸性 + `(4)` 的负号）** |
我上一轮那条「量级比较」的推理 | **`[F]`**（无效） |
`S_slow` 一般非纯两体 | **`[T]`** |
**`(8)>(9)` 的反例** | **`[T]`** |
最坏情形支配 `(10)–(12)` | **`[T]`**（但 `β=2.5` 时 `p≈1`，无用） |
**当前 `S_slow` 满足原始退火 ADP** | **`[U]`** |
**`β=2.5` 下 `Q_e(Y)` 支配混合环境传播** | **`[U]`（关键缺口）** |

--- 

## §5⁗⁗⁗⁗⁗‴ 环境重加权界 `[T]` + 确定性支配 `[T]` + `(10)` 仍 `[U]`（第 275 轮）

### 5⁗⁗⁗⁗⁗‴.1 两个 `[T]`（我逐步核过）

**① 快分量的环境重加权（`[T]`）：`|J_{C,a}| ≤ w_C(a)`（`|R_ex_c|=1` + 三角不等式）。**

$$\mathrm{osc}_{x_C}\log\frac{d\mu_{C,z}}{d\mu_{C,y}}\ \le\ 4\beta_W\!\!\sum_{a\in D(y,z)}\!\!w_C(a)$$
$$\boxed{\ \|\mu_{C,z}-\mu_{C,y}\|_{\rm TV}\ \le\ \tanh\Big(\beta_W\!\!\sum_{a\in D(y,z)}\!\!w_C(a)\Big)\ }\tag{4}$$

$$\beta_W=2.5,\ \text{仅改一个 }w=1\text{ 的边界坐标}:\quad \tanh(2.5)=0.986614$$

**⇒ 这个上界**很大** ⇒ 没有足以控制混合环境的强一致性估计。**

**⚠️ 但它只是上界，**不**证明实际 TV 接近此数；它排除的是「靠这个密度比估计拿到有用小常数」。**

**② 确定性 disagreement percolation（`[T]`）：`p_ab = tanh(Δ̄_ab/2)`，**

$$\mathrm{TV}(\gamma_a,\gamma'_a)\ \le\ 1-\prod_{b\in D_a}(1-p_{ab})\tag{12}$$
$$\boxed{\ |\mathrm{Cov}_{\nu_{\rm slow}}(F_A,F_B)|\ \le\ 2\|F_A\|_\infty\|F_B\|_\infty\!\!\sum_{\gamma:A\leftrightarrow B}\prod_{e\in\gamma}p_e\ }\tag{13}$$

**⚠️ 配 `Δ̄_ab ≤ π²K_ab` 得 `p ≤ tanh(π²K_ab/2)`，`β=2.5` 时 `= 1−3.85\times10^{-11}` ⇒ **无用于目标聚簇**。**

### 5⁗⁗⁗⁗⁗‴.2 【我的定量补充】`Δ` 随环境 `|h_0|` 单调**衰减** ⇒ 环境依赖有方向

**`h_0=0` 时（其余邻居**抵消**），取 `u=v=n`、`u'=v'=−n`：**

$$\Delta_{a\leftarrow b}(0)\ \ge\ 2F(2\beta_W)\Big|_{\beta=2.5}=2F(5)=4.551303\ \ \Longrightarrow\ \ Q=\tanh(1.137826)=0.813681$$

**大 `r=|h_0|` 渐近 `F(t)=t-\tfrac32\log t+O(1)`（我核：`t=20/50/100/200` 时差 `−0.245/−0.233/−0.230/−0.228`）⇒ `Δ(h_0)=O(r^{−1})`。**

**数值（`h_0 = r n'`，`n'⊥n`）：**

| `r` | `Δ` 近似 | `Q = tanh(Δ/4)` |
|---|---|---|
`0` | 最大 | **`≈1`** |
`0.5` | `4.3952` | `0.800067` |
`1.0` | `4.0174` | `0.763412` |
`2.0` | `3.1565` | `0.657911` |
`4.0` | `2.0347` | `0.468921` |
`8.0` | `1.1413` | **`0.277826`** |

$$\boxed{\ \text{铁磁有序环境（}r\text{ 大）}\Rightarrow Q\text{ 小（有利）};\qquad\text{阻挫/无序环境（}h_0\approx0\text{）}\Rightarrow Q\approx0.81\text{（不利）}\ }$$

**⇒ 我上一轮那个 `E[Q]=0.899` 用的是**独立均匀**的 `s,s'`，对应**无序端**。⇒ 真值取决于 `ν_slow` 有多有序，而 `β=2.5` 是交叉区。**

### 5⁗⁗⁗⁗⁗‴.3 `(9)` 的第二不等式（`[T]`）

$$\tanh\Big(\tfrac14\!\!\sum_{b\in D_a}\!\!\Delta_{ab}(Z^{(b)})\Big)\ \le\ \sum_{b\in D_a}Q_{b\to a}(Z^{(b)})\qquad(\text{用 }\tanh\text{ 的次可加性})$$

**我核了 `x=(0.5,0.5)`：`tanh(1)=0.761594 ≤ 0.924234` ✓；`(1,2)`：`0.995055 ≤ 1.725622` ✓；`(0.1,0.2,0.3)`：`0.537050 ≤ 0.588356` ✓**

**⚠️ 但 `(9)` 给出的权重是 `Q_{b→a}(Z^{(b)})`，**不是** `Q_{b→a}(Y)`（同一个 `Y`）。这两个不是一件事。**

### 5⁗⁗⁗⁗⁗‴.4 【`[U]`】唯一缺口 `(10)`

$$\boxed{\ \mathbb E_{\rm coupling}\Big[\prod_{j=1}^nQ_{e_j}(Z^{(j)})\Big]\ \le\ \mathbb E_{\nu_{\rm slow}}\Big[\prod_{j=1}^nQ_{e_j}(Y)\Big]\ }\tag{10}$$

**⚠️ `(3)–(8)` **不**给出它。相反 `(3)` 说明不同混合环境对应的快条件分布可以被**显著重加权**。**

**⇒ 从 Wilson 精确积分到 `(10)`，仍缺一条真正的**概率支配定理**。**

### 5⁗⁗⁗⁗⁗‴.5 【记账】反例的范围必须说准

**⚠️ 上一轮那个一般三变量反例**只能否定一般 Gibbs 模型的 ADP 定理**，**不能**直接作为 Wilson 模型违反 ADP 的反例。**

**⚠️ 同理，Wilson 的因子分解与 `S³` 旋转结构**也还不足以证明 ADP**。**

$$\boxed{\ \text{原始退火 ADP 在 Wilson 中仍是 }[U];\qquad\text{确定性混合环境支配已是 }[T]\ }$$

### 5⁗⁗⁗⁗⁗‴.6 审计（本轮）

| 命题 | 判 |
|---|---|
Wilson 快分量精确条件分布 `(2)` | **`[T]`** |
环境改变导致的重加权界 `(3)–(4)` | **`[T]`** |
**`Δ` 随 `|h_0|` 单调衰减（`h_0=0` 时 `Q≈0.81`，`r=8` 时 `Q≈0.28`）** | **`[T]`（数值，Bessel 结构）** |
Wilson 有效作用量的环境依赖多体性 `(5)–(8)` | **`[T]`** |
混合环境下的单步传播界 `(9)` | **`[T]`** |
确定性最坏环境 disagreement percolation `(11)–(13)` | **`[T]`** |
**原始退火 ADP `(10)`** | **`[U]`** |
**`β=2.5` 下 `Q_e(Y)` 支配混合环境传播** | **`[U]`** |
一般三变量反例的**范围** | **`[T]`（只否一般 Gibbs，不否 Wilson）** |

**⇒ `Hessian` 控制、快测度重加权界、单点 `E[Q_e]` —— 三者**都不够**得到 `(10)`。**

--- 

## §5⁗⁗⁗⁗⁗⁗ 单边退火耦合 `[T]` + 三阶累积量 `[T]` + 其三阶界是**空界**（第 276 轮）

### 5⁗⁗⁗⁗⁗⁗.1 【`[T]`】单边退火耦合引理（且对任意有限 Gibbs 成立）

**条件于 `z`，独立抽 `S,S' ~ π_z = ν(·|z)`，再对 `γ_a(·|S,z)`、`γ_a(·|S',z)` 取最大耦合 `(U,U')`。**

**⇒ `(U,S,z)` 与 `(U',S',z)` 的边缘都恰是 `ν_slow`（因 `S|z ~ ν(·|z)`、`U|S,z ~ γ_a`）。**

$$\boxed{\ \mathbb P(U\ne U')\ \le\ \mathbb E_{\nu_{\rm slow}}Q_{b\to a}(Y)\ }\tag{2}$$

**⭐ 这一步把「固定环境下的影响控制」转成了「Gibbs 平衡环境下的**平均**控制」——而**只需一次坐标扰动**。**

**⚠️ 它对任何满足有限交叉振幅的有限 Gibbs 测度都成立（不只 Wilson）；但**只处理一步**，不给多步路径传播。**

### 5⁗⁗⁗⁗⁗⁗.2 【`[T]`】Wilson 三阶累积量：环境敏感度由快分量边界力三阶矩控制

$$D_aD_b\Psi_C[h_a,h_b]=\beta_W^2\mathrm{Cov}_{\mu_{C,y}}(X_a,X_b)\tag{3}$$
$$D_cD_aD_b\Psi_C[h_c,h_a,h_b]=\beta_W^3\kappa_{\mu_{C,y}}(X_a,X_b,X_c)\tag{4}$$

**`|X_i−EX_i| ≤ 2w_C(i)|h_i|`（因 `|X_i| ≤ w_C(i)|h_i|` 且均值同range）⇒**

$$\boxed{\ \big|D_cD_aD_b\Psi_C[h_c,h_a,h_b]\big|\ \le\ 8\beta_W^3w_C(a)w_C(b)w_C(c)|h_a||h_b||h_c|\ }\tag{5}$$

**我核常数链：`8β³=125`（`β=2.5`）；`32β³=500`；`(9)` 的 `8β³=125` ✓**

### 5⁗⁗⁗⁗⁗⁗.3 【`[T]`】`Q_e` 的环境 Lipschitz 估计

$$\Lambda_{ab,c}=\!\!\sum_{C:\,a,b,c\in\partial C}\!\!w_C(a)w_C(b)w_C(c)\tag{6}$$

$$|\Delta_{a\leftarrow b}(z)-\Delta_{a\leftarrow b}(z')|\ \le\ 32\beta_W^3\!\!\sum_{c\ne a,b}\!\!\Lambda_{ab,c}\|z_c-z'_c\|\tag{8}$$
$$\boxed{\ |Q_{b\to a}(z)-Q_{b\to a}(z')|\ \le\ 8\beta_W^3\!\!\sum_{c\ne a,b}\!\!\Lambda_{ab,c}\|z_c-z'_c\|\ }\tag{9}$$

$$\beta_W=2.5:\quad 8\beta_W^3=125\tag{10}$$

### 5⁗⁗⁗⁗⁗⁗.4 【`[T]`】混合环境传播的三体修正展开

$$Q_{e_j}(Z^{(j)})\ \le\ Q_{e_j}(Y)+R_j,\qquad R_j=8\beta_W^3\!\!\sum_c\!\!\Lambda_{e_j,c}\|Z^{(j)}_c-Y_c\|$$
$$\prod_{j=1}^n\big(Q_{e_j}(Y)+R_j\big)=\!\!\sum_{D\subseteq[n]}\!\!\Big(\prod_{j\notin D}Q_{e_j}(Y)\Big)\Big(\prod_{j\in D}R_j\Big)\tag{12}$$
$$R_j\le16\beta_W^3\!\!\sum_c\!\!\Lambda_{e_j,c}\mathbf 1_{\{Z^{(j)}_c\ne Y_c\}}\ \overset{\beta=2.5}{=}\ 250\!\!\sum_c\!\!\Lambda_{e_j,c}\mathbf 1_{\{\cdot\}}\tag{14}$$

**⭐⭐ **结构**：原始退火路径不是唯一传播对象。环境变化通过**三体重叠** `Λ_{ab,c}` 产生额外传播分支**

$$\boxed{\ \text{严格的多步展开应是带三体超边修正的【路径—超图】展开，不只是普通图上的简单路径}\ }$$

### 5⁗⁗⁗⁗⁗⁗.5 【我的定量判决】`(10)` 是**空界**：`8β³ = 125`，而 `Λ` 求和带 `O(ℓ³)` 重数

**对固定 `(a,b)`，同一块 `∂C` 内的 `c` 有 `|∂C|=O(ℓ³)` 个，每个 `Λ_{ab,c} ≥ 1`，且 `‖z_c−z'_c‖ ≤ 2`：**

| `ℓ` | `Σ_cΛ‖·‖` | 界 `8β³·Σ` | 判 |
|---|---|---|---|
`1` | `2.0` | **`250`** | **空** |
`2` | `16.0` | `2000` | 空 |
`3` | `54.0` | `6750` | 空 |
`4` | `128.0` | `16000` | 空 |

$$\boxed{\ \ell=1\ \text{时界已 }=250\text{，而 }|Q-Q'|\text{ 的平凡界是 }1\ \Longrightarrow\ (10)\ \text{比平凡界还松 }250\times\ }$$

**⇒ 所以三阶累积量路线虽是 `[T]`，**定量上空**。**

**⚠️ 而它与前几轮撞的是**同一堵墙**：**

| 量 | 量级 |
|---|---|
`n_a`（Dobrushin 行和的非零邻居数） | `O(ℓ³)` |
`‖K‖` | `O(ℓ³)` |
`‖S_eff‖` | `βℓ^{10}` |
**`Σ_cΛ_{ab,c}`** | **`O(ℓ³)`** |

$$\boxed{\ \text{二阶（Hessian）与三阶（累积量）路线撞同一堵墙：块边界上的【全连接】}\ }$$

### 5⁗⁗⁗⁗⁗⁗.6 为何这还没证出原始 ADP

$$\text{原始 ADP: }\ \mathbb E_{\rm coupling}\prod_{j=1}^nQ_{e_j}(Z^{(j)})\ \le\ \mathbb E_{\nu_{\rm slow}}\prod_{j=1}^nQ_{e_j}(Y)\tag{15}$$

**但 `(12)` 表明左侧 = 静态路径项 `+` 所有 `D≠∅` 的环境修正项。要证 `(15)` 至少须证这些项不使左侧超过静态项，等价于控**

$$\mathbb E_{\rm coupling}\Big[\prod_{j\in D}\mathbf 1_{\{Z^{(j)}_{c_j}\ne Y_{c_j}\}}\prod_{j\notin D}Q_{e_j}(Y)\Big]\tag{16}$$

**⚠️ 单边引理 `(2)` 不控多步联合事件；三阶界 `(9)` 也不自动提供它。**

**⚠️ **不能**把 `E∏Q(Z^{(j)}) ≤ E∏Q(Y)` 当成未证的交换或比较用——不同时刻的混合环境**不自动**具有同一 Gibbs 构型的联合分布。**

### 5⁗⁗⁗⁗⁗⁗.7 审计（本轮）

| 命题 | 判 |
|---|---|
Wilson 单边退火耦合 `(2)` | **`[T]`** |
快分量三阶累积量公式 `(4)` | **`[T]`** |
三阶响应统一上界 `(5)` | **`[T]`**（常数我核：`8β³=125`） |
`Q_e` 的环境 Lipschitz 估计 `(8)–(10)` | **`[T]`** |
**但 `(10)` 定量上是空界（`≥250` 对平凡界 `1`）** | **`[T]`** |
混合环境传播的三体修正展开 `(11)–(14)` | **`[T]`** |
**原始 ADP 的多步联合支配 `(16)`** | **`[U]`** |
**`β=2.5` 下原始 ADP** | **`[U]`** |

**⚠️ `125`（或 `250`）**本身就不是小参数**；且与 `Λ` 的 `O(ℓ³)` 重数相乘后超出平凡界。**

**⇒ 要完成原始 ADP，须另找：具体块几何的消去 / 快分量协方差的进一步结构 / 或独立的退火联合概率估计。**

--- 

## §5⁗⁗⁗⁗⁗⁗′ 撤回共享环境基例 + 最短因子路径 `[T]` + **路径和发散**（第 277 轮）

### 5⁗⁗⁗⁗⁗⁗′.1 【勘误 32】上一轮的 `(2)` 是**共享环境**版，不是原始 ADP 的基例

$$\pi_{\rm sh}(db,db')=\int\nu_Z(dz)\nu(db\mid z)\nu(db'\mid z)\qquad\text{vs}\qquad\pi_{\rm ind}(db,db')=\nu_B(db)\nu_B(db')$$

$$\pi_{\rm sh}\ne\pi_{\rm ind}\quad\text{(共享 }z\text{ 使 }b,b'\text{ 相关)}$$

$$\boxed{\ \text{`(2)` 证的是共享环境下的单边耦合，【不是】独立边界副本下的原始 ADP 基例}\ }$$

**⇒ 我上一轮把 `(2)` 与 ADP 基例挂在一起的标记**撤回**。`(2)` 本身仍 `[T]`，但**范围**只到 `π_sh`。**

### 5⁗⁗⁗⁗⁗⁗′.2 【`[T]`】最短因子路径的几何事实

$$S_{\rm slow}(y)=\sum_{f\in\mathcal F}\Phi_f(y_{\partial f}),\qquad \Phi_C=-\log Z_C$$

**最短因子路径 `a_0−f_1−a_1−f_2−⋯−f_n−a_n`，`P={a_0,…,a_n}`：**

$$\boxed{\ |\partial f\cap P|\ \le\ 2,\qquad\text{且若为 2 则两个路径变量必【相邻】}\ }\tag{8}$$

**证明：若 `f⊇{a_i,a_j}` 且 `j≥i+2`，则 `a_i−f−a_j` 是长度 2 的路径，可替换更长的片段 ⇒ 与最短矛盾。∎**

**⭐⭐ 这是本轮最重要的结构事实。**

### 5⁗⁗⁗⁗⁗⁗′.3 【`[T]`】固定路径外环境后，路径条件测度是**一维非齐次 Markov 链**

$$\nu(dy_P\mid\xi)\propto\prod_{j=0}^ne^{-U_j(y_{a_j};\xi)}\prod_{j=1}^ne^{-V_j(y_{a_{j-1}},y_{a_j};\xi)}\prod_{j=0}^nd\sigma(y_{a_j})\tag{9}$$

**由 `(8)`，同时含 `a_{j−1},a_j` 的因子**不含其他路径变量** ⇒**

$$\boxed{\ \Delta_j(\xi)=\Delta_{a_{j-1}\leftarrow a_j}(y_{\widehat{a_{j-1},a_j}})\qquad\text{只依赖 }\xi\ \text{（路径外环境）}\ }\tag{10}$$

**⇒ **沿最短因子路径，`Q_j` 对路径内部变量完全不敏感**——先前最棘手的「路径内先前改变会改变当前 `Q`」在此**消失**。**

### 5⁗⁗⁗⁗⁗⁗′.4 【`[T]`】多步同构型乘积控制

$$T_j(dt\mid s;\xi)=\frac{e^{-V_j(s,t;\xi)}h_j(t;\xi)d\sigma(t)}{\int e^{-V_j(s,u;\xi)}h_j(u;\xi)d\sigma(u)},\qquad h_j\ \text{只依赖 }t,\xi\tag{12}$$

$$\boxed{\ \sup_{s,s'}\|T_j(\cdot\mid s;\xi)-T_j(\cdot\mid s';\xi)\|_{\rm TV}\ \le\ Q_j(\xi)\ }\tag{13}$$

**Dobrushin 系数次乘 ⇒**

$$\sup_{s,s'}\big\|\nu(Y_{a_n}\in\cdot\mid Y_{a_0}=s,\xi)-\nu(\cdot\mid Y_{a_0}=s',\xi)\big\|_{\rm TV}\ \le\ \prod_{j=1}^nQ_j(\xi)\tag{14}$$
$$\boxed{\ |\mathrm{Cov}_\nu(F,G\mid\xi)|\ \le\ 2\|F\|_\infty\|G\|_\infty\prod_{j=1}^nQ_j(\xi)\ }\tag{15}$$
$$\boxed{\ \mathbb E_\nu\big[|\mathrm{Cov}_\nu(F,G\mid Y_{P^c})|\big]\ \le\ 2\|F\|_\infty\|G\|_\infty\,\mathbb E_{\nu_{\rm slow}}\Big[\prod_{j=1}^nQ_{a_{j-1}\to a_j}(Y)\Big]\ }\tag{16}$$

**⭐⭐⭐ `(16)` 是**真正的多步同构型退火乘积控制**——不用小耦合展开，也不把不同混合环境强行替换成同一构型。**

### 5⁗⁗⁗⁗⁗⁗′.5 【我的定量判决】即便 `(16)` 成立，路径和**发散**

$$\mathrm{Cov}=\underbrace{\mathbb E[\mathrm{Cov}(F,G\mid Y_{P^c})]}_{(16)\ \text{控住}}+\underbrace{\mathrm{Cov}(\mathbb E[F\mid Y_{P^c}],\mathbb E[G\mid Y_{P^c}])}_{\text{路径外剩余，}[U]}\tag{17}$$

**先看第一项**本身**求和是否收敛。简单路径数 `~Δ^n`，路径权重 `~e^{−kn}` ⇒ 需**

$$\boxed{\ \Delta\,e^{-k}\lt1\iff k\gt\ln\Delta\ }$$

| 度数 `Δ` | 需 `k >` |
|---|---|
`6` | `1.7918` |
`8` | **`2.0794`** |
`10` | `2.3026` |

**而 `k` 有个**上界**，来自 `|h_0| ≤ 2d = 6`（至多 6 个邻居，每个 `|y|=1`）：**

| `r=|h_0|` | `Δ` | `Q` | `k` | `k\gt\ln8`? |
|---|---|---|---|---|
`0` | `4.5513` | `0.813681` | `0.2062` | FAIL |
`2` | `3.1565` | `0.657911` | `0.4187` | FAIL |
`4` | `2.0347` | `0.468921` | `0.7573` | FAIL |
**`6`**（物理上界） | `1.4675` | `0.351252` | **`1.0463`** | **FAIL** |

$$\boxed{\ k_{\max}\approx1.05\ \lt\ \ln8=2.08\ \Longrightarrow\ \text{差 }1.99\times\ \Longrightarrow\ \text{退火路径和【发散】}\ }$$

**⇒ 即使取度数 `10`，仍差 `2.20×`。**

**⇒ 要收敛需 `Δ \lt 0.50`（`Δ=8`），而 `r≤6` 时 `Δ ≥ 1.4675` ⇒ **不可能**。**

$$\boxed{\ \text{所以 }(16)\ \text{这一项【单项】就过不了路径求和——与第二项无关}\ }$$

**⚠️ 这是**新增的独立障碍**：它不是「第二项没控」，而是「第一项控住了但求和仍发散」。**

**⇒ 根因与前面几轮同源：`|h_0| ≤ 2d` 是**图度数**的硬上界 ⇒ `Δ` 有正下界 ⇒ `k` 有上界。**

### 5⁗⁗⁗⁗⁗⁗′.6 【`[U]`】第二项：路径外环境的递归展开

$$\Big|\mathrm{Cov}\big(\mathbb E[F\mid Y_{P^c}],\mathbb E[G\mid Y_{P^c}]\big)\Big|\ \le\ C\|F\|\|G\|\!\!\sum_{\substack{\Gamma:A\leftrightarrow B\\ \Gamma\ \text{经路径外因子}}}\!\!\mathbb E_\nu\prod_{e\in\Gamma}Q_e(Y)\tag{18}$$

**⚠️ 递归并不自动成立：`S_slow` 的因子图**含大量环路**（Wilson 格子不是树），递归会产生**无法控制的重复权重**。**

**⚠️ 且**即使**递归成功，§5⁗⁗⁗⁗⁗⁗′.5 已说明单项求和都发散。**

### 5⁗⁗⁗⁗⁗⁗′.7 审计（本轮）

| 命题 | 判 |
|---|---|
**`(2)` 作为原始 ADP 的基例** | **撤回（勘误 32）** |
`(2)` 在**共享环境** `π_sh` 下成立 | `[T]` |
Wilson `S_slow` 的因子图表示 `(6)` | `[T]` |
**最短因子路径上 `\|∂f∩P\|≤2` 且相邻 `(8)`** | **`[T]`** |
固定 `ξ` 后路径条件测度是一维 Markov 链 `(9)` | `[T]` |
**`Δ_j(ξ)` 只依赖路径外环境 `(10)`** | **`[T]`** |
`(13)`、`(14)`、`(15)`、`(16)` | **`[T]`** |
**`(16)` 单项路径和** | **`[F]`（`k_max≈1.05 < lnΔ≥1.79`）** |
路径外剩余协方差 `(18)` | `[U]` |
完整原始 ADP | `[U]` |
`β=2.5` 下完整退火路径和指数收敛 | `[F]`（同上根因） |

**⚠️ 两条不得混为一谈：**
1. **「第一项控住了但路径和发散」——本轮新增，已判 `[F]`；**
2. **「第二项（路径外剩余）没控」——仍是 `[U]`。**

--- 

## §5⁗⁗⁗⁗⁗⁗″ 路径消元递归：`(29)` `[T]` + `(37)` **失败**（差 `O(ℓ³)`，第 278 轮）

### 5⁗⁗⁗⁗⁗⁗″.1 精确边缘化与递归结构（`[T]`）

$$S(y_P,z)=S_O(z)+S_P(y_P,z),\quad Z_P(z)=\int e^{-S_P}d\sigma_P,\quad \Psi_P(z)=-\log Z_P(z)$$
$$\boxed{\ \nu^P(dz)=\frac1{Z^P}e^{-S_O(z)-\Psi_P(z)}\prod_{c\in O}d\sigma(z_c)\ }\tag{20}$$
$$\mathrm{Cov}_\nu(F,G)=\mathbb E_{\nu^P}[\mathrm{Cov}_\nu(F,G\mid Y_O)]+\mathrm{Cov}_{\nu^P}(\bar F,\bar G)\tag{21}$$

**⚠️ **必须**在**新** Gibbs 测度 `ν^P` 下展开剩余协方差——**不能**把原因子图删掉路径后在剩余图上重复原证明。**

### 5⁗⁗⁗⁗⁗⁗″.2 新有效因子的二阶响应（`[T]`，本轮最实质的结果）

$$D^2_{cd}\Psi_P=\mathbb E_{\mu_{P,z}}[D^2_{cd}S_P]-\mathrm{Cov}_{\mu_{P,z}}(A_c,A_d)\tag{25}$$
$$\boxed{\ |\mathrm{Cov}_{\mu_{P,z}}(A_{c,f},A_{d,g})|\ \le\ \tfrac14L_{c,f}L_{d,g}\Theta_{fg}(z)\ }\tag{27}$$

$$\Theta_{fg}(z)=\begin{cases}\prod_{k=j_f+1}^{i_g}Q_k(z),&j_f\lt i_g\\ \prod_{k=j_g+1}^{i_f}Q_k(z),&j_g\lt i_f\\ 1,&I_f\cap I_g\ne\varnothing\end{cases}\tag{28}$$

$$\boxed{\ \|D^2_{cd}\Psi_P(z)\|_{\rm op}\le B^P_{cd}+\tfrac14\!\!\sum_{f\ni c}\sum_{g\ni d}\!\!L_{c,f}L_{d,g}\Theta_{fg}(z)\ }\tag{29}$$

**⭐⭐ 所以：**路径外剩余协方差不是无关的新问题，而是由路径消元诱导的有效 Gibbs 模型**；间接响应已带**原路径上的 `Q` 乘积**。**

### 5⁗⁗⁗⁗⁗⁗″.3 【我的定量判决】(23) 救了 `L`，但 `(37)` 仍失败 `O(ℓ³)`

**① `(23)` 的真实收益：**`|\partial C\cap I_P|\le2` ⇒**

$$L_{c,C}\le\pi\beta_W^2w_C(c)\!\!\sum_{a\in\partial C\cap I_P}\!\!w_C(a)\ \le\ \pi\beta_W^2w_C(c)\cdot O(1)\ =\ O(\pi\beta_W^2)=39.3$$

$$\boxed{\ L\ \text{是纯数、}O(1)\text{、与 }\ell\text{ 无关（而非 }O(\ell^3)\text{）}\ }$$

**② 但 `B^P` **没有**路径约束：`c,d` 只要**同在某个块的边界**即可，`c` 不必在路径上。**

$$w_C(c)=w_C(d)=1,\ m=0\ \Longrightarrow\ B^P_{cd}\ \ge\ \beta_W^2=6.25$$

$$\Delta^P\ \ge\ \pi^2\beta_W^2=61.69\ \Longrightarrow\ \tanh(\Delta^P/4)=\tanh(15.42)=1.000000000000$$

**⇒ **新边权 `≈1`**。**

**③ 而每个 `c` 有多少这样的 `d`？共享一个块边界 `∂C` 的所有 `d`：`|∂C|=O(ℓ³)` 个：**

| `ℓ` | `\|∂C\|` | `ρ(𝒞̃^P)` 下界 | `<1`? |
|---|---|---|---|
`1` | `1` | `1` | **否** |
`2` | `8` | `8` | **否** |
`3` | `27` | `27` | **否** |
`4` | `64` | `64` | **否** |

$$\boxed{\ \rho(\widetilde{\mathcal C}^P)\ \ge\ O(\ell^3)\ \gg\ 1\ \Longrightarrow\ \text{(37)}\ \text{失败}\ }$$

**④ 根因（本轮说清）：**

$$\text{单个 }-\log Z_C\ \text{把【整个边界 }\partial C\text{】全连接}$$
$$\text{路径消元只限制【路径变量】在 }\partial C\text{ 中的个数（}\le2\text{）}$$
$$\text{而 }B^P\text{ 里的 }c,d\ \text{都是【路径外】变量，\textbf{不受此限制}}$$

$$\boxed{\ \text{同一堵墙：块边界上的全连接 —— 路径消元把它【搬家】而没把它【变小】}\ }$$

### 5⁗⁗⁗⁗⁗⁗″.4 逐轮量级对照：全部同源

| 出现处 | 量级 |
|---|---|
`n_a`（Dobrushin 行和非零邻居数） | `O(ℓ³)` |
`‖K‖` | `O(ℓ³)` |
`Σ_cΛ_{ab,c}` | `O(ℓ³)` |
**`ρ(𝒞̃^P)`（本轮）** | **`O(ℓ³)`** |
`k_max` vs `lnΔ`（路径和） | 差 `1.99×`（**图度数**上界，另一条锁） |

**⇒ 六处、四种不同推导、同一个 `O(ℓ³)`。**

### 5⁗⁗⁗⁗⁗⁗″.5 【`[U]`】两个未闭合条件

**① `(37)` `ρ(𝒞̃^P)<1`：**本轮判 `[F]`**（差 `O(ℓ³)`）。**

**② `(40)` 消元稳定性（退火权重在边缘化下保持支配）：**

$$\mathbb E_{\nu^P}\Big[\prod_{j=1}^m\widetilde Q^P_{c_{j-1}\to c_j}(Z)\Big]\ \le\ \sum_\Gamma K(\Gamma)\,\mathbb E_{\nu_{\rm slow}}\Big[\prod_{e\in\Gamma}Q_e(Y)\Big]$$

**⚠️ 注意 `(29)` 的 `Θ_{fg}(z)` 已是**同一个路径外构型** `z` 上的 `Q` 乘积（这一点已解决）；但 `(36)` 用了**环境上确界** ⇒ 一旦代入 Dobrushin 矩阵就**损失退火平均优势**。**

**⚠️ `(40)` 未证，**不能**当成递归展开的已知性质。**

### 5⁗⁗⁗⁗⁗⁗″.6 审计（本轮）

| 命题 | 判 |
|---|---|
路径消元后的精确边缘 Gibbs 测度 `(19)–(21)` | **`[T]`** |
最短因子路径的因子几何 `(23)` | **`[T]`** |
新有效因子的混合二阶响应 `(25)` | **`[T]`** |
沿原路径的 Markov 响应 `(27)–(29)` | **`[T]`** |
**`(23)` 使 `L=O(1)`（与 `ℓ` 无关）** | **`[T]`（我算）** |
Wilson 局部常数显式化 `(30)–(35)` | **`[T]`** |
**`(37)` 剩余矩阵收缩** | **`[F]`（本轮：`ρ≥O(ℓ³)`，根因是 `B^P` 无路径约束）** |
**`(40)` 消元稳定性** | **`[U]`** |
`β=2.5` 下完整原始 ADP | **`[U]`** |

--- 

## §5⁗⁗⁗⁗⁗⁗‴ 端点响应权重补齐 `[T]` + 游走展开 `[T]` + `(27)` 失败、且 `(26)` **弄坏界**（第 279 轮）

### 5⁗⁗⁗⁗⁗⁗‴.1 【`[T]`】补齐上一轮 `(39)` 的两个端点权重

$$W_0(u,z)=-\log\!\!\int_{(S^3)^n}\!\!e^{-S_P(u,y_{a_1\cdots n},z)}\prod_{j=1}^nd\sigma(y_{a_j})\tag{9}$$
$$\lambda_0(du\mid z)=\frac{e^{-W_0(u,z)}d\sigma(u)}{\int e^{-W_0(v,z)}d\sigma(v)}\tag{10}$$
$$\Omega^L_c(P)=\sup_{\substack{z,z':\,z_{O\setminus c}=z'_{O\setminus c}}}\sup_{u,u'}\big|W_0(u,z)-W_0(u,z')-W_0(u',z)+W_0(u',z')\big|\tag{11}$$
$$\boxed{\ \eta^L_c(P)=\tanh(\Omega^L_c(P)/4),\qquad \delta_c(\bar F)\le2\|F\|_\infty\eta^L_c(P)\ }\tag{12),(14}$$

**⇒ 且 `c ∉ H(P) ⇒ η^L_c=0`（`W_0` 不依赖 `z_c`）⇒ 端点向量只在**路径边界**上有支撑 ✓**

### 5⁗⁗⁗⁗⁗⁗‴.2 【`[T]`】剩余图游走展开（含**环路**）

$$D^P=(I-C^P)^{-1}=\sum_{m\ge0}(C^P)^m\tag{6}$$
$$D^P_{cd}=\mathbf1_{\{c=d\}}+\sum_{m\ge1}\sum_{\substack{c_0=c,\,c_m=d\\ c_1\cdots c_{m-1}\in O}}\prod_{j=1}^mC^P_{c_{j-1}c_j}\tag{7}$$

**⭐ `(7)` **允许重复顶点** ⇒ 它是所有有限**游走**的总权重 ⇒ **环路也被计入**（这正是处理非树因子图所需的）。**

$$|\mathrm{Cov}_{\nu^P}(U,V)|\le\tfrac14\sum_{c,d}\delta_c(U)D^P_{cd}\delta_d(V)\tag{8}$$

**⇒ 完整的有限体积条件定理：**

$$\boxed{\ \begin{aligned}|\mathrm{Cov}_\nu(F,G)|\le{}&2\|F\|\|G\|\,\mathbb E_\nu\Big[\prod_{j=1}^nQ_{a_{j-1}\to a_j}(Y)\Big]\\ &+\|F\|\|G\|\sum_{m\ge0}\sum_{c_0\cdots c_m}\eta^L_{c_0}(P)\Big[\prod_{j=1}^mC^P_{c_{j-1}c_j}\Big]\eta^R_{c_m}(P)\end{aligned}\ }\tag{19}$$

**成立条件 `ρ(C^P)<1`。**⚠️ **这是**条件定理**，不是 `β=2.5` 的无条件证明。**

### 5⁗⁗⁗⁗⁗⁗‴.3 【`[T]`】Wilson 因子还原剩余边权

$$\|D^2_{cd}\Psi_P\|\le B^P_{cd}+\tfrac14\!\!\sum_{f\ni c}\sum_{g\ni d}\!\!L_{c,f}L_{d,g}\Theta_{fg}(z)\tag{22}$$
$$\Delta^P_{c\leftarrow d}\le\Delta^O_{c\leftarrow d}+\pi^2\Big[B^P_{cd}+\tfrac14\!\!\sum_{f,g}\!\!L_{c,f}L_{d,g}\sup_z\Theta_{fg}\Big]\tag{24}$$
$$\boxed{\ L_{c,C}\le\pi\beta_W^2w_C(c)\!\!\sum_{a\in\partial C\cap I_P}\!\!w_C(a)\ }\tag{25}$$
$$\widehat C^P_{cd}=\tanh\Big(\frac{\Delta^O_{c\leftarrow d}+\pi^2\sup_zM^P_{cd}(z)}4\Big)\tag{26}$$

### 5⁗⁗⁗⁗⁗⁗‴.4 【我的定量判决】`(27)` 失败，而 `(26)` 的显式分解**把界弄坏**

**① `L` 的量级（`(23)` 收益后仍很大）：**

$$L_{c,C}\ \le\ \pi\beta_W^2w(c)\cdot2=39.27\qquad(\beta=2.5)$$

**② `(24)` 两项各自都使新边权 `=1`：**

| 项 | 值 | `tanh(π²·/4)` |
|---|---|---|
`B^P ≥ β² = 6.25` | `π²B^P = 61.7` | **`1.0000000000`** |
`¼L² = 385.5` | `π²·381\ldots = 3805` | **`1.0`** |

$$\boxed{\ \tfrac14L^2=385.5\ \text{比}\ B^P=6.25\ \text{大}\ 62\times\ \Longrightarrow\ (26)\ \text{的显式分解把界【弄坏】而不是弄好}\ }$$

**③ `ρ` 下界：**

| `ℓ` | `\|∂C\|` | `ρ(C^P) ≥` |
|---|---|---|
`1` | `1` | `1` |
`2` | `8` | `8` |
`4` | `64` | `64` |

$$\boxed{\ \text{(27)}\ \text{失败，差 }O(\ell^3)\ }$$

**④ 即便保住退火 `Θ`（不放 `sup`），`L` 项仍太大：**

$$\tfrac14L^2\Theta\le0.5\ \text{需}\ L\le1.41,\quad \text{实际}\ L\le39.3\ \Longrightarrow\ \text{差 }27.8\times$$
$$\Theta\le e^{-km}\big|_{k\approx1.05,m\ge1}=0.350\ \Longrightarrow\ \tfrac14L^2\Theta=135\ \gg\ 0.5$$

**⇒ **退火也救不了 `L` 项。****

### 5⁗⁗⁗⁗⁗⁗‴.5 【本轮最有用的横向对照】三条独立路线、三个 `β` 阈值

| 路线 | `β` 阈值 | `2.5/阈值` |
|---|---|---|
`tanh(π²K/4)` 粗界 | `≲0.05` | **`50×`** |
**`L` 项（本轮）** | **`≲0.474`** | **`5.3×`** |
`β_pin`（pinning 障碍） | `1.01092` | `2.47×` |

$$\boxed{\ \text{三条独立路线给出三个阈值，全部}\lt2.5\ }$$

**⚠️ 注意 `L ∝ β²` ⇒ 要 `L≲1.41` 需 `β≲0.474`。**

**⇒ 所以 `(25)` 这个「显式化」虽然让**结构**变清楚了，但**量级**上把最坏的一项（`¼L²`）放大了 `62×`。**

**⇒ 有用的推论：**下一步若要救 `(27)`，必须把 `L` 的**平方**换掉——例如不用 Cauchy–Schwarz 把两个因子得分拆成「振幅乘积」，而直接估它们沿路径的**联合**响应。**

### 5⁗⁗⁗⁗⁗⁗‴.6 审计（本轮）

| 命题 | 判 |
|---|---|
剩余 Gibbs 的精确边缘化 `(1)` | **`[T]`** |
剩余单点影响系数由交叉振幅控制 `(3)–(4)` | **`[T]`** |
Dobrushin 矩阵的**游走**展开 `(6)–(8)`（含环路） | **`[T]`** |
**左右端点响应权重 `(9)–(15)`（补齐 `(39)`）** | **`[T]`** |
剩余图路径和控制剩余协方差 `(17)`（条件 `ρ<1`） | **`[T]`** |
Wilson 因子结构还原剩余边权 `(20)–(26)` | **`[T]`** |
**`(27)` 在 `β=2.5`** | **`[F]`（差 `O(ℓ³)`；`¼L²` 比 `B^P` 大 `62×`）** |
`(24)` 用 `sup_zΘ` 重新引入最坏环境控制 | **`[T]`（记账）** |
**完整原始退火 ADP** | **`[U]`** |

--- 

## §5⁗⁗⁗⁗⁗⁗⁗ 动态退火时空路径：摆脱 Dobrushin 谱半径（`[T]`）+ `(DAPT)` 是唯一缺口（第 280 轮）

### 5⁗⁗⁗⁗⁗⁗⁗.1 【`[T]`】相关边标记：摆脱 Dobrushin 的方式

**关键改变：不再要求不同传播步的边权独立，也不要求它们来自同一静态 Gibbs 构型。**

**逐坐标粘合：`z^(0)` 取 `x'` 全部邻居坐标，逐步把 `b_j` 从 `x'` 换成 `x`，`μ_j=γ_a(·|z^(j))`：**

$$\|\mu_{j-1}-\mu_j\|_{\rm TV}\ \le\ Q_{b_j\to a}(z^{(j-1)})\tag{3}$$

**对每对 `(μ_{j-1},μ_j)` 取最大耦合 `π_j`，由**粘合引理**得联合 `(U_0,…,U_k)`，每对相邻边缘正确：**

$$\mathbb P(U_{j-1}\ne U_j\mid x,x')\le Q_{b_j\to a}(z^{(j-1)})\tag{4}$$
$$\boxed{\ X_a^{\rm new}\ne X_a'^{\rm new}\ \Longrightarrow\ \sum_{b\in D_a}M_{b\to a}\ \ge\ 1\ }\tag{5}$$

**⭐⭐ **这些标记不需要相互独立**——它们来自同一个粘合耦合，允许强相关。**这正是摆脱 Dobrushin 行和条件的关键。****

### 5⁗⁗⁗⁗⁗⁗⁗.2 【`[T]`】时空不一致路径与动态退火权重

$$\boxed{\ X_T(a)\ne X'_T(a)\ \Longrightarrow\ \exists\ \text{从 }D_B(\eta,\eta')\times\{0\}\text{ 到 }(a,T)\text{ 的打开时空路径}\ }\tag{7}$$

$$\mathfrak W_T(\pi)=\mathbb E_{\mathbb P_T^{\rm ann}}\Big[\mathbf1_{\{\eta_{b_0}\ne\eta'_{b_0}\}}\!\!\prod_{(b,t)\to(a,t+1)\in\pi}\!\!M_{b\to a,t+1}\Big]\tag{8}$$

**⇒ 期望同时包含：两份**独立**边界、两份热浴的联合演化、每步的粘合随机性、动态混合环境的边标记。**

$$\boxed{\ |\mathrm{Cov}_\nu(F_A,G_B)|\ \le\ 2\|F_A\|_\infty\|G_B\|_\infty\liminf_{T\to\infty}\!\!\sum_{\pi\in\mathcal P_T(A,B)}\!\!\mathfrak W_T(\pi)\ }\tag{14}$$

$$\boxed{\ \textbf{定理 (14) 不需要任何 Dobrushin 谱半径小于 1 的条件}\ }$$

**⇒ 与静态路线 `E_{ν_{\rm slow}}[ΠQ_e(Y)]` 的本质差别：动态版直接对**传播过程中实际出现的**混合环境、条件耦合与边标记联合平均。**

### 5⁗⁗⁗⁗⁗⁗⁗.3 【我的定量补充】`sup_T` 不是免费的：时空路径计数

**一条 `m` 步空间路径在时间 `T` 内的时刻组合数 `=C(T,m)`，空间分支 `~Δ^m`：**

$$\#\mathcal P_T\ \sim\ \binom{T}{m}\Delta^m$$

| `L` | `T=t_mix~L⁴` | `m=4` 的 `log10#` | `m=8` 的 `log10#` |
|---|---|---|---|
`6` | `1296` | `14.7` | `27.5` |
`10` | `10⁴` | `18.2` | `34.6` |
`16` | `65536` | `21.5` | `41.2` |

$$\boxed{\ \text{路径数随 }T\text{ 与 }L\text{ 迅速爆炸}\ \Longrightarrow\ (DAPT)\ \text{要求每条的权重相应变小}\ }$$

**而权重里确有**时间衰减**：标记打开需要两份构型在 `t_{j-1}` 时**仍不同**，而两份过程在边界不同的前提下内部会部分耦合 ⇒**

$$\mathfrak W_T(\pi)\ \sim\ \prod_j\big[Q\cdot\mathbb P(\text{differ at }t_{j-1})\big]$$

$$\boxed{\ \text{(DAPT) 的全部内容 = 空间分支 }(\Delta^m)\text{ 与时间耦合衰减的竞争}\ }$$

### 5⁗⁗⁗⁗⁗⁗⁗.4 一处精度修正：`(14)` 的 `liminf` 实际等于 `sup_T`

**⚠️ `η≠η'` 时边界**永久不同** ⇒ 两份过程**永不**完全耦合 ⇒ `P_T(disagree)` **不**随 `T→0`（所以 `(14)` 无矛盾）。**

**⚠️ 但 `(14)` 的 `liminf_T` 作用于一个**随 `T` 增长**的和，实际等于 `sup_T`：**

$$\boxed{\ \text{(DAPT) 的 }\sup_T\ \text{与 }(14)\ \text{的 }\liminf\ \text{是同一件事，不是两件事}\ }$$

### 5⁗⁗⁗⁗⁗⁗⁗.5 `[U]` 唯一缺口 (DAPT)

$$\boxed{\ \sup_T\!\!\sum_{\pi\in\mathcal P_T(A,B)}\!\!\mathfrak W_T(\pi)\ \le\ Ce^{-md(A,B)}\ }\tag{DAPT}$$

**⇒ 相比原 Dobrushin 条件更适合目标：**

| 要求 | Dobrushin 路线 | 动态路线 |
|---|---|---|
每个局部系数小 | **需要** | **不需要** |
静态影响矩阵谱半径 `<1` | **需要** | **不需要** |
保留动态环境真实相关 | 丢掉了 | **保留** |
允许单边权接近 1 | 不行 | **可以**（只要长距离时空连通总概率衰减） |

**⚠️ `(14)` 已证；**`(DAPT)` 在 `β=2.5` 下未证**。**

### 5⁗⁗⁗⁗⁗⁗⁗.6 审计（本轮）

| 命题 | 判 |
|---|---|
多邻居逐坐标粘合耦合 `(3)–(6)` | **`[T]`** |
**不一致必有一条标记入边 `(5)`** | **`[T]`** |
**标记无需独立（来自同一粘合耦合）** | **`[T]`** |
不一致沿时空标记路径传播 `(7)` | **`[T]`** |
动态退火路径概率控制协方差 `(14)` | **`[T]`** |
**不需要 Dobrushin 收缩** | **`[T]`** |
**`sup_T` 的实际内容 = 空间分支 vs 时间耦合** | **`[T]`（我的定量补充）** |
**`(DAPT)` 在 `β=2.5`** | **`[U]`** |
`β=2.5` 下指数聚簇 | **`[U]`** |

**⇒ 本轮建立的是一条**新而严格**的路线：**

$$\text{Wilson }S_{\rm slow}\to\text{局部条件核}\to\text{相关边标记}\to\text{时空不一致路径}\to\text{协方差界}$$

**⇒ 它没有把静态 ADP 偷渡回来，也没有假设边标记独立。**

--- 

## §5⁗⁗⁗⁗⁗⁗⁗′ 反向祖先路径 `[T]` + `θ_r<1` 等价于 `β<0.126`（第 281 轮）

### 5⁗⁗⁗⁗⁗⁗⁗′.1 【`[T]`】反向祖先路径：**有限体积确定性因果**，不需任何混合假设

$$D_t(a)=1\ \Longrightarrow\ \sum_{b\in N(a)}M_{b\to a,t}\ge1\tag{1}\qquad M_{b\to a,t}=1\ \Longrightarrow\ D_{t^-}(b)=1\tag{2}$$

**从 `a_0∈A` 向过去追踪，`T>s_0>s_1>⋯≥0`：**

$$\boxed{\ D_T(a)=1\ \Longrightarrow\ \exists\ \text{从 }B\text{ 到 }a\text{ 的反向祖先路径}\ }\tag{4}$$

**⚠️ 这一步处理了「固定边界在整个演化期间**持续**存在」这个易被忽略的问题——边界不是「初始时刻一次性的不一致粒子」，而是**永久的祖先源**。**

$$\mathbb P_T^{\rm ann}(X_T|_A\ne X'_T|_A)\ \le\ \sum_{a_0\in A}\sum_{n\ge N}\sum_{\gamma:|\gamma|=n\atop \gamma_0=a_0}W_T(\gamma),\quad N=\Big\lceil\frac{d(A,B)}{R_0}\Big\rceil\tag{7}$$

### 5⁗⁗⁗⁗⁗⁗⁗′.2 【`[T]`】有限步反向繁殖系数推出 DAPT

$$\Theta_r(T,\pi)=\!\!\sum_{\gamma\in\mathcal W_r(a_k)}\!\!\mathbb P_T^{\rm ann}(E_T(\pi\gamma)\mid E_T(\pi)),\qquad \theta_r=\sup_{T,\Lambda,\pi}\Theta_r(T,\pi)\tag{9),(10)}$$

$$S_{n+r}(T)\le\theta_rS_n(T),\qquad S_{n+s}(T)\le\Delta^sS_n(T)\tag{15),(16)}$$

$$\boxed{\ S_{qr+s}\le|A|\Delta^s\theta_r^q\ \Longrightarrow\ m=\frac{-\log\theta_r}{rR_0},\quad C=\frac{|A|}{1-\theta_r}\Big(\sum_{s=0}^{r-1}\Delta^s\Big)\theta_r^{-1}\ }\tag{19}$$

$$\boxed{\ |\mathrm{Cov}_\nu(F_A,G_B)|\ \le\ \frac{2|A|A_r}{1-\theta_r}\theta_r^{-1}\|F_A\|_\infty\|G_B\|_\infty e^{-md(A,B)}\ }\tag{20}$$

**⭐ 条件只有 `θ_r<1`：**它允许单步影响 `>1/Δ`，也允许不同传播边高度相关。只要若干步后平均祖先延伸数降到 1 以下。****

### 5⁗⁗⁗⁗⁗⁗⁗′.3 【`[T]`】一个必须说准的细节：`(23)` 是**后验**条件

$$\sup_{\Lambda,T,\pi}\sum_{\gamma\in\mathcal W_r(\pi_n)}\frac{\mathbb E[\mathbf1_{E_T(\pi\gamma)}]}{\mathbb E[\mathbf1_{E_T(\pi)}]}\lt1\tag{23}$$

**⚠️ 这是**后验**条件影响——**不能**用正向条件核的 `(22)` 直接代。条件事件本身对过去配置产生**选择偏差**（用户已指出）。**

### 5⁗⁗⁗⁗⁗⁗⁗′.4 【我的定量判决】`θ_r<1` 等价于 `β ≲ 0.126`，且偏差方向**不利**

**① 单条慢—慢 Wilson 因子：`Φ_ab=−βy_a^TRy_b` ⇒ `Δ=4β` ⇒ `Q=tanh(β)`：**

$$\beta=2.5:\qquad Q=\tanh(2.5)=0.986614298$$

**② `θ_1 ≤ Σ_{b∈N(a)}E[Q\mid\text{历史}] = Δ\,q_{\rm eff}$，而 `Δ=8`：**

$$\Delta\tanh(\beta)\Big|_{\beta=2.5}=8\times0.986614=7.8929\ \gg\ 1$$

**③ **偏差方向不利**：条件 `E_T(π)` 说「差异沿 `π` 传过来了」⇒ 祖先**更可能**不同 ⇒ `q_eff ≥ q`：**

$$\theta_r\ \sim\ (\Delta q_{\rm eff})^r,\qquad \Delta q_{\rm eff}\ \ge\ 7.8929$$

| `r` | `(Δq)^r` 下界 |
|---|---|
`1` | `7.9` |
`2` | `62.3` |
`4` | `3881` |
`8` | `1.5\times10^7` |

$$\boxed{\ \theta_r\ \text{随 }r\ \textbf{爆炸}，不是衰减\ \Longrightarrow\ (24)\ \text{对任何 }r\ \text{都失败}\ }$$

**④ 要 `θ_1<1`：**

$$\Delta\tanh(\beta)\lt1\iff\tanh(\beta)\lt\tfrac18\iff\beta\lt\operatorname{artanh}(0.125)=0.12566$$

$$\frac{2.5}{0.12566}=19.9\times$$

### 5⁗⁗⁗⁗⁗⁗⁗′.5 【本轮最重要的横向对照】四条独立路线、四个 `β` 阈值

| 路线 | `β` 阈值 | `2.5/阈值` |
|---|---|---|
`tanh(π²K/4)` 粗界（第 270 轮） | `0.0500` | **`50×`** |
**反向祖先 `Δtanh(β)<1`（本轮）** | **`0.1257`** | **`19.9×`** |
`L` 项（第 279 轮） | `0.4737` | `5.3×` |
`β_pin`（pinning 障碍） | `1.01092` | `2.5×` |

$$\boxed{\ \text{四条独立路线、四个阈值，跨度从 }0.05\text{ 到 }1.01\ \text{——全部}\lt2.5\ }$$

**⭐ 这四条路径**互不相同**（Worst-case `tanh`、倒数二阶导、反向祖先繁殖、pinning 子模型），却都给出 `β` 阈值 `≤1.01`。⇒ **这不是某个常数不够细，而是四条独立机制都指向同一结论。****

### 5⁗⁗⁗⁗⁗⁗⁗′.6 审计（本轮）

| 命题 | 判 |
|---|---|
不一致必有反向祖先路径 `(1)–(7)` | **`[T]`** |
有限步祖先繁殖系数定义 `(9)–(10)` | **`[T]`** |
`θ_r<1` 推出 DAPT `(15)–(19)` | **`[T]`** |
DAPT 推出指数协方差衰减 `(20)` | **`[T]`** |
**`(23)` 是后验条件（不能用正向核代）** | **`[T]`（记账）** |
**`θ_r<1` 等价于 `Δtanh(β)<1` ⇒ `β<0.126`** | **`[T]`（我算；且偏差方向不利）** |
**Wilson 后验祖先估计 `(23)` 在 `β=2.5`** | **`[F]`（`θ_r` 随 `r` 爆炸）** |
`β=2.5` 下存在 `r,θ_r<1` | **`[F]`** |

**⇒ 所以本轮虽然把 DAPT 归结为一个**漂亮的有限步条件**（比单步 Dobrushin 弱、不比它强），但这个条件在 `β=2.5` 处**定量上不成立**——而原因是 `Δ\tanh(β) = 7.89`，与前面几轮同源。**

--- 

## §5⁗⁗⁗⁗⁗⁗⁗″ 【勘误 33】反向 `θ_r<1` **结构性不可能** → 改前向退火繁殖（第 282 轮）

### 5⁗⁗⁗⁗⁗⁗⁗″.1 【`[F]`】**结构性**：`θ_r ≥ 1` 对**任意** `β` 成立

**给定 `E_T(π)`，`a_k` 在 `s_{k-1}` 时不一致。它**不能无来源地**变成不一致 ⇒ 必有一个更早的不一致来源。**

**而 `d(a_k,B) > rR_0` 时前 `r` 步到不了边界；时间够早也到不了 `t=0`：**

$$\bigcup_{\gamma\in\mathcal W_r(a_k)}E_T(\pi\gamma)\ \supseteq\ E_T(\pi)\quad\Longrightarrow\quad \mathbb P_T^{\rm ann}\Big(\bigcup_\gamma E_T(\pi\gamma)\mid E_T(\pi)\Big)=1$$
$$\boxed{\ \Theta_r(T,\pi)\ \ge\ 1\qquad\Longrightarrow\qquad \theta_r\ \ge\ 1\ \text{对每个有限 }r\ }\tag{2),(3)$$

**⇒ **与 `β` 无关**——只要有非平凡传播，反向祖先就不能在内部凭空终止。**

### 5⁗⁗⁗⁗⁗⁗⁗″.2 【勘误 33】我上一轮的 `[F]` **结论对、理由错**

**我说 `θ_r<1` 失败是因为 `Δtanh(β)=7.89>1`（`β=2.5` 的定量障碍）。**

$$\boxed{\ \text{真因是【因果结构】：}\theta_r\ge1\ \text{对一切 }\beta\text{ 成立，与 }2.5\text{ 无关}\ }$$

**⇒ 我那条 `Δtanh(β)<1 ⟺ β<0.126` 算的是**另一个（前向）量**的条件，不是 `θ_r` 的条件。勘误。**

**⚠️ 教训：反向追踪的是「已发生的不一致 → 它过去必定有来源」⇒ 内部必然找到祖先。**不能把祖先平均成灭绝过程。****

### 5⁗⁗⁗⁗⁗⁗⁗″.3 【`[T]`】前向路线：`d_t ≤ A_t d_{t^-}`

$$A_t(i,j)=\begin{cases}M_{j\to a,t},&i=a\ (\text{时刻 }t\text{ 更新 }a)\\ \mathbf 1_{\{i=j\}},&\text{否则}\end{cases}\tag{6}$$
$$\boxed{\ d_t\ \le\ A_t\,d_{t^-}\qquad\Longrightarrow\qquad d_T\ \le\ A_TA_{T-1}\cdots A_1\,d_0\ }\tag{7),(8)}$$

**⭐⭐ **与反向不同，这个数量**可能为零**：所有传播路径都可能在前向更新时因成功耦合成同一值而**灭绝**。**

$$\mathscr K_r\big((\omega,b),d\omega'\times\{a\}\big)=\mathbb E_\omega\big[(\mathcal A_r)_{ab}\mathbf1_{\{\omega_{t+r}\in d\omega'\}}\big]\tag{9}$$
$$\boxed{\ \text{(FRC)}\quad\exists r\ge1,\ \theta\lt1:\ \|\lambda\mathscr K_r\|\le\theta\|\lambda\|\quad\forall\lambda\in\mathfrak M_{\rm ann}\ }\tag{FRC}$$

**⇒ 这个条件**只针对 Gibbs 边界副本生成的可达状态锥**，不对所有人为钉扎状态取最坏上确界 ⇒ **有机会避开 pinning 障碍。****

$$\mathbb P(X|_A\ne X'|_A)\ \le\ \frac{C|\partial B|}{1-\theta}\exp\Big[-\frac{-\log\theta}{rR_0q}R\Big]\tag{12}$$

### 5⁗⁗⁗⁗⁗⁗⁗″.4 【我的定量补充】前向乘数的一次扫描界

**一次更新后，`E[Σ_aD_t(a)] ≤ Σ_b(Σ_aQ_{b→a})D_{t^-}(b)` ⇒ 每次扫描的乘数 `≤ sup_bΣ_aQ_{b→a} ≤ ΔQ_max`。**

**而 `Q_max` 在**最阻挫**环境（`h_0≈0`）取到：`Q_max = tanh(F(2β)/2)`：**

| `β` | `F(2β)` | `Q_max` | `ΔQ_max = 8Q_max` |
|---|---|---|---|
`0.50` | `0.122499` | `0.061173` | **`0.489`** ✓ |
`1.00` | `0.464134` | `0.227989` | `1.824` |
`1.20` | `0.649771` | `0.313918` | `2.511` |
**`2.50`** | `2.275651` | `0.813681` | **`6.509`** |

$$\boxed{\ 8\tanh(F(2\beta)/2)\lt1\iff F(2\beta)\lt2\operatorname{artanh}(1/8)=0.251314\iff\beta\lt0.723654\ }$$

$$\frac{2.5}{0.723654}=3.455\times$$

**⚠️ 而**最佳**环境（`|h_0|=6`，物理上界）给 `Q=0.3513` ⇒ `8Q=2.81>1` ⇒ **即使最有序也过不了**。**

### 5⁗⁗⁗⁗⁗⁗⁗″.5 【本轮最有价值的横向表】五条路线、五个阈值

| 路线 | `β` 阈值 | `2.5/阈值` |
|---|---|---|
`tanh(π²K/4)` 粗界 | `0.0500` | `50.0×` |
`L` 项（倒数二阶导） | `0.4737` | `5.28×` |
**前向退火繁殖（本轮）** | **`0.7237`** | **`3.46×`** |
`β_pin`（pinning 子模型） | `1.01092` | `2.47×` |
反向祖先 | **结构性不可能**（与 `β` 无关） | — |

**⇒ **前向路线是目前**最好**的一条**（除 pinning 外）：阈值 `0.724`，差距 `3.46×`，而不是 `50×`。**

**⇒ 而它与 pinning 的 `1.011` 是同一量级——**两条最紧的路线都在 `β ~ 1`**。**

### 5⁗⁗⁗⁗⁗⁗⁗″.6 审计（本轮）

| 命题 | 判 |
|---|---|
**反向祖先 `θ_r<1`** | **`[F]`（**结构性**：`θ_r≥1` 对一切 `β`）** |
我上一轮给的 `Δtanh(β)<1` 理由 | **`[F]`（勘误 33：算错对象）** |
前向随机影响矩阵不等式 `(7)–(8)` | **`[T]`** |
前向退火繁殖算子定义 `(9)–(10)` | **`[T]`** |
`(FRC)` 推出指数传播尾 `(12)` | **`[T]`** |
**前向一次扫描乘数 `≤ΔQ_max`，阈值 `β<0.724`** | **`[T]`（我算）** |
`β=2.5` 下 `(FRC)` | **`[U]`** |
`β=2.5` 下指数聚簇 | **`[U]`** |

**⚠️ 现有的 `Var(q)`、`χ_q`、`‖K‖`、`Q_{b→a}(y)` **都不足以**给出 `(13)`——它们没控制耦合状态在发生不一致后对后代分布的选择偏差。**

**⇒ 需要一个 Wilson 专属估计，控 `E_ann[Σ_a(𝒜_r)_{ab}]` 在**边界可达的退火状态锥**上随 `r` 的衰减，且必须含**标记的联合相关性**。**

--- 

## §5⁗⁗⁗⁗⁗⁗⁗‴ 回到主线：好快块条件控制（第一阶段）（第 283 轮）

> **⚠️ 记账：此前多轮（§5⁗–§5⁗⁗⁗⁗⁗⁗⁗″）全在**慢变量聚簇**一侧，与本模块（**快变量条件 Poincaré**）是两台不同的机器。本节把主线拉回。**

### 5⁗⁗⁗⁗⁗⁗⁗‴.1 【`[T]`】任意边界下的有限块 Poincaré

$$H_C(x;\eta,\zeta)=\beta_W\!\!\sum_{e\in E_C^\star}\!\!s_e(x,\eta)+H_C^{(0)}(\zeta),\qquad |E_C^\star|\le\Delta n_C$$

$$\boxed{\ \mathrm{osc}\,H_C\ \le\ 2\beta_W\Delta n_C\ }$$

**`σ_C = ⊗dσ(x_v)` 的谱隙 `= 3`（单 `S³` 首个非零 `−Δ` 本征值 `=3`；乘积取最小）。Holley–Stroock：`C_P(μ) ≤ e^{osc H}C_P(σ)`，`C_P(σ_C)=1/3`：**

$$\boxed{\ \mathrm{Var}_{\mu_C^{\eta,\zeta}}(f)\ \le\ \frac{e^{2\beta_W\Delta n_C}}3\!\int\!\sum_{v\in C}|\nabla_vf|^2d\mu\ }\tag{5}$$
$$\boxed{\ \lambda_C(\eta,\zeta)\ \ge\ 3e^{-2\beta_W\Delta n_C}\ \overset{\beta=2.5,\Delta=6}{=}\ 3e^{-30n_C}\ }\tag{7}$$

**我核：`n_C=1` ⇒ `3e^{−30}=2.81e−13`；`n_C=4` ⇒ `3e^{−120}=2.30e−52`。**

**⇒ 对一切 `η`、一切 `ζ` 一致，且不依赖整系统体积 `L`。**

$$\boxed{\ |\mathrm{Cov}_{\mu_C}(f,g)|\ \le\ \frac{e^{30n_C}}3\Big(\!\int|\nabla f|^2\Big)^{1/2}\Big(\!\int|\nabla g|^2\Big)^{1/2}\ }\tag{8}$$

**⚠️ 不足：常数随块体积**指数恶化**——真下界，但不足以用于大块尺度。**

### 5⁗⁗⁗⁗⁗⁗⁗‴.2 【`[T]`】good 区域的点态曲率（`κ_good = 1`）

$$\mathcal G_C(\eta)=\Big\{x\in(S^3)^C:\ s_e(x,\eta)\ge\tfrac{14}{15}\ \forall e\in E_C^\star\Big\}\tag{9}$$

$$\mathrm{Hess}V_e(\xi,\xi)=\beta_W|\xi_u-R_e\xi_v|^2-\beta_W(1-s_e)(|\xi_u|^2+|\xi_v|^2)\tag{10}$$
$$\boxed{\ \mathrm{Ric}+\mathrm{Hess}V_C\ \succeq\ \Big(2-\frac{\beta_W\Delta}{15}\Big)I+\beta_WL_R^C\ }\tag{13}$$

$$\beta_W\Delta/15=2.5\cdot6/15=1\ \Longrightarrow\ \boxed{\ \mathrm{Ric}+\mathrm{Hess}V_C\ \succeq\ I+2.5L_R^C\ \succeq\ I\ }\tag{14}$$

**⇒ 一般式 `κ_good = 2 − β_WΔε`；只要 `ε < 2/(β_WΔ)` 就有严格正曲率。**

### 5⁗⁗⁗⁗⁗⁗⁗‴.3 【`[T]`】单顶点快块的 good 条件 Poincaré

**`C={v}` 时，其余全是固定边界，`r_e=R_eη_{w(e)}`：**

$$\mathcal G_{\{v\}}(\eta)=\bigcap_{e\ni v}\Big\{x_v\in S^3:\ x_v^{\sf T}r_e\ge\tfrac{14}{15}\Big\}\tag{19}$$

**每项是半径 `θ_★=arccos(14/15)<π/2` 的**球面测地球** ⇒ 测地凸；有限交仍凸：**

$$\boxed{\ \mathrm{gap}\big(\mu^{\eta,\zeta,G}_{\{v\}}\big)\ \ge\ 1\ }\tag{20}$$

**⇒ 真正强、且对边界方向统一。⚠️ 交集为空时该配置归入坏块分析。**

### 5⁗⁗⁗⁗⁗⁗⁗‴.4 【`[U]`】多顶点 good 区域**不**自动测地凸（用户反例，我核了）

$$p=(e_1,e_1),\ q=(-e_1,-e_1)\in\mathcal G;\qquad x_1(t)=\cos(\pi t)e_1+\sin(\pi t)e_2,\ x_2(t)=\cos(\pi t)e_1+\sin(\pi t)e_3$$
$$t=\tfrac12:\quad x_1\cdot x_2=e_2\cdot e_3=0\ \lt\ \tfrac{14}{15}\ \Longrightarrow\ \text{中间点不在 }\mathcal G$$

$$\boxed{\ \mathcal G_C(\eta)\ \text{一般不具有自动的测地凸性}\ }\tag{18}$$

### 5⁗⁗⁗⁗⁗⁗⁗‴.5 【我的结构定位】`(18)` 的非凸性**只来自 base**，而 `n_C=2` 的 `(23)` 有具体路线

$$\mathcal G=\{(x,y):\ \angle(x,y)\le\theta_\star\},\qquad\theta_\star=\arccos(14/15)=0.367208$$

**坐标 `(u,\varphi,w)`：`u` = 测地中点（base `S³`），`φ∈[0,θ_★]`，`w∈S²`（方向）：**

$$\boxed{\ \mathcal G\ \text{是 }S^3\text{ 上的【凸纤维丛】：纤维 }\{φ\leθ_\star\}\text{ 是 }R^3\text{ 中的球 ⇒ 凸}\ }$$

**⭐ 而 `(18)` 的反例中：`p` 的 `u=e₁`、`q` 的 `u=−e₁` 是**对径点** ⇒ 只有 base 出问题，纤维没问题。**

**⇒ 纤维谱隙 `~π²/θ_★² = 73.2`；base 谱隙 `= 3`；warp 是 `cosφ`，而**

$$(\log\cos)''=-\sec^2\lt0\ \Longrightarrow\ \text{log 凹}\ \Longrightarrow\ \text{warped product 的 BE 曲率}\ \ge\ \min(\text{base},\text{fiber})$$

$$\boxed{\ \mathrm{gap}(\mathcal G)\ \gtrsim\ \min(3,\,73.2)/O(1)\ =\ O(1)\ \Longrightarrow\ c_G(2)\ \text{是 }O(1)\text{，\textbf{不}指数小}\ }$$

**⇒ 而 base `S³` 的非凸只需 `O(1)` 块覆盖：**

| 覆盖 | 重叠 |
|---|---|
**2 个半球** | **赤道带，`O(1)`** |

**⇒ 这与用户指出的第二条路线一致（划分为有限个局部凸区 + 受控重叠 + 拼接 Poincaré，并注意不同分支间的零能量跳跃模式）。**

### 5⁗⁗⁗⁗⁗⁗⁗‴.6 两层结果对照

| 条件 | 谱隙下界 | 适用 |
|---|---|---|
任意构型、任意边界 | `3e^{−30n_C}` | 任意有限快分量 |
所有边 good 的点态曲率 | `κ=1` | `𝒢_C` 上 |
good 区域测地凸时的条件谱隙 | `1` | 测地凸 good 区域 |
**单顶点快块的 good 条件谱隙** | **`1`** | 任意非空单顶点 good 集合 |
**`n_C=2` 的 `c_G`（我的路线）** | **`O(1)`** | 尚待写成定理 `[U]` |

### 5⁗⁗⁗⁗⁗⁗⁗‴.7 审计（本轮）

| 命题 | 判 |
|---|---|
任意边界的有限块 Poincaré `(5)–(7)` | **`[T]`** |
协方差推论 `(8)` | **`[T]`** |
good 区域点态曲率 `(13)–(14)`（`κ_good=1`） | **`[T]`** |
**单顶点 good 快块条件 Poincaré `(20)–(21)`** | **`[T]`** |
**多顶点 good 区域不自动测地凸 `(18)`** | **`[T]`（反例）** |
**`n_C=2` 的 `c_G` 为 `O(1)`（凸纤维丛 + 半球覆盖）** | **`[U]`（有具体路线）** |
一般多顶点 good 条件 Poincaré `(23)` | **`[U]`** |

**⇒ 下一步按用户建议：从**两个快顶点 + 一条内部连接边 + 至少一个固定边界锚点**的块开始。**

**⇒ 而 §5⁗⁗⁗⁗⁗⁗⁗‴.5 给出的路线是：`𝒢` = 凸纤维丛 → 纤维凸、base 需 `O(1)` 块 → 拼接。**

--- 

## §5⁗⁗⁗⁗⁗⁗⁗⁗ 移除 good 区域依赖：全空间 Poincaré `[T]` —— 及其代价（第 284 轮）

### 5⁗⁗⁗⁗⁗⁗⁗⁗.1 【`[T]`】任意多顶点快块的全空间条件 Poincaré

$$\boxed{\ \mathrm{Var}_{\mu_C^{\eta,\zeta}}(f)\ \le\ \frac13e^{2\beta_W|E_C^\star|}\!\int_{\mathcal M_C}\!\sum_{v\in C}|\nabla_vf|^2d\mu_C\ }\tag{3}$$
$$|E_C^\star|\le\Delta n_C\ \Longrightarrow\ \boxed{\ \mathrm{Var}\le\frac{e^{2\beta_W\Delta n_C}}3\!\int\!\!\sum|\nabla_vf|^2,\qquad \lambda_C\ge3e^{-30n_C}\ }\tag{4),(7)}$$

$$\boxed{\ C_P(C)\ \le\ \frac{e^{5|E_C^\star|}}3\ }\qquad(\beta=2.5,\ 2\beta=5)\tag{8}$$

**⇒ 对一切 `η`、一切 `ζ`、**一切构型（含坏构型）**一致；不依赖全局体积。**

**⚠️ 这是**全空间**结果——**不需要** `𝒢_C` 测地凸。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗.2 【`[T]`】快块边界平均响应

$$|J_{C,a}|\le w_C(a),\qquad \int|\nabla_{x_C}A_a|^2d\mu_C\le w_C(a)^2|h_a|^2\tag{9),(12)}$$

**（`∇_vA_a=P_{x_v^\perp}(Σ_{e=(a,v)}R_e^{\sf T}h_a)`；投影不增范数；`Σ_vm(a,v)^2≤(Σ_vm(a,v))^2=w_C(a)^2` ✓）**

$$\boxed{\ \big|D_a\mathbb E_{\mu_C}[F][h_a]\big|\ \le\ \beta_WC_P(C)\,w_C(a)|h_a|\Big(\!\int|\nabla F|^2d\mu_C\Big)^{1/2}\ }\tag{13}$$
$$\boxed{\ \le\ \frac{\beta_W e^{30n_C}}3w_C(a)|h_a|\,\|\nabla F\|_{L^2(\mu_C)}\ }\tag{14}$$

### 5⁗⁗⁗⁗⁗⁗⁗⁗.3 【`[T]`】但 `(16)` 的 `min{1,C_P}` **恒取 1** —— 全空间 Poincaré 不改善 Hessian

$$\big|D^2_{ab}(-\log Z_C)\big|\ \le\ \beta_W^2w_C(a)w_C(b)\min\{1,C_P(C)\}|h_a||h_b|\tag{16}$$

$$C_P\lt1\iff5|E_C^\star|\lt\ln3=1.0986\iff|E_C^\star|\lt0.2197\ \Longrightarrow\ \text{不可能}\ (|E_C^\star|\ge1)$$

$$\boxed{\ C_P\gt1\ \text{恒成立}\ \Longrightarrow\ \min\{1,C_P\}=1\ \text{恒成立}\ \Longrightarrow\ (16)\ \textbf{永不改善}\ }$$

**⚠️ 这是本轮最重要的诚实结论：全空间 Poincaré 补上了**快块函数响应**控制，但**没有**改善 `−log Z_C` 混合 Hessian 的数值常数。两件事不能混。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗.4 【我的定量审计】移除 `[U]` 的代价：常数是 `10^{12}`

**冻结 `n_0`（块大小上限）使 `C_P^block` **`L`-无关**：**

| `n_0` | `C_P = e^{30n_0}/3` | `(19)` 的系数 `β_WC_P` |
|---|---|---|
`1` | **`3.562×10^{12}`** | **`8.905×10^{12}`** |
`2` | `3.807×10^{25}` | `9.517×10^{25}` |
`4` | `4.347×10^{51}` | `1.087×10^{52}` |

**⇒ `(19)` 在 `n_0=1`、`w=1`、`|h|=1`、`‖∇F‖=1` 时给 `~9×10^{12}`。**

$$\boxed{\ \text{模块【逻辑闭合】（`[U]` 被移除）但【定量不可用】（常数 }10^{12}\text{）}\ }$$

**⇒ 这是一个明确的 **trade-off**（我把它写成正式账）：**

| 路线 | 常数 | 判 |
|---|---|---|
全空间（本轮） | `~e^{30n_C}`（`10^{12}`） | **`[T]` 逻辑闭合，定量不可用** |
good 区域 | `O(1)`（有用） | **`[U]`（`n_C≥2` 测地凸不成立）** |

**⇒ 本轮的选择是**用逻辑闭合换定量强度**——这是合理的工程取舍，但必须记明：这个模块**不能**喂给后续定量链条。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗.5 good 曲率结果的**保留范围**（收窄但保留）

$$\mathrm{Ric}+\mathrm{Hess}V_C\ \succeq\ I+2.5L_R^C\qquad\text{在 }\mathcal G_C\ \text{上}\tag{17}$$

**⚠️ 只用于：good 构型上的局部能量/曲率评估；好坏块的局部几何诊断；对全空间积分的分区贡献。**

**⚠️ **不再**从 `(17)` 直接推多变量 good 条件 LSI（因 `𝒢_C` 一般不必测地凸，见 §5⁗⁗⁗⁗⁗⁗⁗‴.4 的反例）。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗.6 修改后的依赖链

$$\boxed{\ \text{有限快块全空间 Poincaré}\ \Downarrow\ \text{快块边界平均响应}\ \Downarrow\ \text{精确有效慢作用量}\ \Downarrow\ \text{块间联合响应与坏块控制}\ }$$

**⇒ 不再要求中间经过「多顶点 good 区域条件 LSI」。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗.7 审计（本轮）

| 命题 | 判 |
|---|---|
任意多顶点快块全空间 Poincaré `(3)–(4)` | **`[T]`** |
`(8)` 用实际 `|E_C^\star|` | **`[T]`** |
快块边界平均响应 `(13)–(14)` | **`[T]`** |
`∫\|∇A_a\|² ≤ w(a)²\|h_a\|²` `(12)` | **`[T]`** |
**`min{1,C_P}=1` 恒成立 ⇒ `(16)` 永不改善** | **`[T]`（我算）** |
**移除 `[U]` 的代价：常数 `~10^{12}`** | **`[T]`（我算）** |
**多顶点 good 区域条件 Poincaré/LSI** | **`[U]`（但主线不再依赖它）** |
坏块连通簇联合概率指数尾 | **`[U]`** |
块间联合响应（不只最坏环境上界） | **`[U]`** |

**⇒ 关闭的是**快块控制模块的依赖缺口**，不是整个 CLGF。**

--- 

## §5⁗⁗⁗⁗⁗⁗⁗⁗′ D 模块：有效慢块联合响应（第一阶段，第 285 轮）

### 5⁗⁗⁗⁗⁗⁗⁗⁗′.1 【`[T]`】五个结果（我逐步核过）

$$W_{C,i}:=\!\!\sum_{a\in B_i\cap\partial C}\!\!w_C(a)\tag{4}$$

**① 一阶响应**（`D_iF̄_C = βCov(F_C,A_{C,i})`，用 `¼osc·osc`）：

$$\boxed{\ |D_i\overline F_C[h_i]|\ \le\ \beta_WMW_{C,i}\|h_i\|_\infty\ }\tag{10}\qquad \boxed{\ \mathrm{osc}_{B_i}\overline F_C\ \le\ \pi\beta_WMW_{C,i}\ }\tag{11}$$

**⭐⭐ **这里**不用** Poincaré ⇒ **常数里没有 `e^{30|C|}`**。**

**② 二块联合响应**（`D_jD_iF̄_C = β²κ(F_C,A_i,A_j)`）：

$$\boxed{\ \Delta_{ij}\overline F_C\ \le\ 8\pi^2\beta_W^2MW_{C,i}W_{C,j}\ }\tag{15}$$

**③ 有效因子二块混合 Hessian**（`D_jD_iΦ_C = −β²Cov(A_i,A_j)`）：

$$\boxed{\ \Delta_{ij}\Phi_C\ \le\ \pi^2\beta_W^2W_{C,i}W_{C,j}\ }\tag{19}$$

**④ 全部有效慢作用量**：

$$\boxed{\ \Delta^{\rm eff}_{ij}\ \le\ 4\beta_WN^{SS}_{ij}+\pi^2\beta_W^2\!\!\sum_C\!\!W_{C,i}W_{C,j}\ }\tag{20}$$

**⑤ 有限范围**：若两块既无直接慢—慢边、也不共享快分量，则 `Δ^eff_{ij}=0`。**

$$\beta_W=2.5:\quad 4\beta=10.0,\quad \pi^2\beta^2=61.685,\quad 8\pi^2\beta^2=493.48\tag{21),(22)}$$

### 5⁗⁗⁗⁗⁗⁗⁗⁗′.2 【大改善】与上一轮 `10^{12}` 对比：降了 `1.4×10^{11}` 倍

| 轮 | 系数 |
|---|---|
第 284 轮 `(19)`（含 Poincaré） | **`β_WC_P = 8.905×10^{12}`** |
**本轮 `(18)`（不用 Poincaré，`W=1`）** | **`β_W²W_iW_j = 61.69`** |

$$\boxed{\ \text{改善 }1.4\times10^{11}\text{ 倍}\ }$$

**⇒ 本条是本轮的真实收益：用 `¼osc·osc` / 方差界替掉 Poincaré，**指数因子消失**。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗′.3 【我的核对】`(25)` **不是约束**——它被 `(20)` 的有限范围**平凡满足**

$$\text{(20) 已给：不相邻块}\Rightarrow\Delta^{\rm eff}_{ij}=0\Rightarrow\mathcal J_{ij}=0$$
$$\Longrightarrow\quad \mathbb E_{\nu_{\rm slow}}\big[|\mathcal J_{ij}(Y)|\big]\ \le\ C e^{-md(B_i,B_j)}\quad\text{对【任意 }m\gt0\text{】自动成立}$$

$$\boxed{\ \text{(25) 是 (20) 的推论，不是新瓶颈}\ }$$

**⇒ 真正的瓶颈是 `tanh(Δ^eff/4)` 能否**小于 1 且有收缩**。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗′.4 【定量判决】收缩需 `Δ^eff ≲ 0.50`，两条通道给 `β ≲ 0.126`

**块坐标数 `~8` ⇒ 需 `tanh(Δ/4) < 1/8`：**

$$\Delta^{\rm eff}\ \lt\ 4\operatorname{artanh}(0.125)\ =\ 0.5026$$

| 通道 | 条件 | `β` 阈值 |
|---|---|---|
A 共享快分量 | `π²β² < 0.5026` | `β < 0.2257` |
**B 直接慢—慢边** | **`4β < 0.5026`** | **`β < 0.1257`** |

$$\boxed{\ \beta\ \lesssim\ 0.1257\qquad\frac{2.5}{0.1257}=19.9\times\ }$$

**而 `β=2.5` 时两条通道的 `c`：**

| 通道 | `Δ^eff ≥` | `c = tanh(Δ/4)` |
|---|---|---|
**A** | `π²β² = 61.69` | **`1.0000000000`** |
B | `4β = 10.0` | `0.9866143` |

**⇒ 两条都 `c ≈ 1` ⇒ 无收缩。**

**⚠️ 而这个 `β < 0.1257` 与第 282 轮反向祖先那条**完全一致**（`Δtanh(β)<1`）——因为**根是同一个**：一条慢—慢 Wilson 因子的交叉振幅就是 `4β`。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗′.5 你指出的方向对：`(20)` 是安全上界，不是最优

**你已指出：单快顶点 + `R_e=I` 时 `Z_C` 只依赖 `|Σ_ay_a|`，多块相互作用有明显**方向结构**，粗界 `(20)` 不反映它。**

**⇒ 所以下一手的正确对象是：**

$$\mathcal J_{ij}=\sum_C\mathrm{Cov}_{\mu_{C,y}}\big(J_{C,B_i},J_{C,B_j}\big)\tag{24}$$

**而它的**绝对值上界** `|𝒥_{ij}| ≤ Σ_CW_{C,i}W_{C,j}` 仍是最坏情形。**

**⇒ 而我在 §5⁗⁗⁗⁗⁗⁗⁗‴.5 已算出：`Cov` 在**铁磁有序**环境（`|h_0|` 大）远小于振幅乘积；`|h_0| ≤ 2d = 6` 是图度数上界 ⇒ `Δ` 有正下界 ⇒ 收缩不可能。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗′.6 审计（本轮）

| 命题 | 判 |
|---|---|
单快分量一阶块响应 `(9)–(11)` | **`[T]`** |
两块联合响应 `(12)–(15)` | **`[T]`** |
`−log Z_C` 二块混合 Hessian `(17)–(19)` | **`[T]`** |
全部有效慢作用量 `(20)` | **`[T]`** |
有效慢块的**有限范围**结构 | **`[T]`** |
**不用 Poincaré ⇒ 常数从 `10^{12}` 降到 `61.69`** | **`[T]`（我算：改善 `1.4×10^{11}`）** |
**`(25)` 被 `(20)` 平凡满足，不是瓶颈** | **`[T]`（我算）** |
**收缩需 `β ≲ 0.126`（两条通道取更紧）** | **`[T]`（我算）** |
`β=2.5` 下块影响矩阵收缩 | **`[U]`/`[F]`（两通道 `c ≈ 1`）** |
`𝒥_{ij}` 的退火联合衰减（真正的内容） | **`[U]`** |

**⚠️ 但必须记明：`(24)` 的绝对上界仍是**最坏情形**；真正的退火内容不在 `(25)`，而在 `Cov` 的**方向结构**（`|Σ_ay_a|` 型依赖）。**

--- 

## §5⁗⁗⁗⁗⁗⁗⁗⁗″ 两项审计：块收缩 `[F]`（pinning）+ 联合退火 `[U]`；`(17)` 的 Wilson 对应物（第 286 轮）

### 5⁗⁗⁗⁗⁗⁗⁗⁗″.1 【`[F]`】均匀块 Dobrushin 收缩受 pinning 障碍

$$\sup_i\sum_{j\ne i}e^{\alpha d(B_i,B_j)}c_{ij}\le q\lt1\tag{2}\ \Longrightarrow\ |\nu^\eta_\Lambda(F_A)-\nu^{\eta'}_\Lambda(F_A)|\le\frac{2|A|\|F_A\|_\infty}{1-q}e^{-\alpha d(A,\partial\Lambda)}\tag{3}$$

**而 pinning 嵌入给（在它成立的前提下）**

$$\exists F_A:\quad\limsup_{R\to\infty}|\nu_R^\eta(F_A)-\nu_R^{\eta'}(F_A)|\gt0\tag{4}$$

$$\boxed{\ \text{均匀块 Dobrushin 收缩不成立}\ }\tag{5}$$

**⚠️ 这是**条件性 no-go**：正式写入 `[F]` 还须把 pinning 子系统的精确嵌入与耦合常数归一化一并固定。它**不**否定真实退火 Wilson 测度的指数聚簇。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗″.2 【`[T]`】Wilson 联合响应的精确恒等式与有限范围

$$\boxed{\ D_iD_j(-\log Z_C)[h_i,h_j]=-\beta_W^2\mathrm{Cov}_{\mu_{C,y}}\big(A_{C,i}(h_i),A_{C,j}(h_j)\big)\ }\tag{8}$$
$$\boxed{\ \|D_iD_jS_{\rm slow}(y)\|\ \le\ \beta_WN^{SS}_{ij}+\beta_W^2\!\!\sum_C\!\!W_{C,i}W_{C,j}\ }\tag{10}$$

$$\beta=2.5:\quad \|D_iD_jS_{\rm slow}\|\ \le\ 2.5N^{SS}_{ij}+6.25\!\!\sum_C\!\!W_{C,i}W_{C,j}\tag{11}$$

$$\boxed{\ D_iD_jS_{\rm slow}\equiv0\quad\text{若 }d(B_i,B_j)\gt R_*\ }\tag{13}$$

**⚠️ 有限范围是**相互作用支撑**性质，不是 Gibbs 测度的聚簇性质——有限范围作用量仍可有长程关联。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗″.3 【我的核对】`(11)` 比上一轮紧 `4×` 与 `π²×`，但收缩上**无净收益**

| 形式 | 常数 |
|---|---|
四点振幅（第 285 轮 `(21)`） | `4β_WN + π²β²ΣWW = 10.0N + 61.685ΣWW` |
**混合 Hessian（本轮 `(11)`）** | **`β_WN + β²ΣWW = 2.5N + 6.25ΣWW`** |

$$\text{比值：}N\text{ 项 }4\times;\quad ΣWW\text{ 项 }\pi^2=9.87\times$$

**⚠️ 但收缩判据用的是**交叉振幅**：**

$$\Delta\approx\pi^2\|D^2\|\ \Longrightarrow\ \pi^2(2.5N+6.25ΣWW)=24.7N+61.69ΣWW$$

$$\boxed{\ \text{与上一轮【完全相同】}\ \Longrightarrow\ \text{收缩问题上无净收益}\ }$$

### 5⁗⁗⁗⁗⁗⁗⁗⁗″.4 【`[T]`】`(17)` 反例：单边平均小 `\not\Longrightarrow` 路径联合衰减

$$R_1=\cdots=R_n=Z,\quad \mathbb P(Z=1)=p\quad\Longrightarrow\quad \mathbb E R_k=p\lt1\ \text{但}\ \mathbb E\prod_kR_k=p\ \text{与 }n\text{ 无关}$$

$$\boxed{\ \mathbb E R_k\lt1\ \forall k\ \not\Longrightarrow\ \mathbb E\prod_kR_k\le Ce^{-mn}\ }\tag{17}$$

### 5⁗⁗⁗⁗⁗⁗⁗⁗″.5 【我的结构定位】`(17)` 在 Wilson 里的对应物就是 **base 场 `ζ`**

**`ζ` 活在 `Λ_3 ⊂ E_f` 上，维度 `O(L^4)`；**每个**快分量 `C` 的 `μ_{C,y}` 都依赖**同一个** `ζ`。**

$$\boxed{\ ζ\ \text{是【全局公共环境】}\ \Longrightarrow\ (17)\ \text{的机制在 Wilson 中是\textbf{活的}}\ }$$

$$\mathbb E\Big[\prod_k\mathcal R_{i_ki_{k+1}}(Y)\Big]=\mathbb E_\zeta\Big[\mathbb E\big[\textstyle\prod_k\mathcal R\mid\zeta\big]\Big]$$

**⇒ 内层对路径是**可乘**的；外层对 `ζ` 平均**可以把 `n`-衰减抵消掉**。这正是 `(17)`。**

**⭐⭐ 而这解释了一个贯穿本次会话的现象：**

| 路线 | condition 在什么上 | 那个东西就是公共环境 |
|---|---|---|
静态退火 | 同一个 `Y∼ν_slow` | `Y` |
动态退火 | 耦合过程 | **`ζ` 仍是** |
路径消元 | `ν^P` | `ζ` 仍是 |

$$\boxed{\ \text{每一条退火路线都 condition 在某个东西上，而那个东西就是公共环境}\ }$$

**⇒ 要真正摆脱，**必须**：证明 `ζ` 的相关性足够弱，**或**把 `ζ` 的响应也纳入路径权重（即 `ζ` 必须被积分进慢变量图——这与第 275 轮 §5⁗⁗⁗⁗⁗‴.5 的记账是同一条）。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗″.6 【`[U]`】退火多块联合路径估计

$$\mathfrak P_n(i_0,i_n)=\!\!\sum_{i_1\cdots i_{n-1}}\!\!\mathbb E_{\nu_{\rm slow}}\Big[\prod_{k=0}^{n-1}Q_{i_ki_{k+1}}(Y)\Big]\tag{18}$$
$$\mathbb E_{\nu_{\rm slow}}\Big[\prod_{k=0}^{n-1}Q_{i_ki_{k+1}}(Y)\Big]\le C_0\vartheta^n,\quad \vartheta\lt1\tag{20}$$

**⚠️ `(8)–(15)` **不**给出 `(20)`。它需要进一步的 Gibbs 联合概率、反向 Poincaré 或**多尺度好坏块**输入。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗″.7 审计（本轮）

| 命题 | 判 |
|---|---|
所有边界条件下块 Dobrushin 收缩 | **`[F]`（条件性：pinning 有序子系统）** |
Wilson 有效块混合 Hessian 恒等式 `(8)` | **`[T]`** |
显式上界 `(10)–(11)` | **`[T]`** |
响应矩阵有限范围 `(12)–(13)` | **`[T]`** |
单块退火一阶矩 `(15)` | **`[T]`** |
**单块平均小 `\not\Longrightarrow` 路径联合衰减 `(17)`** | **`[T]`（反例）** |
**`(17)` 的 Wilson 对应物 = `ζ`（全局公共环境）** | **`[T]`（我补的结构定位）** |
`(11)` 更紧但收缩上无净收益 | **`[T]`（我算）** |
当前 Wilson 多块退火联合衰减 `(20)` | **`[U]`** |
`β=2.5` 下退火块响应路径收敛 | **`[U]`** |

### 5⁗⁗⁗⁗⁗⁗⁗⁗″.8 架构调整（正式记账）

$$\boxed{\ \text{最坏情形块收缩：放弃}\ \succ\ \text{局部块响应精确公式：完成}\ \succ\ \text{有限范围：完成}\ \succ\ \text{退火多块联合路径：未完成}\ \succ\ \text{全局双块聚簇：未完成}\ }$$

**⇒ 结论：Wilson 有效慢块响应的**精确有限范围**与**单块退火控制**已证；但它们在真实 Gibbs 测度下的**多块联合衰减未证**。**

**⇒ 下一手不是重新估计块矩阵的最大行和，而是直接分析 `(20)` 的联合乘积，利用 `ζ` 的协方差结构排除「公共坏环境」造成的长期相关。**

--- 

## §5⁗⁗⁗⁗⁗⁗⁗⁗‴ 复制 Gibbs 系统：两个 `[U]` 合并为一个对象（第 287 轮）

### 5⁗⁗⁗⁗⁗⁗⁗⁗‴.1 【`[T]`】复制恒等式

$$\boxed{\ \mathrm{Cov}_\mu(A,B)=\tfrac12\mathbb E_{\mu\otimes\mu}\big[(A(X)-A(X'))(B(X)-B(X'))\big]\ }\tag{5}$$

**（展开：`½(2E[AB]−2E[A]E[B])` ✓）**

$$\mathcal J_{C,ij}(y)\ \le\ \tfrac12\mathbb E_{\mu_{C,y}^{\otimes2}}\big[D_{C,i}(X,X')D_{C,j}(X,X')\big]\tag{6}\qquad |\mathcal J_{C,ij}|\le W_{C,i}W_{C,j}\tag{8}$$

### 5⁗⁗⁗⁗⁗⁗⁗⁗‴.2 【`[T]`】多块联合响应的复制表示

$$\widehat\nu_{\boldsymbol C}(dY,dX)=\nu_{\rm slow}(dY)\prod_{k=1}^n\mu_{C_k,Y}(dX_k^+)\mu_{C_k,Y}(dX_k^-)\tag{9}$$

$$\boxed{\ \mathbb E_{\nu_{\rm slow}}\prod_{k=1}^n\big|\mathrm{Cov}_{\mu_{C_k,Y}}(A_{i_{k-1}},A_{i_k})\big|\ \le\ \frac1{2^n}\mathbb E_{\widehat\nu_{\boldsymbol C}}\prod_{k=1}^n\big[D_{C_k,i_{k-1}}D_{C_k,i_k}\big]\ }\tag{10}$$

**⭐⭐ 联合概率空间的**结构**是本轮要点：所有路径边**共享同一个 `Y`**；每条边有自己的复制对；给定 `Y` 后复制对**独立**；它们**通过共同的 `Y` 产生相关性**。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗‴.3 【`[T]`】复制测度的精确结构与**一个可读的等价描述**

$$d\widehat\nu_{\boldsymbol C}\propto e^{-S_{\rm slow}(Y)}\prod_{k=1}^n\frac{e^{H_{C_k}(X_k^+;Y)}e^{H_{C_k}(X_k^-;Y)}}{Z_{C_k}(Y)^2}d\sigma(Y)\prod_kdX_k^+dX_k^-\tag{12}$$

$$\boxed{\ \widehat\nu_{\boldsymbol C}(dY)=\nu_{\rm slow}(dY)\ }\tag{13}$$

**每个复制对与慢变量的耦合为（我核过这个形式）：**

$$\beta_W\!\!\sum_{a\in\partial C}\!\!Y_a^{\sf T}\big[J_{C,a}(X^+)+J_{C,a}(X^-)\big]\ -\ 2\log Z_C(Y)\tag{14}$$

$$\boxed{\ \text{复制系统}\ =\ \text{原 Wilson 模型，但每个快分量的贡献【加倍】}\ }$$

**⇒ 所以 `ρ_ann<1` 问的是：**在这个「加倍模型」里，沿路径的响应是否衰减**。而 `(13)` 说 `Y` 的边缘不变——但**联合**测度是另一个模型。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗‴.4 【`[T]`】`ρ_ann<1` 推出路径衰减（与最坏情形无关）

$$\Xi_n=\sup_i\mathbb E_{\nu_{\rm slow}}\Big[\sum_j(\mathsf R(Y)^n)_{ij}\Big]\overset{(20)}{=}\sup_i\!\!\sum_{i_1\cdots i_n}\!\!\mathbb E_{\nu_{\rm slow}}\Big[\prod_{k=1}^nR_{i_{k-1}i_k}(Y)\Big]$$
$$\rho_{\rm ann}=\limsup_{n\to\infty}\Xi_n^{1/n}\tag{21}$$

$$\boxed{\ \rho_{\rm ann}\lt1\ \Longrightarrow\ \Xi_n\le C\theta^n\ \Longrightarrow\ \text{块响应路径和指数衰减}\ }\tag{24}$$

**⇒ **完全不需要**最坏情形 Dobrushin 行和 `<1`。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗‴.5 【我的核对】`(11)` 的 `Σ_{C_1⋯C_n}` **受共同接触约束**，但 `ρ_ann` 粗界仍 `≫1`

**⚠️ **必须**注意：`Σ_{C_k}` 只跑**同时接触 `i_{k-1}` 与 `i_k`** 的快分量 ⇒ 个数 `O(1)`（有限范围），**不是** `N_C∼L^4`。否则 `(27)` 平凡不成立。**

**⇒ 于是 `(27)` LHS `\le (2\beta_W^2W^2N_{\rm touch})^n$：**

| `N_touch` | `2β²W²N` |
|---|---|
`1` | **`12.50`** |
`2` | `25.00` |

**⇒ 每边增长因子 `≥12.50` ⇒ `ρ_ann` **不可能** `<1`。**

**而粗界直接给（`W=1`、一个共享 `C`）：**

$$R^{\rm fast}_{ij}\ \le\ \beta_W^2\cdot1=6.25\ \Longrightarrow\ \text{行和}\approx6.25\Delta_B$$

| `Δ_B` | 行和 |
|---|---|
`4` | `25.0` |
`8` | `50.0` |

$$\boxed{\ \rho_{\rm ann}\ \gtrsim\ \beta_W^2\Delta_B=50\ \gg\ 1\ }$$

**⇒ 要 `ρ_ann<1` 需 `β²W²<1/Δ_B` ⇒ `W\lt0.141`，而 `W_{C,i}≥1` 只要 `B_i∩∂C≠∅` ⇒ **不可能**。**

**⚠️ 但这是**粗界**。复制公式 `(26)` 的价值正在于它**允许** `E_{\hatν}[∏DD]` 远小于 `(2W)^{2n}`——即允许两个复制副本**趋于一致** ⇒ `D_k→0`。**

$$\boxed{\ \text{真正需要的是：复制的两个快变量副本在共同慢背景下【趋于一致】，使 }D_k\ \text{沿路径衰减}\ }$$

### 5⁗⁗⁗⁗⁗⁗⁗⁗‴.6 剩余的唯一目标 `(27)`

$$\sup_{i_0}\!\!\sum_{\substack{i_1\cdots i_n\\ C_1\cdots C_n}}\!\!\frac{\beta_W^{2n}}{2^n}\mathbb E_{\widehat\nu_{\boldsymbol C}}\prod_{k=1}^nD_{C_k,i_{k-1}}D_{C_k,i_k}\ \le\ C_0\theta^n\tag{27}$$

**⚠️ 现有 Poincaré、Hessian、边界力上界**都不足以**推出它。**

**⚠️ `(17)` 的粗界 `∏[2.5N+6.25ΣWW]` **无任何衰减**——它完全没用到 `(12)` 的共同慢背景结构。**

### 5⁗⁗⁗⁗⁗⁗⁗⁗‴.7 审计（本轮）

| 命题 | 判 |
|---|---|
Wilson 二块响应复制恒等式 `(5)–(6)` | **`[T]`** |
多块路径联合响应复制表示 `(9)–(11)` | **`[T]`** |
复制测度精确 Gibbs 表达 `(12)–(14)` | **`[T]`** |
退火响应矩阵有限范围 `(18)` | **`[T]`** |
`ρ_ann<1` ⇒ 联合路径衰减 `(21)–(24)` | **`[T]`** |
**`Σ_{C_k}` 受共同接触约束（`O(1)` 项）** | **`[T]`（我核）** |
**粗界给 `ρ_ann ≳ β²Δ_B = 50 ≫ 1`** | **`[T]`（我算）** |
所有边界下最坏情形块收缩 | **`[F]`（pinning，前提性）** |
`β=2.5` 下复制响应路径满足 `(27)` | **`[U]`** |
`β=2.5` 下 `ρ_ann<1` | **`[U]`** |

### 5⁗⁗⁗⁗⁗⁗⁗⁗‴.8 主线（正式固定）

$$\boxed{\ \text{精确快慢消元}\succ\text{边界力协方差}\succ\text{多块联合响应复制恒等式}\succ\text{复制 Gibbs 的路径响应尾}\succ\rho_{\rm ann}\lt1\succ\text{有效慢块联合衰减}\ }$$

**⇒ 本轮完成了从「Wilson 有效慢块联合响应」到「复制 Gibbs 系统」的**严格转换**，并证明了 `ρ_ann` 控制联合路径衰减。**

**⇒ 而 `(27)` 在 `β=2.5` 下**不**标 `[T]`：需要的是一条**真正的联合响应估计**（两个复制副本趋于一致），而不是把单块均值相乘。**

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
