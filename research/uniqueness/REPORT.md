# 无限体积极限的唯一性、真空唯一性，与普适性

> 2026-11-22 · 范围：4D 纯 Yang–Mills（紧单群 $G$，主情形 $SU(N)$），Wilson 格点正则化。
> 前置：`docs/yang-mills-prior-work.md`；`lean/YMStatement/Basic.lean`；`research/physics_handoff.md`（弱→强 RG 交接）。
> 本报告只依赖前两份文档对 16 输入的**摘要**，不重读 2611.00492 原文；凡涉该文具体形式处均标注"待核实"。

---

## 0. 给主会话的 TL;DR

1. **"真空唯一"与"极限唯一"是两个不同的要求，而且强度差别很大。**
   - *每个*子列极限是一个 OS 测度、满足 OS4（聚类/遍历性）。**遍历 ⇒ 该极限的真空向量唯一（一维不变子空间）**。所以"极限理论内部真空唯一"几乎免费——它已是输入集的产物（OS4 在 16 输入里）。
   - **真正的缺口是"子列极限本身唯一"**：不同子列给出不同理论，每个各有一个（自己的）唯一真空。把"子列收敛"升级为"全网收敛"需要的是**极限的唯一性**，不是真空的唯一性。

2. **极限唯一性与 MG（指数聚类）不是同一件事，MG 既非必要也非充分（以单态形式）。**
   - **充分且最小**的输入是**边界屏蔽性 W1**：有限体积 DLR 规范下，局部可观测量的条件期望对边界数据的依赖在 $d\to\infty$ 时趋于 $0$（速率任意）。W1 ⇒ 唯一 Gibbs 态 ⇒ 极限唯一 ⇒ 真空唯一。
   - **MG-单态 ⇒ 不推出 W1**：二维 Ising 在 $\beta>\beta_c$ 时，$+$ 相的有限体积关联一致指数衰减（MG 成立），但模型有 $\mu_\pm$ 两个 Gibbs 态，`+bc/−bc` 给不同极限——**唯一性失败**。所以若 hexagonmath 的 MG 只针对单一态，它**不能**推出唯一性（待核对其是否已对边界条件一致量化）。
   - **W1 ⇒ 不推出 MG**：二维 Ising 在 $\beta_c$ 有唯一 Gibbs 态、边界影响 $\to0$（幂律），但关联幂律衰减、无质量隙。所以唯一性**严格弱于**质量隙。
   - **打包的正确目标是 MG-spec**：把 MG 加强为"对一切边界条件、对一切体积**一致**的指数边界衰减"。MG-spec ⇒ W1 ⇒ 唯一性，且 ⇒ MG-inf ⇒ 质量隙。

3. **普适性（跨正则化/跨 blocking 无关）需要的是 RG 稳定性，不是 MG。** 同一正则化内"跨穷尽/跨边界条件无关"由 W1 给出；跨"格距、作用量、blocking 方案"无关需要"无关算符控制 + 相关耦合匹配"（RG 输入），MG 完全不涉及。两者要分开记账。

4. **对物理直觉的校准**：hexagonmath 的反例（E4 聚类成立但无谱隙）恰好是一个**真空唯一但无质量隙**的理论——这与"唯一性 ⇐ W1 ⇐ MG-spec"以及"唯一性 ⇏ 质量隙"完全一致。不要把"唯一真空"当成质量隙的进展。

---

## 1. 精确设定

### 1.1 格点 YM 与 DLR 规范

$G$ 紧连通**单连通**单李群（取单连通覆盖 $\tilde G$；如 $SU(N)$）。记 $N_c:=\dim f$（基本表示维数），$\mathrm dU$ 归一化 Haar 测度。对有限盒 $\Lambda\Subset\mathbb Z^4$，记 $E(\Lambda)$（键）、$P(\Lambda)$（面）。**约定** $\beta:=2N_c/g^2$、$\gamma:=\beta/N_c=2/g^2$（$\gamma$ 是 $\mathrm{Re}\,\mathrm{tr}\,U_p$ 的系数，即共享测度 A0 的不变量）。对 $\beta>0$，Wilson 有限体积测度
$$
\mathrm d\mu_{\Lambda}^{\xi,\beta}(U_\Lambda)
=\frac1{Z_\Lambda^{\xi,\beta}}\;
\exp\!\Big(\frac{\beta}{N_c}\sum_{p\in P(\Lambda)}\mathrm{Re}\,\mathrm{tr}\,U_p\Big)
\prod_{b\in E(\Lambda)}\mathrm dU_b
=\frac1Z\exp\!\Big(\gamma\sum_{p}\mathrm{Re}\,\mathrm{tr}\,U_p\Big)\prod_{b}\mathrm dU_b ,
$$
等价地（差一常数）$\exp\big(\beta\sum_p(1-\tfrac{1}{N_c}\mathrm{Re}\,\mathrm{tr}\,U_p)\big)$。边界条件 $\xi$ 给定 $E(\mathbb Z^4)\setminus E(\Lambda)$ 上的组态（只有与 $P(\Lambda)$ 相邻的边界键进入作用量）。

> **归一化警告（GL1/A0）**：系数是 $\beta/N_c$，**不是** $\beta$。$SU(3)$、$g^2=1$ 给 $\beta=2N_c/g^2=6$、$\gamma=\beta/N_c=2$；若误写 $\exp(6\,\mathrm{Re}\,\mathrm{tr})$，则对应 $\gamma=6=2/g'^2\Rightarrow g'^2=1/3$，是**另一个理论**。

> **重整化轨迹（RT）**：$\beta(a)=2N_c/g^2(a)$ 沿渐近自由轨迹，但**跑动耦合用 $\lambda_G:=h^\vee g^2$**（$h^\vee$ 为双 Coxeter 数，**不是** $\dim f$）：一圈 $1/g^2(\mu)=2b_0\ln(\mu/\Lambda)$、$b_0=11h^\vee/(48\pi^2)$，故 $B:=2b_0=11h^\vee/(24\pi^2)$、$t_{\rm up}=24\pi^2/(11\lambda_G)$、**$C_{\rm RG}:=24\pi^2/11$**（配 $\lambda_G$，与 $\dim f$ 无关；物理归一化 $T(f)=1/2$，`handoff-proof`/`group-uniformity` §2.0 单一真源）。$\Lambda$ 为该方案标度，方案常数 $\Lambda\mapsto\Lambda e^{c/2b_0}$ 此约定内一次钉死。$a\downarrow0\iff\beta(a)\to\infty$；次序：先固定盒 $a\to0$，再 $L\to\infty$。

> **两套变量，勿混**：测度段用 $N_c:=\dim f$（$\beta=2N_c/g^2$、$\gamma=2/g^2$）；RT 段用 $h^\vee$（$b_0\propto h^\vee$、$\lambda_G=h^\vee g^2$）。二者仅在 $G=SU(N)$ 时 $h^\vee=N_c$。最小反例：$SO(3)\cong SU(2)/\mathbb Z_2$ 是同一理论、$b_0$ 必同，但矢量表示 $\dim f=3\neq$ 基本表示 $\dim f=2$，而 $h^\vee=2$ 相同。其余 $SO(2n+1)$: $2n-1$ vs $2n+1$；$Sp(n)$: $n+1$ vs $2n$；$G_2$: $4$ vs $7$。

> **盒 ↔ 环面桥（W1$_0$/W1$_1$）**：本文用盒 $\Lambda\Subset\mathbb Z^4$ + DLR 边界 $\xi$；构造与 GL/glue-c 用 $T^4_L$（周期、无边界）。二者经**边界屏蔽 W1$_0$/$W1_1$** 桥接：$T^4_L$ 对应周期边界数据 $\xi_{\rm per}$，而 W1 说 $\gamma_\Delta(A\mid\xi)$ 对 $\xi$ 的依赖 $\le\varepsilon(d)\to0$，故盒结果对所有 $\xi$（含 $\xi_{\rm per}$）一致，从而转移 $T^4_L$ 及反之。绕非平凡圈的拓扑扇区在 $\mathcal A_{\rm gi}$ 上由 Elitzur + W1 处理（残余缺口记此）。glue-c `I-UNIQ` 行即此桥。

**DLR 规范** $\gamma=\{\gamma_\Delta\}_{\Delta\Subset\mathbb Z^4}$：$\gamma_\Delta(\cdot\mid\xi)$ 是给定外部组态 $\xi$ 后 $E(\Delta)$ 上的条件律。Wilson 作用量是**有限程**（每个面 4 条键），故 $\gamma$ 是标准 DLR 规范。$\mathcal G(\gamma)$ 为满足全部 DLR 方程的 Gibbs 测度集。

- **规范群** $\mathcal G_{\rm loc}=G^{V(\mathbb Z^4)}$ 逐点作用 $U_b\mapsto g_xU_bg_y^{-1}$；$\gamma$ 规范不变；物理可观测代数取规范不变量 $\mathcal A_{\rm gi}$（Wilson 环、场强等）。由 Elitzur 定理，非规范不变局部量期望为 0，所以相结构完全由 $\mathcal A_{\rm gi}$ 承载。
- **平移** $T_a$（$a\in\mathbb Z^4$）协变作用于 $G^{E(\mathbb Z^4)}$，$\gamma$ 平移协变。

### 1.2 关键性质（全部以 $\mathcal A_{\rm gi}$ 表述）

以 $\Lambda_0\Subset\Delta\Subset\mathbb Z^4$ 记局部区域，$A\in\mathcal A_{\rm gi}$ 支撑在 $\Lambda_0$，$\|A\|_\infty\le1$，$d:=\operatorname{dist}(\Lambda_0,\partial\Delta)$。

| 记号 | 含义 |
|---|---|
| **OS4** | 给定测度 $\mu$ 的（弱）聚类/遍历性：$Z_\mu[f+T_ag]\to Z_\mu[f]Z_\mu[g]$（$|a|\to\infty$），或时间平均形式。16 输入已对每个子列极限给出。 |
| **W1**（边界屏蔽） | 存在 $\varepsilon(d)\to0$ 使 $\forall\Delta,\ \forall$ 边界 $\xi,\xi'$（在距 $\Lambda_0\ge d$ 之外相同），$\ |\gamma_\Delta(A\mid\xi)-\gamma_\Delta(A\mid\xi')|\le\varepsilon(d)$。速率任意（幂律亦可）。 |
| **MG-spec** | 同 W1，但 $\varepsilon(d)=C e^{-md}$，$m>0$。 |
| **MG-inf** | 单一无限体积态 $\mu$ 的规范不变量指数聚类：$|\langle A;B\rangle_\mu^c|\le C\|A\|\|B\|e^{-md}$。 |
| **UNI** | $\lvert\mathcal G(\gamma)\rvert=1$。 |
| **UNI-lim** | 对任意穷尽 $\Lambda_n\uparrow\mathbb Z^4$ 与任意边界序列 $\xi_n$，$\mu_{\Lambda_n}^{\xi_n}\Rightarrow\mu_\infty$（弱收敛），且 $\mu_\infty$ 与选择无关。 |
| **VAC** | 平移作用在 $L^2(\mu_\infty)$ 上的不动点空间为一维（⇔ $\mu_\infty$ 遍历 ⇔ OS 重构的 Wightman 理论真空向量唯一，up to 相位）。 |

### 1.3 两个层次，务必分开

- **层次 A（格点，$a$ 固定，$\Lambda\uparrow\mathbb Z^4$）**：$G^{E(\mathbb Z^4)}$ 由 Tychonoff 紧，**紧性/Prokhorov 自动**。子列极限存在性容易；唯一性是难点。
- **层次 B（连续，$a\to0$，$\Lambda\uparrow\mathbb R^4$）**：底空间 $\mathcal S'(\mathbb R^4)$ 非紧，**紧性/胎紧性非平凡**，须由一致 OS1 正则界（两点半范数/生成泛函的指数界，$\Lambda$-一致、$a$-一致）提供。这是 16 输入里的分析部分（AN1 一类）。
  连续极限的顺序为"先固定 $\Lambda$ 取 $a\to0$（UV 极限，已在输入内），再 $\Lambda\to\infty$（子列）"。要把后者升级为全网收敛，**最优做法是把 W1 在格点层证明，且对 $a\in(0,a_0]$ 一致**（定理 B），再下沉到连续层；直接在连续层写"连续规范"更难。

---

## 2. 诊断：现有构造为什么只给子列极限

记现有构造产出：固定 $\Lambda$ 的 UV 极限 $\mu_\Lambda^{\rm cont}$；沿子列 $\Lambda_{n_k}$ 的弱极限 $\mu^{(k)}$。每一 $\mu^{(k)}$ 满足 OS0–OS4（含 E4 聚类）并重构出 Wightman 理论。

**只给子列的唯一原因是：没有"极限点唯一"的论证。** 逻辑上：
1. 胎紧性（层 B 的非平凡输入）⇒ 相对紧 ⇒ 存在子列极限；
2. 极限点是 Gibbs 测度（DLR 取极限）；
3. **缺步**：$\mathcal G(\gamma)$ 是否单点？若无，则不同子列可落到不同 Gibbs 态，得到不同理论。

OS4/E4 是**每个极限点自身**的性质，**不能**排除"另一个子列落到另一个 Gibbs 态"。所以 OS4（已在输入内）**不**填补第 3 步。第 3 步需要的是 $\mathcal G$ 的唯一性——等价于 W1（下节 L2）。

> 一句话：**构造给出"至少一个理论"，缺的是"只有一个理论"。** 唯一性 ⟸ W1；而 W1 不在 16 输入的摘要里（16 输入给紧性 + 子列 + 每点 OS4）。

---

## 3. 精确定理

**定理 A（W1 ⇒ 极限唯一 + 真空唯一）。**
设 (W1) 成立，且网 $\{\mu_\Lambda^\xi\}$ 胎紧（层 A 自动；层 B 由一致 OS1）。则 (UNI)、(UNI-lim)、(VAC) 同时成立。

**定理 B（MG-spec ⇒ 一致速率与 $a$-一致）。**
若 (MG-spec) 的常数 $C,m$ 对 $a\in(0,a_0]$ 一致，且 OS1 正则常数对 $a$ 一致，则：(i) 每个 $a$ 的无限体积极限唯一；(ii) 当 $a\to0$ 时 $\mu_\infty^{(a)}$ 收敛到唯一的 OS 测度 $\mu_\infty$；(iii) $\mu_\infty$ 遍历（VAC）且以速率 $m$ 指数聚类（MG-inf）。——MG-spec 是"打包目标"。

**定理 C（普适性——需 RG，独立于 MG）。**
设正则化 $R,R'$ 的裸耦合沿同一物理轨迹匹配（同一物理方案不变量，如 $m/\Lambda$），且**额外的 RG 输入 (RG)**：轨迹处线性化 RG 映射在无关方向上是压缩的（谱半径 $<1$），相关方向可逆匹配，两有效作用量之差为 $O(a^\theta)$（$\theta>0$）的无关算符。则在定理 B 的假设下，$R,R'$ 的 $a\to0$ 极限在 $\mathcal A_{\rm gi}$ 上一致。**MG/MG-spec 本身不蕴含 (RG)。**

---

## 4. 引理链（含证明骨架）

**L1（紧性 ⇒ Gibbs 极限点）。**
层 A：$G^{E(\mathbb Z^4)}$ 紧 ⇒ 任意族胎紧 ⇒ Prokhorov ⇒ 弱极限存在。局部函数连续有界 ⇒ DLR 在极限下保持 ⇒ 任一极限点 $\in\mathcal G$。
层 B：一致 OS1（$\exp(c(\|f\|_1+\|f\|_p^p))$ 界，$\Lambda,a$-一致）⇒ Fernique/Milnos 型胎紧 ⇒ 同上。**这是唯一需要分析输入的地方，已在 16 输入内（AN1）。**

**L2（W1 ⇒ UNI）。** 设 $\mu,\nu\in\mathcal G$，$A$ 支撑于 $\Lambda_0$。取 $\Delta\supset\Lambda_0$ 充分大，距 $\Lambda_0$ 为 $d$。DLR：$\mu(A)=\int\gamma_\Delta(A\mid\xi)\,\mu(\mathrm d\xi)$。由 W1，$\gamma_\Delta(A\mid\xi)$ 与边界数据 $\xi$ 在距 $\Lambda_0$ 为 $d$ 之外的取值**无关**到误差 $\varepsilon(d)$；特别地对所有 $\xi$ 与某参考 $\xi_0$，
$$
\bigl|\gamma_\Delta(A\mid\xi)-\gamma_\Delta(A\mid\xi_0)\bigr|\le\varepsilon(d).
$$
于是 $|\mu(A)-\gamma_\Delta(A\mid\xi_0)|\le\varepsilon(d)$，对 $\nu$ 同式，故 $|\mu(A)-\nu(A)|\le2\varepsilon(d)$；令 $\Delta\uparrow$（$d\to\infty$，$\varepsilon(d)\to0$）：$\mu(A)=\nu(A)$。对全部局部 $A$ ⇒ $\mu=\nu$。∎
（这正是 Dobrushin–Föllmer 唯一性定理：影响矩阵/边界依赖 $\to0$ ⇒ 唯一 Gibbs 态。文献：Dobrushin 1968；Georgii《Gibbs Measures and Phase Transitions》第 8 章。）

**L3（唯一极限点 ⇒ 收敛）。** 相对紧的网若有唯一聚点，则收敛到该点。结合 L1+L2，得定理 A 的 (UNI-lim)。∎

**L4（UNI ⇒ VAC）。** 对平移协变规范，$\mathcal G(\gamma)$ 是 Choquet 单形，其极点 = 遍历（尾 σ-代数平凡）Gibbs 态；Gibbs 测度的遍历分量仍是 Gibbs 测度（Georgii 第 7 章）。$\mathcal G$ 单点 ⇒ 该点即唯一极点 ⇒ 遍历 ⇒ $L^2(\mu)$ 中平移不动子空间 $=\mathbb C\cdot\mathbf1$ ⇒ 真空向量唯一。∎
**要点：L4 不需要指数速率。** 只需唯一性 + 遍历性。

**L5（MG-spec ⇒ W1 ⇒ …，并 ⇒ MG-inf）。** MG-spec 是 W1 的指数特例，故定理 A 全部结论成立。另外对 $\Delta\uparrow$ 取极限把边界依赖界转成无限体积的指数聚类界（MG-inf）：给定 $A,B$，用条件期望做"远距因子化" $|\langle A;B\rangle|\le\sup_\xi|\gamma_\Delta(A\mid\xi)-\mu(A)|\cdot\|B\|$ 得 $Ce^{-md}$。这给出质量隙的速率层（接 M1 的等价性定理）。

**L6（$a$-一致性 ⇒ 连续层唯一）。** 若 L2/L5 的常数对 $a$ 一致，则 $\{\mu_\infty^{(a)}\}$ 的局部边缘被同一组界控制 ⇒ 胎紧 ⇒ 任聚点为 OS 测度；再需 UV 段（固定体积 $a\to0$ 收敛，在输入内）+ 一致性把体积边缘钉住 ⇒ 唯一 $a\to0$ 极限。定理 B。∎

**L7（普适性）。** 对 $R,R'$，由 (RG) 把两有效作用量之差表为相关方向（沿轨迹相消）+ 无关算符（$O(a^\theta)$）。两 Schwindger 函数之差 $\le Ca^\theta$（用稳定性/有界性估计），$a\to0$ 相消。**该步只用 (RG)，不用 MG。** 定理 C。∎

---

## 5. 唯一性的充分条件层级（核心）

```
        MG-spec ──⇒── W1 ──⇒── UNI ──⇒── UNI-lim
          │                              │
          │                              └──⇒── VAC
          └────⇒── MG-inf ──(需 W1/UNI)──⇒── VAC + 速率(=质量隙层)

W1 ──⇒── VAC（经 UNI ⇒ 遍历；不需指数速率）
OS4 ──✗──⇒── UNI          （单态聚类不排除他相）
MG-inf ──✗──⇒── UNI        （同上，Ising 两相）
UNI ──✗──⇒── MG-inf        （唯一但幂律衰减，临界 Ising）
```

> 注：上图以"单一聚类半径"表述；按 N-ality **分扇区**后的精确版本见 §5.5（中性 $+$ 基本两扇区，中心序参量由 OS4 免费）。W1/MG 须在每个扇区分别成立。

两个**尖锐反例**（把层级钉死）：

**例 1（MG-单态 ⇏ UNI）。** $\mathbb Z^2$ Ising，$\beta>\beta_c$。纯相 $\mu_+$（+bc）在体积一致意义下指数聚类，速率 $m(\beta)>0$——即"单一态满足 MG"。但 $\mu_+,\mu_-$ 皆 Gibbs 态；`+bc` 与 `−bc` 给不同无限体积极限。**故单态指数聚类不能钉住极限。**

**例 2（W1/UNI ⇏ MG）。** $\mathbb Z^2$ Ising，$\beta=\beta_c$。Gibbs 态唯一（无自发磁化），边界影响 $\to0$（幂律）⇒ W1、UNI、VAC 均成立；但关联 $\sim|x|^{-1/4}$，无指数衰减、无质量隙。**故唯一性严格弱于质量隙。**

**推论（与 hexagonmath 反例一致）。** 其"E4 聚类但无谱隙"的模型正是一个**真空唯一、无质量隙**的理论：满足 UNI/VAC，破坏 MG。这说明"唯一真空"不是质量隙的进展，恰是另一条轴。

---

## 5.5 扇区分解（N-ality）与 't Hooft 代数：对唯一性论证的修正

> 对接 `research/confinement/REPORT.md` §2.1、§4.5。
> 修正后的结论：**唯一性需要两个扇区，但电荷扇区的"序参量"部分已由 OS4 免费给出**；剩余的电荷输入是"弦扇区唯一"，属禁闭任务的对象——但只需其**唯一性**，不需 $\sigma>0$。

### 5.5.1 分次

规范不变代数 $\mathcal A_{\rm gi}$ 按中心荷（N-ality $n\in\mathbb Z_N$）分次：$\mathcal A_{\rm gi}=\bigoplus_n\mathcal A_n$。
- $\mathcal A_0$：**中性**（局域胶球 $\mathrm{tr}\,F^2$、胶球算子等）；
- $\mathcal A_{\pm1}$：**基本**（Wilson 环、Polyakov 环——非局域但规范不变）；
- 中心（全局）作用：$A\in\mathcal A_n\mapsto\omega^n A$，$\omega=e^{2\pi i/N}$。周期/自由边界下有限体积测度**中心不变**。

### 5.5.2 事实 1：OS4（全强度）⇒ 混合 (MIX) ⇒ 中心不破（电荷序参量已由输入给出）

**假设 (MIX)**：极限测度 $\mu_\infty$ 在时空平移下**混合**，即对一切有界可测 $F,G$，
$$
\int F\,(G\circ T_a)\,d\mu_\infty \xrightarrow{\ |a|\to\infty\ }\int F\,d\mu_\infty\int G\,d\mu_\infty .
$$
(MIX) 蕴含**任意**有界可测观测量的聚类——**含非局域 Wilson/Polyakov 环**。

**(MIX) 由 16 输入的 OS4 给出（需全强度）。** `OS4_Clustering` 的定义是对**一切** $f,g\in\mathcal S$：$Z[f+T_ag]\to Z[f]Z[g]$。而 $Z[f+T_ag]=\mathbb E[e^{i(\langle\omega,f\rangle+\langle\omega,T_ag\rangle)}]$ 正是 $(\langle\omega,f\rangle,\langle\omega,T_ag\rangle)$ 的**联合特征函数**。
1. 对一切 $f,g$：联合特征函数 $\to$ 边缘之积；由 OS2（平移不变）边缘即 $Z[f],Z[g]$。
2. Lévy 连续性（$\mathbb R^2$）⇒ 联合律 $\to$ 积律。
3. Cramér–Wold（对 $\sum c_if_i,\ \sum d_jg_j$ 的一切线性组合）⇒ 任意有限柱集对的联合律 $\to$ 积律；Portmanteau ⇒ 柱集混合。
4. 柱集是生成 $\sigma$-代数的 $\pi$-系 ⇒（$\pi$-$\lambda$）全可测集混合 ⇒ (MIX)。
故**若构造的 OS4 是全强度（对一切 $f,g$、一切方向 $a$）**，则 (MIX) **已成立**，§5.5.2 不超出假设。**若 OS4 仅对受限类（如只局域算子、或只时间方向）成立**，则须把 (MIX) 列为**显式独立假设**（等价于转移矩阵/谱聚类）。为对接 gluing，下文契约**一律显式列入 (MIX)**。

**由 (MIX) 到中心不破。** 设 (MIX) 且 bc 中心不变。若中心自发破缺，则有 $N$ 个真空 $\mu_k$（Polyakov 环取 $|P|\omega^k$）；中心对称平均 $\nu=\tfrac1N\sum_k\mu_k$ 中心不变，但
$$
\langle P^*P(x)\rangle_\nu\to|P|^2>0,\qquad \langle P^*\rangle_\nu\langle P\rangle_\nu=0\;(\text{因 }\langle P\rangle_\nu=0),
$$
与 (MIX) 矛盾（$P$ 有界可测；取 $F=P^*,\ G=P$）。故 (MIX) 排除中心破缺 ⇒ 中心不破。

- **修正 §8 障碍 1**：退禁闭是**有限温度**（$L_t$ 有限）现象；$T=0$ 极限下 (MIX) **已强制**中心不破。取 $L_t\to\infty$ 即可，无需额外假设。
- **前提**：(i) 构造用的 bc 中心不变（周期/自由均可；中心扭曲 bc 须另查）；(ii) 连续层 $P$ 为**重整化后** $L^2$ 量（$P$ 本身有周长发散，聚类作用于重整化关联；格点层 $|P|\le N$ 有界，无碍）。

### 5.5.3 事实 2：'t Hooft 代数给出二择一与 $|Z_N|$ 多重度

$$
\langle W(C)V(C')\rangle=\omega^{\,\mathrm{link}(C,C')}\langle V(C')W(C)\rangle .
$$
- ⇒ $W$ 与 $V$ 的**面积律不能同时**成立（至多一个面积、另一个周长）。
- ⇒ 中心破缺时，极值物理 Gibbs 态数 $\ge N$（'t Hooft 环在 $N$ 个真空上置换）。

由此得组织关系（'t Hooft 代数**不证明**任一因子，只把它们组织起来）：
$$
\boxed{\text{物理唯一性}\iff \text{中心不破}\;\wedge\;\text{中性扇区唯一}\;\wedge\;\text{弦扇区（}\mathcal A_{\pm1}\text{关联）唯一}.}
$$

### 5.5.4 事实 3：中性唯一是否蕴含弦扇区唯一？

- **(a) 未必，需 W1_1**：逻辑上允许两个中心不破纯相在中性上一致而在弦张力 $\sigma$ 上不同。4D $SU(N)$ 无已知此类共存，但未排除。
- **(b) 可能由 Wightman 侧蕴含（可省掉电荷输入）**：局域规范不变代数 $\mathcal A_0$ 对真空循环（Reeh–Schlieder/Borchers），Wilson 环为非局域量、疑由 $\mathcal A_0$ 经 OPE 决定；若成立，则中性唯一 ⇒ 弦扇区含 $\sigma$ 唯一。物理上自洽（$\sigma=c_\sigma\Lambda^2$，中性物理钉住 $\Lambda$ 即钉住 $\sigma$）。**格点 Gibbs 框架中未严格化。**

### 5.5.5 对 §5/§6 的修正

- **W1 须按扇区分层**：W1$_0$（中性）+ W1$_1$（基本）。§1.2 的"W1"应读作二者之合。
- **MG 的可观测类至关重要**：若 MG（单态）只含局域中性算子（典型），则 MG ⇒ 至多 W1$_0$，**不含电荷扇区**。故"MG-spec ⇒ 全部"仅在 MG 对**所有规范不变观测（含 Wilson/Polyakov 环）**一致量化时成立——**待核实其可观测类**（与"是否对 bc 一致"并列的第二项核对）。
- **修正后的唯一性充分条件**：
$$
\text{唯一性}\;\Longleftarrow\;\text{W1}_0\;\wedge\;(\text{W1}_1\ \text{或 (b)}),
$$
  中心不破因子由 **OS4 + 中心不变 bc 免费**给出（事实 1）。
- **唯一性仍严格弱于禁闭**：只需"弦扇区唯一"，**不需 $\sigma>0$**。故唯一性弱于质量隙（例 2）**且**弱于禁闭（$\sigma>0$）。
- **与禁闭任务的耦合点**：唯一性的剩余电荷输入 = "弦扇区唯一"，其证明很可能与禁闭/面积律共用同一 RG 墙，但**目标更弱**（唯一性，而非张力正性）。

---

## 5.6 可消费的 typed 契约（供 gluing 直接使用）

记号：$G$ 紧单；$\Lambda_n\uparrow\mathbb Z^4$ 穷尽；格距 $a\in(0,a_0]$；$\gamma$ 为 Wilson 模型的 DLR 规范（有限程）；$\mathcal A_{\rm gi}=\bigoplus_{n\in\mathbb Z_N}\mathcal A_n$ 按 N-ality 分次；$\Lambda_0\Subset\Delta$，$d=\operatorname{dist}(\Lambda_0,\partial\Delta)$，$A$ 支撑于 $\Lambda_0$，$\|A\|_\infty\le1$；$\xi,\xi'$ 为在距 $\Lambda_0\ge d$ 之外相同的边界数据。

```text
Contract UniquenessInput

  # --- 边界假设（显式）---
  BC : 构造所用边界条件中心不变（周期或自由）。
       若为中心扭曲 bc，则中心不破结论失效，须另证。

  # --- 混合（电荷扇区的关键；由全强度 OS4 蕴含，否则独立假设）---
  MIX : 极限测度 μ∞ 在时空平移下混合，即
        ∀ F,G ∈ L∞(μ∞):  ∫ F·(G∘T_a) dμ∞ → ∫F dμ∞ · ∫G dμ∞  (|a|→∞).
        蕴含：一切有界可测观测量聚类，含非局域 Wilson/Polyakov 环。
        来源：OS4_Clustering(∀ f,g, ∀ 方向) ⇒ MIX（§5.5.2 的四步论证）。
        若 OS4 仅受限，则 MIX 为独立假设（=转移矩阵/谱聚类）。

  # --- 中性扇区边界屏蔽 ---
  W1_0 : ∀ Λ0⋐Δ, ∀ A ∈ 𝒜0 支于 Λ0, ‖A‖∞≤1, ∀ ξ,ξ' :
           |γΔ(A|ξ) − γΔ(A|ξ')| ≤ ε0(d),   ε0(d) → 0.   (速率任意)

  # --- 基本扇区边界屏蔽（弦扇区）---
  W1_1 : 同上，A ∈ 𝒜±1（有界的重整化 Wilson/Polyakov 环）。
         注：若接受 §5.5.4(b)（局域中性代数经 Reeh–Schlieder/OPE 决定非局域环），
             可省 W1_1，但格点侧未严格化。

  # --- 紧性（层 A 自动；层 B 由一致 OS1）---
  TIGHT : {μΛ^ξ : Λ,a,ξ} 胎紧（⇐ 一致 OS1 正则界）。

  # --- 待核实的 MG 可观测类（不由本会话审计）---
  MG-class : hexagonmath MG 是否 (i) 对 bc 一致；(ii) 量化 Wilson/Polyakov 环。
             默认（保守）：若只含局域中性算子，则 MG 只服务 W1_0，不含 W1_1。

  # --- 结论（定理 A/B 的扇区版）---
  MIX ∧ BC                         ⇒ 每个子列极限中心不破。
  TIGHT ∧ W1_0 ∧ W1_1              ⇒ UNI 与 UNI-lim（唯一物理 Gibbs 态 + 全网收敛）。
  UNI（+ 平移不变）                 ⇒ VAC（唯一真空）。
  W1_0 单独                         ⇒ 仅中性扇区唯一（不足）。
  W1_0 ∧ (W1_1 ∨ §5.5.4(b)) ∧ MIX  ⇒ 全物理唯一性（中心因子免费）。

  # --- 明确不在本契约内 ---
  σ>0（禁闭）；质量隙速率（MG）；跨正则化普适性（需 RG，(RG) 假设）。
```

**给 gluing 的一句话**：本契约把"唯一性"分解为 **(MIX) + BC + W1_0 + W1_1 + TIGHT**；其中 (MIX) 若 OS4 全强度则免费，BC 为构造事实，TIGHT 属 16 输入，**真正的新增数学输入是 W1_0 与 W1_1**（弦扇区）。MG 只（至多）覆盖 W1_0，故 **"MG for Wilson/Polyakov" 应作为独立需求**，与本会话未审计的 (MG-class) 一致。

---

## 6. 超出 MG 到底还需要什么（逐项）

现有（16 输入摘要）：胎紧 + 子列极限 + 每极限 OS0–OS4。在此基础上（按扇区，见 §5.5）：

| 目标 | 所需**额外**输入 | 是否由 MG-单态推出 |
|---|---|---|
| **UNI-lim**（子列 → 全网收敛） | **W1$_0$ ∧ W1$_1$**（中性 + 基本扇区边界屏蔽，速率任意）；中心不破因子由 OS4 免费 | **否**（例 1；且扇区见 §5.5） |
| **UNI**（唯一物理 Gibbs 态） | W1$_0$ ∧（W1$_1$ 或 §5.5.4(b)）；或独立的"中性唯一 ∧ 无弦共存" | 否 |
| **VAC**（每个极限真空唯一） | 已由 OS4/遍历性给出（核实 OS4 在 $\mathcal A_{\rm gi}$ 上成立）；跨极限"同一真空"回到 UNI-lim | 部分 |
| **MG-inf 速率 / 质量隙** | MG 本身（主攻方向，另一会话） | — |
| **MG-spec（打包）** | W1$_0$ 的指数版（+ $a$-一致）；若 MG 仅含局域中性算子，则**不含** W1$_1$ | 是（若已对 bc **且**对 Wilson/Polyakov 一致量化，待核对原文） |

**结论**
- **对唯一性，最小充分输入是 W1$_0$ ∧ W1$_1$，而非 MG。** 每扇区 W1 逻辑上弱于该扇区的 MG-spec（MG-spec ⇒ W1），且 W1 不蕴含 MG（例 2）。中心不破因子由 OS4 免费（§5.5.2）。
- **MG 若只以"单一态/极限态"且"仅局域中性算子"形式陈述，则不足以给唯一性**（例 1 + §5.5）。主会话在攻 MG 时应**同时**证明：对 bc 一致（MG-spec）、**且**对 Wilson/Polyakov 等电荷算符成立，否则拿不到唯一性；或另证"无相共存 / 弦扇区唯一"。**这是本报告最可操作的一条。**
- **逻辑上 W1 弱于 MG，但证明上未必更易**：4D 非阿贝尔 YM 无 FKG/关联不等式（非阿贝尔无正关联结构），W1 很可能仍需同一套 RG 技术。目前**无已知捷径**。

---

## 7. 普适性：独立于 MG 的 RG 输入

分两类，别混：

- **普适性 (i)：同一正则化内，与穷尽顺序/边界条件无关。** 由 W1（定理 A）给出。属于"唯一性"范畴。
- **普适性 (ii)：与正则化（格距 $a$、作用量 $S^a$、blocking/RG 方案）无关。** 需 (RG)：
  - **相关/边缘算符匹配**：两方案的裸耦合沿同一物理轨迹（同一 $\Lambda_{\rm QCD}$、同一 $m/\Lambda$、同一 $\sigma/\Lambda^2$ 等方案不变量）。方案变换 $\beta'=\beta+cg^4+\dots$ 只平移 $\beta$ 的地址，不改物理；见 `research/physics_handoff.md` 的 scheme-无关量讨论。
  - **无关算符控制**：$S-S'=O(a^\theta)\times$（无关算符），且 RG 线性化在无关方向压缩。此即 Balaban 型 UV 稳定性 + 流动可微性的定量版，部分或已含于 RG1–RG5/UV1，但"跨正则化一致"作为**独立命题**须显式陈述与证明。
  - 结论：$\mu_\infty^{R}=\mu_\infty^{R'}$ 在 $\mathcal A_{\rm gi}$ 上。
- **MG 在此不出现。** 质量隙与普适性是两条独立轴：前者是 IR 速率，后者是 UV/中尺度算符匹配。

---

## 8. 障碍清单

1. **相共存 / 退禁闭——已大部分被 OS4 消解，残余在弦扇区。** 若 bc 中心不变且取 $T=0$（$L_t\to\infty$），OS4 已强制中心不破（§5.5.2），故"退禁闭型 $N$ 重真空"**不是**障碍。**残余障碍**：两个中心不破纯相可在中性一致而在弦张力 $\sigma$ 上不同（§5.5.4(a)），需"弦扇区唯一"（W1$_1$ 或 §5.5.4(b)）。**注意前提**：构造必须用中心不变 bc；若用中心扭曲 bc，则须另证其极限中心不破。
2. **规范商。** 唯一性只能对**规范轨道**（或 $\mathcal A_{\rm gi}$）陈述；链接测度本身在规范群下不是"唯一"的对象（规范拷贝被 Haar 积分掉）。Elitzur 保证只有规范不变量承载相结构，故这不损失内容，但**陈述必须落在 $\mathcal A_{\rm gi}$**。
3. **$a$-一致性（层 B）。** W1/MG-spec 须对格距 $a$ 一致，否则只能得到"每个 $a$ 唯一"，跨 $a$ 仍可能分叉。这是相对层 A 的**额外**技术负担。
4. **连续层规范。** 连续层没有逐点 DLR；要么走"格点一致 W1 + UV 段"的路线（推荐），要么建立连续规范（更难）。
5. **无关联不等式。** 非阿贝尔 YM 无 FKG/正关联，排除"单调性"捷径；W1 只能靠 RG/展开。
6. **W1 与 MG 的证明成本相近。** 逻辑上 W1 弱，但现有工具都指向同一瓶颈（弱耦合 RG 交接）；不能指望"只证 W1 就绕开 MG 难题"。

---

## 9. Lean 形式化计划（本目录内，不重构建）

现有 `YMStatement/Basic.lean` 是**标量场**版（`FieldConfiguration d`）。YM 需换成**规范场组态 + 规范不变可观测代数**；"真空"落在物理（规范不变）Hilbert 空间。以下按可落地性排序：

1. **`exponentialClustering_implies_ergodic`**（高价值、可做）。
   从 `ExponentialClustering μ m`（生成泛函指数聚类）⇒ 平移作用在 $L^2(\mu)$ 上**混合/遍历**。技术点：指数聚类代数稠密 + Mathlib `Mathlib.Dynamics.Ergodic.*` 的均值遍历定理/遍历等价刻画。
2. **`unique_vacuum_of_ergodic`**（可做）。
   遍历 ⇒ 平移不动子空间 $=\mathbb C\cdot\mathbf1$ ⇒ 真空向量一维。本质是遍历性的不动集刻画。
   合并 1+2 ⇒ **定理 A 的 (VAC) 部分**，在概率层；不依赖 OS 重构的 GNS（当前 Mathlib 无 GNS）。
3. **抽象 DLR 唯一性 `Specification.unique_of_boundaryDecay`**（中）。
   在一般可测乘积空间上把 L2 形式化：边界依赖 $\to0$ ⇒ 唯一 Gibbs 态。避开 YM 具体作用量，只留 $(\varepsilon(d))_{d}$ 的抽象条件。
4. **`unique_clusterPoint_implies_tendsto`**（中，取决于 Mathlib Prokhorov）。
   紧空间 + 唯一聚点 ⇒ 收敛。查 `MeasureTheory` 是否有 Prokhorov/紧性支撑；无则用有限维/乘积紧 + Tychonoff 手搭。
5. **（后续）MG-spec ⇒ MG-inf**：把边界依赖界转成无限体积聚类界（L5），接 M1 的"指数聚类 ⇔ 谱隙"。

**注意**：以上 1–2 是"真空唯一"的形式化，**不是**"极限唯一"。真正的极限唯一（定理 A 的 (UNI-lim)）依赖 3–4，且其**数学输入（W1）本身未证**——形式化只能固化"W1 ⇒ 唯一"的**蕴含**，不能替代 W1 的证明。

---

## 10. 下一步（按性价比）

1. **【立即·核对】确认 hexagonmath MG 的量化形式**：是"单一极限态"还是"对边界条件/体积一致"？若为前者，主会话应把目标改为 **MG-spec**（或另立"无相共存"）。本报告 §5 例 1 表明这一步不可省。
2. **【近期·形式化】VAC 的概率层**（§9 的 1–2）：给出"指数聚类 ⇒ 真空一维"的 0-sorry 定理，接 `YMStatement` 陈述层。这是"真空唯一"的可交付件。
3. **【中期·数学】W1 的弱耦合路线**：把"边界影响衰减"作为独立引理，尝试用 RG 步进 + 大场估计证明，且对 $a$ 一致；与主会话的 RG 交接共享技术，但**目标更弱**（不需速率）。列为比 MG 更可达的中间里程碑。
4. **【中期·数学】普适性 (ii)**：显式陈述 (RG) 假设与"跨正则化一致"定理；核对是否已在 RG1–RG5/UV1 中隐含。
5. **【长期】W1/MG-spec 的 4D 正面证明**：即弱→强 RG 交接的本体，与主会话合流。

---

## 附录 A：关键文献（按本文用途）

- **DLR / Gibbs 测度 / 唯一性**：Dobrushin, *Theory Probab. Appl.* 13 (1968) 197–224；Georgii, *Gibbs Measures and Phase Transitions*, 2nd ed., de Gruyter (2011), 第 7–8 章（Choquet 单形、唯一性定理）。
- **格点规范唯一性/质量隙（强耦合）**：Osterwalder–Seiler, *Ann. Physics* 110 (1978) 440–471；Seiler, *Gauge Theories as a Problem of CQFT and Stat. Mech.*, LNP 159 (1982)。
- **Dobrushin 在规范理论**：Chatterjee, *Gauge Theory*（书/预印本），第 15–16 章（Dobrushin 天花板：仅强耦合）。
- **Langevin / Bakry–Émery 强耦合唯一性 + 质量隙**：Shen–Zhu–Zhu, arXiv:2204.12737，*CMP* 400 (2023) 805–851。
- **OS 重构与"聚类 ⇒ 唯一真空"**：Glimm–Jaffe, *Quantum Physics* (1987)；Simon, *The $P(\varphi)_2$ Euclidean (Quantum) Field Theory* (1974)。
- **Elitzur 定理**：Elitzur, *Phys. Rev. D* 12 (1975) 3978。
- **中心对称/退禁闭**：Polyakov (1978)；Svetitsky–Yaffe, *Nucl. Phys. B* 210 (1982) 423。
- **Prokhorov**：Prokhorov, *Theory Probab. Appl.* 1 (1956) 157–214。
- **临界 Ising（反例 2）**：标准结果（如 McCoy–Wu；Georgii 书）；$\beta_c$ 唯一 Gibbs 态、幂律关联。
- **16 输入与 MG**：hexagonmath, 2611.00492v1 (2026-11-15)；摘要见 `docs/yang-mills-prior-work.md` §3.4（**具体量化形式待核对**）。

> 诚实说明：本文未重读 2611.00492 原文，对 16 输入的描述取自项目摘要；定理 A/B/C 的假设是**分析性**的（W1/MG-spec/(RG)），其证明不属本报告。反例 1、2 是标准事实，用以确定层级；它们说明的是"单态 MG 不足"与"唯一性弱于质量隙"，不构成对 YM 的任何断言。
