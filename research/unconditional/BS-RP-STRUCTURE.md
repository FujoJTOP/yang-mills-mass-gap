# BS-RP-STRUCTURE — 两胶子复合体的反射正性能否给 $0^{++}$ 质量一个**下**界？

> 2026-11-23 · 4D 纯 SU(N) YM · 角度：$0^{++}$ = J=0 两胶子 BS 束缚态，试探"复合结构 $+$ 反射正性 (RP)"是否给出比单模 KL 正性**更强**的东西。
> **判词 (Verdict)：NOT 5/5 — 2/5。** 复合结构确有单模 KL 没有的三条精确事实：(i) 作用量密度 $O=-\operatorname{tr}F^2\ge0$ **逐点非负**；(ii) 阈值恒等式 $\inf\operatorname{supp}\nu_{F^2}=(\sqrt{s_1}+\sqrt{s_2})^2\big|_{s_i=\inf}$ $=4\inf\operatorname{supp}\rho_\phi=(2m_g)^2$；(iii) 连续谱的复合仍为连续谱（**复合体不制造原子**）。但这三条**全部是上界或条件式**：阈值恒等式把 $0^{++}$ 卡在 $(2m_g)^2=$ 两胶子阈值，方向给出 $m_{0^{++}}\le 2m_g$（**上界**，反方向），其前提"两胶子有正阈值"正是 $m_g>0=$ H1；逐点非负性被**自由麦克斯韦完整饱和**（$O\ge0$、$\nu\ge0$、支撑自 $0$ 起、无原子、却无胶球）；连续谱⟶连续谱说明束缚态原子必须来自**另行的吸引输入**（BS-KERNEL-SIGN 的符号不定核，$=$ H1/$f_1^2>0$）。故复合 RP 未越过"正性只给上界"这堵墙。精确卡点：**两胶子阈值 $=(2m_g)^2$ 中的 $m_g>0$（正范数扇区的胶子质量隙 $=$ H1）**。

## 1 · 对象：复合体的 RP 究竟是什么

令 $O(x)=-\operatorname{tr}F_{\mu\nu}F_{\mu\nu}(x)$（反厄米生成元，Euclidean）。RP/Osterwalder–Seiler 只对**规范不变**的 $O$ 成立，给出 KL 测度
$$
W_c(t)=\langle O(t)O(0)\rangle_c=\int_0^\infty e^{-st}\,\nu(ds),\qquad \nu\ge0,\qquad m_{0^{++}}^2=\inf\operatorname{supp}\nu .
$$
这一步 $=$ BS-BLIND §(c) 的单测度 KL 正性：只约束 $W_c$ 为完全单调（对 $t$ 的各阶导交替、上界），**允许 $\nu$ 在 $s=0$ 堆积**。任务问的是：$O$ 是**两胶子复合**（$O\sim F\cdot F$，$1^{-{-}}\otimes1^{-{-}}\supset0^{++}$），这个复合/两粒子结构能否补上单测度缺的**下**界？下面认真钢人化（steelman），再攻击。

## 2 · Steelman：复合特有的三条精确事实

**(i) 逐点非负。** $O=-\operatorname{tr}F^2\ge0$ 作为函数（Euclidean 作用量密度）。故 $O$ 落于不变函数正锥内，$\langle O\rangle\ge0$（凝聚为正），$W_c$ 是逐点 $\ge0$ 的正定函数。这是复合体最接近"正性改良"的性质。
**(ii) 阈值恒等式（Källén）。** Gaussian（自由）复合 $:\!\phi^2\!:$ 的连通谱测度是两体相空间
$$
\nu_{F^2}(s)=\!\int\!\!\!\int\rho(ds_1)\rho(ds_2)\,\Pi_0(s;s_1,s_2),\quad \operatorname{supp}\subseteq\big[(\sqrt{s_1}+\sqrt{s_2})^2,\infty\big),
$$
$\Pi_0\propto\sqrt{\lambda(s,s_1,s_2)}/s$ 为 s-波相空间。于是 $\inf\operatorname{supp}\nu_{F^2}=(\!2\sqrt{\inf\operatorname{supp}\rho})^2=(2m_g)^2$：复合底**恰为两胶子阈值**，且**无原子**（相空间在阈处 $\sim\sqrt{s-4m_g^2}\to0$）。**但**这是 BS-THRESHOLD 的 (T1)——**微扰/OPE 两胶子切**，仅对 Gaussian 精确；禁闭相 (T2) 里胶子非渐近态，精确测度**无 (T1) 切**，物理连续谱是**多胶球**，阈 $\theta_c\ge(2m_1)^2\ge(2m_g)^2$。两阈皆条件于 $m_g>0$，差别只是"把哪个 $m$ 当构成质量"，不改下界。
**(iii) 复合不造原子。** 连续（无原子）构成测度的复合仍为连续谱；原子（束缚态）只能来自构成测度本身有原子或核有吸引——复合结构自身不产生 $0^{++}$ 原子。
**(iv) OPE。** $O(x)O(0)$ 的短距展开给出受保护系数（单位算子、$T_{\mu\nu}$、迹反常 $c>0$）——但见 GAP-G2-edge：这是**总/UV**数据，只把总量变**上界**，dim-4 无有限受保护矩。

## 3 · Attack：三条全部上界/条件，无一给下界

**(a) 方向错。** 阈值恒等式说 $0^{++}$ 位于 $(2m_g)^2$，但稳定的胶球应是 $m_B^2<(2m_g)^2$ 的**极**（低于阈值），"低于阈值"是绑能，**反方向**于下界。能从严证明的只有 $m_{0^{++}}\le2m_g$（**上界**），对间隙无用——与 BS-KERNEL-SIGN §4、G2-edge"变分给上界"一致。更糟：任何变分/Temple 方案的"下一能级" $E_2$ 其实是**径向激发** $0^{++*}$（BS-VARIATIONAL-LOWER），一个动力学 BS 束缚态而非运动学两体阈，故连 $E_2$ 的锐度也是动力学的、不可用。
**(b) 循环。** 阈值 $(2m_g)^2$ 需要 $m_g>0$，即**胶子**在正范数扇区的质量隙。这正是 H1（$G''>0$/$f_1^2>0$），是间隙本身；K0-TRANSVERSE-DECOUPLING 的负模、Savvidy 不定背景让它非微扰不可得。复合体把缺口**搬**到"构成胶子质量隙"，未移除。
**(c) 构成分解本身不 RP-合法。** 严格说，$O=\operatorname{tr}F^2$ 的两胶子分解是**规范固定/构成图像**的启发式；胶子被禁闭、不在物理 Hilbert 空间，其传播子无 KL 正性（Oehme–Zimmermann 变号，$\int\rho_A=0$）。故"两胶子阈值"在 YM 里**不是一个严格 RP 对象**。诚实二难（正是 BS-THRESHOLD 的 (T1)/(T2) 二分）：**若**承认构成分解 ⟹ (T1) 阈值 $=(2m_g)^2$，但 $m_g>0=$ H1 且方向为上界；**若**不承认 ⟹ 精确测度 (T2) 无两胶子切、连续谱为多胶球 $\theta_c\ge(2m_1)^2$，**唯一的严格复合 RP 就是 §1 的单测度 KL 正性**（BS-BLIND(c)，上界、允许 $\delta=0$），且 $\nu((0,m_1^2))=0$ **字面上就是** Clay 的支撑陈述 (A)。"禁闭抬升"是把 $m_g$ 改名为 $m_1$ 的**重贴标签**，非动力学抬升。两读皆无下界。
**(d) 逐点非负被自由场饱和（酸测）。** 自由麦克斯韦：$O\ge0$ 逐点、$\nu\ge0$、$W_c>0\ \forall x$、阈值 $0$、无原子——**完整满足** (i)(ii)(iii) 而**没有任何胶球**、谱自 $0$ 起。故逐点非负 $+$ 复合 $+$ RP $+$ 正性改良，全被无质量情形饱和；它们**无法**判别 gapped 与 gapless。

## 4 · 数值酸测（numpy 本地复核，全部通过）

- 自由有质量 $:\!\phi^2\!:$：阈 $(2m)^2=4$，密度 $\sim(s-4)^{1/2}$（log-log 斜率 $p=0.4999$，**无 $\delta$**）。
- 自由无质量：支撑自 $0$ 起，$\rho\equiv1>0$，RP 成立而 gap $=0$。
- 阈值 $(\sqrt{s_1}+\sqrt{s_2})^2$：$(m,m)\to4$、$(m,0)\to1$、$(0,0)\to0$——**有质量构成$\Rightarrow$阈值 $2m$；无质量构成$\Rightarrow$复合无隙**；连续$\to$连续（无原子）。
- 矩退化：$\delta_0+\delta_2$ 与 $2\delta_1$ 共享 $M_0=M_1=2$，$\inf\operatorname{supp}$ 为 $0$ vs $1$——前两矩**不能**卡边（GAP-G2-edge）。

## 5 · 记分与最锐一句

1. **非微扰 — PARTIAL**（RP/KL 非微扰；阈值恒等式仅在 Gaussian 精确；逐点非负精确）。
2. **符号/下界方向 — FAIL**（阈值给 $m\le2m_g$ 上界；KL/PF 皆上界；RP 是 PSD⟹Cauchy–Schwarz⟹上界）。
3. **非循环 — FAIL**（需 $m_g>0$=H1，或原子 $=f_1^2>0$ 绑能）。
4. **非 SUSY$+$自由场失败 — 部分/空洞**（自由麦克斯韦完整满足结构却无隙；判别力被完全外包给 $m_g$）。
5. **一致 $a\to0$、$\Delta\ge\lambda$ — FAIL**（同墙）。

**2/5；Clay 0%。** 真东西：Källén 阈值恒等式、逐点非负的酸测、连续谱不造原子引理。**非真**：任何 $m_{0^{++}}>0$ 的下界。

**最锐一句。** "复合体 $+$ RP"比单模 KL 只多三样，且每一样都朝**错方向**：阈值恒等式把质量**上**卡在 $2m_g$（条件式）；逐点非负被自由麦克斯韦饱和（判别力为零）；连续$\to$连续明示原子（束缚态）须另求（BS 核的吸引，符号不定）。单模 RP 是"上界/PSD"，复合 RP **仍是**"上界/PSD"——两粒子结构不改 RP 的**单边性**。

**精确卡点。** 需要**严格强于单测度 KL 正性**的两胶子 RP 不等式；但任何 RP 不等式是 PSD/Cauchy–Schwarz（上界），而唯一的严格复合内容（Källén 阈值 $(2m_g)^2$）条件于**构成胶子质量隙 $m_g>0$（$=$ H1）**——正是间隙本身。复合 RP 未越过墙；缺口在 BS 坐标下重述为 "$m_g>0$ 的正范数扇区下界 $+$ 核吸引给出 $0<m_B<2m_g$"，与 BS-KERNEL-SIGN / GAP-G2-edge / SUPPORT-NOT-ATOM 的卡点同一。

*Refs.* Källén, Helv. Phys. Acta **25** (1952) 417; Lehmann, Nuovo Cim. **11** (1954) 342; Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; Osterwalder–Schrader, CMP **31** (1973) 83; Salpeter–Bethe, PR **84** (1951) 1232; Maris–Roberts–Tandy, nucl-th/0009064; Oehme–Zimmermann, PRD **21** (1980) 1661（胶子谱变号 $\int\rho_A=0$）。Internal: BS-BLIND §(c)（单测度 KL、$\delta=0$ 允许）; BS-KERNEL-SIGN（核符号不定、阈在 $0$、变分上界）; BS-THRESHOLD（(T1)/(T2) 两阈、禁闭无两胶子切、多胶球阈 $\theta_c\ge(2m_1)^2$、抬升 $=$ 重贴标签）; BS-VARIATIONAL-LOWER（$E_2=$ 径向激发 $0^{++*}$、非运动学阈）; BS-NUMERICAL-RIGOROUS（RP 只上界、下界需 Simon–Spencer 包络）; K0-BLIND（KL 倒数 $=$ 正密度、负裸模不入）; GAP-G2-edge（矩只上界）; SUPPORT-NOT-ATOM（唯一右号仪器 $=G''=$ H1）; PROOF-GAUSS-INVARIANT-PF（PF 钉扇区顶 $=$ 真空）; K0-TRANSVERSE-DECOUPLING（负模与立方顶点）; FORGE-FESHBACH-ANGLE（$C\prec0\Rightarrow$ 只加正、$=H1$）.
