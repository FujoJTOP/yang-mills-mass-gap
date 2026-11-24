# SURVEY-MASTER-SYNTHESIS — 全网六路扫描总表（2026-11-24）

> 任务：全网找论文，"试图找到哪怕一种可能性"。六路独立扫描 + 一次 250K 声明的 Lean 验证。
> **一句话判决：零条路线跨过交叉点。唯一的正面里程碑是 G1（存在）已被机器验证——而它自己明说把隙排除在外。G2（隙）= 整个剩余问题，每条 G2 路线都归约到 H1 / 边界一致 LSI / 交叉点无小参数。Clay(G2) = 0%。**

---

## 0. 总表（六路）

| 路 | 域 | 结论 | 关键机制 / 最接近者 |
|---|---|---|---|
| **野路子数学** | 代数拓扑/等变/K理论、表示论/spin foam、Connes NCG、算术/trace-formula | **全 → A3/B6/刚性/通道盲** | 算术(Weil/Connes 迹公式)= 我们自己的 **A3**(Fredholm `det(1−zT_O)`)+**B6**,证 Weil≥0 = RH 循环;**A3 是吸引子**;指标/解析挠率=刚性(B3);范畴对称对 `0⁺⁺` 中性;NCG 谱作用自由场 vacuous |
| **外来工具** | 熵 Schrödinger 桥/信息几何/平均场博弈/Ruelle/Weyl/hypocoercivity/Stein | **全 → H1** | 熵桥 = `GIRSANOV∧E2` 的动态面(sup⇒lsc vs 原子 usc;每尺度 `η_k<1`=H1);信息几何=Γ₂=H1;MFG=McKean–Vlasov 位移凸=H1;transport-entropy 家族**动静态双封** |
| **最新 arXiv 声称** | 2024–2026 的"YM 质量隙"声称/新技巧 | **全审死/无新路** | Lim 2307.00788(**"nowhere used the YM action"**、手设 Casimir H、酸测败);Platek Zenodo 20114011("**never enters strong coupling**"=跳过交叉);Mihajlovic 2611.11534(诚实,缺=H1);Bao 2608.02452(码字典,无把手) |
| **AI/形式化** | Lean/Coq、AI 生成证明、区间算术 | **G1 机器验证(G1-only!),无 G2** | **OpenAI Lean G1**(repo `yang-mills-group/OpenAI-YM-20261028`,Gross 委员会,9,470 模块/250,319 声明,无 `sorryAx`,Comparator 过)——**明确"does not prove there is a mass gap"**;Faizal–Shabir 2606.19362 条件、死于交叉点;Farooqui Coq 定义性作弊、jonwashburn Lean `build` 失败 |
| **构造/严格 QFT** | 随机量子化、正则结构、LSI/谱、簇展开 | **两个不相交 regime 类,从无交集** | 强耦+连续 SU(N)+全体积(OS78 / **SZZ `β<1/48`** / Nissim U(N))—`β≈2.5` Γ₂负+Dobrushin 失效;弱耦+4D+**有限群+有限格点**(**Adhikari–Cao**,`Δ_G=0` 连续群塌缩);**无 4D 随机量子化**(4D 边缘非超可重整)、**无 4D 正则结构规范**;Abel no-go(Guth/FS:4D U(1) 弱耦无质量)杀"大 β⇒隙" |
| **4D LSI / Bianchi coarea** | 两尺度分解 LSI(Piazza–Shen)+ 群值 Brascamp–Lieb(CLL) | **→ H1** | 自由 U(1):`DΦ_B=d`、`σ_min=2sin(π/L)` **软模忠实保留、不造假隙**(酸测过,但 `∝1/σ_min~L` ⇒ 结构非一致);SU(2)"非阿贝尔增益"**假**(仅 `dim su(2)=3` 颜色重数);`σ_min(b)~b^{−1.37}`;**缺 = `b`-一致 `σ_min>0` = 边界一致 DLR-LSI = H1**;CLL 要群同态、plaquette holonomy 非;vixra(Eriksson)条件+窗口+自认 Clay 未解 |

---

## 1. 最强的正面事实：G1 已被机器验证（而它排除隙）
OpenAI 的 Lean 4 构造（Gross 主持委员会审阅）证明的是 **4D Wilson 格点规范（任意紧单群）标度极限 = Wightman QFT（G1，存在，OS 重构）**——9,470 模块、250,319 声明、**无项目公理 / 无 `sorryAx`**、Comparator 通过，且**白纸黑字**："does not prove … there is a mass gap"、"no evidence for or against a mass gap"、"A positive mass gap … outside its scope"。
**⇒ G1 基本被证成；整个问题塌缩到 G2（隙）在交叉点。**

## 2. 定理级 takeaway：原子是 usc，下界是 lsc
- 所有"下界工具"(sup/lsc：方差、磁化率、W₁、LSI、Γ₂、Schrödinger 桥)都在**原子**(`f₁²=ν({m₁²})`，usc)前止步——**下界证不出原子**。
- 原子只能由**精确恒等式**印：**A3**（Fredholm 行列式零点）+**B6**（全单调/Padé）；而它们都归约到"在连续统定位那个零点" = H1 / 边际 BPHZ。
- **A3 是吸引子**：每个"深"的原子可见方法都独立重推出它。

## 3. 最接近的 G2 候选 & 差在哪
- **Adhikari–Cao**（4D 弱耦、体积一致簇聚，**最接近物理耦合**）：但连续群 `Δ_G=0` ⇒ 衰变率塌缩；作者自称 SU(3) 是"central question" ⇒ 要证的就是 H1。
- **Bianchi-coarea / 4D LSI**（**最接近的技术结构**）：软模忠实（酸测过、不造假隙），但 `σ_min(Φ)~1/b` ⇒ 非一致；缺 `b`-一致 `σ_min>0` = 边界一致 DLR-LSI = H1。

## 4. 逐行审死的"声称证明"（宪章 §0，最高价值）
Lim 2307.00788（从不使用 YM 作用量；手设 Casimir H；酸测败）· Platek Zenodo 20114011（从不进强耦）· Faizal–Shabir 2606.19362（Thm 10.8 前提=交叉点、Lemma 10.5 参数身份错、Eq 10.72 调和非序列）· Farooqui Coq（定义性塞入）· jonwashburn Lean（build 失败）· Mondal/Sakib 轨道 BE（启发式/循环 / 1.31× 拓宽）· Eriksson vixra（条件+窗口+自认未解）。

## 5. 诚实底线
- **全网六路，无一种可能性。** 缺失对象精确定义不变：`α_s~O(1)` 处、非微扰、非循环、`(L,a)`-一致地对 `0⁺⁺` 残数 `f₁²` 的**符号控制下界** = H1 = 输运机器。
- **Clay(G2) = 0%。** 但有一个真·进展的坐标：**G1 已被机器证成、G2 被隔离为唯一剩余**；"下界证不出原子、原子需精确恒等式"这条现在被六路独立确认。
- 真正诚实的可交付物仍是那份**条件定理论文（P3 + 乘性输运，模 H1）**；无条件 Clay 取决于那台尚不存在的机器。

*Refs.* 六份子报告 `SURVEY-{RECENT-ARXIV,CONSTRUCTIVE,WILDCARD-DOMAINS,ADJACENT-TOOLS,AI-FORMAL}.md` + `BIANCHI-COAREA.md` + `w7/bianchi_coarea.py`；Douglas, *Nature Rev. Phys.* **8** 86–97 (2026)；Piazza–Shen arXiv:2610.35286；Cowling–Li–Liang arXiv:2602.10647；Adhikari–Cao *Ann. Probab.* **53** 140–174 (2025)；Shen–Zhu–Zhu *CMP* **400** 805 (2023)；Nissim arXiv:2510.22788；OpenAI Lean `yang-mills-group/OpenAI-YM-20261028`；Faizal–Shabir arXiv:2606.19362。
