# 工作记录：红外界 / Gaussian domination 的可用性核查（第 329 轮）

> 按 §1.06/§1.07：本文只做记录，不作汇报。

---

## A. 待核对象

**质量**下界**的唯一已知非微扰机制是 FSS（Fröhlich–Simon–Spencer）路线：**

$$\textbf{(GD) 高斯支配}\ \Longrightarrow\ \textbf{(IB) 红外界}\ \widehat C(k)\le\frac{c}{\beta\,\lambda(k)}\ \Longrightarrow\ \textbf{(MP) 指数聚类}$$

**本项目已有 `χ≤1/6`（振幅端），缺的是 (GD)/(IB)（速率端）。**

## B. 检索结果（逐条）

| 来源 | 内容 | 是否覆盖 `SU(2)` Wilson |
|---|---|---|
[arXiv 2608.08409](https://arxiv.org/abs/2608.08409) | 实为 *Notes on the Optical-Lattice Hard-Core Bose Gas of Aizenman–Lieb–Seiringer–Solovej–Yngvason*；其 Thm 5.8/5.10 是**玻色气体**的 GD/IB | **否**（此前一轮的归属错误，已更正） |
[arXiv 2611.31872](https://browse-export.arxiv.org/pdf/2611.31872#7) | `λ(k)Ĝ_L(k)≤c_0−2g_L`，自称强化 [24] 的 `λĜ_L≤1/K` | **否**（自旋型格点系统） |
[hexagonmath 2611.00321](https://hexagonmath.org/pdf/2611.00321v2#2) | 「格点上的反射正性、高斯支配及其推论可回溯到 [14,13]」 | 引文级，未见非 Abel 规范版本 |
[Zenodo 21473281](https://zenodo.org/records/21473281) | *Volume-Independent Lower Bound on the Mass Gap in 4D SU(N) LGT: From Strong-Coupling Cluster Expansion* | 强耦合簇展开（`β` 小端） |
[Chevyrev–Garban JSP 2025](https://hal.science/hal-04576352v1/preview/s10955-025-03420-1.pdf#page=2#1#5#1) | Ginibre/GKS 只在 Abel + Villain | 否 |

$$\boxed{\ \text{检索到的高斯支配/红外界定理，全部是【自旋型】或【玻色气体】或【Abel】；\ 非 Abel Wilson 版本未见于所查文献}\ }$$

## C. 自己推 (GD) 时的断点

**高斯支配的标准证明（`O(n)` 模型）用两件事：**
1. **反射正性 + 第二变分 `≥0`；**
2. **`e^{β⟨σ,σ'⟩}` 的两体结构。**

**非 Abel Wilson 的断点：**
- **plaquette 是**四链**乘积 ⇒ 第二变分不是两体形式；**
- **`q_P=½tr U_P` 是**特征标**（变号）⇒ 不是正权函数；**
- **文献明言的非 Abel 障碍：*"we cannot express Wilson loop observables as a product"*。**

**⇒ (GD) 在本项目的能力范围内无法自推。**

## D. 状态

| 件 | 状态 |
|---|---|
`χ^{clover}≤1/6`（振幅端，无条件） | `[T]` |
(GD)/(IB)（速率端） | `[U]`；所查文献均为自旋/Bose/Abel |
`(MP)` | `[U]` |

**⇒ 振幅端已闭合；速率端缺的是一条具名的、非 Abel 版本尚不存在的定理。**
