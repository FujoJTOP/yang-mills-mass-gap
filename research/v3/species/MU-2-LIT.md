# MU-2-LIT · 环路方程接触项的文献处理

**判定（第一段）**：**(b)**——文献里接触项**不是**被"吸收进质量项并给出系数"的；Makeenko 本人给出的处理是**显式 UV 正则化**（Halpern–Makeenko 连续正则化：把 `x=y` 的重合点用长度 `Λ⁻¹` 的 Brownian 桥"抹开"），正则化后的方程在 `Λ → ∞` 时回到原式。因此接触项 = **cutoff 尺度上的量**，其系数只能是 `Λ_cut` 的量纲幂，**不能**作为"物理质量/蝌蚪"的独立证据；`σ = 2g²N·μ²` 里的 `μ²` 在文献口径下就是**短程截断的倒数平方**（`μ² ≃ 1/ℓ² = Λ²`），而不是任意处方，也**不是**已证实的物理质量。

**对广播问题的直接回答（`1/a²` vs `1/Λ²`）**：
1. 系数是 `1/a` 型（UV 发散）——文献里它的名字是 cutoff `Λ`，不是有限物理量。
2. 有明确"正则化"说法：正则化方程 (2.36) "manifestly recovers Eq. (2.18) when Λ → ∞"。
3. 判：**`1/Λ²`**。理由：`a`（格距）与 `Λ⁻¹`（抹开长度）在文献里是**同一个短程尺度**；`δ³_⊥(0) ~ 1/ℓ²`，`ℓ = Λ⁻¹` ⇒ `μ² = Λ²`。所以 `ε = a` 只在 `a = 1/Λ_UV` 的意义下成立；把 `a` 当"物理长度"读才会得到 `σ ∝ 1/a²` 的错形状。

---

## 1. Makeenko–Migdal 原文：方程里有没有接触项？

- **Yu.M. Makeenko, A.A. Migdal, "Exact Equation for the Loop Average in Multicolor QCD", Phys. Lett. B88 (1979) 135**；勘误 Phys. Lett. B89 (1980) 437。DOI `10.1016/0370-2693(79)90131-X`（INSPIRE record 143896，已核）。
- **状态：正文未逐字核到（Elsevier 付费墙；检索式：INSPIRE `t "Exact equation for the loop average in multicolor QCD"`、arXiv 无该文）。** 因此"原文是否写出 `δ³_⊥(0)` 接触项"这一条，本文件**不作断言**。
- `Nucl. Phys. B188 (1981) 269`（"QCD as a theory of strings" 一系）：**未核到书目记录**（检索式：INSPIRE `a Makeenko and j Nucl.Phys.B and de 1981`、`a Migdal and t "strings" and de 1981`）。**引用前必须再核。**

## 2. Makeenko 的后续：接触项怎么正规化（唯一逐字核到的权威文本）

**Yuri Makeenko, "Notes on the Loop Equation in Loop Space", arXiv:2508.09705 [hep-th]（2025-08-13 提交；正文自注 "YM-7-94, December, 1994"）**，§2.6 "Nonperturbative regularization"。逐字：

> "A convenient way of ultraviolet regularization of the loop equation (2.18) was proposed by Halpern and the author [17] using the relation to the stochastic quantization discussed in the previous subsection."

> （Eq. 2.36）"Averaging over the gauge field, A, one arrives at the regularized loop equation
> ΔW(C) = λ ∮_C dx_μ ∮_C dy_μ ∫_{r(0)=x}^{r(Λ^{-2})=y} Dr e^{−(1/2)∫_0^{Λ^{-2}} dτ ṙ²(τ)} W(C_xy r_yx) W(C_yx r_xy)
> which **manifestly recovers Eq. (2.18) when Λ → ∞**."

> "The advantage of this regularization of the loop equation is that the contours C_xy r_yx and C_yx r_xy on the r.h.s. of Eq. (2.36) both are closed and **do not have marked points** if C does not have. Therefore, Eq. (2.36) is written entirely on loop space."

**读法**：`x = y` 的重合不是被"吸收"，而是被**替换**为一个长度 `Λ⁻¹` 的开链 `r`（proper time `Λ⁻²`）。接触项的全部内容 = `Λ` 依赖；`Λ → ∞` 才回到含 `δ(x−y)` 的原方程。**没有任何"吸收进质量项/周长项"的句子**，也没有给出"吸收系数 = 某物理质量"的表述。

参考文献（该文 list 已逐字核）：
- [17] M.B. Halpern, Yu.M. Makeenko, "Continuum regularized loop-space equation", **Phys. Lett. B218 (1989) 230**。
- [21] Z. Bern, M.B. Halpern, L. Sadun, C. Taubes, "Continuum regularization of quantum field theory. 2. Gauge theory", **Nucl. Phys. B284 (1987) 35**。
- [22] Yu.M. Makeenko, "Polygon discretization of the loop-space equation", **Phys. Lett. B212 (1988) 221** ← 对应"多边形/格点"版本，即 `a` 型截断。

## 3. "带质量项的环路方程"：有没有那一句？

**未查到**（检索式：arXiv API `abs:"loop equation" AND abs:"mass term"` → totalResults = 0；`abs:"loop equation" AND abs:"contact term"` 无对应权威原文；Google 检索 "loop equation" "mass term" Migdal perimeter 只回落到非权威二手页面）。因此**本文件不引用任何"质量项由接触项生成"的逐字句子**——这样的句子我没能核到。凡在本项目文档中写"文献说接触项吸收进质量项"，都必须先补上这一条引用，否则属于未核断言。

## 4. Polyakov / 综述 / NPB 线：周长律

**未逐字核到**（Migdal, Phys. Rept. 102 (1983) 199 与 Makeenko 的 "Large N gauge theories" (hep-th/0001047) 均为付费/无 HTML 版本；ar5iv 对该 id 返回 404）。唯一能确定的口径仍是 §2 的 `Λ` 截断：**周长律的系数在文献口径下属于"要由截断/重整化条件固定"的量**，而不是由方程自身算出的数。

---

## 对 `rv-2` 记账的含义（诚实裁决）

1. 文献把接触项当**正则化对象**（`Λ`），不当**物理质量**。所以 `μ²` **有物理意义**这一条**没有**文献支持；`[C]` 不能靠"文献说吸收进质量项"活。
2. 但 `μ² = δ³_⊥(0)` 也**不是**任意处方：在 `ℓ = Λ⁻¹` 的短程截断下 `μ² = Λ²`（唯一自洽读法），故 `σ = 2g²N Λ²`——**"σ ∝ 截断平方"是文献级别的标准形状**，`σ ∝ 1/a²`（把 `a` 当物理长度）才是读错。
3. 缺口精确化：**"哪个重整化条件把 `Λ²` 减掉、剩下什么"——文献未给**。这一步是我们自己的账，`(b)` 的缺口即此。
