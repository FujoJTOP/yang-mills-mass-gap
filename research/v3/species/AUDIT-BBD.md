# Yang--Mills mass gap：Bauerschmidt--Bodineau--Dagallier 框架核查

## 核查结论

**理论入口真实存在，但目前还不能据此宣称已经得到 Yang--Mills
物理质量隙。**

新版第 15.5 节所使用的核心思路------允许重整化势具有负 Hessian
方向，再通过紧李群的 Ricci
曲率补偿，并将负曲率代价压缩为一个可积函数------与
Bauerschmidt--Bodineau--Dagallier 的黎曼几何推广框架有实质对应。

关键文献不是单独的 Theorem 3.6，而是综述中的 **Theorem 3.20**。Theorem
3.6 是欧氏高斯协方差框架；Theorem 3.20 才包含非平坦黎曼流形上的 Ricci
曲率和时间变化度量项。

------------------------------------------------------------------------

## 1. 原始文献

相关作者姓名应写为：

-   Roland Bauerschmidt
-   Thierry Bodineau
-   Benoit Dagallier

主要综述：

-   Roland Bauerschmidt, Thierry Bodineau, Benoit Dagallier, *Stochastic
    dynamics and the Polchinski equation: an introduction*,
    **Probability Surveys**, Vol. 21, pp. 200--290 (2024).
-   DOI: [10.1214/24-PS27](https://doi.org/10.1214/24-PS27)

相关应用论文：

-   [arXiv:2202.02301](https://arxiv.org/abs/2202.02301)
-   [arXiv:2202.02295](https://arxiv.org/abs/2202.02295)
-   [arXiv:2307.07619](https://arxiv.org/abs/2307.07619)

注意：此前草稿中的 "Bachmann--Dagallier" 是作者姓名误写，应改为
**Bauerschmidt--Bodineau--Dagallier**，并根据具体命题使用准确的作者组合。

------------------------------------------------------------------------

## 2. Theorem 3.6：欧氏高斯版本

Theorem 3.6
的基本环境是有限维欧氏空间中的线性空间，以及随时间变化的高斯协方差分解
(C_t)。

其关键矩阵不等式为

\[
`\dot `{=tex}C_t,`\operatorname{Hess}`{=tex}V_t(`\varphi`{=tex}),`\dot `{=tex}C_t
-`\frac12`{=tex}`\ddot `{=tex}C_t `\geq`{=tex}
`\dot`{=tex}`\lambda`{=tex}\_t`\dot `{=tex}C_t. `\tag{3.32}`{=tex} \]

这里：

-   (C_t)：高斯协方差；
-   (V_t)：由高斯卷积/重整化构造的有效势；
-   (`\dot`{=tex}`\lambda`{=tex}\_t)：允许为负的时间函数；
-   (`\operatorname{Hess}`{=tex}V_t)：有效势的 Hessian。

定义

\[
`\lambda`{=tex}\_t=`\int`{=tex}\_0\^t`\dot`{=tex}`\lambda`{=tex}\_s,ds,
`\qquad`{=tex}
`\frac1`{=tex}`\gamma`{=tex}=`\int`{=tex}\_0^`\infty `{=tex}e^{-2`\lambda`{=tex}\_t},dt.
\]

在原文列出的假设及连续性条件下，可以推出相应的 Log-Sobolev 不等式：

\[ `\operatorname{Ent}`{=tex}\_{`\nu`{=tex}*0}(F) `\leq`{=tex}
`\frac2`{=tex}`\gamma`{=tex} `\mathbb `{=tex}E*{`\nu`{=tex}\_0} `\left[
(\nabla\sqrt F)^2_{\dot C_0}
\right]`{=tex}. \]

重要之处在于：原始势不必凸，沿重整化流的 Hessian
也可以出现负方向；关键是相应的时间累计控制条件成立。

这与新版第 15.5 节"允许局部负曲率、控制累计负曲率代价"的思路一致，但
Theorem 3.6 本身不是紧李群上的几何定理。

------------------------------------------------------------------------

## 3. Theorem 3.20：黎曼几何版本

对紧李群设置，更直接相关的是综述中的 Theorem 3.20。

原文的几何条件为

\[ `\boxed{
\operatorname{Hess}_{g_t}V_t
+\frac12\operatorname{Ric}_{g_t}
+\frac12\dot g_t
\geq
\dot\lambda_t g_t.
}`{=tex} `\tag{3.100}`{=tex} \]

其中：

-   (g_t)：时间依赖的黎曼度量；
-   (`\operatorname{Ric}`{=tex}\_{g_t})：该度量的 Ricci 张量；
-   (`\dot `{=tex}g_t)：度量对时间的导数；
-   (V_t)：重整化势；
-   (`\dot`{=tex}`\lambda`{=tex}\_t)：允许为负的控制函数。

在原文要求的连续性等条件成立时，该框架给出相应的 Log-Sobolev 不等式。

**这才是把 Polchinski
多尺度方法迁移到非平坦紧李群时，应当优先核对的定理。**

------------------------------------------------------------------------

## 4. 代入固定双不变度量的紧李群

考虑固定双不变度量：

\[ g_t=g,`\qquad `{=tex}`\dot `{=tex}g_t=0. \]

若参考体积取为固定的 Haar 体积，并且热半群与该度量相匹配，则 Theorem
3.20 的几何条件简化为

\[ `\operatorname{Hess}`{=tex}\_gV_t
+`\frac12`{=tex}`\operatorname{Ric}`{=tex}\_g `\geq`{=tex}
`\dot`{=tex}`\lambda`{=tex}\_tg. \]

若所采用的度量归一化满足

\[ `\operatorname{Ric}`{=tex}\_g=`\frac `{=tex}N4g, \]

那么

\[ `\boxed{
\dot\lambda_t
\leq
\frac N8+
\inf_\phi
\lambda_{\min}
\left(\operatorname{Hess}_gV_t(\phi)\right).
}`{=tex} \]

这解释了新版第 15.5 节中 (N/8) 的 Ricci 补偿项。

定义

\[ q(t)= `\sup`{=tex}*`\phi`{=tex} `\left`{=tex}\|
`\left`{=tex}(`\operatorname{Hess}`{=tex}*gV_t(`\phi`{=tex})`\right`{=tex})*-
`\right`{=tex}\|*{`\mathrm{op}`{=tex}}, \]

其中 ((A)\_-) 表示对称算子 (A) 的负部分。若

\[ `\operatorname{Hess}`{=tex}\_gV_t`\succeq`{=tex}-q(t)g, \]

则可以选择

\[ `\dot`{=tex}`\lambda`{=tex}\_t=`\frac `{=tex}N8-q(t). \]

因此

\[ `\lambda`{=tex}\_t=`\frac `{=tex}N8t-Q(t), `\qquad`{=tex}
Q(t)=`\int`{=tex}\_0\^tq(s),ds. \]

如果

\[ Q(`\infty`{=tex})\<`\infty`{=tex}, \]

则

\[
```{=tex}
\begin{aligned}
\frac1\gamma
&=\int_0^\infty e^{-2\lambda_t}\,dt\\
&=\int_0^\infty e^{-Nt/4+2Q(t)}\,dt\\
&\leq e^{2Q(\infty)}
\int_0^\infty e^{-Nt/4}\,dt\\
&=\frac4N e^{2Q(\infty)}.
\end{aligned}
```
\]

于是

\[ `\boxed{
\gamma\geq\frac N4e^{-2Q(\infty)}.
}`{=tex} \]

这与新版第 15.5 节给出的核心估计链条一致。

### 必须保留的假设

上述推导并不自动验证所有应用条件。至少还需要核实：

1.  所选紧李群、度量和热半群确实符合 Theorem 3.20 的几何设定；
2.  参考体积、热核及重整化势的定义与原文约定一致；
3.  Theorem 3.20 所要求的连续性和可积性条件成立；
4.  Ricci 张量的归一化 (`\operatorname{Ric}`{=tex}\_g=(N/4)g)
    与实际使用的 Lie algebra 内积、生成元归一化一致；
5.  (q(t))
    是对整个配置空间取上确界的严格控制，而不是有限样本或单一格点的数值估计。

------------------------------------------------------------------------

## 5. 群卷积下的 Hessian---协方差恒等式

考虑

\[ P_te\^{-S_W}(x) = `\int `{=tex}K_t(z)e\^{-S_W(xz)},dz, \]

并定义条件概率测度

\[ `\mu`{=tex}\_{t,x}(dz) `\propto`{=tex} K_t(z)e\^{-S_W(xz)},dz. \]

令

\[ V_t(x)=-`\log`{=tex}(P_te\^{-S_W})(x). \]

在将不同群点的切向量通过左平移识别到共同切空间，并使用与双不变度量相容的联络后，对数积分求导给出

\[ `\nabla `{=tex}V_t =
`\mathbb `{=tex}E\_{`\mu`{=tex}\_{t,x}}\[`\nabla `{=tex}S_W\], \]

以及

\[ `\boxed{
\operatorname{Hess}V_t
=
\mathbb E_{\mu_{t,x}}
[\operatorname{Hess}S_W]
-
\operatorname{Cov}_{\mu_{t,x}}
(\nabla S_W).
}`{=tex} \]

协方差项前面的负号来自对数配分积分的二阶导数结构，并非数值拟合出来的修正项。

不过，正式证明仍应明确：

-   梯度如何在不同切空间之间平行移动或通过左平移识别；
-   Hessian 使用哪一种联络；
-   群乘法 (xz) 的微分如何进入二阶导数；
-   热核的变量约定及左、右不变性；
-   是否存在额外的联络或曲率项。

数值误差达到 (10\^{-8}) 可以验证程序实现，但不能替代这些解析步骤。

------------------------------------------------------------------------

## 6. 对新版第 15.5 节的逐项修订建议

  --------------------------------------------------------------------------------------------------------------------------------
  当前写法或主张          核查结果                                   建议
  ----------------------- ------------------------------------------ -------------------------------------------------------------
  Bachmann--Dagallier     作者姓名误写                               改为
                                                                     Bauerschmidt--Bodineau--Dagallier，并按具体论文核对作者列表

  以 Theorem 3.6          Theorem 3.6 是欧氏高斯框架                 将 Theorem 3.20 作为黎曼几何主引用
  作为群几何版本                                                     

  群上的条件只写 Hessian  几何版本包含 Ricci 项和 (`\dot `{=tex}g_t) 明确列出 Theorem 3.20 的完整条件
                          项                                         

  (C_t=tg\^{-1})          这是协方差式记法，不等于完整的群热流定义   明确区分欧氏协方差参数与群上的热半群

  引入 (N/8)              在 (`\operatorname{Ric}`{=tex}\_g=(N/4)g)  给出度量和 Lie algebra 归一化
                          时与半 Ricci 项一致                        

  "只剩一个标量"          仅在其余分析环节已严格闭合时成立           明确 (Q(`\infty`{=tex}))
                                                                     的统一性是当前关键缺口，同时列出尚待核验的适用条件
  --------------------------------------------------------------------------------------------------------------------------------

建议把第 15.5 节分成三个独立层次：

1.  **原文定理层**：完整陈述 Theorem 3.20；
2.  **群几何应用层**：证明固定双不变度量下的简化条件及 Hessian 恒等式；
3.  **Yang--Mills 统一性层**：证明 (Q(`\infty`{=tex}))
    对格点尺度和体积的一致控制。

不要把这三个层次合并成一句"定理直接适用"。

------------------------------------------------------------------------

## 7. 最关键的逻辑边界：配置空间谱隙不自动等于物理质量隙

Theorem 3.20 所给出的 Log-Sobolev 不等式，控制的是相应概率测度及其
Markov 半群的函数空间不等式。它可以推出相应 Dirichlet 生成元的 Poincaré
型谱隙。

但这并不自动等于 Yang--Mills 物理 Hilbert
空间上的转移矩阵谱隙，更不自动给出 (0\^{++}) 通道的胶球质量隙。

两者之间需要独立建立严格的谱论桥梁，包括：

-   配置空间 Dirichlet 生成元与物理转移矩阵之间的关系；
-   规范不变物理子空间的处理；
-   连续极限和无限体积极限中的统一性；
-   (0\^{++}) 谱测度上的正谱原子；
-   原子的正残数及其隔离。

另一个重要提醒：固定紧李群的 Ricci 曲率可以在 (V_t`\equiv0`{=tex})
时仍然提供几何 Log-Sobolev
控制。因此，不能把这个几何常数直接认作一个必须在无质量自由场极限消失的物理质量隙。需要区分配置空间几何谱隙与物理质量隙。

------------------------------------------------------------------------

## 8. 最终判断

### 已经有根据的结论

1.  Bauerschmidt--Bodineau--Dagallier
    的多尺度框架确实允许沿重整化流出现负 Hessian 方向。
2.  原综述 Theorem 3.20 提供了包含 Ricci
    曲率和时间变化度量项的黎曼几何版本。
3.  在固定双不变度量下，若 Ricci 归一化为 ((N/4)g)，则出现 (N/8)
    补偿项。
4.  若负 Hessian 的统一积分代价 (Q(`\infty`{=tex}))
    有限，就能得到上述形式的配置空间 Log-Sobolev 常数下界。

### 仍未由这条定理单独解决的事项

1.  (Q(`\infty`{=tex})) 对 ((L,a)) 的一致控制；
2.  相关连续性与热流适用条件的完整验证；
3.  配置空间不等式到物理 (0\^{++}) 谱隙的严格传递；
4.  正残数和谱隔离。

**因此，准确的定位是：发现了一个有原始定理支撑的紧李群几何分析入口，并把一个重要的统一性问题集中到
(Q(`\infty`{=tex}))；但这条入口本身尚不是 Yang--Mills 物理质量隙证明。**

------------------------------------------------------------------------

## 参考资料

1.  Bauerschmidt, Bodineau, Dagallier, *Stochastic dynamics and the
    Polchinski equation: an introduction*, Probability Surveys 21
    (2024), 200--290. [DOI](https://doi.org/10.1214/24-PS27)
2.  [arXiv:2202.02301](https://arxiv.org/abs/2202.02301)
3.  [arXiv:2202.02295](https://arxiv.org/abs/2202.02295)
4.  [arXiv:2307.07619](https://arxiv.org/abs/2307.07619)
5.  Yang--Mills 项目主论文：第 15.5
    节及相关附录。本文对项目论文的描述以本轮已核查的上传版本为准。
