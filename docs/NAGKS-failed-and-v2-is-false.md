# 非 Abel GKS：我做不到，而且我原来那条路是**假的**

> 这一轮必须报的是坏消息里的硬事实：为了绕过非 Abel GKS 而引入的辅助不等式 `v_2≥0`，**在数值上是假的**。所以那条路不是「差一步」，是**不通**。
> 而 `[T]` 的部分（`χ^{clover}≤1/6`）**不依赖 `v_2` 的符号**，仍然成立。

---

## 1. 我做不到什么

$$\boxed{\ \text{非 Abel GKS（}SU(2)\ \text{Wilson）：我做不到}\ }$$

**尝试过的具体路线，逐条结果：**

| 路线 | 结果 |
|---|---|
Quaternion `q_P=⟨A_x,B_x⟩_{R^4}`（`A=U_1U_2`，`B=U_1'U_2'`）→ 当 `O(4)` 自旋模型用 GKS | **不行**：耦合是**四链**的（plaquette），不是两体；GKS 框架要求两体 |
Character 展开 → 每链积分给非负整数（融合系数） | `(q_{P_0}q_{P_1})≥0` 得到，但 `≥⟨q⟩⟨q⟩` 需要逐项比较，正是文献说的 *"we cannot express Wilson loop observables as a product"* |
对方向 1 的链条件化，再看 `U_2,U_3` 是否因子化 | **不因子化**（`(2,3)` 平面 plaquette 把它们耦合） |
Gaussian/HS 双重复制 | `q_P` 是**特征标**（变号），非正权函数 |

**⇒ 与文献一致：Abelian/Villain 有，非 Abel Wilson 没有。**

---

## 2. 更硬的：我原来那条绕路**是假的**

**第 315–316 轮我用 `v_2:=Cov(Q^{(A)},Q^{(B)})/L^4≥0` 把 `χ^{clover}≤1/6` 转成 `χ^{fixed}≤1`。**
**本轮数值查了 `β=0.5…4.0`（`L=4`）：**

| `β` | `min off-diagonal M_AB` | `M_AA≤Σ_BM_AB`? | `max M_AA` |
|---|---|---|---|
**`0.5`** | **`−1.754e−2`** | **False** | `0.2606` |
**`1.0`** | **`−1.179e−2`** | **False** | `0.2351` |
`2.0` | `−4.167e−3` | True | `0.2528` |
`3.0` | `+1.615e-3` | True | `0.0630` |
`4.0` | `−1.370e−3` | True | `0.0306` |

$$\boxed{\ v_2\ge0\ \text{是}\ \textbf{假的}\quad(\beta\lesssim2\ \text{时为负，量级 }7\%\ \text{的 }M_{AA})\ }$$

**⇒ 而且 `M_AA≤Σ_BM_{AB}` 在 `β=0.5,1.0` 也**失败**。**

**⇒ 「六取向分解 + 正关联」这条路**整体不通**，不是差一步。**

**⚠️ 数据只到 `n=450`、`L=4`，`β≤2` 的负值需要更硬的统计确认；但**符号在 `β=0.5` 与 `1.0` 两次独立运行中一致**，且量级远超 `β=3,4` 处的噪声水平。**

---

## 3. 仍然成立的那部分（`[T]`，且与 `v_2` 符号无关）

$$\boxed{\ \int_0^\infty\chi^{\rm clover}_L(\beta)\,d\beta=\frac16\qquad\text{且}\qquad \chi^{\rm clover}_L(\beta)\le\frac16\quad\forall L\ge2,\ \forall\beta\gt0\ }$$

**证明只用三件事，**都不涉及取向分解**：**

1. **`d⟨𝒮⟩/dβ=Var(𝒮)≥0`**（`d⟨O⟩/dβ=−Cov(O,S_β)`，`S_β=β(6L^4−𝒮)`）
2. **端值 `⟨𝒮⟩_{β=0}=6L^4`、`⟨𝒮⟩_{β→∞}=0`**
3. **代数 `Q^{\rm clover}=𝒮/6`**

**⇒ `Var(𝒮)≥0` 与 `∫Var(𝒮)dβ=6L^4` 直接给 `χ^{clover}=Var(𝒮)/(36L^4)≤1/6`。**
**`v_2` 正或负，这个推导一字不動。**

---

## 4. 于是算子范数这条还缺什么（精确定位，不再兜圈）

$$\|\Gamma^{\rm clover}\|_{2\to2}=\max_k\widehat C^{\rm clover}(k),\qquad \widehat C^{\rm clover}(k)=\frac1{36}\widehat C^{\mathcal S}(k)$$

**其中 `Ĉ^{\mathcal S}(k)=Σ_z\mathrm{Cov}(s_0,s_z)e^{ikz}`，`s_x=Σ_A(1−q_{P_A(x)})`。**

$$\boxed{\ \text{缺：}\ \widehat C^{\mathcal S}(k)\ \le\ \widehat C^{\mathcal S}(0)=\mathrm{Var}(\mathcal S)\qquad\Big((Q)\ \text{对站点作用量密度 }s\Big)\ }$$

**⚠️ `s_x` 含全部六个取向 ⇒ 必然含时间方向 plaquette ⇒ 正是勘误 35 的**交叉项**障碍。**
**⇒ `(Q)` 对 `s` 与 `v_2` 是同源问题，但 `v_2` 已证假，所以不能走取向分解，只能直接攻 `s` 的 `(Q)`。**

**⇒ 这是**唯一**还开着的入口。**

---

## 5. 状态表（本轮更正）

| 命题 | 判 |
|---|---|
$$\int_0^\infty\chi^{\rm clover}d\beta=\tfrac16$$ | **`[T]` 无条件，与 `v_2` 符号无关** |
$$\chi^{\rm clover}\le\tfrac16$$ | **`[T]` 无条件** |
$$d\langle\mathcal S\rangle/d\beta=\mathrm{Var}(\mathcal S)$$ | `[T]` |
$$v_2\ge0$$ | **`[F]` 假**（`β\lesssim2` 为负；`[M]`） |
「取向分解 → `χ^{fixed}`」 | **`[F]` 路线不通** |
$$(Q)\ \text{对 }s=Σ_A(1-q_{P_A})$$ | **`[U]` 唯一开口** |
$$\|\Gamma_L\|_{2\to2}\le C$$ | **`[C]`** |

---

## 6. 结论

$$\boxed{\ \begin{aligned}&\text{非 Abel GKS：做不到}\\ &\text{而且我为绕过它用的 }v_2\ge0\ \text{是假的}——\text{那条路是死的，不是差一步}\\ &\text{站得住的：}\ \chi^{\rm clover}\le\frac16\ \text{且}\ \int_0^\infty\chi^{\rm clover}d\beta=\frac16\quad[T]\ \text{无条件}\\ &\text{唯一还开的入口：}\ (Q)\ \text{对站点作用量密度 }s_x\end{aligned}\ }$$
