# A Hodge–Helffer–Sjöstrand reduction for the plaquette-scalar covariance in $d=4$ lattice $\mathrm{SU}(2)$ gauge theory

**Abstract.** Let $\mu_{\beta,L}$ be the Wilson measure on $\mathrm{SU}(2)^{E(\Lambda_L)}$, $\Lambda_L$ the $4$-torus of side $L$, and $q_x=\tfrac12\operatorname{tr}U_{P(x)}$ the scalar plaquette observable. Writing $\Gamma_L(x,y)=\operatorname{Cov}_\mu(q_x,q_y)$, we prove the exact identity
$$\langle f,\Gamma_L f\rangle=\big\langle dF_f,\ \mathscr H_{\beta,L}^{-1}\,dF_f\big\rangle,$$
where $\mathscr H_{\beta,L}=dd_\mu^*+d_\mu^*d$ is the Witten–Hodge Laplacian on $1$-forms and $F_f=\sum_x f(x)(q_x-\langle q_x\rangle)$. We prove $\|d\|^2\le 24$ uniformly in $L$ and $\beta$, and we reduce the operator-norm bound $\|\Gamma_L\|_{2\to2}\le C\Lambda_L$ ($\Lambda_L\asymp 1+\log L$ the critical-kernel row sum) to a single, momentum-resolved **infrared spectral-leakage hypothesis (M)**.

We then analyse (M) honestly and prove two negative results that delimit what this machine can achieve: (i) (M) is **strictly stronger** than the target — it forbids any atom at spectral parameter $0$, whereas the target permits one of mass $o(\varepsilon/\log L)$ at $\lambda=\varepsilon$; (ii) the natural weakening of (M) which is equivalent to the target coincides with the target at zero momentum. Hence the machine is a correct **reformulation**, not a reduction in difficulty.

The machine's value is diagnostic. It shows that a local Hessian eigenvalue exceeding $1$ — for the present model, the transverse eigenvalue $2.6975$ at $\beta=2.5$ — is **not** an obstruction to the target, because the target constrains only $\langle D_k,\mathscr H^{-1}D_k\rangle$ on the **scalar** channel. Numerical data at $L=6,8,10,12$, $\beta=2.5$ is reported which pinpoints the remaining content as a single quantitative statement: *negative (Savvidy-type) modes of $\operatorname{Hess}S$ do not couple to the plaquette-scalar channel at low $\mathscr H$-eigenvalue.*

---

## 1. Setting

Let $\Lambda_L=\{0,\dots,L-1\}^4$ with periodic boundary conditions, $L$ even. Let $E_L$ be the set of **directed** edges, with $U_{e^{-1}}=U_e^*$. Let $dU_e$ denote Haar probability measure on $\mathrm{SU}(2)$.

For a plaquette $P=(e_1,e_2,e_3,e_4)$ (boundary in cyclic order), put $U_P=U_{e_1}U_{e_2}U_{e_3}U_{e_4}$ and

$$S_{\beta,L}(U)=\beta\sum_{P}\Big(1-\tfrac12\operatorname{tr}U_P\Big),\qquad
d\mu_{\beta,L}(U)=\frac{1}{Z_{\beta,L}}e^{-S_{\beta,L}}\prod_{e\in E_L}dU_e .$$

Here $\beta>0$; the range of interest is $\beta\in[1/48,\,2.5]$.

**Observable.** For each plaquette $P$ in the positive $12$-plane based at a lattice site, let $x$ denote its centre, so that $x\mapsto P(x)$ is a bijection onto the set of such plaquettes. Put

$$\boxed{\ q_x:=\tfrac12\operatorname{tr}U_{P(x)}\in[-1,1]\ }$$

This is a real, scalar, **gauge-invariant** function of $U$.

**Target operator.**

$$\Gamma_L(x,y):=\operatorname{Cov}_{\mu_{\beta,L}}(q_x,q_y),\qquad \|\Gamma_L\|_{2\to2}:=\sup_{\|f\|_{\ell^2}\le1}\|\Gamma_L f\|_{\ell^2}.$$

**Critical kernel.** Fix $g_c(z)=(1+\|z\|_2)^{-4}$ and

$$\boxed{\ \Lambda_L:=1+\sup_{x\in\Lambda_L}\sum_{y\in\Lambda_L}g_c(x-y)\ }$$

with $\|x-y\|_2$ the minimum-image distance. One has $\Lambda_L=\alpha\ln L+O(1)$ with $\alpha=2\pi^2$ for the Euclidean kernel ($\alpha=8/3$ for $\ell^1$; the ratio is $3\pi^2/4=7.402203$).

**Note on normalization.** With $g_c(0)=(1+0)^{-4}=1$ included in the sum, one computes
$$\Lambda_L\big|_{L=6,8,10,12}=7.7269,\ 10.2988,\ 12.6480,\ 14.7840 .$$

---

## 2. The problem

$$\boxed{\ \textbf{(T)}\qquad \|\Gamma_L\|_{2\to2}\ \le\ C\,\Lambda_L\qquad\text{uniformly in }L\text{ and }\beta\in[1/48,2.5]\ }\tag{T}$$

with $C$ independent of $L,\beta$.

Sections 3–5 reduce (T) to a spectral hypothesis. Sections 6–7 delimit that reduction. Section 8 gives numerical evidence. Section 9 isolates the remaining content.

---

## 3. Exact Hodge identity

Let $\mathcal M_L=\mathrm{SU}(2)^{E_L}$ with the product of bi-invariant metrics. Let $\mathcal L_\mu=\Delta-\nabla S_{\beta,L}\cdot\nabla$ be the (negative semidefinite) weighted Laplacian, $d$ the exterior derivative, $d_\mu^*$ its $L^2(\mu)$-adjoint, and

$$\boxed{\ \mathscr H_{\beta,L}:=dd_\mu^*+d_\mu^*d\ }\tag{3.1}$$

the Witten–Hodge (Helffer–Sjöstrand) operator on $1$-forms. On $1$-forms the Weitzenböck formula reads

$$\mathscr H_{\beta,L}=\nabla_\mu^*\nabla+\operatorname{Ric}+\operatorname{Hess}S_{\beta,L}.\tag{3.2}$$

Since $\mathcal M_L$ is a compact Lie group with bi-invariant metric, $\operatorname{Ric}\succeq0$; the term $\operatorname{Hess}S_{\beta,L}$ is **indefinite**, and we do not bound it from below.

For $f\in\ell^2(\Lambda_L)$ write $F_f:=\sum_x f(x)\,(q_x-\langle q_x\rangle_\mu)$, so that

$$\langle f,\Gamma_L f\rangle=\int|F_f|^2\,d\mu_{\beta,L}.\tag{3.3}$$

**Theorem 3.1 (exact covariance identity).** Let $F$ be a centred real function on $\mathcal M_L$ (i.e. $\int F\,d\mu=0$) and let $u$ be the unique centred solution of $-\mathcal L_\mu u=F$. Then $du=\mathscr H_{\beta,L}^{-1}dF$ on the space of exact $1$-forms, and

$$\boxed{\ \int F^2\,d\mu=\big\langle dF,\ \mathscr H_{\beta,L}^{-1}\,dF\big\rangle\ }\tag{3.4}$$

Consequently

$$\Gamma_L=D^*\mathscr H_{\beta,L}^{-1}D,\qquad Df:=dF_f.\tag{3.5}$$

*Proof.* Ellipticity on a compact manifold makes $-\mathcal L_\mu$ invertible on centred functions, so $u$ exists and is unique. By definition of $d_\mu^*$, $d_\mu^*du=-\mathcal L_\mu u=F$. Since $d^2=0$,
$$\mathscr H du=dd_\mu^*du+d_\mu^*d^2u=dF .$$
Hence $du=\mathscr H^{-1}dF$. Integrating by parts,
$$\int F^2d\mu=\langle F,d_\mu^*du\rangle=\langle dF,du\rangle=\langle dF,\mathscr H^{-1}dF\rangle .$$
For complex $f$ apply the identity to real and imaginary parts. $\square$

**Remark 3.2.** (3.4) is an identity, not an estimate. No spectral gap of $\mathscr H$ is assumed; when $\mathscr H$ has zero modes, the inverse is taken on the subspace orthogonal to them, which is legitimate because $dF_f$ has no component there whenever $\int F_f\,d\mu=0$. $\square$

---

## 4. The source operator

**Theorem 4.1 (uniform source norm).** For every $L,\beta$,

$$\boxed{\ \|D\|^2=\sup_{\|f\|_2\le1}\|dF_f\|^2_{L^2(\mu)}\ \le\ 24\ }\tag{4.1}$$

*Proof.* Write $U_P=a_0I+i\sum_{j=1}^3a_j\sigma_j$ with $a_0^2+|a|^2=1$. For a unit tangent perturbation $U_e(t)=e^{tX}U_e$, $X=i\sum_j\xi_j\sigma_j$, $|\xi|=1$, a direct differentiation gives $|d_eq_P|\le1$.

Every edge belongs to at most $6$ plaquettes, and every plaquette has $4$ edges. By Cauchy–Schwarz over the at most $6$ summands,
$$\|dF_f\|^2=\sum_{e}\Big|\sum_{x:\,e\in\partial P(x)}f(x)\,d_eq_x\Big|^2\le 6\sum_e\sum_{x:\,e\in\partial P(x)}|f(x)|^2 .$$
Now $\sum_e\sum_{x:e\in\partial P(x)}|f(x)|^2=\sum_x|f(x)|^2\cdot 4=4\|f\|_2^2$, whence $\|dF_f\|^2\le 24\|f\|_2^2$. $\square$

**Proposition 4.2 (Ward).** Let $\xi$ be the infinitesimal generator of a vertex gauge transformation, $\delta_\xi U_{x,\mu}=\xi_xU_{x,\mu}-U_{x,\mu}\xi_{x+e_\mu}$. Then $dq_P(\delta_\xi U)=0$, hence

$$dF_f(\delta_\xi U)=0.\tag{4.2}$$

*Proof.* $\operatorname{tr}U_P$ is gauge invariant. $\square$

**Remark 4.3.** (4.2) places the source in the orthogonal complement of the gauge orbits. It does **not** by itself show that the physical transverse projection is small; that requires the full plaquette-block integration of §9. $\square$

---

## 5. Spectral reduction

By translation invariance $\Gamma_L(x,y)=C_L(y-x)$; plane waves diagonalize $\Gamma_L$ and

$$\|\Gamma_L\|_{2\to2}=\max_{k}\widehat C_L(k).\tag{5.1}$$

Let $f_k(x)=L^{-2}e^{ik\cdot x}$ ($\|f_k\|_2=1$), $k\in\tfrac{2\pi}{L}\mathbb Z_L^4$, $D_k:=dF_{f_k}$, and let $\nu_{k,L,\beta}$ be the spectral measure of $\mathscr H_{\beta,L}$ against $D_k$:

$$\nu_{k,L,\beta}(A):=\big\langle D_k,\ \mathsf E_{\mathscr H}(A)\,D_k\big\rangle_{L^2(\mu)},\qquad \nu_{k,L,\beta}\big((0,\infty)\big)=\|D_k\|^2\le24 .$$

Put $\omega_L(k):=L^{-2}+4\sum_{\mu=1}^4\sin^2(k_\mu/2)$.

**Theorem 5.1 (spectral representation).**

$$\boxed{\ \widehat C_L(k)=\int_{(0,\infty)}\frac{1}{\lambda}\,d\nu_{k,L,\beta}(\lambda)\ }\tag{5.2}$$

*Proof.* Combine (3.4) with the spectral theorem. $\square$

**Lemma 5.2 (layer-cake).** For a positive measure $\nu$ on $(0,\infty)$ with $\int\lambda^{-1}d\nu<\infty$,

$$\int_{(0,\infty)}\frac{d\nu(\lambda)}{\lambda}=\int_0^\infty\frac{\nu([0,s])}{s^2}\,ds.\tag{5.3}$$

*Proof.* Fubini: $\int_0^\infty\nu([0,s])s^{-2}ds=\int\big(\int_\lambda^\infty s^{-2}ds\big)d\nu(\lambda)=\int\lambda^{-1}d\nu(\lambda)$. $\square$

**Hypothesis (M) (infrared spectral leakage).** There exist $A<\infty$, $\delta>0$, independent of $L,\beta,k$, such that for all $\beta\in[1/48,2.5]$,

$$\boxed{\ \nu_{k,L,\beta}\big([0,s]\big)\ \le\ A\begin{cases} s\,\big(s/\omega_L(k)\big)^{\delta}, & 0<s\le\omega_L(k),\\[2pt] s, & \omega_L(k)<s\le1,\\[2pt] 1, & s>1 .\end{cases}\ }\tag{M}$$

**Theorem 5.3 ((M) implies (T)).** Under (M),

$$\widehat C_L(k)\le A\Big(1+\frac1\delta+\log\frac1{\omega_L(k)}\Big)\le A\Big(1+\frac1\delta+2\log L\Big),\tag{5.4}$$

and hence $\|\Gamma_L\|_{2\to2}\le C\Lambda_L$, i.e. (T) holds.

*Proof.* Insert (M) into (5.3) and split at $\omega_L(k)$ and $1$:
$$\int_0^{\omega}s\cdot(s/\omega)^\delta s^{-2}ds=\omega^{-\delta}\int_0^\omega s^{\delta-1}ds=\tfrac1\delta,\quad
\int_\omega^1 s\cdot s^{-2}ds=\log\tfrac1\omega,\quad
\int_1^\infty 1\cdot s^{-2}ds=1 .$$
Summing gives (5.4). Since $\omega_L(k)\ge L^{-2}$ (attained at $k=0$), $\log(1/\omega_L(k))\le 2\log L$. Finally $\Lambda_L\asymp 1+\log L$. $\square$

---

## 6. (M) is strictly stronger than (T)

We record honestly what the reduction of §5 does **not** achieve.

**Proposition 6.1.** (M) implies (T). The converse fails.

*Proof of the second claim.* Suppose $\nu_{k}$ has an atom of mass $m$ at $\lambda=\varepsilon>0$. Then $\nu_k([0,s])\ge m$ for all $s\ge\varepsilon$, while the right side of (M), evaluated at $s=\varepsilon$, equals $A\varepsilon(\varepsilon/\omega)^\delta\to0$ as $\varepsilon\to0$. Hence (M) forbids **any** atom at the origin.

On the other hand,
$$\widehat C_L(k)=\int\lambda^{-1}d\nu_k\ \ge\ \frac{m}{\varepsilon},$$
so (T) is compatible with an atom of mass $m\lesssim C\varepsilon\log L$ placed at $\lambda=\varepsilon$. Such a configuration satisfies (T) and violates (M). $\square$

**Proposition 6.2.** Let (M$'$) be the weakened hypothesis
$$\int_0^1\frac{\nu_{k}([0,s])}{s^2}\,ds\ \le\ C\log\frac{1}{\omega_L(k)}\qquad\text{for all }k .\tag{M$'$}$$
Then (M$'$) $\Longleftrightarrow$ (T) (up to the constant in front of $\log L$).

*Proof.* (M$'$) $\Rightarrow$ (T) is Theorem 5.3 with the pointwise bound replaced by the integral bound. Conversely, at $k=0$ one has $\omega_L(0)=L^{-2}$, and (5.3) gives $\widehat C_L(0)=\int_0^1\nu_0([0,s])s^{-2}ds+\int_1^\infty\nu_0([0,s])s^{-2}ds$, whose second term is $\le\|D_0\|^2\le24$ by (4.1). So (T) at $k=0$ is exactly (M$'$) at $k=0$. $\square$

**Corollary 6.3.** Between (M) and (T) there is **no room for improvement**: the sharpest hypothesis of this shape which implies (T) at $k=0$ coincides with (T) at $k=0$.

> **Interpretation.** The machine of §3–5 is a correct **reformulation** of (T) as a momentum-resolved spectral statement. It is not a reduction of the difficulty of (T). Its content is diagnostic, and is summarised in §7.

---

## 7. What the reduction does achieve: the local-Hessian red herring

Consider the local model consisting of a single fast variable $x\in S^3$ coupled to two boundary variables $y_1,y_2\in S^3$ through $H(x;y_1,y_2)=\beta_Wx^{\mathsf T}(y_1+y_2)$. With $y_1=y_2=n$, $|n|=1$, the conditional law is the four-dimensional von Mises–Fisher distribution $\mu\propto e^{tn^{\mathsf T}x}d\sigma(x)$, $t=2\beta_W$, and

$$Z(t)=\frac{2I_1(t)}{t},\qquad \mathbb E_\mu|P_{n^\perp}X|^2=3\,\frac{I_2(t)}{tI_1(t)} .$$

**Proposition 7.1.** For $\beta_W=2.5$ ($t=5$), the "one-step replicated weight"
$$g:=\frac{\beta_W^2}{2}\mathbb E_{\mu\otimes\mu}\big|P_{n^\perp}(X-X')\big|^2$$
equals
$$\boxed{\ g(n,n)=\frac{3\beta_W^2}{t}\cdot\frac{I_2(t)}{I_1(t)}=3\cdot6.25\cdot0.14386812=2.697527>1\ }\tag{7.1}$$
Moreover $\operatorname{Hess}(-\log Z)$ splits exactly as
$$\underbrace{\beta_W^2F''(t)(\hat h\,\cdot)(\hat h\,\cdot)}_{\text{along }\hat h\ (\text{rank }1),\ =0.3184}\ +\ \underbrace{\frac{\beta_WF'(t)}{r}\big[I-\hat h\hat h^{\mathsf T}\big]}_{\perp\hat h\ (\text{rank }2),\ =2.6975},\qquad F(t)=\log\frac{2I_1(t)}{t},\tag{7.2}$$
and the value (7.1) is produced entirely by the **second** term. $\square$

**Consequence (the red herring).** (7.1) shows that the *absolute-value envelope* of the covariance is not a subcritical activity, and that no path-breeding scheme built on it can converge. It does **not** show that (T) fails, because (T) constrains only $\langle D_k,\mathscr H^{-1}D_k\rangle$ on the scalar channel — a quantity that can remain small even when a local Hessian eigenvalue exceeds $1$.

**Methodological point.** The order of operations in §5 is:
$$\text{project onto the scalar source } D_k\ \longrightarrow\ \text{then take the spectral bound},$$
and **not** $\text{bound the worst local eigenvalue}\longrightarrow\text{then sum}$. This is the essential difference between this machine and the ones it replaces.

---

## 8. Numerical evidence

Monte Carlo at $\beta=2.5$, $\mathrm{SU}(2)$ Wilson action, $L=6,8,10,12$. Here $\chi:=\widehat C_L(0)=L^{-4}\sum_{x,y}\Gamma_L(x,y)$.

| $L$ | $\chi$ | $\Lambda_L$ | $\chi/\Lambda_L$ | $1/(\chi/\Lambda_L)$ |
|---|---|---|---|---|
| 6 | 0.03192553 | 7.7269 | 0.004132 | $242\times$ |
| 8 | 0.03295298 | 10.2988 | 0.003200 | $313\times$ |
| 10 | 0.03428803 | 12.6480 | 0.002711 | $369\times$ |
| 12 | 0.03269092 | 14.7840 | 0.002211 | $452\times$ |

A second independent run at $L=6$ (different seed, shorter) gives $\chi=0.03059447$.

Shell sums $S_j:=\sum_{\|z\|=j}|\Gamma_L(0,z)|\,g_c(z)$ satisfy
$$\sup_{L,j}S_j\approx7.43\times10^{-3},\qquad \text{ratios for }L=6,8,10,12:\ 1.006,\ 1.010,\ 0.982,\ 0.925,$$
i.e. $S_j$ is $L$-independent. Shells with $j>L/2$ are contaminated by the minimum-image convention and must be excluded from fits; on clean shells $j\le L/2$ the fitted effective power is $p_{\rm eff}=4.5\text{–}5.1$.

**Two quantitative readings.**

1. **$\chi$ is $L$-independent while $\Lambda_L\sim\ln L$**, so $\chi/\Lambda_L\searrow0$. The margin **grows** with $L$. The target (T) is therefore satisfied with a wide and increasing margin.

2. **Effective spectral scale.** By (4.1), $\nu_0$ has total mass $\le24$; by (5.2), $\chi=\int\lambda^{-1}d\nu_0$. Hence the $\nu_0$-average of $\lambda^{-1}$ is
$$\frac{\chi}{\|D_0\|^2}=\frac{0.032}{24}=1.33\times10^{-3}\quad\Longrightarrow\quad \lambda_{\rm eff}\approx 750\gg1 .$$
$$\boxed{\ \text{The }\mathscr H\text{-spectrum against the scalar source lies far outside the infrared; effective scale}\ \sim750\ }$$

**Consequence for (M).** Since $\chi$ is *constant* rather than $\asymp(1/\delta)\log L$, the true spectral leakage is far faster than (M) permits: the data is consistent with an $O(1)$ spectral gap in the scalar channel, i.e. with
$$\nu_0([0,s])\le A s^{1+\delta}\qquad\text{with }\delta\text{ effectively large},$$
and **without** the volume-dependent factor $\omega_L(0)^{-\delta}=L^{2\delta}$.

---

## 9. The remaining problem

Sections 3–5 are complete; §6 shows the hypothesis cannot be weakened into a strictly weaker statement; §8 shows the hypothesis holds with a wide margin. What remains is a single model-specific statement.

**Problem 9.1.** Prove that there are $A,\delta>0$ with

$$\boxed{\ \big\langle D_0,\ \mathsf E_{\mathscr H}\big([0,s]\big)D_0\big\rangle\ \le\ A\,s^{1+\delta}\qquad \text{for all }s\in(0,1],\ \text{uniformly in }L\text{ and }\beta\in[1/48,2.5]\ }\tag{9.1}$$

i.e. the version of (M) with the volume factor $\omega_L(0)^{-\delta}$ removed. Equivalently: the scalar source sees no infrared spectral mass.

**Why this is the sharp form.** By (3.2),
$$\mathscr H=\nabla^*_\mu\nabla+\underbrace{\operatorname{Ric}}_{\succeq0}+\underbrace{\operatorname{Hess}S_{\beta,L}}_{\text{indefinite}} .$$
The only possible source of low $\mathscr H$-spectrum is the negative part of $\operatorname{Hess}S_{\beta,L}$ (the unstable/Savvidy-type modes). Numerically these modes exist (a measured $7.1\%$ negative fraction in a related diagnostic) but do not produce low $\mathscr H$-spectrum against $D_0$ — that is precisely what $\lambda_{\rm eff}\approx750$ expresses.

$$\boxed{\ \textbf{Problem 9.1} \equiv \text{negative modes of }\operatorname{Hess}S\text{ do not couple to the plaquette-scalar channel in the infrared.}\ }$$

**Proposed technical route** (each step is a concrete computation, not an estimate):

1. Fix a gauge-invariant **maximal-tree (forest) gauge**, removing the redundant gauge directions of Remark 4.3;
2. group three adjacent free edge variables into a **plaquette block**, so that no single edge is bounded in absolute value;
3. integrate out the block variable explicitly, retaining the full four-dimensional von Mises–Fisher Hessian (as in (7.2));
4. **contract with $D_k$ first**, then take the spectral bound — in that order;
5. use the plaquette **Bianchi constraint** together with $\mathrm{SU}(2)$ **intertwiner** identities to convert the result into (9.1).

**Remark 9.2.** (4.2) alone is insufficient for step 1: the Ward identity removes pure-gauge directions only, and the transverse projection must be quantified in step 3. $\square$

**Remark 9.3 (falsifiability).** If for some sequence $k_L$ and $s_L\to0$ the ratio
$$\frac{\nu_{k_L,L,\beta}([0,s_L])}{s_L^{1+\delta}}$$
diverges as $L\to\infty$, Problem 9.1 fails. Whether that also forces (T) to fail is a separate question, decided by (5.3). The numerical value $\lambda_{\rm eff}\approx750$ is the first test and is favourable. $\square$

---

## 10. What is not claimed

1. **Not claimed:** that (M), or Problem 9.1, holds at $\beta=2.5$; §8 is numerical evidence, not proof.
2. **Not claimed:** that (T) implies any statement about the physical spectrum of a transfer matrix, or that a spectral gap of $\mathscr H$ on the full $1$-form space holds. Only the scalar channel, against the specific source $D_k$, is addressed.
3. **Not claimed:** that the reduction of §5 lowers the difficulty. By Proposition 6.2 it does not.
4. **Not claimed:** that the numerics at $L\le12$ extrapolate. The $\ln L$ growth of $\chi$ is not resolvable against the $4.26\%$ seed spread, and $L=12$ lies below $L=10$.
5. **Not claimed:** that $\lambda_\perp=2.6975$ is harmless in every scheme. It is harmless **for the quantity (T)**; it correctly rules out absolute-value path-breeding schemes.
6. **Not claimed:** that the local model of §7 is a subsystem of the $d=4$ theory in a sense that transfers estimates. It is used only to exhibit the mechanism of (7.2).
7. **Not claimed:** that $\chi_L=O(1)$ implies $\|\Gamma_L\|_{2\to2}=O(1)$. See §11.1 — it does not.
8. **Not claimed:** that the effective spectral scale of §8.2 is a lower bound on the low-lying spectrum. See §11.2 — it is an upper bound on a harmonic mean.

---

## 11. Three logical gaps closed

### 11.1 Zero-momentum data does not bound the operator norm

Recall $\|\Gamma_L\|_{2\to2}=\max_k\widehat C_L(k)$, whereas §8 reports only $\chi_L=\widehat C_L(0)$.

**Counterexample (within the range $q_x\in[-1,1]$).** Let $\xi$ and $\{\eta_x\}$ be mutually independent Rademacher variables and set
$$q_x=a\,\xi(-1)^{x_1}+b\,\eta_x,\qquad a+b\le1 .$$
Then $q_x\in[-1,1]$ and
$$\Gamma(x,y)=a^2(-1)^{x_1+y_1}+b^2\mathbf 1_{\{x=y\}} .$$
For $L$ even, $\widehat C_L(0)=b^2$ is independent of $L$, while at $k_*=(\pi,0,0,0)$
$$\widehat C_L(k_*)=a^2L^4+b^2 .$$
So $\chi_L$ may be held constant while $\|\Gamma_L\|$ grows like $L^4$.

**Consequence.** The numerics of §8 support $\chi_L=O(1)$ and nothing more. The operator-norm statement requires either the full Fourier data (Problem 11.4) or the pointwise decay of Problem 11.5.

### 11.2 The effective spectral scale is a harmonic mean

Let $M_0:=\nu_0((0,\infty))=\|D_0\|^2\le24$ and $\lambda_{\rm harm}:=M_0/\chi_L$. From (5.2), $\lambda_{\rm harm}$ is the $\nu_0$-mean of $\lambda^{-1}$, inverted; with only $M_0\le24$ one gets the **upper** bound
$$\lambda_{\rm harm}=\frac{M_0}{\chi_L}\le\frac{24}{0.032}=750 .$$
That is not the location of the low spectrum. Indeed, for any $\varepsilon>0$ put
$$\nu_\varepsilon=\frac{\chi\varepsilon}{2}\,\delta_\varepsilon+\Big(24-\frac{\chi\varepsilon}{2}\Big)\delta_{\lambda_h},\qquad \lambda_h=\frac{2}{\chi}\Big(24-\frac{\chi\varepsilon}{2}\Big).$$
Then $\nu_\varepsilon((0,\infty))=24$ and $\int\lambda^{-1}d\nu_\varepsilon=\chi$ for every $\varepsilon$, so $\lambda_{\rm harm}=750$ throughout, yet
$$\nu_\varepsilon([0,\varepsilon])=\frac{\chi}{2}\,\varepsilon ,$$
which is $O(\varepsilon)$, not $O(\varepsilon^{1+\delta})$. Hence a large $\lambda_{\rm harm}$ is compatible with violation of (9.1).

### 11.3 The zero- versus small-positive-spectrum distinction

Since $\mathscr H\ge0$ (from $\langle\omega,\mathscr H\omega\rangle=\|d_\mu^*\omega\|^2+\|d\omega\|^2$) and $H^1(\mathrm{SU}(2)^{|E_L|})=0$ in finite volume, $\mathscr H$ has **no zero modes**. So the difference between (T) and (M) is not about atoms at $0$; it is:

> (T) permits mass $O(\varepsilon\log L)$ at a small positive $\varepsilon$; (M) requires mass $O(\varepsilon^{1+\delta}\omega_L^{-\delta})$, which is strictly smaller.

Proposition 6.1 stands with this correction.

### 11.4 Problem (Fourier data)

**Problem A.** For $L=6,8,10,12$ and $\beta=2.5$, compute **all** $\widehat C_L(k)$ and report $\arg\max_k\widehat C_L(k)$ and $\max_k\widehat C_L(k)/\widehat C_L(0)$.

If the correlator array $\Gamma_L(0,\cdot)$ is available, $\widehat C_L=\mathrm{FFT}_4(\Gamma_L(0,\cdot))$. From Monte Carlo configurations $q^{(n)}$,
$$\widehat C_L(k)=\frac{1}{L^4}\mathbb E\big|Q_n(k)-\mathbb E Q_n(k)\big|^2,\qquad Q_n(k)=\sum_xq^{(n)}_xe^{ik\cdot x},$$
estimated in practice by $\frac{1}{L^4(N-1)}\sum_n|Q_n(k)-\overline Q(k)|^2$. **The zero mode must be centred across configurations, not within each configuration**, otherwise the genuine zero-momentum fluctuation is removed.

### 11.5 Problem (pointwise decay; the exponent $4$ is insufficient)

$$\|\Gamma_L\|_{2\to2}\le\|\Gamma_L\|_{1\to1}=\sum_{z\in\Lambda_L}|\Gamma_L(0,z)| .$$

**But $|\Gamma_L(0,z)|\le C(1+\|z\|_2)^{-4}$ gives only $O(\log L)$.** One needs
$$\boxed{\ \textbf{(C2)}\qquad |\Gamma_L(0,z)|\le C\,(1+\|z\|_2)^{-4-\varepsilon},\quad\varepsilon>0\ }\tag{C2}$$
which does give absolute summability.

**Why the shell data of §8 does not suffice.** $S_j\le M$ implies only $\sum_{\|z\|=j}|\Gamma_L(0,z)|\le M(1+j)^4$, hence $\sum_z|\Gamma_L(0,z)|\le M\sum_{j\le O(L)}(1+j)^4=O(L^5)$. Moreover, with the counterexample of §11.1, $\Gamma_L(0,z)=a^2(-1)^{z_1}+b^2\mathbf 1_{\{z=0\}}$ gives $\chi_L=b^2$ constant, $S_j\asymp a^2/j$ bounded, and yet $\sum_z|\Gamma_L(0,z)|\asymp a^2L^4$. So (SB) boundedness carries no implication of absolute summability.

**The fitted $p_{\rm eff}=4.5$–$5.1$ does support (C2)** (with $\varepsilon\approx0.5$–$1.1$), but it comes from fitting $\Gamma_L(0,z)$ against $\|z\|$, not from $S_j$-boundedness.

---

## 12. Negative part of $\operatorname{Hess}S$: exact local computation

This section records what can be proved unconditionally about the negative part, and identifies precisely why it does not close Problem 9.1.

### 12.1 Exact single-plaquette Hessian

Use the bi-invariant metric $\langle X,Y\rangle=-\tfrac12\operatorname{tr}(XY)$ on $\mathfrak{su}(2)$, so that $X(x)=ix\cdot\sigma$ has $|X(x)|=|x|$, and $\operatorname{Ric}=2I$.

For a plaquette $U_P=qI+ia\cdot\sigma$ with $q^2+|a|^2=1$, parallel-transport the edge tangents to a common base point, obtaining $y_1,\dots,y_4\in\mathbb R^3$ (reversed edges carry a sign and an adjoint). Then $U_P(t)=e^{tY_1}\cdots e^{tY_4}U_P$ with $Y_i=iy_i\cdot\sigma$, and using $(ix\cdot\sigma)(iy\cdot\sigma)=-(x\cdot y)I-i(x\times y)\cdot\sigma$,

$$\frac{d^2}{dt^2}\Big|_{t=0}\tfrac12\operatorname{tr}U_P(t)=-q\Big|\sum_{i=1}^4y_i\Big|^2+2\sum_{i<j}(y_i\times y_j)\cdot a .$$

Since $S_{\beta,L}=\beta\sum_P(1-q_P)$,

$$\boxed{\ Q_P(Y)=\beta q\Big|\sum_{i=1}^4y_i\Big|^2-2\beta\sum_{i<j}(y_i\times y_j)\cdot a\ }\tag{12.1}$$

an identity, not an estimate. It exhibits the indefiniteness directly.

### 12.2 An explicit configuration with $V_-\ne0$

**Proposition 12.1.** Take the configuration with $U_e=I$ for all $e$ except a single edge $e_0$ with $U_{e_0}=-I$. Let $X$ be a unit tangent vector at $e_0$ and $0$ elsewhere. Then
$$\langle X,\operatorname{Hess}S\,X\rangle=-6\beta|X|^2 .\tag{12.2}$$
Hence $\|V_-\|\ge6\beta$, and $\|V_-\|\ge15$ at $\beta=2.5$.

*Proof.* The six plaquettes containing $e_0$ have holonomy $-I$, so $q=-1$, $a=0$; the other plaquettes do not move. Since $\operatorname{Ad}(-I)=I$, all six contribute with the same sign. By (12.1), each gives $Q_P=-\beta|X|^2$. $\square$

**Local criterion.** If a plaquette has $q_P<0$ and one takes $y_1=y_2=y_3=y_4=v$ (all cross products vanish), then $Q_P=16\beta q_P|v|^2<0$. So negative curvature occurs away from the $-I$ configuration as well.

### 12.3 Uniform bounds

From $|\sum y_i|^2\le4\sum|y_i|^2$ and $\sum_{i<j}|y_i||y_j|\le\tfrac32\sum|y_i|^2$,
$$|Q_P(Y)|\le\beta\big(4|q|+3|a|\big)\sum_i|y_i|^2\le5\beta\sum_i|y_i|^2 ,$$
using $4|q|+3|a|\le\sqrt{16+9}=5$ on $q^2+|a|^2=1$. Summing over the six plaquettes per edge,

$$\boxed{\ \|V\|\le30\beta,\qquad\|V_-\|\le30\beta\ }\tag{12.3}$$

and, with $\dim T_U\mathcal M_L=3\cdot4L^4=12L^4$,

$$\boxed{\ \operatorname{tr}_{\rm fiber}V_-\ \le\ 30\beta\cdot12L^4\ =\ 360\beta L^4\ }\tag{12.4}$$

At $\beta=2.5$: $\|V_-\|\le75$ and $\operatorname{tr}_{\rm fiber}V_-\le900L^4$. Combining (12.2) and (12.3):
$$6\beta\ \le\ \|V_-\|\ \le\ 30\beta .$$

**Warning.** $\operatorname{tr}_{\rm fiber}$ is the trace on the finite-dimensional tangent space at a single configuration; it is **not** the trace of the multiplication operator on $L^2$.

### 12.4 Source-weighted Birman–Schwinger form

Write $W:=V_-$, $V=V_+-W$, and $A_+:=\nabla_\mu^*\nabla+\operatorname{Ric}+V_+\ge2I$, so that

$$\mathscr H=A_+-W,\qquad A_+\ge2I .\tag{12.5}$$

For a normalized eigenvector $\psi_j$ with $0\le\lambda_j<2$, put $\phi_j:=W^{1/2}\psi_j$. Then $(A_+-\lambda_j)\psi_j=W\psi_j$, hence

$$\boxed{\ W^{1/2}(A_+-\lambda_j)^{-1}W^{1/2}\phi_j=\phi_j\ }\tag{12.6}$$

and for the zero-momentum source

$$\boxed{\ \nu_0([0,s])=\sum_{\lambda_j\le s}\Big|\big\langle W^{1/2}(A_+-\lambda_j)^{-1}D_0,\ \phi_j\big\rangle\Big|^2\ }\tag{12.7}$$

Since $A_+-\lambda_j\ge(2-s)I$ for $\lambda_j\le s$, one obtains the bound

$$\boxed{\ \nu_0([0,s])\ \le\ \frac{\|W\|\,\|D_0\|^2}{(2-s)^2}\sum_{\lambda_j\le s}\|\phi_j\|^2\ \le\ \frac{720\beta}{(2-s)^2}\sum_{\lambda_j\le s}\langle\psi_j,W\psi_j\rangle\ }\tag{12.8}$$

and at $\beta=2.5$, $s\le1$:
$$\nu_0([0,s])\ \le\ 1800\sum_{\lambda_j\le s}\langle\psi_j,W\psi_j\rangle .\tag{12.9}$$

**This converts the low-spectral mass exactly into the negative-Hessian overlap of the low eigenfunctions.** It does **not** close the problem, for two independent reasons.

### 12.5 Why the Birman–Schwinger / CLR route is closed

**(a) Fiber trace does not control the global trace.** The right side of (12.9) is $\operatorname{Tr}_{L^2}(P_sW)$ with $P_s$ the low spectral projector, a quantity depending on how eigenfunctions are distributed over the whole configuration manifold. The bound $\operatorname{tr}_{\rm fiber}W(U)=O(L^4)$ is a statement about a finite-dimensional matrix at each $U$. There is no inequality
$$\operatorname{Tr}_{L^2}(P_sW)\ \le\ C\,\operatorname{tr}_{\rm fiber}W ,$$
and none can hold: a single highly concentrated eigenfunction contributes $O(1)$ to the left side while the right side is averaged over configurations.

**(b) The spectral-count bounds diverge.** Cwikel–Lieb–Rozenblum / Lieb–Thirring on $\mathcal M_L$ would give $\#\{\lambda_j<2\}\le C\int\operatorname{tr}(V_-^{d/2})$ with $d=\dim\mathcal M_L=12L^4\to\infty$; and the Birman–Schwinger trace $\operatorname{Tr}_{L^2}\big(W^{1/2}(A_+-s)^{-1}W^{1/2}\big)$ is bounded only by $\|W\|(2-s)^{-1}\operatorname{Tr}_{L^2}(\mathbf 1)=\infty$. Both are vacuous.

**(c) The compensation requested in §9 cannot be obtained.** What would suffice is $\sum_{\lambda_j\le s}\langle\psi_j,W\psi_j\rangle=O(s^{1+\delta})$ uniformly, equivalently a vanishing multiplicity of the low $\mathscr H$-spectrum. Even assuming delocalization of low eigenfunctions gives only $O(1)$ per $j$ and hence $O(\#\{\lambda_j\le s\})$ in total, so one needs $\#\{\lambda_j\le s\}=O(s^{1+\delta})$ — a **gap statement on the full $1$-form space**, which is stronger than (T).

> **Conclusion of §12.** The local geometry of the negative part of $\operatorname{Hess}S$ is now closed: (12.1)–(12.4) are exact, and they show that the negative part is a first-order effect ($\|V_-\|\ge6\beta$, i.e. $\ge15$ at $\beta=2.5$), not a small perturbation. But pushing this route to its end requires control of the low-$\mathscr H$ eigenvalue **multiplicity**, which is strictly stronger than the plaquette-scalar target. **The negative-Hessian route is therefore a detour, not an approach.**

---

## 13. Revised problem list

| label | problem | relation to (T) |
|---|---|---|
**A** | compute all $\widehat C_L(k)$, decide $\arg\max_k$ | **decisive** test; needs raw data |
**B** | $\nu_0([0,s])\le As^{1+\delta}$ | implies $\chi_L=O(1)$ only |
**C** | $|\Gamma_L(0,z)|\le C(1+\|z\|)^{-4-\varepsilon}$ | **implies (T) with $O(1)$**, stronger than needed |
**D** | low-$\mathscr H$ multiplicity | **stronger than (T)**; route closed in §12.5 |
**E** | $\sup_{L,\beta,k}\int_0^1\nu_{k}([0,s])s^{-2}ds<\infty$ | **equivalent** to $\|\Gamma\|=O(1)$ |

**Recommended order.** A (cheapest, decisive) → C (direct route to $O(1)$, skips all spectral machinery) → B (weakest layer) → E (the exact spectral form) → D (closed).

---

## Appendix A. Why the four earlier approaches fail (with numbers)

Four schemes were tested against (T) and are excluded. None of them is a counterexample to (T).

| scheme | resulting threshold | reason for failure |
|---|---|---|
| worst-case $\sup$ of the local Hessian, one-edge coefficient $\tanh(\pi^2\beta/4)$ | $\beta<0.050$ | uses $\sup_y$ worst case; (T) is a $2\to2$ norm |
| absolute-value envelope as path weight | §7: $g=2.697527>1$ | envelope is not an activity; confirmed by (7.1) |
| backward "ancestry" breeding coefficient $\theta_r$ | **structurally impossible:** $\theta_r\ge1$ for all parameters | an inconsistent state cannot appear without a source, so the conditional probability of at least one $r$-step extension is $1$ |
| signed vector-valued response operator | $\beta_*=2.244977$ | this measures the transverse (angular) response, a different object from (T) |

The last line is the closest: $\beta_*=2.244977$ gives $2.5/\beta_*=1.11$. But it is the transverse channel; the obstruction there does not transfer to the scalar channel of (T).

## Appendix B. Free-field and $d=2$ checks

- Free massless scalar, $q\sim(\partial\phi)^2$: $\operatorname{Cov}(q_0,q_z)\sim g^4\|z\|_2^{-8}$, which is integrable in $d=4$, so (T) holds.
- In $d=2$ the analogous $q_x$ does not propagate ($\Gamma\equiv0$), so (T) holds trivially.

Both models satisfy (T), so neither discriminates (T). This is a known weakness of the target and is recorded here rather than hidden.

---

*All identities in §3–§5 and Propositions 6.1, 6.2, 7.1 are proved above. Problem 9.1 is open.*
