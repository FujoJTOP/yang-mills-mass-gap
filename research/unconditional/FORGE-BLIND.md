# FORGE-BLIND
### A blind, first-principles attempt at a cutoff-uniform gap-lower-bound engine across an indefinite Hessian

## 1. The obstruction: three faces, one hole in the toolbox

A mass-gap lower bound $m\ge m_0$ is the statement that the self-adjoint excitation operator $\mathcal A = H-E_0$ has **no spectrum in $(0,m_0)$**. This is *spectral stability*, and every robust stability certificate rests on a positive quadratic form: completing the square, Cauchy-Schwarz, reflection, min-max on a semibounded form. At the Savvidy point the fluctuation form is indefinite — $A=-D^2+2\,\mathrm{ad}F$ has negative eigenvalues — so the toolbox splits and both halves fail:

- **Direction (right hypothesis, wrong sign).** Tools that survive indefiniteness conclude the wrong thing. Reflection positivity yields only $H\ge 0$: a *zero* gap guaranteed, value $0$. Moment / Chebyshev-Markov / Stieltjes bounds need a *positive* spectral measure; the saddle carries a *signed* (Krein) measure, so every bound reverses or goes vacuous. Index / spectral-flow bounds need a spectral gap merely to *define* the Fredholm index (the gap condition) and to protect it topologically; across the Savvidy crossing the relevant index vanishes, so they certify a **no-go**, never a lower bound.
- **Positivity (right direction, wrong hypothesis).** Log-Sobolev, Brascamp-Lieb, and cluster expansions *do* point the right way ($\lambda_1\ge 1/C_{LS}$; decay rate $\ge$ rate), but each requires a convex potential / positive Hessian: $\Gamma_2\ge\rho\Gamma$, Hess $>0$, a positive covariance. That hypothesis is false *exactly* in the indefinite region, and convergence radii shrink precisely as one approaches it. These engines prove gaps only where the quadratic reference is already positive.
- **Locality.** A certificate uniform in cutoff and volume must be local and gauge-covariant. But the negative mode is *delocalized* (a plane / lowest-Landau mode), and its negative curvature is *extensive*, $\propto VF$. A sum of local positive terms cannot dominate an extensive negative term without a global density input.

The indefinite region is the no-man's-land between "right direction, wrong hypothesis" and "right hypothesis, wrong direction."

## 2. The tool: a Feshbach-Bogoliubov (FB) gap certificate

The opening is **two-block positivity**: we never need $\mathcal A\ge0$ from a positive Hessian — only positivity on the *small* tachyonic block, plus a gap on its complement. Split $\mathcal H=\mathcal H_-\oplus\mathcal H_+$, $\mathcal A=\begin{pmatrix}A&B\\B^*&C\end{pmatrix}$, where $\mathcal H_-$ is the isolated negative/near-zero spectral subspace of the quadratic Savvidy operator (the lowest Landau band, $\dim\mathcal H_-\propto VF$) and $\mathcal H_+$ its complement, on which the quadratic operator has a genuine gap $\gamma_0\sim F>0$ (the negative band is isolated from the $\gtrsim F$ continuum at every field strength). $A$ may be arbitrarily negative; this is allowed.

**Lemma (two-block gap; rigorous).** If $C\ge\gamma I$ and the *Schur complement* $S:=A-BC^{-1}B^*\ge\mu I$ on $\mathcal H_-$, then
$$\mathcal A\ \ge\ \min\!\Big\{\tfrac{\mu}{1+2\|B\|^2/\gamma^2},\ \tfrac{\gamma}{2}\Big\}\,I\ >\ 0 .$$
*Proof.* $\langle u,\mathcal Au\rangle=\langle u_1,Su_1\rangle+\gamma\,\|u_2+\gamma^{-1}B^*u_1\|^2$; put $v_2=u_2+\gamma^{-1}B^*u_1$ and use $\|u_2\|^2\le 2\|v_2\|^2+2\|B\|^2\gamma^{-2}\|u_1\|^2$. $\square$

So the full indefinite operator is positive with an *explicit* retained gap **iff the Schur complement on the tachyonic block is positive.** In Bogoliubov (squeezing) language, $S$ is the curvature of the squeezed vacuum, and the lemma says that vacuum has no negative-frequency mode. The indefinite Hessian is not hidden — it is *converted* into the demand $S\ge\mu>0$.

## 3. The single missing inequality

With $A|_{_-}\approx-cF<0$, positivity of $S$ requires the **interacting** contributions to overcompensate:
$$\underbrace{\langle\phi,\mathcal Q\phi\rangle}_{\text{quartic 2nd variation}}+\underbrace{\langle\phi,\mathcal K\phi\rangle}_{\text{Feshbach lift}}\ -\ cF\ \ge\ \mu>0\qquad\text{for all }\|\phi\|=1,\ \phi\in\mathcal H_-. \tag{TCI}$$
This is the **tachyonic-compensation inequality (TCI)**. For delocalized $\phi$ it reduces, by translation averaging, to a *density* inequality; the covariant split keeps $\gamma,\|B\|$ UV-finite. **The obstacle is localization:** because $\mathcal H_-$ is extensive, TCI must hold for *localized wave packets probing the worst cell* — equivalently, that the gauge-invariant effective potential is **convex at its minimum, pointwise and nonperturbatively**, and that the Savvidy vacuum is stable (no imaginary part). That, and only that, is where the certificate can fail.

## 4. Ledger: proved / conditional / open

- **Proved (unconditional, abstract).** The two-block lemma: gap of $\mathcal A$ $\ge$ an explicit function of (Schur gap $\mu$, complement gap $\gamma$, coupling $\|B\|$), with **no hypothesis of a positive Hessian**. Uniform in cutoff and volume *provided* $\gamma,\|B\|,S$ are.
- **Conditional.** For 4D pure SU(N): $m\ge\min\{\mu/(1+2\|B\|^2/\gamma^2),\,\gamma/2\}$, conditional on **TCI** holding uniformly (worst cell, all $N$, all field strengths $F$).
- **Acid test (free massless field).** No $F$, no quartic $\Rightarrow\mathcal Q=0$, $A|_-$ has $0$ at its bottom, $S\not\ge\mu$; TCI fails ($0\not\ge cF+\mu$) and the certificate returns $0$. Correct: gapless. The tool is therefore *not vacuous* (it must fail on free fields and does) and is *not* reflection positivity (it never assumes a positive Hessian).
- **Open / precise stuck point.** Prove TCI pointwise and uniformly. What the tool *achieves* is a sharp reduction: the entire gap-lower-bound is equivalent to **convexity-at-the-vacuum of the nonperturbative effective potential plus absence of the Savvidy imaginary part**, paid for by one explicit leakage constant $\|B\|^2/\gamma^2$. Reflection positivity, LSI, and cluster expansions cannot even *state* this reduction; the FB certificate can — its only debt is TCI.
