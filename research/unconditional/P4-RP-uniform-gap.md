# P4 — Certificate-free G2: reflection positivity + transfer matrix, volume-uniform gap

**Verdict.** Lattice RP + transfer matrix make the spectral reformulation rigorous and certificate-free, but RP yields no gap, and every known lower-bound technique needs the one uniform input that *is* the mass gap. The route is valid and non-circular but relocates the IR problem without weakening it.

## 1. What OS'78 gives, exactly
Osterwalder–Seiler (Ann. Phys. **110** (1978) 440) prove the Wilson measure is reflection positive and reconstruct $(H,T,\Omega)$: $T=e^{-aH}$, $0\le T\le1$, $H\ge0$, $T\Omega=\Omega$; on a finite box, discrete spectrum $1=\lambda_0>\lambda_1\ge\cdots$; KL holds for gauge invariants; $\Delta_V=-\tfrac1a\ln\lambda_1=\min\operatorname{supp}\nu$. RP is a *sign* inequality $\langle\theta F,F\rangle\ge0$; it fixes signs, **not** the size of $1-\lambda_1$. OS prove no gap.

## 2. Finite volume upper-bounds only
$\Delta_{\mathbb R^4}=\inf_V\Delta_V$ ($V\to\infty$ limit, free BC; min–max lets larger boxes hold more trial states). Each finite-$V$ gap is an *upper* bound; a positive $\Delta_{\mathbb R^4}$ needs $\inf_V\Delta_V\ge m$, a volume-uniform **lower** bound — and that infimum is the gap itself. "$\Delta_L>0\ \forall L$" is insufficient: it would "prove" uniform clustering for subcritical 2D Ising, where $\Delta_L\sim e^{-\sigma L}\to0$ (F9-1).

## 3. Techniques for a uniform lower bound
- **Variational / Wilson-loop:** min–max gives *upper* bounds on excited energies; loop states bound $\Delta$ above, never below.
- **Gap monotonicity:** compression and RG coarse-graining ($\Delta_{\rm coarse}\ge\Delta_{\rm fine}$) run the wrong way — upper bounds.
- **Perron–Frobenius:** finite-$V$ primitivity gives $1-\lambda_1>0$, but quantitative PF/Dobrushin constants scale like $e^{-cV}$.
- **Martingale / finite-size criterion (Nachtergaele 1996; survey Nachtergaele–Sims, *Markov Proc. Rel. Fields* **12** (2006)):** the one genuine uniform lower-bound machine — finite-block gap $\ge\gamma$ plus local gap stability implies a bulk gap. Nachtergaele–Sims, *CMP* **265** (2006) 119 is the converse (gap ⇒ uniform Lieb–Robinson clustering), one-way. Known successes: strong-coupling spin models; Shen–Zhu–Zhu (*arXiv:2204.12737*), the sole 4D gauge uniform gap, also strong-coupling.

## 4. Direct route, no certificate?
Yes logically: the uniform Poincaré bound $\operatorname{Var}_\mu(F)\le\tfrac{a}{m}\sum_x\mu(|\partial_xF|^2)$ is exactly $\lambda_1\le e^{-am}$, i.e. $\Delta_{\mathbb R^4}\ge m$ — no expansion, no $\beta\approx2.47$ crossover. But this inequality *is* the mass gap in measure language; not weaker.

## 5. The single missing lemma
> **L\*.** $\exists m>0$: uniform Poincaré constant $a/m$ for all $a,V$ at critical coupling — equivalently a finite-size criterion at physical coupling.

**Not provable with present tools.** Volume-uniform gaps exist only at strong coupling (polymer/Dobrushin; OS'78; SZZ), in $d=2$ (Glimm–Jaffe–Spencer), or solvable/large-$N$ corners. At 4D weak coupling the uniform constant is blocked by the same IR gluon-mass / Savvidy-concavity wall (W3-D). L\* holds iff YM is gapped — it *is* the Clay problem, certificate-free.

**Sharpest line.** RP+TM convert a uniform lower bound into the mass gap crossing no crossover; producing that bound equals the mass gap and is open.
