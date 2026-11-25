# X2-FACTORIZATION — the mass gap as (the failure of) local constancy of the observable factorization algebra: a **sign audit**, an **exact solvable-model obstruction**, and the **kill verdict** (sixth species, or the transfer matrix renamed?)

> 2026-11-24 · First step of the **factorization-algebra / observables** paradigm. The proposed
> definition under test is *“the mass gap is the failure of local constancy of the factorization
> algebra of observables.”* This memo does three things and stops honestly: **(1)** it pins the
> **sign** of the definition against the actual Costello–Gwilliam theorem and the project's own
> acid-test convention; **(2)** it **computes the obstruction exactly** in the solvable cases —
> free scalar in $d=1,2,3,4$, 2D Yang–Mills, 3D topologically massive gauge theory; **(3)** it renders
> the kill verdict: is this a new object, or is it weak-$*$/H1/W/I/the transfer matrix in homotopical
> dress?
>
> **Labels:** `[T]` 真定理 · `[C]` 条件定理 · `[NG]` no-go for a *named class* · `[O]` 精确缺失步 ·
> `[V]` verified-here (numerical witness run this session) / bibliographic-verified in
> `TOOL-FACTORIZATION-ALGEBRA.md`. **No QED; Clay 0%.** This memo proves no gap and no
> global impossibility; `∀class ≠ ∀possible`.
>
> **Do not double-count.** The scouting memo `research/unconditional/TOOL-FACTORIZATION-ALGEBRA.md`
> already established the *two-tier blindness lemma* (deformation-invariant derived invariants are
> mass-blind; the only jump is the IR limit $=$ restatement). X2 is the **focused first step**:
> the **sign audit** (missing there, decisive here), the **explicit 2D-YM decoupling computation**,
> and the **tier-mapping onto H1/W/I**.

---

## §0 · Verdict up front

> The literal Costello–Gwilliam (CG) theorem is
> $$\textbf{mass gap}\iff \textbf{the observable factorization algebra is locally constant}$$
> (free *massive* scalar $=$ locally constant; free *massless* $=$ **not**). The task's
> “mass gap $:=$ **failure** of local constancy” is the **logical negation** of this theorem, and it
> is falsified by the free-field anchor and, sharply, by 2D Yang–Mills.

Four findings, each labeled:

1. **Sign audit `[T]` `[V]`.** Under the precise CG/Ayala–Francis notion, *local constancy* (LC) holds
   for the free **massive** field and fails for the free **massless** one. The obstruction to LC is the
   **propagator's IR decay rate**, which equals the mass. Hence “mass gap $=$ failure of LC” has the
   wrong sign; “mass gap $=$ **LC**” is the theorem. The only reading that makes the task's stated
   acid test (*“free massless locally constant; massive/YM not”*) literally true is
   *LC $:=$ scale-invariance / vanishing trace anomaly* — a **different** object (§1.3), which
   reduces to the **gluon condensate** $\langle\tfrac14\mathrm{tr}F^2\rangle$ (already a banked name of
   the same open target, `CLAUDE.md` §5), hence circular.

2. **2D Yang–Mills obstructions it `[T]/[C]` `[V]`.** 2D YM has **no local propagating degrees of
   freedom** (all curvature is a top form; observables are holonomy data), so its observable FA is
   **locally constant for topological reasons** (defect $\equiv 0$) — *yet* its transfer matrix is
   **gapped**, gap $=\tfrac{e^2}{2}\min_{R\neq0}C_2(R)>0$ ($3e^2/8$ for $SU(2)$, $2e^2/3$ for
   $SU(3)$ fundamental). So **LC and the mass gap are independent axes**; the task's sign predicts 2D
   YM is gapless — **contradicted**.

3. **“Obstruction” $=$ the analytic kernel `[T]` `[V]`.** In every solvable case the object that
   decides LC is literally the **propagator / transfer-matrix eigenvalue**: free scalar $\to m$;
   3D topologically massive $\to\mu$ (the Chern–Simons mass, an **input**). The criterion **renames
   the transfer matrix**; it does not compute a new number.

4. **Kill criterion `[NG]`/`[O]`.** Correctly signed (gap $\iff$ LC), the criterion is the **gap
   restated qualititatively** — strictly **stronger** than H1 (H1 is a quantitative ceiling
   $\chi\le1/\lambda$; LC $\Rightarrow$ clustering $\Rightarrow$ gap $\Rightarrow$ H1-type), and it
   **does not** reach W (atom weight) or I (isolation), because *LC $\iff$ the full gap* and the gap
   does not force an atom/isolation (R-wave: gapped-atomless measures exist). It is **not** a sixth
   species: it sits at the **gap tier** the campaign already knows does not close Clay; it does **not**
   reduce to H1/W/I (different, lower/orthogonal tiers), but it **does** reduce to the transfer
   spectrum. Kill clause “**only renames the transfer matrix**” is **triggered.**

---

## §1 · The sign audit

### §1.1 The precise object and the precise theorem

**Prefactorization algebra `[V]`.** A prefactorization algebra on a manifold $M$ assigns a cochain
complex $\mathcal F(U)$ to each open $U\subset M$ with multiplication maps
$\mathcal F(U_1)\otimes\cdots\otimes\mathcal F(U_n)\to\mathcal F(V)$ for pairwise-disjoint
$U_i\subset V$; it is a **factorization algebra** if it satisfies the local-to-global (Weiss/Čech)
condition. (CG Vol 1 §2–3.) `[V]` K. Costello & O. Gwilliam, *Factorization Algebras in QFT*, Vol 1
(New Math. Monographs **31**, CUP 2017), Vol 2 (**41**, CUP 2021).

**Local constancy `[V]`.** $\mathcal F$ is **locally constant** if for every nested-ball inclusion
$B_r(x)\subset B_R(x)$ ($0<r<R$) the structure map $\mathcal F(B_r(x))\to\mathcal F(B_R(x))$ is a
**quasi-isomorphism**. Equivalently (Ayala–Francis): on $\mathbb R^d$, locally constant
$\iff$ $\mathcal F$ is (up to equivalence) the factorization algebra of an $E_d$-algebra. `[V]`
D. Ayala & J. Francis, *Factorization homology of topological manifolds*, J. Topology **8** (2015)
1045–1084, arXiv:1206.5522.

**The CG free-scalar theorem `[T]`.** For the free scalar of mass $m$ on $\mathbb R^d$:
- $m>0$: $\mathrm{Obs}_m$ **is** locally constant.
- $m=0$ ($d\ge3$): $\mathrm{Obs}_0$ is **not** locally constant.

*(This is the standard CG free-scalar example. For the precise **quantum** observables the IR statement
carries the usual massless IR-divergence subtlety — see Kawahira–Shigemura, who handle it by removing
IR divergences and treat the mass as an input; the *classical/perturbative* statement is clean. The
sign, massive $=$ LC, is not in doubt.)*

*Mechanism (why this is the whole story).* The structure map that extends an observable from $B_r$ to
$B_R$ is built from the **propagator** $P_m=(-\Delta+m^2)^{-1}$. Because $P_m$ decays like
$e^{-m|x|}$, the extra shell $B_R\setminus B_r$ is **cohomologically trivial** — one homotopes it away
along the decaying kernel, so $H^\bullet(\iota)$ is an isomorphism. The massless $P_0$ has a long-range
tail ($|x|^{-(d-2)}$, $\log|x|$, or linear in $d=2,1$); the outer shell carries genuine cohomology (the
harmonic / zero mode), so $\iota$ is **not** a quasi-isomorphism. **The mass is the decay rate of the
kernel, and LC is exactly the statement that this rate is positive.**

### §1.2 The sign, confronted with the free-field anchor `[T]` `[V]`

Witness `w7/x2_factorization.py` (output `w7/x2_factorization_results.json`) extracts the decay rate
$\mathrm{rate}(r)=-\tfrac{d}{dr}\log|G(r)|$ of the exact Euclidean propagators. For $m=1$:

| $d$ | massive rate at $r{=}6,10,14$ | $\to m$ | massless rate at $r{=}6,10,14$ | $\to 0$ |
|---|---|---|---|---|
| 1 | $1.00,\,1.00,\,1.00$ | $1$ | $-0.167,-0.100,-0.071$ ($1/r$) | $0$ |
| 2 | $1.06$ @ $r{=}8$ | $1$ | $-0.060$ @ $r{=}8$ ($\log$) | $0$ |
| 3 | $1.167,\,1.100,\,1.071$ | $1$ | $0.167,\,0.100,\,0.071$ ($1/r$) | $0$ |
| 4 | $1.259,\,1.153,\,1.109$ | $1$ | $0.333,\,0.200,\,0.143$ ($2/r$) | $0$ |

So the exact free propagator gives: **massive $\Rightarrow$ exponential decay $\Rightarrow$ locally
constant $\Rightarrow$ gapped**; **massless $\Rightarrow$ no decay $\Rightarrow$ not locally constant
$\Rightarrow$ gapless.** The obstruction to LC is the mass, literally. **Therefore the assertion
“mass gap $=$ failure of local constancy” is the negation of the truth**: it would declare the free
massless field (no decay, *not* locally constant) to be **gapped** and the massive field (decay,
locally constant) to be **gapless** — backwards in both directions.

**The project's own acid test confirms the same sign `[T]`.** `CLAUDE.md` §1: the tool must **fail
(be $0$/vacuous) on the free massless field** and be **positive on Yang–Mills**. Take *LC itself* as
the detector ($\mathrm{LC}=1\!\Rightarrow$ gapped):

| detector | free massless | YM / gapped | passes `CLAUDE.md` acid test? |
|---|---|---|---|
| **LC** (gap $\iff$ LC, standard sign) | $\mathrm{LC}=0$ (vacuous ✓) | $\mathrm{LC}=1$ (positive ✓) | **passes both ends** |
| **failure of LC** (task's sign) | defect $>0$ (positive on free ✗) | defect $=0$ (zero on YM ✗) | **fails both ends** |

The task's stated acid test (*“free massless locally constant; massive/YM not”*) is exactly the
**negation** of the project's standard 酸测 and of the CG theorem. This is the decisive sign finding:
**use LC, not its failure.**

### §1.3 The one charitable reading that matches the task's acid test — and why it is circular `[C]`

There **is** a reading under which *“free massless locally constant; massive/YM not”* is literally
true: take
$$\text{“local constancy”} \;:=\; \text{scale-invariance on flat space},\qquad
\text{“obstruction”} \;:=\; \langle T^\mu_{\ \mu}\rangle=\frac{\beta(g)}{2g}\langle \mathrm{tr}\,F^2\rangle.$$
Then: free massless (flat) $\to$ anomaly $0$ $\to$ LC $\to$ “gapless” ✓; massive scalar $\to m$-term
anomaly $\to$ not-LC $\to$ “gapped” ✓; 4D YM $\to$ dimensional transmutation / trace anomaly $\to$
not-LC $\to$ “gapped” ✓; 2D YM $\to$ dimensionful $e^2$ $\to$ not-LC. The acid tests line up. But
this reading is **not** the CG object, and it fails as a tool on three counts:

- **(a) Circular in the value.** The anomaly *coefficient* $\beta/2g$ is perturbative (one-loop), but
  the nonzero **value** $\langle T^\mu_{\ \mu}\rangle\ne0$ requires the nonperturbative **gluon
  condensate** $\langle\tfrac14\mathrm{tr}F^2\rangle\ne0$ — which `CLAUDE.md` §5 lists as *another name*
  of the same open object. Computing the obstruction presupposes the condensate.
- **(b) Wrong direction.** A mass **needs** a scale, so gap $\Rightarrow$ scale-breaking; but
  scale-breaking $\not\Rightarrow$ gap in general (a running/marginally-deformed theory can flow to a
  gapless IR fixed point with $\beta\neq0$ along the way). It is a **one-way screen**, not a gap
  detector — exactly the `[NG]` “trace anomaly is the wrong direction” obstruction already recorded.
- **(c) The exact flat-space free-massless anomaly is forced to zero**, so “LC $=$ scale-invariance”
  is **perturbatively anchored to the free CFT**, i.e. **mass-blind at fixed scale** — the two-tier
  blindness lemma (`TOOL-FACTORIZATION-ALGEBRA.md` §2), the rigidity no-go B3 in conformal dress.

So the rescuable reading collapses onto the **gluon condensate / trace anomaly**, an already-banked
name of the target, not a new engine.

---

## §2 · The obstruction, computed exactly in the solvable models `[T]` `[C]` `[V]`

### §2.1 Free scalar ($d=1,2,3,4$) — the obstruction is the mass `[T]` `[V]`

§1.2's table is the exact statement: LC-defect $=$ Kernel decay rate; it is $m$ (massive, LC) versus
$0$ (massless, not-LC). No free-field subtlety; the propagators are standard Bessel kernels
($K_0,K_1$ in $d=2,4$; elementary in $d=1,3$). **The obstruction and the mass are the same number.**

### §2.2 2D Yang–Mills — LC and the gap are **independent axes** `[T]/[C]` `[V]`

`[T]` *(standard 2D YM on a circle).* The transfer matrix acts on holonomy states labelled by
irreps $R$ of the gauge group, with $H|R\rangle=\tfrac{e^2}{2}C_2(R)|R\rangle$ (convention-dependent
normalization; the structure is standard). The vacuum is $R=\mathbf 1$ and the gap is
$$\Delta_{2D}=\tfrac{e^2}{2}\min_{R\ne\mathbf 1}C_2(R).$$
Witness values ($e^2=1$): $SU(2)$: $C_2(1/2)=\tfrac34$ (lightest) $\Rightarrow\Delta=0.375$;
$SU(3)$: $C_2(\mathbf 3)=\tfrac43\Rightarrow\Delta=0.667$. **The theory is gapped.**

`[C]` *(the observable FA is locally constant).* 2D YM has **no local propagating degrees of
freedom**: $F$ is a 2-form $=$ a top form on the surface, so all local curvature data is
topological, and the only gauge-invariant local observables are holonomy/Wilson data; the area
dependence is a **scalar weight** $e^{-\frac{\text{area}}{2e^2}C_2}$ on the TQFT sectors, not a
scale-dependent *local algebra*. Hence the observable FA depends only on the homotopy type of regions
— **locally constant**, defect $\equiv 0$ (witness). This is the standard “2D YM is a (area-deformed)
TQFT” identification; conditioned on it, LC holds.

**The obstruction computation, and the kill it supplies.** 2D YM is simultaneously **LC** (defect $0$,
topological) and **gapped** ($\Delta=3e^2/8$). Therefore:
- Under the **task's sign** (*gap $=$ failure of LC*): defect $=0$ $\Rightarrow$ predicted **gapless**,
  but 2D YM **is** gapped — a **direct counterexample**, `[NG]` for “failure of LC $=$ gap.”
- Under the **standard sign** (*gap $\iff$ LC*): LC ✓ and gap ✓ are consistent — but LC is satisfied
  for a **topological** reason (no local dof) that **never computes the number** $\Delta$; it is a
  Boolean “gapped/not,” blind to the value, and satisfied even by theories (pure topological) where
  there is **no** propagating gap at all (§2.3). So even correctly signed, LC does not track the
  **value** or furnish a **lower bound** on $\Delta$.

**The structural lesson:** *LC is a **local/algebraic** property* (“are there scale-dependent local
degrees of freedom adjacent to the vacuum?”), while *the mass gap is a **spectral** property of the
global transfer matrix.* In free massive fields both are controlled by $m$ and coincide; **2D YM
separates them.**

### §2.3 3D solvable models — where LC and gap split cleanly `[T]` `[V]`

- **Topologically massive gauge theory (Maxwell–Chern–Simons, Deser–Jackiw–Templeton) `[T]` `[V]`.**
  Propagator $\sim e^{-\mu r}/\sqrt r$ with $\mu$ the CS mass; witness rate $\to1.0$ for $\mu=1$
  ($r=6,10,14$: $1.167,1.100,1.071$). LC holds iff $\mu>0$, and the “obstruction” **is** $\mu$ — an
  **input parameter** inserted into the quadratic operator. Local constancy and the gap coincide here,
  and the gap is never **derived**, only **read off** the kernel.
- **Pure 3D Chern–Simons `[T]`.** Topological: the FA is locally constant (defect $0$), and the
  Hamiltonian vanishes on the physical Hilbert space — **no propagating gap**. Under the task's sign
  this “works” (LC $\to$ no gap) but **vacuously**: it is a TQFT, with no mass scale to detect. Under
  the standard sign it is a degenerate “gapped” with no dynamics.
- **3D free scalar `[T]` `[V]`.** Massive LC / massless not (same as §1.2 at $d=3$).

**Synthesis of §2.** Across free, 2D YM, and 3D topological-massive, the “obstruction to local
constancy” evaluates to the **analytic kernel** — the propagator decay rate or transfer eigenvalue
($m$, $\Delta$, $\mu$) — or is $0$ for **topological** reasons unrelated to any gap. It is never a
new lower-bound engine; it is the transfer matrix, re-labeled.

---

## §3 · Genuinely different, or a repackaging? — the tier map onto H1 / W / I / weak-$*$

The campaign's banked coordinate/tier structure (R-/W-wave): **weak-$*$ atom** $f_1^2>0$ (order-1
Perron jump of the measure, `ch09`), **H1** ceiling $\chi\le1/\lambda$, **W** atom weight
$w_1=f_1^2/M_0\ge w_0$, **I** isolation $\rho^2=\mu_2/\mu_1\le C^2$ (window theorem); the gap needs
$W\wedge H1$ or $I\wedge(\rho^2\le C^2)$.

| object | tier / strength | what it states | relation to LC |
|---|---|---|---|
| **weak-$*$ atom** $f_1^2$ | order-1 (Perron), bottom of Källén–Lehmann | discontinuity of the spectral measure | **orthogonal**: LC $\iff$ gap, and a gap does **not** force an atom (gapped-atomless measures, R6); LC says nothing about $f_1^2$ |
| **H1** $\chi\le1/\lambda$ | quantitative **ceiling** | susceptibility bound of the zero mode | **LC $\Rightarrow$ H1-type**: LC ($=$ clustering) $\Rightarrow$ gap $\Rightarrow$ finite $\chi$; LC is **strictly stronger**, but supplies **no** value of $\chi$ |
| **W** $w_1\ge w_0$ | atom-weight **input** to the gap | normalized first-excited KL mass | **not implied**: LC $\iff$ gap $\ne$ “$\exists$ atom”; W is a *sufficient ingredient* the campaign cannot yet prove, LC the *conclusion* restated |
| **I** $\rho^2\le C^2$ | order-2 isolation **input** | separation $m_{0^{++*}}/m_{0^{++}}$ | **not implied** for the same reason; I $\Rightarrow$ gap $\Rightarrow$ LC, but not conversely as a *tool* |
| **LC** | $=$ the **gap** (qualitative) | observable FA clusters exponentially / transfer operator has a spectral gap | — |

**Verdict of the map.**
- It does **not reduce to H1**: LC is strictly stronger (H1 is one quantitative face). 
- It does **not reduce to W or I**: those are the *sufficient inputs* the campaign is trying to
  establish; LC is the *conclusion* (“the gap”) restated, so LC is what W$\wedge$H1 / I$\wedge$window
  are **trying to prove** — proving LC is proving the gap.
- It **does** reduce to the **transfer matrix**: LC $\iff$ exponential clustering $\iff$ a positive
  spectral gap of the transfer operator (the Tauberian identity `gap $\iff$ decay`,
  `HISTORICAL-TOOL-SWEEP` #12, in homotopical dress). This is exactly the kill clause **“only renames
  the transfer matrix.”**

**The in-principle escape, and why it is foreclosed `[NG]`.** The one feature of the FA not present in
the two-point function is the **homotopy-coherent / derived ($E_d$, higher-operation) structure** —
in principle it could see an obstruction no correlator sees. But (i) every **deformation-invariant**
derived invariant is **mass-blind** (two-tier blindness lemma: the family $m\mapsto\mathcal F_m$ is
continuous through $0$, so any invariant frozen on a connected family returns the free value — the
rigidity no-go B3); (ii) the **only** invariant that jumps is the **IR / large-scale limit**, which
**is** “gapped,” i.e. the restatement, whose nonperturbative computation is the missing machine; (iii)
the derived candidate functionals (a scale defect `δ_ℓ`, factorization-homology/Hochschild traces) are
themselves either the clustering identity or a determinant (A3) — all three already-closed bins
(`TOOL-FACTORIZATION-ALGEBRA.md` §3). **Two horns, no third** — the regress is *factored*, not broken.

---

## §4 · Kill verdict, exact missing step, and flip condition

**Kill verdict (`[NG]`-for-the-class + `[O]`).** The paradigm, **correctly signed**
(gap $\iff$ LC), is the **gap restated qualitatively**, computing (in solvable cases) nothing but the
analytic propagator / transfer eigenvalue. It is **not** a sixth species:
- reduces to **H1**? **No** — strictly stronger (H1 is the ceiling; LC is the gap). 
- reduces to **W / I**? **No** — those are the sufficient *inputs*; LC is the *conclusion*.
- **renames the transfer matrix**? **Yes** — triggered. LC $\iff$ a positive transfer gap; the
  “obstruction” is the kernel decay rate / eigenvalue ($m$, $\Delta$, $\mu$).
The gap tier is exactly what the R-/Q-waves established **does not close Clay** (the atom/isolation,
not the gap, is irreducible), and the task's **failure-of-LC** sign additionally contradicts the free
anchor and 2D YM. Neither reading gives a nonperturbative engine for $\beta=2.5$.

**Exact missing step `[O]`** (identical shape to `TOOL-FACTORIZATION-ALGEBRA.md` §5, so as not to
double-count). The object that would flip the verdict is an invariant $\mathcal I$ of the observable
factorization algebra with:
1. $\mathcal I$ **jumps** at the gap (not deformation-invariant) **without** being handed the
   propagator $P$ or the IR limit;
2. $\mathcal I$ is **UV-computable** (from the micro-local / CG UV data), not requiring the
   nonperturbative IR / scaling limit;
3. $\mathcal I$ yields a **uniform-in-$(L,a)$ lower bound** $\ge\lambda>0$;
4. $\mathcal I$ **fails the free massless field non-vacuously** ($=0$ there, $>0$ on YM).

No such invariant is known, and (1)$\wedge$(2) are in structural tension (a derived invariant is built
to be deformation-invariant; forcing it to jump forces it onto the analytic kernel, which is the IR
data). **Flip condition:** exhibit one and the verdict flips to BROKEN.

---

## §5 · Witness and references

**Witness (this session, `[V]`).**
- `D:\Dev\math\research\v3\w7\x2_factorization.py` — exact Euclidean free propagators ($d=1,2,3,4$),
  decay-rate extraction (massive $\to m$, massless $\to 0$); 2D YM Casimir spectrum and gap
  ($SU(2)$: $3e^2/8$; $SU(3)$: $2e^2/3$) with LC-defect $=0$; 3D topological-massive rate $\to\mu$.
  Run `py -3.12 -I x2_factorization.py` from `w7/`.
- `D:\Dev\math\research\v3\w7\x2_factorization_results.json` — the tabulated output.

**Internal.**
- `D:\Dev\math\CLAUDE.md` §1 (acid test: fail on free massless, positive on YM), §5 (the wall; the
  gluon condensate $\langle\tfrac14\mathrm{tr}F^2\rangle$ as a name of the key object; the missing
  machine).
- `D:\Dev\math\research\unconditional\TOOL-FACTORIZATION-ALGEBRA.md` — CG scouting memo: the
  two-tier blindness lemma (§2), the `δ_ℓ` build-and-kill (§3), the exact stuck step (§5). X2 does
  **not** re-derive these; it supplies the **sign audit** and the **2D-YM decoupling**.
- `D:\Dev\math\research\v3\Q-WAVE-SYNTHESIS.md` §1/§4 (H1 ceiling, the atom as irreducible target),
  `R-WAVE-SYNTHESIS.md` §1/§8 (isolation $I$ is the parent; $W$ and $I$ independent; gapped-atomless),
  `W-WAVE-CHARTER.md` §0 ($W$: $w_1=f_1^2/M_0\ge w_0$; $I$: $\rho^2=\mu_2/\mu_1\le C^2$),
  `docs\paper-v2\ch09-weakstar.tex` (weak-$*$ atom, Theorem 4: no continuous forcing functional).
- `HISTORICAL-TOOL-SWEEP.md` #12 (Tauberian identity gap $\iff$ decay); `INVENT-NEW-INDEX.md` (rigidity
  no-go B3); `INVENT-NEW-IDENTITY.md`/A3 (Fredholm determinant).

**External (bibliographic; `[V]`-checked in `TOOL-FACTORIZATION-ALGEBRA.md`, not re-fetched here).**
- K. Costello & O. Gwilliam, *Factorization Algebras in Quantum Field Theory*, Vol 1 (New Math.
  Monographs **31**, CUP 2017), Vol 2 (**41**, CUP 2021) — LC of the massive free field; the free-scalar
  examples.
- K. Costello, *Renormalization and Effective Field Theory*, AMS Surveys **170** (2011) — perturbative,
  formal $\hbar$.
- O. Gwilliam & K. Rejzner, *Relating nets and factorization algebras of observables: free field
  theories*, arXiv:1711.06674, *CMP* **373** (2020) 107–174; *The observables of a pAQFT form a
  factorization algebra*, arXiv:2212.08175.
- M. Kawahira & T. Shigemura, *States and IR divergences in factorization algebras*, arXiv:2412.08183
  (2024; v3 2025) — IR-divergence handling in the **massless** theory; mass is an input.
- D. Ayala & J. Francis, *Factorization homology of topological manifolds*, J. Topology **8** (2015)
  1045–1084, arXiv:1206.5522 — LC $\iff$ $E_d$.
- J. Lurie, *On the classification of topological field theories*, arXiv:0905.0465 (2009) —
  cobordism hypothesis (deformation-invariant classification).
- 2D YM on a circle (Casimir spectrum / heat-kernel transfer matrix): the Migdal–Witten exact solution
  (standard; reported, not re-derived here).

---

*Status: **Sign audit `[T]`** — the CG theorem is gap $\iff$ **LC**, the task's “gap $=$ failure of LC”
is its negation, and the sole acid-test-matching reading collapses to the trace anomaly $=$ gluon
condensate `[C]`. **Solvable obstruction `[T]/[C]/[V]`** — the obstruction is the propagator /
transfer eigenvalue ($m$, $\Delta=3e^2/8$, $\mu$); 2D YM is locally constant (topological, defect $0$)
**and** gapped, so LC and the mass gap are independent axes, and the task's sign is a direct
counterexample. **Kill verdict `[NG]`/`[O]`** — not H1, not W, not I; it **renames the transfer
matrix**; not a sixth species. The escape (a non-deformation-invariant, UV-computable, free-massless-
failing derived invariant, §4) is the exact missing step and is structurally obstructed. **The wall is
relabeled, not crossed; no fabrication; Clay stays open.***
