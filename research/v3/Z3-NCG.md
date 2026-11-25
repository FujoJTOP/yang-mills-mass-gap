# Z3-NCG — the noncommutative-geometry / spectral-triple paradigm: the **spectral action is the classical YM action re-expressed** (finite and positive on the free massless field ⇒ acid-vacuous); the **Dirac spectral gap is either an input (finite triple), a UV/compact artifact (classical background), or exactly H1 (quantum transfer generator)**; **index-theoretic NCG is the rigidity no-go**. No 4D YM *quantum* spectral triple exists.

> **2026-11-24 · Wave Z3 — the noncommutative-geometry / spectral-triple paradigm.** The task: test whether
> Connes' spectral action / Dirac "spectral gap" gives a **rigorous mass-gap statement** for **4D pure
> Yang–Mills**; test on **finite spectral triples** and the **Moyal plane**; apply the free-massless acid
> test. Kill criteria: (1) no 4D YM spectral triple exists; (2) the spectral gap is **assumed**; (3) it
> **reduces to H1**.
> **Witness:** `research/v3/w7/z3_ncg.py` → `z3_ncg_results.json` (`py -3.12 -I`, numpy only, deterministic;
> **21/21 PASS**, exit 0).
> **Reads.** `docs/paper-v2/ch07-nogos.tex` (§rigidity, §wrong-half), `Q6-H1-STATUS.md` (precise H1),
> `X1-MODULAR.md` (the modular $=\beta\cdot$H1 reduction pattern), `docs/paper-v2/ch02-squeeze.tex` (H1).
> **Labels (charter §2).** **[T]** theorem/identity · **[C]** conditional · **[NG]** no-go for a named
> class · **[O]** exact missing step · **[V]** numerical witness. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED)

A spectral triple is $({\mathcal A},{\mathcal H},D)$: a $\ast$-algebra on a Hilbert space ${\mathcal H}$ and
a self-adjoint **Dirac operator** $D$. The paradigm offers three things, each of which could in principle
carry "the mass gap." **None does, and each dies for a named reason.**

1. **[NG] The spectral action $\operatorname{Tr} f(D/\Lambda)$ is the classical Yang–Mills action rewritten,
   not a quantization.** Its heat-kernel expansion reproduces $\int \tfrac14\operatorname{tr}F^2$ at leading
   gauge order — the same classical functional written down by Yang & Mills in 1954. It is *not* a path-integral
   measure, a Hamiltonian, or a vacuum state of the quantum theory. **Decisive acid datum [V]:** the spectral
   action is **finite and positive on the free massless field** (and on the free massive field alike, differing
   only by the smooth Gaussian factor $e^{-m^2}$). By the charter's acid rule — *a gap tool must fail on the
   free massless field* — the spectral action is **acid-vacuous**: it fires on the gapless free theory, so it
   cannot be a gap order parameter. **(Witness P5.)**

2. **[NG] The Dirac "spectral gap" $\min|\operatorname{spec}D\setminus\{0\}|$ is never the dynamical quantum
   mass.** It splits by construction:
   - **Finite spectral triple:** the internal $D_F$ is a finite matrix whose entries are the **input**
     masses/Yukawas; its gap *is* the smallest input mass, exactly, and scales linearly with the scaling of
     the input. A finite-triple "mass gap" is **among the assumptions**. **(Witness P2.)** Kill criterion (2)
     **triggered.**
   - **Classical background:** the Dirac operator of a fixed (flat or curved) geometry has a spectrum that is
     the *single-particle* spectrum of a classical field, not the quantum excitation spectrum. On the free
     massive field it has gap $=m$ — but only because $m$ is the **input mass parameter** [(P1e)]. On a box the
     gap is either an exact zero mode (periodic, massless) or $\sim 1/L\to0$ (thermodynamic); and for the
     almost-commutative $D=D_M+\gamma_5 D_F$ the internal gap $m_F$ is **washed out** by the continuum of
     $D_M$ momenta (apparent gap $\to0$ as the box grows, $(P3)$). It is a **UV/"wrong-half" scale**, not an
     IR mass.
   - **Quantum:** to make "$D$" an operator of the *quantum* theory is to take the **generator of the transfer
     operator**; its gap is the **transfer-matrix spectral gap $=H1$** (ch02 `eq:h1`). **(Sec. 4.)** Kill
     criterion (3) **triggered** — exactly the reduction X1 found for the modular operator ($=\beta\cdot H1$).

3. **[NG] Index-theoretic NCG (local index theorem, spectral flow, JLO cocycle) is a discrete invariant and
   dies by the rigidity no-go of ch07.** Any integer locally constant under the connected $\beta$-family
   equals its free value $=0$; the gap must be a $\beta$-non-analytic *jump*, not a discrete invariant.
   `ch07-nogos.tex` already records this for "noncommutative geometry" explicitly. **[NG]**

4. **[NG]/[V] The Moyal plane $\mathbb R^4_\theta$ gives nothing.** Its **free** Dirac spectrum is *identical*
   to the commutative $\mathbb R^4$ (the star-product phase is unit-modulus on free plane waves), hence
   **gapless**; and $\operatorname{Tr}|D|^{-p}$ diverges at both ends for every $p$, so the free massless
   Moyal plane is not even a genuine spectral triple with compact resolvent / finite summability. **(Witness
   P4.)**

**Answer to the question.** *No* — the NCG / spectral-triple paradigm does not give a rigorous mass-gap
statement for 4D pure YM. There is **no quantum 4D YM spectral triple** (the existing object is the *classical*
geometry whose spectral action is the classical YM action); the apparent "Dirac gap" is an **input**, a
**UV/compact artifact**, or exactly **H1**; and the index sector is the **rigidity** no-go. **All three kill
criteria are triggered. Clay 0%.**

---

## 1 · What the paradigm actually is, precisely

| Object | Definition | What it is |
|---|---|---|
| spectral triple | $({\mathcal A},{\mathcal H},D)$; $D=D^\ast$, $[D,a]$ bounded, $|D|^{-1}$ compact (summable) | a *differential-geometric* (classical) datum |
| spectral action | $\operatorname{Tr} f(D/\Lambda)$; heat-kernel $\sim\sum_k f_k\Lambda^{4-k}a_k(D^2)$ | the **classical** bosonic action at scale $\Lambda$ |
| Dirac "spectral gap" | $\min|\operatorname{spec}D\setminus\{0\}|$ (or gap of $D^2$, $=m^2$) | a **single-particle** spectrum of a fixed operator |
| index-theoretic NCG | local index / spectral flow / JLO | a **discrete, homotopy** invariant |
| almost-commutative | $D=D_M\otimes 1+\gamma_5\otimes D_F$; gauge fields $=$ inner fluctuations of $D$ | the geometry whose action is YM$(+$gravity$)$ |

The **crux** is the mismatch of categories: the spectral action and the Dirac operator are **classical**
(functional of a fixed background $D$), whereas "the mass gap" is a **quantum** statement about the spectrum of
the Hamiltonian / the transfer operator of the *quantized* theory. The paradigm never crosses that line. For
pure YM the internal space is trivial ($D_F=0$): the only content is $D_M$ and its inner fluctuations, whose
spectral action is $\int\tfrac14\operatorname{tr}F^2+\cdots$ — the classical action.

---

## 2 · Does a 4D YM spectral triple exist? (kill criterion 1) — [NG]

**A classical one, yes; a quantum one, no.** Precisely:

- **(exists, [T])** An almost-commutative spectral triple whose **inner fluctuations** of $D$ are the YM
  connection, and whose **spectral action** expands to $\int\tfrac14\operatorname{tr}F^2$ (plus Einstein–Hilbert
  and cosmological terms), exists. This is Chamseddine–Connes / Chamseddine–Connes–Marcolli. It is a
  reformulation of the **classical** action.
- **(does not exist, [NG]/[O])** A spectral triple $({\mathcal A}_{\rm YM},{\mathcal H}_{\rm YM},D_{\rm YM})$
  encoding the **quantized** 4D pure-YM theory — a Hilbert space of *physical states* with a Dirac operator
  whose spectral gap is the *dynamical* mass — does **not** exist. Constructing ${\mathcal H}_{\rm YM}$
  presupposes the existence of the 4D continuum theory (the corpus's **G1**, itself open: `ch03-existence.tex`,
  `PROOF-SKELETON.md`). And once ${\mathcal H}_{\rm YM}$ is the physical Hilbert space, the natural "$D$" is the
  transfer generator, whose gap is **H1** (Sec. 4).

So kill criterion (1) is **triggered in the sense that matters**: the paradigm supplies a *classical
differential-geometric reformulation*, not a *quantum spectral triple* for 4D YM. The only "existence" it has is
the existence of the classical action, which was known before NCG.

---

## 3 · The spectral action is the classical action — and is acid-vacuous — [NG]/[V]

The heat-kernel expansion of $\operatorname{Tr}f(D/\Lambda)$ gives, at leading gauge order,
$a_4(D^2)\supset c_f\int\operatorname{tr}F_{\mu\nu}F^{\mu\nu}$ for a universal constant $c_f$ (the
Chamseddine–Connes coefficient, convention-dependent) — the **classical** Yang–Mills term. It is a functional of
the *fixed classical* field strength $F$; it contains no path integral and no quantum fluctuation. (The exact
$c_f$ is not load-bearing here; what matters is that the term is *the classical action*, not a measure.)

**The acid test (the sharpest form of the no-go).** On a 4-torus of side $L$, with $D^2$ having eigenvalues
$(2\pi|\mathbf n|/L)^2+m^2$ (spinor degeneracy 4), the spectral action is

$$\operatorname{Tr}f(D^2/\Lambda^2)=4\sum_{\mathbf n\in\mathbb Z^4} f\!\big(((2\pi|\mathbf n|/L)^2+m^2)/\Lambda^2\big),\qquad f(x)=e^{-x}.$$

Witness **P5** evaluates it: $\operatorname{Tr}f(D^2)=1660.05$ at $m=0$ (**gapless** spectrum) and
$610.70$ at $m=1$ (**gapped**), with ratio exactly $e^{-1}$ and cutoff stable to $10^{-16}$. Consequences:

- The spectral action is **finite and positive on the free massless field** $\Rightarrow$ it **does not fail**
  on the gapless free theory $\Rightarrow$ it is **acid-vacuous** as a gap detector (it fires on a gapless
  spectrum). **[NG]**
- Its $m$-dependence is the **smooth Gaussian factor** $e^{-m^2}$ — no threshold, no jump, no order-parameter
  structure. It cannot separate "gapped" from "gapless"; it is a smooth classical functional. **[V]**

This is the cleanest statement of why the spectral action cannot carry the mass gap: **it is blind to
gappedness because it is the classical action, and the classical action is the same functional for the free
massless and the free massive fields (up to the explicit mass term).**

---

## 4 · The Dirac "spectral gap": assumed, or UV, or $=H1$ — [NG]

Three readings, three failure modes.

**(a) Finite spectral triple — the gap is an input [V], kill (2) triggered.** The internal $D_F$ is a finite
matrix. For $D_F=\operatorname{diag}(m_1,\dots,m_n)$ the gap $\min_{m_i\ne0}|m_i|$ **is** the smallest input
mass, exactly, and scales linearly in a uniform rescaling of the input (P2: masses $\{0,0.05,0.5,2\}\to$ gap
$0.05$; scale $\alpha\to$ gap $\alpha\cdot0.05$). An off-diagonal ("Yukawa") $2\times2$ block has smallest
$|\text{eigenvalue}|$ given exactly by the input formula (P2c). There is **no dynamics**: the finite-triple gap
is a **parameter**, i.e., it presupposes the answer.

**(b) Classical background — the gap is UV / tunable [V], "wrong-half" no-go.** On the free field the
single-particle Dirac operator has gap $=m$ — but $m$ is the **input mass** (P1e: gap $=1.3000$, the input).
For the almost-commutative $D=D_M+\gamma_5 D_F$ the total spectrum is $\{\lambda\pm m_F:\lambda\in
\operatorname{spec}D_M\}$; on a box the apparent IR gap is **washed out** by the $D_M$ momenta as $L$ grows
(P3: apparent gap $0.285\to0.055\to0.019\to\dots\to7.0\times10^{-6}$ for $L=8\to256$, while the internal
$m_F=0.5$ stays fixed). The "$D_F$ gap" is a **compact/Kaluza–Klein (UV) artifact**, not an IR mass. On a
non-compact space the spectrum is a union of branches $\lambda\pm m_F$ with $\lambda$ filling $\mathbb R$ — the
total $D$ is **gapless for every finite $m_F$**.

**(c) Quantum — the gap is exactly H1 [T], kill (3) triggered.** The only reading in which "$D$" is an
operator of the *quantum* theory is $D=$ generator of the **transfer operator** on the physical Hilbert space.
Its spectral gap is, by definition, the **transfer-matrix gap**, i.e. the corpus's **H1**
(`ch02-squeeze.tex` `eq:h1`: $\chi=\operatorname{Var}_\mu(\xi_{k\to0})=\int s^{-1}d\nu\le 1/\lambda$). This
is the exact same reduction X1 established for the modular operator: *"gap in the transfer generator"* $=$
$\beta^{-1}\times$*"gap in the modular operator"* $=$ H1. The Dirac spectrum of the transfer generator **is not
new information**; it *is* the spectral gap one is trying to bound.

```
   finite D_F gap        ==  min input mass                         ASSUMED (parameter)      [P2]
   classical Dirac gap   ==  single-particle / KK scale,           WRONG-HALF (UV/compact)  [P1e,P3]
                            tunable, washed out by D_M
   quantum "D" gap       ==  transfer-operator spectral gap        = H1                     [T]
   index-theoretic NCG   ==  discrete homotopy invariant            RIGIDITY no-go (ch07)    [NG]
```

---

## 5 · The Moyal plane $\mathbb R^4_\theta$ — [NG]/[V]

The Moyal plane, $[x^\mu,x^\nu]=i\theta^{\mu\nu}$, is the canonical noncommutative test space. Its spectral
data:

- **Free Dirac spectrum is $\theta$-independent.** Free plane waves obey $e^{ikx}\star e^{ipx}
  =e^{i(k+p)x}\,e^{\frac i2 k\theta p}$, a **unit-modulus** phase (P4a: $|{\rm phase}|=1$ to $10^{-12}$ for a
  random sample). So $\operatorname{spec}D=\{\pm|\mathbf k|:\mathbf k\in\mathbb R^4\}$ — identical to the
  commutative plane: **continuous, gapless** (zero mode at $\mathbf k=0$). The deformation does not move the
  free spectrum.
- **Not a genuine spectral triple.** $\operatorname{Tr}|D|^{-p}=C\!\int_0^\infty k^{3-p}\,dk$ diverges at
  both ends for every $p$ (P4c/P4d: $\operatorname{Tr}|D|^{-4}$ grows as $\varepsilon\to0$ *and* as
  $\Lambda\to\infty$). Hence the free massless Moyal plane has **no compact resolvent / finite summability** —
  it fails the spectral-triple axioms, and the "dimension" only appears with both IR and UV cutoffs.
- **NCYM has no rigor.** The $\theta$-deformed YM theory (Seiberg–Witten map from ordinary YM) suffers UV/IR
  mixing; no rigorous construction, no proven mass gap. The Moyal plane is a deformation of a theory whose gap
  is already the open problem — it neither proves nor transports one.

So the Moyal plane gives **zero positive content** on the gap, and its free case is the acid test's "must
fail" (gapless) — a vacuous pass.

---

## 6 · Index-theoretic NCG $=$ the rigidity no-go — [NG]

The genuinely *algebraic* NCG invariants — the **local index theorem**, **spectral flow**, the **JLO
(co)cycle representation of $K$-homology / cyclic cohomology**, and the **Witten index** — are **discrete** and
**locally constant** under continuous deformation. `ch07-nogos.tex` §`rigidity` proves: any locally constant
integer index along the connected $\beta$-family (which reaches the free Gaussian at $\beta=\infty$) equals its
free value, hence $0$ for YM; and states explicitly that *"noncommutative geometry, localization, BPS, and
Witten-spinor indices all return $0$ at the free point; none can carry the gap."* The gap must be a
**$\beta$-non-analytic jump**, not any discrete invariant. This sector is therefore **closed**, not merely
unpromising. **[NG]**

---

## 7 · Acid-test table and kill ledger

| Acid test | Spectral action | Dirac gap (finite $D_F$) | Dirac gap (classical bkgd) | Quantum "D" gap | Index NCG |
|---|---|---|---|---|---|
| free massless → no gap / fails | **fires positively** (1660.05) ⇒ **VACUOUS** | n/a (input) | gap $=0$ ✓ but reads no input | $=0$ (H1 with $\lambda=0$) | $=$ free value $0$ |
| free massive → gap | also positive (smooth $e^{-m^2}$) | gap $=$ input $m$ | gap $=m$ (**input**) | $=m\,/\,$H1 | $0$ |
| Moyal plane | finite, $\theta$-indep. | n/a | gapless | n/a | $0$ |
| **a gap tool for YM?** | **no** (acid-vacuous) | **no** (assumed) | **no** (UV/tunable) | **no** ($=H1$) | **no** (rigidity) |
| **reduces to** | classical YM action | an input | wrong-half (UV) | **H1** | **rigidity (ch07)** |

**Kill ledger.** The kill criteria, and how each fires:

```
 (1) no 4D YM spectral triple       -> only a CLASSICAL triple exists; the QUANTUM one
                                        ({H_YM, D_YM}) presupposes G1 and is absent.   TRIGGERED
 (2) spectral gap is ASSUMED        -> finite-triple gap = input masses (P2); free
                                        massive Dirac gap = input m (P1e).            TRIGGERED
 (3) reduces to H1                  -> quantum "D" is the transfer generator; its gap
                                        = transfer spectral gap = H1 (Sec. 4c).        TRIGGERED
     (+ index NCG -> rigidity no-go; + spectral action -> acid-vacuous; + Moyal -> gapless/non-ST.)
```

---

## 8 · What would upgrade it — the exact missing step [O]

The one object that would turn NCG into a positive tool is a **canonically associated quantum spectral triple**
$({\mathcal A}_{\rm YM},{\mathcal H}_{\rm YM},D_{\rm YM})$ with:

- ${\mathcal H}_{\rm YM}$ the physical Hilbert space of 4D pure YM (so its construction $\equiv$ **G1**), and
- $D_{\rm YM}$ an operator whose **spectral gap is the dynamical mass gap** — i.e., $D_{\rm YM}$ is *not* the
  transfer generator (else the gap is H1), and *not* a fixed classical background (else the gap is a UV
  artifact), but a canonically attached operator whose spectrum encodes the *quantum* excitation gap.

No such $D_{\rm YM}$ exists in the literature or the corpus. Constructing it and proving its gap $>0$ is either
**equivalent to proving the mass gap** (circular) or, if $D_{\rm YM}$ is the transfer generator, an exact
**reduction to H1**. This is the **exact missing step** — the same "wrong-half / UV⊥gap" and "$=H1$" walls the
corpus hits everywhere (`ch07-nogos.tex` §wrong-half; `Q6-H1-STATUS.md`). NCG does not supply the missing
machine; it supplies a classical reformulation of the action the machine must act on.

---

## 9 · References (standard theory; not re-verified here)

- A. Connes, *Noncommutative Geometry* (Academic Press, 1994) — spectral triples, reconstruction, spectral
  distance $d(p,q)=\sup\{|a(p)-a(q)|:\|[D,a]\|\le1\}$.
- A. H. Chamseddine, A. Connes, "The spectral action principle," *Commun. Math. Phys.* **186** (1997) 731,
  arXiv:hep-th/9606001 — $\operatorname{Tr}f(D/\Lambda)$ and its heat-kernel expansion.
- A. H. Chamseddine, A. Connes, M. Marcolli, "Gravity and the standard model with neutrino mixing,"
  *Adv. Theor. Math. Phys.* **11** (2007) 991, arXiv:hep-th/0610241 — almost-commutative SM.
- W. D. van Suijlekom, *Noncommutative Geometry and Particle Physics* (Springer, 2015) — spectral action,
  finite triples, inner fluctuations; the gauge term $\propto\int\tfrac14\operatorname{tr}F^2$.
- J. M. Gracia-Bondía, J. C. Várilly, H. Figueroa, *Elements of Noncommutative Geometry* (Birkhäuser, 2001).
- A. Connes, M. Marcolli, *Noncommutative Geometry, Quantum Fields and Motives* (AMS, 2008) — NCG and QFT,
  renormalization, the spectral action.
- V. Gayral, J. M. Gracia-Bondía, J. C. Várilly, "Moyal planes are spectral triples," *Commun. Math. Phys.*
  **246** (2004) 569, arXiv:hep-th/0307241 — the *regulated* Moyal spectral triple (dimension-spectrum, with
  cutoffs; the free massless case is non-summable without them).
- N. Seiberg, E. Witten, "String theory and noncommutative geometry," *JHEP* **09** (1999) 032,
  arXiv:hep-th/9908142 — the Seiberg–Witten map and UV/IR mixing in NCYM.
- Corpus: `docs/paper-v2/ch07-nogos.tex` (§rigidity, §wrong-half), `ch02-squeeze.tex` (`eq:h1`),
  `research/v3/Q6-H1-STATUS.md`, `research/v3/X1-MODULAR.md`. `PROOF-SKELETON.md` (G1 existence).

*One line.* **NCG/spectral triples do not give a 4D pure-YM mass-gap statement: the spectral action is the
classical YM action (finite and positive on the free massless field ⇒ acid-vacuous); the Dirac spectral gap is
an input (finite triple), a UV/compact artifact (classical background, tunable/washed-out), or exactly H1
(transfer generator); index-theoretic NCG is the rigidity no-go; the Moyal plane is free-gapless and not a
genuine spectral triple. No quantum 4D YM spectral triple exists (constructing one presupposes G1). All three
kill criteria triggered: assumed / UV-wrong-half / =H1. No QED, Clay 0%.**
