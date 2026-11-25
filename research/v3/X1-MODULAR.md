# X1-MODULAR — the mass gap as a "modular / Connes-spectrum gap": the **Connes-spectrum reading is a no-go** (local vacuum algebras are type III₁, Connes spectrum = ℝ); the **global thermomodular reading is an exact repackaging of H1** (`= β·` transfer gap); the **local entanglement reading passes the gappedness acid test but is a non-linear, state-derived *read* (U5-class), not the mass value and not a forcing tool**

> **2026-11-24 · Wave X1 — first step of the modular / algebraic-QFT paradigm.** The task: define the
> mass gap as a gap in the **modular operator** (Tomita–Takesaki) / **Connes spectrum**, test massive
> vs massless free fields in 1D/2D solvable models, and decide whether the modular gap is **equivalent**
> to the mass gap or only a **repackaging**. Acid tests: free massless → continuous/0; free massive →
> positive gap; 2D solvable. Kill criterion: if it reduces to H1/W/I, say so.
> **Witness:** `research/v3/w7/x1_modular.py` (`py -3.12 -I`, numpy 2.4.6 + scipy 1.18.1, deterministic;
> **21/21 PASS**, exit 0).
> **Reads.** `docs/paper-v2/ch09-weakstar.tex`, `T-WAVE-SYNTHESIS.md`, `U-WAVE-SYNTHESIS.md`,
> `U5-ENTANGLEMENT-BINDING.md`, `U6-BINDING-LOGIC.md`, `Q6-H1-STATUS.md`, `R2-ISOLATION.md`.
> **Labels (charter §2).** **[T]** theorem/identity · **[C]** conditional · **[NG]** no-go for a named
> class · **[O]** exact missing step · **[V]** numerical witness. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED)

There are **two inequivalent** operators that the phrase "the modular operator" could name, and **three
readings** of "mass gap = a modular/Connes-spectrum gap." All three fail to give a new tool.

1. **[NG] Connes-spectrum reading — impossible in principle.** The genuine algebraic object is the
   modular operator of the **vacuum** on a **local algebra** $\mathcal A(O)$. By the Haag–Kastler /
   Araki theory, $\mathcal A(O)$ is a **type III₁ factor** (Araki 1964 for the free field; the standard
   result for interacting QFT). A type III₁ factor has **Connes spectrum all of $\mathbb R$** — its
   modular spectrum is *by definition* gapless and mass-independent. So "the mass gap is a gap in the
   Connes spectrum of the vacuum local algebra" is **not false, it is unsatisfiable**: there is no
   Connes-spectrum gap for a type III₁ factor, for *any* mass. The Bisognano–Wichmann theorem makes
   this concrete: on a wedge the modular operator is $e^{-2\pi K_{\text{boost}}}$, the boost generator,
   with continuous spectrum $\mathbb R$ — for massive and massless alike. **[NG]**

2. **[T] Global thermomodular reading — exact, but a pure repackaging of H1.** For the finite-volume
   equilibrium (Gibbs) state $\rho_\beta=e^{-\beta H}/Z$ the modular Hamiltonian is
   $K_{\text{mod}}=\beta(H-E_0)$, so the **modular gap** $\inf(\operatorname{spec}K_{\text{mod}}\setminus\{0\})
   =\beta\cdot(\text{mass gap})$ — an **exact information-preserving bijection** whose only extra datum
   is the free scale $\beta$. The modular operator *is* the Boltzmann weight; "gap in the modular
   operator" is the **transfer-matrix / spectral gap** ($=$ H1) rescaled by $\beta$. To know it you must
   know the spectrum of $H$, i.e. the original problem. **Repackaging of H1.** **[T]**

3. **[T]+[V] Local-entanglement reading — a genuine gappedness order parameter, but a *read*.** The
   lowest eigenvalue $\varepsilon_{\min}$ of the subregion modular Hamiltonian $K_A=-\log\rho_A$
   ("entanglement gap", Li–Haldane) is **positive for every $m>0$** (a saturating plateau, $N$- and
   $\ell$-independent) and **$\to0$ only at exactly $m=0$**. So it **passes the gapped/critical acid
   dichotomy**. But $\varepsilon_{\min}$ is (i) a **non-linear, cutoff- and region-shape-dependent**
   functional of $m$ — *not* a bijection to the value ($\varepsilon_{\min}/m$ runs from $35$ at
   $m{=}0.05$ to $3.1$ at $m{=}2$); (ii) a **read of the reduced density matrix** — it *presupposes the
   ground state*, so computing it for YM needs the YM vacuum's entanglement spectrum; and (iii) in the
   continuum the local algebra is type III₁ and the modular spectrum is dense (W4), so the entanglement
   gap is a **lattice-scale feature**. It is the **U5 / ch09 state-derived read-out class, gated by
   H1** — a **read, not a force**. **[T]+[V]**

**Answer to the question.** The modular gap is **not an independent equivalent** of the mass gap. The
**Connes-spectrum** form is a **no-go** (type III₁: no such gap exists). The **global** form is an
**exact repackaging of H1** ($\beta\times$ transfer gap). The **local** form is a **faithful gappedness
*detector*** but a **non-linear, state-derived read** (U5-class), not the value and not a forcing tool.
**Kill criterion TRIGGERED: it reduces to H1 (global, exactly) and to the U5/W state read (local).**

---

## 1 · The two modular operators, precisely

For a von Neumann algebra $\mathcal M$ with a cyclic separating vector $\Omega$, Tomita–Takesaki gives
the modular operator $\Delta$ and modular conjugation $J$; $\Delta^{it}$ generates $\sigma_t$, and
$K=-\log\Delta$ is the modular Hamiltonian. The phrase "modular operator" is used two ways:

| | GLOBAL (thermomodular) | LOCAL (subregion / algebraic) |
|---|---|---|
| state | equilibrium $\rho_\beta=e^{-\beta H}/Z$ on all of space | **vacuum** $\Omega$ restricted to $\mathcal A(O)$ (or $\rho_A=\operatorname{tr}_{\bar A}|\Omega\rangle\langle\Omega|$) |
| operator | $K_{\text{mod}}=\beta(H-E_0)$ | $K_A=-\log\rho_A$ (entanglement Hamiltonian) |
| gap | $\beta\times(\text{transfer gap})$ | "entanglement gap" $\varepsilon_{\min}$ |
| massive vs massless | $\beta m$ vs $0$ | plateau $>0$ vs $\to0$ |
| status | **$=H1$ repackaging** [T] | **gappedness read** [T]+[V] |

The algebraic-QFT paradigm (Bisognano–Wichmann, Connes, Haag–Kastler) is about the **local** one; the
"mass gap" is a property of the **spectrum of translations**, not of the modular operator. The two live
on different objects, and the rest of this note makes that exact.

---

## 2 · [NG] The Connes-spectrum form is unsatisfiable (type III₁)

A finite-volume lattice theory has a **type I** algebra (finite-dimensional), so its modular spectrum
is discrete and *does* have a gap — this is the **finite-volume artifact** the ch09 corollary flags as
vacuous for Clay (a compact transfer matrix is trivially gapped at finite $V$). In the **continuum**:

- Local algebras $\mathcal A(O)$ for $O$ a bounded region (or a wedge) are **type III₁** — proven by
  Araki (1964) for the free field and standard for interacting QFT satisfying the Haag–Kastler axioms.
- A type III₁ factor has **Connes spectrum $\Gamma=\mathbb R$** (equivalently $S(\mathcal M)=[0,\infty)$,
  the "flow of weights" with full spectrum). The modular spectrum is the **whole real line** — there is
  **no gap**, by definition, for any mass.
- Bisognano–Wichmann: for a wedge, $K_A=2\pi K_{\text{boost}}$, the boost has spectrum $\mathbb R$.

So "$m_1>0$ $\iff$ a gap in the Connes spectrum of the vacuum local algebra" is impossible: the RHS is
*always* gapless. This is the sharpest result of the wave, and it is a clean **no-go** for the named
class ("Connes-spectrum gap of a local QFT vacuum algebra"). **[NG]**

---

## 3 · [T] The global thermomodular gap is $=\beta\times$ the transfer gap

The one setting where "modular gap" is a *number* is the finite-volume Gibbs state. Then

```
   K_mod = beta (H - E0),     gap(K_mod) = beta * gap(H) = beta * (mass gap).     [T]
```

In the free theory the mass gap is exactly $m$ (the lightest one-particle mode $k=0$ has $\omega_0=m$;
massless has $\omega_0=0$ exactly, the constant zero mode). Witness W1/W3/W5 confirms, in 1D and 2D:
gap $=m$ for the massive field, $=0$ for the massless, **for every $N$**; and the ratio
`gap(K_mod)/gap(H)` is identically $\beta$. The map is an **exact bijection** — faithful, and therefore
**empty**: it adds only the free scale $\beta$ and requires the spectrum of $H$, which is the original
problem. This is the transfer-matrix / spectral-gap content, i.e. **H1**. Witness: W1, W3, W5. **[T]**

---

## 4 · [T]+[V] The local entanglement gap: a genuine gappedness *read*

The subregion modular Hamiltonian $K_A=-\log\rho_A$ of the vacuum has lowest eigenvalue
$\varepsilon_{\min}$, the **entanglement gap**. Exact Gaussian computation (W2), $N=512$–$768$:

| $m$ | 0 | 0.05 | 0.1 | 0.3 | 0.5 | 1.0 | 2.0 |
|---|---|---|---|---|---|---|---|
| $\varepsilon_{\min}$ ($\ell{=}64$) | $\to0$ (log-slow) | 1.77 | 2.14 | 2.95 | 3.51 | 4.61 | 6.28 |

- **Massive:** $\varepsilon_{\min}$ **saturates** to a positive plateau, independent of region size
  $\ell$ (for $\ell\ge32$) and of chain length $N\in\{512,768\}$ — a thermodynamic value. **[V]**
- **Massless:** $\varepsilon_{\min}\to0$ (logarithmically slowly, $\propto c-b\log\ell$), the hallmark of
  a critical system.

So the local modular operator **passes the gapped/critical acid dichotomy**: entanglement gap open
$\iff$ gapped ($m>0$); closed $\iff$ critical ($m=0$). This is a real, computable order parameter
(physically the Li–Haldane entanglement spectrum). **[T]+[V]**

**But it is not the mass gap value, and not a forcing tool (W2b):**
1. **Non-linear.** $\varepsilon_{\min}/m$ falls monotonically from $35$ ($m{=}0.05$) to $3.1$
   ($m{=}2$). There is no bijection $\varepsilon_{\min}\leftrightarrow m$; it is a monotone *scale/
   gappedness* readout, cutoff- and region-shape dependent.
2. **State-derived read.** $\varepsilon_{\min}$ is obtained by diagonalising $\rho_A$ — it
   **presupposes the ground state**. For YM this needs the YM vacuum's entanglement spectrum, which is
   exactly what one does not have. It is the **U5 / ch09 read-out class**, and it is **gated by H1**.
3. **Continuum.** The local algebra is type III₁ (Sec. 2); the genuine modular spectrum is dense (W4),
   so the entanglement gap is a lattice-scale feature, not a continuum invariant without a cutoff.

---

## 5 · [T]+[V] The continuum modular spectrum is dense (type III₁) for every mass

A **fixed** one-particle energy window holds an **extensive** number of modes: in the massless field,
$\#\{\omega_k\in[0,1]\}/N\to1/3$; in the massive field, $\#\{\omega_k\in[m,m+0.6]\}/N\to0.43$ — both
**stable** as $N{:}\,64\to1024$ (W4). Hence the one-particle spectrum is **dense**, $\operatorname{spec}
(\log\Delta)$ is a **continuum** for **every** mass, and no operator gap survives the thermodynamic
limit. This is the lattice shadow of type III₁ / Connes spectrum $\mathbb R$ (Sec. 2). **[T]+[V]**

---

## 6 · Acid-test table (both operators)

| Acid test | GLOBAL thermomodular | LOCAL (subregion) |
|---|---|---|
| free massless → continuous/0 | yes: gap $=\beta\cdot0=0$ | yes: $\varepsilon_{\min}\to0$ (critical) |
| free massive → positive gap | yes: gap $=\beta m>0$ | yes: $\varepsilon_{\min}$ plateau $>0$ |
| 2D solvable (W5) | yes | same structure |
| **equal to the mass gap?** | **no**: $=\beta\times$ transfer gap | **no**: non-linear, $\ne m$ |
| **a forcing tool?** | **no**: it *is* $H$, a repackaging | **no**: a *read* of the state |
| reduces to | **H1** (exactly, $\times\beta$) | **U5** (entanglement read), gated by H1 |

---

## 7 · Kill ledger — the reduction map

The kill criterion: *if it reduces to H1/W/I, say so.* It does, in every reading.

```
   global modular gap  ==  beta * (transfer / H1 spectral gap)      EXACT bijection -> H1
   local modular gap   ==  entanglement gap of rho_A                U5 / ch09 read-out, gated by H1
   first-excited state psi_1, residue  ==  the U-wave w_1 master identity   (W)
   Connes-spectrum gap of the vacuum algebra  ==  empty (type III1, spectrum R)   [NG]
```

- **To H1**: the global thermomodular gap is LITERALLY $\beta\times$H1's spectral gap (W1/W3/W5).
- **To W**: the relevant modular/Liouvillian first-excited state $\psi_1$ and its overlap with the
  source is the U-wave master identity $w_1=|\langle\psi_1|o-1\rangle|^2/\mathrm{Var}(o)$; the atom
  fraction is not produced by the modular spectrum.
- **To U5 (entanglement)**: the local modular gap is the entanglement spectrum — operator-blind and
  anti-correlated with $w_1$ in the atom regime, exactly as U5 established.
- **No-go**: the Connes-spectrum form is unsatisfiable for a type III₁ factor [NG].

**What is *not* proved / would be a genuinely new tool [O].** The one thing that would upgrade the
local reading from a "read" to a "force" is a **state-independent lower bound** on the YM vacuum's
entanglement gap, derivable from modular/Connes structure alone — $\varepsilon_{\min}(\text{YM})\ge c>0$
*without* knowing the state. No such bound exists in the corpus, and obtaining it would require
controlling the YM ground state's entanglement, which is the binding/residue problem again (it collapses
to H1). This is the **exact missing step**, not a solved one.

---

## 8 · What the wave proves and does not prove

**Proved.**
- **[NG]** The Connes-spectrum form of "mass gap" is unsatisfiable for the vacuum local algebra: local
  QFT algebras are type III₁, Connes spectrum $\mathbb R$; the BW wedge modular operator is the boost,
  gapless. (Standard theory; the no-go is its direct consequence.)
- **[T]** The global thermomodular gap is exactly $\beta\times$(transfer gap) — a faithful bijection and
  hence an exact repackaging of H1 (witness W1/W3/W5).
- **[T]+[V]** The local entanglement gap is a genuine gapped/critical order parameter (massive plateau,
  massless $\to0$, $N$- and $\ell$-independent) but a non-linear, state-derived *read*, not the mass
  value (witness W2/W2b).

**Not proved / open.** No new forcing functional. No state-independent lower bound on the YM
entanglement gap (the [O] above). No reduction of Clay: as everywhere in the corpus, the positive
content (gap or binding) remains gated by H1.

---

## 9 · References (standard theory; not re-verified here)

- J. J. Bisognano, E. H. Wichmann, *J. Math. Phys.* **16** (1975) 985; **17** (1976) 303 — modular
  operator of the wedge $=$ boost; $K_A=2\pi K_{\text{boost}}$.
- M. Takesaki, *Tomita's Theory of Modular Hilbert Algebras* (1970) — modular operator / group.
- A. Connes, *Ann. Sci. ÉNS (4)* **6** (1973) 133 — classification of type III factors; Connes
  spectrum $\Gamma$, invariant $S(\mathcal M)$; type III₁ has $\Gamma=\mathbb R$, $S=[0,\infty)$.
- H. Araki, *Prog. Theor. Phys.* **32** (1964) 956 — free-field local algebras are type III.
- H. Li, F. D. M. Haldane, *Phys. Rev. Lett.* **101**, 010504 (2008), arXiv:0805.0332 —
  entanglement spectrum; gap as a gapped/critical diagnostic.
- R. Haag, *Local Quantum Physics* (2nd ed., 1996) — Haag–Kastler net, type III₁ of $\mathcal A(O)$,
  Reeh–Schlieder.

*One line.* **The mass gap is not a modular/Connes-spectrum gap: the Connes form is a no-go (type III₁,
spectrum $\mathbb R$); the global thermomodular form is exactly $\beta\times$H1; the local entanglement
form is a genuine gappedness *read* but non-linear, state-derived, and gated by H1 — a read, not a force.
Kill criterion triggered: reduces to H1 / U5 / W. No QED, Clay 0%.**
