# U5-ENTANGLEMENT-BINDING — entanglement does **not** give a rigorous binding bound: the one exact link (single-mode: `w₁=Z`, `S=h(Z)`) is **symmetric, hence anti-correlated with the atom**; the free massless chain is a **counterexample to any universal floor** `w₁≥f(S)` (`S` large, `w₁=0`); MI is **never negative** (SSA); and `w₁` is **not a function of `ρ_A`** (basis-independence) — so the missing object is an operator-specific **coherence**, which is exactly the residue `= isolation = H*`

> **2026-11-24 · Attack surface U5 — the nonperturbative binding theorem attacked through entanglement.**
> The T wave reduced the binding criterion to a uniform floor `w₁=f₁²/M₀≥w₀>0` on the atom fraction of
> the minimal four-wall tube at `β=2.5` (`T2-TUBE-OVERLAP`, `T-WAVE-SYNTHESIS`). **This note asks whether an
> entanglement / strong-subadditivity / relative-entropy functional can supply that floor** — the tube is a
> correlated region, and binding might look like a mutual-information deficit or an entanglement deficit.
> **Reads:** `T-WAVE-SYNTHESIS.md`, `T2-TUBE-OVERLAP.md`, `T6-BINDING-ADVERSARIAL.md`, `SM5-TCI-HEAD.md`,
> `H6-HEAD-SPECTRAL.md`; `docs/paper-v2/ch06-tools.tex` (§acid, §toolwave), `ch07-nogos.tex` (§schur-master,
> §impossible, wrong-direction); `unconditional/IMPOSSIBLE-ENTANGLEMENT.md`, `unconditional/H1-T4-entanglement.md`;
> `R5-EXACT-ANCHORS.md` (§1.1 2D `trF²=H`), `GAP-G2-edge.md` (residue–level-gap L2).
> **Witness:** `research/v3/w7/u5_entanglement_binding.py` (`py -3.12 -I`, numpy only, deterministic, exit 0;
> output quoted §6). No third-party downloads.
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/model · **[NG]** no-go for a named
> class · **[O]** exact missing step · **[V]** numerical witness here · **[std]** standard tool used as textbook.
> **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

**Entanglement gives one exact relation and no floor, and the exact relation is turned against the
program by its own symmetry.** Six findings:

1. **The single-mode exact identity `[T]`.** In a genuine two-level (single-excitation / single-mode)
   subsystem the reduced state is `ρ_A=diag(Z,1−Z)`, the atom fraction is `w₁=Z`, and the subsystem
   entropy is the binary entropy `S=h(Z)`. Hence the *exact* lower bound `w₁ ≥ h⁻¹_<(S)` (small root).
   **But `h(Z)=h(1−Z)` is symmetric**, so `S=0` is compatible with *both* `w₁=0` (free field) *and*
   `w₁=1` (perfect atom); the bound is **vacuous at both ends** and strongest only at `w₁=½`. This is the
   **anti-correlation lemma**: *entanglement is maximal exactly where the atom is smallest.*

2. **No universal floor `w₁≥f(S)` `[NG]`.** The free massless chain is the counterexample: its bipartite
   entropy grows `∝ log l` (large, §6·B), yet the free `0⁺⁺` composite has `w₁=0` **exactly** (two-particle
   cut). So *large entropy with zero atom* is realized; no decreasing function of an entanglement monotone
   can lower-bound `w₁`. (The naive move — "plug the area-law `S` into `h⁻¹_<` ⟹ `w₁≥½`" — is exactly what
   the free chain falsifies.)

3. **Strong subadditivity forbids "negative mutual information" `[T]`.** `I(A:B)=S(A)+S(B)−S(AB)≥0`,
   `S(ABC)+S(BCD)≥S(BC)`, and the Araki–Lieb triangle `|S(A)−S(B)|≤S(AB)` hold identically (§6·C, 0
   violations in 4000 random pure 4-qubit states). A "mutual-information deficit" must therefore be a
   **relative entropy** `S(ρ_A‖σ)` against a *reference* tube state `σ` — and that reference *is* the atom,
   so the bound is **circular**.

4. **`w₁` is not a function of `ρ_A`; entanglement is operator-blind `[T]`.** Entanglement monotones are
   functions of `ρ_A` alone. But `w₁` is the overlap of the source with the **operator-chosen** tube ray:
   for the *same* `ρ_A` (same `S`) two operators give two different `w₁` (§6·D: `0.200` vs `0.824`). This is
   the entanglement-language form of the **Schur-convex master obstruction** (`ch07 §schur-master`): entropy
   is a basis-independent *aggregate*; the atom is a basis-dependent *coherence* (a vacuum↔tube off-diagonal,
   `f₁²=|⟨G|O|Ω⟩|²`). **An aggregate of `ρ_A` cannot lower-bound a coherence in an operator basis.**

5. **The MI/confinement diagnostic lives in the wrong (1-form/center) channel `[NG]`.** Mutual information
   across a confining string measures `σ>0`, but `σ>0` does not imply the `0⁺⁺` atom: **2D YM has
   `σ=g²C₂/2>0` (confinement) and `f₁²(trF²)=0` exactly** (`trF²=H` conserved, `R5 §1.1`). The witness
   §6·E reproduces the structure: a highly correlated state (`MI>0`) whose *conserved* density has `w₁=0`.

6. **What does survive, and the exact missing estimate `[O]`.** The one non-vacuous entanglement bound
   requires knowing `w₁≥½` already (the majority branch: `w₁=1−h⁻¹_<(S)`), i.e. it presupposes a floor
   *stronger* than the target — circular. The honest missing datum is an **operator-specific coherence /
   residue lower bound**: `w₁=|⟨G|O|Ω⟩|²/‖QO Ω‖² ≥ w₀`, the off-diagonal weight of the isolated pole in
   the source's spectral measure. This is **exactly** the polymer four-wall residue (`T2`) `=` isolation
   `I₂` (`L2`) `=` `H*`. The entanglement route **re-expresses the wall and adds a sharp no-go (basis
   independence); it does not lower it.**

**One line.** *Entanglement gives one exact identity `w₁=Z`, `S=h(Z)`, but it is symmetric (`h(Z)=h(1−Z)`),
so it is vacuous where the atom lives (anti-correlated); the free massless chain (large `S`, `w₁=0`) kills any
universal floor `w₁≥f(S)`; MI is never negative (SSA) and a deficit must be a circular relative entropy; and
`w₁` is not even a function of `ρ_A` (basis-independence) — the missing object is a coherence, which is the
residue `=` isolation `=` `H*`. No QED. Clay 0%.*

---

## 1 · Definitions: entropy / MI / relative entropy of the tube region

Lattice 4D pure SU(N) Wilson YM on `T³_L×_aℤ`, vacuum `Ω`, gauge-invariant measure `μ` on `𝒜/𝒢`. The
flux tube is a correlated region; take a spatial bipartition `A|A^c` (`A` = a slab or the transverse
tube neighborhood) and the reduced state `ρ_A=Tr_{A^c}|Ω⟩⟨Ω|`.

| object | definition | role |
|---|---|---|
| von Neumann entropy | `S(A) = −Tr ρ_A log ρ_A = −Σ λ_i log λ_i` | entanglement of the region |
| mutual information | `I(A:B) = S(A)+S(B)−S(A∪B) ≥ 0` | correlation across the cut |
| Rényi-2 / purity | `S₂(A) = −log Tr ρ_A²` | `e^{−S₂}=Tr ρ_A²` |
| relative entropy | `S(ρ‖σ) = Tr ρ(log ρ − log σ) ≥ 0` | "deficit" vs a reference |
| atom fraction (`T2`) | `w₁ = f₁²/M₀`, `f₁²=|⟨G|O|Ω⟩|²`, `M₀=‖QO Ω‖²`, `O=¼trF²` | the target |

The Källén–Lehmann measure is `dν=Σ_n w_n δ_{E_n} + (cut)` with `C(t)=∫e^{−√s t}dν`, `f₁²=dν({E₁})`,
`M₀=∫dν`; `w₁=dν({E₁})/M₀∈[0,1]` (`T2 §1.3`). The tube state is `|G⟩`, the lightest isolated `0⁺⁺`; the
source direction is `φ̂=QO Ω/‖·‖`, so `w₁=|⟨G|φ̂⟩|²` is the squared overlap of the source with the tube ray.

**The three candidate routes, stated exactly.** (i) *Entropy/MI*: bound `w₁` by a function of `S(A)` or
`I(A:B)` using SSA/subadditivity. (ii) *Area law*: use the leading area term (or its deficit). (iii)
*Relative entropy*: bound `w₁` (or leakage `1−w₁`) by `S(ρ_A^{tube}‖ρ_A^{ref})` via Pinsker / fidelity. §2–§4
show each either returns nothing, returns the wrong sign, or is circular.

---

## 2 · The one exact link — and its symmetry (the anti-correlation lemma)

The only place entanglement *equals* a function of the atom is the **single-mode / two-level** reduction.
Take a subsystem that carries a single excitation (one flux-tube mode, or one impurity coupled to a bath).
In the state `|Ψ⟩ = √Z |1⟩_A|0⟩_{A^c} + √(1−Z) |0⟩_A|1⟩_{A^c}` (one quantum shared between `A` and the
bath), the reduced state is genuinely two-level, `ρ_A=diag(Z,1−Z)`.

> **Lemma U5-1 (single-mode exact identity) `[T]`.**
> For the two-level reduction with `ρ_A=diag(w₁,1−w₁)`,
> ```
>   w₁ = Z ,      S(A) = h(w₁) := −w₁ log w₁ − (1−w₁) log(1−w₁) ,
>   hence     w₁ ≥ h⁻¹_<(S)      (h⁻¹_< = the smaller root in [0,½]) ,     with equality iff w₁≤½.
> ```
> *Proof.* `ρ_A` has eigenvalues `Z,1−Z`; `S=−Z log Z−(1−Z)log Z`. `h` is strictly increasing on `[0,½]`,
> strictly decreasing on `[½,1]`, `h(0)=h(1)=0`, `h(½)=log2`; `h(w₁)=S` has the two roots
> `h⁻¹_<(S)≤½` and `1−h⁻¹_<(S)≥½`, and `w₁` is one of them, so `w₁≥h⁻¹_<(S)`. ∎

Witness §6·A confirms `w₁=Z`, `S=h(Z)`, `w₁≥h⁻¹_<(S)` to machine precision over a full sweep (both
branches `w₁≥½` and `w₁<½`).

> **Lemma U5-2 (anti-correlation) `[T]`.** `h` is **symmetric**, `h(Z)=h(1−Z)`. Consequently:
> (a) `S(A)=0` is compatible with `w₁=0` **and** with `w₁=1` — the single-mode entropy **cannot distinguish
> "no atom" from "perfect atom."**
> (b) the bound `w₁≥h⁻¹_<(S)` is **vacuous (`w₁≥0`) at both `w₁→0` and `w₁→1`** and non-trivial only near
> `w₁=½`, the *least* atomic configuration;
> (c) equivalently, **the bound strengthens as `w₁→½` (atom dilutes) and weakens as `w₁→1` (atom forms): it
> is anti-correlated with the object it is supposed to certify.**

Witness §6·A makes this quantitative: as the source decouples from the bath (`V↓`), `w₁→0.985` but the bound
falls to `0.015`; as the mixing grows (`V↑`), `w₁→0.502` while the bound rises to `0.498`. *The atom is largest
exactly where the entanglement bound is smallest.*

**The precise positive conditional — and why it is circular `[C]`.**
> If one *knows* the atom is the majority branch (`w₁≥½`), Lemma U5-1 gives the **equality** `w₁=1−h⁻¹_<(S)`:
> a genuine non-vacuous floor `w₁≥1−h⁻¹_<(S)` that increases with `S`. But the hypothesis `w₁≥½` is *stronger*
> than any target floor `w₁≥w₀` with `w₀≤½`, so it **already implies the conclusion** and cannot be used to
> prove it.

So the single-mode route is **either vacuous (no majority hypothesis) or circular (with it)**. There is no
free entanglement floor even in the solvable two-level model.

---

## 3 · The universal-floor no-go: free massless chain (large `S`, `w₁=0`)

Could one at least bound `w₁` from below by an entanglement monotone `f(S)` that is *not* `h⁻¹_<`? No.

> **Theorem U5-3 (no universal entanglement floor) `[T]`+`[V]`.**
> There is no decreasing function `f≥0` with `w₁ ≥ f(S(A))` for every theory, because the free massless
> field violates it for every `f` with `f(S)>0`: its bipartite entropy grows `∝ log l` (arbitrarily large
> with region size) while its `0⁺⁺` composite `¼trF²` has `w₁=f₁²=0` **exactly** (a two-particle cut, no
> atom; `T2 §4.3`, `ch06 §acid`).

Witness §6·B: free-fermion half-filled chain, `S(l)=1.419, 1.649, 1.878, 2.098, 2.286` for
`l=8,16,32,64,128` — **monotonically growing** — while `w₁(free 0⁺⁺)=0`. *Large entropy coexists with a
zero atom.* Any "entanglement → atom" inference must therefore invoke extra YM structure, i.e. the measure
`μ` — which is `H1`.

**Why the naive area-law move is exactly the trap.** The tempting argument "the `β=2.5` region is
area-law-entangled, `S` is large, so `h⁻¹_<(S)≈½` and `w₁≥½`" is *invalid*: `h⁻¹_<` is the **two-level**
(single-mode) inverse, and the gauge region is a **many-mode field** whose entropy is dominated by the UV
area term — identical for gapped and gapless theories at short distance (`H1-T4 §3`). The free chain is the
counterexample: maximal-looking `S`, zero atom. The area term is **gap-blind** (`H1-T4 §3`); only the `O(1)`
subleading piece carries any IR information, and no theorem converts subleading entropy into an atom floor.

**Area law status `[std]`.** A rigorous area law for interacting 4D YM is unproven; even the *gap-blind*
entropy is not available at `β=2.5` (Borgs–Seiler `β<0.706` vs `β≈2.5`; `H1-T4 §2`). Hastings' theorem is
the **converse** (gap ⟹ area law), so it cannot be run backwards (`IMPOSSIBLE-ENTANGLEMENT §1(ii)`).

---

## 4 · SSA forbids negative MI; the deficit is circular; relative entropy sees only diagonals

**(a) Negative mutual information is impossible `[T]`.** Mutual information is defined by SSA to be
non-negative; §6·C (4000 random pure 4-qubit states) finds **0** violations of `I(A:B)≥0`,
`S(ABC)+S(BCD)≥S(BC)`, and `|S(A)−S(B)|≤S(AB)`. So the task's "negative mutual information" *cannot occur*;
the meaningful object would be a **deficit relative to a reference**, `S(ρ_A‖σ)`.

**(b) The relative-entropy ("deficit") route `[NG]`/`[O]`.** Pinsker gives
`S(ρ_A‖σ) ≥ ½‖ρ_A−σ‖₁²`, so `‖ρ_A−σ‖₁ ≤ √(2S(ρ_A‖σ))`, bounding **diagonal** observable differences
`|Tr(ρ_A−σ)X| ≤ ‖ρ_A−σ‖₁‖X‖∞`. But:
- `w₁=f₁²/M₀` is an **off-diagonal coherence** between the vacuum `Ω` and the one-tube state `|G⟩`
  (different energy/particle sectors), *not* a diagonal population; the trace distance between two
  nearly-diagonal reductions does **not** control it.
- the reference `σ` must be the **tube state**, which is the very object `w₁` measures — so using it is
  **circular** (the "one route that would work — and it is circular" of `T2 §2.4`, in entanglement language).

**(c) The sign/direction obstruction `[NG]`.** Entropies, MI, and Rényi entropies are **concave, aggregate,
upper-bounding** functionals: they bound how *few* Schmidt states / how *small* the correlation can be. A
lower bound on a single spectral atom needs a **constructive overlap** (a positive recoupling / residue), not
an information aggregate. This is the **wrong-direction** no-go (`ch07 §wrong-direction`: "variational and
entropy methods give *upper* bounds; a lower bound needs a threshold oracle") and the **aggregate-vs-edge-
atom** sign wall (`H1-T4 §4`, `GAP-G2-edge`), now in entanglement coordinates.

---

## 5 · The decisive structural theorem: `w₁` is not a function of `ρ_A`

The deepest reason entanglement cannot carry the atom is a **basis-independence** obstruction.

> **Theorem U5-4 (operator-blindness of entanglement) `[T]`+`[V]`.**
> Every entanglement monotone is a function of `ρ_A` alone. But `w₁=w₁(ρ_A,O)=|⟨G_O|φ̂_O⟩|²` depends on the
> **operator** `O` through the tube ray `G_O` and the source `φ_O=QO|Ω⟩`. For a fixed `ρ_A` (fixed `S`),
> two operators `O,O'` give different `w₁`. Hence **no functional `F` with `w₁=F(S(A))` exists.**
> *Witness §6·D.* Fix `ρ_sys=diag(0.1,0.2,0.7)` on `{v,G,⊥}` (`S=0.80182`). Operator `O₁` with tube `=|G⟩`
> gives `w₁=0.200`; operator `O₂` with tube `=(|G⟩+|⊥⟩)/√2` gives `w₁=0.824` — **same `ρ_A`, same `S`,
> `w₁` differs by `0.62`.** ∎

**Reading.** The atom is a **coherence** in the operator-selected basis — the vacuum↔tube off-diagonal
`⟨G|O|Ω⟩` — while entanglement is a **basis-independent aggregate** of `ρ_A` that is invariant under
coarse-graining/mixing in that basis. A decoherence-robust quantity (entropy) cannot certify a coherence- fragile
one (the residue). This is exactly the **Schur-convex master obstruction** (`ch07 §schur-master`): "positivity,
reflection positivity, spectral positivity give upper bounds and the vacuum, never the edge." In YM the
operator `O=¼trF²` is *fixed*, so the lesson is sharper: the missing information is the **`O`-matrix-element
itself** (`f₁²`), an off-diagonal datum no `ρ_A`-functional contains. Bounding it from `ρ_A` alone is
impossible; bounding it from `ρ_A` **plus the `O`-alignment** re-introduces `f₁²` = circular.

**The 1-form/center wrong channel `[NG]`.** A natural entanglement diagnostic of a *confining* tube is the
MI between its two ends (or across a cut of the string). That measures the **1-form/center** sector (`σ>0`),
which is **Schur-orthogonal** to the `0⁺⁺` bulk scalar `trF²` (`IMPOSSIBLE-ENTANGLEMENT §3`). The decisive
countermodel is **2D YM**: `σ=g²C₂/2>0` (confining, area law, large end-to-end MI) yet `trF²=H` is conserved,
`⟨0|trF²|G⟩=0`, and `w₁=0` **exactly** (`R5 §1.1`, `T6 §4`). Witness §6·E reproduces the structure in a
solvable chain: a state with `MI=0.34>0` whose conserved density `N=Σc†c` has a `t`-constant correlator
(KL `=δ₀`, no atom at `m>0`), `w₁=0`. *Entanglement fires where the `0⁺⁺` atom is zero.*

---

## 6 · Witness — `research/v3/w7/u5_entanglement_binding.py` (`py -3.12 -I`, exit 0)

Deterministic (seeded), numpy only, ~0.7 s. Quoted output (abridged):

```
(A) single-mode exact identity w1=Z, S=h(Z), w1>=h^{-1}_<(S):
      V=0.02  Z=0.98473  S=0.07901  bound=0.01527   (atom strong, bound WEAK)
      V=0.70  Z=0.55611  S=0.68684  bound=0.44389
      V=20.0  Z=0.50198  S=0.69314  bound=0.49802   (atom diluted, bound STRONG)
    other branch: eps=5,V=.02 -> Z=0.00064 S=0.00534 ; eps=1,V=1 -> Z=0.46059 S=0.69004
    => S->0 at BOTH w1->0 and w1->1 ; bound vacuous at the atom; anti-correlated with the atom.
(B) free massless chain S(l): l=8:1.4187 16:1.6493 32:1.8778 64:2.0983 128:2.2858  (grows)
    while w1(free 0++) = 0 exactly  => NO universal floor w1>=f(S)  [counterexample]
(C) 4000 random pure 4-qubit states: I(A:B)<0 viol=0 ; SSA viol=0 ; Araki-Lieb viol=0
    (min I = 3.6e-2 ; worst SAB-|SA-SB| = -2.8e-1 ; worst SSA slack -1.4e-2 => all strictly satisfied)
(D) same rho_sys (S=0.80182): O1(tube=|G>) -> w1=0.20000 ; O2(tube=(|G>+|perp>)/v2) -> w1=0.82417
    => w1 is NOT a function of rho_A / S  (operator / basis dependence)
(E) 2D-analogue: correlated state MI=0.34175 (gap 4) / 0.09736 (gap 16) with a CONSERVED density
    (N=sum c^dag c commutes with H => C(t)=const => KL=delta_0 => w1=0)
```

Readings: **(A)** the identity is exact but **symmetric** — the bound is weakest exactly at the atom
(`w₁→1`) and strongest at dilution (`w₁→½`); **(B)** large `S` with `w₁=0` falsifies any `w₁≥f(S)`;
**(C)** SSA/subadditivity hold identically, so "negative MI" is impossible; **(D)** `w₁` is not a function
of `ρ_A` (two operators, same `S`, `w₁` `0.20` vs `0.82`); **(E)** `MI>0` coexists with `w₁=0` for a conserved
density (the 2D `trF²=H` analogue). **Witness, not proof; the theorems are §2–§5.**

---

## 7 · Ledger

| item | statement | label |
|---|---|---|
| single-mode identity | `w₁=Z`, `S=h(Z)`, `w₁≥h⁻¹_<(S)` (equality iff `w₁≤½`) | **[T]** |
| anti-correlation lemma | `h(Z)=h(1−Z)`: `S=0` for `w₁=0` **and** `w₁=1`; bound vacuous at the atom | **[T]** |
| conditional floor | `w₁≥½ ⟹ w₁=1−h⁻¹_<(S)`; hypothesis `w₁≥½` already implies any target `w₀≤½` | **[C]** (circular) |
| universal-floor no-go | free massless chain: `S∝log l` large, `w₁=0` ⟹ no `w₁≥f(S)` | **[NG]**+**[V]** |
| area law at `β=2.5` | unproven (Borgs–Seiler `β<0.706`); area term gap-blind | **[std]** |
| MI ≥ 0 / SSA / Araki-Lieb | identity; "negative MI" impossible | **[T]**+**[V]** |
| relative-entropy deficit | Pinsker bounds only *diagonal* differences; reference = tube = circular | **[NG]**/**[O]** |
| wrong-direction / aggregate-vs-edge | entropies upper-bound; atom needs constructive overlap | **[NG]** |
| operator-blindness | `w₁` not a function of `ρ_A`; same `S`, `w₁∈{0.20,0.82}` | **[T]**+**[V]** |
| 1-form/center channel | MI/confinement `=σ>0`; 2D YM `σ>0` yet `w₁(trF²)=0` | **[NG]** (`R5 §1.1`) |
| conserved-density analogue | `MI>0` with `w₁=0` (chain); 2D `trF²=H` (exact) | **[T]**+**[V]** |
| **uniform `w₁≥w₀` at `β=2.5`** | the binding target — not supplied by entanglement | **[O]** |
| **missing datum** | operator-specific coherence `f₁²/M₀` = residue = isolation `I₂` = `H*` | **[O]** |
| **Clay** | unbuilt | **[O] · 0%** |

---

## 8 · The exact missing estimate and the falsification test

**Missing estimate (stated once, exactly).** An **operator-specific coherence floor**:
> There exist `w₀>0`, `δ>0` such that, uniformly in `(L,a)` through `β=2.5`,
> ```
>   w₁ = f₁²/M₀ = |⟨G|O|Ω⟩|² / ‖QO Ω‖² ≥ w₀ ,    equivalently  Σ_{n≥2} w_n ≤ 1−w₀ ,
> ```
> the **off-diagonal weight of the isolated pole in the source's spectral measure** — a *coherence* in the
> operator-`O` basis, **not** a function of `ρ_A` and hence invisible to every entanglement monotone.
> By `T2 §8` / `L2` this is the four-wall polymer residue `=` isolation `I₂` `=` `H*`.
> Entanglement contributes only the **symmetric** part `h(w₁)` (one bit: `w₁` vs `1−w₁`) and its
> **anti-correlation** with the atom; the missing datum is the **phase/selection of the `O`-basis** = the
> residue sign/amplitude, which is basis-dependent and information-theoretically **lost** under the
> coarse-graining that defines `ρ_A`.

**Falsification test (cheap, pre-registered, decisive for the no-go).**
Measure on SU(2) `β=2.5`, `L=8,12`, the same `0⁺⁺` `trF²` correlator as `R4` (GEVP + two-state fit) to
extract `w₁(L,T)` (as in `T2 §8`), **and** the region-complement entanglement entropy `S(A)` (a flux-tube-
region slab; proxy: the Rényi-2 purity `Tr ρ_A²` from a replica / free-field-limit estimator, subtracting the
UV area term `S_area`).

- **This note predicts (falsifiable):** `w₁` and the **non-area** entropy `S(A)−S_area` are **anti-correlated**
  (negative correlation coefficient). The bound-atom regime (`w₁→1`) has `S(A)−S_area→0` (a single mode
  decouples ⟹ its reduced state purifies; Lemma U5-1/2); the entangled regime (`S(A)−S_area` large) has `w₁`
  near the symmetric value. *A positive correlation, or a config with `w₁≥½` **and** `S(A)−S_area ≥ ½log2`
  simultaneously, **falsifies** the anti-correlation/no-floor picture and would mean YM's multi-mode structure
  beats Lemma U5-2.*
- **Acid controls:** the same pair `(w₁, S(A))` must return `w₁=0` at **any** `S(A)` on (i) a free-field /
  `β=0`-Haar ensemble and (ii) 2D YM (`trF²=H` conserved). If the estimator returns `w₁>0` there, it is
  reading a carrier/aggregate artifact, not the atom (`ch06 §acid`, `ch07 §schur-master`).
- **What a positive result would require:** a functional of the **joint** `(ρ_A,O)` — i.e. a coherence
  measure anchored to the physical `O=¼trF²` — not of `ρ_A` alone (Theorem U5-4). Constructing one *is*
  constructing the residue, i.e. the polymer four-wall recoupling, i.e. the binding theorem.

---

## 9 · Sharpest line

**Entanglement is the wrong species for the binding theorem, and the reason is now a theorem, not a slogan.**
The one exact link — the single-mode identity `w₁=Z`, `S=h(Z)` — is **symmetric** (`h(Z)=h(1−Z)`), so the
single-mode entropy cannot tell a perfect atom (`w₁=1`) from the free field (`w₁=0`); its lower bound
`w₁≥h⁻¹_<(S)` is **vacuous exactly in the atom regime** and strengthens only as the atom dilutes toward `½`
— it is **anti-correlated with the object it should certify**. No universal floor `w₁≥f(S)` can exist, because
the free massless chain realizes *large entropy with a zero atom* (`S∝log l`, `w₁=0`). Strong subadditivity
makes mutual information **non-negative**, so "negative MI" is impossible and an "entanglement deficit" must be
a **relative entropy** against a reference tube state — the atom itself, hence circular. And the deepest
obstruction is **operator-blindness**: `w₁` is the **coherence** `|⟨G|O|Ω⟩|²/‖QO Ω‖²` in the operator-`O`
basis, and for a fixed `ρ_A` (fixed `S`) two operators give different `w₁` (`0.20` vs `0.82`) — an entropic
aggregate of `ρ_A` cannot lower-bound a basis-dependent coherence (the Schur-convex master obstruction, in
entanglement coordinates). The natural confining-tube MI diagnostic lives in the **1-form/center** sector and
is killed by **2D YM** (`σ>0`, `f₁²(trF²)=0`). The exact missing estimate is therefore an **operator-specific
coherence floor** `w₁≥w₀` — the polymer four-wall residue `=` isolation `I₂` `=` `H*`. The entanglement route
**re-expresses the wall and adds a sharp no-go; it does not lower it.** **No QED. Clay 0%.**

---

*Refs.* **This wave:** `T-WAVE-SYNTHESIS.md`, `T2-TUBE-OVERLAP.md` (§1.3 `w₁`, §2.4 circular L2), `T6-BINDING-
ADVERSARIAL.md`, `SM5-TCI-HEAD.md`, `H6-HEAD-SPECTRAL.md`. **Paper:** `docs/paper-v2/ch06-tools.tex` (§acid,
§toolwave — aggregate upper bounds), `ch07-nogos.tex` (§schur-master, §impossible `IMPOSSIBLE-ENTANGLEMENT`,
wrong-direction/wrong-regime). **Entanglement base:** `unconditional/IMPOSSIBLE-ENTANGLEMENT.md` (edge vs bulk,
Li–Haldane; area law not ⇒ gap; functional of `μ`), `unconditional/H1-T4-entanglement.md` (gap ⇒ area law one-
way; area term gap-blind; aggregate-vs-edge-atom). **Corpus anchors:** `R5-EXACT-ANCHORS.md` (§1.1 2D
`trF²=H` ⇒ atom `0`; §5 scaling), `GAP-G2-edge.md` (residue–level-gap L2), `R2-ISOLATION.md`, `Q6-H1-STATUS.md`.
**Standard (textbook, used as such, not re-fetched):** Lieb–Ruskai (SSA, 1973); Araki–Lieb (1970); Pinsker /
Fannes–Audenaert continuity; Fuchs–van de Graaf (fidelity–trace); Hastings JSTAT (2007) P08024 (gap ⇒ area law);
Borgs–Seiler (area law `β<0.706`); Bombelli–Koul–Lee–Sorkin (1986), Srednicki (1993), Casini–Huerta (2009)
(free-field area laws); Li–Haldane PRL **101** (2008) (entanglement spectrum = edge). **Witness:**
`research/v3/w7/u5_entanglement_binding.py` (`py -3.12 -I`, numpy only, deterministic, exit 0). **No third-party
downloads; no `[C]` upgraded to `[T]`; no QED; Clay 0%.**
