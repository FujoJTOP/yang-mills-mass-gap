# BASE-INDEX — the gap is an **order parameter, not an integer index**: homotopy-rigidity + atom/gap type mismatch + acid failure ⇒ the index class is floor-blind

> 2026-11-24 · Base-round piece: does the missing **one sign** (SLOT-SYNTHESIS) package as a **global integer index** — spectral flow of the adjoint Dirac operator, the eta invariant / APS index, or the Callias index — evaluated in the presence of the would-be mass? Witness `research/unconditional/w7/base_index_witness.py` (`py -3.12 -I`, tiny finite matrices, logical structure only).
> Labels (`CLAUDE.md §2`): **[T]** true theorem / banked · **[C]** conditional · **[O]** exact missing step · **[NG]** no-go for a *named class* · **[V-here]** new this round (incl. witness) · **[R]** inherited/corpus.
> **Does NOT prove Clay; does NOT prove no index can ever exist (`∀class ≠ ∀possible`). No QED.**

---

## VERDICT (up front)

> **NO — the gap is not an integer index.** Four independent [T]-level obstructions, any one fatal:
> **(I1) Rigidity/continuity.** An index is a **homotopy invariant**: integer-valued and **constant on each connected component** of the operator-family parameter space. The coupling `β` is a continuous parameter and the gap is a **non-analytic `β`-jump** (corpus **B3**: "discrete indices frozen at the free value; the gap is a `β`-non-analytic jump"). A continuous ℤ-valued invariant on the connected interval `β∈[0.5,4]` is **constant** — it **cannot** jump at `β_c≈2.5`. It is blind to the crossover by construction.
> **(I2) Atom/gap type mismatch.** An index counts **zero modes / vacua** (ground-state, *atom-side* data — Witten index, Callias zero modes, APS obstruction). The gap is the **first excited level above** the ground state (*isolation*, `H*`). The corpus already proved **atom ⇏ gap** and **gap ⇏ atom**; an index certifies at most atom-side structure and is **blind to the gap**.
> **(I3) Zero-mode ≠ positive-gap.** Spectral flow changes **only when the gap closes** (a level crosses zero). A *positive* gap `g>0` carries the **same** index for every `g` (witness W2: flow `=0` for `g∈[10⁻³,10]`). So an index can certify **gaplessness** (gap `=0`), never **`gap>0`** — the wrong direction for us.
> **(I4) Acid failure.** A **nonzero** topological index (instanton charge `∫tr F∧F`, winding) **coexists with a massless mode** (gap `0`) — it fails the discriminating acid test (witness W4). On the free massless field the index is `0` (vacuous, as required — but `0`-index is compatible with *both* gapped and gapless, so it carries nothing).
> **⇒ The gap is an _order parameter_ (real, continuous within a phase, jumps at the transition), NOT a _topological invariant_ (integer, constant on components).** The index class is floor-blind. **Exact failing step:** the only index-type route that *can* force a gap is **'t Hooft anomaly matching** (an 't Hooft anomaly forbids a unique gapped symmetric vacuum), and 4D pure YM has **no such anomaly** — the confined vacuum is unique+gapped with the 1-form center symmetry unbroken (area law), so anomaly matching neither forces nor forbids, and *assuming* it would assume the structure to be proven. **Clay 0%.**

---

## §1 · What the three candidate indices actually count

- **Spectral flow** (Atiyah 1985, *Global theory of elliptic operators*; Phillips 1996, *Self-adjoint Fredholm operators and spectral flow*). For a continuous path of self-adjoint (Fredholm) operators, the spectral flow is the net number of eigenvalues crossing `0`; it is an **integer**, a **homotopy invariant** of the path, and equals the index (in the APS sense) of the corresponding Dirac operator on the mapping cylinder. `[T]` (witness W1: flow `=1` for one crossing; unchanged under a continuous off-diagonal perturbation).
- **Eta invariant / APS index** (Atiyah–Patodi–Singer 1975, *Spectral asymmetry and Riemannian geometry I*, Math. Proc. Cambridge Philos. Soc. **77**, 43). `η(H)` is the spectral asymmetry (real, regularized `Σ sgn λ`); `η(0) mod 2 ∈ ℤ`; along a family `∂η/∂t = 2·spectral flow`; the APS index `ind(D⁺) = ∫ ch − (η+h)/2` is an **integer**. It is a property of a **Dirac-type operator with a given mass/background**.
- **Callias index** (Callias 1978, *Axial anomalies and index theorems on open spaces*, CMP **62**, 213). On odd-dimensional (e.g. 3D) open space with a mass/Higgs matrix `Φ(x)` **invertible at spatial infinity** and interpolating between asymptotic limits, the index counts zero modes of the Dirac operator; **integer**, topological. It counts the **spectral flow between two asymptotic "vacua"** of a *given* mass profile.

**Common structure.** All three are **integer-valued homotopy invariants of an operator family**, and all three count **zero modes** of a **Dirac-type** operator with a **given** mass/background. None is a functional of a **bosonic** `0⁺⁺` correlation, and none is defined absent a given mass.

---

## §2 · (I1) Rigidity: an integer invariant cannot jump at the crossover

The gap in this campaign is a **`β`-non-analytic object** (corpus **B3**; `PROFILE` P2/P4: the species is `usc`, supported on the closed-nowhere-dense gapped set `{m₁≥δ}`, vanishing on the dense gapless class). An index is the **opposite kind of object**:

> **A continuous ℤ-valued function on a connected interval is constant** (ℤ discrete, interval connected). Witness W3: any candidate integer built from the gap (`round Δ`) is **discontinuous** — a jump proves it is **not** a homotopy invariant; a genuine invariant is forced **constant** on `β∈[0.5,4]`.

So a `β`-index `I(β)` that is one integer at `β=0` (strong coupling, P3-gapped) and the *same* integer at `β→∞` (weak coupling) is **constant in between** and **cannot** change at `β_c≈2.5`. It is the **B3 rigidity no-go** exactly: the discrete index is frozen at the free value, while the gap is a jump no invariant tracks.

---

## §3 · (I2) Atom/gap type mismatch

An index counts **zero modes / vacua** — the `Tr(−1)^F`-weighted **ground-state** sector (Witten 1982, *Constraints on supersymmetry breaking*, NPB **202**, 253: robust, `β`-independent — precisely because it is **gap-blind**). The gap is the **first excitation above** the ground state. The corpus established both directions are independent: **gap ⇏ atom** (free massive `:φ²:` is gapped with `f₁²=0`) and **atom ⇏ gap** (a thin atom can sit at a near-zero mass without isolation). An index lives on the **atom side** and inherits **atom ⇏ gap**: it can be nonzero while the gap is `0` (W4) and zero while the gap is `>0`.

**Consequence.** Even a *successful* index computation would certify, at best, atom/ground-state structure — not `H*` (atom **+ isolation**). The campaign's `H*` needs the **isolation**, which is a first-excitation datum no ground-state index sees.

---

## §4 · (I3)+(I4) Zero-mode blindness and the acid test — the decisive failures

**W2 (positive-gap blindness).** For the gapped family `H(s)=[[s,g],[g,−s]]`, eigenvalues `±√(s²+g²)`, the spectral flow over `s∈[−1,1]` is `0` for **every** `g∈{10⁻³,…,10}` while the gap `2g` varies over **four decades**. The flow changes **only** at `g=0` (the crossing). ⇒ **An index detects `gap=0` (zero modes), never `gap>0`.** For the Clay problem (certify `gap>0`) this is the **wrong half** of the dichotomy.

**W4 (charge ≠ gap).** A winding-`1` background has topological charge `Q=1` (an integer, computed exactly) **independent of the mass parameter `m∈{0,0.1,1,5}`**, and coexists with a **massless mode** (Dirac zero mode, gap `0`). Real instanton number `∫tr F∧F /8π² ≠ 0` is therefore **gap-blind** — it is a UV/topological datum present in gapless-compatible configurations (dilute-instanton gas is photon-like/gapless at weak coupling). **Acid FAIL**: the nonzero index does not vanish on the gapless adversary.

**Eta/Callias add nothing and worsen two ways.** (i) `η` is **real**, an integer only `mod 2` — even weaker than an integer index; on the symmetric free field `η=0` (vacuous). (ii) **Callias is circular**: it needs a **mass matrix `Φ(x)` invertible at infinity** as *input* to count the flow that would certify the mass — but the **existence of the mass is the unknown** (this is the A3/BS circularity: *read, not force*).

---

## §5 · The one honest loophole, and why it fails for pure YM

There **is** a legitimate way a *topological*/index argument can force a gap: **'t Hooft anomaly matching** — an 't Hooft anomaly in the global symmetry forbids a **unique, trivially gapped, symmetric** vacuum; the theory must either break the symmetry spontaneously or have a massless mode. This is an index-type statement that *does* bear on the gap.

**For 4D pure `SU(N)` Yang–Mills this does not apply (at `θ=0`, the mass-gap setting).**
- Pure YM has **no chiral symmetry** and hence no standard 't Hooft anomaly. The relevant symmetry is the **1-form center symmetry `ℤ_N`** (lattice), and in the **confined** phase it is **unbroken** (area law) with a **unique gapped vacuum** — *consistent* with there being no anomaly obstructing that vacuum. The 1-form symmetry is the **order parameter** for confinement, not an anomaly that forbids a gapped vacuum.
- Generalized anomalies of the 1-form symmetry **do exist** but are about the **`θ`-dependence of the phase diagram** (e.g. a mixed `ℤ_N^{(1)}`–`CP` anomaly constraining `θ=π` to be non-trivially gapped or a transition point). They **do not forbid** a unique gapped symmetric vacuum at **`θ=0`** — so they neither force nor forbid the gap we seek; invoking them for the `θ=0` mass gap would be a category error.
- Using anomaly matching to *force* the gap at `θ=0` would require assuming the confined/center-unbroken structure — i.e., assuming the very thing to be proven.

So the anomaly loophole is **empty for pure YM** (it is the correct tool for QCD-like theories with chiral symmetry; not for the pure gauge theory).

---

## §6 · What this means for the base

The index language, made precise, **isolates the structural obstruction**: **the gap is an order parameter, not an invariant.** This is a **positive** base contribution in the negative-space sense — it **rules out the entire topological-invariant class** and tells the base *what kind of object the supporting structure must be*:

- The supporting object is **not** integer/topological/homotopy-invariant (that class is constant on `β` — B3).
- It **must** be a **continuous, phase-distinguishing functional** of the **field measure** (the P1 door): real-valued, `usc`, supported on the gapped set, non-constant within the
 coupling axis. Equivalently: a **non-topological order parameter**.
- The base therefore belongs to the **field-measure-functional** family (`PROFILE` P1–P4), not the index family. This **redirects base construction away from topology** (instantons, spectral flow, eta, Callias, Witten) — all now provably floor-blind — toward **measure-ordering** functionals.

This is consistent with, and sharpens, `SLOT-SYNTHESIS §3`: the irreducible core is **one sign**, and that sign is a **phase/order-parameter** sign, not a **topological charge**. A topological charge is **invariant**; a phase sign is **order parameter** — and the gap is the latter.

---

## §7 · Verdict, labels, exact failing step

- **[T]** Spectral flow is integer + homotopy-invariant (Atiyah 1985; Phillips 1996; APS 1975); Callias index integer (Callias 1978); Witten index `β`-independent (Witten 1982). Witness W1.
- **[T]** Continuous ℤ-valued function on a connected interval is constant ⇒ any `β`-index is constant across the crossover (B3 rigidity). Witness W3.
- **[T/V-num]** Spectral flow is blind to a **positive** gap, changes only at gap-closure (W2); topological charge is independent of the mass and coexists with a zero mode (W4).
- **[NG]** *No-go for the index class:* no integer index of a Dirac-type operator family certifies `gap>0` for a bosonic `0⁺⁺` excitation of 4D pure YM.
- **[O]** *Exact failing step:* the gap-sign is an **order parameter**, not an invariant; the only index route that bears on gappedness is 't Hooft anomaly matching, **absent** for pure YM; a *non-local* index of a *family of measures* would need a closed family (the `β`-line is not) and a qualitative detector = **circular = A3**.
- **Acid test: FAIL** (nonzero charge coexists with gap `0`; free-field index `=0` is vacuous). **IL-i** fail (blind to positive gap / gappedness); **IL-ii** fail (Witten/Callias indices exist for `U(1)`); **IL-iii** the rigidity obstruction is independent of `α_s` — the wrong tool.

**Clay = 0%. No QED.** `∀class ≠ ∀possible`: this rules out the *index class*, not a hypothetical non-local, non-integer, measure-family invariant-*like* object (which would be the P1 door, not an index).

*Refs:* `SLOT-SYNTHESIS.md` (the one sign), `PROFILE.md` P1–P4, `docs/new-tools-registry.md` **B3** (rigidity), **A3** (Fredholm read-not-force); corpus `INVENT-VARIANCE-INEQ`, `MS-GAP-EQUATION`. Witness `w7/base_index_witness.py`. External: Atiyah–Patodi–Singer 1975; Atiyah 1985; Phillips 1996; Callias 1978; Witten 1982; 't Hooft 1980 (anomaly matching).
