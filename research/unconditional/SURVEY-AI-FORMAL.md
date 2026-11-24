# SURVEY-AI-FORMAL — Sweep & audit of AI-generated / formalized / computer-assisted Yang–Mills-gap attempts

> **2026-11-24.** Scope: (a) Lean/Coq/Isabelle formalizations of the 4D YM mass gap or its prerequisites (G1 existence, OS axioms, area law); (b) AI-generated / AI-assisted proof attempts (2025–2026 arXiv, Zenodo, GitHub, blogs); (c) computer-assisted / interval-arithmetic / verified-numerics gap extractions for 4D gauge.
>
> **Labels (charter §2).** `[T]` 真定理 · `[C]` 条件定理 · `[NG]` no-go（仅对具名类） · `[O]` 精确缺失步 / open · `[R]` 外部据报 / 未独立核实 · `[J]` 判断（非定理） · `[V]` 本会话亲自核实（读取一手 PDF / 仓库 / 委员会声明原文） · `[V-prior]` 本项目前序审计已核实（见 `CHATGPT-METHOD-RGENTRY.md`）。
>
> **Iron rules.** No "QED" on anything unproven. Every claim typed. Free-massless acid test applied where a "new tool" is claimed. Unverifiable / self-published-only sources flagged explicitly.

---

## 0. TOP-LINE（先读这一段）

**Q1 — Is there any AI-generated or formalized 4D YM-gap claim that survives audit?**
**No.** After sweeping the 2025–2026 landscape and reading the primary artifacts, **no** AI-generated or machine-checked source establishes a 4D Yang–Mills **mass gap** (G2). The landscape splits cleanly:

- The one *serious, machine-checked* AI result — OpenAI's Lean 4 construction, reviewed by the Gross-chaired committee (Douglas, Jaffe, Klebanov, Nekrasov) — is **G1-only (existence)**. It **explicitly excludes the mass gap, uniqueness, and the area law**: "It does not prove that the infinite-volume limit is unique or that there is a mass gap"; "The construction provides no evidence for or against a mass gap in the theory." **[V]** (assessment PDF, KITP statement, repo README). It is a real, substantial **partial** — not a gap claim.
- The one *serious 4D mass-gap claim* — **Faizal–Shabir 2026** (arXiv:2606.19362, *Fortschr. Phys.* 74 e70097) — is **human-written, published, and conditional**. It is not AI-generated, and it dies at the same crossover wall (**H1**): the AF↔IR "common polydisc" is *assumed*, not proved. **[V-prior] + [V] this session.**
- Every explicit "full solution" Lean/Coq claim (jonwashburn, Farooqui, "Synthesis Reasoning/ASI", merchantmoh, the Zenodo "Rigged-Hilbert scaffolding") is either **non-compiling**, **definitionally vacuous**, **conditional on an open hypothesis**, or **unverifiable/self-published overclaim**. **[V]/[R]**, itemised below.

**Q2 — Is there a reusable *new technique* buried inside?**
**Yes, at the infrastructure / 2D level; No, at the 4D-gap-crossing level.**
- **Reusable (real):** (i) a verified Lean infrastructure for OS axioms, reflection positivity, chessboard estimates, and the forest formula (Douglas–Mei–Yin `OSreconstruct`; Hoback–Douglas; `mrdouglasny/lgt`); (ii) the **2D decomposition template** — "the RG shows a gap *test* is cutoff-independent; one crude estimate evaluates it once" — with the O(N) **twist → stiffness ≤ C²/log L → exponential decay** theorem actually formalized in Lean 4 (117 Mathlib files) **[V]**; (iii) **verified-numerics** gadgets: the 4D **box-doubling criterion Δ₄** and an interval-arithmetic **Hessian certificate** (h_F ≤ 14.420001). These are genuine and transferable.
- **Not reusable:** no source supplies a **4D lower-bound tool that crosses the `α_s~O(1)` crossover** and is `(L,a)`-uniform. Douglas states the gap in exactly the project's words **[V]**: *"For YM4 the crude estimate must come from the strong coupling regime, so the RG has to flow from weak coupling to couplings of order one and hand off — a regime the sigma model papers never treat."*

**Q3 — The most serious claim and its precise missing step.**
**Faizal–Shabir 2026.** Its proof kernel hinges on Theorem 10.8's premise *"at a common intermediate scale k large enough that both trajectories lie in the contracting polydisc"* **[V, quoted]** — i.e. a common UV/AF and IR/strong-coupling polydisc. That is **O2-Bridge = the crossover = H1**, assumed not proved; compounded by a **parameter-identity error** in Lemma 10.5 (`g_bare→0` is *weak* coupling, mislabelled "strong-coupling sector") and a **non-sequitur** at Eq (10.72) (`Σ|ε_k|<∞` while `ε_k ⊇ O(g_k²) ~ 1/k`, harmonic). **Verdict: [C] conditional, not a proof; Clay still 0%.**

**One-line verdict.** *No AI/formal 4D mass-gap claim survives audit; the only clean machine-checked AI result is G1 (existence) and says so, and it places the remaining gap exactly where the project does — across the crossover (H1).*

---

## 1. Audit frame: G1 vs G2 vs H1 (so the table can be read)

- **G1 = existence**: a Euclidean measure / correlation functions satisfying OS (⇒ Wightman) axioms on `ℝ⁴`, non-trivial, with the right UV scaling. A **subsequential** infinite-volume limit via compactness is *not* the Jaffe–Witten bar unless additional properties of the limit are established.
- **G2 = gap**: `spec(H) ⊆ {0} ∪ [Δ,∞)`, `Δ>0`; equivalently uniform exponential clustering of all local gauge-invariant observables; in the project's residue form `f₁²=|⟨0|trF²|G⟩|²>0` **and** isolated.
- **H1 = the wall**: the `(L,a)`-uniform, sign-controlling lower bound across the Savvidy crossover `α_s~O(1)`. Every route that claims the gap must, at some equation, produce this. A claim that "both banks are the same" or "the strong-coupling basin is reached" is assuming H1.
- **Audit test used below:** does the claimed artifact *establish the uniform `(L,a)` gap*, or does it hide (i) G1-only / (ii) strong-coupling only / (iii) a non-sequitur / (iv) the H1 wall / (v) definitional smuggling?

---

## 2. Master table

| # | Source (verified) | What is ACTUALLY done | 4D? | Audit verdict | Missing step |
|---|---|---|---|---|---|
| 1 | **OpenAI YM construction, Lean 4** — repo `yang-mills-group/OpenAI-YM-20261028`; committee **assessment 10 Nov 2026**; KITP 3rd statement **[V]** | Constructs continuum correlators as scaling limit of 4D Wilson LGT; UV convergence at fixed volume; **subsequential** infinite-volume via Arzelà–Ascoli; phase decomposition; OS reconstruction ⇒ Wightman; non-free; stress tensor/Ward/OPE | **4D yes** | **Real G1 partial [V]/[C]**. **Explicitly no gap, no uniqueness, no area law.** | G2 entirely; uniqueness/clustering; area law; universality of spacing sequences |
| 2 | **Gross-chaired validation panel** (Douglas, Gross, Jaffe, Klebanov, Nekrasov) **[V]** | Kernel-checks (1) above; 9,470 modules, 250,319 declarations, no project axioms / `sorryAx`, only `propext/choice/Quot.sound`; Comparator passed; independent OS cross-check | 4D | **Confirms G1-only**; "We do not claim to fully understand the proof yet" **[V]** | Validation ≠ truth of definitions; gap excluded |
| 3 | **Sussman–Alpöge (Claude-assisted)** `ethanwsussman/Yang--Mills-LLM-experiment` **[R]** | YM4 "zeroth draft" existence assuming Balaban's inputs; QED3/YM3 OS manuscripts; "slop, extremely long, basically not understandable" (Alpöge) | 4D (draft) | **G1-only draft [R]**; explicitly no mass gap/uniqueness/clustering | Same as #1; unvetted |
| 4 | **Faizal–Shabir 2026** arXiv:2606.19362 / *Fortschr. Phys.* 74 e70097 **[V]** | Reflection-positive lattice → transfer gap → area law → multiscale transport → OS ⇒ Minkowski gap + confinement; claims continuum + universality | **4D yes** | **Serious but [C] conditional.** Dies at crossover. Deep audit §3 | **Thm 10.8 common-polydisc (H1)**; Lemma 10.5 parameter-identity; Eq 10.72 non-sequitur |
| 5 | **`mrdouglasny/lgt` Lean 4** (Lean Reservoir) **[V]** | `ym_mass_gap_exponential_decay`, `ym_mass_gap_rate_exists` for U(n) Wilson LGT on `(ℤ/Nℤ)^d`, **Dobrushin**, `β<1/(32n(d−1))`; d≥2; 0 sorry | 4D **at strong coupling** | **True strong-coupling theorem [T]**; **not continuum, not crossover** | Weak coupling / continuum / crossover |
| 6 | **Farooqui "Machine-Verified Proof"** Coq, repo `Shariq81/yang-mills-mass-gap` **[V]** | Claims "mass gap at any β>0"; but explicit `m=β/10−4 for β>50`; README admits conditional on **Balaban pointwise convergence, open for YM_4**; CHANGELOG "0 hypotheses" achieved by **encoding Wilson suppression in a Coq `Record`** | 4D claimed | **Overclaim + definitional smuggling.** §5.1 | Continuum limit; correct scaling; un-smuggled hypothesis |
| 7 | **`jonwashburn/ym-proof` Lean 4** (Zulip) **[V]** | Claims sorry/axiom-free full Clay proof ("Recognition Science") | 4D claimed | **Does not compile** — `lake build` "errors about things not existing" (Wieser); "still does not compile at all" (Rothgang) | Everything |
| 8 | **Zenodo 23788087 v4.9.0**, Bhattacharjee–Bhattacharya–Lockwood **[V]** | Lattice transfer-matrix gap: SU(2) to `β_W≤0.345`, SU(3) to `β_W≤0.76`, general G below curvature threshold; 86 Lean theorems, 0 sorry; conditional continuum via **unproved** Dobrushin–Shlosman cube estimate | 4D lattice | **Honest partial [T]; explicitly open** ("continuum … remains open", "crossover remain open") | Cube estimate; continuum; crossover |
| 9 | **Zenodo 24450219**, Rayn-McLean **[V]** | MC box-doubling `Δ₄` on 2⁴–16⁴; constants audit of Alpöge–Sussman; Lean bookkeeping; strong-coupling anchor `Δ₄(16)≤2⁻¹⁰` for `β_W≤0.01273` | 4D lattice | **Not a mass-gap proof** ("It does not solve the mass gap"); self-published; neutral-sector open | Continuum; neutral sector |
| 10 | **Zenodo 22904674 "Rigged-Hilbert scaffolding"**, Norris et al. **[V]** | 108 Lean theorems, 0 sorry: Gel'fand triplet, Riesz embedding, "generalized spectral decomposition", asserts `Δ>0` on vacuum complement | 4D claimed | **Definitional/structural**, no dynamics shown; "scaffolding", "reduction lemmas" | Actual dynamics; gap derived, not encoded |
| 11 | **Zenodo 22767474 "Lean 4 Formalization Program" v54** **[R]** | ~114 modules, ~35k lines, small-β finite-volume Lipschitz stability; **states thermodynamic/continuum/mass-gap remain open** | 4D | **Partial [R]**, self-described incomplete | Continuum; crossover |
| 12 | **OpenAI 2D O(3) five preprints** (Sep–Oct 2026) **[V via Douglas slides]** | Problem A: O(N), N≥3, all β, exponential decay (twist/commutator trick); Problem B: O(3) continuum via Balaban-style RG (~22k lines) | **2D** | **2D, not 4D.** A = formalized [V]; B = "**Plausible, unaudited, and not formalized**" **[V]** | (2D) independent audit; (4D) whole crossover hand-off |
| 13 | **OpenAI O(N) decay theorem, Lean 4** (Problem A) **[V]** | `Z(θ)/Z(0)→1` ⇒ stiffness `≤C²/log L` ⇒ Peierls ⇒ exp. decay; 117 Mathlib files, std axioms; rebuilt/checked (Douglas) | 2D lattice | **Genuine new theorem, formalized [T]** | Nothing for 2D-A; irrelevant to 4D gap |
| 14 | **Zenodo 17083731** Ancri (Claude Code/Cursor) **[V]** | 97-page AI "exploration"; author "not a mathematician"; "does not constitute a peer-reviewed proof"; "a central argument that assumes what it needs to prove" | 4D claim | **Self-caveated, not a proof** | Everything; circular by admission |
| 15 | **"Synthesis Reasoning/ASI" Zenodo 19432416; merchantmoh "Ark Research"; hexagonmath "every compact simple group β≈0.76"** **[R]** | "Millennium #2 SOLVED", "10/10 verification", etc. | 4D claimed | **Overclaims / unverified** [R]; hexagonmath finite-lattice cover only | Continuum; crossover; or non-existent |

> **Isabelle:** no Isabelle formalization of any part of YM / OS / area law was found. Negative result.

---

## 3. Deepest audit A — Faizal–Shabir 2026 (the serious 4D mass-gap claim)

> Provenance: **arXiv:2606.19362v1** [hep-lat], 9 Jun 2026, cross-list hep-th; **Fortschr. Phys. 74 (2026) 4, e70097**, DOI 10.1002/prop.70097; companion papers *Int. J. Geom. Meth. Mod. Phys.* (2026) 2650111–2650114. **[V]** (abstract + full PDF re-extracted this session). 593 pp. Published but contested; **not** peer-certified as a solution. This is a **human** paper, not AI-generated — it is the serious 4D *gap* claim in the field, included because the task is to audit any 4D claim.

### 3.1 What is genuine (do not throw it away)
- §3 strong-coupling character/polymer convergence at `β<β₀` (fixed `a`) — essentially Osterwalder–Seiler, **[T]** in its domain.
- §5 interlacing inequality `T_{k+1}=V_k T_k^b V_k − D_k + E_k`, `D_k≥0`, `‖E_k‖≤ε_k` — a genuine **conditional** persistence-of-gap theorem: *if* `Σε_k<∞`, the gap survives the multiscale step. **[C]**
- Thm 8.6 Tauberian gap from exponential clustering, using **`H-RS(vac)`** (vacuum-representation Reeh–Schlieder) — correctly scoped; `H-RS(all)` is stronger and not required by Clay. **[V]**
- The paper is **honest about scope** in its companion Part 1 (§3.8 below).

### 3.2 Thm 10.4 — local, not global (paper admits it) **[V, quoted]**
> "Theorem (10.4) does not by itself guarantee that an arbitrary admissible trajectory will enter the polydisc `D(ε,ρ)`."

So the contraction is **local** near the Gaussian fixed point. Entry is a separate lemma.

### 3.3 Lemma 10.5 — parameter-identity error **[V, quoted]**
The entry lemma requires `0 < g_bare ≤ ḡ` (small bare coupling) and concludes `(g₀,K₀) ∈ D(ε/2,ρ/2)`. Its proof says:
> "By construction **in the strong-coupling sector**, `(g₀,K₀)` is finite and reflection-positive; moreover, **as `g_bare → 0`** one has `g₀ → 0` and `K₀ → 0` …"

**Flag.** `g_bare → 0` is the **weak-coupling / UV** limit (in Wilson normalization `β_W = 2N/g²`, so `g_bare→0 ⇔ β_W→∞`). Section 10 is explicitly "Weak-Coupling Extension and Asymptotic Freedom"; its polydisc is the **Gaussian/UV** polydisc. The strong-coupling base of §3 is `β<β₀` (large `g`), a **different, disjoint** region. So Lemma 10.5 puts a **UV** trajectory into the Gaussian polydisc; it does **not** put the effective action at the IR/crossover scale into a strong-coupling KP basin. The phrase "strong-coupling sector" is a **parameter-identity error**, and the paper uses "bare coupling" inconsistently across §3 and §10.

### 3.4 Theorem 10.8 — the crossover assumed **[V, quoted]**
> "We compare the asymptotically free trajectory … with any admissible trajectory that realizes the infrared construction of Sections 5–8 … **At a common intermediate scale `k` large enough that both trajectories lie in the contracting polydisc**, we write the difference of their one-step kernels as `T_k` …"

**This premise is assumed, not proved.** The AF trajectory (small `g`, weak coupling) and the IR construction (strong coupling, `β<β₀`) are **disjoint** in the known rigorous regimes. To make them share a polydisc one must run the strong-coupling trajectory to the UV — exactly the unproven transport. The "common intermediate scale" **is** the crossover **is** H1.

### 3.5 Lemma F.10 — the "same bare action" requirement **[V, quoted]**
> "Lemma F.10 (UV alignment). Let `Ω₀, Ω₁` be admissible regulators **based on the same local bare action** (for instance the Wilson plaquette action at the same bare coupling) …"

Theorem F.6 gives a bound proportional to `d_adm(Ω₀,Ω₁)`; to conclude equality one needs `d_adm→0`, supplied by F.10 **only** for the *same* local bare action at the *same* bare coupling. The AF trajectory (small `g_bare`) and the IR construction (strong coupling) are **not** based on the same bare action, so the use of F.10 to identify them is **unjustified**. If one reinterprets the IR base as a *blocked* version of the same weak-coupling action, proving the blocked action is KP-small **is** O2-Bridge.

### 3.6 Eq (10.72) — a non-sequitur, harmonic divergence **[V, quoted]**
> "`h_{k+1} − h_k = 2β₀ log b + ε_k`, where `ε_k = O(g_k²) + O(K_k/g_k) + O(K_k²/g_k²)`. … The previously established decay `K_k→0` together with `g_k→0` implies `Σ_{k≥k₀}|ε_k| < ∞`. Consequently, `h_k = h_{k₀} + 2β₀(log b)(k−k₀) + O(1)`, and the `O(1)` remainder has a finite `k→∞` limit."

**False as written.** AF gives `g_k ~ (2β₀ log b · k)^{-1/2}`, hence `g_k² ~ 1/(2β₀ log b · k)` and
`Σ_k |ε_k| ⊇ Σ_k c g_k² ~ (c/(2β₀ log b)) Σ_k 1/k = ∞` — **harmonic**, not summable. (`g_k→0` does *not* imply `Σg_k²<∞`; only `Σg_k^{2+δ}` converges.) The correct remainder is `O(log k)`, not `O(1)`. The leading asymptotic `g_k ~ (2β₀ log b·k)^{-1/2}` survives, so this is **repairable**, but the claimed summability and the telescoping in Thm 10.8 are **not valid as written**. (Also: `ε_k` is reused as a *different* symbol in Eq (10.78)'s interlacing defect — expositional trap.)

### 3.7 Companion Part 1 admits the gap is fixed-`a` only **[V, quoted]**
> "… confined to the strong-coupling domain `0<β<β(N)` at fixed lattice spacing `a`, and … the critical question of whether this nonzero gap survives as `β` increases (approaching the continuum limit) is **deliberately left** [unanswered here]."

This is a direct admission that the strong-coupling gap is **fixed-`a`**, and the continuum transport is the open question. Sections 5–10 claim to supply it; §3.4–§3.6 show the bridge is **assumed**.

### 3.8 FS verdict / exact missing step
The FS paper is a **well-organized conditional construction**, not an unconditional proof. Its single exact open obstacle, in the project's language:

```
sup_{a,L}  𝔄_κ( R_b^J S_{a,L} )  <  ε_KP     at ℓ_* = 1/(cΛ),  J = ⌈log_b(M/c)⌉,
```

i.e. **(`a,L`)-uniform KP-smallness of the gauge-invariant effective action at the crossover scale**, equivalently a **common contracting/KP polydisc for the UV (AF) and IR (strong-coupling) trajectories** — without assuming the gap. In FS: the common-polydisc premise of Thm 10.8 + the same-bare-action requirement of Lemma F.10 + the Eq (10.72) bookkeeping. **Label: `[C]` conditional; `[O]` the missing step; Clay 0%.** (This refines `CHATGPT-METHOD-RGENTRY.md §5`, which reached the same verdict; the quotes above are re-verified against the arXiv text this session.)

---

## 4. Deepest audit B — OpenAI 4D existence (G1-only, machine-checked)

### 4.1 What is actually established **[V]**
From the committee's own **assessment (10 Nov 2026)** and the KITP **3rd statement (9 Nov 2026)**, and the public repo:
- Continuum correlation functions as a **scaling limit of 4D Wilson lattice gauge theory** at **fixed physical volume**; UV convergence ("Cauchy in N") along a chosen spacing sequence.
- **Subsequential** infinite-volume limit by **Arzelà–Ascoli compactness**; one subsequence works for all correlators simultaneously.
- **Phase decomposition** (Gelfand + Riesz–Markov) into component theories that satisfy **clustering**; **OS reconstruction** (their own theorem, cross-checked against Douglas–Mei–Yin `OSreconstruct`) ⇒ Wightman axioms; **non-free** (trace anomaly is the "fingerprint").
- Gauge-invariant curvature polynomials, conserved stress tensor, Ward identity, an OPE with error bounds, asymptotic-freedom scaling.
- Lean: clean rebuild **Lean 4.32.0**, 9,470 modules, **250,319 declarations**, **no project axioms**, **no `sorryAx` dependency** (sole source `sorry` is the reference `Challenge.lean` stub), principal results use only `propext, Classical.choice, Quot.sound`; **Comparator** passed (30 Oct public target; 3 Nov combined target); independent OS adapter passed.

### 4.2 What is explicitly NOT proved **[V, quoted]**
> "It does **not** prove that the infinite-volume limit is unique or that there is a **mass gap**." (assessment, p.1)

> "no mass gap is claimed. Uniqueness across volume subsequences and agreement of the component theories are also unresolved." (assessment, p.1)

> "The construction provides **no evidence for or against a mass gap** in the theory." (assessment, §2, item 4)

> "A positive mass gap … are **outside its scope**." (repo README)

Open-questions list (assessment §8): (1) convergence + uniqueness of the infinite-volume limit; (2) a probability measure on fields; (3) **a uniform mass gap**; (4) Wilson-loop observables + **area law**; (5) single-particle spectrum. Also: **universality** of spacing sequences is not proved; convergence only for test functions away from coincidences; contact-term OS positivity on unrestricted tests not established.

### 4.3 Why it is NOT a gap claim, and where G2 would enter
The construction is a **G1** result: it builds *a* Wightman theory (along a subsequence / a phase) as a scaling limit. **Nothing in the chain produces a mass gap lower bound.** The gap would enter at exactly two places, neither touched:
1. **Clustering / uniqueness of the vacuum** — the subsequential limit is not shown to cluster; clustering would imply vacuum uniqueness but is not proved. No clustering ⇒ no spectral gap.
2. **An `(a,L)`-uniform lower bound** — the compactness (Arzelà–Ascoli) gives *subsequence* existence and *local boundedness/equicontinuity*, not *uniform exponential decay* and not a positive `m₁`. There is no `G″`-type or residue-type sign-controlling estimate anywhere.

So OpenAI **confirms the project's G1/G2 split**: existence is separable and now (modulo independent peer review of the 1,147-page manuscript) plausibly done; **the gap is untouched**, and the "hard core" is exactly where `PROOF-SKELETON.md` puts it (L3/G1-UV vs H1/G2-IR, same wall).

### 4.4 The validation risk (the real caveat)
The committee is explicit: *"We do not claim to fully understand the proof yet"*; the assessment is based on the **Lean code**, **not** the manuscript; they used OpenAI's "Astra" to digest the code; human review of definitions was "extensive (but not exhaustive)." Douglas's own framing **[V]**: *"A computer can check that a proof is correct, but not that the definitions written in the code mean what we think they mean."* The honest label: **the existence of this Lean-checked claim and its G1-only scope are `[V]`; the truth of the theorem is `[C]/[R]` pending independent peer review**, and in any case it is **not a gap**.

---

## 5. Formalization failure modes — definitional smuggling vs real formalization

The most instructive audits are the ones that *look* like machine-verified full solutions. They fail in a recurring way: **Lean/Coq verifies the term, not the meaning.** This is the single most important lesson of the whole sweep (cf. Douglas, §4.4).

### 5.1 Farooqui Coq — the textbook case **[V]**
Repo `Shariq81/yang-mills-mass-gap` (Coq 8.18; "1,306 Qed"). README claims *"for lattice gauge theory at any coupling constant beta > 0, there exists a positive mass gap m > 0."* But:
- The explicit theorem is `ym_explicit_mass_gap`: **`m = β/10 − 4 for β > 50`.** `β>50` is the **weak-coupling / continuum** regime; a gap **growing linearly with β** is the exact **dimensional error** the project recorded in `PROOF-SKELETON.md §0` (`m_lat=√G″_a = a·m_phys`; "fixed-`a` `G″_a≥λ` constant ⇒ `m_phys=√λ/a→∞`, a lattice/UV artifact, not a particle"). So the "gap" is a UV artifact, not a physical continuum mass.
- The README **admits** the proof is *"conditional on one research-level hypothesis: Balaban Pointwise Convergence … open for YM_4,"* plus *"Five additional textbook hypotheses."*
- The CHANGELOG claims "711 Qed theorems. 0 hypotheses" by **"encoding the Wilson bound definitionally using a Coq `Record` type"** — i.e. moving the needed physical hypothesis (`β>50 ⇒ |activity| ≤ e^{-(β/10-4)|size|}`) into the **definition** of "PhysicalPolymer." Only polymers satisfying the bound count, so the "gap" theorem is **vacuous relative to the actual YM problem**: the premise was assumed by fiat.
**Verdict:** `[R]` self-published overclaim; **definitional smuggling + wrong scaling + open Balaban hypothesis.** Does **not** establish the uniform `(L,a)` gap.

### 5.2 jonwashburn/ym-proof — does not compile **[V]**
Lean Zulip (Eric Wieser, Michael Rothgang): *"it is rejected by Lean; `lake build YangMillsProof` gives errors about things not existing"*; *"still does not compile at all"* — "looking-plausible-at-a-distance-but-actually-not-working LLM output." Claims sorry/axiom-free full Clay proof; **fails to build**.

### 5.3 "Synthesis Reasoning/ASI" (Zenodo 19432416), merchantmoh "Ark Research", hexagonmath **[R]**
"Millennium #2 SOLVED," "10/10 verification," "positive Ricci curvature forces the gap" — implausible/incoherent; **unverified overclaims**; the hexagonmath "every compact simple group β≈0.76" is a finite-lattice certified cover only (continuum open). Treat as `[R]`, not contenders.

**Lesson.** A machine-checked `Qed`/`#print axioms` line is **necessary but not sufficient**: the audit must read *what the definition smuggles* and *what scaling the statement encodes*. The three failure modes — (i) definitional smuggling, (ii) fixed-`a`/UV artifact mistaken for a physical gap, (iii) conditional-on-open-hypothesis relabelled as theorem — are precisely the "G1-only / strong-coupling / non-sequitur / H1" masks this survey is built to detect.

---

## 6. Verified-numerics / interval-arithmetic (4D gauge)

These are **real** computer-assisted results, but **finite-volume / strong-coupling**, **not** continuum gaps:
- **Box-doubling criterion `Δ₄`** (4D form of OpenAI's finite-volume criterion; formalized first by I. J. Sessford, used by Bhattacharjee et al.): measured on 2⁴–16⁴; strong-coupling anchor `Δ₄(16) ≤ 2⁻¹⁰` for `β_W ≤ 0.01273`. **[V]** (Rayn-McLean, Zenodo 24450219, self-published; explicitly "does not solve the mass gap"; neutral sector open).
- **Interval-arithmetic Hessian bound** `h_F ≤ 14.420001` (forest-gauge Hessian), with outward-rounded interval elimination (least eigenvalue `> 559/716800`) and a **Lean 4 kernel-checked certificate** (Romanov, `github.com/romanov360/yang-mills-4d-forest-certificate`). **[R]** (reported; a genuine verified-numerics gadget, but a **finite-lattice Hessian bound**, not a continuum gap).
- **Zenodo 23788087 v4.9.0** (Bhattacharjee–Bhattacharya–Lockwood): lattice transfer-matrix gap up to `β_W≤0.345` (SU(2)), `β_W≤0.76` (SU(3)); 86 Lean theorems, 0 sorry; **conditional** continuum via a **not-proved** Dobrushin–Shlosman cube estimate; paper states *"the continuum Yang–Mills existence and mass-gap problem remains open"* and *"the crossover remain[s] open."* **[V]** — an honest partial.

**Verdict:** `[R]/[T-in-domain]` — verified-numerics is maturing (interval arithmetic + Lean certificates are real and reusable), but every 4D certificate to date certifies a **finite-lattice/strong-coupling** quantity; none certifies the `(L,a)`-uniform continuum gap across the crossover.

---

## 7. The 2D O(3) case — contested, and what it does/doesn't give 4D

- **Five AI-generated preprints** (late Sep–early Oct 2026) on the nearest-neighbor O(N) model on `ℤ²`. **Problem A** (O(N), N≥3, all β, exponential decay) is proved by the **twist/commutator** trick (~1,500 lines) and **formalized in Lean 4** (117 Mathlib files, standard axioms; Douglas: *"we rebuilt and checked it"*). **Problem B** (O(3) continuum via Balaban-style RG, ~22,000 lines) is **"Plausible, unaudited, and not formalized."** **[V]** (Douglas slides).
- **Status: contested `[R]`.** Douglas calls it "(apparently) solved"; Bauerschmidt "seems very reasonable"/"horrible" on different papers; a community digest (`bauerscr/nlsm`) is in progress; reconstruction notes (e.g. Zenodo 23827755, self-described "a guess, not OpenAI's work") argue the infinite-volume limit is merely a uniform finite-box bound and the mass gap does not follow. **Not independently certified.**
- **What it gives 4D — and what it does not.** Douglas's own 4D caveat is decisive **[V]**: the sigma-model closer is a **strong-coupling** estimate; *"For YM4 the crude estimate must come from the strong coupling regime, so the RG has to flow from weak coupling to couplings of order one and hand off — a regime the sigma model papers never treat. Gauge fixing and renormalizability … are the other YM-specific difficulties."* In other words the 2D template **relocates** the wall; it does **not** cross it. Gauge fixing + infinitely many divergences + the `O(1)`-coupling hand-off are exactly the `H1`/L3 wall.

**Correction to the `[R]` context.** The prompt's "Lean formalization reportedly covers only a classical lemma" is **not accurate**: the formalized portion is Douglas's **Problem A decay theorem** — a genuine new research-level result (twist → stiffness bound → exponential decay for all N≥3, all β), formalized. The *continuum mass gap* (Problem B) is the part that is **unformalized and unaudited**. So: "a real partial formalization (2D-A) + an unformalized contested 2D-B claim," not "only a classical lemma."

---

## 8. Reusable-technique inventory (honest)

**Yes — reusable now:**
1. **Lean OS/RP infrastructure**: Douglas–Mei–Yin `xiyin137/OSreconstruct` (independently formulated OS axioms, used as the OpenAI cross-check); Hoback–Douglas OS axioms + reflection positivity + chessboard estimates + forest formula; `mrdouglasny/lgt` (lattice gauge theory, Dobrushin strong-coupling mass-gap theorems). These can formalize **future G1/G2 partials**.
2. **The 2D decomposition template**: "gap *test* is invariant under exact blocking; prove it cutoff-independent; evaluate it once." This is a clean methodological idea.
3. **The O(N) twist→decay theorem** (formalized): a real reusable 2D technique.
4. **Verified-numerics**: `Δ₄` box-doubling criterion + interval-arithmetic Hessian certificate with Lean check.

**No — not found / not reusable for the 4D gap:**
- Any **4D non-perturbative, sign-controlling, `(L,a)`-uniform lower bound** that crosses `α_s~O(1)`. This is the project's `H1`/O2-Bridge; nothing in the swept literature supplies it. The OpenAI G1 construction and the Sussman–Alpöge draft are **UV/existence**, not gap; LGT and Zenodo 23788087 are **strong-coupling**; Farooqui is **definitionally smuggled**; the 2D template explicitly **stops before** the 4D hand-off.
- **Free-massless acid test** on the reusable items: the Lean/OS infra is *structural* (not a gap tool, so not vacuous in the acid sense); the `Δ₄` criterion and forest-Hessian bound are *finite-volume* and do not claim the continuum atom; **none** is a candidate "new tool" for the atom `f₁²>0`. So none is a route to G2.

**Bottom line on Q2:** the *infrastructure* is a real gain and lowers the cost of future formal G1/G2 partials; the *mathematics of the gap* is unchanged.

---

## 9. Provenance, verification, corrections

- **Verified this session `[V]` (primary):** KITP 3rd committee statement (Nov 9 2026); committee assessment PDF (Nov 10 2026); OpenAI repo `yang-mills-group/OpenAI-YM-20261028` (README + file list); Douglas Yang Symposium slides (Oct 13 2026); Faizal–Shabir arXiv:2606.19362 full text (Thm 10.4/10.5/10.8, Lemma F.10, Eq 10.72–10.74, companion Part 1 admission); Lean Reservoir `mrdouglasny/lgt`; Zenodo 23788087 / 24450219 / 22904674 / 17083731; Farooqui repo `README.md`/`CHANGELOG.md`/tree; Lean Zulip thread (jonwashburn).
- **As-reported `[R]` (secondary / self-published / not independently fetched):** Sussman–Alpöge draft details (GitHub is WIP); Zenodo 22767474 v54; "Synthesis Reasoning/ASI" 19432416; merchantmoh; hexagonmath; Zenodo 23827755 reconstruction; parts of the 2D-O(3) community chatter.
- **Secondary press** (e.g. dataconomy "Nobel panel validates …") was **not** used as evidence; the committee's own PDFs and repo were. Note the press headline overstates: the panel validated **existence**, explicitly **not** the gap, and is "not a ruling on the Clay Prize."
- **Citation correction carried forward:** Zenodo **22840696** is **not** a Faizal–Shabir audit (it is a different G-field program, V. Liashkevych); FS-specific criticisms previously attributed to it are `[R]` and unsupported. (This corrects `RG-ENTRY-VERIFICATION.md §4`, per `CHATGPT-METHOD-RGENTRY.md §5.7`.)
- **Naming note:** first author is **Mir Faizal** (not "Faizal–Shabir" as a surname pair); co-author Arshid Shabir.
- **Open question for the field:** whether the OpenAI 1,147-page manuscript survives independent mathematical peer review. Even if it does, it is **G1-only**. **The gap remains open.**

---

## 10. One-line verdict & scorecard

> **No AI-generated or formalized 4D Yang–Mills *mass-gap* claim survives audit.** The only clean machine-checked AI result is **G1 (existence), explicitly G1-only**, and it places the remaining gap exactly where this project does — across the `α_s~O(1)` crossover (`H1`). The only serious 4D *gap* claim (Faizal–Shabir) is **human, published, and conditional**, failing at the same crossover (common-polydisc premise + parameter-identity + harmonic-sum non-sequitur). Reusable gains are **infrastructure, the 2D template, and verified-numerics** — real, but **not** a route to G2.

| Question | Answer | Label |
|---|---|---|
| Any audited-valid AI/formal **4D mass-gap** claim? | **No** | `[NG]` for the surveyed class; `[V]` |
| Any audited-valid AI/formal **4D existence (G1)** claim? | **Yes, plausibly** (OpenAI, Lean-checked, panel-reviewed) — **explicitly no gap** | `[V]` scope / `[C]/[R]` truth pending peer review |
| Most serious 4D *gap* claim? | **Faizal–Shabir 2026** | `[C]` conditional |
| Its precise missing step? | `sup_{a,L} 𝔄_κ(R_b^J S_{a,L}) < ε_KP` at the crossover = **common UV/IR polydisc = H1** | `[O]` |
| Reusable **new 4D-gap-crossing** tool found? | **No** | `[NG]`/`[O]` |
| Reusable **infrastructure / 2D / verified-numerics**? | **Yes** (OS/RP/forest Lean; 2D decomposition + twist theorem; `Δ₄`, interval Hessian) | `[T]`/`[R]` |
| Does this change Clay's status? | **No.** Gap open. | `[V]` |

*Status of this note: a verified survey + audit, not a proof and not a no-go theorem for the problem. It records what is actually formalized/proved vs asserted, and the exact missing step of the most serious 4D claim. No QED on anything unproven.*

---

### References (external, with verification status)

- **OpenAI / committee:** KITP, *Proof by OpenAI of the existence of Quantum Yang–Mills Theory — 3rd Statement by the Shepherding Committee*, Nov 9 2026 `[V]`; *Validation of the OpenAI Yang–Mills construction*, Nov 10 2026 (repo `yang-mills-group/YM`, `docs/assessment-1.pdf`) `[V]`; repo `yang-mills-group/OpenAI-YM-20261028` `[V]`; M. R. Douglas, *The Yang–Mills Millennium Prize Problem*, Yang Symposium, Stony Brook, Oct 13 2026 (slides) `[V]`; M. R. Douglas, *Nature Rev. Phys.* **8** (2026) 86–97, DOI 10.1038/s42254-025-00909-2.
- **Faizal–Shabir:** arXiv:2606.19362 `[V]`; *Fortschr. Phys.* **74** (2026) e70097, DOI 10.1002/prop.70097 `[V]`; companions *Int. J. Geom. Meth. Mod. Phys.* (2026) 2650111–2650114.
- **Lean packages:** `mrdouglasny/lgt` (Lean Reservoir) `[V]`; `xiyin137/OSreconstruct` (Douglas–Mei–Yin) `[V]` (cited in assessment); Zenodo 22767474 (Lean 4 Formalization Program v54) `[R]`; Zenodo 22904674 (Rigged-Hilbert scaffolding, Norris et al.) `[V]`.
- **Verified numerics:** Zenodo 24450219 (Rayn-McLean) `[V]`; Zenodo 23788087 v4.9.0 (Bhattacharjee–Bhattacharya–Lockwood) `[V]`; `github.com/romanov360/yang-mills-4d-forest-certificate` (interval Hessian, Lean) `[R]`.
- **Outliers / overclaims:** `Shariq81/yang-mills-mass-gap` (Coq) `[V]`; `jonwashburn/ym-proof` (Lean, fails to build) `[V]`; Zenodo 19432416 ("Synthesis Reasoning/ASI") `[R]`; Zenodo 17083731 (Ancri, AI exploration, self-caveated) `[V]`; merchantmoh/"Ark Research", hexagonmath `[R]`.
- **2D O(3):** OpenAI five preprints, status per Douglas slides `[V]`; community digest `bauerscr/nlsm` and reconstruction Zenodo 23827755 `[R]`.
- **Internal cross-refs:** `D:\Dev\math\docs\PROOF-SKELETON.md`; `D:\Dev\math\docs\BRIEF-FOR-MODEL.md`; `D:\Dev\math\research\unconditional\CHATGPT-METHOD-RGENTRY.md` §5; `D:\Dev\math\research\unconditional\SOTA-VERIFIED.md`; `D:\Dev\math\research\unconditional\TRANSPORT-FIVE-WAY-VERDICT.md`; `D:\Dev\math\docs\LEAN-STATUS.md`.
