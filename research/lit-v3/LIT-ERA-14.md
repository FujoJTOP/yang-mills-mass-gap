# LIT-ERA-14 — 2024–2026: formalization, AI-assisted mathematics, current frontier

Prepared 2026-11-24 for paper v3. Window: January 2024 – November 2026.
Entries: 28. Bibliographic records checked: 22. No entry is upgraded from announcement to theorem.

Labels (per `research/lit-v3/README.md`, with one local extension):

- `[T]` theorem, including machine-checked theorems (scope stated per entry)
- `[C]` conditional theorem / conditional construction
- `[NG]` no-go for a named class
- `[O]` open, unestablished, or **announced-but-not-peer-reviewed** (the status tag says which)
- `[V]` numerical or computational evidence, including prover benchmark/capability results

Honesty rules applied: a preprint is not peer review; a committee validation is not peer review; a Lean-checked proof is not a checked *statement*; existence is not mass gap. Publisher pages at Nature, Wiley and Science were paywalled or blocked in this environment, so those records were checked through search metadata and are flagged accordingly.

---

## §0 · Status in one paragraph

The Clay Mathematics Institute still lists Yang–Mills existence and mass gap as **Unsolved** (C1), and no peer-reviewed proof of either clause appeared in 2024–2026. The era's real advances are four: (i) the Osterwalder–Schrader/Glimm–Jaffe axioms for the **4D free bosonic field** are now machine-checked in Lean 4 (A1), together with a physics library and a physics benchmark (A3–A8); (ii) AI provers reached olympiad-medal and large-scale collaborative-mathematics performance (B1–B6); (iii) rigorous progress continued in **3D** stochastic quantization and in lattice area-law expansions (D1–D2), plus scalar 4D RG/triviality (D4–D5); (iv) in late 2026 an OpenAI-generated proof of the **existence half** of 4D pure Yang–Mills was submitted to a five-member committee, Lean-formalized, and reported validated on 10 November 2026 — explicitly **not** the mass gap, not peer-reviewed, with the 1,147-page natural-language proof unreleased (C4–C5). Item (iv) is an announcement, not a theorem in the literature. The campaign's actual target — a cutoff-consistent, nonperturbative positive mass gap — is untouched by any verified result in this window.

---

## §1 · Formalization: statements versus proofs

Three distinctions govern this whole era.

1. **Formalizing a statement** means encoding definitions (gauge field, Wilson action, OS/Wightman axioms, "mass gap") and the claim, then checking that the encoding is faithful. This is a definitional and semantic task; Lean cannot do the faithfulness check by itself.
2. **Formalizing a proof** means every inference is checked by the kernel. This is the strong, reliable part.
3. **Existence versus mass gap** are separate clauses of the Clay problem. A construction of the measure/field theory does not give a gap; a gap proof presupposes an object to have a gap.

A fourth distinction matters for this campaign: **free field versus interacting field**. The only OS formalization completed in this window is for the Gaussian free field (A1). The Gaussian theory is massless, so it passes the axioms while having no mass gap. This is the acid-test failure in formalized form: infrastructure, not evidence.

### A1. [T — free field only] Douglas–Hoback–Mei–Nissim (2026): OS/Glimm–Jaffe axioms for the 4D free field in Lean 4

M. R. Douglas, S. Hoback, A. Mei, R. Nissim, "Formalization of QFT," arXiv:2603.15770 [hep-th; cs.LO, math-ph], submitted 16 March 2026, 35 pp. No journal reference. Companion ICML 2026 workshop poster: A. Mei, M. R. Douglas, "Formalization of QFT: Osterwalder–Schrader Axioms for the Free Field in Lean 4," 3rd AI for Math Workshop, ICML 2026, 11 July 2026.

What is done: the paper constructs the free bosonic field on 4D Euclidean spacetime in Lean 4 and proves it satisfies the **Glimm–Jaffe axioms**, a variant of the Osterwalder–Schrader axioms. The first release assumed Minlos' theorem, nuclearity of Schwartz space, and Goursat's theorem; the authors report that later releases, by their group and Lean-community contributors, either proved or avoided these, so the axioms are now established using Lean and Mathlib alone. Project reports describe a companion `OSforGFF` library of roughly 32,000 lines with OS0–OS4 and zero `sorry`s.

Scope and status: 4D Euclidean, **free (Gaussian) field only**. The paper makes no Yang–Mills or mass-gap claim. It is an arXiv preprint; the ICML item is a workshop poster. `[T]` for the formalized free-field theorem; **not** evidence for an interacting 4D theory and **not** evidence for a gap.

Why it matters to v3: it demonstrates that the OS axiom layer and the Gaussian-measure infrastructure (Minlos, nuclear spaces, distributions) are formalizable now. Any future formalization of a YM construction will need exactly this layer. It also fixes the benchmark for what "formalized OS axioms" means: a variant, for a free field, with statement-faithfulness checked by humans, not by Lean.

### A3. [V] Tooby-Smith (2024–2025): HepLean / PhysLean

J. Tooby-Smith, "HepLean: Digitalising high energy physics," arXiv:2405.08863 [hep-th], 14 May 2024, 16 pp. Journal: *Computer Physics Communications* **308** (2025) 109457, DOI 10.1016/j.cpc.2024.109457 (journal metadata from search; arXiv page not updated).

Formalizes in Lean 4: CKM/flavour matrices, local anomaly cancellation, Higgs-potential physics, with a broader library (Lorentz group, spinors, Standard Model, BSM models) developed under the names HepLean then PhysLean/PhysLib. `[V]` as infrastructure rather than a theorem about YM. Relevant because it shows the "definition and calculation" layer of particle physics is now machine-checkable, while the measure-theoretic QFT layer is only beginning (A1).

### A4. [V] Tooby-Smith (2024): index notation in Lean 4

J. Tooby-Smith, "Formalization of physics index notation in Lean 4," arXiv:2411.07667 [cs.LO; hep-ph, hep-th], 12 November 2024, 24 pp.

First formally verified implementation of physics tensor index notation in Lean 4, using a category-theoretic encoding. `[V]` infrastructure. It lowers the cost of stating tensor identities (field strengths, Ward identities) in a proof assistant.

### A5. [T] Tooby-Smith (2025): Wick's theorem in Lean 4

J. Tooby-Smith, "Digitalizing Wick's theorem," arXiv:2505.07939 [hep-th; cs.LO], 12 May 2025, 10 pp.

Formalizes Wick's theorem and its proof, including static and normal-ordered versions, in Lean 4 within PhysLean. `[T]` for perturbative QFT combinatorics. Directly useful to this campaign: the variance decomposition `G″ = S − D` and the residue-transport identities are Gaussian/Wick-level statements; having a checked Wick theorem is the right substrate for formalizing those finite-dimensional pieces.

### A6. [T] Tooby-Smith (2026): formalization finds an error in a published physics paper

J. Tooby-Smith, "Formalizing the stability of the two Higgs doublet model potential into Lean: identifying an error in the literature," arXiv:2603.08139 [hep-ph], v1 9 March 2026, v2 20 April 2026.

Formalizing the 2006 Maniatis–von Manteuffel–Nachtmann–Nagel analysis of two-Higgs-doublet-model vacuum stability exposed an error in their arguments, invalidating their main stability theorem. The author calls it the first non-trivial error in a physics paper found through formalization and asks how many physics papers would survive this scrutiny. `[T]` about the *argument*, not about the physical model: the formalization proves that the published proof does not establish its claim. Important precedent for v3: formalization can police published proofs, but it cannot by itself decide whether the encoded physical statement is the right one.

### A7. [O] Tooby-Smith (2025): perspective on interactive theorem provers in physics

J. Tooby-Smith, "A Perspective on Interactive Theorem Provers in Physics," *Advanced Science* (2025), DOI 10.1002/advs.202517294. UNVERIFIED: Wiley returned HTTP 403; title/venue/DOI from search only.

Perspective article arguing that ITPs are becoming practical for physics. `[O]`; not load-bearing. Recorded because it marks the transition of Lean-in-physics from isolated projects to a community program.

### A8. [V] Li et al. (2025): Lean4PHYS / LeanPhysBench

Y. Li, M. Liu, R. Wang, W. Ji, Z. He, R. Pan, J. Huang, T. Zhang, Y. R. Fung, "Lean4Physics: Comprehensive Reasoning Framework for College-level Physics in Lean4," arXiv:2510.26094 [cs.AI; cs.LG], 30 October 2025.

A 200-problem hand-written benchmark of college-level physics statements to be proved in Lean 4 (the individual problems are described as peer-reviewed; the paper itself is an arXiv preprint), plus PhysLib. Baselines: DeepSeek-Prover-V2-7B 16%, Claude-Sonnet-4 35%; adding PhysLib improves models by 11.75% on average. `[V]`. The low absolute numbers are the useful fact: current provers are far from research physics, let alone research QFT.

### A9. [O — announcement, not credible] Washburn (2025): claimed sorry-free YM mass-gap formalization

J. Washburn, `github.com/jonwashburn/ym-proof`, announced on the Lean Zulip "New Members" stream, 11 July 2025; roughly 9,250 lines, claimed axiom- and `sorry`-free formalization of Yang–Mills existence and mass gap.

Lean-community review reported that `lake build` fails and that the OS-to-Wightman reconstruction step is not established; the argument reduces the conclusion to four primitive constants of the author's own "Recognition Science" framework. `[O]`, UNVERIFIED, and not evidence. Recorded as the clearest 2025 example of the claim pattern that the 2026 committee process was designed to filter.

### A10. [O — withdrawn by arXiv] Jacobsen (2025): claimed constructive SU(3) existence and mass gap

D. C. Jacobsen, "A Constructive Proof of Existence and Mass Gap for Pure SU(3) Yang-Mills in Four-Dimensional Space-Time," arXiv:2506.00284 [physics.gen-ph], v1 30 May 2025; v2 13 June 2025 marked withdrawn. Withdrawal notice: "This paper has been withdrawn by arXiv Admin," reason "does not meet arXiv's research content quality standards."

Claimed a 5D orbifold regulator, convergent polymer expansion, OS reconstruction, nonperturbative BRST, and a Sturm–Liouville glueball mass. `[O]`. This is the cleanest documented case in the window of an existence-plus-gap claim that failed editorial screening. Cite it only as a cautionary data point.

### A11. [O — unrefereed claim] Brechtel (2025): claimed machine-verified SU(3) mass gap

Brechtel, "Machine-Verified Proof and High-Precision Numerical Validation of the SU(3) Yang–Mills Mass Gap," Zenodo, 2025, DOI 10.5281/zenodo.16267341. UNVERIFIED: metadata from search; record not independently inspected.

Claims a fully machine-verified Lean 4 proof of a strictly positive spectral gap plus numerical matching to lattice-QCD glueball masses. `[O]`. No arXiv or peer-reviewed version was found. Do not cite as a theorem.

### A12. [C/O — unrefereed cluster] Zenodo "constructive 4D Yang–Mills" claims (2025–2026)

A cluster of self-published claims, recorded to document the 2025–2026 announcement environment:

- Santos Borom, "Constructive Existence of Quantum Yang–Mills Theory on R⁴ with a Strictly Positive Mass Gap for Each Compact Simple Gauge Group G," Zenodo, 8 December 2025, `zenodo.org/records/17851647`. Claimed OS/BRST construction via a "Fourfold Algebra/GNS" scheme. UNVERIFIED.
- "The Eriksson Programme," 64-paper Zenodo volume (2026), assembling a conditional 4D SU(N) construction from Balaban RG, log-Sobolev/Poincaré inequalities, and polymer expansions, with Lean 4 lattice results (Wilson-loop area law, 2D SU(2) heat kernel, reflection positivity for Z_N gauge theory, transfer-operator gaps, thermodynamic limit). The volume is explicitly labelled **conditional**; 11 papers carry retraction/erratum/supersession notices, including one full retraction; the author withdraws any claim of solving the Clay problem. `[C]`. UNVERIFIED.
- "The Vacuum Has Steps," Zenodo, 14 December 2025, `zenodo.org/records/17929933`; "A Constructive Einstein–Cartan–Yang–Mills Theory with Positive Mass Gap in Four Dimensions," Zenodo, October 2025, `zenodo.org/records/17721926`. UNVERIFIED.

`[C]/[O]`. The cluster is useful only as evidence that the barrier is real and that self-publication does not move it. The Eriksson volume's explicitly conditional status and retractions are the honest part; the rest should not enter a citation base as results.

---

## §2 · AI-assisted theorem proving

None of B1–B5 proves a new theorem in analysis or QFT. They are capability results: they establish that neural provers can search Lean proof states at olympiad and benchmark scale. The first claimed research-level QFT result from an AI system is the 2026 Yang–Mills existence claim (C5), and it is an announcement.

### B1. [V — peer-reviewed] AlphaProof: olympiad-level formal reasoning

T. Hubert, R. Mehta, L. Sartran, et al. (Google DeepMind), "Olympiad-level formal mathematical reasoning with reinforcement learning," *Nature* (2025), DOI 10.1038/s41586-025-09833-y, published 12 November 2025 (publisher page blocked; metadata from Nature search result).

AlphaProof is an AlphaZero-style RL agent that learns to find Lean proofs, pretrained on code/math text, fine-tuned on about 300,000 Lean state–tactic pairs from Mathlib, with a Gemini-based auto-formalizer turning about one million natural-language problems into about 80 million formal problems. At IMO 2024 it solved 4 of 6 problems, including P6, reaching silver-medal level (about 28–29/42). Test-time RL gives 10–15% gains. Limitations: high compute, slow inference, weak on open-ended problems. `[V]`. This is the strongest peer-reviewed evidence that the Lean+AI loop works at olympiad scale; it says nothing about research analysis.

### B2. [V — preprint] DeepSeek-Prover-V2

Z. Z. Ren, Z. Shao, J. Song, H. Xin, H. Wang, W. Zhao, L. Zhang, Z. Fu, Q. Zhu, D. Yang, Z. F. Wu, Z. Gou, S. Ma, H. Tang, Y. Liu, W. Gao, D. Guo, C. Ruan, "DeepSeek-Prover-V2: Advancing Formal Mathematical Reasoning via Reinforcement Learning for Subgoal Decomposition," arXiv:2504.21801 [cs.AI], 30 April 2025 (v2 18 July 2025).

Open-source Lean 4 prover with recursive subgoal decomposition; 88.9% on MiniF2F-test, 49/658 on PutnamBench; ProverBench introduced. `[V]`, preprint.

### B3. [V — preprint] Kimina-Prover

H. Wang, M. Unsal, X. Lin, M. Baksys, J. Liu, M. Dos Santos, F. Sung, M. Vinyes, Z. Ying, Z. Zhu, et al. (40 authors), "Kimina-Prover Preview: Towards Large Formal Reasoning Models with Reinforcement Learning," arXiv:2504.11354 [cs.AI], 15 April 2025.

Qwen2.5-72B fine-tuned with large-scale RL; 80.7% on MiniF2F at pass@8192; open 1.5B and 7B distillations. `[V]`, preprint.

### B4. [V — preprint] Goedel-Prover

Y. Lin, S. Tang, B. Lyu, J. Wu, H. Lin, K. Yang, J. Li, M. Xia, D. Chen, S. Arora, C. Jin, "Goedel-Prover: A Frontier Model for Open-Source Automated Theorem Proving," arXiv:2502.07640 [cs.LG; cs.AI], 11 February 2025 (v3 19 April 2025).

Open data pipeline (1.64M formal statements; 800K+ with proofs), expert iteration, no RL for the SFT model: 57.6% MiniF2F pass@32, over 60% with RL/DPO, 7 PutnamBench problems. `[V]`, preprint.

### B5. [V — preprint] Goedel-Prover-V2

Y. Lin, S. Tang, B. Lyu, Z. Yang, J.-H. Chung, H. Zhao, L. Jiang, Y. Geng, J. Ge, J. Sun, J. Wu, J. Gesi, X. Lu, D. Acuna, K. Yang, H. Lin, Y. Choi, D. Chen, S. Arora, C. Jin, "Goedel-Prover-V2: Scaling Formal Theorem Proving with Scaffolded Data Synthesis and Self-Correction," arXiv:2508.03613 [cs.LG; cs.AI], 5 August 2025.

Scaffolded data synthesis, verifier-guided self-correction, model averaging. 8B: 84.6% MiniF2F pass@32; 32B: 88.1% standard and 90.4% with self-correction; 86 PutnamBench problems at pass@184. `[V]`, preprint.

### B6. [T — finite combinatorics] The Equational Theories Project

M. Bolan, J. Breitner, J. Brox, N. Carlini, M. Carneiro, F. van Doorn, M. Dvorak, A. Goens, A. Hill, H. Husum, H. Ibarra Mejia, Z. A. Kocsis, B. Le Floch, A. Livne Bar-on, L. Luccioli, D. McNeil, A. Meiburg, P. Monticone, P. P. Nielsen, E. O. Osazuwa, G. Paolini, M. Petracci, B. Reinke, D. Renshaw, M. Rossel, C. Roux, J. Scanvic, S. Srinivas, A. R. Tadipatri, T. Tao, V. Tsyrklevich, F. Vaquerizo-Villar, D. Weber, F. Zheng, "The Equational Theories Project: Advancing Collaborative Mathematical Research at Scale," arXiv:2512.07087 [math.RA; cs.LO], 8 December 2025 (v2 16 December 2025), 74 pp.

Determined all 22,028,942 edges of the implication graph between the 4,694 simplest equational laws on magmas; all proofs, human and automated, validated in Lean. `[T]` for the finite implication graph. This is the best available model for a mass-scale, machine-checked collaboration, not a result about QFT. It shows what "AI-assisted mathematics" looks like when the object is finite and the statement layer is unambiguous — precisely the conditions absent in the YM mass-gap problem.

### B7. [V — unrefereed audit] ProofGate

"ProofGate: A Reproducible Audit of Faithfulness, Alignment, and Vacuity in State-of-the-Art Lean Theorem Provers," ICML 2026 listing, `icml.cc/virtual/2026/82433`. Authors not retrieved; UNVERIFIED.

Reported kernel-level audit of 635 DeepSeek-Prover-V2 and Kimina-Prover-72B miniF2F proofs: 100% kernel-faithful under a trust set excluding `sorry`/extra axioms; about 15.2% statement-misalignment on DeepSeek miniF2F-test; Goedel-Prover-V2's released "proofs" for miniF2F alleged to be input statements with `sorry` placeholders. `[V]` as a caution, not a theorem. The methodological point is solid regardless of the audit's own review status: benchmark solve rates do not measure statement fidelity or vacuity.

---

## §3 · Clay status, surveys, and the November 2026 existence claim

### C1. [O] The official Clay statement

A. Jaffe, E. Witten, "Quantum Yang–Mills Theory," official problem description, Clay Mathematics Institute, 2000. Problem page: "Yang-Mills & the Mass Gap," Clay Mathematics Institute, `claymath.org/millennium/yang-mills-the-maths-gap/`, accessed 24 November 2026.

The page states the requirement: prove that the quantum Yang–Mills theory exists **and** that it has a mass gap; the prize is for a proof, not a disproof. Status label on the page: **Unsolved**. The page says experiment and simulations suggest a mass gap, but "no proof of this property is known," and that fundamental new ideas in both physics and mathematics are needed. `[O]`. This is the canonical statement against which any claimed formalization must be compared.

### C2. [O — peer-reviewed survey] Douglas (2026): The Yang–Mills Millennium problem

M. R. Douglas, "The Yang–Mills Millennium problem," *Nature Reviews Physics* **8**, 86–97 (2026); published online 2025; DOI 10.1038/s42254-025-00909-2. VERIFIED via search metadata; publisher page paywalled/redirected.

A peer-reviewed review of 25 years of progress from physics and mathematics. It reports that recent mathematical approaches center on rigorous stochastic quantization and rigorous strong-coupling expansions (probability theory), and that computational renormalization-group work by physicists would help. Its conclusion is that no clear path to the mass gap is yet known. `[O]`. This is the authoritative survey to cite for the state of the field at the start of the era's final year.

### C3. [O] Chatterjee (2025): Clay/Harvard Millennium Prize lecture

S. Chatterjee, "Yang–Mills and the foundations of quantum field theory," Millennium Prize Problems Lecture, Harvard University, 15 October 2025, 5:00–6:00 pm, Science Center Hall D; organizers M. Bridson (CMI), D. Freed, M. Hopkins.

Abstract: Yang–Mills theories underlie the Standard Model but lack a rigorous mathematical basis; perturbative calculations are thought to miss phenomena; constructive field theory, launched in the 1960s, has not reached its original goal; finishing that program is the Clay problem. The talk surveys "exciting recent progress" from the past decade that has "rejuvenated the quest." `[O]`. A program survey, not a solution; useful for the claim that rigorous stochastic quantization and strong-coupling probability are the live mathematical routes.

### C4. [O — primary document] The 2026 Shepherding Committee, second statement

Yang–Mills Working Group / Shepherding Committee: M. R. Douglas, D. J. Gross (Chair), A. Jaffe, I. R. Klebanov, N. Nekrasov, "Proof by OpenAI of the existence of Quantum Yang-Mills Theory — 2nd Statement by the Shepherding Committee," Kavli Institute for Theoretical Physics, 3 November 2026, `kitp.ucsb.edu/sites/default/files/2nd-Statement-About-Quantum-YangMills-Theory_Nov-3.pdf`.

Verbatim from the statement: "Currently, we are focused on checking the Lean proof of the existence of pure Yang-Mills theory, specifically that the definitions and problem statements correspond with those that are standard in physics and mathematics, and that the formalization aligns with the informal treatment. This effort is proceeding well. We hope to complete it very soon. Once this validation is done, we will release the Lean proof, hopefully next week." It also says the committee will announce any release of LLM-produced material a full day in advance.

`[O]`. The primary document establishes the existence of the process and names the committee; it does **not** itself assert that the proof is correct. As of 3 November 2026, statement-fidelity checking was still in progress.

### C5. [O — announcement; existence only] The OpenAI Yang–Mills existence claim and the 10 November 2026 validation report

Reported facts (Dataconomy, 12 November 2026; MaddMaths!/SIMAI, 18 November 2026; Science news, "AI rivals tackle one of physics' deepest mathematical challenges"; India Today, 11 November 2026; committee statement C4):

- OpenAI says an internal model produced a written proof of the **existence** of quantum Yang–Mills theory on 9 September 2026; a QCD proof followed on 18 September.
- Lean formalization was completed on 16 October 2026; the Lean code was submitted on 28 October 2026; OpenAI asked David Gross to take custody of the proof and manage release.
- On 10 November 2026 the committee released the Lean code and a 35-page validation report. The report is said to conclude that the code proves existence of pure Yang–Mills theory as a scaling limit of four-dimensional Wilson lattice gauge theory satisfying the Wightman axioms, for any simple compact gauge group.
- Explicitly **not** established, per the report as reported: the mass gap; uniqueness of the infinite-volume limit; the Wilson-loop area law. The group says the finding "is not a ruling on the Clay Mathematics Institute's $1 million Millennium Prize."
- The 1,147-page natural-language proof has **not** been released. The committee did not rely on OpenAI's written manuscript. It used OpenAI's "Astra" model to inspect the code and states it "did not claim to fully understand the proof yet." A warnings list appears on page 13 of the report; Matt von Hippel is reported skeptical of at least one required element.
- The code reportedly passed a clean rebuild, an axiom audit, and a Comparator check. An outside project, Humanise (associated with Martin Hairer), is importing the code and mapping dependencies; it had not completed an independent full compilation, axiom check, or line-by-line verification at the time of the reports.
- Status: **not peer-reviewed**. The committee's own framing is that Lean checks logical steps but cannot confirm that the encoded definitions reflect the intended physics; that is why the committee audited definitions by hand.

`[O]`. If the formalization faithfully encodes the standard statement — the committee's stated task — this would be a proof of the **first clause** of the Clay problem (existence) in 4D pure YM. It would not be the mass gap, and it is not yet a theorem in the literature. The correct citation phrase for v3 is "announced and committee-validated Lean formalization of existence, November 2026; not peer-reviewed; mass gap not addressed."

### C6. [O — press only] The rival AI-assisted effort

E. Sussman and L. Alpöge (Anthropic) are reported to have released unfinished, independently developed AI-assisted work on the same problem on 6 October 2026 (MaddMaths!/SIMAI, 18 November 2026). No paper, preprint, or formal artifact was retrieved. UNVERIFIED. Recorded because it corroborates that the 2026 activity is a race, not a single result, and because it reinforces the need for a standing verification protocol.

---

## §4 · Rigorous 2024–2026 progress on construction and the lattice

No peer-reviewed construction of 4D Yang–Mills and no mass-gap proof appeared in the window. The rigorous progress is 3D, lattice, and scalar.

### D1. [T — 3D, local] Chandra–Chevyrev–Hairer–Shen (2024): stochastic quantization of Yang–Mills–Higgs in 3D

A. Chandra, I. Chevyrev, M. Hairer, H. Shen, "Stochastic quantisation of Yang–Mills–Higgs in 3D," *Inventiones Mathematicae* **237** (2024) 541–696, DOI 10.1007/s00222-024-01264-2; arXiv:2201.03487.

Constructs a nonlinear metric space of distributions as initial conditions for the deterministic and stochastic YMH flow, a quotient by gauge orbits, local-in-time solutions of the renormalized stochastic YMH flow via regularity structures, and a canonical Markov process on the gauge-orbit quotient modulo a possible finite-time blow-up. `[T]` for 3D local theory. Not a global 3D construction, not 4D, and not a mass gap. It is the strongest peer-reviewed rigorous gauge-theory construction in the window.

### D2. [T — lattice] Cao–Nissim–Sheffield (2025): expanded area-law regimes for lattice Yang–Mills

S. Cao, R. Nissim, S. Sheffield, "Expanded regimes of area law for lattice Yang–Mills theories," arXiv:2505.16585 [math.PR; math-ph], v1 22 May 2025, v3 26 September 2025, 35 pp.

Expands the set of parameters for which the area law is proved for pure U(N) lattice Yang–Mills, especially at large N, beyond the 1978 Osterwalder–Seiler result. The method treats the master loop equation as a linear inhomogeneous equation for Wilson-string expectations, proves an a priori bound on solutions, and controls the merger term by a truncated model plus an approximation theorem. `[T]` for the lattice area law in the stated regime. Not a continuum limit and not a mass gap. It is the main 2025 rigorous lattice advance relevant to confinement.

### D4. [T — mean-field scope] Wang–Kopper (2024–2025): triviality for mean-field φ⁴₄

P. Wang, C. Kopper, "Triviality proof for mean-field φ⁴₄-theories," arXiv:2407.01309 [math-ph], v1 1 July 2024, v3 26 February 2025; *Journal of Mathematical Physics* **66**, 062302 (2025) per search metadata (arXiv page not updated).

Within the mean-field approximation to the Wilson RG equations, proves triviality of one-component φ⁴₄ for arbitrary (including infinite) UV-regularized bare coupling, and extends the result to O(N) vector models, replacing the unphysical IR cutoff by a physical mass. `[T]` for the mean-field model, not the full 4D φ⁴ theory, not gauge theory, and not a mass gap. Useful as a benchmark for how far rigorous 4D RG has reached: even scalar 4D triviality is proved only in a mean-field truncation here, while the Bauerschmidt–Brydges–Slade program supplies the non-mean-field method being extended elsewhere.

### D5. [T — scalar RG] Park (2025): an RG map for |φ|⁴ in d ≥ 4

J. Park, "A Renormalisation Group Map for Short- and Long-ranged Weakly Coupled |φ|⁴ Models in d ≥ 4 at and Above the Critical Point," arXiv:2511.03495 [math.PR; math-ph], 5 November 2025, 94 pp.

Constructs and analyzes an RG map for the weakly coupled n-component |φ|⁴ model on a periodic box in d ≥ 4, for short- and long-range interactions with upper critical dimension four. Extends and refines the Bauerschmidt–Brydges–Slade map; aims at exact correlation decay and gives evidence for a finite-volume plateau when combined with Michta–Park–Slade finite-size scaling. `[T]` for the RG map. Scalar, not gauge; infrastructure for rigorous 4D QFT, not a YM mass gap.

### N1 (derived; not a bibliographic entry). [NG — structural] The free-field route cannot give the gap

Derived no-go, not a bibliographic entry: the Gaussian free field in 4D is massless; its spectrum is continuous. A formal proof that it satisfies OS/Glimm–Jaffe (A1) therefore cannot imply a positive mass gap, and any argument that proves the gap only for a Gaussian/free model fails the campaign's acid test (`f₁² = 0`, `G″ = 0` for the free field). Label `[NG]` for the class "free-field or Gaussian formalization used as evidence for the YM gap." This is the formalized version of the campaign's standing acid test.

---

## §5 · What formalization can and cannot do for this problem

**Can do:**

- Check every inference of a claimed proof, including long technical estimates, once the statement is encoded (A1, B6).
- Catch real errors in published arguments (A6).
- Make reductions and dependencies explicit: a formalized chain from assumptions to the Clay statement would show exactly which step is missing.
- Certify finite-lattice or finite-volume statements, where the statement layer is unambiguous (B6; the Eriksson volume's lattice fragments, though unrefereed).
- Audit statement fidelity mechanically at the level of definitional equality, while leaving the semantic question to humans (C4).

**Cannot do:**

- Certify that the encoded definitions are the intended physical/mathematical ones. This is the committee's stated central task in C4 and the reason a Lean proof of "existence" is not automatically a proof of the Clay existence clause.
- Discover the missing nonperturbative idea. Formalization is downstream of construction; it cannot supply the machine that transports the strong-coupling gap across the Savvidy crossing.
- Replace peer review or community understanding. The 2026 committee explicitly says it does not yet fully understand the proof (C5); a kernel check is not understanding.
- Handle physics semantics: A6 proves an argument invalid, not a physical model wrong.
- Currently state the full Clay problem without new infrastructure: Mathlib before A1 lacked infinite-dimensional Gaussian measure, Minlos, nuclear spaces, and OS reconstruction, and it has no Yang–Mills theory to state the problem against. That absence, not any single missing lemma, is the baseline.

**Existence versus mass gap:** the 2026 claim is the first clause only. Even a complete, peer-reviewed existence proof would leave the campaign's target untouched. Conversely, the free-field OS theorem is an existence proof for a gapless theory; it is the precise formal counterpart of the acid-test failure.

---

## §6 · Community open-problem lists currently in use

1. **Clay official list (C1).** Existence plus mass gap; proof, not disproof; publication and community-acceptance requirements. The canonical statement.
2. **Douglas review (C2).** The peer-reviewed map of approaches: rigorous stochastic quantization, rigorous strong-coupling expansions, computational RG; conclusion that no clear path to the gap is known. Use this as the survey-level problem list.
3. **2026 Shepherding Committee's remaining list (C4–C5).** After the claimed existence proof: (a) uniqueness of the infinite-volume limit; (b) a uniform mass gap; (c) the Wilson-loop area law. This is now the de facto short list for the existence program.
4. **Lean-QFT roadmap (A1 and Douglas's 2025–2026 talks, reported).** OS reconstruction theorem; P(φ)₂ via lattice construction; Yang–Mills at strong coupling; Minlos/nuclear-space infrastructure. The talk-level items are reported, not all documented in A1; treat the roadmap as a program statement, not a published list.
5. **ETP-style collaboration model (B6).** Not a problem list, but the current model for a mass-scale verified attack: decompose into finite, unambiguous claims; formalize each; publish the graph. The analogue for YM would be a formalized reduction graph from axioms to the gap statement.
6. **The campaign's own list (CLAUDE.md, current frontier).** Savvidy concavity; cutoff-consistent `G″(0) > 0`; the missing nonperturbative transport machine. This list is consistent with C2's conclusion and is not contradicted by any verified 2024–2026 result.

---

## §7 · Position in 14 eras

One line per era, from the formalization/AI angle.

| # | Window / theme | Position seen from era 14 |
|---|---|---|
| 01 | 1954–1965, origin | Problem posed in physics prose; no axioms, hence nothing to formalize. |
| 02 | 1963–1973, quantization/renormalization/AF | Perturbative series and renormalization; not proof-grade, not formalized. |
| 03 | 1974–1979, instantons/Savvidy/confinement pictures | Classical solutions are formalizable in principle; the unstable-vacuum argument remains prose. |
| 04 | 1974–1982, lattice gauge theory/strong coupling | The one layer that is formalizable today: finite lattice, Wilson action, strong-coupling expansions. |
| 05 | 1982–1990, constructive RG/gluon condensate | Balaban's programme is the technical ancestor of the 2026 existence claim; it is the object a formalization would have to encode. |
| 06 | 1985–1995, confinement mechanisms/lattice phenomenology | Area law is a formalizable statement (D2); confinement mechanisms are not proofs. |
| 07 | 1993–2000, DSE/exact RG | Truncation schemes; not proof-grade and not formalized. |
| 08 | 1998–2005, holography/large-N/AdS-QCD | Conjecture-level; outside the current formalization stack. |
| 09 | 2002–2010, lattice precision/step scaling/conformal window | Numerical `[V]`; now benchmark data for formal statements, not proofs. |
| 10 | 2008–2014, CFT constraints/a-theorem/bootstrap | Analytic theorems (e.g. a-theorem) are formalizable in principle; none has been formalized. |
| 11 | 2013–2018, rigorous probabilistic/constructive QFT | OS axioms and 2D/3D constructions become the natural Lean targets; A1 is the first 4D fragment. |
| 12 | 2018–2021, resurgence/transseries | Asymptotic analysis; not formalized, and not a route to the gap. |
| 13 | 2021–2024, rigorous frontier/large-scale lattice | 3D stochastic quantization (D1) and lattice area-law expansions (D2) are the immediate mathematical input to era 14. |
| 14 | 2024–2026, formalization/AI/current frontier | OS axioms for the 4D free field in Lean (A1); AI provers at scale (B1–B6); claimed AI existence proof with committee-validated Lean (C4–C5); mass gap still open. |

---

## §8 · Gaps and attack surfaces for v3

- **Statement fidelity is the weak link.** The 2026 committee's own central task is checking that the Lean definitions match the standard physics/math statements. v3 should treat any formalization claim as two claims: (i) the encoding is faithful; (ii) the proof checks. Only (ii) is mechanical. Attack surface: build a public, versioned "Clay statement" Lean file with explicit axioms and a documented translation, and require every claimed proof to import it.
- **Existence is not the mass gap.** The 2026 claim is existence only; it explicitly leaves the mass gap, uniqueness of the infinite-volume limit, and the area law open. v3 must not let a validated existence formalization be reported as progress on the gap. Attack surface: state the gap clause separately in Lean, with the Hamiltonian/spectral definitions and the continuum limit, and keep it as the headline `sorry`.
- **Free-field formalization is infrastructure, not evidence.** A1 passes OS/Glimm–Jaffe for a massless Gaussian theory. That is the acid-test failure in machine-checked form. Any v3 argument that reduces to Gaussian/Wick manipulations must be flagged `[NG]` for the gap (N1). Attack surface: formalize the campaign's acid test itself — a Lean theorem that the free field has no gap — so that any proposed tool can be tested against it.
- **Mathlib gaps are tractable and worth filling.** Infinite-dimensional Gaussian measure, Minlos, nuclear spaces, distributions, OS reconstruction. A1 shows the layer is within reach. Filling it is useful regardless of the gap; it is not the gap. Attack surface: contribute the OS reconstruction theorem and P(φ)₂ construction as reusable Mathlib/PhysLib modules.
- **The missing machine is not a formalization problem.** The campaign's key object — a cutoff-consistent nonperturbative transport across the Savvidy crossing — has not been built by anyone in this window. Formalization can certify it once found; it cannot find it. Attack surface: keep the "造机器" line separate from the "formalize" line, and do not let formalization activity substitute for the missing construction.
- **AI provers are benchmark-strong, research-weak.** B1–B5 are olympiad/benchmark results; B7 shows solve rates overstate statement fidelity and can hide vacuity. No peer-reviewed new theorem in analysis/QFT has come from these systems. Attack surface: use AI as a generator of candidate inequalities and counterexamples (acid tests), not as an oracle; require Lean build, axiom audit, statement diff, and independent rebuild for any AI-produced claim.
- **Announcement hygiene is now a research-infrastructure problem.** The 2025–2026 record includes a withdrawn arXiv claim (A10), a failing Lean build (A9), multiple Zenodo claims (A11–A12), and the 2026 committee process (C4–C5). The committee process is a template, not a substitute for peer review. Attack surface: adopt a standing protocol — clean rebuild, `#print axioms`, no `sorry`, statement diff against the canonical file, independent second rebuild — and apply it to the campaign's own artifacts.
- **Formalize the campaign's conditional chain.** P3 strong-coupling residue theorem, the strong-coupling gap, OS reconstruction, and the continuum limit can be stated and partially proved in Lean now. A formalized reduction graph would expose the exact missing step and make the conditional theorem citable as a machine-checkable artifact. Attack surface: start with the finite-lattice/strong-coupling layer, where statements are unambiguous.
- **Open-problem list drift.** The community now has several lists (C1, C2, C4–C5, the Lean-QFT roadmap). v3 should adopt a single tracked list: Clay mass gap; unique infinite-volume limit; uniform gap; area law; and the campaign's transport machine. Attack surface: maintain it in `research/lit-v3/` with status tags and re-verify each item per revision.

---

## §9 · Verification ledger

Status definitions: **VERIFIED** = bibliographic record checked at the cited source (arXiv abstract page, DOI landing page, Clay page, KITP PDF, or conference page). **VERIFIED (event)** = the event/document exists and is quoted, but the mathematical claim is not verified. **UNVERIFIED** = metadata from secondary search only; do not cite as established.

| # | Short citation | Source actually checked | Status |
|---|---|---|---|
| A1 | Douglas–Hoback–Mei–Nissim, arXiv:2603.15770 (2026) | arXiv abs page; ICML 2026 poster page | VERIFIED (preprint) |
| A3 | Tooby-Smith, HepLean, arXiv:2405.08863; CPC 308 (2025) 109457 | arXiv abs page; journal metadata via search | VERIFIED (arXiv); journal metadata search-only |
| A4 | Tooby-Smith, arXiv:2411.07667 (2024) | arXiv abs page | VERIFIED (preprint) |
| A5 | Tooby-Smith, arXiv:2505.07939 (2025) | arXiv abs page | VERIFIED (preprint) |
| A6 | Tooby-Smith, arXiv:2603.08139 (2026) | arXiv abs page | VERIFIED (preprint) |
| A7 | Tooby-Smith, Advanced Science (2025), DOI 10.1002/advs.202517294 | Search metadata only; Wiley 403 | UNVERIFIED |
| A8 | Li et al., arXiv:2510.26094 (2025) | arXiv abs page | VERIFIED (preprint) |
| A9 | Washburn, ym-proof / Lean Zulip (2025) | Search summary of Zulip thread | UNVERIFIED |
| A10 | Jacobsen, arXiv:2506.00284 (withdrawn 2025) | arXiv abs page, withdrawal notice | VERIFIED (withdrawn) |
| A11 | Brechtel, Zenodo DOI 10.5281/zenodo.16267341 (2025) | Search metadata only | UNVERIFIED |
| A12 | Borom; Eriksson volume; Zenodo cluster (2025–2026) | Search metadata only | UNVERIFIED |
| B1 | Hubert–Mehta–Sartran et al., Nature (2025), DOI 10.1038/s41586-025-09833-y | Nature search result; publisher page blocked | VERIFIED (metadata); publisher blocked |
| B2 | Ren et al., arXiv:2504.21801 (2025) | arXiv abs page | VERIFIED (preprint) |
| B3 | Wang et al., arXiv:2504.11354 (2025) | arXiv abs page | VERIFIED (preprint) |
| B4 | Lin et al., arXiv:2502.07640 (2025) | arXiv abs page | VERIFIED (preprint) |
| B5 | Lin et al., arXiv:2508.03613 (2025) | arXiv abs page | VERIFIED (preprint) |
| B6 | Bolan–Breitner–…–Tao–Zheng, arXiv:2512.07087 (2025) | arXiv abs page | VERIFIED (preprint) |
| B7 | ProofGate, ICML 2026 listing | Conference page metadata via search | UNVERIFIED |
| C1 | Jaffe–Witten, Clay official description; CMI problem page | Clay problem page | VERIFIED |
| C2 | Douglas, Nature Rev. Phys. 8, 86–97 (2026), DOI 10.1038/s42254-025-00909-2 | Search metadata; publisher paywalled | VERIFIED (metadata); publisher blocked |
| C3 | Chatterjee, Harvard Clay lecture, 15 Oct 2025 | Harvard event page | VERIFIED (event) |
| C4 | Shepherding Committee 2nd Statement, 3 Nov 2026 | KITP PDF, text extracted | VERIFIED (document) |
| C5 | OpenAI existence claim + 10 Nov 2026 validation report | Committee statement; Dataconomy; MaddMaths!; Science title | VERIFIED (event/announcement); report UNVERIFIED |
| C6 | Sussman–Alpöge, reported 6 Oct 2026 | MaddMaths! report only | UNVERIFIED |
| D1 | Chandra–Chevyrev–Hairer–Shen, Invent. Math. 237 (2024) 541–696 | arXiv abs page (journal ref and DOI listed) | VERIFIED |
| D2 | Cao–Nissim–Sheffield, arXiv:2505.16585 (2025) | arXiv abs page | VERIFIED (preprint) |
| D4 | Wang–Kopper, arXiv:2407.01309; J. Math. Phys. 66 (2025) 062302 | arXiv abs page; journal metadata via search | VERIFIED (arXiv); journal metadata search-only |
| D5 | Park, arXiv:2511.03495 (2025) | arXiv abs page | VERIFIED (preprint) |

Count: 28 entries; 22 with bibliographic records checked; 6 unverified (A7, A9, A11, A12, B7, C6), with C5's claim flagged separately as an announcement.

**Do not cite as a theorem:** A9, A10 (withdrawn), A11, A12 (except the explicitly conditional/retracted status), C5 (existence only, not peer-reviewed), C6.

**Do not cite as peer-reviewed:** A1, A3–A6, A8, B2–B7, C4–C5, D2, D4 (arXiv version), D5.

**Peer-reviewed in this window:** B1 (Nature), C2 (Nature Reviews Physics), D1 (Inventiones), D4 (J. Math. Phys. per search), A3 (CPC per search). No peer-reviewed 4D YM construction or mass-gap proof.

---

## §10 · Sources consulted

Primary/records:

- https://arxiv.org/abs/2603.15770 — Formalization of QFT
- https://icml.cc/virtual/2026/82544 — OS axioms for the free field (ICML 2026 poster)
- https://arxiv.org/abs/2405.08863 — HepLean
- https://doi.org/10.1016/j.cpc.2024.109457 — HepLean journal DOI
- https://arxiv.org/abs/2411.07667 — index notation
- https://arxiv.org/abs/2505.07939 — Wick's theorem
- https://arxiv.org/abs/2603.08139 — 2HDM formalization error
- https://doi.org/10.1002/advs.202517294 — Perspective on ITPs in physics
- https://arxiv.org/abs/2510.26094 — Lean4Physics
- https://arxiv.org/abs/2506.00284 — withdrawn SU(3) claim
- https://doi.org/10.5281/zenodo.16267341 — Brechtel claim
- https://zenodo.org/records/17851647 — Borom claim
- https://zenodo.org/records/17929933 — The Vacuum Has Steps
- https://zenodo.org/records/17721926 — Einstein–Cartan–YM claim
- https://www.nature.com/articles/s41586-025-09833-y — AlphaProof
- https://arxiv.org/abs/2504.21801 — DeepSeek-Prover-V2
- https://arxiv.org/abs/2504.11354 — Kimina-Prover
- https://arxiv.org/abs/2502.07640 — Goedel-Prover
- https://arxiv.org/abs/2508.03613 — Goedel-Prover-V2
- https://arxiv.org/abs/2512.07087 — Equational Theories Project
- https://icml.cc/virtual/2026/82433 — ProofGate
- https://www.claymath.org/millennium/yang-mills-the-maths-gap/ — Clay problem page
- https://www.nature.com/articles/s42254-025-00909-2 — Douglas review
- https://www.math.harvard.edu/event/millennium-prize-problems-lecture-sourav-chatterjee-yang-mills-existence-and-mass-gap/ — Chatterjee lecture
- https://www.kitp.ucsb.edu/sites/default/files/2nd-Statement-About-Quantum-YangMills-Theory_Nov-3.pdf — committee statement
- https://dataconomy.com/2026/11/12/nobel-panel-validates-openais-yang-mills-existence-proof/ — press report
- https://maddmaths.simai.eu/news-2/sembra-che-lai-abbia-dimostrato-un-altro-problema-del-millennio-o-meglio-mezzo/ — press report
- https://www.science.org/content/article/ai-rivals-tackle-one-physics-deepest-mathematical-challenges — Science news (403; title only)
- https://arxiv.org/abs/2201.03487 — 3D YMH
- https://arxiv.org/abs/2505.16585 — lattice area law
- https://arxiv.org/abs/2407.01309 — mean-field φ⁴ triviality
- https://arxiv.org/abs/2511.03495 — RG map d ≥ 4

Search-only (no stable URL retained): Zulip thread on Washburn's claim; the Eriksson Programme volume landing page; India Today report; Spektrum report.
