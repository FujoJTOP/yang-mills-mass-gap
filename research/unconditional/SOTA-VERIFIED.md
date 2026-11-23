# SOTA-VERIFIED — Literature check of the campaign's premises

**As of 2026-11-23.** Every item is labeled **Verified / Unverified / Corrected**. Citations give author, year, venue, arXiv where known. Claims I could not confirm are stated as such, not glossed.

---

## 1. Is the Clay problem still open? (Verified: YES)

- **Official statement (Jaffe–Witten, Clay):** "Prove that for any compact simple gauge group G, a non-trivial quantum Yang–Mills theory exists on ℝ⁴ and has a mass gap Δ > 0," with rigor at least that of Streater–Wightman / Osterwalder–Schrader. Clay still lists it **unsolved**.
- **Most recent authoritative overview:** M. R. Douglas, "The Yang–Mills Millennium problem," *Nature Reviews Physics* **8**, 86–97 (2026), DOI 10.1038/s42254-025-00909-2 (12 Jan 2026). Douglas states existence needs renormalization made rigorous "which has never been made rigorous in the needed generality," and on the gap that "no clear path is known."
- **Major Nov 2026 development (Verified via primary repo + press + Wikipedia):** OpenAI produced a **Lean 4 machine-checked construction of the *existence* half** — a Wightman-axiom QFT as a scaling limit of 4D Wilson lattice gauge theory, any simple compact group. Reviewed by a 5-member group chaired by D. Gross (Douglas, Jaffe, Klebanov, Nekrasov); assessment dated 10 Nov 2026, repo `github.com/yang-mills-group/YM`. It **explicitly does not prove** infinite-volume uniqueness, Wilson-loop area law, or the mass gap, and is "not a ruling on the Clay Prize"; not yet journal-peer-reviewed. Wikipedia (accessed today): "The mass gap problem remains unproven."

**Bottom line:** *Existence* now has a credible machine-verified claim (caveated, non-prize); ***the mass gap is still open***, which is the Clay-hard core.

## 2. Strongest rigorous results toward the problem

- **SZZ (Verified).** H. Shen, R. Zhu, X. Zhu, "A stochastic analysis approach to lattice Yang–Mills at strong coupling," *Comm. Math. Phys.* **400**(2), 805–851 (2023), arXiv:2204.12737. Lattice YM at strong coupling in any d; for SU(N) at **β < 1/(16(d−1)) = 1/48** in d = 4; unique infinite-volume measure, exponential ergodicity, LSI/Poincaré, Wilson-loop clustering ⇒ positive mass gap in infinite volume.
- **CNS 2025 (Verified, preprints).** S. Cao, R. Nissim, S. Sheffield, "Dynamical approach to area law for lattice Yang–Mills," arXiv:2509.04688 (Sep 2025); "Expanded regimes of area law," arXiv:2505.16585. Area law in the 't Hooft regime (U(N), SU(N), SO(2N)). R. Nissim, "U(N) lattice Yang–Mills in the 't Hooft regime," arXiv:2510.22788 (Oct 2025): mass gap + unique infinite-volume limit for U(N) at β < β̃(d).
- **Chatterjee (Verified).** S. Chatterjee, "Rigorous solution of SO(N) lattice gauge theory at strong coupling" (2015); Chatterjee–Jafarov, 1/N expansion / gauge–string duality (2016); survey "Yang-Mills for probabilists," arXiv:1803.01950 (2018; Springer Proc. Math. Stat. 283, 2019). **Large-N / strong-coupling only; no continuum mass gap.**
- **Balaban (Verified).** T. Bałaban, multiscale RG program, *Comm. Math. Phys.* series ~1983–1989 (RG approach I–II; large-field renormalization I–II; UV stability). Established UV stability and partial continuum control; **program incomplete, no mass gap**.
- **2024–2026:** **no rigorous continuum 4D pure YM mass gap.** (One 2026 *Fortschr. Phys.* paper — Faizal–Shabir, 74, e70097 — claims a full reflection-positive continuum construction with mass gap + confinement; it is **not community-accepted**, is written in first person, and contradicts the working-group/Clay consensus. Treat as **Unverified / disputed**, not established.)

## 3. Campaign-specific claims

| Campaign claim | Verdict |
|---|---|
| "CPS/CNS 2025 λ < 1/24 strong-coupling threshold" | **UNVERIFIED.** No source uses 1/24. Closest rigorous constants: SZZ **1/48** (d=4, SU(N)); Osterwalder–Seiler β < c_d/N². Likely a normalization/doubling error or a non-peer-reviewed source. |
| "All rigorous mass-gap results are strong-coupling / large-N only" | **TRUE / Verified** (SZZ, Chatterjee, CNS, Nissim). |
| "No rigorous continuum 4D pure YM mass gap exists" | **TRUE / Verified** per consensus (gap open), the disputed Faizal–Shabir aside. |
| "Existence + gap both open" | **PARTIALLY OUTDATED (Corrected).** As of Nov 2026 *existence* has a machine-verified claim (caveated); the *gap* is open. |

## 4. Does the campaign's framing match the field?

**Partly.** The field's live routes are: (a) constructive/RG continuum (Balaban-style, plus new stochastic-quantization/SPDE methods); (b) probabilistic strong-coupling + large-N (SZZ / CNS / Nissim); (c) tensor-RG / effective-action control (Douglas's suggested route). The mass gap is the acknowledged hard core, with "no clear path known."

The campaign's reduction to a **single variance inequality `G″(0) > 0`** is **idiosyncratic** — it is the campaign's own program, *not* a recognized reduction in the literature. The field does not presently reduce the problem to one variance inequality. So while the campaign's *conclusion* that the gap is open holds, its *internal framing* (one-inequality reduction, "existence+gap both open") does not match the external state of the art and should be revised: existence is now plausibly resolved modulo publication/acceptance, and the gap remains the sole open Clay-critical step.

---

### Corrections to the campaign premises
1. **Drop "λ < 1/24."** Use the verified SZZ threshold **β < 1/48** (d=4, SU(N)) or Osterwalder–Seiler β < c_d/N².
2. **Update "existence open."** As of 10 Nov 2026, *existence* has a Gross-panel-reviewed, Lean-machine-checked claim (mass gap/uniqueness/area law excluded; not Clay-awarded; not yet journal-peer-reviewed). The **mass gap** is the live open target.
3. **Qualify "all results strong-coupling/large-N":** true for accepted results, but note the *disputed* Faizal–Shabir continuum claim so the statement is not silently overreaching.
4. **Flag the variance-inequality reduction as internal, not field-standard.**

*Verified against: Clay/Jaffe–Witten; Douglas, Nat. Rev. Phys. 8, 86–97 (2026); Shen–Zhu–Zhu CMP 400 (2023); Cao–Nissim–Sheffield arXiv:2509.04688 / 2505.16585; Nissim arXiv:2510.22788; Chatterjee arXiv:1803.01950; Bałaban CMP 1983–1989; Yang–Mills Working Group repo `yang-mills-group/YM` (10 Nov 2026); Wikipedia "Yang–Mills existence and mass gap" (accessed 2026-11-23).*
