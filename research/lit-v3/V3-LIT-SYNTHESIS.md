# V3-LIT-SYNTHESIS — 14 eras, ~443 entries: the literature contains a strong-coupling theorem and a finite-group weak-coupling theorem, but nothing for continuous SU(N) at weak coupling in 4D; no monotone connects them

> **VERDICT (no QED, Clay 0%).** Fourteen era reports, ~443 entries, ~299 bibliographically
> checked. The literature does **not** contain a proof of the 4D pure-YM mass gap, and it does not
> contain the transport machine. What it does contain is a precise **shape of the hole**, sharper
> than the campaign's own previous statement:
>
> **(LIT-GAP)** *In 4D there is a rigorous mass-gap theorem at strong coupling for continuous
> SU(N) (Osterwalder–Seiler 1978; Shen–Zhu–Zhu 2023, β < 1/48), and a rigorous weak-coupling
> theorem only for **finite** gauge groups (Adhikari–Cao), but **no theorem for continuous SU(N)
> at weak coupling**. The physical coupling sits between the two, and no known monotone quantity
> connects them.*
>
> Three attack surfaces survive the survey, ranked in §5. The rest of this file is the evidence
> map, the v2→v3 correction list, and the list of what is dead.

---

## §0 · Labels

`[T]` theorem · `[C]` conditional · `[NG]` no-go for a named class · `[O]` exact missing step ·
`[V]` numerical. `[LIT]` marks a statement about the literature map, not about Yang–Mills.

---

## §1 · The strongest result in each era, and what it does not give

| Era | Strongest rigorous statement | Regime | What it does **not** give |
|---|---|---|---|
| 01 (1954–65) | gauge invariance forbids a local `m²A²` (algebraic fact) | all | no spectral statement; "gauge invariance forbids mass" was never a theorem |
| 02 (1963–73) | renormalizability order by order; AF (`b₀>0`) | UV | no IR statement; scale generation ≠ mass gap |
| 03 (1974–79) | exact Hessian spectrum in a constant background; one-loop `V_eff` | fixed background | no true-vacuum stability, no condensate, no gap |
| 04 (1974–82) | strong-coupling gap + area law (OS 1978; SZZ 2023 at `β<1/48`) | `β<1/48` | no transport to `β~2.5–6`, no `a→0`, no `L→∞` isolation |
| 05 (1982–90) | Balaban UV stability, fixed volume | UV, weak coupling | no OS reconstruction, no volume-uniform, no IR |
| 06 (1985–95) | Mack–Petkova Peierls condition (model) | model | no gauge-invariant spectral bound; KO criterion fails via Gribov |
| 07 (1993–2000) | truncated Landau-gauge DSE scaling solution | gauge-fixed | no gauge-invariant mass; truncation error uncontrolled |
| 08 (1998–2005) | exact dualities for SUSY theories | N=4/SUSY | no pure-YM dual; hard wall not a solution; soft wall NEC-violating |
| 09 (2002–10) | continuum-controlled SF running coupling | UV/`1/L` | no infinite-volume `0⁺⁺` residue; window only perturbative |
| 10 (2008–14) | `a`-theorem; `1/3≤a/c≤31/18`; bootstrap | CFT existence | phase-blind; no bridge from `a` to `f₁²`; no flow statement |
| 11 (2013–18) | 3D stochastic quantization; `φ⁴₃` | d=3 | 4D gauge: nothing (marginal, infinite counterterms) |
| 12 (2018–21) | resurgence fixes Borel exponent/Stokes | transseries | does not fix the transseries parameter `C`; no `Δ/Λ` |
| 13 (2021–24) | SZZ `β<1/48`; CNS area law in 't Hooft regime | strong coupling / large-N | physical `λ≈1–3`; 16–48× away; no continuous weak coupling |
| 14 (2024–26) | Lean 4 OS axioms for the 4D free boson; formalized existence (external) | free/existence | mass gap explicitly excluded; free field is gapless |

**Every row's "does not give" is the same object in a different coordinate.** This is the
literature-level confirmation of the campaign's one-species claim.

---

## §2 · Three facts the whole literature agrees on

### 2.1 The coupling gap is exact and structural

- SZZ's 4D SU(N) theorem holds for `β < 1/48` (`λ > 48`). The physical SU(3) coupling is
  `λ ≈ 1–3`. The rigorous window stops **16–48× short** (Era 13).
- At weak coupling, the only 4D area-law/mass-gap theorem for continuous groups is absent; the
  available weak-coupling decay theorem is for **finite groups** (Adhikari–Cao, Era 13).
- The 2025 CNS results (area law, and a companion mass gap + unique infinite-volume limit) live in
  the **'t Hooft / large-N regime**, not at physical SU(3) (independently re-verified 2026-11-24:
  arXiv:2505.16585, arXiv:2509.04688).
- **So (LIT-GAP) is not rhetoric:** there is a theorem at strong coupling, a theorem for finite
  groups at weak coupling, and a hole in between for continuous SU(N). `[LIT]`

### 2.2 No monotone along `β`

The transport machine needs a quantity that is (i) monotone along the flow, (ii) positive at
strong coupling, (iii) a lower bound on the physical mass. The survey found none:

- `f₁²(β)` is not monotone (it is a residue; multiplicative transport gives sign, not size);
- `m_{0⁺⁺}(β)` is not monotone (no theorem; the gap can in principle close);
- `χ(β)` is not monotone;
- `a` (the anomaly coefficient) **is** monotone (Komargodski–Schwimmer) but **phase-blind**: at a
  gapped endpoint `a_IR=0` is allowed, so it is quantitatively silent on `Δ` (Era 10);
- entropy/`c`-functions are monotone but blind to the atom (Era 10, and the campaign's TOOL wave).

**The closest existing bridges** (all unproved, all candidates for v3):
1. **Condensate–residue sum rule** (Era 05): `⟨G²⟩ ≈ ∫dν`; if the atom dominates the IR end, the
   leading IR piece is `f₁²`. The renormalon ambiguity is the same order, so the bridge must be
   renormalon-subtracted.
2. **`a`–residue sum rule** (Era 10): prove `a_IR ≥ F(f₁², m₁)` or `Δa ≥ c·f₁²`; this would turn
   the exact `a`-monotonicity into a lower bound on the residue.
3. **Polymer-weight–residue dictionary** (Era 04): SZZ's gap is the spacing `m₁=4ln(1/u)/a`; P3's
   residue is the weight. They are complementary, and the same `u` controls both.

### 2.3 Every "mechanism" is a picture, a gauge-fixed computation, or a model

The survey confirms, era by era, that the only gauge-invariant rigorous statements are the
strong-coupling polymer results, the area law, and the finite-group weak-coupling theorem. Vortex
and monopole condensation, abelian projection, Kugo–Ojima, DSE scaling, Gribov–Zwanziger, hard/soft
wall holography, and FRG finite-`k` masses are **not** gauge-invariant spectral theorems. This is
not a criticism of the literature; it is the precise reason the gap is open.

---

## §3 · Era verdicts (one line each, with the correction that matters)

- **01.** The "no mass term" worry is algebraic, not spectral; Stueckelberg 1938 and the 2D
  Schwinger model 1962 already show gauge theories can be massive. The only no-go is Goldstone's,
  for global symmetries.
- **02.** Renormalizability is UV order-by-order and conditional on a gauge-invariant regulator;
  Coleman–Weinberg generates a scale but no gap. FP is globally incomplete (Gribov); BRST
  cohomology is the right physical-subspace handle.
- **03.** Savvidy/Nielsen–Olesen/Leutwyler compute the operator and the one-loop potential;
  **true-vacuum stability and the condensate are not theorems**, and "barrier-top vs vacuum" is a
  modern reinterpretation. The instanton gas is uncontrolled at large `ρ`.
- **04.** `1/48` is **SZZ 2023**, not OS 1978 (OS give a group/dimension-dependent small-`β` gap).
  The mechanism is the four-wall tube; P3 (weight) and OS/SZZ (spacing) are complementary.
- **05.** Balaban proves **UV stability at weak coupling in fixed volume** — no OS, no volume
  uniformity, no IR. The condensate has an `O(Λ⁴)` IR-renormalon ambiguity of its own size.
- **06.** No 1985–95 mechanism proves a gauge-invariant gap; Mack–Petkova is a conditional Peierls
  bound; Kugo–Ojima's premise is broken by Gribov copies; GZ has horizon geometry but no physical
  spectrum.
- **07.** The DSE scaling solution is a truncated Landau-gauge solution (gluon→0, ghost→∞,
  `α_c≈9.5`); finite-`k` ERG masses are regulator/background artifacts; scaling vs decoupling is a
  dispute about gauge-fixed correlators.
- **08.** Holography has exact duals only for SUSY theories; AdS/QCD inserts the IR scale by hand;
  hard wall is not a solution, soft wall violates the NEC or needs a tachyon; the "gap" is often
  the KK scale `1/R`.
- **09.** Step scaling is UV-only; `N_f=0` is outside the window only perturbatively; the
  `N_f=12` classification is **contested** (LSD vs Hasenfratz–Peterson); all lattice window claims
  are finite `(a,L,m)`.
- **10.** Exact logic: `¬(S) ∧ ¬(X) ⇔ gap`, where S = gapless CFT IR, X = non-conformal gapless
  phase. Every 2008–14 result constrains S only; `a` is silent on `Δ`; a bridge from `a` to `f₁²`
  does not exist.
- **11.** 3D results are weaker than often stated: CCHS 2024 is YM–**Higgs**, local-in-time, allows
  finite-time blow-up — no global measure, no OS, no gap. Pure 3D YM still has no construction.
  4D is marginal (`[g²]=0`), so regularity structures do not sum.
- **12.** The gap is `Δ=CΛ`, `Λ=μe^{−1/(2β₀g²)}`; all `g²`-derivatives vanish. Resurgence fixes
  the exponent and Stokes jump, **not** the transseries parameter `C` or the real part across the
  Savvidy saddle.
- **13.** SZZ `β<1/48`; weak coupling only finite groups; `N_f=12` contested; the 2025 CNS/Nissim
  results are 't Hooft-regime, not physical SU(3); no `N`-uniform gap constant is stated.
- **14.** Lean 4 now formalizes the OS/Glimm–Jaffe axioms for the **free** 4D boson
  (Douglas–Hoback–Mei–Nissim, arXiv:2603.15770; independently re-verified 2026-11-24). The
  external 2026 OpenAI Lean **existence** proof was committee-validated (Gross, Douglas, Jaffe,
  Klebanov, Nekrasov; independently re-verified 2026-11-24) but explicitly **excludes the mass
  gap, uniqueness of the infinite-volume limit, and the area law**; the 1,147-page proof was
  unreleased at review, the committee states it does not fully understand it, and the statement
  fidelity is not peer-reviewed.

---

## §4 · Corrections to apply in v3 (must-fix list)

1. **`1/48` attribution.** SZZ 2023, not OS 1978. OS give a small-`β` gap with a group- and
   dimension-dependent constant; SZZ give `|β| < 1/(16(d−1)) = 1/48`, uniform in `N`.
2. **Conformal-window wording.** v2 says `SU(3)` with `9≤N_f≤16` is "genuinely single-coupling
   AF gapless-conformal". Correct to: the Banks–Zaks zero exists perturbatively for
   `N_f > 8.05`; lattice evidence supports IR conformality for sufficiently large `N_f`, but
   **`N_f=12` is contested** (LSD: no IRFP; Hasenfratz–Peterson: IRFP), and all statements are
   finite-`(a,L,m)`. The use of "window is real" must be restricted to "BZ fixed points exist", not
   "the whole 9–16 range is established conformal".
3. **Savvidy historical status.** The barrier-top/vacuum re-reading is a modern reinterpretation,
   not a 1970s theorem; the 1970s results are the exact Hessian and the one-loop potential. Keep
   the campaign's removal claims labelled as the campaign's, with this history stated.
4. **DKSST venue.** JHEP **10** (2015) 171, not 06.
5. **SZZ reference.** arXiv:2204.12737 (2022), CMP **400** (2023) 805–851.
6. **Small bibliographic fixes caught by the wave:** Gribov NPB 139 (1978) 1–19; Leutwyler NPB 179
   (1981) 129; Zwanziger NPB 412 (1994) 657 for the critical-limit chronology; 't Hooft's 1972 AF
   cited via the 1998 retrospective; Ito 1985 is hierarchical Migdal–Kadanoff, **not** 4D lattice
   SU(2); Fodor 2016 not 2017; Lucini–Panero; Phys. Rev. E not PRL; Abbott 2023 not Kanwar 2020.
7. **3D status.** CCHS 2024 is YM–Higgs, local-in-time, blow-up-allowed; do not cite it as a 3D
   pure-YM construction.
8. **2026 existence formalization.** Cite as **committee-validated Lean check of an external
   announcement, not peer-reviewed, statement fidelity open, gap explicitly excluded**. Never as a
   theorem "G1 solved" in the Clay sense without those qualifiers.
9. **CNS 2025.** Area law in the 't Hooft regime; the companion mass gap + unique infinite-volume
   limit is large-N/'t Hooft, not physical SU(3).
10. **Balaban.** UV stability, fixed volume, weak coupling; no OS reconstruction, no IR.
11. **DF80 mass-gap condition.** The 2025 dynamical area-law proof verifies the Durhuus–Fröhlich
    mass-gap condition; state that this is the mechanism, not a new gap theorem at physical
    coupling.

---

## §5 · The three attack surfaces that survive (ranked)

### A. A monotone-to-residue sum rule (highest value)
Two independent eras point to the same bridge: the condensate is `∫dν` up to renormalons (Era 05),
and the `a`-anomaly is monotone but phase-blind (Era 10). A proof of either
`Δa ≥ c·f₁²` or a renormalon-subtracted `⟨G²⟩ = F(f₁², m₁)` would convert an exact monotone into a
**lower bound on the residue** — a genuinely new species, and the only route that gives transport a
*rail* rather than a ruler.
First concrete step: derive `a` as an integral over the Källén–Lehmann spectral density in the
`0⁺⁺` channel, check the free-field and `:φ²:` acid tests, and see whether the atom term is
isolated with a positive coefficient.

### B. An `N_f`-monotonicity theorem (cleanest new target)
Prove that a universal IR observable (existence of an IRFP, `γ_m`, the `0⁺⁺` mass ratio) is
monotone in `N_f`, together with a rigorous Banks–Zaks IRFP at large `N_f`. Then `N_f=0` is on the
gapped side without computing the strong-coupling β. Era 09 flags this as the cleanest new-tool
target: an inequality in `N_f`, not another β-function fit.

### C. Extend SZZ/CNS in the two directions that matter
(i) **`N`-uniformity:** SZZ's hypothesis is `N`-uniform but no explicit `N`-scaling of the gap
constant is stated; if it does not degrade exponentially, the large-`N` route opens. (ii) **The
`(a,Λ)`-uniform RG step:** Balaban's constants grow with volume; redo the polymer expansion with
`Λ`-tempered weights and a block-dependent small-field domain (Era 05), then transport the
strong-coupling residue/gap across the crossover with a uniform RG map (Era 04). This is the
closest rigorous precedent to the missing machine.

**Infrastructure, not progress on the gap:** formalize the acid test (free field has no gap) and
the conditional chain (Era 14). Useful, citable, and explicitly *not* a gap proof.

---

## §6 · What is dead (literature-verified)

- **Instanton gas/liquid as a gap derivation** — uncontrolled large-`ρ`; Witten 1979.
- **Kugo–Ojima as a gap theorem** — premise broken by Gribov copies; confinement ≠ gap.
- **DSE scaling as a physical mass** — gauge-fixed, truncated; decoupling favored post-era.
- **FRG finite-`k` masses** — regulator/background artifacts.
- **Holographic hard/soft wall as proof** — model; non-solution or NEC violation; KK scale.
- **Bootstrap as a gap proof** — tests CFT existence, not RG flow.
- **Resurgence determining `C`** — the transseries parameter is not computed in 4D.
- **`a/c` bounds as exclusion** — phase-blind; free vector sits at the inclusive edge `31/18`.
- **Continuous SU(N) weak coupling in 4D** — no theorem exists (only finite groups).
- **2026 existence formalization as gap progress** — explicitly excluded by its reviewers.

---

## §7 · The v3 framing sentence

> **The literature contains a mass-gap theorem at strong coupling for continuous SU(N) in 4D, and a
> weak-coupling theorem only for finite gauge groups in 4D, but nothing for continuous SU(N) at weak
> coupling. The physical coupling lies between them, and no monotone quantity is known that
> connects the two. The missing machine must therefore be a *transport object* across that gap, and
> the only exact monotones available (`a`) are phase-blind unless bridged to the residue `f₁²`.**

This is `[LIT]` — a statement about the map, not a theorem about Yang–Mills. It is sharper than
"the gap is hard," and it names the three attack surfaces above.

---

## §8 · Verification status of this wave

- 14 reports, **~443 entries**, **~299 bibliographically verified**, ~143 flagged `UNVERIFIED`
  (search-budget exhaustion, not guesses). Each report carries a per-entry ledger.
- Independently re-verified by the synthesizer on 2026-11-24: Douglas–Hoback–Mei–Nissim
  arXiv:2603.15770; the 2026 Gross-committee validation of the OpenAI Lean existence proof
  (gap excluded); Cao–Nissim–Sheffield arXiv:2505.16585 and 2509.04688.
- No claim in this synthesis upgrades an announcement to a theorem. `Clay = 0%`; no QED.

*Refs:* `LIT-ERA-01.md` … `LIT-ERA-14.md`; `README.md`.
