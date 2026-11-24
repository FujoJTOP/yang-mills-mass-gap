# SIGNED-FOURIER-BLIND — independent blind red-team on the "representation-ring / signed-Fourier (Bochner-type)" species for `f₁²`

> 2026-11-23 · **判词：SPECIES DEAD. All four kill targets land; survival score 1/5.** The proposed tool claims a *lower* bound on the atom `f₁²=|⟨0|tr F²|G⟩|²` from the character/Weingarten **sign structure**, asserted to be *not* a measure-positivity statement and therefore able to certify an atom (coherence) rather than disorder. It dies on the project's own sharpest line: **`f₁²` is a modulus-square, and a modulus-square annihilates every phase** — the Weingarten/character "sign" is a phase in the amplitude `A₁=⟨0|trF²|G⟩`, and `|A₁|²` is blind to it (`docs/yang-mills-frontier.md` line 68: "符号交替翻不动平方交叠；Weingarten 符号管的是面积律求和，触不到定 `m₁` 的逐片 `c_R`"). Making the claim precise forces it onto one of two dead rails: (i) the **positive** character expansion `c_r(β)≥0`, which **is** Kotecký–Preiss / strong coupling = P3, only at `β<1/48`; or (ii) **Bochner positivity** of the Källén–Lehmann measure `ν≥0`, which **is** OS reflection positivity (the dead family), and which is sign-blind to atoms and gives **upper** bounds. The only place a sign genuinely lives — the signed area-law **surface** sum — cannot reach the per-slice coefficient that sets `m₁`. There is no "third object". Kills: **COLLAPSE ✓, FREE-FIELD LEAK ✓, CROSSOVER ⟂ ✓, CIRCULARITY ✓.** The exact stuck step: *a signed sum that lower-bounds the modulus-square of a non-vanishing per-slice character coefficient* — which is the gap.

---

## §0. The tool, stated precisely (kill the real thing, not a strawman)

The proposal: decompose the connected `tr F²` correlator / the KL spectral measure in the SU(N) representation ring; write

  `f₁² = |A₁|²`,  `A₁ = ⟨0| tr F² |G⟩ = Σ_r (character coefficient of the lightest `0⁺⁺` channel)`,

and extract a **lower** bound on `f₁²` from the sign/orthogonality structure of the characters `χ_r` and the Weingarten function `Wg(σ)` (the `(−1)`-signed inversion of the character convolution). The claim of novelty is that this is *sign* information (Bochner/Fourier on the group), not a variance/positivity lower bound — so it could certify the **atom** `ν({m₁²})>0` directly, escaping the wall that all positivity tools (RP/LSI/FKG/BL/Prékopa) only certify disorder and are blind to atoms.

To be a *tool* it must do all of: (1) be nonperturbative; (2) deliver a sign-controlling **lower** bound on `f₁²` (strict `>0`, at `m₁>0`); (3) be non-circular (not assume the gap / `m_g>0`); (4) fail on the free field (`f₁²=0`); (5) be uniform in `L, a` and reach `β≈2.5`. I show it fails (2),(3),(4),(5) and passes only (1) — collapsing onto already-banked P3 or dead RP.

The two horns it must choose between, *both* already named in the repo:

- **Horn P (positive).** On the Wilson lattice, `e^{(β/N)Re tr U_p}=Σ_r c_r(β)χ_r(U_p)` with **`c_r(β)≥0`** for all `β>0` (`PROOF-P3-UNCONDITIONAL.md` §1; `PROOF-FISHER-ZERO-FREE.md` §3: "at β>0 the Wilson character coefficients are nonnegative (no sign cancellation)"). The Weingarten signs live in the **inverse** area-law series, **not** in this forward one (`PROOF-P3-UNCONDITIONAL.md` §1, citing W2). So at strong coupling the coefficient `A₁` is a **positive sum** and `|A₁|²>0` is a modulus-square of a positive number — this is exactly Kotecký–Preiss / Strong coupling = P3.
- **Horn B (Bochner).** On the spectral side, requiring the KL measure `ν` to be a positive measure — which is what makes `f₁²` a "squared coefficient" at all — **is** OS reflection positivity (`BS-BLIND.md` §(c); `PROOF-RKHS-GAUSSIAN-PROCESS.md` §1: "`ν≥0` (OS reflection positivity)"). Bochner's theorem (positive-definite ⟺ positive spectral measure) is a re-writing of RP on the gauge group. That is precisely the "measure-positivity statement" the tool claims to escape.

The species has **no third home**. §1–§4 make each death exact.

---

## §1. KILL 1 — COLLAPSE: the "sign" is a phase killed by the modulus-square; what remains is KP (strong coupling) or RP (dead family)

**(a) The modulus-square obstruction (fatal, elementary).** `f₁²=|A₁|²` (`PROOF-P3-UNCONDITIONAL.md` §3: "`f₁²=|⟨Ω|O|G⟩|²` … is a modulus-square, so no sign can undo it"; `docs/yang-mills-frontier.md` line 68). The Weingarten function's alternating sign `(−1)^{|p∨q|+k}` (`PROOF-QUANTUM-GROUP-DEFORM.md` §2, Collins–Śniady) is a **phase** carried by `A₁`. Take any sign change `A₁→−A₁` (or any `U(1)` phase from a center twist / a summed signed surface): `|A₁|²` is **invariant**. The proposed mechanism — "the sign structure lower-bounds the *squared* coefficient" — is therefore a category error: the sign lives *before* the square and the square destroys it. The project states this as a verified structural theorem: `C8-SPECIES-NOGO.md` §(a) finds the minimal frustrated SU(2) plaquette has `Cov>0` *twist-invariant* — "frustration is invisible" — i.e. the character/center sign is provably decoupled from the quantity being bounded. **A sign cannot lower-bound a modulus-square.** The only sign-robust statement is `|A₁|≠0`, and that is *not* a sign argument (see (b)).

**(b) `|A₁|≠0` has exactly two sources, both dead rails.** To make `f₁²>0` you must show `A₁≠0`. On Horn P this is automatic and positive: `A₁ = a₀ · u^{p}` is a product of fundamental character coefficients `c_f>0` and the singlet recoupling `f⊗f̄→1` coefficient `1>0` (`PROOF-P3-UNCONDITIONAL.md` §3, "`a₀(β)>0` … `f⊗f̄→1`, whose coefficient is `1>0`"). That is **P3 / KP**, not a new signed-Fourier species — and it lives only at `β<β₁=O(1/N)` where the character/KP expansion converges absolutely. On Horn B, Bochner positivity gives only `ν≥0` hence `f₁²≥0` (non-strict); getting `>0` requires an atom, which is what is being sought (§4). So *every* route to strictness is "KP at `β<1/48`" or "assume the atom". **The species collapses to KP/strong-coupling** — exactly kill target 1.

**(c) The Weingarten sign is in the wrong sum.** `PROOF-P1P2.md` §2 and its *Sharpest line* state verbatim: "SU(N) Weingarten sign-alternation lives in the area-law **surface** sum, not in the per-tube character expansion (positive for `β>0`)." `AREALAW-GAP.md` line 11: the inverse (area-law) series "is Weingarten, whose local weights are **signed** …; positive `c_r≥0` live only in the forward series." The per-tube coefficient that fixes `m₁=ln(1/u)/a` is a **forward, positive** object; the sign sits in a *different* sum (the surface sum) which is **not positive** and is tamed only by absolute convergence = KP, or by CPS–CNS at small `λ`/large `N`, which the file itself flags as missing the positive regrouping. So the sign is **both** irrelevant (killed by the square) **and** located where it cannot bound the amplitude. **COLLAPSE lands.**

---

## §2. KILL 2 — FREE-FIELD LEAK: Bochner/character positivity is satisfied by the free field, so it returns `≥0` (vacuous) or `>0` (false)

The acid test is `f₁²(free Maxwell)=0`; a tool must **fail** there.

**(a) Horn B leaks outright.** Bochner positivity ⟺ `ν≥0`. The free field's KL measure `ν_free` is a **positive** measure (pure massless two-gluon continuum, threshold 0) — so *every* Bochner/positive-definiteness inequality holds on the free field with `f₁²=0` (`PROOF-RKHS-GAUSSIAN-PROCESS.md` §3: "free field … `ν_free` purely continuous … any *correct* detector returns `f₁²=0` — passes; a tool returning `f₁²>0` for `K_free` is vacuous"). Concretely, the Mercer/KL eigenvalues are positive functionals `λ_k=⟨w_k,ν⟩≥f₁²w_k(m₁²)`, giving `f₁²≤λ_k/w_k` — an **upper** bound (`PROOF-RKHS-GAUSSIAN-PROCESS.md` §2). A positivity/Bochner *lower* bound that returns `f₁²≥(positive)` on `ν_free` is **false**; one that returns `f₁²≥0` is **vacuous**. Either way the acid test is not passed by sign information — it is passed (if at all) by the non-Gaussian dynamics.

**(b) The finite-regulator leak (REDTEAM cross-cutting flaw).** `REDTEAM-FORGE.md` (Cross-cutting flaw): "the acid test `f₁²(free)=0` is true only in the `L→∞, a→0` limit; at any finite regulator the free lattice field has a positive lowest mode and nonzero residue." A signed-Fourier/character bound is a **finite-`L`, finite-`β`** object (it sums finitely many character coefficients). Evaluated at finite regulator it returns `>0` on the *free* lattice field ⇒ **free-field leak**. To make it vanish on the free field you must take `L→∞, a→0` — "the very limit that is the theorem" ⇒ circular (§4). **LEAK lands.**

---

## §3. KILL 3 — CROSSOVER IMPOSSIBILITY: at `β≈2.5` the character expansion is asymptotic/renormalon with no controllable sign

The bound must reach the Savvidy crossover `β≈2.5` where the physical atom sits.

- The character/positive expansion converges only inside the KP/Pringsheim radius `β<1/48`; the disk is **maximal** and ends on the positive real activity axis (`PROOF-FISHER-ZERO-FREE.md` §3). Beyond it, resummation meets a positive-axis singularity — the IR renormalon (`CONTINUUM-ASYMPTOTIC-BEYOND.md` §1–3). `β≈2.5` is uncovered **from both sides** (`PROOF-FISHER-ZERO-FREE.md` §2: "`β≈2.5` is uncovered from both sides"; large-`β` open for continuous SU(N), Balaban incomplete).
- The Borel/renormalon structure has **no definite sign**: lateral sums differ by `ΔS∝2πi Σ_n C_n e^{−2n/(β₀g²)}`, `C_n∝γ_n⟨O_{4n}⟩`; on ℝ⁴ the Stokes parameter `σ` is free (`CONTINUUM-ASYMPTOTIC-BEYOND.md` §2). Worse for a **lower** bound: the median sum is `g`-continuous and vanishes at `g=0` (`0↦0`), so "Borel summation … cannot manufacture `f₁²≥r₀>0` from a `g`-trivial start" (`CONTINUUM-ASYMPTOTIC-BEYOND.md` §4). A trans-series does not give a sign-definite coefficient lower bound; it gives an expansion fixed only up to the physics it omits (`CONTINUUM-ASYMPTOTIC-BEYOND.md` Sharpest line).
- The blocked measure itself turns **signed** at the crossover: RG blocking one e-fold "meets the same sign one block in … the *blocked* measure at the next scale is *signed* (Weingarten=Savvidy)" (`PROOF-MULTISCALE-4D.md` §3), so the positive-weight character control "is lost, not transported."

Hence any coefficient **lower** bound at `β≈2.5` is either (a) an uncontrolled truncation of a factorial-growth series (the `O(u^{p+1})` slop is the same order as the next term — `PROOF-P1P2.md` §3), or (b) circular on the condensate `⟨Tr F²⟩` itself. **CROSSOVER ⟂ lands.**

---

## §4. KILL 4 — CIRCULARITY: a lower bound on the *atom* already consumes the second level `E₂` / the gap

The tool promises to certify an **atom** (coherence). The project separates three objects that the claim conflates:

1. **Residue `I1`** (`f₁²>0`) — proved at `β<1/48` by P3.
2. **Isolation `I2`** (`inf supp(ν−f₁²δ_{m₁²})−m₁²≥δ`) — **not** obtained; "the cost grading is a statement about polymer *amplitudes*, not about transfer-matrix eigenvalue spacing, and cannot be converted into a uniform isolation gap" (`PROOF-P1P2.md` §3; `PROOF-P3-UNCONDITIONAL.md` §5: "an earlier draft claimed the grading yields `I2` — that is **wrong**").
3. **Atom at `m₁>0`** — needs both: `f₁²>0` **and** an isolated pole at `m₁>0`. `MS-BLIND.md` §3: "`f₁²>0` does **not** imply `m>0`: a massless excitation carries `f₁²>0` at a pole `s=0` … the residue is the intensity; the gap is its lower support edge."

So "certify the atom" `=` `I1 ∧ I2 ∧ m₁>0`. The only *proven* lower bound on `f₁²` beyond the trivial, and the only way to localize it, **consumes `E₂`**: the residue–level-gap lemma `f₁² ≥ (ΛM₀−M₁)/(Λ−m₁²)`, `Λ=E₂` "consumes the second level, the one non-moment input" (`INEQ-R10-sumrule-positive.md`, last line; `GAP-G2-edge`). And to bound the *binding coefficient* `A₁` as a genuine glueball overlap you must define the two-gluon threshold `≥2m_g`, which presupposes the **constituent-gluon mass `m_g>0`** — H1 itself. Equivalently: `INVENT-TRANSPORT-INEQ` transports the residue **sign** with `f₁²=Z_O²f₁²(β_eff)` but `Z_O²` "multiplies a base value that is the spectral-edge atom" (`CONTINUUM-ASYMPTOTIC-BEYOND.md` §3) — it never *produces* the positive base point. **To lower-bound the atom you must already assume `E₂` / `m_g>0` = the gap.** **CIRCULARITY lands.**

---

## §5. The class is provably empty (this is why it is dead, not merely hard)

`GOLDS-POS-MASS.md` proves (uniqueness lemma) that within "positivity + analyticity + model-independent symmetry": the **only** bridge forcing a specific atom to nonzero weight is the **Goldstone** theorem, whose conservation law pins the atom at **`m=0`** (`Qu|0⟩` has exact zero energy, weight `‖Qu|0⟩‖²>0`). A **positive-mass** atom has no symmetry anchor — the only exception is a **SUSY BPS** state (`{Q,Q̄}∼P`, central charge), which pure SU(N) lacks. Pure YM has neither anchor: center/1-form symmetry is discrete and adjoint-blind, gluon condensate is a Lorentz-scalar vev that does not break Poincaré, and scale symmetry is explicitly broken by the anomaly. `EDGE-f1sq-attack.md` and `INEQ-T5-sumrule-hunt.md` show the sign is structural: every OS/Cauchy–Schwarz/Hölder/positivity inequality runs `(atom)≤(aggregate)` or `(cross)²≤(same)(same)` — an **upper** bound — and forcing `f₁²>0` needs the *Schur-convex reversal* `total≤C·atom`, which "never appears" (`INEQ-R10-sumrule-positive.md` §the sign lemma). The peer-reviewed confirmation is exact: Caron-Huot–Pokraka–Zahraee, JHEP **01** (2024) 195: "positivity alone forces the minimal residue to zero"; a genuine lower bound "requires extra input" (= the gap).

The proposed tool asserts it delivers exactly the one forbidden thing — a positivity-class lower bound on a **positive-mass** atom — *without* being positivity. But §0–§1 show it has no non-positivity home: Horn P is positivity (`c_r≥0`, KP), Horn B is positivity (`ν≥0`, RP). It is not in a new class; it is RP/KP re-labeled, and the object it claims (positive-mass atom from sign structure) sits in the class `GOLDS-POS-MASS.md` proves empty.

---

## §6. Scorecard (five hard constraints) — survival 1/5

| # | Constraint | Verdict | Exact failure |
|---|---|---|---|
| 1 | Nonperturbative | **PASS** | measure/representation-level; no loops. (Fails at crossover — §3.) |
| 2 | Sign-controlling **lower** bound on `f₁²` | **FAIL** | sign is a phase killed by `|A₁|²` (§1a); Bochner/positive gives **upper** (`f₁²≤λ_k/w_k`, RKHS §2; INEQ-R10) |
| 3 | Non-circular | **FAIL** | atom needs `I2`/`E₂`/`m_g>0` = gap (§4); residue–level-gap consumes `E₂` |
| 4 | Non-SUSY + **fails free field** | **FAIL** | Bochner/character positivity holds on `ν_free` with `f₁²=0` (≥0 vacuous or >0 false); finite-regulator free field has spurious residue (§2) |
| 5 | Uniform `≥λ` (reach `β≈2.5`) | **FAIL** | character series asymptotic/renormalon, no controllable sign at 2.5; 0↦0 under Borel (§3) |

**Survival score: 1/5** (only nonperturbativity). Compare `REDTEAM-FORGE.md`: this joins the dominated field with PROOF-QUANTUM-GROUP-DEFORM (positivity aggregate, Weingarten sign returns at `q→1`) and PROOF-RKHS-GAUSSIAN-PROCESS (kernel positivity = RP, upper bound) — the two closest prior candidates, both dead for the same reason.

---

## §7. What, if anything, survives (honest residue of the attack)

- **Nothing new.** The only true statement inside the proposal is: at `β<1/48`, `A₁` is a product of **positive** fundamental recouplings, so `f₁²=|A₁|²=|a₀|²u^{2p}>0` — which is exactly `PROOF-P3-UNCONDITIONAL.md` (I1, residue), already banked, and it is KP/strong-coupling, not signed-Fourier. Note even this is `I1`, not the atom (no `I2`).
- The negative-space content — "Bochner/character positivity is sign-blind to the atom and gives upper bounds" — is already the banked statement of `PROOF-RKHS-GAUSSIAN-PROCESS.md`, `INEQ-R10`, `INEQ-T5`. It is not a tool; it is another witness to the wall.
- **The one genuinely-open object the attack clarifies:** a *signed sum that lower-bounds the modulus-square of a non-vanishing per-slice character coefficient `c_R` without consuming `E₂`*. The sign lives (only) in the area-law **surface** sum, which is the wrong sum; the relevant per-slice coefficient is **forward-positive** (so it is KP, not signed); and the square erases any residual phase. Constructing the missing signed-sum coercion is, by `SCHUR-CONVEX.md`/`EDGE-f1sq-attack.md`, the edge theorem itself.

---

## §8. Sharpest line

**The proposed species dies on an elementary algebra fact the whole campaign already banked: `f₁²=|A₁|²` is a modulus-square, and a sign is a phase — the Weingarten/character sign lives in the amplitude, `|·|²` deletes it, and the sign that *is* real (the signed area-law surface sum) is in the wrong sum and cannot reach the per-slice `c_R` that fixes `m₁`. What is left, when made precise, is either the *positive* character expansion `c_r≥0` (which is Kotecký–Preiss / P3, only at `β<1/48`) or *Bochner positivity* `ν≥0` (which is OS reflection positivity — the exact "measure-positivity" family the tool claimed to escape), and both are sign-blind to the atom and bound it from **above**. There is no third object.**

**Exact stuck step.** A signed sum `S` with `f₁² ≥ Φ(S) > 0` in which the sign of `S` is controlled without (a) absolute convergence (=KP, `β<1/48`), (b) `ν≥0` (=RP, upper bound), or (c) an assumed `E₂`/`m_g>0` (=the gap). Producing it is equivalent to the Schur-convex coercivity `total ≤ C·atom` — the edge theorem, i.e. the Clay gap. The species does not exist.

---

*Refs.* Internal (read/cited as ground, not invented): `docs/yang-mills-frontier.md` line 68 (`f₁²` modulus-square, Weingarten sign in area-law sum); `PROOF-P1P2.md` §2,§3, Sharpest line (positive `c_r`, Weingarten in surface sum; I1 not I2); `PROOF-P3-UNCONDITIONAL.md` §1,§3,§5 (`c_r≥0`, `a₀>0` positive recoupling, modulus-square, `I2` correction); `PROOF-FISHER-ZERO-FREE.md` §2,§3 (character coefficients nonnegative for `β>0`; `β≈2.5` uncovered both sides; Pringsheim); `PROOF-MULTISCALE-4D.md` §3 (blocked measure signed Weingarten=Savvidy); `CONTINUUM-ASYMPTOTIC-BEYOND.md` §1–4 (`0↦0`, renormalon no sign, `⟨TrF²⟩`); `MS-BLIND.md` §3 (residue ≠ gap; massless pole carries `f₁²>0`); `PROOF-RKHS-GAUSSIAN-PROCESS.md` §1–3 (Bochner/`ν≥0`=RP; Mercer gives **upper** `f₁²≤λ_k/w_k`; free field passes with `f₁²=0`); `INEQ-R10-sumrule-positive.md` (positive weight ⇒ upper bound; residue–level-gap consumes `E₂`); `INEQ-T5-sumrule-hunt.md` (Weingarten = mass ordering via RP+CS, upper on coupling; positivity ⇒ `f₁²≥0` trivial); `EDGE-f1sq-attack.md` (`(atom)≤(aggregate)` sign; needs Schur-convex reversal); `GOLDS-POS-MASS.md` (only positivity→atom bridge is Goldstone at `m=0`; massive atom needs SUSY BPS); `SCHUR-CONVEX.md`, `INEQ-T3-schur-convex-nonprot.md`; `C5-NEG-ASSOC.md`, `C6-DEPENDENCE-CLASS.md`, `C8-SPECIES-NOGO.md` (torsion ⇒ no order; minimal plaquette `Cov>0`, twist-invariant); `PROOF-QUANTUM-GROUP-DEFORM.md` §2 (Weingarten sign `(−1)^{k+|p∨q|}` returns at `q→1`); `AREALAW-GAP.md` line 11 (Weingarten signed, forward positive only); `REDTEAM-FORGE.md` (cross-cutting finite-regulator acid-test flaw; RKHS/quantum-group precedents). External (pre-vetted in-repo): Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737; Kotecký–Preiss, CMP **103** (1986) 19; Källén, Helv. Phys. Acta **25** (1952); Lehmann, Nuovo Cim. **11** (1954); Weingarten, PRL **51**,1830 ('83); Collins–Śniady, CMP **264** (2006) 773; Caron-Huot–Pokraka–Zahraee, JHEP **01** (2024) 195, arXiv:2309.04472; Bochner, Math. Ann. **108** (1932) (positive-definite ⟺ positive spectral measure); Aronszajn (1950), Trans. AMS **68**, 337 (RKHS/Mercer).
