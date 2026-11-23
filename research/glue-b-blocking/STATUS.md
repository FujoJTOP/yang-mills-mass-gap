# STATUS — glue-b-blocking (GL2a, GL2b, GL2b′, GL5)

> 2026-11-22 · 4D pure YM · full report: `REPORT.md` · checks: `python -I verify.py`, `bash research/glue-b-blocking/check.sh`
> Scope: only `research/glue-b-blocking/`. No other directory touched.

## 0. Verified Lean artefacts

**One command:** `bash research/glue-b-blocking/check.sh` → `GL2a + GL2bPrime COMPILE (no sorry, standard axioms frozen)`, **exit 0**.

| file | theorems | axioms |
|---|---|---|
| `lean/GL2a.lean` | `rp_descends`, `compression_form_nonneg` | `[propext, Classical.choice, Quot.sound]` |
| `lean/GL2bPrime.lean` | `gap_lipschitz`, `gap_tendsto`, `gap_ge_liminf` | `[propext, Classical.choice, Quot.sound]` |

Axioms are frozen in-file by `#guard_msgs`. **Both guards were tested negatively**: injecting a `sorry` makes the `#guard_msgs` fail with `sorryAx` in the printed list; the comment-stripping source scan flags a real `sorry` and ignores mentions inside docstrings/line comments.

**Declared scope, not `sorry`** (REPORT §2.8, §3.8.1):

1. GL2a is formalised as the abstract theorem. The instantiation check — path-level (B1) for the face-average equivariant blocking of gluing §1.7 — is **not** machine-checked.
2. GL2b′ uses `SPEC-RES` in its **norm-convergence** form `‖H_k − H‖ → 0` (constrains the operators, not `gap`/`liminf`; stronger than the paper's strong-resolvent hypothesis). The strong-resolvent-only version needs Reed–Simon I Thm VIII.24(b), absent from Mathlib, and remains a paper remark.
3. `gap` is the **variational** gap; identifying it with `inf(σ H \ {0})` needs `ker H = ℝΩ` (vacuum simplicity) and is **not** formalised — nor needed downstream.

## 1. Ledger (honest)

| Lemma | Verdict | Where |
|---|---|---|
| **GL2a** | **✅ PROVEN** (abstract) — 5-line proof; all content is (B1)+(B2); criterion `A∈W∩RWR⁻¹` (+ offset `v∈Λ'`) verified in exact arithmetic | REPORT §2 |
| **GL2b** | **⚖️ DECIDED** — interlacing direction proven; the *needed* lower-bound transfer is **false** without a visibility hypothesis; truncation error cannot flip the sign; replaced by GL2b′ | REPORT §3.1–3.6 |
| **GL2b′** | **✅-IMPLICATION** — `SPEC-RES ⟹ Δ_cont ≥ liminf_{a→0} Δ_a` (4 lines via `σ(H) ⊆ limsup σ(H_a)`). Necessity counterexample: **GL3 cannot supply it** (same limit, arbitrary gap). It is `≥`, not `=`. | REPORT §3.8 |
| **GL5-KIN** | **✅ STRICT** — `SYM-DIM4` (Lean'd) **+ `OP-CLASS`** (added here, machine-checked `S₄` family exchange). Note `SYMDim4.lean` explicitly scopes `OP-CLASS` out, so the two must be booked together. | REPORT §4.7.1 |
| **GL5-DYN** | **✅-IMPLICATION** — `IN-2 (≡ U2 ≡ sym) ⟹ OS2`; translations free, all content is `IN-2`. This is bookkeeping, not a theorem. | REPORT §4.7.2 |

Count toward the team goal: **GL2a proven; GL5-KIN strict; two ✅-implications (GL2b′, GL5-DYN) whose hypotheses (`SPEC-RES`, `IN-2`) are open and not among the 16 inputs; GL2b decided as a negative result.** I do not count hypotheses as proved.

## 2. The three things the coordinator should act on

1. **GL2a is not where the work is.** No estimates enter. The RG difficulty is Balaban's class preservation, not RP descent. Remaining item is a finite check: (B1) path-level for the face-average equivariant blocking (§1.7). Hand to `lean-infra` or `gauge-invariance`.
2. **GL2b must be re-worded.** `Δ_coarse ≥ Δ_fine` holds in *every* reading (compression, exact marginal, truncated+ε). A coarse-scale lower bound is an **upper** bound on the physical gap. Error control does not fix it. Actionable: state `CPB-W`/`CPB-S` on the physical theory at its own cutoff, and carry `GL2b′` (`Δ_cont ≥ liminf_{a→0} Δ_a`) as the bridge — a limit statement, strictly weaker than truncation norm bounds. `GL9` is unaffected.
3. **GL5 route A has a tautology dilemma.** If "inclined lattice" means the spacing-5 alternative regularization, then `5Rℤ⁴ ⊂ ℤ⁴` (index 625) makes the inclined limit automatically `R_*μ`, so invariance under `RWR⁻¹` is *equivalent* to the known `W`-invariance — no new information. Route A is therefore either tautological or it is BLK-IND. Only the **density lemma (Prop D)** is done. Recommend route B (NP-SYM), whose gap is local, quantitative and lattice-falsifiable.

## 3. Falsification tests (all in REPORT, two machine-checked)

- GL2a: non-θ-covariant blocking crossing `t=0` must give a negative `⟨θA,A⟩`; `B₁₂` must fail the criterion. ✅ (`verify.py`)
- GL2b: `Δ_c = +∞` counterexample must be reproduced by anyone claiming `Δ_c ≤ Δ_f + C`; lattice: blocked `ξ_eff(t)` must approach the plateau from **above**. ✅ structural
- GL5: anisotropic lattice YM (`β₀≠βᵢ`, symmetry `(ℤ/2)⁴⋊S₃`) must have a preferred frame — proves OS2 is non-vacuous; measure anisotropy scaling `∝ a²` **and** its `L`-independence at fixed `a`. ✅ structural

## 4. Lean status (DONE — supersedes the earlier "Lean-able core" sketch)

Both targets are formalised and hang together in one script (see §0). The earlier
sketch is obsolete; what was actually proved differs in two instructive ways:

- **GL2a** went through exactly as sketched, but *without* `condExp`: the only used
  property is `integral_map`. The "θ-covariant conditional expectation" is a
  consequence of (B1), not a hypothesis.
- **GL2b′** is *not* formalised as interlacing (that is GL2b, and false in the needed
  direction). It is formalised as the **Lipschitz estimate for the variational gap**,
  `|gap H Ω − gap K Ω| ≤ ‖H − K‖`, from which `gap (H k) → gap H` and the `≥ liminf`
  form follow with no spectral theory at all. This is the machine-checked reason the
  GL2b′ route is cheaper than GL2b.
- Still **not** formalisable today (declared, not `sorry`): Reed–Simon VIII.24(b)
  (strong resolvent ⇒ spectrum inclusion) is absent from Mathlib; the variational-gap
  ⟷ spectral-gap identification needs `ker H = ℝΩ` and is skipped.

Do **not** attempt `Δ_coarse ≥ Δ_fine` as a theorem about measures — it is not true
without visibility (Prop 3); formalize only the finite-dimensional compression
statement, which is `compression_form_nonneg` in `lean/GL2a.lean`.

## 5. Blockers / confidence

- **Blockers:** `SPEC-RES` and `IN-2` are both open; neither is in the 16 inputs. (Old `BLK-IND` archived — route A withdrawn, R28. Old `NP-SYM` refuted, R33.)
- **High confidence:** GL2a (short proof, independently consistent with `os3` §4.6/§5); GL2b's direction and the failure of truncation to repair it (elementary, machine-checked); GL2b′'s implication (one-line σ-inclusion); GL5-KIN; the bookkeeping `IN-2 ≡ U2 ⟹ OS2`.
- **Medium:** whether `(Q-CONS)` closes `IN-2` (route D, `REVERSAL-OPEN`); whether `SYM-4D + SYM-IR` can be upgraded from approximate U2 to U2.
- **Unknown, not guessed:** truth of `SPEC-RES`; truth of `IN-2` / `(Q-CONS)`.
- **Would change verdicts:** a construction proving strong resolvent convergence of `H_a` (closes GL2b′/`hcert`); a counterexample satisfying `Conserved.symmetric + StressLaws + MG` but breaking `U2` (reverses the route-D reduction); lattice evidence of `a`-linear (not `a²`) anisotropy at fixed `L` (falsifies `SYM-4D`).
