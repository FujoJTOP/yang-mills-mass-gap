# STATUS — glue-c (GL1, GL6, GL8)

> 2026-11-22 · full report: `REPORT.md` · scope: 4D pure Yang–Mills.
> **One-command check:** `bash research/glue-c-gauge-group/check-glue-c.sh` → last line `ALL GLUE-C CHECKS PASS`.

## 0. Verdict

| Lemma | Status | Evidence | Named gap / caveat |
|---|---|---|---|
| **GL6** gauge-fixing independence | ✅ **proven and machine-checked** (finite lattice + orbit-space RG) | `GL6Fubini.lean` (`orbit_sum_identity`, `gauge_fixing_independence`, S₃ `native_decide` example; axioms `[propext, choice, Quot.sound]`); `gl6_check.py` (1296-config cross-check) | Balaban quantitative class preservation (L2.5) is **not** GL6; reflection covariance is GL2a (glue-b); the orbit-space measure `μ̄` is not constructed in Lean (not needed for the finite layer). |
| **GL8** group-uniform constants | **✅ assembly (unconditional) + ✅-implication (supply)** | `GroupUniform.lean` (`uniform_over_groups`, `supply_min_pos`; axioms `[propext, Classical.choice, Quot.sound]`); `gl8_check.py` | GL8-A: group-blind assembly. GL8-S: `IN-3 ∧ IN-4` supply `(δ_W,δ_S)` per `G` ⟹ `min>0`; `IN-3`/`IN-4` are **inputs**, not a GL8 gap. GL8-U (class-uniform) needs `IN-8`. |
| **GL1** single measure family | ✅ **9/9 (measure + declarations)** | `REPORT.md` §2.7 re-audit; `gl1_check.py` | None remaining. E1–E4 plus the final four declarations are closed. |

**Honest count:** 3 of my 3 lemmas are ✅ at the level stated (GL1 is a compliance audit, now 9/9; GL6/GL8 have proofs/artifacts). Team strict count 3/9 → 4/9 per the peer's ledger.

## 1. GL1 — final re-audit (9/9)

The shared measure depends only on `γ := β/N_c` (coefficient of `Σ_p Re tr U_p`).
Full table in `REPORT.md` §2.7: **measure 9/9, declarations 9/9.**

| Interface | Written weight (post-fix) | `γ` | Measure | Declaration |
|---|---|---|---|---|
| I-OS2 (`os2` §1.1) | `exp(β Σ(1 − (1/N_c)Re tr))`, `β=2N_c/g²` | `β/N_c` | ✅ | ✅ same RT/C |
| I-OS3 (`os3` §2.1) | `exp((β/N_c)Re tr_f)`, `β=2N_c/g²` | `2/g²` | ✅ | ✅ (T⁴ bridge declared) |
| I-UNIQ (`uniqueness` §1.1) | `exp((β/N_c)Σ Re tr)` | `β/N_c` | ✅ | ✅ (RT + W1 bridge) |
| I-NT (`nontriviality`) | explicit A0 reference | `β/N_c` | ✅ | ✅ |
| I-GAUGE (`gauge-invariance` §2.1) | `β Σ(1 − (1/N_c)Re tr)`, `β=2N_c/g²` | `β/N_c` | ✅ | ✅ `G̃`, `T⁴_L`, `θ=0`, RT/C |
| I-GROUP (`group-uniformity` §2.0) | `β Σ(1 − (1/dim f)Re tr_f)`, `β=2dim f/g²` | `β/dim f` | ✅ | ✅ `T⁴_L`, `θ=0`, RT/C |
| I-HANDOFF (`handoff-proof` §2) | `T⁴_L`, `β=2N_c/g²` | `β/N_c` | ✅ | ✅ |
| I-CONF (`confinement` §1.1) | `exp((β/d_f)Re tr_f)`, `β=2d_f/g²` | `2/g²` | ✅ | ✅ L7 `2d_f`, `C₂(adj)`, RT/C |
| I-SPEC (`mass-spectrum` §1.1) | explicit A0 reference | `β/N_c` | ✅ | ✅ |

Concrete: `SU(3)`, `g²=1`, `β_A0=6`; a session writing the bare form `exp(k Σ Re tr)` must use `k=2`, not `6` (factor 3).

**Remaining declarations:** none. A0 §1.2's E4 is now the two-form text (normalized `β=2N_c/g²`; bare `k=2/g²=β/N_c`), consistent with its `6 vs 2` example; the leak into §3.5 and the handover note is cleaned. The "4 families + 5 exceptions" enumeration is fixed in A0 C6/R6/§2.6; `gl8_check.py` (G6) verifies it by rank (`1,3,3,5,4,5,5,5`).

## 2. GL6 — what is proven and what is not

**Machine-checked (`GL6Fubini.lean`, no sorry, axioms frozen):**
- `orbit_sum_identity`: general finite orbit Fubini
  `∑_x F x μ x = ∑_{a∈A} F a μ a · |orbit a| / |A ∩ orbit a|`;
- `gauge_fixing_independence`: two slices meeting every orbit give the same expectation;
- S₃ on `Fin 3` (non-free, multiplicity 1 vs 2): exact rational identity by `native_decide`.

Proven on paper (elementary, no open problem):
- finite-lattice orbit identity `Σ_X F dμ = Σ_{a∈A} F(a)μ(a)|Γa|/m(a)`, `m(a)=|A∩Γa|`;
- compact-group orbit-space identity `∫_X F dμ = ∫_{X/Γ} F̄ dμ̄` (no section);
- RG compatibility: exact-equivariant `Φ` ⇒ `ν` invariant, `Φ̄` descends, `ν̄ = Φ̄_*μ̄`, iteration, defect bound;
- Balaban's internal gauge fixing is removable (GL6-FP + GL6-RG(iii)).

Not proven (and not part of GL6): quantitative class preservation (L2.5); exact equivariance is proven in `gauge-invariance` §2.6, not re-proved here; reflection covariance is GL2a.

Executable check (`S3`, 2×1 torus, 1296 configs, 49 orbits):
`E[F_inv] = 0.430665463009` exactly; slice with `|orbit|/m` gives the same; dropping `1/m` gives `0.250303751724`; multiplicity `m ∈ {1,2,3,4,6,12}`; non-invariant `F` is slice-dependent (exact 0, slice 0.52).

## 3. GL8 — split as GL9: assembly (✅) + supply (✅-implication)

- **GL8-A (assembly, unconditional ✅).** `GroupUniform.lean`: `handoff_gap`, `uniform_over_groups`, `min_pos`, `class_uniform`. The statement never mentions the group index type. Axiom footprint frozen.
- **GL8-S (supply, ✅-implication).** `GroupUniform.lean`: `WeakSideInput` (IN-3), `StrongSideInput` (IN-4), `supply_min_pos`. Type: `IN-3 ∧ IN-4` per `G` ⟹ `min(δW,δS) > 0` per `G`. The `bound` fields are the named inputs; GL8 does not prove them. `overlap : q1 ≤ q0` is `IN-5`, needed by GL8-A's covering but not by GL8-S's positivity.

Numbers (`C_RG = 24π²/11`, `m/Λ=8`, `q₁=1`, fixed `g²=1`): `q₀(SU(2))=1.69e−4`, `SU(3)=6.11e−3`, `G2=0.037` — **all fail at one loop**. At fixed `λ*=10.4` all groups pass (`q₀=1.009`), but the required `g*²` is `5.2` (SU(2)), `3.47` (SU(3)), `2.60` (G2). The hard inputs are `IN-3` (CPB-W, ≥2-loop `t_up`, `G`-uniform), `IN-4`, `IN-5`, and `(H_m)` for GL8-U.

GL8-A/S ⇏ **GL8-U** (class-uniform): `δ_n=1/(n+1)` is a counterexample. A0's one-line "min>0 per G" covers only S1.

**Consumption list.** GL8-A consumes per-`G` `HandoffData` (four constants, monotonicity, `overlap`=`IN-5`). GL8-S consumes `IN-3` (weak side) and `IN-4` (strong side); it produces `(δW,δS)`, which feed GL8-A. `IN-3` is the physical wall and belongs to `handoff-proof`/`IN-8`, not to GL8.

## 4. What I need next

1. **GL1: closed.** E1–E4 and the four declaration items are all in the source; A0 §1.2's E4 two-form text is correct. The only maintenance item is to keep `N_c := dim f` out of `b₀` (use `C₂(adj)`), which `group-uniformity` §2.0 and A0 §1.4 now warn about.
2. **`IN-3`/`IN-4`'s `G`-uniformity**: does the weak-side (`CPB-W`) analysis have ledger `(0,0,0)` in `(dim G, h^∨, dim f)`? Even a conjecture-level statement would let GL8-U be stated conditionally.
3. **`IN-5` (`q₀≥q₁`) and `(H_m)`**: the remaining inputs for GL8-U.
4. **Optional Lean**: the finite GL6 orbit identity is still a Python check; a ~150-line Lean file would upgrade it from "reproducible computation" to "compiled artifact".

## 5. Confidence

- **High:** GL6 finite/compact/RG statements (proofs self-contained; executable check). GL8 assembly (Lean, axiom-clean). The GL1 dictionary and the compliance verdicts (line-referenced).
- **Medium:** GL1's audit is text-line-referenced (the reports can drift); re-run §2.7's row checks before the conditional theorem document is frozen.
- **Low / not claimed:** `IN-3` (`CPB-W`) truth; `IN-4`; `IN-5` (`q₀≥q₁`) for any group; GL8-U.
