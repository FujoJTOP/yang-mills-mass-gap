# TRANSPORT-PACKING-AUDIT — Is the physical instanton vacuum DENSE (𝔡 = n ρ̄⁴ ≳ 1) or DILUTE (≪ 1)? Audit of the "packing fraction = O(10)" premise in `E2-FLIP-TRANSPORT.md`

> 2026-11-24 · Independent audit of the linchpin numerical claim in `research/unconditional/E2-FLIP-TRANSPORT.md` §4.1: *"the one-loop instanton gas is DENSE at the mass-generating scale, 𝔡 := n ρ⁴ ≈ O(10) for SU(3)."* That claim is the stated premise of the `[NG]` no-go ("diluteness fails exactly where the mass is generated"). This note settles the number from scratch, cross-checks ILM and lattice, and reports what it does to the no-go.
>
> **Labels (`CLAUDE.md §2`):** **[T]** 真定理 / **[C]** 条件定理 / **[O]** 精确缺失步 / **[V-here]** fresh computation this session (`research/unconditional/w7/transport_packing.py`) / **[R]** external as-reported, not re-verified here. **No QED on anything unproven.**

---

## VERDICT UP FRONT

> ### **DILUTE-to-BORDERLINE — the O(10) is an ARTIFACT; the no-go's DILUTENESS premise is FALSE.**
>
> The physical instanton vacuum is **dilute**: the best empirical numbers give **𝔡 = n ρ̄⁴ ≈ 0.01** (ILM, Schäfer–Shuryak `[R]`) to **≈ 0.06** (literal one-loop full integral, `[V-here]`; and cooling-based larger sizes, `[R]`). Honest range **𝔡 ∈ [0.008, 0.07]** — always **≪ 1**. There is **no controlled computation** that yields 𝔡 ≳ 1.
>
> The corpus's **𝔡 ≈ O(10)** arises from **two compounding errors** in `E2-FLIP-TRANSPORT §4.1`:
> 1. **Wrong coefficient `[T-computable]`.** It assumes the one-loop determinant prefactor is `C = O(10⁻¹)`. The actual 't Hooft/Callan–Dashen–Gross coefficient for SU(3) is **`C_3 = 0.466 e^{−1.679·3}/(2!·1!) = 1.51×10⁻³`** — **~66× smaller**. This alone turns a would-be 𝔡 = 3.7 into 𝔡 = 0.056.
> 2. **Wrong quantity `[T-computable]`.** The number `(11)⁶e^{−11} = 29.6` that the note multiplies by `C` is **`D(ρ*)ρ*⁵ = C_N S*^{2N}e^{−S*}`, a density-per-unit-size with dimension `1/length`** — it is **not** the packing fraction `n ρ̄⁴` (which needs `n = ∫D(ρ)dρ`, and has **no** `ρ*⁵`). Even at face value the note's "nρ*⁴ = C·29.6" is off by a factor `ρ*` and mislabels the object.
>
> **What is RIGHT in the note:** the peak is genuinely at **ρ* = e^{−2N/(b−5)}/Λ = 0.368/Λ ≈ 0.37 fm**, i.e. the mass-generating instantons really are **large = IR** (this matches the measured ILM size ρ ≈ 0.33 fm — see §3). So "mass is generated at ρ ~ 1/Λ" is correct. It is only the second half — "and the gas is **dense** there" — that is false by 2–3 orders of magnitude.
>
> **Consequence for the no-go.** The `[NG]` in `E2-FLIP-TRANSPORT` as **stated** (premise "𝔡 = O(10) ⇒ diluteness fails ⇒ KP does not converge") **does not survive**: the diluteness premise is numerically false. **However, this does NOT reopen the dilute route as a theorem.** Kotecký–Preiss convergence requires control of the *activity × range* in the continuum, `Σ_p |z_p| e^{a|p|} < ∞` — a statement about `z_p ~ e^{−S_p}` and geometric range at ρ ~ 1/Λ, **not** about `n ρ⁴`. Small 𝔡 is *necessary but not sufficient*. So the no-go's **conclusion** (no rigorous cluster expansion at the physical scale) may still stand — but it must be **re-derived on correct grounds** (continuum range/overlap control = the confinement problem), **not** asserted from a spurious O(10) density.

---

## §1 The exact one-loop instanton size distribution, from scratch

**The formula `[R]` ('t Hooft 1976; Callan–Dashen–Gross 1978; Bernard 1979; Schäfer–Shuryak RMP 1998 Eq. 2.13).** For 4D SU(N), one-loop instanton density:

```
   D(ρ) dρ = C_N · [ 8π²/g²(ρ) ]^{2N} · e^{−8π²/g²(ρ)} · dρ / ρ⁵        (1)

   C_N = 0.466 · e^{−1.679 N} / [ (N−1)! (N−2)! ]                        (2)
```

Every factor accounted for: the `ρ^{−5}` is the size measure (four translational + one dilatation zero modes → `ρ⁵` in the measure); `[8π²/g²]^{2N}` are the `2N` nonzero-mode fluctuation determinants; `e^{−8π²/g²}` is the classical action weight; `C_N` is the collective zero-mode/GH determinant (`0.466 = 4.6/π²`, `1.679 = 16/... ` scheme constants). One-loop running (pure YM, `N_f=0`):

```
   8π²/g²(ρ) = S(ρ) = b · ln(1/(ρΛ)),     b = 11N/3.                       (3)
```

**Definitions (explicit, as demanded).** The density per unit 4-volume per unit size is `D(ρ)` with measure `dρ` (not `dρ/ρ`); the total number density is

```
   n := ∫_0^{1/Λ} D(ρ) dρ ,          ρ̄ := ∫ρ D dρ / ∫D dρ ,          𝔡 := n · ρ̄⁴.       (4)
```

`𝔡` is the packing fraction; it is dimensionless and (by construction, since `n ∝ Λ⁴`, `ρ̄ ∝ 1/Λ`) **Λ-independent** — a pure function of `N` and of the one-loop convention.

**The exact coefficient `[V-here]`.** Eq. (2) gives (script §1):

| group | `b=11N/3` | `C_N = 0.466 e^{−1.679N}/((N−1)!(N−2)!)` |
|---|---|---|
| SU(2) | 7.333 | 1.62×10⁻² |
| **SU(3)** | **11.000** | **1.513×10⁻³** |
| SU(4) | 14.667 | 4.70×10⁻⁵ |
| SU(5) | 18.333 | 7.31×10⁻⁷ |

**The note's assumed `C = O(10⁻¹)` is wrong for SU(3) by a factor ≈ 66** (and by ≈ 200 if one takes the note's "O(10)" as needing `C ≈ 0.3`). This is the single largest source of the O(10).

---

## §2 Peak, and the full-integral packing fraction

**Peak `[V-here]`.** Extremising Eq. (1) with Eq. (3): `d/dρ[S^{2N}e^{−S}ρ^{−5}] = 0` gives `S* = 2Nb/(b−5)`, hence

```
   ρ* = e^{−2N/(b−5)}/Λ ,        S* = 8π²/g²(ρ*) = 2Nb/(b−5).              (5)
```

For SU(3): `S* = 2·3·11/6 = 11.0` **exactly** (which is why `e^{−11}`), `ρ* = e^{−1}/Λ = 0.368/Λ`. This reproduces the note's `ρ* = 0.368/Λ` and `e^{−S*} = e^{−11} = 1.67×10⁻⁵`. ✅ *The peak location is correct.*

**Where the note's 29.6 comes from `[V-here]`.** At the peak,

```
   D(ρ*)ρ*⁵ = C_N · S*^{2N} e^{−S*} = C_N · (11)⁶ e^{−11} = C_N · 29.588 .   (6)
```

The number **29.588 is `D(ρ*)ρ*⁵`** — it carries dimension `L⁻¹` and, crucially, it is the **density per unit size at a point**, not an integrated density. The note writes it as `"n ρ*⁴ = C · 29.6"`, i.e. it (a) substitutes the *peak value* `D(ρ*)` for the *integral* `n = ∫D dρ`, and (b) keeps a `ρ*⁵` where a packing fraction has `ρ̄⁴`. With the **true** `C_3 = 1.513×10⁻³`, even this quantity is `C_3 · 29.588 = 4.48×10⁻²` — not O(10).

**The actual packing fraction `[V-here]`.** Because `D = C_N ρ^{b−5}(b ln(1/ρ))^{2N}` (Λ=1), the full integral is elementary:

```
   n   = C_N · b^{2N} · (2N)! / (b−4)^{2N+1}                        (7)
   ρ̄   = (b−4)^{2N+1}/(b−3)^{2N+1},      𝔡 = n ρ̄⁴.                 (8)
```

Evaluated (script §3; analytic = numeric to 4 digits):

| group | n (Λ⁴ units) | ρ̄ (1/Λ) | **𝔡 = n ρ̄⁴** |
|---|---|---|---|
| SU(2) | 2.736 | 0.269 | **0.0144** |
| **SU(3)** | **2.343** | **0.393** | **0.0557** |

> **The literal one-loop full integral gives 𝔡(SU(3)) = 0.056 — BORDERLINE, ~180× below the note's O(10).**

Sanity: this `𝔡` is the same object across units (Λ cancels). It is not a peak artifact: `n` in Eq. (7) integrates the whole distribution; the mean size `ρ̄ = 0.393/Λ` sits right at the peak `ρ* = 0.368/Λ` (distribution is narrow in `ln ρ`, as expected).

**Sensitivity `[V-here]` (script §5).** 𝔡 scales **linearly** in the prefactor (the only knob), and mildly in the (unphysical) IR cutoff `ρ_max`:
- `ρ_max = 0.3, 0.5, 0.7, 1.0 /Λ` → `𝔡 = 0.002, 0.025, 0.052, 0.056` (the integral is dominated by `ρ ~ 1/Λ`, where one loop is least reliable — a caution, §5).
- `C = C₃, 0.01, 0.03, 0.1, 0.3` → `𝔡 = 0.056, 0.37, 1.10, 3.68, 11.05`.

So **O(1–10) is reached only if one inflates `C_N` from 1.5×10⁻³ to ≳ 3×10⁻²** — i.e. only via the wrong prefactor.

---

## §3 Literature cross-check: ILM and lattice

**Instanton-liquid model `[R]` (Schäfer–Shuryak, Rev. Mod. Phys. 70 (1998) 323; arXiv:hep-ph/9610451; Shuryak, Nucl. Phys. B 203 (1982) 93/116/140; Schäfer–Shuryak, Phys. Rev. D 53 (1996) 6522, arXiv:hep-ph/9509337).** The ILM fits to lattice correlators fix:

```
   n ≈ 1 fm⁻⁴  (quenched ~1.3–1.6 fm⁻⁴),     ρ ≈ 0.33 fm ,
   ⇒ 𝔡 = n ρ⁴ ≈ 1 × (0.33)⁴ = 0.012 .                                        (9)
```

The review itself quotes the packing fraction as `π² n̄ ρ̄⁴ ≈ 0.12`, i.e. **`n ρ⁴ ≈ 0.12/π² = 0.012`** — *exactly* the dilute number. (The factor `π²` is the 4-ball volume normalization; the 4D **volume** fraction is `n·(π²/2)ρ⁴ ≈ 0.06` — still ≪ 1.)

**Lattice topological charge `[R]`.** In a gas of unit-charge instantons and anti-instantons, `χ_top = ⟨Q²⟩/V = n_{I+Ā} = n` (independent Poisson). Direct measurements of pure SU(3):

- Teper, Phys. Lett. B 202 (1988) 553: `χ_t = (179 ± 4 MeV)⁴`.
- Smith & Teper, Phys. Rev. D 58 (1998) 014505: `χ_t^{1/4} ≈ 187 MeV`, `ρ̄ ≈ 0.5 fm`.
- de Forcrand, García Pérez, Hetrick, Stamatescu, arXiv:hep-lat/9802017: `χ_t ≈ (185 MeV)⁴`, `ρ̄ ≈ 0.5–0.6 fm`.
- Teper, "Topology in QCD", arXiv:hep-lat/9909124: `χ_t ≈ (200 MeV)⁴ ≈ 1 fm⁻⁴`.

Converting (`1 fm⁻¹ = 197.33 MeV`): `χ_t = (179–200 MeV)⁴ = 0.68–1.06 fm⁻⁴`, so **`n = χ_top ≈ 0.7–1.1 fm⁻⁴`** — *independently* confirming the ILM density with no semiclassical input.

**Packing fraction from the lattice-consistent numbers `[V-here]` (script §6):**

| n (fm⁻⁴) | ρ (fm) | 𝔡 = nρ⁴ | source |
|---|---|---|---|
| 0.68 | 0.33 | **0.008** | χ_t(179), ILM size |
| 1.0 | 0.33 | **0.012** | ILM central |
| 1.0 | 0.50 | **0.062** | χ_t, cooling-based size (Smith–Teper) |
| 1.06 | 0.50 | **0.066** | χ_t(200), cooling size |

Note the striking consistency of **scales**: the one-loop peak `ρ* = 0.368/Λ`; with `Λ ≈ 0.2 GeV`, `1/Λ = 0.99 fm`, so **`ρ* ≈ 0.36 fm` — matching the measured ILM `ρ ≈ 0.33 fm`.** The semiclassical formula gets the *size* right. Its *density* at that size is dilute (𝔡 ≈ 0.01–0.06), not dense.

---

## §4 Resolving the discrepancy — is the discrepancy real, and which convention matters?

**Yes, the 3-orders-of-magnitude gap is real and the note is on the wrong side of it.** The two "answers" in the task prompt are:
- one-loop peak, **as the note computes it** (wrong `C`, wrong object): 𝔡 = O(10) ✗;
- ILM fit: 𝔡 ≈ 0.01 ✓.

The resolution is **not** a convention ambiguity that could rescue O(10). It is a **computational error in the note**:
1. `C = O(10⁻¹)` vs true `C_3 = 1.5×10⁻³` — factor 66.
2. `D(ρ*)ρ*⁵` (dim 1/length) reported as `n ρ̄⁴` (dimensionless packing fraction) — the "29.6" is a *per-unit-size density at a point*, and the true full integral Eq. (7) gives 𝔡 = 0.056 with the right `C`.

**IS the one-loop prefactor reliable? Honest answer: NO, not to better than an order of magnitude, and it does not need to be.** `C_N` is convention- and scheme-dependent (the `0.466`, `1.679` carry scheme constants; the `[8π²/g²]^{2N}` exponent arrangement varies by reference). But **two independent routes agree on diluteness**: (i) the *literal* standard formula, full integral, true `C_3` → 𝔡 = 0.056; (ii) the *empirical* ILM+lattice → 𝔡 = 0.01. Both are ≪ 1. The honest range is **𝔡 ∈ [0.008, 0.07]**. To get 𝔡 ≳ 1 from the one-loop formula requires `C_N ≳ 3×10⁻²`, i.e. ~20× the true value — not a convention shift but an error. **The physical instanton vacuum is dilute; the "dense at the mass-generating scale" premise is false.** `[T-computable]` that the note's arithmetic is wrong; `[R]` that the vacuum is dilute (ILM/lattice); the *combination* is a robust DILUTE verdict.

**Caveat on "liquid" vs "gas" (honest).** Schäfer–Shuryak deliberately call it an *instanton liquid*, not a gas, because **interactions are O(1) and non-negligible** — but "liquid" here is a statement about *correlation strength*, **not** about packing fraction. The packing fraction is 𝔡 ≈ 0.01. A dilute-but-interacting system (like a dilute electrolyte) is fully consistent with 𝔡 ≪ 1.

**What Kotecký–Preiss actually requires `[R]` (Kotecký–Preiss, Commun. Math. Phys. 103 (1986) 491–498).** KP convergence is the summability `Σ_p |z_p| e^{a|p|} < ∞` over polymers `p`, with activity `z_p` and geometric size `|p|` (range), for some `a > 0`. **It is a statement about the activity and the interaction range, not about `n ρ⁴`.** For the instanton gas:
- the *activity* `z ~ e^{−S}` at the peak is `e^{−11} = 1.7×10⁻⁵` — genuinely small;
- but the *range* `|p| ~ ρ ~ 1/Λ ≈ 0.4–1 fm` is large and the mutual overlap at separation `~ρ` is O(1), so `e^{a|p|}` can grow faster than `z` decays for any `a` one actually needs;
- therefore **`𝔡 ≪ 1` is consistent with diluteness but does NOT by itself prove KP.** The KP obstruction at the physical scale — if there is one — lives in the **continuum range/activity control** at `ρ ~ 1/Λ`, which is the confinement problem itself, **not** in an O(10) density.

So: the note's *reason* (𝔡 = O(10)) is wrong; the note's *conclusion* (no rigorous cluster expansion at the physical scale) is a separate open question that this audit neither proves nor refutes. Correcting the premise forces the conclusion to be re-argued on different grounds.

---

## §5 Free-field / scale-invariant acid check `[V-here]` (script §7)

Remove the running: hold `g` fixed, so `S = 8π²/g²` is a constant and `D(ρ) = C_N S^{2N} e^{−S}/ρ⁵ ∝ ρ^{−5}`. Then:

- `D(λρ) = λ^{−5}D(ρ)` — **exactly scale invariant**; no preferred size, no `ρ̄`, no `Λ`.
- The density `n = ∫dρ/ρ⁵` **diverges** (at `ρ→0`) — no finite `n` without introducing an external scale.
- Hence `ρ̄` and `𝔡 = n ρ̄⁴` are **undefined**; there is no packing fraction and no mass.

**Acid test PASSED.** The diluteness statement (like the mass it is supposed to certify) **vanishes in the scale-invariant limit** — it exists only because the running coupling introduces `Λ` and a preferred `ρ*`. So the corrected diluteness number 𝔡 ≈ 0.01–0.06 is a *bona fide* dimensionally-transmuted statement, not a free-field artifact. ✅ This is consistent with `CLAUDE.md §1`'s acid-test rule: any new tool must fail on the free field; here the tool (semiclassical diluteness) correctly returns "undefined / no mass" when the running is switched off.

---

## §6 Verdict and what it does to the no-go

**VERDICT: DILUTE-to-BORDERLINE, NOT DENSE.** `𝔡 = n ρ̄⁴`:
- ILM fit (`[R]`): **0.012**;
- lattice χ_top (independent, `[R]`): **n ≈ 0.7–1.1 fm⁻⁴** ⇒ 𝔡 = **0.008–0.066**;
- literal one-loop full integral (`[V-here]`): **0.056**;
- honest range: **𝔡 ∈ [0.008, ~0.07]** — always ≪ 1.

**What survives of the note's §4.1:**
- ✅ The peak is at `ρ* = 0.368/Λ ≈ 0.37 fm` (large = IR; matches measured size). "Mass generates at `ρ ~ 1/Λ`" is **correct** `[V-here]`.
- ✅ `S* = 11` at the peak — marginal but adequate for one-loop semiclassics.
- ❌ `𝔡 = O(10)` — **FALSE**. Driven by `C = O(0.1)` (true `C_3 = 1.5×10⁻³`) and by using `D(ρ*)ρ*⁵` (a per-unit-size density, dim `1/length`) as if it were the dimensionless packing fraction. **The correct value is 𝔡 ≈ 0.056 (one-loop) / 0.012 (ILM) — the note is wrong by a factor ~200–800.**

**Effect on the `E2-FLIP-TRANSPORT` no-go `[NG]`:**

1. **As stated, the no-go's diluteness argument FAILS.** The premise "diluteness fails at the mass-generating scale because 𝔡 = O(10)" is numerically false; the scale is genuinely diluteness-compatible (𝔡 ~ 0.01–0.06). The `E2-FLIP-TRANSPORT §4.1` "two decisive facts" reduce to **one** true fact (large/IR peak) plus **one** false fact (dense).
2. **But the conclusion is NOT reopened as a theorem.** `𝔡 ≪ 1` is necessary, not sufficient, for Kotecký–Preiss: KP needs continuum **activity × range** control at `ρ ~ 1/Λ`, which is untouched by the (corrected) small density. Whether a convergent cluster expansion exists at the physical scale remains `[O]` — the same step the corpus identifies elsewhere as `H1`/`I2`, the confining-RG problem.
3. **Required correction.** `E2-FLIP-TRANSPORT §4.1–§4.2` and its headline claim "packing O(1) (SU(3): O(10))" must be **corrected**: replace with `𝔡 ≈ 0.01–0.06 (dilute)`, and re-derive the `[NG]` (if it is to stand) from the **range/overlap/continuum-control** obstruction — *not* from a spurious high density. The "size-modulus ⟺ wrong-scale" dichotomy is **weakened**: the modulus is still un-fixed (a real obstruction), but the gas is *not* parametrically dense, so "diluteness fails" cannot be the mechanism.

**Net:** the highest-value outcome is a **correction of a corpus error**. The physical instanton vacuum is dilute (𝔡 ≈ 0.01–0.06). The `E2-FLIP-TRANSPORT` no-go **in its stated form does not stand**; its conclusion may survive only on the different (and harder) ground of continuum cluster-expansion control. **This does not prove Clay, does not prove the dilute route closes, and does not prove any tool can exist.** The size modulus remains unfixed; the crossover remains open.

---

## §7 References

**Internal (this audit):**
- `D:\Dev\math\research\unconditional\E2-FLIP-TRANSPORT.md` — the audited note (§4.1 "n ρ*⁴ = O(10)", §4.2 obstruction, `[NG]`).
- `D:\Dev\math\research\unconditional\w7\transport_packing.py` — **the `[V-here]` witness** (one-loop D(ρ), exact `C_N`, peak, full-integral `n`, `ρ̄`, `𝔡`; ILM/lattice cross-check; free-field acid test). Run: `py -3.12 transport_packing.py` (single job).
- `D:\Dev\math\CLAUDE.md` (§1 acid-test / build-don't-stop; §2 labels; §5 wall).

**External `[R]` (as-reported; formula cross-checked against two sources; numbers NOT re-derived):**
- G. 't Hooft, "Computation of the quantum effects due to a four-dimensional pseudoparticle," *Phys. Rev. D* **14** (1976) 3432; "Symmetry breaking through Bell–Jackiw anomalies," *Phys. Rev. D* **14** (1976) 3432 (one-loop measure, prefactor `C_N`).
- C. Bernard, "Gauge zero modes, instanton determinants, and QCD," *Phys. Rev. D* **19** (1979) 3013 (SU(N) normalization).
- C. Callan, R. Dashen, D. Gross, "The structure of the gauge theory vacuum," *Phys. Rev. D* **17** (1978) 2717.
- T. Schäfer, E. Shuryak, "Instantons in QCD," *Rev. Mod. Phys.* **70** (1998) 323, arXiv:hep-ph/9610451 (Eq. 2.13: `C_N = 0.466 e^{−1.679N}/((N−1)!(N−2)!)`; ILM `n≈1 fm⁻⁴, ρ≈0.33 fm`; packing `π²n̄ρ̄⁴≈0.12`).
- E. Shuryak, *Nucl. Phys. B* **203** (1982) 93, 116, 140 (instanton liquid).
- T. Schäfer, E. Shuryak, "Interacting instanton liquid in QCD at zero and finite temperatures," *Phys. Rev. D* **53** (1996) 6522, arXiv:hep-ph/9509337.
- M. Teper, "The SU(3) topological susceptibility at zero and finite temperature," *Phys. Lett. B* **202** (1988) 553 (`χ_t=(179±4 MeV)⁴`); "Topology in QCD," arXiv:hep-lat/9909124 (`χ_t≈(200 MeV)⁴`).
- D. Smith, M. Teper, "Topological structure of the SU(3) vacuum," *Phys. Rev. D* **58** (1998) 014505 (arXiv:hep-lat/9802017 area), `χ_t^{1/4}≈187 MeV, ρ̄≈0.5 fm`.
- P. de Forcrand, M. García Pérez, J. Hetrick, I.-O. Stamatescu, arXiv:hep-lat/9802017, *Nucl. Phys. B Proc. Suppl.* (`χ_t≈(185 MeV)⁴` SU(3)).
- R. Kotecký, D. Preiss, "Cluster expansion for abstract polymer models," *Commun. Math. Phys.* **103** (1986) 19 (KP criterion `Σ_p |z_p|e^{a|p|}<∞`).

---

*Status: **AUDIT COMPLETE — CORPUS ERROR CORRECTED.** The `E2-FLIP-TRANSPORT §4.1` claim `𝔡 = nρ⁴ = O(10)` is **wrong** `[V-here]`: the SU(3) one-loop prefactor is `C_3 = 1.51×10⁻³` (not O(0.1)), and the packing fraction is the **full integral** `𝔡 = 0.056` (one-loop), `0.012` (ILM), `0.008–0.066` (lattice χ_top). The physical instanton vacuum is **DILUTE/ BORDERLINE** (`𝔡 ≈ 0.01–0.06`), **not dense**. The `[NG]`'s diluteness premise fails; its conclusion (no rigorous cluster expansion at the physical scale) may survive only via the separate, open continuum range/activity control (= `H1`/confinement), not via an O(10) density. **No QED; Clay stays open.***
