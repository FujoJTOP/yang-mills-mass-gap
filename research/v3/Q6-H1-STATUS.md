# Q6-H1-STATUS — the exact logical status of H1: **H1 is strictly weaker than the mass gap** (machine-checked non-converse); the three claimed equivalences (clustering / transfer gap / `G″≥λ`) are **false as biconditionals** and hold **iff the coarse marginal is Gaussian**; H1 is a theorem in seven regimes and **fails (correctly) on the free massless field**; proving H1 in 4D at `β=2.5` is **not circular** but leaves the atom as the one remaining input

> **2026-11-24 · Wave Q6 — the status question, not a new inequality.**
> H1 is `χ := Var_μ(ξ_{k→0}) = ∫s⁻¹dν ≤ 1/λ` uniformly in `(L,a)` at `β=2.5` (ch02 `eq:h1`).
> The campaign claims `Clay = H1 + atom`, with H1 necessary-not-sufficient and the atom
> `f₁² = ν({m₁²}) ≥ r₀` (with isolation `I₂`) the missing half. This note fixes the **exact**
> logical content: which equivalences are true, in which regime H1 is already a theorem, whether a
> proof of H1 would be circular, and the cleanest target + falsification battery.
>
> **Reads.** `P-WAVE-SYNTHESIS.md`, `STAR2-CLASSIFICATION.md`, `D6-DECISION.md`, `H6-HEAD-SPECTRAL.md`,
> `P6-PERIODIC-TO-H1.md`, `H1-RG-DESCENT.md`, `H2-DEFECT-BOUND.md`; `docs/paper-v2/ch02-squeeze.tex`
> (squeeze, H1 necessary-not-sufficient, thin-atom Prop., `§schur`), `ch03-existence.tex` (§L3,
> dimensionless correction), `app.tex` (`app:sharp`, `app:postaudit`); `lean/YMStatement/H1Sharpening.lean`;
> `C1-SZZ-NSCALING.md`, `C2-BALABAN-UNIFORM.md`, `V3-LIT-SYNTHESIS.md`;
> `unconditional/H1-IS-IT-THE-GAP.md`, `D3-SUSCEPTIBILITY-VARIATIONAL.md`, `REDUCTION-AUDIT.md`,
> `K0-WARD-LOWER.md`.
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness (cited, not re-run here). No QED.
> **No witness run** — every check quoted is already exact algebra or already machine-checked in the
> corpus (Lean `H1Sharpening.lean`; H6 witness §6·B double well; thin-atom Prop.). No third-party sources.

---

## 0 · Verdict (up front; no QED, Clay 0%)

1. **H1 is STRICTLY WEAKER than the mass gap `[NG]`.** Gap `⇒` H1 is a theorem (`[T]`, squeeze:
   `m₁²≥λ ⇒ χ=∫s⁻¹dν≤M₀/λ`). The converse is **false and machine-checked false**
   (`H1_does_not_imply_atom`, countermodel `(f,m,M,χ,λ)=(0,1,1,1,1)`; thin-atom Prop.
   `ν=wδ_ε+(1−w)δ_{s₂}`, `w=ε/2λ`: `χ=1/λ` yet `inf supp ν=ε→0`). So **a proof of H1 at `β=2.5`
   would NOT close Clay** — it is the necessary ceiling, not the edge.

2. **The three claimed equivalences are FALSE as exact biconditionals `[NG]`, and TRUE only in the
   Gaussian-marginal limit `[T]`.**
   - H1 vs **exponential clustering**: `clustering ⇒ H1` `[T]`; `H1 ⇒ clustering` `[NG]` (thin atom /
     massive `:φ²:`: `χ<∞`, `f₁²=0`, no isolated edge).
   - H1 vs **transfer-matrix gap in the k=0 sector**: `gap ⇒ H1` `[T]` (test function `φ=a`); `H1 ⇒ gap`
     `[NG]` (double well `V_eff=(g/2)(y²−1)²`: `Var=O(1)` but `λ₁∼e^{−g/2}→0`).
   - H1 vs **`G″_{k=0}≥λ`**: `ess inf G″≥λ ⇒ H1` `[T]` (Brascamp–Lieb); `H1 ⇒ G″≥λ` `[NG]` (non-Gaussian
     defect `η`; Savvidy flank `ess inf G″<0`).
   - **Unifying reason `[T]`:** H1 is the variance of the **single** coordinate `a=⟨v,X⟩`; clustering
     controls the **edge atom**, the transfer gap controls the **full inf-ratio over all test functions
     `φ`**, and `G″≥λ` controls the **essential infimum of the curvature**. Each is strictly richer than
     one coordinate's variance. They collapse to H1 **iff** `ν` is Gaussian (`G″·χ=1` exactly).

3. **H1 is a theorem in seven regimes and open in the eighth** (the physical 4D `β=2.5`): strong
   coupling `β<1/48` `[T]` (SZZ); `β=0` pure Haar `[T]` (`χ=3.7198` exact); 1D periodic chain `[T]`
   (1D exponential mixing); finite volume `[T]` (trivially, only with `λ=c_L→0`); 2D `[T-standard]`
   (Migdal solvability) with an abelian/KT caveat; large-`N` strong coupling `[T]` (`K_S∝N`); UV face
   (L3) `[C]` (true, certification missing — a different object). The free massless `U(1)` field
   **violates H1** (`χ=∞`) — the discriminating acid test. 4D `SU(N)` at `β=2.5` is `[O]` (H1 itself).

4. **Circularity: NOT circular.** H1 is strictly weaker than the gap, so proving it cannot assume the
   conclusion; but it also cannot deliver it. The sharpest separation: the feasible set
   `{χ≤1/λ, M₀ fixed}` contains, for every `ε>0`, a measure with `inf supp ν=ε` (thin atom) — so H1 fixes
   the **residue-weighted harmonic mean** and is silent on the **bottom of support**. Proving H1 leaves
   exactly one strictly-harder input: the atom floor `f₁²≥r₀` (with isolation `I₂`) — the corpus's `H*`.

5. **Carrier hazard (must be stated before "prove H1"):** on the gauge-invariant composite `tr F²`,
   `Var(q̄)=χ/V→0` `[V]`, so H1 is **vacuously true** and carries no gap content; on the gluon
   `a=⟨v,log U⟩` H1 holds but the carrier is gauge-dependent and fails the free-field acid test. Both
   branches are reported in `app:postaudit` / `H1-IS-IT-THE-GAP §5`. Hence **"H1 ⟺ gap" presupposes a
   gauge-invariant, non-self-averaging carrier whose linear existence is exactly the structural theorem
   that does not hold.** This is part of why H1 is weaker than the gap, and why the informative target is
   the atom, not H1.

6. **Cleanest target theorem + falsification battery** in §4. The campaign's minimal sufficient target is
   the **pair** `H* = (f₁²≥r₀ ∧ I₂) ∧ (χ≤1/λ)` on the correct carrier `[T given the pair]`; the irreducible
   single object is the **atom** on the `0⁺⁺` Källén–Lehmann measure. Falsification battery F1–F5 below.

---

## 1 · H1, stated unambiguously (the definition the rest of this note uses)

**The ceiling.**
```
H1:    χ := Var_μ(ξ_{k→0}) = ∫_0^∞ s^{-1} dν(s) ≤ 1/λ ,   uniformly in (L,a),  β=2.5.     (1.1)
```
`ξ_{k→0}` is the **total-normalised** coarse mode, `a^{(tot)}_k = (1/√V)Σ_x e^{ikx}O(x)`, for which
`Var(a^{(tot)}_k)=C̃(k)=∫e^{ikx}C(x)→χ` at `k→0` (H1-IS-IT-THE-GAP §1). The **mean**-normalised mode
`Õ(k)=(1/V)Σ` instead has `Var=χ/V→0`, making (1.1) trivially true and **vacuous** — so the
normalisation is load-bearing, not cosmetic.

**The four faces (all one object).** `G″(k=0) = V_eff″(0) = Γ″(0) = 1/χ = 1/C̃(0) = 1/∫s⁻¹dν`
`[T]`, the **inverse zero-momentum susceptibility**. It equals the pole mass `m₁²=inf supp ν` **iff**
`ν` is a single atom (Gaussian marginal). `G″·Var=1−η` with `η` the non-Gaussianity, **measured `2–5%`**,
not an identity (`[C]`; ch02 §G″ vs 1/χ, Ward `E[af]=1` is exact but `G″·Var=1` is Gaussian-only).

**The squeeze (the whole logical content).** For any positive `ν` on `[m₁²,∞)` `[T]`:
```
f₁²/χ  ≤  m₁²  ≤  M₀/χ ,      M₀=∫dν,   f₁²=ν({m₁²})=ν(edge atom).                          (1.2)
```
The **lower** rung is the atom contribution; the **upper** is the harmonic mean. Hence
`G″=1/χ` is an **upper** bound on `m₁²` (`1/χ ≥ m₁²/M₀` because harmonic mean `≥` min), and a
**gap** needs the **lower** rung with a **residue floor** `f₁²≥r₀`:
```
( f₁²≥r₀ ) ∧ ( χ≤1/λ )   ⇒   m₁² ≥ f₁²/χ ≥ r₀λ > 0 .                                        (1.3)
```
This is `H*`; it is `[T]` **given the pair**, and the **pair is the campaign's central content**.

---

## 2 · (1) Do the three claimed equivalences hold? — two are one-way, all die off the Gaussian

### 2.1 H1 vs exponential clustering of the k=0 channel — **FALSE as a biconditional**

"Exponential clustering of the k=0 channel" means `C(τ) ≤ C₀ e^{−mτ}` for a rate `m>0` bounded below,
equivalently (by Källén–Lehmann) an **isolated atom** of `ν` at `m₁²>0` with `f₁²>0` (`ν` has no support
in `(0,m₁²+δ)`).

- **`clustering ⇒ H1` `[T]`.** An isolated atom `f₁²δ_{m₁²}` + continuum `≥Λ★>m₁²` gives
  `χ = f₁²/m₁² + ∫_{≥Λ★}s⁻¹dν ≤ M₀/m₁² < ∞`. So a uniform mass gap forces `χ ≤ M₀/λ` (H1 up to `M₀`). ✓
- **`H1 ⇒ clustering` `[NG]`.** `χ = ∫s⁻¹dν < ∞` only says `C` is **integrable** (`∫₀^∞C dτ < ∞`), not
  that it decays **exponentially**. Two rigorous same-χ counterexamples already in the corpus:
  - **Thin atom** `ν = wδ_ε + (1−w)δ_{s₂}`, `w=ε/(2λ)` `[T]`: `χ=1/λ` (H1 holds) yet `inf supp ν=ε→0`;
    the mode `w e^{−ετ}` forces a correlation length `∼1/ε→∞` — no uniform rate. (ch02 Prop. thin atom;
    machine-checked non-converse.)
  - **Massive free composite `:φ²:`** `[T]`: gapped (`m₁²=4m²>0`), but `ν` is a **cut** on `[4m²,∞)`,
    `f₁²=0`; `χ_{:φ²:}` finite, so H1 holds, **yet there is no isolated atom** — the channel clusters only
    into a continuum (H1-IS-IT-THE-GAP §2.3-iii, the discriminating acid model).
  So H1 holds in a theory with no exponential clustering rate **bounded below** (thin atom) and in one with
  a gap but no atom (`:φ²:`). **H1 strictly weaker.** ∎

### 2.2 H1 vs a positive transfer-matrix gap in the k=0 sector — **FALSE as a biconditional**

The k=0 sector transfer matrix is the **coarse marginal transfer** `e^{−H_eff}` with ground state
`ν=|ψ₀|²`; its gap is the Poincaré (Rayleigh) constant `λ₁(ν)=inf_φ ℰ_ν(φ)/Var_ν(φ)`, `ℰ_ν=∫|∇φ|²dν`
(H6 §1 object (ii)).

- **`gap ⇒ H1` `[T]`.** Take the single test function `φ(a)=a` (the coordinate). Then
  `ℰ(a)=E_ν[|∂_a a|²]=1`, `Var(a)=χ`, so `λ₁ ≤ 1/χ`. If `λ₁ ≥ λ` then `χ ≤ 1/λ₁ ≤ 1/λ` — H1. ✓ (This is
  the **Rayleigh upper bound on the gap**, ch02 `G″≥λ ⇒ Var(a)≤1/λ` one way.)
- **`H1 ⇒ gap` `[NG]`.** `λ₁ ≤ 1/χ` is an **upper** bound on the gap; `χ ≤ 1/λ` gives `λ₁ ≤ 1/χ` and
  `1/χ ≥ λ`, which does **not** order `λ₁` against `λ`. Decisive counterexample `[T]+[V]`: the **double
  well** `V_eff=(g/2)(y²−1)²` (H6 witness §6·B): `Var=O(1)` (bounded, H1 holds with a fixed `λ`) while
  `λ₁ ∼ e^{−g/2} → 0` as `g→∞` (Kramers) — **bounded variance, vanishing transfer gap**. So H1 holds in
  every regime the double well is metastable, and the transfer gap collapses. **H1 strictly weaker.** ∎

  *The exact equivalence `λ₁(ν)=1/χ` holds iff `φ=a` is the minimiser iff `ν` is Gaussian* (then the
  ground state is Gaussian and the first excited mode is the coordinate with the same curvature).
  Non-Gaussian `⇒ λ₁ < 1/χ` strictly; this is the "faithful rename" of H6 reading (ii) being
  **`[C]` (Gaussian only)**, not `[T]` — a correction this note records.

### 2.3 H1 vs `G″_{k=0}≥λ` — **FALSE as a biconditional (holds only in the Gaussian limit)**

- **`ess inf G″≥λ ⇒ H1` `[T]`** (Brascamp–Lieb 1-D): `V_eff″≥λ` **a.e.** `⇒ Var(a)≤1/λ`. (This is the
  only correctly-signed *rigorous* route, D3 face (i).)
- **`G″≥λ ⇒ H1` `[NG]`.** A **point** value `G″(a₀*)≥λ` (typical/minimiser) is not an `ess inf`: a
  density sharp at `a₀*` but heavy-tailed has **both** large `G″(a₀*)` and large `Var`. The Savvidy
  flank makes `G″` negative on ~half the support (`H2 §0.3`; `H6 §3` `ess inf G″=m²<0`), so the BL
  hypothesis fails exactly where the floor is quoted. `[NG]` for the pointwise class.
- **`H1 ⇒ G″≥λ` `[NG]`.** `G″=(1−η)/χ`; `χ≤1/λ` gives `G″≥(1−η)λ` **only**, and for `η>0` (the
  measured `2–5%`) this can fall below `λ`. Exact equality `G″·χ=1` needs `η=0` (Gaussian). So H1 and
  `G″≥λ` differ by the non-Gaussianity defect everywhere `ν` is non-Gaussian.
- **Equivalence `[T]` iff Gaussian marginal.** Under `η=0`, `G″=1/χ` exactly, and all of H1, `G″≥λ`,
  and `λ₁≥λ` coincide. ∎

### 2.4 The unifying theorem `[T]/[NG]`

> **`[T]`** For a Gaussian coarse marginal `ν`, `G″ = 1/χ` exactly, and
> `H1 ⟺ G″_{k=0}≥λ ⟺ λ₁(ν)≥λ ⟺` (isolated atom at `m₁²=χ⁻¹`) exponential clustering — all one object.
> **`[NG]`** For a non-Gaussian `ν`, the chain splits into **strict one-way implications converging on
> H1**: `{clustering, transfer gap, ess-inf G″≥λ} ⇒ H1`, and each converse fails (thin atom / `:φ²:`;
> double well; Savvidy flank). **H1 is the projection of all richer objects onto the variance of the single
> coordinate `a`.**

This is the precise answer to (1): **none of the three equivalences is a theorem; each is a strict
implication into H1, and all collapse to equivalence exactly in the Gaussian-marginal limit** — which is
*not* the generic Yang–Mills case (every gauge-invariant carrier is a composite; `η≠0` measured).

---

## 3 · (2) Regimes where H1 is already a theorem — and the precise hypothesis of each

| regime | H1 status | precise hypothesis / object | source |
|---|---|---|---|
| **strong coupling `β<1/48`** | **`[T]`** | SU(N) 't Hooft-scaled Wilson action; Bakry–Émery `K_S=N(½−24|β|)>0` (radius `1/48` **N-independent**, SU(N)); volume-uniform LSI pushed through the block map (Lemma SC) ⇒ coarse Poincaré ⇒ `χ≤C/K_S`. **All `N`, all `L`, fixed `a`.** | C1 (SZZ CMP 400 (2023) 805), H6 §5C, STAR2 §5 |
| **`β=0` pure Haar** | **`[T]`** | i.i.d. Haar links; marginal = `L`-fold convolution; `χ=Var_Haar(ξ₃)=(4π²−6)/9=3.7198` **exactly**, `G″→9/(4π²−6)=0.26883`, `L`-independent; `Φ` convex ⇒ `G″≥N/12` pointwise ⇒ `χ≤E[1/G″]≤12/N`. **All `L`, all modes.** | D1/D2/D6 §3.1, D3 §2, H6 §5B |
| **1D periodic chain (`1×L`)** | **`[T]`** | 1D finite-range Gibbs measure at finite `β` ⇒ unique transfer-matrix ground state, exponential mixing ⇒ `Σ_d C(d)<∞` ⇒ `χ_L` bounded ⇒ H1. (Route A `1/L→0` is `[NG]` here.) | D6 §4.1, D1, P4 |
| **finite volume `L<∞`** | **`[T]` (trivially)** | finite-dimensional compactly-supported measure ⇒ `χ_L=Var(a)<∞`; but only with `λ(L)=c_L=4π²/L²→0`. **The uniformity in `L` is the open part; finite-`V` alone gives no uniform `λ`.** | D6 §1; ch03 §dim |
| **2D** | **`[T-standard]`** with caveat | 2D pure Yang–Mills is exactly solvable (Migdal 1975): area law + mass gap ⇒ `χ<∞` ⇒ H1. **Caveat:** for a 2D *abelian/scalar* proxy the constant-mode `χ` can diverge (KT/Mermin–Wagner Goldstone) — so "2D" is model-dependent, not a single theorem. No dedicated 2D campaign theorem beyond the exact 1D-transfer chain. | standard (Migdal); P6 §1.1 caveat |
| **large `N` (strong coupling)** | **`[T]`** | SZZ at fixed `β<1/48`: `K_S=N(½−24|β|)∝N` **improves** with `N`; per-colour `K_S/N` N-independent. **This is the strong-coupling 't Hooft limit `N→∞` at fixed `β<1/48`, NOT physical large-`N` at `λ≈1–3`.** | C1 §2–§3; V3-LIT §5C |
| **UV face / Balaban–Lemma B (L3)** | **`[C]`** (true, needs certification) | Deep UV `g_j→0`: Hessian `⪰0`, cone empty, so `Var≲E[H]` holds up to one-loop equipartition; **Lemma B `[T]`** closes the large-field entropy (volume-free per unit volume, const inequality `2cβ₀>d ln L`); Balaban gives action-level UV stability uniform in `a`, **fixed volume**. **But this is L3 (force, 1 derivative), NOT H1 (Hessian, `k→0`); Balaban does not reach the force Hessian, and is not volume-uniform.** | C2 §1–§2; ch03 §L3 |
| **free massless `U(1)`** | **`[NG]` — FAILS (correctly)** | `χ_L≍L²` (box gap), `G″≍4π²/L²→0`; `χ=∞`, so H1 **fails** — the discriminating acid test. A tool that "proves" H1 here is vacuous. | D5/D6 §6, H6 §5A |
| **4D `SU(N)` at `β=2.5`** | **`[O]` = H1 itself** | The open step: `χ` bounded uniformly in `(L,a)` at physical coupling. **Not a theorem anywhere.** | ch02 `eq:h1`; STAR2 §3 |

**Honest summary of the UV row.** The campaign does **not** have H1 at the UV end as a theorem — it has the
UV face (L3) as `[C]`: the *inequality* is essentially true (positive Hessian, cone empty) but the
*uniform non-perturbative certification* of the **force** variance is open, and force ≠ Hessian ≠ `χ`.
So "H1 is a theorem in the UV" is **not** correct as stated; the correct statement is "the UV face is
`[C]`, on a different derivative order, and its entropic half is `[T]` (Lemma B)." This note records
the correction (charter §2: no upgrade of `[C]` to `[T]`).

---

## 4 · (3) Circularity check: does a proof of H1 at `β=2.5` give the gap? — **NO; strictly weaker, not circular**

### 4.1 The two directions, labelled

- **`gap ⇒ H1` `[T]`** (necessity). `m₁²≥λ ⇒ χ=∫s⁻¹dν≤M₀/λ`. (Squeeze upper rung.) This is the only
  direction the corpus proves, and it is exact.
- **`H1 ⇒ gap` `[NG]`** (non-converse). Machine-checked: `H1_does_not_imply_atom` (Lean, zero `sorry`,
  axioms `[propext, Classical.choice, Quot.sound]`) — the reduction's hypotheses
  `{0<λ,0<χ,χ≤1/λ,0<m,f≤χm,χm≤M}` admit the model `(f,m,M,χ,λ)=(0,1,1,1,1)` with `f=0` (no atom).
  Analytic witness: thin-atom Prop. `ν=wδ_ε+(1−w)δ_{s₂}`, `w=ε/2λ` ⇒ `χ=1/λ`, `inf supp ν=ε→0` (ch02).

### 4.2 The sharpest separation statement

> **`[T]+[NG]` (sharp separation).** Let `𝒮_λ={ν≥0 : χ(ν)=∫s⁻¹dν≤1/λ, M₀(ν)=M₀ fixed}`. Then
> ```
>   inf_{ν∈𝒮_λ}  inf supp ν  =  0 ,        but       {ν : inf supp ν ≥ √λ}  ⊊  𝒮_λ   strictly.
> ```
> Equivalently, H1 controls the **residue-weighted harmonic mean** `M₀/χ` (an **upper** rung `m₁²≤M₀/χ`)
> and is **silent on the bottom of support**: the deficiency between H1 and the gap is **exactly the edge
> residue `f₁²`** (and its isolation `I₂`), because the lower rung `m₁²≥f₁²/χ` degenerates to `m₁²≥0` as
> `f₁²→0`. **H1 and the gap differ by one datum `f₁²∈[0,M₀]`, and nothing in H1 bounds `f₁²` below**
> (every moment/positivity/convex-order functional is Schur-concave — an **upper** bound on `f₁²`, never a
> lower bound; ch02 `§schur`; the free field is a feasible atomless point).

### 4.3 Is it circular? — **No; and the distinction is load-bearing**

- **Not circular.** H1 is strictly weaker than the gap, hence a proof of H1 does **not** presuppose the
  conclusion. It is a genuine **necessary** partial result: it would establish the *ceiling* face of the
  pair `H*` and would fix the **susceptibility** (screening length / inverse zero-momentum curvature).
- **But also not sufficient.** It leaves the atom floor `f₁²≥r₀` + isolation — the corpus's `H*` and
  "the hard half." The logical shape is a **diamond converging on H1**, with `H*` the sufficient apex
  (REDUCTION-AUDIT §"what survives"; ch02 `§verdict`). Written as `⟺`, the reduction is **unsound**;
  H1 alone is one conjunct of the sufficient pair.

### 4.4 The conjecture that would change this — and its current status

The corpus *once* claimed `H*=H1` as a structural identification (`G2-ASSEMBLE-100 §7`). That would make
H1 **equivalent** to the gap (and a proof of H1 would close Clay). **But the corpus has downgraded it to
a doubted conjecture** (`ch02 §verdict`: "now doubted"; H1 is the inverse susceptibility, equal to the
pole only for elementary single-atom fields; the free massive `:φ²:` has `G″=1/χ≠m₁²`). Under the corpus's
**current** stance the identification is **false in general**, so the safe and honest logical position is:

> **`[NG]` H1 is STRICTLY WEAKER than the mass gap.** A proof of H1 at `β=2.5` is neither circular nor
> Clay-closing; it proves the **necessary** ceiling `χ≤1/λ` and leaves the **atom** as the single
> remaining, strictly-harder input. (Only the unproven `H*=H1` identification would upgrade it.)

### 4.5 Carrier caveat (why "prove H1" is under-specified)

Before any of 4.1–4.4 can be run, the **object** must be pinned:
- **Composite `tr F²` (the natural gauge-invariant `0⁺⁺` carrier):** `Var(q̄)∝V^{−1.0}` `[V]`, so
  `Var≤1/λ` holds for **every** `λ` at large `V` — **H1 vacuous**, no gap content (`app:postaudit`;
  H1-IS-IT-THE-GAP §5).
- **Gluon `a=⟨v,log U⟩`:** `Var` has the right behaviour, but it is **gauge-dependent**, sits *below* the
  free field at every measured `β` (`G″/free=0.18–0.63`), and **fails the free-field acid test**.
- **No third option** (no linear gauge-invariant carrier of `Var` with `Var↛0`; H1-IS-IT-THE-GAP §5
  structural theorem). So the meaningful H1 lives on the **total-normalised propagator branch**, whose
  gauge-invariant carrier is exactly the structural object that does not exist. **This is independent
  confirmation that H1 is weaker than (and partially orthogonal to) the gauge-invariant gap.**

**Conclusion of §4.** H1 is **strictly weaker** than the gap, **not circular**, and **not sufficient**;
its content under the natural carrier is either vacuous (composite) or gauge-dependent (gluon). Proving
it would relocate the problem to the atom — it does not remove it.

---

## 5 · (4) The cleanest target theorem and the falsification battery

### 5.1 Cleanest target theorem (the campaign's minimal sufficient statement)

> **`(T★)` [the cleanest sufficient target].** There exist constants `λ>0`, `r₀>0`, `δ>0` and a
> **gauge-invariant `0⁺⁺` carrier `O`** (total-normalised `k→0` mode, `Var↛0`) such that, uniformly in
> `(L,a)` at `β=2.5`:
> ```
>   (C1)  χ = ∫s⁻¹dν ≤ 1/λ                         [the ceiling; H1 on the correct carrier]
>   (C2)  f₁² = ν({m₁²}) ≥ r₀   and   Λ★ − m₁² ≥ δ   [the edge atom, with isolation I₂]
> ```
> Then `m₁² ≥ f₁²/χ ≥ r₀λ > 0` — the mass gap `[T given (C1),(C2)]` (ch02 `eq:sufficient`/`eq:hstar`).

**Three equivalent single-statement forms** (the campaign should pick one and commit):
1. **Pair form (above).** Most faithful to the logic `Clay = H1 + atom`.
2. **Conditioned-variance form (`app:sharp`, the sharpest no-slack single inequality):** there is
   `Δ₀>0` with, uniformly in `(L,a)` and every coarse transverse `v` at `k→0`,
   ```
   Var_μ(⟨v,∇V⟩ | ξ=0)  ≤  E_μ[H_vv | ξ=0] − a²Δ₀² ,      V = S_W + Φ ,        (5.1)
   ```
   equivalently `D ≤ S − λ` **together with** `S ≥ λ`, equivalently the dimensionless ratio
   `G″_a/a² ≥ Δ₀²` (dimension-correct under `Z_O` transport; a constant `λ` is a lattice artifact —
   ch03 §dim). **This is the gap in curvature coordinates; proving it closes Clay, but it is not weaker
   than the gap** — it is the target, not a half-step.
3. **Atom form (the irreducible single object):** `f₁² = ν({m₁²}) ≥ r₀` with isolation `I₂` on the KL
   measure of the gauge-invariant `0⁺⁺` operator. Given the carrier vacuity of (C1) on the composite, this
   is the **operative** content: on the natural carrier H1 is free, so **the entire Clay gap reduces to the
   atom** `[O]` (`H1-IS-IT-THE-GAP §6`, `K0-WARD-LOWER §5`).

**Recommended cleanest target.** State `(T★)` as the pair, and pursue the **atom form** as the operative
lemma: it is minimal, carrier-unambiguous, and the one unprotected number. H1 is the *matching ceiling*,
not the target.

### 5.2 Falsification battery for any proposed proof of H1 (or of `(T★)`)

A proposed proof must survive all five; failing any **refutes** it (each is an exact/machine-checked
countermodel in the corpus, so no new witness is needed):

- **F1 — free massless field (must FAIL, i.e. must not manufacture a floor) `[T]+[V]`.** Free massless
  `U(1)`/Gaussian: `χ=∞`, `G″=0` (`λ₁=G″=1/Var→0`). If the proof "establishes H1" on the free massless
  field (finite `χ`), it is **vacuous** and refuted. (Acid test; `H6 §5A`, `D6 §6.1`.)
- **F2 — thin atom (must NOT certify the gap from `χ` alone) `[T]` + Lean.** Feed the mechanism
  `ν=wδ_ε+(1−w)δ_{s₂}`, `w=ε/2λ`: it has `χ=1/λ` (H1) but `m₁=ε→0`. If the proof claims "H1 alone ⇒ gap,"
  the thin atom refutes it. (ch02 Prop.; `H1_does_not_imply_atom`.) **This is the decisive falsifier for
  any claim that H1 closes Clay.**
- **F3 — deep double well (must not claim a barrier-uniform `λ`) `[T]+[V]`.** On
  `V_eff=(g/2)(y²−1)²`, `Var=O(1)` but `λ₁∼e^{−g/2}→0` as `g→∞`; the Savvidy normal form is exactly this
  (H6 §3 double well). A proof claiming a **uniform-positive** transfer gap (or `ess inf G″≥λ`) at fixed
  positive typical curvature is refuted as `g→∞`.
- **F4 — carrier / vacuity (must be non-vacuous) `[V]`.** On the gauge-invariant composite `tr F²`,
  `Var(q̄)=χ/V→0`, so `Var≤1/λ` is free. A proof that only works because `Var→0` (composite) proves nothing;
  it must be checked on the **total-normalised** mode and on a carrier with `Var↛0`. If the only such
  carrier is gauge-dependent, the proof has not reached the gauge-invariant gap. (`app:postaudit`;
  H1-IS-IT-THE-GAP §5 structural theorem.)
- **F5 — dimension / `a→0` (must be a ratio) `[T]`.** `G″_a=a²m_phys²`; a constant `λ` (not
  `λ_a∼a²Δ₀²`) would give `m_phys=√λ/a→∞`, a UV lattice artifact, not a particle. A proof whose bound
  does not scale as `G″_a/a²` (or is not multiplicatively `Z_O²`-transported) is dimensionally wrong
  (ch03 §dim).

**What a valid proof would therefore look like.** A gauge-invariant, total-normalised, `(L,a)`-uniform
bound that (i) **returns `∞`/`0` on the free massless field** (F1), (ii) **does not derive the gap from
`χ` alone** (F2) but supplies the atom or the conditioned-variance `G″_a/a²≥Δ₀²` (5.1), (iii) **survives
the Savvidy double well uniformly** (F3) — which forces it to beat the *essential infimum*, not a typical
point (H6 §2–§3, the precise obstruction), and (iv) is dimensionally a ratio (F5). No such object exists
in the present toolbox (H6 §4 `[NG]` for the spectral/PF/BE/Cheeger/block class; STAR2 §5 pincer).

---

## 6 · Status ledger

| item | statement | label |
|---|---|---|
| H1 definition | `χ=Var_μ(ξ_{k→0})=∫s⁻¹dν≤1/λ`, total-normalised mode, uniform `(L,a)`, `β=2.5` | ch02 `eq:h1` |
| `G″(k=0)=1/χ=V_eff″=Γ″=1/C̃(0)=1/∫s⁻¹dν` | inverse zero-momentum susceptibility | **[T]** |
| `G″·χ=1` | Gaussian marginal only; measured `2–5%` defect | **[T] (Gaussian)**, not identity |
| squeeze `f₁²/χ≤m₁²≤M₀/χ` | any positive measure | **[T]** |
| `gap ⇒ H1` (`m₁²≥λ ⇒ χ≤M₀/λ`) | necessity | **[T]** |
| `H1 ⇒ gap` | non-converse (thin atom; Lean countermodel) | **[NG]** |
| `clustering ⇒ H1` / `H1 ⇒ clustering` | one-way / false (thin atom, `:φ²:`) | **[T] / [NG]** |
| `transfer gap ⇒ H1` / `H1 ⇒ transfer gap` | one-way (Rayleigh `φ=a`) / false (double well) | **[T] / [NG]** |
| `ess inf G″≥λ ⇒ H1` / `H1 ⇒ G″≥λ` | one-way (Brascamp–Lieb) / false (η>0, Savvidy) | **[T] / [NG]** |
| all four coincide | **iff Gaussian marginal** | **[T]** (Gaussian) / **[NG]** (general) |
| H1 at `β<1/48` | SZZ `K_S=N(½−24|β|)>0`, all N,L, fixed a | **[T]** |
| H1 at `β=0` | `χ=3.7198` exact, `G″→0.26883`, `L`-indep. | **[T]+[V]** |
| H1 on `1×L` chain | 1D exponential mixing | **[T]** |
| H1 finite volume | trivial `χ_L<∞`, `λ=c_L→0`; uniform-in-`L` open | **[T] (non-uniform) / [O] (uniform)** |
| H1 in 2D | Migdal solvability (gapped); abelian/KT caveat | **[T-standard]** |
| H1 large N (strong coupling) | `K_S∝N`, per-colour N-indep.; not physical large-N | **[T] / [O]** |
| H1 at UV end (L3/Balaban) | true modulo certification; force≠Hessian; entropy half Lemma B | **[C] (L3) / [T] (entropy)** |
| H1 on free massless `U(1)` | `χ=∞`, fails (acid test) | **[NG] (correct)** |
| H1 4D `SU(N)` `β=2.5` | the open step | **[O]** |
| carrier: composite `tr F²` | `Var(q̄)∝V⁻¹`, H1 vacuous | **[V]** |
| carrier: gluon | gauge-dependent, acid-fails | **[V]/[O]** |
| `H*=H1` (one object) | doubted conjecture; if true, H1⟺gap | **[C] (doubted)** |
| `(T★)` pair `(f₁²≥r₀,I₂)∧(χ≤1/λ) ⇒ m₁²≥r₀λ` | minimal sufficient | **[T] given the pair; [O] to prove** |
| Clay gap | unbuilt | **[O]** · **0%** |

---

## 7 · Sharpest line

**H1 (`χ=Var_μ(ξ_{k→0})≤1/λ`) is the inverse zero-momentum susceptibility — the variance of the single
`k=0` coordinate — and it is *strictly weaker* than the mass gap.** The gap implies it (`[T]`, squeeze);
it does not imply the gap (`[NG]`, machine-checked: thin atom `ν=wδ_ε+(1−w)δ_{s₂}` with `χ=1/λ`,
`m₁=ε→0`). The three claimed equivalences — to exponential clustering, to a transfer-matrix gap, to
`G″≥λ` — are **false as biconditionals**: each is a strict implication *into* H1 (`clustering`/`transfer
gap`/`ess inf G″≥λ` `⇒` H1) and the converses fail on the thin atom / the double well / the Savvidy
flank; **all collapse to equivalence iff the coarse marginal is Gaussian**, which the generic composite
carrier is not (`G″·χ=1−η`, `η≈2–5%`). H1 is a theorem at `β<1/48` (SZZ), `β=0` (exact Haar), the 1D
chain, finite volume (non-uniformly), 2D (solvability), and large-`N` strong coupling; the UV face is
`[C]` on a *different* object (force, L3); the free massless `U(1)` field correctly **fails** it; and 4D
`SU(N)` at `β=2.5` is `[O]`. Proving H1 in 4D is **not circular** — it is strictly weaker than the
conclusion — but it is **not sufficient**: it proves the ceiling and leaves the **atom** `f₁²≥r₀` (with
isolation) as the one remaining, strictly-harder input, because H1 fixes a residue-weighted harmonic mean
and is blind to the bottom of support (convex-order / Schur-concavity). On the natural gauge-invariant
carrier the ceiling is even vacuous (`Var(q̄)∝V⁻¹`), so the Clay gap lives in the **atom**, not in H1. The
cleanest sufficient target is the pair `(T★) = (atom+isolation) ∧ (ceiling on the correct carrier)`; any
proposed proof must pass the falsification battery F1–F5 (fail on the free field, not derive the gap from
`χ` alone, survive the deep double well uniformly, be non-vacuous on a total-normalised carrier, and be
dimensionally a ratio). **The wall is unmoved and now unambiguously located: H1 is the necessary ceiling;
the gap is the edge atom.** **No QED. Clay 0%.**

---

*Refs.* **The squeeze / H1 not sufficient:** `docs/paper-v2/ch02-squeeze.tex` (`eq:chi`, `eq:h1`,
`eq:squeeze`, thin-atom Prop., `§h1-notsuff`, `§schur`, `eq:sufficient`, `eq:hstar`, `§verdict`);
`lean/YMStatement/H1Sharpening.lean` (`lemmaG_twoSided`, `cert_necessity`, `cert_sufficiency`,
`H1_blind_to_thin_atom`, `H1_does_not_imply_atom`; 0 `sorry`, axioms `[propext, Classical.choice,
Quot.sound]`). **Susceptibility ≠ pole:** `unconditional/H1-IS-IT-THE-GAP.md` (§2.3 `:φ²:`, §3 direction,
§5 carrier); `unconditional/D3-SUSCEPTIBILITY-VARIATIONAL.md` (BL/BL-defect/Legendre/KL faces, `β=0
χ≤12/N`, 2.5 log-concavity fails); `unconditional/K0-WARD-LOWER.md` (`E[af]=1`, tower = `p″(0)<0` = local
H1, no residue floor); `unconditional/REDUCTION-AUDIT.md` (diamond, one-way, false `⟺`). **Spectral /
equivalence failures:** `research/v3/H6-HEAD-SPECTRAL.md` (coarse Poincaré = H1 [Gaussian]; typical point
vs `ess inf`; double well `λ₁∼e^{−ΔV}`; spectral no-go `[NG]`); `H2-DEFECT-BOUND.md` (rank-one split,
Fisher identity `E_ρ[G″]=I(ρ)`, `ess inf G″<0` at head). **Dichotomy / decompactification:**
`D6-DECISION.md` (`χ<∞⇒saturation`, `m₁>0⇒χ<∞` necessary-not-sufficient, `β=0` exact);
`P-WAVE-SYNTHESIS.md`, `P6-PERIODIC-TO-H1.md` (periodic floor = `k=0` ceiling, never atom). **Regimes:**
`C1-SZZ-NSCALING.md` (SZZ `1/48`, `K_S=N(½−24|β|)`); `C2-BALABAN-UNIFORM.md` (Balaban action-level, Lemma B
entropy `[T]`, L3 `[C]`); `V3-LIT-SYNTHESIS.md` (LIT-GAP: strong coupling + finite-group weak coupling,
no continuous weak coupling). **Dimension / tower:** `docs/paper-v2/ch03-existence.tex` (`eq:l3`,
`eq:dim`), `app.tex` (`app:sharp` `(5.1)`, `app:postaudit` carrier audit); `research/v3/H1-RG-DESCENT.md`
(soft-mode failure, curvature ≠ atom). **Standard inputs (textbook, not re-fetched):** Källén 1952;
Lehmann 1954; Osterwalder–Schrader 1973–75; Brascamp–Lieb, Adv. Math. **22** (1976) 366; Prékopa,
Acta Sci. Math. **32** (1971) 301 / **34** (1973) 335; Bakry–Émery 1985; Shen–Zhu–Zhu, CMP **400** (2023)
805, arXiv:2204.12737; Migdal 1975 (2D YM solvability). **No third-party sources downloaded; no witness
re-run** — every quantitative claim quoted is exact algebra, already-Lean-checked, or an H6/D1/D6 witness
cited as such. No `[C]` upgraded to `[T]`; no QED; Clay 0%.
