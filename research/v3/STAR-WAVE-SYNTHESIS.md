# STAR-WAVE-SYNTHESIS — (★) survives the boundary falsification, is exact in the one-plaquette model, and decomposes into a finite list: good-set log-concavity, Schur loss, S≥δ (=H1), and H2-global for (★dec)

> **VERDICT (no QED, Clay 0%).** Six angles attacked the single live inequality
> `(★) Var_{μ(·|U′)}(⟨Pᵀ∇S_j,v⟩) ≤ E⟨PᵀH_jP v,v⟩ + W_j(v) − δ_W`, uniform in `j,U′,a`, plus its
> cross-shell companion `(★dec)`. **`(★)` was not proved and not falsified.** What the wave did:
> (i) it killed the most natural counterexample (boundary aligned with a Savvidy negative mode) in
> an exactly solvable model, where `(★)` reduces to `2W_geo ≥ δ_W` with `W_geo>0` `[T]`;
> (ii) it showed that the **strong form "uniform in all `U′`" is false** — uniform log-concavity
> fails on chromomagnetic boundaries — so the correct statement needs the small-field set
> (excision), an average over `U′`, or a boundary-dependent `δ_W`;
> (iii) it proved that `(★)` and `(★dec)` are **independent** and must be assumed separately;
> (iv) it removed the `N`- and `a`-obstructions (linear `N`, `a`-uniform), and (v) it showed that
> even a perfect conditional Poincaré controls the defect `D` but **not** the coarse `S`, so the UV
> step still contains the IR condition `S≥δ`. Clay remains 0%.

---

## §0 · Labels

`[T]` theorem · `[C]` conditional · `[NG]` no-go for a named class · `[O]` exact missing step ·
`[V]` numerical. All six witnesses run under `py -3.12` (numpy 2.4.6), exit 0.

---

## §1 · The six verdicts

| Angle | Verdict | Load-bearing result |
|---|---|---|
| **S1** conditional Poincaré | `[T]`+`[NG]`+`[O]` | On the small-field boundary set `𝒮_t`, the modulus `ρ(U′) ≥ min{c_□+N/12−C₁βF_∂(U′), N/12+R₀}` gives a uniform LSI/BL `[T]` (corollary `[C]`). Uniform log-concavity for **all** `U′` fails: any chromomagnetic boundary makes `K=Π₀HessVΠ₀` indefinite on the soft mode (witness `λ_min: 0→−1.43→…→−8`) `[NG]` for the BE/BL/PI/Prékopa class. Even a perfect conditional PI controls `D`, **not** `S`; the free zero mode has `K⪰0` yet `S=1/Var(a)→0`, so `(★)` still needs `S≥δ` `[O]`. |
| **S2** boundary counterexample | `[T]`+`[V]` (falsification failed) | Exact SU(2) one-plaquette/one-link model with boundary holonomy `B=e^{iθσ₃/2}`. The flat-chart Stein identity gives `Var(G_v+Φ′_v)=E[H_vv+Φ″_vv]`, hence `Cov(G_v,Φ′_v)=W_geo(v)` and **`G″_model(v)=2W_geo(v)` exactly** (residual ≤5e-7). Thus `(★) ⟺ 2W_geo(v|U′)≥δ_W`. `W_geo>0` for all `β>0`, all `θ∈[0,π]`, both `v`; the worst case (max holonomy `θ=π`, transverse `v=e₁`) gives `δ_W≈0.132` at `β=2.5`, growing to `≈1/3` along the UV tower. Acid test: flat/free group and `β=0` give `G″=0`, so `(★)` fails there `[T]`. |
| **S3** UV Feshbach | `[C]` reduction | Exact two-block split `G″=E[Σ]−D`, `Σ=A−⟨K⁻¹c,c⟩`; if `K⪰γ_j` then Prékopa+BL give `D≤0`, so **`(★) ⟸ [fibre gap] ∧ [E[Σ]≥δ_W]`** `[T]`. `W_j=W_geo+X`; its only provable content `Φ″⪰N/12` must pay the Schur loss `⟨K⁻¹c,c⟩` and **cannot absorb `X`**. TCI is the same Schur algebra but the opposite sign structure (TCI retains a tachyonic `A≺0`; `(★)` has a positive coarse mode), confirming "same form, different object". |
| **S4** decimation | `[T]` independence | `(★)` and `(★dec)` are the diagonal and off-diagonal of one covariance matrix. `(★)⇏(★dec)` (diagonal bound gives only `|Cov|≤const`; power-law counterexample) and `(★dec)⇏(★)` (boundedness ≠ curvature; Savvidy-negative `E[H]`). The conjunction gives `G″⪰½δ_W` by Gershgorin with interface `κ>(2c_dC/δ_W)^{1/d}` (`c₄=12π²`; witness sign-flip `κ*≈2.52` in `d=2`). The only bridge `(★)⇒(★dec)` needs uniform global fiber LSI (`H2-global`, `C(β)<1`), which fails at the head (`C≈15`) `[C]`. |
| **S5** `N/a` uniformity | `[T]`+`[C]` | The block-conditional `(★)` preserves **linear `N`** (the adjoint dimension is multiplicity, not magnitude; SZZ Lemma 4.1) and **`a`-uniformity**. `δ_W=cN−O(1)`, `c≈1/12–0.13`, from the Haar reservoir `E[Φ'']≥N/12` and `W_geo≈0.125N`. Conditioning adds no factor; the only exponential is the measure of the large-field boundary set, `e^{−cNj}`, which **shrinks with `N`**. `[NG]` only for "all `U′`, no excision" uniform conditional LSI. SZZ is off-regime in the tower (`β_{SZZ,j}=11j/(12π²)>1/48`), so the block constant comes independently from `W_geo`. |
| **S6** Dobrushin/decoupling | `[T]` at strong coupling + `[NG]`+`[O]` | Dobrushin coefficient `α(β_j)=κβ_j`, `κ=O(d)≈72`; `α<1` gives `(★)∧(★dec)` via Dobrushin–Shlosman complete analyticity `[T]`, but only for `β_j<β_Dobr≈1/72`. The strong-coupling window `[0,1/72]` and the physical tower `[2.5,∞)` are **disjoint** — it proves `(★)` for zero tower shells. At `β≈2.5`, `α≈180`, percolation is supercritical `[NG]` for the local-influence class. At the deep UV end, BL gives `Var≤E[H]` for log-concave fibres, reducing `(★)` to `W_j≥δ_W` `[C]`. Normalization-honest: SZZ's `β<1/48` maps to standard Wilson `β_W≈0.19` (`λ_Hooft≈96`), so the honest head-to-threshold distance is `≈13×`, not `≈120×`. |

---

## §2 · The corrected statement of `(★)`

The wave's main structural output is that the strong form, uniform in **all** `U′`, is false
(S1: chromomagnetic boundaries; S5: boundary coexistence). The corrected forms, in increasing
strength:

- **`(★_sf)` (small-field boundary / excision).** On the set `𝒮_t` of boundaries with
  `F_∂(U′)≤t`, `Var_{μ(·|U′)}(⟨Pᵀ∇S_j,v⟩) ≤ E⟨PᵀH_jP v,v⟩+W_j(v)−δ_W`, with `δ_W>0` uniform in
  `U′∈𝒮_t`, `v`, `a`, and depending on `β_j`. This is the form compatible with Balaban's excision.
- **`(★_avg)` (averaged boundary).** The same inequality after averaging over `U′` with the
  equilibrium boundary measure; weaker, but enough if the RG map averages.
- **`(★_{δ(U′)})` (boundary-dependent margin).** Allow `δ_W(U′)>0` with a lower bound only on
  `𝒮_t`; `δ_W` may vanish on the excised set.

`(★dec)` is separate: `|Cov_μ(s_B,s_B′)| ≤ C e^{−κ dist(B,B′)/a}`, tempered per unit volume, with
`κ>(2c_dC/δ_W)^{1/d}` for the Gershgorin conjunction. Its bridge from `(★)` needs `H2-global`
(`C(β)<1`), which fails at the head.

**The exact toy-model theorem (S2).** In the one-plaquette/one-link SU(2) model,
`G″=2W_geo` exactly, so `(★_sf)` is *equivalent* to `2W_geo≥δ_W`, and `W_geo>0` gives it with an
explicit constant. The remaining hiding place is the **interior coarse-projection boundary
term**, not boundary alignment with the negative mode.

---

## §3 · What was removed

- **"Boundary aligned with Savvidy kills `(★)`"** — refuted in the exact model (S2): the worst
  case is the targeted one and still gives `δ_W≈0.132` at `β=2.5`.
- **"Volume growth is intrinsic"** — already refuted by Lemma B (C2); S5 confirms the boundary
  measure shrinks with `N`.
- **"`N`-degradation"** — removed: linear `N` survives blocking (S5).
- **"`a`-collapse"** — `δ_W=Θ(a⁰)`, `j*~log(1/a)` (S5).
- **"TCI closes L3"** — refuted (S3, C2): same form, different object.
- **"`(★)` implies `(★dec)`"** — false in general (S4); they are independent.

## §4 · What remains (the finite list)

1. **Good-set log-concavity** `K⪰0` on `𝒮_t` — `[T]` on small-field boundaries (S1), `[NG]` for all
   `U′`; the excision/measure estimate is the standard Balaban step.
2. **Schur loss** `E[Σ]≥δ_W` — `[O]` (S3); the geometric term `Φ″⪰N/12` pays part of it, `X`
   remains.
3. **Coarse `S≥δ`** — `[O]` and **identical to H1** (S1). This is the crucial structural point: the
   UV step is not closed by UV analysis alone; the IR condition reappears inside it.
4. **`(★dec)` with `H2-global`** — `[C]`, fails at the head (S4).
5. **Interior coarse-projection term** — the only surviving hiding place in the exact model (S2).
6. **Strong-coupling-only theorem** — S6's Dobrushin result is real but disjoint from `β≈2.5`.

## §5 · Honest bottom line

The `(★)` wave did not prove the gap. It did something more useful than another failed route: it
turned a black-box inequality into a **finite, labelled system** with
- one exact theorem in a solvable model (S2),
- one theorem on the good set (S1),
- one exact decomposition (S3),
- one independence/conjunction theorem (S4),
- `N/a` bookkeeping that removes two potential obstructions (S5),
- one strong-coupling theorem that is honest about being disjoint from the physical point (S6).

The remaining list is short and specific. The most striking item is (3): even a perfect
conditional Poincaré inequality controls only the defect `D`, not the coarse `S`, and `S≥δ` is
H1. So "one form, two instances" is now precise at the level of the **condition**: both faces need
the same coarse floor. Proving the UV instance does not give the IR instance for free; it gives
the same inequality again.

**Clay = 0%. No QED.** The wall is still there; the UV side of it now has a map with six labelled
rooms, one of which is the IR wall itself.

*Refs:* `S1-CONDITIONAL-POINCARE.md`, `S2-BOUNDARY-COUNTEREXAMPLE.md`, `S3-UV-FESHBACH.md`,
`S4-DECIMATION.md`, `S5-N-UNIFORM.md`, `S6-DOBRUSHIN.md`; witnesses in `w7/`;
`C1-SZZ-NSCALING.md`, `C2-BALABAN-UNIFORM.md`, `V3-WAVE-SYNTHESIS.md`.
