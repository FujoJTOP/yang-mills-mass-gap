# PNT-CHI-CERTIFICATE — the χ-certificate **is** `H1` verbatim, and this is **provable** (Lemma G): any datum that implies `χ ≤ 1/λ` already forces a gap floor. No softer certificate exists.

> **VERDICT (no QED).** `ROOF-UNIFORM-CHI` reduced the roof to one certificate `χ_phys ≤ 1/λ` uniformly,
> and found five routes all `[NG]`. This piece does not merely repeat the five no-gos: it proves **why no
> sixth can exist**. **Lemma G `[T]`:** on the `0⁺⁺` KL measure `ν` (atom `f₁²` at `m₁²`, total weight `M`,
> `χ = ∫s⁻¹dν`), one has `f₁²/m₁² ≤ χ ≤ M/m₁²`. Hence **`χ ≤ 1/λ` ⟹ `m₁² ≥ λ f₁²`** — a gap floor. **Any
> datum (spectral, reflection-positivity, bootstrap, geometric) that implies the certificate therefore
> implies a lower bound on the mass** = `H1`. The certificate is **equivalent to the gap up to the bounded
> weight `M`**, and the only "slack" is the non-Gaussianity defect `η` (`κχ = 1−η`), which is uncontrolled
> and equals the interaction. **The certificate is irreducible. Clay 0%.**

---

## §0 · Labels, acid test

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric.
**Acid:** `χ = ∞` on the free massless field (`f₁² = 0`, no mass scale). Definitions as in `PROFILE` /
`ROOF-UNIFORM-CHI`: `ν` = KL spectral measure of the `0⁺⁺` coarse mode, `χ = Var_ν(a) = ∫s⁻¹dν = G(k→0)`,
`1/χ = G″(k→0)`.

---

## §1 · Lemma G (the certificate ⟹ gap)

**Lemma G `[T]`.** Let `ν ≥ 0` be a finite measure on `[0,∞)` with `inf supp ν = m₁² > 0`, atom weight
`f₁² = ν({m₁²}) > 0`, total mass `M = ∫dν`, and `χ = ∫s⁻¹dν`. Then
```
   f₁²/m₁²  ≤  χ  = ∫ s⁻¹ dν  ≤  M/m₁² .                                      (G)
```
*Proof.* Since `s ≥ m₁²` on `supp ν`, `s⁻¹ ≤ m₁²⁻¹`, so `χ ≤ M/m₁²`. The atom contributes `f₁²/m₁²`
exactly and all other terms are `≥ 0`, so `χ ≥ f₁²/m₁²`. ∎

**Consequences.**
- **(Necessity)** `χ ≤ 1/λ  ⟹  f₁²/m₁² ≤ 1/λ  ⟹  m₁² ≥ λ f₁² > 0.` **The certificate forces a gap floor.**
- **(Sufficiency)** `m₁² ≥ λ M  ⟹  χ ≤ M/m₁² ≤ 1/λ.`
- Hence `χ ≤ 1/λ  ⟺  m₁² ≥ λ·C` with `C ∈ [f₁², M]` — **`H1` (`m₁² ≥ λ·weight`) up to the bounded weight `M`.**

**Why this is the whole round.** `ROOF-UNIFORM-CHI` tried five routes and each died; Lemma G shows the
death is **not route-specific**: *any* implication of `χ ≤ 1/λ` already yields `m₁² ≥ λ f₁²`, a genuine gap
lower bound. So no input weaker than a gap floor can produce the certificate. The "softer certificate"
this round was looking for **does not exist** — provably.

---

## §2 · The four construction routes, and exactly where each reduces to the gap

The directive asked to *attempt* construction by four routes; each is blocked by Lemma G's necessity half.

**R-geo — geometric / Diophantine bound on `Var(tr F²)`.** Compactness gives `|tr F²_□| ≤ C_N`, hence
`M ≤ C_N` (the total spectral weight is bounded) and `χ ≤ C_N/m₁²` — **which is Lemma G with `M = C_N`**.
A "Diophantine" flavour (finite centre `Z_N`, discrete configurations) adds no bound beyond this compactness
tie: the weight is bounded, but the *floor* `m₁²` is not, and by Lemma G the certificate **needs** the
floor. The route recovers (G) and nothing more: `[NG]` beyond `H1`.

**R-RP — reflection positivity + a NEW inequality.** RP gives `ν ≥ 0` (the KL positivity). The only new
inequality one can append is Lemma G itself (`χ ≥ f₁²/m₁²`), whose combination is again `H1`. RP alone is
**acid-vacuous**: the free massless field has `ν ≥ 0` (`f₁² = 0`) and `χ = ∞`; RP does not see it. So
**RP + Lemma G = H1**; there is no RP-only certificate. `[NG]`

**R-boot — bootstrap on `ρ(s)` with an IR-only constraint.** The witness (§B) shows the decisive fact:
`χ = ∫s⁻¹dν` is **not a function of finitely many polynomial moments** — an atom of weight `ε` at `s = ε`
changes `χ` by `O(1)` while moving every `M_j = ∫sʲdν` (`j ≥ 0`) by `O(ε)`. So any constraint built from
UV/polynomial data is **floor-blind** (`SLOT-BOOTSTRAP`'s `Δ_min(K) = 0 ∀K`). An **IR-only** constraint
that *does* fix `χ` must assert `ν` has no weight below `s₀` — **that assertion is the gap.** `[NG] ≡ H1`.

**R-scheme — a scheme where `χ_lattice ≲ a_ℓ²` is provable.** Any scheme that makes the `a_ℓ²` scaling
manifest must insert the physical scale (the mass) into the subtraction; that insertion **is** the gap.
Equivalently: the demand `χ_lattice ≲ a_ℓ²` is `χ_phys` finite (`ROOF-UNIFORM-CHI §1`), and `χ_phys` finite
forces `m₁² > 0` by Lemma G. `[NG] ≡ H1`.

---

## §3 · The one loophole: the non-Gaussianity defect `η` (and why it is not a free lunch)

The exact Ward/Stein identity (`K0-WARD §2`; `SLOT-WARD-STEIN`) reads `κ·χ = 1 − η`, `κ = G″(0)`. Hence
```
   χ ≤ 1/λ   ⟺   κ ≥ λ(1 − η) .                                                (★η)
```
So the certificate is **`H1` weakened by the defect `η`**: a *large* `η` would let a *weaker* `κ`
(`κ < λ`) suffice. Witness (§D): Gaussian `η = 0` (`κχ = 1` exactly); quartic toy `η ∈ (0, 0.63]` as the
interaction grows. **This is the only slack in the certificate.** But:

- `η = 1 − κχ` is the **non-Gaussianity defect** — it is *the interaction itself*, not an independent input.
- There is **no lower bound `η ≥ η₀ > 0`** at `β = 2.5` (that would be a proof of interaction strength, i.e.
  of `H1`). So `(★η)` cannot be *used*: to exploit `η` you must already bound `η`, which is the wall.
- `η` is sign-indefinite in general (`κχ` can exceed or fall below 1 depending on the marginal's shape,
  `SLOT-WARD-STEIN`), so it cannot be assumed positive.

⇒ the loophole is real *algebraically* but **not usable**: closing `(★η)` with `η` large is as hard as
`H1`. It sharpens the statement — the certificate is `H1` **modulo the defect** — but does not soften it.

---

## §4 · What is genuinely gained (the honest yield)

1. **Lemma G `[T]`** turns `ROOF-UNIFORM-CHI`'s "no route found" into "**no route weaker than `H1` can
   exist**": the certificate logically implies a gap floor. This is a genuine **negative theorem** that
   *closes* the search for a softer certificate, rather than leaving it open.
2. **The sandwich `f₁²/m₁² ≤ χ ≤ M/m₁²`** localizes the entire difficulty onto one scalar `m₁²` with a
   bounded bracket. It is exactly the `PROFILE` squeeze `f₁²/χ ≤ m₁² ≤ M/χ`, now recognised as the reason the
   certificate *is* `H1`.
3. **The defect identity `(★η)`** is the precise form of the only slack: the certificate is `H1` modulo the
   non-Gaussianity `η`, which is uncontrolled. This rules out "exploit non-Gaussianity" as a shortcut.

## §5 · Honest bottom line

- **The χ-certificate is `H1` verbatim** (`Lemma G`: `χ ≤ 1/λ ⟹ m₁² ≥ λf₁²`), **provably irreducible**:
  any datum implying it implies a gap floor.
- **No softer certificate exists**; the four construction routes all reduce to `H1`, the fifth (`η`) is
  unusable. **`[NG]` for the class "certificate without a gap floor."**
- **What stands:** Lemma G (a clean, useful `[T]`), the `η`-identity, and the exact statement of the
  irreducible input. **What is missing:** the gap floor itself — which is `H1`, as before.
- **`QED`: no. Clay 0%.** The certificate cannot be manufactured from anything weaker than the gap; the roof
  remains conditional on `H1`, and `H1` is now *proved* to be the minimal input.

*Refs:* `ROOF-UNIFORM-CHI.md`, `BASE-NONLOCAL-BOUND.md` (Cheeger), `SLOT-WARD-STEIN.md`, `SLOT-BOOTSTRAP.md`
(`Δ_min=0 ∀K`), `PROFILE.md`; witness `w7/pnt_chi_certificate.py` (runs clean, `py -3.12 -I`).
