# WI-P3-RIGOR — the approach law as an isolation detector: rigorous exp-vs-power-law in the infinite-volume limit, but at fixed finite L every discrete spectrum is exponential (content = the L-scaling of the approach rate = isolation I2 = H1)

> **VERDICT (no QED, Clay 0%).** WI-PLATEAU-RIGOR's P3 (the effective-mass **approach law** separates an
> isolated atom — exponential `e₁+Ce^{−Δτ}` — from a continuum edge — power law `e₁+c/τ`) is made precise and
> lattice-testable. **In the infinite-volume (`L→∞`) limit it is a rigorous theorem `[T]`** (the atom gives an
> exponential approach with rate `Δ=e₂−e₁`; a continuum edge of density `ρ∼(s−E)^p` gives `m(τ)−E ∼ (p+1)/τ`,
> a power law). **But at a single finite `L` every transfer matrix is compact, so the spectrum is discrete and
> the lightest level is *always* isolated by a finite gap `δ_L`** — the approach is exponential for **both** an
> atom **and** a cut (rate `δ_L`). So **a single finite L cannot discriminate.** The discriminating content is
> the **`L`-scaling of `δ_L`**: an isolated atom ⇒ `δ_L→Δ>0`; a two-particle cut ⇒ `δ_L→0` (levels crowd). A
> uniform `δ_L≥Δ>0` is **isolation = I2 = H1**. The approach law is a genuine `[C]` diagnostic, not a
> certificate.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric.
Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · Setup and exact decomposition

Connected 0⁺⁺ correlator on a finite lattice, `C_O(τ)=Σ_n c_n e^{−e_nτ}` (`c_n=|⟨n|O|Ω⟩|²≥0`, `e_n` discrete).
Effective mass `m_O(τ) = −∂_τ log C_O(τ) = Σ_n e_n c_n e^{−e_nτ} / Σ_n c_n e^{−e_nτ}` (the `τ`-weighted mean of levels). Split off the lightest level `e₁`:

```
C = c₁ e^{−e₁τ} (1 + ε(τ)),   ε(τ) = (1/c₁) Σ_{n≥2} c_n e^{−(e_n−e₁)τ},
m − e₁ = −ε′/(1+ε).
```

---

## §2 · The infinite-volume theorem `[T]`

- **Isolated atom** (next level at `e₂`, gap `Δ=e₂−e₁`): `ε=(c₂/c₁)e^{−Δτ}`, so
  `m(τ)−e₁ = Δ·(c₂/c₁)e^{−Δτ} + O(e^{−2Δτ})` — **exponential**, rate `Δ`.
- **Continuum edge** (density `ρ(s)=ρ₀(s−E)^p+…` for `s≥E=m₁²`, no floor atom): by Watson's lemma
  `∫_E (s−E)ρ(s)e^{−(s−E)τ}ds = ρ₀ Γ(p+2)/τ^{p+2}`, and `∫ρ(s)e^{−(s−E)τ}ds = ρ₀ Γ(p+1)/τ^{p+1}`, so
  `m(τ)−E = (p+1)/τ + O(1/τ²)` — **power law**.

**Discriminator:** the local log-derivative `T(τ) := −d/dτ log|m(τ)−ê₁|`:
- atom: `T(τ) → Δ` (a **constant**, `L`-independent);
- continuum: `T(τ) = 1/τ → 0` (decays like `1/τ`).

Equivalently: `log|m−ê₁|` is **linear in `τ`** (atom) vs **linear in `log τ`** (cut). (Witness §1/§2: atom
`T=0.88,0.98,0.9997`; continuum `m−E=2/τ` and `T=1/τ` to 4 digits at `τ=15,30`.) `[T]`

---

## §3 · The finite-`L` subtlety `[T]`/`[V-num]` (the honest core)

On a finite lattice the transfer matrix is **compact** — the spectrum is discrete and the lightest level has a
**finite gap `δ_L = e₂^{(L)} − e₁^{(L)}`** to the next level. For `τ ≫ 1/δ_L` the sum is dominated by the two
lowest levels, so
```
m(τ) − e₁  =  δ_L (c₂/c₁) e^{−δ_L τ} (1+…)   —  EXPONENTIAL, rate δ_L,
```
**for an atom AND for a cut alike.** A single finite `L` therefore **cannot** distinguish them: both look
"exponential" (any discrete spectrum does). (Witness §3: `N` levels of spacing `δ` give `T ≈ δ` — `0.105/0.053/
0.021/0.0105/0.0053/0.0026` for `δ=0.1…0.0025` — and `δ→0` as `N` grows.) **The content is the `L`-scaling
of `δ_L`:** atom ⇒ `δ_L→Δ>0` (a genuinely separated state); cut ⇒ `δ_L→0` (level spacing `∝ L^{−d}·` density).
This is the same "finite `L` everything is discrete" caveat as `WI-THEOREM`/`CONSTRUCT-IR-ANALYTIC-NUM`.

---

## §4 · Estimator, test statistic, protocol

- **Estimator (model-independent):** `ê₁` = late-`τ` effective mass; then `T(τ)=−d/dτ log|m(τ)−ê₁|` on a
  `τ`-window. **Test:** `T(τ)` constant (flat in `τ`) ⇒ atom; `T(τ)∝1/τ` (log-log slope −1) ⇒ cut.
- **Estimator (finite-`L`, practical):** fit `C_O(τ)` to 2–4 exponentials; extract `δ_L = e₂−e₁` and check its
  **operator-independence** (atom: `δ_L` robust under reweighting; cut: `δ_L` tracks overlaps). Repeat for
  `L=4,6,8` (SU(2)) / `16–24` (SU(3)) and **extrapolate `δ_L(L)`**.
- **Protocol:** 0⁺⁺ operators (`tr F²`, plaquette, APE/`HYP`-smeared, clover) at `β≈2.5` (SU(2)) or `β≈6/3.75`
  (SU(3)); `τ` up to the SNR limit; **`τ ≫ 1/δ_L`** so the two-level regime is reached. The SU(3) protocol of
  `WI-SU3` applies verbatim (`τ`-fit + spread over smearing).
- **What it certifies:** an `L`-extrapolation with `δ_L ≥ Δ > 0` **bounded below** ⇒ isolation (`I2`);
  `δ_L→0` ⇒ no isolation (a cut). Either is a **`[C]` diagnostic**; a single `L` certifies nothing.

---

## §5 · Honest bottom line

- **Rigorous `[T]`:** the infinite-volume exp-vs-power-law separation (Watson's lemma) — `T(τ)→Δ` (atom) vs
  `1/τ→0` (cut). Genuine, acid-clean (a free/massless field: `m=0`, no approach).
- **Honest `[T]/[V-num]`:** at fixed finite `L` every spectrum is discrete ⇒ exponential for **both**; the
  functional form does **not** discriminate — **only the `L`-scaling of `δ_L` does.**
- **The wall `[O]`:** certifying `δ_L ≥ Δ>0` **uniformly in `L`** (and `a→0`) is **isolation `I2` = `H1`** —
  the same wall, now in the approach-law coordinate. The `δ_L(L)` extrapolation is gated by statistics
  (`WI-STATS`' `Δ_cut`) and is `[C]`, never a proof.
- **This is a real diagnostic refinement of P3, not a gap certificate.** It makes P3 lattice-testable *and*
  states honestly that a single lattice cannot decide it — the decision is the `L→∞` isolation = `H1`.

*Refs:* `WI-PLATEAU-RIGOR.md` (P3), `WI-SU3.md`, `WI-STATS.md`, `WI-THEOREM.md`, `WI-EMBEDDED.md`,
`CONSTRUCT-IR-ANALYTIC-NUM.md`, `ONLY-EXIT-SYNTHESIS.md`; witness `w7/wi_p3_rigor.py`.
