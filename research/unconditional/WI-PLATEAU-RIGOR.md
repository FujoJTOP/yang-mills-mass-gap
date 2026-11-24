# WI-PLATEAU-RIGOR — the effective-mass plateau is a theorem GIVEN isolation; its uniformity in (L,a) reduces to `H1` by *both* the spectral and the complex-analytic routes. New: the approach LAW (exponential vs 1/τ) is a rigorous atom/cut discriminator, and the exponential RATE = the isolation gap.

## Setup

Connected correlator `C_O(τ) = Σ_{n≥1} c_n e^{−e_n τ}` (VEV subtracted), `c_n=|⟨n|O|Ω⟩|² ≥ 0`,
`0 < e₁ ≤ e₂ ≤ …`, effective mass `m_O(τ) := −∂_τ log C_O(τ) = (Σ_n c_n e_n e^{−e_nτ})/(Σ_n c_n e^{−e_nτ})`
— the weighted average of the levels.

## [T] P1 — operator-independent limit, given isolation + coupling

If `e₁` is **isolated** (`e₂ > e₁`) and **coupled** (`c₁ > 0`):
`m_O(τ) = e₁ + (c₂/c₁)(e₂−e₁) e^{−(e₂−e₁)τ} + O(e^{−(e₃−e₁)τ})`.
Hence `lim_{τ→∞} m_O(τ) = e₁` for **every** O — the plateau value is the spectral gap, **operator-independent**.
But the finite-τ value is **operator-dependent** through `R_O := (c₂/c₁)(e₂−e₁)`.

## [T] P2 — the spread decays at the isolation gap

For two operators: `|m_O(τ) − m_{O′}(τ)| = |R_O − R_{O′}| e^{−Δτ} + O(…)`, `Δ = e₂−e₁`.
So the WI spread shrinks **at rate Δ = the gap to the second level** (the isolation gap). Reading the spread's
τ-decay gives Δ — but only if Δ is resolvable, i.e. isolation.

## [T] P3 — the approach LAW distinguishes atom from cut

- **Isolated atom**: exponential approach `e₁ + C e^{−Δτ}`.
- **Continuum at the floor** (`ρ(e₁)>0`, smooth): `C(τ) = ∫ρ e^{−eτ}de ∼ ρ(e₁) e^{−e₁τ}/τ` ⇒
  `m(τ) = e₁ + 1/τ + o(1/τ)` — **power-law**, not exponential.
So the *shape* of the approach is a rigorous atom/cut discriminator (sharper than comparing values).

## [T] P4 — complex-analytic route, given a uniform gap

`m_O(τ)` ↑ `e₁` ⟺ `e₁` is the nearest pole of the Stieltjes transform `∫ dν/(s+·)` on the real axis.
**Hadamard/Padé convergence** (and matrix-pencil/ESPRIT) recovers `e₁` from a finite window **iff** the
singularity is an isolated pole and `e₂−e₁ ≥ Δ > 0`. So the complex-analytic mechanism makes the readout
rigorous — *given* a uniform gap.

## [O]/[NG] P5 — uniformity *is* `H1`; no mechanism escapes

Uniform-in-`(L,a)` operator-independence at `β≈2.5` requires `e₂−e₁ ≥ Δ > 0` (a uniform isolation gap) and
`c₁ > 0` — that is **I2 = `H1`**. The spectral route (P1–P3) and the complex route (P4) **both reduce to it**;
neither removes it. The Padé convergence theorem is *conditional on* the uniform gap, which is the wall.

**Verdict.** The plateau's operator-independence is a `[T]` given isolation (τ→∞ limit = gap, O-independent;
spread decays at the isolation rate; approach law discriminates atom/cut). A genuinely new rigorous
discriminator is the **exponential-vs-1/τ approach law**. But its **uniformity in (L,a) is exactly `H1`** —
no escape via complex analysis: Padé convergence *presupposes* the uniform gap. Clay 0%; no QED.

*Witness:* `w7/wi_plateau_rigor.py` (two operators → e₁, spread ~ e^{−2τ}; atom exp vs cut 1/τ).
