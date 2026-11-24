# WI-P3-CONTINUUM — the P3 approach-law SURVIVES a background continuum (atom+continuum ⇒ exponential, pure continuum ⇒ power-law); it is the robust, acid-clean isolator the complex-zero `𝒜` is not — but it discriminates only in the `L→∞` limit

> **VERDICT.** The approach-law discriminator (`WI-PLATEAU-RIGOR` P3: `m_O(τ)` approaches `e₁` **exponentially** for an isolated atom, **power-law** for a continuum edge) **does survive a background continuum**. Concretely, with an atom at `e₁` **below** a continuum threshold `E_thr` (`Δ=E_thr−e₁>0`, coupling `c₁>0`), `m_O(τ)−e₁` decays **exponentially** at rate `Δ`; a **pure continuum** (no atom below) decays **power-law** `~1/τ`. This is exactly where the WI complex-zero `𝒜` fails (`WI-EMBEDDED`: any continuum, even 0.01%, destroys it). So the approach-law **is** the robust, acid-clean isolator the complex-zero is not. **But** it is an **infinite-volume** discriminator: on a finite lattice the spectrum is discrete, both are exponential, and the power-law is the `L→∞`/`a→0` limit — so using it to certify the YM atom still needs the `L→∞` **survival** of the exponential = isolation = `I2` = `H1`.

---

## §1 · The calculation (Watson's lemma)

`C_O(τ)=c₁e^{−e₁τ}+∫_{E_thr}^∞ρ(e)e^{−eτ}de`. For smooth `ρ`, `∫_{E_thr}^∞ρ(e)e^{−(e−e₁)τ}de ≈ ρ(E_thr)e^{−Δτ}/τ` (`Δ=E_thr−e₁`), so

`m_O(τ)−e₁ = (ρ(E_thr)/c₁)(Δ/τ + 1/τ²)e^{−Δτ} + O(e^{−2Δτ})`  → **exponential, rate `Δ`.** `[T]`

Pure continuum (`c₁=0`): `C≈ρe^{−E_thrτ}/τ` ⇒ `m_O−E_thr = −1/τ` → **power law `1/τ`.** `[T]`

**Exact condition:** `Δ=E_thr−e₁>0` (the atom is strictly below the cut = isolation `I2`) and `c₁>0` (atom coupled). **Robust to the operator:** any `O` with `c₁>0` gives the exponential (no fine-tuned overlap needed — unlike `WI-CONTINUUM-SUBTRACT`, which needed `⟨0|O−O′|cont⟩=0`).

## §2 · Witness (`w7/wi_p3_continuum.py`, exact quadrature)

| case | `m_O(τ=30)` | `d log|m−floor|/dτ` | reading |
|---|---|---|---|
| **atom + continuum** | `1.000000` (dev `2.0e−15`) | **`−1.066 ≈ −Δ`** | **EXPONENTIAL** (isolated atom) |
| **pure continuum** | `2.033` (dev `3.3e−2`) | `−0.062 ≈ 0` (`dev·τ≈const`) | **POWER LAW** |
| **free massless** (acid) | `0.033` (floor 0, dev·τ≈1) | `−0.062` | **POWER LAW ⇒ says NO atom** |

**Separation `−1.066` vs `−0.062`** — clean. And **acid-clean**: free massless / `:φ²:` (continuum) give power-law = "no atom", correct — unlike the complex-zero `𝒜` (acid-dirty, fires on both).

## §3 · Contrast with the complex-zero `𝒜`

- **Complex-zero `𝒜=π/|Im τ₀|`**: destroyed by **any** continuum, even 0.01% (`WI-EMBEDDED`, `WI-ADMISSIBLE`); acid-dirty. **Not** robust to the physical atom+continuum.
- **Approach-law (P3)**: **survives** a background continuum (atom below the cut dominates `C_O` exponentially); acid-clean. **This is the robust isolator.**

## §4 · The honest caveat — infinite-volume only

The power-law is a property of a **continuous** spectrum (`L→∞`, `a→0`). On a **finite lattice** the spectrum is discrete ⇒ `C_O` is a finite exponential sum ⇒ `m_O→e₁` **exponentially** for **any** finite `L`, whether or not the `L→∞` limit is an atom. So at fixed `L`, both atom and continuum give exponential; the discriminator is a **limit statement**. Using it to certify the physical atom needs the **`L→∞` survival** of the exponential approach = **isolation `I2` = `H1`**.

## §5 · Verdict

- `[T]` the approach-law distinguishes atom+continuum (exponential, rate `Δ=E_thr−e₁`) from pure continuum (power-law `1/τ`), **with the continuum present**; acid-clean on free/`:φ²:`.
- `[T]` exact condition: `Δ>0` (isolation) ∧ `c₁>0` (coupling); robust to operator choice.
- `[V-num]` witness separation: atom slope `−1.066` vs continuum `−0.062`.
- `[O]` it is an **infinite-volume** discriminator; its finite-`L` realization reads the finite gap, and the `L→∞` survival = isolation = `I2` = `H1`.
- **Net:** the approach-law is a **genuinely robust, acid-clean isolator** (better than the complex-zero `𝒜`) — but, like every good WI route, its decisive use is gated by the isolation `L→∞` limit = `H1`. **A real diagnostic; not a gap certifier.** Clay 0%; no QED.
