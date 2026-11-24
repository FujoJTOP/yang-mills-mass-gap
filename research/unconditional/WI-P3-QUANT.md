# WI-P3-QUANT — the approach law gives a QUANTITATIVE one-sided lower bound on the isolation: `Δ ≥ κ(τ)`, `κ(τ)→Δ`; but it is read-not-force (needs the mass), window-limited, and not a uniform `Δ≥δ` certificate

> WI LINE 1 (P3 approach-law), quantitative. **Result:** the P3 approach law is strictly more than a
> qualitative atom-vs-cut: it yields a **one-sided certified lower bound on the isolation** `Δ=e₂−e₁`.
> In the effective two-level regime the measured approach rate `κ(τ) = −d/dτ log(m(τ)−e₁)` satisfies
> **`κ(τ) ≤ Δ` for every `τ`** and **`κ(τ) → Δ` as `τ→∞`**; a threshold (continuum edge) gives `κ(τ) → 0`.
> So `κ` quantitatively lower-bounds the isolation and discriminates gap from cut. **But** measuring `κ`
> needs `e₁` (the plateau mass), the bound is loose for thin `Δ` (window-limited, `Δτ_max ≳ 1` for ≥90%),
> and multi-level spectra make `κ(τ)` unreliable at small `τ`. So it is a **quantitative diagnostic for
> isolation (`I2`), read-not-force, not a uniform certificate.** `[T]` exact `κ` formulae; `[O]` `Δ≥δ`
> certification is window-limited. Clay 0%; no QED.

---

## §0 · Labels
**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric.

## §1 · Exact approach rate

Connected correlator `C(τ)=Σ cₙe^{−eₙτ}`, effective mass `m(τ)=−∂_τ log C`. Define the approach rate
`κ(τ) = −∂_τ log(m(τ)−e₁)`. Two exact closed forms (witness `w7/wi_p3_quant.py`):
- **atom + gap** (2 level, `r=c₂/c₁`): `m−e₁ = Δ w₂`, `w₂=r e^{−Δτ}/(1+r e^{−Δτ})`, so
  **`κ(τ) = Δ/(1 + r e^{−Δτ}) ∈ (0, Δ)`**, increasing to `Δ`. ⇒ **`Δ ≥ κ(τ)`** (one-sided lower bound).
- **atom + flat threshold**: `C = e^{−e₁τ}(c₁+ρ₀/τ)`, **`κ(τ) = (2c₁τ+ρ₀)/(τ(c₁τ+ρ₀)) → 0`** (`~2/τ`).

Witness (`e₁=0`, `r=0.3`, `ρ₀=0.3`):

| model | κ(5) | κ(10) | κ(20) | κ(30) | κ(50) | κ(100) |
|---|---|---|---|---|---|---|
| gap Δ=0.5 | 0.488 | 0.499 | 0.500 | 0.500 | 0.500 | 0.500 |
| gap Δ=0.2 | 0.180 | 0.192 | 0.199 | 0.200 | 0.200 | 0.200 |
| gap Δ=0.05 | 0.0405 | 0.0423 | 0.0450 | 0.0469 | 0.0488 | 0.0499 |
| threshold | 0.389 | 0.197 | 0.099 | 0.066 | 0.040 | 0.020 |

**Gap ⇒ `κ` rises to `Δ`; threshold ⇒ `κ` falls to 0.** The limit `κ→Δ` reads `Δ`; the bound `Δ≥κ(τ)` holds at every `τ`.

## §2 · Certifiability — window-limited

`κ(τ)/Δ = 1/(1+r e^{−Δτ})`. Tightness table (`r=0.3`):

| Δ \ τ_max | 5 | 20 | 50 | 100 |
|---|---|---|---|---|
| 0.5 | 0.98 | 1.00 | 1.00 | 1.00 |
| 0.2 | 0.90 | 0.99 | 1.00 | 1.00 |
| 0.1 | 0.85 | 0.96 | 1.00 | 1.00 |
| 0.05 | 0.81 | 0.90 | 0.98 | 1.00 |
| 0.02 | 0.79 | 0.83 | 0.90 | 0.96 |

- **A tight (≥90%) lower bound needs `Δ·τ_max ≳ 1`.**
- **Thin isolation** (`Δτ_max≪1`) yields only `≈1/(1+r)≈0.77` of `Δ`: a **conservative** bound — **no false
  positives** (if measured `κ≥δ` then `Δ≥δ`), but **false negatives** for small `Δ` (fails to certify
  `Δ≥δ` even when true). So `Δ≥δ` is certifiable only up to `δ ≈ O(1)/τ_max`.

## §3 · Multi-level caveat

A far heavy level (`e₃=100`, weight `1e-3`) inflates `κ(τ)` **only for `τ≲1/e₃`** (`κ≈89, 86, 79` at
`τ=0.002,0.005,0.01`) and is **fully decayed in the `Δ`-window** (`κ→1.0` by `τ=0.2`). So `Δ≥κ(τ)` is
**regime-conditional**: it holds over the window where the `e₂`-transient dominates (fast transients
decayed), not at arbitrarily small `τ`. The large-`τ` asymptote is the clean read of `Δ`.

## §4 · Verdict

- **`[T]`** the approach law gives a **quantitative one-sided lower bound on the isolation**: `Δ ≥ κ(τ)`,
  `κ(τ)→Δ`, with a quantitative gap-vs-threshold discriminant (`κ→Δ` vs `κ→0`). This is **strictly stronger**
  than the qualitative P3 atom-vs-cut.
- **`[O]/[NG]`** it is **not a uniform `Δ≥δ` certificate**: measuring `κ` needs `e₁` (the mass — read-not-force);
  the bound is **window-limited** (thin `Δ`, `Δτ_max≲1`, under-certified by the `1/(1+r)` factor); and
  multi-level at small `τ` breaks `Δ≥κ`.
- **Deep status:** to use `κ` you must already know `e₁` (the plateau = the mass). So it lower-bounds the
  **isolation given the mass** — a genuine quantitative `I2` read, still **read-not-force** = `H1`-adjacent
  (needs `e₁`, and needs the isolation to exist). It quantifies how **isolated** a detected mass is; it does
  not certify the isolation exists uniformly.

*Files:* `research/unconditional/WI-P3-QUANT.md`, `research/unconditional/w7/wi_p3_quant.py`.
