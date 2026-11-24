# KAPPA-BREATHING-MATCH — the rigid-string model's second, independent prediction (0⁺⁺*) CONTRADICTS lattice: κ is NOT pinned by two-state consistency; the model is falsified at the excited level

> **VERDICT (no QED).** The rigid-string model (MASTERFIELD/FLUXTUBE-BREATHING) makes two predictions in terms of ONE κ: the lightest `0⁺⁺` and the breathing `0⁺⁺*`. Computing κ from each gives **κ_light = 0.193** vs **κ_exc = 0.111** (factor 1.73 apart). The model's **internal excitation ratio** `m(0⁺⁺*)/m(0⁺⁺) = 1/(π√2 κ) = 1.168` **contradicts** the lattice ratio **1.538** (Morningstar–Peardon) by a factor 1.32 — and this ratio is **normalisation-independent**, so no choice of `m(0⁺⁺)/√σ` rescues it. **The rigid string is a one-state fit, not a two-state prediction; its breathing mode is NOT the physical `0⁺⁺*`. κ remains underived = `H1`, and now additionally the model's excited-state relation is falsified.**

---

## §0 · Setup — two relations, one κ

Rigid string `E(R)=2πσR+πκ/R`, minimum `R*=√(κ/2σ)`:
- **lightest:** `m₀ = 2π√(2σκ)` ⇒ `m₀/√σ = 2π√(2κ)` ⇒ **`κ = (m₀/√σ)²/(8π²)`**.
- **breathing (radial osc):** `m_br = √(E''/M) = 2√(σ/κ)` ⇒ `m_br/√σ = 2/√κ` ⇒ **`κ = 4/(m_br/√σ)²`**.
- **model-internal ratio:** `m_br/m₀ = [2/√κ]/[2π√(2κ)] = 1/(π√2·κ)` — a function of κ alone, **normalisation-free**.

## §1 · The numbers (witness `w7/kappa_breathing_match.py`)

| quantity | value |
|---|---|
| `m(0⁺⁺)/√σ` (input) | 3.9 |
| **κ from lightest** = (3.9)²/(8π²) | **0.19264** |
| model-predicted `m(0⁺⁺*)/√σ` = 2/√κ | **4.5568** |
| lattice `m(0⁺⁺*)/m(0⁺⁺)` (Morningstar–Peardon) = 5.46/3.55 | **1.538** |
| lattice `m(0⁺⁺*)/√σ` = 3.9·1.538 | **5.998** |
| **κ from excited** = 4/(5.998)² | **0.11117** |
| model-internal ratio 1/(π√2·κ) | **1.168** vs lattice **1.538** (factor **1.316**) |
| **κ_light / κ_exc** | **1.733** |

**Verdict: DISAGREE.** Two κ from two states differ by 73%; the model's own predicted excitation ratio (1.168) missed the lattice (1.538) by 32%.

## §2 · Normalisation robustness (why this is not a `√σ`-convention artefact)

The model's excitation ratio `1/(π√2 κ)` and the lattice ratio `1.538` are **both normalisation-independent** — they do not depend on the value of `m(0⁺⁺)/√σ` or `√σ`. So the contradiction stands regardless. (For completeness: with Morningstar's own `r₀√σ≈1.16`, `m(0⁺⁺)/√σ≈3.06`, the excited lands at ≈4.71 ≈ the model's 4.56 — *close* — but then κ_light=0.119 vs κ_exc=0.180 **still disagree**, factor 1.52. No convention makes both states consistent.)

## §3 · What it means

- **κ is NOT empirically pinned by two-state consistency.** The one-state fit κ=0.193 does not survive the second state. So the earlier hope "two independent states pin κ" is **refuted**.
- **The breathing-mode identification is wrong**: the radial breathing excitation of the rigid closed string (at 1.17× the ground state) is **not** the physical first-excited scalar glueball (at 1.54×). The rigid-string model, as built, **under-predicts the excitation gap**.
- **κ remains underived = `H1`** (CONSTRUCT-KAPPA-DERIVE: `κ = (m/√σ)²/(8π²)`, so deriving κ ≡ deriving the gap). The two-state route does not deliver it.
- **A genuine, falsifiable negative:** the rigid-string model is now **falsified at the excited level**, not merely "missing an input." This sharpens the string route: a correct flux-tube model must reproduce the **ratio** 1.538, which the simple `2πσR+πκ/R` ansatz does not.

*Citation:* Morningstar & Peardon, *Phys. Rev. D* **60**, 034509 (1999), arXiv:hep-lat/9901004 (`m(0⁺⁺)=3.55(7)`, `m(0⁺⁺*)=5.46(22)` in `r₀⁻¹`). ⚠ the specific table values should be re-verified against the paper; the **ratio** 1.538 is the scale-independent content and is robust to the ± errors.

*Labels:* `[V-num]` the two κ and the ratio comparison; `[T]` the model-internal ratio formula; `[NG]` for the rigid-string-breathing identification; `[O]` κ underived = `H1`. **No QED; Clay 0%.**
