# WI-GENERAL — the "operator-reweighting invariance = atom signature" principle, generalized

> **VERDICT (no QED).** The reweighting-invariance (WI) idea of `CONSTRUCT-IR-ANALYTIC-NUM` generalizes
> cleanly, but with one important correction. **Reweighting = changing the residues at fixed poles**
> (the spectrum is operator-independent; the weights are not). A functional of the spectral data is a
> **candidate atom-certifier iff it is (a) reweighting-invariant and (b) acid-clean** (distinguishes the
> atomless `:φ²:` edge). **[T]** the *nearest complex-`τ` zero* `π/|Im τ|` is WI **only** for a pure
> two-term atom; in the presence of a continuum it grabs the *continuum's* zero (`π/2 < π/1`) and is NOT
> atom-valued — a correction to `PNT-IR-ANALYTIC`. **[T]/[V-num]** the *robust* WI atom signature is the
> **effective-mass plateau** `m_eff(τ)→E_min` (equivalently the **Padé-pole position**), which is the
> standard lattice mass extraction and whose WI property *is* operator-independence of the mass. The
> existence of a stable plateau across `(L,a)` at `β≈2.5` is isolation `H*` = `H1`. **Clay 0%.**

---

## §0 · Setup and the reweighting equivalence

A finite system gives `Z(τ) = w₀ + Σᵢ aᵢ e^{−Eᵢ τ}`, where `Eᵢ` = spectrum of `H` (**operator-independent** = the **support** of the KL measure `ν`), `aᵢ = |⟨i|O|0⟩|²`, `w₀ = |⟨0|O|0⟩|²` (**operator-dependent** = the **weights**).

**Definition (reweighting).** Change the residues `{aᵢ}, w₀` at **fixed** poles `{Eᵢ}` (mix `O → O+tO'`). A functional `𝒜(Z)` is **reweighting-invariant (WI)** if `𝒜` depends only on `{Eᵢ}` (supp `ν`), not `{aᵢ, w₀}`.

**Definition (atom).** `ν` has an atom at the edge `m₁=E₀` iff `E₀` is an **isolated** point of the support (uniform gap `δ>0` to the rest). Continuum-at-edge (dense poles at `E₀`, no gap) = no atom.

---

## §1 · The functional landscape

| functional | WI iff atom? | lattice-measurable? | verdict |
|---|---|---|---|
| **Stieltjes-transform pole structure** `Π(z)=∫dν/(z−s)` | **YES** — poles = atoms = support; branch cut = continuum | reading a pole from finite moments = ill-posed analytic continuation (`SLOT-BOOTSTRAP`, `PNT-OPE-RHO`) | **ideal certifier, un-extractable** `[O]` |
| **Padé-pole positions / effective-mass plateau** | **YES** — pole position = spectrum (operator-indep. mass) | **YES** — standard multi-exponential / Prony / GEVP / variational, from **real-`τ`** data | **most robust AND measurable** `[C]` |
| **Zero-ladder ratios** `|Im τ_k|=π/E₀·(2k+1)` | **YES for isolated atom** (arithmetic, odd multiples); cut = irregular | needs complex-`τ` (ill-posed) or Padé-zero count | clean but **fragile** (continuum zeros interleave) `[V-num]` |
| **Zero-density profile** (count / window) | atom sparse+periodic, cut dense+crowding | complex-`τ` | fragile, hard to measure |
| **Christoffel-function zeros** (from moments) | zeros → supp `ν` as `N→∞`; finite-`N` not invariant, **stability ⟺ atom** | moment stability only | moderate; needs high-precision moments |

---

## §2 · The general theorem (schematic)

**Theorem (WI atom-certifier principle).** Fix the pole set `P={Eᵢ}` (supp `ν`, operator-independent). A functional `𝒜` of the spectral sum is an **atom-certifier** — returns the atom mass and is nonzero/nontrivial on atomic vs atomless measures — **iff**:
- **(WI)** `𝒜` is invariant under all residue changes at fixed `P` (depends only on `P`), **and**
- **(AC)** `𝒜` is **acid-clean**: it takes a distinguished value on atomless measures (the `:φ²:` continuum edge, `f₁²=0`).

Both are necessary: a constant is WI but not AC; the raw nearest-zero is AC-ish but not WI.

**Proof sketch (the structural content).**
- *WI ⟹ depends only on `P`.* In `Z = Σ aᵢ e^{−Eᵢτ}`, the `Eᵢ` are fixed by `H`; any `{aᵢ}` is reachable by operator mixing. A WI functional factors through `P`.
- *`P` isolated edge ⇒ atom-valued (T1, T2).* Pure two-term `Z=w₀+a e^{−mτ}`: zero at `|Im τ|=π/m`, independent of `w₀,a` — **𝒜=π/|Im τ|=m exactly**. With higher states at gap `δ`: the atom's zero is corrected by `O(e^{−δ/m})`. The **effective-mass plateau** `m_eff(τ)=[ΣaᵢEᵢe^{−Eᵢτ}]/[Σaᵢe^{−Eᵢτ}]→E_min` at large `τ`, **WI** (ratio → `E_min` for any residues).
- *`P` dense at edge ⇒ not WI (T3).* Many terms contribute to the zeros / to `m_eff(τ)`; both depend on `{aᵢ}` — reweighting moves them (spread `0.2`).
- *AC excludes trivial invariants.* The atomless `:φ²:` edge (`f₁²=0`) must give a different value; a constant fails.

---

## §3 · Witness (`w7/wi_general_witness.py`, numpy, runs clean)

- **(T1)** single atom `E₀=1`: `𝒜=π/|Im τ|=1.00000000` for `(w₀,a)∈{(1,.5),(1,.2),(1,2),(.3,5),(7,.05)}` — **exact, weight-independent**.
- **(T2)** atom `E₀=1` + continuum from `2.0` (gap 1): `m_eff(τ)=2.43,2.13,1.65,1.10,1.001` at `τ=0.3,1,2,4,8` — **plateaus at `E₀=1`**, reweighting-invariant (spread `3×10⁻³`); meanwhile the **nearest-zero** `𝒜` is `≈2.2` (grabs the continuum zero `≈π/2`) — **NOT atom-valued** (this corrects `PNT-IR-ANALYTIC`'s "raw 𝒜 is acid-clean").
- **(T3)** continuum at edge from `1.0`: nearest `𝒜` moves under reweighting (`spread 0.176`); no plateau.
- **(T4)** atom zero ladder `|Im|/π = 1,1,3,3,5,5` (odd multiples); cut `0.67,0.67,0.78…` (irregular).
- **(T5)** `m_eff(τ)` atom flat `1.0000`; cut drifts `1.68→1.52→1.33→1.18`.

---

## §4 · Ranking and the honest caveat

**Most robust + lattice-measurable WI atom signature: the effective-mass plateau / Padé-pole position.** It is *the* standard lattice glueball-mass method, and its WI property *is* the physical statement "the mass is operator-independent." Measurable from real-`τ` data by varying the operator and looking for a `τ`-plateau.

**Why this is still not a free win (the bite).** The plateau exists **iff** the edge is isolated (gap above `E_min`); proving a **uniform-in-`(L,a)`** plateau at `β≈2.5` *is* isolation `H*` = `H1`. So the WI criterion is a **diagnostic/consistency check** (reweight the operator, test invariance) that a candidate atom passes, not a proof that it is one — the uniformity is the wall.

**Net.**
- [T] WI ⟺ support-only dependence; isolated atom ⟹ WI + atom-valued; dense edge ⟹ not WI.
- [T] The **nearest-zero** `π/|Im τ|` is WI **only** for the pure 2-term atom; **not** in the presence of a continuum (correction to `PNT-IR-ANALYTIC`).
- [V-num] The **effective-mass plateau** is the robust, real-`τ`-measurable WI signature; the **zero ladder** is a clean complex-`τ` signature.
- [O] Uniform isolation (the plateau's existence across `L,a` at `β≈2.5`) is `H1`.

*Labels:* [T] structural (WI ⟺ support-only; pure-atom exactness; dense-edge non-invariance); [V-num] witness; [O] uniform isolation = H1. No QED; Clay 0%.

*Witness:* `research/unconditional/w7/wi_general_witness.py`.
